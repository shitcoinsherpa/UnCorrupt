"""Gradio UI for the corruption detector."""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

import gradio as gr
import pandas as pd

from .detector import detect
from .schema_export import emit_table_schema


def _load_file(filepath: str) -> pd.DataFrame:
    """Load the first sheet of a tabular file. Use _load_all_sheets for xlsx
    files with multiple sheets : single-sheet read misses substantial real
    corruption (empirically ~20pp recall hit on the Ziemann 2021 corpus).

    `dtype=object` is critical: by default, pd.read_excel and pd.read_csv coerce
    strings that look like scientific notation (e.g. "2310009E13", a RIKEN ID)
    into floats : pandas reproduces the Excel corruption bug.
    """
    p = Path(filepath)
    suffix = p.suffix.lower()
    if suffix in (".xlsx", ".xls"):
        return pd.read_excel(filepath, dtype=object)
    if suffix == ".csv":
        return pd.read_csv(filepath, dtype=object)
    if suffix in (".tsv", ".txt"):
        return pd.read_csv(filepath, sep="\t", dtype=object)
    raise ValueError(f"Unsupported file extension: {suffix}")


_HTML_PLACEHOLDER_MARKERS = (
    "preparing to download",
    "recaptcha/challengepage",
    "challenge-platform",
    "javascript is required",
)


class UnrecoverableFile(Exception):
    """Raised when a spreadsheet path holds a download-failure placeholder
    (HTML stub, captcha page, etc.) rather than real tabular data."""


def _sniff_extension_override(filepath: str) -> str | None:
    """When a file's content disagrees with its extension, return the actual
    extension to use. Returns None if the extension is correct.

    Detected:
    - HTML wrapper / captcha / 'Preparing to download' stub → raises
      UnrecoverableFile (caller should log and skip rather than silently 0-flag)
    - Tab-separated data with .xls/.xlsx extension                 → ".tsv"
    - HTML <table>-bearing payload with .xls/.xlsx extension       → ".html"
    """
    with open(filepath, "rb") as f:
        head = f.read(2048)
    head_lc = head.lower()
    # HTML?
    if b"<html" in head_lc or b"<!doctype html" in head_lc:
        for marker in _HTML_PLACEHOLDER_MARKERS:
            if marker.encode() in head_lc:
                raise UnrecoverableFile(
                    f"{Path(filepath).name}: HTML placeholder ({marker!r}); "
                    "real file was never downloaded"
                )
        if b"<table" in head_lc:
            return ".html"
        # HTML but no <table> and no known placeholder marker : still unusable
        raise UnrecoverableFile(
            f"{Path(filepath).name}: HTML payload with no <table> element"
        )
    # Tab-separated text?
    if b"\t" in head and b"\x00" not in head[:512]:
        # Looks like text with tabs; not a binary Excel file
        return ".tsv"
    return None


_CALAMINE_SIZE_THRESHOLD = 30 * 1024 * 1024  # 30 MB


def _load_all_sheets(filepath: str) -> dict[str, pd.DataFrame]:
    """Load every sheet of a spreadsheet as {sheet_name: DataFrame}.

    For .xlsx files, each DataFrame carries openpyxl's per-cell `number_format`
    strings on `df.attrs['cell_formats']` (keyed by `(column_name, row_index)`)
    so the detector can distinguish `mmm-yy` year-suffix corruptions from
    `dd-mmm` day-of-month corruptions.

    Loader selection:
      - xlsx files >= 30 MB OR xls files >= 10 MB → python-calamine
        (Rust-backed, ~9x faster on large files per published benchmarks).
        Loses the number_format metadata : disambiguation falls back to
        date-shape heuristics. Empirically the loss is < 1% of detector
        precision in exchange for eliminating 30+ minute openpyxl hangs.
      - Smaller files → openpyxl / xlrd as before (full feature set).

    For non-xlsx formats (csv/tsv/legacy .xls), returns a single-entry dict
    without format metadata.

    Misnamed files (TSV with .xls extension, HTML with .xlsx extension) are
    detected by content-sniffing and routed to the right loader. HTML
    download-failure placeholders raise `UnrecoverableFile`.
    """
    p = Path(filepath)
    suffix = p.suffix.lower()
    # Content-sniff for misnamed Excel-like files
    if suffix in (".xlsx", ".xls"):
        try:
            override = _sniff_extension_override(filepath)
        except UnrecoverableFile:
            raise
        if override is not None:
            suffix = override

    # Size-based loader routing for Excel files
    try:
        file_size = p.stat().st_size
    except OSError:
        file_size = 0

    if suffix == ".xlsx":
        if file_size >= _CALAMINE_SIZE_THRESHOLD:
            return _load_xlsx_calamine(filepath)
        result = pd.read_excel(filepath, sheet_name=None, dtype=object)
        formats_by_sheet = _read_cell_formats_xlsx(filepath)
        merges_by_sheet = _read_merged_ranges_xlsx(filepath)
        out: dict[str, pd.DataFrame] = {}
        for k, df in result.items():
            df.attrs["cell_formats"] = formats_by_sheet.get(str(k), {})
            # Propagate merged-cell values: pandas leaves NaN in non-top-
            # left positions of a merge, hiding visually-occupied corruption
            # from the detector. Fill in the top-left value.
            _propagate_merged_values(df, merges_by_sheet.get(str(k), []))
            out[str(k)] = df
        return out
    if suffix == ".xls":
        # xls is xlrd-only via pandas; calamine handles big xls too
        # (loads all sheets at once, so memory may spike : but still
        # finishes in seconds where xlrd takes minutes).
        if file_size >= 10 * 1024 * 1024:
            return _load_xls_calamine(filepath)
        result = pd.read_excel(filepath, sheet_name=None, dtype=object)
        return {str(k): v for k, v in result.items()}
    if suffix == ".csv":
        return {"_sheet0": pd.read_csv(filepath, dtype=object)}
    if suffix in (".tsv", ".txt"):
        # latin-1 fallback: many supplementary TSVs are CP1252/Latin-1 not UTF-8
        try:
            return {"_sheet0": pd.read_csv(filepath, sep="\t", dtype=object)}
        except UnicodeDecodeError:
            return {
                "_sheet0": pd.read_csv(filepath, sep="\t", dtype=object,
                                        encoding="latin-1"),
            }
    if suffix == ".html":
        # pandas.read_html returns a list[DataFrame], one per <table>
        tables = pd.read_html(filepath)
        return {f"_table{i}": t for i, t in enumerate(tables)}
    raise ValueError(f"Unsupported file extension: {suffix}")


def _read_merged_ranges_xlsx(
    filepath: str,
) -> dict[str, list[tuple[int, int, int, int]]]:
    """Return `{sheet_name: [(r_top, c_left, r_bottom, c_right), ...]}` for
    every merged cell range, in 0-based row/col coords matching pandas
    iterrows() (i.e., header row = 0, body rows start at 1).

    openpyxl's `merged_cells.ranges` is 1-based with header row included;
    we shift to pandas semantics.
    """
    import warnings
    import openpyxl

    out: dict[str, list[tuple[int, int, int, int]]] = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            wb = openpyxl.load_workbook(filepath, data_only=True, read_only=False)
        except Exception:
            return out
        try:
            for sn in wb.sheetnames:
                ws = wb[sn]
                ranges = []
                for rng in ws.merged_cells.ranges:
                    # openpyxl is 1-based, row 1 = header; subtract 2 to get
                    # body-row index for first data row.
                    ranges.append((
                        rng.min_row - 2, rng.min_col - 1,
                        rng.max_row - 2, rng.max_col - 1,
                    ))
                out[str(sn)] = ranges
        finally:
            wb.close()
    return out


def _propagate_merged_values(
    df: pd.DataFrame,
    ranges: list[tuple[int, int, int, int]],
) -> None:
    """In-place: for each merged range, copy the top-left cell's value to
    every other cell in the range. Pandas leaves NaN there; the detector
    needs the visually-occupied value to flag corruption uniformly."""
    if not ranges:
        return
    columns = list(df.columns)
    for r_top, c_left, r_bottom, c_right in ranges:
        if r_top < 0 or r_top >= len(df) or c_left < 0 or c_left >= len(columns):
            continue
        top_col = columns[c_left]
        top_value = df.iloc[r_top][top_col]
        if pd.isna(top_value):
            continue
        for r in range(r_top, min(r_bottom + 1, len(df))):
            for c in range(c_left, min(c_right + 1, len(columns))):
                col = columns[c]
                if r == r_top and c == c_left:
                    continue
                # Only overwrite NaN : never clobber an existing value
                if pd.isna(df.iloc[r][col]):
                    df.iloc[r, df.columns.get_loc(col)] = top_value


def _calamine_rows_to_dataframe(rows: list[list]) -> pd.DataFrame:
    """Convert a list-of-lists from CalamineSheet.to_python() to a
    pandas DataFrame with header in row 0 (mirrors pd.read_excel).

    Duplicate column names are suffixed `.1`, `.2`, ... to match pandas'
    convention. Without this, `df[col]` returns a DataFrame instead of a
    Series for duplicated columns, and downstream `int(idx)` calls fail
    because `.items()` iterates over column names instead of row indices.
    Empirically caught: `PMC11150474/41467_2024_49202_MOESM4_ESM.xlsx`
    has seven `DayN` columns each repeated twice; without the suffix,
    detect() crashed with `ValueError: invalid literal for int(): 'Day1'`.
    """
    if not rows:
        return pd.DataFrame()
    raw_header = [
        str(c) if c is not None and c != "" else f"Unnamed: {i}"
        for i, c in enumerate(rows[0])
    ]
    # Deduplicate column names (pandas-style: `Day1`, `Day1.1`, `Day1.2`)
    seen: dict[str, int] = {}
    header: list[str] = []
    for name in raw_header:
        if name in seen:
            seen[name] += 1
            header.append(f"{name}.{seen[name]}")
        else:
            seen[name] = 0
            header.append(name)
    body = rows[1:]
    n_cols = len(header)
    body = [(r + [None] * (n_cols - len(r)))[:n_cols] for r in body]
    df = pd.DataFrame(body, columns=header, dtype=object)
    return df


def _load_xlsx_calamine(filepath: str) -> dict[str, pd.DataFrame]:
    """Fast xlsx loader for large files via python-calamine.

    Trade-off: number_format metadata is not available (calamine doesn't
    expose it). Date cells still return as datetime.date / datetime.datetime
    objects, so the detector's main `gene-date` path is unaffected.
    """
    from python_calamine import CalamineWorkbook
    wb = CalamineWorkbook.from_path(filepath)
    out: dict[str, pd.DataFrame] = {}
    for name in wb.sheet_names:
        sheet = wb.get_sheet_by_name(name)
        df = _calamine_rows_to_dataframe(sheet.to_python())
        df.attrs["cell_formats"] = {}  # calamine doesn't expose formats
        df.attrs["loader"] = "calamine"
        out[name] = df
    return out


def _load_xls_calamine(filepath: str) -> dict[str, pd.DataFrame]:
    """Fast xls loader for large BIFF8 files. xls date cells come back as
    raw integer serials (calamine 0.6.2 doesn't auto-decode for xls).
    The detector's Pass-1.5 numeric-decode path handles this correctly : 
    integers in [_SERIAL_MIN, _SERIAL_MAX] are tested as date serials.
    """
    from python_calamine import CalamineWorkbook
    wb = CalamineWorkbook.from_path(filepath)
    out: dict[str, pd.DataFrame] = {}
    for name in wb.sheet_names:
        sheet = wb.get_sheet_by_name(name)
        df = _calamine_rows_to_dataframe(sheet.to_python())
        df.attrs["cell_formats"] = {}
        df.attrs["loader"] = "calamine-xls"
        out[name] = df
    return out


def _read_cell_formats_xlsx(
    filepath: str,
) -> dict[str, dict[tuple[str, int], str]]:
    """Extract `number_format` for every date-typed cell in every sheet.

    Returns `{sheet_name: {(column_name, row_index): number_format}}`. Keys
    align with what pandas exposes : `column_name` is the pandas column label
    (header row 1, or `Unnamed: N` when missing), `row_index` is the 0-based
    data row matching `df.iterrows()`.

    Only date-typed cells are kept; everything else stays out of the dict
    (the detector only consults format strings on date values).
    """
    import warnings

    import openpyxl

    out: dict[str, dict[tuple[str, int], str]] = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        wb = openpyxl.load_workbook(filepath, data_only=True, read_only=True)
        for sn in wb.sheetnames:
            ws = wb[sn]
            sheet_map: dict[tuple[str, int], str] = {}
            header: list[str] = []
            for row_i, row in enumerate(ws.iter_rows(values_only=False)):
                if row_i == 0:
                    # Match pandas' rule: empty header cells become "Unnamed: N"
                    header = [
                        (str(c.value) if c.value is not None else f"Unnamed: {i}")
                        for i, c in enumerate(row)
                    ]
                    continue
                pandas_row = row_i - 1  # pandas iterrows is 0-based, post-header
                for col_i, c in enumerate(row):
                    if c.is_date and c.value is not None and c.number_format:
                        col_name = header[col_i] if col_i < len(header) else f"Unnamed: {col_i}"
                        sheet_map[(col_name, pandas_row)] = c.number_format
            out[str(sn)] = sheet_map
        wb.close()
    return out


def analyze(file: gr.utils.NamedString | None) -> tuple[pd.DataFrame, str, str | None]:
    if file is None:
        return pd.DataFrame(), "No file uploaded.", None
    df = _load_file(file.name)
    report = detect(df)
    summary = (
        f"Scanned {report.rows_scanned} rows x {report.columns_scanned} columns.\n"
        f"Identifier-shaped columns: {report.identifier_columns or 'none'}\n"
        f"Suspicions found: {len(report.suspicions)}"
    )
    rows = [
        {
            "column": s.column,
            "row": s.row,
            "value": repr(s.value),
            "kind": s.kind,
            "suggestion": s.suggestion or "",
            "confidence": s.confidence,
            "reason": s.reason,
        }
        for s in report.suspicions
    ]

    # Emit a Frictionless Table Schema sidecar : the prevention side. If the
    # scientist commits this JSON next to their data, downstream readers will
    # catch corruption at import time instead of after.
    schema = emit_table_schema(df, report, name=Path(file.name).stem)
    sidecar_path = Path(tempfile.mkdtemp(prefix="gc_schema_")) / (
        Path(file.name).name + ".schema.json"
    )
    sidecar_path.write_text(json.dumps(schema, indent=2))

    return pd.DataFrame(rows), summary, str(sidecar_path)


def classify_with_local_llm(
    file: "gr.utils.NamedString | None",
    progress=None,
) -> pd.DataFrame:
    """Optional second step: for each detected suspicion, run a local LLM
    over a structured review packet and return per-cell verdicts.

    Uses llama-cpp-python with Qwen2.5-3B-Instruct-GGUF (model downloads
    on first run, ~2.1 GB cached after).
    """
    if file is None:
        return pd.DataFrame()
    from .corpus import load_ziemann_2021_corpus
    from .detector import detect_file
    from .local_classifier import classify_packet
    from .subagent_packet import build_packet_for_suspicion, format_packet_as_markdown

    report = detect_file(file.name)
    if not report.suspicions:
        return pd.DataFrame([{"info": "No suspicions to classify."}])

    # Best-effort Ziemann context (will be empty for non-Ziemann files)
    try:
        z_entries = {e.pmc_id: e for e in load_ziemann_2021_corpus()}
    except Exception:
        z_entries = {}

    rows = []
    for s in report.suspicions[:50]: # cap at 50 per UI invocation
        susp_dict = {
            "pmc_id": "(uploaded)", "sheet": s.sheet, "column": s.column,
            "row": s.row, "flagged_value": repr(s.value), "kind": s.kind,
            "suggestion": s.suggestion, "confidence": s.confidence,
            "reason": s.reason, "file_name": Path(file.name).name,
        }
        pkt = build_packet_for_suspicion(susp_dict, Path(file.name), None)
        if pkt is None:
            continue
        md = format_packet_as_markdown(pkt)
        try:
            result = classify_packet(md)
        except Exception as exc:
            rows.append({
                "column": s.column, "row": s.row, "value": repr(s.value)[:40],
                "verdict": "ERROR", "reason": str(exc)[:80],
            })
            continue
        rows.append({
            "column": s.column, "row": s.row, "value": repr(s.value)[:40],
            "suggestion": s.suggestion, "verdict": result.verdict,
            "reason": result.reason[:120],
        })
    return pd.DataFrame(rows)


def _build_corrections_table(file: "gr.utils.NamedString | None"):
    """Detector pass + propose corrections, return the editable table + analytics.

    Output shape:
      (proposals_dataframe, message_md, original_corrections_state,
       summary_md, by_kind_df, by_family_df, by_confidence_df, by_sheet_col_df).
    """
    from .corrections import (
        compute_analytics, corrections_to_dataframe, propose_corrections,
    )
    if file is None:
        empty = pd.DataFrame({"bucket": [], "count": []})
        return (
            pd.DataFrame([{"info": "Upload a file first."}]),
            "Upload a file to compute proposed corrections.",
            [],
            "Upload a file to see batch analytics.",
            empty, empty, empty, empty,
        )
    corrs = propose_corrections(file.name)
    if not corrs:
        empty = pd.DataFrame({"bucket": [], "count": []})
        return (
            pd.DataFrame([{"info": "No corruption detected : nothing to correct."}]),
            "0 proposed corrections (nothing flagged).",
            [],
            "**No corruptions found** : file is clean.",
            empty, empty, empty, empty,
        )
    df = corrections_to_dataframe(corrs)
    analytics = compute_analytics(corrs)
    n_high = sum(1 for c in corrs if c.confidence >= 0.85)
    msg = (
        f"{len(corrs)} proposed corrections. {n_high} pre-accepted "
        f"(confidence >= 0.85); the rest are unchecked : review and tick "
        f"`accept` on the ones you want applied. You can also edit the "
        f"`proposed` column to choose between alternate candidates "
        f"(shown in `alt_candidates`)."
    )
    return (
        df, msg, corrs,
        analytics["summary_md"],
        analytics["by_kind_df"],
        analytics["by_family_df"],
        analytics["by_confidence_df"],
        analytics["by_sheet_col_df"],
    )


def _toggle_all_accept(select_all: bool, current_df):
    """Set every row's `accept` column to `select_all`."""
    if current_df is None or len(current_df) == 0:
        return current_df
    df = current_df.copy()
    if "accept" in df.columns:
        df["accept"] = bool(select_all)
    return df


def _apply_accepted_corrections(
    file: "gr.utils.NamedString | None",
    edited_df,
    original_corrections,
):
    """Take the user-edited Dataframe and the original proposals, apply only
    the accepted rows, return (corrected_xlsx_path, audit_log_path,
    summary_text)."""
    from .corrections import apply_corrections, dataframe_to_decisions
    if file is None or original_corrections is None or len(original_corrections) == 0:
        return None, None, "Nothing to apply : re-run detection first."
    decisions = dataframe_to_decisions(edited_df, original_corrections)
    if not decisions:
        return None, None, (
            "0 rows accepted. Tick `accept` for the corrections you want, "
            "then click *Apply selected* again."
        )
    if Path(file.name).suffix.lower() != ".xlsx":
        return None, None, (
            "Corrections currently support .xlsx only. Re-save your file as "
            "xlsx and re-upload, or apply corrections via the Python API."
        )
    corrected, audit, summary = apply_corrections(file.name, decisions)
    msg = (
        f"Applied {summary['n_applied']} / {summary['n_decisions']} corrections "
        f"({summary['n_skipped']} skipped). "
        f"Corrected file + audit log ready for download."
    )
    return corrected, audit, msg


def build_app() -> gr.Blocks:
    with gr.Blocks(title="UnCorrupt : Excel-corruption detector for scientific spreadsheets") as app:
        gr.Markdown(
            "# UnCorrupt : Excel-corruption detector for scientific spreadsheets\n"
            "Upload a spreadsheet. Detects Excel-style auto-corruption of gene "
            "symbols and alphanumeric IDs (Ziemann 2016, 2021).\n\n"
            "After analysis, a **Frictionless Table Schema sidecar** is generated "
            "for download : commit it next to your data file and downstream "
            "readers will catch corruption at import time.\n\n"
            "**Optional steps:**\n"
            "- *Propose corrections* : every flagged cell becomes an editable row "
            "with the proposed gene symbol; review per-row diffs, tick `accept`, "
            "and download the corrected xlsx.\n"
            "- *Classify with local LLM* : run a Qwen2.5-3B model over each "
            "suspicion for TP / FP / INCONCLUSIVE verdicts (model auto-downloads, "
            "~2.1 GB on first use)."
        )
        file_in = gr.File(label="Spreadsheet (.xlsx, .csv, .tsv)")

        with gr.Tab("1. Detect"):
            summary_out = gr.Textbox(label="Summary", lines=4)
            table_out = gr.Dataframe(label="Suspicions", interactive=False)
            schema_out = gr.File(
                label="Schema sidecar (Frictionless Table Schema)",
                interactive=False,
            )

        with gr.Tab("2. Review & batch-correct"):
            gr.Markdown(
                "Each row is one flagged cell. The **proposed** column is the "
                "first-candidate suggestion (e.g. SEPT2 from `SEPT2 | SEP2`); "
                "alternate candidates are in **alt_candidates**. Edit `proposed` "
                "to pick a different candidate per row. Click *Apply selected* to "
                "download a new xlsx with only the ticked rows fixed."
            )
            propose_btn = gr.Button("Propose corrections (re-runs detection)")

            # Broad-view batch analytics : shown immediately after proposing
            with gr.Accordion("Batch analytics", open=True):
                analytics_summary = gr.Markdown(
                    "Click **Propose corrections** above to compute analytics."
                )
                with gr.Row():
                    plot_by_kind = gr.BarPlot(
                        pd.DataFrame({"kind": [], "count": []}),
                        x="kind", y="count",
                        title="Suspicions by detection kind",
                        x_title="kind", y_title="count",
                        height=260,
                    )
                    plot_by_family = gr.BarPlot(
                        pd.DataFrame({"family": [], "count": []}),
                        x="family", y="count",
                        title="Proposed corrections by gene family",
                        x_title="family", y_title="count",
                        height=260,
                    )
                with gr.Row():
                    plot_by_confidence = gr.BarPlot(
                        pd.DataFrame({"confidence": [], "count": []}),
                        x="confidence", y="count",
                        title="Confidence distribution",
                        x_title="confidence bucket", y_title="count",
                        height=260,
                    )
                    plot_by_sheet_col = gr.BarPlot(
                        pd.DataFrame({"sheet_column": [], "count": []}),
                        x="sheet_column", y="count",
                        title="Top sheets×columns (max 15)",
                        x_title="sheet :: column", y_title="count",
                        height=260,
                    )

            corrections_msg = gr.Markdown()
            select_all = gr.Checkbox(label="Accept ALL proposed corrections",
                                       value=False)
            corrections_table = gr.Dataframe(
                label="Proposed corrections (editable)",
                datatype=[
                    "bool", "str", "str", "number", "str", "str",
                    "str", "number", "str", "str",
                ],
                interactive=True,
                wrap=True,
                column_widths=["90px", "100px", "150px", "60px", "150px",
                                "150px", "150px", "90px", "100px", "300px"],
            )
            with gr.Row():
                apply_btn = gr.Button("Apply selected → download corrected xlsx",
                                       variant="primary")
                apply_all_btn = gr.Button("Apply ALL (one-click batch correct)",
                                          variant="primary")
            apply_msg = gr.Markdown()
            corrected_file = gr.File(label="Corrected spreadsheet (download)",
                                       interactive=False)
            audit_file = gr.File(label="Corrections audit log (JSON)",
                                  interactive=False)
            # Hidden state to remember the original Correction objects so we
            # can map back from the edited Dataframe to the right (sheet,
            # column, row) for writeback.
            original_corrections_state = gr.State([])

        with gr.Tab("3. LLM verdicts (optional)"):
            classify_btn = gr.Button("Classify each suspicion with local LLM "
                                       "(Qwen2.5-3B)")
            verdicts_out = gr.Dataframe(label="LLM verdicts", interactive=False)

        # Detection flow : fires automatically on upload
        file_in.change(
            fn=analyze, inputs=file_in,
            outputs=[table_out, summary_out, schema_out],
        )
        # Corrections flow
        propose_btn.click(
            fn=_build_corrections_table, inputs=file_in,
            outputs=[corrections_table, corrections_msg,
                     original_corrections_state,
                     analytics_summary,
                     plot_by_kind, plot_by_family,
                     plot_by_confidence, plot_by_sheet_col],
        )
        select_all.change(
            fn=_toggle_all_accept, inputs=[select_all, corrections_table],
            outputs=corrections_table,
        )
        apply_btn.click(
            fn=_apply_accepted_corrections,
            inputs=[file_in, corrections_table, original_corrections_state],
            outputs=[corrected_file, audit_file, apply_msg],
        )
        # Apply-ALL: toggle the select-all checkbox + re-run apply in one click
        def _apply_all_in_one(file, corrections):
            from .corrections import (
                apply_corrections, corrections_to_dataframe,
                dataframe_to_decisions,
            )
            if not corrections:
                # User clicked Apply ALL without first clicking Propose
                if file is None:
                    return None, None, "Upload a file first."
                from .corrections import propose_corrections
                corrections = propose_corrections(file.name)
                if not corrections:
                    return None, None, "No corruption found : nothing to apply."
            df = corrections_to_dataframe(corrections)
            df["accept"] = True
            decisions = dataframe_to_decisions(df, corrections)
            if Path(file.name).suffix.lower() != ".xlsx":
                return None, None, (
                    "Corrections currently support .xlsx only."
                )
            corrected, audit, summary = apply_corrections(file.name, decisions)
            msg = (
                f"Auto-applied **all {summary['n_applied']}** proposed "
                f"corrections (one-click batch). Corrected file + audit log "
                f"ready for download."
            )
            return corrected, audit, msg
        apply_all_btn.click(
            fn=_apply_all_in_one,
            inputs=[file_in, original_corrections_state],
            outputs=[corrected_file, audit_file, apply_msg],
        )
        # LLM verdicts
        classify_btn.click(
            fn=classify_with_local_llm, inputs=file_in, outputs=verdicts_out,
        )
    return app


if __name__ == "__main__":
    import os

    port = int(os.environ.get("GRADIO_SERVER_PORT", "7860"))
    build_app().launch(server_name="127.0.0.1", server_port=port)
