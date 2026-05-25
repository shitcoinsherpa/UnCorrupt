"""Gradio UI for the corruption detector.

Single-screen drag-drop-download. Accepts one file or many. Runs the
multi-sheet detector (matches the CLI), shows a plain-English report,
and bundles every output (cleaned spreadsheet, audit log, JSON report,
schema sidecar) into one downloadable ZIP.
"""
from __future__ import annotations

import json
import tempfile
import zipfile
from pathlib import Path
from typing import Any

import gradio as gr
import pandas as pd

from .detector import detect_file
from .schema_export import emit_table_schema

# --------------------------------------------------------------------------
# Loader plumbing
# --------------------------------------------------------------------------

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
    extension to use. Returns None if the extension is correct."""
    with open(filepath, "rb") as f:
        head = f.read(2048)
    head_lc = head.lower()
    if b"<html" in head_lc or b"<!doctype html" in head_lc:
        for marker in _HTML_PLACEHOLDER_MARKERS:
            if marker.encode() in head_lc:
                raise UnrecoverableFile(
                    f"{Path(filepath).name}: HTML placeholder ({marker!r}); "
                    "real file was never downloaded"
                )
        if b"<table" in head_lc:
            return ".html"
        raise UnrecoverableFile(
            f"{Path(filepath).name}: HTML payload with no <table> element"
        )
    if b"\t" in head and b"\x00" not in head[:512]:
        return ".tsv"
    return None


_CALAMINE_SIZE_THRESHOLD = 5 * 1024 * 1024  # 5 MB


def _load_all_sheets(filepath: str) -> dict[str, pd.DataFrame]:
    """Load every sheet of a spreadsheet as {sheet_name: DataFrame}.

    Multi-sheet read is critical: single-sheet readers miss ~20pp recall
    on the Ziemann 2021 corpus. Cell-format metadata is attached on
    `df.attrs['cell_formats']` for xlsx so the detector can distinguish
    `mmm-yy` vs `dd-mmm` corruptions.
    """
    p = Path(filepath)
    suffix = p.suffix.lower()
    if suffix in (".xlsx", ".xls"):
        override = _sniff_extension_override(filepath)
        if override is not None:
            suffix = override

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
            _propagate_merged_values(df, merges_by_sheet.get(str(k), []))
            out[str(k)] = df
        return out
    if suffix == ".xls":
        if file_size >= 10 * 1024 * 1024:
            return _load_xls_calamine(filepath)
        result = pd.read_excel(filepath, sheet_name=None, dtype=object)
        return {str(k): v for k, v in result.items()}
    if suffix == ".csv":
        return {"_sheet0": pd.read_csv(filepath, dtype=object)}
    if suffix in (".tsv", ".txt"):
        try:
            return {"_sheet0": pd.read_csv(filepath, sep="\t", dtype=object)}
        except UnicodeDecodeError:
            return {
                "_sheet0": pd.read_csv(filepath, sep="\t", dtype=object,
                                        encoding="latin-1"),
            }
    if suffix == ".html":
        tables = pd.read_html(filepath)
        return {f"_table{i}": t for i, t in enumerate(tables)}
    raise ValueError(f"Unsupported file extension: {suffix}")


def _read_merged_ranges_xlsx(
    filepath: str,
) -> dict[str, list[tuple[int, int, int, int]]]:
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
                if pd.isna(df.iloc[r][col]):
                    df.iloc[r, df.columns.get_loc(col)] = top_value


def _calamine_rows_to_dataframe(rows: list[list]) -> pd.DataFrame:
    if not rows:
        return pd.DataFrame()
    raw_header = [
        str(c) if c is not None and c != "" else f"Unnamed: {i}"
        for i, c in enumerate(rows[0])
    ]
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
    return pd.DataFrame(body, columns=header, dtype=object)


def _load_xlsx_calamine(filepath: str) -> dict[str, pd.DataFrame]:
    from python_calamine import CalamineWorkbook

    wb = CalamineWorkbook.from_path(filepath)
    out: dict[str, pd.DataFrame] = {}
    for name in wb.sheet_names:
        sheet = wb.get_sheet_by_name(name)
        df = _calamine_rows_to_dataframe(sheet.to_python())
        df.attrs["cell_formats"] = {}
        df.attrs["loader"] = "calamine"
        out[name] = df
    return out


def _load_xls_calamine(filepath: str) -> dict[str, pd.DataFrame]:
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
                    header = [
                        (str(c.value) if c.value is not None else f"Unnamed: {i}")
                        for i, c in enumerate(row)
                    ]
                    continue
                pandas_row = row_i - 1
                for col_i, c in enumerate(row):
                    if c.is_date and c.value is not None and c.number_format:
                        col_name = header[col_i] if col_i < len(header) else f"Unnamed: {col_i}"
                        sheet_map[(col_name, pandas_row)] = c.number_format
            out[str(sn)] = sheet_map
        wb.close()
    return out


# --------------------------------------------------------------------------
# UI-facing helpers
# --------------------------------------------------------------------------

_KIND_LABELS = {
    "gene-date":                "Date mistaken for gene name",
    "gene-date-string":         "Date text in gene column",
    "gene-date-serial":         "Excel-serial integer in gene column",
    "id-float":                 "Float (gene ID lost to scientific notation)",
    "long-int-precision-loss":  "Integer too big for Excel (trailing digits lost)",
    "leading-zero-stripped":    "Leading zero stripped",
    "autofill-sequence":        "Autofill drag (Excel filled a sequence)",
    "cas-registry":             "Chemical CAS number read as a date",
    "decimal-comma":            "Decimal-comma locale collision",
    "time-coercion":            "Time string (plate well IDs)",
    "homoglyph":                "Unicode look-alike (Cyrillic / Greek / fullwidth)",
    "unrecognized-symbol":      "Symbol not in any registry",
    "file-too-large":           "File too large to scan safely",
}


def _confidence_label(c: float) -> str:
    if c >= 0.85:
        return "High"
    if c >= 0.50:
        return "Medium"
    if c >= 0.30:
        return "Low"
    return "Info"


def _scan_one(path: Path) -> dict[str, Any]:
    """Run detect_file on a single path. Return a dict with the report rows,
    a summary fragment, the JSON-serialisable report, and the schema sidecar
    path. Errors are returned as `error=...` instead of raised."""
    name = path.name
    try:
        report = detect_file(str(path))
    except UnrecoverableFile as exc:
        return {"file": name, "error": (
            "Looks like a download-failure placeholder (HTML or captcha page) "
            "instead of a real spreadsheet. Re-download and re-upload."
        ), "detail": str(exc)}
    except ValueError as exc:
        return {"file": name, "error": (
            f"Unsupported file type: {exc}. Supported: .xlsx, .xls, .csv, .tsv."
        ), "detail": str(exc)}
    except Exception as exc:
        return {"file": name, "error": (
            f"Could not read this file: {type(exc).__name__}: {exc}. "
            f"Try re-saving as .xlsx or .csv and uploading again."
        ), "detail": f"{type(exc).__name__}: {exc}"}

    rows = [
        {
            "File": name,
            "Sheet": s.sheet or "",
            "Column": s.column,
            "Row": s.row,
            "Original": repr(s.value),
            "Proposed fix": s.suggestion or "(no fix; ID was lost)",
            "What happened": _KIND_LABELS.get(s.kind, s.kind),
            "Confidence": _confidence_label(s.confidence),
            "Why": s.reason,
        }
        for s in report.suspicions
    ]

    bands = {"High": 0, "Medium": 0, "Low": 0, "Info": 0}
    for s in report.suspicions:
        bands[_confidence_label(s.confidence)] += 1

    return {
        "file": name,
        "rows_scanned": report.rows_scanned,
        "identifier_columns": list(report.identifier_columns),
        "suspicions": report.suspicions,
        "report_rows": rows,
        "bands": bands,
        "report_json": {
            "file": name,
            "rows_scanned": report.rows_scanned,
            "columns_scanned": report.columns_scanned,
            "identifier_columns": list(report.identifier_columns),
            "suspicions": [
                {
                    "sheet": s.sheet,
                    "column": s.column,
                    "row": s.row,
                    "value": repr(s.value),
                    "kind": s.kind,
                    "suggestion": s.suggestion,
                    "confidence": s.confidence,
                    "reason": s.reason,
                }
                for s in report.suspicions
            ],
        },
    }


def _build_schema_sidecar(path: Path) -> tuple[Path, dict] | None:
    """Re-derive a Frictionless Table Schema for the first sheet of `path`.
    Returns (sidecar_path, schema_dict) or None on failure."""
    try:
        sheets = _load_all_sheets(str(path))
        first = next(iter(sheets.values()))
        report = detect_file(str(path))
        schema = emit_table_schema(first, report, name=path.stem)
    except Exception:
        return None
    sidecar = Path(tempfile.mkdtemp(prefix="uncorrupt_schema_")) / (
        path.name + ".schema.json"
    )
    sidecar.write_text(json.dumps(schema, indent=2))
    return sidecar, schema


def _summarise_batch(scans: list[dict]) -> str:
    """Top-level markdown summary across all uploaded files."""
    total_files = len(scans)
    failed = [s for s in scans if s.get("error")]
    ok = [s for s in scans if not s.get("error")]
    total_susp = sum(len(s["suspicions"]) for s in ok)
    total_bands = {"High": 0, "Medium": 0, "Low": 0, "Info": 0}
    for s in ok:
        for k, v in s["bands"].items():
            total_bands[k] += v

    head_parts: list[str] = []
    if total_files == 1 and ok:
        head_parts.append(f"Scanned **{ok[0]['rows_scanned']:,} rows** of `{ok[0]['file']}`.")
    elif total_files > 1:
        head_parts.append(f"Scanned **{total_files} files**.")
    if failed:
        head_parts.append(f"{len(failed)} could not be read (see below).")

    if total_susp == 0 and ok:
        head_parts.append("**Nothing looks corrupted.** "
                          "You can still download the schema sidecar(s) to "
                          "prevent future corruption on re-import.")
    elif total_susp > 0:
        parts = [f"{n} {label.lower()}"
                 for label, n in total_bands.items() if n > 0]
        head_parts.append(
            f"**Found {total_susp:,} likely corrupted cells** "
            f"({', '.join(parts)}). The table below shows each one and the "
            f"fix UnCorrupt proposes. Click *Download cleaned files* below "
            f"to apply the high-confidence fixes."
        )

    if failed:
        bullets = "\n".join(f"  - `{s['file']}`: {s['error']}" for s in failed)
        head_parts.append("**Files that could not be scanned:**\n" + bullets)

    return "\n\n".join(head_parts) if head_parts else "Upload a spreadsheet to scan it."


def _aggregate_rows(scans: list[dict]) -> pd.DataFrame:
    """Flatten per-file rows into one table; hide the File column when only
    one file was uploaded."""
    rows: list[dict] = []
    for s in scans:
        rows.extend(s.get("report_rows", []))
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    single_file = len({r["File"] for r in rows}) == 1
    if single_file and "File" in df.columns:
        df = df.drop(columns=["File"])
    return df


# --------------------------------------------------------------------------
# Top-level UI callbacks
# --------------------------------------------------------------------------


def _analyze_upload(files):
    """Run on every file upload. Returns:
      (summary_md, report_table, has_xlsx_flag, scan_state)."""
    if not files:
        return (
            "Upload one or more spreadsheets to scan them.",
            gr.update(value=pd.DataFrame(), visible=False),
            False,
            [],
        )
    if not isinstance(files, list):
        files = [files]

    scans = [_scan_one(Path(f.name)) for f in files]
    summary = _summarise_batch(scans)
    table = _aggregate_rows(scans)
    has_xlsx = any(
        Path(f.name).suffix.lower() == ".xlsx" for f in files
    )
    return (
        summary,
        gr.update(value=table, visible=not table.empty),
        has_xlsx,
        scans,
    )


def _clean_label(has_xlsx: bool, has_files: bool) -> str:
    if not has_files:
        return "Download cleaned files (upload a spreadsheet first)"
    if not has_xlsx:
        return "Auto-clean supports .xlsx only — re-save as .xlsx"
    return "Download cleaned files (.zip)"


def _build_outputs_zip(files, scans: list[dict]) -> tuple[str | None, str]:
    """Apply auto-clean to every uploaded .xlsx and bundle:
      <basename>.cleaned.xlsx
      <basename>.audit.csv
      <basename>.report.json
      <basename>.schema.json
    plus a top-level `report.csv` aggregating findings across files.

    Returns (zip_path, status_markdown).
    """
    if not files:
        return None, "Upload a file first."
    if not isinstance(files, list):
        files = [files]

    from .corrections import (
        apply_corrections,
        corrections_to_dataframe,
        dataframe_to_decisions,
        propose_corrections,
    )

    tmpdir = Path(tempfile.mkdtemp(prefix="uncorrupt_outputs_"))
    zip_path = tmpdir / "uncorrupt_outputs.zip"
    n_total_corrupt = 0
    n_total_applied = 0
    n_files_processed = 0
    n_files_skipped: list[str] = []

    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        # Per-file payloads
        for f, scan in zip(files, scans, strict=False):
            src = Path(f.name)
            base = src.stem
            # Skip files that failed to scan
            if scan.get("error"):
                n_files_skipped.append(src.name)
                continue
            # Report JSON: always
            zf.writestr(
                f"{base}.report.json",
                json.dumps(scan["report_json"], indent=2),
            )
            # Schema sidecar: best-effort
            sidecar = _build_schema_sidecar(src)
            if sidecar:
                _, schema = sidecar
                zf.writestr(
                    f"{base}.schema.json",
                    json.dumps(schema, indent=2),
                )
            # Cleaned xlsx + audit: xlsx only
            if src.suffix.lower() != ".xlsx":
                n_files_skipped.append(src.name)
                continue
            try:
                corrections = propose_corrections(str(src))
            except Exception as exc:
                n_files_skipped.append(f"{src.name} ({type(exc).__name__})")
                continue
            n_total_corrupt += len(corrections)
            if not corrections:
                continue
            df = corrections_to_dataframe(corrections)
            df["accept"] = df["confidence"] >= 0.85
            if not df["accept"].any():
                continue
            decisions = dataframe_to_decisions(df, corrections)
            cleaned, audit, summary = apply_corrections(str(src), decisions)
            n_total_applied += summary["n_applied"]
            n_files_processed += 1
            zf.write(cleaned, f"{base}.cleaned.xlsx")
            zf.write(audit, f"{base}.audit.csv")

        # Aggregated top-level report.csv across all files
        agg = _aggregate_rows(scans)
        if not agg.empty:
            zf.writestr("report.csv", agg.to_csv(index=False))

    parts = [
        f"**Bundled {n_files_processed} cleaned file(s) into a ZIP.**",
        f"Applied **{n_total_applied:,} high-confidence fixes** across "
        f"{n_total_corrupt:,} total flagged cells.",
    ]
    if n_files_skipped:
        parts.append(
            "Skipped (not .xlsx, or failed to load): "
            + ", ".join(f"`{n}`" for n in n_files_skipped[:6])
            + (f" and {len(n_files_skipped) - 6} more" if len(n_files_skipped) > 6 else "")
            + "."
        )
    parts.append("Unzip the download to get one cleaned spreadsheet, audit "
                 "log, JSON report, and schema sidecar per input file, plus a "
                 "top-level `report.csv` aggregating every flag.")
    return str(zip_path), "\n\n".join(parts)


def _build_report_json(scans: list[dict]) -> str | None:
    """Write the aggregated JSON report to a tempfile and return its path."""
    if not scans:
        return None
    payload = {
        "version": 1,
        "tool": "uncorrupt",
        "files": [s.get("report_json") for s in scans if s.get("report_json")],
        "errors": [
            {"file": s["file"], "error": s["error"]}
            for s in scans if s.get("error")
        ],
    }
    path = Path(tempfile.mkdtemp(prefix="uncorrupt_report_")) / "uncorrupt_report.json"
    path.write_text(json.dumps(payload, indent=2))
    return str(path)


def _build_schema_sidecars_zip(files) -> str | None:
    """Bundle one schema sidecar per uploaded file into a single ZIP."""
    if not files:
        return None
    if not isinstance(files, list):
        files = [files]
    tmpdir = Path(tempfile.mkdtemp(prefix="uncorrupt_schemas_"))
    zip_path = tmpdir / "uncorrupt_schemas.zip"
    n_written = 0
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for f in files:
            src = Path(f.name)
            sidecar = _build_schema_sidecar(src)
            if not sidecar:
                continue
            _, schema = sidecar
            zf.writestr(src.name + ".schema.json",
                        json.dumps(schema, indent=2))
            n_written += 1
    if n_written == 0:
        return None
    return str(zip_path)


# --------------------------------------------------------------------------
# build_app
# --------------------------------------------------------------------------


def build_app() -> gr.Blocks:
    """Single screen. Drop one file or many. Get back a ZIP with cleaned
    spreadsheets, audit logs, JSON reports, and schema sidecars."""
    with gr.Blocks(title="UnCorrupt: find and fix Excel-mangled gene names") as app:
        gr.Markdown(
            "# UnCorrupt\n"
            "**Excel keeps turning your gene names into dates.** Drop one or "
            "more spreadsheets here. Get back a clean copy and an audit log.\n\n"
            "Catches the corruption families documented in the published "
            "literature: `SEPT2` becoming `2-Sep`, `OCT4` becoming `4-Oct`, "
            "RIKEN identifiers turning into floats like `2.31E+19`, leading "
            "zeros vanishing, Cyrillic look-alikes hiding inside ASCII names. "
            "Accepts `.xlsx`, `.xls`, `.csv`, `.tsv`. Files are processed in "
            "this container and discarded immediately; nothing is logged."
        )

        file_in = gr.File(
            label="Your spreadsheet(s) — drop one or many",
            file_count="multiple",
            file_types=[".xlsx", ".xls", ".csv", ".tsv"],
            height=160,
        )

        summary_out = gr.Markdown()
        table_out = gr.Dataframe(
            label="Every cell UnCorrupt flagged",
            interactive=False,
            wrap=True,
            visible=False,
        )

        # Shared state: list of per-file scan dicts; xlsx-presence flag
        scan_state = gr.State([])
        has_xlsx_state = gr.State(False)

        with gr.Row():
            clean_btn = gr.Button(
                "Download cleaned files (upload a spreadsheet first)",
                variant="primary",
                size="lg",
                interactive=False,
            )
            report_btn = gr.Button(
                "Download JSON report",
                variant="secondary",
                interactive=False,
            )
            schemas_btn = gr.Button(
                "Download schema sidecars",
                variant="secondary",
                interactive=False,
            )

        clean_msg = gr.Markdown()
        outputs_zip = gr.File(
            label="Cleaned files + audit + reports + schemas (.zip)",
            interactive=False,
            visible=False,
        )
        report_file = gr.File(
            label="JSON report (for methods sections and lab notebooks)",
            interactive=False,
            visible=False,
        )
        schemas_zip = gr.File(
            label="Frictionless schema sidecars (.zip)",
            interactive=False,
            visible=False,
        )

        # --- behaviour wiring ---

        def _on_upload(files):
            summary, table_update, has_xlsx, scans = _analyze_upload(files)
            has_files = bool(files)
            return (
                summary,
                table_update,
                scans,
                has_xlsx,
                gr.update(
                    value=_clean_label(has_xlsx, has_files),
                    interactive=has_xlsx and has_files,
                ),
                gr.update(interactive=has_files),
                gr.update(interactive=has_files),
                gr.update(value=None, visible=False),
                gr.update(value=None, visible=False),
                gr.update(value=None, visible=False),
                "",
            )

        file_in.change(
            fn=_on_upload,
            inputs=file_in,
            outputs=[
                summary_out, table_out, scan_state, has_xlsx_state,
                clean_btn, report_btn, schemas_btn,
                outputs_zip, report_file, schemas_zip, clean_msg,
            ],
        )

        def _on_clean(files, scans):
            zip_path, msg = _build_outputs_zip(files, scans)
            return (
                msg,
                gr.update(value=zip_path, visible=zip_path is not None),
            )

        clean_btn.click(
            fn=_on_clean,
            inputs=[file_in, scan_state],
            outputs=[clean_msg, outputs_zip],
        )

        def _on_report(scans):
            path = _build_report_json(scans)
            return gr.update(value=path, visible=path is not None)

        report_btn.click(
            fn=_on_report,
            inputs=scan_state,
            outputs=report_file,
        )

        def _on_schemas(files):
            path = _build_schema_sidecars_zip(files)
            return gr.update(value=path, visible=path is not None)

        schemas_btn.click(
            fn=_on_schemas,
            inputs=file_in,
            outputs=schemas_zip,
        )

        gr.Markdown(
            "### What the confidence labels mean\n\n"
            "- **High**: an independent source (a database identifier in the "
            "same row, or two or more similar corruptions in the same column) "
            "agrees. Auto-fixed in the cleaned download.\n"
            "- **Medium**: the pattern matches a known corruption family but "
            "there is no second source to confirm. Worth reviewing.\n"
            "- **Low**: a single suspicious-looking cell with no corroboration. "
            "Often a false alarm. Review.\n"
            "- **Info**: hidden by default; visible in the CLI for completeness."
        )

        gr.Markdown(
            "### What this tool does NOT do\n\n"
            "- **Float-coerced IDs are gone for good.** Once Excel turned "
            "`2310009E13` into the number `2.31E+19`, the original digits are "
            "mathematically unrecoverable. UnCorrupt flags the cell so you "
            "know to re-fetch from your source. It cannot reconstruct.\n"
            "- **Real dates in a date column are left alone.** Column context "
            "tells the difference between a date that belongs there and a "
            "date that used to be `SEPT2`.\n"
            "- **Brand-new corruption patterns we have never seen.** Knows "
            "what Ziemann 2016, Abeysooriya 2021, and Koh 2022 catalogued. "
            "If your file did something genuinely novel, file an issue with "
            "the file attached and we add it.\n"
            "- **Auto-clean writes `.xlsx` only.** CSV/TSV uploads still get "
            "a JSON report and a schema sidecar; the cleaned-spreadsheet "
            "rewrite path is xlsx-specific."
        )

    return app


def _launch() -> None:
    """Entry point for the `uncorrupt-app` console script (pip install)."""
    import os
    port = int(os.environ.get("GRADIO_SERVER_PORT", "7860"))
    host = os.environ.get("GRADIO_SERVER_NAME", "127.0.0.1")
    build_app().launch(server_name=host, server_port=port, theme=gr.themes.Soft())


if __name__ == "__main__":
    _launch()
