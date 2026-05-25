"""Batch-correction support for the Gradio UI.

Workflow per Mark Ziemann's "never silently repair data" principle:

    1. `propose_corrections(file_path)` runs the detector and proposes a
       replacement value for each flagged cell. Multi-candidate suggestions
       (e.g. "SEPT2 | SEP2") default to the FIRST candidate (the modern
       symbol) but the user can edit the proposed value before applying.

    2. `apply_corrections(file_path, decisions)` writes a NEW xlsx file with
       only the user-accepted corrections applied. The original file is
       never modified in place. Each applied correction is recorded in a
       sidecar audit log so the change history is recoverable.

Output paths are returned as absolute strings so Gradio's `gr.File` widget
can serve them for download.
"""
from __future__ import annotations

import json
import shutil
import tempfile
import warnings
from dataclasses import asdict, dataclass
from datetime import UTC, date, datetime
from pathlib import Path
from typing import Any

import openpyxl

from .detector import Suspicion, detect_file


@dataclass
class Correction:
    """A proposed replacement for one flagged cell."""
    sheet: str | None
    column: str               # pandas column label
    row: int                  # pandas 0-based row index
    original_value: Any       # the cell value as the detector saw it
    proposed_value: str       # the first-candidate suggestion (user-editable)
    all_candidates: list[str] # full suggestion split
    confidence: float
    kind: str
    reason: str


def _first_candidate(suggestion: str | None) -> str:
    """Return the first '|'-separated candidate from a suggestion string.

    Convention: the modern HGNC symbol comes first (SEPT2 | SEP2 → SEPT2).
    """
    if not suggestion:
        return ""
    return suggestion.split("|")[0].strip()


def _all_candidates(suggestion: str | None) -> list[str]:
    if not suggestion:
        return []
    return [c.strip() for c in suggestion.split("|") if c.strip()]


def propose_corrections(file_path: str) -> list[Correction]:
    """Run the detector and propose a replacement for each suspicion."""
    report = detect_file(file_path)
    out: list[Correction] = []
    for s in report.suspicions:
        out.append(Correction(
            sheet=s.sheet,
            column=str(s.column),
            row=int(s.row),
            original_value=s.value,
            proposed_value=_first_candidate(s.suggestion),
            all_candidates=_all_candidates(s.suggestion),
            confidence=float(s.confidence),
            kind=s.kind,
            reason=s.reason,
        ))
    return out


def _resolve_openpyxl_cell(ws, pandas_column: str, pandas_row: int) -> "openpyxl.cell.Cell | None":
    """Map pandas (column_label, 0-based row) → openpyxl cell.

    Convention used elsewhere in this codebase:
      - openpyxl row_i=0 (1-based row 1) is the header
      - pandas row 0 = openpyxl 1-based row 2 = row_i=1
      - column name = header cell text, OR "Unnamed: N" when header was empty
    """
    target_row = pandas_row + 2  # +1 for header, +1 for openpyxl 1-indexing

    # Find the column letter by walking the header row
    header_row = ws[1]
    target_col_idx: int | None = None
    for col_i, c in enumerate(header_row):
        header_value = str(c.value) if c.value is not None else f"Unnamed: {col_i}"
        if header_value == pandas_column:
            target_col_idx = col_i + 1  # openpyxl is 1-indexed
            break
    if target_col_idx is None:
        return None
    return ws.cell(row=target_row, column=target_col_idx)


def apply_corrections(file_path: str,
                       decisions: list[Correction],
                       output_dir: Path | None = None,
                       ) -> tuple[str, str, dict]:
    """Write a NEW xlsx file with the accepted corrections applied.

    Returns (corrected_xlsx_path, audit_log_path, summary_dict).

    `decisions` is the LIST of corrections to apply (already filtered to the
    user-accepted ones with the user's final `proposed_value` per row).
    """
    src = Path(file_path)
    if output_dir is None:
        output_dir = Path(tempfile.mkdtemp(prefix="gc_corrected_"))
    dst = output_dir / src.name.replace(src.suffix, f".corrected{src.suffix}")
    audit_path = output_dir / f"{src.name}.corrections-audit.json"

    # Copy source -> destination, then edit in place. openpyxl can save
    # only files it loaded with load_workbook (read_only=False).
    shutil.copy2(src, dst)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        wb = openpyxl.load_workbook(dst, data_only=False)

    applied: list[dict] = []
    skipped: list[dict] = []
    for d in decisions:
        if not d.proposed_value:
            skipped.append({**asdict(d), "skip_reason": "empty proposed_value"})
            continue
        sheet_name = d.sheet
        if sheet_name is None or sheet_name not in wb.sheetnames:
            # Single-sheet file: pandas read produces sheet=None
            if len(wb.sheetnames) == 1:
                sheet_name = wb.sheetnames[0]
            else:
                skipped.append({**asdict(d), "skip_reason": f"sheet {d.sheet!r} not in workbook"})
                continue
        ws = wb[sheet_name]
        cell = _resolve_openpyxl_cell(ws, d.column, d.row)
        if cell is None:
            skipped.append({**asdict(d), "skip_reason": f"could not resolve cell {d.column!r}/{d.row}"})
            continue
        # Preserve the original value for the audit trail
        original_repr = repr(cell.value)
        # Overwrite with the user-approved proposed value (always a string : 
        # gene symbols are always strings).
        cell.value = d.proposed_value
        # Clear the date number_format so the cell stops displaying as a date.
        if cell.number_format and ("y" in cell.number_format.lower()
                                    or "d" in cell.number_format.lower()
                                    or "m" in cell.number_format.lower()):
            cell.number_format = "General"
        applied.append({
            "sheet": sheet_name,
            "openpyxl_coord": cell.coordinate,
            "column": d.column,
            "row": d.row,
            "original_value": original_repr,
            "corrected_to": d.proposed_value,
            "all_candidates_offered": d.all_candidates,
            "confidence": d.confidence,
            "kind": d.kind,
        })

    wb.save(dst)
    wb.close()

    summary = {
        "source_file": str(src),
        "corrected_file": str(dst),
        "audit_log": str(audit_path),
        "n_decisions": len(decisions),
        "n_applied": len(applied),
        "n_skipped": len(skipped),
        "applied_at": datetime.now(UTC).isoformat(),
        "applied": applied,
        "skipped": skipped,
    }
    audit_path.write_text(json.dumps(summary, indent=2, default=str))
    return str(dst), str(audit_path), summary


def corrections_to_dataframe(corrections: list[Correction]):
    """Render the proposals as a pandas DataFrame for Gradio.

    Columns: [accept, sheet, column, row, original, proposed, alt_candidates,
             confidence, kind, reason]. The `accept` column is bool : Gradio
    will render it as a checkbox when `datatype=["bool", ...]` is set on the
    Dataframe component.

    Default accept = True for high-confidence (>= 0.85) corrections,
    False for ambiguous ones (the human should disambiguate).
    """
    import pandas as pd
    return pd.DataFrame([
        {
            "accept": c.confidence >= 0.85,
            "sheet": c.sheet or "_",
            "column": c.column,
            "row": c.row,
            "original": str(c.original_value)[:50],
            "proposed": c.proposed_value,
            "alt_candidates": " | ".join(c.all_candidates[1:]) if len(c.all_candidates) > 1 else "",
            "confidence": c.confidence,
            "kind": c.kind,
            "reason": c.reason[:120],
        }
        for c in corrections
    ])


_GENE_FAMILY_PREFIXES = ("SEPT", "SEP", "MARCH", "MARC", "OCT", "DEC",
                          "NOV", "APR", "FEB")


def _family_of(proposed: str) -> str:
    """Bucket a proposed gene symbol into its Excel-corruption family."""
    s = proposed.strip().upper()
    for prefix in ("SEPTIN", "MARCHF", "POU2F"): # post-rename modern symbols
        if s.startswith(prefix):
            return prefix
    for prefix in _GENE_FAMILY_PREFIXES:
        if s.startswith(prefix) and s[len(prefix):].lstrip("-").isdigit():
            return prefix
    return "OTHER"


def _confidence_bucket(conf: float) -> str:
    if conf >= 0.85:
        return "high (>=0.85)"
    if conf >= 0.5:
        return "medium (0.5-0.85)"
    return "low (<0.5)"


def compute_analytics(corrections: list[Correction]) -> dict:
    """Broad-view batch analytics for the proposals.

    Returns:
        {
          "summary_md": markdown text summarizing the batch
          "by_kind_df": DataFrame for the per-detection-kind bar chart
          "by_family_df": DataFrame for the per-gene-family bar chart
          "by_confidence_df": DataFrame for the per-confidence-bucket bar chart
          "by_sheet_col_df": DataFrame for the per-(sheet,column) bar chart
          "totals": raw counts dict for tests
        }
    """
    import pandas as pd
    from collections import Counter

    if not corrections:
        empty = pd.DataFrame({"bucket": [], "count": []})
        return {
            "summary_md": "**No corruptions found** : file is clean.",
            "by_kind_df": empty,
            "by_family_df": empty,
            "by_confidence_df": empty,
            "by_sheet_col_df": empty,
            "totals": {"n_corrections": 0},
        }

    n = len(corrections)
    kinds = Counter(c.kind for c in corrections)
    families = Counter(_family_of(c.proposed_value) for c in corrections)
    conf_buckets = Counter(_confidence_bucket(c.confidence) for c in corrections)
    sheet_cols = Counter(
        f"{c.sheet or '_'}::{c.column}" for c in corrections
    )

    n_sheets = len({c.sheet or "_" for c in corrections})
    n_columns = len({(c.sheet or "_", c.column) for c in corrections})
    n_high = sum(1 for c in corrections if c.confidence >= 0.85)
    n_med = sum(1 for c in corrections if 0.5 <= c.confidence < 0.85)
    n_low = sum(1 for c in corrections if c.confidence < 0.5)

    summary_md = (
        f"### Batch summary\n"
        f"- **{n:,}** cells flagged for correction across "
        f"**{n_sheets}** sheet(s), **{n_columns}** distinct identifier "
        f"column(s)\n"
        f"- By confidence: **{n_high}** high (auto-accept default), "
        f"**{n_med}** medium (review), **{n_low}** low (review)\n"
        f"- By detection kind: " +
        ", ".join(f"**{k}** = {v}" for k, v in kinds.most_common()) + "\n"
        f"- By gene family (proposed): " +
        ", ".join(f"**{k}** = {v}"
                   for k, v in families.most_common() if k != "OTHER") +
        (f", OTHER = {families['OTHER']}" if families.get("OTHER") else "")
    )

    return {
        "summary_md": summary_md,
        "by_kind_df": pd.DataFrame(
            sorted(kinds.items(), key=lambda x: -x[1]),
            columns=["kind", "count"],
        ),
        "by_family_df": pd.DataFrame(
            sorted(families.items(), key=lambda x: -x[1]),
            columns=["family", "count"],
        ),
        "by_confidence_df": pd.DataFrame(
            [(b, conf_buckets.get(b, 0)) for b in
             ("high (>=0.85)", "medium (0.5-0.85)", "low (<0.5)")],
            columns=["confidence", "count"],
        ),
        "by_sheet_col_df": pd.DataFrame(
            sorted(sheet_cols.items(), key=lambda x: -x[1])[:15],  # top 15
            columns=["sheet_column", "count"],
        ),
        "totals": {
            "n_corrections": n,
            "n_sheets": n_sheets,
            "n_columns": n_columns,
            "n_high_confidence": n_high,
            "n_medium_confidence": n_med,
            "n_low_confidence": n_low,
            "n_kinds": dict(kinds),
            "n_families": dict(families),
        },
    }


def dataframe_to_decisions(df, original_corrections: list[Correction]) -> list[Correction]:
    """Reverse mapping: take the user-edited Dataframe + the original
    Correction list, return only the rows the user accepted, with their
    edited `proposed_value` honored."""
    out: list[Correction] = []
    if df is None or len(df) == 0:
        return out
    for i, row in df.iterrows():
        if not row.get("accept"):
            continue
        if i >= len(original_corrections):
            continue
        orig = original_corrections[i]
        # Allow the user to override the proposed value text
        user_proposed = str(row.get("proposed", orig.proposed_value)).strip()
        out.append(Correction(
            sheet=orig.sheet,
            column=orig.column,
            row=orig.row,
            original_value=orig.original_value,
            proposed_value=user_proposed,
            all_candidates=orig.all_candidates,
            confidence=orig.confidence,
            kind=orig.kind,
            reason=orig.reason,
        ))
    return out
