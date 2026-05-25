"""Sub-agent review packets for inconclusive detector flags.

For each cell our detector flagged but neither row-context xref nor
article-body-text could classify, build a structured review packet
containing the full context a sub-agent (LLM with tool access) needs to
decide: is this cell a real uncorrupt, or is it legitimate data?

Each packet is self-contained : no ambiguity, no missing context : so
the sub-agent's classification is grounded in actual evidence rather
than guess.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

import pandas as pd

from .app import _load_all_sheets
from .column_classifier import classify_column
from .corpus import load_ziemann_2021_corpus
from .detector import detect

PROJECT_ROOT = Path(__file__).resolve().parents[2]


@dataclass
class ReviewPacket:
    pmc_id: str
    article_journal: str
    article_year: int | None
    article_species: str
    ziemann_note: str

    file_name: str
    sheet_name: str | None
    column_name: str
    row_index: int

    cell_value: str  # repr
    cell_kind: str  # detector's kind classification
    cell_suggestion: str | None
    cell_confidence: float
    cell_reason: str

    column_type: str  # column_classifier verdict
    column_type_reason: str
    column_sample_values: list[str]

    row_context: dict[str, str]  # column_name → repr(cell value) for same row

    review_question: str  # prompt for the sub-agent


def build_packet_for_suspicion(
    suspicion_dict: dict,
    file_path: Path,
    ziemann_entry,
) -> ReviewPacket | None:
    """Build a structured review packet from a row-context-validator result."""
    sheets = _load_all_sheets(str(file_path))
    sheet = suspicion_dict.get("sheet")
    column = suspicion_dict["column"]
    row = suspicion_dict["row"]

    # Resolve the right sheet's dataframe
    df = sheets.get(sheet) if sheet else next(iter(sheets.values()))
    if df is None or column not in df.columns:
        return None

    series = df[column]
    cls = classify_column(series, column)

    # Sample values from the column (de-duplicated)
    sample = []
    for v in series.dropna().head(15):
        s = repr(v) if not isinstance(v, str) else v
        if s not in sample:
            sample.append(s[:60])

    # Row context: other columns at the same row index
    row_context: dict[str, str] = {}
    try:
        row_data = df.iloc[row]
        for col, val in row_data.items():
            if str(col) == column:
                continue
            if pd.isna(val):
                continue
            row_context[str(col)] = repr(val) if not isinstance(val, str) else val
            if len(row_context) >= 12:
                break
    except (IndexError, KeyError):
        pass

    return ReviewPacket(
        pmc_id=suspicion_dict["pmc_id"],
        article_journal=ziemann_entry.journal if ziemann_entry else "",
        article_year=ziemann_entry.year if ziemann_entry else None,
        article_species=ziemann_entry.species if ziemann_entry else "",
        ziemann_note=ziemann_entry.notes if ziemann_entry else "",
        file_name=suspicion_dict.get("file_name", file_path.name),
        sheet_name=sheet,
        column_name=column,
        row_index=row,
        cell_value=suspicion_dict.get("flagged_value", ""),
        cell_kind=suspicion_dict.get("kind", ""),
        cell_suggestion=suspicion_dict.get("suggestion"),
        cell_confidence=suspicion_dict.get("confidence", 0.0),
        cell_reason=suspicion_dict.get("reason", ""),
        column_type=cls.column_type,
        column_type_reason=cls.reason,
        column_sample_values=sample,
        row_context=row_context,
        review_question=(
            f"Is the cell value in column '{column}' at row {row} (sheet '{sheet}') "
            f"a REAL Excel-style gene-name corruption (originally a gene symbol "
            f"that got auto-converted by Excel to its current form), or is it "
            f"LEGITIMATE data (the value as stored is what the original authors "
            f"intended)?\n\n"
            f"Our detector flagged this cell and suggested it should be the gene: "
            f"{suspicion_dict.get('suggestion')!r}. "
            f"Use the column type, the column's other sample values, and the "
            f"surrounding row context to decide."
        ),
    )


def format_packet_as_markdown(packet: ReviewPacket) -> str:
    """Render a packet as markdown for sub-agent consumption."""
    lines = [
        f"# Cell Review Packet : {packet.pmc_id}",
        f"",
        f"## Article",
        f"- Journal: {packet.article_journal}",
        f"- Year: {packet.article_year}",
        f"- Species: {packet.article_species}",
        f"- Ziemann's recorded corruption note: `{packet.ziemann_note}`",
        f"",
        f"## File",
        f"- File: `{packet.file_name}`",
        f"- Sheet: `{packet.sheet_name}`",
        f"- Column: `{packet.column_name}`",
        f"- Row: {packet.row_index}",
        f"",
        f"## Flagged cell",
        f"- Value: `{packet.cell_value}`",
        f"- Detector kind: `{packet.cell_kind}`",
        f"- Detector suggestion: `{packet.cell_suggestion}`",
        f"- Detector confidence: {packet.cell_confidence}",
        f"- Detector reason: {packet.cell_reason}",
        f"",
        f"## Column-type analysis",
        f"- Classified as: **{packet.column_type}**",
        f"- Reason: {packet.column_type_reason}",
        f"- Sample values from the column: {packet.column_sample_values}",
        f"",
        f"## Row context (other columns at the same row)",
    ]
    if packet.row_context:
        for c, v in packet.row_context.items():
            lines.append(f"- `{c}`: `{v}`")
    else:
        lines.append("- (row had no other non-empty cells)")
    lines.extend([
        f"",
        f"## Review question",
        f"{packet.review_question}",
        f"",
        f"## Required output",
        f"Respond with exactly one of: `TP`, `FP`, or `INCONCLUSIVE`",
        f"followed by a one-sentence rationale.",
    ])
    return "\n".join(lines)
