"""Row-context validation: ground-truth check for detector flags.

For each cell our detector flags as a uncorrupt, look at all OTHER
cells in the same row for external gene IDs (UniProt, RefSeq, Ensembl,
Entrez, HGNC). Resolve those IDs to canonical gene symbols. Compare with
our detector's suggestion.

Possible outcomes:
- TP-confirmed: row contains an external ID that resolves to a gene
  matching our suggestion. Strong evidence we flagged a real corruption.
- FP-likely: row contains external IDs that resolve to a DIFFERENT gene
  than our suggestion. Evidence our suggestion is wrong (or the cell
  isn't actually corruption).
- inconclusive: row has no external IDs we can use.

This is automatic, deterministic ground truth : no manual review or model
opinions required. Sub-agent review can then focus on the inconclusive
subset.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

from .app import _load_all_sheets
from .corpus import load_hgnc
from .detector import detect, Suspicion
from .xref_lookup import XrefIndex, load_xref_index

PROJECT_ROOT = Path(__file__).resolve().parents[2]
LEDGER_PATH = PROJECT_ROOT / "results" / "row_context_validation.jsonl"


@dataclass
class CellGroundTruth:
    pmc_id: str
    file_name: str
    sheet: str | None
    column: str
    row: int
    flagged_value: str  # repr of the cell value
    suggestion: str | None
    confidence: float
    row_xrefs: dict[str, str]  # {external_id: resolved_gene}
    outcome: str  # "TP-confirmed", "FP-likely", "inconclusive"


@dataclass
class RowContextRun:
    timestamp: str
    n_articles: int
    n_files_processed: int
    n_flags_total: int
    n_tp_confirmed: int
    n_fp_likely: int
    n_inconclusive: int
    cell_precision: float  # TP-confirmed / (TP-confirmed + FP-likely)
    runtime_seconds: float
    results: list[CellGroundTruth] = field(default_factory=list)


def _canonical_set(suggestions: set[str]) -> set[str]:
    """Expand a set of suggestions to include their HGNC-rename equivalents.

    Our detector emits old gene-symbol names (SEPT2, MARCH1, OCT2, etc.)
    because those are what Excel actually corrupted. The row xref resolves
    to CURRENT HGNC symbols (SEPTIN2, MARCHF1, POU2F2, etc.). To compare
    fairly, expand the suggestion set to include both old and new forms.
    """
    hgnc = load_hgnc()
    expanded: set[str] = set()
    for s in suggestions:
        expanded.add(s)
        # If old → new map exists, include the new form
        if s in hgnc.prev_symbol_to_current:
            expanded.add(hgnc.prev_symbol_to_current[s])
        # If new → also include all old forms that map to it
        for old, new in hgnc.prev_symbol_to_current.items():
            if new == s:
                expanded.add(old)
        # Aliases too
        if s in hgnc.alias_to_current:
            expanded.add(hgnc.alias_to_current[s])
        for alias, new in hgnc.alias_to_current.items():
            if new == s:
                expanded.add(alias)
    return expanded


def _validate_suspicion_against_row(
    s: Suspicion, df: pd.DataFrame, xref: XrefIndex,
) -> CellGroundTruth:
    """Look at the row of a flagged cell and try to ground-truth via xrefs."""
    # Get the row containing the flagged cell
    try:
        row_series = df.iloc[s.row]
    except (IndexError, KeyError):
        return CellGroundTruth(
            pmc_id="", file_name="", sheet=s.sheet, column=s.column,
            row=s.row, flagged_value=repr(s.value),
            suggestion=s.suggestion, confidence=s.confidence,
            row_xrefs={}, outcome="row-unreadable",
        )

    # Find any xref-shaped cells in the same row (excluding the flagged cell)
    xrefs: dict[str, str] = {}
    for col, val in row_series.items():
        if str(col) == s.column:
            continue
        if not isinstance(val, str):
            continue
        resolved = xref.lookup(val.strip())
        if resolved is not None:
            xrefs[val.strip()] = resolved

    if not xrefs:
        outcome = "inconclusive"
    else:
        # Did any resolved gene match our suggestion? Compare AFTER
        # HGNC-rename expansion so SEPT2 (our suggestion) matches SEPTIN2
        # (the modern symbol the row xref resolves to).
        suggestion_genes: set[str] = set()
        if s.suggestion:
            for part in s.suggestion.split("|"):
                suggestion_genes.add(part.strip())
        expanded_suggestions = _canonical_set(suggestion_genes)
        resolved_set = set(xrefs.values())
        if expanded_suggestions & resolved_set:
            outcome = "TP-confirmed"
        else:
            outcome = "FP-likely"

    return CellGroundTruth(
        pmc_id="", file_name="", sheet=s.sheet, column=s.column,
        row=s.row, flagged_value=repr(s.value),
        suggestion=s.suggestion, confidence=s.confidence,
        row_xrefs=xrefs, outcome=outcome,
    )


def validate_file(file_path: Path, xref: XrefIndex,
                  pmc_id: str = "") -> list[CellGroundTruth]:
    """Run detector on a file; validate each flag against row-context xrefs."""
    try:
        sheets = _load_all_sheets(str(file_path))
    except Exception:
        return []

    results: list[CellGroundTruth] = []
    for sheet_name, df in sheets.items():
        rep = detect(df, sheet=sheet_name)
        for s in rep.suspicions:
            gt = _validate_suspicion_against_row(s, df, xref)
            gt.pmc_id = pmc_id
            gt.file_name = file_path.name
            gt.sheet = sheet_name
            results.append(gt)
    return results


def run_row_context_validation(
    article_dirs: list[Path],
    xref: XrefIndex | None = None,
) -> RowContextRun:
    """Run row-context validation across a set of cached articles."""
    import time
    if xref is None:
        xref = load_xref_index()
    t0 = time.monotonic()

    all_results: list[CellGroundTruth] = []
    n_files = 0
    for article_dir in article_dirs:
        pmc_id = article_dir.name
        # find xlsx/xls/csv files in the article dir
        for f in article_dir.iterdir():
            if f.suffix.lower() not in (".xlsx", ".xls", ".csv", ".tsv"):
                continue
            results = validate_file(f, xref, pmc_id=pmc_id)
            all_results.extend(results)
            n_files += 1

    n_tp = sum(1 for r in all_results if r.outcome == "TP-confirmed")
    n_fp = sum(1 for r in all_results if r.outcome == "FP-likely")
    n_inc = sum(1 for r in all_results if r.outcome == "inconclusive")
    n_total = len(all_results)
    precision = n_tp / max(1, n_tp + n_fp)

    return RowContextRun(
        timestamp=datetime.now(UTC).isoformat(),
        n_articles=len(article_dirs),
        n_files_processed=n_files,
        n_flags_total=n_total,
        n_tp_confirmed=n_tp,
        n_fp_likely=n_fp,
        n_inconclusive=n_inc,
        cell_precision=precision,
        runtime_seconds=time.monotonic() - t0,
        results=all_results,
    )


def append_to_ledger(run: RowContextRun, path: Path = LEDGER_PATH) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(json.dumps(asdict(run)) + "\n")


def format_report(run: RowContextRun) -> str:
    n_grounded = run.n_tp_confirmed + run.n_fp_likely
    lines = [
        f"=== ROW-CONTEXT GROUND-TRUTH VALIDATION ===",
        f"timestamp: {run.timestamp}",
        f"articles processed: {run.n_articles}",
        f"files processed: {run.n_files_processed}",
        f"flags total: {run.n_flags_total}",
        f"---",
        f"TP-confirmed (row xref matches our suggestion): {run.n_tp_confirmed:>6}",
        f"FP-likely (row xref disagrees with our suggestion): {run.n_fp_likely:>6}",
        f"Inconclusive (no row xref to check): {run.n_inconclusive:>6}",
        f"---",
        f"Cell-level precision: {run.cell_precision:.4f}  (TP / (TP+FP), n={n_grounded})",
        f"runtime: {run.runtime_seconds:.1f}s",
    ]
    return "\n".join(lines)
