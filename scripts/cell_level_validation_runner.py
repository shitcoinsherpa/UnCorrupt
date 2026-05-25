"""Cell-level validation against Ziemann S2 - multiprocess edition.

Same logic as the single-process runner but uses a worker pool so the
walk completes in N-cores time instead of single-core time. Each worker
loads its own HGNC tables once (~1.3 GB) - keep WORKERS modest on
memory-constrained boxes.
"""
from __future__ import annotations

import json
import multiprocessing as mp
import os
import re
import time
from dataclasses import asdict, dataclass, field
from datetime import UTC, date, datetime
from pathlib import Path

S2_PATH = Path("/local/S2_articles_2026-05-14.xlsx")
CACHE_DIR = Path("/local/files")
LEDGER_PATH = Path("/work/cell_level_validation_fresh.jsonl")
WORKERS = int(os.environ.get("UNCORRUPT_WORKERS", "6"))


@dataclass
class CellValidation:
    pmc_id: str
    notes: str
    expected_genes: list[str]
    detector_suggestions: list[str]
    outcome: str
    n_detector_suspicions: int = 0


@dataclass
class CellValidationRun:
    timestamp: str
    n_entries: int
    n_files_cached: int
    n_tp: int
    n_fn_missed: int
    n_fn_suggestion: int
    n_ambiguous: int
    n_file_missing: int
    cell_level_recall: float
    cell_level_precision: float
    runtime_seconds: float
    workers: int
    results: list[CellValidation] = field(default_factory=list)


def _parse_ziemann_note(note: str) -> list[str]:
    """Mirror of uncorrupt_research.cell_level_validation._parse_ziemann_note.
    Uses the detector's authoritative `_reverse_gene_date` so the expected-
    gene set matches what the detector would emit for the same date."""
    from uncorrupt.detector import (
        _parse_date_string, _reverse_gene_date, _serial_to_date,
    )

    note = (note or "").strip()
    if not note or note.lower() in ("nan", "none"):
        return []

    candidates: set[str] = set()

    # Try ISO-format datetime "2014-03-06 00:00:00"
    iso_match = re.match(r"^(\d{4})-(\d{1,2})-(\d{1,2})", note)
    if iso_match:
        try:
            d = date(int(iso_match.group(1)), int(iso_match.group(2)),
                     int(iso_match.group(3)))
            candidates.update(_reverse_gene_date(d))
        except (ValueError, TypeError):
            pass

    for token in note.split(","):
        token = token.strip()
        if not token:
            continue
        if token.isdigit():
            n = int(token)
            d = _serial_to_date(n)
            if d is not None:
                candidates.update(_reverse_gene_date(d))
        for parsed in _parse_date_string(token):
            candidates.update(_reverse_gene_date(parsed))
        if re.fullmatch(r"[A-Z][A-Z0-9-]{1,14}", token):
            candidates.add(token)

    return sorted(candidates)


def _expected_cache_filename(pmc_id: str, affected_file_url: str) -> str:
    """Cache filename convention from the original validation script."""
    if not affected_file_url:
        return ""
    fname = affected_file_url.rsplit("/", 1)[-1]
    return f"{pmc_id}_{fname}"


def _process_one(entry_tuple: tuple) -> dict:
    """Worker: detect on the SPECIFIC file Ziemann's S2 row pinned the
    annotation to (per `affected_file_url`). Falls back to globbing all
    PMC files only when the URL is missing - mirrors the original
    validation methodology so the recall number is comparable."""
    from uncorrupt.detector import detect_file

    pmc_id, notes, affected_file_url = entry_tuple
    expected = _parse_ziemann_note(notes)
    if not expected:
        return {
            "pmc_id": pmc_id, "notes": notes,
            "expected_genes": [], "detector_suggestions": [],
            "outcome": "ambiguous-note", "n_detector_suspicions": 0,
        }
    # Pin to the specific affected file when possible
    if affected_file_url:
        target_name = _expected_cache_filename(pmc_id, affected_file_url)
        candidate = CACHE_DIR / target_name
        if candidate.exists():
            matches = [candidate]
        else:
            return {
                "pmc_id": pmc_id, "notes": notes,
                "expected_genes": expected, "detector_suggestions": [],
                "outcome": "file-missing", "n_detector_suspicions": 0,
            }
    else:
        matches = list(CACHE_DIR.glob(f"{pmc_id}_*"))
        matches = [m for m in matches if m.suffix.lower() in (".xls", ".xlsx")]
        if not matches:
            return {
                "pmc_id": pmc_id, "notes": notes,
                "expected_genes": expected, "detector_suggestions": [],
                "outcome": "file-missing", "n_detector_suspicions": 0,
            }
    suggestions: list[str] = []
    n_susp = 0
    for f in matches:
        try:
            report = detect_file(str(f), row_context_boost=True)
        except Exception:
            continue
        n_susp += len(report.suspicions)
        for s in report.suspicions:
            if s.suggestion:
                # Original splits only on '|'; match exactly.
                for tok in s.suggestion.split("|"):
                    tok = tok.strip().upper()
                    if tok:
                        suggestions.append(tok)
    match = any(e in suggestions for e in expected)
    if match:
        outcome = "TP"
    elif n_susp == 0:
        outcome = "FN-missed"
    else:
        outcome = "FN-suggestion-mismatch"
    return {
        "pmc_id": pmc_id, "notes": notes,
        "expected_genes": expected,
        "detector_suggestions": suggestions[:30],
        "outcome": outcome, "n_detector_suspicions": n_susp,
    }


def main() -> None:
    from uncorrupt.corpus import load_ziemann_2021_corpus

    t0 = time.monotonic()
    entries = load_ziemann_2021_corpus(S2_PATH)
    cached_pmcs = {p.name.split("_", 1)[0] for p in CACHE_DIR.glob("PMC*")}
    cached_entries = [e for e in entries if e.pmc_id in cached_pmcs]
    print(f"[mp-validate] {len(cached_entries)} entries with cached files "
          f"(of {len(entries)} total)", flush=True)
    print(f"[mp-validate] launching {WORKERS} worker processes", flush=True)

    # Worker payload: (pmc_id, notes, affected_file_url) tuples
    payloads = [(e.pmc_id, e.notes, e.affected_file_url) for e in cached_entries]

    results: list[dict] = []
    n_done = 0
    n_total = len(payloads)
    last_print = 0

    # `spawn` to avoid copying the (large) parent state; each worker
    # imports the package fresh.
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=WORKERS) as pool:
        # imap_unordered streams results as they complete -> live progress
        for r in pool.imap_unordered(_process_one, payloads, chunksize=4):
            results.append(r)
            n_done += 1
            if n_done - last_print >= 50:
                last_print = n_done
                n_tp = sum(1 for x in results if x["outcome"] == "TP")
                n_fnm = sum(1 for x in results if x["outcome"] == "FN-missed")
                n_fns = sum(1 for x in results
                            if x["outcome"] == "FN-suggestion-mismatch")
                elapsed = time.monotonic() - t0
                rate = n_done / elapsed
                eta = (n_total - n_done) / rate / 60 if rate > 0 else -1
                print(f"[mp-validate]   {n_done}/{n_total} "
                      f"(TP={n_tp}, FN-miss={n_fnm}, FN-sugg={n_fns}) "
                      f"rate={rate:.2f}/s eta={eta:.1f}min", flush=True)

    n_tp = sum(1 for x in results if x["outcome"] == "TP")
    n_fnm = sum(1 for x in results if x["outcome"] == "FN-missed")
    n_fns = sum(1 for x in results if x["outcome"] == "FN-suggestion-mismatch")
    n_amb = sum(1 for x in results if x["outcome"] == "ambiguous-note")
    n_miss = sum(1 for x in results if x["outcome"] == "file-missing")
    total_validatable = n_tp + n_fnm + n_fns
    recall = n_tp / total_validatable if total_validatable else 0.0

    run = CellValidationRun(
        timestamp=datetime.now(UTC).isoformat(),
        n_entries=len(cached_entries),
        n_files_cached=len(cached_entries) - n_miss,
        n_tp=n_tp, n_fn_missed=n_fnm, n_fn_suggestion=n_fns,
        n_ambiguous=n_amb, n_file_missing=n_miss,
        cell_level_recall=recall, cell_level_precision=0.0,
        runtime_seconds=time.monotonic() - t0,
        workers=WORKERS,
        results=[CellValidation(**r) for r in results],
    )

    LEDGER_PATH.parent.mkdir(parents=True, exist_ok=True)
    with LEDGER_PATH.open("w") as fh:
        fh.write(json.dumps(asdict(run)) + "\n")

    print()
    print(f"=== FRESH CELL-LEVEL VALIDATION RESULT ({WORKERS} workers) ===")
    print(f"  entries with notes: {run.n_entries}")
    print(f"  files cached:       {run.n_files_cached}")
    print(f"  TP (correctly fixed):      {run.n_tp}")
    print(f"  FN-missed (no flag):       {run.n_fn_missed}")
    print(f"  FN-suggestion-mismatch:    {run.n_fn_suggestion}")
    print(f"  ambiguous-note:            {run.n_ambiguous}")
    print(f"  file-missing:              {run.n_file_missing}")
    print(f"  CELL-LEVEL RECALL: {run.cell_level_recall:.4f} "
          f"({run.n_tp} / {total_validatable})")
    print(f"  runtime: {run.runtime_seconds:.1f}s "
          f"({run.runtime_seconds/60:.1f} min)")
    print(f"  ledger written: {LEDGER_PATH}")


if __name__ == "__main__":
    main()
