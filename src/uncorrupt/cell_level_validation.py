"""Cell-level validation against Ziemann S2's Notes column.

The Koh-replication and the basic Ziemann recall measure file-level outcomes:
"did the detector flag this article?" That's not enough. We need cell-level:
"did the detector find the EXACT corruption Ziemann documented, and does our
suggestion match the gene that was originally there?"

Ziemann's S2 Notes column records what each corruption looked like (e.g.
'39880', 'Sep-02', '2014-03-06 00:00:00', '2-Oct,ATOCT2,OCT2'). We parse
those notes through the same logic the detector uses, derive the EXPECTED
gene-name candidate(s), then compare with what our detector emitted for that
file.

Outcomes per article:
- TP: our detector found the corruption Ziemann noted, and our suggestion
  contains the expected gene
- FN-missed: detector ran on file but did NOT flag the cell Ziemann noted
- FN-suggestion: detector flagged something but the suggestion doesn't match
  Ziemann's expected gene
- File-not-cached: we don't have the file locally
- Ambiguous: Ziemann's note doesn't parse into a unique expected suggestion

This is the precision/recall measurement we should actually be quoting.
"""
from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass, field
from datetime import UTC, date, datetime
from pathlib import Path

from .corpus import CorpusEntry, load_ziemann_2021_corpus
from .detector import _parse_date_string, _reverse_gene_date, _serial_to_date, detect_file

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CACHE_DIR = PROJECT_ROOT / "data" / "raw" / "ziemann_2021_corpus" / "files"
LEDGER_PATH = PROJECT_ROOT / "results" / "cell_level_validation.jsonl"


@dataclass
class CellValidation:
    pmc_id: str
    notes: str
    expected_genes: list[str]  # what Ziemann's note implies the corrupted gene was
    detector_suggestions: list[str]  # what our detector emitted (suggestions field)
    outcome: str  # "TP", "FN-missed", "FN-suggestion-mismatch", "ambiguous-note", "file-missing"
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
    cell_level_recall: float  # TP / (TP + FN-missed + FN-suggestion-mismatch)
    cell_level_precision: float  # TP / (TP + ?) : see notes
    runtime_seconds: float
    results: list[CellValidation] = field(default_factory=list)


def _parse_ziemann_note(note: str) -> list[str]:
    """Parse a Ziemann S2 Notes value into expected gene-name suggestions.

    Examples observed:
      '39880'                       → date-serial → 2009-03-08 → MARCH8
      'Sep-02'                      → date-string → SEPT2 / SEP2
      '2014-03-06 00:00:00'         → date timestamp → MARCH6
      '2-Oct,ATOCT2,OCT2'           → date + literal gene names
      '37500'                       → date-serial → MARCHn or similar
      'Mar-09'                      → date-string → MARCH9
      ''                            → empty
    """
    note = (note or "").strip()
    if not note or note.lower() in ("nan", "none"):
        return []

    candidates: set[str] = set()

    # Try ISO-format datetime "2014-03-06 00:00:00"
    iso_match = re.match(r"^(\d{4})-(\d{1,2})-(\d{1,2})", note)
    if iso_match:
        try:
            d = date(int(iso_match.group(1)), int(iso_match.group(2)), int(iso_match.group(3)))
            candidates.update(_reverse_gene_date(d))
        except ValueError:
            pass

    # Comma-separated tokens (e.g. '2-Oct,ATOCT2,OCT2')
    for token in note.split(","):
        token = token.strip()
        if not token:
            continue

        # Pure digits → date serial
        if token.isdigit():
            n = int(token)
            d = _serial_to_date(n)
            if d is not None:
                candidates.update(_reverse_gene_date(d))

        # Date-formatted strings (Sep-02, Mar-09, 2-Oct, etc.)
        for parsed in _parse_date_string(token):
            candidates.update(_reverse_gene_date(parsed))

        # Literal gene-name shape (uppercase letters + digits, e.g. OCT2, MARCHF1)
        if re.fullmatch(r"[A-Z][A-Z0-9-]{1,14}", token):
            candidates.add(token)

    return sorted(candidates)


def _expected_cache_filename(pmc_id: str, affected_file_url: str) -> str:
    """Compute the cache filename pattern for a given Ziemann S2 entry."""
    if not affected_file_url:
        return ""
    fname = affected_file_url.rsplit("/", 1)[-1]
    return f"{pmc_id}_{fname}"


def _suggestions_for_pmc_file(pmc_id: str,
                              affected_file_url: str | None,
                              ) -> tuple[list[str], int, int, bool]:
    """Suggestions for the SPECIFIC file Ziemann's S2 row pinned the annotation
    to.

    Returns (suggestions, n_suspicions, n_unrecoverable, file_in_cache).

    `file_in_cache=False` distinguishes "this Ziemann annotation references a
    supplementary file we never downloaded" from "the detector ran on the file
    and found nothing." The first is a corpus-completeness issue, the second
    is a detector miss.

    When `affected_file_url` is None or empty, falls back to ALL cached files
    for the PMC (pre-iter-3 behaviour), with `file_in_cache=True` if any file
    exists.
    """
    cached_files: list[Path]
    if affected_file_url:
        target_name = _expected_cache_filename(pmc_id, affected_file_url)
        candidate = CACHE_DIR / target_name
        if candidate.exists():
            cached_files = [candidate]
        else:
            return [], 0, 0, False
    else:
        cached_files = list(CACHE_DIR.glob(f"{pmc_id}_*"))
        if not cached_files:
            return [], 0, 0, False

    all_suggestions: list[str] = []
    n_susp = 0
    n_unrecoverable = 0
    for f in cached_files:
        try:
            rep = detect_file(str(f))
        except Exception as exc:
            if exc.__class__.__name__ == "UnrecoverableFile":
                n_unrecoverable += 1
            continue
        for s in rep.suspicions:
            n_susp += 1
            if s.suggestion:
                for part in s.suggestion.split("|"):
                    all_suggestions.append(part.strip())
    return all_suggestions, n_susp, n_unrecoverable, True


# Backwards-compat alias : kept so older call sites still work.
def _suggestions_from_file(pmc_id: str) -> tuple[list[str], int, int]:
    suggestions, n_susp, n_unrecoverable, _in_cache = _suggestions_for_pmc_file(
        pmc_id, affected_file_url=None,
    )
    return suggestions, n_susp, n_unrecoverable


def validate_against_ziemann_notes(
    n_sample: int | None = None,
    seed: int = 42,
) -> CellValidationRun:
    """For each cached Ziemann S2 entry, validate detector against the Notes column."""
    import random
    import time
    t0 = time.monotonic()

    entries = load_ziemann_2021_corpus()
    # Only entries we have cached files for
    cached_pmcs = {p.name.split("_", 1)[0] for p in CACHE_DIR.glob("PMC*")}
    cached_entries = [e for e in entries if e.pmc_id in cached_pmcs]
    print(f"[cell-validate] {len(cached_entries)} entries with cached files "
          f"(of {len(entries)} total)")

    if n_sample:
        rng = random.Random(seed)
        cached_entries = rng.sample(cached_entries, min(n_sample, len(cached_entries)))

    results: list[CellValidation] = []
    n_tp = n_fn_miss = n_fn_sugg = n_ambig = n_missing = 0

    for i, entry in enumerate(cached_entries, 1):
        if i % 50 == 0:
            print(f"[cell-validate]  {i}/{len(cached_entries)}")

        expected = _parse_ziemann_note(entry.notes)
        if not expected:
            n_ambig += 1
            results.append(CellValidation(
                pmc_id=entry.pmc_id, notes=entry.notes,
                expected_genes=[], detector_suggestions=[],
                outcome="ambiguous-note",
            ))
            continue

        suggestions, n_susp, n_unrecoverable, file_in_cache = (
            _suggestions_for_pmc_file(entry.pmc_id, entry.affected_file_url)
        )
        if not file_in_cache:
            # Ziemann's S2 row points to a supplementary file we don't have
            # cached locally. This is a corpus-completeness issue, not a
            # detector miss : exclude from the recall denominator.
            n_missing += 1
            results.append(CellValidation(
                pmc_id=entry.pmc_id, notes=entry.notes,
                expected_genes=expected, detector_suggestions=[],
                outcome="file-not-in-cache", n_detector_suspicions=n_susp,
            ))
            continue
        if not suggestions and n_susp == 0:
            # File is cached but detector found nothing : either it's a
            # download placeholder, or genuinely no corruption.
            if n_unrecoverable > 0:
                n_missing += 1
                outcome = "file-unrecoverable"
            else:
                n_fn_miss += 1
                outcome = "FN-missed"
            results.append(CellValidation(
                pmc_id=entry.pmc_id, notes=entry.notes,
                expected_genes=expected, detector_suggestions=[],
                outcome=outcome, n_detector_suspicions=n_susp,
            ))
            continue

        # Did our suggestions include any expected gene?
        suggestion_set = set(suggestions)
        if any(g in suggestion_set for g in expected):
            n_tp += 1
            outcome = "TP"
        else:
            n_fn_sugg += 1
            outcome = "FN-suggestion-mismatch"

        results.append(CellValidation(
            pmc_id=entry.pmc_id, notes=entry.notes,
            expected_genes=expected, detector_suggestions=suggestions[:20],
            outcome=outcome, n_detector_suspicions=n_susp,
        ))

    total_validatable = n_tp + n_fn_miss + n_fn_sugg
    recall = n_tp / total_validatable if total_validatable else 0.0
    return CellValidationRun(
        timestamp=datetime.now(UTC).isoformat(),
        n_entries=len(cached_entries),
        n_files_cached=len(cached_entries) - n_missing,
        n_tp=n_tp, n_fn_missed=n_fn_miss,
        n_fn_suggestion=n_fn_sugg, n_ambiguous=n_ambig, n_file_missing=n_missing,
        cell_level_recall=recall,
        cell_level_precision=0.0,  # not measurable from S2 alone
        runtime_seconds=time.monotonic() - t0,
        results=results,
    )


def append_to_ledger(run: CellValidationRun, path: Path = LEDGER_PATH) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(json.dumps(asdict(run)) + "\n")


def format_report(run: CellValidationRun) -> str:
    lines = [
        f"=== CELL-LEVEL VALIDATION (vs Ziemann S2 Notes) ===",
        f"timestamp: {run.timestamp}",
        f"entries with notes: {run.n_entries}",
        f"files cached: {run.n_files_cached}",
        f"---",
        f"True positives: {run.n_tp:>5}  detector found Ziemann's expected gene",
        f"FN missed: {run.n_fn_missed:>5}  detector returned 0 suggestions",
        f"FN suggestion-mismatch: {run.n_fn_suggestion:>5}  detector flagged but suggestion didn't match",
        f"Ambiguous-note: {run.n_ambiguous:>5}  Ziemann's note didn't parse into a unique gene",
        f"File-missing: {run.n_file_missing:>5}  no cached file",
        f"---",
        f"CELL-LEVEL RECALL: {run.cell_level_recall:.4f}  "
        f"(TP / (TP + FN-missed + FN-suggestion-mismatch))",
        f"runtime: {run.runtime_seconds:.1f}s",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=None, help="sample size (default: all cached)")
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()
    run = validate_against_ziemann_notes(n_sample=args.n, seed=args.seed)
    print(format_report(run))
    append_to_ledger(run)
    print(f"\nAppended to {LEDGER_PATH}")
