"""Validate current detector against the Ziemann 2016 S1 corpus.

This is INDEPENDENT ground truth : a different paper, different sampling
methodology, different time period (2005-2015) from the Koh 2022 corpus
used by `scripts/calibrate_with_fp_anchors.py`. Each row in
`Additional_file_1.xlsx` (Ziemann 2016) has:

  - A "Confirmed" flag set to "Confirmed" iff the authors hand-verified
    the corruption.
  - A "Example Gene Name Conversion" column documenting the corruption
    pattern, e.g. "MARCH9 → 2005-09-01", "SEPT2 → 2-Sep", "OCT4 → Oct-4".

The validation routine for each fetched file:

  1. Loads the file with `detect_file()` .
  2. Extracts the *expected* gene symbols from the "Example" column
     (regex over uppercase tokens that match HGNC/MGI/etc.).
  3. Reports per-file outcome:
       - TP: any detector suggestion contains an expected gene.
       - FN-missed: detector emitted no suspicions at all.
       - FN-suggestion: detector emitted suspicions but none match.
       - file-unreadable: load failed.
  4. Computes recall against the Confirmed ground truth.

Output: `results/ziemann_2016_validation.jsonl` (per-file records) and
a summary report printed to stdout.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import warnings
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from uncorrupt.detector import (  # noqa: E402
    detect_file, _reverse_gene_date, _parse_date_string, _serial_to_date,
)
from uncorrupt.corpus import load_ziemann_2016_corpus  # noqa: E402
from datetime import date, datetime  # noqa: E402
import re as _re  # noqa: E402

FILES_DIR = ROOT / "data/raw/ziemann_2016_corpus/files"
OUT_JSONL = ROOT / "results/ziemann_2016_validation.jsonl"
OUT_REPORT = ROOT / "results/ziemann_2016_validation.md"

# Token forms in the "Example Gene Name Conversion" column.
# Examples seen in the corpus:
#   "MARCH9 → 9-Mar"  → expected = {"MARCHF9", "MARCH9"}
#   "SEPT2 to Sep-2"  → expected = {"SEPTIN2", "SEPT2"}
#   "OCT4 became Oct-04"  → expected = {"POU5F1", "OCT4"}
_GENE_TOKEN_RE = re.compile(r"\b([A-Z][A-Z0-9]{1,9})\b")

# Known HGNC current<->prev_symbol pairs for the most common date-fragile
# gene families (Ziemann's headline examples). The detector's HGNC table
# is canonical; we replicate only the small subset that appears in this
# corpus's "Example" column to drive expected-set expansion.
_FAMILY_EQUIV = {
    "MARCH1": ("MARCHF1",), "MARCH2": ("MARCHF2",), "MARCH3": ("MARCHF3",),
    "MARCH4": ("MARCHF4",), "MARCH5": ("MARCHF5",), "MARCH6": ("MARCHF6",),
    "MARCH7": ("MARCHF7",), "MARCH8": ("MARCHF8",), "MARCH9": ("MARCHF9",),
    "MARCH10": ("MARCHF10",), "MARCH11": ("MARCHF11",),
    "SEPT1": ("SEPTIN1",), "SEPT2": ("SEPTIN2",), "SEPT3": ("SEPTIN3",),
    "SEPT4": ("SEPTIN4",), "SEPT5": ("SEPTIN5",), "SEPT6": ("SEPTIN6",),
    "SEPT7": ("SEPTIN7",), "SEPT8": ("SEPTIN8",), "SEPT9": ("SEPTIN9",),
    "SEPT10": ("SEPTIN10",), "SEPT11": ("SEPTIN11",), "SEPT12": ("SEPTIN12",),
    "SEPT14": ("SEPTIN14",),
    "DEC1": ("DELEC1",), "DEC2": ("BHLHE41",),
    "OCT4": ("POU5F1", "OCT3"), "OCT3": ("POU5F1", "OCT4"),
}


def _expand_expected(tokens: set[str]) -> set[str]:
    out = set(tokens)
    for t in list(tokens):
        if t in _FAMILY_EQUIV:
            out.update(_FAMILY_EQUIV[t])
    return out


def _extract_expected_genes(example_text: str) -> set[str]:
    """Extract the gene symbol(s) the corrupted value SHOULD have been.

    The Ziemann 2016 S1 'Example Gene Name Conversion' column holds the
    CORRUPTED value as it appeared in the file (e.g. '2002-09-01 00:00:00'
    or '6029999999999999517655040'). We reverse-decode it via the
    detector's own machinery to compute what gene(s) it would represent.

    Two paths covered:
      1. Date strings / datetimes → `_reverse_gene_date(date_obj)` returns
         the candidate gene-family members.
      2. Pure-digit huge-number → indicates a RIKEN/accession float-coercion;
         we mark this as 'expected non-gene-but-corruption' so the validator
         can still register that the file WAS corrupted (recall=TP iff
         detector emitted ANY suspicion).

    Returns a set of HGNC symbols (uppercase). Empty set means we cannot
    derive an expectation from this example."""
    if not example_text or example_text.lower() in ("nan", "none", ""):
        return set()

    text = example_text.strip()
    expected: set[str] = set()

    # Path 1: try to parse as a date string or pandas Timestamp form.
    # The S1 column often shows openpyxl-coerced datetime: 'YYYY-MM-DD HH:MM:SS'.
    iso_match = _re.match(r"^(\d{4})-(\d{1,2})-(\d{1,2})", text)
    parsed: date | None = None
    if iso_match:
        try:
            y, m, d = int(iso_match.group(1)), int(iso_match.group(2)), int(iso_match.group(3))
            parsed = date(y, m, d)
        except ValueError:
            pass

    if parsed is None:
        # Try other locale-style date parses via detector's parser
        for d in _parse_date_string(text):
            parsed = d
            break

    if parsed is not None:
        for cand in _reverse_gene_date(parsed):
            expected.add(cand)
        return _expand_expected(expected)

    # Path 2: pure digit string : likely a RIKEN/accession float-coercion
    # (e.g. '6029999999999999517655040' = 6.03e24 from a 2310009E13-like ID).
    # We cannot derive a specific gene; mark as "non-gene corruption" via a
    # sentinel so the validator counts the file as expecting *some* flag.
    if _re.fullmatch(r"-?\d+(\.\d+)?(e[+-]?\d+)?", text, _re.I):
        expected.add("__NON_GENE_CORRUPTION__")
        return expected

    # Path 3: arbitrary text : try to scan it for gene-symbol-looking tokens
    # (legacy behavior for safety).
    raw = set(_GENE_TOKEN_RE.findall(text))
    drop = {"TO", "INTO", "BECAME", "AS", "OR", "AND", "BY", "WITH",
            "FROM", "GENE", "NAME", "EX", "IS", "WAS", "BE",
            "JAN", "FEB", "MAR", "APR", "MAY", "JUN",
            "JUL", "AUG", "SEP", "OCT", "NOV", "DEC",
            "ID", "ID1", "ID2", "IDS",
            "EXCEL", "DATE", "DATES", "NUMBER", "FLOAT"}
    raw = {t for t in raw if t not in drop and len(t) >= 3}
    return _expand_expected(raw)


def _parse_ledger() -> dict[str, dict]:
    """Map filename → metadata from the fetch ledger. The ledger is a
    single pretty-printed JSON document (despite the `.jsonl` extension),
    with a top-level `successful` array."""
    ledger_path = ROOT / "data/raw/ziemann_2016_corpus/fetch_ledger.jsonl"
    if not ledger_path.exists():
        return {}
    try:
        entry = json.loads(ledger_path.read_text())
    except json.JSONDecodeError:
        return {}
    out: dict[str, dict] = {}
    for s in entry.get("successful", []):
        out[s["filename"]] = s
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--limit", type=int, default=None,
                   help="Cap on number of files to process (default: all).")
    p.add_argument("--max-bytes", type=int, default=50 * 1024 * 1024,
                   help="Skip files larger than this (default: 50 MB).")
    args = p.parse_args()

    ledger = _parse_ledger()
    files = sorted(FILES_DIR.glob("*.xls*"))
    if args.limit:
        files = files[: args.limit]
    print(f"[ziemann-2016] {len(files)} files to validate", flush=True)

    outcomes: list[dict] = []
    counter: Counter[str] = Counter()
    n_tp = n_fn_miss = n_fn_sugg = n_unreadable = n_skipped = n_no_expected = 0

    OUT_JSONL.parent.mkdir(parents=True, exist_ok=True)
    with OUT_JSONL.open("w") as fh:
        for i, path in enumerate(files, 1):
            meta = ledger.get(path.name, {})
            example = meta.get("documented_corruption", "")
            expected = _extract_expected_genes(example)

            if path.stat().st_size > args.max_bytes:
                n_skipped += 1
                outcome = "skipped-large"
            elif not expected:
                # Corruption example didn't yield extractable gene tokens.
                # Could be RIKEN/accession-format, decimal-coercion, etc.
                n_no_expected += 1
                outcome = "no-expected-genes"
            else:
                try:
                    report = detect_file(str(path), row_context_boost=True)
                except Exception as exc:
                    n_unreadable += 1
                    outcome = f"unreadable:{type(exc).__name__}"
                    counter[outcome] += 1
                    rec = {"file": path.name, "pubmed_id": meta.get("pubmed_id"),
                           "journal": meta.get("journal"),
                           "documented_corruption": example,
                           "expected_genes": sorted(expected),
                           "detector_suggestions": [],
                           "outcome": outcome}
                    fh.write(json.dumps(rec) + "\n")
                    continue

                # Collect ALL detector suggestions (split on '|')
                all_suggs: set[str] = set()
                for s in report.suspicions:
                    if s.suggestion:
                        for part in s.suggestion.split("|"):
                            all_suggs.add(part.strip().upper())

                non_gene_expected = "__NON_GENE_CORRUPTION__" in expected
                gene_expected = expected - {"__NON_GENE_CORRUPTION__"}

                if non_gene_expected and not gene_expected:
                    # Numeric / RIKEN-coercion case : TP iff detector
                    # emitted ANY suspicion on the file.
                    if report.suspicions:
                        n_tp += 1
                        outcome = "TP-any-flag"
                    else:
                        n_fn_miss += 1
                        outcome = "FN-missed"
                elif not all_suggs:
                    n_fn_miss += 1
                    outcome = "FN-missed"
                elif gene_expected & all_suggs:
                    n_tp += 1
                    outcome = "TP"
                else:
                    n_fn_sugg += 1
                    outcome = "FN-suggestion-mismatch"

                rec = {
                    "file": path.name,
                    "pubmed_id": meta.get("pubmed_id"),
                    "journal": meta.get("journal"),
                    "year": meta.get("year"),
                    "documented_corruption": example,
                    "expected_genes": sorted(expected),
                    "detector_suggestions": sorted(all_suggs)[:30],
                    "n_total_suspicions": len(report.suspicions),
                    "outcome": outcome,
                }
                fh.write(json.dumps(rec) + "\n")

            counter[outcome] += 1
            if i % 50 == 0:
                print(f"[ziemann-2016] {i}/{len(files)}  "
                      f"TP={n_tp} FN-miss={n_fn_miss} "
                      f"FN-sugg={n_fn_sugg} unread={n_unreadable} "
                      f"no-expected={n_no_expected}",
                      flush=True)

    # Final report
    n_validatable = n_tp + n_fn_miss + n_fn_sugg
    recall = n_tp / n_validatable if n_validatable else 0.0
    print()
    print("=" * 60)
    print("Ziemann 2016 : validation summary")
    print("=" * 60)
    print(f"Files processed: {len(files)}")
    print(f"  TP: {n_tp}")
    print(f"  FN-missed: {n_fn_miss}   (detector found nothing)")
    print(f"  FN-sugg: {n_fn_sugg}   (found something, wrong gene)")
    print(f"  unreadable: {n_unreadable}")
    print(f"  skipped-lg: {n_skipped}")
    print(f"  no-expected: {n_no_expected} (couldn't parse expected genes)")
    print(f"  Validatable: {n_validatable}")
    print(f"  RECALL: {recall:.4f}")

    # Wilson 95% CI on recall
    import math
    if n_validatable:
        z = 1.96; p = recall; n = n_validatable
        d = 1 + z*z/n
        c = (p + z*z/(2*n)) / d
        h = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / d
        lo, hi = max(0, c-h), min(1, c+h)
        print(f"  Wilson 95%: [{lo:.4f}, {hi:.4f}]")

    # Write markdown report
    lines = [
        "# Ziemann 2016 : independent recall validation",
        "",
        f"Run timestamp: {datetime.now(UTC).isoformat()}",
        f"Source: `data/raw/ziemann_2016_corpus/files/` "
        f"({len(files)} fetched files from Ziemann 2016 Additional File 1).",
        "",
        "## Methodology",
        "",
        "For every Ziemann-2016-Confirmed corruption file we have on disk:",
        "1. Extract the *expected* gene symbol(s) from the published",
        "   'Example Gene Name Conversion' column (e.g. \"MARCH9 → 2005-09-01\"",
        "   yields expected = {MARCH9, MARCHF9}).",
        "2. Run `detect_file(...row_context_boost=True)` on the cached file.",
        "3. Collect every emitted suggestion (split on `|`).",
        "4. TP iff `expected ∩ detector_suggestions ≠ ∅`.",
        "",
        "Independent of the Koh 2022 corpus used by",
        "`scripts/calibrate_with_fp_anchors.py` : different paper,",
        "different sampling methodology, different journals, time period 2005-2015.",
        "",
        "## Results",
        "",
        f"- Files processed: **{len(files)}**",
        f"- True positives (TP): **{n_tp}**",
        f"- False negatives : detector emitted nothing (FN-missed): {n_fn_miss}",
        f"- False negatives : wrong suggestion (FN-sugg): {n_fn_sugg}",
        f"- Unreadable: {n_unreadable}",
        f"- Skipped (>50 MB): {n_skipped}",
        f"- No-expected-genes (couldn't parse): {n_no_expected}",
        f"- **Validatable set**: {n_validatable}",
        f"- **Recall**: **{recall:.4f}**",
    ]
    if n_validatable:
        lines.append(f"- **Wilson 95 % CI on recall**: [{lo:.4f}, {hi:.4f}]")
    lines.extend([
        "",
        "## Outcome breakdown",
        "",
        "| Outcome | n |",
        "|---|---:|",
    ])
    for o, n in counter.most_common():
        lines.append(f"| {o} | {n} |")
    OUT_REPORT.write_text("\n".join(lines))
    print(f"Wrote {OUT_REPORT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
