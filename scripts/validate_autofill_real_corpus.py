"""Validate the autofill-sequence detector against the real Ziemann 2021 corpus.

This is *real-data validation*, not a synthetic unit test. The autofill flag
fires on a column when three or more row-adjacent cells form a monotonic
date series AND each maps to a gene-family member (Ziemann 2021 PMC6330011
documents this Excel propagation mode).

We run detect_file across every supplementary file in:
  data/raw/ziemann_2021_corpus/files/

and report:
  - autofill_files       : files with at least one autofill-sequence flag
  - total_autofill_flags : total runs flagged (a file can have multiple)
  - per_file_top_5       : the five largest runs across the corpus
  - sanity_count         : total gene-date suspicions (sanity: the detector
                           is still seeing the older corruption modes)

Run:
    python scripts/validate_autofill_real_corpus.py
"""
from __future__ import annotations

import json
import sys
import time
import traceback
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from uncorrupt.detector import detect_file  # noqa: E402

CORPUS = ROOT / "data/raw/ziemann_2021_corpus/files"
OUT = ROOT / "results/autofill_real_corpus_validation.jsonl"


def main() -> int:
    files = sorted(p for p in CORPUS.iterdir() if p.suffix.lower() in (".xls", ".xlsx"))
    OUT.parent.mkdir(parents=True, exist_ok=True)

    autofill_files = 0
    total_autofill_flags = 0
    sanity_gene_date = 0
    largest_runs: list[tuple[int, str, str, int, str]] = []  # (len, file, col, row, reason)
    failed: list[tuple[str, str]] = []

    t0 = time.time()
    with OUT.open("w") as fh:
        for i, f in enumerate(files):
            try:
                report = detect_file(str(f), row_context_boost=False)
            except Exception as exc:
                failed.append((f.name, f"{type(exc).__name__}: {exc}"))
                fh.write(json.dumps({"file": f.name, "error": str(exc)}) + "\n")
                continue

            autofill = [s for s in report.suspicions if s.kind == "autofill-sequence"]
            gene_dates = [s for s in report.suspicions if s.kind == "gene-date"]
            sanity_gene_date += len(gene_dates)

            if autofill:
                autofill_files += 1
                total_autofill_flags += len(autofill)
                for s in autofill:
                    # The run length is embedded in the value field
                    # as "rows X..Y (N cells)" : parse N.
                    n_cells = 0
                    val = str(s.value)
                    if "(" in val and "cells)" in val:
                        try:
                            n_cells = int(val.split("(")[1].split()[0])
                        except (ValueError, IndexError):
                            pass
                    largest_runs.append((n_cells, f.name, s.column, s.row, s.reason))
                    fh.write(json.dumps({
                        "file": f.name,
                        "sheet": s.sheet,
                        "column": s.column,
                        "row": s.row,
                        "value": val,
                        "confidence": s.confidence,
                        "reason": s.reason,
                    }) + "\n")

            if (i + 1) % 50 == 0:
                elapsed = time.time() - t0
                print(f"[{i+1}/{len(files)}] {autofill_files} autofill files, "
                      f"{total_autofill_flags} flags, {sanity_gene_date} gene-date "
                      f"sanity hits ({elapsed:.1f}s)", flush=True)

    largest_runs.sort(reverse=True)
    print()
    print(f"Files scanned: {len(files)}")
    print(f"Failed to load: {len(failed)}")
    print(f"Files with autofill: {autofill_files} ({100*autofill_files/max(1,len(files)):.1f}%)")
    print(f"Total autofill flags: {total_autofill_flags}")
    print(f"Sanity gene-date: {sanity_gene_date}")
    print()
    print("Top 5 longest autofill runs:")
    for length, fname, col, row, reason in largest_runs[:5]:
        print(f"  {length} cells | {fname} | col={col!r} row={row}")
        print(f"     {reason[:140]}")
    if failed:
        print()
        print(f"Failures ({len(failed)}):")
        for name, err in failed[:5]:
            print(f"  {name}: {err}")
    print()
    print(f"Detail ledger: {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
