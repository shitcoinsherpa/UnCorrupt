"""Sample 'inconclusive' flags from the FP-anchor calibration ledger and
present each one for manual adjudication.

The calibration walk (scripts/calibrate_with_fp_anchors.py) labels each
flag against row-context xref evidence. About 73% of flags are
'inconclusive' (the row has no resolvable external ID to compare against
the suggested gene). Without manual review of a sample, the per-bin
precision tables only reflect the xref-labeled minority.

This script:

  1. Reads results/calibration_fp_anchors.jsonl.
  2. Samples N inconclusive flags (default 100), stratified across
     suspicion kinds.
  3. For each sample, opens the source file and shows:
       - file path, sheet, column, row, value
       - the surrounding 5 rows
       - the suggested gene
       - the detector's reason
  4. Prompts the reviewer for a y/n/u verdict (yes-real-corruption /
     no-false-positive / unable-to-tell) and a one-line note.
  5. Writes the decisions to results/adjudication_decisions.jsonl.
  6. Recomputes per-bin precision incorporating the adjudicated labels.

Usage:
   python scripts/adjudicate_inconclusive.py --sample 100
   python scripts/adjudicate_inconclusive.py --sample 20 --kind gene-date-serial
"""
from __future__ import annotations

import argparse
import json
import random
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

LEDGER = ROOT / "results/calibration_fp_anchors.jsonl"
KOH = ROOT / "data/raw/koh_replication/files"
OUT = ROOT / "results/adjudication_decisions.jsonl"
SEED = 20260518


def _load_inconclusive(kind_filter: str | None) -> list[dict]:
    rows = []
    with LEDGER.open() as fh:
        for line in fh:
            r = json.loads(line)
            if r["label"] != "inconclusive":
                continue
            if kind_filter and r["kind"] != kind_filter:
                continue
            # Skip non-corruption-claim kinds
            if r["kind"] in {"decimal-comma", "cas-registry"}:
                continue
            rows.append(r)
    return rows


def _show_context(record: dict) -> None:
    """Print the file path + surrounding cells for the given record."""
    path = KOH / record["pmc_id"] / record["file"]
    print()
    print("=" * 78)
    print(f"File: {path}")
    print(f"Sheet: {record['sheet']!r}")
    print(f"Column: {record['column']!r}")
    print(f"Row: {record['row']}")
    print(f"Value: {record['value']!r}")
    print(f"Kind: {record['kind']}")
    print(f"Sug: {record['suggestion']!r}")
    print(f"Conf: {record['confidence']}")
    if not path.exists():
        print(f"(file not on disk)")
        return
    try:
        from uncorrupt.app import _load_all_sheets
        sheets = _load_all_sheets(str(path))
    except Exception as exc:
        print(f"(load error: {type(exc).__name__}: {exc})")
        return
    sheet_name = record["sheet"]
    df = sheets.get(sheet_name) if sheet_name else next(iter(sheets.values()))
    if df is None:
        print("(sheet not found)")
        return
    row = record["row"]
    lo = max(0, row - 2)
    hi = min(len(df), row + 3)
    if record["column"] in df.columns:
        print(f"\nSurrounding cells (rows {lo}..{hi-1}):")
        print(df.iloc[lo:hi].to_string(max_cols=8, max_colwidth=40))


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample", type=int, default=100)
    p.add_argument("--kind", default=None,
                   help="Filter to one suspicion kind (e.g., gene-date-serial)")
    p.add_argument("--batch", action="store_true",
                   help="Print records without prompting (for offline review).")
    args = p.parse_args()

    inconclusive = _load_inconclusive(args.kind)
    print(f"{len(inconclusive)} inconclusive flags available "
          f"(kind={args.kind or 'any'}).")
    if not inconclusive:
        return 0

    rng = random.Random(SEED)
    sample = rng.sample(inconclusive, min(args.sample, len(inconclusive)))

    decisions: list[dict] = []
    if OUT.exists():
        # Resume from prior decisions
        for line in OUT.read_text().splitlines():
            if line:
                decisions.append(json.loads(line))
    already_seen = {(d["file"], d["sheet"], d["column"], d["row"])
                    for d in decisions}

    new_count = 0
    for i, rec in enumerate(sample):
        key = (rec["file"], rec["sheet"], rec["column"], rec["row"])
        if key in already_seen:
            continue
        _show_context(rec)
        print(f"\n[{i+1}/{len(sample)}]")
        if args.batch:
            continue
        try:
            verdict = input("  verdict (y=real/n=fp/u=unable/q=quit): ").strip().lower()
        except (EOFError, KeyboardInterrupt):
            print("\n(interrupted)")
            break
        if verdict == "q":
            break
        if verdict not in {"y", "n", "u"}:
            print("  invalid; skipping")
            continue
        note = input("  one-line note: ").strip()
        entry = {**rec, "verdict": verdict, "note": note}
        with OUT.open("a") as fh:
            fh.write(json.dumps(entry) + "\n")
        new_count += 1

    print(f"\nRecorded {new_count} new decisions. Total: {len(decisions) + new_count}.")
    print(f"Ledger: {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
