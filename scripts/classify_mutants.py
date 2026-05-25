"""Manual classification of cosmic-ray surviving mutants per
Schuler & Zeller (2013).

For each surviving mutant, show:
  - The mutation diff
  - Coverage profile (does the mutated line get exercised by the test
    suite at all?)
  - A suggested classification: equivalent / test-gap / live-bug

The reviewer enters a verdict. Outputs `results/mutant_classifications.jsonl`.

Recomputes the adjusted kill rate excluding confirmed equivalents:

    adjusted_kill_rate = killed / (total - n_equivalent)

Per Schuler & Zeller's reported ~45% equivalent rate in similar codebases,
our 39 surviving / 80 total may have ~18 equivalents → adjusted ~65%.

Run:
    python scripts/classify_mutants.py
    python scripts/classify_mutants.py --batch      # print without prompts
"""
from __future__ import annotations

import argparse
import json
import sqlite3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DB = ROOT / "results/cosmic-ray-detector.sqlite"
OUT = ROOT / "results/mutant_classifications.jsonl"


def _suggest_classification(diff: str) -> str:
    """Quick heuristic for a suggested verdict : reviewer overrides."""
    d = diff.lower()
    if "==" in d and "!=" in d:
        return "boundary"
    if "and " in d and "or " in d:
        return "boundary"
    if "true" in d or "false" in d:
        return "boundary"
    if "remove" in d and "decorator" in d:
        return "test-gap"
    if "return " in d:
        return "live-bug"
    return "needs-review"


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--batch", action="store_true",
                   help="Print survivors without interactive prompting.")
    args = p.parse_args()

    if not DB.exists():
        print(f"FATAL: cosmic-ray DB missing at {DB}", file=sys.stderr)
        return 1

    con = sqlite3.connect(DB)
    cur = con.cursor()
    cur.execute(
        "SELECT job_id, diff FROM work_results "
        "WHERE test_outcome = 'SURVIVED'"
    )
    survivors = list(cur.fetchall())
    print(f"Cosmic-ray survivors: {len(survivors)}")

    # Load prior decisions if any
    already_seen: set[str] = set()
    if OUT.exists():
        for line in OUT.read_text().splitlines():
            if line:
                already_seen.add(json.loads(line)["job_id"])

    new_count = 0
    n_equiv = 0
    n_gap = 0
    n_bug = 0
    for job_id, diff in survivors:
        if job_id in already_seen:
            continue
        suggestion = _suggest_classification(diff)
        print()
        print("=" * 78)
        print(f"job_id: {job_id}")
        print(f"suggested: {suggestion}")
        print()
        # Show only the meaningful lines of the diff
        for line in diff.split("\n"):
            if line.startswith(("-    ", "+    ", "@@")):
                print(f"  {line.rstrip()}")
        if args.batch:
            continue
        try:
            v = input("  verdict (e=equivalent / g=test-gap / b=live-bug / s=skip / q=quit): ").strip().lower()
        except (EOFError, KeyboardInterrupt):
            print("\n(interrupted)")
            break
        if v == "q":
            break
        if v == "s":
            continue
        if v == "e":
            label = "equivalent"
            n_equiv += 1
        elif v == "g":
            label = "test-gap"
            n_gap += 1
        elif v == "b":
            label = "live-bug"
            n_bug += 1
        else:
            print("  invalid; skipping")
            continue
        note = input("  one-line note (optional): ").strip()
        with OUT.open("a") as fh:
            fh.write(json.dumps({
                "job_id": job_id,
                "verdict": label,
                "note": note,
                "suggestion": suggestion,
            }) + "\n")
        new_count += 1

    print()
    if new_count:
        print(f"Recorded {new_count} new decisions "
              f"(equivalent={n_equiv}, gap={n_gap}, bug={n_bug}).")

    # Recompute adjusted kill rate
    total_jobs = cur.execute("SELECT count(*) FROM work_items").fetchone()[0]
    killed = cur.execute(
        "SELECT count(*) FROM work_results WHERE test_outcome = 'KILLED'"
    ).fetchone()[0]
    n_total_equivalent = 0
    if OUT.exists():
        for line in OUT.read_text().splitlines():
            if line and json.loads(line).get("verdict") == "equivalent":
                n_total_equivalent += 1
    raw_kill = killed / total_jobs if total_jobs else 0.0
    adjusted_kill = (
        killed / (total_jobs - n_total_equivalent)
        if (total_jobs - n_total_equivalent) > 0 else 0.0
    )
    print()
    print(f"Mutation testing summary:")
    print(f"  total mutants: {total_jobs}")
    print(f"  killed: {killed}")
    print(f"  surviving: {len(survivors)}")
    print(f"  classified equivalent: {n_total_equivalent}")
    print(f"  raw kill rate: {raw_kill:.3f}")
    print(f"  adjusted kill rate: {adjusted_kill:.3f}")
    print()
    print(f"Decisions ledger: {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
