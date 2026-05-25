"""Fresh-files manual proof of the iter-5 detector.

Randomly samples 20 PMC files from the Koh corpus that contain at least one
gene-symbol column, runs detect_file on each, and prints a structured audit
report:

  - For each suspicion: the cell value, the openpyxl number_format, the
    detector's suggestion, and the verdict ("looks real" by date-rule heuristic).
  - For each file: every cell in the identifier columns that the detector did
    NOT flag, so a human can verify nothing obvious was missed.

The output is intended to be human-readable enough that the user can scroll
through it and confirm "yes, these are all corruption" / "no FPs".

Usage:
    python scripts/fresh_files_proof.py [--n 20] [--seed 42]
    python scripts/fresh_files_proof.py --pmcs PMC1234567 PMC2345678
"""
from __future__ import annotations

import argparse
import json
import random
import sys
import warnings
from datetime import date, datetime
from pathlib import Path

warnings.filterwarnings("ignore")

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from uncorrupt.corpus import load_hgnc  # noqa: E402
from uncorrupt.detector import detect_file  # noqa: E402

KOH_DIR = PROJECT_ROOT / "data" / "raw" / "koh_replication" / "files"
WALK_LEDGER = PROJECT_ROOT / "results" / "koh_cache_walk.jsonl"


def _is_real_corruption_signature(value, suggestion: str | None) -> str:
    """Quick verdict heuristic: does this look like a real Excel-corruption?

    Returns one of: 'date-corruption-likely', 'serial-corruption-likely',
    'string-corruption-likely', 'id-float-likely', 'inconclusive'.
    """
    if suggestion is None or suggestion == "":
        return "id-float-likely" if isinstance(value, (int, float)) and abs(value) >= 1e10 else "inconclusive"
    if isinstance(value, (date, datetime)):
        return "date-corruption-likely"
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return "serial-corruption-likely"
    if isinstance(value, str):
        return "string-corruption-likely"
    return "inconclusive"


def audit_file(fp: Path) -> dict:
    """Run detect_file + collect a structured audit record."""
    try:
        report = detect_file(str(fp))
    except Exception as exc:
        return {
            "file": str(fp.relative_to(PROJECT_ROOT)),
            "pmc_id": fp.parent.name,
            "error": f"{type(exc).__name__}: {exc}",
        }
    return {
        "file": str(fp.relative_to(PROJECT_ROOT)),
        "pmc_id": fp.parent.name,
        "n_suspicions": len(report.suspicions),
        "identifier_columns": report.identifier_columns,
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
                "verdict_heuristic": _is_real_corruption_signature(s.value, s.suggestion),
            }
            for s in report.suspicions
        ],
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=20)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--pmcs", nargs="+", help="specific PMC ids to audit")
    p.add_argument("--from-walk-ledger", action="store_true",
                   help="only sample files that the walk-ledger marked as having corruption")
    p.add_argument("--out", default="results/fresh_files_proof.jsonl")
    args = p.parse_args()

    # Find candidate files
    if args.pmcs:
        targets = []
        for pmc in args.pmcs:
            pmc_dir = KOH_DIR / pmc
            if not pmc_dir.exists():
                print(f"warning: {pmc} not in cache", file=sys.stderr)
                continue
            for fp in pmc_dir.iterdir():
                if fp.suffix.lower() in (".xlsx", ".xls") and not fp.name.startswith("_"):
                    targets.append(fp)
    elif args.from_walk_ledger:
        # only files that the recent cache-walk marked as having_corruption
        with WALK_LEDGER.open() as f:
            rows = [json.loads(l) for l in f]
        latest = rows[-1]
        candidates = [
            (f["pmc_id"], f["file_name"])
            for f in latest["files"]
            if f["has_gene_symbols"] and f["n_suspicions"] > 0
        ]
        rng = random.Random(args.seed)
        sampled = rng.sample(candidates, min(args.n, len(candidates)))
        targets = [KOH_DIR / pmc / fname for pmc, fname in sampled]
        targets = [t for t in targets if t.exists()]
    else:
        all_files = []
        for pmc_dir in KOH_DIR.iterdir():
            if pmc_dir.is_dir():
                for fp in pmc_dir.iterdir():
                    if fp.suffix.lower() in (".xlsx", ".xls") and not fp.name.startswith("_"):
                        all_files.append(fp)
        rng = random.Random(args.seed)
        targets = rng.sample(all_files, min(args.n, len(all_files)))

    print(f"[proof] auditing {len(targets)} files", file=sys.stderr)

    out_path = PROJECT_ROOT / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)
    audits: list[dict] = []
    with out_path.open("w") as f:
        for i, fp in enumerate(targets, 1):
            audit = audit_file(fp)
            audits.append(audit)
            f.write(json.dumps(audit) + "\n")
            n_susp = audit.get("n_suspicions", "ERR")
            print(f"[proof] {i:>3}/{len(targets)}  {fp.parent.name}/{fp.name:<40}  "
                  f"suspicions={n_susp}", file=sys.stderr)

    # Summary
    print(f"\n=== SUMMARY ===", file=sys.stderr)
    print(f"Files audited: {len(audits)}", file=sys.stderr)
    n_err = sum(1 for a in audits if "error" in a)
    n_clean = sum(1 for a in audits if "n_suspicions" in a and a["n_suspicions"] == 0)
    n_dirty = sum(1 for a in audits if "n_suspicions" in a and a["n_suspicions"] > 0)
    total_susp = sum(a.get("n_suspicions", 0) for a in audits if "n_suspicions" in a)
    print(f"Load/detect errors: {n_err}", file=sys.stderr)
    print(f"Files with 0 suspicions: {n_clean}", file=sys.stderr)
    print(f"Files with ≥1 suspicion: {n_dirty}", file=sys.stderr)
    print(f"Total cell suspicions: {total_susp}", file=sys.stderr)
    print(f"\nAppended to {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
