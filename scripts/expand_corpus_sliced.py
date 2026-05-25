"""Sliced EPMC corpus expansion : splits PMC IDs of a single journal across
N workers. Used to parallelize the remaining work in Scientific Reports
and eLife after the initial 4-way split completed the other journals.

Each worker calls run_koh_replication with a single journal name and
(pmc_id_offset, pmc_id_limit) : the script slices that journal's IDs
to the requested range and processes only those.

Run:
    python scripts/expand_corpus_sliced.py \
        --journal "Scientific Reports" \
        --offset 0 --limit 2500 \
        --worker-id sr0 \
        --rate 1.0
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from uncorrupt.koh_replication import run_koh_replication, append_to_ledger


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--journal", type=str, required=True)
    p.add_argument("--offset", type=int, default=0)
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--worker-id", type=str, required=True)
    p.add_argument("--rate", type=float, default=1.0)
    p.add_argument("--year-min", type=int, default=2019)
    p.add_argument("--year-max", type=int, default=2026)
    args = p.parse_args()

    print(f"[sliced-{args.worker_id}] journal='{args.journal}' "
          f"slice=[{args.offset}:{(args.offset + args.limit) if args.limit else 'end'}]")
    run = run_koh_replication(
        year_min=args.year_min,
        year_max=args.year_max,
        journals=[args.journal],
        rate_limit_per_second=args.rate,
        checkpoint_path=ROOT / f"results/koh_cache_walk_sliced_{args.worker_id}.checkpoint.json",
        pmc_id_offset=args.offset,
        pmc_id_limit=args.limit,
    )
    append_to_ledger(run)
    print(f"[sliced-{args.worker_id}] DONE files={run.n_supplementary_files} "
          f"corruption={run.n_with_corruption} runtime={run.runtime_seconds:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
