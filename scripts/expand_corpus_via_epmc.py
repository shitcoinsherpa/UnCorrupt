"""Expand the validation corpus by re-running the Koh-style scrape against a
broader journal set + a wider time window.

The original Koh 2022 list (11 journals) gave us 7,945 cached files. To push
labeled-cell n from 1,630 toward 6-10K (and Wilson lower bound 0.998 → 0.9995)
we widen the journal set and back-fill 2019-2021.

Added journals beyond the Koh 11:
  - Cell                                 (high-impact, supplementary-rich)
  - Cell Reports                         (open-access, weekly)
  - Bioinformatics                       (computational papers, many tables)
  - Briefings in Bioinformatics          (Ziemann's own outlet)
  - eLife                                (open-access)
  - Scientific Reports                   (Koh 2022 itself was published here)

Time window: 2019-01 through 2026-05 (covers Koh's June-2022 sample + 3 years
each side).
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from uncorrupt.koh_replication import run_koh_replication, KOH_JOURNALS  # noqa: E402

# Wider journal list : Koh 11 + 6 more
EXTENDED_JOURNALS = KOH_JOURNALS + [
    "Cell",
    "Cell Reports",
    "Bioinformatics",
    "Briefings in Bioinformatics",
    "eLife",
    "Scientific Reports",
]

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--year-min", type=int, default=2019)
    p.add_argument("--year-max", type=int, default=2026)
    p.add_argument("--rate", type=float, default=2.5,
                   help="requests/sec for NCBI esearch + EuropePMC")
    args = p.parse_args()

    run = run_koh_replication(
        year_min=args.year_min,
        year_max=args.year_max,
        journals=EXTENDED_JOURNALS,
        rate_limit_per_second=args.rate,
        checkpoint_path=ROOT / "results/koh_cache_walk_extended.checkpoint.json",
    )
    print(f"[epmc-expand] DONE: {run.n_supplementary_files} files scanned, "
          f"{run.n_with_gene_symbols} with gene symbols, "
          f"{run.n_with_corruption} with corruption "
          f"in {run.runtime_seconds:.0f}s")
