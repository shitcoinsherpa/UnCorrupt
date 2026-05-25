"""Inspect the detector's flags on the 'clean with genes' pool.

For each flag found by null_result_clean_pool.py, fetch the actual cell
value, the column header, and the surrounding rows. This lets a human
adjudicate: is it (a) detector regression, (b) cache-classifier miss
(real corruption never previously caught), or (c) genuine false positive?
"""
from __future__ import annotations

import json
import random
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from uncorrupt.detector import detect_file  # noqa: E402

WALK_LEDGER = ROOT / "results/koh_cache_walk.jsonl"
KOH = ROOT / "data/raw/koh_replication/files"
SAMPLE_N = 100
SEED = 20260518


def main() -> int:
    with WALK_LEDGER.open() as fh:
        walk = json.load(fh)
    clean_with_genes = [
        f for f in walk.get("files", [])
        if f.get("has_gene_symbols") and f.get("n_suspicions", 0) == 0
    ]
    rng = random.Random(SEED)
    sample = rng.sample(clean_with_genes, min(SAMPLE_N, len(clean_with_genes)))

    for entry in sample:
        pmc = entry["pmc_id"]
        fname = entry["file_name"]
        path = KOH / pmc / fname
        if not path.exists():
            continue
        try:
            r = detect_file(str(path), row_context_boost=False)
        except Exception:
            continue
        if not r.suspicions:
            continue
        # Skip decimal-comma : informational, not a corruption claim
        real = [s for s in r.suspicions if s.kind != "decimal-comma"]
        if not real:
            continue
        print(f"\n=== {pmc}/{fname} ({len(real)} flags) ===")
        for s in real[:6]:
            print(f"  sheet={s.sheet!r} col={s.column!r} row={s.row} "
                  f"kind={s.kind!r} conf={s.confidence:.2f}")
            print(f"     value={s.value!r}")
            print(f"     suggestion={s.suggestion!r}")
            print(f"     reason: {s.reason[:160]}")
        if len(real) > 6:
            print(f"  ... +{len(real)-6} more")
    return 0


if __name__ == "__main__":
    sys.exit(main())
