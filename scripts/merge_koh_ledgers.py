"""Merge the original Koh walk ledger and the EPMC-expanded walk ledger
into a single unified ledger that `calibrate_with_fp_anchors.py` can
sample from.

The original ledger (`results/koh_cache_walk.jsonl`) carries
`n_suspicions`; the EPMC checkpoint
(`results/koh_cache_walk_extended.checkpoint.json`) carries
`n_corruption_suspicions` plus extra fields (`journal`, `year`,
`identifier_columns`).

This script:
  - reads both ledgers
  - normalises the schema (drops extras, renames the suspicion field)
  - de-duplicates by (pmc_id, file_name)
  - writes the merged ledger to `results/koh_cache_walk_merged.jsonl`

Run:
    python scripts/merge_koh_ledgers.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT / "results/koh_cache_walk.jsonl"
EPMC = ROOT / "results/koh_cache_walk_extended.checkpoint.json"
EPMC_FINAL = ROOT / "results/koh_replication.jsonl"
OUT = ROOT / "results/koh_cache_walk_merged.jsonl"


def _normalize(f: dict) -> dict:
    """Project a file record to the old schema."""
    return {
        "pmc_id": f["pmc_id"],
        "file_name": f["file_name"],
        "has_gene_symbols": f.get("has_gene_symbols", False),
        "n_suspicions": f.get("n_suspicions", f.get("n_corruption_suspicions", 0)),
        "error": f.get("error"),
    }


def main() -> int:
    if not OLD.exists():
        print(f"FATAL: {OLD} missing", file=sys.stderr)
        return 1
    with OLD.open() as fh:
        old = json.load(fh)
    old_files = old.get("files", [])
    print(f"[merge] old ledger: {len(old_files)} files")

    # Combine ALL KohRun records from koh_replication.jsonl : multiple
    # parallel workers each append their own KohRun, and we want files
    # from every run that touched this session's cache.
    epmc_files: list[dict] = []
    if EPMC_FINAL.exists():
        runs = [json.loads(line) for line in EPMC_FINAL.open()]
        for run in runs:
            epmc_files.extend(run.get("files", []))
        print(f"[merge] epmc final ledger: {len(epmc_files)} files combined "
              f"across {len(runs)} runs")
    # Also pick up in-progress worker checkpoints (defensive : if a worker
    # crashed before append, its checkpoint may have files the ledger doesn't)
    for cp in sorted(ROOT.glob("results/koh_cache_walk_worker_*.checkpoint.json")):
        try:
            with cp.open() as fh:
                cpdata = json.load(fh)
            epmc_files.extend(cpdata.get("files", []))
            print(f"[merge] picked up {cp.name}: {len(cpdata.get('files', []))} files")
        except Exception as exc:
            print(f"[merge] skip {cp.name}: {exc}", file=sys.stderr)
    # Also include the original single-threaded checkpoint if present
    if EPMC.exists():
        try:
            with EPMC.open() as fh:
                cpdata = json.load(fh)
            epmc_files.extend(cpdata.get("files", []))
            print(f"[merge] picked up {EPMC.name}: {len(cpdata.get('files', []))} files")
        except Exception as exc:
            print(f"[merge] skip {EPMC.name}: {exc}", file=sys.stderr)
    if not epmc_files:
        print(f"[merge] no EPMC ledgers found; merged == old", file=sys.stderr)

    # De-duplicate by (pmc_id, file_name). Prefer EPMC records when present
    # since the detector that produced them is the more recent build.
    by_key: dict[tuple[str, str], dict] = {}
    for f in old_files:
        by_key[(f["pmc_id"], f["file_name"])] = _normalize(f)
    n_old_only = len(by_key)
    for f in epmc_files:
        by_key[(f["pmc_id"], f["file_name"])] = _normalize(f)
    n_merged = len(by_key)
    n_new = n_merged - n_old_only

    merged = {
        "timestamp": "merged",
        "n_files_scanned": n_merged,
        "n_files_failed_load": sum(1 for f in by_key.values() if f.get("error")),
        "n_with_gene_symbols": sum(1 for f in by_key.values() if f.get("has_gene_symbols")),
        "n_with_corruption": sum(1 for f in by_key.values()
                                 if f.get("has_gene_symbols") and f.get("n_suspicions", 0) > 0),
        "runtime_seconds": 0.0,
        "per_pmc_with_corruption": {},
        "files": list(by_key.values()),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w") as fh:
        json.dump(merged, fh)
    print(f"[merge] wrote {OUT} with {n_merged} files "
          f"({n_old_only} from old, {n_new} new from EPMC)")
    print(f"[merge] n_with_gene_symbols={merged['n_with_gene_symbols']} "
          f"n_with_corruption={merged['n_with_corruption']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
