"""Parallel EPMC corpus expansion : splits journals across N worker
processes for ~Nx speedup on multi-core machines.

The single-threaded `expand_corpus_via_epmc.py` was using one core
out of eight. The detect_file() calls per cached XLSX are CPU-bound
and embarrassingly parallel across journals; this script splits the
17-journal list into N chunks and runs each chunk in a subprocess.

Each worker:
  - calls `run_koh_replication(journals=chunk, checkpoint_path=worker_N.json)`
  - appends a KohRun to results/koh_replication.jsonl on exit
  - logs to results/epmc_worker_N.log

When all workers exit, the merger can combine all KohRun records
from koh_replication.jsonl into a unified file list.

Run:
    python scripts/expand_corpus_parallel.py --workers 4 --rate 1.5
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from uncorrupt.koh_replication import KOH_JOURNALS  # noqa: E402

EXTENDED_JOURNALS = KOH_JOURNALS + [
    "Cell",
    "Cell Reports",
    "Bioinformatics",
    "Briefings in Bioinformatics",
    "eLife",
    "Scientific Reports",
]


WORKER_BODY = '''
import sys
from pathlib import Path
ROOT = Path("{root}")
sys.path.insert(0, str(ROOT / "src"))
from uncorrupt.koh_replication import run_koh_replication, append_to_ledger
journals = {journals!r}
run = run_koh_replication(
    year_min={year_min},
    year_max={year_max},
    journals=journals,
    rate_limit_per_second={rate},
    checkpoint_path=ROOT / "results/koh_cache_walk_worker_{wid}.checkpoint.json",
)
append_to_ledger(run)
print(f"[worker-{wid}] DONE journals={{journals}} files={{run.n_supplementary_files}} corruption={{run.n_with_corruption}} runtime={{run.runtime_seconds:.0f}}s")
'''


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--rate", type=float, default=1.5,
                   help="requests/sec per worker (total = workers * rate)")
    p.add_argument("--year-min", type=int, default=2019)
    p.add_argument("--year-max", type=int, default=2026)
    args = p.parse_args()

    journals = EXTENDED_JOURNALS
    n = len(journals)
    w = args.workers
    chunks = [journals[i::w] for i in range(w)]

    print(f"[parallel] launching {w} workers, rate {args.rate}/sec each "
          f"({w * args.rate:.1f}/sec total)")
    for i, chunk in enumerate(chunks):
        print(f"[parallel] worker {i}: {chunk}")

    log_dir = ROOT / "results"
    log_dir.mkdir(parents=True, exist_ok=True)

    pids: list[tuple[int, subprocess.Popen]] = []
    for wid, chunk in enumerate(chunks):
        if not chunk:
            continue
        body = WORKER_BODY.format(
            root=str(ROOT), journals=chunk,
            year_min=args.year_min, year_max=args.year_max,
            rate=args.rate, wid=wid,
        )
        log_path = log_dir / f"epmc_worker_{wid}.log"
        log_fh = log_path.open("w")
        env = os.environ.copy()
        env["PYTHONUNBUFFERED"] = "1"
        proc = subprocess.Popen(
            [sys.executable, "-u", "-c", body],
            stdout=log_fh, stderr=subprocess.STDOUT,
            env=env, cwd=str(ROOT),
        )
        print(f"[parallel] worker {wid} pid={proc.pid} log={log_path.name}")
        pids.append((wid, proc))

    # Wait for all workers
    print(f"[parallel] waiting for {len(pids)} workers")
    t0 = time.monotonic()
    for wid, proc in pids:
        rc = proc.wait()
        elapsed = time.monotonic() - t0
        print(f"[parallel] worker {wid} exited rc={rc} at {elapsed:.0f}s")

    print(f"[parallel] all workers done at {time.monotonic() - t0:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
