"""Multi-process Koh-corpus walk.

Walks every .xlsx/.xls under /research/koh, runs the gene-symbol prefilter,
then detect_file. Records per-file results and aggregates totals.

Output: /work/koh_walk_fresh.jsonl  (single JSON line)
        /work/koh_walk_fresh.checkpoint.json  (live counters)
"""
from __future__ import annotations

import json
import multiprocessing as mp
import os
import time
import warnings
from collections import Counter
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

warnings.filterwarnings("ignore")

CACHE_ROOT = Path("/research")
LEDGER = Path("/work/koh_walk_fresh.jsonl")
CHECKPOINT = Path("/work/koh_walk_fresh.checkpoint.json")
WORKERS = int(os.environ.get("UNCORRUPT_WORKERS", "6"))


@dataclass
class FileResult:
    pmc_id: str
    file_name: str
    has_gene_symbols: bool
    n_suspicions: int
    n_high: int = 0
    n_medium: int = 0
    n_low: int = 0
    n_info: int = 0
    error: str | None = None


@dataclass
class RunSummary:
    timestamp: str
    n_files_on_disk: int
    n_files_scanned: int
    n_files_failed_load: int
    n_with_gene_symbols: int
    n_with_corruption: int
    n_total_suspicions: int
    n_high_total: int
    n_medium_total: int
    n_low_total: int
    n_info_total: int
    runtime_seconds: float
    workers: int
    per_pmc_with_corruption: int = 0
    files: list[FileResult] = field(default_factory=list)


def _has_gene_symbols(df: pd.DataFrame, hgnc_set: set[str],
                       min_distinct: int = 5,
                       scan_limit: int = 5000) -> bool:
    for col in df.columns:
        seen: set[str] = set()
        for i, val in enumerate(df[col]):
            if i >= scan_limit:
                break
            if isinstance(val, str) and val in hgnc_set:
                seen.add(val)
                if len(seen) >= min_distinct:
                    return True
    return False


# Worker-local cached HGNC set (loaded once per worker process)
_HGNC_SET: set[str] | None = None


def _worker_init() -> None:
    global _HGNC_SET
    from uncorrupt.corpus import load_hgnc
    h = load_hgnc()
    _HGNC_SET = h.current_symbols | set(h.prev_symbol_to_current.keys())


def _process_one(fp_str: str) -> dict:
    """Run prefilter + detect_file on a single path."""
    from uncorrupt.app import _load_all_sheets
    from uncorrupt.detector import detect_file

    fp = Path(fp_str)
    pmc_id = fp.parent.name
    try:
        sheets = _load_all_sheets(str(fp))
    except Exception as exc:
        return {
            "pmc_id": pmc_id, "file_name": fp.name,
            "has_gene_symbols": False, "n_suspicions": 0,
            "n_high": 0, "n_medium": 0, "n_low": 0, "n_info": 0,
            "error": f"load: {type(exc).__name__}: {exc}",
        }
    has_genes = any(
        _has_gene_symbols(df, _HGNC_SET) for df in sheets.values()
    )
    if not has_genes:
        return {
            "pmc_id": pmc_id, "file_name": fp.name,
            "has_gene_symbols": False, "n_suspicions": 0,
            "n_high": 0, "n_medium": 0, "n_low": 0, "n_info": 0,
            "error": None,
        }
    try:
        report = detect_file(str(fp))
    except Exception as exc:
        return {
            "pmc_id": pmc_id, "file_name": fp.name,
            "has_gene_symbols": True, "n_suspicions": 0,
            "n_high": 0, "n_medium": 0, "n_low": 0, "n_info": 0,
            "error": f"detect: {type(exc).__name__}: {exc}",
        }
    n_high = sum(1 for s in report.suspicions if s.confidence >= 0.85)
    n_med = sum(1 for s in report.suspicions
                if 0.50 <= s.confidence < 0.85)
    n_low = sum(1 for s in report.suspicions
                if 0.30 <= s.confidence < 0.50)
    n_info = sum(1 for s in report.suspicions if s.confidence < 0.30)
    return {
        "pmc_id": pmc_id, "file_name": fp.name,
        "has_gene_symbols": True,
        "n_suspicions": len(report.suspicions),
        "n_high": n_high, "n_medium": n_med, "n_low": n_low, "n_info": n_info,
        "error": None,
    }


def main() -> None:
    t0 = time.monotonic()
    targets = sorted(
        str(p) for p in CACHE_ROOT.rglob("*")
        if p.is_file() and p.suffix.lower() in (".xlsx", ".xls")
        and not p.name.startswith("_")
    )
    n_on_disk = len(targets)
    print(f"[koh-walk] {n_on_disk:,} files to scan, {WORKERS} workers",
          flush=True)

    results: list[dict] = []
    n_done = 0
    band_totals = Counter()
    n_failed = 0
    pmcs_with_corruption: set[str] = set()
    last_print = 0
    last_checkpoint = 0

    def _checkpoint() -> None:
        CHECKPOINT.parent.mkdir(parents=True, exist_ok=True)
        with CHECKPOINT.open("w") as f:
            json.dump({
                "n_done": n_done,
                "n_on_disk": n_on_disk,
                "elapsed": time.monotonic() - t0,
                "band_totals": dict(band_totals),
                "n_failed": n_failed,
                "n_with_corruption_pmcs": len(pmcs_with_corruption),
            }, f)

    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=WORKERS, initializer=_worker_init) as pool:
        for r in pool.imap_unordered(_process_one, targets, chunksize=8):
            results.append(r)
            n_done += 1
            if r.get("error"):
                n_failed += 1
            else:
                band_totals["high"] += r["n_high"]
                band_totals["medium"] += r["n_medium"]
                band_totals["low"] += r["n_low"]
                band_totals["info"] += r["n_info"]
                if r["has_gene_symbols"] and r["n_suspicions"] > 0:
                    pmcs_with_corruption.add(r["pmc_id"])
            if n_done - last_print >= 200:
                last_print = n_done
                elapsed = time.monotonic() - t0
                rate = n_done / elapsed
                eta = (n_on_disk - n_done) / rate / 60 if rate > 0 else -1
                pct = 100.0 * n_done / n_on_disk
                print(f"[koh-walk] {n_done:>6}/{n_on_disk} ({pct:4.1f}%) "
                      f"rate={rate:.2f}/s eta={eta:.1f}min "
                      f"flagged_pmcs={len(pmcs_with_corruption)} "
                      f"high={band_totals['high']} med={band_totals['medium']}",
                      flush=True)
            if n_done - last_checkpoint >= 500:
                last_checkpoint = n_done
                _checkpoint()

    runtime = time.monotonic() - t0
    n_with_genes = sum(1 for r in results if r.get("has_gene_symbols"))
    n_with_corr = sum(
        1 for r in results
        if r.get("has_gene_symbols") and r["n_suspicions"] > 0
    )
    n_total_susp = sum(r.get("n_suspicions", 0) for r in results)

    summary = RunSummary(
        timestamp=datetime.now(UTC).isoformat(),
        n_files_on_disk=n_on_disk,
        n_files_scanned=len(results),
        n_files_failed_load=n_failed,
        n_with_gene_symbols=n_with_genes,
        n_with_corruption=n_with_corr,
        n_total_suspicions=n_total_susp,
        n_high_total=band_totals["high"],
        n_medium_total=band_totals["medium"],
        n_low_total=band_totals["low"],
        n_info_total=band_totals["info"],
        per_pmc_with_corruption=len(pmcs_with_corruption),
        runtime_seconds=runtime,
        workers=WORKERS,
        files=[FileResult(**r) for r in results],
    )

    LEDGER.parent.mkdir(parents=True, exist_ok=True)
    with LEDGER.open("w") as f:
        f.write(json.dumps(asdict(summary)) + "\n")

    print()
    print(f"=== FRESH KOH WALK ({WORKERS} workers) ===")
    print(f"  files on disk        : {summary.n_files_on_disk:,}")
    print(f"  files scanned        : {summary.n_files_scanned:,}")
    print(f"  failed to load       : {summary.n_files_failed_load:,}")
    print(f"  with gene symbols    : {summary.n_with_gene_symbols:,}")
    print(f"  with corruption flag : {summary.n_with_corruption:,}")
    print(f"  total suspicions     : {summary.n_total_suspicions:,}")
    print(f"    high (>=0.85)      : {summary.n_high_total:,}")
    print(f"    medium (0.50-0.85) : {summary.n_medium_total:,}")
    print(f"    low (0.30-0.50)    : {summary.n_low_total:,}")
    print(f"    info (<0.30)       : {summary.n_info_total:,}")
    print(f"  PMCs with corruption : {summary.per_pmc_with_corruption:,}")
    print(f"  runtime              : {summary.runtime_seconds/60:.1f} min "
          f"({summary.runtime_seconds:.0f}s)")
    print(f"  ledger written       : {LEDGER}")


if __name__ == "__main__":
    main()
