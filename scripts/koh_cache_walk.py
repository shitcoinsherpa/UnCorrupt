"""Run the iter-5 detector against every cached Koh corpus file on disk.

Unlike `uncorrupt.koh_replication`, this script bypasses the per-article
manifest cache layer. It walks `data/raw/koh_replication/files/PMC*/` for
every `.xlsx`/`.xls`, runs the gene-symbol screen + detector, and writes a
fresh ledger entry.

This is the apples-to-apples comparison we need:
- the pre-iter-5 run captured 8,441 files (with a healthy EPMC manifest cache)
- the iter-5 manifest-based run captured only 2,079 (stale manifests skipped
  several journals)
- this cache-walk captures *every cached file on disk* and runs iter-5 detect

Output: `results/koh_cache_walk.jsonl` (one line per run).
"""
from __future__ import annotations

import json
import sys
import time
import warnings
from collections import Counter
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path

warnings.filterwarnings("ignore")

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from uncorrupt.corpus import load_hgnc  # noqa: E402
from uncorrupt.detector import detect_file  # noqa: E402
from uncorrupt.koh_replication import has_gene_symbols_in_file  # noqa: E402

CACHE_ROOT = PROJECT_ROOT / "data" / "raw" / "koh_replication" / "files"
LEDGER = PROJECT_ROOT / "results" / "koh_cache_walk.jsonl"


@dataclass
class FileResult:
    pmc_id: str
    file_name: str
    has_gene_symbols: bool
    n_suspicions: int
    error: str | None = None


@dataclass
class RunSummary:
    timestamp: str
    n_files_scanned: int
    n_files_failed_load: int
    n_with_gene_symbols: int
    n_with_corruption: int
    runtime_seconds: float
    per_pmc_with_corruption: int = 0  # how many distinct PMCs have ≥1 corrupted file
    files: list[FileResult] = field(default_factory=list)


def main() -> None:
    hgnc = load_hgnc()
    hgnc_set = hgnc.current_symbols | set(hgnc.prev_symbol_to_current.keys())
    print(f"[walk] HGNC pool: {len(hgnc_set):,} symbols", file=sys.stderr)

    targets = sorted(
        p for p in CACHE_ROOT.rglob("*")
        if p.is_file() and p.suffix.lower() in (".xlsx", ".xls", ".XLSX", ".XLS")
        and not p.name.startswith("_")
    )
    print(f"[walk] {len(targets):,} files to scan", file=sys.stderr)

    results: list[FileResult] = []
    pmcs_with_corruption: set[str] = set()
    t0 = time.monotonic()
    n_failed = 0

    for i, fp in enumerate(targets, start=1):
        pmc_id = fp.parent.name  # parent dir is PMC{ID}
        # Gene-symbol screen first (cheap)
        try:
            has_genes = has_gene_symbols_in_file(fp, hgnc_set)
        except Exception as exc:
            results.append(FileResult(
                pmc_id=pmc_id, file_name=fp.name,
                has_gene_symbols=False, n_suspicions=0,
                error=f"gene-screen: {type(exc).__name__}: {exc}",
            ))
            n_failed += 1
            continue
        if not has_genes:
            results.append(FileResult(
                pmc_id=pmc_id, file_name=fp.name,
                has_gene_symbols=False, n_suspicions=0,
            ))
            if i % 200 == 0:
                _progress(i, len(targets), results, t0, pmcs_with_corruption)
            continue
        # Full detect
        try:
            report = detect_file(str(fp))
            n_susp = len(report.suspicions)
            results.append(FileResult(
                pmc_id=pmc_id, file_name=fp.name,
                has_gene_symbols=True, n_suspicions=n_susp,
            ))
            if n_susp > 0:
                pmcs_with_corruption.add(pmc_id)
        except Exception as exc:
            results.append(FileResult(
                pmc_id=pmc_id, file_name=fp.name,
                has_gene_symbols=True, n_suspicions=0,
                error=f"detect: {type(exc).__name__}: {exc}",
            ))
            n_failed += 1
        if i % 200 == 0:
            _progress(i, len(targets), results, t0, pmcs_with_corruption)

    runtime = time.monotonic() - t0
    n_with_genes = sum(1 for r in results if r.has_gene_symbols)
    n_with_corr = sum(1 for r in results if r.has_gene_symbols and r.n_suspicions > 0)

    summary = RunSummary(
        timestamp=datetime.now(UTC).isoformat(),
        n_files_scanned=len(results),
        n_files_failed_load=n_failed,
        n_with_gene_symbols=n_with_genes,
        n_with_corruption=n_with_corr,
        per_pmc_with_corruption=len(pmcs_with_corruption),
        runtime_seconds=runtime,
        files=results,
    )
    LEDGER.parent.mkdir(parents=True, exist_ok=True)
    with LEDGER.open("a") as f:
        f.write(json.dumps(asdict(summary)) + "\n")

    print(
        f"\n=== DONE ({runtime/60:.1f} min) ===\n"
        f"  files scanned          : {summary.n_files_scanned:,}\n"
        f"  with gene symbols      : {summary.n_with_gene_symbols:,}\n"
        f"  with corruption        : {summary.n_with_corruption:,}\n"
        f"  PMCs with ≥1 corrupted : {summary.per_pmc_with_corruption:,}\n"
        f"  load/detect failures   : {summary.n_files_failed_load:,}",
        file=sys.stderr,
    )
    rate = (n_with_corr / n_with_genes * 100) if n_with_genes else 0.0
    print(f"  corruption rate (of gene-files): {rate:.1f}%", file=sys.stderr)
    print(f"\nAppended to {LEDGER}", file=sys.stderr)


def _progress(i: int, total: int, results: list[FileResult], t0: float,
              pmcs_corr: set[str]) -> None:
    n_g = sum(1 for r in results if r.has_gene_symbols)
    n_c = sum(1 for r in results if r.has_gene_symbols and r.n_suspicions > 0)
    elapsed = time.monotonic() - t0
    print(
        f"[walk] {i:,}/{total:,}  with_genes={n_g:,}  with_corr={n_c:,}  "
        f"pmcs_corr={len(pmcs_corr):,}  elapsed={elapsed/60:.1f}min",
        file=sys.stderr, flush=True,
    )


if __name__ == "__main__":
    main()
