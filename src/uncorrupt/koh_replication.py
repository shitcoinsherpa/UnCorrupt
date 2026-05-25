"""Koh-2022-style benchmark replicated at substantially larger scale.

Koh et al. (Sci Rep 2022) sampled supplementary files from 11 journals during
June 2022 (1 month, 356 files harvested, 81 with gene symbols, 28 with date
errors, 22 auto-corrected). We replicate the same 11 journals over a longer
window using the same three headline metrics so the comparison is apples-to-
apples, then add per-journal and per-year breakdowns.

Source channels (legitimate per NCBI policy):
- NCBI E-utilities (esearch) for article discovery by journal + date
- Europe PMC supplementaryFiles ZIP for the actual files
- NCBI OA tarball as fallback (already in validate.py)
"""
from __future__ import annotations

import io
import json
import re
import time
import zipfile
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path

import httpx
import pandas as pd

from .corpus import load_hgnc
from .detector import detect_file

# Koh 2022 list of journals, with NCBI-friendly forms.
# Confirmed via NCBI's catalog (e.g., "PLoS One" not "PLOS ONE" for older content).
KOH_JOURNALS: list[str] = [
    "BMC Genomics",
    "Nature",
    "Genome Biology",
    "Nucleic Acids Research",
    "Human Molecular Genetics",
    "BMC Bioinformatics",
    "Nature Communications",
    "PLoS One",
    "Genome Research",
    "Genes & Development",
    "RNA",
]

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EPMC_SUPPL_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/{full_pmc_id}/supplementaryFiles"

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data" / "raw" / "koh_replication"
LEDGER_PATH = PROJECT_ROOT / "results" / "koh_replication.jsonl"

USER_AGENT = "uncorrupt-research/0.1 (https://github.com/shitcoinsherpa)"


@dataclass
class KohFileResult:
    pmc_id: str
    journal: str
    year: int | None
    file_name: str
    has_gene_symbols: bool
    n_corruption_suspicions: int
    identifier_columns: list[str] = field(default_factory=list)
    error: str | None = None


@dataclass
class KohRun:
    timestamp: str
    journals: list[str]
    year_min: int
    year_max: int
    runtime_seconds: float
    n_articles_queried: int
    n_supplementary_files: int
    n_with_gene_symbols: int
    n_with_corruption: int
    files: list[KohFileResult] = field(default_factory=list)


def _esearch_pmc(journal: str, year_min: int, year_max: int, client: httpx.Client,
                 month: int | None = None) -> list[str]:
    """Return PMC IDs for articles in the given journal/date window.

    If `month` is given, constrains to that specific month in each year of the
    range (use year_min == year_max for a single month).
    """
    if month is not None:
        day_max = 31 if month in (1, 3, 5, 7, 8, 10, 12) else (29 if month == 2 else 30)
        date_clause = (
            f'"{year_min}/{month:02d}/01"[PubDate] : '
            f'"{year_max}/{month:02d}/{day_max:02d}"[PubDate]'
        )
    else:
        date_clause = f'"{year_min}/01/01"[PubDate] : "{year_max}/12/31"[PubDate]'
    params = {
        "db": "pmc",
        "term": f'"{journal}"[Journal] AND ({date_clause})',
        "retmode": "json",
        "retmax": "10000",
        "tool": "uncorrupt-research",
        "contact": "https://github.com/shitcoinsherpa",
    }
    resp = client.get(f"{EUTILS_BASE}/esearch.fcgi", params=params, timeout=60)
    resp.raise_for_status()
    data = resp.json()
    return data.get("esearchresult", {}).get("idlist", [])


def _esummary_pmc(pmc_ids: list[str], client: httpx.Client) -> dict[str, dict]:
    """Bulk-fetch article metadata for PMC IDs. Returns {pmc_id: summary}."""
    if not pmc_ids:
        return {}
    params = {
        "db": "pmc",
        "id": ",".join(pmc_ids),
        "retmode": "json",
        "tool": "uncorrupt-research",
        "contact": "https://github.com/shitcoinsherpa",
    }
    resp = client.get(f"{EUTILS_BASE}/esummary.fcgi", params=params, timeout=120)
    resp.raise_for_status()
    data = resp.json()
    result = data.get("result", {})
    return {k: v for k, v in result.items() if k != "uids"}


def _fetch_epmc_zip(full_pmc_id: str, client: httpx.Client) -> bytes | None:
    """Fetch the Europe PMC supplementary ZIP. `full_pmc_id` is the form
    `PMC12345` (i.e. with the prefix). Returns None on 404, HTTPError, or
    body too small to be a real ZIP (EPMC returns a ~165-byte stub for
    articles with no supplementaries)."""
    try:
        resp = client.get(EPMC_SUPPL_URL.format(full_pmc_id=full_pmc_id),
                          follow_redirects=True, timeout=120)
        if resp.status_code == 404:
            return None
        resp.raise_for_status()
        # EPMC sometimes returns 200 with a tiny no-content payload
        if len(resp.content) < 1000:
            return None
        return resp.content
    except httpx.HTTPError:
        return None


def _iter_xlsx_in_zip(zip_bytes: bytes):
    """Yield (filename, bytes) for each xlsx/xls inside the ZIP. Skips non-tabular."""
    try:
        with zipfile.ZipFile(io.BytesIO(zip_bytes)) as z:
            for info in z.infolist():
                if info.is_dir():
                    continue
                name = info.filename.rsplit("/", 1)[-1]
                if name.lower().endswith((".xlsx", ".xls", ".csv", ".tsv")):
                    yield name, z.read(info)
    except zipfile.BadZipFile:
        return


def has_gene_symbols(df: pd.DataFrame, hgnc_set: set[str], min_distinct: int = 5,
                    scan_limit: int = 5000) -> bool:
    """Match Koh's criterion: file 'contains gene symbols'.

    Heuristic: at least `min_distinct` distinct HGNC current-or-prev symbols
    appear in any single column. Bounded by `scan_limit` rows per column to
    keep large files cheap.
    """
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


def has_gene_symbols_in_file(path: Path, hgnc_set: set[str]) -> bool:
    """Check across all sheets of a file."""
    from .app import _load_all_sheets
    try:
        sheets = _load_all_sheets(str(path))
    except Exception:
        return False
    for df in sheets.values():
        if has_gene_symbols(df, hgnc_set):
            return True
    return False


def run_koh_replication(
    year_min: int = 2022,
    year_max: int = 2026,
    journals: list[str] | None = None,
    rate_limit_per_second: float = 3.0,
    checkpoint_path: Path | None = None,
    month_year_pairs: list[tuple[int, int]] | None = None,
    pmc_id_offset: int = 0,
    pmc_id_limit: int | None = None,
) -> KohRun:
    """Run the Koh-2022-style benchmark.

    If `month_year_pairs` is provided, esearch is constrained to those specific
    (year, month) windows : replicates Koh's single-month methodology at the
    requested temporal scale.

    `pmc_id_offset` and `pmc_id_limit` slice the esearch-returned PMC IDs per
    journal (after de-dup): processes IDs [offset : offset+limit). Used by
    the parallel runner to split a single journal's IDs across N workers.
    """
    import sys

    journals = list(journals) if journals else list(KOH_JOURNALS)
    t0 = time.monotonic()
    delay = 1.0 / rate_limit_per_second

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    headers = {"User-Agent": USER_AGENT}

    hgnc = load_hgnc()
    hgnc_set = hgnc.current_symbols | set(hgnc.prev_symbol_to_current.keys())
    print(f"[koh] HGNC pool: {len(hgnc_set)} symbols", file=sys.stderr)
    if month_year_pairs:
        print(f"[koh] month-year windows: {month_year_pairs}", file=sys.stderr)

    results: list[KohFileResult] = []
    n_articles_queried = 0

    with httpx.Client(headers=headers) as client:
        for j_idx, journal in enumerate(journals, start=1):
            print(f"[koh] journal {j_idx}/{len(journals)}: {journal}", file=sys.stderr)
            # Build PMC ID list across all configured month-year windows
            pmc_ids: list[str] = []
            try:
                if month_year_pairs:
                    for year, month in month_year_pairs:
                        ids = _esearch_pmc(journal, year, year, client, month=month)
                        pmc_ids.extend(ids)
                        time.sleep(delay)
                        print(f"[koh]   {journal} {year}-{month:02d}: {len(ids)} IDs",
                              file=sys.stderr)
                    # De-duplicate
                    pmc_ids = list(dict.fromkeys(pmc_ids))
                else:
                    pmc_ids = _esearch_pmc(journal, year_min, year_max, client)
                    time.sleep(delay)
            except Exception as exc:
                print(f"[koh]   esearch failed: {exc}", file=sys.stderr)
                continue
            print(f"[koh]   {len(pmc_ids)} PMC IDs total for {journal}", file=sys.stderr)
            n_articles_queried += len(pmc_ids)
            if pmc_id_offset or pmc_id_limit is not None:
                end = (pmc_id_offset + pmc_id_limit) if pmc_id_limit is not None else len(pmc_ids)
                pmc_ids = pmc_ids[pmc_id_offset:end]
                print(f"[koh]   sliced to [{pmc_id_offset}:{end}] -> {len(pmc_ids)} PMC IDs",
                      file=sys.stderr)

            # Bulk-fetch summaries in batches of 200
            id_to_year: dict[str, int | None] = {}
            for batch_start in range(0, len(pmc_ids), 200):
                batch = pmc_ids[batch_start:batch_start + 200]
                try:
                    sums = _esummary_pmc(batch, client)
                    time.sleep(delay)
                    for pid, s in sums.items():
                        pub = s.get("pubdate", "") or s.get("epubdate", "")
                        m = re.match(r"(\d{4})", pub)
                        id_to_year[pid] = int(m.group(1)) if m else None
                except Exception as exc:
                    print(f"[koh]   esummary batch failed: {exc}", file=sys.stderr)

            for i, pmc_id in enumerate(pmc_ids, start=1):
                # Progress check at top so it fires even when articles continue early
                if i % 50 == 0 or i == 1:
                    n_so_far = len(results)
                    n_w_genes = sum(1 for r in results if r.has_gene_symbols)
                    n_w_corr = sum(1 for r in results
                                   if r.has_gene_symbols and r.n_corruption_suspicions > 0)
                    elapsed = time.monotonic() - t0
                    print(
                        f"[koh]   {journal} {i}/{len(pmc_ids)}  "
                        f"files={n_so_far} with_genes={n_w_genes} with_corr={n_w_corr}  "
                        f"elapsed={elapsed:.0f}s",
                        file=sys.stderr,
                    )
                    sys.stderr.flush()
                    if checkpoint_path is not None and i % 50 == 0:
                        partial = KohRun(
                            timestamp=datetime.now(UTC).isoformat(),
                            journals=journals, year_min=year_min, year_max=year_max,
                            runtime_seconds=elapsed,
                            n_articles_queried=n_articles_queried,
                            n_supplementary_files=n_so_far,
                            n_with_gene_symbols=n_w_genes,
                            n_with_corruption=n_w_corr,
                            files=list(results),
                        )
                        # WSL2 DrvFs on /mnt/ drives sporadically raises
                        # OSError(EINVAL) on rewrite. A skipped checkpoint
                        # is recoverable; a crashed multi-hour walk is not.
                        try:
                            with checkpoint_path.open("w") as fh:
                                json.dump(asdict(partial), fh)
                        except OSError as exc:
                            print(f"[koh]   checkpoint write skipped: {exc}",
                                  file=sys.stderr)

                full_pmc = f"PMC{pmc_id}"
                article_cache_dir = DATA_DIR / "files" / full_pmc
                article_cache_dir.mkdir(parents=True, exist_ok=True)
                year = id_to_year.get(pmc_id)

                # Get supplementary files, caching the ZIP extract list
                manifest_path = article_cache_dir / "_manifest.json"
                if not manifest_path.exists():
                    zip_bytes = _fetch_epmc_zip(full_pmc, client)
                    time.sleep(delay)
                    if zip_bytes is None:
                        manifest_path.write_text(json.dumps({"files": [], "error": "no epmc"}))
                        continue
                    file_names: list[str] = []
                    for name, content in _iter_xlsx_in_zip(zip_bytes):
                        out = article_cache_dir / name
                        out.write_bytes(content)
                        file_names.append(name)
                    manifest_path.write_text(json.dumps({"files": file_names}))
                manifest = json.loads(manifest_path.read_text())
                file_names = manifest.get("files", [])

                for name in file_names:
                    local = article_cache_dir / name
                    if not local.exists():
                        continue
                    try:
                        has_genes = has_gene_symbols_in_file(local, hgnc_set)
                    except Exception as exc:
                        results.append(KohFileResult(
                            pmc_id=full_pmc, journal=journal, year=year,
                            file_name=name, has_gene_symbols=False,
                            n_corruption_suspicions=0,
                            error=f"gene-screen failed: {type(exc).__name__}",
                        ))
                        continue
                    if not has_genes:
                        results.append(KohFileResult(
                            pmc_id=full_pmc, journal=journal, year=year,
                            file_name=name, has_gene_symbols=False,
                            n_corruption_suspicions=0,
                        ))
                        continue
                    try:
                        report = detect_file(str(local))
                        results.append(KohFileResult(
                            pmc_id=full_pmc, journal=journal, year=year,
                            file_name=name, has_gene_symbols=True,
                            n_corruption_suspicions=len(report.suspicions),
                            identifier_columns=report.identifier_columns,
                        ))
                    except Exception as exc:
                        results.append(KohFileResult(
                            pmc_id=full_pmc, journal=journal, year=year,
                            file_name=name, has_gene_symbols=True,
                            n_corruption_suspicions=0,
                            error=f"detect failed: {type(exc).__name__}",
                        ))

                # Per-journal progress now handled at top of loop

    elapsed = time.monotonic() - t0
    n_with_genes = sum(1 for r in results if r.has_gene_symbols)
    n_with_corr = sum(1 for r in results
                      if r.has_gene_symbols and r.n_corruption_suspicions > 0)

    return KohRun(
        timestamp=datetime.now(UTC).isoformat(),
        journals=journals,
        year_min=year_min,
        year_max=year_max,
        runtime_seconds=elapsed,
        n_articles_queried=n_articles_queried,
        n_supplementary_files=len(results),
        n_with_gene_symbols=n_with_genes,
        n_with_corruption=n_with_corr,
        files=results,
    )


def format_koh_report(run: KohRun) -> str:
    from collections import Counter
    lines = [
        f"Koh-2022-replication run @ {run.timestamp}",
        f"  Window: {run.year_min}-{run.year_max}",
        f"  Journals: {len(run.journals)}",
        f"  Articles queried: {run.n_articles_queried}",
        f"  Suppl. files: {run.n_supplementary_files}",
        f"  With gene syms: {run.n_with_gene_symbols}  "
        f"({run.n_with_gene_symbols / max(1, run.n_supplementary_files) * 100:.1f}%)",
        f"  With corruption: {run.n_with_corruption}  "
        f"({run.n_with_corruption / max(1, run.n_with_gene_symbols) * 100:.1f}% "
        f"of gene-symbol files)",
        f"  Runtime: {run.runtime_seconds / 60:.1f} min",
        "",
        "Per-journal:",
    ]
    by_j = Counter()
    by_j_genes = Counter()
    by_j_corr = Counter()
    for r in run.files:
        by_j[r.journal] += 1
        if r.has_gene_symbols:
            by_j_genes[r.journal] += 1
            if r.n_corruption_suspicions > 0:
                by_j_corr[r.journal] += 1
    for journal in run.journals:
        n = by_j.get(journal, 0)
        ng = by_j_genes.get(journal, 0)
        nc = by_j_corr.get(journal, 0)
        lines.append(
            f"  {journal:<32}  files={n:>5}  with_genes={ng:>5}  "
            f"with_corruption={nc:>5}  "
            f"({nc / max(1, ng) * 100:>5.1f}% of gene-symbol files)"
        )
    return "\n".join(lines)


def append_to_ledger(run: KohRun, path: Path = LEDGER_PATH) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(json.dumps(asdict(run)) + "\n")


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--year-min", type=int, default=2022)
    p.add_argument("--year-max", type=int, default=2026)
    p.add_argument("--journals", nargs="+", help="Subset of journals (default: all 11)")
    p.add_argument("--rate", type=float, default=3.0)
    p.add_argument("--months", nargs="+", help="Specific YYYY-MM windows, e.g. 2022-06 2023-06")
    args = p.parse_args()

    month_pairs = None
    if args.months:
        month_pairs = []
        for ym in args.months:
            y, m = ym.split("-")
            month_pairs.append((int(y), int(m)))

    checkpoint = PROJECT_ROOT / "results" / "koh_replication_partial.json"
    run = run_koh_replication(
        year_min=args.year_min,
        year_max=args.year_max,
        journals=args.journals,
        rate_limit_per_second=args.rate,
        checkpoint_path=checkpoint,
        month_year_pairs=month_pairs,
    )
    print(format_koh_report(run))
    append_to_ledger(run)
    print(f"\nAppended to {LEDGER_PATH}")
