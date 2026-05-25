"""Precision measurement against random non-Ziemann PMC articles.

The Ziemann 2021 corpus gives us a recall number (we catch 99.2% of confirmed
corrupted articles). It does NOT give us a precision number : what fraction
of articles we flag are actually corrupted vs false positives.

This module samples random PMC genetics articles that are NOT in Ziemann's S2
list, runs the detector, and surfaces every flagged article for manual triage.
True positives are findings Ziemann's tool missed (real corruption in the
wild). False positives reveal detector over-flagging that needs tightening.

Output goes to a separate ledger so the precision corpus stays distinct from
the Koh-replication corpus.
"""
from __future__ import annotations

import json
import random
import time
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path

import httpx

from .corpus import load_ziemann_2021_corpus
from .detector import detect_file
from .koh_replication import (
    DATA_DIR as KOH_DATA_DIR,
    USER_AGENT,
    _esearch_pmc,
    _esummary_pmc,
    _fetch_epmc_zip,
    _iter_xlsx_in_zip,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data" / "raw" / "precision_check"
LEDGER_PATH = PROJECT_ROOT / "results" / "precision_check.jsonl"

# A broad PMC search expression for genetics-adjacent papers.
GENETICS_QUERY_TERM = (
    '("gene expression"[Title/Abstract] OR "RNA-seq"[Title/Abstract] '
    'OR "genomics"[Title/Abstract] OR "transcriptome"[Title/Abstract] '
    'OR "differentially expressed"[Title/Abstract]) '
    'AND open access[filter]'
)


@dataclass
class PrecisionResult:
    pmc_id: str
    journal: str
    year: int | None
    file_name: str
    n_supplementary_files: int
    n_corruption_suspicions: int
    flagged: bool
    sample_suggestions: list[str] = field(default_factory=list)
    error: str | None = None


@dataclass
class PrecisionRun:
    timestamp: str
    seed: int
    year_min: int
    year_max: int
    n_sampled: int
    n_processed: int
    n_with_supplementaries: int
    n_flagged: int
    runtime_seconds: float
    results: list[PrecisionResult] = field(default_factory=list)


def _esearch_genetics(year_min: int, year_max: int, client: httpx.Client) -> list[str]:
    """Find PMC IDs for genetics-adjacent open-access articles."""
    params = {
        "db": "pmc",
        "term": (
            f"{GENETICS_QUERY_TERM} AND "
            f'("{year_min}/01/01"[PubDate] : "{year_max}/12/31"[PubDate])'
        ),
        "retmode": "json",
        "retmax": "100000",
        "tool": "uncorrupt-research",
        "contact": "https://github.com/shitcoinsherpa",
    }
    resp = client.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        params=params, timeout=120,
    )
    resp.raise_for_status()
    return resp.json().get("esearchresult", {}).get("idlist", [])


def _exclude_ziemann_set() -> set[str]:
    """Numeric PMC IDs (without 'PMC' prefix) that appear in Ziemann's S2."""
    entries = load_ziemann_2021_corpus()
    return {e.pmc_id.removeprefix("PMC") for e in entries}


def run_precision_check(
    n_sample: int = 200,
    year_min: int = 2018,
    year_max: int = 2020,
    seed: int = 42,
    rate_limit_per_second: float = 1.5,
    checkpoint_path: Path | None = None,
) -> PrecisionRun:
    """Sample N random non-Ziemann genetics articles, run detector, record flags.

    Defaults to 2018-2020 because (a) Ziemann's corpus covers 2014-2020 so we
    can exclude it cleanly, and (b) older-than-2021 articles have better EPMC
    coverage than recent ones (per our earlier empirical finding).
    """
    import sys

    rng = random.Random(seed)
    t0 = time.monotonic()
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    delay = 1.0 / rate_limit_per_second

    headers = {"User-Agent": USER_AGENT}
    with httpx.Client(headers=headers) as client:
        all_ids = _esearch_genetics(year_min, year_max, client)
        time.sleep(delay)
        ziemann = _exclude_ziemann_set()
        candidates = [pid for pid in all_ids if pid not in ziemann]
        print(
            f"[precision] discovery: {len(all_ids)} genetics IDs, "
            f"{len(ziemann)} excluded, {len(candidates)} candidates",
            file=sys.stderr,
        )
        sampled = rng.sample(candidates, min(n_sample, len(candidates)))
        print(f"[precision] sampling {len(sampled)} for fetch+detect", file=sys.stderr)

        results: list[PrecisionResult] = []
        n_with_suppl = 0

        for i, pmc_id in enumerate(sampled, start=1):
            full_pmc = f"PMC{pmc_id}"
            cache_dir = DATA_DIR / full_pmc
            cache_dir.mkdir(parents=True, exist_ok=True)
            manifest_path = cache_dir / "_manifest.json"

            if not manifest_path.exists():
                zip_bytes = _fetch_epmc_zip(full_pmc, client)
                time.sleep(delay)
                if zip_bytes is None:
                    manifest_path.write_text(json.dumps({"files": [], "error": "no epmc"}))
                    results.append(PrecisionResult(
                        pmc_id=full_pmc, journal="", year=None,
                        file_name="", n_supplementary_files=0,
                        n_corruption_suspicions=0, flagged=False,
                        error="no epmc",
                    ))
                    continue
                file_names: list[str] = []
                for name, content in _iter_xlsx_in_zip(zip_bytes):
                    (cache_dir / name).write_bytes(content)
                    file_names.append(name)
                manifest_path.write_text(json.dumps({"files": file_names}))

            file_names = json.loads(manifest_path.read_text()).get("files", [])
            if not file_names:
                results.append(PrecisionResult(
                    pmc_id=full_pmc, journal="", year=None,
                    file_name="", n_supplementary_files=0,
                    n_corruption_suspicions=0, flagged=False,
                ))
                continue

            n_with_suppl += 1
            for name in file_names:
                local = cache_dir / name
                if not local.exists():
                    continue
                try:
                    rep = detect_file(str(local))
                    sample_suggestions = [
                        s.suggestion for s in rep.suspicions[:5]
                        if s.suggestion
                    ]
                    results.append(PrecisionResult(
                        pmc_id=full_pmc, journal="", year=None,
                        file_name=name,
                        n_supplementary_files=len(file_names),
                        n_corruption_suspicions=len(rep.suspicions),
                        flagged=len(rep.suspicions) > 0,
                        sample_suggestions=sample_suggestions,
                    ))
                except Exception as exc:
                    results.append(PrecisionResult(
                        pmc_id=full_pmc, journal="", year=None,
                        file_name=name,
                        n_supplementary_files=len(file_names),
                        n_corruption_suspicions=0, flagged=False,
                        error=f"detect failed: {type(exc).__name__}",
                    ))

            if i % 25 == 0 or i == 1 or i == len(sampled):
                n_proc = sum(1 for r in results if r.error is None)
                n_flagged = sum(1 for r in results if r.flagged)
                elapsed = time.monotonic() - t0
                print(
                    f"[precision]   {i}/{len(sampled)}  with_suppl={n_with_suppl}  "
                    f"flagged={n_flagged}  elapsed={elapsed:.0f}s",
                    file=sys.stderr,
                )
                sys.stderr.flush()
                if checkpoint_path is not None and i % 25 == 0:
                    partial = PrecisionRun(
                        timestamp=datetime.now(UTC).isoformat(),
                        seed=seed, year_min=year_min, year_max=year_max,
                        n_sampled=i, n_processed=n_proc,
                        n_with_supplementaries=n_with_suppl,
                        n_flagged=n_flagged,
                        runtime_seconds=elapsed,
                        results=list(results),
                    )
                    with checkpoint_path.open("w") as fh:
                        json.dump(asdict(partial), fh)

    elapsed = time.monotonic() - t0
    n_proc = sum(1 for r in results if r.error is None)
    n_flagged = sum(1 for r in results if r.flagged)
    return PrecisionRun(
        timestamp=datetime.now(UTC).isoformat(),
        seed=seed, year_min=year_min, year_max=year_max,
        n_sampled=len(sampled), n_processed=n_proc,
        n_with_supplementaries=n_with_suppl,
        n_flagged=n_flagged,
        runtime_seconds=elapsed,
        results=results,
    )


def append_to_ledger(run: PrecisionRun, path: Path = LEDGER_PATH) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(json.dumps(asdict(run)) + "\n")


def format_precision_report(run: PrecisionRun) -> str:
    """Show flagged articles for manual triage. Each is either a Ziemann-missed
    finding or a false positive."""
    lines = [
        f"Precision check @ {run.timestamp} (seed={run.seed})",
        f"  Window: {run.year_min}-{run.year_max}",
        f"  Sampled: {run.n_sampled}",
        f"  Processed (no error): {run.n_processed}",
        f"  With supplementaries: {run.n_with_supplementaries}",
        f"  Flagged by detector: {run.n_flagged}",
        f"  Runtime: {run.runtime_seconds / 60:.1f} min",
        "",
        f"Flagged articles for manual triage:",
        f"{'PMC ID':<12} {'file':<40} {'n_susp':>7} {'sample suggestions'}",
    ]
    for r in run.results:
        if r.flagged:
            sugg = ", ".join(r.sample_suggestions[:3]) if r.sample_suggestions else ""
            lines.append(
                f"{r.pmc_id:<12} {r.file_name[:40]:<40} {r.n_corruption_suspicions:>7} {sugg}"
            )
    return "\n".join(lines)


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=200)
    p.add_argument("--year-min", type=int, default=2018)
    p.add_argument("--year-max", type=int, default=2020)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--rate", type=float, default=1.5)
    args = p.parse_args()

    checkpoint = PROJECT_ROOT / "results" / "precision_check_partial.json"
    run = run_precision_check(
        n_sample=args.n, year_min=args.year_min, year_max=args.year_max,
        seed=args.seed, rate_limit_per_second=args.rate,
        checkpoint_path=checkpoint,
    )
    print(format_precision_report(run))
    append_to_ledger(run)
    print(f"\nAppended to {LEDGER_PATH}")
