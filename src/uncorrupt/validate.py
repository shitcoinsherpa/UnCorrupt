"""Real-corpus validation of the corruption detector.

Downloads a random sample of confirmed-corrupted supplementary files from the
Ziemann 2021 corpus, runs the detector, and writes results to a JSONL ledger
under results/.

Per the real-data rule: this : not the synthetic smoke test : produces the
canonical recall number for the detector.
"""
from __future__ import annotations

import hashlib
import io
import json
import random
import time
import zipfile
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path

import httpx

# Two legitimate automated channels for PMC content:
#   1. Europe PMC REST: supplementaryFiles per article (preferred : flat ZIP).
#   2. NCBI OA service: oa.fcgi returns a tarball URL with the full article.
# Direct scraping of pmc.ncbi.nlm.nih.gov is rate-limited by PoW and is also
# explicitly disallowed by NCBI policy ("Download failed" response cites OA,
# OAI-PMH, FTP, and E-Utilities as the only permitted automated channels).
EPMC_SUPPL_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc_id}/supplementaryFiles"
NCBI_OA_FCGI = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id={pmc_id}"


def _fetch_supplementary_zip(pmc_id: str, client: httpx.Client) -> bytes:
    url = EPMC_SUPPL_URL.format(pmc_id=pmc_id)
    resp = client.get(url, follow_redirects=True, timeout=120)
    resp.raise_for_status()
    return resp.content


def _fetch_oa_tarball_url(pmc_id: str, client: httpx.Client) -> str | None:
    """Query NCBI OA service for the article's tarball URL. None if not in OA set."""
    import xml.etree.ElementTree as ET

    resp = client.get(NCBI_OA_FCGI.format(pmc_id=pmc_id), timeout=30)
    resp.raise_for_status()
    try:
        root = ET.fromstring(resp.text)
    except ET.ParseError:
        return None
    link = root.find(".//link[@format='tgz']")
    if link is None:
        return None
    url = link.get("href")
    if url and url.startswith("ftp://"):
        # Same path is served via HTTPS at the same host
        url = "https://" + url[len("ftp://"):]
    return url


def _extract_from_zip(zip_bytes: bytes, target_filename: str) -> bytes | None:
    """Return the bytes of `target_filename` inside the ZIP, or None if absent."""
    with zipfile.ZipFile(io.BytesIO(zip_bytes)) as z:
        for name in z.namelist():
            if name == target_filename or name.endswith("/" + target_filename):
                return z.read(name)
    return None


def _extract_from_tarball(tarball_bytes: bytes, target_filename: str) -> bytes | None:
    """Find target inside a (gzipped) tar archive."""
    import tarfile
    with tarfile.open(fileobj=io.BytesIO(tarball_bytes), mode="r:gz") as tf:
        for member in tf.getmembers():
            if not member.isfile():
                continue
            name = member.name.rsplit("/", 1)[-1]
            if name == target_filename:
                f = tf.extractfile(member)
                if f is not None:
                    return f.read()
    return None


# Backward-compat alias
_extract_target = _extract_from_zip

from .corpus import load_ziemann_2021_corpus
from .detector import detect_file

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_FILES_DIR = PROJECT_ROOT / "data" / "raw" / "ziemann_2021_corpus" / "files"
LEDGER_PATH = PROJECT_ROOT / "results" / "real_corpus_validation.jsonl"


@dataclass
class FileResult:
    pmc_id: str
    file_name: str
    file_size: int
    file_sha256: str
    rows: int
    columns: int
    identifier_columns: list[str]
    n_suspicions: int
    detected: bool
    error: str | None = None


@dataclass
class ValidationRun:
    timestamp: str
    seed: int
    n_sampled: int
    n_downloaded: int
    n_processed: int
    n_detected: int
    recall: float
    runtime_seconds: float
    results: list[FileResult] = field(default_factory=list)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def run_validation(
    n_sample: int | None = 100,
    seed: int = 42,
    rate_limit_per_second: float = 1.5,
    only_confirmed: bool = True,
    checkpoint_path: Path | None = None,
    progress_every: int = 20,
) -> ValidationRun:
    """Run real-corpus validation.

    Args:
        n_sample: how many entries to process; None means all.
        seed: RNG seed for sampling.
        rate_limit_per_second: 1.5 by default (well under EPMC's 10/s ceiling).
        only_confirmed: if True, only entries Ziemann manually confirmed.
            False also includes the tool-flagged-but-unconfirmed entries
            (useful as a soft precision proxy: high agreement = high precision).
        checkpoint_path: if provided, write results JSON every `progress_every`
            entries so a crash doesn't lose progress.
        progress_every: stderr progress + checkpoint frequency.
    """
    import sys

    rng = random.Random(seed)
    t0 = time.monotonic()

    entries = load_ziemann_2021_corpus()
    if only_confirmed:
        pool = [e for e in entries if e.confirmed_code and e.confirmed_code != "nan"]
    else:
        pool = entries
    if n_sample is None:
        sampled = pool
    else:
        sampled = rng.sample(pool, min(n_sample, len(pool)))
    total = len(sampled)
    print(f"[validate] corpus pool={len(pool)}, sampled={total}", file=sys.stderr)

    DATA_FILES_DIR.mkdir(parents=True, exist_ok=True)
    results: list[FileResult] = []
    delay = 1.0 / rate_limit_per_second

    headers = {
        "User-Agent": "uncorrupt-research/0.1 (https://github.com/shitcoinsherpa)",
    }
    with httpx.Client(headers=headers) as client:
        for i, entry in enumerate(sampled, start=1):
            if i % progress_every == 0 or i == 1 or i == total:
                elapsed = time.monotonic() - t0
                eta_s = elapsed / i * (total - i) if i > 0 else 0
                detected_so_far = sum(1 for r in results if not r.error and r.detected)
                processed_so_far = sum(1 for r in results if not r.error)
                rec = detected_so_far / processed_so_far if processed_so_far else 0
                print(
                    f"[validate] {i:>5}/{total}  elapsed={elapsed:.0f}s  eta={eta_s:.0f}s  "
                    f"detected={detected_so_far}/{processed_so_far} recall_so_far={rec:.3f}",
                    file=sys.stderr,
                )
                if checkpoint_path is not None and i % progress_every == 0:
                    partial = ValidationRun(
                        timestamp=datetime.now(UTC).isoformat(),
                        seed=seed,
                        n_sampled=i,
                        n_downloaded=sum(1 for r in results if r.file_size > 0),
                        n_processed=processed_so_far,
                        n_detected=detected_so_far,
                        recall=rec,
                        runtime_seconds=elapsed,
                        results=list(results),
                    )
                    with checkpoint_path.open("w") as fh:
                        json.dump(asdict(partial), fh)

            target_name = entry.affected_file_url.rsplit("/", 1)[-1]
            local_path = DATA_FILES_DIR / f"{entry.pmc_id}_{target_name}"

            if not local_path.exists():
                file_bytes: bytes | None = None
                fetch_errors: list[str] = []

                # Path 1: Europe PMC supplementaryFiles ZIP
                try:
                    zip_bytes = _fetch_supplementary_zip(entry.pmc_id, client)
                    time.sleep(delay)
                    try:
                        file_bytes = _extract_from_zip(zip_bytes, target_name)
                        if file_bytes is None:
                            fetch_errors.append(f"target not in epmc zip")
                    except zipfile.BadZipFile:
                        fetch_errors.append("epmc returned non-zip")
                except Exception as exc:
                    fetch_errors.append(f"epmc: {type(exc).__name__}")

                # Path 2: NCBI OA tarball (only if Path 1 didn't produce file)
                if file_bytes is None:
                    try:
                        tarball_url = _fetch_oa_tarball_url(entry.pmc_id, client)
                        if tarball_url is None:
                            fetch_errors.append("article not in OA set")
                        else:
                            tr = client.get(tarball_url, follow_redirects=True, timeout=180)
                            tr.raise_for_status()
                            file_bytes = _extract_from_tarball(tr.content, target_name)
                            if file_bytes is None:
                                fetch_errors.append(f"target not in OA tarball")
                            time.sleep(delay)
                    except Exception as exc:
                        fetch_errors.append(f"oa: {type(exc).__name__}: {exc}")

                if file_bytes is None:
                    results.append(FileResult(
                        pmc_id=entry.pmc_id, file_name=target_name,
                        file_size=0, file_sha256="", rows=0, columns=0,
                        identifier_columns=[], n_suspicions=0, detected=False,
                        error=" ; ".join(fetch_errors),
                    ))
                    continue
                local_path.write_bytes(file_bytes)

            try:
                report = detect_file(str(local_path))
                results.append(FileResult(
                    pmc_id=entry.pmc_id, file_name=target_name,
                    file_size=local_path.stat().st_size,
                    file_sha256=_sha256(local_path),
                    rows=report.rows_scanned,
                    columns=report.columns_scanned,
                    identifier_columns=report.identifier_columns,
                    n_suspicions=len(report.suspicions),
                    detected=len(report.suspicions) > 0,
                ))
            except Exception as exc:
                results.append(FileResult(
                    pmc_id=entry.pmc_id, file_name=target_name,
                    file_size=local_path.stat().st_size,
                    file_sha256=_sha256(local_path),
                    rows=0, columns=0,
                    identifier_columns=[], n_suspicions=0, detected=False,
                    error=f"detect failed: {type(exc).__name__}: {exc}",
                ))

    n_processed = sum(1 for r in results if r.error is None)
    n_detected = sum(1 for r in results if r.detected)
    recall = n_detected / n_processed if n_processed > 0 else 0.0

    return ValidationRun(
        timestamp=datetime.now(UTC).isoformat(),
        seed=seed,
        n_sampled=len(sampled),
        n_downloaded=sum(1 for r in results if r.file_size > 0),
        n_processed=n_processed,
        n_detected=n_detected,
        recall=recall,
        runtime_seconds=time.monotonic() - t0,
        results=results,
    )


def _format_report(run: ValidationRun) -> str:
    lines = [
        f"Real-corpus validation @ {run.timestamp} (seed={run.seed})",
        "-" * 100,
        f"  Sampled: {run.n_sampled}",
        f"  Downloaded: {run.n_downloaded}",
        f"  Processed: {run.n_processed}",
        f"  Detected: {run.n_detected}",
        f"  Recall: {run.recall:.3f}  (detected / processed; ground truth = all confirmed corrupted)",
        f"  Runtime: {run.runtime_seconds:.1f}s",
        "",
        f"{'PMC ID':<12} {'size':>7} {'rows':>5} {'cols':>4} {'ID cols':<28} {'susp':>4} {'flagged'}",
    ]
    for r in run.results:
        if r.error:
            lines.append(f"{r.pmc_id:<12} ERROR: {r.error}")
        else:
            id_cols_s = (",".join(r.identifier_columns)[:26]) or "-"
            lines.append(
                f"{r.pmc_id:<12} {r.file_size:>7} {r.rows:>5} {r.columns:>4} "
                f"{id_cols_s:<28} {r.n_suspicions:>4} {'YES' if r.detected else 'no'}"
            )
    return "\n".join(lines)


def append_to_ledger(run: ValidationRun, ledger_path: Path = LEDGER_PATH) -> None:
    ledger_path.parent.mkdir(parents=True, exist_ok=True)
    with ledger_path.open("a") as f:
        f.write(json.dumps(asdict(run)) + "\n")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=100,
                   help="number to sample (use 0 or 'all' for full corpus)")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--all", action="store_true", help="process the full corpus")
    p.add_argument("--include-unconfirmed", action="store_true",
                   help="also include tool-flagged-but-unconfirmed entries")
    p.add_argument("--rate", type=float, default=1.5, help="requests per second")
    args = p.parse_args()

    checkpoint = PROJECT_ROOT / "results" / "real_corpus_validation_partial.json"
    run = run_validation(
        n_sample=None if args.all else (None if args.n == 0 else args.n),
        seed=args.seed,
        rate_limit_per_second=args.rate,
        only_confirmed=not args.include_unconfirmed,
        checkpoint_path=checkpoint,
    )
    print(_format_report(run))
    append_to_ledger(run)
    print(f"\nAppended to {LEDGER_PATH}")
