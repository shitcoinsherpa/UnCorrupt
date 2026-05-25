"""Fetch the Ziemann 2021 S2 supplementary files we don't have cached yet.

The S2 article list has 5,086 confirmed-corruption entries; we have 840
files cached locally (1 file per row of the S2 sample-set used by
`uncorrupt.cell_level_validation`). This script fetches the rest via:

  1. Europe PMC supplementary-files ZIP endpoint (when the PMC ID is OA)
  2. Direct publisher URL from the `Affected_file` column when EPMC fails

The fetched files go to `data/raw/ziemann_2021_corpus/files/` with the
naming convention `{pmc_id}_{filename}` matching the existing cache.

Run:
   python scripts/fetch_missing_ziemann_2021.py --limit 500 --rate 1.5
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import sys
import time
import warnings
import zipfile
from datetime import UTC, datetime
from pathlib import Path
from urllib.parse import urlparse

warnings.filterwarnings("ignore")

import httpx  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from uncorrupt.corpus import load_ziemann_2021_corpus  # noqa: E402

CACHE = ROOT / "data/raw/ziemann_2021_corpus/files"
LEDGER = ROOT / "data/raw/ziemann_2021_corpus/fetch_ledger_missing.json"

EPMC_SUPPL = "https://www.ebi.ac.uk/europepmc/webservices/rest/{full_pmc_id}/supplementaryFiles"

USER_AGENT = "UnCorrupt/0.7 (research; github.com/shitcoinsherpa)"


def _safe_filename(name_or_url: str) -> str:
    raw = name_or_url.rsplit("/", 1)[-1] or "supplementary.xlsx"
    raw = raw.split("?", 1)[0]
    return re.sub(r"[^A-Za-z0-9._-]", "_", raw)[:140]


def _try_epmc_zip(full_pmc_id: str, client: httpx.Client) -> bytes | None:
    try:
        r = client.get(EPMC_SUPPL.format(full_pmc_id=full_pmc_id), timeout=60)
    except Exception:
        return None
    if r.status_code != 200 or len(r.content) < 500:
        return None
    return r.content


def _try_publisher(url: str, client: httpx.Client) -> bytes | None:
    if not url or not url.startswith("http"):
        return None
    try:
        r = client.get(url, timeout=60)
    except Exception:
        return None
    if r.status_code != 200 or len(r.content) < 256:
        return None
    return r.content


def _iter_xlsx_in_zip(zip_bytes: bytes):
    try:
        z = zipfile.ZipFile(io.BytesIO(zip_bytes))
    except zipfile.BadZipFile:
        return
    for info in z.infolist():
        name_lower = info.filename.lower()
        if not (name_lower.endswith(".xlsx") or name_lower.endswith(".xls")):
            continue
        try:
            yield info.filename, z.read(info.filename)
        except Exception:
            continue


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--limit", type=int, default=None,
                   help="Cap on number of new fetches per run (default: all missing).")
    p.add_argument("--rate", type=float, default=1.5)
    args = p.parse_args()

    entries = load_ziemann_2021_corpus()
    cached_pmcs = {p.name.split("_", 1)[0] for p in CACHE.glob("PMC*")}
    missing = [e for e in entries if e.pmc_id not in cached_pmcs]
    print(f"[fetch-z21] {len(missing)} missing (of {len(entries)} total)",
          flush=True)

    if args.limit:
        missing = missing[: args.limit]

    delay = 1.0 / max(args.rate, 0.1)
    headers = {"User-Agent": USER_AGENT, "Accept": "*/*"}
    success: list[dict] = []
    failures: list[dict] = []

    with httpx.Client(headers=headers, follow_redirects=True, timeout=60) as client:
        for i, e in enumerate(missing, 1):
            full_pmc = e.pmc_id
            if not full_pmc.startswith("PMC"):
                full_pmc = "PMC" + full_pmc

            # Try Europe PMC first (returns a ZIP of all supplementary files)
            zip_bytes = _try_epmc_zip(full_pmc, client)
            time.sleep(delay)
            saved = False
            if zip_bytes:
                for name, content in _iter_xlsx_in_zip(zip_bytes):
                    safe = _safe_filename(name)
                    out_path = CACHE / f"{full_pmc}_{safe}"
                    if not out_path.exists():
                        out_path.write_bytes(content)
                    saved = True
                    success.append({
                        "pmc_id": full_pmc,
                        "source": "epmc",
                        "filename": out_path.name,
                        "size_bytes": out_path.stat().st_size,
                        "sha256": hashlib.sha256(content).hexdigest(),
                    })

            if not saved and e.affected_file_url:
                # Fall back to direct publisher URL
                pub_bytes = _try_publisher(e.affected_file_url, client)
                time.sleep(delay)
                if pub_bytes:
                    safe = _safe_filename(urlparse(e.affected_file_url).path)
                    out_path = CACHE / f"{full_pmc}_{safe}"
                    if not out_path.exists():
                        out_path.write_bytes(pub_bytes)
                    saved = True
                    success.append({
                        "pmc_id": full_pmc,
                        "source": "publisher",
                        "filename": out_path.name,
                        "url": e.affected_file_url,
                        "size_bytes": out_path.stat().st_size,
                        "sha256": hashlib.sha256(pub_bytes).hexdigest(),
                    })

            if not saved:
                failures.append({"pmc_id": full_pmc,
                                  "url": e.affected_file_url})

            if i % 50 == 0:
                print(f"[fetch-z21] {i}/{len(missing)} "
                      f"{len(success)} ok, {len(failures)} fail",
                      flush=True)

    LEDGER.parent.mkdir(parents=True, exist_ok=True)
    LEDGER.write_text(json.dumps({
        "timestamp": datetime.now(UTC).isoformat(),
        "n_attempted": len(missing),
        "n_success": len(success),
        "n_failures": len(failures),
        "successful": success,
        "failed": failures,
    }, indent=2))

    print()
    print(f"[fetch-z21] DONE: {len(success)} success, {len(failures)} failures")
    print(f"           ledger: {LEDGER}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
