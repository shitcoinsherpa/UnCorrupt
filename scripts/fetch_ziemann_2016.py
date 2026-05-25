"""Fetch a fresh sample of files from Ziemann 2016 S1 (Additional file 1).

This is independent ground truth : a corpus from a different paper than
Ziemann 2021, never touched by this project. Each row in S1 contains:

  - Supplementary file URL (publisher-direct, not PMC)
  - Example Gene Name Conversion (the documented corruption: "MARCH9", "2005-09-01", etc.)
  - Journal / Year / Confirmed status / PubmedID

We sample `--n` random `Confirmed` rows from publishers we can reliably
fetch (Springer static-content, PLoS journals, CSH journals). Files cached
to `data/raw/ziemann_2016_corpus/files/<pubmed_id>__<safe-filename>`.

Usage:
    python scripts/fetch_ziemann_2016.py --n 30 --seed 42
"""
from __future__ import annotations

import argparse
import hashlib
import json
import random
import re
import sys
import time
import warnings
from datetime import UTC, datetime
from pathlib import Path
from urllib.parse import urlparse

warnings.filterwarnings("ignore")

import httpx
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
S1_PATH = PROJECT_ROOT / "data/raw/ziemann_2016_corpus/Additional_file_1.xlsx"
OUT_DIR = PROJECT_ROOT / "data/raw/ziemann_2016_corpus/files"
LEDGER = PROJECT_ROOT / "data/raw/ziemann_2016_corpus/fetch_ledger.jsonl"

# Hosts we can reliably fetch (open-access publishers)
RELIABLE_HOSTS = {
    "static-content.springer.com",
    "journals.plos.org",
    "dx.plos.org",
    "genome.cshlp.org",
    "genesdev.cshlp.org",
    "rnajournal.cshlp.org",
    "hmg.oxfordjournals.org",
    "nar.oxfordjournals.org",
    "mbe.oxfordjournals.org",
    "bioinformatics.oxfordjournals.org",
    "gbe.oxfordjournals.org",
    "dnaresearch.oxfordjournals.org",
}

USER_AGENT = (
    "UnCorrupt/0.4 (research; github.com/shitcoinsherpa)"
)


def _safe_filename(url: str) -> str:
    parsed = urlparse(url)
    raw = parsed.path.rsplit("/", 1)[-1] or "supplementary.xlsx"
    return re.sub(r"[^A-Za-z0-9._-]", "_", raw)[:120]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=30)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--rate", type=float, default=2.0, help="requests/sec")
    args = p.parse_args()

    df = pd.read_excel(S1_PATH, dtype=object)
    # Filter to confirmed corruption on reliable hosts
    df = df[df["Confirmed"] == "Confirmed"].copy()
    df["host"] = df["Supplementary file URL"].apply(
        lambda u: urlparse(str(u)).netloc
    )
    df = df[df["host"].isin(RELIABLE_HOSTS)]
    print(f"[fetch] {len(df)} confirmed entries on reliable hosts", file=sys.stderr)

    rng = random.Random(args.seed)
    sampled = df.sample(args.n, random_state=args.seed).reset_index(drop=True)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    LEDGER.parent.mkdir(parents=True, exist_ok=True)
    delay = 1.0 / max(args.rate, 0.1)

    headers = {
        "User-Agent": USER_AGENT,
        "Accept": "*/*",
    }
    success: list[dict] = []
    failures: list[dict] = []
    with httpx.Client(headers=headers, follow_redirects=True, timeout=30) as client:
        for i, row in sampled.iterrows():
            url = str(row["Supplementary file URL"])
            pubmed = str(row["PubmedID"])
            safe = _safe_filename(url)
            out_path = OUT_DIR / f"{pubmed}__{safe}"
            if out_path.exists():
                print(f"[fetch] {i+1}/{len(sampled)} CACHED  {pubmed}  {safe}",
                      file=sys.stderr)
                size = out_path.stat().st_size
            else:
                try:
                    r = client.get(url)
                    if r.status_code != 200:
                        failures.append({"pubmed": pubmed, "url": url,
                                          "status": r.status_code, "len": 0})
                        print(f"[fetch] {i+1}/{len(sampled)} FAIL {r.status_code}  "
                              f"{pubmed}", file=sys.stderr)
                        time.sleep(delay)
                        continue
                    out_path.write_bytes(r.content)
                    size = len(r.content)
                except Exception as exc:
                    failures.append({"pubmed": pubmed, "url": url,
                                      "status": "error", "error": str(exc)[:100]})
                    print(f"[fetch] {i+1}/{len(sampled)} ERROR  {pubmed}  "
                          f"{type(exc).__name__}", file=sys.stderr)
                    time.sleep(delay)
                    continue
                time.sleep(delay)
            # Hash the file
            h = hashlib.sha256(out_path.read_bytes()).hexdigest()
            success.append({
                "pubmed_id": pubmed,
                "url": url,
                "filename": out_path.name,
                "size_bytes": size,
                "sha256": h,
                "journal": str(row["Journal"]),
                "year": int(row["Year Published"]) if pd.notna(row["Year Published"]) else None,
                "documented_corruption": str(row["Example Gene Name Conversion"]),
                "other_identifiers": str(row["Other accession numbers/identifying info"]),
            })

    # Write ledger
    with LEDGER.open("w") as f:
        f.write(json.dumps({
            "timestamp": datetime.now(UTC).isoformat(),
            "n_requested": len(sampled),
            "n_success": len(success),
            "n_failures": len(failures),
            "successful": success,
            "failed": failures,
        }, indent=2))
    print(f"\n[fetch] DONE: {len(success)}/{len(sampled)} success, "
          f"{len(failures)} failures")
    print(f"  files in: {OUT_DIR}")
    print(f"  ledger: {LEDGER}")


if __name__ == "__main__":
    main()
