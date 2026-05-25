"""Refresh the HGNC bulk-download snapshot used by registries.py.

HGNC publishes the complete-set bulk download daily at:
    https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt

UnCorrupt depends on this file for current-symbol membership, prev-symbol
rename map, and alias-to-current map. The detector is only as accurate as
its HGNC snapshot : symbols added or renamed *after* our snapshot will
silently miss. This script:

  1. Computes SHA256 of the current local snapshot (if any).
  2. Fetches the upstream file to a date-stamped path.
  3. Compares SHA256s. Drift → emit "stale" notice with the date delta.
  4. Records a manifest entry under data/raw/registries/MANIFEST.jsonl.

Run weekly as a cron (Pillar 5 : continuous validation) or before any
external claim. The detector imports the most recently dated file in
data/raw/registries/ matching hgnc_complete_set_*.tsv; pointing imports
at the new file is a single-line update to registries.py once the fetch
succeeds.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from datetime import UTC, date, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REG_DIR = ROOT / "data/raw/registries"
MANIFEST = REG_DIR / "MANIFEST.jsonl"
HGNC_URL = (
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/"
    "hgnc_complete_set.txt"
)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def _current_snapshot() -> Path | None:
    candidates = sorted(REG_DIR.glob("hgnc_complete_set_*.tsv"))
    return candidates[-1] if candidates else None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Report drift but do not download.",
    )
    parser.add_argument(
        "--check-only", action="store_true",
        help="Print whether the local snapshot is older than 30 days; exit non-zero if stale.",
    )
    args = parser.parse_args()

    REG_DIR.mkdir(parents=True, exist_ok=True)
    current = _current_snapshot()
    if current is None:
        print("No local HGNC snapshot found.", file=sys.stderr)
        if args.check_only:
            return 2
    else:
        # Parse YYYY-MM-DD out of the filename.
        date_part = current.stem.replace("hgnc_complete_set_", "")
        try:
            snapshot_date = date.fromisoformat(date_part)
        except ValueError:
            snapshot_date = None
        if snapshot_date is not None:
            age = (date.today() - snapshot_date).days
            print(f"Local snapshot: {current.name} ({age} days old)")
            if args.check_only:
                if age > 30:
                    print(f"STALE: local snapshot is {age} days old.")
                    return 1
                print("Fresh.")
                return 0
        else:
            print(f"Local snapshot: {current.name} (date unparseable)")

    if args.dry_run:
        print("Dry-run: skipping download.")
        return 0

    today_iso = date.today().isoformat()
    target = REG_DIR / f"hgnc_complete_set_{today_iso}.tsv"
    print(f"Fetching {HGNC_URL} → {target}")
    try:
        urllib.request.urlretrieve(HGNC_URL, target)  # noqa: S310
    except Exception as exc:
        print(f"FETCH FAILED: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 3

    new_sha = _sha256(target)
    old_sha = _sha256(current) if current is not None else None
    if old_sha == new_sha:
        print(f"SHA256 unchanged ({new_sha[:12]}…). HGNC upstream not yet refreshed.")
        # Remove the redundant date-stamped copy.
        target.unlink()
        return 0
    print(f"NEW HGNC snapshot : SHA256: {new_sha[:12]}…")
    if old_sha is not None:
        print(f"  Previous: {old_sha[:12]}…")

    # Append to MANIFEST.jsonl
    entry = {
        "filename": target.name,
        "downloaded_at": datetime.now(UTC).isoformat(),
        "sha256": new_sha,
        "url": HGNC_URL,
        "previous_sha256": old_sha,
    }
    with MANIFEST.open("a") as fh:
        fh.write(json.dumps(entry) + "\n")
    print(f"Manifest updated: {MANIFEST}")
    print()
    print("Next: update `data/raw/registries/` import path in registries.py")
    print(f"if needed (it should auto-discover the newest dated file).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
