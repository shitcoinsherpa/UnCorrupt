"""Build a SHA256 manifest of every cached corpus file.

Writes `data/raw/CORPUS_MANIFEST.json`. The manifest pins every file we used
to derive our headline numbers : anyone re-fetching the same files from
Europe PMC + NCBI can verify their cache bit-matches ours by replaying the
same hashes.

Output schema:
    {
        "generated": "ISO-8601 UTC timestamp",
        "ziemann_2021": {
            "<filename>": {
                "sha256": "...64 hex chars...",
                "size_bytes": int,
                "pmc_id": "PMC1234567"
            },
            ...
        },
        "koh_replication": {
            "<PMC_ID>/<filename>": {
                "sha256": "...",
                "size_bytes": int,
                "journal": "<from EPMC manifest>",
                "source_url": "...EPMC URL if available..."
            },
            ...
        }
    }

Usage:
    python scripts/build_corpus_manifest.py
    python scripts/build_corpus_manifest.py --verify   # check current cache matches manifest

Pillar 4 (FAIR data sharing): without this manifest, a reproduction attempt
can claim "I downloaded the corpus" without proving the bits match.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import UTC, datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ZIEMANN_DIR = PROJECT_ROOT / "data" / "raw" / "ziemann_2021_corpus" / "files"
KOH_DIR = PROJECT_ROOT / "data" / "raw" / "koh_replication" / "files"
MANIFEST_PATH = PROJECT_ROOT / "data" / "raw" / "CORPUS_MANIFEST.json"


def _sha256_of(path: Path, chunk: int = 1024 * 1024) -> tuple[str, int]:
    h = hashlib.sha256()
    size = 0
    with path.open("rb") as f:
        while True:
            block = f.read(chunk)
            if not block:
                break
            h.update(block)
            size += len(block)
    return h.hexdigest(), size


def _hash_one(arg: tuple[str, Path]) -> tuple[str, Path, str, int]:
    """Worker: (label, path) -> (label, path, sha256, size)."""
    label, path = arg
    digest, size = _sha256_of(path)
    return label, path, digest, size


def _ziemann_jobs() -> list[tuple[str, Path]]:
    if not ZIEMANN_DIR.exists():
        return []
    return [("ziemann_2021", p) for p in sorted(ZIEMANN_DIR.iterdir()) if p.is_file()]


def _koh_jobs() -> list[tuple[str, Path]]:
    if not KOH_DIR.exists():
        return []
    out: list[tuple[str, Path]] = []
    for pmc_dir in sorted(KOH_DIR.iterdir()):
        if not pmc_dir.is_dir():
            continue
        for p in sorted(pmc_dir.iterdir()):
            if p.name == "_manifest.json":
                continue
            if p.is_file():
                out.append(("koh_replication", p))
    return out


def _load_koh_url_index() -> dict[str, str]:
    """Map (pmc_id, filename) -> source URL from per-article manifest jsons.

    Falls back gracefully when manifests are missing the URL field.
    """
    out: dict[str, str] = {}
    if not KOH_DIR.exists():
        return out
    for pmc_dir in KOH_DIR.iterdir():
        if not pmc_dir.is_dir():
            continue
        mp = pmc_dir / "_manifest.json"
        if not mp.exists():
            continue
        try:
            data = json.loads(mp.read_text())
        except Exception:
            continue
        # Manifest schema: {"files": [name, ...], "error": optional}
        # URLs aren't captured in the current manifest, but we can reconstruct
        # the EPMC supplementaryFiles URL from the PMC ID.
        pmc_id = pmc_dir.name
        for fname in data.get("files", []):
            key = f"{pmc_id}/{fname}"
            out[key] = (
                f"https://europepmc.org/article/PMC/{pmc_id.removeprefix('PMC')}"
                f"/supplementaryFiles/{fname}"
            )
    return out


def build_manifest(workers: int | None = None) -> dict:
    """Compute SHA256 for every cached file. Returns the manifest dict."""
    jobs = _ziemann_jobs() + _koh_jobs()
    if not jobs:
        print("[manifest] no corpus files found : nothing to manifest",
              file=sys.stderr)
        return {"generated": datetime.now(UTC).isoformat(),
                "ziemann_2021": {}, "koh_replication": {}}

    print(f"[manifest] hashing {len(jobs):,} files "
          f"(ziemann_2021={len(_ziemann_jobs()):,}, koh_replication={len(_koh_jobs()):,})",
          file=sys.stderr)
    koh_urls = _load_koh_url_index()

    out: dict = {
        "generated": datetime.now(UTC).isoformat(),
        "ziemann_2021": {},
        "koh_replication": {},
    }
    done = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for future in as_completed([ex.submit(_hash_one, j) for j in jobs]):
            label, path, digest, size = future.result()
            done += 1
            if done % 500 == 0 or done == len(jobs):
                print(f"[manifest]   {done:,}/{len(jobs):,}", file=sys.stderr,
                      flush=True)
            if label == "ziemann_2021":
                pmc_prefix = path.name.split("_", 1)[0]
                out["ziemann_2021"][path.name] = {
                    "sha256": digest,
                    "size_bytes": size,
                    "pmc_id": pmc_prefix,
                }
            else:
                pmc_id = path.parent.name
                rel = f"{pmc_id}/{path.name}"
                out["koh_replication"][rel] = {
                    "sha256": digest,
                    "size_bytes": size,
                    "source_url": koh_urls.get(rel, ""),
                }

    return out


def verify_manifest() -> int:
    """Recompute SHA256s and compare against the stored manifest.

    Returns the number of mismatches (0 = perfect verification).
    """
    if not MANIFEST_PATH.exists():
        print(f"[verify] no manifest at {MANIFEST_PATH}", file=sys.stderr)
        return -1
    manifest = json.loads(MANIFEST_PATH.read_text())
    mismatches = 0
    missing = 0
    matched = 0

    for fname, entry in manifest.get("ziemann_2021", {}).items():
        p = ZIEMANN_DIR / fname
        if not p.exists():
            missing += 1
            continue
        actual, _ = _sha256_of(p)
        if actual != entry["sha256"]:
            mismatches += 1
            print(f"  MISMATCH ziemann_2021/{fname}", file=sys.stderr)
        else:
            matched += 1

    for rel, entry in manifest.get("koh_replication", {}).items():
        p = KOH_DIR / rel
        if not p.exists():
            missing += 1
            continue
        actual, _ = _sha256_of(p)
        if actual != entry["sha256"]:
            mismatches += 1
            print(f"  MISMATCH koh_replication/{rel}", file=sys.stderr)
        else:
            matched += 1

    print(f"[verify] matched={matched:,}  mismatched={mismatches:,}  "
          f"missing={missing:,}", file=sys.stderr)
    return mismatches


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--verify", action="store_true",
                   help="recompute hashes and compare to stored manifest")
    p.add_argument("--workers", type=int, default=None,
                   help="parallel worker count (default: CPU count)")
    args = p.parse_args()

    if args.verify:
        rc = verify_manifest()
        sys.exit(0 if rc == 0 else 1)

    manifest = build_manifest(workers=args.workers)
    MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2, sort_keys=True))
    n_z = len(manifest["ziemann_2021"])
    n_k = len(manifest["koh_replication"])
    total_size = sum(e["size_bytes"] for e in manifest["ziemann_2021"].values()) + \
                 sum(e["size_bytes"] for e in manifest["koh_replication"].values())
    print(f"\n[manifest] wrote {MANIFEST_PATH}", file=sys.stderr)
    print(f"  ziemann_2021     : {n_z:>6,} files", file=sys.stderr)
    print(f"  koh_replication  : {n_k:>6,} files", file=sys.stderr)
    print(f"  total size       : {total_size / 1024 / 1024 / 1024:>6.2f} GB",
          file=sys.stderr)


if __name__ == "__main__":
    main()
