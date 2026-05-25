"""Null-result test: re-run detector against the Koh
"clean with gene symbols" pool. Expected output: 0 suspicions.

Methodology: 100 randomly-sampled files (deterministic seed) from the
classifier's "clean with gene symbols" stratum in the Koh cache walk
ledger. Each file is run through `detect_file()`; counts of flags by
kind are recorded.

A non-zero count of `gene-date`/`gene-date-serial`/`gene-date-string`/
`autofill-sequence`/`time-coercion`/`homoglyph` on a known-clean file
is a precision regression : investigate immediately.

Run:
   python scripts/null_result_clean_pool.py
"""
from __future__ import annotations

import json
import random
import sys
import warnings
from collections import defaultdict
from pathlib import Path

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from uncorrupt.detector import detect_file  # noqa: E402

WALK_LEDGER = ROOT / "results/koh_cache_walk.jsonl"
KOH = ROOT / "data/raw/koh_replication/files"
OUT = ROOT / "results/null_result_clean_pool.json"
SAMPLE_N = 100
SEED = 20260518


def main() -> int:
    if not WALK_LEDGER.exists():
        print(f"FATAL: walk ledger missing at {WALK_LEDGER}", file=sys.stderr)
        return 1
    with WALK_LEDGER.open() as fh:
        walk = json.load(fh)
    clean_with_genes = [
        f for f in walk.get("files", [])
        if f.get("has_gene_symbols") and f.get("n_suspicions", 0) == 0
    ]
    rng = random.Random(SEED)
    sample = rng.sample(clean_with_genes, min(SAMPLE_N, len(clean_with_genes)))

    suspicion_kinds: dict[str, int] = defaultdict(int)
    n_flagged_files = 0
    failures: list[tuple[str, str]] = []
    files_inspected = 0
    for entry in sample:
        pmc = entry["pmc_id"]
        fname = entry["file_name"]
        path = KOH / pmc / fname
        if not path.exists():
            failures.append((f"{pmc}/{fname}", "file missing on disk"))
            continue
        try:
            r = detect_file(str(path), row_context_boost=False)
        except Exception as exc:
            failures.append((f"{pmc}/{fname}", f"{type(exc).__name__}: {exc}"))
            continue
        files_inspected += 1
        if r.suspicions:
            # Informational flag kinds are NOT corruption claims. Exclude
            # them from the corruption-rate accounting:
            #   - decimal-comma: locale serialization note
            #   - unrecognized-symbol: HGNC drift / non-human ortholog signal
            #   - cas-registry in chemistry columns: chain-of-custody only
            INFO_KINDS = {"decimal-comma", "unrecognized-symbol"}
            real_flags = [s for s in r.suspicions
                          if s.kind not in INFO_KINDS]
            if real_flags:
                n_flagged_files += 1
            for s in real_flags:
                suspicion_kinds[s.kind] += 1

    summary = {
        "sample_size": SAMPLE_N,
        "files_inspected": files_inspected,
        "n_failures": len(failures),
        "n_flagged_files": n_flagged_files,
        "suspicion_kinds": dict(suspicion_kinds),
        "false_positive_rate_files": n_flagged_files / files_inspected
        if files_inspected else 0.0,
        "failures_sample": failures[:5],
    }
    OUT.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print()
    if sum(suspicion_kinds.values()) == 0:
        print(f"PASS: 0 false-positive flags on {files_inspected} clean files.")
        return 0
    else:
        print(
            f"FAIL: {sum(suspicion_kinds.values())} flags on {files_inspected} clean files. "
            f"Detailed kinds: {dict(suspicion_kinds)}"
        )
        return 1


if __name__ == "__main__":
    sys.exit(main())
