"""Augmented confidence calibration using xref-cross-validation to
recover both TP and FP anchors from a large heterogeneous corpus.

Problem with `calibrate_confidence.py`: the labeled audit corpus
(478 files, 55,622 cells) had 0 false positives by construction : 
every flagged cell was confirmed via openpyxl. Without FP examples,
the calibration data shows 100% precision per bin but tells us nothing
about how stated confidence relates to *actual* false-positive rates.

This script builds the missing FP-anchor stratum:

  1. Walk a stratified random sample of the Koh cache (4,000 files
     by default : covers all three classifier strata: with-corruption,
     clean-with-genes, no-genes).
  2. Run the detector on each file.
  3. For each flag, cross-validate against the row-context xref:
        - CORROBORATED: row contains an Entrez/Ensembl/RefSeq/UniProt/HGNC
          identifier that resolves to the suggested gene (or a HGNC alias
          equivalent) → confirmed TP.
        - CONTRADICTED: row contains a resolvable external ID that resolves
          to a *different* gene → candidate FP.
        - INCONCLUSIVE: no resolvable external ID in the row → cannot
          adjudicate from xref alone.
  4. The CORROBORATED + CONTRADICTED stratum is the labeled subset
     (xref ground truth). Compute Brier/ECE/Wilson/Platt on that subset
     and on the broader labeled audit corpus combined.
  5. Emit per-bin precision tables and a reliability diagram with
     two series: audit-only (existing) and xref-augmented.

This is "real-data calibration" : the labels come from independent
external identifiers in the same row, not from a re-run of the detector.

Run:
   python scripts/calibrate_with_fp_anchors.py
   python scripts/calibrate_with_fp_anchors.py --sample 1000   # smaller
"""
from __future__ import annotations

import argparse
import json
import math
import random
import sys
import time
import warnings
from collections import defaultdict
from pathlib import Path

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from uncorrupt.detector import detect_file  # noqa: E402
from uncorrupt.xref_lookup import load_xref_index  # noqa: E402
from uncorrupt.corpus import load_hgnc  # noqa: E402

WALK_LEDGER = ROOT / "results/koh_cache_walk.jsonl"
KOH = ROOT / "data/raw/koh_replication/files"
OUT_JSONL = ROOT / "results/calibration_fp_anchors.jsonl"
OUT_REPORT = ROOT / "results/calibration_with_fp_anchors.md"

SEED = 20260518


def _expand_with_hgnc_aliases(hgnc, suggestions: set[str]) -> set[str]:
    """Expand a set of gene-name suggestions to include every HGNC
    equivalent (prev_symbol, alias) : case-insensitively. So `OCT4`
    expands to include `POU5F1` because POU5F1 has alias `Oct4` (in
    HGNC's canonical mixed case)."""
    expanded = set(suggestions)
    for s in list(suggestions):
        # Forward: stale or alias → current canonical
        canonical = hgnc.resolve(s)
        if canonical:
            expanded.add(canonical)
        # Reverse: current → all known prev_symbols + aliases pointing to it
        for old, new in hgnc.prev_symbol_to_current.items():
            if new == s or new == canonical:
                expanded.add(old)
        for alias, new in hgnc.alias_to_current.items():
            if new == s or new == canonical:
                expanded.add(alias)
    return expanded


def _wilson_ci(n_correct: int, n_total: int) -> tuple[float, float]:
    if n_total == 0:
        return (0.0, 1.0)
    p = n_correct / n_total
    z = 1.96
    denom = 1 + z * z / n_total
    centre = (p + z * z / (2 * n_total)) / denom
    half = z * math.sqrt(p * (1 - p) / n_total + z * z / (4 * n_total * n_total)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def _brier(scores: list[float], outcomes: list[int]) -> float:
    if not scores:
        return 0.0
    return sum((s - o) ** 2 for s, o in zip(scores, outcomes)) / len(scores)


def _ece(scores: list[float], outcomes: list[int]) -> float:
    if not scores:
        return 0.0
    bins: dict[float, list[int]] = defaultdict(list)
    for s, o in zip(scores, outcomes):
        bins[s].append(o)
    n = len(scores)
    return sum(
        (len(outs) / n) * abs(conf - sum(outs) / len(outs))
        for conf, outs in bins.items()
    )


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample", type=int, default=4000)
    p.add_argument("--per-pool", action="store_true",
                   help="Stratify the sample equally across the three Koh pools.")
    p.add_argument("--ledger", type=str, default=str(WALK_LEDGER),
                   help="Path to the walk ledger to sample from. "
                        "Defaults to results/koh_cache_walk.jsonl.")
    p.add_argument("--sample-offset", type=int, default=0,
                   help="Skip the first N files of the sampled list "
                        "(for parallel sharding).")
    p.add_argument("--sample-limit", type=int, default=None,
                   help="Process at most N files starting from --sample-offset.")
    p.add_argument("--out", type=str, default=str(OUT_JSONL),
                   help="Output JSONL path (override for sharded runs).")
    args = p.parse_args()

    ledger_path = Path(args.ledger)
    if not ledger_path.exists():
        print(f"FATAL: walk ledger missing at {ledger_path}", file=sys.stderr)
        return 1

    print(f"Loading xref + HGNC...", flush=True)
    xref = load_xref_index()
    hgnc = load_hgnc()

    with ledger_path.open() as fh:
        walk = json.load(fh)
    files = walk.get("files", [])

    # Stratify
    by_pool: dict[str, list] = defaultdict(list)
    for f in files:
        if f.get("error"):
            continue
        if f.get("has_gene_symbols") and f.get("n_suspicions", 0) > 0:
            pool = "with_corruption"
        elif f.get("has_gene_symbols"):
            pool = "clean_with_genes"
        else:
            pool = "no_genes"
        by_pool[pool].append(f)

    rng = random.Random(SEED)
    if args.per_pool:
        per = args.sample // 3
        sample = (
            rng.sample(by_pool["with_corruption"], min(per, len(by_pool["with_corruption"])))
            + rng.sample(by_pool["clean_with_genes"], min(per, len(by_pool["clean_with_genes"])))
            + rng.sample(by_pool["no_genes"], min(per, len(by_pool["no_genes"])))
        )
    else:
        sample = rng.sample(files, min(args.sample, len(files)))

    print(f"Pool sizes: {dict((k, len(v)) for k, v in by_pool.items())}", flush=True)
    print(f"Sampled {len(sample)} files for cross-validation", flush=True)

    # Slice for parallel sharding : random.Random(SEED) is deterministic so
    # all workers produce the same `sample` list; each worker picks its slice.
    if args.sample_offset or args.sample_limit is not None:
        end = (args.sample_offset + args.sample_limit) if args.sample_limit is not None else len(sample)
        sample = sample[args.sample_offset:end]
        print(f"Slicing sample to [{args.sample_offset}:{end}] -> {len(sample)} files",
              flush=True)
    OUT_JSONL_LOCAL = Path(args.out)

    # Pathological-file guard: openpyxl can hang for hours on certain
    # xlsx files even at 16 MB (empirically caught: PMC9246204/
    # pone.0250137.s001.xlsx, 16 MB, openpyxl stuck > 1 hour). The
    # python-calamine fast loader kicks in at 30 MB; below that we use
    # openpyxl and risk these hangs. Lower the guard to 10 MB for
    # calibration walks where we'd rather sample more files than wait
    # on one pathological one.
    MAX_BYTES = 10 * 1024 * 1024
    OUT_JSONL_LOCAL.parent.mkdir(parents=True, exist_ok=True)
    n_corroborated = 0
    n_contradicted = 0
    n_inconclusive = 0
    n_failures = 0
    n_skipped_large = 0
    files_processed = 0
    t0 = time.time()

    with OUT_JSONL_LOCAL.open("w") as fh:
        for i, entry in enumerate(sample):
            pmc = entry["pmc_id"]
            fname = entry["file_name"]
            path = KOH / pmc / fname
            if not path.exists():
                continue
            try:
                if path.stat().st_size > MAX_BYTES:
                    n_skipped_large += 1
                    continue
            except OSError:
                continue
            try:
                # row_context_boost=False so we get the *unboosted* confidence
                report = detect_file(str(path), row_context_boost=False)
            except Exception:
                n_failures += 1
                continue
            if not report.suspicions:
                files_processed += 1
                continue
            # Lazy-load the sheets so we can do row-context lookups
            from uncorrupt.app import _load_all_sheets
            try:
                sheets = _load_all_sheets(str(path))
            except Exception:
                n_failures += 1
                continue
            for s in report.suspicions:
                if not s.suggestion:
                    continue
                df = sheets.get(s.sheet) if s.sheet else (
                    next(iter(sheets.values())) if len(sheets) == 1 else None
                )
                if df is None or s.row not in df.index:
                    continue
                row = df.loc[s.row]
                resolved: list[tuple[str, str]] = []
                for col, val in row.items():
                    if str(col) == s.column or not isinstance(val, str):
                        continue
                    gene = xref.lookup(val.strip())
                    if gene is not None:
                        resolved.append((val.strip(), gene))
                suggestion_set = {p.strip() for p in s.suggestion.split("|") if p.strip()}
                match_set = _expand_with_hgnc_aliases(hgnc, suggestion_set)
                # Case-insensitive matching : cross-species ortholog
                # symbols are mixed-case (mouse Marchf10, fly mei-9)
                # while HGNC suggestions are uppercase.
                match_set_lower = {s.lower() for s in match_set}
                corroborated = [
                    (eid, g) for eid, g in resolved
                    if g.lower() in match_set_lower
                ]
                if not resolved:
                    n_inconclusive += 1
                    label = "inconclusive"
                elif corroborated:
                    n_corroborated += 1
                    label = "corroborated"
                else:
                    n_contradicted += 1
                    label = "contradicted"
                # Tag each record with the Koh classifier stratum so the
                # downstream analysis can break down precision per pool.
                # pool_tag is file-level (`entry`) but we compute it inside
                # the loop because the write itself is per-suspicion : every
                # flag in the file inherits the same pool.
                if entry.get("has_gene_symbols") and entry.get("n_suspicions", 0) > 0:
                    pool_tag = "A"  # with-corruption
                elif entry.get("has_gene_symbols"):
                    pool_tag = "B"  # clean-with-genes
                else:
                    pool_tag = "C"  # no-genes
                fh.write(json.dumps({
                    "pmc_id": pmc,
                    "file": fname,
                    "pool": pool_tag,
                    "sheet": s.sheet,
                    "column": s.column,
                    "row": s.row,
                    "value": str(s.value),
                    "kind": s.kind,
                    "suggestion": s.suggestion,
                    "confidence": s.confidence,
                    "label": label,
                    "resolved_in_row": resolved[:6],
                    "corroborated_by": corroborated[:6],
                }) + "\n")
            files_processed += 1
            if (i + 1) % 100 == 0:
                elapsed = time.time() - t0
                print(
                    f"[{i+1}/{len(sample)}] {n_corroborated} corroborated, "
                    f"{n_contradicted} contradicted, {n_inconclusive} inconclusive "
                    f"({elapsed:.0f}s)", flush=True
                )

    print()
    print(f"Files processed: {files_processed} ({n_failures} failed, "
          f"{n_skipped_large} skipped as >50MB)")
    print(f"Corroborated: {n_corroborated}")
    print(f"Contradicted: {n_contradicted}")
    print(f"Inconclusive: {n_inconclusive}")

    # Calibration metrics on the corroborated + contradicted subset
    scores: list[float] = []
    outcomes: list[int] = []
    per_bin_data: dict[float, dict[str, int]] = defaultdict(
        lambda: {"corroborated": 0, "contradicted": 0}
    )
    per_pool_per_bin: dict[str, dict[float, dict[str, int]]] = defaultdict(
        lambda: defaultdict(lambda: {"corroborated": 0, "contradicted": 0})
    )
    with OUT_JSONL.open() as fh:
        for line in fh:
            r = json.loads(line)
            if r["label"] == "inconclusive":
                continue
            o = 1 if r["label"] == "corroborated" else 0
            scores.append(r["confidence"])
            outcomes.append(o)
            per_bin_data[r["confidence"]][r["label"]] += 1
            pool = r.get("pool", "?")
            per_pool_per_bin[pool][r["confidence"]][r["label"]] += 1

    # Markdown report
    lines = [
        "# UnCorrupt confidence calibration with FP anchors",
        "",
        f"Source: stratified sample of {files_processed} Koh cache files "
        f"({n_failures} failed to load).",
        "",
        f"- Corroborated (xref → same gene): {n_corroborated}",
        f"- Contradicted (xref → different gene): {n_contradicted}",
        f"- Inconclusive (no resolvable row xref): {n_inconclusive}",
        "",
        "## Per-bin precision on the labeled subset",
        "",
        "| Stated | n_total | n_corroborated | n_contradicted | Precision | Wilson 95% CI |",
        "|--------|---------|----------------|----------------|-----------|---------------|",
    ]
    for conf in sorted(per_bin_data):
        b = per_bin_data[conf]
        total = b["corroborated"] + b["contradicted"]
        if total == 0:
            continue
        prec = b["corroborated"] / total
        lo, hi = _wilson_ci(b["corroborated"], total)
        lines.append(
            f"| {conf:.2f} | {total} | {b['corroborated']} | "
            f"{b['contradicted']} | {prec:.3f} | "
            f"[{lo:.3f}, {hi:.3f}] |"
        )
    bs = _brier(scores, outcomes)
    ece = _ece(scores, outcomes)

    # Per-pool tables (Pool A=with-corruption, B=clean-with-genes, C=no-genes)
    pool_titles = {"A": "with-corruption", "B": "clean-with-genes",
                   "C": "no-genes"}
    for pool in ("A", "B", "C"):
        if not per_pool_per_bin[pool]:
            continue
        lines.append(f"")
        lines.append(f"### Pool {pool} ({pool_titles[pool]})")
        lines.append("")
        lines.append("| Stated | n_labeled | corroborated | contradicted | Precision |")
        lines.append("|--------|-----------|--------------|--------------|-----------|")
        for conf in sorted(per_pool_per_bin[pool]):
            b = per_pool_per_bin[pool][conf]
            t = b["corroborated"] + b["contradicted"]
            if t == 0:
                continue
            p = b["corroborated"] / t
            lines.append(
                f"| {conf:.2f} | {t} | {b['corroborated']} | "
                f"{b['contradicted']} | {p:.3f} |"
            )

    lines += [
        "",
        "## Scalar metrics on the xref-labeled subset",
        "",
        f"- Brier score: {bs:.4f}",
        f"- Expected Calibration Error (ECE): {ece:.4f}",
        f"- Total labeled cells: {len(scores)} "
        f"({sum(outcomes)} positive, {len(outcomes)-sum(outcomes)} negative)",
        "",
        "## Methodology",
        "",
        "Each flagged cell is labeled by independent xref evidence in the",
        "same row. The detector's *unboosted* confidence (`row_context_boost=False`)",
        "is the raw score; the xref decision is the ground-truth label.",
        "Corroborated rows are TP; contradicted rows are FP candidates; rows",
        "without any resolvable external ID are inconclusive and excluded.",
        "",
        "HGNC equivalence: a suggested gene matches when any element of",
        "`prev_symbol → current_symbol` and `alias → current_symbol`",
        "transitively resolves to the corroborating ID's gene.",
        "",
        "Wilson 95% CIs (Wilson 1927) over Beta posterior estimates",
        "(Laplace prior). Brier score & ECE follow Brier 1950, Murphy 1973,",
        "Naeini et al. 2015, Guo et al. 2017.",
        "",
    ]
    OUT_REPORT.write_text("\n".join(lines))
    print(f"\nWrote {OUT_REPORT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
