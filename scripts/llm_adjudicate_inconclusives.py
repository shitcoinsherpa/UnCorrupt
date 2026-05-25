"""LLM-assisted adjudication of inconclusive cells from the calibration walk.

The xref-anchored calibration walk labels cells as `corroborated` (TP),
`contradicted` (FP), or `inconclusive` (no resolvable external ID in the
same row). The inconclusive set is the largest population (~27 K cells in
.

This script samples inconclusive cells stratified by confidence band,
re-builds the review packet for each (column header + neighbouring cells
+ detector's reasoning), and feeds it to the local Qwen2.5-3B-Instruct-GGUF
model used by the Gradio app. Each LLM verdict (TP / FP / INCONCLUSIVE)
becomes a weakly-supervised label.

The resulting labels expand the validatable set so we can compute a
post-boost precision number that includes the inconclusive-by-xref subset.

Run:
   python scripts/llm_adjudicate_inconclusives.py --sample 200 --seed 20260521
"""
from __future__ import annotations

import argparse
import json
import random
import sys
import warnings
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from uncorrupt.local_classifier import classify_packet  # noqa: E402
from uncorrupt.subagent_packet import (  # noqa: E402
    build_packet_for_suspicion,
    format_packet_as_markdown,
)
from uncorrupt.detector import Suspicion  # noqa: E402

LEDGER = ROOT / "results/calibration_fp_anchors.jsonl"
KOH = ROOT / "data/raw/koh_replication/files"
OUT_JSONL = ROOT / "results/llm_adjudication.jsonl"
OUT_REPORT = ROOT / "results/llm_adjudication.md"


def _stratify_sample(records: list[dict], n: int, seed: int) -> list[dict]:
    """Stratified sample: keep at least 30 cells per confidence band in
    [0.30, 0.50, 0.55, 0.60] where the post-boost layer leaves
    inconclusive cells in the user-visible range."""
    rng = random.Random(seed)
    by_conf: dict[float, list[dict]] = {}
    for r in records:
        c = round(r["confidence"], 2)
        by_conf.setdefault(c, []).append(r)
    sample: list[dict] = []
    per_bucket = max(30, n // len(by_conf))
    for c, rs in by_conf.items():
        rng.shuffle(rs)
        sample.extend(rs[:per_bucket])
    rng.shuffle(sample)
    return sample[:n]


def _suspicion_from_record(r: dict) -> Suspicion:
    return Suspicion(
        column=r["column"], row=r["row"], value=r["value"],
        kind=r["kind"], suggestion=r.get("suggestion"),
        confidence=r["confidence"], reason=r.get("reason"),
        sheet=r.get("sheet"),
    )


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample", type=int, default=200)
    p.add_argument("--seed", type=int, default=20260521)
    p.add_argument("--only-kinds", type=str, default="gene-date,gene-date-string,gene-date-serial",
                   help="Comma-separated list of kinds to adjudicate.")
    p.add_argument("--sample-offset", type=int, default=0,
                   help="Skip first N stratified-sampled cells (for parallel sharding).")
    p.add_argument("--sample-limit", type=int, default=None,
                   help="Process at most N cells starting from --sample-offset.")
    p.add_argument("--out", type=str, default=str(OUT_JSONL),
                   help="Output JSONL path (override for sharded runs).")
    args = p.parse_args()

    if not LEDGER.exists():
        print(f"FATAL: {LEDGER} missing", file=sys.stderr)
        return 1

    keep_kinds = set(args.only_kinds.split(","))
    inconclusive = [
        json.loads(line)
        for line in LEDGER.open()
        if (r := json.loads(line))["label"] == "inconclusive"
        and r["kind"] in keep_kinds
    ]
    # Re-read because the walrus shadowed the iterator
    inconclusive = []
    for line in LEDGER.open():
        r = json.loads(line)
        if r["label"] == "inconclusive" and r["kind"] in keep_kinds:
            inconclusive.append(r)

    print(f"[llm-adj] {len(inconclusive)} inconclusive cells in scope")
    sample = _stratify_sample(inconclusive, args.sample, args.seed)
    print(f"[llm-adj] sampled {len(sample)} (stratified by confidence)")

    # Slice for parallel sharding : _stratify_sample is deterministic on
    # (records, n, seed), so all workers compute the same sample list.
    if args.sample_offset or args.sample_limit is not None:
        end = (args.sample_offset + args.sample_limit) if args.sample_limit is not None else len(sample)
        sample = sample[args.sample_offset:end]
        print(f"[llm-adj] sliced to [{args.sample_offset}:{end}] -> {len(sample)} cells",
              flush=True)

    out_path = Path(args.out)
    results: list[dict] = []
    verdicts: Counter[str] = Counter()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        for i, r in enumerate(sample, 1):
            pmc = r["pmc_id"]; fname = r["file"]
            path = KOH / pmc / fname
            if not path.exists():
                continue
            try:
                susp = _suspicion_from_record(r)
                susp_dict = {
                    "column": susp.column, "row": susp.row, "value": susp.value,
                    "kind": susp.kind, "suggestion": susp.suggestion,
                    "confidence": susp.confidence, "reason": susp.reason,
                    "sheet": susp.sheet,
                    "pmc_id": r["pmc_id"],
                    "file_name": r["file"],
                }
                packet = build_packet_for_suspicion(susp_dict, path, None)
                if packet is None:
                    continue
                packet_md = format_packet_as_markdown(packet)
                cls = classify_packet(packet_md)
            except Exception as exc:
                print(f"[llm-adj] {i}: {type(exc).__name__} {exc!s:.80}",
                      flush=True)
                continue
            verdicts[cls.verdict] += 1
            rec = {
                **r,
                "llm_verdict": cls.verdict,
                "llm_reason": cls.reason[:240],
                "llm_inference_seconds": cls.inference_seconds,
            }
            fh.write(json.dumps(rec) + "\n")
            results.append(rec)
            if i % 25 == 0:
                print(f"[llm-adj] {i}/{len(sample)}  TP={verdicts['TP']} "
                      f"FP={verdicts['FP']} INC={verdicts['INCONCLUSIVE']}",
                      flush=True)

    n = sum(verdicts.values())
    if n == 0:
        print("[llm-adj] no records produced", file=sys.stderr)
        return 2
    n_tp = verdicts["TP"]; n_fp = verdicts["FP"]; n_inc = verdicts["INCONCLUSIVE"]
    n_decided = n_tp + n_fp
    precision = n_tp / n_decided if n_decided else 0.0
    import math
    if n_decided:
        z = 1.96; p = precision
        d = 1 + z*z/n_decided
        c = (p + z*z/(2*n_decided)) / d
        h = z * math.sqrt(p*(1-p)/n_decided + z*z/(4*n_decided*n_decided)) / d
        lo, hi = max(0, c-h), min(1, c+h)
    else:
        lo, hi = 0.0, 1.0

    OUT_REPORT.write_text("\n".join([
        "# LLM-assisted adjudication : inconclusive subset",
        "",
        f"Timestamp: {datetime.now(UTC).isoformat()}",
        f"Model: Qwen2.5-3B-Instruct-GGUF (local; temperature=0.0)",
        f"Ledger: `results/calibration_fp_anchors.jsonl`",
        f"Sample size: {n} (stratified by confidence band)",
        "",
        "## Outcomes",
        "",
        "| Verdict | n | Fraction |",
        "|---|---:|---:|",
        f"| TP | {n_tp} | {n_tp/n:.3f} |",
        f"| FP | {n_fp} | {n_fp/n:.3f} |",
        f"| INCONCLUSIVE | {n_inc} | {n_inc/n:.3f} |",
        "",
        "## Inconclusive-subset precision (LLM-anchored)",
        "",
        f"Among the {n_decided} cells the LLM could decide (TP+FP), precision is "
        f"**{precision:.4f}** Wilson 95 % [{lo:.4f}, {hi:.4f}].",
        "",
        "This expands the validatable set beyond the xref-anchored ground truth.",
    ]))
    print()
    print(f"LLM-adjudicated {n} cells. TP={n_tp} FP={n_fp} INC={n_inc}")
    print(f"Inconclusive-subset precision: {precision:.4f} "
          f"Wilson 95% [{lo:.4f}, {hi:.4f}]")
    print(f"Wrote {OUT_REPORT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
