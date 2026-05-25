"""Aggregate local Qwen classification results into a summary report.

Reads results/local_classification.jsonl produced by classify_packets_locally.py
and the corresponding packet_*.md files (for detector kind, column type,
PMC id, journal, year). Emits:

  - tally by verdict
  - tally by detector kind (gene-date / gene-date-serial / id-float)
  - tally by column type (gene_symbol / measurement / free_text / ...)
  - per-PMC verdict summary
  - inference timing stats

Usage:
    python scripts/aggregate_local_classifications.py
"""
from __future__ import annotations

import json
import re
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
LEDGER = PROJECT_ROOT / "results/local_classification.jsonl"
PACKET_DIR = PROJECT_ROOT / "data/derived/subagent_packets"


def parse_packet(md: str) -> dict:
    """Extract the structured fields we need from a packet_*.md file."""
    out: dict = {}
    for line in md.splitlines():
        if line.startswith("# Cell Review Packet"):
            m = re.search(r"PMC\d+", line)
            if m:
                out["pmc_id"] = m.group(0)
        elif line.startswith("- Journal:"):
            out["journal"] = line.split(":", 1)[1].strip()
        elif line.startswith("- Year:"):
            out["year"] = line.split(":", 1)[1].strip()
        elif line.startswith("- Species:"):
            out["species"] = line.split(":", 1)[1].strip()
        elif line.startswith("- Detector kind:"):
            out["kind"] = line.split(":", 1)[1].strip().strip("`")
        elif line.startswith("- Detector suggestion:"):
            out["suggestion"] = line.split(":", 1)[1].strip().strip("`")
        elif line.startswith("- Detector confidence:"):
            try:
                out["confidence"] = float(line.split(":", 1)[1].strip())
            except ValueError:
                out["confidence"] = None
        elif line.startswith("- Classified as:"):
            out["column_type"] = (
                line.split(":", 1)[1].strip().lstrip("*").rstrip("*").strip()
            )
    return out


def main() -> None:
    if not LEDGER.exists():
        print(f"missing ledger: {LEDGER}", file=sys.stderr)
        sys.exit(1)

    rows: list[dict] = []
    with LEDGER.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rows.append(json.loads(line))

    if not rows:
        print("no rows in ledger", file=sys.stderr)
        sys.exit(1)

    # Enrich each row with packet metadata
    enriched: list[dict] = []
    for r in rows:
        packet_path = PACKET_DIR / r["packet"]
        meta = parse_packet(packet_path.read_text()) if packet_path.exists() else {}
        enriched.append({**r, **meta})

    n = len(enriched)
    print(f"=== local Qwen classification summary ({n} packets) ===\n")

    # Verdict tally
    verdicts = Counter(r["verdict"] for r in enriched)
    print("Verdict distribution:")
    for v in ("TP", "FP", "INCONCLUSIVE"):
        c = verdicts.get(v, 0)
        print(f"  {v:<14} {c:>4}  ({c / n * 100:5.1f}%)")
    print()

    # By detector kind
    by_kind: dict[str, Counter] = defaultdict(Counter)
    for r in enriched:
        by_kind[r.get("kind", "?")][r["verdict"]] += 1
    print("Verdict by detector kind:")
    print(f"  {'kind':<22} {'TP':>4} {'FP':>4} {'INC':>4} {'total':>6}")
    for kind, c in sorted(by_kind.items()):
        total = sum(c.values())
        print(
            f"  {kind:<22} {c.get('TP', 0):>4} {c.get('FP', 0):>4} "
            f"{c.get('INCONCLUSIVE', 0):>4} {total:>6}"
        )
    print()

    # By column type
    by_col: dict[str, Counter] = defaultdict(Counter)
    for r in enriched:
        by_col[r.get("column_type", "?")][r["verdict"]] += 1
    print("Verdict by column-classifier type:")
    print(f"  {'column_type':<16} {'TP':>4} {'FP':>4} {'INC':>4} {'total':>6}")
    for col, c in sorted(by_col.items()):
        total = sum(c.values())
        print(
            f"  {col:<16} {c.get('TP', 0):>4} {c.get('FP', 0):>4} "
            f"{c.get('INCONCLUSIVE', 0):>4} {total:>6}"
        )
    print()

    # Inference timing
    times = [r["inference_seconds"] for r in enriched if "inference_seconds" in r]
    if times:
        print(
            f"Inference timing (n={len(times)}): "
            f"mean={statistics.mean(times):.1f}s  "
            f"median={statistics.median(times):.1f}s  "
            f"p95={sorted(times)[int(0.95 * len(times))]:.1f}s  "
            f"max={max(times):.1f}s  "
            f"total={sum(times) / 60:.1f}min"
        )
        print()

    # Per-PMC summary (1 line per PMC, across all flagged cells)
    by_pmc: dict[str, Counter] = defaultdict(Counter)
    for r in enriched:
        by_pmc[r.get("pmc_id", "?")][r["verdict"]] += 1
    n_pmcs = len(by_pmc)
    pmcs_any_tp = sum(1 for c in by_pmc.values() if c.get("TP", 0) > 0)
    print(f"Per-PMC corruption signal: {pmcs_any_tp}/{n_pmcs} PMCs have ≥1 TP cell "
          f"({pmcs_any_tp / n_pmcs * 100:.1f}%)")


if __name__ == "__main__":
    main()
