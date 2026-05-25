"""Batch-classify all generated review packets using the local Qwen2.5-3B GGUF.

Reads packet_*.md files from data/derived/subagent_packets/, runs each
through the local LLM, writes results to results/local_classification.jsonl.

Usage:
    python scripts/classify_packets_locally.py [--n 10]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from uncorrupt.local_classifier import classify_packet  # noqa: E402

PACKET_DIR = PROJECT_ROOT / "data/derived/subagent_packets"
LEDGER = PROJECT_ROOT / "results/local_classification.jsonl"


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=None,
                   help="number of packets to classify (default: all)")
    args = p.parse_args()

    packets = sorted(PACKET_DIR.glob("packet_*.md"))
    if args.n:
        packets = packets[:args.n]
    print(f"Classifying {len(packets)} packets with local Qwen2.5-3B...",
          file=sys.stderr, flush=True)

    LEDGER.parent.mkdir(parents=True, exist_ok=True)

    t0 = time.monotonic()
    counts = {"TP": 0, "FP": 0, "INCONCLUSIVE": 0}
    with LEDGER.open("a") as ledger_f:
        for i, pkt_path in enumerate(packets, 1):
            md = pkt_path.read_text()
            try:
                result = classify_packet(md)
            except Exception as exc:
                print(f"  [{i}/{len(packets)}] {pkt_path.name}: ERROR {type(exc).__name__}: {exc}",
                      file=sys.stderr, flush=True)
                continue
            counts[result.verdict] = counts.get(result.verdict, 0) + 1
            entry = {
                "timestamp": datetime.now(UTC).isoformat(),
                "packet": pkt_path.name,
                "verdict": result.verdict,
                "reason": result.reason,
                "inference_seconds": result.inference_seconds,
            }
            ledger_f.write(json.dumps(entry) + "\n")
            ledger_f.flush()
            if i <= 5 or i % 10 == 0:
                print(f"  [{i}/{len(packets)}] {pkt_path.name}: {result.verdict}  "
                      f"({result.inference_seconds:.1f}s) : {result.reason[:60]}",
                      file=sys.stderr, flush=True)

    elapsed = time.monotonic() - t0
    total = sum(counts.values())
    print(f"\n=== DONE ({elapsed/60:.1f} min, {len(packets)} packets) ===",
          file=sys.stderr)
    for k, v in counts.items():
        print(f"  {k:<14} {v:>4}  ({v/max(1,total)*100:.1f}%)", file=sys.stderr)


if __name__ == "__main__":
    main()
