"""Synthetic smoke test for the corruption detector : NOT scientific evidence.

This module runs the detector against synthetically-constructed columns whose
corruption was applied by our own simulator. It exists to catch regressions in
CI: if the detector stops finding obvious dates-in-an-identifier-column or
floats-in-an-identifier-column, this fails fast.

It is **not** the canonical performance measurement. The whole reason the
gene-name corruption problem persisted for ten years is that synthetic
benchmarks said everything was fine while real supplementary files were
30% corrupted. Treat the output of this module as a regression guard, not as
evidence. Canonical evidence lives in `src/uncorrupt/corpus.py`
(Ziemann 2021 PMC corpus, real registry bulk downloads).
"""
from __future__ import annotations

import json
import random
import time
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

from .corpus import derive_date_prone_map, load_hgnc
from .detector import detect
from .simulator import excel_autoconvert


@dataclass
class RegistryScore:
    registry: str
    n_total: int
    n_corrupted: int
    true_positives: int
    false_positives: int
    false_negatives: int
    precision: float
    recall: float
    f1: float


@dataclass
class BenchmarkRun:
    timestamp: str
    seed: int
    column_size: int
    corruption_fraction: float
    runtime_seconds: float
    scores: list[RegistryScore] = field(default_factory=list)


def _make_riken_id(rng: random.Random) -> str:
    digits = [str(rng.randint(0, 9)) for _ in range(7)]
    letter = chr(rng.randint(ord("A"), ord("Z")))
    suffix = [str(rng.randint(0, 9)) for _ in range(2)]
    return "".join(digits + [letter] + suffix)


def _score_registry(
    registry: str,
    clean_pool: list[str],
    column_size: int,
    corruption_fraction: float,
    rng: random.Random,
) -> RegistryScore:
    chosen = rng.choices(clean_pool, k=column_size)
    target_corrupt = rng.sample(range(column_size), int(column_size * corruption_fraction))
    corrupt_set = set(target_corrupt)

    column: list[object] = []
    truth: list[bool] = []
    for i, gid in enumerate(chosen):
        if i in corrupt_set:
            corrupted = excel_autoconvert(gid)
            if corrupted != gid:
                column.append(corrupted)
                truth.append(True)
                continue
        column.append(gid)
        truth.append(False)

    df = pd.DataFrame({
        "symbol": pd.array(column, dtype=object),
        "value": [rng.random() for _ in range(column_size)],
    })

    report = detect(df)
    flagged = {s.row for s in report.suspicions if s.column == "symbol"}

    tp = sum(1 for i, t in enumerate(truth) if t and i in flagged)
    fp = sum(1 for i, t in enumerate(truth) if not t and i in flagged)
    fn = sum(1 for i, t in enumerate(truth) if t and i not in flagged)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 1.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 1.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    return RegistryScore(
        registry=registry,
        n_total=column_size,
        n_corrupted=sum(truth),
        true_positives=tp,
        false_positives=fp,
        false_negatives=fn,
        precision=precision,
        recall=recall,
        f1=f1,
    )


def run_benchmark(
    seed: int = 42,
    column_size: int = 200,
    corruption_fraction: float = 0.2,
) -> BenchmarkRun:
    rng = random.Random(seed)
    t0 = time.monotonic()

    riken_pool = [_make_riken_id(rng) for _ in range(2000)]
    hgnc_date_prone = list(derive_date_prone_map(load_hgnc()).keys())

    scores = [
        _score_registry("HGNC_DATE_PRONE", hgnc_date_prone, column_size,
                        corruption_fraction, rng),
        _score_registry("RIKEN", riken_pool, column_size, corruption_fraction, rng),
    ]

    return BenchmarkRun(
        timestamp=datetime.now(UTC).isoformat(),
        seed=seed,
        column_size=column_size,
        corruption_fraction=corruption_fraction,
        runtime_seconds=time.monotonic() - t0,
        scores=scores,
    )


def append_to_ledger(run: BenchmarkRun, ledger_path: Path) -> None:
    """Append one benchmark run to a JSONL ledger."""
    ledger_path.parent.mkdir(parents=True, exist_ok=True)
    with ledger_path.open("a") as f:
        f.write(json.dumps(asdict(run)) + "\n")


def _format_report(run: BenchmarkRun) -> str:
    lines = [
        f"Benchmark run @ {run.timestamp} (seed={run.seed}, n={run.column_size}, "
        f"corrupt={run.corruption_fraction:.0%}, runtime={run.runtime_seconds:.3f}s)",
        "-" * 80,
        f"{'registry':<20} {'n':>5} {'corrupt':>8} {'TP':>4} {'FP':>4} {'FN':>4} "
        f"{'P':>6} {'R':>6} {'F1':>6}",
    ]
    for s in run.scores:
        lines.append(
            f"{s.registry:<20} {s.n_total:>5} {s.n_corrupted:>8} "
            f"{s.true_positives:>4} {s.false_positives:>4} {s.false_negatives:>4} "
            f"{s.precision:>6.3f} {s.recall:>6.3f} {s.f1:>6.3f}"
        )
    return "\n".join(lines)


if __name__ == "__main__":
    run = run_benchmark()
    print(_format_report(run))
    ledger = Path("results") / "benchmark_ledger.jsonl"
    append_to_ledger(run, ledger)
    print(f"\nAppended to {ledger}")
