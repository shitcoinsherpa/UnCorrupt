"""Regression tests for the synthetic benchmark.

These pin minimum precision/recall thresholds per registry. If the detector
regresses, these tests fail before any user sees the bad grade.
"""
from __future__ import annotations

from uncorrupt.benchmark import run_benchmark


def test_benchmark_runs_deterministically() -> None:
    a = run_benchmark(seed=42, column_size=200)
    b = run_benchmark(seed=42, column_size=200)
    assert [s.f1 for s in a.scores] == [s.f1 for s in b.scores]


def test_hgnc_date_prone_recall_above_threshold() -> None:
    run = run_benchmark(seed=42, column_size=400)
    hgnc = next(s for s in run.scores if s.registry == "HGNC_DATE_PRONE")
    assert hgnc.recall >= 0.95, f"HGNC recall regressed: {hgnc.recall:.3f}"
    assert hgnc.precision >= 0.95, f"HGNC precision regressed: {hgnc.precision:.3f}"


def test_riken_recall_above_threshold() -> None:
    run = run_benchmark(seed=42, column_size=400)
    riken = next(s for s in run.scores if s.registry == "RIKEN")
    assert riken.recall >= 0.95, f"RIKEN recall regressed: {riken.recall:.3f}"
    assert riken.precision >= 0.90, f"RIKEN precision regressed: {riken.precision:.3f}"


def test_benchmark_produces_some_corrupted_samples() -> None:
    run = run_benchmark(seed=42, column_size=400, corruption_fraction=0.2)
    for s in run.scores:
        assert s.n_corrupted > 0, f"{s.registry} produced 0 corrupted samples"
