"""Wiring smoke tests for the Gradio "Classify with local LLM" button.

These verify the upload → detect → packet → verdict-dataframe pipeline
without forcing CI to download/load the 2 GB Qwen model. The slow
end-to-end variant runs only when UNCORRUPT_RUN_LLM_SMOKE=1.
"""
from __future__ import annotations

import os
from datetime import date
from pathlib import Path

import pandas as pd
import pytest

from uncorrupt import app as app_module
from uncorrupt.local_classifier import ClassificationResult


def _write_xlsx(path: Path, data: dict[str, list]) -> Path:
    df = pd.DataFrame({k: pd.array(v, dtype=object) for k, v in data.items()})
    df.to_excel(path, index=False)
    return path


def _file_obj(p: Path):
    """Mimic the gradio NamedString shape : only the .name attr is read."""
    return type("F", (), {"name": str(p)})()


def _gene_date_xlsx(path: Path) -> Path:
    return _write_xlsx(path, {
        "symbol": ["BRCA1", "TP53", date(2024, 9, 2), "MARCHF1",
                   "SEPTIN2", date(2024, 3, 1), "PTEN", "EGFR", "KRAS"],
        "fold_change": [1.2, 3.4, 0.8, 2.1, 1.7, 4.0, 0.9, 1.5, 2.0],
    })


def test_classify_with_local_llm_returns_empty_for_no_file() -> None:
    df = app_module.classify_with_local_llm(None)
    assert df.empty


def test_classify_with_local_llm_handles_clean_file_with_no_suspicions(
    monkeypatch, tmp_path: Path,
) -> None:
    """If the detector finds nothing, the function returns the 'no suspicions' marker
    and never invokes the LLM."""
    calls = {"n": 0}

    def boom(_md: str, **_kw) -> ClassificationResult:
        calls["n"] += 1
        raise AssertionError("LLM must not be called when there are no suspicions")

    monkeypatch.setattr("uncorrupt.local_classifier.classify_packet", boom)

    f = _write_xlsx(tmp_path / "clean.xlsx", {
        "symbol": ["BRCA1", "TP53", "EGFR", "KRAS", "PTEN"],
        "fold_change": [1.2, 3.4, 0.8, 2.1, 1.7],
    })
    df = app_module.classify_with_local_llm(_file_obj(f))
    assert calls["n"] == 0
    assert "info" in df.columns
    assert df.iloc[0]["info"] == "No suspicions to classify."


def test_classify_with_local_llm_wires_verdicts_into_dataframe(
    monkeypatch, tmp_path: Path,
) -> None:
    """End-to-end wiring: upload → suspicions → packets → stub-LLM → dataframe.
    Stubs the model so the test runs in <1 s."""
    invocations: list[str] = []

    def fake_classify(md: str, **_kw) -> ClassificationResult:
        invocations.append(md)
        return ClassificationResult(
            verdict="TP", reason="stubbed verdict for wiring smoke",
            raw_response="VERDICT: TP\nREASON: stubbed verdict for wiring smoke",
            inference_seconds=0.0,
        )

    monkeypatch.setattr("uncorrupt.local_classifier.classify_packet", fake_classify)

    f = _gene_date_xlsx(tmp_path / "gene-date.xlsx")
    df = app_module.classify_with_local_llm(_file_obj(f))

    # the file holds two date-corruption cells; both should reach the LLM
    assert len(invocations) >= 1, "expected at least one packet to reach the LLM"
    # verdict column is populated with the stub
    assert "verdict" in df.columns
    assert (df["verdict"] == "TP").all()
    # the per-row context columns are present
    for col in ("column", "row", "value", "suggestion", "reason"):
        assert col in df.columns


def test_classify_with_local_llm_records_errors_per_row(
    monkeypatch, tmp_path: Path,
) -> None:
    """If the LLM raises, the row gets an ERROR verdict : no crash, no skip."""
    def fake_classify(_md: str, **_kw) -> ClassificationResult:
        raise RuntimeError("simulated model load failure")

    monkeypatch.setattr("uncorrupt.local_classifier.classify_packet", fake_classify)

    f = _gene_date_xlsx(tmp_path / "gene-date.xlsx")
    df = app_module.classify_with_local_llm(_file_obj(f))

    assert len(df) >= 1
    assert (df["verdict"] == "ERROR").all()
    assert df.iloc[0]["reason"].startswith("simulated model load failure")


def test_classify_with_local_llm_caps_at_50_invocations(
    monkeypatch, tmp_path: Path,
) -> None:
    """The UI path caps per-invocation work at 50 cells to keep latency bounded.
    Build a file guaranteed to flag >50 corruptions and verify the cap holds."""
    invocations = {"n": 0}

    def fake_classify(_md: str, **_kw) -> ClassificationResult:
        invocations["n"] += 1
        return ClassificationResult(
            verdict="TP", reason="stub",
            raw_response="VERDICT: TP\nREASON: stub",
            inference_seconds=0.0,
        )

    monkeypatch.setattr("uncorrupt.local_classifier.classify_packet", fake_classify)

    # 60 distinct gene-date corruptions : each is a real Excel auto-convert pattern.
    # Use distinct years so the values dedupe to separate cells.
    bulk_dates = [date(2014 + (i % 10), 9, 2 + (i % 25)) for i in range(60)]
    symbols = ["BRCA1"] * 5 + bulk_dates + ["TP53"] * 5
    f = _write_xlsx(tmp_path / "bulk.xlsx", {
        "symbol": symbols,
        "x": list(range(len(symbols))),
    })
    df = app_module.classify_with_local_llm(_file_obj(f))

    assert invocations["n"] <= 50, (
        f"expected ≤50 LLM calls per UI invocation, got {invocations['n']}"
    )
    assert len(df) <= 50


@pytest.mark.skipif(
    os.environ.get("UNCORRUPT_RUN_LLM_SMOKE") != "1",
    reason="real-model smoke; set UNCORRUPT_RUN_LLM_SMOKE=1 to run "
           "(downloads ~2 GB on first call, ~30 s/cell on CPU)",
)
def test_classify_with_local_llm_real_model_end_to_end(tmp_path: Path) -> None:
    """Gated end-to-end with the real Qwen model. Runs only when explicitly opted in."""
    f = _gene_date_xlsx(tmp_path / "gene-date.xlsx")
    df = app_module.classify_with_local_llm(_file_obj(f))

    assert len(df) >= 1
    assert "verdict" in df.columns
    assert df["verdict"].iloc[0] in {"TP", "FP", "INCONCLUSIVE", "ERROR"}
