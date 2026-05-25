"""Unit tests for the local LLM classifier : model-load-free.

The actual model download + inference is exercised by a smoke script; these
tests check the parsing/dispatch logic so CI doesn't depend on a 2GB
GGUF download.
"""
from __future__ import annotations

from uncorrupt.local_classifier import classify_packet  # noqa: F401
from uncorrupt.local_classifier import ClassificationResult, SYSTEM_PROMPT
import re


def test_system_prompt_contains_required_definitions() -> None:
    """Sanity: the system prompt explicitly defines TP/FP/INCONCLUSIVE."""
    assert "TP:" in SYSTEM_PROMPT
    assert "FP:" in SYSTEM_PROMPT
    assert "INCONCLUSIVE:" in SYSTEM_PROMPT
    assert "Excel" in SYSTEM_PROMPT


def test_system_prompt_includes_excel_serial_correction() -> None:
    """The system prompt must teach the model the correct Excel-serial math
    (sub-agents got this wrong on initial test)."""
    assert "41897" in SYSTEM_PROMPT
    assert "2014-09-15" in SYSTEM_PROMPT


def test_response_parser_handles_strict_format() -> None:
    from uncorrupt.local_classifier import _parse_response
    cases = [
        ("VERDICT: TP\nREASON: gene-symbol column with date.", "TP", "gene-symbol column with date."),
        ("verdict: fp\nreason: legitimate measurement.", "FP", "legitimate measurement."),
        ("VERDICT:INCONCLUSIVE\nREASON:insufficient context.", "INCONCLUSIVE", "insufficient context."),
        ("Some preamble.\n\nVERDICT: TP\nREASON: clear corruption.", "TP", "clear corruption."),
    ]
    for raw, exp_v, exp_r in cases:
        v, r = _parse_response(raw)
        assert v == exp_v, f"verdict mismatch on {raw!r}"
        assert r == exp_r, f"reason mismatch on {raw!r}"


def test_response_parser_handles_loose_format() -> None:
    """Models often skip the 'VERDICT:' prefix. Parser must handle 'TP: ...'
    and similar forms."""
    from uncorrupt.local_classifier import _parse_response
    cases = [
        ("TP: gene-symbol column with date corruption.", "TP"),
        ("FP - legitimate measurement value.", "FP"),
        ("INCONCLUSIVE\nThe context is unclear.", "INCONCLUSIVE"),
        ("**TP**: clear corruption pattern.", "TP"),
        ("Based on the evidence: TP. The cell is corrupted.", "TP"),
    ]
    for raw, exp_v in cases:
        v, _ = _parse_response(raw)
        assert v == exp_v, f"verdict mismatch on {raw!r} → got {v}"


def test_classification_result_dataclass() -> None:
    r = ClassificationResult(
        verdict="TP", reason="strong evidence",
        raw_response="VERDICT: TP\nREASON: strong evidence",
        inference_seconds=1.42,
    )
    assert r.verdict == "TP"
    assert r.reason == "strong evidence"
    assert r.inference_seconds == 1.42
