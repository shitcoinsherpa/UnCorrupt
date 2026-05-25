"""Smoke tests over the public test fixtures.

Each fixture in tests/fixtures/ demonstrates one corruption class. These tests
let external users run `pytest tests/test_fixtures.py` and see the detector
work on real (tiny) xlsx files without needing the 22 GB validation corpus.

The expected-outputs JSON keeps this test data-driven : adding a new fixture
means appending one entry to expected.json, not adding test functions.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from uncorrupt.detector import detect_file

FIXTURE_DIR = Path(__file__).parent / "fixtures"
EXPECTED = json.loads((FIXTURE_DIR / "expected.json").read_text())


@pytest.mark.parametrize("fixture_name", sorted(EXPECTED.keys()))
def test_fixture_meets_expectations(fixture_name: str) -> None:
    """Every fixture honors its expected.json contract."""
    fp = FIXTURE_DIR / fixture_name
    assert fp.exists(), f"fixture file {fp} missing"

    spec = EXPECTED[fixture_name]
    report = detect_file(str(fp))
    n_susp = len(report.suspicions)

    # Count constraints
    if "min_suspicions" in spec:
        assert n_susp >= spec["min_suspicions"], (
            f"{fixture_name}: expected >= {spec['min_suspicions']} suspicions, "
            f"got {n_susp}; suspicions: {[(s.column, s.row, s.value, s.suggestion) for s in report.suspicions]}"
        )
    if "max_suspicions" in spec:
        assert n_susp <= spec["max_suspicions"], (
            f"{fixture_name}: expected <= {spec['max_suspicions']} suspicions, "
            f"got {n_susp}; suspicions: {[(s.column, s.row, s.value, s.suggestion) for s in report.suspicions]}"
        )

    # Suggestion constraints: each expected gene-name family must appear in
    # at least one suggestion. We compare against the candidate set, allowing
    # for the detector to emit either form (SEPT2 vs SEPT2 | SEP2 etc.).
    if "expected_suggestions" in spec:
        all_suggestion_tokens: set[str] = set()
        for s in report.suspicions:
            if s.suggestion:
                for part in s.suggestion.split("|"):
                    all_suggestion_tokens.add(part.strip())
        for expected_family in spec["expected_suggestions"]:
            assert any(g in all_suggestion_tokens for g in expected_family), (
                f"{fixture_name}: expected one of {expected_family} in suggestions; "
                f"got {sorted(all_suggestion_tokens)}"
            )

    # Kind constraint (e.g. fixture_06 expects id-float)
    if "expected_kind" in spec:
        kinds = {s.kind for s in report.suspicions}
        assert spec["expected_kind"] in kinds, (
            f"{fixture_name}: expected kind {spec['expected_kind']!r} not found; "
            f"got {sorted(kinds)}"
        )

    # Xref-corroboration constraint (Proposal B): for fixtures with paired
    # RefSeq/Ensembl IDs, confirm the row-context boost fired AND the
    # corroborating ID appears in the reason field for audit.
    if "min_xref_corroborated" in spec:
        n_corr = sum(1 for s in report.suspicions
                      if "CORROBORATED" in (s.reason or ""))
        assert n_corr >= spec["min_xref_corroborated"], (
            f"{fixture_name}: expected >= {spec['min_xref_corroborated']} "
            f"row-xref-corroborated suspicions, got {n_corr}"
        )
        # Each xref-boosted suspicion has confidence >= 0.85 (auto-accept band)
        for s in report.suspicions:
            if "CORROBORATED" in (s.reason or ""):
                assert s.confidence >= 0.85, (
                    f"{fixture_name}: xref-corroborated cell had "
                    f"confidence {s.confidence} (expected >= 0.85)"
                )
    if "expected_xref_corroboration_keywords" in spec:
        all_reasons = " ".join(s.reason or "" for s in report.suspicions)
        for kw in spec["expected_xref_corroboration_keywords"]:
            assert kw in all_reasons, (
                f"{fixture_name}: expected {kw!r} to appear in some "
                f"suspicion reason as xref corroboration"
            )


def test_fixture_directory_has_readme() -> None:
    """Self-documenting fixture directory : meta-test."""
    # Ensure expected.json describes every actual fixture (no orphans)
    actual_fixtures = {
        p.name for p in FIXTURE_DIR.iterdir()
        if p.suffix.lower() in (".xlsx", ".xls")
    }
    documented_fixtures = set(EXPECTED.keys())
    orphan = actual_fixtures - documented_fixtures
    missing = documented_fixtures - actual_fixtures
    assert not orphan, f"fixtures without expected.json entries: {orphan}"
    assert not missing, f"expected.json entries without fixture files: {missing}"
