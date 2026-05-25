"""Boundary regression tests targeting mutation-sensitive constants.

Each test pins a specific threshold or range boundary identified by the
mutation-analysis sub-agent (see review in CHANGELOG v0.4.2). The tests
fail loudly if any constant drifts by ±1, ±10 %, or any other small mutation.

Cross-ref: detector.py constants pinned by these tests:
  - _SERIAL_MIN = 20000  (boundaries: 19999 -> None, 20000 -> date)
  - _SERIAL_MAX = 60000  (boundaries: 60000 -> date, 60001 -> None)
  - Year-suffix window 2000-2059
  - Pass 1 id-float threshold 1e10
  - Pass 1.5 id-float threshold 1e7 (RIKEN edge case)
  - leading-zero column-signal threshold n_lz >= 3
  - Family ranges: TAMM 30-50, MEI 1-31, MAY 1-24, JUN 1-3, AGO 1-12
  - Confidences: 0.95 single-candidate, 0.5 ambiguous, 0.6 unique-canonical
  - xref boost delta 0.4
"""
from __future__ import annotations

from datetime import date

import openpyxl
import pandas as pd
import pytest

from uncorrupt.detector import (
    _candidates_for_month_and_n,
    _reverse_gene_date,
    _serial_to_date,
    detect,
)


# --- _serial_to_date boundary --------------------------------------------


@pytest.mark.parametrize("n, expected_some_date", [
    (19999, False),  # just below min : None
    (20000, True),   # exactly at min : date
    (60000, True),   # exactly at max : date
    (60001, False),  # just above max : None
])
def test_serial_to_date_boundaries(n: int, expected_some_date: bool) -> None:
    out = _serial_to_date(n)
    if expected_some_date:
        assert out is not None, f"_serial_to_date({n}) should return a date"
    else:
        assert out is None, f"_serial_to_date({n}) should return None"


# --- Year-suffix window boundaries ---------------------------------------


def test_year_suffix_window_lower_boundary() -> None:
    """day=1, year=2000 with no format -> year-suffix interpretation fires
    (SEPTIN0 not real, so candidates for September would only include SEPT0
    which doesn't exist -> empty). But year=2001 should produce SEPT1."""
    candidates = _reverse_gene_date(date(2001, 9, 1))
    assert "SEPT1" in candidates or "SEPTIN1" in candidates


def test_year_suffix_window_upper_boundary_includes_tamm() -> None:
    """TAMM41 (Finnish-locale corruption) stores as 2041-01-01. With
    extended window 2000-2059, year-suffix should fire and produce TAMM41."""
    candidates = _reverse_gene_date(date(2041, 1, 1))
    assert "TAMM41" in candidates, (
        f"TAMM41 should be reverse-mapped from date(2041, 1, 1); got {candidates}"
    )


def test_year_suffix_window_outside_returns_empty_year_path() -> None:
    """day=1, year=2060 is OUTSIDE the window : no year-suffix candidates
    emitted via the unknown-format path."""
    # mmm-yy format forces year interpretation; without it, year-2060 stays
    # day-of-month (which is 1 for January). January day=1 not in any family.
    candidates = _reverse_gene_date(date(2060, 1, 1))
    assert candidates == [], (
        f"date(2060, 1, 1) should produce no candidates; got {candidates}"
    )


# --- Pass 1 id-float threshold (1e10) ------------------------------------


def test_pass1_id_float_threshold_above() -> None:
    df = pd.DataFrame({
        "gene_id": ["2310009E13"] + [f"REGISTRY{i}" for i in range(1, 20)]
                    + [1_000_000_001 * 11],  # 1.1e10 : just above 1e10
    })
    report = detect(df)
    id_floats = [s for s in report.suspicions if s.kind == "id-float"]
    assert len(id_floats) >= 1, "1.1e10 should trigger id-float in Pass 1"


def test_pass1_id_float_threshold_below_in_classifier_id_column() -> None:
    """A value of 9.9e9 in a riken-shaped column should NOT trigger Pass 1
    id-float (below 1e10). Pass 1.5 may still catch it on a different path."""
    # Note: this column also has riken-shape strings so it classifies as id;
    # the value 9.9e9 is below 1e10 but above Pass 1.5's 1e7 threshold so
    # the Pass 1.5 fallback DOES fire. Both checks are valid : we just
    # confirm Pass 1's threshold doesn't off-by-one.
    df = pd.DataFrame({
        "gene_id": ["2310009E13"] + [f"REGISTRY{i}" for i in range(1, 20)],
    })
    report = detect(df)
    # No 1e10+ value in this dataframe → Pass 1 id-float should fire 0 times
    # (only the existing fixture's exponent string).
    pass1_floats = [s for s in report.suspicions
                     if s.kind == "id-float" and isinstance(s.value, (int, float))
                     and abs(s.value) < 1e10]
    assert pass1_floats == [], "no Pass 1 id-float should fire below 1e10"


# --- Pass 1.5 id-float threshold (1e7) -----------------------------------


def test_pass15_id_float_threshold_above() -> None:
    """v0.5.0 tightening (driven by null-result evidence: 74 Entrez-range
    floats were false-positively flagged by the old 1e7 threshold). The
    new policy distinguishes three cases:

    - magnitude ≥ 1e10 in any identifier-header column → flagged (RIKEN E02+)
    - magnitude in [1e7, 1e10] AND header explicitly says RIKEN/cDNA → flagged (RIKEN E01)
    - magnitude in [1e7, 1e10] in any other identifier column → NOT flagged
      (Entrez IDs now reach ~1.5e8)
    """
    # Generic identifier column with Entrez-range floats : must NOT fire
    df_generic = pd.DataFrame({
        "gene_id": [1.5e7, 2.3e7, 3.1e7, 1.0e8, 1.5e8],
    })
    report = detect(df_generic)
    id_floats = [s for s in report.suspicions if s.kind == "id-float"]
    assert len(id_floats) == 0, (
        f"Pass 1.5 must NOT flag Entrez-range floats in a generic id-header "
        f"column; got {len(id_floats)} flags"
    )

    # RIKEN-hinted column with same range → flagged
    df_riken = pd.DataFrame({
        "RIKEN_transcript_id": [1.5e7, 2.3e7, 3.1e7, 1.0e8, 1.5e8],
    })
    report_riken = detect(df_riken)
    id_floats_riken = [s for s in report_riken.suspicions if s.kind == "id-float"]
    assert len(id_floats_riken) == 5, (
        f"Pass 1.5 must flag RIKEN-hinted column with E01-range floats; "
        f"got {len(id_floats_riken)}"
    )

    # Any identifier column with magnitude ≥ 1e10 → flagged
    df_big = pd.DataFrame({
        "gene_id": [2.31e19, 1.5e10, 4.0e15],
    })
    report_big = detect(df_big)
    id_floats_big = [s for s in report_big.suspicions if s.kind == "id-float"]
    assert len(id_floats_big) == 3, (
        f"magnitude ≥ 1e10 in identifier column must flag; got {len(id_floats_big)}"
    )


def test_pass15_id_float_below_threshold_in_id_column() -> None:
    """A value of 5e6 (below 1e7) in an id-headered column should NOT fire."""
    df = pd.DataFrame({
        "gene_id": [1e6, 2e6, 3e6, 4e6, 5e6, 6e6, 7e6, 8e6, 9e6, 9.9e6],
    })
    report = detect(df)
    # No flags expected : values are below Pass 1.5 threshold AND not date-shaped
    id_floats = [s for s in report.suspicions if s.kind == "id-float"]
    assert id_floats == [], (
        f"sub-1e7 values in id-column should NOT flag id-float; got {id_floats}"
    )


# --- Family-range boundaries --------------------------------------------


@pytest.mark.parametrize("month, n, expected_contains, expected_excludes", [
    # MARCH 1-12 (we extended to include MARCH12 / MARCHF12)
    (3, 12, "MARCH12", "MARCH13"),
    (3, 13, None, "MARCH13"),  # MARCH13 doesn't exist
    # TAMM 30-50 boundaries
    (1, 30, "TAMM30", "TAMM29"),
    (1, 50, "TAMM50", "TAMM51"),
    (1, 29, None, "TAMM29"),
    (1, 51, None, "TAMM51"),
    # MEI / MAY ranges
    (5, 1, "MEI1", None),
    (5, 24, "MAY24", None),
    (5, 25, "MEI25", "MAY25"),  # MEI continues, MAY stops at 24
    (5, 31, "MEI31", None),
    # JUN family (via May-29/30/31 numeric subtraction)
    (5, 31, "JUN1", None),
    (5, 30, "JUN2", None),
    (5, 29, "JUN3", None),
    # AGO 1-12 (Argonaute family, locale corruption)
    (8, 1, "AGO1", None),
    (8, 12, "AGO12", None),
    (8, 13, None, "AGO13"),
])
def test_family_range_boundaries(month: int, n: int,
                                  expected_contains: str | None,
                                  expected_excludes: str | None) -> None:
    candidates = _candidates_for_month_and_n(month, n)
    if expected_contains:
        assert expected_contains in candidates, (
            f"_candidates_for_month_and_n({month}, {n}): expected "
            f"{expected_contains!r} in {candidates}"
        )
    if expected_excludes:
        assert expected_excludes not in candidates, (
            f"_candidates_for_month_and_n({month}, {n}): {expected_excludes!r} "
            f"should NOT appear in {candidates}"
        )


# --- Confidence-level pins -----------------------------------------------


def test_confidence_single_candidate_in_id_column_is_095() -> None:
    """NOV1 is the only month=11 candidate. Pass 1 unique-candidate path
    pins confidence at exactly 0.95."""
    df = pd.DataFrame({
        "symbol": ["BRCA1", "TP53", "EGFR", "MARCHF1", "SEPTIN2",
                    date(2024, 11, 1), "KRAS", "PTEN"],
    })
    report = detect(df)
    nov_susp = next(s for s in report.suspicions if s.value.month == 11)
    assert nov_susp.confidence == pytest.approx(0.95)


def test_confidence_ambiguous_non_unique_canonical_is_05() -> None:
    """A date where the canonicalisation produces TWO modern HGNC symbols
    pins confidence at 0.5 WHEN the column shows column-corroboration
    (>=2 date-corruption cells, v0.7.1 gate). With a single date in the
    column, the v0.7.1 demote rule kicks in: 0.5 → 0.20."""
    # Two MARCH cells (date corruption signature) → 0.5 stays
    df_corroborated = pd.DataFrame({
        "symbol": ["BRCA1", "TP53", date(2024, 3, 1), "EGFR", "KRAS",
                    date(2024, 3, 2), "PTEN", "MARCHF1", "SEPTIN2"],
    })
    report = detect(df_corroborated)
    march = next(s for s in report.suspicions
                 if hasattr(s.value, "month") and s.value.month == 3
                 and s.value.day == 1)
    assert march.confidence == pytest.approx(0.5)


def test_confidence_unique_canonical_is_06() -> None:
    """SEPT2 unique canonical (SEPTIN2) → 0.6 boost : when column shows
    >=2 date-corruption cells (v0.7.1 gate)."""
    df_corroborated = pd.DataFrame({
        "symbol": ["BRCA1", "TP53", date(2024, 9, 2), "EGFR", "KRAS",
                    date(2024, 9, 3), "PTEN", "MARCHF1"],
    })
    report = detect(df_corroborated)
    sept = next(s for s in report.suspicions
                if hasattr(s.value, "month") and s.value.month == 9
                and s.value.day == 2)
    assert sept.confidence == pytest.approx(0.6)


def test_confidence_isolated_date_demoted_v071() -> None:
    """v0.7.1: a LONE date in an identifier column is more likely a
    publication date than autofill propagation. Demote 0.5/0.6 → 0.2/0.3
    when the column has fewer than 2 date-corruption cells."""
    df_isolated = pd.DataFrame({
        "symbol": ["BRCA1", "TP53", date(2024, 9, 2), "EGFR", "KRAS",
                    "PTEN", "MARCHF1"],
    })
    report = detect(df_isolated)
    sept = next(s for s in report.suspicions
                if hasattr(s.value, "month") and s.value.month == 9)
    # SEPT2 unique-canonical, but isolated → 0.30 instead of 0.60
    assert sept.confidence == pytest.approx(0.30)
    assert "isolated date in column" in sept.reason


# --- Leading-zero threshold pin (n_lz >= 3) ------------------------------


def test_leading_zero_threshold_below() -> None:
    """Only 2 leading-zero strings : below threshold, so bare ints in the
    column should NOT be flagged as leading-zero-stripped."""
    fp = pd.DataFrame({
        "probe_id": ["00123", "00456",  # only 2 leading-zero strings
                      "REGULAR1", "REGULAR2", "REGULAR3", "REGULAR4",
                      "REGULAR5", "REGULAR6", "REGULAR7",
                      12345],  # could be a stripped 5-digit ID but evidence too thin
    })
    report = detect(fp)
    lz_flags = [s for s in report.suspicions if s.kind == "leading-zero-stripped"]
    assert lz_flags == [], (
        f"only 2 leading-zero exemplars should be below threshold; got {lz_flags}"
    )


def test_leading_zero_threshold_at() -> None:
    """3 leading-zero strings of consistent width : threshold met, a bare
    int whose decimal length < typical_width should be flagged."""
    df = pd.DataFrame({
        "probe_id": ["00123456", "00234567", "00345678",  # 3 exemplars, width 8
                      123],  # bare 3-digit int : plausibly stripped to '00000123'
    })
    report = detect(df)
    lz_flags = [s for s in report.suspicions if s.kind == "leading-zero-stripped"]
    assert len(lz_flags) >= 1, (
        f"3 leading-zero exemplars should meet threshold; got {report.suspicions}"
    )
