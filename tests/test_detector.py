"""Unit tests for the corruption detector."""
from __future__ import annotations

from datetime import date

import pandas as pd
import pytest

from uncorrupt.detector import detect


def test_detector_flags_gene_date_in_named_column() -> None:
    df = pd.DataFrame({
        "symbol": ["BRCA1", "TP53", date(2024, 9, 2), date(2024, 3, 1), "MARCHF1"],
    })
    report = detect(df)
    assert "symbol" in report.identifier_columns
    assert len(report.suspicions) == 2
    sept_susp = next(s for s in report.suspicions if s.value.month == 9)
    # Canonical leads: SEPTIN2 (modern HGNC current), then SEPT2 (deprecated),
    # then SEP2 (historical alias for SELENOF : kept for completeness, NOT
    # silently dropped per the "no silent repair" rule).
    assert "SEPTIN2" in (sept_susp.suggestion or "")
    assert "SEPT2" in (sept_susp.suggestion or "")
    assert "SEP2" in (sept_susp.suggestion or "")
    # +0.1 confidence boost because exactly one canonical resolved
    assert sept_susp.confidence == pytest.approx(0.6)


def test_detector_flags_ambiguous_march_date_with_lower_confidence() -> None:
    df = pd.DataFrame({
        "symbol": ["BRCA1", "TP53", date(2024, 3, 1), "MARCHF1"],
    })
    report = detect(df)
    march_susp = next(s for s in report.suspicions if s.value.month == 3)
    assert "MARCH1" in (march_susp.suggestion or "")
    assert "MARC1" in (march_susp.suggestion or "")
    assert march_susp.confidence < 0.9


def test_detector_flags_float_in_identifier_column() -> None:
    df = pd.DataFrame({
        "gene_id": ["2310009E13"] + [f"REGISTRY{i}" for i in range(1, 20)] + [2.310009e19],
    })
    report = detect(df)
    assert "gene_id" in report.identifier_columns
    float_suspicions = [s for s in report.suspicions if s.kind == "id-float"]
    assert len(float_suspicions) == 1
    assert "precision lost" in float_suspicions[0].reason


def test_detector_ignores_columns_without_identifier_pattern() -> None:
    df = pd.DataFrame({
        "patient_age": [42, 51, 38, 29, 64],
        "diagnosis_date": [date(2023, 1, 5), date(2023, 6, 12), date(2024, 3, 1), date(2024, 9, 2), date(2024, 11, 30)],
    })
    report = detect(df)
    assert report.identifier_columns == []
    assert report.suspicions == []


def test_detector_uses_header_hint() -> None:
    df = pd.DataFrame({
        "gene": [date(2024, 9, 2)],
    })
    report = detect(df)
    assert "gene" in report.identifier_columns
    assert len(report.suspicions) == 1


def test_detector_decodes_mmm_yy_format_via_year_suffix() -> None:
    """When openpyxl supplies a `mmm-yy` format, the date(2007, 9, 1) cell came
    from Excel auto-converting "SEPT7" → year 2007. Decode via year-suffix,
    not day-of-month."""
    df = pd.DataFrame({
        "symbol": ["BRCA1", "TP53", date(2007, 9, 1), "MARCHF1", "EGFR", "KRAS"],
    })
    df.attrs["cell_formats"] = {("symbol", 2): "mmm-yy"}
    report = detect(df)
    assert "symbol" in report.identifier_columns
    sept = next(s for s in report.suspicions if s.value.month == 9)
    assert "SEPT7" in (sept.suggestion or "")
    assert "SEP7" in (sept.suggestion or "")
    # Should NOT have offered SEPT1 : the format unambiguously says year-suffix
    assert "SEPT1" not in (sept.suggestion or "")
    assert "mmm-yy" in sept.reason


def test_detector_decodes_dd_mmm_format_via_day_of_month() -> None:
    """`d-mmm` format means day shown / year implicit. Stay with day-of-month."""
    df = pd.DataFrame({
        "symbol": ["BRCA1", "TP53", date(2024, 9, 7), "MARCHF1", "EGFR", "KRAS"],
    })
    df.attrs["cell_formats"] = {("symbol", 2): "d-mmm"}
    report = detect(df)
    sept = next(s for s in report.suspicions if s.value.month == 9)
    assert "SEPT7" in (sept.suggestion or "")
    assert "SEP7" in (sept.suggestion or "")


def test_detector_unknown_format_emits_both_when_day_is_one_and_year_in_window() -> None:
    """Without format info, day=1 + year in 2000-2029 is the year-suffix
    signature. Emit both interpretations so the human can disambiguate."""
    df = pd.DataFrame({
        "symbol": ["BRCA1", "TP53", date(2007, 9, 1), "MARCHF1", "EGFR", "KRAS"],
    })
    # No df.attrs["cell_formats"] set
    report = detect(df)
    sept = next(s for s in report.suspicions if s.value.month == 9)
    sugg = sept.suggestion or ""
    # year-suffix candidate
    assert "SEPT7" in sugg and "SEP7" in sugg
    # day-of-month candidate (the older interpretation)
    assert "SEPT1" in sugg and "SEP1" in sugg


def test_detector_unknown_format_outside_year_window_uses_day_only() -> None:
    """day=1 + year=2024 is too recent to be a year-suffix corruption (no
    SEPT24 gene). Stay with day-of-month interpretation only."""
    df = pd.DataFrame({
        "symbol": ["BRCA1", "TP53", date(2024, 9, 1), "MARCHF1", "EGFR", "KRAS"],
    })
    report = detect(df)
    sept = next(s for s in report.suspicions if s.value.month == 9)
    sugg = sept.suggestion or ""
    assert "SEPT1" in sugg
    assert "SEPT24" not in sugg  # SEPT24 doesn't exist as a gene


def test_detector_parses_iso_date_string_in_gene_column() -> None:
    """ISO `'2012-03-03 00:00:00'` strings appear when an Excel date round-trips
    through CSV/TSV : the date-string parser must recognise them."""
    df = pd.DataFrame({
        "symbol": ["BRCA1", "TP53", "2012-03-03 00:00:00", "MARCHF1", "EGFR", "KRAS"],
    })
    report = detect(df)
    march = next(s for s in report.suspicions
                 if s.kind == "gene-date-string" and "2012-03-03" in str(s.value))
    assert "MARCH3" in (march.suggestion or "")


def test_detector_parses_iso_date_string_no_time_component() -> None:
    df = pd.DataFrame({
        "symbol": ["BRCA1", "TP53", "2007-09-07", "MARCHF1", "EGFR", "KRAS"],
    })
    report = detect(df)
    sept = next(s for s in report.suspicions if s.kind == "gene-date-string")
    # day=7 (day-of-month interpretation since unambiguous ISO date)
    assert "SEPT7" in (sept.suggestion or "")
