"""Unit tests for the Excel-corruption simulator."""
from __future__ import annotations

from datetime import date

import pytest

from uncorrupt.simulator import corrupt_column, excel_autoconvert


@pytest.mark.parametrize(
    "input_value,expected_month,expected_day",
    [
        ("SEPT1", 9, 1),
        ("SEPT2", 9, 2),
        ("SEPT15", 9, 15),
        ("MARCH1", 3, 1),
        ("MARCH11", 3, 11),
        ("MARC1", 3, 1),
        ("DEC1", 12, 1),
    ],
)
def test_known_genes_become_dates(input_value: str, expected_month: int, expected_day: int) -> None:
    result = excel_autoconvert(input_value)
    assert isinstance(result, date), f"{input_value!r} should be a date, got {type(result).__name__}"
    assert result.month == expected_month
    assert result.day == expected_day


def test_riken_id_becomes_float() -> None:
    result = excel_autoconvert("2310009E13")
    assert isinstance(result, float)
    assert result == pytest.approx(2.310009e19)


def test_simple_exponent_string_becomes_float() -> None:
    assert excel_autoconvert("1E10") == 1e10
    assert excel_autoconvert("12E5") == 12e5


def test_normal_gene_symbols_unchanged() -> None:
    assert excel_autoconvert("BRCA1") == "BRCA1"
    assert excel_autoconvert("TP53") == "TP53"
    assert excel_autoconvert("MARCHF1") == "MARCHF1"
    assert excel_autoconvert("SEPTIN2") == "SEPTIN2"


def test_non_string_passthrough() -> None:
    assert excel_autoconvert(42) == 42
    assert excel_autoconvert(None) is None
    assert excel_autoconvert(3.14) == 3.14


def test_corrupt_column_round_trip() -> None:
    inputs = ["BRCA1", "SEPT2", "TP53", "2310009E13", "DEC1"]
    result = corrupt_column(inputs)
    assert result[0] == "BRCA1"
    assert isinstance(result[1], date) and result[1].month == 9 and result[1].day == 2
    assert result[2] == "TP53"
    assert isinstance(result[3], float)
    assert isinstance(result[4], date) and result[4].month == 12 and result[4].day == 1


# --- Locale-specific simulator coverage ---
# Each row pins a corruption mode that the detector now reverses but the
# simulator previously could not produce, leaving the simulator/detector
# pair untestable as a round trip.


@pytest.mark.parametrize(
    "gene,expected_month,expected_day",
    [
        # AGO family (Italian/Spanish/Portuguese locale, Ziemann 2021)
        ("AGO2", 8, 2),
        ("AGO3", 8, 3),
        ("AGO4", 8, 4),
        # MEI family (Dutch locale)
        ("MEI1", 5, 1),
        ("MEI15", 5, 15),
        # MAY family (May plant species)
        ("MAY1", 5, 1),
        ("MAY24", 5, 24),
        # English MARCHF12 (off-by-one fix, v0.4.1)
        ("MARCH12", 3, 12),
        # OCT and APR (English)
        ("OCT4", 10, 4),
        ("APR1", 4, 1),
        ("APR2", 4, 2),
        # FEB family : HGNC has only FEB3 and FEB4
        ("FEB3", 2, 3),
        ("FEB4", 2, 4),
        # NOV / DEC2
        ("NOV1", 11, 1),
        ("DEC2", 12, 2),
    ],
)
def test_locale_genes_become_dates(
    gene: str, expected_month: int, expected_day: int,
) -> None:
    result = excel_autoconvert(gene)
    assert isinstance(result, date), f"{gene!r} should be a date, got {type(result).__name__}"
    assert result.month == expected_month
    assert result.day == expected_day


@pytest.mark.parametrize("gene,expected_year", [
    ("TAMM30", 2030),
    ("TAMM47", 2047),
    ("TAMM50", 2050),
])
def test_tamm_finnish_locale_becomes_year_dated(gene: str, expected_year: int) -> None:
    """TAMM30..TAMM50 in Finnish locale renders as January 1 of years
    2030-2050 (day-of-month encoded as year suffix)."""
    result = excel_autoconvert(gene)
    assert isinstance(result, date)
    assert result.year == expected_year
    assert result.month == 1
    assert result.day == 1


@pytest.mark.parametrize("gene,expected_month,expected_day", [
    # JUN family: numeric-subtraction quirk. "jun-N" = June - N days
    ("JUN1", 5, 31),
    ("JUN2", 5, 30),
    ("JUN3", 5, 29),
])
def test_jun_numeric_subtraction(
    gene: str, expected_month: int, expected_day: int,
) -> None:
    result = excel_autoconvert(gene)
    assert isinstance(result, date)
    assert result.month == expected_month
    assert result.day == expected_day


def test_simulator_detector_round_trip_locale_families() -> None:
    """Round-trip property: simulator corrupts a known gene → detector
    reverses corrupted date back to the gene candidate.

    Without this end-to-end test, simulator and detector can drift apart
    silently : the simulator's locale coverage caught up to the detector's
    in v0.4.3."""
    from uncorrupt.detector import _reverse_gene_date

    for gene in ["AGO2", "AGO3", "MEI1", "MAY15", "TAMM30",
                 "JUN1", "MARCH12", "OCT4", "APR2", "FEB3", "FEB4"]:
        corrupted = excel_autoconvert(gene)
        assert isinstance(corrupted, date), f"{gene!r} did not corrupt"
        candidates = _reverse_gene_date(corrupted)
        assert gene in candidates, (
            f"round trip broke: {gene!r} → {corrupted} → {candidates} "
            f"(does not contain {gene!r})"
        )
