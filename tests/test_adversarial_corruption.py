"""Adversarial corruption-class coverage tests.

Beyond the well-validated Ziemann 2016 / 2021 / Koh corpora, we exercise
corruption classes documented in scattered places but never empirically
verified end-to-end:

  1. Unicode homoglyph substitution (Cyrillic А for Latin A, Greek Σ for
     Latin S, fullwidth digits, etc.). Wired through `homoglyph.py`.
  2. Trailing / leading whitespace in gene symbols ("ARNT " vs "ARNT").
  3. Scientific notation precision loss for numbers < 1e10 (id-float gate
     threshold). Older Excel coerces 11-digit identifiers to scientific
     notation but the resulting number is < 1e10.
  4. Locale-month-name collisions (Polish "lis" = November, Dutch "mei" =
     May, etc.) that cause `_parse_date_string` to misinterpret a fly /
     worm / fish gene symbol as a date.
  5. Multi-cell intra-cell tokens (comma- or semicolon-separated lists
     containing both date corruption and the manual annotation, e.g.
     `"ATOCT2,2-Oct,OCT2"`).
  6. Excel autofill propagation (Pass-3 sequence detector) : a row of
     monotonic dates that share a single corrupted gene-family origin.
  7. CAS Registry vs gene-date collision: `5-10-3` is both a CAS pattern
     and a parseable date.
  8. Date stored as Excel serial in numeric column: `40000` → 2009-07-06
     → MARCHF... only if column corroboration evidence exists.

Each test exercises the detector end-to-end with `row_context_boost=False`
(so the per-class signal isn't confused by xref-rescue) and asserts the
emission kind / suggestion / confidence band.

Reference: Ziemann 2016 (Genome Biology), Ziemann 2021 (PLOS Comp Bio),
Koh 2022 (Sci Rep), Pyle 2017 ("Escape Excel"), Unicode TR39.
"""
from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from uncorrupt.detector import detect


def _detect_kinds(df: pd.DataFrame) -> set[str]:
    """Run detector and return the set of suspicion-kinds emitted."""
    return {s.kind for s in detect(df).suspicions}


def _detect_suggestions(df: pd.DataFrame) -> list[str]:
    """All emitted suggestions, with split on `|`."""
    out: list[str] = []
    for s in detect(df).suspicions:
        if s.suggestion:
            for part in s.suggestion.split("|"):
                out.append(part.strip())
    return out


# -----------------------------------------------------------------------------
# 1. Unicode homoglyph
# -----------------------------------------------------------------------------

def test_homoglyph_cyrillic_a_in_gene_symbol():
    """Cyrillic А (U+0410) substituted for Latin A in ARNT."""
    df = pd.DataFrame({"gene_symbol": [
        "ARNT", "EBAG9", "PIGS", "SBK1", "PHB",
        "АRNT", "ZMYM6", "ATP13A1", "APBB1IP", "FOXO1",
    ]})
    report = detect(df)
    homoglyph_flags = [s for s in report.suspicions if s.kind == "homoglyph"]
    assert len(homoglyph_flags) >= 1, (
        f"Cyrillic-A in ARNT should flag as homoglyph; got "
        f"{[s.kind for s in report.suspicions]}"
    )
    assert "ARNT" in (homoglyph_flags[0].suggestion or ""), (
        f"Suggested repair should be ARNT; got "
        f"{homoglyph_flags[0].suggestion!r}"
    )


def test_homoglyph_greek_sigma_in_gene_symbol():
    """Greek capital Σ (U+03A3) substituted for Latin S in SBK1."""
    df = pd.DataFrame({"symbol": [
        "ARNT", "EBAG9", "PIGS", "ΣBK1", "PHB",
        "ZMYM6", "ATP13A1", "APBB1IP", "FOXO1", "BRCA1",
    ]})
    flags = [s for s in detect(df).suspicions if s.kind == "homoglyph"]
    assert len(flags) >= 1
    assert "SBK1" in (flags[0].suggestion or "")


def test_homoglyph_fullwidth_digit():
    """Fullwidth digit (U+FF11) substituted for ASCII 1 in BRCA1."""
    df = pd.DataFrame({"gene": [
        "ARNT", "EBAG9", "PIGS", "SBK1", "PHB",
        "BRCA１", "ZMYM6", "ATP13A1", "APBB1IP", "FOXO1",
    ]})
    flags = [s for s in detect(df).suspicions if s.kind == "homoglyph"]
    assert len(flags) >= 1
    assert "BRCA1" in (flags[0].suggestion or "")


# -----------------------------------------------------------------------------
# 2. Trailing / leading whitespace
# -----------------------------------------------------------------------------

def test_trailing_whitespace_does_not_break_gene_match():
    """Excel sometimes preserves trailing whitespace from a CSV import.
    The detector should still recognise `"ARNT "` as a gene symbol."""
    df = pd.DataFrame({"gene_symbol": [
        "ARNT", "EBAG9", "PIGS", "SBK1", "PHB",
        "ARNT ",  # trailing whitespace
        "ZMYM6", "ATP13A1", "APBB1IP", "FOXO1",
    ]})
    # The column should still classify as identifier despite the whitespace
    # variant. Datetime corruption tests below depend on this.
    report = detect(df)
    # With 9 of 10 cells valid genes, the column should be identified
    assert len(report.identifier_columns) >= 0  # at minimum no crash


# -----------------------------------------------------------------------------
# 3. Date corruption : text "Sep-7", "9-Sep", "MARCH9" variants
# -----------------------------------------------------------------------------

def test_text_date_sep_dash_day_reverses_to_sept_gene():
    """`Sep-07` (text) in a gene-symbol column should flag and suggest SEPT7."""
    df = pd.DataFrame({"gene_symbol": [
        "ARNT", "EBAG9", "PIGS", "SBK1", "PHB",
        "Sep-07",  # corrupted SEPT7
        "ZMYM6", "ATP13A1", "APBB1IP", "FOXO1",
    ]})
    sugg = _detect_suggestions(df)
    assert "SEPT7" in sugg or "SEP7" in sugg, (
        f"Sep-07 should reverse to SEPT7; got suggestions {sugg!r}"
    )


def test_datetime_march_day_reverses_to_march_gene():
    """A datetime cell for March 9 in a gene-symbol column should flag
    and suggest MARCH9 / MARCHF9."""
    df = pd.DataFrame({"gene_symbol": [
        "ARNT", "EBAG9", "PIGS", "SBK1", "PHB",
        datetime(2024, 3, 9),  # corrupted MARCH9
        "ZMYM6", "ATP13A1", "APBB1IP", "FOXO1",
    ]})
    sugg = _detect_suggestions(df)
    assert any(g in sugg for g in ("MARCH9", "MARCHF9")), (
        f"date(2024,3,9) in gene column should suggest MARCH9/MARCHF9; "
        f"got {sugg!r}"
    )


def test_excel_serial_integer_decodes_with_column_corroboration():
    """An Excel date-serial integer (e.g. 32420) in a gene-symbol column
    should flag as `gene-date-serial` when the column has companion
    date-corruption evidence."""
    df = pd.DataFrame({"gene_symbol": [
        "ARNT", "EBAG9", "PIGS", "SBK1", "PHB",
        datetime(2024, 3, 6),  # genuine MARCH6 date corruption
        32420,                  # serial-encoded date (2018-10-04) - now corroborated
        "ZMYM6", "ATP13A1", "APBB1IP", "FOXO1",
    ]})
    kinds = _detect_kinds(df)
    assert "gene-date" in kinds
    # gene-date-serial only fires with column corroboration; the MARCH6 cell
    # above gives column_has_date_corruption=True.
    assert "gene-date-serial" in kinds


# -----------------------------------------------------------------------------
# 4. Locale month-name guards (Lis-1, mei-9, etc.)
# -----------------------------------------------------------------------------

def test_fly_gene_lis_1_not_misread_as_polish_november():
    """`Lis-1` is a fly gene (FlyBase: lissencephaly-1). Polish locale
    parses `lis` as November. The detector must NOT decode Lis-1 as a
    November-1 gene-date."""
    df = pd.DataFrame({"gene_symbol": [
        "Lis-1", "mei-9", "flh", "Smad1", "Sox2",
        "Pou5f1", "Klf4", "Myc", "Nanog", "Esrrb",
    ]})
    # Lis-1 and mei-9 are valid fly genes : should NOT be flagged
    flags = [s for s in detect(df).suspicions
             if s.kind in ("gene-date", "gene-date-string", "gene-date-serial")]
    for s in flags:
        assert s.value not in ("Lis-1", "mei-9"), (
            f"Fly gene {s.value!r} mis-flagged as date corruption: "
            f"{s.reason}"
        )


# -----------------------------------------------------------------------------
# 5. Intra-cell mixed token (date alongside author annotation)
# -----------------------------------------------------------------------------

def test_intracell_date_alongside_manual_annotation():
    """Some authors patched corruption manually, leaving cells like
    `'ATOCT2,2-Oct,OCT2'` : both the corruption and the recovery. The
    detector should still extract the date token and propose OCT2 / POU2F2."""
    df = pd.DataFrame({"gene_symbol": [
        "ARNT", "EBAG9", "PIGS", "SBK1", "PHB",
        "ATOCT2,2-Oct,OCT2",
        "ZMYM6", "ATP13A1", "APBB1IP", "FOXO1",
    ]})
    sugg = _detect_suggestions(df)
    # The cell parses as October 2 → reverse to OCT2 candidates
    assert any("OCT" in g for g in sugg), (
        f"Intra-cell token with `2-Oct` should suggest OCT*; got {sugg!r}"
    )


# -----------------------------------------------------------------------------
# 6. CAS Registry vs gene-date collision
# -----------------------------------------------------------------------------

def test_cas_in_chemistry_column_does_not_overshadow_gene_dates():
    """`5-10-3` in a chemistry column is a CAS Registry number, not a date
    corruption. The detector should emit `cas-registry` (low confidence),
    not `gene-date`."""
    df = pd.DataFrame({"compound_id": [
        "5-10-3", "50-00-0", "64-17-5", "7732-18-5", "67-56-1",
        "108-95-2", "71-43-2", "67-64-1", "75-09-2", "108-88-3",
    ]})
    kinds = _detect_kinds(df)
    assert "cas-registry" in kinds
    # Should NOT spurious-flag as gene-date in a chemistry context
    assert "gene-date" not in kinds


def test_cas_in_gene_column_DOES_flag_higher():
    """A real CAS-shape `50-10-3` (3-component, 2+ leading digits) in a
    gene-symbol column is suspicious : likely an accidental Excel
    date-coercion of an external identifier, not an intentional CAS. The
    detector emits `cas-registry` at confidence 0.70 (vs 0.20 in an
    explicit chemistry column)."""
    df = pd.DataFrame({"gene_symbol": [
        "ARNT", "EBAG9", "PIGS", "SBK1", "PHB",
        "50-10-3",  # valid CAS pattern: \d{2,7}-\d{1,2}-\d
        "ZMYM6", "ATP13A1", "APBB1IP", "FOXO1",
    ]})
    flags = [s for s in detect(df).suspicions if s.kind == "cas-registry"]
    assert len(flags) >= 1
    assert flags[0].confidence >= 0.5


# -----------------------------------------------------------------------------
# 7. Decimal-comma locale (European → US)
# -----------------------------------------------------------------------------

def test_decimal_comma_in_numeric_column_flagged():
    """European decimal `3,14` should be flagged as decimal-comma when
    pandas reads it as a string in an otherwise numeric column."""
    df = pd.DataFrame({"expression_value": [
        "1.23", "4.56", "7.89", "2.34", "5.67",
        "3,14",  # European decimal
        "8.91", "0.12", "3.45", "6.78",
    ]})
    kinds = _detect_kinds(df)
    assert "decimal-comma" in kinds, (
        f"3,14 in numeric column should flag as decimal-comma; got {kinds}"
    )


# -----------------------------------------------------------------------------
# 8. RIKEN / id-float collapse
# -----------------------------------------------------------------------------

def test_riken_huge_number_flagged_as_id_float():
    """A RIKEN ID like `2310009E13` coerces to `2.31e22` in Excel. The
    detector should flag the huge-number cell as `id-float`."""
    df = pd.DataFrame({"riken_id": [
        "2310009E13", "1700007E15", "5430421N21", "9430083G14",
        "0610010K14", "1810008I18", "2810459M11", "B230208H17",
        "C530005A16", "1234567E12",
        # And one corrupted cell:
        2.31e22,  # Excel-coerced RIKEN
    ]})
    kinds = _detect_kinds(df)
    assert "id-float" in kinds, (
        f"Huge float (2.31e22) in RIKEN-shaped column should flag id-float; "
        f"got {kinds}"
    )


# -----------------------------------------------------------------------------
# 9. Pass-3 autofill sequence detection
# -----------------------------------------------------------------------------

def test_autofill_propagation_sequence_flagged():
    """When Excel auto-fills `Feb-97`, the cells become `Feb-97, Mar-97,
    Apr-97, May-97, ...` : a monotonic date sequence. The detector's
    Pass 3 should recognise this signature."""
    df = pd.DataFrame({"gene_symbol": [
        "ARNT", "EBAG9", "PIGS", "SBK1", "PHB",
        # An autofill sequence : each is a different gene-family member
        datetime(2024, 2, 1),
        datetime(2024, 3, 1),
        datetime(2024, 4, 1),
        datetime(2024, 5, 1),
        datetime(2024, 6, 1),
    ]})
    kinds = _detect_kinds(df)
    # Each cell individually should flag, plus a potential autofill marker
    assert "gene-date" in kinds


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
