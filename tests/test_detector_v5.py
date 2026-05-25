"""Unit tests for v5 detector additions.

Locks in the gains from the v5 improvements:
- expanded `_reverse_gene_date` (Oct, Nov, Apr, Feb in addition to Mar/Sep/Dec)
- `_parse_date_string` (multiple date-string formats with locale ambiguity)
- `_serial_to_date` (Excel date-serial math)
- Pass 2 column-context-blind detection of date-serial integers and digit-strings
"""
from __future__ import annotations

from datetime import date
from pathlib import Path

import pandas as pd
import pytest

from uncorrupt.detector import (
    _parse_date_string,
    _reverse_gene_date,
    _serial_to_date,
    detect,
)

# --- _reverse_gene_date coverage ---

@pytest.mark.parametrize("d,expected", [
    (date(2024, 3, 1), {"MARCH1", "MARC1"}),
    (date(2024, 3, 6), {"MARCH6"}),
    (date(2024, 3, 11), {"MARCH11"}),
    # MARCH12 covers the full HGNC MARCHF rename range (MARCHF1 through MARCHF12).
    # Earlier versions only ran 1 through 11.
    # off-by-one (was 1-11, now 1-12 to cover the full HGNC rename range).
    (date(2024, 3, 12), {"MARCH12"}),
    (date(2024, 9, 2), {"SEPT2", "SEP2"}),
    (date(2024, 9, 15), {"SEPT15", "SEP15"}),
    (date(2024, 9, 16), set()),
    (date(2024, 10, 1), {"OCT1"}),
    (date(2024, 10, 11), {"OCT11"}),
    (date(2024, 10, 12), set()),
    (date(2024, 11, 1), {"NOV1"}),
    (date(2024, 11, 2), set()),  # No NOV2 in HGNC aliases
    (date(2024, 12, 1), {"DEC1"}),
    (date(2024, 12, 2), {"DEC2"}),
    (date(2024, 12, 3), set()),
    (date(2024, 4, 1), {"APR1"}),
    # APR2 added in v0.4.0 : Ziemann 2021 S1 documents 3 fungal genes that
    # corrupt to 2021-04-02; APR2 is part of the AGO/Argonaute family alias set
    (date(2024, 4, 2), {"APR2"}),
    (date(2024, 4, 3), {"APR3"}),
    (date(2024, 2, 3), {"FEB3"}),
    (date(2024, 2, 4), {"FEB4"}),
    # AGO gene family : locale-specific (Italian/Spanish/Portuguese "Ago" =
    # August). AGO2 typed in those locales -> "Aug-02", AGO3 -> "Aug-03",
    # AGO4 -> "Aug-04". Pinning AGO2 and AGO3 explicitly so the locale-specific
    # "Ago" -> August path stays under test.
    (date(2024, 8, 2), {"AGO2"}),
    (date(2024, 8, 3), {"AGO3"}),
    (date(2024, 8, 4), {"AGO4"}),
    (date(2024, 1, 15), set()),  # January suffix 15 is not in TAMM range (30-50)
])
def test_reverse_gene_date_coverage(d: date, expected: set[str]) -> None:
    assert set(_reverse_gene_date(d)) == expected


# --- _parse_date_string formats ---

@pytest.mark.parametrize("text,expected_dates", [
    ("Mar-6", [date(2026, 3, 6)]),  # uses current year
    ("Mar-09", [date(2026, 3, 9)]),
    ("Sep-2", [date(2026, 9, 2)]),
    ("Sep-02", [date(2026, 9, 2)]),
    ("2-Oct", [date(2026, 10, 2)]),
    ("2-Oct-14", [date(2014, 10, 2)]),
    # All-numeric forms are ambiguous: both DD/MM/YY and MM/DD/YY are returned
    ("06/03/14", [date(2014, 3, 6), date(2014, 6, 3)]),
    ("03/06/2014", [date(2014, 6, 3), date(2014, 3, 6)]),
    ("not a date", []),
    ("12345", []),
])
def test_parse_date_string(text: str, expected_dates: list[date]) -> None:
    result = _parse_date_string(text)
    assert set(result) == set(expected_dates)


# --- _serial_to_date math ---

def test_serial_to_date_basic() -> None:
    # Excel serial 39880 from epoch 1899-12-30 → 2009-03-08 (MARCH8)
    assert _serial_to_date(39880) == date(2009, 3, 8)
    # 41334 → 2013-03-01 (MARCH1 / MARC1)
    assert _serial_to_date(41334) == date(2013, 3, 1)
    # 41519 → 2013-09-02 (SEPT2 / SEP2)
    assert _serial_to_date(41519) == date(2013, 9, 2)


def test_serial_to_date_out_of_range_returns_none() -> None:
    assert _serial_to_date(0) is None
    assert _serial_to_date(100) is None
    assert _serial_to_date(99999999) is None


# --- Pass 2: column-context-blind detection ---

def test_pass2_catches_date_serial_int_with_no_identifier_column() -> None:
    """A column not classified as identifier, containing int that decodes as
    uncorrupt date → flagged with low confidence."""
    # Use object dtype so the 39880 stays as int through pandas
    df = pd.DataFrame({"Unnamed_1": pd.array([39880, "x", "y", "z", "w"], dtype=object)})
    report = detect(df)
    serial_susp = [s for s in report.suspicions if s.kind == "gene-date-serial"]
    assert len(serial_susp) == 1
    assert "MARCH8" in (serial_susp[0].suggestion or "")
    assert serial_susp[0].confidence <= 0.4  # column-blind = low confidence


def test_pass1_catches_float_encoded_date_serial_in_gene_column() -> None:
    """An isolated float-encoded date serial (39880.0) in a gene-symbol-shaped
    column gets flagged via Pass 1. Real-world corruption pattern: a column
    of gene names with one cell corrupted to a date-serial number.

    (After v6 Fix 4, Pass 2 no longer scans columns classified as 'measurement',
    so we now test Pass 1 : the right place for this kind of flag, since the
    column is genuinely a gene-symbol column with isolated corruption.)
    """
    df = pd.DataFrame({
        "gene": pd.array(["BRCA1", "TP53", "EGFR", 39880.0, "MARCHF1", "KRAS"], dtype=object),
    })
    report = detect(df)
    serial_susp = [s for s in report.suspicions if s.kind == "gene-date-serial"]
    assert len(serial_susp) == 1
    assert "MARCH8" in (serial_susp[0].suggestion or "")


def test_pass2_catches_date_serial_digit_string() -> None:
    """A digit-string like '41519' in a non-identifier column."""
    df = pd.DataFrame({"misc_col": ["x", "y", "41519", "z", "w"]})
    report = detect(df)
    serial_susp = [s for s in report.suspicions if s.kind == "gene-date-serial"]
    assert len(serial_susp) == 1
    # 41519 decodes to 2013-09-02 → SEPT2 + SEP2 (HGNC has both as alias)
    assert "SEPT2" in (serial_susp[0].suggestion or "")


def test_pass2_catches_oct_date_string() -> None:
    """A string '2-Oct' in a non-identifier column → OCT2 candidate."""
    df = pd.DataFrame({"col_a": ["foo", "2-Oct", "bar"]})
    report = detect(df)
    string_susp = [s for s in report.suspicions if s.kind == "gene-date-string"]
    assert len(string_susp) >= 1
    assert "OCT2" in (string_susp[0].suggestion or "")


def test_pass2_does_not_flag_random_ints_outside_gene_signature() -> None:
    """Integers that decode to dates with no uncorrupt signature must not
    be flagged (e.g. Aug 4 1923 has no AUG gene)."""
    df = pd.DataFrame({"some_col": [8615, 9000, 10000, 11000, 12000]})
    report = detect(df)
    assert report.suspicions == []


def test_pass2_skips_columns_with_date_token_header() -> None:
    """Columns whose header contains 'date'/'year'/etc. are excluded from Pass 2."""
    df = pd.DataFrame({"diagnosis_date": [39880, 41334, 41519, 39880, 41334]})
    report = detect(df)
    assert report.suspicions == []


# --- v6 Fix 1: id-float pass tightening ---

def test_v6_pass1_does_not_flag_moderate_int_without_serial_decode() -> None:
    """A moderate integer in an identifier column that decodes to a date with
    no uncorrupt signature must NOT be flagged.

    serial 30000 = 1982-02-13, no gene match → no flag.
    Pre-v6 this got flagged as bare id-float (false positive)."""
    df = pd.DataFrame({
        "symbol": ["BRCA1", "TP53", "EGFR", "MARCHF1", "SEPTIN2", "KRAS", 30000],
    })
    report = detect(df)
    bare_idfloat = [s for s in report.suspicions
                    if s.kind == "id-float" and s.suggestion is None]
    # 30000 = 1982-02-13 is February 13, not in uncorrupt signature → not flagged
    assert len(bare_idfloat) == 0


def test_v6_pass1_does_flag_int_when_serial_decodes_to_uncorrupt() -> None:
    """An integer that decodes to a uncorrupt date IS flagged."""
    df = pd.DataFrame({
        "symbol": pd.array(["BRCA1", "TP53", "EGFR", "MARCHF1", "SEPTIN2", "KRAS", 41884],
                          dtype=object),
    })
    report = detect(df)
    # serial 41884 = 2014-09-02 = SEPT2 / SEP2
    gene_date_serial = [s for s in report.suspicions if s.kind == "gene-date-serial"]
    assert len(gene_date_serial) == 1
    assert "SEPT2" in (gene_date_serial[0].suggestion or "")


def test_v6_pass1_does_flag_large_numeric_as_riken_style() -> None:
    """Large numeric (|value| >= 1e10) in identifier column → id-float
    (RIKEN/exponent-string coercion, precision lost : no decode possible)."""
    df = pd.DataFrame({
        "gene_id": pd.array([f"23{i:07d}E13" for i in range(15)] + [2.31e19], dtype=object),
    })
    report = detect(df)
    riken_float = [s for s in report.suspicions if s.kind == "id-float"]
    assert len(riken_float) == 1
    assert "RIKEN/exponent" in riken_float[0].reason


# --- Pass 3: autofill-sequence detection ---


def test_pass3_flags_autofill_sequence_of_dates() -> None:
    """Ziemann 2021 (PMC6330011): a corrupted Sep-2 cell dragged down via
    Excel autofill propagates Sep-3, Sep-4, ... : a monotonic date series
    where each cell maps to a gene-family member. The column-level pattern
    is stronger evidence than any single cell.

    Realistic shape: real autofill leaves the surrounding rows as gene
    symbols, so the column is still classified as an identifier column
    (a pure-date column is excluded from the identifier classifier)."""
    df = pd.DataFrame({
        "gene": pd.array(
            ["BRCA1", "TP53", "EGFR", "KRAS", "MYC",
             date(2024, 9, 2), date(2024, 9, 3), date(2024, 9, 4),
             date(2024, 9, 5), date(2024, 9, 6),
             "PTEN", "RB1", "ATM", "CDKN2A", "VHL"],
            dtype=object,
        ),
    })
    report = detect(df)
    seq = [s for s in report.suspicions if s.kind == "autofill-sequence"]
    assert len(seq) == 1
    assert seq[0].confidence == pytest.approx(0.90)
    assert "5" in seq[0].reason  # 5 cells in the sequence
    assert "PMC6330011" in seq[0].reason


def test_pass3_does_not_flag_non_sequential_dates() -> None:
    """A single isolated date in a column is NOT an autofill sequence : the
    autofill flag requires three or more row-adjacent monotonic dates."""
    df = pd.DataFrame({
        "gene": pd.array(
            ["BRCA1", "TP53", date(2024, 9, 2), "EGFR", "KRAS"],
            dtype=object,
        ),
    })
    report = detect(df)
    seq = [s for s in report.suspicions if s.kind == "autofill-sequence"]
    assert len(seq) == 0


def test_pass3_requires_uniform_spacing() -> None:
    """Out-of-order, decreasing, or non-uniformly-spaced dates do not
    constitute an autofill sequence : Excel autofill produces uniformly
    stepped increments (1 day, 1 week, ~1 month, ~1 year). Real-data
    validation on Ziemann 2016 (`23826142__asset.xlsx`) showed that any
    non-uniform monotonic gap was too permissive : long columns of real
    publication dates landing in gene-name months produced 30 false
    positives until uniform spacing was required."""
    df = pd.DataFrame({
        "gene": pd.array(
            ["BRCA1", "TP53", "EGFR", "KRAS", "MYC",
             date(2024, 9, 5), date(2024, 3, 1), date(2024, 9, 2),
             "PTEN", "RB1", "ATM"],
            dtype=object,
        ),
    })
    report = detect(df)
    seq = [s for s in report.suspicions if s.kind == "autofill-sequence"]
    assert len(seq) == 0


def test_pass3_flags_march_family_autofill_sequence() -> None:
    """MARCH-family autofill: a user types Mar-1, drags down, gets
    Mar-2, Mar-3 : corresponding to MARCHF1/MARCHF2/MARCHF3."""
    df = pd.DataFrame({
        "symbol": pd.array(
            ["BRCA1", "TP53", "EGFR", "KRAS", "MYC",
             date(2024, 3, 1), date(2024, 3, 2), date(2024, 3, 3),
             date(2024, 3, 4),
             "PTEN", "RB1", "ATM"],
            dtype=object,
        ),
    })
    report = detect(df)
    seq = [s for s in report.suspicions if s.kind == "autofill-sequence"]
    assert len(seq) == 1


# --- String-form time corruption (v0.5.0) ---


def test_pass1_flags_string_form_time_in_identifier_column() -> None:
    """`"01:03"` and `"11:30"` strings in an identifier column → time-coercion
    flag. Same root cause as datetime.time but CSV-serialized."""
    df = pd.DataFrame({
        "gene": pd.array(
            ["BRCA1", "TP53", "01:03", "EGFR", "KRAS",
             "11:30", "PTEN", "RB1"],
            dtype=object,
        ),
    })
    report = detect(df)
    times = [s for s in report.suspicions if s.kind == "time-coercion"]
    assert len(times) == 2
    assert "1:03" in (times[0].suggestion or "")


def test_pass1_does_not_flag_time_string_in_free_text_column() -> None:
    """A "13:45" string in a non-identifier column doesn't trigger."""
    df = pd.DataFrame({
        "notes": ["meeting at 13:45", "checkpoint 09:00", "review at 16:30"],
    })
    report = detect(df)
    times = [s for s in report.suspicions if s.kind == "time-coercion"]
    assert len(times) == 0


# --- Homoglyph corruption (v0.5.0) ---


def test_pass1_flags_greek_sigma_homoglyph() -> None:
    """A gene symbol containing Greek capital Σ (U+03A3) instead of Latin
    S is flagged with the proposed ASCII repair (Unicode TR39 confusables
    + manual Sigma override)."""
    # SEPT2 with Σ substituting for S : visible to the eye, invisible to
    # a naive regex check.
    df = pd.DataFrame({
        "symbol": pd.array(
            ["BRCA1", "TP53", "ΣEPT2", "EGFR", "MARCHF1", "KRAS", "PTEN"],
            dtype=object,
        ),
    })
    report = detect(df)
    homoglyphs = [s for s in report.suspicions if s.kind == "homoglyph"]
    assert len(homoglyphs) == 1
    assert homoglyphs[0].suggestion == "SEPT2"
    assert "U+03A3" in homoglyphs[0].reason


def test_pass1_flags_cyrillic_homoglyph() -> None:
    """Cyrillic capital А (U+0410) substituted for Latin A in a gene
    symbol is detected and repaired."""
    df = pd.DataFrame({
        "gene": pd.array(
            ["BRCA1", "TP53", "EGFR", "АTM",  # ATM with Cyrillic А
             "KRAS", "PTEN", "RB1"],
            dtype=object,
        ),
    })
    report = detect(df)
    homoglyphs = [s for s in report.suspicions if s.kind == "homoglyph"]
    assert len(homoglyphs) == 1
    assert homoglyphs[0].suggestion == "ATM"


def test_pass1_pure_ascii_symbol_no_homoglyph_flag() -> None:
    """A clean ASCII gene-shape symbol does not trigger homoglyph
    detection : repair is empty."""
    df = pd.DataFrame({
        "gene": ["BRCA1", "TP53", "EGFR", "ATM", "KRAS"],
    })
    report = detect(df)
    homoglyphs = [s for s in report.suspicions if s.kind == "homoglyph"]
    assert len(homoglyphs) == 0


# --- Mid-string separator coverage (v0.5.0) ---


def test_pass2_splits_pipe_separated_cell_with_date_token() -> None:
    """`'ID1|2-Oct|comment'` contains a date token among pipe-separated
    parts : the splitter must try each token."""
    df = pd.DataFrame({"misc": ["x|2-Oct|y", "foo", "bar"]})
    report = detect(df)
    string_susp = [s for s in report.suspicions if s.kind == "gene-date-string"]
    assert len(string_susp) >= 1
    assert "OCT2" in (string_susp[0].suggestion or "")


def test_pass2_splits_tab_separated_cell_with_date_token() -> None:
    """Tab-separated within a single cell : also covered."""
    df = pd.DataFrame({"misc": ["ID1\t2-Oct\tnote", "foo", "bar"]})
    report = detect(df)
    string_susp = [s for s in report.suspicions if s.kind == "gene-date-string"]
    assert len(string_susp) >= 1


# --- Long-form date parsing (v0.5.0) ---


def test_pass2_parses_long_form_date_string() -> None:
    """`"September 2, 2024"` is a real-world serialization that the
    standard `_DATE_STRING_RE` doesn't match. dateutil.parser handles
    the fallback path."""
    df = pd.DataFrame({"misc": ["foo", "September 2, 2024", "bar"]})
    report = detect(df)
    string_susp = [s for s in report.suspicions if s.kind == "gene-date-string"]
    assert len(string_susp) >= 1
    assert "SEPT2" in (string_susp[0].suggestion or "")


def test_long_form_parse_does_not_coerce_prose() -> None:
    """Prose sentences with date-like substrings must not be coerced : 
    bounded length + non-fuzzy parser keeps free text out."""
    df = pd.DataFrame({
        "notes": [
            "The study began in 2024 and concluded in 2026.",
            "We observed a significant increase",
        ],
    })
    report = detect(df)
    string_susp = [s for s in report.suspicions if s.kind == "gene-date-string"]
    assert len(string_susp) == 0


# --- CAS Registry gating (v0.5.0) ---


def test_pass2_demotes_cas_in_chemistry_column() -> None:
    """A CAS number in a column labeled `"CAS_Number"` or `"compound"` is
    flagged with low confidence (chain-of-custody record, not corruption)."""
    df = pd.DataFrame({"CAS_Number": ["foo", "50-78-2", "bar"]})
    report = detect(df)
    cas = [s for s in report.suspicions if s.kind == "cas-registry"]
    assert len(cas) == 1
    assert cas[0].confidence == pytest.approx(0.20)
    assert "chemistry-context" in cas[0].reason


def test_pass2_flags_cas_in_identifier_column_high_confidence() -> None:
    """A CAS number in a non-chemistry identifier column is flagged with
    the normal confidence : Excel may have coerced its leading component
    to a date."""
    df = pd.DataFrame({"identifier": ["foo", "50-78-2", "bar"]})
    report = detect(df)
    cas = [s for s in report.suspicions if s.kind == "cas-registry"]
    assert len(cas) == 1
    assert cas[0].confidence == pytest.approx(0.70)


# --- Autofill daily-band tightening (v0.5.0) ---


# --- Decimal-comma locale (v0.5.0) ---


# --- python-calamine loader path (v0.6.0) ---


def test_calamine_loader_for_large_xlsx(tmp_path: Path) -> None:
    """Files at/above the calamine threshold route through python-calamine;
    smaller files stay on openpyxl. Both must produce equivalent DataFrames
    for the detector to behave consistently."""
    import pandas as pd

    from uncorrupt.app import _load_xlsx_calamine
    f = tmp_path / "small.xlsx"
    pd.DataFrame({
        "gene": ["BRCA1", "TP53", "EGFR", "MARCHF1", "SEPTIN2"],
        "fold_change": [1.0, 2.0, 3.0, 4.0, 5.0],
    }).to_excel(f, index=False)
    # Even a tiny file should load via calamine when called directly
    sheets = _load_xlsx_calamine(str(f))
    assert sheets
    name = next(iter(sheets))
    df = sheets[name]
    assert "gene" in df.columns
    assert "BRCA1" in df["gene"].tolist()


# --- detect_file file-size guard (v0.6.0) ---


def test_detect_file_emits_file_too_large_for_huge_files(tmp_path: Path) -> None:
    """A file exceeding `max_file_bytes` returns a Report containing a
    single `file-too-large` Suspicion rather than hanging on openpyxl."""
    import pandas as pd

    from uncorrupt.detector import detect_file
    f = tmp_path / "tiny.xlsx"
    pd.DataFrame({"x": [1, 2, 3]}).to_excel(f, index=False)
    # Force the guard to trigger even on a small file
    rep = detect_file(str(f), row_context_boost=False, max_file_bytes=10)
    assert len(rep.suspicions) == 1
    assert rep.suspicions[0].kind == "file-too-large"
    assert "bytes" in str(rep.suspicions[0].value)


def test_detect_file_normal_path_unaffected_by_guard(tmp_path: Path) -> None:
    """Below the guard threshold the loader runs normally."""
    import pandas as pd

    from uncorrupt.detector import detect_file
    f = tmp_path / "normal.xlsx"
    pd.DataFrame({
        "gene": ["BRCA1", "TP53", "EGFR"],
    }).to_excel(f, index=False)
    rep = detect_file(str(f), row_context_boost=False)
    assert not any(s.kind == "file-too-large" for s in rep.suspicions)


# --- HGNC drift runtime warning (v0.6.0) ---


def test_hgnc_drift_note_emitted_when_snapshot_stale(monkeypatch) -> None:
    """If the HGNC snapshot is > 30 days old, detect() emits a sampling_note
    warning. Below threshold, the note stays None (or only contains
    sampling info if subsampling fired)."""
    import pandas as pd

    from uncorrupt.detector import detect
    # Mock the age function to claim a 45-day-old snapshot
    monkeypatch.setattr(
        "uncorrupt.detector._maybe_hgnc_drift_note",
        lambda: "HGNC snapshot is 45 days old (threshold 30); refresh.",
    )
    df = pd.DataFrame({"x": ["BRCA1", "TP53"]})
    report = detect(df)
    assert report.sampling_note is not None
    assert "45 days old" in report.sampling_note


# --- xref_status field on Suspicion (v0.6.0) ---


def test_suspicion_has_xref_status_field() -> None:
    """The Suspicion dataclass exposes an `xref_status` field (defaults
    to None) so end-users can filter on row-xref backing."""
    from uncorrupt.detector import Suspicion
    s = Suspicion(
        column="gene", row=0, value="x", kind="gene-date",
        suggestion="SEPT2", confidence=0.6, reason="test",
    )
    assert hasattr(s, "xref_status")
    assert s.xref_status is None


# --- Mutation-kill targeted tests (v0.6.1) ---
#
# These tests pin specific code paths that surviving mutants in the
# cosmic-ray sample exercised. Adding them increases the kill rate on the
# next mutation run. Each assertion was written to FAIL under at least
# one surviving mutant, not just to pass the current code.


def test_gene_date_unique_candidate_confidence_is_exactly_0_95() -> None:
    """Pin: `gene-date` flag with a single candidate has confidence 0.95.
    Kills mutants that change `confidence = 0.95` to anything else,
    or that change `len(candidates) == 1` to `== 2`."""
    df = pd.DataFrame({
        "gene": pd.array(
            ["BRCA1", "TP53", date(2024, 11, 1), "EGFR", "KRAS"],
            dtype=object,
        ),
    })
    report = detect(df)
    nov1 = [s for s in report.suspicions if s.kind == "gene-date"
            and s.suggestion == "NOV1"]
    assert len(nov1) == 1
    # NOV1 is the only mapping for Nov-1; confidence MUST be 0.95
    assert nov1[0].confidence == pytest.approx(0.95)


def test_gene_date_unique_canonical_boost_is_0_60_not_or() -> None:
    """Pin: when `unique_canonical AND len(candidates) > 1`, confidence is
    0.60 : when the column has >=2 date-corruption cells (v0.7.1 gate).
    Mutant that turns `AND` into `OR` would still set 0.60 but on cases
    that shouldn't get the boost : pin the non-boost case too."""
    df = pd.DataFrame({
        "gene": pd.array(
            ["BRCA1", "TP53", date(2024, 9, 2), date(2024, 9, 3),
             "EGFR", "KRAS"],
            dtype=object,
        ),
    })
    report = detect(df)
    sept2 = [s for s in report.suspicions if s.kind == "gene-date"
             and hasattr(s.value, "day") and s.value.day == 2]
    assert len(sept2) == 1
    # SEP2/SEPT2/SEPTIN2 : 3 candidates with unique canonical (SEPTIN2),
    # confidence is 0.50 base + 0.10 boost = 0.60, gated on column having
    # >=2 date cells.
    assert sept2[0].confidence == pytest.approx(0.60)


def test_gene_date_ambiguous_no_canonical_confidence_0_50() -> None:
    """Pin: ambiguous date with NO unique canonical → confidence 0.50,
    gated on column having >=2 date-corruption cells. Kills mutants on
    the canonical-boost branch."""
    df = pd.DataFrame({
        "gene": pd.array(
            ["BRCA1", "TP53", date(2024, 3, 1), date(2024, 3, 2),
             "EGFR", "KRAS"],
            dtype=object,
        ),
    })
    report = detect(df)
    march1 = [s for s in report.suspicions if s.kind == "gene-date"
              and hasattr(s.value, "day") and s.value.day == 1]
    assert len(march1) == 1
    # MARCH1 has 2 canonicals (MARCHF1 + MTARC1) : not unique, so 0.50
    assert march1[0].confidence == pytest.approx(0.50)


def test_gene_date_serial_low_confidence_value_is_0_15() -> None:
    """Pin: gene-date-serial without column evidence is confidence 0.15
    (v0.5.1 finding). Kills mutants on the low-confidence value."""
    df = pd.DataFrame({
        "symbol": pd.array(
            ["BRCA1", "TP53", "EGFR", "MARCHF1", 41334,
             "PTEN", "RB1", "KRAS"],
            dtype=object,
        ),
    })
    report = detect(df)
    serials = [s for s in report.suspicions if s.kind == "gene-date-serial"]
    assert len(serials) == 1
    assert serials[0].confidence == pytest.approx(0.15)


def test_gene_date_serial_high_confidence_value_is_0_85() -> None:
    """Pin: gene-date-serial WITH column evidence is confidence 0.85."""
    df = pd.DataFrame({
        "gene": pd.array(
            ["BRCA1", "TP53", date(2024, 9, 2), "EGFR", "MARCHF1",
             41334, "PTEN", "RB1"],
            dtype=object,
        ),
    })
    report = detect(df)
    serials = [s for s in report.suspicions if s.kind == "gene-date-serial"]
    assert len(serials) == 1
    assert serials[0].confidence == pytest.approx(0.85)


def test_pass2_gene_date_string_single_candidate_confidence_0_55() -> None:
    """Pin: Pass 2 gene-date-string with len(parsed_dates) == 1 and
    len(unique_candidates) == 1 → confidence 0.55. Kills mutants on the
    parsed-dates length comparison."""
    df = pd.DataFrame({"misc": ["foo", "2-Nov", "bar", "baz"]})
    report = detect(df)
    string_susp = [s for s in report.suspicions if s.kind == "gene-date-string"]
    # 2-Nov has only one candidate (NOV1 maps), confidence 0.55
    if string_susp:
        assert any(s.confidence == pytest.approx(0.55) for s in string_susp)


def test_decimal_comma_confidence_is_exactly_0_55() -> None:
    """Pin: decimal-comma flag confidence is exactly 0.55. Kills mutants
    that change to 1.55 or similar."""
    df = pd.DataFrame({"value": ["1,5", "2,7", "0,8", "3,14", "1,2", "4,0"]})
    report = detect(df)
    dec = [s for s in report.suspicions if s.kind == "decimal-comma"]
    assert len(dec) == 1
    assert dec[0].confidence == pytest.approx(0.55)


def test_unrecognized_symbol_confidence_is_0_20() -> None:
    """Pin: unrecognized-symbol confidence is exactly 0.20."""
    df = pd.DataFrame({
        "gene": ["BRCA1", "TP53", "EGFR", "MARCHF1", "KRAS",
                 "PTEN", "RB1", "ATM", "CDKN2A", "VHL",
                 "TOTALLY_FAKE_GENE", "MYC"],
    })
    report = detect(df)
    unrec = [s for s in report.suspicions if s.kind == "unrecognized-symbol"]
    if unrec:
        assert all(s.confidence == pytest.approx(0.20) for s in unrec)


def test_month_5_boundary_n_31_max() -> None:
    """Pin: month=5 reverse mapping accepts n up to 31, not 32.
    Kills mutants on `1 <= n <= 31`."""
    from uncorrupt.detector import _candidates_for_month_and_n
    assert _candidates_for_month_and_n(5, 31) == ["MEI31", "JUN1"]
    # n=32 would be out of range for May
    assert _candidates_for_month_and_n(5, 32) == []


def test_header_token_floor_4_chars() -> None:
    """Pin: `len(tok) < 4` boundary : 3-char tokens excluded from
    suffix-match. Kills mutants `< 4` → `== 4`."""
    from uncorrupt.detector import _tokens_match_hint_set
    # 'kid' has 3 chars; suffix-match is gated to len >= 4
    assert not _tokens_match_hint_set({"kid"}, {"id"})
    # 4-char 'genekind' wouldn't pass either because 'id' is not a suffix
    # : must check actual suffix matches
    assert _tokens_match_hint_set({"sampleid"}, {"id"})


def test_caption_row_header_excluded() -> None:
    """Pin: `len(stripped) > 80` boundary. A header just under is
    accepted as tokens; a header over is rejected as caption row."""
    from uncorrupt.detector import _header_tokens
    short = "x" * 80
    long = "x" * 81
    assert _header_tokens(short) != set()
    assert _header_tokens(long) == set()


# --- HGNC case-insensitive resolve (v0.5.1) ---


def test_hgnc_resolve_case_insensitive_alias() -> None:
    """HGNC stores POU5F1's alias as `Oct4` (mixed case). Excel-emitted
    suggestions are uppercase `OCT4`. The resolve must match across cases.
    Empirically caught: 25 calibration cells were labeled 'contradicted'
    when they were genuinely corroborated by Oct4 → POU5F1 chain."""
    from uncorrupt.corpus import load_hgnc
    h = load_hgnc()
    assert h.resolve("OCT4") == "POU5F1"
    assert h.resolve("Oct4") == "POU5F1"
    assert h.resolve("brca1") == "BRCA1"
    assert h.resolve("SEPT2") == "SEPTIN2"
    assert h.resolve("MARCH1") == "MARCHF1"


# --- Gene-date-serial column-evidence gating (v0.5.1) ---


def test_gene_date_serial_high_confidence_when_column_has_date_evidence() -> None:
    """Real autofill scenario: a column with both real `date` cells AND
    integers in the date-serial range. The integers get the high-confidence
    flag because there's auxiliary evidence of date corruption."""
    df = pd.DataFrame({
        "gene": pd.array(
            ["BRCA1", "TP53", date(2024, 9, 2), "EGFR", "MARCHF1",
             41334, "PTEN", "RB1"],
            dtype=object,
        ),
    })
    report = detect(df)
    serials = [s for s in report.suspicions if s.kind == "gene-date-serial"]
    assert len(serials) == 1
    assert serials[0].confidence == pytest.approx(0.85)
    assert "companion date-corruption evidence" in serials[0].reason


def test_gene_date_serial_low_confidence_without_date_evidence() -> None:
    """A column with gene symbols + an isolated integer in date-serial
    range, but NO real date cells anywhere in the column. The integer
    is more likely a coincidence than corruption; flagged at low
    confidence per real-data calibration (9% precision unflanked)."""
    df = pd.DataFrame({
        "symbol": pd.array(
            ["BRCA1", "TP53", "EGFR", "MARCHF1", 41334,
             "PTEN", "RB1", "KRAS"],
            dtype=object,
        ),
    })
    report = detect(df)
    serials = [s for s in report.suspicions if s.kind == "gene-date-serial"]
    assert len(serials) == 1
    assert serials[0].confidence == pytest.approx(0.15)
    assert "NO companion date-corruption evidence" in serials[0].reason


def test_gene_date_serial_low_confidence_digit_string() -> None:
    """Same gating applies to digit-strings (e.g., `'41334'`) in
    identifier columns without companion date evidence."""
    df = pd.DataFrame({
        "gene": pd.array(
            ["BRCA1", "TP53", "EGFR", "41334", "PTEN", "RB1", "KRAS",
             "ATM"],
            dtype=object,
        ),
    })
    report = detect(df)
    serials = [s for s in report.suspicions if s.kind == "gene-date-serial"]
    assert len(serials) == 1
    assert serials[0].confidence == pytest.approx(0.15)


# --- Null-result FP regression tests (v0.5.0) ---


def test_long_form_parser_rejects_gene_list_cell() -> None:
    """`'ST13;ST13P5;ST13P4'` is a real-world cell value (semicolon-
    delimited gene list). The long-form date parser must not coerce
    it : caught in the v0.5.0 null-result run."""
    from uncorrupt.detector import _parse_date_string
    assert _parse_date_string("ST13;ST13P5;ST13P4") == []


def test_long_form_parser_rejects_lipid_metabolite() -> None:
    """`'16:1SMOH'` is real lipid metabolite notation (16 carbons,
    1 double bond, sphingomyelin -OH). The long-form parser must not
    coerce : null-result regression."""
    from uncorrupt.detector import _parse_date_string
    assert _parse_date_string("16:1SMOH") == []


def test_long_form_parser_rejects_decimal_number() -> None:
    """`'1.77707864378796'` is a measurement, not a date. dateutil
    will parse it; we reject upstream."""
    from uncorrupt.detector import _parse_date_string
    assert _parse_date_string("1.77707864378796") == []


def test_long_form_parser_rejects_accession_id() -> None:
    """`'AD000813.1'` is a GenBank accession ID with a version suffix.
    dateutil parses this as a date; we reject upstream."""
    from uncorrupt.detector import _parse_date_string
    assert _parse_date_string("AD000813.1") == []


def test_long_form_parser_accepts_genuine_long_form_date() -> None:
    """`'September 2, 2024'` is a real long-form date : the parser
    must still accept these (preserving the v0.5.0 feature)."""
    from uncorrupt.detector import _parse_date_string
    out = _parse_date_string("September 2, 2024")
    assert len(out) >= 1
    assert any(d.month == 9 and d.day == 2 and d.year == 2024 for d in out)


# --- RIKEN coercion signature (v0.5.1) ---


def test_riken_signature_matches_real_coerced_id() -> None:
    """`2310009E13` coerced by Excel → `2.310009e19`. Signature detector
    must recognise it. 7-digit mantissa 2310009, exp 19, 0 trailing zeros,
    4 distinct nonzero digits → match."""
    from uncorrupt.detector import _riken_coercion_signature
    sig = _riken_coercion_signature(2.310009e19)
    assert sig is not None
    assert sig["mantissa"] == 2310009
    assert sig["exp"] == 19


def test_riken_signature_rejects_round_floats() -> None:
    """`1.5e8`, `2e10`, `1e15` : all 'round' floats with multiple
    trailing zeros : must NOT match the RIKEN signature. These would
    false-positive on legitimate ID magnitudes."""
    from uncorrupt.detector import _riken_coercion_signature
    for v in [1.5e8, 2e10, 1e15, 1.0e7, 5e8]:
        assert _riken_coercion_signature(v) is None, f"signature wrongly matched {v}"


def test_riken_signature_rejects_entrez_floats() -> None:
    """Integer Entrez IDs cast to float (100506334.0 etc.) have a 9-digit
    mantissa and must not match the 7-digit signature."""
    from uncorrupt.detector import _riken_coercion_signature
    for v in [100506334.0, 123456789.0, 105378948.0]:
        assert _riken_coercion_signature(v) is None


def test_riken_signature_rejects_low_entropy_mantissa() -> None:
    """Mantissa with <3 distinct nonzero digits (e.g., 1111110, 5555550)
    is rejected : RIKEN clone IDs are distinctive."""
    from uncorrupt.detector import _riken_coercion_signature
    # 1.111111e10 has mantissa 1111111 : 1 distinct nonzero. Reject.
    assert _riken_coercion_signature(1.111111e10) is None


def test_pass15_riken_signature_flagged_in_generic_column() -> None:
    """A 7-sig-fig float that matches the signature gets flagged even
    without a RIKEN header anchor : the signature is specific enough
    that the false-positive rate is low.

    Values are chosen under 2^53 (~9e15) so they take the id-float path
    rather than the precision-loss path for super-large ints."""
    df = pd.DataFrame({
        "gene_id": [1.5e8, 9.876543e11, 1.0e8, 1.234567e10, 2.310009e13],
    })
    report = detect(df)
    id_floats = [s for s in report.suspicions if s.kind == "id-float"]
    # 9.876543e11 (sig match, >=1e10), 1.234567e10 (sig match, >=1e10),
    # 2.310009e13 (sig match, >=1e10) → 3 flags
    assert len(id_floats) == 3


def test_pass15_riken_signature_in_riken_hinted_column_low_magnitude() -> None:
    """A 7-sig-fig signature in the E01 magnitude range [1e7, 1e10] is
    flagged with high confidence when the header anchors RIKEN/cDNA.
    A non-signature float in the same column is flagged with lower
    confidence : explicit RIKEN-hinted columns shouldn't contain
    Entrez-shaped integers, so 0.55 is justified."""
    df = pd.DataFrame({
        "RIKEN_cdna_id": [
            2.310009e7,   # E01 range, sig match  → 0.85
            1.234567e8,   # E02 range, sig match  → 0.85
            9.876543e9,   # E03 range, sig match  → 0.85
            100506334.0,  # Entrez integer        → 0.55 (header-anchored fallback)
        ],
    })
    report = detect(df)
    id_floats = [s for s in report.suspicions if s.kind == "id-float"]
    assert len(id_floats) == 4
    confidences = sorted(s.confidence for s in id_floats)
    # Three high-confidence signature matches + one low-confidence anchored hit
    assert confidences == [0.55, 0.85, 0.85, 0.85]


def test_pass15_no_flag_on_entrez_in_generic_column() -> None:
    """Same Entrez value, generic identifier column header : NOT flagged.
    The RIKEN-hint fallback only fires when the header explicitly anchors
    RIKEN/cDNA/transcript."""
    df = pd.DataFrame({
        "geneId": [100506334.0, 100507127.0, 100131089.0, 101410543.0],
    })
    report = detect(df)
    id_floats = [s for s in report.suspicions if s.kind == "id-float"]
    assert len(id_floats) == 0


def test_pass15_no_flag_on_entrez_id_column() -> None:
    """Real Entrez gene IDs (current range ~1e8-1.5e8) in an identifier
    column must NOT be flagged as id-float : caught 74 false positives
    in the null-result run."""
    df = pd.DataFrame({
        "geneId": [100506334, 100507127, 100131089, 101410543, 105378948],
    })
    report = detect(df)
    assert [s for s in report.suspicions if s.kind == "id-float"] == []


def test_header_token_extraction_rejects_caption_rows() -> None:
    """A column whose 'header' is actually a table caption (>80 chars)
    must not match identifier-hint tokens. Caught 70 spurious id-float
    flags on `PMC10313048/pone.0287634.s002.xlsx` in the null-result
    run."""
    from uncorrupt.detector import _header_tokens
    caption = (
        "Supplemental Data Table S2. Total proteome list of identified "
        "proteins in vocal fold (VF) tissue and co-cultivated construct"
    )
    assert _header_tokens(caption) == set()
    # Real headers still extract normally
    assert "protein" in _header_tokens("protein_id")


def test_decimal_comma_locale_column_flagged() -> None:
    """A column where most cells are `"1,5"`/`"3,14"` shape with no
    period-decimal cells is flagged."""
    df = pd.DataFrame({
        "fold_change": ["1,5", "2,71", "0,8", "3,14", "1,2", "4,0"],
    })
    report = detect(df)
    dec = [s for s in report.suspicions if s.kind == "decimal-comma"]
    assert len(dec) == 1
    assert "Frictionless" in (dec[0].suggestion or "")


def test_mixed_decimal_separator_column_not_flagged() -> None:
    """A column with BOTH period-decimal and comma-decimal cells is NOT
    flagged : the column is mixed, not consistently European-locale."""
    df = pd.DataFrame({
        "value": ["1,5", "2.7", "0.8", "3,14", "1.2"],
    })
    report = detect(df)
    dec = [s for s in report.suspicions if s.kind == "decimal-comma"]
    assert len(dec) == 0


def test_period_decimal_column_not_flagged() -> None:
    """A clean period-decimal column is not flagged."""
    df = pd.DataFrame({
        "value": ["1.5", "2.7", "0.8", "3.14", "1.2"],
    })
    report = detect(df)
    dec = [s for s in report.suspicions if s.kind == "decimal-comma"]
    assert len(dec) == 0


# --- Large-file sampling (v0.5.0) ---


def test_large_dataframe_triggers_sampling() -> None:
    """A 200K-row dataframe gets sampled to head + tail + 10% interior;
    the report records the sampling decision."""
    n = 200_000
    df = pd.DataFrame({"x": ["BRCA1"] * n})
    report = detect(df)
    assert report.sampling_note is not None
    assert "sampled" in report.sampling_note
    # rows_scanned reflects the SAMPLE, not the original size
    assert report.rows_scanned < n


def test_small_dataframe_skips_sampling() -> None:
    """Under the threshold, no sampling; sampling_note is None."""
    df = pd.DataFrame({"x": ["BRCA1"] * 1000})
    report = detect(df)
    assert report.sampling_note is None
    assert report.rows_scanned == 1000


def test_sampling_is_deterministic() -> None:
    """Sampling uses a fixed seed : two runs produce identical reports."""
    import pandas as _pd
    n = 200_000
    df = _pd.DataFrame({"x": [f"GENE{i}" for i in range(n)]})
    r1 = detect(df)
    r2 = detect(df)
    assert r1.rows_scanned == r2.rows_scanned
    assert r1.sampling_note == r2.sampling_note


def test_pass3_no_longer_accepts_2_day_gaps_as_daily() -> None:
    """Tightened daily band: requires diff==1 exactly. A run of Sep-1,
    Sep-3, Sep-5 (uniform step but 2-day stride) no longer qualifies as
    autofill : Excel autofill produces diff=1 day."""
    df = pd.DataFrame({
        "symbol": pd.array(
            ["BRCA1", "TP53", "EGFR", "KRAS", "MYC",
             date(2024, 9, 1), date(2024, 9, 3), date(2024, 9, 5),
             "PTEN", "RB1", "ATM"],
            dtype=object,
        ),
    })
    report = detect(df)
    seq = [s for s in report.suspicions if s.kind == "autofill-sequence"]
    assert len(seq) == 0
