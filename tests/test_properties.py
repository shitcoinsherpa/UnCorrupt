"""Property-based tests using `hypothesis`.

These complement the example-based unit tests by stress-testing invariants
across the full input space. Each property is a claim about the detector
that should hold for ANY input matching the given strategy.

Run targeted: `pytest tests/test_properties.py -q`

Per Mark Ziemann's "tests grow with the code, not after it" principle, every
detector decision rule deserves at least one property a future reader could
verify is invariant.
"""
from __future__ import annotations

from datetime import date

from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st

from uncorrupt.detector import (
    _candidates_for_month_and_n,
    _canonicalize_candidates,
    _parse_date_string,
    _reverse_gene_date,
    _serial_to_date,
)

# --- Strategies --------------------------------------------------------------

# Excel-corruption-relevant month + day-of-month combinations. These are the
# months Excel auto-converts gene symbols into.
# Includes August (AGO family : locale-specific Italian/Spanish/Portuguese
# corruption per Ziemann 2021), May (MEI/MAY/JUN families), and January
# (TAMM Finnish-locale corruption).
_CORRUPTION_MONTHS = st.sampled_from([1, 2, 3, 4, 5, 8, 9, 10, 11, 12])

# Day-of-month bounded to plausible gene-suffix range (15 covers SEPT15 : the
# highest-numbered Excel-corruption gene). Hypothesis will also explore the
# boundary values 1 and 31.
_GENE_SUFFIX_RANGE = st.integers(min_value=1, max_value=31)


# --- _candidates_for_month_and_n invariants ---------------------------------


@given(month=_CORRUPTION_MONTHS, n=_GENE_SUFFIX_RANGE)
def test_candidates_for_month_returns_list_of_strings(month: int, n: int) -> None:
    out = _candidates_for_month_and_n(month, n)
    assert isinstance(out, list)
    for c in out:
        assert isinstance(c, str)
        assert len(c) > 0


@given(month=_CORRUPTION_MONTHS, n=_GENE_SUFFIX_RANGE)
def test_candidates_for_month_uppercase_only(month: int, n: int) -> None:
    """Every candidate is uppercase A-Z + digits + optional hyphen : matches
    the HGNC current/prev-symbol convention."""
    for c in _candidates_for_month_and_n(month, n):
        assert c.upper() == c, f"{c!r} contains lowercase characters"


@given(month=_CORRUPTION_MONTHS, n=_GENE_SUFFIX_RANGE)
def test_candidates_match_n_in_suffix(month: int, n: int) -> None:
    """Each candidate's numeric suffix equals `n` : proving the reversal
    is faithful and doesn't off-by-one.

    EXCEPTION: the JUN family is *intentionally* not a faithful reversal : 
    Ziemann 2021 documents that "jun-1" → "May-31" via numeric subtraction
    (June minus 1 = May 31), so n=31 produces JUN1, n=30 produces JUN2,
    n=29 produces JUN3. This is a documented Excel quirk, not an off-by-one.
    """
    import re
    for c in _candidates_for_month_and_n(month, n):
        if c.startswith("JUN"):
            # JUN family is a documented exception : see docstring.
            continue
        m = re.search(r"(\d+)$", c)
        assert m is not None, f"candidate {c!r} has no trailing digits"
        assert int(m.group(1)) == n, (
            f"candidate {c!r} suffix != input n={n}"
        )


# --- _canonicalize_candidates invariants -----------------------------------


@given(candidates=st.lists(st.sampled_from([
    "SEPT2", "SEP2", "SEPTIN2", "MARCH1", "MARC1", "MARCHF1", "MTARC1",
    "DEC1", "DELEC1", "OCT4", "NOV1", "APR1", "FEB3", "BRCA1", "TP53",
    "EGFR", "KRAS",
]), max_size=8))
def test_canonicalize_never_drops_inputs(candidates: list[str]) -> None:
    """Per the 'no silent repair' rule, every input candidate must survive
    canonicalization (possibly reordered, never dropped)."""
    out, _ = _canonicalize_candidates(candidates)
    for c in candidates:
        assert c in out, f"canonicalize dropped {c!r}"


@given(candidates=st.lists(st.sampled_from([
    "SEPT2", "SEP2", "MARCH1", "MARC1", "DEC1", "BRCA1",
]), max_size=8))
def test_canonicalize_no_duplicates(candidates: list[str]) -> None:
    out, _ = _canonicalize_candidates(candidates)
    assert len(out) == len(set(out)), f"duplicates in {out}"


@given(candidates=st.lists(st.sampled_from([
    "SEPT2", "SEP2", "MARCH1", "MARC1", "DEC1", "BRCA1",
]), max_size=8))
def test_canonicalize_idempotent(candidates: list[str]) -> None:
    """canonicalize(canonicalize(x)) == canonicalize(x). Future re-application
    must be a no-op."""
    out1, unique1 = _canonicalize_candidates(candidates)
    out2, unique2 = _canonicalize_candidates(out1)
    assert out2 == out1
    # unique-canonical bit may flip from True to False if applying twice
    # because modern symbols already lead : but should never flip from False
    # to True
    if unique1:
        assert unique2 or "SEPTIN" in (out1[0] if out1 else "") or \
               "MARCHF" in (out1[0] if out1 else "") or \
               "DELEC" in (out1[0] if out1 else "") or \
               "MTARC" in (out1[0] if out1 else "")


# --- _serial_to_date round-trip ---------------------------------------------


@given(n=st.integers(min_value=20000, max_value=60000))
def test_serial_to_date_round_trip(n: int) -> None:
    """Every serial in our valid range maps to a date; converting back via
    the Excel-epoch math yields the same integer."""
    d = _serial_to_date(n)
    assert d is not None
    # Round-trip
    from datetime import date as _date
    epoch = _date(1899, 12, 30)
    recovered = (d - epoch).days
    assert recovered == n


@given(n=st.integers(min_value=-10_000, max_value=10_000))
def test_serial_to_date_out_of_range_returns_none(n: int) -> None:
    """Values outside [20000, 60000] never decode."""
    assume(n < 20000 or n > 60000)
    assert _serial_to_date(n) is None


# --- _reverse_gene_date symmetry -------------------------------------------


@settings(suppress_health_check=[HealthCheck.too_slow], deadline=None)
@given(
    year=st.integers(min_value=1990, max_value=2030),
    month=_CORRUPTION_MONTHS,
    day=st.integers(min_value=1, max_value=28),  # 28 covers every month safely
)
def test_reverse_gene_date_returns_only_uppercase_candidates(
    year: int, month: int, day: int,
) -> None:
    """For any plausible date, every candidate emitted is uppercase shape."""
    try:
        d = date(year, month, day)
    except ValueError:
        assume(False)
    for c in _reverse_gene_date(d):
        assert c.upper() == c
        assert c[0].isalpha()


@given(
    year=st.integers(min_value=2001, max_value=2015),
    day=st.integers(min_value=1, max_value=1),  # only day=1 triggers year-suffix
)
def test_year_suffix_mode_when_format_is_mmm_yy(year: int, day: int) -> None:
    """With format `mmm-yy` AND day=1, the year-suffix interpretation fires.
    The emitted candidate's numeric suffix must equal year mod 100."""
    d = date(year, 9, day)
    candidates = _reverse_gene_date(d, number_format="mmm-yy")
    yy = year % 100
    if 1 <= yy <= 15:
        # Must include SEPT{yy} or SEP{yy}
        assert any(c.endswith(str(yy)) for c in candidates), (
            f"year={year} -> candidates {candidates} contain no suffix {yy}"
        )


@given(
    year=st.integers(min_value=2020, max_value=2024),
    day=st.integers(min_value=1, max_value=11),
)
def test_day_of_month_mode_when_format_is_d_mmm(year: int, day: int) -> None:
    """With format `d-mmm` (no year), day-of-month interpretation wins."""
    d = date(year, 3, day)  # March
    candidates = _reverse_gene_date(d, number_format="d-mmm")
    # Must include MARCH{day}
    assert any(c == f"MARCH{day}" for c in candidates), (
        f"day={day} -> candidates {candidates} missing MARCH{day}"
    )


# --- _parse_date_string round-trip -----------------------------------------


@given(
    day=st.integers(min_value=1, max_value=28),
    month_abbrev=st.sampled_from(["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]),
)
def test_parse_date_string_handles_dd_mmm(day: int, month_abbrev: str) -> None:
    """`12-Sep`-style strings parse to a unique date."""
    s = f"{day}-{month_abbrev}"
    out = _parse_date_string(s)
    assert len(out) == 1
    assert out[0].day == day
    month_num = {"Jan": 1, "Feb": 2, "Mar": 3, "Apr": 4, "May": 5, "Jun": 6,
                  "Jul": 7, "Aug": 8, "Sep": 9, "Oct": 10, "Nov": 11, "Dec": 12,
                  }[month_abbrev]
    assert out[0].month == month_num


@given(
    year=st.integers(min_value=2000, max_value=2030),
    month=st.integers(min_value=1, max_value=12),
    day=st.integers(min_value=1, max_value=28),
)
def test_parse_iso_date_string_round_trip(year: int, month: int, day: int) -> None:
    """ISO `YYYY-MM-DD` strings parse to the exact corresponding date."""
    s = f"{year}-{month:02d}-{day:02d}"
    out = _parse_date_string(s)
    assert len(out) == 1
    assert out[0] == date(year, month, day)


@given(
    year=st.integers(min_value=2000, max_value=2030),
    month=st.integers(min_value=1, max_value=12),
    day=st.integers(min_value=1, max_value=28),
)
def test_parse_iso_datetime_string_round_trip(year: int, month: int, day: int) -> None:
    """ISO `YYYY-MM-DD HH:MM:SS` strings also parse correctly."""
    s = f"{year}-{month:02d}-{day:02d} 00:00:00"
    out = _parse_date_string(s)
    assert len(out) == 1
    assert out[0] == date(year, month, day)
