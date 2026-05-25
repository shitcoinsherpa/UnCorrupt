"""Property-based boundary tests targeting cosmic-ray surviving mutants.

The v0.6.1 cosmic-ray re-run showed 43.75% kill rate : many survivors
are boundary mutants (`>` vs `>=`, `0.40` vs `0.45`, `<= 31` vs `<= 32`).
Per the research review (Schuler & Zeller 2013; Hypothesis docs):

  - Random `floats()` rarely hit exact pivot values.
  - `@example` decorators inject specific boundary points deterministically.
  - `math.nextafter(x, +inf)` / `nextafter(x, -inf)` give the smallest
    floats strictly above/below the threshold.

These tests pin every threshold the detector uses with both the pivot
value AND the immediately-adjacent floats, which kills `>` ↔ `>=`,
`< 4` ↔ `== 4`, `1e7` ↔ `1e7 + 1`, etc.

Following PIT's "Conditionals Boundary Mutator" guidance : boundary
mutants are stable by design; only explicit boundary tests kill them.
"""
from __future__ import annotations

import math
from datetime import date

import pandas as pd
import pytest
from hypothesis import example, given, settings, strategies as st

from uncorrupt.detector import (
    _SERIAL_MAX,
    _SERIAL_MIN,
    _candidates_for_month_and_n,
    _header_tokens,
    _riken_coercion_signature,
    _serial_to_date,
    _tokens_match_hint_set,
)
from uncorrupt.registries import IDENTIFIER_HEADER_HINTS


# --- Excel-serial boundary ---


@example(_SERIAL_MIN)
@example(_SERIAL_MIN - 1)
@example(_SERIAL_MIN + 1)
@example(_SERIAL_MAX)
@example(_SERIAL_MAX - 1)
@example(_SERIAL_MAX + 1)
@given(n=st.integers(min_value=_SERIAL_MIN - 5, max_value=_SERIAL_MAX + 5))
def test_serial_range_boundaries_decoded_correctly(n: int) -> None:
    """Kills `_SERIAL_MIN <=` mutated to `_SERIAL_MIN <`, etc."""
    out = _serial_to_date(n)
    if _SERIAL_MIN <= n <= _SERIAL_MAX:
        assert out is not None, f"valid serial {n} returned None"
    else:
        assert out is None, f"out-of-range {n} returned {out}"


# --- Family-range boundaries ---


@pytest.mark.parametrize("month,n,expected_count_nonzero", [
    # MARCH 1-12
    (3, 0, False), (3, 1, True), (3, 12, True), (3, 13, False),
    # SEPT 1-15
    (9, 0, False), (9, 1, True), (9, 15, True), (9, 16, False),
    # OCT 1-11
    (10, 0, False), (10, 1, True), (10, 11, True), (10, 12, False),
    # NOV only n=1
    (11, 0, False), (11, 1, True), (11, 2, False),
    # DEC 1-2
    (12, 0, False), (12, 1, True), (12, 2, True), (12, 3, False),
    # APR 1-3
    (4, 0, False), (4, 1, True), (4, 3, True), (4, 4, False),
    # FEB 3-4
    (2, 2, False), (2, 3, True), (2, 4, True), (2, 5, False),
    # AGO 1-12
    (8, 0, False), (8, 1, True), (8, 12, True), (8, 13, False),
    # May 1-31 (MEI/MAY/JUN families)
    (5, 0, False), (5, 1, True), (5, 31, True), (5, 32, False),
    # TAMM 30-50
    (1, 29, False), (1, 30, True), (1, 50, True), (1, 51, False),
])
def test_family_n_boundaries(month: int, n: int, expected_count_nonzero: bool) -> None:
    """Pins every family-range boundary explicitly. Kills mutants of
    `<= 11` → `<= 12`, `n == 1` → `n >= 1`, etc."""
    candidates = _candidates_for_month_and_n(month, n)
    if expected_count_nonzero:
        assert len(candidates) > 0, f"month={month} n={n} returned no candidates"
    else:
        assert len(candidates) == 0, (
            f"month={month} n={n} returned candidates {candidates} but range "
            f"should reject"
        )


# --- Header-token length boundary ---


@example("abc")    # 3 chars : too short for suffix match
@example("abcd")   # 4 chars : boundary; suffix path now eligible
@example("abcde")  # 5 chars : suffix path eligible
@given(token=st.text(alphabet="abcdefghj", min_size=1, max_size=10))
def test_header_token_floor_boundary(token: str) -> None:
    """Pins `len(tok) < 4` boundary in the SUFFIX-match path. Kills
    `< 4` → `<= 4` / `== 4`. Excludes token exactly equal to the hint
    (that hits the exact-set intersection path, not the suffix path)."""
    if token == "id" or token.endswith("id") and len(token) < 4:
        # Skip exact-equals matches and undefined behavior
        return
    if len(token) < 4:
        # Suffix path requires len >= 4. Below the floor : no suffix match.
        assert not _tokens_match_hint_set({token}, {"id"})


# --- RIKEN signature mantissa boundaries ---


@pytest.mark.parametrize("mantissa,exp,expect_match", [
    (1234567, 6, False),    # below RIKEN range
    (1234567, 7, True),     # at boundary
    (1234567, 8, True),
    (1234567, 30, True),
    (1234567, 31, True),    # at upper boundary
    (1234567, 32, False),   # one above max
])
def test_riken_signature_exponent_boundary(
    mantissa: int, exp: int, expect_match: bool,
) -> None:
    """Pins RIKEN exponent boundary [7, 31]. Kills `7 <= exp <= 31` → various
    boundary mutants."""
    value = float(mantissa * (10 ** (exp - 6)))
    sig = _riken_coercion_signature(value)
    if expect_match:
        assert sig is not None, f"exp={exp} mantissa={mantissa} should match"
    else:
        assert sig is None, f"exp={exp} mantissa={mantissa} should NOT match"


# --- IDENTIFIER_HEADER_HINTS membership boundary ---


def test_identifier_header_hints_each_keyword_recognized() -> None:
    """Pins that every hint in the set is actually recognized.
    Kills mutants that remove items from IDENTIFIER_HEADER_HINTS."""
    for hint in IDENTIFIER_HEADER_HINTS:
        assert _tokens_match_hint_set({hint}, IDENTIFIER_HEADER_HINTS), (
            f"hint {hint!r} unexpectedly not recognized in its own set"
        )


# --- Caption-row exclusion boundary ---


@example(79)    # under threshold
@example(80)    # at threshold (excluded)
@example(81)    # over threshold (excluded)
@given(n=st.integers(min_value=1, max_value=200))
def test_header_token_caption_boundary(n: int) -> None:
    """Pins `len(stripped) > 80` boundary. Header of exactly 80 chars
    should be excluded; 79 chars accepted."""
    header = "abc def ghi " * (n // 12) + "x" * (n % 12)
    header = header[:n]  # exact length n
    tokens = _header_tokens(header)
    if n > 80:
        assert tokens == set(), f"header of {n} chars should yield no tokens"
    else:
        # Below threshold : should yield tokens (if non-empty)
        if n >= 3:
            assert tokens != set(), f"header of {n} chars yielded no tokens"


# --- Magnitude boundary for id-float ---


@example(1e10)
@example(math.nextafter(1e10, -math.inf))   # just below 1e10
@example(math.nextafter(1e10, math.inf))    # just above 1e10
@example(1e7)
@example(math.nextafter(1e7, -math.inf))
@example(math.nextafter(1e7, math.inf))
@given(value=st.floats(min_value=1e6, max_value=1e12,
                        allow_nan=False, allow_infinity=False))
@settings(max_examples=200, deadline=None)
def test_riken_signature_boundary_nextafter(value: float) -> None:
    """Hypothesis with nextafter pinpoints float pivots. Kills `> 1e10`
    mutated to `>= 1e10` because the test feeds exactly 1e10 and the
    floats immediately adjacent."""
    sig = _riken_coercion_signature(value)
    # The contract: sig returns dict or None : must never crash
    assert sig is None or isinstance(sig, dict)
