"""Deterministic simulator of Microsoft Excel's import-time auto-conversion.

Empirically confirmed on 2026-05-14 against LibreOffice 7.3.7.2 headless:
- "2310009E13" -> float 2.310009e+19 (matches Excel behavior; LibreOffice agrees)
- "SEPT2" / "MARCH1" / "DEC1" -> string in LibreOffice (LibreOffice does NOT
  corrupt these); Excel DOES corrupt to dates per Ziemann 2016. The simulator
  encodes Excel's behavior, not LibreOffice's.
"""
from __future__ import annotations

from datetime import date
from typing import Any

from .registries import EXPONENT_PATTERN

_REFERENCE_YEAR = 2024

_GENE_TO_DATE: dict[str, date] = {}
# English-locale conversions : Excel's default behavior with US/UK regional
# settings. MARCH extended to 12 per the HGNC MARCHF1-12 alias range.
for _i in range(1, 13):
    _GENE_TO_DATE[f"MARCH{_i}"] = date(_REFERENCE_YEAR, 3, _i)
for _i in range(1, 16):
    _GENE_TO_DATE[f"SEPT{_i}"] = date(_REFERENCE_YEAR, 9, _i)
    _GENE_TO_DATE[f"SEP{_i}"] = date(_REFERENCE_YEAR, 9, _i)
for _i in range(1, 3):
    _GENE_TO_DATE[f"MARC{_i}"] = date(_REFERENCE_YEAR, 3, _i)
for _i in range(1, 12):
    _GENE_TO_DATE[f"OCT{_i}"] = date(_REFERENCE_YEAR, 10, _i)
for _i in range(1, 4):
    _GENE_TO_DATE[f"APR{_i}"] = date(_REFERENCE_YEAR, 4, _i)
for _i in (3, 4):
    # Only FEB3 and FEB4 are real HGNC gene aliases; FEB1, FEB2, FEB5+
    # have no HGNC mapping (so the detector intentionally doesn't reverse
    # those, and the simulator must not invent corruption for non-genes).
    _GENE_TO_DATE[f"FEB{_i}"] = date(_REFERENCE_YEAR, 2, _i)
_GENE_TO_DATE["NOV1"] = date(_REFERENCE_YEAR, 11, 1)
_GENE_TO_DATE["DEC1"] = date(_REFERENCE_YEAR, 12, 1)
_GENE_TO_DATE["DEC2"] = date(_REFERENCE_YEAR, 12, 2)

# Locale-specific corruptions documented in Ziemann 2021. These do NOT fire
# in the default en-US Excel; they require regional settings that interpret
# the abbreviation as a month name. AGO = "agosto" (August) in
# Italian/Spanish/Portuguese; MEI = May in Dutch; TAMM = "tammikuu"
# (January) in Finnish, where the day-of-month is encoded in the year
# component (Tamm-30 → 1930-01-01 by Excel's two-digit-year rules).
for _i in range(1, 13):
    _GENE_TO_DATE[f"AGO{_i}"] = date(_REFERENCE_YEAR, 8, _i)
for _i in range(1, 32):
    _GENE_TO_DATE[f"MEI{_i}"] = date(_REFERENCE_YEAR, 5, _i)
    if _i <= 24:
        _GENE_TO_DATE[f"MAY{_i}"] = date(_REFERENCE_YEAR, 5, _i)
for _i in range(30, 51):
    _GENE_TO_DATE[f"TAMM{_i}"] = date(2000 + _i, 1, 1)
# JUN family: documented Excel numeric-subtraction quirk. "jun-1" parses
# as "June minus 1" → 2024-05-31; "jun-2" → 2024-05-30; "jun-3" → 2024-05-29.
_GENE_TO_DATE["JUN1"] = date(_REFERENCE_YEAR, 5, 31)
_GENE_TO_DATE["JUN2"] = date(_REFERENCE_YEAR, 5, 30)
_GENE_TO_DATE["JUN3"] = date(_REFERENCE_YEAR, 5, 29)


def excel_autoconvert(value: Any) -> Any:
    """Apply Excel's import-time autoconversion to a single cell value.

    Inputs that aren't strings pass through unchanged.
    Returns: date for date-converted gene names, float for scientific-notation
    IDs, str otherwise.
    """
    if not isinstance(value, str):
        return value
    if value in _GENE_TO_DATE:
        return _GENE_TO_DATE[value]
    if EXPONENT_PATTERN.fullmatch(value):
        try:
            return float(value)
        except ValueError:
            pass
    return value


def corrupt_column(values: list[Any]) -> list[Any]:
    """Apply Excel autoconversion to every value in a column."""
    return [excel_autoconvert(v) for v in values]
