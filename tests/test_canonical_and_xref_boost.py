"""Tests for the HGNC-canonical-ranking + row-context xref-boost mechanisms.

These are the two confidence-improving mechanisms added per Proposals A + B:

  A) `_canonicalize_candidates`: modern HGNC symbol leads, plus a +0.1
     confidence bump when exactly one canonical resolves. Never drops
     ambiguity.
  B) `_apply_row_context_boost`: when another column in the same row holds
     an Ensembl/RefSeq/UniProt/Entrez/HGNC ID that resolves (via xref) to
     the suggested gene, confidence is boosted to min(0.99, conf + 0.4).
"""
from __future__ import annotations

from datetime import date

import pandas as pd
import pytest

from uncorrupt.detector import (
    _apply_row_context_boost,
    _canonicalize_candidates,
)

# === Proposal A : HGNC canonical ranking =====================================


def test_canonicalize_puts_modern_hgnc_first() -> None:
    out, unique = _canonicalize_candidates(["SEPT2", "SEP2"])
    assert out[0] == "SEPTIN2"  # modern HGNC current symbol leads
    assert "SEPT2" in out
    assert "SEP2" in out
    assert unique is True  # exactly one canonical resolution


def test_canonicalize_marchf_rename() -> None:
    out, unique = _canonicalize_candidates(["MARCH1", "MARC1"])
    assert "MARCHF1" in out  # MARCH1 -> MARCHF1
    assert "MTARC1" in out   # MARC1 -> MTARC1
    assert unique is False   # TWO canonicals resolved : not the disambiguating case


def test_canonicalize_preserves_singletons() -> None:
    out, unique = _canonicalize_candidates(["BRCA1"])
    assert out == ["BRCA1"]  # no rename for BRCA1
    assert unique is False


def test_canonicalize_empty_input() -> None:
    out, unique = _canonicalize_candidates([])
    assert out == []
    assert unique is False


def test_canonicalize_dec1_rename() -> None:
    out, unique = _canonicalize_candidates(["DEC1"])
    assert out[0] == "DELEC1"
    assert "DEC1" in out
    assert unique is True


# === Proposal B : row-context xref boost =====================================


class _FakeXrefIndex:
    """Test-friendly xref stub mapping known external IDs to gene symbols.
    Bypasses the 150 MB HGNC-backed loader so tests stay fast and isolated."""

    def __init__(self, mapping: dict[str, str]):
        self._map = mapping

    def lookup(self, identifier: str) -> str | None:
        return self._map.get(identifier.strip())


def _build_fake_report(suspicion_kwargs: dict) -> object:
    """Wrap one Suspicion in a Report-shaped namespace for the boost function."""
    from uncorrupt.detector import Report, Suspicion
    rep = Report(rows_scanned=0, columns_scanned=0)
    rep.suspicions.append(Suspicion(**suspicion_kwargs))
    return rep


def test_row_xref_boost_corroborates_via_refseq() -> None:
    """A RefSeq ID in the same row that resolves to the suggested gene → boost."""
    df = pd.DataFrame({
        "symbol": ["BRCA1", "TP53", date(2024, 9, 2), "EGFR"],
        "refseq": ["NM_007294", "NM_000546", "NM_001008491", "NM_005228"],
    })
    rep = _build_fake_report({
        "column": "symbol", "row": 2, "value": date(2024, 9, 2),
        "kind": "gene-date",
        "suggestion": "SEPTIN2 | SEPT2 | SEP2",
        "confidence": 0.6,
        "reason": "date in identifier column",
        "sheet": "data",
    })
    xref = _FakeXrefIndex({"NM_001008491": "SEPTIN2"})
    n_boosted = _apply_row_context_boost(rep, {"data": df}, xref)
    assert n_boosted == 1
    assert rep.suspicions[0].confidence == pytest.approx(1.0) or \
           rep.suspicions[0].confidence == pytest.approx(0.99)
    assert "row-xref CORROBORATED" in rep.suspicions[0].reason
    assert "SEPTIN2" in rep.suspicions[0].reason


def test_row_xref_boost_downgrades_when_contradicted() -> None:
    """If the row's xrefs resolve to a DIFFERENT gene than the suggestion,
    v0.6.4 downgrades confidence to 0.20 to reflect the contradicting
    evidence (was: leave confidence alone). Empirically caught: the
    `supplementary_table_8_ddac017.xlsx` row-279 case where a publication-
    date column was mis-classified as identifier, the detector decoded
    the date to MARCHF5, but the row's Ensembl ID resolved to GARS1."""
    df = pd.DataFrame({
        "symbol": ["BRCA1", "TP53", date(2024, 9, 2), "EGFR"],
        "refseq": ["NM_007294", "NM_000546", "NM_005228", "NM_001008491"],
    })
    rep = _build_fake_report({
        "column": "symbol", "row": 2, "value": date(2024, 9, 2),
        "kind": "gene-date",
        "suggestion": "SEPTIN2 | SEPT2 | SEP2",
        "confidence": 0.6,
        "reason": "date in identifier column",
        "sheet": "data",
    })
    # NM_005228 resolves to EGFR, not SEPTIN2 : contradicting evidence
    xref = _FakeXrefIndex({"NM_005228": "EGFR"})
    n_boosted = _apply_row_context_boost(rep, {"data": df}, xref)
    assert n_boosted == 0
    assert rep.suspicions[0].confidence == pytest.approx(0.20)
    assert rep.suspicions[0].xref_status == "contradicted"
    assert "CONTRADICTED" in rep.suspicions[0].reason
    assert "EGFR" in rep.suspicions[0].reason


def test_row_xref_boost_silent_on_no_row_xrefs() -> None:
    """No external IDs in the row → no boost, no annotation. Honest silence."""
    df = pd.DataFrame({
        "symbol": ["BRCA1", date(2024, 9, 2)],
        "other": ["foo",   "bar"],
    })
    rep = _build_fake_report({
        "column": "symbol", "row": 1, "value": date(2024, 9, 2),
        "kind": "gene-date", "suggestion": "SEPTIN2 | SEPT2 | SEP2",
        "confidence": 0.6, "reason": "...", "sheet": "data",
    })
    xref = _FakeXrefIndex({})
    n_boosted = _apply_row_context_boost(rep, {"data": df}, xref)
    assert n_boosted == 0


def test_row_xref_boost_caps_confidence_at_099() -> None:
    df = pd.DataFrame({
        "symbol": [date(2024, 9, 2), "EGFR"],
        "ensembl": ["ENSG00000168385", "ENSG00000146648"],
    })
    rep = _build_fake_report({
        "column": "symbol", "row": 0, "value": date(2024, 9, 2),
        "kind": "gene-date", "suggestion": "SEPTIN2 | SEPT2 | SEP2",
        "confidence": 0.95,  # already high
        "reason": "...", "sheet": "data",
    })
    xref = _FakeXrefIndex({"ENSG00000168385": "SEPTIN2"})
    _apply_row_context_boost(rep, {"data": df}, xref)
    assert rep.suspicions[0].confidence == pytest.approx(0.99)


def test_row_xref_boost_handles_renamed_gene_match() -> None:
    """Row xref resolves to SEPTIN2 (modern), suggestion has SEPT2 (deprecated).
    HGNC rename map should still match them."""
    df = pd.DataFrame({
        "symbol": ["BRCA1", date(2024, 9, 2)],
        "refseq": ["NM_007294", "NM_001008491"],
    })
    rep = _build_fake_report({
        "column": "symbol", "row": 1, "value": date(2024, 9, 2),
        "kind": "gene-date",
        "suggestion": "SEPT2",  # ONLY the deprecated form (canonical not run)
        "confidence": 0.5,
        "reason": "...", "sheet": "data",
    })
    # Xref resolves to the MODERN symbol
    xref = _FakeXrefIndex({"NM_001008491": "SEPTIN2"})
    n_boosted = _apply_row_context_boost(rep, {"data": df}, xref)
    assert n_boosted == 1  # rename expansion catches SEPT2 == SEPTIN2
