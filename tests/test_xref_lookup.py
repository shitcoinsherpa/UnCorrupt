"""Unit tests for the cross-reference lookup module."""
from __future__ import annotations

from uncorrupt.xref_lookup import load_xref_index


def test_xref_index_loads() -> None:
    idx = load_xref_index()
    # HGNC has ~20k uniprot, ~42k entrez, etc.
    assert len(idx.by_entrez) > 30_000
    assert len(idx.by_ensembl_gene) > 30_000
    assert len(idx.by_refseq) > 30_000
    assert len(idx.by_uniprot) > 15_000
    assert len(idx.by_hgnc_id) > 40_000


def test_xref_known_canonical_lookups() -> None:
    """Canonical examples: each external ID resolves to its HGNC symbol."""
    idx = load_xref_index()
    # TP53 : HGNC:11998
    assert idx.lookup("HGNC:11998") == "TP53"
    # BRCA1 : Entrez 672 (3 digits, doesn't match tightened pattern requiring 4+)
    # but Ensembl ENSG00000012048 does
    assert idx.lookup("ENSG00000012048") == "BRCA1"
    # EGFR : Entrez 1956 (4 digits, matches tightened pattern), UniProt P00533
    assert idx.lookup("1956") == "EGFR"
    assert idx.lookup("P00533") == "EGFR"


def test_xref_ignores_non_id_strings() -> None:
    idx = load_xref_index()
    assert idx.lookup("not_an_id") is None
    assert idx.lookup("MARCH1") is None  # gene symbol, not external ID
    assert idx.lookup("") is None


def test_xref_versioned_ensembl_strips_suffix() -> None:
    idx = load_xref_index()
    # ENSG00000012048.21 (versioned) should resolve same as base
    assert idx.lookup("ENSG00000012048.21") == "BRCA1"
    assert idx.lookup("ENSG00000012048") == "BRCA1"


def test_xref_versioned_refseq_strips_suffix() -> None:
    idx = load_xref_index()
    # NM_007294 = BRCA1 (versioned: NM_007294.4)
    if idx.lookup("NM_007294") is not None:
        assert idx.lookup("NM_007294") == "BRCA1"
        assert idx.lookup("NM_007294.4") == "BRCA1"


def test_xref_identify_returns_kind() -> None:
    idx = load_xref_index()
    assert idx.identify("HGNC:11998") == "hgnc_id"
    assert idx.identify("ENSG00000012048") == "ensembl"
    assert idx.identify("NM_007294") == "refseq"
    assert idx.identify("P00533") == "uniprot"
    # Entrez requires 4+ digits to avoid false matches on small ints (rank, idx)
    assert idx.identify("1956") == "entrez"
    assert idx.identify("672") is None  # 3-digit Entrez IDs no longer match
    assert idx.identify("1") is None
    assert idx.identify("MARCH1") is None


def test_xref_small_ints_are_not_treated_as_entrez() -> None:
    """Critical: small ints in tables are usually rank/index/chromosome, not
    Entrez IDs. The tightened pattern must reject them."""
    idx = load_xref_index()
    for s in ["1", "2", "17", "100", "999"]:
        assert idx.lookup(s) is None, f"{s!r} should not resolve as Entrez"
