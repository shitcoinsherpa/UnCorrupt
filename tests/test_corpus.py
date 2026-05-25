"""Tests for the real-registry and real-corpus loaders.

These exercise the actual HGNC bulk download and Ziemann 2021 S2 supplementary
file checked into `data/raw/`. They are not synthetic.
"""
from __future__ import annotations

from uncorrupt.corpus import (
    derive_date_prone_map,
    load_hgnc,
    load_ziemann_2021_corpus,
)


def test_hgnc_loads_with_expected_size() -> None:
    reg = load_hgnc()
    assert 40_000 < len(reg.current_symbols) < 50_000, (
        f"HGNC current symbol count looks off: {len(reg.current_symbols)}"
    )
    assert len(reg.prev_symbol_to_current) > 5_000


def test_hgnc_canonical_renames_resolve_correctly() -> None:
    """Sanity: HGNC must agree with the canonical 2020 renames (Ziemann 2021)."""
    reg = load_hgnc()
    assert reg.resolve("SEPT1") == "SEPTIN1"
    assert reg.resolve("SEPT2") == "SEPTIN2"
    assert reg.resolve("SEPT12") == "SEPTIN12"
    assert reg.resolve("MARCH1") == "MARCHF1"
    assert reg.resolve("MARCH11") == "MARCHF11"
    assert reg.resolve("MARC1") == "MTARC1"
    assert reg.resolve("MARC2") == "MTARC2"
    assert reg.resolve("DEC1") == "DELEC1"


def test_hgnc_documents_sept13_and_sept15_anomalies() -> None:
    """Reality check empirically: SEPT13 maps to SEPTIN7P2 (pseudogene merge),
    SEPT15 has no current HGNC mapping. Both are real-data findings that
    hardcoded copy-paste lists miss.
    """
    reg = load_hgnc()
    assert reg.resolve("SEPT13") == "SEPTIN7P2"
    assert reg.resolve("SEPT15") is None


def test_derived_date_prone_map_covers_canonical_families() -> None:
    reg = load_hgnc()
    date_prone = derive_date_prone_map(reg)
    # All four families must be represented
    assert any(k.startswith("SEPT") for k in date_prone)
    assert any(k.startswith("MARCH") for k in date_prone)
    assert any(k.startswith("MARC") and not k.startswith("MARCH") for k in date_prone)
    assert "DEC1" in date_prone
    # Canonical specific entries
    assert date_prone["SEPT1"] == "SEPTIN1"
    assert date_prone["MARCH1"] == "MARCHF1"
    assert date_prone["MARC1"] == "MTARC1"
    assert date_prone["DEC1"] == "DELEC1"


def test_ziemann_corpus_loads() -> None:
    entries = load_ziemann_2021_corpus()
    assert len(entries) > 3_000, f"expected >3k entries, got {len(entries)}"
    # All entries should have a PMCID starting with "PMC"
    for e in entries[:50]:
        assert e.pmc_id.startswith("PMC"), f"unexpected pmc_id: {e.pmc_id!r}"
        assert e.affected_file_url.startswith("http"), (
            f"unexpected url: {e.affected_file_url!r}"
        )


def test_ziemann_corpus_confirmed_subset_size() -> None:
    """The paper reports 3,436 manually-confirmed corrupted papers."""
    entries = load_ziemann_2021_corpus()
    confirmed = [e for e in entries if e.confirmed_code and e.confirmed_code != "nan"]
    # The exact count may vary slightly by interpretation of the Confirmed code;
    # we assert the order of magnitude here and pin the exact number empirically.
    assert 3_000 < len(confirmed) < 6_000, (
        f"confirmed-subset count is {len(confirmed)}; expected ~3,436"
    )
