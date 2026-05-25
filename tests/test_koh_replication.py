"""Unit tests for the Koh-2022-replication benchmark module.

These exercise the local logic. Network-dependent paths (esearch, EPMC fetch,
OA tarball) are covered by the smoke-test script and the full-corpus run; we
keep unit tests offline so CI doesn't depend on NCBI uptime.
"""
from __future__ import annotations

import io
import zipfile
from pathlib import Path

import pandas as pd

from uncorrupt.corpus import load_hgnc
from uncorrupt.koh_replication import (
    KOH_JOURNALS,
    has_gene_symbols,
    _iter_xlsx_in_zip,
)


def test_koh_journals_count_matches_paper() -> None:
    """Koh 2022 explicitly named 11 journals."""
    assert len(KOH_JOURNALS) == 11


def test_has_gene_symbols_finds_real_hgnc_in_dataframe() -> None:
    hgnc = load_hgnc()
    hgnc_set = hgnc.current_symbols | set(hgnc.prev_symbol_to_current.keys())
    df = pd.DataFrame({
        "name": ["BRCA1", "TP53", "EGFR", "MARCHF1", "SEPTIN2", "KRAS", "PTEN"],
        "other": [1, 2, 3, 4, 5, 6, 7],
    })
    assert has_gene_symbols(df, hgnc_set, min_distinct=5)


def test_has_gene_symbols_rejects_non_gene_data() -> None:
    hgnc = load_hgnc()
    hgnc_set = hgnc.current_symbols | set(hgnc.prev_symbol_to_current.keys())
    df = pd.DataFrame({
        "patient_id": ["P001", "P002", "P003", "P004", "P005"],
        "age": [42, 51, 38, 29, 64],
    })
    assert not has_gene_symbols(df, hgnc_set)


def test_has_gene_symbols_requires_minimum_distinct_count() -> None:
    """A single accidental gene-shape match must not flip the file to 'gene file'."""
    hgnc = load_hgnc()
    hgnc_set = hgnc.current_symbols | set(hgnc.prev_symbol_to_current.keys())
    df = pd.DataFrame({
        "col": ["random", "values", "BRCA1", "more", "stuff"],
    })
    assert not has_gene_symbols(df, hgnc_set, min_distinct=5)


def test_iter_xlsx_in_zip_filters_to_tabular_extensions() -> None:
    """Build a synthetic ZIP with mixed file types; only xlsx/xls/csv/tsv emerge."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr("table1.xlsx", b"PK fake-content")
        z.writestr("table2.csv", b"a,b,c")
        z.writestr("image.png", b"\x89PNG fake")
        z.writestr("doc.pdf", b"%PDF fake")
        z.writestr("nested/legacy.xls", b"\xd0\xcf fake")
        z.writestr("note.txt", b"hello")

    names = [name for name, _ in _iter_xlsx_in_zip(buf.getvalue())]
    assert "table1.xlsx" in names
    assert "table2.csv" in names
    assert "legacy.xls" in names
    assert "image.png" not in names
    assert "doc.pdf" not in names
    assert "note.txt" not in names


def test_iter_xlsx_in_zip_handles_bad_zip_gracefully() -> None:
    """Non-ZIP bytes return empty iterator, no exception."""
    assert list(_iter_xlsx_in_zip(b"<html>not a zip</html>")) == []
