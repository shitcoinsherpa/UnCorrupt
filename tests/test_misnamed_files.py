"""Tests for content-sniffing of misnamed Excel files."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from uncorrupt.app import (
    UnrecoverableFile,
    _load_all_sheets,
    _sniff_extension_override,
)
from uncorrupt.column_classifier import classify_column


def test_sniff_detects_html_placeholder(tmp_path: Path) -> None:
    p = tmp_path / "fake.xlsx"
    p.write_bytes(
        b"<html><head><title>Preparing to download ...</title></head><body></body></html>"
    )
    with pytest.raises(UnrecoverableFile, match="HTML placeholder"):
        _sniff_extension_override(str(p))


def test_sniff_detects_recaptcha_placeholder(tmp_path: Path) -> None:
    p = tmp_path / "fake.xls"
    p.write_bytes(
        b"<!doctype html><html><head><base href=\"https://www.google.com/recaptcha/challengepage/\"></html>"
    )
    with pytest.raises(UnrecoverableFile, match="recaptcha"):
        _sniff_extension_override(str(p))


def test_sniff_detects_tsv_with_xls_extension(tmp_path: Path) -> None:
    p = tmp_path / "fake.xls"
    p.write_bytes(b"col1\tcol2\tcol3\nA\t1\t2\nB\t3\t4\n")
    assert _sniff_extension_override(str(p)) == ".tsv"


def test_load_all_sheets_routes_tsv_in_xls_extension(tmp_path: Path) -> None:
    p = tmp_path / "fake.xls"
    p.write_text("symbol\tvalue\nBRCA1\t1.5\nTP53\t2.7\n")
    sheets = _load_all_sheets(str(p))
    assert "_sheet0" in sheets
    df = sheets["_sheet0"]
    assert list(df.columns) == ["symbol", "value"]
    assert len(df) == 2


def test_load_all_sheets_raises_for_unrecoverable_html(tmp_path: Path) -> None:
    p = tmp_path / "fake.XLSX"
    p.write_bytes(
        b"<html><head><title>Preparing to download ...</title></head></html>"
    )
    with pytest.raises(UnrecoverableFile):
        _load_all_sheets(str(p))


def test_column_classifier_recognises_numeric_strings_as_measurement() -> None:
    """A TSV column of integer strings (`'-96'`, `'-148592'`) should classify
    as `measurement`, not `free_text` : otherwise Pass 2 leaks into the
    column and decodes the integers as uncorrupt serials."""
    series = pd.Series(["-96", "-148592", "-528", "-564", "-5124", "12.5", "1e6", "1234"])
    cls = classify_column(series, "FEATURE_TO_PEAK_DISTANCE")
    assert cls.column_type == "measurement"


def test_column_classifier_treats_zero_as_placeholder_not_measurement() -> None:
    """30 gene symbols + 170 zero placeholders should classify as gene_symbol : 
    bare 0 is a 'no measurement' sentinel in omics columns, not a measurement
    that should outvote gene-symbol evidence (PMC4079602 / Unnamed: 16 case)."""
    series = pd.Series(
        ["INSR", "RBM39", "SIM1", "TXNIP", "TFR2", "SNX27", "DDEF1",
         "BRCA1", "TP53", "EGFR"] * 3 + [0] * 170
    )
    cls = classify_column(series, "Unnamed: 16")
    assert cls.column_type == "gene_symbol", (
        f"expected gene_symbol, got {cls.column_type} ({cls.reason})"
    )


def test_column_classifier_treats_string_placeholders_as_missing() -> None:
    """`'NA'`, `'-'`, `'.'` are missing-value sentinels : they should not
    contribute to free_text/string-other classification."""
    series = pd.Series(
        ["BRCA1", "TP53", "EGFR", "MARCHF1", "SEPTIN2", "PTEN", "KRAS"] * 3 +
        ["NA", "-", ".", "n/a"] * 25
    )
    cls = classify_column(series, "gene_id")
    assert cls.column_type == "gene_symbol"
