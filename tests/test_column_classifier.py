"""Unit tests for the value-pattern column-type classifier."""
from __future__ import annotations

from datetime import date

import pandas as pd

from uncorrupt.column_classifier import classify_column


def test_gene_symbol_column() -> None:
    series = pd.Series(["BRCA1", "TP53", "EGFR", "KRAS", "PTEN", "MARCHF1"])
    result = classify_column(series, "symbol")
    assert result.column_type == "gene_symbol"


def test_entrez_column_correctly_classified() -> None:
    """The previously-misclassified EntrezGeneID case. Must classify as 'entrez',
    not 'gene_symbol'."""
    # Entrez IDs are 4+ digit integers
    series = pd.Series(["55801", "1956", "672", "7157", "5290"])
    result = classify_column(series, "EntrezGeneID")
    assert result.column_type == "entrez"


def test_refseq_column() -> None:
    series = pd.Series(["NM_007294", "NM_001354609", "NM_005228", "NM_004985", "NM_000314"])
    result = classify_column(series, "RefSeq")
    assert result.column_type == "refseq"


def test_uniprot_column() -> None:
    series = pd.Series(["P38398", "P04637", "P00533", "P01116", "P60484"])
    result = classify_column(series, "UniProt")
    assert result.column_type == "uniprot"


def test_ensembl_column() -> None:
    series = pd.Series([
        "ENSG00000012048", "ENSG00000141510", "ENSG00000146648",
        "ENSG00000133703", "ENSG00000171862",
    ])
    result = classify_column(series, "ensembl_gene_id")
    assert result.column_type == "ensembl"


def test_measurement_column_with_id_header() -> None:
    """Critical case: column header looks identifier-like (has 'id'), but the
    values are clearly numeric measurements. Must NOT classify as gene-symbol."""
    series = pd.Series([48567.0, 44701.5, 34554.2, 41888.7, 49685.3, 38714.1])
    result = classify_column(series, "luciferase_value_id")
    assert result.column_type == "measurement"


def test_date_column() -> None:
    series = pd.Series([
        date(2023, 1, 5), date(2023, 6, 12), date(2024, 3, 1),
        date(2024, 9, 2), date(2024, 11, 30),
    ])
    result = classify_column(series, "diagnosis_date")
    assert result.column_type == "date"


def test_empty_column() -> None:
    series = pd.Series([None, None, None])
    result = classify_column(series, "x")
    assert result.column_type == "empty"


def test_free_text_column() -> None:
    series = pd.Series(["lorem ipsum", "dolor sit", "amet consectetur", "elit sed do"])
    result = classify_column(series, "description")
    assert result.column_type == "free_text"


def test_mixed_column_majority_gene_symbols_minority_measurements() -> None:
    """Realistic real-world shape: a column whose intent is gene symbols
    but the file author accidentally interspersed numeric measurements
    (e.g., fold-change values). The classifier must commit to
    `gene_symbol` because the majority of values are gene-like : the
    minority measurements remain anomalies the detector can flag.

    Regression check: the classifier previously had no test for
    mixed-content columns; pure cases hid a decision-rule weakness."""
    series = pd.Series([
        "BRCA1", "TP53", "EGFR", "KRAS", "MYC", "PTEN", "RB1", "ATM",
        2.34, 5.67,  # 20% numeric measurements
    ])
    result = classify_column(series, "Gene")
    assert result.column_type == "gene_symbol"
    # Confidence should reflect that not every cell is gene-like.
    assert result.confidence < 1.0


def test_mixed_column_majority_measurements_minority_genes() -> None:
    """The inverse : mostly measurements with a couple of stray gene-like
    strings. Must classify as `measurement` (the dominant signal). This
    pins the rule that the classifier follows the majority pattern."""
    series = pd.Series([
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
        "BRCA1", "TP53",  # 20% gene-like
    ])
    result = classify_column(series, "fold_change")
    assert result.column_type == "measurement"
