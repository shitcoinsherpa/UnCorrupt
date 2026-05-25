"""Tests for the Frictionless Table Schema sidecar emitter."""
from __future__ import annotations

import json
from datetime import date
from pathlib import Path

import pandas as pd
import pytest
from frictionless import Schema

from uncorrupt.detector import detect
from uncorrupt.schema_export import emit_table_schema, write_schema_sidecar


def _has_field(schema: dict, name: str) -> dict:
    """Find a field in a Frictionless Table Schema dict by name."""
    return next(f for f in schema["fields"] if f["name"] == name)


def test_identifier_column_pinned_to_string_with_pattern() -> None:
    df = pd.DataFrame({
        "symbol": ["BRCA1", "TP53", "EGFR", "MARCHF1", "SEPTIN2", "KRAS", "PTEN"],
        "fold_change": [1.2, 3.4, 0.8, 2.1, 1.7, 4.0, 0.9],
    })
    report = detect(df)
    schema = emit_table_schema(df, report)
    sym = _has_field(schema, "symbol")
    assert sym["type"] == "string"
    assert sym["constraints"]["pattern"] == "^[A-Z][A-Z0-9-]{1,14}$"
    assert "HGNC" in sym["description"]
    fc = _has_field(schema, "fold_change")
    assert fc["type"] == "number"


def test_riken_column_pinned_to_string_with_riken_pattern() -> None:
    riken_ids = [f"23{10000 + i:05d}E13" for i in range(8)]
    df = pd.DataFrame({"gene_id": pd.array(riken_ids, dtype=object)})
    report = detect(df)
    schema = emit_table_schema(df, report)
    f = _has_field(schema, "gene_id")
    assert f["type"] == "string"
    assert f["constraints"]["pattern"] == "^\\d{7}[A-Z]\\d{2}$"
    assert "RIKEN" in f["description"]


def test_non_identifier_columns_get_inferred_types() -> None:
    df = pd.DataFrame({
        "patient_age": [42, 51, 38, 29, 64],
        "diagnosis_date": [
            date(2023, 1, 5), date(2023, 6, 12), date(2024, 3, 1),
            date(2024, 9, 2), date(2024, 11, 30),
        ],
    })
    report = detect(df)
    schema = emit_table_schema(df, report)
    age = _has_field(schema, "patient_age")
    assert age["type"] == "integer"
    date_field = _has_field(schema, "diagnosis_date")
    assert date_field["type"] == "date"
    # No identifier-pinning on these columns
    assert "constraints" not in age or not age.get("constraints")


def test_emitted_schema_validates_against_frictionless_spec() -> None:
    df = pd.DataFrame({
        "symbol": ["BRCA1", "TP53", "EGFR", "MARCHF1", "SEPTIN2"],
        "value": [1.0, 2.0, 3.0, 4.0, 5.0],
    })
    report = detect(df)
    schema_dict = emit_table_schema(df, report)
    # Frictionless can parse it
    parsed = Schema.from_descriptor(schema_dict)
    assert parsed is not None
    assert len(parsed.fields) == 2


def test_sidecar_round_trips_through_disk(tmp_path: Path) -> None:
    df = pd.DataFrame({
        "symbol": ["BRCA1", "TP53", "EGFR", "MARCHF1", "SEPTIN2"],
        "value": [1.0, 2.0, 3.0, 4.0, 5.0],
    })
    report = detect(df)
    sidecar = tmp_path / "mydata.csv.schema.json"
    write_schema_sidecar(df, report, sidecar)
    assert sidecar.exists()
    loaded = json.loads(sidecar.read_text())
    assert loaded["fields"][0]["name"] == "symbol"
    assert loaded["fields"][0]["type"] == "string"


def test_schema_handles_special_character_column_names() -> None:
    """Real-world supplementary files use column headers with spaces,
    parentheses, dots, slashes, and quotes ("Gene Symbol", "log2(FC)",
    "p-value", "% identity"). The schema emitter must preserve these
    verbatim, Frictionless must accept them, and the round-trip through
    JSON must not corrupt them.

    Regression check: every schema test previously used clean snake_case
    names, hiding any quoting bugs."""
    df = pd.DataFrame({
        "Gene Symbol": ["BRCA1", "TP53", "EGFR", "MARCHF1", "SEPTIN2"],
        "log2(FC)": [1.2, 3.4, 0.8, 2.1, 1.7],
        "p-value": [0.001, 0.01, 0.05, 0.005, 0.02],
        "% identity": [99.5, 88.2, 76.4, 95.1, 80.0],
        "gene/ortholog (mouse)": ["Brca1", "Trp53", "Egfr", "March1", "Septin2"],
    })
    report = detect(df)
    schema = emit_table_schema(df, report)
    # Every header is preserved verbatim
    field_names = [f["name"] for f in schema["fields"]]
    assert "Gene Symbol" in field_names
    assert "log2(FC)" in field_names
    assert "p-value" in field_names
    assert "% identity" in field_names
    assert "gene/ortholog (mouse)" in field_names
    # The identifier-shaped columns are pinned to string
    sym = _has_field(schema, "Gene Symbol")
    assert sym["type"] == "string"
    assert "constraints" in sym
    # Frictionless accepts the special-character names
    parsed = Schema.from_descriptor(schema)
    assert parsed is not None
    assert len(parsed.fields) == 5


def test_multi_sheet_identifier_column_strips_sheet_prefix() -> None:
    """When report.identifier_columns has 'Sheet!Col' form, the schema field
    should be matched against the bare column name."""
    from uncorrupt.detector import Report
    df = pd.DataFrame({"symbol": ["BRCA1", "TP53"]})
    # Manually construct a report as if from detect_file (multi-sheet)
    rep = Report(
        rows_scanned=2,
        columns_scanned=1,
        identifier_columns=["Sheet1!symbol"],
        suspicions=[],
    )
    schema = emit_table_schema(df, rep)
    sym = _has_field(schema, "symbol")
    assert sym["type"] == "string"
    assert "constraints" in sym
