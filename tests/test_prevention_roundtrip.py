"""Prevention round-trip: schema sidecar prevents re-corruption.

The CLI `uncorrupt schema <file>` emits a Frictionless Table Schema sidecar
that pins identifier columns to type="string" with regex constraints. The
contract: any later tool that re-imports the data (Excel, locale-aware CSV
loader, etc.) will fail validation when it tries to coerce gene symbols to
dates / floats / etc.

These tests verify:
 1. Schema emission produces a valid Frictionless schema (parseable).
 2. Identifier columns get type="string" + a regex pattern.
 3. A CLEAN copy of the data validates against the schema.
 4. A CORRUPTED copy (gene symbol coerced to date) FAILS schema validation.
"""
from __future__ import annotations

import json
import subprocess
import sys
from datetime import date
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from uncorrupt.detector import detect_file
from uncorrupt.schema_export import emit_table_schema, write_schema_sidecar


def _write_xlsx(tmp_path: Path, df: pd.DataFrame, name: str = "data.xlsx") -> Path:
    path = tmp_path / name
    df.to_excel(path, index=False)
    return path


def test_schema_emission_pins_gene_column_to_string(tmp_path):
    """Identifier column with gene symbols gets type=string + pattern."""
    df = pd.DataFrame({"gene_symbol": [
        "ARNT", "EBAG9", "PIGS", "SBK1", "PHB",
        "ZMYM6", "ATP13A1", "APBB1IP", "FOXO1", "BRCA1",
    ]})
    path = _write_xlsx(tmp_path, df)
    report = detect_file(str(path), row_context_boost=False)
    schema = emit_table_schema(df, report, name="test")

    gene_field = next((f for f in schema["fields"]
                       if f["name"] == "gene_symbol"), None)
    assert gene_field is not None
    assert gene_field["type"] == "string", (
        f"Gene column should be pinned to type=string; got {gene_field!r}"
    )
    constraints = gene_field.get("constraints", {})
    assert "pattern" in constraints, (
        f"Gene column should have regex pattern constraint; got {gene_field!r}"
    )


def test_cli_schema_subcommand_writes_sidecar(tmp_path):
    """Invoke the CLI directly and verify it writes a sidecar file."""
    df = pd.DataFrame({"gene_symbol": [
        "ARNT", "EBAG9", "PIGS", "SBK1", "PHB",
        "ZMYM6", "ATP13A1", "APBB1IP", "FOXO1", "BRCA1",
    ]})
    path = _write_xlsx(tmp_path, df)
    out = subprocess.run(
        [sys.executable, "-m", "uncorrupt.cli", "schema", str(path),
         "-o", str(tmp_path)],
        capture_output=True, text=True, timeout=600,
    )
    assert out.returncode == 0, f"CLI failed: {out.stderr}"
    sidecars = list(tmp_path.glob("*.schema.json"))
    assert len(sidecars) >= 1
    schema = json.loads(sidecars[0].read_text())
    assert schema["fields"][0]["name"] == "gene_symbol"


def _validate_with_frictionless(csv_path: Path, sidecar_path: Path,
                                tmp_path: Path):
    """Run frictionless validate from inside tmp_path so relative
    paths are 'safe' (frictionless's default policy rejects absolute
    paths outside the cwd)."""
    import os

    from frictionless import Resource
    cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        resource = Resource(csv_path.name, schema=sidecar_path.name)
        return resource.validate()
    finally:
        os.chdir(cwd)


def test_corrupted_copy_fails_schema_validation(tmp_path):
    """The schema rejects a corrupted version of the data : a gene symbol
    coerced to a date string violates the regex pattern."""
    df_clean = pd.DataFrame({"gene_symbol": [
        "ARNT", "EBAG9", "PIGS", "SBK1", "PHB",
        "ZMYM6", "ATP13A1", "APBB1IP", "FOXO1", "BRCA1",
    ]})
    path_clean = _write_xlsx(tmp_path, df_clean, "clean.xlsx")
    report = detect_file(str(path_clean), row_context_boost=False)
    sidecar = write_schema_sidecar(
        df_clean, report, tmp_path / "clean.schema.json", name="clean"
    )

    df_corrupt = df_clean.copy()
    df_corrupt.loc[0, "gene_symbol"] = "2023-09-07"
    csv_corrupt = tmp_path / "corrupt.csv"
    df_corrupt.to_csv(csv_corrupt, index=False)

    rpt = _validate_with_frictionless(csv_corrupt, sidecar, tmp_path)
    assert not rpt.valid, (
        "Corrupted copy should fail schema validation (gene_symbol regex). "
        f"Report: {rpt!r}"
    )


def test_clean_copy_validates_against_schema(tmp_path):
    """A clean copy of the data validates against the same schema."""
    df_clean = pd.DataFrame({"gene_symbol": [
        "ARNT", "EBAG9", "PIGS", "SBK1", "PHB",
        "ZMYM6", "ATP13A1", "APBB1IP", "FOXO1", "BRCA1",
    ]})
    path_clean = _write_xlsx(tmp_path, df_clean, "clean.xlsx")
    report = detect_file(str(path_clean), row_context_boost=False)
    sidecar = write_schema_sidecar(
        df_clean, report, tmp_path / "clean.schema.json", name="clean"
    )

    csv_path = tmp_path / "clean.csv"
    df_clean.to_csv(csv_path, index=False)
    rpt = _validate_with_frictionless(csv_path, sidecar, tmp_path)
    assert rpt.valid, (
        f"Clean copy should pass schema validation; got {rpt!r}"
    )


def test_cli_detect_returns_correct_exit_code(tmp_path):
    """Exit codes per CLI contract:
    0 = no corruption | 1 = high-confidence flags | 2 = mid only."""
    df_clean = pd.DataFrame({"gene_symbol": ["BRCA1", "ARNT", "TP53"] * 4})
    path_clean = _write_xlsx(tmp_path, df_clean, "clean.xlsx")
    out = subprocess.run(
        [sys.executable, "-m", "uncorrupt.cli", "detect", str(path_clean)],
        capture_output=True, text=True, timeout=600,
    )
    assert out.returncode == 0, (
        f"Clean file should exit 0; got {out.returncode}\n"
        f"stdout: {out.stdout}\nstderr: {out.stderr}"
    )

    # Two corrupted date cells in a gene column so the column-corroboration
    # gate kicks in and emits at conf >= 0.30 (a single isolated date is
    # correctly demoted by the lone-date rule).
    df_corrupt = pd.DataFrame({"gene_symbol": [
        "ARNT", "EBAG9", "PIGS", "SBK1", "PHB",
        date(2024, 3, 1),  # MARCH1 corruption
        date(2024, 3, 9),  # MARCH9 corruption : gives column corroboration
        "ATP13A1", "APBB1IP", "FOXO1",
    ]})
    path_corrupt = _write_xlsx(tmp_path, df_corrupt, "corrupt.xlsx")
    out = subprocess.run(
        [sys.executable, "-m", "uncorrupt.cli", "detect", str(path_corrupt),
         "--no-boost"],
        capture_output=True, text=True, timeout=600,
    )
    assert out.returncode in (1, 2), (
        f"Corrupted file should exit 1 or 2; got {out.returncode}\n"
        f"stdout: {out.stdout[:500]}"
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
