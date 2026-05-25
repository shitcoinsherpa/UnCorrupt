"""Integration tests driving the Gradio app via gradio_client.

Exercise the full upload -> detection -> report pipeline without a browser.
The Playwright/browser layer is a separate, slower test.

Table response shape (from gradio_client):
    {"headers": [...], "data": [[col, row, value, kind, suggestion, confidence, reason], ...]}
Column indices used below:
    0=column  1=row  2=value  3=kind  4=suggestion  5=confidence  6=reason
"""
from __future__ import annotations

from datetime import date
from pathlib import Path

import pandas as pd
import pytest
from gradio_client import Client, handle_file


def _write_xlsx(path: Path, data: dict[str, list]) -> Path:
    df = pd.DataFrame({k: pd.array(v, dtype=object) for k, v in data.items()})
    df.to_excel(path, index=False)
    return path


@pytest.fixture(scope="session")
def client(gradio_url: str) -> Client:
    return Client(gradio_url)


def test_clean_file_yields_no_suspicions(client: Client, tmp_path: Path) -> None:
    f = _write_xlsx(tmp_path / "clean.xlsx", {
        "symbol": ["BRCA1", "TP53", "EGFR", "MARCHF1", "SEPTIN2", "KRAS", "PTEN"],
        "fold_change": [1.2, 3.4, 0.8, 2.1, 1.7, 4.0, 0.9],
    })
    table, summary, _schema = client.predict(file=handle_file(str(f)), api_name="/analyze")
    assert "Suspicions found: 0" in summary
    assert table["data"] == []


def test_gene_date_corruption_is_detected_with_correct_confidence(
    client: Client, tmp_path: Path
) -> None:
    mixed = ["BRCA1", "TP53", date(2024, 9, 2), "MARCHF1", "SEPTIN2", date(2024, 3, 1), "PTEN"]
    f = _write_xlsx(tmp_path / "gene-date.xlsx", {
        "symbol": mixed,
        "fold_change": [1.2, 3.4, 0.8, 2.1, 1.7, 4.0, 0.9],
    })
    table, summary, _schema = client.predict(file=handle_file(str(f)), api_name="/analyze")
    assert "Identifier-shaped columns: ['symbol']" in summary
    assert "Suspicions found: 2" in summary

    rows = table["data"]
    assert len(rows) == 2
    sept_row = next(r for r in rows if "SEPT2" in r[4])
    # SEPT2 + SEP2 + canonical SEPTIN2: unique modern canonical leads the
    # suggestion list, so the 0.5 base gets the +0.1 canonical boost → 0.6.
    assert "SEP2" in sept_row[4]
    assert sept_row[5] == pytest.approx(0.6)
    march_row = next(r for r in rows if "MARCH1" in r[4] and "MARC1" in r[4])
    # MARCH1 has TWO canonical successors (MARCHF1 and MTARC1) : not unique,
    # so confidence stays at the 0.5 ambiguity base.
    assert march_row[5] == pytest.approx(0.5)


def test_riken_float_corruption_is_detected(client: Client, tmp_path: Path) -> None:
    riken_ids = [f"23{10000 + i:05d}E13" for i in range(15)]
    column = riken_ids[:7] + [2.310009e19] + riken_ids[8:]
    f = _write_xlsx(tmp_path / "riken-float.xlsx", {
        "gene_id": column,
        "expression": [i * 0.5 for i in range(15)],
    })
    table, summary, _schema = client.predict(file=handle_file(str(f)), api_name="/analyze")
    assert "gene_id" in summary

    rows = table["data"]
    id_float_rows = [r for r in rows if r[3] == "id-float"]
    assert len(id_float_rows) == 1
    assert "precision lost" in id_float_rows[0][6]


def test_non_identifier_columns_are_not_flagged(client: Client, tmp_path: Path) -> None:
    """Negative control: dates and floats in clearly non-identifier columns must not trigger."""
    f = _write_xlsx(tmp_path / "clinical.xlsx", {
        "patient_age": [42, 51, 38, 29, 64],
        "diagnosis_date": [
            date(2023, 1, 5), date(2023, 6, 12), date(2024, 3, 1),
            date(2024, 9, 2), date(2024, 11, 30),
        ],
    })
    table, summary, _schema = client.predict(file=handle_file(str(f)), api_name="/analyze")
    assert "Suspicions found: 0" in summary
    assert table["data"] == []


# --- Server-free integration path -------------------
# These tests bypass gradio_client / Gradio server entirely by calling
# `analyze()` directly. The full client tests above require the live server
# fixture; this path lets the same wiring be exercised in environments where
# binding to a port is undesirable (e.g., locked-down CI runners).


def test_analyze_serverfree_clean_file(tmp_path: Path) -> None:
    from types import SimpleNamespace
    from uncorrupt.app import analyze
    f = _write_xlsx(tmp_path / "clean.xlsx", {
        "symbol": ["BRCA1", "TP53", "EGFR", "MARCHF1", "SEPTIN2"],
        "value": [1.0, 2.0, 3.0, 4.0, 5.0],
    })
    df, summary, schema_path = analyze(SimpleNamespace(name=str(f)))
    assert "Suspicions found: 0" in summary
    assert len(df) == 0
    assert schema_path is not None
    assert Path(schema_path).exists()


def test_analyze_serverfree_gene_date_corruption(tmp_path: Path) -> None:
    from types import SimpleNamespace
    from uncorrupt.app import analyze
    f = _write_xlsx(tmp_path / "corrupted.xlsx", {
        "symbol": ["BRCA1", "TP53", date(2024, 9, 2), "MARCHF1", "SEPTIN2",
                   date(2024, 3, 1), "PTEN"],
        "value": [1.2, 3.4, 0.8, 2.1, 1.7, 4.0, 0.9],
    })
    df, summary, schema_path = analyze(SimpleNamespace(name=str(f)))
    assert "Suspicions found: 2" in summary
    suggestions = df["suggestion"].tolist()
    assert any("SEPT2" in s for s in suggestions)
    assert any("MARCH1" in s for s in suggestions)
    # Schema sidecar pins the identifier column to string + HGNC pattern
    import json
    schema = json.loads(Path(schema_path).read_text())
    sym = next(f for f in schema["fields"] if f["name"] == "symbol")
    assert sym["type"] == "string"
    assert "pattern" in sym["constraints"]


def test_analyze_serverfree_no_file_returns_empty() -> None:
    from uncorrupt.app import analyze
    df, summary, schema_path = analyze(None)
    assert "No file uploaded" in summary
    assert len(df) == 0
    assert schema_path is None
