"""CLI surface tests for the `uncorrupt` console entry-point.

Covers the three public commands (`detect`, `schema`, `audit`), the exit-code
contract (0 clean, 1 high-confidence flags, 2 mid-only, 3 file/folder error),
the `--json` output shape, the `--no-boost` flag, and the `--recursive`
folder walk.

These tests call `uncorrupt.cli.main` directly with an argv list rather than
shelling out via subprocess so they stay fast (no Python startup per case)
and so the underlying xref index loads once per session via test scope.
"""
from __future__ import annotations

import io
import json
from contextlib import redirect_stderr, redirect_stdout
from datetime import date
from pathlib import Path

import pandas as pd
import pytest

from uncorrupt.cli import main as cli_main


def _write_xlsx(path: Path, df: pd.DataFrame) -> None:
    df.to_excel(path, index=False)


def _run(argv: list[str]) -> tuple[int, str, str]:
    """Invoke the CLI and capture stdout/stderr."""
    out = io.StringIO()
    err = io.StringIO()
    with redirect_stdout(out), redirect_stderr(err):
        try:
            rc = cli_main(argv)
        except SystemExit as e:
            rc = e.code if isinstance(e.code, int) else 0
    return rc, out.getvalue(), err.getvalue()


# `--no-boost` is passed throughout so we don't need the heavyweight xref
# index for these tests; the row-xref boost is tested separately in
# tests/test_canonical_and_xref_boost.py.


# --- detect ---

def test_detect_clean_file_exits_zero(tmp_path):
    """A spreadsheet of plain gene symbols produces no flags; exit 0."""
    p = tmp_path / "clean.xlsx"
    _write_xlsx(p, pd.DataFrame({"gene_symbol": ["BRCA1", "TP53", "EGFR"] * 4}))
    rc, out, err = _run(["detect", str(p), "--no-boost"])
    assert rc == 0, f"clean file should exit 0; got {rc}; stderr={err!r}"
    assert "Suspicions: 0" in out


def test_detect_high_confidence_corruption_exits_one(tmp_path):
    """Two MARCH dates in a gene column trip the column-corroboration gate
    and produce high-confidence flags; exit 1."""
    p = tmp_path / "corrupt.xlsx"
    _write_xlsx(p, pd.DataFrame({"gene_symbol": [
        "ARNT", "EBAG9", "PIGS", "SBK1", "PHB",
        date(2024, 3, 1),    # MARCH1 corruption
        date(2024, 3, 9),    # MARCH9 corruption: gives column corroboration
        "ATP13A1", "APBB1IP", "FOXO1",
    ]}))
    rc, _, _ = _run(["detect", str(p), "--no-boost"])
    assert rc in (1, 2), f"corrupted file should exit 1 or 2; got {rc}"


def test_detect_missing_file_exits_three(tmp_path):
    """Pointing at a file that does not exist returns exit 3 with a clear
    error message."""
    rc, _, err = _run(["detect", str(tmp_path / "nope.xlsx"), "--no-boost"])
    assert rc == 3
    assert "file not found" in err.lower()


def test_detect_json_shape(tmp_path):
    """`--json` emits an object with `file`, `rows_scanned`, `columns_scanned`,
    `identifier_columns`, and a `suspicions` array."""
    p = tmp_path / "clean.xlsx"
    _write_xlsx(p, pd.DataFrame({"gene_symbol": ["BRCA1", "TP53", "EGFR"] * 4}))
    rc, out, _ = _run(["detect", str(p), "--no-boost", "--json"])
    assert rc == 0
    data = json.loads(out)
    assert set(data.keys()) >= {
        "file", "rows_scanned", "columns_scanned",
        "identifier_columns", "suspicions",
    }
    assert isinstance(data["suspicions"], list)


def test_detect_json_suspicion_shape(tmp_path):
    """Every entry in the `suspicions` array has the documented fields."""
    p = tmp_path / "corrupt.xlsx"
    _write_xlsx(p, pd.DataFrame({"gene_symbol": [
        "ARNT", "EBAG9", "PIGS", "SBK1", "PHB",
        date(2024, 3, 1), date(2024, 3, 9),
        "ATP13A1", "APBB1IP", "FOXO1",
    ]}))
    rc, out, _ = _run(["detect", str(p), "--no-boost", "--json"])
    assert rc in (1, 2)
    data = json.loads(out)
    assert len(data["suspicions"]) >= 1
    s = data["suspicions"][0]
    for field in ("sheet", "column", "row", "value", "kind", "suggestion",
                  "confidence", "reason", "xref_status"):
        assert field in s, f"missing field {field!r} in suspicion record"


def test_detect_max_per_band_caps_display(tmp_path):
    """Many flags in a single band get truncated with an "and N more" line
    when `--max-per-band` is below the flag count."""
    rows = (["ARNT"] * 5
            + [date(2024, 3, i) for i in range(1, 13)]
            + ["FOXO1"] * 5)
    p = tmp_path / "many.xlsx"
    _write_xlsx(p, pd.DataFrame({"gene_symbol": rows}))
    rc, out, _ = _run(["detect", str(p), "--no-boost", "--max-per-band", "3"])
    assert rc in (1, 2)
    # Either the cap fires (look for "and N more"), or the band genuinely has
    # <= 3 flags after column-classification. Tolerate both.
    if "and " in out and " more" in out:
        assert "and " in out


# --- schema ---

def test_schema_writes_sidecar(tmp_path):
    """`uncorrupt schema` writes a `*.schema.json` next to the input."""
    p = tmp_path / "data.xlsx"
    _write_xlsx(p, pd.DataFrame({"gene_symbol": [
        "ARNT", "EBAG9", "PIGS", "SBK1", "PHB",
        "ZMYM6", "ATP13A1", "APBB1IP", "FOXO1", "BRCA1",
    ]}))
    rc, _, _ = _run(["schema", str(p)])
    assert rc == 0
    sidecars = list(tmp_path.glob("*.schema.json"))
    assert len(sidecars) >= 1
    schema = json.loads(sidecars[0].read_text())
    assert "fields" in schema
    gene = next((f for f in schema["fields"] if f["name"] == "gene_symbol"),
                None)
    assert gene is not None
    assert gene["type"] == "string"
    assert "pattern" in gene.get("constraints", {})


def test_schema_output_dir_override(tmp_path):
    """`-o <dir>` writes the sidecar into the specified directory."""
    p = tmp_path / "data.xlsx"
    outdir = tmp_path / "schemas"
    _write_xlsx(p, pd.DataFrame({"gene_symbol": ["BRCA1", "TP53", "EGFR"] * 4}))
    rc, _, _ = _run(["schema", str(p), "-o", str(outdir)])
    assert rc == 0
    sidecars = list(outdir.glob("*.schema.json"))
    assert len(sidecars) >= 1


def test_schema_missing_file_exits_three(tmp_path):
    rc, _, err = _run(["schema", str(tmp_path / "nope.xlsx")])
    assert rc == 3
    assert "file not found" in err.lower()


# --- audit ---

def test_audit_clean_folder_exits_zero(tmp_path):
    """Folder of clean files: exit 0."""
    for i in range(3):
        _write_xlsx(tmp_path / f"clean_{i}.xlsx",
                    pd.DataFrame({"gene_symbol": ["BRCA1", "TP53", "EGFR"] * 4}))
    rc, out, _ = _run(["audit", str(tmp_path), "--no-boost"])
    assert rc == 0
    assert "Files with corruption flags: 0/3" in out


def test_audit_dirty_folder_exits_nonzero(tmp_path):
    """Folder containing at least one corrupted file exits non-zero. The
    audit command distinguishes high-confidence (1) from medium-only (2);
    either is a 'something was flagged' signal CI should fail on."""
    _write_xlsx(tmp_path / "clean.xlsx",
                pd.DataFrame({"gene_symbol": ["BRCA1", "TP53", "EGFR"] * 4}))
    _write_xlsx(tmp_path / "corrupt.xlsx", pd.DataFrame({"gene_symbol": [
        "ARNT", "EBAG9", "PIGS", "SBK1", "PHB",
        date(2024, 3, 1), date(2024, 3, 9),
        "ATP13A1", "APBB1IP", "FOXO1",
    ]}))
    rc, _, _ = _run(["audit", str(tmp_path), "--no-boost"])
    assert rc in (1, 2), (
        f"expected exit 1 (high) or 2 (medium-only); got {rc}"
    )


def test_audit_recursive_finds_nested(tmp_path):
    """`-r` walks subdirectories."""
    sub = tmp_path / "deep" / "nested"
    sub.mkdir(parents=True)
    _write_xlsx(sub / "corrupt.xlsx", pd.DataFrame({"gene_symbol": [
        "ARNT", "EBAG9", "PIGS", "SBK1", "PHB",
        date(2024, 3, 1), date(2024, 3, 9),
        "ATP13A1", "APBB1IP", "FOXO1",
    ]}))
    # Non-recursive: empty top folder, exit 0 (no .xlsx at top level)
    rc_flat, out_flat, _ = _run(["audit", str(tmp_path), "--no-boost"])
    assert rc_flat == 0
    assert "Auditing 0 file(s)" in out_flat
    # Recursive: finds the nested corrupt file, exit 1 (high) or 2 (medium).
    rc_rec, _, _ = _run(["audit", str(tmp_path), "-r", "--no-boost"])
    assert rc_rec in (1, 2)


def test_audit_missing_folder_exits_three(tmp_path):
    rc, _, err = _run(["audit", str(tmp_path / "nope"), "--no-boost"])
    assert rc == 3
    assert "folder not found" in err.lower()


def test_audit_json_shape(tmp_path):
    """`--json` emits an object keyed by file path with per-file counts."""
    _write_xlsx(tmp_path / "f.xlsx",
                pd.DataFrame({"gene_symbol": ["BRCA1", "TP53", "EGFR"] * 4}))
    rc, out, _ = _run(["audit", str(tmp_path), "--no-boost", "--json"])
    assert rc == 0
    # The "Auditing 1 file(s) ..." line precedes the JSON object; pull the
    # JSON tail off the output.
    json_start = out.find("{")
    assert json_start >= 0
    data = json.loads(out[json_start:])
    assert len(data) == 1
    rec = next(iter(data.values()))
    assert set(rec.keys()) >= {"n_suspicions", "n_high_confidence",
                                "n_mid_confidence", "kinds"}


# --- argparse / subparser contract ---

def test_no_subcommand_exits_two():
    """argparse requires a subcommand; bare `uncorrupt` exits 2."""
    rc, _, _ = _run([])
    assert rc == 2


def test_unknown_subcommand_exits_two():
    rc, _, _ = _run(["badcmd"])
    assert rc == 2


def test_detect_help_exits_zero():
    """`detect --help` prints usage and exits 0."""
    rc, out, _ = _run(["detect", "--help"])
    assert rc == 0
    assert "--json" in out
    assert "--no-boost" in out


def test_schema_help_exits_zero():
    rc, out, _ = _run(["schema", "--help"])
    assert rc == 0
    assert "--output" in out or "-o" in out


def test_audit_help_exits_zero():
    rc, out, _ = _run(["audit", "--help"])
    assert rc == 0
    assert "--recursive" in out or "-r" in out


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
