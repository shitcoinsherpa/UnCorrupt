"""Tests for the batch-correction module."""
from __future__ import annotations

from datetime import date
from pathlib import Path

import openpyxl

from uncorrupt.corrections import (
    Correction,
    _all_candidates,
    _family_of,
    _first_candidate,
    apply_corrections,
    compute_analytics,
    corrections_to_dataframe,
    dataframe_to_decisions,
    propose_corrections,
)


def _write_xlsx(path: Path, rows: list[list], fmts: dict[str, str] | None = None) -> Path:
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "data"
    for r in rows:
        ws.append(r)
    if fmts:
        for coord, fmt in fmts.items():
            ws[coord].number_format = fmt
    wb.save(path)
    return path


def test_first_candidate_picks_modern_symbol() -> None:
    assert _first_candidate("SEPT2 | SEP2") == "SEPT2"
    assert _first_candidate("MARCH1 | MARC1") == "MARCH1"
    assert _first_candidate("") == ""
    assert _first_candidate(None) == ""


def test_all_candidates_splits_pipe_separated() -> None:
    assert _all_candidates("SEPT2 | SEP2") == ["SEPT2", "SEP2"]
    assert _all_candidates("OCT4") == ["OCT4"]
    assert _all_candidates(None) == []


def test_propose_corrections_returns_one_per_suspicion(tmp_path: Path) -> None:
    fp = _write_xlsx(
        tmp_path / "f.xlsx",
        [
            ["symbol"],
            ["BRCA1"], ["TP53"],
            [date(2024, 9, 2)],   # SEPT2 corruption
            ["EGFR"], ["KRAS"], ["PTEN"],
        ],
        fmts={"A4": "d-mmm"},
    )
    corrs = propose_corrections(str(fp))
    assert len(corrs) == 1
    c = corrs[0]
    # Canonical leads: SEPTIN2 is the modern HGNC symbol
    assert c.proposed_value == "SEPTIN2"
    assert "SEPTIN2" in c.all_candidates
    assert "SEPT2" in c.all_candidates  # historical form preserved as alt
    assert c.column == "symbol"
    assert c.row == 2  # pandas 0-based: header row 0 = "symbol", data starts row 0... actually row 2


def test_apply_corrections_audit_log_records_full_provenance(tmp_path: Path) -> None:
    """The audit log JSON must record original value, corrected value, confidence,
    sheet, column, row, and all candidates offered : for every applied correction.
    This is the chain-of-custody record; missing fields would break reproducibility.
    """
    import json
    fp = _write_xlsx(
        tmp_path / "f.xlsx",
        [["symbol"], ["BRCA1"], ["TP53"], [date(2024, 9, 2)], ["EGFR"]],
        fmts={"A4": "d-mmm"},
    )
    corrs = propose_corrections(str(fp))
    _, audit, summary = apply_corrections(str(fp), corrs, output_dir=tmp_path)
    audit_data = json.loads(Path(audit).read_text())

    # Top-level chain-of-custody fields
    assert "source_file" in audit_data
    assert "corrected_file" in audit_data
    assert "applied_at" in audit_data
    assert audit_data["n_applied"] >= 1
    assert "applied" in audit_data and audit_data["applied"]

    # Each applied correction has full provenance
    for a in audit_data["applied"]:
        assert "sheet" in a
        assert "openpyxl_coord" in a
        assert "column" in a
        assert "row" in a
        assert "original_value" in a, f"audit row missing original_value: {a}"
        assert "corrected_to" in a, f"audit row missing corrected_to: {a}"
        assert "confidence" in a, f"audit row missing confidence: {a}"
        assert "all_candidates_offered" in a, (
            f"audit row missing all_candidates_offered: {a}"
        )
        assert "kind" in a, f"audit row missing kind: {a}"
        # The corrected_to value is non-empty
        assert a["corrected_to"], f"corrected_to is empty: {a}"


def test_apply_corrections_writes_corrected_xlsx(tmp_path: Path) -> None:
    fp = _write_xlsx(
        tmp_path / "f.xlsx",
        [
            ["symbol"],
            ["BRCA1"], ["TP53"],
            [date(2024, 9, 2)],   # SEPT2
            [date(2024, 3, 1)],   # MARCH1
            ["EGFR"],
        ],
        fmts={"A4": "d-mmm", "A5": "d-mmm"},
    )
    corrs = propose_corrections(str(fp))
    assert len(corrs) >= 2
    corrected, audit, summary = apply_corrections(str(fp), corrs, output_dir=tmp_path)
    assert Path(corrected).exists()
    assert Path(audit).exists()
    assert summary["n_applied"] >= 2

    # Verify the corrected xlsx actually has the suggested gene symbol
    wb = openpyxl.load_workbook(corrected, data_only=False)
    ws = wb.active
    sept_cell = ws["A4"]
    march_cell = ws["A5"]
    # Canonical leads : writeback uses the modern HGNC symbol
    assert sept_cell.value == "SEPTIN2"
    assert march_cell.value == "MARCHF1"
    # Number format was cleared (no longer a date display)
    assert "d" not in (sept_cell.number_format or "").lower()


def test_apply_corrections_skips_unaccepted_rows(tmp_path: Path) -> None:
    """Only rows in `decisions` get applied : anything missing is untouched."""
    fp = _write_xlsx(
        tmp_path / "f.xlsx",
        [
            ["symbol"],
            ["BRCA1"], ["TP53"],
            [date(2024, 9, 2)],   # SEPT2 : will be fixed
            [date(2024, 3, 1)],   # MARCH1 : will NOT be fixed
            ["EGFR"],
        ],
        fmts={"A4": "d-mmm", "A5": "d-mmm"},
    )
    corrs = propose_corrections(str(fp))
    sept_only = [c for c in corrs if c.original_value.month == 9]
    assert len(sept_only) == 1
    corrected, audit, summary = apply_corrections(str(fp), sept_only, output_dir=tmp_path)
    wb = openpyxl.load_workbook(corrected, data_only=False)
    ws = wb.active
    # SEPT2 cell corrected to modern canonical SEPTIN2
    assert ws["A4"].value == "SEPTIN2"
    # MARCH1 cell untouched (still a date)
    from datetime import datetime
    assert isinstance(ws["A5"].value, (date, datetime))


def test_corrections_to_dataframe_marks_high_confidence_accepted(tmp_path: Path) -> None:
    """The default-accept logic ticks the box on confidence >= 0.85."""
    c_high = Correction(
        sheet=None, column="symbol", row=2,
        original_value="something", proposed_value="MARCH3",
        all_candidates=["MARCH3"], confidence=0.95,  # >= 0.85
        kind="gene-date", reason="...",
    )
    c_low = Correction(
        sheet=None, column="symbol", row=4,
        original_value="something", proposed_value="SEPT2",
        all_candidates=["SEPT2", "SEP2"], confidence=0.5,  # < 0.85
        kind="gene-date", reason="...",
    )
    df = corrections_to_dataframe([c_high, c_low])
    assert df.iloc[0]["accept"] is True or df.iloc[0]["accept"] == True  # noqa: E712
    assert df.iloc[1]["accept"] is False or df.iloc[1]["accept"] == False  # noqa: E712


def test_dataframe_to_decisions_respects_user_edits(tmp_path: Path) -> None:
    """If the user edits `proposed` and `accept`, the resulting decisions
    reflect those edits."""
    corrs = [
        Correction(sheet=None, column="symbol", row=2,
                    original_value="x", proposed_value="SEPT2",
                    all_candidates=["SEPT2", "SEP2"], confidence=0.5,
                    kind="gene-date", reason="..."),
    ]
    df = corrections_to_dataframe(corrs)
    # User unchecks and re-checks with a different candidate
    df.at[0, "accept"] = True
    df.at[0, "proposed"] = "SEP2"  # picked the older alias
    decisions = dataframe_to_decisions(df, corrs)
    assert len(decisions) == 1
    assert decisions[0].proposed_value == "SEP2"


def test_family_of_buckets_known_excel_corruption_families() -> None:
    assert _family_of("SEPT2") == "SEPT"
    assert _family_of("SEP2") == "SEP"
    assert _family_of("MARCH1") == "MARCH"
    assert _family_of("MARC1") == "MARC"
    assert _family_of("OCT4") == "OCT"
    assert _family_of("DEC1") == "DEC"
    assert _family_of("NOV1") == "NOV"
    assert _family_of("APR3") == "APR"
    assert _family_of("FEB4") == "FEB"
    # Modern post-rename symbols
    assert _family_of("SEPTIN2") == "SEPTIN"
    assert _family_of("MARCHF1") == "MARCHF"
    assert _family_of("POU2F2") == "POU2F"
    # Non-corruption symbols fall through
    assert _family_of("BRCA1") == "OTHER"


def test_compute_analytics_summarises_batch() -> None:
    """Smoke test the analytics returns the expected shape + sane counts."""
    corrs = [
        Correction(sheet="s1", column="Gene", row=0, original_value="x",
                    proposed_value="SEPT2",
                    all_candidates=["SEPT2", "SEP2"], confidence=0.95,
                    kind="gene-date", reason="..."),
        Correction(sheet="s1", column="Gene", row=1, original_value="x",
                    proposed_value="SEPT7",
                    all_candidates=["SEPT7", "SEP7"], confidence=0.5,
                    kind="gene-date", reason="..."),
        Correction(sheet="s1", column="Gene", row=2, original_value="x",
                    proposed_value="MARCH1",
                    all_candidates=["MARCH1", "MARC1"], confidence=0.95,
                    kind="gene-date-serial", reason="..."),
        Correction(sheet="s2", column="Gene", row=10, original_value="x",
                    proposed_value="OCT4",
                    all_candidates=["OCT4"], confidence=0.45,
                    kind="gene-date-string", reason="..."),
    ]
    analytics = compute_analytics(corrs)
    assert analytics["totals"]["n_corrections"] == 4
    assert analytics["totals"]["n_sheets"] == 2
    assert analytics["totals"]["n_columns"] == 2
    assert analytics["totals"]["n_high_confidence"] == 2
    assert analytics["totals"]["n_medium_confidence"] == 1
    assert analytics["totals"]["n_low_confidence"] == 1
    # Kind breakdown
    kind_df = analytics["by_kind_df"]
    assert set(kind_df["kind"]) == {"gene-date", "gene-date-serial", "gene-date-string"}
    # Family breakdown
    fam_df = analytics["by_family_df"]
    assert "SEPT" in set(fam_df["family"])
    assert "MARCH" in set(fam_df["family"])
    assert "OCT" in set(fam_df["family"])
    # Summary markdown is non-empty
    assert "4" in analytics["summary_md"]
    assert "2 sheet" in analytics["summary_md"] or "2** sheet" in analytics["summary_md"]


def test_compute_analytics_empty_corrections_returns_clean_message() -> None:
    analytics = compute_analytics([])
    assert analytics["totals"]["n_corrections"] == 0
    assert "No corruptions" in analytics["summary_md"]
    assert len(analytics["by_kind_df"]) == 0


def test_dataframe_to_decisions_excludes_unaccepted(tmp_path: Path) -> None:
    corrs = [
        Correction(sheet=None, column="symbol", row=2,
                    original_value="x", proposed_value="SEPT2",
                    all_candidates=["SEPT2"], confidence=0.95,
                    kind="gene-date", reason="..."),
        Correction(sheet=None, column="symbol", row=4,
                    original_value="y", proposed_value="MARCH1",
                    all_candidates=["MARCH1"], confidence=0.95,
                    kind="gene-date", reason="..."),
    ]
    df = corrections_to_dataframe(corrs)
    df.at[0, "accept"] = True
    df.at[1, "accept"] = False
    decisions = dataframe_to_decisions(df, corrs)
    assert len(decisions) == 1
    assert decisions[0].proposed_value == "SEPT2"
