"""Extensive fresh-files audit of the iter-5 detector.

For every file in three sampling pools, this script independently inspects the
file with openpyxl and compares against the detector's verdict, looking for:

  - True positives  : detector flagged a cell AND openpyxl confirms it's a date
                      (or date-like string, or integer-in-serial-range)
  - False positives : detector flagged a cell but the cell is clean
  - False negatives : detector did NOT flag a cell, but openpyxl finds a
                      date / date-string / serial-shaped integer in a gene-symbol column
  - True negatives  : detector + openpyxl agree the cell is clean

Three pools (all stratified-random from the Koh cache):
  A) all files marked "with corruption" by the cache-walk (278)
  B) 100 random files marked "clean with gene symbols"
  C) 100 random files marked "no gene symbols"

Outputs:
  results/extensive_audit.jsonl       full per-cell evidence
  results/extensive_audit_summary.csv  one row per file: pool, pmc, file, n_flags,
                                       n_tp, n_fp, n_fn, n_potential_fn

Usage:
    python scripts/extensive_fresh_audit.py
    python scripts/extensive_fresh_audit.py --pool A           # just the dirty pool
    python scripts/extensive_fresh_audit.py --n-clean 200 --n-empty 200
"""
from __future__ import annotations

import argparse
import csv
import json
import random
import sys
import time
import warnings
from datetime import UTC, date, datetime, timedelta
from pathlib import Path

warnings.filterwarnings("ignore")

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import openpyxl  # noqa: E402

import pandas as pd

from uncorrupt.app import _load_all_sheets, UnrecoverableFile  # noqa: E402
from uncorrupt.column_classifier import classify_column  # noqa: E402
from uncorrupt.detector import (  # noqa: E402
    _parse_date_string,
    _reverse_gene_date,
    _serial_to_date,
    detect_file,
)

KOH = PROJECT_ROOT / "data/raw/koh_replication/files"
WALK_LEDGER = PROJECT_ROOT / "results/koh_cache_walk.jsonl"
OUT_LEDGER = PROJECT_ROOT / "results/extensive_audit.jsonl"
OUT_SUMMARY = PROJECT_ROOT / "results/extensive_audit_summary.csv"

EXCEL_EPOCH = date(1899, 12, 30)
SERIAL_MIN, SERIAL_MAX = 20000, 60000


def is_date_serial(v) -> bool:
    if isinstance(v, bool):
        return False
    if not isinstance(v, (int, float)):
        return False
    if isinstance(v, float) and not v.is_integer():
        return False
    n = int(v)
    return SERIAL_MIN <= n <= SERIAL_MAX


def has_date_signature(v) -> bool:
    """True if cell looks like a date/serial/date-string in any form we would
    want to flag if it appeared in an identifier column."""
    if isinstance(v, (date, datetime)):
        return True
    if is_date_serial(v):
        return True
    if isinstance(v, str):
        s = v.strip()
        if _parse_date_string(s):
            return True
        # Digit-string Excel-serial form (e.g. "41897") : what TSV/CSV round-
        # tripped serials look like
        if s.isdigit():
            try:
                return SERIAL_MIN <= int(s) <= SERIAL_MAX
            except ValueError:
                return False
        # Also catch the multi-token cells (e.g. "2-Oct,ATOCT2,OCT2")
        for delim in (",", ";"):
            if delim in s:
                first = s.split(delim, 1)[0].strip()
                if _parse_date_string(first):
                    return True
    return False


def is_id_float_signature(v) -> bool:
    """True if cell looks like an alphanumeric ID coerced to scientific-notation
    float (the RIKEN/2310009E13 case). Detector emits these as kind='id-float'."""
    if isinstance(v, bool):
        return False
    if isinstance(v, (int, float)) and abs(v) >= 1e10:
        return True
    if isinstance(v, str):
        s = v.strip()
        # E-notation strings that didn't survive pandas read
        import re as _re
        m = _re.fullmatch(r"^[+-]?\d*\.?\d+[eE]([+-]?\d+)$", s)
        if m and abs(int(m.group(1))) >= 10:
            return True
    return False


def decodes_to_known_corruption(v) -> bool:
    """True if a cell value matches a real Excel-gene-corruption pattern (i.e.
    `_reverse_gene_date` returns a non-empty candidate list for the decoded
    date). This is the tightest "potential FN" criterion : it excludes integer
    IDs in date-serial range that don't actually encode gene symbols (donor IDs,
    supplier IDs, sample IDs, etc.)."""
    # Date object directly
    if isinstance(v, (date, datetime)):
        d = v.date() if isinstance(v, datetime) else v
        return bool(_reverse_gene_date(d))
    # Integer in date-serial range
    if isinstance(v, (int, float)) and not isinstance(v, bool):
        if isinstance(v, float) and not v.is_integer():
            return False
        n = int(v)
        if SERIAL_MIN <= n <= SERIAL_MAX:
            d = _serial_to_date(n)
            return d is not None and bool(_reverse_gene_date(d))
    # String forms
    if isinstance(v, str):
        s = v.strip()
        # Digit-string serial
        if s.isdigit():
            try:
                n = int(s)
                if SERIAL_MIN <= n <= SERIAL_MAX:
                    d = _serial_to_date(n)
                    return d is not None and bool(_reverse_gene_date(d))
            except ValueError:
                pass
        # Parseable date-string
        parsed = _parse_date_string(s)
        if parsed:
            return any(_reverse_gene_date(d) for d in parsed)
        # Multi-token cells (first comma-separated date)
        for delim in (",", ";"):
            if delim in s:
                first = s.split(delim, 1)[0].strip()
                parsed = _parse_date_string(first)
                if parsed:
                    return any(_reverse_gene_date(d) for d in parsed)
    return False


def cell_matches_detection_kind(cell_value, kind: str) -> bool:
    """Per-kind verification : does the cell value match what the detector
    claims it found? Critical for honest precision measurement: a flag is a
    TP only if the cell type matches the kind, not just any date-signature."""
    if kind == "id-float":
        return is_id_float_signature(cell_value)
    if kind in ("gene-date", "gene-date-serial", "gene-date-string"):
        return has_date_signature(cell_value)
    # Unknown kind : fall back to either signature
    return has_date_signature(cell_value) or is_id_float_signature(cell_value)


def col_letter_to_index(letter: str) -> int:
    """Convert 'A' -> 0, 'AA' -> 26 (matches pandas 0-based column index)."""
    n = 0
    for ch in letter:
        n = n * 26 + (ord(ch.upper()) - ord("A") + 1)
    return n - 1


def _normalize_value_for_signature(v) -> "object":
    """Return a comparable form. openpyxl reads dates as `datetime`; detector
    stores them as `date`. Coerce both to `date` (no time component) for
    signature comparison."""
    if isinstance(v, datetime):
        return v.date()
    return v


def _identifier_columns_for_sheet(df) -> dict[str, "object"]:
    """Return `{stringified-column-name: original-pandas-column-key}` for every
    column that classifies as an identifier-shaped column (gene_symbol or
    riken). Some files have integer or float column headers : preserving the
    original key lets `df[original_key]` index correctly while the string form
    matches detector-emitted positions."""
    out: dict[str, object] = {}
    for col in df.columns:
        try:
            cls = classify_column(df[col], str(col))
        except Exception:
            continue
        if cls.column_type in ("gene_symbol", "riken"):
            out[str(col)] = col
    return out


def audit_file(fp: Path) -> dict:
    """For each detector flag, verify by *position* (sheet/column/row) that
    openpyxl agrees the cell has a date signature → TP. Cells where the
    detector flagged something with no date signature → real FP. For every
    cell in an identifier-shaped column that openpyxl says has a date
    signature but the detector didn't flag → potential FN.

    Position-based matching (not value-repr) avoids the `datetime` vs `date`
    repr-mismatch artefact that the original audit script suffered from.
    Column-aware FN filtering (only count gene-symbol-shaped columns) avoids
    counting legitimate date columns (`diagnosis_date`, etc.) as FN candidates.
    """
    record: dict = {
        "pmc_id": fp.parent.name,
        "file": fp.name,
        "detect_error": None,
        "openpyxl_error": None,
        "n_suspicions": 0,
        "flagged_cells": [],
        "tp_cells": [],
        "fp_cells": [],
        "missed_date_candidates": [],
        "n_identifier_columns": 0,
        "n_total_date_cells_scanned": 0,
    }

    # --- detector pass ---
    try:
        report = detect_file(str(fp))
    except UnrecoverableFile as exc:
        record["detect_error"] = f"UnrecoverableFile: {exc}"
        return record
    except Exception as exc:
        record["detect_error"] = f"{type(exc).__name__}: {exc}"
        return record
    record["n_suspicions"] = len(report.suspicions)
    record["identifier_columns"] = [str(c) for c in report.identifier_columns]
    for s in report.suspicions:
        record["flagged_cells"].append({
            "sheet": s.sheet, "column": s.column, "row": s.row,
            "value": repr(s.value), "kind": s.kind,
            "suggestion": s.suggestion, "confidence": s.confidence,
        })

    # --- pandas pass to know which columns are identifier-shaped per sheet ---
    try:
        sheets_dfs = _load_all_sheets(str(fp))
    except Exception as exc:
        record["openpyxl_error"] = f"pandas-load: {type(exc).__name__}: {exc}"
        return record
    # {sheet_name: {stringified-col-name: original-pandas-col-key}}
    id_cols_by_sheet: dict[str, dict[str, object]] = {}
    for sn, df in sheets_dfs.items():
        id_cols_by_sheet[sn] = _identifier_columns_for_sheet(df)
    record["n_identifier_columns"] = sum(len(v) for v in id_cols_by_sheet.values())

    # Build the FLAG INDEX by (sheet, column, row) position
    flagged_positions: set[tuple[str, str, int]] = {
        (f["sheet"] or "_sheet0", f["column"], f["row"])
        for f in record["flagged_cells"]
    }

    # --- check each detector flag against pandas-loaded value at that position ---
    for f in record["flagged_cells"]:
        sn = f["sheet"] or "_sheet0"
        df = sheets_dfs.get(sn)
        if df is None and len(sheets_dfs) == 1:
            df = next(iter(sheets_dfs.values()))
        # Match by stringified column name to handle integer/float headers
        col_key = None
        if df is not None:
            for c in df.columns:
                if str(c) == f["column"]:
                    col_key = c
                    break
        if df is None or col_key is None or f["row"] >= len(df):
            record["fp_cells"].append({**f, "verify_reason": "cell out of bounds"})
            continue
        cell_value = df.iloc[f["row"]][col_key]
        norm = _normalize_value_for_signature(cell_value)
        # Verify by KIND: id-float looks for large numerics, gene-date* looks
        # for date signatures. The previous version only checked dates and so
        # falsely marked every id-float as FP.
        if (cell_matches_detection_kind(norm, f["kind"])
                or isinstance(cell_value, (date, datetime))):
            record["tp_cells"].append({
                "sheet": sn, "column": f["column"], "row": f["row"],
                "value": repr(cell_value), "kind": f["kind"],
            })
        else:
            record["fp_cells"].append({
                **f, "verify_reason": (
                    f"cell value {cell_value!r} does not match kind={f['kind']!r}"
                )
            })

    # --- scan identifier-shaped columns for un-flagged date cells (FN candidates) ---
    # Tight FN criterion: the cell must decode to a KNOWN gene-corruption
    # pattern (not just any date or any integer in the date-serial range).
    # This excludes the donor-ID / sample-ID class : integers that happen to
    # fall in [20000, 60000] but don't encode a real gene symbol.
    for sn, df in sheets_dfs.items():
        id_cols = id_cols_by_sheet.get(sn, {})
        if not id_cols:
            continue
        for col_str, col_key in id_cols.items():
            col = col_str
            for idx, val in df[col_key].items():
                if val is None or (isinstance(val, float) and pd.isna(val)):
                    continue
                norm = _normalize_value_for_signature(val)
                # Cheap pre-filter: must at least look like a date signature
                if not (has_date_signature(norm)
                         or isinstance(val, (date, datetime))):
                    continue
                record["n_total_date_cells_scanned"] += 1
                if (sn, col, int(idx)) in flagged_positions:
                    continue  # already TP
                # Tight check: the cell must actually decode to a gene-corruption
                if not decodes_to_known_corruption(norm):
                    continue
                record["missed_date_candidates"].append({
                    "sheet": sn, "column": col, "row": int(idx),
                    "value": repr(val),
                })

    return record


def pick_samples(args) -> dict[str, list[Path]]:
    with WALK_LEDGER.open() as f:
        rows = [json.loads(l) for l in f]
    latest = rows[-1]
    files = latest["files"]

    pool_A = [(f["pmc_id"], f["file_name"]) for f in files
              if f["has_gene_symbols"] and f["n_suspicions"] > 0]
    pool_B = [(f["pmc_id"], f["file_name"]) for f in files
              if f["has_gene_symbols"] and f["n_suspicions"] == 0]
    pool_C = [(f["pmc_id"], f["file_name"]) for f in files
              if not f["has_gene_symbols"]]

    rng = random.Random(args.seed)
    sampled: dict[str, list[Path]] = {}
    if args.pool in (None, "A", "all"):
        # Pool A is small enough to take all
        chosen = pool_A if not args.n_dirty else rng.sample(pool_A, min(args.n_dirty, len(pool_A)))
        sampled["A"] = [KOH / pmc / fn for pmc, fn in chosen]
    if args.pool in (None, "B", "all"):
        chosen = rng.sample(pool_B, min(args.n_clean, len(pool_B)))
        sampled["B"] = [KOH / pmc / fn for pmc, fn in chosen]
    if args.pool in (None, "C", "all"):
        chosen = rng.sample(pool_C, min(args.n_empty, len(pool_C)))
        sampled["C"] = [KOH / pmc / fn for pmc, fn in chosen]
    return sampled


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--pool", choices=["A", "B", "C", "all"], default=None,
                   help="only audit one pool (default: all)")
    p.add_argument("--n-dirty", type=int, default=0,
                   help="cap pool A at this size (0 = all 278)")
    p.add_argument("--n-clean", type=int, default=100)
    p.add_argument("--n-empty", type=int, default=100)
    args = p.parse_args()

    samples = pick_samples(args)
    total = sum(len(v) for v in samples.values())
    print(f"[audit] {total} files: "
          f"A(dirty)={len(samples.get('A', []))} "
          f"B(clean)={len(samples.get('B', []))} "
          f"C(empty)={len(samples.get('C', []))}",
          file=sys.stderr)

    OUT_LEDGER.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    summary_rows: list[dict] = []
    with OUT_LEDGER.open("w") as f_ledger:
        done = 0
        for pool_name, files in samples.items():
            for fp in files:
                done += 1
                if not fp.exists():
                    print(f"[audit] {done}/{total} MISSING {fp}", file=sys.stderr)
                    continue
                rec = audit_file(fp)
                rec["pool"] = pool_name
                f_ledger.write(json.dumps(rec) + "\n")
                # Summary row
                n_flags = rec.get("n_suspicions", 0)
                n_tp = len(rec["tp_cells"])
                n_fp = len(rec["fp_cells"])
                n_potential_fn = len(rec["missed_date_candidates"])
                summary_rows.append({
                    "pool": pool_name,
                    "pmc_id": rec["pmc_id"],
                    "file": rec["file"],
                    "n_flags": n_flags,
                    "n_tp": n_tp,
                    "n_fp": n_fp,
                    "n_potential_fn": n_potential_fn,
                    "n_total_date_cells": rec["n_total_date_cells_scanned"],
                    "detect_error": rec["detect_error"] or "",
                    "openpyxl_error": rec["openpyxl_error"] or "",
                })
                if done % 25 == 0 or done == total:
                    el = time.monotonic() - t0
                    print(f"[audit] {done:>4}/{total}  ({el/60:.1f} min elapsed)",
                          file=sys.stderr, flush=True)

    # Write CSV summary
    OUT_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with OUT_SUMMARY.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        w.writeheader()
        w.writerows(summary_rows)

    # Aggregate report
    elapsed = time.monotonic() - t0
    print(f"\n=== DONE ({elapsed/60:.1f} min) ===", file=sys.stderr)
    for pool in ("A", "B", "C"):
        pool_rows = [r for r in summary_rows if r["pool"] == pool]
        if not pool_rows:
            continue
        total_flags = sum(r["n_flags"] for r in pool_rows)
        total_tp = sum(r["n_tp"] for r in pool_rows)
        total_fp = sum(r["n_fp"] for r in pool_rows)
        total_fn = sum(r["n_potential_fn"] for r in pool_rows)
        precision = total_tp / total_flags * 100 if total_flags else 0.0
        files_w_fn = sum(1 for r in pool_rows if r["n_potential_fn"] > 0)
        print(f"\nPool {pool}: {len(pool_rows)} files", file=sys.stderr)
        print(f"  flags: {total_flags:,}", file=sys.stderr)
        print(f"  cross-verified TP: {total_tp:,}  ({precision:.1f}% precision)", file=sys.stderr)
        print(f"  unmatched flags (rare): {total_fp:,}", file=sys.stderr)
        print(f"  potential missed dates: {total_fn:,}  (across {files_w_fn} files)", file=sys.stderr)
    print(f"\nLedger: {OUT_LEDGER}", file=sys.stderr)
    print(f"Summary: {OUT_SUMMARY}", file=sys.stderr)


if __name__ == "__main__":
    main()
