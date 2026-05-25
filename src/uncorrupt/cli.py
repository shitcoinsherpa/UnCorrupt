"""Console entry-point for UnCorrupt.

Researchers invoke this directly:

    uncorrupt detect path/to/spreadsheet.xlsx
    uncorrupt schema path/to/spreadsheet.xlsx
    uncorrupt audit  path/to/folder/

Designed to be self-sufficient : once installed, no further docs or
human-in-the-loop required. The `detect` command emits a structured
report; `schema` emits a Frictionless Table Schema sidecar that
prevents future corruption when the file is re-opened in a different
locale or Excel version; `audit` walks a folder and produces a
summary report.

Exit codes:
   0 : nothing flagged
   1 : at least one high-confidence flag (conf >= 0.85); block the build
   2 : medium-confidence flags only (0.50 <= conf < 0.85); warn
   3 : file unreadable or other error
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .detector import detect_file
from .schema_export import write_schema_sidecar


def _cmd_detect(args) -> int:
    path = Path(args.file)
    if not path.exists():
        print(f"ERROR: file not found: {path}", file=sys.stderr)
        return 3
    try:
        report = detect_file(str(path), row_context_boost=not args.no_boost)
    except Exception as exc:
        print(f"ERROR: detect_file failed: {type(exc).__name__}: {exc}",
              file=sys.stderr)
        return 3

    if args.json:
        out = {
            "file": str(path),
            "rows_scanned": report.rows_scanned,
            "columns_scanned": report.columns_scanned,
            "identifier_columns": report.identifier_columns,
            "suspicions": [
                {
                    "sheet": s.sheet, "column": s.column, "row": s.row,
                    "value": str(s.value), "kind": s.kind,
                    "suggestion": s.suggestion, "confidence": s.confidence,
                    "reason": s.reason, "xref_status": s.xref_status,
                }
                for s in report.suspicions
            ],
        }
        json.dump(out, sys.stdout, indent=2)
        print()
    else:
        # Plain-English labels for the internal `kind` enum so the CLI output
        # reads as something a bench scientist understands, not an enum dump.
        # Keep in sync with src/uncorrupt/app.py:_KIND_LABELS.
        KIND_LABELS = {
            "gene-date":                "Date mistaken for gene name",
            "gene-date-string":         "Date text in gene column ('2-Sep')",
            "gene-date-serial":         "Excel-serial integer in gene column",
            "id-float":                 "Float (gene ID lost to scientific notation)",
            "long-int-precision-loss":  "Integer too big for Excel (trailing digits lost)",
            "leading-zero-stripped":    "Leading zero stripped",
            "autofill-sequence":        "Autofill drag (Excel filled a sequence)",
            "cas-registry":             "Chemical CAS number Excel read as a date",
            "decimal-comma":            "Decimal-comma locale collision",
            "time-coercion":            "Time string (plate well IDs)",
            "homoglyph":                "Unicode look-alike (Cyrillic/Greek/fullwidth)",
            "unrecognized-symbol":      "Symbol not in any registry",
            "file-too-large":           "File too large to scan safely",
        }
        # Confidence bands must match src/uncorrupt/app.py:_confidence_label.
        BANDS = [
            (0.85, "High",   "what shows is virtually always real corruption"),
            (0.50, "Medium", "review before accepting; mostly correct"),
            (0.30, "Low",    "isolated flags; pattern not corroborated"),
            (0.00, "Info",   "informational only (suppressed by default)"),
        ]
        print(f"File: {path}")
        print(f"Rows scanned: {report.rows_scanned}")
        print(f"Columns scanned: {report.columns_scanned}")
        print(f"Identifier columns: {len(report.identifier_columns)}")
        print(f"Suspicions: {len(report.suspicions)}")

        # Bucket by confidence band
        from .detector import Suspicion
        bucket: dict[str, list[Suspicion]] = {label: [] for _, label, _ in BANDS}
        for s in report.suspicions:
            for threshold, label, _ in BANDS:
                if s.confidence >= threshold:
                    bucket[label].append(s)
                    break

        for threshold, label, hint in BANDS:
            flags = bucket[label]
            if not flags:
                continue
            print(f"\n--- {label} confidence ({len(flags)} flags) -- {hint}")
            for s in flags[: args.max_per_band]:
                loc = f"{s.sheet}!{s.column}" if s.sheet else s.column
                kind_label = KIND_LABELS.get(s.kind, s.kind)
                print(f"  {kind_label}  {loc!r} row {s.row}: {s.value!r}")
                if s.suggestion:
                    print(f"    suggest: {s.suggestion}")
            if len(flags) > args.max_per_band:
                print(f"    ... and {len(flags) - args.max_per_band} more")

    # Exit code by highest band (bands match UI: High >= 0.85, Medium >= 0.50).
    high = sum(1 for s in report.suspicions if s.confidence >= 0.85)
    mid = sum(1 for s in report.suspicions if 0.30 <= s.confidence < 0.85)
    if high:
        return 1
    if mid:
        return 2
    return 0


def _cmd_schema(args) -> int:
    path = Path(args.file)
    if not path.exists():
        print(f"ERROR: file not found: {path}", file=sys.stderr)
        return 3
    try:
        report = detect_file(str(path), row_context_boost=False)
    except Exception as exc:
        print(f"ERROR: detect_file failed: {type(exc).__name__}: {exc}",
              file=sys.stderr)
        return 3

    # Load the underlying df(s) and pick the first/largest sheet
    from .app import _load_all_sheets
    sheets = _load_all_sheets(str(path))
    # Schema is emitted per sheet; for multi-sheet files, write one sidecar
    # per non-trivial sheet.
    out_dir = Path(args.output) if args.output else path.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for sheet_name, df in sheets.items():
        if df.empty: continue
        suffix = f".{sheet_name}" if len(sheets) > 1 else ""
        sidecar = out_dir / f"{path.stem}{suffix}.schema.json"
        write_schema_sidecar(df, report, sidecar,
                             name=f"{path.stem}{suffix}")
        written.append(sidecar)
    print(f"Wrote {len(written)} schema sidecar(s):")
    for p in written:
        print(f"  {p}")
    print()
    print("# How to use the sidecar to prevent re-corruption")
    print("# 1. Validate any future copy of this file against the schema:")
    print(f"#       frictionless validate --schema {written[0].name} {path.name}")
    print("# 2. Or programmatically (Python):")
    print("#       from frictionless import Schema, Resource")
    print(f"#       Schema('{written[0].name}').validate(Resource('{path.name}'))")
    print("# 3. The schema pins identifier columns to type='string' with regex")
    print("#    constraints : Excel re-import or locale-aware tools that")
    print("#    coerce will fail validation instead of silently corrupting.")
    return 0


def _cmd_audit(args) -> int:
    root = Path(args.folder)
    if not root.exists() or not root.is_dir():
        print(f"ERROR: folder not found: {root}", file=sys.stderr)
        return 3
    pattern = "**/*.xlsx" if args.recursive else "*.xlsx"
    files = list(root.glob(pattern))
    if args.recursive:
        files += list(root.glob("**/*.xls"))
    else:
        files += list(root.glob("*.xls"))
    print(f"Auditing {len(files)} file(s) under {root}")
    # Plain-English kind labels (mirror detect / app). Coverage of every
    # detector `kind=` is enforced by tests/test_kind_labels.py; a missing
    # entry raises KeyError below instead of leaking the raw enum.
    KIND_LABELS = {
        "gene-date":                "date-mistaken-for-gene",
        "gene-date-string":         "date-text-in-gene-column",
        "gene-date-serial":         "excel-serial-in-gene-column",
        "id-float":                 "id-lost-to-scientific-notation",
        "long-int-precision-loss":  "integer-too-big-for-excel",
        "leading-zero-stripped":    "leading-zero-stripped",
        "autofill-sequence":        "autofill-drag",
        "cas-registry":             "cas-number-read-as-date",
        "decimal-comma":            "decimal-comma-collision",
        "time-coercion":            "time-coerced-id",
        "homoglyph":                "unicode-look-alike",
        "unrecognized-symbol":      "symbol-not-in-registry",
        "file-too-large":           "file-too-large",
    }
    summary: dict[str, dict] = {}
    n_high_total = 0
    n_corrupted_files = 0
    for f in files:
        try:
            report = detect_file(str(f), row_context_boost=not args.no_boost)
        except Exception as exc:
            summary[str(f)] = {"error": f"{type(exc).__name__}: {exc!s}"}
            continue
        # Band thresholds match the rest of the tool: High >= 0.85.
        high = [s for s in report.suspicions if s.confidence >= 0.85]
        mid = [s for s in report.suspicions if 0.50 <= s.confidence < 0.85]
        low = [s for s in report.suspicions if 0.30 <= s.confidence < 0.50]
        summary[str(f)] = {
            "n_suspicions": len(report.suspicions),
            "n_high_confidence": len(high),
            "n_mid_confidence": len(mid),
            "n_low_confidence": len(low),
            "kinds": sorted({s.kind for s in report.suspicions}),
        }
        n_high_total += len(high)
        if high or mid:
            n_corrupted_files += 1
    if args.json:
        json.dump(summary, sys.stdout, indent=2)
        print()
    else:
        print(f"\nFiles with corruption flags: {n_corrupted_files}/{len(files)}")
        for fpath, s in summary.items():
            if "error" in s:
                print(f"  [error] {Path(fpath).name}: {s['error']}")
            elif s.get("n_high_confidence") or s.get("n_mid_confidence"):
                kinds = ", ".join(KIND_LABELS.get(k, k) for k in s["kinds"])
                print(f"  {Path(fpath).name}: "
                      f"{s['n_high_confidence']} high + "
                      f"{s['n_mid_confidence']} mid "
                      f"({kinds})")
    # Exit-code semantics (CI-friendly):
    #   0 : nothing flagged at medium-or-above
    #   1 : at least one high-confidence flag (block the build)
    #   2 : only medium-confidence flags (warn)
    if n_high_total > 0:
        return 1
    if n_corrupted_files > 0:
        return 2
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="uncorrupt",
        description="Detect and prevent Excel gene-symbol corruption.",
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_detect = sub.add_parser("detect",
                              help="Scan a spreadsheet for corruption.")
    p_detect.add_argument("file")
    p_detect.add_argument("--json", action="store_true",
                           help="JSON output instead of human-readable.")
    p_detect.add_argument("--no-boost", action="store_true",
                           help="Skip row-xref boost (faster but noisier).")
    p_detect.add_argument("--max-per-band", type=int, default=20,
                           help="Cap on flags shown per confidence band.")
    p_detect.set_defaults(func=_cmd_detect)

    p_schema = sub.add_parser("schema",
                               help="Emit a Frictionless Table Schema sidecar.")
    p_schema.add_argument("file")
    p_schema.add_argument("-o", "--output", default=None,
                           help="Output directory (default: alongside the file).")
    p_schema.set_defaults(func=_cmd_schema)

    p_audit = sub.add_parser("audit",
                              help="Scan a folder of spreadsheets.")
    p_audit.add_argument("folder")
    p_audit.add_argument("-r", "--recursive", action="store_true")
    p_audit.add_argument("--json", action="store_true")
    p_audit.add_argument("--no-boost", action="store_true")
    p_audit.set_defaults(func=_cmd_audit)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
