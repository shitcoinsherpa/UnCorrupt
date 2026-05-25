"""Automatically adjudicate inconclusive flags using the multi-species
xref directly (NOT just same-row external IDs).

The vanilla calibration walk labels a flag "inconclusive" when the row
has no resolvable external ID. But many cells in those rows ARE the
gene symbol directly (e.g., a string cell with "BRCA1" or "Sept2"). The
multi-species xref now indexes >1.6M symbols : we can check whether the
row CONTAINS the suggested gene name as a string, not just as an
external ID.

Three labels:
  - auto_corroborated: another cell in the row IS the suggested gene
    symbol (case-insensitive, in any species). Likely TP.
  - auto_contradicted: another cell in the row IS a DIFFERENT gene
    symbol. Likely FP (the row is about that other gene).
  - still_inconclusive: no gene-symbol-shape strings in any other column.

Output: results/auto_adjudicated.jsonl with the new labels merged.
       results/calibration_ci_auto_report.md with Wilson + BCa on the
       combined xref-labeled + auto-labeled set.

This is "weak supervision" applied to labeled-set expansion. Per the
active-learning guidance:
auto-labels are noisy but cheaper than human review; pair with the
margin-sampling queue for the hardest cases.

Run:
    python scripts/auto_adjudicate_inconclusives.py
"""
from __future__ import annotations

import json
import math
import sys
import warnings
from collections import defaultdict
from pathlib import Path

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

LEDGER = ROOT / "results/calibration_fp_anchors.jsonl"
OUT_LEDGER = ROOT / "results/auto_adjudicated.jsonl"
OUT_REPORT = ROOT / "results/calibration_ci_auto_report.md"

EXCLUDED = {"decimal-comma", "cas-registry", "unrecognized-symbol"}


def _wilson_ci(k: int, n: int) -> tuple[float, float]:
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    z = 1.96
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def _bca(outcomes: list[int]) -> tuple[float, float]:
    n = len(outcomes)
    if n == 0:
        return (0.0, 1.0)
    try:
        import numpy as np
        from scipy.stats import bootstrap
    except ImportError:
        return _wilson_ci(sum(outcomes), n)
    arr = np.asarray(outcomes, dtype=float)
    if arr.sum() == 0 or arr.sum() == n:
        return _wilson_ci(int(arr.sum()), n)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = bootstrap((arr,), np.mean, n_resamples=10_000,
                         confidence_level=0.95, method="BCa",
                         random_state=20260520)
    return (float(res.confidence_interval.low),
            float(res.confidence_interval.high))


def _load_symbol_set() -> frozenset[str]:
    """Lower-cased set of every known gene symbol across all species
    registries + HGNC current symbols."""
    syms: set[str] = set()
    reg_dir = ROOT / "data/raw/registries"
    for p in reg_dir.glob("*.xref.tsv"):
        with p.open() as fh:
            next(fh, None)  # header
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                _sp, _ext, cur, alias = parts
                if cur:
                    syms.add(cur.strip().lower())
                if alias:
                    syms.add(alias.strip().lower())
    # HGNC current symbols
    from uncorrupt.corpus import load_hgnc
    h = load_hgnc()
    for s in h.current_symbols:
        syms.add(s.lower())
    for s in h.prev_symbol_to_current:
        syms.add(s.lower())
    for s in h.alias_to_current:
        syms.add(s.lower())
    return frozenset(syms)


def _expand_with_hgnc(suggestions: set[str]) -> set[str]:
    from uncorrupt.corpus import load_hgnc
    h = load_hgnc()
    expanded = set(suggestions)
    for s in list(suggestions):
        canonical = h.resolve(s)
        if canonical:
            expanded.add(canonical)
        for old, new in h.prev_symbol_to_current.items():
            if new == s or new == canonical:
                expanded.add(old)
        for alias, new in h.alias_to_current.items():
            if new == s or new == canonical:
                expanded.add(alias)
    return expanded


def main() -> int:
    if not LEDGER.exists():
        print(f"FATAL: {LEDGER} missing", file=sys.stderr)
        return 1

    print("Loading symbol set...")
    symbol_set = _load_symbol_set()
    print(f"  {len(symbol_set):,} symbols indexed")

    # Re-open every row, scan the source file's row, look for any string
    # cell that IS a gene symbol.
    sys.path.insert(0, str(ROOT / "src"))
    from uncorrupt.app import _load_all_sheets
    from pathlib import Path as _Path
    KOH = ROOT / "data/raw/koh_replication/files"

    inconclusive_records = []
    other_records = []
    for line in LEDGER.open():
        r = json.loads(line)
        if r["kind"] in EXCLUDED or not r.get("suggestion"):
            continue
        # Use the same case-insensitive match the bootstrap script uses
        suggestions = {p.strip() for p in r["suggestion"].split("|") if p.strip()}
        ms = _expand_with_hgnc(suggestions)
        msl = {s.lower() for s in ms}
        resolved = [g for _, g in r.get("resolved_in_row") or []]
        if resolved:
            is_corr = any(g.lower() in msl for g in resolved)
            r["_xref_label"] = "corroborated" if is_corr else "contradicted"
            other_records.append(r)
        else:
            r["_xref_label"] = "inconclusive"
            inconclusive_records.append(r)

    print(f"\nInconclusive flags to auto-adjudicate: {len(inconclusive_records)}")
    auto_corr = 0
    auto_contr = 0
    still_inc = 0
    adjudicated = []

    file_sheets_cache: dict[str, dict] = {}
    for i, r in enumerate(inconclusive_records):
        pmc = r["pmc_id"]
        fname = r["file"]
        path = KOH / pmc / fname
        if not path.exists():
            still_inc += 1
            continue
        key = str(path)
        try:
            if key not in file_sheets_cache:
                file_sheets_cache[key] = _load_all_sheets(key)
            sheets = file_sheets_cache[key]
        except Exception:
            still_inc += 1
            continue
        sheet_name = r["sheet"]
        df = sheets.get(sheet_name)
        if df is None and len(sheets) == 1:
            df = next(iter(sheets.values()))
        if df is None or r["row"] not in df.index:
            still_inc += 1
            continue
        row = df.loc[r["row"]]
        suggestions = {p.strip() for p in r["suggestion"].split("|") if p.strip()}
        ms = _expand_with_hgnc(suggestions)
        msl = {s.lower() for s in ms}
        other_genes_lower: set[str] = set()
        suggestion_match = False
        for col, val in row.items():
            if str(col) == r["column"] or not isinstance(val, str):
                continue
            v = val.strip().lower()
            if not (2 <= len(v) <= 20):
                continue
            if v in msl:
                suggestion_match = True
            elif v in symbol_set:
                other_genes_lower.add(v)
        if suggestion_match:
            r["_auto_label"] = "auto_corroborated"
            auto_corr += 1
        elif other_genes_lower:
            r["_auto_label"] = "auto_contradicted"
            auto_contr += 1
        else:
            r["_auto_label"] = "still_inconclusive"
            still_inc += 1
        adjudicated.append(r)
        if (i + 1) % 100 == 0:
            print(f"  {i+1}/{len(inconclusive_records)} auto_corr={auto_corr} "
                  f"auto_contr={auto_contr} still_inc={still_inc}")
            # Free cache periodically to avoid memory bloat
            if len(file_sheets_cache) > 200:
                file_sheets_cache.clear()

    print(f"\nAuto-adjudicated {len(adjudicated)}:")
    print(f"  auto_corroborated: {auto_corr}")
    print(f"  auto_contradicted: {auto_contr}")
    print(f"  still_inconclusive: {still_inc}")

    OUT_LEDGER.write_text(
        "\n".join(json.dumps(r) for r in (other_records + adjudicated)) + "\n"
    )

    # Combined CI
    outcomes_xref: list[int] = []
    outcomes_combined: list[int] = []
    for r in other_records:
        v = 1 if r["_xref_label"] == "corroborated" else 0
        outcomes_xref.append(v)
        outcomes_combined.append(v)
    for r in adjudicated:
        if r["_auto_label"] == "still_inconclusive":
            continue
        v = 1 if r["_auto_label"] == "auto_corroborated" else 0
        outcomes_combined.append(v)

    n_xref = len(outcomes_xref)
    n_combined = len(outcomes_combined)
    p_xref = sum(outcomes_xref) / n_xref if n_xref else 0.0
    p_combined = sum(outcomes_combined) / n_combined if n_combined else 0.0
    w_xref = _wilson_ci(sum(outcomes_xref), n_xref)
    w_combined = _wilson_ci(sum(outcomes_combined), n_combined)
    b_combined = _bca(outcomes_combined)

    lines = [
        "# Calibration CI with auto-adjudicated inconclusives",
        "",
        f"Xref-labeled only: {n_xref} cells, precision {p_xref:.3f} "
        f"Wilson [{w_xref[0]:.3f}, {w_xref[1]:.3f}]",
        "",
        f"Combined (xref + auto): **{n_combined} cells**, "
        f"**precision {p_combined:.3f}** "
        f"Wilson [{w_combined[0]:.3f}, {w_combined[1]:.3f}] "
        f"BCa [{b_combined[0]:.3f}, {b_combined[1]:.3f}]",
        "",
        "## Methodology",
        "",
        "Auto-adjudication checks every other cell in the same row of the",
        "source xlsx file. If any cell IS a known gene symbol (across all",
        "species: HGNC + MGI + ZFIN + FlyBase + WormBase + rat, ~1.8M",
        "symbols total), the flag is labeled `auto_corroborated`/`auto_contradicted`",
        "based on whether that symbol matches the suggestion. Cells with no",
        "in-row symbol stay inconclusive.",
        "",
        "This is weak supervision: faster than human review, noisier than",
        "manual adjudication. Use as a labeled-set expansion lever, NOT as",
        "the primary precision claim. Sample the auto-labeled set for",
        "human verification before publication.",
    ]
    OUT_REPORT.write_text("\n".join(lines))
    print(f"\nWrote {OUT_REPORT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
