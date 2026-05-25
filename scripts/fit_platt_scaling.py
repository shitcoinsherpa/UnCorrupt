"""Fit Platt scaling on the FP-anchor labeled calibration corpus.

Methodology (per the literature review): Platt 1999 logistic regression
P(y=1 | s) = sigmoid(coef * s + intercept), fit on the corroborated +
contradicted subset. This is the right calibration method below the
~1000-sample isotonic threshold (Niculescu-Mizil & Caruana 2005).

Inputs: results/calibration_fp_anchors.jsonl
Outputs: src/uncorrupt/calibration_params.json  (coef, intercept, n_samples)
         results/platt_scaling_report.md         (human-readable diagnostic)

The detector loads `calibration_params.json` at import time; if present,
the Suspicion dataclass carries a `calibrated_probability` alongside the
legacy `confidence`. The legacy value stays for backward compatibility.

Run:
   python scripts/fit_platt_scaling.py
"""
from __future__ import annotations

import json
import math
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "results/calibration_fp_anchors.jsonl"
OUT_PARAMS = ROOT / "src/uncorrupt/calibration_params.json"
OUT_REPORT = ROOT / "results/platt_scaling_report.md"

EXCLUDED = {"decimal-comma", "cas-registry", "unrecognized-symbol"}


def _load_hgnc_resolver():
    sys.path.insert(0, str(ROOT / "src"))
    from uncorrupt.corpus import load_hgnc
    h = load_hgnc()

    def expand(suggestions: set[str]) -> set[str]:
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
    return expand


def _label_records():
    """Yield (score, outcome) pairs from the ledger, excluding info-only kinds."""
    expand = _load_hgnc_resolver()
    for line in LEDGER.open():
        r = json.loads(line)
        if r["kind"] in EXCLUDED:
            continue
        if not r.get("suggestion"):
            continue
        suggestions = {p.strip() for p in r["suggestion"].split("|") if p.strip()}
        match_set = expand(suggestions)
        # Case-insensitive comparison : cross-species ortholog symbols
        # are mixed-case (mouse Marchf10, fly mei-9) while HGNC
        # suggestions are uppercase. Without case-folding, cross-species
        # corroborations get mis-labeled as contradicted.
        match_set_lower = {s.lower() for s in match_set}
        resolved = [g for _, g in r["resolved_in_row"]]
        if not resolved:
            continue  # inconclusive : excluded from fit
        is_corroborated = any(g.lower() in match_set_lower for g in resolved)
        yield r["confidence"], (1 if is_corroborated else 0)


def main() -> int:
    if not LEDGER.exists():
        print(f"FATAL: {LEDGER} missing", file=sys.stderr)
        return 1

    scores, outcomes = [], []
    for s, o in _label_records():
        scores.append(float(s))
        outcomes.append(int(o))

    n = len(scores)
    n_pos = sum(outcomes)
    n_neg = n - n_pos
    print(f"Labeled samples: {n} ({n_pos} positive, {n_neg} negative)")

    if n < 5:
        print("Too few samples to fit. Run a fuller calibration walk first.", file=sys.stderr)
        return 2
    if n_pos == 0 or n_neg == 0:
        print("Only one class observed : Platt scaling degenerate.", file=sys.stderr)
        return 3

    try:
        import numpy as np
        from sklearn.linear_model import LogisticRegression
    except ImportError:
        print("scikit-learn + numpy required: pip install scikit-learn",
              file=sys.stderr)
        return 4

    Xs = np.array(scores).reshape(-1, 1)
    ys = np.array(outcomes)
    lr = LogisticRegression(C=1e6)  # near-MLE
    lr.fit(Xs, ys)
    intercept = float(lr.intercept_[0])
    coef = float(lr.coef_[0][0])
    print(f"Platt: P(y=1|s) = sigmoid({coef:.4f} * s + {intercept:.4f})")

    # Per-bin recalibrated values
    bins = defaultdict(list)
    for s, o in zip(scores, outcomes):
        bins[s].append(o)
    per_bin = []
    for s in sorted(bins):
        outs = bins[s]
        n_b = len(outs)
        n_correct = sum(outs)
        empirical = n_correct / n_b
        platt_p = 1.0 / (1.0 + math.exp(-(coef * s + intercept)))
        per_bin.append({
            "stated": s, "n": n_b, "correct": n_correct,
            "empirical": empirical, "platt": platt_p,
        })

    # Brier on fit
    brier_raw = sum((s - o) ** 2 for s, o in zip(scores, outcomes)) / n
    brier_platt = sum(
        (1.0 / (1.0 + math.exp(-(coef * s + intercept))) - o) ** 2
        for s, o in zip(scores, outcomes)
    ) / n

    params = {
        "method": "platt-scaling",
        "intercept": intercept,
        "coef": coef,
        "n_samples": n,
        "n_positive": n_pos,
        "n_negative": n_neg,
        "brier_raw": brier_raw,
        "brier_platt": brier_platt,
        "fit_date_iso": "2026-05-19",
        "reference": (
            "Platt 1999, Probabilistic outputs for SVMs; "
            "Niculescu-Mizil & Caruana 2005, ICML."
        ),
        "per_bin": per_bin,
    }
    OUT_PARAMS.parent.mkdir(parents=True, exist_ok=True)
    OUT_PARAMS.write_text(json.dumps(params, indent=2))
    print(f"Wrote {OUT_PARAMS}")

    lines = [
        "# Platt scaling calibration",
        "",
        f"Labeled samples: {n} ({n_pos} positive, {n_neg} negative).",
        f"Platt fit: P(y=1|s) = sigmoid({coef:.4f} * s + {intercept:.4f}).",
        "",
        "## Per-bin recalibration",
        "",
        "| Stated | n | correct | Empirical | Platt-calibrated |",
        "|--------|---|---------|-----------|-------------------|",
    ]
    for b in per_bin:
        lines.append(
            f"| {b['stated']:.2f} | {b['n']} | {b['correct']} | "
            f"{b['empirical']:.3f} | {b['platt']:.3f} |"
        )
    lines += [
        "",
        f"## Brier scores",
        f"- raw (stated confidence): {brier_raw:.4f}",
        f"- Platt-calibrated: {brier_platt:.4f}",
        f"- improvement: {(brier_raw - brier_platt):.4f}",
        "",
        "## Use",
        "",
        "The detector loads `src/uncorrupt/calibration_params.json` at import",
        "time if present. Each Suspicion's `calibrated_probability` field is",
        "set to the Platt-mapped value. Legacy `confidence` stays for",
        "backwards compatibility.",
    ]
    OUT_REPORT.write_text("\n".join(lines))
    print(f"Wrote {OUT_REPORT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
