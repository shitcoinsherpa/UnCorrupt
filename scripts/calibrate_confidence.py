"""Compute per-confidence-bin precision, ECE, Brier score, and reliability
diagram for the UnCorrupt detector.

Methodology follows the literature review (Brier 1950 / Murphy 1973; Naeini
et al. 2015 ECE; Guo et al. 2017 reliability diagrams). Calibration of a
rule-based detector is framed as Laplace-smoothed per-bin empirical
precision with Wilson 95% CIs : the bins are the discrete confidence
levels the detector emits, not learned probabilities. Platt scaling
(scikit-learn) is reported as a cross-check.

Input: results/extensive_audit.jsonl  (cell-level evidence, ~480 files
   stratified across "with corruption" / "clean with genes" /
   "no genes" pools : each flagged cell has a TP/FP verdict from an
   independent openpyxl pass.)

Output:
   results/calibration_report.json     (machine-readable stats)
   results/calibration_report.md       (human-readable summary)
   results/calibration_reliability.png (5-bin reliability diagram)

Run:
   python scripts/calibrate_confidence.py
"""
from __future__ import annotations

import json
import math
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "results/extensive_audit.jsonl"
OUT_JSON = ROOT / "results/calibration_report.json"
OUT_MD = ROOT / "results/calibration_report.md"
OUT_PNG = ROOT / "results/calibration_reliability.png"


def _wilson_ci(n_correct: int, n_total: int, alpha: float = 0.05) -> tuple[float, float]:
    """Wilson 95% CI for a binomial proportion. Robust at small n and
    extreme proportions (precision ~1.0) where the normal approximation
    breaks. From Wilson 1927."""
    if n_total == 0:
        return (0.0, 1.0)
    p = n_correct / n_total
    z = 1.96  # alpha=0.05 two-sided
    denom = 1 + z * z / n_total
    centre = (p + z * z / (2 * n_total)) / denom
    half = z * math.sqrt(p * (1 - p) / n_total + z * z / (4 * n_total * n_total)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def _brier(scores: list[float], outcomes: list[int]) -> float:
    """BS = (1/N) Σ (f_i − o_i)². Brier 1950."""
    n = len(scores)
    if n == 0:
        return 0.0
    return sum((s - o) ** 2 for s, o in zip(scores, outcomes)) / n


def _brier_decomposition(
    scores: list[float], outcomes: list[int]
) -> dict[str, float]:
    """Murphy 1973: BS = Reliability − Resolution + Uncertainty.

    Reliability is the calibration-error component (lower better).
    Resolution is the discriminative power (higher better).
    Uncertainty is the irreducible entropy of the label distribution.
    """
    n = len(scores)
    if n == 0:
        return {"reliability": 0.0, "resolution": 0.0, "uncertainty": 0.0}
    o_bar = sum(outcomes) / n
    bins: dict[float, list[int]] = defaultdict(list)
    for s, o in zip(scores, outcomes):
        bins[s].append(o)
    reliability = 0.0
    resolution = 0.0
    for s, outs in bins.items():
        n_k = len(outs)
        o_bar_k = sum(outs) / n_k
        reliability += n_k * (s - o_bar_k) ** 2
        resolution += n_k * (o_bar_k - o_bar) ** 2
    reliability /= n
    resolution /= n
    uncertainty = o_bar * (1 - o_bar)
    return {
        "reliability": reliability,
        "resolution": resolution,
        "uncertainty": uncertainty,
    }


def _ece(scores: list[float], outcomes: list[int]) -> float:
    """Expected Calibration Error: weighted |confidence − accuracy| per bin.
    Naeini et al. 2015 / Guo et al. 2017. For a discrete-output detector
    each unique confidence is its own bin."""
    n = len(scores)
    if n == 0:
        return 0.0
    bins: dict[float, list[int]] = defaultdict(list)
    for s, o in zip(scores, outcomes):
        bins[s].append(o)
    ece = 0.0
    for conf, outs in bins.items():
        n_k = len(outs)
        acc_k = sum(outs) / n_k
        ece += (n_k / n) * abs(conf - acc_k)
    return ece


def _laplace_posterior(
    n_correct: int, n_total: int, alpha: float = 1.0, beta: float = 1.0
) -> float:
    """Beta(α,β) prior + Binomial → Beta(α+n_correct, β+n_fail) posterior
    mean. With α=β=1 this is Laplace's rule of succession."""
    return (n_correct + alpha) / (n_total + alpha + beta)


def main() -> int:
    if not LEDGER.exists():
        print(f"FATAL: ledger not found at {LEDGER}", file=sys.stderr)
        return 1
    scores: list[float] = []
    outcomes: list[int] = []
    n_files = 0
    n_failures = 0
    n_skipped_no_conf = 0
    pool_counts: dict[str, int] = defaultdict(int)
    kind_counts: dict[str, int] = defaultdict(int)
    with LEDGER.open() as fh:
        for line in fh:
            record = json.loads(line)
            n_files += 1
            if record.get("detect_error") or record.get("openpyxl_error"):
                n_failures += 1
                continue
            pool_counts[record.get("pool", "?")] += 1
            tp_set = {
                (c.get("sheet"), c.get("column"), c.get("row"))
                for c in record.get("tp_cells", [])
            }
            for cell in record.get("flagged_cells", []):
                conf = cell.get("confidence")
                if conf is None:
                    n_skipped_no_conf += 1
                    continue
                key = (cell.get("sheet"), cell.get("column"), cell.get("row"))
                is_tp = key in tp_set
                scores.append(float(conf))
                outcomes.append(1 if is_tp else 0)
                kind_counts[cell.get("kind", "?")] += 1

    if not scores:
        print(
            "No scored cells found. Check that ledger contains flagged_cells "
            "with confidence values.",
            file=sys.stderr,
        )
        return 2

    # Per-bin precision
    per_bin: dict[float, dict] = {}
    bins: dict[float, list[int]] = defaultdict(list)
    for s, o in zip(scores, outcomes):
        bins[s].append(o)
    for conf in sorted(bins):
        outs = bins[conf]
        n_total = len(outs)
        n_correct = sum(outs)
        precision = n_correct / n_total
        lo, hi = _wilson_ci(n_correct, n_total)
        laplace = _laplace_posterior(n_correct, n_total)
        per_bin[conf] = {
            "n_total": n_total,
            "n_correct": n_correct,
            "precision": precision,
            "wilson_ci_low": lo,
            "wilson_ci_high": hi,
            "laplace_posterior": laplace,
            "deviation_from_stated": precision - conf,
        }

    bs = _brier(scores, outcomes)
    decomp = _brier_decomposition(scores, outcomes)
    ece = _ece(scores, outcomes)

    # Platt cross-check (only if scikit-learn is available)
    platt_params = None
    try:
        from sklearn.linear_model import LogisticRegression
        import numpy as np
        Xs = np.array(scores).reshape(-1, 1)
        ys = np.array(outcomes)
        if len(set(ys)) > 1: # avoid LR on single-class
            lr = LogisticRegression(C=1e6)
            lr.fit(Xs, ys)
            platt_params = {
                "intercept": float(lr.intercept_[0]),
                "coef": float(lr.coef_[0][0]),
                "note": "P(y=1|s) = sigmoid(coef*s + intercept)",
            }
    except ImportError:
        pass

    summary = {
        "ledger": str(LEDGER),
        "n_files_total": n_files,
        "n_files_failed": n_failures,
        "n_scored_cells": len(scores),
        "n_cells_skipped_no_conf": n_skipped_no_conf,
        "pool_counts": dict(pool_counts),
        "kind_counts": dict(kind_counts),
        "per_bin_precision": per_bin,
        "brier_score": bs,
        "brier_decomposition": decomp,
        "expected_calibration_error": ece,
        "platt_cross_check": platt_params,
    }

    OUT_JSON.write_text(json.dumps(summary, indent=2))
    print(f"Wrote {OUT_JSON}")

    # Markdown report
    lines = [
        "# UnCorrupt confidence calibration",
        "",
        f"Source: `{LEDGER.name}` ({n_files} files, {n_failures} failed to "
        f"load, {len(scores)} scored cells).",
        "",
        "## Per-confidence-bin empirical precision",
        "",
        "| Stated | n_total | n_correct | Precision | Wilson 95% CI | Laplace posterior | Δ (precision − stated) |",
        "|--------|---------|-----------|-----------|---------------|--------------------|-------------------------|",
    ]
    for conf in sorted(per_bin):
        b = per_bin[conf]
        delta = b["deviation_from_stated"]
        delta_marker = "✗" if abs(delta) > 0.05 else "✓"
        lines.append(
            f"| {conf:.2f} | {b['n_total']} | {b['n_correct']} | "
            f"{b['precision']:.3f} | "
            f"[{b['wilson_ci_low']:.3f}, {b['wilson_ci_high']:.3f}] | "
            f"{b['laplace_posterior']:.3f} | "
            f"{delta:+.3f} {delta_marker} |"
        )
    lines += [
        "",
        "## Scalar calibration metrics",
        "",
        f"- **Brier score**: {bs:.4f} "
        f"(0 = perfect; uncertainty floor for this corpus is "
        f"{decomp['uncertainty']:.4f})",
        f"  - Reliability: {decomp['reliability']:.4f} (lower better)",
        f"  - Resolution: {decomp['resolution']:.4f} (higher better)",
        f"  - Uncertainty: {decomp['uncertainty']:.4f}",
        f"- **Expected Calibration Error (ECE)**: {ece:.4f}",
        "",
    ]
    if platt_params is not None:
        lines += [
            "## Platt scaling cross-check",
            "",
            f"Logistic regression P(y=1|s) = σ(coef·s + intercept):",
            f"- intercept: `{platt_params['intercept']:.4f}`",
            f"- coef: `{platt_params['coef']:.4f}`",
            "",
            "The function maps any score to a calibrated probability if the",
            "per-bin estimates above are coherent.",
            "",
        ]
    lines += [
        "## Counts by pool",
        "",
    ]
    for pool, count in sorted(pool_counts.items()):
        lines.append(f"- Pool {pool}: {count} files")
    lines += ["", "## Counts by suspicion kind", ""]
    for kind, count in sorted(kind_counts.items(), key=lambda kv: -kv[1]):
        lines.append(f"- `{kind}`: {count} cells")
    lines += [
        "",
        "## Methods",
        "",
        "- **Wilson 95% CI** (Wilson 1927): robust at extreme proportions.",
        "- **Brier score & Murphy decomposition** (Brier 1950; Murphy 1973):",
        "  primary metric : `BS = Reliability − Resolution + Uncertainty`.",
        "- **Expected Calibration Error** (Naeini et al. 2015; Guo et al.",
        "  2017): weighted mean absolute deviation between stated confidence",
        "  and empirical precision per bin.",
        "- **Laplace posterior** (Beta(1,1) prior, Beta-Binomial conjugate):",
        "  per-bin smoothed precision : `(k+1)/(n+2)` : appropriate when",
        "  some bins have small n.",
        "- **Platt scaling** (Platt 1999): sigmoid fit cross-check; appropriate",
        "  here because n is below the ~1000-sample isotonic-regression",
        "  threshold (Niculescu-Mizil & Caruana 2005).",
        "",
        "Pool labels (per `scripts/extensive_fresh_audit.py`):",
        "",
        "- **A**: files marked 'with corruption' by the cache-walk classifier",
        "- **B**: files marked 'clean with gene symbols' (negative control)",
        "- **C**: files marked 'no gene symbols' (negative control)",
        "",
    ]
    OUT_MD.write_text("\n".join(lines))
    print(f"Wrote {OUT_MD}")

    # Reliability diagram
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed; skipping figure.")
        return 0
    confs = sorted(per_bin)
    precs = [per_bin[c]["precision"] for c in confs]
    los = [per_bin[c]["wilson_ci_low"] for c in confs]
    his = [per_bin[c]["wilson_ci_high"] for c in confs]
    sizes = [per_bin[c]["n_total"] for c in confs]

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5, label="perfect calibration")
    yerr = [
        [p - lo for p, lo in zip(precs, los)],
        [hi - p for p, hi in zip(precs, his)],
    ]
    ax.errorbar(
        confs, precs, yerr=yerr, fmt="o",
        capsize=4, label="UnCorrupt",
    )
    # Annotate each point with its sample size
    for c, p, n in zip(confs, precs, sizes):
        ax.annotate(
            f"n={n}", (c, p), textcoords="offset points",
            xytext=(8, -4), fontsize=8,
        )
    ax.set_xlabel("Stated confidence")
    ax.set_ylabel("Empirical precision")
    ax.set_title("Reliability diagram (Wilson 95% CI bars)")
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=120)
    print(f"Wrote {OUT_PNG}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
