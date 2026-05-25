"""Active-learning adjudication queue for the inconclusive calibration
flags.

Per Settles 2009 (Active Learning Literature Survey): when the labeled
corpus is small and most flags are inconclusive (no row xref evidence),
the cheapest way to tighten the precision interval is NOT to label
everything : it's to label the *most informative* unlabeled flags first.

Method: train a lightweight calibrated classifier on the labeled subset
using the detector's own features (confidence, kind, suggestion length,
column header tokens), score every inconclusive flag by margin-sampling
uncertainty (Settles eq. 4.4 : binary case: distance to 0.5), rank
descending. The top-k get manual review.

For UnCorrupt at n=27, adding 100 actively-selected labels can drop the
Wilson half-width from ±15pp to ±8pp (per the research review's
arithmetic).

This is the framework; the manual-review step is human-in-the-loop and
not automated.

Run:
    python scripts/active_learning_queue.py --top 100
    python scripts/active_learning_queue.py --top 50 --kind gene-date
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "results/calibration_fp_anchors.jsonl"
OUT_QUEUE = ROOT / "results/active_learning_queue.jsonl"

EXCLUDED = {"decimal-comma", "cas-registry", "unrecognized-symbol"}


def _featurize(record: dict) -> list[float]:
    """Extract a simple feature vector for the classifier.

    Features chosen to be cheap and informative without leaking the
    label:
      - confidence
      - number of candidate suggestions (more candidates → more ambiguous)
      - kind is gene-date (1) vs gene-date-serial (2) vs gene-date-string (3)
      - column header looks identifier-like (heuristic: contains
        'gene'/'symbol'/'id') vs date-like ('date'/'year'/'time') vs neither
    """
    conf = float(record["confidence"])
    sugg = record.get("suggestion") or ""
    n_candidates = sum(1 for p in sugg.split("|") if p.strip())
    kind_map = {"gene-date": 1, "gene-date-serial": 2, "gene-date-string": 3}
    kind_n = kind_map.get(record["kind"], 0)
    col = (record.get("column") or "").lower()
    id_like = int(any(t in col for t in ("gene", "symbol", "probe", "id")))
    date_like = int(any(t in col for t in ("date", "year", "time", "day")))
    return [conf, n_candidates, kind_n, id_like, date_like]


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--top", type=int, default=100,
                   help="How many inconclusive flags to queue for review.")
    p.add_argument("--kind", default=None,
                   help="Filter to one suspicion kind (e.g., gene-date).")
    args = p.parse_args()

    if not LEDGER.exists():
        print(f"FATAL: {LEDGER} missing", file=sys.stderr)
        return 1

    sys.path.insert(0, str(ROOT / "src"))
    from uncorrupt.corpus import load_hgnc
    hgnc = load_hgnc()

    def expand(suggs):
        e = set(suggs)
        for s in list(suggs):
            c = hgnc.resolve(s)
            if c:
                e.add(c)
            for o, n in hgnc.prev_symbol_to_current.items():
                if n == s or n == c:
                    e.add(o)
            for a, n in hgnc.alias_to_current.items():
                if n == s or n == c:
                    e.add(a)
        return {x.lower() for x in e}

    # Partition into labeled (corroborated/contradicted) and unlabeled
    # (inconclusive).
    labeled_X: list[list[float]] = []
    labeled_y: list[int] = []
    unlabeled: list[dict] = []
    for line in LEDGER.open():
        r = json.loads(line)
        if r["kind"] in EXCLUDED:
            continue
        if args.kind and r["kind"] != args.kind:
            continue
        if not r.get("suggestion"):
            continue
        suggs = {p.strip() for p in r["suggestion"].split("|") if p.strip()}
        ms = expand(suggs)
        resolved = [g for _, g in r.get("resolved_in_row") or []]
        if not resolved:
            unlabeled.append(r)
            continue
        labeled_X.append(_featurize(r))
        labeled_y.append(1 if any(g.lower() in ms for g in resolved) else 0)

    print(f"Labeled: {len(labeled_X)} cells "
          f"({sum(labeled_y)} positive, {len(labeled_y) - sum(labeled_y)} negative)")
    print(f"Inconclusive (candidates for adjudication): {len(unlabeled)}")

    if len(labeled_X) < 5 or len(set(labeled_y)) < 2:
        print("Too few labeled cells (need >=5 with both classes). "
              "Run the calibration walk first.", file=sys.stderr)
        return 2

    try:
        import numpy as np
        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import StandardScaler
    except ImportError:
        print("scikit-learn required: pip install scikit-learn", file=sys.stderr)
        return 3

    X = np.asarray(labeled_X, dtype=float)
    y = np.asarray(labeled_y, dtype=int)
    scaler = StandardScaler()
    Xs = scaler.fit_transform(X)
    clf = LogisticRegression(C=1.0, max_iter=1000)
    clf.fit(Xs, y)

    # Score unlabeled by uncertainty (margin sampling for binary = |p - 0.5|)
    U = np.asarray([_featurize(r) for r in unlabeled], dtype=float)
    Us = scaler.transform(U)
    probs = clf.predict_proba(Us)[:, 1]  # P(corroborated)
    uncertainty = 1.0 - np.abs(probs - 0.5) * 2  # 0 to 1, higher = more uncertain

    # Pair, sort, write top-k
    ranked = sorted(
        zip(uncertainty, probs, unlabeled),
        key=lambda t: -t[0],  # descending uncertainty
    )
    top_k = ranked[:args.top]

    OUT_QUEUE.parent.mkdir(parents=True, exist_ok=True)
    with OUT_QUEUE.open("w") as fh:
        for unc, prob, rec in top_k:
            entry = {
                "uncertainty": float(unc),
                "predicted_p_corroborated": float(prob),
                "pmc_id": rec.get("pmc_id"),
                "file": rec.get("file"),
                "sheet": rec.get("sheet"),
                "column": rec.get("column"),
                "row": rec.get("row"),
                "value": rec.get("value"),
                "kind": rec.get("kind"),
                "suggestion": rec.get("suggestion"),
                "confidence": rec.get("confidence"),
            }
            fh.write(json.dumps(entry) + "\n")
    print(f"Wrote top-{len(top_k)} active-learning queue: {OUT_QUEUE}")
    print()
    print("Highest-uncertainty samples (first 5):")
    for unc, prob, rec in top_k[:5]:
        print(f"  uncertainty={unc:.3f} p_corr={prob:.3f} "
              f"kind={rec['kind']} val={rec.get('value')!r}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
