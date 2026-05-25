# Calibration confidence intervals (Wilson + BCa bootstrap)

Source: `calibration_fp_anchors.jsonl` — 3258 xref-labeled cells (2180 positive, 1078 negative).

## Overall precision

- **Empirical precision**: 0.669
- **Wilson 95% CI**: [0.653, 0.685]
- **BCa bootstrap 95% CI** (B=10,000): [0.653, 0.685]

## Per pool

| Pool | n | k | Precision | Wilson 95% | BCa 95% |
|------|---|---|-----------|------------|---------|
| A | 1897 | 1730 | 0.912 | [0.898, 0.924] | [0.899, 0.924] |
| C | 1361 | 450 | 0.331 | [0.306, 0.356] | [0.306, 0.356] |

## Per kind

| Kind | n | k | Precision | Wilson 95% | BCa 95% |
|------|---|---|-----------|------------|---------|
| gene-date | 1728 | 1676 | 0.970 | [0.961, 0.977] | [0.961, 0.977] |
| gene-date-serial | 1119 | 326 | 0.291 | [0.265, 0.319] | [0.265, 0.318] |
| gene-date-string | 285 | 178 | 0.625 | [0.567, 0.679] | [0.568, 0.677] |
| leading-zero-stripped | 118 | 0 | 0.000 | [0.000, 0.032] | [0.000, 0.032] |
| time-coercion | 8 | 0 | 0.000 | [0.000, 0.324] | [0.000, 0.324] |

## Per confidence bin — PRE-BOOST (column-corroboration layer only)

This is the raw output of the column-internal corroboration gate, with
`row_context_boost=False`. Useful for calibrating the bottom-layer signal.

| Stated | n | k | Precision | Wilson 95% | BCa 95% |
|--------|---|---|-----------|------------|---------|
| 0.15 | 1109 | 326 | 0.294 | [0.268, 0.321] | [0.268, 0.321] |
| 0.20 | 12 | 11 | 0.917 | [0.646, 0.985] | [0.583, 1.000] |
| 0.25 | 138 | 96 | 0.696 | [0.614, 0.766] | [0.616, 0.768] |
| 0.30 | 42 | 27 | 0.643 | [0.492, 0.770] | [0.476, 0.786] |
| 0.40 | 1 | 0 | 0.000 | [0.000, 0.793] | [0.000, 0.793] |
| 0.50 | 264 | 255 | 0.966 | [0.936, 0.982] | [0.936, 0.985] |
| 0.55 | 146 | 82 | 0.562 | [0.481, 0.640] | [0.479, 0.644] |
| 0.60 | 1405 | 1377 | 0.980 | [0.971, 0.986] | [0.972, 0.986] |
| 0.65 | 112 | 0 | 0.000 | [0.000, 0.033] | [0.000, 0.033] |
| 0.80 | 8 | 0 | 0.000 | [0.000, 0.324] | [0.000, 0.324] |
| 0.85 | 6 | 0 | 0.000 | [0.000, 0.390] | [0.000, 0.390] |
| 0.95 | 15 | 6 | 0.400 | [0.198, 0.643] | [0.200, 0.667] |

## Per confidence bin — POST-BOOST (user-facing default)

This is the actual `detect_file(... row_context_boost=True)` output —
what the user sees. The boost moves corroborated cells to ≥0.95 and
demotes contradicted cells to ≤0.20.

| Stated | n | k | Precision | Wilson 95% | BCa 95% |
|--------|---|---|-----------|------------|---------|
| 0.15 | 783 | 0 | 0.000 | [0.000, 0.005] | [0.000, 0.005] |
| 0.20 | 295 | 0 | 0.000 | [0.000, 0.013] | [0.000, 0.013] |
| 0.55 | 326 | 326 | 1.000 | [0.988, 1.000] | [0.988, 1.000] |
| 0.60 | 11 | 11 | 1.000 | [0.741, 1.000] | [0.741, 1.000] |
| 0.65 | 96 | 96 | 1.000 | [0.962, 1.000] | [0.962, 1.000] |
| 0.70 | 27 | 27 | 1.000 | [0.875, 1.000] | [0.875, 1.000] |
| 0.90 | 255 | 255 | 1.000 | [0.985, 1.000] | [0.985, 1.000] |
| 0.95 | 82 | 82 | 1.000 | [0.955, 1.000] | [0.955, 1.000] |
| 0.99 | 1383 | 1383 | 1.000 | [0.997, 1.000] | [0.997, 1.000] |

## Post-boost precision at user-facing thresholds

Cumulative precision: the precision of every flag whose post-boost
confidence is ≥ the threshold. This is what a user reading the
report at confidence T actually sees.

| Threshold | n flags | k correct | Precision | Wilson 95% | BCa 95% |
|-----------|---------|-----------|-----------|------------|---------|
| conf>=0.20 | 2475 | 2180 | 0.881 | [0.867, 0.893] | [0.867, 0.893] |
| conf>=0.30 | 2180 | 2180 | 1.000 | [0.998, 1.000] | [0.998, 1.000] |
| conf>=0.50 | 2180 | 2180 | 1.000 | [0.998, 1.000] | [0.998, 1.000] |
| conf>=0.60 | 1854 | 1854 | 1.000 | [0.998, 1.000] | [0.998, 1.000] |
| conf>=0.95 | 1465 | 1465 | 1.000 | [0.997, 1.000] | [0.997, 1.000] |

## Methodology

- **Wilson 95% CI** (Wilson 1927; Agresti & Coull 1998): standard
  precision/recall reporting CI. Robust at boundary proportions.
- **BCa bootstrap 95% CI** (scipy.stats.bootstrap default): bias-
  corrected accelerated bootstrap with 10,000 resamples. Used as
  robustness check at small n and near-boundary p.
- **Excluded kinds**: decimal-comma, cas-registry, unrecognized-symbol
  (informational flags, not corruption claims).
- **Excluded labels**: inconclusive (no xref evidence in the row).

Sample-size target per the validation-methodology research: n=457 for
Wilson 95% CI half-width ≤ ±2pp at p=0.95. Current n above reflects
the multi-species xref-expanded labeled set (HGNC + MGI + ZFIN +
FlyBase + WormBase, 1.6M total cross-references).

### Pre-boost vs post-boost — what they measure

- **Pre-boost** isolates the column-internal corroboration layer.
  Useful as a research diagnostic when redesigning the gate (e.g.
  added the column-corruption-count split inside this layer).
- **Post-boost** is the contract with the user. `detect_file(...)`
  defaults to `row_context_boost=True`; the boost step uses
  independent external IDs in the same row to corroborate or
  contradict, and the reported confidence reflects that adjudication.

The two numbers diverging by ~20pp is not a regression — it is the
boost layer doing its job. Cells that the column layer scored 0.60
but which the boost layer adjudicates as contradicted are demoted
to 0.20, i.e. moved OUT of the reported-flag band, not silently
kept at high confidence.