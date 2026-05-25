# Platt scaling calibration

Labeled samples: 25 (14 positive, 11 negative).
Platt fit: P(y=1|s) = sigmoid(0.0363 * s + 0.2252).

## Per-bin recalibration

| Stated | n | correct | Empirical | Platt-calibrated |
|--------|---|---------|-----------|-------------------|
| 0.15 | 3 | 3 | 1.000 | 0.557 |
| 0.25 | 2 | 0 | 0.000 | 0.558 |
| 0.30 | 4 | 2 | 0.500 | 0.559 |
| 0.45 | 1 | 0 | 0.000 | 0.560 |
| 0.50 | 2 | 2 | 1.000 | 0.561 |
| 0.55 | 5 | 0 | 0.000 | 0.561 |
| 0.60 | 8 | 7 | 0.875 | 0.561 |

## Brier scores
- raw (stated confidence): 0.2859
- Platt-calibrated: 0.2464
- improvement: 0.0395

## Use

The detector loads `src/uncorrupt/calibration_params.json` at import
time if present. Each Suspicion's `calibrated_probability` field is
set to the Platt-mapped value. Legacy `confidence` stays for
backwards compatibility.