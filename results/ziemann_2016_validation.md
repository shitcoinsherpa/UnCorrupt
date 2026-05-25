# Ziemann 2016 — independent recall validation (v0.7.5)

Run timestamp: 2026-05-21T06:01:37.082930+00:00 (against v0.7.5 detector)
Source: `data/raw/ziemann_2016_corpus/files/` (450 fetched files from Ziemann 2016 Additional File 1).

## Methodology

For every Ziemann-2016-Confirmed corruption file we have on disk:
1. Extract the *expected* gene symbol(s) from the published
   'Example Gene Name Conversion' column (e.g. "MARCH9 → 2005-09-01"
   yields expected = {MARCH9, MARCHF9}).
2. Run `detect_file(...row_context_boost=True)` on the cached file.
3. Collect every emitted suggestion (split on `|`).
4. TP iff `expected ∩ detector_suggestions ≠ ∅`.

Independent of the Koh 2022 corpus used by
`scripts/calibrate_with_fp_anchors.py` — different paper,
different sampling methodology, different journals, time period 2005-2015.

## Results

- Files processed: **450**
- True positives (TP): **279**
- False negatives — detector emitted nothing (FN-missed): 0
- False negatives — wrong suggestion (FN-sugg): 0
- Unreadable: 152
- Skipped (>50 MB): 0
- No-expected-genes (couldn't parse): 19
- **Validatable set**: 279
- **Recall**: **1.0000**
- **Wilson 95 % CI on recall**: [0.9864, 1.0000]

## Outcome breakdown

| Outcome | n |
|---|---:|
| TP | 275 |
| unreadable:UnrecoverableFile | 152 |
| no-expected-genes | 19 |
| TP-any-flag | 4 |