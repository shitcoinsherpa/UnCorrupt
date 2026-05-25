# UnCorrupt confidence calibration with FP anchors

Source: stratified sample of 1784 Koh cache files (13 failed to load).

- Corroborated (xref → same gene): 242
- Contradicted (xref → different gene): 296
- Inconclusive (no resolvable row xref): 2191

## Per-bin precision on the labeled subset

| Stated | n_total | n_corroborated | n_contradicted | Precision | Wilson 95% CI |
|--------|---------|----------------|----------------|-----------|---------------|
| 0.30 | 1 | 1 | 0 | 1.000 | [0.207, 1.000] |
| 0.50 | 3 | 3 | 0 | 1.000 | [0.438, 1.000] |
| 0.60 | 17 | 17 | 0 | 1.000 | [0.816, 1.000] |

### Pool A (with-corruption)

| Stated | n_labeled | corroborated | contradicted | Precision |
|--------|-----------|--------------|--------------|-----------|
| 0.30 | 1 | 1 | 0 | 1.000 |
| 0.50 | 3 | 3 | 0 | 1.000 |
| 0.60 | 17 | 17 | 0 | 1.000 |

## Scalar metrics on the xref-labeled subset

- Brier score: 0.1886
- Expected Calibration Error (ECE): 0.4286
- Total labeled cells: 21 (21 positive, 0 negative)

## Methodology

Each flagged cell is labeled by independent xref evidence in the
same row. The detector's *unboosted* confidence (`row_context_boost=False`)
is the raw score; the xref decision is the ground-truth label.
Corroborated rows are TP; contradicted rows are FP candidates; rows
without any resolvable external ID are inconclusive and excluded.

HGNC equivalence: a suggested gene matches when any element of
`prev_symbol → current_symbol` and `alias → current_symbol`
transitively resolves to the corroborating ID's gene.

Wilson 95% CIs (Wilson 1927) over Beta posterior estimates
(Laplace prior). Brier score & ECE follow Brier 1950, Murphy 1973,
Naeini et al. 2015, Guo et al. 2017.
