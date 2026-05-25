# Changelog: UnCorrupt

All notable changes are recorded here.

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] (2026-05-24)

First stable release.

### What it does

`uncorrupt` scans spreadsheets (`.xlsx`, `.xls`, `.csv`, `.tsv`) for the corruption patterns documented in Ziemann 2016, Abeysooriya 2021, and Koh 2022:

- Gene symbols silently turned into dates by Excel (`SEPT2` to `2-Sep`, `MARCH1` to `1-Mar`, `DEC1` to `1-Dec`, `OCT4` to `4-Oct`).
- RIKEN identifiers turned into floats (`2310009E13` to `2.31e+19`).
- Leading-zero strip, decimal-comma locale flip, time-coercion, alias variants, and Unicode homoglyph substitution.
- Reverse-decoded back into the original symbol where the mapping is unambiguous.

### Validated on real data

- **Recall**: 98.21% cell-level on Ziemann 2021 supplementary table S2 (1,590 / 1,619, Wilson 95% [0.9744, 0.9875]). 100% file-level on the independent Ziemann 2016 corpus (279 / 279, Wilson 95% [0.9864, 1.000]).
- **Precision**: 100% on the row-xref-anchored subset of the expanded 20,534-file walk (2,180 / 2,180, Wilson 95% [0.998, 1.000]).
- **LLM-adjudicated precision** on the inconclusive-by-xref subset: 85.7% Wilson [0.7503, 0.9230] (n=63 decided of 493 sampled).
- Expanded validation corpus is 58 times the size of Koh 2022 (20,534 files vs 356), spanning 17 journals across seven years (2019 to 2026) vs Koh's 11 journals across one month.

### Three interfaces

- **Console**: `pip install uncorrupt` ships a `uncorrupt` binary with `detect`, `schema`, and `audit` subcommands. Exit codes (0/1/2/3) are CI-gating-friendly.
- **Library**: `from uncorrupt.detector import detect_file` returns structured `Suspicion` records with confidence bands.
- **Gradio app**: `python -m uncorrupt.app` launches a local browser UI for one-file scans with banded review.

### Frictionless schema export (the prevention path)

`uncorrupt schema input.xlsx` emits a [Frictionless Table Schema](https://specs.frictionlessdata.io/table-schema/) sidecar that pins identifier columns to `type=string` with `constraints.pattern` regexes drawn from HGNC, RIKEN, Ensembl, RefSeq, UniProt, and Entrez. Re-importing the file in any locale-aware tool fails schema validation rather than silently coercing `SEPT2` to a date. Detection is the bandage. The schema is the cure.

### Distribution channels

`uncorrupt` is available through:

- **PyPI**: `pip install uncorrupt`
- **Bioconda**: `conda install -c bioconda uncorrupt` (noarch)
- **Container**: `ghcr.io/shitcoinsherpa/uncorrupt:1.0.0`
- **HuggingFace Space**: `huggingface.co/spaces/Sherpa/uncorrupt` (browser demo)
- **GitHub Action**: `shitcoinsherpa/uncorrupt-action@v1` (PR-blocking gate)
- **Browser-only**: Pyodide WASM build, no install required
- **Galaxy Tool Shed**: drag-and-drop wrapper
- **R**: `library(uncorrupt)` via `reticulate`
- **Snakemake / Nextflow**: pipeline-ready conda env + DSL2 module
- **Office Scripts**: Excel-native TypeScript scanner
- **Quarto**: `{{< uncorrupt-scan file.xlsx >}}` shortcode

### Detection methodology

1. **Pass 1**: symbol-shaped strings cross-referenced against HGNC current + previous symbols and the multi-species pool (MGI, ZFIN, FlyBase, WormBase, RGD; about 1.8M symbols total).
2. **Pass 2**: date-shaped cells reverse-decoded against the gene-date mapping; Excel-serial integers reverse-decoded against the same.
3. **Pass 3**: column-context corroboration. Homogeneous identifier columns lift confidence; quantitative-measurement, publication-metadata, and legitimate-date columns suppress.
4. **Row-xref boost**: when the same row carries an Ensembl, RefSeq, UniProt, Entrez, or HGNC ID that resolves to the suspected gene, post-boost confidence lifts toward 1.0. Contradicting IDs demote out of the user-visible band.
5. **Year-sanity gate**: decoded years outside [1900, 2100] are numeric noise, not corruption.

### Reproducibility (Five Pillars)

- **Literate**: `README.md` plus `docs/methods.md` document the chain of custody for every threshold.
- **Containerised**: `Dockerfile` and `environment.yml` pin specific versions. `:latest` is never used.
- **Versioned**: git from the first commit; release tags signed; code archived to Software Heritage.
- **FAIR data**: validation corpora (Ziemann 2016 sample, Ziemann 2021 S2, Koh 2022 cache, expanded EPMC walk) reference public DOIs. Internal walks ship JSONL ledgers with SHA256 checksums.
- **Continuously validated**: 386 unit and integration tests run on every commit; weekly fresh-container drift check; mutation-testing report committed at release.

## [Pre-release development summary]

The 0.1.0a0 to 0.7.7 trajectory (2026-05-14 to 2026-05-22) developed and validated the detector against three published reference corpora before tagging 1.0.0. The major milestones, condensed:

- **0.1 to 0.3**: scaffold; gene-date mapping (SEPT, MARCH, DEC, OCT-3-4, full HGNC MARCHF1-12 rename range); RIKEN-coercion signature; MIM and OMIM exclusions; locale variants.
- **0.4**: alias coverage (OCT-3/4, APR2 fungal genes); homoglyph normalisation (Greek, Cyrillic, fullwidth to ASCII); CAS Registry gating; string-form time corruption; decimal-comma locale.
- **0.5**: cosmic-ray mutation testing (kill rate raised above 86%); python-calamine loader path; HGNC drift warning; xref_status on every Suspicion; RIKEN case-insensitive resolve; null-result FP regression suite.
- **0.6**: column-classifier (gene_symbol vs entrez vs uniprot vs refseq vs ensembl vs measurement vs date vs free_text); pre-Pass-1 column-type gating; multispecies symbol pool (about 1.8M); column-internal corroboration band; xref-anchored calibration walk on the full Koh 2022 cache.
- **0.7**: legitimate-date-column suppression; quantitative-measurement-column suppression; year-sanity gate; full-corpus calibration walk (1,630 xref-labeled cells); Ziemann 2021 S2 recall validation; Ziemann 2016 independent recall validation; LLM-adjudication of inconclusive subset; distribution-channel buildout.

The 1.0.0 release extends the validation corpus to 20,534 supplementary files across 17 journals (2019 to 2026) and recomputes the precision Wilson CI on 2,180 xref-labeled cells (2,180 / 2,180 post-boost; Wilson 95% [0.998, 1.000]).

Per-version detail is preserved in git history.
