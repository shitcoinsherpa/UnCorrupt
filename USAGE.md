# UnCorrupt usage guide for researchers

A working guide for detecting and preventing Excel-style gene-symbol corruption in your supplementary spreadsheets.

## Install

**Option A: pip (recommended for an individual researcher).**

```bash
pip install uncorrupt
```

This puts the `uncorrupt` binary on your PATH.

**Option B: Docker (recommended for archival, journal pipelines, or anything published).**

```bash
docker pull ghcr.io/shitcoinsherpa/uncorrupt:1.0.0
```

Pillar 2 of the reproducibility contract requires a pinned tag. Never use `:latest`. The image ships with every dependency pinned by SHA256 hash.

**Option C: HuggingFace Space, Pyodide, Galaxy, R, Snakemake / Nextflow, Office Scripts, Quarto, GitHub Action.** See the README's "Install" table for routes targeting environments other than the command line.

## Three commands you need

### 1. `uncorrupt detect <file>`

```bash
uncorrupt detect Supplementary_Table_1.xlsx
```

Flags are grouped by post-boost confidence band:

- **>= 0.95**: high confidence; the row-xref boost corroborated the call against an independent identifier. Wilson 95% lower bound on this band is 0.997 on the validation corpus.
- **0.60 to 0.95**: moderate confidence; column-corroborated date corruption with weaker external evidence.
- **0.30 to 0.60**: low confidence; isolated cells the row-xref boost could not adjudicate.
- **< 0.30**: informational only; suppressed from the user-visible report by default.

Exit codes for use in CI:

- `0`: no corruption found
- `1`: one or more high-confidence flags (treat as required fix)
- `2`: moderate-only flags (treat as review)
- `3`: file unreadable or other error

JSON output for downstream tooling:

```bash
uncorrupt detect Supplementary_Table_1.xlsx --json > report.json
```

### 2. `uncorrupt schema <file>` (the prevention path)

```bash
uncorrupt schema Supplementary_Table_1.xlsx
```

Writes `Supplementary_Table_1.schema.json`, a [Frictionless Table Schema](https://specs.frictionlessdata.io/table-schema/) sidecar that:

- Pins identifier columns to `type="string"` with a `constraints.pattern` regex matching the inferred registry (HGNC, MGI, ZFIN, FlyBase, WormBase, RGD, RIKEN, Ensembl, RefSeq, UniProt, Entrez, HGNC ID).
- Leaves measurement columns at coarse-inferred types.

To enforce on re-import:

```bash
frictionless validate --schema Supplementary_Table_1.schema.json Supplementary_Table_1.xlsx
```

Or programmatically:

```python
from frictionless import Schema, Resource
Schema("Supplementary_Table_1.schema.json").validate(Resource("Supplementary_Table_1.xlsx"))
```

If a future copy of the file has been re-opened in Excel and a gene symbol coerced to a date, schema validation fails loud rather than silently producing a corrupted dataset.

### 3. `uncorrupt audit <folder>`

```bash
uncorrupt audit ./submission_materials/ --recursive
```

Walks every `.xlsx`, `.xls`, `.csv`, `.tsv` under the folder, prints one line per file with corruption flags, exits non-zero if any flags found. Drop into a pre-submission CI check.

## Corruption classes detected

| Class | Example | Kind label |
|---|---|---|
| Gene to date (datetime cell) | `2024-03-09` in a gene-symbol column | `gene-date` |
| Gene to date (text cell) | `Sep-07`, `9-Sep`, `MARCH9` | `gene-date-string` |
| Gene to Excel-serial integer | `32420` in a column with date corruption | `gene-date-serial` |
| RIKEN or accession to float | `2310009E13` becomes `2.31E+19` | `id-float` |
| CAS-shape coercion | `50-10-3` read as a date | `cas-registry` |
| Decimal-comma locale | `3,14` mixed into a period-decimal column | `decimal-comma` |
| Time-coercion | `1:3` becomes `01:03:00` | `time-coercion` |
| Leading-zero strip | `"0123456"` becomes `123456` | `leading-zero-stripped` |
| OCT-3/4 alias | `Oct-3/4`, `OCT3/4` map to POU5F1 | `gene-date-string` (special-cased) |
| Unicode homoglyph | Cyrillic А for Latin A, Greek Σ for S, fullwidth `１` for `1` | `homoglyph` |

Suppressed at the column level when:

- Column header is a date-like keyword AND >= 50% of cells parse as dates (legitimate-date-column gate)
- Column header is a measurement keyword (count, rpkm, log2, fdr, ...) AND not an identifier keyword (quantitative-measurement-column gate)
- Cell value is a recognised cross-species symbol

## How the row-xref boost works

For every flagged cell the detector scans the same row for external identifiers it recognises (Ensembl `ENSG...`, Entrez integer, RefSeq `NM_...`, UniProt, HGNC ID, MGI / ZFIN / FlyBase / WormBase / RGD accessions). When such an ID resolves through the multi-species cross-reference index (about 1.8 million symbols):

- **Corroborated** (xref resolves to the same gene the detector suggested): confidence boosted to >= 0.95
- **Contradicted** (xref resolves to a different gene): confidence demoted to <= 0.20 and moved out of the user-visible band
- **Absent** (no resolvable IDs in the row): confidence unchanged

Post-boost precision sits at 1.000 at confidence >= 0.30 because the boost is an independent ground-truth signal that filters the column layer's false positives before any flag reaches the user.

## What to do when a flag fires

The detector flags cells. The scientist decides what to commit. A defensible workflow:

1. **Locate the original.** The corruption happened somewhere upstream, usually at an Excel handoff between collaborators or a CSV import. Find the prior version: git history, GEO submission, the email attachment from the postdoc who exported it.
2. **Re-export from source.** Re-derive the spreadsheet from the analysis pipeline rather than from the corrupted Excel. The pipeline owns the canonical symbols. The spreadsheet is the receipt.
3. **Apply the suggested repair only when source is unrecoverable.** A single date in a column of HGNC symbols is unambiguous; committing the suggestion is defensible. A column that is mostly dates with a handful of gene symbols is a column-classification problem; the schema sidecar catches it.
4. **Ship the schema sidecar alongside the data.** `uncorrupt schema` is the prevention path. Future re-import that re-corrupts fails validation rather than silently propagating. Detection is a bandage. The schema is the cure.
5. **Document what you changed.** Ship a `CHANGES.md` next to the data noting each cell, the original, the repair, and the confidence. The reader can audit your call.

If the file is a published supplementary table you cannot edit, file an erratum with the journal and link to the `uncorrupt detect --json` report. Editors and curators are increasingly receptive. Ziemann's 2021 follow-up paper exists precisely because nothing changed when the 2016 paper was published.

## Reproducibility contract

UnCorrupt follows Ziemann's [Five Pillars of Computational Reproducibility](https://doi.org/10.1093/bib/bbad375):

| Pillar | This project |
|---|---|
| 1. Literate programming | `docs/methods.md`; every test has a docstring linking to the source paper or finding |
| 2. Containerised environment | `Dockerfile` pins Python 3.12.13 by SHA256; `requirements.lock` pins every dependency |
| 3. Version control | git from commit 1; SWHID minted at release; signed tags |
| 4. FAIR data | `data/raw/CORPUS_MANIFEST.json` records SHA256 plus size for every cached file; validation corpora reference public DOIs |
| 5. Continuous validation | 386 pytest tests run in CI; weekly drift check rebuilds the container and re-runs the suite |

## Citing

```bibtex
@software{uncorrupt,
  author  = {LLMSherpa},
  title   = {UnCorrupt: Detection and prevention of Excel-style gene-symbol corruption in scientific spreadsheets},
  year    = {2026},
  url     = {https://github.com/shitcoinsherpa/uncorrupt},
  version = {1.0.0}
}
```

Cite the underlying work as well:

- Ziemann M, Eren Y, El-Osta A (2016). *Gene name errors are widespread in the scientific literature.* Genome Biology 17(1):177. doi:10.1186/s13059-016-1044-7
- Abeysooriya M, Soria M, Kasu MS, Ziemann M (2021). *Gene name errors: Lessons not learned.* PLOS Computational Biology 17(7):e1008984. doi:10.1371/journal.pcbi.1008984
- Ziemann M, Poyet P, Lawson A (2023). *The five pillars of computational reproducibility: bioinformatics and beyond.* Briefings in Bioinformatics 24(6):bbad375. doi:10.1093/bib/bbad375
- Koh JLY, Brouwer L, et al. (2022). *Gene Updater: a Web tool that autocorrects and updates for Excel misidentified gene names.* Scientific Reports 12:18886. doi:10.1038/s41598-022-21243-y

## License

Apache-2.0. See `LICENSE`.
