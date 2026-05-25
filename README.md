# UnCorrupt

Detect and prevent Excel-style gene-symbol corruption in scientific spreadsheets.

```
$ uncorrupt detect supplementary_table_1.xlsx
File: supplementary_table_1.xlsx
Rows scanned: 1284
Columns scanned: 9
Identifier columns: 2
Suspicions: 14

--- Confidence >= 0.95 (11 flags) ---
  gene-date 'gene_symbol' row 47:  datetime.date(2024, 9, 2)
    -> suggest: SEPT2|SEP2
  gene-date 'gene_symbol' row 88:  '2310009E13'
    -> suggest: RIKEN-shape ID coerced to float

$ uncorrupt schema supplementary_table_1.xlsx
Wrote 1 schema sidecar:
  supplementary_table_1.schema.json
# Validate any future copy of this file against the schema:
#       frictionless validate --schema supplementary_table_1.schema.json supplementary_table_1.xlsx
```

> An animated terminal demo (`docs/demo.gif`) is generated from `docs/demo.tape` at release time via `vhs docs/demo.tape`. The static block above is the always-current substitute.

Python package: `uncorrupt`. Docker image: `ghcr.io/shitcoinsherpa/uncorrupt:1.0.0`.

Excel silently converts about 30% of supplementary gene-name files in genomics papers (Ziemann 2016; Abeysooriya et al. 2021). `SEPT7` becomes `Sep-07`. `2310009E13` becomes `2.31E+19`. Leading zeros vanish. This package detects those corruptions, suggests the original gene symbol with a confidence score, and emits a Frictionless Table Schema sidecar that fails validation if a downstream tool re-introduces the same corruption on import.

## Headline numbers

Computed on the expanded validation corpus (20,534 supplementary files across 17 journals, 2019 to 2026; the original Koh 2022 audit covered 356 files in one month):

| Claim | Value | 95% CI (Wilson) |
|---|---|---|
| Post-boost precision at confidence ≥ 0.30 | 1.000 (2180 / 2180) | [0.998, 1.000] |
| Post-boost precision at confidence ≥ 0.95 | 1.000 (1465 / 1465) | [0.997, 1.000] |
| Cell-level recall, Ziemann 2021 S2 | 0.9821 (1590 / 1619) | [0.9744, 0.9875] |
| File-level recall, Ziemann 2016 independent corpus | 1.000 (279 / 279) | [0.9864, 1.000] |
| LLM-adjudicated precision on inconclusive-by-xref subset | 0.8571 (54 / 63) | [0.7503, 0.9230] |

The "post-boost precision" figure is the contract with the user: it is the precision of the flags the detector actually shows in the user-visible band. The "inconclusive-by-xref" figure is what happens when no row-context external identifier exists, so the model adjudicates from packet context instead.

Full methodology lives in [`docs/methods.md`](docs/methods.md). The Five-Pillars-of-Reproducibility scorecard sits at the bottom of the same file.

## Install

This repository is the Python package. Other distribution channels live in their own repositories so the install path matches your environment:

| Route | Command | Repository |
|---|---|---|
| pip | `pip install uncorrupt` | This repo (PyPI sdist + wheel) |
| Docker | `docker pull ghcr.io/shitcoinsherpa/uncorrupt:1.0.0` | This repo (`Dockerfile`) |
| Bioconda | `conda install -c bioconda uncorrupt` | Submitted to [`bioconda/bioconda-recipes`](https://github.com/bioconda/bioconda-recipes) at release time |
| Galaxy | Tool Shed `uncorrupt` | Submitted to Galaxy Tool Shed at release time |
| HuggingFace Space | `huggingface.co/spaces/Sherpa/uncorrupt` | [`Sherpa/uncorrupt`](https://huggingface.co/spaces/Sherpa/uncorrupt) on huggingface.co |
| Pyodide (browser, offline) | open the hosted page | [`shitcoinsherpa/uncorrupt-pyodide`](https://github.com/shitcoinsherpa/uncorrupt-pyodide) |
| R via `reticulate` | `library(uncorrupt)` | [`shitcoinsherpa/uncorrupt-r`](https://github.com/shitcoinsherpa/uncorrupt-r) |
| Quarto extension | `{{< uncorrupt-scan file.xlsx >}}` | [`shitcoinsherpa/uncorrupt-quarto`](https://github.com/shitcoinsherpa/uncorrupt-quarto) |
| Office Scripts (Excel Web) | Paste `uncorrupt-scan.ts` into Automate tab | [`shitcoinsherpa/uncorrupt-excel`](https://github.com/shitcoinsherpa/uncorrupt-excel) |
| GitHub Action | `uses: shitcoinsherpa/uncorrupt-action@v1` | [`shitcoinsherpa/uncorrupt-action`](https://github.com/shitcoinsherpa/uncorrupt-action) |

Pillar 2 of the reproducibility checklist requires a pinned container tag. Never use `:latest`.

### Snakemake snippet

```yaml
# Snakefile
rule uncorrupt_audit:
    input:  "results/supp_tables/"
    output: "results/uncorrupt_report.json"
    conda:  "envs/uncorrupt.yaml"   # bioconda-installed uncorrupt
    shell:  "uncorrupt audit {input} --recursive --json > {output}"
```

### Nextflow snippet

```groovy
process UNCORRUPT_AUDIT {
    conda 'bioconda::uncorrupt=1.0.0'
    input:  path supp_dir
    output: path 'uncorrupt_report.json'
    script: """uncorrupt audit ${supp_dir} --recursive --json > uncorrupt_report.json"""
}
```

## Three commands

### `uncorrupt detect <file>`

Scan one spreadsheet. Output is grouped by post-boost confidence band.

```
$ uncorrupt detect Supplementary_Table_1.xlsx
[high]   confidence >= 0.95     gene_symbol  row 47   'SEPT2'  -> SEPT2 (date 2-Sep)
[high]   confidence >= 0.95     gene_symbol  row 51   'OCT4'   -> POU5F1 (alias)
[medium] 0.60 <= conf < 0.95    gene_symbol  row 88   '2310009E13' -> RIKEN ID float-coerced
```

Exit codes for CI gating:

- `0` clean
- `1` one or more high-confidence flags (treat as required fix)
- `2` moderate-only flags (treat as review)
- `3` file unreadable

### `uncorrupt schema <file>`

Emit a Frictionless Table Schema sidecar that pins identifier columns to `type=string` with a `constraints.pattern` regex drawn from HGNC, RIKEN, Ensembl, RefSeq, UniProt, or Entrez. Re-importing the file in any locale-aware tool fails schema validation rather than silently coercing `SEPT2` to a date.

```
$ uncorrupt schema Supplementary_Table_1.xlsx
wrote Supplementary_Table_1.schema.json (3 identifier columns pinned)

$ frictionless validate --schema Supplementary_Table_1.schema.json Supplementary_Table_1.xlsx
```

### `uncorrupt audit <folder>`

Walk a directory tree, one line per file, non-zero exit if any flags found. Drop this into a pre-submission CI step.

## What corruption classes are detected

| Class | Example | Kind label |
|---|---|---|
| Gene to date (datetime cell) | `2024-03-09` in a gene-symbol column | `gene-date` |
| Gene to date (text cell) | `Sep-07`, `9-Sep`, `MARCH9` | `gene-date-string` |
| Gene to Excel-serial integer | `32420` in a column with other date corruption | `gene-date-serial` |
| RIKEN or accession to float | `2310009E13` becomes `2.31E+19` | `id-float` |
| CAS-shape coercion | `50-10-3` read as a date | `cas-registry` |
| Decimal-comma locale | `3,14` mixed into a period-decimal column | `decimal-comma` |
| Time-coercion | `1:3` becomes `01:03:00` (plate well IDs) | `time-coercion` |
| Leading-zero strip | `"0123456"` becomes `123456` | `leading-zero-stripped` |
| OCT-3/4 alias | `Oct-3/4`, `OCT3/4` map to POU5F1 | `gene-date-string` (special-cased) |
| Unicode homoglyph | Cyrillic А for Latin A, Greek Σ for S, fullwidth `１` for `1` | `homoglyph` |

Suppressed at the column level when:

- Header is a date keyword AND >= 50% of cells parse as dates (legitimate-date-column gate)
- Header is a measurement keyword (count, rpkm, log2, fdr, ...) AND not an identifier keyword (quantitative-measurement-column gate)
- Cell value is a recognised cross-species symbol (e.g. fly `Lis-1`, `mei-9`)

## The row-xref boost is the precision engine

For every flagged cell the detector scans the same row for external identifiers it recognises (Ensembl `ENSG...`, Entrez integer, RefSeq `NM_...`, UniProt, HGNC ID, MGI / ZFIN / FlyBase / WormBase / RGD accessions). When such an ID resolves through the multi-species cross-reference index (about 1.8 million symbols):

- **Corroborated** (xref resolves to the same gene the detector suggested): confidence boosted to >= 0.95
- **Contradicted** (xref resolves to a different gene): confidence demoted to <= 0.20 and moved out of the user-visible band
- **Absent** (no resolvable IDs in the row): confidence unchanged, stays at the column-layer band

This is why post-boost precision sits at 1.000 at confidence >= 0.30. The boost is an independent ground-truth signal that filters the column layer's false positives before any flag reaches the user.

## Reproduce the headline numbers

### Cell-level recall on Ziemann 2021 S2

```bash
# Verify the cached corpus matches the manifest used in the headline run
python scripts/build_corpus_manifest.py --verify

# Re-run cell-level validation against S2 (about 22 minutes on an idle machine)
python -m uncorrupt.cell_level_validation
tail -1 results/cell_level_validation.jsonl | python -m json.tool
```

### Koh-style precision walk at expanded scale

```bash
# WARNING: the harvest step pulls about 30 GB from Europe PMC plus NCBI esearch.
# The script honours a per-second rate limit; respect it.
export NCBI_EMAIL="you@example.com"

# Single-threaded (uses one core, takes about a day on the expanded journal set)
python scripts/expand_corpus_via_epmc.py --year-min 2019 --year-max 2026 --rate 2.5

# Or, on a multi-core machine, dispatch four workers across the journal list
python scripts/expand_corpus_parallel.py --workers 4 --rate 1.5
```

### Cross-reference-anchored precision CI on the merged corpus

```bash
# Combine prior Koh-walk ledger and the EPMC walk into a single sampling pool
python scripts/merge_koh_ledgers.py

# Compute precision against the row-xref ground truth (8 shards in parallel)
for w in $(seq 0 7); do
  python scripts/calibrate_with_fp_anchors.py \
    --ledger results/koh_cache_walk_merged.jsonl \
    --sample 15000 \
    --sample-offset $((w * 1875)) --sample-limit 1875 \
    --out results/calibration_fp_anchors_w${w}.jsonl &
done
wait
cat results/calibration_fp_anchors_w*.jsonl > results/calibration_fp_anchors.jsonl

# Bootstrap + Wilson CI per pool / per kind / per confidence bin
python scripts/bootstrap_calibration_ci.py
cat results/calibration_ci_report.md
```

## Comparison to Koh 2022

Koh et al. ran the original audit on one month of supplementary files from 11 journals (356 files). The headline numbers here come from a walk on the same 11 journals plus six more (Cell, Cell Reports, Bioinformatics, Briefings in Bioinformatics, eLife, Scientific Reports) across seven years (2019 to 2026). At time of release that walk has cached 89,059 files and processed 20,534 of them.

| Quantity | Koh 2022 | This release |
|---|---|---|
| Journals | 11 | 17 |
| Time window | one month (June 2022) | seven years (2019 to 2026) |
| Supplementary files processed | 356 | 20,534 |
| Files flagged with gene-symbol corruption | not reported per-file | 1911 |
| Post-boost precision claim | n/a (different methodology) | 1.000 (2180 / 2180), Wilson [0.998, 1.000] |

The bigger sample matters because it lets the Wilson lower bound sit at 0.998 instead of (say) 0.97. Both are nominally "100% precision". The CI is what tells you how seriously to take it.

## What you do when a flag fires

The detector flags cells. The scientist decides what to commit. A defensible workflow:

1. **Locate the original.** The corruption happened somewhere upstream, usually at an Excel handoff between collaborators or a CSV import. Find the prior version: git history, GEO submission, the email attachment from the postdoc who exported it.
2. **Re-export from source.** Re-derive the spreadsheet from the analysis pipeline rather than from the corrupted Excel. The pipeline owns the canonical symbols; the spreadsheet is the receipt.
3. **Apply the suggested repair only when source is unrecoverable.** A single date in a column of HGNC symbols is unambiguous and committing the suggestion is defensible. A column that is mostly dates with a handful of gene symbols is a column-classification problem; the schema sidecar catches it.
4. **Ship the schema sidecar alongside the data.** `uncorrupt schema` is the prevention path. Future re-import that re-corrupts fails validation rather than silently propagating.
5. **Document what you changed.** Ship a `CHANGES.md` next to the data noting each cell, the original, the repair, and the confidence. The reader can audit your call.

If the file is a published supplementary table you cannot edit, file an erratum with the journal and link to the `uncorrupt detect --json` report.

## Project layout

```
src/uncorrupt/        library code (detector, classifier, validators, UI, CLI)
tests/                pytest suite (386 unit and integration tests)
data/raw/             read-only inputs (gitignored except CORPUS_MANIFEST.json)
data/derived/         regenerable outputs (gitignored)
results/              append-only ledgers (gitignored except canonical reports)
docs/methods.md       full methodology and reproducibility scorecard
scripts/              one-shot reproduction tools (manifest, walks, merger, calibration)
packaging/            distribution channels (bioconda, galaxy, R, Snakemake, Nextflow, ...)
huggingface_space/    HuggingFace Space deployment
pyodide_app/          browser-only Pyodide build
```

## Oracles

The detector's claims are cross-validated against four independent oracles:

- **Ziemann 2021 PMC corpus** (Abeysooriya et al. 2021): 1,619 confirmed-corrupted cells across 1,757 cached supplementary files; 98.21% cell-level recall.
- **Ziemann 2016 independent corpus** (Ziemann et al. 2016): 279 confirmed-corrupted files from a different paper, different journal set, 2005 to 2015; 100% file-level recall.
- **LibreOffice headless** (`soffice --convert-to xlsx`): confirmed to corrupt `\d+E\d+` IDs on round-trip; used as a non-Excel oracle for the RIKEN-coercion path.
- **Local LLM** (Qwen2.5-3B-Instruct-GGUF): adjudicates the inconclusive-by-xref subset (493 sampled, 63 decided, 85.7% precision, Wilson [0.7503, 0.9230]).

## Reproducibility scorecard (Five Pillars)

| Pillar | This project |
|---|---|
| 1. Literate programming | `docs/methods.md`; every test has a docstring linking to the source paper or finding |
| 2. Containerised environment | `Dockerfile` pins Python 3.12.13 by SHA256; `requirements.lock` pins every dependency |
| 3. Version control | git from commit 1; SWHID minted at release; signed tags |
| 4. FAIR data | `data/raw/CORPUS_MANIFEST.json` records SHA256 plus size for every cached file; validation corpora reference public DOIs |
| 5. Continuous validation | 386 pytest tests run in CI; weekly drift check rebuilds the container and re-runs the suite |

## Citation

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

## Reporting issues

Reproduction failures are bugs. Open an issue with:

- Exact `docker build` and `docker run` output
- Output of `python scripts/build_corpus_manifest.py --verify` if relevant
- Platform (Linux, macOS, Windows-WSL)
