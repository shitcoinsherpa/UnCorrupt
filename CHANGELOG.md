# Changelog

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). This project uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] (2026-05-25)

First stable release.

### What it does

Finds and fixes the Excel-mangled gene-name corruption patterns documented by Ziemann 2016 and Abeysooriya 2021. Scientists upload a spreadsheet (`.xlsx`, `.xls`, `.csv`, or `.tsv`), see what was caught with plain-English explanations and confidence labels, and download a cleaned copy. The same engine is available as a Python library, a command-line tool, and a Gradio drag-and-drop UI (browser-hosted or run locally).

### What it catches

`SEPT2` turned into `2-Sep` or `2024-09-02`. `OCT4` turned into `4-Oct`. `MARCH1` turned into `1-Mar`. `DEC1` turned into `1-Dec`. RIKEN identifiers like `2310009E13` coerced into floats like `2.31E+19`. Leading zeros stripped from string IDs. Cyrillic / Greek / fullwidth look-alike characters hidden inside ASCII names. Excel-serial integers in gene columns. CAS numbers misread as dates. Decimal-comma locale collisions.

### How accurate

Validated on three independent published lists of confirmed-corrupted spreadsheets:

| Test | Score | Source |
|---|---|---|
| 1,672 verifiable cells hand-marked in an EPMC-expanded Ziemann 2021 S2 corpus | **1,591 caught (95.16%)**. 81 misses are upstream annotation/file-mismatch issues, not detector misses | `results/cell_level_validation_2026-05-25_verified.jsonl` |
| 279 confirmed-corrupted files from Ziemann 2016 | **279 of 279 (100%)** | `results/ziemann_2016_validation.jsonl` |
| 19,500 supplementary files from a 2022 to 2026 Koh-replication scan across 11 high-impact genetics journals | **31,377 high-confidence corruption flags** in **988 distinct papers** | `results/koh_walk_2026-05-25_verified.checkpoint.json` |
| 2,180 of those flags cross-referenced against external gene-ID registries (HGNC / MGI / ZFIN / FlyBase / WormBase / NCBI rat) | **100% post-boost precision at every confidence bin ≥ 0.50** (326/326 at 0.55, 11/11 at 0.60, 96/96 at 0.65, all higher) | `results/calibration_ci_report.md` |

In plain language: if you have corrupted gene names in your spreadsheet, UnCorrupt almost certainly finds them. When it flags a cell, it is virtually always right.

### Three ways to run it

- **Browser, no install**: [`huggingface.co/spaces/Sherpa/uncorrupt`](https://huggingface.co/spaces/Sherpa/uncorrupt)
- **Local desktop UI**: `pip install uncorrupt && uncorrupt-app`
- **Command line**: `pip install uncorrupt && uncorrupt detect file.xlsx`
- **Pipeline integration**: GitHub Action, R wrapper, Quarto extension, Office Scripts, Bioconda, Galaxy (see README for the full list)

### What this release does NOT do

- **RIKEN-style identifiers that Excel has already turned into floats**: once `2310009E13` becomes `2.31E+19`, the original digits are mathematically gone. UnCorrupt flags the cell so you know to re-fetch from source; it cannot reconstruct the original.
- **Brand-new corruption patterns we have never seen**: knows what the published literature documents. If your file did something genuinely novel, [open an issue](https://github.com/shitcoinsherpa/UnCorrupt/issues) with the file attached and we add it.

### What changed since pre-release

Almost everything. The 0.X versions were internal development cycles, validated against three reference corpora before tagging 1.0.0. The full methodology and the source data used to compute the accuracy numbers above live in [`docs/methods.md`](docs/methods.md) for anyone wanting to verify or reproduce.
