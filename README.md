# UnCorrupt

[![CI](https://github.com/shitcoinsherpa/UnCorrupt/actions/workflows/ci.yml/badge.svg)](https://github.com/shitcoinsherpa/UnCorrupt/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/uncorrupt.svg)](https://pypi.org/project/uncorrupt/)
[![Python](https://img.shields.io/pypi/pyversions/uncorrupt.svg)](https://pypi.org/project/uncorrupt/)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![HuggingFace Space](https://img.shields.io/badge/%F0%9F%A4%97%20Try%20it-Space-yellow.svg)](https://huggingface.co/spaces/Sherpa/uncorrupt)

**Excel keeps turning your gene names into dates. This unfucks them.**

![UnCorrupt CLI demo](docs/demo.gif)

`SEPT2` becomes `2-Sep`. `MARCH1` becomes `1-Mar`. `OCT4` becomes `4-Oct`. RIKEN identifiers get turned into floats like `2.31E+19`. Leading zeros disappear. About one in three supplementary gene-name spreadsheets in published genomics papers has at least one of these (Ziemann 2016, 2021). It is your data being silently corrupted by your own spreadsheet program. UnCorrupt finds the damage, tells you what it was originally, and gives you a clean file back.

## Just use it (no install)

[**huggingface.co/spaces/Sherpa/uncorrupt**](https://huggingface.co/spaces/Sherpa/uncorrupt)

Drop your `.xlsx` or `.csv` in. Get a report. Download the fixed file. Your data is processed once and discarded; nothing is logged.

For data that can't leave your machine: the [browser-only Pyodide version](https://shitcoinsherpa.github.io/uncorrupt-pyodide/) runs entirely in your own browser tab via WebAssembly. Your file never gets uploaded anywhere; the detector runs locally on your CPU. ([source code](https://github.com/shitcoinsherpa/uncorrupt-pyodide))

## How well does it work

We tested it against three independently-published lists of confirmed-corrupted spreadsheets:

| The test | Score |
|---|---|
| 1,672 cells hand-marked as corrupted in an EPMC-expanded extension of the Ziemann 2021 supplementary table | UnCorrupt caught **1,591 of them** (about 95 out of every 100). Of the 81 it missed, prior analysis showed every one is an upstream annotation/file-mismatch issue, not a detector miss. (Re-run today 2026-05-25; see [`results/cell_level_validation_2026-05-25_verified.jsonl`](results/cell_level_validation_2026-05-25_verified.jsonl).) |
| 279 files from the older Ziemann 2016 study, every one of which contained corruption | UnCorrupt caught **every single one** (279 of 279) |
| **19,500 supplementary files** from a recent (2022 to 2026) replication of the Koh 2022 scan, across 11 high-impact genetics journals | **31,377 high-confidence corruption flags** in **988 distinct papers.** Of the xref-validated subset (877 files, 2,180 positive cells cross-referenced against HGNC / MGI / ZFIN / FlyBase / WormBase / NCBI rat), **post-boost precision was 100% at every confidence band at or above 0.50** (326/326 at 0.55, 11/11 at 0.60, 96/96 at 0.65, and all higher). Re-run today 2026-05-25; see [`results/koh_walk_2026-05-25_verified.checkpoint.json`](results/koh_walk_2026-05-25_verified.checkpoint.json) and [`results/calibration_ci_report.md`](results/calibration_ci_report.md). |

In plain English: if you have corrupted gene names in your spreadsheet, UnCorrupt almost certainly finds them. And when it flags a cell, it is virtually always right.

## What corruption looks like, with examples

| What you typed in | What Excel saved | What UnCorrupt tells you it was |
|---|---|---|
| `SEPT2` | `2-Sep` (a date) | `SEPT2` (the septin-2 gene) |
| `SEPT2` | `2024-09-02` (full datetime) | `SEPT2` |
| `MARCH1` | `1-Mar` | `MARCH1` (now renamed MARCHF1 by HGNC) |
| `OCT4` | `4-Oct` | `OCT4` (now formally POU5F1) |
| `DEC1` | `1-Dec` | `DEC1` |
| `2310009E13` | `2.31E+19` (huge float) | flags it, says the original digits can't be recovered from a float (re-upload the original if you have it) |
| `0123456` | `123456` (leading zero stripped) | flags it |
| `BRCA1` written with a Cyrillic `А` instead of Latin `A` | looks identical to humans, breaks every database lookup | normalises to plain ASCII |

## Other ways to run it

**As a desktop app (same drop-the-file UI, runs on your machine):**

```bash
pip install uncorrupt
uncorrupt-app
```

That opens the UI at http://localhost:7860 in your browser.

**Or via Docker** if you do not want to touch your local Python:

```bash
docker run --rm -p 7860:7860 ghcr.io/shitcoinsherpa/uncorrupt:1.0.2
```

**As a command line tool** for pipelines and CI:

```bash
uncorrupt detect supplementary_table_1.xlsx          # scan one file
uncorrupt audit ./submission_materials/ --recursive  # scan a folder
```

Exit code `0` is clean; `1` means corruption found.

**As a plugin in your existing tool**:

- R: [`uncorrupt-r`](https://github.com/shitcoinsherpa/uncorrupt-r)
- Quarto / R Markdown reports: [`uncorrupt-quarto`](https://github.com/shitcoinsherpa/uncorrupt-quarto)
- Excel on the web (paste this script into the Automate tab): [`uncorrupt-excel`](https://github.com/shitcoinsherpa/uncorrupt-excel)
- GitHub Actions (block a pull request if a corrupted supplementary file got committed): [`uncorrupt-action`](https://github.com/shitcoinsherpa/uncorrupt-action)
- Bioconda, Galaxy: shipping at the next release

## What it cannot do

We are honest about the limits:

- **RIKEN-style identifiers that Excel turned into floats.** Once `2310009E13` becomes the number `2.31E+19`, the original digits are mathematically gone. UnCorrupt flags the cell, but it cannot reconstruct the ID from the float. The fix is to find a non-corrupted source of the same data and reupload that.
- **Brand-new corruption patterns we have never seen.** UnCorrupt knows the corruption families catalogued by Ziemann 2016, Abeysooriya 2021, and Koh 2022. If your Excel did something genuinely novel, file an issue with the file attached and we will add it.
- **Cells that look like they could be either a gene or a real date.** A cell that says `2024-03-09` in a column that is clearly publication dates is a real date and we leave it alone. A cell that says `2024-03-09` in a column of gene symbols is `MARCH9` and we flag it. We do this by looking at the column as a whole, but if your file is weird enough we might guess wrong. The fix is to send us the file.

## Found something wrong?

Open an issue: https://github.com/shitcoinsherpa/UnCorrupt/issues. The fastest fix for any bug is a copy of the file that broke it.

## Citing this in a paper

```bibtex
@software{uncorrupt,
  author  = {LLMSherpa},
  title   = {UnCorrupt: Repair Excel-mangled gene symbols in genomics spreadsheets},
  year    = {2026},
  url     = {https://github.com/shitcoinsherpa/UnCorrupt},
  version = {1.0.2}
}
```

Please also cite the scientists whose work this exists to address:

- Ziemann M, Eren Y, El-Osta A (2016). *Gene name errors are widespread in the scientific literature.* Genome Biology 17(1):177.
- Abeysooriya M, Soria M, Kasu MS, Ziemann M (2021). *Gene name errors: Lessons not learned.* PLOS Computational Biology 17(7):e1008984.

## License

Apache-2.0. See [LICENSE](LICENSE).

## How accurate is it really

Skeptical readers and anyone implementing it into a journal submission pipeline: the full methodology, every validation walk, exact confidence intervals, and the source data used to compute the scores in the table above all live in [`docs/methods.md`](docs/methods.md). We tried to put nothing in this README that we cannot show our work on.

---

*Maintainer: [@LLMSherpa](https://x.com/LLMSherpa) / [bt6.gg](https://bt6.gg). Issues and bug reports: [github.com/shitcoinsherpa/UnCorrupt/issues](https://github.com/shitcoinsherpa/UnCorrupt/issues).*
