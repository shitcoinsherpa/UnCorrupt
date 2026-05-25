# How to use UnCorrupt

The fastest way is the browser: drop your spreadsheet at [huggingface.co/spaces/Sherpa/uncorrupt](https://huggingface.co/spaces/Sherpa/uncorrupt) and download a cleaned copy. This file is for everyone else: people running it on their own machine, putting it in a pipeline, or wrapping it from another tool.

## Install

The recommended path for a scientist:

```bash
pip install uncorrupt
```

That gives you both the desktop drag-and-drop UI and the command line. To open the UI:

```bash
uncorrupt-app
```

A browser tab opens at `http://localhost:7860`. Drop your file in, get a cleaned copy back.

If you do not want to touch your local Python, run the published Docker image:

```bash
docker run --rm -p 7860:7860 ghcr.io/shitcoinsherpa/uncorrupt:1.0.0
```

Same UI, available at `http://localhost:7860`.

For a lab managed with conda:

```bash
conda install -c bioconda uncorrupt    # (submitted at release time)
```

## The three commands

### 1. `uncorrupt detect <file>`

Scan a single spreadsheet.

```
$ uncorrupt detect Supplementary_Table_1.xlsx
File: Supplementary_Table_1.xlsx
Rows scanned: 1284
Columns scanned: 9
Identifier columns: 2
Suspicions: 14

--- High confidence (11 flags) -- what shows is virtually always real corruption
  Date mistaken for gene name  'gene_symbol' row 47: datetime.date(2024, 9, 2)
    suggest: SEPT2|SEP2
  Float (gene ID lost to scientific notation)  'gene_symbol' row 88: '2.31E+19'

--- Medium confidence (3 flags) -- review before accepting; mostly correct
  Excel-serial integer in gene column  'gene_symbol' row 142: 41897
    suggest: SEPT15
```

The confidence bands mean what they say:

- **High** (>= 0.95): an independent source in the same row, or several similar corruptions in the same column, agree the cell is mangled.
- **Medium** (0.60 to 0.95): the pattern matches a known corruption family but only one example in the column, so worth a glance.
- **Low** (0.30 to 0.60): isolated cells that look suspicious. Often a false alarm.
- **Info** (< 0.30): hidden by default; mostly noise.

Exit codes for CI:

- `0`: nothing flagged
- `1`: at least one high-confidence flag (block the merge / fail the build)
- `2`: medium-only flags (warn)
- `3`: file unreadable or other error

JSON output for scripts:

```bash
uncorrupt detect Supplementary_Table_1.xlsx --json > report.json
```

### 2. `uncorrupt schema <file>`

Write a [Frictionless](https://specs.frictionlessdata.io/table-schema/) JSON sidecar that pins the identifier columns to type `string` with a regex pattern. Commit this next to your spreadsheet. If anyone later re-imports the spreadsheet in a tool that would re-introduce the corruption, validation fails loud instead of silently rewriting your data.

```bash
uncorrupt schema Supplementary_Table_1.xlsx
# wrote Supplementary_Table_1.schema.json

frictionless validate --schema Supplementary_Table_1.schema.json Supplementary_Table_1.xlsx
# PASS
```

This is the prevention path. Use it on any file you publish.

### 3. `uncorrupt audit <folder>`

Walk a folder and scan every spreadsheet under it. Exits non-zero if anything is flagged. Drop into a pre-submission check.

```bash
uncorrupt audit ./submission_materials/ --recursive
```

## When a flag fires, what to do

The tool flags cells. The decision to commit a change is yours. A defensible workflow:

1. **Find the original.** The corruption happened somewhere upstream, usually an Excel handoff between collaborators or a CSV import. The git history, the GEO submission, the attachment from the postdoc who exported it: somewhere there is a version of the file that has not yet been mangled.
2. **Re-export from source.** Re-derive the spreadsheet from the analysis pipeline rather than from the corrupted Excel. The pipeline owns the canonical symbols. The spreadsheet is the receipt.
3. **Apply the suggested fix only when the source is unrecoverable.** A single date in a column of HGNC symbols is unambiguous and committing the suggestion is defensible. A column that is mostly dates with a handful of gene symbols is a column-classification problem and the schema sidecar catches it.
4. **Ship the schema sidecar.** Whatever you submit alongside your paper, attach the JSON sidecar from `uncorrupt schema`. Any future re-import that re-corrupts will fail validation rather than silently propagating.
5. **Document what you changed.** If you edited values, ship a one-page `CHANGES.md` next to the data noting each cell, the original, the repair, and the confidence. Anyone reading the paper later can audit your call.

If the file is a published supplementary table you cannot edit, file an erratum with the journal and attach the `uncorrupt detect --json` report. The 2021 follow-up paper by Ziemann's group exists precisely because nothing changed when the 2016 paper documenting this was published.

## Citation

```bibtex
@software{uncorrupt,
  author  = {LLMSherpa},
  title   = {UnCorrupt: Repair Excel-mangled gene symbols in genomics spreadsheets},
  year    = {2026},
  url     = {https://github.com/shitcoinsherpa/UnCorrupt},
  version = {1.0.0}
}
```

Please also cite the scientists whose work this implements:

- Ziemann M, Eren Y, El-Osta A (2016). *Gene name errors are widespread in the scientific literature.* Genome Biology 17(1):177. doi:10.1186/s13059-016-1044-7
- Abeysooriya M, Soria M, Kasu MS, Ziemann M (2021). *Gene name errors: Lessons not learned.* PLOS Computational Biology 17(7):e1008984. doi:10.1371/journal.pcbi.1008984

## License

Apache-2.0. See [LICENSE](LICENSE).
