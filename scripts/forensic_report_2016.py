"""Generate forensic per-file reports for the Ziemann 2016 sample.

For each successfully-fetched file:
  1. Run UnCorrupt detector
  2. Compare detector flags against Ziemann's documented corruption note
  3. Emit a markdown one-pager: original cells, proposed fixes, ground-truth match

Outputs:
  results/ziemann_2016_forensic/<pubmed>.md     one report per file
  results/ziemann_2016_forensic/SUMMARY.md      aggregate (precision/recall)
"""
from __future__ import annotations

import json
import sys
import warnings
from datetime import UTC, datetime
from pathlib import Path

warnings.filterwarnings("ignore")

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from uncorrupt.app import UnrecoverableFile  # noqa: E402
from uncorrupt.detector import detect_file  # noqa: E402

LEDGER = PROJECT_ROOT / "data/raw/ziemann_2016_corpus/fetch_ledger.jsonl"
OUT_DIR = PROJECT_ROOT / "results/ziemann_2016_forensic"


def _parse_documented_corruption(note: str) -> set[str]:
    """Ziemann 2016 documented corruption can be e.g. 'CDK4-MARCH9',
    '2004-10-01 00:00:00', 'SEPT2'. Extract expected gene candidates."""
    import re
    out: set[str] = set()
    # Gene-name token
    for m in re.finditer(r"\b((?:SEPT|SEP|MARCH|MARC|DEC|DELEC|OCT|POU2F|"
                          r"NOV|APR|FEB|SEPTIN|MARCHF|MTARC)\d{1,2})\b", note,
                          re.IGNORECASE):
        out.add(m.group(1).upper())
    # Date string → reverse-decode
    m = re.match(r"(\d{4})-(\d{2})-(\d{2})", note)
    if m:
        from uncorrupt.detector import _reverse_gene_date
        from datetime import date
        try:
            d = date(int(m.group(1)), int(m.group(2)), int(m.group(3)))
            out.update(_reverse_gene_date(d))
        except (ValueError, TypeError):
            pass
    return out


def main() -> None:
    ledger = json.loads(LEDGER.read_text())
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    aggregate = {
        "n_files_tested": 0,
        "n_with_flags": 0,
        "n_with_ground_truth_match": 0,
        "n_load_errors": 0,
        "per_file": [],
    }

    for entry in ledger["successful"]:
        fp = PROJECT_ROOT / "data/raw/ziemann_2016_corpus/files" / entry["filename"]
        pubmed = entry["pubmed_id"]
        out_md = OUT_DIR / f"{pubmed}__{entry['filename']}.md"
        aggregate["n_files_tested"] += 1

        try:
            report = detect_file(str(fp))
        except UnrecoverableFile as exc:
            out_md.write_text(
                f"# {entry['filename']} (PubMed {pubmed})\n\n"
                f"**STATUS: Unrecoverable file** : {exc}\n"
            )
            aggregate["n_load_errors"] += 1
            continue
        except Exception as exc:
            out_md.write_text(
                f"# {entry['filename']} (PubMed {pubmed})\n\n"
                f"**STATUS: Detector error** : {type(exc).__name__}: {exc}\n"
            )
            aggregate["n_load_errors"] += 1
            continue

        # Ground-truth check
        expected = _parse_documented_corruption(entry["documented_corruption"])
        emitted: set[str] = set()
        for s in report.suspicions:
            if s.suggestion:
                for part in s.suggestion.split("|"):
                    emitted.add(part.strip().upper())
        ground_truth_match = bool(expected & emitted) if expected else None
        if report.suspicions:
            aggregate["n_with_flags"] += 1
        if ground_truth_match:
            aggregate["n_with_ground_truth_match"] += 1

        # Render markdown
        lines = [
            f"# UnCorrupt forensic report : PubMed {pubmed}",
            "",
            f"**Source file:** `{entry['filename']}`  ({entry['size_bytes']:,} bytes)",
            f"**SHA256:** `{entry['sha256']}`",
            f"**Publisher URL:** {entry['url']}",
            f"**Journal:** {entry['journal']} ({entry['year']})",
            f"**Ziemann 2016 documented corruption:** `{entry['documented_corruption']}`",
            "",
            f"## Detector verdict",
            "",
            f"- **Total flagged cells:** {len(report.suspicions)}",
            f"- **Sheets scanned:** {report.rows_scanned} rows × {report.columns_scanned} cols",
            f"- **Identifier columns:** `{report.identifier_columns}`",
            "",
        ]
        if expected:
            match_emoji = "✅" if ground_truth_match else "⚠️"
            lines.append(
                f"## Ground-truth match {match_emoji}\n"
                f"\n"
                f"- **Ziemann expected:** `{sorted(expected)}`\n"
                f"- **UnCorrupt emitted:** `{sorted(emitted) if emitted else 'none'}`\n"
                f"- **Match:** {'YES : at least one expected gene appears in UnCorrupt suggestions' if ground_truth_match else 'No exact match (may still be real corruption; see flags below)'}\n"
            )
        else:
            lines.append(
                f"## Ground-truth note unparseable\n\n"
                f"Ziemann's documented note (`{entry['documented_corruption']}`) "
                f"does not parse to an explicit gene-symbol expectation in this script. "
                f"Manual inspection required.\n"
            )
        lines.append("")
        if report.suspicions:
            lines.append("## Flagged cells\n")
            lines.append(
                "| Sheet | Column | Row | Original | Proposed | Kind | Confidence |\n"
                "|---|---|---:|---|---|---|---:|"
            )
            for s in report.suspicions[:20]:
                col = str(s.column)[:30]
                lines.append(
                    f"| {s.sheet or '_'} | `{col}` | {s.row} | "
                    f"`{str(s.value)[:30]}` | `{s.suggestion or ''}` | "
                    f"{s.kind} | {s.confidence:.2f} |"
                )
            if len(report.suspicions) > 20:
                lines.append(f"\n*…+{len(report.suspicions)-20} more*")
        else:
            lines.append("## No flags emitted\n\n"
                          "Either the file has been silently cleaned since "
                          "Ziemann's 2016 review, or the corruption is in a "
                          "different sheet/column the detector did not classify.\n")

        out_md.write_text("\n".join(lines))
        aggregate["per_file"].append({
            "pubmed_id": pubmed,
            "filename": entry["filename"],
            "n_flags": len(report.suspicions),
            "expected": sorted(expected),
            "emitted": sorted(emitted),
            "ground_truth_match": ground_truth_match,
            "report_md": str(out_md.relative_to(PROJECT_ROOT)),
        })

    # Summary
    total = aggregate["n_files_tested"] - aggregate["n_load_errors"]
    n_with_truth = sum(1 for r in aggregate["per_file"]
                        if r["ground_truth_match"] is not None)
    n_matched = aggregate["n_with_ground_truth_match"]

    summary_lines = [
        "# UnCorrupt vs Ziemann 2016 ground truth : fresh-corpus validation\n",
        f"**Generated:** {datetime.now(UTC).isoformat()}",
        f"**Sample size:** {aggregate['n_files_tested']} files (random sample from Ziemann 2016 Additional file 1)",
        f"**Load/detect errors:** {aggregate['n_load_errors']}",
        f"**Files with detector flags:** {aggregate['n_with_flags']} / {total}",
        "",
        "## Ground-truth recall (independent corpus)\n",
        f"- **Files with parseable Ziemann ground truth:** {n_with_truth}",
        f"- **UnCorrupt suggestion matches Ziemann expected gene:** {n_matched} / {n_with_truth} "
        f"({n_matched/max(1,n_with_truth)*100:.1f}%)",
        "",
        "## Per-file results\n",
        "| PubMed | File | Ziemann expected | UnCorrupt emitted | Match | n_flags | Report |",
        "|---|---|---|---|---|---:|---|",
    ]
    for r in aggregate["per_file"]:
        m = "✅" if r["ground_truth_match"] else ("⚠️" if r["ground_truth_match"] is False else ": ")
        e_short = ", ".join(r["expected"]) if r["expected"] else "(unparseable)"
        em_short = ", ".join(r["emitted"][:5]) if r["emitted"] else ": "
        if len(r["emitted"]) > 5:
            em_short += f", …+{len(r['emitted'])-5}"
        summary_lines.append(
            f"| {r['pubmed_id']} | `{r['filename'][:40]}` | "
            f"`{e_short}` | `{em_short}` | {m} | {r['n_flags']} | "
            f"[md]({Path(r['report_md']).name}) |"
        )
    (OUT_DIR / "SUMMARY.md").write_text("\n".join(summary_lines))
    print(f"\nWrote {len(aggregate['per_file'])} forensic reports + SUMMARY.md to {OUT_DIR}")
    print(f"\nHeadline: {n_matched}/{n_with_truth} = "
          f"{n_matched/max(1,n_with_truth)*100:.1f}% ground-truth match on this fresh corpus.")


if __name__ == "__main__":
    main()
