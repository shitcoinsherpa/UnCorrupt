"""Value-pattern-based column-type classification.

Inspired by Fortemi's document-type auto-detection (filename pattern →
extension → content magic). We apply the same hierarchical approach to
columns inside spreadsheets:

- Inspect the column's actual values
- Match against known identifier patterns (HGNC symbols, Entrez IDs,
  UniProt accessions, RefSeq, Ensembl, etc.)
- Classify the column by dominant pattern

Detector behaviour then becomes column-type-aware: only run
uncorrupt detection on `gene_symbol` columns. Skip external-DB-ID
columns (where integers/accessions are the canonical format, not
corruption). This addresses the EntrezGeneID-false-positive class
identified by row-context validation.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import date, datetime
from typing import Literal

import pandas as pd

from .registries import GENE_LIKE_PATTERN, RIKEN_PATTERN
from .xref_lookup import XrefIndex, _ID_PATTERNS

ColumnType = Literal[
    "gene_symbol",      # gene symbol candidates (date-corruption plausible)
    "entrez",           # Entrez Gene IDs (numeric, NOT uncorrupt)
    "uniprot",          # UniProt accessions
    "refseq",           # RefSeq accessions (NM_/NR_/NP_)
    "ensembl",          # Ensembl IDs (ENSG/ENST/ENSP)
    "hgnc_id",          # HGNC IDs
    "riken",            # RIKEN-style identifiers (date-corruption plausible via E notation)
    "measurement",      # numeric measurement data
    "date",             # date/time values
    "free_text",        # mixed text, no obvious pattern
    "empty",            # no usable data
]


@dataclass
class ColumnClassification:
    column_type: ColumnType
    confidence: float
    n_samples: int
    n_matches: int  # cells matching the winning pattern
    reason: str


_NUMERIC_STRING_RE = __import__("re").compile(
    r"^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$"
)


_MISSING_SENTINELS: frozenset[str] = frozenset({
    "", "-", "--", "n/a", "na", "null", "none", "nan", ".", "?", "n.a.",
    "n/d", "n.d.", "nd",
})


def _classify_value(v, xref: XrefIndex | None) -> str:
    """Return the kind of value (matches one of the patterns, or 'numeric',
    'date', 'string-other', 'empty', 'placeholder')."""
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return "empty"
    if isinstance(v, (date, datetime)):
        return "date"
    if isinstance(v, bool):
        return "bool"
    if isinstance(v, (int, float)):
        # Bare 0 is an "absent measurement" sentinel in most omics columns : 
        # it should not count as evidence of a measurement column when other
        # cells in the same column carry gene-symbol strings.
        if v == 0:
            return "placeholder"
        return "numeric"
    if not isinstance(v, str):
        return "other"
    s = v.strip()
    if not s:
        return "empty"
    if s.lower() in _MISSING_SENTINELS:
        return "placeholder"
    if s in ("0", "0.0", "-0"):
        return "placeholder"
    # Identifier patterns
    for kind, pat in _ID_PATTERNS.items():
        if pat.fullmatch(s):
            return kind
    # Gene-symbol shape (high-quality match using HGNC if available)
    if xref is not None:
        # Check membership in HGNC current symbols by quick string comparison
        # : done via xref's HGNC index (built indirectly)
        pass
    if GENE_LIKE_PATTERN.fullmatch(s):
        return "gene_like"
    if RIKEN_PATTERN.fullmatch(s):
        return "riken_pattern"
    # Numeric strings (signed/unsigned, integer/float, scientific-notation).
    # Necessary for files loaded via pd.read_csv with dtype=object : every
    # measurement comes through as a string and would otherwise misclassify
    # as `string-other` → `free_text`, leaking measurement columns into Pass 2.
    #
    # But: large-exponent scientific-notation strings like `'2310009E13'` are
    # the RIKEN-coercion sentinel : those are corrupted IDs, not measurements.
    # Reject |exponent| >= 10 to keep id-float detection working.
    if _NUMERIC_STRING_RE.fullmatch(s):
        import re as _re
        sci = _re.fullmatch(r"^[+-]?\d*\.?\d+[eE]([+-]?\d+)$", s)
        if sci and abs(int(sci.group(1))) >= 10:
            return "string-other"
        return "numeric"
    return "string-other"


from functools import cache as _cache


@_cache
def _hgnc_symbol_set() -> set[str]:
    """Cached union of HGNC current + prev symbols (built once per process)."""
    from .corpus import load_hgnc
    hgnc = load_hgnc()
    return hgnc.current_symbols | set(hgnc.prev_symbol_to_current.keys())


def classify_column(
    series: pd.Series,
    column_name: str,
    xref: XrefIndex | None = None,
    threshold: float = 0.30,
    scan_limit: int = 200,
    min_samples: int = 5,
) -> ColumnClassification:
    """Decide the column's type from value patterns.

    Returns the dominant column type if any pattern reaches `threshold`
    fraction of non-empty cells. Otherwise returns 'free_text'.

    `min_samples` is the floor for confident classification. Columns with
    fewer non-empty values fall back to 'free_text' regardless of pattern
    distribution : too small a sample for value-pattern verdicts to be
    reliable.
    """
    from collections import Counter

    hgnc_symbols = _hgnc_symbol_set()

    non_null = series.dropna().head(scan_limit)
    if len(non_null) == 0:
        return ColumnClassification(
            column_type="empty", confidence=0.0,
            n_samples=0, n_matches=0,
            reason="empty column",
        )
    if len(non_null) < min_samples:
        return ColumnClassification(
            column_type="free_text", confidence=0.0,
            n_samples=len(non_null), n_matches=0,
            reason=f"too few samples ({len(non_null)} < {min_samples}) "
                   f"for confident pattern classification",
        )

    # also tally cross-species symbol membership so that columns of
    # mouse `Arnt, Ebag9, Pigs`, fly `mei-9`, zebrafish `flh`, etc. classify
    # as gene_symbol : not free_text. Without this, HGNC-only matching
    # misses every non-human supplementary file (and Ziemann's 2016
    # corpus is dominated by mouse + fly + worm + fish).
    # _multispecies_symbols_lower() returns a frozenset cached at module
    # level : reference it directly, do NOT copy (1.8M elements).
    multispecies_lower: frozenset[str] | set[str] = frozenset()
    try:
        from .detector import _multispecies_symbols_lower
        multispecies_lower = _multispecies_symbols_lower()
    except Exception:
        pass

    kinds: Counter = Counter()
    n_hgnc = 0
    n_multispecies = 0
    for v in non_null:
        k = _classify_value(v, xref)
        kinds[k] += 1
        # Also tally explicit HGNC membership (case-insensitive : Excel
        # often round-trips uppercase symbols as mixed-case after autoFit)
        if isinstance(v, str):
            vs = v.strip()
            if vs in hgnc_symbols or vs.upper() in hgnc_symbols:
                n_hgnc += 1
            elif multispecies_lower and vs.lower() in multispecies_lower:
                n_multispecies += 1

    # Placeholders ('-', 'NA', '0', empty string, ...) are missing-value
    # sentinels in scientific data : exclude them from the denominator so a
    # column of {30 gene symbols, 170 placeholder zeros} classifies as
    # gene_symbol, not measurement.
    n_placeholder = kinds.pop("placeholder", 0)
    kinds.pop("empty", 0)  # also drop empties from the denominator
    total = sum(kinds.values())
    if total == 0:
        return ColumnClassification(
            column_type="empty", confidence=0.0,
            n_samples=n_placeholder, n_matches=0,
            reason=f"only {n_placeholder} placeholder/empty cells; no signal",
        )
    frac_hgnc = n_hgnc / total
    frac_combined = (n_hgnc + n_multispecies) / total

    # Prioritized classification
    # Highest specificity first: explicit HGNC matches override others
    if frac_hgnc >= threshold:
        return ColumnClassification(
            column_type="gene_symbol", confidence=frac_hgnc,
            n_samples=total, n_matches=n_hgnc,
            reason=f"{n_hgnc}/{total} cells are HGNC current/prev symbols",
        )
    # cross-species symbol match (mouse / fly / worm / fish).
    # Same threshold as HGNC; non-human supplementary files (Ziemann 2016
    # Drosophila / mouse files) classify correctly without losing HGNC
    # precision (the HGNC check ran first and would have already matched).
    if frac_combined >= threshold:
        return ColumnClassification(
            column_type="gene_symbol", confidence=frac_combined,
            n_samples=total, n_matches=n_hgnc + n_multispecies,
            reason=f"{n_hgnc} HGNC + {n_multispecies} multi-species symbols "
                   f"out of {total} cells",
        )

    # External-DB-ID columns: integer/accession patterns dominate
    for db in ("hgnc_id", "ensembl", "refseq", "uniprot", "entrez"):
        n = kinds[db]
        if n / total >= threshold:
            return ColumnClassification(
                column_type=db, confidence=n / total,
                n_samples=total, n_matches=n,
                reason=f"{n}/{total} cells match {db} pattern",
            )

    n_riken_pat = kinds["riken_pattern"]
    if n_riken_pat / total >= threshold:
        return ColumnClassification(
            column_type="riken", confidence=n_riken_pat / total,
            n_samples=total, n_matches=n_riken_pat,
            reason=f"{n_riken_pat}/{total} cells match RIKEN pattern",
        )

    n_gene_like = kinds["gene_like"]
    if n_gene_like / total >= threshold:
        return ColumnClassification(
            column_type="gene_symbol", confidence=n_gene_like / total,
            n_samples=total, n_matches=n_gene_like,
            reason=f"{n_gene_like}/{total} cells match gene-name shape "
                   f"(uppercase letters + digits)",
        )

    # structural date-column detection (LOWER threshold).
    # Before the numeric / measurement check, recognise columns where
    # a minority but non-trivial fraction of cells are real date/datetime
    # objects. Empirically caught in the calibration walk:
    # supplementary_table_8_ddac017.xlsx row 279 had a publication-date
    # column with only ~5% of rows populated as dates; the rest were
    # strings or NaN. At 5% the column wasn't classified as "date"
    # under the 0.5 threshold, so gene-date FPs fired. The fix: at
    # >= 20% date cells with NO conflicting gene-like cells, classify
    # as date : these are publication-date / sample-date columns where
    # gene-date detection is always wrong.
    n_date = kinds["date"]
    if (n_date / total >= 0.20
            and kinds["gene_like"] == 0
            and kinds["riken_pattern"] == 0):
        return ColumnClassification(
            column_type="date", confidence=n_date / total,
            n_samples=total, n_matches=n_date,
            reason=(
                f"{n_date}/{total} cells are dates and no gene-like "
                f"strings present : sample/publication date column"
            ),
        )

    # Numeric / date dominant (legacy higher threshold)
    n_numeric = kinds["numeric"]
    if n_numeric / total >= 0.5:
        return ColumnClassification(
            column_type="measurement", confidence=n_numeric / total,
            n_samples=total, n_matches=n_numeric,
            reason=f"{n_numeric}/{total} cells are numeric : measurement data, not identifiers",
        )

    if n_date / total >= 0.5:
        return ColumnClassification(
            column_type="date", confidence=n_date / total,
            n_samples=total, n_matches=n_date,
            reason=f"{n_date}/{total} cells are dates",
        )

    return ColumnClassification(
        column_type="free_text", confidence=kinds["string-other"] / total,
        n_samples=total, n_matches=kinds["string-other"],
        reason="no dominant pattern found",
    )
