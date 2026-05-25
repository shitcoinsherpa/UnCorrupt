"""Identifier registries: patterns and live-loaded data from HGNC.

Symbol sets and rename maps are derived from the real HGNC complete-set bulk
download cached under `data/raw/registries/`. Regex patterns are hand-coded
because they describe identifier *shape*, not *identity*.

Per the project's Pillar 4 (FAIR data) and the real-data rule: the hardcoded
SEPT/MARCH lists that previously lived here were copied from secondary
literature and disagreed with HGNC in non-trivial ways (e.g. SEPT13 maps to
SEPTIN7P2, not SEPTIN13; SEPT15 has no current mapping). HGNC is the truth.

If the bulk download is missing, import fails with a clear error pointing at
the fetch script.
"""
from __future__ import annotations

import re
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .corpus import HGNCRegistry

GENE_LIKE_PATTERN = re.compile(r"^[A-Z][A-Z0-9-]{1,14}$")
# Case-insensitive: TSV/CSV exports may write lowercase `2310009e13`.
RIKEN_PATTERN = re.compile(r"^\d{4}\d{3}[A-Z]\d{2}$", re.IGNORECASE)
EXPONENT_PATTERN = re.compile(r"^\d+[eE]\d+$")

IDENTIFIER_HEADER_HINTS: set[str] = {
    "gene", "genes", "symbol", "symbols", "gene_symbol", "gene_name",
    "id", "identifier", "accession", "name", "gene_id",
    "protein", "proteins", "transcript", "mrna", "locus",
    "entrez", "entrezid", "uniprot", "ensembl", "refseq", "probe",
    "ncbi", "hgnc", "mgi", "ortholog", "geneid",
}

NON_IDENTIFIER_HEADER_HINTS: set[str] = {
    "date", "time", "datetime", "timestamp", "year", "month", "day",
    "age", "count", "score", "value", "level", "duration",
}


def _load() -> HGNCRegistry:
    from .corpus import load_hgnc
    try:
        return load_hgnc()
    except FileNotFoundError as exc:
        raise RuntimeError(
            "HGNC bulk download not found under data/raw/registries/. "
            "Run `scripts/fetch_data.sh` to download it, or place "
            "`hgnc_complete_set_<date>.tsv` there manually."
        ) from exc


_HGNC = _load()

HGNC_NEW_SYMBOLS: set[str] = _HGNC.current_symbols
HGNC_RENAME_MAP: dict[str, str] = _HGNC.prev_symbol_to_current
HGNC_ALIAS_MAP: dict[str, str] = _HGNC.alias_to_current


def get_hgnc() -> HGNCRegistry:
    """Access the loaded HGNC registry."""
    return _HGNC
