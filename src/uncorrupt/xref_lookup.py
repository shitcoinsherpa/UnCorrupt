"""Cross-reference lookup: given any external gene ID, return the HGNC symbol.

Most gene tables include external identifiers (UniProt, RefSeq, Ensembl,
Entrez, HGNC) in cells adjacent to the gene-symbol cell. When the
gene-symbol cell looks corrupted, we can use the external ID in the same
row to look up what the gene SHOULD be : that's our automatic ground truth.

Build the index from the HGNC bulk download (already cached). Patterns
covered:
- Entrez Gene ID (digits, e.g. "1", "503538")
- Ensembl gene/transcript/protein (ENSG/ENST/ENSP + digits)
- RefSeq mRNA/protein (NM_/NR_/NP_ + digits)
- UniProt accession (e.g. "P04217", "Q9NQ94")
- HGNC ID (HGNC:N)
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from functools import cache
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    pass

# Pre-compiled patterns for ID-shape detection in arbitrary cell strings
_ID_PATTERNS = {
    # Ensembl: human (ENSG/ENST/ENSP/ENSR) + mouse (ENSMUSG/ENSMUST/...) +
    # zebrafish (ENSDARG/ENSDART/...). The G/T/P/R letter encodes
    # gene/transcript/protein/regulatory and is sometimes followed by
    # additional species letters before the digits.
    "ensembl": re.compile(r"^ENS[A-Z]{0,4}[GTPR]\d{6,}\.?\d*$"),
    "refseq": re.compile(r"^N[MRPC]_\d+\.?\d*$"),
    "uniprot": re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"),
    "hgnc_id": re.compile(r"^HGNC:\d+$"),
    "mgi_id": re.compile(r"^MGI:\d+$"),
    "wormbase_id": re.compile(r"^WBGene\d+$"),
    "flybase_id": re.compile(r"^FBgn\d+$"),
    "zfin_id": re.compile(r"^ZDB-GENE-\d+-\d+$"),
    # Entrez IDs are integers, but small integers (1-999) are too commonly
    # row numbers, ranks, chromosome numbers, etc. to be safe Entrez matches.
    # Require 4+ digits : covers the vast majority of real Entrez IDs while
    # avoiding false positives on rank/index columns.
    "entrez": re.compile(r"^\d{4,9}$"),
}


@dataclass
class XrefIndex:
    """Reverse-lookup index: external ID → HGNC current symbol."""
    by_entrez: dict[str, str]
    by_ensembl_gene: dict[str, str]
    by_refseq: dict[str, str]
    by_uniprot: dict[str, str]
    by_hgnc_id: dict[str, str]

    def lookup(self, value: str) -> str | None:
        """Return the resolved gene symbol for any matching external ID, or None.

        Multi-species IDs (MGI/WBGene/FBgn/ZDB-GENE-...) resolve to their
        own organism's current symbol (mouse Mecp2, fly sevenless, etc.),
        not to a forced human ortholog. This is correct for the
        calibration walk: if a row has both a mouse Entrez ID and a
        symbol like "Sept2", the correct resolution is Septin2 (mouse),
        which would equivalent-match SEPT2 (human) under HGNC alias
        expansion in the calibration script.
        """
        v = value.strip()
        if not v:
            return None
        for kind, pat in _ID_PATTERNS.items():
            if pat.fullmatch(v):
                if kind == "ensembl":
                    base = v.split(".", 1)[0]
                    return self.by_ensembl_gene.get(base)
                if kind == "refseq":
                    base = v.split(".", 1)[0]
                    return self.by_refseq.get(base)
                if kind == "uniprot":
                    return self.by_uniprot.get(v)
                if kind == "hgnc_id":
                    return self.by_hgnc_id.get(v)
                if kind == "entrez":
                    return self.by_entrez.get(v)
                if kind in ("mgi_id", "wormbase_id", "flybase_id", "zfin_id"):
                    # Multi-species IDs share the by_ensembl_gene dict
                    # (it's used as a general "non-bare gene ID" bucket
                    # populated by _ingest_multispecies).
                    return self.by_ensembl_gene.get(v)
        return None

    def identify(self, value: str) -> str | None:
        """Return the kind of ID (entrez/ensembl/refseq/uniprot/hgnc_id) or None."""
        v = value.strip()
        for kind, pat in _ID_PATTERNS.items():
            if pat.fullmatch(v):
                return kind
        return None


def _parse_pipe_separated(s: str) -> list[str]:
    """HGNC stores some xrefs as pipe-separated (e.g. multiple RefSeq IDs)."""
    if pd.isna(s):
        return []
    return [t.strip() for t in str(s).split("|") if t.strip()]


@cache
def load_xref_index(path: Path | None = None) -> XrefIndex:
    """Load the HGNC bulk file and build reverse-lookup dicts."""
    if path is None:
        from .corpus import _latest
        path = _latest("registries/hgnc_complete_set_*.tsv")
    df = pd.read_csv(path, sep="\t", dtype=object, low_memory=False)
    df = df[df["status"] == "Approved"]

    by_entrez: dict[str, str] = {}
    by_ensembl_gene: dict[str, str] = {}
    by_refseq: dict[str, str] = {}
    by_uniprot: dict[str, str] = {}
    by_hgnc_id: dict[str, str] = {}

    for row in df.itertuples(index=False):
        sym = row.symbol
        if pd.isna(sym):
            continue
        sym = str(sym)

        ent = getattr(row, "entrez_id", None)
        if pd.notna(ent):
            by_entrez.setdefault(str(ent).strip(), sym)
        eg = getattr(row, "ensembl_gene_id", None)
        if pd.notna(eg):
            by_ensembl_gene.setdefault(str(eg).strip(), sym)
        rs = getattr(row, "refseq_accession", None)
        if pd.notna(rs):
            for r in _parse_pipe_separated(rs):
                by_refseq.setdefault(r.split(".", 1)[0], sym)
        up = getattr(row, "uniprot_ids", None)
        if pd.notna(up):
            for u in _parse_pipe_separated(up):
                by_uniprot.setdefault(u, sym)
        hid = getattr(row, "hgnc_id", None)
        if pd.notna(hid):
            by_hgnc_id.setdefault(str(hid).strip(), sym)

    # multi-species xref expansion (MGI mouse, ZFIN zebrafish,
    # FlyBase Drosophila, WormBase C. elegans). All CC-BY 4.0; populated by
    # `scripts/refresh_multispecies_xref.py`. Each species file lives under
    # `data/raw/registries/{species}.xref.tsv` with columns:
    #     species  external_id  current_symbol  prev_or_alias
    # We treat any `external_id` as a resolvable lookup key. The species's
    # current_symbol is the resolution target (we do NOT collapse to HGNC
    # because non-human IDs don't always have human orthologs).
    _ingest_multispecies(
        by_entrez=by_entrez,
        by_ensembl_gene=by_ensembl_gene,
        by_hgnc_id=by_hgnc_id,  # noop : multi-species don't add HGNC IDs
        registries_dir=path.parent,
    )

    return XrefIndex(
        by_entrez=by_entrez,
        by_ensembl_gene=by_ensembl_gene,
        by_refseq=by_refseq,
        by_uniprot=by_uniprot,
        by_hgnc_id=by_hgnc_id,
    )


def _ingest_multispecies(
    *,
    by_entrez: dict[str, str],
    by_ensembl_gene: dict[str, str],
    by_hgnc_id: dict[str, str],  # currently unused : placeholder
    registries_dir: Path,
) -> None:
    """Merge `*.xref.tsv` species files into the in-memory dicts.

    Conflict policy: HGNC entries already in `by_entrez` / `by_ensembl_gene`
    are preserved (HGNC is the human canonical). Multi-species IDs add
    coverage; they never overwrite an existing key.

    Each row is `species\\texternal_id\\tcurrent_symbol\\tprev_or_alias`.
    External IDs may carry prefixes (`MGI:`, `NCBI_Gene:`, `ENSMUSG`,
    `ENSDARG`, `FBgn`, `WBGene`, `ZDB-GENE-`). We strip `NCBI_Gene:` to
    canonicalize Entrez integer keys; other prefixes are kept verbatim.
    """
    for xref_path in sorted(registries_dir.glob("*.xref.tsv")):
        try:
            with xref_path.open() as fh:
                next(fh, None)  # header
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 4:
                        continue
                    _species, ext_id, symbol, _alias = parts
                    ext_id = ext_id.strip()
                    symbol = symbol.strip()
                    if not ext_id or not symbol:
                        continue
                    # NCBI_Gene:N → bare integer, goes into by_entrez
                    if ext_id.startswith("NCBI_Gene:"):
                        key = ext_id[len("NCBI_Gene:"):]
                        by_entrez.setdefault(key, symbol)
                        continue
                    # Ensembl gene IDs (ENSG/ENSMUSG/ENSDARG/etc.)
                    if ext_id.startswith("ENS"):
                        by_ensembl_gene.setdefault(ext_id, symbol)
                        continue
                    # Species-specific bare IDs (MGI:..., WBGene..., FBgn..., ZDB-GENE-...)
                    # : stash under by_ensembl_gene for now (it serves as a
                    # general "non-human gene ID → symbol" bucket).
                    by_ensembl_gene.setdefault(ext_id, symbol)
        except Exception:
            continue  # malformed file → skip; don't break the whole load


def lookup_row_genes(row_values: list, xref: XrefIndex) -> dict[str, str]:
    """For a row, return {cell_value: resolved_symbol} for every cell that
    looks like a known external ID and resolves to a gene."""
    found: dict[str, str] = {}
    for v in row_values:
        if not isinstance(v, str):
            continue
        sym = xref.lookup(v)
        if sym is not None:
            found[v] = sym
    return found
