"""Refresh the multi-species xref tables (MGI, ZFIN, FlyBase, WormBase).

Downloads each registry's bulk file, parses the columns documented in the
research review, writes one TSV per species under
`data/raw/registries/` with a unified schema:

    species  external_id  current_symbol  prev_or_alias

External ID prefix indicates source (`MGI:`, `ZDB-GENE-`, `FBgn`, `WB:`,
`NCBI_Gene:`, `ENSMUSG`, `ENSDARG`, `FBgn`, etc.). The xref index loader
in `src/uncorrupt/xref_lookup.py` ingests all of `*.xref.tsv` under
`data/raw/registries/` at startup and merges them into one mapping.

All four sources are CC-BY 4.0 (or equivalent : RGD pending).

Run:
   python scripts/refresh_multispecies_xref.py              # all species
   python scripts/refresh_multispecies_xref.py --only mgi   # one species
   python scripts/refresh_multispecies_xref.py --check-only # offline check
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import sys
import urllib.error
import urllib.request
from datetime import UTC, date, datetime
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[1]
REG_DIR = ROOT / "data/raw/registries"
MANIFEST = REG_DIR / "MULTISPECIES_MANIFEST.jsonl"

SOURCES: dict[str, dict] = {
    "mgi": {
        "name": "MGI (Mus musculus)",
        "license": "CC-BY 4.0",
        "files": [
            ("MGI_EntrezGene.rpt",
             "https://www.informatics.jax.org/downloads/reports/MGI_EntrezGene.rpt"),
            ("MRK_ENSEMBL.rpt",
             "https://www.informatics.jax.org/downloads/reports/MRK_ENSEMBL.rpt"),
        ],
    },
    "zfin": {
        "name": "ZFIN (Danio rerio)",
        "license": "CC-BY 4.0",
        "files": [
            ("zfin_gene.txt", "https://zfin.org/downloads/gene.txt"),
            ("zfin_aliases.txt", "https://zfin.org/downloads/aliases.txt"),
            ("zfin_ensembl_1_to_1.txt",
             "https://zfin.org/downloads/ensembl_1_to_1.txt"),
        ],
    },
    "flybase": {
        "name": "FlyBase (Drosophila melanogaster)",
        "license": "CC-BY 4.0",
        "files": [
            # Alliance BGI fallback (FlyBase FTP often returns 404 from
            # automated mirrors).
            ("BGI_FB.json.gz",
             "https://fms.alliancegenome.org/download/BGI_FB.json.gz"),
        ],
    },
    "wormbase": {
        "name": "WormBase (Caenorhabditis elegans)",
        "license": "CC-BY 4.0",
        "files": [
            ("BGI_WB.json.gz",
             "https://fms.alliancegenome.org/download/BGI_WB.json.gz"),
        ],
    },
    "rat": {
        # NCBI Rattus_norvegicus.gene_info.gz : US public domain.
        # Cleanest license for redistribution alongside Apache-2.0 code.
        # RGD direct is CC-BY 4.0 (verified at
        # https://rgd.mcw.edu/wg/home/disclaimer) but requires
        # attribution; NCBI requires only courtesy attribution.
        "name": "NCBI Rattus norvegicus gene_info (rat)",
        "license": "US public domain (NCBI)",
        "files": [
            ("Rattus_norvegicus.gene_info.gz",
             "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/"
             "Mammalia/Rattus_norvegicus.gene_info.gz"),
        ],
    },
}


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def _download(url: str, dest: Path) -> None:
    print(f"  fetching {url}", flush=True)
    try:
        with urllib.request.urlopen(url, timeout=60) as r: # noqa: S310
            dest.write_bytes(r.read())
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"HTTP {e.code} on {url}") from e
    except urllib.error.URLError as e:
        raise RuntimeError(f"URL error on {url}: {e.reason}") from e


def _parse_mgi(reg_dir: Path) -> list[dict]:
    """Yield {species, external_id, current_symbol, prev_or_alias} rows."""
    rows: list[dict] = []
    entrez = reg_dir / "MGI_EntrezGene.rpt"
    ensembl = reg_dir / "MRK_ENSEMBL.rpt"
    if not entrez.exists():
        print(f"  SKIP MGI: {entrez} missing")
        return rows
    # MGI_EntrezGene.rpt: tab-separated, no header
    # col0=MGI_ID, col1=symbol, col2=status (O/W), col8=Entrez, col9=aliases
    with entrez.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10 or parts[2] != "O":
                continue
            mgi_id, symbol = parts[0], parts[1]
            entrez_id = parts[8].strip()
            aliases_raw = parts[9].strip()
            rows.append({
                "species": "mouse", "external_id": mgi_id,
                "current_symbol": symbol, "prev_or_alias": symbol,
            })
            if entrez_id:
                rows.append({
                    "species": "mouse", "external_id": f"NCBI_Gene:{entrez_id}",
                    "current_symbol": symbol, "prev_or_alias": symbol,
                })
            for alias in aliases_raw.split("|"):
                alias = alias.strip()
                if alias and alias != symbol:
                    rows.append({
                        "species": "mouse", "external_id": mgi_id,
                        "current_symbol": symbol, "prev_or_alias": alias,
                    })
    # MRK_ENSEMBL.rpt: col0=MGI_ID, col1=symbol, col5=Ensembl gene ID
    if ensembl.exists():
        with ensembl.open() as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue
                mgi_id, symbol, ens_id = parts[0], parts[1], parts[5].strip()
                if ens_id:
                    rows.append({
                        "species": "mouse", "external_id": ens_id,
                        "current_symbol": symbol, "prev_or_alias": symbol,
                    })
    return rows


def _parse_zfin(reg_dir: Path) -> list[dict]:
    rows: list[dict] = []
    gene = reg_dir / "zfin_gene.txt"
    aliases = reg_dir / "zfin_aliases.txt"
    ens = reg_dir / "zfin_ensembl_1_to_1.txt"
    if not gene.exists():
        return rows
    with gene.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4 or not parts[0].startswith("ZDB-GENE"):
                continue
            zdb_id, _so, symbol, entrez = parts[0], parts[1], parts[2], parts[3]
            rows.append({
                "species": "zebrafish", "external_id": zdb_id,
                "current_symbol": symbol, "prev_or_alias": symbol,
            })
            if entrez.strip():
                rows.append({
                    "species": "zebrafish",
                    "external_id": f"NCBI_Gene:{entrez.strip()}",
                    "current_symbol": symbol, "prev_or_alias": symbol,
                })
    if aliases.exists():
        with aliases.open() as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5 or not parts[0].startswith("ZDB-GENE"):
                    continue
                zdb_id, _name, symbol, alias = parts[0], parts[1], parts[2], parts[3]
                if alias and alias != symbol:
                    rows.append({
                        "species": "zebrafish", "external_id": zdb_id,
                        "current_symbol": symbol, "prev_or_alias": alias,
                    })
    if ens.exists():
        with ens.open() as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4 or not parts[0].startswith("ZDB-GENE"):
                    continue
                _zdb_id, _so, symbol, ens_id = parts
                if ens_id.strip():
                    rows.append({
                        "species": "zebrafish",
                        "external_id": ens_id.strip(),
                        "current_symbol": symbol, "prev_or_alias": symbol,
                    })
    return rows


def _parse_ncbi_gene_info(reg_dir: Path, fname: str, species: str) -> list[dict]:
    """Parse NCBI gene_info.gz. Columns (per
    https://ftp.ncbi.nlm.nih.gov/gene/DATA/README):
      0  tax_id
      1  GeneID (Entrez)
      2  Symbol (current)
      3  LocusTag
      4  Synonyms (pipe-separated, '-' if absent)
      5  dbXrefs (pipe-separated, e.g. 'RGD:1306236|HGNC:HGNC:...')
      ...
    Header is `#tax_id\\tGeneID\\t...` : strip the leading `#`.
    """
    path = reg_dir / fname
    rows: list[dict] = []
    if not path.exists():
        return rows
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        header = next(fh, None)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            _tax, entrez, symbol, _locus, synonyms, dbxrefs = parts[:6]
            symbol = symbol.strip()
            if not symbol:
                continue
            rows.append({
                "species": species,
                "external_id": f"NCBI_Gene:{entrez.strip()}",
                "current_symbol": symbol,
                "prev_or_alias": symbol,
            })
            if synonyms and synonyms != "-":
                for alias in synonyms.split("|"):
                    alias = alias.strip()
                    if alias and alias != symbol:
                        rows.append({
                            "species": species,
                            "external_id": f"NCBI_Gene:{entrez.strip()}",
                            "current_symbol": symbol,
                            "prev_or_alias": alias,
                        })
            if dbxrefs and dbxrefs != "-":
                for xref in dbxrefs.split("|"):
                    xref = xref.strip()
                    # RGD:N → RGD reference
                    if xref.startswith("RGD:") or xref.startswith("Ensembl:"):
                        rows.append({
                            "species": species,
                            "external_id": xref,
                            "current_symbol": symbol,
                            "prev_or_alias": symbol,
                        })
    return rows


def _parse_alliance_bgi(reg_dir: Path, fname: str, species: str) -> list[dict]:
    """Parse Alliance of Genome Resources BGI JSON for FlyBase or WormBase."""
    path = reg_dir / fname
    rows: list[dict] = []
    if not path.exists():
        return rows
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        data = json.load(fh)
    for entry in data.get("data", []):
        symbol = entry.get("symbol")
        if not symbol:
            continue
        bge = entry.get("basicGeneticEntity", {})
        primary = bge.get("primaryId", "").replace("WB:", "").replace("FB:", "")
        if primary:
            rows.append({
                "species": species, "external_id": primary,
                "current_symbol": symbol, "prev_or_alias": symbol,
            })
        for alias in bge.get("synonyms") or []:
            if alias and alias != symbol:
                rows.append({
                    "species": species, "external_id": primary,
                    "current_symbol": symbol, "prev_or_alias": alias,
                })
        for sec in bge.get("secondaryIds") or []:
            sid = sec.replace("WB:", "").replace("FB:", "")
            if sid:
                rows.append({
                    "species": species, "external_id": sid,
                    "current_symbol": symbol, "prev_or_alias": symbol,
                })
        for xref in bge.get("crossReferences") or []:
            xid = xref.get("id", "") if isinstance(xref, dict) else ""
            if xid.startswith("NCBI_Gene:") or xid.startswith("ENSEMBL:"):
                rows.append({
                    "species": species, "external_id": xid,
                    "current_symbol": symbol, "prev_or_alias": symbol,
                })
    return rows


def _write_xref(rows: Iterable[dict], path: Path) -> int:
    """Write deduplicated rows to a unified TSV. Returns row count."""
    seen: set[tuple] = set()
    written = 0
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["species", "external_id", "current_symbol",
                             "prev_or_alias"],
            delimiter="\t",
        )
        writer.writeheader()
        for r in rows:
            key = (r["species"], r["external_id"], r["prev_or_alias"])
            if key in seen:
                continue
            seen.add(key)
            writer.writerow(r)
            written += 1
    return written


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--only", choices=sorted(SOURCES) + ["none"], default=None)
    p.add_argument("--check-only", action="store_true",
                   help="Don't download; just parse what's already on disk.")
    args = p.parse_args()

    REG_DIR.mkdir(parents=True, exist_ok=True)
    chosen = [args.only] if args.only else list(SOURCES)
    total = 0
    for species in chosen:
        meta = SOURCES[species]
        print(f"\n== {meta['name']} ({meta['license']}) ==")
        if not args.check_only:
            for fname, url in meta["files"]:
                dest = REG_DIR / fname
                try:
                    _download(url, dest)
                except Exception as exc:
                    print(f"  WARN: {exc}", file=sys.stderr)
        # Parse
        if species == "mgi":
            rows = _parse_mgi(REG_DIR)
        elif species == "zfin":
            rows = _parse_zfin(REG_DIR)
        elif species == "flybase":
            rows = _parse_alliance_bgi(REG_DIR, "BGI_FB.json.gz", "fly")
        elif species == "wormbase":
            rows = _parse_alliance_bgi(REG_DIR, "BGI_WB.json.gz", "worm")
        elif species == "rat":
            rows = _parse_ncbi_gene_info(
                REG_DIR, "Rattus_norvegicus.gene_info.gz", "rat",
            )
        else:
            rows = []
        out_path = REG_DIR / f"{species}.xref.tsv"
        n = _write_xref(rows, out_path)
        total += n
        print(f"  wrote {n} rows to {out_path}")
        if not args.check_only:
            entry = {
                "species": species,
                "name": meta["name"],
                "rows": n,
                "downloaded_at": datetime.now(UTC).isoformat(),
                "license": meta["license"],
            }
            with MANIFEST.open("a") as fh:
                fh.write(json.dumps(entry) + "\n")
    print(f"\nTotal rows across species: {total}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
