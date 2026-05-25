"""Loaders for real-world registries and corpora.

Per the Ziemann principle: synthetic data hides the bug. The detector's
canonical evidence comes from real registries (HGNC bulk download) and the
real-world corpus of confirmed-corrupted papers (Ziemann 2021 S2 table).

Files expected under `data/raw/`:
- `registries/hgnc_complete_set_<date>.tsv` : HGNC complete set, downloaded
  from https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
- `ziemann_2021_corpus/S2_articles_<date>.xlsx` : supplementary table from
  Ziemann et al. 2021 PLOS Comp Bio, journal.pcbi.1008984.s002
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

_PACKAGE_DIR = Path(__file__).resolve().parent
# Registry data files (HGNC + multispecies xref) ship bundled inside the
# package at src/uncorrupt/_data/. Fall back to the dev-checkout layout
# if running from a source tree where the data has been split out.
DATA_RAW = _PACKAGE_DIR / "_data"
_DEV_FALLBACK = _PACKAGE_DIR.parents[1] / "data" / "raw"
if not DATA_RAW.exists() and _DEV_FALLBACK.exists():
    DATA_RAW = _DEV_FALLBACK


def _latest(glob: str) -> Path:
    matches = sorted(DATA_RAW.glob(glob))
    if not matches:
        raise FileNotFoundError(f"no file matches {glob} under {DATA_RAW}")
    return matches[-1]


@dataclass
class HGNCRegistry:
    """Real HGNC registry view: current symbols + reverse maps from old / alias.

    Lookup is case-insensitive on the input (the maps themselves preserve
    HGNC's canonical case for outputs). HGNC stores aliases in mixed case
    (e.g., POU5F1's alias is `Oct4`, not `OCT4`) but Excel-corrupted
    suggestions are uppercase. Case-folding the lookup recovers the alias
    relationship : empirically caught when calibration treated `OCT4` →
    `POU5F1` rows as 'contradicted' instead of 'corroborated'.
    """
    current_symbols: set[str]
    prev_symbol_to_current: dict[str, str]
    alias_to_current: dict[str, str]
    _current_lower: set[str] = None  # type: ignore[assignment]
    _prev_lower: dict[str, str] = None  # type: ignore[assignment]
    _alias_lower: dict[str, str] = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        # Lower-cased lookup mirrors (built once at construction).
        self._current_lower = {s.lower() for s in self.current_symbols}
        self._prev_lower = {k.lower(): v for k, v in self.prev_symbol_to_current.items()}
        self._alias_lower = {k.lower(): v for k, v in self.alias_to_current.items()}

    def resolve(self, symbol: str) -> str | None:
        """Return current canonical symbol for a possibly-stale input, or None.
        Case-insensitive on the input : HGNC mixed-case aliases (`Oct4`,
        `MGC22487`) match uppercase Excel-output forms (`OCT4`)."""
        if symbol in self.current_symbols:
            return symbol
        if symbol in self.prev_symbol_to_current:
            return self.prev_symbol_to_current[symbol]
        if symbol in self.alias_to_current:
            return self.alias_to_current[symbol]
        # Case-insensitive fallback
        s_lower = symbol.lower()
        if s_lower in self._current_lower:
            # Need to recover canonical case
            for s in self.current_symbols:
                if s.lower() == s_lower:
                    return s
        if s_lower in self._prev_lower:
            return self._prev_lower[s_lower]
        if s_lower in self._alias_lower:
            return self._alias_lower[s_lower]
        return None


def hgnc_snapshot_age_days(path: Path | None = None) -> int | None:
    """Return the age in days of the local HGNC snapshot, or None if
    unparseable. Used to emit drift warnings in Report.sampling_note."""
    if path is None:
        try:
            path = _latest("registries/hgnc_complete_set_*.tsv")
        except FileNotFoundError:
            return None
    from datetime import date as _date
    date_part = path.stem.replace("hgnc_complete_set_", "")
    try:
        snapshot_date = _date.fromisoformat(date_part)
    except ValueError:
        return None
    return (_date.today() - snapshot_date).days


from functools import lru_cache as _lru_cache


@_lru_cache(maxsize=4)
def _load_hgnc_cached(path: Path) -> HGNCRegistry:
    """Internal: parse and cache the HGNC registry per file path. Without
    this the 75 MB TSV gets re-parsed on every detect() call, which made
    multi-sheet xlsx walks 60x slower than necessary (caught in v0.6 real-
    data validation: a 66-sheet file took 5+ minutes when each sheet
    triggered a full HGNC reload)."""
    return _load_hgnc_uncached(path)


def load_hgnc(path: Path | None = None) -> HGNCRegistry:
    """Load the HGNC complete set as a registry view (cached per path)."""
    if path is None:
        path = _latest("registries/hgnc_complete_set_*.tsv")
    return _load_hgnc_cached(path)


def _load_hgnc_uncached(path: Path) -> HGNCRegistry:
    """The actual parser. Separated so the cache decorator stays simple."""
    df = pd.read_csv(path, sep="\t", dtype=object, low_memory=False)
    df = df[df["status"] == "Approved"]

    current_symbols = set(df["symbol"].dropna())
    prev_map: dict[str, str] = {}
    alias_map: dict[str, str] = {}

    for sym, prev, alias in zip(df["symbol"], df["prev_symbol"], df["alias_symbol"], strict=True):
        if pd.isna(sym):
            continue
        if pd.notna(prev):
            for p in str(prev).split("|"):
                p = p.strip()
                if p and p not in current_symbols:
                    prev_map.setdefault(p, sym)
        if pd.notna(alias):
            for a in str(alias).split("|"):
                a = a.strip()
                if a and a not in current_symbols and a not in prev_map:
                    alias_map.setdefault(a, sym)

    return HGNCRegistry(
        current_symbols=current_symbols,
        prev_symbol_to_current=prev_map,
        alias_to_current=alias_map,
    )


DATE_PRONE_PATTERN = re.compile(r"^(MARCH|MARC|SEPT|DEC)(\d+)$")


def derive_date_prone_map(registry: HGNCRegistry) -> dict[str, str]:
    """Return {old_symbol: current_symbol} for genes Excel converts to dates.

    Pattern (Ziemann 2016, 2021; HGNC 2020 rename announcement):
        MARCH<n>, MARC<n>, SEPT<n>, DEC<n>
    """
    return {
        old: new
        for old, new in registry.prev_symbol_to_current.items()
        if DATE_PRONE_PATTERN.fullmatch(old)
    }


@dataclass
class CorpusEntry:
    pmc_id: str
    journal: str
    year: int | None
    species: str
    affected_file_url: str
    confirmed_code: str
    notes: str


def load_ziemann_2021_corpus(path: Path | None = None) -> list[CorpusEntry]:
    """Load the Ziemann 2021 S2 article list as structured entries."""
    if path is None:
        path = _latest("ziemann_2021_corpus/S2_articles_*.xlsx")
    df = pd.read_excel(path, dtype=object)
    entries: list[CorpusEntry] = []
    for row in df.itertuples(index=False):
        year = None
        try:
            year = int(row.Year) if pd.notna(row.Year) else None
        except (ValueError, TypeError):
            pass
        entries.append(CorpusEntry(
            pmc_id=str(row.PMCID),
            journal=str(row.Journal_Name) if pd.notna(row.Journal_Name) else "",
            year=year,
            species=str(row.Species) if pd.notna(row.Species) else "",
            affected_file_url=str(row.Affected_file) if pd.notna(row.Affected_file) else "",
            confirmed_code=str(getattr(row, "Confirmed", "")),
            notes=str(getattr(row, "_6", "") or getattr(row, "Notes_Example", "")),
        ))
    return entries


@dataclass
class Ziemann2016Entry:
    """One row from Ziemann 2016 S1: covers journal articles 2005-2015 + GEO."""
    source: str  # "journal" or "geo"
    file_url: str
    example_corruption: str
    journal_or_db: str
    year: int | None
    confirmed: str
    other_id: str
    pubmed_id: str = ""
    geo_accession: str = ""


def load_ziemann_2016_corpus(path: Path | None = None) -> list[Ziemann2016Entry]:
    """Load Ziemann 2016 S1: 999 journal articles + 228 GEO datasets."""
    if path is None:
        path = _latest("ziemann_2016_corpus/S1_articles_*.xlsx")

    entries: list[Ziemann2016Entry] = []

    journal_df = pd.read_excel(path, sheet_name="JournalArticles", dtype=object)
    for r in journal_df.itertuples(index=False):
        year = None
        try:
            year = int(getattr(r, "Year_Published", None)) if pd.notna(
                getattr(r, "Year_Published", None)) else None
        except (ValueError, TypeError):
            pass
        url = getattr(r, "_0", "") or getattr(r, "Supplementary_file_URL", "")
        entries.append(Ziemann2016Entry(
            source="journal",
            file_url=str(url) if pd.notna(url) else "",
            example_corruption=str(getattr(r, "_1", "") or
                                   getattr(r, "Example_Gene_Name_Conversion", "")),
            journal_or_db=str(getattr(r, "Journal", "")),
            year=year,
            confirmed=str(getattr(r, "Confirmed", "")),
            other_id=str(getattr(r, "_5", "") or
                         getattr(r, "Other_accession_numbers_identifying_info", "")),
            pubmed_id=str(getattr(r, "PubmedID", "")),
        ))

    geo_df = pd.read_excel(path, sheet_name="GEOdatasets", dtype=object)
    for r in geo_df.itertuples(index=False):
        year = None
        try:
            year = int(getattr(r, "Year_of_Release", None)) if pd.notna(
                getattr(r, "Year_of_Release", None)) else None
        except (ValueError, TypeError):
            pass
        url = getattr(r, "_0", "") or getattr(r, "Supplementary_file_URL", "")
        entries.append(Ziemann2016Entry(
            source="geo",
            file_url=str(url) if pd.notna(url) else "",
            example_corruption=str(getattr(r, "_1", "") or
                                   getattr(r, "Example_Gene_Name_Conversion", "")),
            journal_or_db=str(getattr(r, "Database", "")),
            year=year,
            confirmed=str(getattr(r, "Confirmed", "")),
            other_id=str(getattr(r, "_5", "") or
                         getattr(r, "Other_accession_numbers_identifying_info", "")),
            geo_accession=str(getattr(r, "Dataset_accession_number", "")),
        ))

    return entries
