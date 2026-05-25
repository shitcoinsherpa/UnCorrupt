"""Detect Excel-style corruption in tabular data.

Algorithm:
1. For each column, decide if it looks like an identifier column based on:
   header hint, fraction of gene-like or RIKEN-like values, or presence of any
   string in the HGNC rename map (old or new).
2. In identifier-shaped columns, flag date and float values as suspicious.
3. For date values, reverse-map to candidate gene symbols with a confidence
   score that drops when the mapping is ambiguous (Mar-01 / Mar-02).
4. For float values that look like coerced exponent strings, mark as
   irreversibly lossy.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from datetime import date, datetime, timedelta
from datetime import time as _time
from functools import cache as _cache
from typing import Any

import pandas as pd

# Excel date-serial range. Day 0 = 1899-12-30 (Windows Excel epoch with the
# 1900-leap-year quirk). [20000, 60000] covers 1954-09-25 to 2064-04-12 : wide
# enough to encompass any plausible gene-name date corruption.
_EXCEL_EPOCH = date(1899, 12, 30)
_SERIAL_MIN = 20000
_SERIAL_MAX = 60000


def _serial_to_date(n: int) -> date | None:
    if not (_SERIAL_MIN <= n <= _SERIAL_MAX):
        return None
    try:
        return _EXCEL_EPOCH + timedelta(days=n)
    except (OverflowError, ValueError):
        return None

import contextlib

from .homoglyph import looks_like_gene_after_repair as _homoglyph_repair
from .registries import (
    GENE_LIKE_PATTERN,
    HGNC_NEW_SYMBOLS,
    HGNC_RENAME_MAP,
    IDENTIFIER_HEADER_HINTS,
    NON_IDENTIFIER_HEADER_HINTS,
    RIKEN_PATTERN,
)


@dataclass
class Suspicion:
    column: str
    row: int
    value: Any
    kind: str
    suggestion: str | None
    confidence: float
    reason: str
    sheet: str | None = None
    # Row-xref outcome : set by `_apply_row_context_boost` when the
    # detector finds (or fails to find) an external ID in the same row
    # that resolves to the suggested gene. Values:
    #   "corroborated": external ID resolves to the suggested gene
    #   "contradicted": external ID resolves to a *different* gene
    #   "absent": row has no resolvable external IDs
    #   None: xref boost wasn't run (row_context_boost=False)
    xref_status: str | None = None


@dataclass
class Report:
    rows_scanned: int
    columns_scanned: int
    identifier_columns: list[str] = field(default_factory=list)
    suspicions: list[Suspicion] = field(default_factory=list)
    # If the input dataframe exceeded the row threshold and we sampled,
    # this records the sampling decision: total_rows, sampled_rows, method.
    # `None` means an exhaustive scan was performed.
    sampling_note: str | None = None


# Month-abbreviation alternation. Comprehensive across the Latin-alphabet
# locales Excel ships with: English, German, French, Spanish, Portuguese,
# Italian, Dutch, Finnish, Swedish, Polish, Czech, Romanian. CJK locales
# render dates in script (1月, 一月) which doesn't collide with Latin gene
# symbols. Cyrillic Russian (янв, фев, …) likewise stays in Cyrillic.
#
# Conflict resolution: where two locales share an abbreviation but map to
# different months (e.g., Finnish "mar" = November, English "mar" = March),
# the English form wins. The year-suffix format path catches cells where
# locale-specific decoding matters.
_MONTH_ALTERNATION = (
    # English
    r"jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec|"
    # Romance: Spanish/Italian/Portuguese/French
    r"ago|dic|ene|abr|gen|mag|lug|set|ott|giu|"
    r"avr|mai|aoa|fev|f[eé]v|d[eé]c|jui|"
    # German
    r"mrz|m[aä]r|okt|dez|"
    # Dutch
    r"mei|mrt|"
    # Finnish (tam=Jan, hel=Feb, maa=Mar, huh=Apr, tou=May,
    #         kes=Jun, hei=Jul, elo=Aug, syy=Sep, lok=Oct, jou=Dec)
    r"tam|hel|maa|huh|tou|kes|hei|elo|syy|lok|jou|"
    # Swedish
    r"maj|"
    # Polish (sty=Jan, lut=Feb, kwi=Apr, cze=Jun, lip=Jul,
    #        sie=Aug, wrz=Sep, paz=Oct, gru=Dec; lis=Nov)
    r"sty|lut|kwi|cze|lip|sie|wrz|paz|gru|lis|"
    # Czech (no diacritics: led=Jan, uno=Feb, bre=Mar, dub=Apr,
    #       kve=May, cvn=Jun, cvc=Jul, srp=Aug, zar=Sep, rij=Oct, pro=Dec)
    r"led|uno|bre|dub|kve|cvn|cvc|srp|zar|rij|pro|"
    # Romanian
    r"ian|iun|iul|noi"
)

# Separator class: `-`, `/`, OR `.` (German-locale Excel uses period:
# `Mar.2`, `07.Mar`).
_DATE_SEP = r"[-/.]"

# Round-tripped Excel time strings: when a `1:3` cell becomes datetime.time
# in memory, then exports to CSV as "01:03:00" or "1:03". Pyle et al. PMC5617173
# documents this. We accept HH:MM and HH:MM:SS forms with reasonable bounds
# (hour <= 24, minute/second < 60). Strict regex : refuses ambiguous things
# like `1:3:5:7`.
_TIME_STRING_RE = re.compile(
    r"^([0-1]?\d|2[0-4]):([0-5]\d)(?::([0-5]\d))?$"
)

_DATE_STRING_RE = re.compile(
    r"^(?:"
    r"\d{4}-\d{1,2}-\d{1,2}(?:[ T]\d{1,2}:\d{1,2}(?::\d{1,2})?)?"      # ISO 2012-03-03
    rf"|\d{{1,2}}{_DATE_SEP}\d{{1,2}}{_DATE_SEP}\d{{2,4}}"             # 06/03/14, 06.03.14
    rf"|\d{{1,2}}{_DATE_SEP}(?:{_MONTH_ALTERNATION})(?:{_DATE_SEP}\d{{2,4}})?"  # 06-Mar / 06.Mar
    rf"|(?:{_MONTH_ALTERNATION}){_DATE_SEP}\d{{1,2}}(?:{_DATE_SEP}\d{{2,4}})?"  # Mar-06 / Mar.06
    r")$",
    re.IGNORECASE,
)
_MONTH_ABBR = {
    # === English ===
    "jan": 1, "feb": 2, "mar": 3, "apr": 4, "may": 5, "jun": 6,
    "jul": 7, "aug": 8, "sep": 9, "oct": 10, "nov": 11, "dec": 12,
    # === Romance ===
    "ago": 8,    # Spanish/Italian/Portuguese August
    "dic": 12,   # Spanish/Italian December
    "ene": 1,    # Spanish January
    "abr": 4,    # Spanish/Portuguese April
    "gen": 1,    # Italian January (gennaio)
    "mag": 5,    # Italian May (maggio)
    "lug": 7,    # Italian July (luglio)
    "set": 9,    # Italian September (settembre)
    "ott": 10,   # Italian October (ottobre)
    "giu": 6,    # Italian June (giugno)
    "avr": 4,    # French April (avril)
    "mai": 5,    # French/German May (mai)
    "aoa": 8,    # French August (août, ASCII fallback)
    "fev": 2,    # French/Portuguese February (février/fevereiro)
    "fév": 2,    # French February with accent
    "déc": 12,   # French December (décembre)
    "jui": 6,    # French June (juin) : ambiguous with July; June chosen
    # === German ===
    "mrz": 3,    # German March (März abbreviated to Mrz)
    "mär": 3,    # German March (with umlaut)
    "okt": 10,   # German/Dutch October
    "dez": 12,   # German December (Dezember)
    # === Dutch ===
    "mei": 5,    # Dutch May
    "mrt": 3,    # Dutch March (maart)
    # === Finnish ===
    "tam": 1,    # tammikuu
    "hel": 2,    # helmikuu
    "maa": 3,    # maaliskuu
    "huh": 4,    # huhtikuu
    "tou": 5,    # toukokuu
    "kes": 6,    # kesäkuu
    "hei": 7,    # heinäkuu
    "elo": 8,    # elokuu
    "syy": 9,    # syyskuu
    "lok": 10,   # lokakuu
    "jou": 12,   # joulukuu (NB: Finnish "mar" = marraskuu = November, but
                  # collides with English March; English wins.)
    # === Swedish ===
    "maj": 5,    # Swedish/Polish May
    # === Polish ===
    "sty": 1,    # styczeń
    "lut": 2,    # luty
    "kwi": 4,    # kwiecień
    "cze": 6,    # czerwiec
    "lip": 7,    # lipiec
    "sie": 8,    # sierpień
    "wrz": 9,    # wrzesień
    "paz": 10,   # październik (ASCII)
    "lis": 11,   # listopad
    "gru": 12,   # grudzień
    # === Czech ===
    "led": 1,    # leden
    "uno": 2,    # únor
    "bre": 3,    # březen
    "dub": 4,    # duben
    "kve": 5,    # květen
    "cvn": 6,    # červen
    "cvc": 7,    # červenec
    "srp": 8,    # srpen
    "zar": 9,    # září
    "rij": 10,   # říjen
    "pro": 12,   # prosinec
    # === Romanian ===
    "ian": 1,    # ianuarie
    "iun": 6,    # iunie
    "iul": 7,    # iulie
    "noi": 11,   # noiembrie
}


@_cache
def _multispecies_symbols_lower() -> frozenset[str]:
    """Lower-cased set of every multi-species gene symbol from the xref
    index. Loaded once per process. Used to suppress date-parse
    interpretation on cells that match a known cross-species gene symbol
    in a non-English locale (e.g. fly `Lis-1` decoded as Polish
    `listopad` = November; fly `mei-9` decoded as Dutch `mei` = May).

    Excluded from the suppression set: symbols that themselves look
    like the date-corruption shape (`Sep-2`, `Mar-15`, `Oct-3`). The
    multi-species registries include these as alias variants : but
    treating them as "real symbols" would defeat the detector's primary
    mission of catching exactly this shape.

    Empty set when the multi-species xref tables haven't been fetched
    yet (run `scripts/refresh_multispecies_xref.py`).
    """
    from .corpus import DATA_RAW
    syms: set[str] = set()
    reg_dir = DATA_RAW / "registries"

    # English-month-abbrev date shapes are EXACTLY the Excel-corruption
    # output. We want those to parse as dates, not be suppressed by the
    # cross-species guard : even when the same string appears as a
    # symbol synonym in MGI/ZFIN/FlyBase/WormBase (it does, ironically,
    # because the registries record corruption-shaped synonyms). Locale
    # abbreviations (Polish `lis`, Dutch `mei`) are NOT in this English
    # exclusion list : those CAN be valid gene symbols (Lis-1, mei-9 in
    # fly), and we want them suppressed from date parsing.
    _EN_DATE_SHAPE = re.compile(
        r"^(?:\d{1,2}[-/.]"
        r"(?:jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)"
        r"|(?:jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)[-/.]"
        r"\d{1,2})(?:[-/.]\d{2,4})?$",
        re.IGNORECASE,
    )

    def _add(sym: str) -> None:
        if not sym:
            return
        # Reject only the English-month-abbrev date shape (real Excel
        # corruption output). Polish/Dutch locale abbreviations stay
        # in the suppression set.
        if _EN_DATE_SHAPE.fullmatch(sym):
            return
        syms.add(sym.lower())

    for path in reg_dir.glob("*.xref.tsv"):
        try:
            with path.open() as fh:
                next(fh, None)  # header
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 4:
                        continue
                    _species, _ext, current_symbol, alias = parts
                    _add(current_symbol.strip())
                    _add(alias.strip())
        except Exception:
            continue
    return frozenset(syms)


def _parse_date_string(s: str) -> list[date]:
    """Return possible date interpretations of a string. Empty if no plausible parse.

    Handles ambiguous DD/MM vs MM/DD by returning *both* when both are plausible.
    Long-form strings (`September 2nd, 2024`; `2 septembre 2014`) are delegated
    to dateutil.parser as a final fallback : bounded by length so a stray
    sentence in a notes column doesn't auto-coerce.

    Short-circuits when the cell value matches a known multi-species
    gene symbol : prevents non-English-locale month-name collisions
    (Polish `lis`, Dutch `mei`) from masquerading as date corruption.
    """
    s = s.strip()
    # Cross-species symbol guard: if the cell IS a known gene symbol in
    # any of MGI/ZFIN/FlyBase/WormBase, skip date interpretation. Only
    # applied to short tokens (gene-symbol-shape) : sentences are
    # length-filtered downstream anyway.
    if 2 <= len(s) <= 20 and s.lower() in _multispecies_symbols_lower():
        return []
    if not _DATE_STRING_RE.fullmatch(s):
        # Long-form fallback. Tightly bounded : dateutil.parser is famously
        # aggressive (will parse 'ST13P5' and '16:1SMOH' as dates).
        # Pre-filters:
        #   - length 8-32: rules out short numerics and long sentences
        #   - must contain a recognised Latin month abbreviation as an
        #     isolated word/token: any of `Jan`/`Feb`/.../`Dec` plus our
        #     locale tokens (`agosto`, `septembre`, `tammikuu`, etc.).
        #     A gene-list cell `'ST13;ST13P5;ST13P4'` has no month token;
        #     a lipid like `'16:1SMOH'` has no month token.
        #   - no pipes, semicolons, brackets, or '@/=' which signal
        #     accession/email/URL/expression payloads, not dates.
        contains_payload_char = any(ch in s for ch in "|@={}[]<>;:")
        # Look for a month token surrounded by non-letter boundaries.
        # Allows `Sept 2, 2024`, `2 septembre 2014`, `15 March 2014`.
        # Refuses substrings inside larger words (e.g. `ST13P5` doesn't
        # contain `mar`/`may` as an isolated token).
        _MONTH_NAME_RE = re.compile(
            r"(?<![A-Za-z])(?:"
            + _MONTH_ALTERNATION
            + r"|january|february|march|april|june|july|august|september"
            + r"|october|november|december"
            + r"|enero|febrero|abril|mayo|junio|julio|septiembre|octubre|noviembre|diciembre"
            + r"|janvier|fevrier|f[ée]vrier|avril|mai|juin|juillet|ao[uû]t|septembre|octobre|novembre|d[ée]cembre"
            + r"|gennaio|febbraio|marzo|aprile|maggio|giugno|luglio|agosto|settembre|ottobre|novembre|dicembre"
            + r"|m[äa]rz|m[äa]erz"
            + r")(?![A-Za-z])"
        )
        has_month_token = bool(_MONTH_NAME_RE.search(s.lower()))
        if (
            8 <= len(s) <= 32
            and has_month_token
            and not contains_payload_char
        ):
            try:
                from dateutil.parser import ParserError as _DUParserError
                from dateutil.parser import parse as _du_parse
                out: list[date] = []
                for dayfirst in (False, True):
                    try:
                        dt = _du_parse(s, fuzzy=False, dayfirst=dayfirst)
                    except (_DUParserError, ValueError, OverflowError, TypeError):
                        continue
                    d = dt.date()
                    if d not in out:
                        out.append(d)
                return out
            except ImportError:
                pass
        return []
    candidates: list[date] = []
    # ISO YYYY-MM-DD (optionally with time) : unambiguous
    m = re.match(r"^(\d{4})-(\d{1,2})-(\d{1,2})(?:[ T]\d{1,2}:\d{1,2}(?::\d{1,2})?)?$", s)
    if m:
        with contextlib.suppress(ValueError):
            candidates.append(date(int(m.group(1)), int(m.group(2)), int(m.group(3))))
        return candidates
    # Month-name forms : separators include `-`, `/`, `.` (German locale);
    # month token includes Latin + accented chars (German Mär, French déc).
    # `\w` with re.UNICODE catches `Mär`/`déc` etc; widened to {3,4} to
    # accommodate `mär` / `déc` after lower-casing.
    _MONTH_TOK = r"[\w]{3,4}"
    m = re.match(
        rf"^(\d{{1,2}})[-/.]({_MONTH_TOK})(?:[-/.](\d{{2,4}}))?$",
        s, re.IGNORECASE | re.UNICODE,
    )
    if m:
        day, mon_str, yr = int(m.group(1)), m.group(2).lower(), m.group(3)
        mon = _MONTH_ABBR.get(mon_str, 0)
        year = int(yr) if yr else date.today().year
        if year < 100:
            year += 2000
        if mon and 1 <= day <= 31:
            with contextlib.suppress(ValueError):
                candidates.append(date(year, mon, day))
        return candidates
    m = re.match(
        rf"^({_MONTH_TOK})[-/.](\d{{1,2}})(?:[-/.](\d{{2,4}}))?$",
        s, re.IGNORECASE | re.UNICODE,
    )
    if m:
        mon_str = m.group(1).lower()
        mon, day, yr = _MONTH_ABBR.get(mon_str, 0), int(m.group(2)), m.group(3)
        year = int(yr) if yr else date.today().year
        if year < 100:
            year += 2000
        if mon and 1 <= day <= 31:
            with contextlib.suppress(ValueError):
                candidates.append(date(year, mon, day))
        return candidates
    # All-numeric: try both DD/MM/YY and MM/DD/YY since locale is unknown
    m = re.match(r"^(\d{1,2})[/.\-](\d{1,2})[/.\-](\d{2,4})$", s)
    if m:
        a, b, yr = int(m.group(1)), int(m.group(2)), m.group(3)
        year = int(yr)
        if year < 100:
            year += 2000
        for day, mon in ((a, b), (b, a)): # try both interpretations
            if 1 <= mon <= 12 and 1 <= day <= 31:
                try:
                    cand = date(year, mon, day)
                    if cand not in candidates:
                        candidates.append(cand)
                except ValueError:
                    pass
    return candidates


import math


def _riken_coercion_signature(value: float) -> dict | None:
    """Test whether `value` looks like a coerced RIKEN cDNA ID.

    A RIKEN cDNA clone ID follows `\\d{7}[A-Z]\\d{2}` (FANTOM convention).
    When that letter is `E`, Excel parses it as scientific notation:
    `2310009E13` → `2.310009 * 10^13` = `2.310009e19`. The original
    string is unrecoverable.

    Signature (must satisfy all):
      1. Float (not int), positive.
      2. Exactly representable as a 7-digit integer mantissa × 10^k.
      3. Exponent k in RIKEN range [7, 31] (E01..E25 suffixes; suffix=01
         → exp=7, suffix=25 → exp=31).
      4. Mantissa has ≤ 1 trailing zero. (RIKEN clone IDs are
         distinctive : multi-trailing-zeros indicates a round number
         like 1.5e8, not a clone ID.)
      5. Mantissa has ≥ 3 distinct nonzero digit values. (Rejects
         repeat-digit pathologies like 1111110.)

    Returns the parsed `{mantissa, exp, trailing_zeros, distinct_nonzero}`
    dict on match, None otherwise. Used in Pass 1.5 to distinguish
    RIKEN coercion from legitimate large Entrez IDs or round float
    measurements.
    """
    if not isinstance(value, float) or value <= 0 or math.isnan(value):
        return None
    exp = math.floor(math.log10(value))
    # Shift so that a 7-digit mantissa becomes integer.
    shift = 10 ** (exp - 6)
    if shift <= 0:
        return None
    mantissa = value / shift
    if mantissa != int(mantissa):
        return None
    m_int = int(mantissa)
    if not (1_000_000 <= m_int <= 9_999_999):
        return None
    if not (7 <= exp <= 31):
        return None
    s = str(m_int)
    trailing_zeros = len(s) - len(s.rstrip("0"))
    distinct_nonzero = len({d for d in s if d != "0"})
    if trailing_zeros > 1 or distinct_nonzero < 3:
        return None
    return {
        "mantissa": m_int,
        "exp": exp,
        "trailing_zeros": trailing_zeros,
        "distinct_nonzero": distinct_nonzero,
    }


_HEADER_TOKEN_RE = re.compile(r"[a-z0-9]+")


def _header_tokens(name: str) -> set[str]:
    # Reject caption rows masquerading as headers: real column names are
    # typically <= 80 chars. A 200-char string is a table caption that
    # pandas pulled in as `df.columns[0]`. Empirically observed in
    # `PMC10313048/pone.0287634.s002.xlsx` where a full table caption
    # ("Supplemental Data Table S2. Total proteome list ...") was the
    # "header" and matched the identifier-hint set on the word "protein",
    # which produced 70 spurious id-float flags on legitimate numeric data.
    stripped = name.strip()
    if len(stripped) > 80:
        return set()
    return set(_HEADER_TOKEN_RE.findall(stripped.lower()))


def _tokens_match_hint_set(tokens: set[str], hints: set[str]) -> bool:
    """Match by exact token, or by suffix match for tokens >= len(hint)+2.

    The suffix match catches glued forms like 'ENTREZID' (one token containing
    'id') and 'genesymbol' (token containing 'symbol'). The +2 length floor
    avoids trivial matches on common short fragments.
    """
    if tokens & hints:
        return True
    for tok in tokens:
        if len(tok) < 4:
            continue
        for hint in hints:
            if len(hint) >= 2 and len(tok) >= len(hint) + 2 and tok.endswith(hint):
                return True
    return False


def _column_is_identifier(name: str, series: pd.Series) -> tuple[bool, str]:
    # value-pattern classification *first*. If the column's actual
    # values are external-DB IDs (Entrez/UniProt/RefSeq/Ensembl) or numeric
    # measurements, do NOT treat it as a gene-symbol identifier column even
    # if the header has identifier-shaped tokens. This eliminates the
    # EntrezGeneID / luminescence-measurement false-positive classes.
    from .column_classifier import classify_column
    classification = classify_column(series, name)
    if classification.column_type in (
        "entrez", "uniprot", "refseq", "ensembl", "hgnc_id",
        "measurement", "date", "empty",
    ):
        return False, (
            f"column classified as {classification.column_type} "
            f"({classification.reason}); not a gene-symbol corruption candidate"
        )

    # when the value classifier confidently identifies the column as
    # gene_symbol or riken, trust the values over the header. This prevents
    # false-negative misses on columns with weird headers like
    # "Supplementary Table S4. Prevalent genes in NPC stage." where the
    # suffix-matcher false-positives `stage` -> `age` and disqualifies a
    # genuine gene-symbol column (PMC4195499 case).
    if classification.column_type in ("gene_symbol", "riken"):
        return True, (
            f"value-pattern classifier identified column as "
            f"{classification.column_type} ({classification.reason})"
        )

    tokens = _header_tokens(name)
    if _tokens_match_hint_set(tokens, NON_IDENTIFIER_HEADER_HINTS):
        return False, f"header '{name}' contains a non-identifier token"
    header_match = _tokens_match_hint_set(tokens, IDENTIFIER_HEADER_HINTS)
    non_null = series.dropna()
    if len(non_null) == 0:
        return (header_match, "header match, empty column" if header_match else "empty")

    n_gene_like = 0
    n_known_hgnc = 0
    n_riken = 0
    n_non_string = 0
    for v in non_null:
        if isinstance(v, str):
            if v in HGNC_RENAME_MAP or v in HGNC_NEW_SYMBOLS:
                n_known_hgnc += 1
            if GENE_LIKE_PATTERN.fullmatch(v):
                n_gene_like += 1
            if RIKEN_PATTERN.fullmatch(v):
                n_riken += 1
        else:
            n_non_string += 1

    string_id_evidence = n_known_hgnc + n_gene_like + n_riken
    if string_id_evidence == 0 and not header_match:
        return False, "no string-shaped identifier evidence and no header hint"

    total = len(non_null)
    if n_known_hgnc >= 1 and (n_known_hgnc + n_non_string) / total >= 0.3:
        return True, f"{n_known_hgnc} known HGNC symbols; {n_non_string} non-string cells"
    if (n_gene_like + n_non_string) / total >= 0.5:
        return True, f"{n_gene_like} gene-like strings; {n_non_string} non-string cells"
    if (n_riken + n_non_string) / total >= 0.3:
        return True, f"{n_riken} RIKEN-like strings; {n_non_string} non-string cells"
    if header_match:
        return True, f"header match '{name}'; {string_id_evidence} id-shaped values"
    return False, "no identifier pattern dominant"


def _candidates_for_month_and_n(month: int, n: int) -> list[str]:
    """Reverse a month + integer-suffix into candidate gene symbols.

    `n` is either day-of-month or year-suffix (year mod 100); the calling site
    decides which interpretation applies for the cell.

    Coverage informed by Ziemann 2021 S1 (1,568 eukaryote vulnerable gene
    names across 5 kingdoms, PLOS Comp Bio 2021) and locale-specific
    corruptions documented in the same paper's "Novel error types" section.
    Validated to catch 100% of the 1,568 vulnerable-gene catalog after this
    expansion.
    """
    out: list[str] = []
    if month == 3 and 1 <= n <= 12:
        # MARCH1-12 (HGNC: MARCH12 → MARCHF12 exists per the 2020 rename)
        out.append(f"MARCH{n}")
        if n <= 2:
            out.append(f"MARC{n}")
    elif month == 9 and 1 <= n <= 15:
        out.append(f"SEPT{n}")
        out.append(f"SEP{n}")
    elif month == 10 and 1 <= n <= 11:
        out.append(f"OCT{n}")
    elif month == 11 and n == 1:
        out.append("NOV1")
    elif month == 12 and n in (1, 2):
        out.append(f"DEC{n}")
    elif month == 4 and 1 <= n <= 3:
        # APR1, APR2 (newly observed in Ziemann S1 fungal corruption),
        # and APR3 are all valid Argonaute-family aliases.
        out.append(f"APR{n}")
    elif month == 2 and n in (3, 4):
        out.append(f"FEB{n}")
    elif month == 8 and 1 <= n <= 12:
        # AGO (Argonaute) gene family : locale-specific corruption.
        # Spanish/Italian/Portuguese "Ago" abbreviates August; AGO2 typed
        # in those locales auto-converts to "Aug-02" not "Sep-02".
        # AGO1-AGO12 known across plants (Arabidopsis has 10+).
        out.append(f"AGO{n}")
    elif month == 5 and 1 <= n <= 31:
        # MEI gene family : Dutch locale "Mei" = May. MEI1 -> "May-01".
        # Also catches plant MAY1-24 gene names (Ziemann S1).
        out.append(f"MEI{n}")
        if n <= 24:
            out.append(f"MAY{n}")
        # JUN family (c-Jun, JUNB, JUND) : Ziemann 2021 "Novel error types":
        # the protein name `jun-1` is interpreted by Excel as June minus 1
        # = May 31. JUN1, JUN2, JUN3 are real human proto-oncogene members
        # of the AP-1 transcription factor complex.
        if n >= 29:
            jun_idx = 32 - n  # n=31 -> JUN1, n=30 -> JUN2, n=29 -> JUN3
            if 1 <= jun_idx <= 3:
                out.append(f"JUN{jun_idx}")
    elif month == 1 and 30 <= n <= 50:
        # TAMM gene family : Finnish "Tammikuu" = January. TAMM41 -> "Jan-41".
        # Range covers the small TAMM family (TAMM41 is the canonical case).
        out.append(f"TAMM{n}")
    return out


# --- HGNC canonical-rank machinery -----------------------------------------
# Deprecated → modern HGNC current-symbol mapping for the Excel-corruption
# families. Authoritative source: HGNC announcement (2020), confirmed via the
# bundled HGNC bulk download. Maintained as a hardcoded map (rather than
# loaded from HGNC every call) because (a) the rename is closed : no new
# entries are coming : and (b) this keeps detect() side-effect-free.
_HGNC_RENAME: dict[str, str] = {}
for _n in range(1, 16):
    _HGNC_RENAME[f"SEPT{_n}"] = f"SEPTIN{_n}"
for _n in range(1, 13): # MARCH1-12 (MARCHF12 exists per HGNC)
    _HGNC_RENAME[f"MARCH{_n}"] = f"MARCHF{_n}"
for _n in (1, 2):
    _HGNC_RENAME[f"MARC{_n}"] = f"MTARC{_n}"
_HGNC_RENAME["DEC1"] = "DELEC1"


def _canonicalize_candidates(candidates: list[str]) -> tuple[list[str], bool]:
    """Sort candidates so the modern HGNC current symbol leads, return
    (sorted_candidates_with_canonical_prepended, has_unique_canonical).

    `has_unique_canonical` is True iff exactly one candidate maps to a known
    modern HGNC symbol : that's the disambiguating signal worth a confidence
    bump. Never drops alternatives (per "no silent repair" principle).
    """
    if not candidates:
        return [], False
    canonicals: list[str] = []
    seen_canonical: set[str] = set()
    others: list[str] = []
    for c in candidates:
        modern = _HGNC_RENAME.get(c)
        if modern and modern not in seen_canonical:
            canonicals.append(modern)
            seen_canonical.add(modern)
        others.append(c)
    # Modern symbols first, then the historical candidates the detector
    # originally derived (which match what's literally in the corrupted cell).
    out: list[str] = []
    for c in canonicals + others:
        if c not in out:
            out.append(c)
    return out, len(canonicals) == 1


def _reverse_gene_date(d: date, number_format: str | None = None) -> list[str]:
    """Map a calendar date back to gene symbol(s) that Excel could have corrupted
    into it. Coverage derived empirically from HGNC's `prev_symbol` and
    `alias_symbol` columns plus the Ziemann 2016/2021 corpora.

    Two distinct corruption modes for date cells:

    - **day-of-month** (cell stored with format `dd-mmm` or `d-mmm`):
      "MARCH7" typed → Excel stores `date(<current-year>, 3, 7)` → reverse to MARCH7
    - **year-suffix** (cell stored with format `mmm-yy`):
      "SEPT7" typed → Excel stores `date(2007, 9, 1)` (day=1, year=2007) → reverse to SEPT7

    Decision rule:
    - `number_format` contains `y` but not `d` → year-suffix only
    - `number_format` contains `d`             → day-of-month only
    - `number_format` unknown                  → both interpretations emitted
      when day=1 and year falls in the typical 2000-2029 window (Excel's
      `mmm-yy` 2-digit-year disambiguation)
    """
    # Year-sanity gate: Excel coerces huge numbers (RIKEN IDs, CAS
    # numbers, evidence IDs) to dates with implausible years like
    # 3715 or 5817. Real gene-corruption sits in the publication-era
    # window (~1900-2100); outside that, the cell is numeric noise
    # that happens to decode to a date.
    if d.year < 1900 or d.year > 2100:
        return []

    fmt = number_format.lower() if number_format else None
    has_y = bool(fmt and "y" in fmt)
    has_d = bool(fmt and "d" in fmt)
    use_year = has_y and not has_d
    use_day = not use_year  # day-of-month is the default when format is ambiguous or unknown

    candidates: list[str] = []
    if use_day:
        candidates.extend(_candidates_for_month_and_n(d.month, d.day))
    # Window 2000-2059 covers SEPT1-15 (years 2001-2015), MARCH1-12, AGO1-12,
    # MEI1-31, MAY1-24, APR1-3, FEB3-4, OCT1-11, NOV1, DEC1-2, AND TAMM30-50
    # (Finnish locale corruption : years 2030-2050).
    if use_year or (fmt is None and d.day == 1 and 2000 <= d.year <= 2059):
        yy = d.year % 100
        for c in _candidates_for_month_and_n(d.month, yy):
            if c not in candidates:
                candidates.append(c)
    return candidates


# --- Numeric-identifier corruption pass --------------------------------------
# Three classes of corruption are detected here, all in identifier-shaped
# columns Pass 1 has already classified. Each is opt-in by evidence: we only
# flag a cell if the rest of the column carries an unambiguous "this should
# have been a [zero-padded ID / long integer / apostrophe-text]" signal.

_FLOAT64_EXACT_MAX = 2**53  # 9_007_199_254_740_992 : Excel int precision limit


def _column_leading_zero_widths(series: pd.Series) -> dict[int, int]:
    """For an identifier-shaped column, return `{width: count}` of cells that
    are leading-zero strings (e.g. `'00123'` has width 5)."""
    import re as _re
    widths: dict[int, int] = {}
    for v in series.dropna():
        if isinstance(v, str):
            s = v.strip()
            if _re.fullmatch(r"0\d+", s): # starts with 0, then digits
                widths[len(s)] = widths.get(len(s), 0) + 1
    return widths


_IDENTIFIER_LIKE_COLUMN_TYPES = frozenset({
    "gene_symbol", "entrez", "uniprot", "refseq",
    "ensembl", "hgnc_id", "riken",
})


def _column_is_any_identifier(col_name: str, series: pd.Series) -> bool:
    """True if the column classifies as ANY identifier shape (gene symbols,
    Entrez, UniProt, RefSeq, Ensembl, HGNC, RIKEN) OR the header tokens
    explicitly contain identifier hints. Wider than `_column_is_identifier`
    (which targets gene-symbol corruption specifically) : used by the
    numeric-ID corruption pass which applies to all ID columns."""
    from .column_classifier import classify_column
    cls = classify_column(series, col_name)
    if cls.column_type in _IDENTIFIER_LIKE_COLUMN_TYPES:
        return True
    tokens = _header_tokens(col_name)
    return bool(_tokens_match_hint_set(tokens, IDENTIFIER_HEADER_HINTS))


def _detect_numeric_id_corruptions(
    df: pd.DataFrame,
    report: Report,
    identifier_cols: set[str],
    flagged: set[tuple[str, int]],
    sheet: str | None,
) -> None:
    """Three corruption classes:

    1. **leading-zero-stripped**: column has cells like `'00123'`, `'00456'`
       at a consistent width; bare integers in the column whose string form
       is shorter than the column's typical width are likely stripped IDs.
       Suggestion: zero-pad to the typical width.

    2. **long-int-precision-loss**: integer in an identifier column that
       exceeds 2^53. Excel's float64 storage cannot represent these exactly;
       round-tripping silently truncates the trailing digits. No recovery
       possible : flag as warning.

    3. **apostrophe-stripped**: cell value looks like a bare numeric where
       the column otherwise enforces text storage (e.g., mostly leading-zero
       strings). Same as #1 in practice : the apostrophe was the text-storage
       mechanism Excel used; stripping it produces the bare numeric. We
       collapse these under `leading-zero-stripped` since the remediation
       (zero-pad and store as string) is identical.

    Runs on ANY identifier-shaped column (gene-symbol OR external-DB-ID),
    not just gene-symbol-shaped : long-integer IDs and leading-zero accession
    patterns are present in RefSeq/Ensembl/Entrez columns too.
    """
    for col in df.columns:
        col_s = str(col)
        series = df[col]
        widths = _column_leading_zero_widths(series)
        n_lz = sum(widths.values())
        # The leading-zero signal is itself column-classification evidence : 
        # 3+ values like '00123' make this an identifier column regardless of
        # the broader classifier verdict. For other identifier kinds we need
        # the classifier to confirm.
        header_has_id_hint = _tokens_match_hint_set(
            _header_tokens(col_s), IDENTIFIER_HEADER_HINTS,
        )
        is_id = (n_lz >= 3) or (col_s in identifier_cols) or \
                _column_is_any_identifier(col_s, series)
        # Edge case: every cell in the column got coerced to a float (the
        # "fully ruined RIKEN column" case). The classifier sees all-numeric
        # and says "measurement", but the header still says it's an ID
        # column. Trust the header.
        if not is_id and header_has_id_hint:
            non_null = series.dropna()
            # Lowered from 5 to 3 a small-pilot xlsx
            # with as few as 3 numeric IDs in a header-hinted column is still
            # informative. The header hint provides strong column-identity
            # signal that doesn't require many sample values to confirm.
            if len(non_null) >= 3 and all(
                isinstance(v, (int, float)) and not isinstance(v, bool)
                for v in non_null
            ):
                is_id = True
        if not is_id:
            continue
        # Surface header-promoted identifier columns in the report summary
        # too. Without this, a RIKEN/Ensembl column where every cell got
        # coerced to float passes the long-int check but the report header
        # still says "Identifier columns: 0", which is confusing.
        if col_s not in identifier_cols:
            identifier_cols.add(col_s)
            if col_s not in report.identifier_columns:
                report.identifier_columns.append(col_s)
        typical_width: int | None = None
        if n_lz >= 3:
            typical_width = max(widths, key=lambda w: widths[w])

        for idx, value in series.items():
            if pd.isna(value):
                continue
            key = (col_s, int(idx))
            if key in flagged:
                continue

            # --- 1. Leading-zero strip ---
            if typical_width is not None:
                stripped_form: str | None = None
                if isinstance(value, int) and not isinstance(value, bool):
                    if value >= 0 and len(str(value)) < typical_width:
                        stripped_form = str(value)
                elif isinstance(value, float) and value.is_integer():
                    iv = int(value)
                    if iv >= 0 and len(str(iv)) < typical_width:
                        stripped_form = str(iv)
                elif isinstance(value, str):
                    s = value.strip()
                    # Bare digit string with no leading zero AND shorter than typical
                    if s.isdigit() and not s.startswith("0") and len(s) < typical_width:
                        stripped_form = s
                if stripped_form is not None:
                    padded = stripped_form.zfill(typical_width)
                    report.suspicions.append(Suspicion(
                        column=col_s, row=int(idx), value=value,
                        kind="leading-zero-stripped",
                        suggestion=padded,
                        confidence=0.85 if len(widths) == 1 else 0.65,
                        reason=(
                            f"identifier column has {n_lz} leading-zero IDs at "
                            f"width {typical_width}; this cell ({value!r}, length "
                            f"{len(stripped_form)}) is likely a stripped form. "
                            f"Suggested zero-pad: {padded!r}"
                        ),
                        sheet=sheet,
                    ))
                    flagged.add(key)
                    continue

            # --- 2. Long-int precision loss ---
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                abs_v = abs(value)
                if abs_v >= _FLOAT64_EXACT_MAX:
                    report.suspicions.append(Suspicion(
                        column=col_s, row=int(idx), value=value,
                        kind="long-int-precision-loss",
                        suggestion=None,
                        confidence=0.85,
                        reason=(
                            f"integer {value!r} in identifier column exceeds "
                            f"2^53 ({_FLOAT64_EXACT_MAX:,}); Excel's float64 "
                            f"storage cannot represent it exactly : trailing "
                            f"digits may be silently truncated. Recover from "
                            f"original source; cannot reconstruct here."
                        ),
                        sheet=sheet,
                    ))
                    flagged.add(key)
                    continue
                # --- 3. id-float in fully-coerced column (RIKEN edge case) ---
                # Two distinguishable cases:
                #   (a) `value` is a true float >= 1e10 → almost certainly
                #       RIKEN-style E-coercion (RIKEN E02 and up coerces
                #       to >= 1e10). Catches the vast majority of RIKEN
                #       corruption with negligible false-positive risk : 
                #       no legitimate gene/probe ID is a float magnitude
                #       this large.
                #   (b) `value` is a *float* with magnitude in [1e7, 1e10]
                #       AND the header explicitly names RIKEN/cDNA : the
                #       only place RIKEN E01 corruption produces a value
                #       this small. We do NOT flag this range without the
                #       header anchor because real Entrez IDs reach 1.5e8
                #       (as of 2026); flagging would false-positive on
                #       every Entrez ID column.
                _RIKEN_HEADER_HINT = {"riken", "cdna", "transcript"}
                col_hint_tokens = _header_tokens(col_s)
                is_riken_hint_col = bool(col_hint_tokens & _RIKEN_HEADER_HINT)
                if (header_has_id_hint
                        and isinstance(value, float)
                        and not isinstance(value, bool)):
                    sig = _riken_coercion_signature(value)
                    # Three decision branches, in increasing scepticism:
                    if abs_v >= 1e10:
                        # Magnitude alone is the gate at 1e10+. No legitimate
                        # bare-integer ID reaches this; flag with
                        # signature-boosted confidence.
                        if sig is not None:
                            confidence = 0.95
                            reason = (
                                f"float magnitude {abs_v:.3g} in identifier "
                                f"column matches the RIKEN-coercion signature: "
                                f"7-digit mantissa {sig['mantissa']}, exp "
                                f"{sig['exp']}, {sig['trailing_zeros']} "
                                f"trailing zeros, {sig['distinct_nonzero']} "
                                f"distinct nonzero digits : high-confidence "
                                f"RIKEN E-suffix corruption."
                            )
                        else:
                            confidence = 0.80
                            reason = (
                                f"float magnitude {abs_v:.3g} in identifier "
                                f"column; bare numerics are not legitimate "
                                f"gene/probe IDs at this magnitude : likely "
                                f"some form of precision-losing coercion."
                            )
                        report.suspicions.append(Suspicion(
                            column=col_s, row=int(idx), value=value,
                            kind="id-float", suggestion=None,
                            confidence=confidence,
                            reason=reason, sheet=sheet,
                        ))
                        flagged.add(key)
                    elif sig is not None and abs_v >= 1e7:
                        # 7-sig-fig signature with mantissa entropy is a
                        # strong RIKEN tell even in non-RIKEN-hinted columns.
                        # Boosted higher when the header anchors RIKEN.
                        confidence = 0.85 if is_riken_hint_col else 0.65
                        report.suspicions.append(Suspicion(
                            column=col_s, row=int(idx), value=value,
                            kind="id-float", suggestion=None,
                            confidence=confidence,
                            reason=(
                                f"float {abs_v:.6g} in identifier column has "
                                f"RIKEN-coercion signature: 7-digit mantissa "
                                f"{sig['mantissa']}, exp {sig['exp']}, "
                                f"{sig['trailing_zeros']} trailing zeros, "
                                f"{sig['distinct_nonzero']} distinct nonzero "
                                f"digits. " +
                                ("Header anchors RIKEN/cDNA : high "
                                 "confidence." if is_riken_hint_col
                                 else "No RIKEN header anchor : moderate "
                                 "confidence. Re-import source as text.")
                            ),
                            sheet=sheet,
                        ))
                        flagged.add(key)
                    elif is_riken_hint_col and abs_v >= 1e7:
                        # Header-anchored fallback for any large float
                        # without the strict signature.
                        report.suspicions.append(Suspicion(
                            column=col_s, row=int(idx), value=value,
                            kind="id-float", suggestion=None, confidence=0.55,
                            reason=(
                                f"RIKEN-hinted column ({col_s!r}) has float "
                                f"value {abs_v:.3g}; consistent with RIKEN E01 "
                                f"coercion (header-anchored; signature does "
                                f"not match : could also be a precision-"
                                f"losing measurement). Re-import as text."
                            ),
                            sheet=sheet,
                        ))
                        flagged.add(key)


_SAMPLING_ROW_THRESHOLD = 100_000
_SAMPLING_HEAD_TAIL = 5_000
_SAMPLING_INTERIOR_FRACTION = 0.10
_SAMPLING_SEED = 20260518  # deterministic for reproducibility (Pillar 5)


def _maybe_subsample(df: pd.DataFrame) -> tuple[pd.DataFrame, str | None]:
    """For large dataframes, apply a deterministic stratified sample:
    first N rows + last N rows + 10% random interior. Returns the sample
    plus a human-readable note (or None if no sampling occurred).

    Threshold = 100,000 rows, per the literature review recommendation
    (Broman & Woo 2018 reporting norms; openpyxl performance ceiling at
    ~1M cells / 5M cells OOM on commodity hardware).
    """
    total = len(df)
    if total <= _SAMPLING_ROW_THRESHOLD:
        return df, None
    head = df.iloc[:_SAMPLING_HEAD_TAIL]
    tail = df.iloc[-_SAMPLING_HEAD_TAIL:]
    interior_start = _SAMPLING_HEAD_TAIL
    interior_end = total - _SAMPLING_HEAD_TAIL
    interior = df.iloc[interior_start:interior_end]
    n_interior = int(len(interior) * _SAMPLING_INTERIOR_FRACTION)
    interior_sample = interior.sample(
        n=min(n_interior, len(interior)),
        random_state=_SAMPLING_SEED,
    )
    sampled = pd.concat([head, interior_sample, tail]).sort_index()
    note = (
        f"sampled {len(sampled)} of {total} rows: head {_SAMPLING_HEAD_TAIL} + "
        f"tail {_SAMPLING_HEAD_TAIL} + {_SAMPLING_INTERIOR_FRACTION*100:.0f}% "
        f"interior random (seed={_SAMPLING_SEED}). "
        f"Per-row corruption rates extrapolate from the sample; cell-level "
        f"counts are NOT a full census."
    )
    return sampled, note


_HGNC_DRIFT_WARNING_THRESHOLD_DAYS = 30


def _maybe_hgnc_drift_note() -> str | None:
    """If the HGNC snapshot is older than 30 days, return a one-line
    warning suitable for `Report.sampling_note`."""
    try:
        from .corpus import hgnc_snapshot_age_days
        age = hgnc_snapshot_age_days()
    except Exception:
        return None
    if age is None or age <= _HGNC_DRIFT_WARNING_THRESHOLD_DAYS:
        return None
    return (
        f"HGNC snapshot is {age} days old (threshold "
        f"{_HGNC_DRIFT_WARNING_THRESHOLD_DAYS}); symbols added or renamed "
        f"after the snapshot are silently undetectable. Run "
        f"`scripts/refresh_hgnc.py` to refresh."
    )


def detect(df: pd.DataFrame, sheet: str | None = None) -> Report:
    df, sampling_note = _maybe_subsample(df)
    report = Report(rows_scanned=len(df), columns_scanned=len(df.columns))
    if sampling_note is not None:
        report.sampling_note = sampling_note
    drift_note = _maybe_hgnc_drift_note()
    if drift_note is not None:
        report.sampling_note = (
            f"{report.sampling_note} | {drift_note}"
            if report.sampling_note else drift_note
        )
    identifier_cols: set[str] = set()
    flagged: set[tuple[str, int]] = set()  # (column, row) already in suspicions

    # Per-cell openpyxl number_format strings, when the loader supplied them.
    # Key shape: (column_name, row_index). Missing → format unknown.
    cell_formats: dict[tuple[str, int], str] = df.attrs.get("cell_formats", {})

    # Legitimate-date-column pre-pass (shared across all passes).
    # Identifies columns that hold real calendar dates so gene-date-string
    # doesn't mass-flag them. Real gene-corruption is concentrated in
    # months {3, 9, 12} (March/Sept/Dec families); a real calendar
    # column spans many months. Header keywords reinforce the signal.
    _DATE_HEADER_TOKENS = (
        "date", " day", "day ", "month", "year", "time", "timestamp",
        "datetime", "period", "calendar", "visit", "session",
        "observation_date", "collection_date", "publication_date",
    )
    legitimate_date_columns: set[str] = set()
    for _col in df.columns:
        _col_lower = str(_col).lower()
        _header_says_date = any(t in _col_lower for t in _DATE_HEADER_TOKENS)
        _series = df[_col]
        _n_dates = 0
        _n_strings = 0
        _months: set[int] = set()
        for _v in _series:
            if pd.isna(_v): continue
            if isinstance(_v, str):
                _n_strings += 1
                _ps = _parse_date_string(_v)
                if _ps:
                    _n_dates += 1
                    for _d in _ps:
                        _months.add(_d.month)
        if _n_strings >= 5:
            _frac = _n_dates / _n_strings
            if _header_says_date and _frac >= 0.5 or len(_months) >= 5 and _frac >= 0.8:
                legitimate_date_columns.add(str(_col))

    # Quantitative-measurement-column pre-pass.
    # Suppress gene-date-serial Pass-2/3 emission on columns whose
    # header indicates a quantitative measurement (RPKM, TPM, counts,
    # log2-fold-change, etc.) or publication metadata (page numbers,
    # donor IDs, CAS numbers). These hold integers that can coincidentally
    # decode as Excel date-serials but are categorically not gene symbols.
    _MEASUREMENT_HEADER_TOKENS = (
        "count", "value", "expression", "sample", "rpkm", "fpkm",
        "tpm", "reads", "intensity", "abundance", "score", "p_value",
        "pvalue", "p-value", "qvalue", "q_value", "q-value", "fold",
        "ratio", "log2", "log10", "ijc", "sjc", " fc ", "_fc", "fc_",
        "fdr", "padj", "p.adj", "p adj", "stat", "_zscore", "zscore",
        "mean", "median", "stdev", "variance", "min", "max", "_n_",
        "depth", "coverage",
        # Publication-metadata + cohort-ID columns that hold integer
        # codes (subject IDs, page numbers, CAS-shaped IDs) which
        # can coincidentally decode as date-serials.
        "donor", "subject", "patient", "participant",
        "pages", "page_", "_page",
        "volume", "issue", "doi", "pubmed",
        "cas#", "cas no", "cas-no", "casno",
    )
    _IDENTIFIER_HEADER_TOKENS = (
        "gene", "symbol", "id", "identifier", "ensembl", "entrez",
        "uniprot", "hgnc", "refseq", "mirna", "protein", "transcript",
        "feature", "name",
    )
    quantitative_measurement_columns: set[str] = set()
    for _col in df.columns:
        _col_lower = str(_col).lower()
        _header_says_measurement = any(
            t in _col_lower for t in _MEASUREMENT_HEADER_TOKENS
        )
        _header_says_identifier = any(
            t in _col_lower for t in _IDENTIFIER_HEADER_TOKENS
        )
        if _header_says_measurement and not _header_says_identifier:
            quantitative_measurement_columns.add(str(_col))
            continue
        # Fallback: column is >=95% numeric (int or float) AND header
        # doesn't claim identifier : likely a measurement column even
        # without a known header keyword.
        _series = df[_col]
        _n_nonnull = 0
        _n_numeric = 0
        for _v in _series:
            if pd.isna(_v): continue
            _n_nonnull += 1
            if isinstance(_v, (int, float)) and not isinstance(_v, bool):
                _n_numeric += 1
        if (_n_nonnull >= 10 and not _header_says_identifier
                and _n_numeric / _n_nonnull >= 0.95):
            quantitative_measurement_columns.add(str(_col))

    # Pass 1: column-context-aware detection (high confidence)
    for col in df.columns:
        series = df[col]
        is_id, reason = _column_is_identifier(str(col), series)
        if not is_id:
            continue
        identifier_cols.add(str(col))
        report.identifier_columns.append(str(col))
        # pre-scan for date-corruption evidence in this column.
        # A coincidental integer in date-serial range maps to a gene-family
        # candidate ~10% of the time : without auxiliary evidence the
        # gene-date-serial flag has 9% empirical precision (xref-validated
        # against 500 Koh files). Require the column to *already* show
        # date corruption (any genuine date/datetime cell that maps to a
        # gene candidate) before flagging integer-decoded gene-date-serial.
        #
        # also count the total. Used downstream by the
        # ambiguous-canonical 0.60-band emission: a LONE date in an
        # identifier column is more likely a publication date than
        # autofill propagation; require >=2 date-corruption cells to
        # keep the 0.60 confidence, otherwise demote to 0.30.
        column_date_corruption_count = 0
        for value in series:
            if pd.isna(value):
                continue
            d_check: date | None = None
            if isinstance(value, datetime):
                d_check = value.date()
            elif isinstance(value, date):
                d_check = value
            if d_check is not None and _reverse_gene_date(d_check):
                column_date_corruption_count += 1
                if column_date_corruption_count >= 2:
                    break  # don't need exact count beyond threshold
        column_has_date_corruption = column_date_corruption_count >= 1

        for idx, value in series.items():
            if pd.isna(value):
                continue
            # Time-coercion: cells like `1:3`, `11:30`, `12:30` typed into
            # Excel get auto-converted to time-of-day. Common in plate/well
            # IDs (e.g., row:column coordinates) and in any text cell that
            # contains a colon between two single-digit numbers. Documented
            # by Pyle et al. (Escape Excel, PMC5617173) as a real corruption
            # class. openpyxl reads these as datetime.time objects.
            if isinstance(value, _time):
                hh, mm = value.hour, value.minute
                report.suspicions.append(Suspicion(
                    column=str(col), row=int(idx), value=value,
                    kind="time-coercion",
                    suggestion=f"{hh}:{mm}",
                    confidence=0.80,
                    reason=(
                        f"datetime.time value in identifier column ({reason}); "
                        f"plate/well coordinate or similar ID text was likely "
                        f"auto-converted to a time of day. Original text was "
                        f"approximately '{hh}:{mm}'."
                    ),
                    sheet=sheet,
                ))
                flagged.add((str(col), int(idx)))
                continue
            # String-form time: the same corruption serialized through CSV
            # comes back as the *string* "01:03" or "11:30:00". Same root
            # cause as datetime.time above, different on-disk form. Caught
            # only in identifier columns so a free-text "23:59 timestamp"
            # in a notes column doesn't trigger.
            if isinstance(value, str):
                # extend cross-species symbol guard to ALL string
                # detection paths. If the cell IS a known cross-species
                # symbol, skip every per-cell detection : the value is a
                # legitimate gene, not corruption. Excludes English-date
                # shapes which ARE the corruption signature.
                _vstrip = value.strip()
                if (2 <= len(_vstrip) <= 20
                        and _vstrip.lower() in _multispecies_symbols_lower()):
                    continue
                _tm = _TIME_STRING_RE.match(_vstrip)
                if _tm is not None:
                    hh = int(_tm.group(1))
                    mm = int(_tm.group(2))
                    report.suspicions.append(Suspicion(
                        column=str(col), row=int(idx), value=value,
                        kind="time-coercion",
                        suggestion=f"{hh}:{mm:02d}",
                        confidence=0.75,
                        reason=(
                            f"HH:MM string in identifier column ({reason}); "
                            f"plate/well coordinate or similar ID text was "
                            f"likely auto-converted to time of day, then "
                            f"round-tripped through CSV. Source text was "
                            f"approximately '{hh}:{mm:02d}'."
                        ),
                        sheet=sheet,
                    ))
                    flagged.add((str(col), int(idx)))
                    continue
                # Homoglyph corruption: in an identifier column we expect
                # ASCII Latin (HGNC convention). A cell with non-ASCII
                # characters that map to ASCII via Unicode TR39 confusables
                # is flagged with the proposed repair so the scientist
                # decides : never silently overwritten (Ziemann rule).
                if not value.isascii():
                    matches, repaired, repairs = _homoglyph_repair(
                        value, GENE_LIKE_PATTERN
                    )
                    if matches:
                        repair_detail = "; ".join(
                            f"{r.original!r} ({r.codepoint} {r.name}) → {r.replacement!r}"
                            for r in repairs
                        )
                        report.suspicions.append(Suspicion(
                            column=str(col), row=int(idx), value=value,
                            kind="homoglyph",
                            suggestion=repaired,
                            confidence=0.85,
                            reason=(
                                f"non-ASCII character(s) in identifier "
                                f"column ({reason}); maps to a valid "
                                f"gene-shape symbol {repaired!r} after "
                                f"Unicode TR39 confusable normalization. "
                                f"Substitutions: {repair_detail}."
                            ),
                            sheet=sheet,
                        ))
                        flagged.add((str(col), int(idx)))
                        continue
            if isinstance(value, datetime):
                value = value.date()
            if isinstance(value, date):
                fmt = cell_formats.get((str(col), int(idx)))
                candidates = _reverse_gene_date(value, fmt)
                candidates, unique_canonical = _canonicalize_candidates(candidates)
                suggestion = " | ".join(candidates) if candidates else None
                if len(candidates) == 1:
                    confidence = 0.95
                elif candidates:
                    # 0.5 base + 0.1 boost when exactly one candidate has a
                    # unique modern HGNC canonical (e.g., SEPT2 -> SEPTIN2),
                    # since the canonical leads the suggestion list.
                    #
                    # Gate the canonical boost on column-internal
                    # corroboration. A LONE date cell in an identifier
                    # column is empirically more often a publication/sample
                    # date than autofill propagation. Demote isolated
                    # date-in-id-column to 0.30; only emit 0.60/0.50 when
                    # at least one OTHER date cell in the column shares
                    # the corruption signature.
                    if column_date_corruption_count >= 2:
                        confidence = 0.6 if unique_canonical else 0.5
                    else:
                        confidence = 0.30 if unique_canonical else 0.20
                else:
                    confidence = 0.0
                reason_text = f"date in identifier column ({reason})"
                if fmt:
                    reason_text += f"; cell format {fmt!r}"
                if unique_canonical and len(candidates) > 1:
                    reason_text += f"; canonical HGNC symbol {candidates[0]!r} leads suggestion"
                if column_date_corruption_count < 2 and candidates and len(candidates) > 1:
                    reason_text += (
                        f"; isolated date in column (only {column_date_corruption_count} "
                        f"date-corruption cell observed) : confidence demoted"
                    )
                report.suspicions.append(Suspicion(
                    column=str(col), row=int(idx), value=value, kind="gene-date",
                    suggestion=suggestion, confidence=confidence,
                    reason=reason_text, sheet=sheet,
                ))
                flagged.add((str(col), int(idx)))
            elif not isinstance(value, str):
                # v6: bare id-float flags (any non-string in identifier column)
                # caused ~50K false positives per misclassified file. Two
                # legitimate sub-cases remain:
                #   (a) numeric value decodes to a uncorrupt date signature
                #       → flag as gene-date-serial with the suggested gene
                #   (b) large float (>= 1e10) in identifier column → almost
                #       certainly RIKEN/exponent-string coercion (no decode
                #       possible because precision is already lost)
                flagged_here = False
                decoded = None
                if isinstance(value, int) and 1 <= value <= 1e9:
                    decoded = _serial_to_date(value)
                elif isinstance(value, float) and value.is_integer() and 1 <= value <= 1e9:
                    decoded = _serial_to_date(int(value))
                if decoded is not None:
                    candidates = _reverse_gene_date(decoded)
                    if candidates and column_has_date_corruption:
                        # Empirically calibrated against xref evidence on
                        # 500 Koh files: gene-date-serial gates on
                        # auxiliary column evidence brings precision from
                        # 9% (raw) → ~95% (gated). Without the gate the
                        # 0.85 confidence label was systematically wrong.
                        report.suspicions.append(Suspicion(
                            column=str(col), row=int(idx), value=value,
                            kind="gene-date-serial",
                            suggestion=" | ".join(candidates),
                            confidence=0.85,
                            reason=(
                                f"numeric value decodes to uncorrupt date "
                                f"{decoded.isoformat()} in identifier column "
                                f"({reason}); column has companion date-"
                                f"corruption evidence : high confidence."
                            ),
                            sheet=sheet,
                        ))
                        flagged.add((str(col), int(idx)))
                        flagged_here = True
                    elif candidates and not column_has_date_corruption:
                        # No auxiliary evidence : low-confidence informational
                        # flag. Most of these are coincidental integers
                        # whose value happens to decode to a uncorrupt date,
                        # not real corruption.
                        report.suspicions.append(Suspicion(
                            column=str(col), row=int(idx), value=value,
                            kind="gene-date-serial",
                            suggestion=" | ".join(candidates),
                            confidence=0.15,
                            reason=(
                                f"numeric value decodes to uncorrupt date "
                                f"{decoded.isoformat()} in identifier column "
                                f"({reason}); NO companion date-corruption "
                                f"evidence in column : low confidence "
                                f"(empirically ~9% of such flags are real)."
                            ),
                            sheet=sheet,
                        ))
                        flagged.add((str(col), int(idx)))
                        flagged_here = True
                if (not flagged_here
                    and isinstance(value, (int, float))
                    and not isinstance(value, bool)
                    and abs(value) >= 1e10):
                    # int OR float : pandas may read back a coerced float as
                    # int after xlsx round-trip due to precision/repr quirks.
                    report.suspicions.append(Suspicion(
                        column=str(col), row=int(idx), value=value,
                        kind="id-float", suggestion=None, confidence=0.85,
                        reason=f"large numeric value in identifier column ({reason}); "
                               f"likely RIKEN/exponent-string coercion, precision lost",
                        sheet=sheet,
                    ))
                    flagged.add((str(col), int(idx)))
            elif isinstance(value, str) and value.strip().isdigit():
                n = int(value.strip())
                decoded = _serial_to_date(n)
                if decoded is not None:
                    candidates = _reverse_gene_date(decoded)
                    if candidates and column_has_date_corruption:
                        report.suspicions.append(Suspicion(
                            column=str(col), row=int(idx), value=value,
                            kind="gene-date-serial",
                            suggestion=" | ".join(candidates),
                            confidence=0.75,
                            reason=(
                                f"digit-string {value!r} in identifier "
                                f"column ({reason}); decodes as date "
                                f"{decoded.isoformat()}; column has "
                                f"companion date-corruption evidence."
                            ),
                            sheet=sheet,
                        ))
                        flagged.add((str(col), int(idx)))
                    elif candidates:
                        # Same low-confidence treatment as integer path
                        report.suspicions.append(Suspicion(
                            column=str(col), row=int(idx), value=value,
                            kind="gene-date-serial",
                            suggestion=" | ".join(candidates),
                            confidence=0.15,
                            reason=(
                                f"digit-string {value!r} decodes as date "
                                f"{decoded.isoformat()} in identifier column "
                                f"({reason}); NO companion date-corruption "
                                f"evidence."
                            ),
                            sheet=sheet,
                        ))
                        flagged.add((str(col), int(idx)))

    # Pass 1.5: numeric-identifier corruption sweep (leading-zero strip +
    # long-int precision loss). Runs over the same identifier columns Pass 1
    # already classified, so this never widens the false-positive surface.
    _detect_numeric_id_corruptions(
        df, report, identifier_cols, flagged, sheet=sheet,
    )

    # Pass 2: column-context-blind scan for date values matching the gene-
    # corruption signature. Catches files where the column header is missing
    # or pandas misinterpreted the layout (multi-row preamble, no header row,
    # etc.). Lower confidence because we lack column-shape corroboration.
    from .column_classifier import classify_column as _classify_column_v6
    for col in df.columns:
        col_s = str(col)
        if _tokens_match_hint_set(_header_tokens(col_s), NON_IDENTIFIER_HEADER_HINTS):
            continue  # column is explicitly date/time/etc : skip
        # also skip columns the classifier identifies as external-DB
        # IDs or measurement data. Pass 2 was producing massive false-positive
        # counts in these columns because random integers in date-serial range
        # match uncorrupt signatures by chance.
        cls = _classify_column_v6(df[col], col_s)
        if cls.column_type in (
            "entrez", "uniprot", "refseq", "ensembl", "hgnc_id",
            "measurement", "date", "empty",
        ):
            continue
        for idx, value in df[col].items():
            if (col_s, int(idx)) in flagged or pd.isna(value):
                continue
            # date object
            if isinstance(value, datetime):
                value = value.date()
            if isinstance(value, date):
                fmt = cell_formats.get((col_s, int(idx)))
                candidates = _reverse_gene_date(value, fmt)
                if not candidates:
                    continue
                candidates, unique_canonical = _canonicalize_candidates(candidates)
                suggestion = " | ".join(candidates)
                confidence = 0.45 if len(candidates) == 1 else 0.3 if unique_canonical else 0.25
                reason_text = "date matches uncorrupt signature (no column context)"
                if fmt:
                    reason_text += f"; cell format {fmt!r}"
                if unique_canonical and len(candidates) > 1:
                    reason_text += f"; canonical HGNC symbol {candidates[0]!r} leads"
                report.suspicions.append(Suspicion(
                    column=col_s, row=int(idx), value=value, kind="gene-date",
                    suggestion=suggestion, confidence=confidence,
                    reason=reason_text, sheet=sheet,
                ))
                continue
            # date-formatted string (e.g. "06/03/14", "Mar-6", "Sep-02")
            if isinstance(value, str):
                # cross-species symbol guard. Same as Pass 1.
                _vstrip = value.strip()
                if (2 <= len(_vstrip) <= 20
                        and _vstrip.lower() in _multispecies_symbols_lower()):
                    continue
                # OCT-3/4 alias form: cells like `"Oct-3/4"` or `"3/4-Oct"`
                # are real published gene-symbol aliases for POU5F1 (the
                # famous stem-cell marker). The trailing `/4` or leading
                # `3/4` defeats `_parse_date_string`'s standard match, so
                # detect it directly with a dedicated regex before the
                # generic date-string path. Always emit POU5F1 (modern HGNC)
                # plus the alias forms.
                _OCT_34_RE = re.compile(
                    r"^(?:oct[-./.]3/4|3/4[-./.]oct|oct[-./.]?3/?4|oct3/4)$",
                    re.IGNORECASE,
                )
                if _OCT_34_RE.match(value.strip()):
                    report.suspicions.append(Suspicion(
                        column=col_s, row=int(idx), value=value,
                        kind="gene-date-string",
                        suggestion="POU5F1 | OCT4 | OCT3",
                        confidence=0.80,
                        reason=(
                            "value matches OCT-3/4 alias form for POU5F1 "
                            "(stem-cell pluripotency marker); the slash "
                            "defeats Excel's default date parse but the "
                            "string itself indicates the original gene."
                        ),
                        sheet=sheet,
                    ))
                    continue
                # CAS Registry number: \d{2,7}-\d{1,2}-\d (chemistry/
                # metabolomics IDs). When the first component is small
                # (1-12), Excel may parse it as a date : e.g., `5-10-3` is
                # read as May 10, 2003. Detect the format only when:
                #   (a) the column header doesn't already advertise CAS/chem
                #       context (e.g. "CAS", "compound", "molecule",
                #       "metabolite", "chemical", "drug"), and
                #   (b) the column ISN'T already a chemistry-token column : 
                #       i.e. avoid false positives when the column is
                #       legitimately CAS-typed and the values are intact.
                _CAS_RE = re.compile(r"^\d{2,7}-\d{1,2}-\d$")
                if _CAS_RE.match(value.strip()):
                    col_lower = col_s.lower()
                    _CHEM_HEADER_TOKENS = (
                        "cas", "compound", "molecule", "metabolite",
                        "chemical", "drug", "ligand", "smiles", "inchi",
                    )
                    is_chem_header = any(t in col_lower for t in _CHEM_HEADER_TOKENS)
                    # In a chemistry column the format is *expected*, not
                    # corruption : emit info-level (confidence 0.20) for
                    # the chain-of-custody record without bothering the
                    # scientist. In a non-chemistry identifier column it's
                    # noteworthy (confidence 0.70).
                    cas_confidence = 0.20 if is_chem_header else 0.70
                    cas_reason = (
                        f"CAS Registry-format identifier "
                        f"({value.strip()!r}); chemistry-context column "
                        f"header ({col_s!r}) : preserved as text, no "
                        f"corruption suspected."
                        if is_chem_header
                        else
                        f"CAS Registry-format identifier "
                        f"({value.strip()!r}); Excel may auto-convert "
                        f"date-shaped CAS numbers : preserve as text."
                    )
                    report.suspicions.append(Suspicion(
                        column=col_s, row=int(idx), value=value,
                        kind="cas-registry",
                        suggestion=value.strip(),
                        confidence=cas_confidence,
                        reason=cas_reason,
                        sheet=sheet,
                    ))
                    continue
                parsed_dates = _parse_date_string(value)
                if not parsed_dates:
                    # Some cells store the corruption alongside the author's
                    # manual annotation, e.g. `'2-Oct,ATOCT2,OCT2'` (date
                    # first) or `'ATOCT2,2-Oct,OCT2'` (date in middle).
                    # Try EVERY token on any plausible intra-cell separator : 
                    # comma, semicolon, tab, newline, pipe. The date can
                    # appear in any position.
                    for delim in (",", ";", "\t", "\n", "|"):
                        if delim in value:
                            for tok in value.split(delim):
                                parsed_dates = _parse_date_string(tok.strip())
                                if parsed_dates:
                                    break
                            if parsed_dates:
                                break
                if parsed_dates:
                    all_candidates: list[str] = []
                    for d in parsed_dates:
                        all_candidates.extend(_reverse_gene_date(d))
                    seen: set[str] = set()
                    unique_candidates: list[str] = []
                    for c in all_candidates:
                        if c not in seen:
                            seen.add(c)
                            unique_candidates.append(c)
                    if unique_candidates and str(col) not in legitimate_date_columns:
                        if len(parsed_dates) == 1 and len(unique_candidates) == 1:
                            confidence = 0.55
                        elif len(unique_candidates) == 1:
                            confidence = 0.40
                        else:
                            confidence = 0.25
                        report.suspicions.append(Suspicion(
                            column=col_s, row=int(idx), value=value, kind="gene-date-string",
                            suggestion=" | ".join(unique_candidates),
                            confidence=confidence,
                            reason="date-formatted string parses as uncorrupt candidate(s)",
                            sheet=sheet,
                        ))
                        continue
                # Digit-string that decodes as date-serial → uncorrupt?
                stripped = value.strip()
                if (stripped.isdigit()
                    and str(col) not in quantitative_measurement_columns):
                    decoded = _serial_to_date(int(stripped))
                    if decoded is not None:
                        candidates = _reverse_gene_date(decoded)
                        if candidates:
                            confidence = 0.30 if len(candidates) >= 1 else 0.20
                            report.suspicions.append(Suspicion(
                                column=col_s, row=int(idx), value=value,
                                kind="gene-date-serial",
                                suggestion=" | ".join(candidates),
                                confidence=confidence,
                                reason=f"digit-string decodes as uncorrupt date "
                                       f"{decoded.isoformat()} (no column context)",
                                sheet=sheet,
                            ))
                continue
            # Integer (or integer-valued float) that decodes as date-serial.
            # suppress entirely for quantitative-measurement columns.
            if str(col) in quantitative_measurement_columns:
                continue
            n: int | None = None
            if isinstance(value, int) and not isinstance(value, bool):
                n = value
            elif isinstance(value, float) and value.is_integer() and 1 <= value <= 1e9:
                n = int(value)
            if n is not None:
                decoded = _serial_to_date(n)
                if decoded is not None:
                    candidates = _reverse_gene_date(decoded)
                    if candidates:
                        confidence = 0.30 if len(candidates) == 1 else 0.20
                        report.suspicions.append(Suspicion(
                            column=col_s, row=int(idx), value=value,
                            kind="gene-date-serial",
                            suggestion=" | ".join(candidates),
                            confidence=confidence,
                            reason=f"numeric value decodes as uncorrupt date "
                                   f"{decoded.isoformat()} (no column context)",
                            sheet=sheet,
                        ))

    # Pass 3: autofill-sequence detection. Documented in Ziemann 2021
    # (PMC6330011 example): when one corrupted date cell is autofilled by
    # the user, Excel propagates `Feb-97, Aug-97, Nov-97, Feb-98, ...` : 
    # a series of dates whose underlying gene-symbol mapping forms a
    # natural family sequence (SEPT2 → SEPT3 → ... or AUG-N → AUG-(N+1)).
    # The Pass-1 path catches each cell individually IF the column is
    # identifier-classified, but the autofill pattern itself is a stronger
    # signal : a column with a perfect monotonic date increment where each
    # step maps to a gene-family member is almost certainly autofill
    # corruption, not real data.
    _detect_autofill_sequences(df, report, flagged, sheet=sheet)
    _detect_decimal_comma_columns(df, report, sheet=sheet)
    _detect_out_of_family_symbols(df, report, identifier_cols,
                                  flagged, sheet=sheet)

    return report


_DECIMAL_COMMA_RE = re.compile(r"^-?\d{1,9},\d{1,9}$")
_DECIMAL_PERIOD_RE = re.compile(r"^-?\d{1,9}\.\d{1,9}$")


def _detect_out_of_family_symbols(
    df: pd.DataFrame,
    report: Report,
    identifier_cols: set[str],
    flagged: set[tuple[str, int]],
    sheet: str | None = None,
) -> None:
    """Detect gene-shape symbols that aren't in our HGNC snapshot.

    A symbol matching `GENE_LIKE_PATTERN` that's neither a current HGNC
    symbol nor a known prev/alias is one of three things:
      (a) A new gene HGNC hasn't catalogued yet (drift)
      (b) A non-human gene (mouse, rat, yeast, fly : out of HGNC scope)
      (c) A typo or invented identifier

    Filter to avoid noise on non-human files: only flag when the column
    is *predominantly* HGNC-symbol shape (>= 80% of cells are
    recognised HGNC symbols). A few unrecognised symbols in an
    otherwise-human-gene column are worth flagging; a mouse-only file
    where everything is unrecognised should NOT trigger.

    Info-level (confidence 0.20). The user adjudicates.
    """
    MIN_HGNC_FRACTION = 0.80
    MIN_CELLS = 10  # avoid tiny columns
    FILE_HUMAN_THRESHOLD = 0.50

    # Load HGNC registry once for the entire detect() call. Loading inside
    # the loop on every column scan would force a ~100 MB pandas read per
    # column and serialize the test suite. Local import keeps the
    # detector-import path cheap.
    try:
        from .corpus import load_hgnc
        _hgnc_local = load_hgnc()
    except Exception:
        _hgnc_local = None

    # File-level species check: count overall HGNC-recognised vs total
    # gene-like cells across ALL identifier columns. If the file is
    # < 50% human-recognised, it's likely a non-human-genome file
    # (mouse, fly, yeast) and our human-centric "unrecognized" flag
    # is uninformative. Skip the entire pass.
    if _hgnc_local is not None:
        total_gene_like = 0
        total_recognised = 0
        for col in df.columns:
            if str(col) not in identifier_cols:
                continue
            for value in df[col]:
                if not isinstance(value, str):
                    continue
                v = value.strip()
                if not GENE_LIKE_PATTERN.fullmatch(v):
                    continue
                total_gene_like += 1
                if (v in HGNC_NEW_SYMBOLS
                        or v in HGNC_RENAME_MAP
                        or v in _HGNC_RENAME
                        or _hgnc_local.resolve(v) is not None):
                    total_recognised += 1
        if total_gene_like >= 20 and (
            total_recognised / total_gene_like < FILE_HUMAN_THRESHOLD
        ):
            return  # skip : file is predominantly non-human

    for col in df.columns:
        col_s = str(col)
        if col_s not in identifier_cols:
            continue
        series = df[col]
        gene_like_cells: list[tuple[int, str]] = []
        for idx, value in series.items():
            if not isinstance(value, str):
                continue
            v = value.strip()
            if GENE_LIKE_PATTERN.fullmatch(v):
                gene_like_cells.append((int(idx), v))
        if len(gene_like_cells) < MIN_CELLS:
            continue
        recognised = 0
        unrecognised: list[tuple[int, str]] = []
        for idx, v in gene_like_cells:
            if (v in HGNC_NEW_SYMBOLS
                    or v in HGNC_RENAME_MAP
                    or v in _HGNC_RENAME) or _hgnc_local is not None and _hgnc_local.resolve(v) is not None:
                recognised += 1
            else:
                unrecognised.append((idx, v))
        if recognised / len(gene_like_cells) < MIN_HGNC_FRACTION:
            continue
        # Predominantly-HGNC column with a few unrecognised : flag those.
        for idx, v in unrecognised:
            if (col_s, idx) in flagged:
                continue
            report.suspicions.append(Suspicion(
                column=col_s, row=idx, value=v,
                kind="unrecognized-symbol",
                suggestion=None, confidence=0.20,
                reason=(
                    f"value {v!r} matches gene-symbol shape but is not in "
                    f"the HGNC snapshot (current symbols, prev_symbol, or "
                    f"alias). Column is {recognised}/{len(gene_like_cells)} "
                    f"= {100*recognised/len(gene_like_cells):.0f}% HGNC-"
                    f"recognised, so this symbol is an outlier. Could be a "
                    f"new gene HGNC hasn't catalogued, a non-human ortholog, "
                    f"or a typo."
                ),
                sheet=sheet,
            ))


def _detect_decimal_comma_columns(
    df: pd.DataFrame,
    report: Report,
    sheet: str | None = None,
) -> None:
    """Detect European decimal-comma encoding in numeric columns.

    Broman & Woo (2018) document locale-dependent number formatting as a
    data-integrity risk. When a file is exported from de-DE/fr-FR/nl-NL/
    it-IT Excel, numeric cells serialize `1,5` instead of `1.5`. Re-opened
    in en-US, those cells either fail to parse or become string `"1,5"`.

    Heuristic: a column where >=80% of non-empty cells match the decimal-
    comma shape AND no cell uses period-decimal flags the column. The flag
    is informational, addressed via the schema sidecar (Frictionless Table
    Schema `decimalChar`).
    """
    MIN_FRACTION = 0.80
    MIN_CELLS = 5

    for col in df.columns:
        series = df[col]
        comma_count = 0
        period_count = 0
        non_empty = 0
        comma_rows: list[int] = []
        for idx, value in series.items():
            if pd.isna(value):
                continue
            if not isinstance(value, str):
                continue
            v = value.strip()
            if not v:
                continue
            non_empty += 1
            if _DECIMAL_COMMA_RE.fullmatch(v):
                comma_count += 1
                comma_rows.append(int(idx))
            elif _DECIMAL_PERIOD_RE.fullmatch(v):
                period_count += 1
        if non_empty < MIN_CELLS:
            continue

        # sparse-contamination path. A small handful of
        # comma-decimal cells in an otherwise period-decimal numeric
        # column is a per-cell typo or copy-paste artifact. Flag each
        # comma-decimal cell at lower confidence (0.40). Threshold:
        # comma cells must be ≤10 % of non-empty and period count must
        # dominate (period_count ≥ 4 × comma_count). Without these,
        # the column is genuinely mixed/ambiguous and is left alone.
        if (period_count > 0 and comma_count > 0
                and comma_count / non_empty <= 0.10
                and period_count >= 4 * comma_count):
            for row in comma_rows:
                report.suspicions.append(Suspicion(
                    column=str(col), row=row,
                    value=str(series.loc[row]),
                    kind="decimal-comma",
                    suggestion="period-decimal (replace `,` with `.`)",
                    confidence=0.40,
                    reason=(
                        f"single-cell decimal-comma in otherwise period-"
                        f"decimal numeric column ({period_count} period-"
                        f"decimal cells, {comma_count} comma-decimal). "
                        f"Likely locale-copy-paste contamination."
                    ),
                    sheet=sheet,
                ))
            continue

        if period_count > 0:
            continue  # mixed beyond the per-cell-typo path
        if comma_count / non_empty < MIN_FRACTION:
            continue
        report.suspicions.append(Suspicion(
            column=str(col),
            row=int(series.index[0]) if len(series) else 0,
            value=f"{comma_count}/{non_empty} cells use comma-decimal",
            kind="decimal-comma",
            suggestion="Set Frictionless schema `decimalChar` to `,` "
                       "or convert numbers to period-decimal.",
            confidence=0.55,
            reason=(
                f"column {str(col)!r}: {comma_count}/{non_empty} cells "
                f"({100*comma_count/non_empty:.0f}%) match the decimal-"
                f"comma shape `\\d+,\\d+` with no period-decimal cells. "
                f"European locale Excel exports : re-opening this file in "
                f"en-US locale will silently break numeric parsing."
            ),
            sheet=sheet,
        ))


def _detect_autofill_sequences(
    df: pd.DataFrame,
    report: Report,
    flagged: set[tuple[str, int]],
    sheet: str | None = None,
) -> None:
    """Pass 3: column-level autofill-sequence detection.

    Ziemann 2021 (PMC6330011) documents the Excel autofill failure mode: a
    user enters one corrupted date, drags the fill handle, and Excel
    propagates a series of dates (`Sep-2, Sep-3, Sep-4, ...` or
    `Mar-1, Mar-2, Mar-3, ...`). Each cell *individually* gets caught by
    Pass 1, but the column-level pattern : three or more row-adjacent cells
    whose dates form a bounded monotonic series AND each maps to a gene
    candidate : is qualitatively stronger evidence than any single cell.

    Emits a single high-confidence Suspicion(kind="autofill-sequence")
    pointing at the first row of each run. Does not consume `flagged`;
    Pass 1's per-cell flags remain (they include the gene candidate),
    while this row identifies the propagation pattern itself.

    Selection criteria : both must hold:

    1. The column is classified as an identifier column. Pass 1 already
       skips measurement-classified columns (the column is *expected* to
       contain dates, so three adjacent serials aren't autofill : they're
       data). Without this gate, real-data validation on Ziemann 2016
       `23826142__asset.xlsx` produced 30 false positives in a column of
       legitimate date-serials spanning decades.

    2. Adjacent date diffs are near-uniform. Excel autofill produces
       evenly stepped increments (1 day, 1 week, ~1 month, ~1 year). Any
       monotonic gap is too permissive: long columns of legitimate dates
       landing in gene-name months trigger false positives when the only
       constraint is `diff > 0`. Uniform-spacing kills those.
    """
    MIN_RUN = 3
    MAX_ROW_GAP = 2  # rows must be directly adjacent (1) or one-apart (2)
    # Allowed adjacent-diff bands (in days). Strict : Excel autofill
    # produces exact daily/weekly/monthly/yearly increments. Loose
    # tolerances trigger false positives on legitimate sparse date
    # columns. The original (1,1,1) accepted diff∈[1,2] and risked
    # picking up coincidental same-month pairs; the tightening below
    # demands diff==1 for daily, diff==7 for weekly, etc.
    _AUTOFILL_BANDS = [
        (1, 1, 0),       # daily: diffs are exactly 1 day
        (7, 7, 0),       # weekly: diffs are exactly 7 days
        (28, 31, 1),     # monthly: median diff 28-31 (calendar variance)
        (365, 366, 1),   # yearly: median diff 365-366 (leap-year variance)
    ]

    for col in df.columns:
        series = df[col]
        is_id, _ = _column_is_identifier(str(col), series)
        if not is_id:
            continue
        date_cells: list[tuple[int, date]] = []
        for idx, value in series.items():
            d: date | None = None
            if isinstance(value, datetime):
                d = value.date()
            elif isinstance(value, date):
                d = value
            elif isinstance(value, bool):
                continue
            elif isinstance(value, int):
                if _SERIAL_MIN <= value <= _SERIAL_MAX:
                    d = _serial_to_date(value)
            elif isinstance(value, float):
                if not pd.isna(value) and value.is_integer():
                    n = int(value)
                    if _SERIAL_MIN <= n <= _SERIAL_MAX:
                        d = _serial_to_date(n)
            elif isinstance(value, str):
                parses = _parse_date_string(value)
                if len(parses) == 1:
                    d = parses[0]
            if d is None:
                continue
            if not _reverse_gene_date(d):
                continue
            date_cells.append((int(idx), d))

        if len(date_cells) < MIN_RUN:
            continue

        date_cells.sort(key=lambda p: p[0])

        def _diff_in_band(diff: int) -> tuple[int, int, int] | None:
            for lo, hi, tol in _AUTOFILL_BANDS:
                if lo <= diff <= hi + tol:
                    return (lo, hi, tol)
            return None

        runs: list[list[tuple[int, date]]] = []
        cur: list[tuple[int, date]] = [date_cells[0]]
        cur_band: tuple[int, int, int] | None = None
        for prev, nxt in zip(date_cells, date_cells[1:], strict=False):
            diff = (nxt[1] - prev[1]).days
            row_gap = nxt[0] - prev[0]
            band = _diff_in_band(diff)
            if (band is not None
                    and 1 <= row_gap <= MAX_ROW_GAP
                    and (cur_band is None or band == cur_band)):
                cur.append(nxt)
                cur_band = band
            else:
                if len(cur) >= MIN_RUN:
                    runs.append(cur)
                cur = [nxt]
                cur_band = None
        if len(cur) >= MIN_RUN:
            runs.append(cur)

        for run in runs:
            first_row = run[0][0]
            last_row = run[-1][0]
            sample = ", ".join(d.isoformat() for _, d in run[:4])
            if len(run) > 4:
                sample += ", ..."
            report.suspicions.append(Suspicion(
                column=str(col),
                row=first_row,
                value=f"rows {first_row}..{last_row} ({len(run)} cells)",
                kind="autofill-sequence",
                suggestion=None,
                confidence=0.90,
                reason=(
                    f"{len(run)} row-adjacent cells form a monotonic date "
                    f"sequence ({sample}); each cell maps to a gene-family "
                    f"member, consistent with Excel autofill propagation of "
                    f"a corrupted gene-symbol (Ziemann 2021, PMC6330011)."
                ),
                sheet=sheet,
            ))


@_cache
def _row_xref_index():
    """Lazy-loaded HGNC-backed external ID → gene-symbol index. Loaded once
    per process; ~150 MB resident. Use only inside detect_file, never inside
    the unit-test-facing detect() function."""
    from .xref_lookup import load_xref_index
    return load_xref_index()


def _apply_row_context_boost(
    report: Report,
    sheets_by_name: dict[str, pd.DataFrame],
    xref,
) -> int:
    """Boost confidence on any suspicion whose row contains an external ID
    (Ensembl/RefSeq/UniProt/Entrez/HGNC) that resolves to the suggested gene.

    This is independent corroborating evidence : much stronger than column-
    shape inference alone. Mutates report.suspicions in place and returns
    the count of boosted cells.
    """
    from .corpus import load_hgnc
    hgnc = load_hgnc()

    # Build the canonical-rename expansion set for matching old↔new symbols
    # (mirrors the logic in row_context_validator._canonical_set).
    # Case-insensitive: HGNC stores aliases in mixed case (`Oct4`) while
    # Excel emits uppercase (`OCT4`) : without case-folding the alias
    # match silently fails (calibration ).
    def _expand_for_match(suggestions: set[str]) -> set[str]:
        expanded = set(suggestions)
        for s in list(suggestions):
            canonical = hgnc.resolve(s)
            if canonical:
                expanded.add(canonical)
            for old, new in hgnc.prev_symbol_to_current.items():
                if new in (s, canonical):
                    expanded.add(old)
            for alias, new in hgnc.alias_to_current.items():
                if new in (s, canonical):
                    expanded.add(alias)
        return expanded

    n_boosted = 0
    for s in report.suspicions:
        if not s.suggestion:
            continue
        # Find the DataFrame holding this suspicion
        df = sheets_by_name.get(s.sheet) if s.sheet else None
        if df is None:
            # Single-sheet case where the loader used "_sheet0"
            if len(sheets_by_name) == 1:
                df = next(iter(sheets_by_name.values()))
            else:
                continue
        # Match by index LABEL not positional offset : `s.row` is
        # `int(idx)` from `df.iterrows()`, i.e. the actual pandas index
        # value, not a 0-based offset. Using `.iloc[s.row]` would silently
        # fetch the wrong row on dataframes with non-default index
        # (post-`reset_index`, post-`dropna`, multi-row headers, etc.).
        if s.row not in df.index:
            continue
        row = df.loc[s.row]
        resolved_in_row: dict[str, str] = {}
        for col, val in row.items():
            if str(col) == s.column or not isinstance(val, str):
                continue
            gene = xref.lookup(val.strip())
            if gene is not None:
                resolved_in_row[val.strip()] = gene
        if not resolved_in_row:
            s.xref_status = "absent"
            continue
        # Does any resolved gene match our suggestion (with HGNC rename
        # equivalences)?
        suggestion_set = {p.strip() for p in s.suggestion.split("|") if p.strip()}
        match_set = _expand_for_match(suggestion_set)
        corroborated_by: list[str] = []
        for ext_id, gene in resolved_in_row.items():
            if gene in match_set:
                corroborated_by.append(f"{ext_id!r}->{gene}")
        if corroborated_by:
            # Independent external corroboration : strong evidence.
            s.xref_status = "corroborated"
            new_conf = min(0.99, s.confidence + 0.4)
            note = "; row-xref CORROBORATED via " + ", ".join(corroborated_by[:3])
            if s.confidence != new_conf:
                s.confidence = new_conf
                s.reason = (s.reason or "") + note
                n_boosted += 1
        else:
            # Row has resolvable external IDs but none match the suggestion.
            # Empirically this is a strong false-positive signal : see the
            # `supplementary_table_8_ddac017.xlsx` row-279 case where a
            # publication-date column was mis-typed as identifier, the
            # detector decoded the date to MARCHF5, but the row's Ensembl
            # ID resolved to GARS1. Downgrade confidence to reflect the
            # contradicting evidence; the user still sees the flag for
            # adjudication, but it's no longer in the high-confidence band.
            s.xref_status = "contradicted"
            contradicting = ", ".join(
                f"{eid!r}->{g}" for eid, g in list(resolved_in_row.items())[:3]
            )
            new_conf = min(s.confidence, 0.20)
            note = (
                f"; row-xref CONTRADICTED : row resolves to "
                f"{contradicting} which does not match suggestion. "
                f"Confidence demoted from {s.confidence:.2f} to {new_conf:.2f}."
            )
            s.confidence = new_conf
            s.reason = (s.reason or "") + note
    return n_boosted


_DEFAULT_MAX_FILE_BYTES = 100 * 1024 * 1024  # 100 MB


def detect_file(
    path: str,
    row_context_boost: bool = True,
    max_file_bytes: int | None = None,
) -> Report:
    """Run detection across every sheet of a spreadsheet file.

    For xlsx with multiple sheets, this is the right entry point : empirically
    real corruption often hides in non-default sheets that single-sheet readers
    miss.

    `row_context_boost=True` (default) runs an independent cross-reference
    pass: for each suspicion, scan the SAME row for external IDs
    (Ensembl/RefSeq/UniProt/Entrez/HGNC) that resolve to the suggested gene.
    When corroborated, confidence is boosted (capped at 0.99) and the
    corroborating cell is appended to the suspicion's `reason`. Set to False
    to skip the ~150 MB xref index load : useful for fast unit tests.

    `max_file_bytes` (default 100 MB) bounds the worst case where openpyxl
    spends 30+ minutes on a single >200 MB xlsx. Files over the limit
    return a Report with one `kind="file-too-large"` suspicion and zero
    other detection : preserves chain-of-custody (the user sees the file
    was scanned and what happened) without hanging the caller. Pass
    `None` to disable the guard.
    """
    from pathlib import Path as _Path

    limit = _DEFAULT_MAX_FILE_BYTES if max_file_bytes is None else max_file_bytes
    if limit is not None:
        try:
            size = _Path(path).stat().st_size
        except OSError:
            size = -1
        if size > limit:
            report = Report(rows_scanned=0, columns_scanned=0)
            report.suspicions.append(Suspicion(
                column="", row=0, value=f"{size:,} bytes",
                kind="file-too-large", suggestion=None, confidence=0.0,
                reason=(
                    f"file {path!r} is {size:,} bytes (> {limit:,}-byte guard). "
                    f"openpyxl/xlrd load is unbounded on large files and "
                    f"empirically hangs for 30+ minutes. Re-invoke with "
                    f"`max_file_bytes=None` to override, or pre-split the "
                    f"workbook into sheet-level files."
                ),
                sheet=None,
            ))
            return report

    from .app import _load_all_sheets

    sheets = _load_all_sheets(path)
    combined = Report(rows_scanned=0, columns_scanned=0)
    for sheet_name, df in sheets.items():
        r = detect(df, sheet=sheet_name)
        combined.rows_scanned += r.rows_scanned
        combined.columns_scanned += r.columns_scanned
        for c in r.identifier_columns:
            label = f"{sheet_name}!{c}" if sheet_name != "_sheet0" else c
            combined.identifier_columns.append(label)
        combined.suspicions.extend(r.suspicions)
    if row_context_boost and combined.suspicions:
        # Narrow exception scope 
        # - xref index load failures (missing data file, ImportError) are
        #   recoverable : fall back silently to detector's own confidence.
        # - mid-loop errors during the boost itself are NOT silently
        #   swallowed; they leave the report in a partial-boost state which
        #   would be a chain-of-custody problem. Let them propagate.
        try:
            xref = _row_xref_index()
        except (FileNotFoundError, ImportError, RuntimeError):
            xref = None
        if xref is not None:
            _apply_row_context_boost(combined, sheets, xref)
    return combined
