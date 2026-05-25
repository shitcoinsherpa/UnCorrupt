"""Homoglyph detection for ASCII-targeting identifier text.

HGNC stores gene symbols in ASCII Latin. Any non-ASCII character in a
gene-symbol field is an error by construction: it could be a Greek capital
Σ (U+03A3) that looks like Latin S, a Cyrillic А (U+0410) that looks like
Latin A, a fullwidth digit (U+FF11), and so on. In every such case the
correct response is to flag the cell : never silently overwrite : and emit
the proposed ASCII repair so the scientist can adjudicate.

The lookup table is built once at import time from Unicode Technical
Standard 39's `confusables.txt` (via the `confusable_homoglyphs` package),
with one manual addition: TR39 maps Greek capital Sigma to Latin ESH
(U+01A9), which is not ASCII; we override to Latin capital S since that's
what the printed glyph looks like and what HGNC uses.
"""
from __future__ import annotations

import unicodedata
from typing import NamedTuple

from confusable_homoglyphs import confusables


class HomoglyphRepair(NamedTuple):
    """A single character substitution in a repaired symbol."""

    original: str
    codepoint: str  # `U+XXXX` form
    name: str
    replacement: str


_HOMOGLYPH_MAP: dict[int, str] = {}


def _build_table() -> None:
    """Populate the lookup table from TR39. Runs once at import time."""
    # Stay inside the BMP : confusable_homoglyphs queries above 0x10000
    # are expensive and rarely productive for ASCII-target lookups.
    for codepoint in range(0x0080, 0x10000):
        try:
            ch = chr(codepoint)
            result = confusables.is_confusable(ch, preferred_aliases=["latin"])
        except Exception:
            continue
        if not result:
            continue
        # Prefer ASCII Latin uppercase letter; accept ASCII digit too
        # (gene symbols are letters + digits, e.g. fullwidth `１` in `BRCA１`
        # confuses with ASCII `1`).
        for h in result[0]["homoglyphs"]:
            target = h["c"]
            if (len(target) == 1 and ord(target) < 128
                    and (target.isupper() or target.isdigit())):
                _HOMOGLYPH_MAP[codepoint] = target
                break
    # Manual gap (per UTS #39): Greek capital Sigma (U+03A3) maps to Latin
    # capital ESH (U+01A9), which is not ASCII. The printed glyph and
    # HGNC convention is plain Latin S : override.
    _HOMOGLYPH_MAP[0x03A3] = "S"
    # Same family: capital Eta (U+0397) → Latin H, capital Iota (U+0399) → I,
    # capital Kappa (U+039A) → K, capital Mu (U+039C) → M, capital Nu (U+039D) → N,
    # capital Omicron (U+039F) → O, capital Rho (U+03A1) → P. Most are present
    # already; these explicit overrides defend against TR39 data drift.
    _HOMOGLYPH_MAP.setdefault(0x0397, "H")
    _HOMOGLYPH_MAP.setdefault(0x0399, "I")
    _HOMOGLYPH_MAP.setdefault(0x039A, "K")
    _HOMOGLYPH_MAP.setdefault(0x039C, "M")
    _HOMOGLYPH_MAP.setdefault(0x039D, "N")
    _HOMOGLYPH_MAP.setdefault(0x039F, "O")
    _HOMOGLYPH_MAP.setdefault(0x03A1, "P")

    # Fullwidth digits and fullwidth Latin letters (CJK forms,
    # U+FF10..U+FF19 and U+FF21..U+FF3A) are visually identical to ASCII
    # 0-9 and A-Z but live in the CJK Compatibility Forms block. TR39's
    # confusables.txt happens to omit fullwidth digits, but they DO appear
    # in scientific data (CJK Excel inputs that round-trip through a
    # locale-aware tool). Manual override.
    for digit in range(10):
        _HOMOGLYPH_MAP.setdefault(0xFF10 + digit, str(digit))
    for i in range(26):
        _HOMOGLYPH_MAP.setdefault(0xFF21 + i, chr(ord("A") + i))


_build_table()


def repair(text: str) -> tuple[str, list[HomoglyphRepair]]:
    """Return (repaired, repairs).

    Repairs lists each per-character substitution that was applied. Empty
    list means the input was already pure ASCII and no work was done.
    """
    if not text or text.isascii():
        return text, []
    out: list[str] = []
    repairs: list[HomoglyphRepair] = []
    for ch in text:
        cp = ord(ch)
        if cp > 127 and cp in _HOMOGLYPH_MAP:
            replacement = _HOMOGLYPH_MAP[cp]
            repairs.append(HomoglyphRepair(
                original=ch,
                codepoint=f"U+{cp:04X}",
                name=unicodedata.name(ch, "<unnamed>"),
                replacement=replacement,
            ))
            out.append(replacement)
        else:
            out.append(ch)
    return "".join(out), repairs


_RELAXED_GENE_SHAPE = __import__("re").compile(
    # HGNC symbols are *typically* uppercase A-Z + digits + hyphen, but
    # the canonical convention also admits a handful of dotted forms
    # (`RP-117C8.1`, `MIR205HG`, lncRNA-style names) and lowercase
    # `orf`/`as`/`os` suffixes. Relax the target check beyond
    # GENE_LIKE_PATTERN so homoglyph corruption of those is also caught.
    r"^[A-Z][A-Za-z0-9.\-]{1,18}$"
)


def looks_like_gene_after_repair(text: str, gene_like_pattern) -> tuple[bool, str, list[HomoglyphRepair]]:
    """Check whether `text` matches the gene-like shape after homoglyph
    normalization. Returns (matches, repaired_text, repairs).

    Two passes:
      1. Strict match against the caller's GENE_LIKE_PATTERN (uppercase +
         digits + hyphen). High-confidence "this is a normal gene symbol
         corrupted by a Greek/Cyrillic letter."
      2. Fallback against a *relaxed* pattern that also admits dotted
         lncRNA-style names and mixed-case suffixes. Reported back to the
         caller so it can apply a lower confidence.
    """
    if text.isascii():
        return False, text, []  # nothing to do
    repaired, repairs = repair(text)
    if not repairs:
        return False, text, []
    if gene_like_pattern.fullmatch(repaired):
        return True, repaired, repairs
    # Fallback to relaxed pattern. Still returns matches=True so the
    # caller emits a flag; the caller can downgrade confidence based on
    # whether the strict pattern matched.
    if _RELAXED_GENE_SHAPE.fullmatch(repaired):
        return True, repaired, repairs
    return False, repaired, repairs
