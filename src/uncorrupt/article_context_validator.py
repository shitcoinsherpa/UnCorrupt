"""Article-context validation for inconclusive detector flags.

When row-context xref can't classify a flag (the 92% of flags whose row
doesn't contain external IDs to cross-check), we fall back to a second
ground-truth signal: does the suggested gene appear anywhere in the
article's full text?

This is a NECESSARY but not SUFFICIENT signal:
- Gene NOT in article text → almost certainly a false positive
- Gene IS in article text → consistent with real corruption, doesn't prove it

We can therefore tighten precision without false-negative risk on the
"gene NOT mentioned" side. The "gene IS mentioned" cases remain
inconclusive for now (sub-agent task next).

Channel: Europe PMC `fullTextXML` endpoint. Returns JATS XML.
"""
from __future__ import annotations

import json
import re
import time
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from xml.etree import ElementTree as ET

import httpx

from .corpus import load_hgnc

PROJECT_ROOT = Path(__file__).resolve().parents[2]
ARTICLE_TEXT_DIR = PROJECT_ROOT / "data" / "raw" / "ziemann_2021_corpus" / "article_texts"
LEDGER_PATH = PROJECT_ROOT / "results" / "article_context_validation.jsonl"

EPMC_FULLTEXT_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc_id}/fullTextXML"
USER_AGENT = "uncorrupt-research/0.1 (https://github.com/shitcoinsherpa)"


@dataclass
class ArticleContextResult:
    pmc_id: str
    suggestion: str
    expanded_genes: list[str]
    article_text_chars: int
    matched_genes: list[str]  # which expanded forms appeared in article text
    outcome: str  # "TP-article", "FP-article", "fetch-failed", "no-genes-in-suggestion"


@dataclass
class ArticleContextRun:
    timestamp: str
    n_inconclusive_input: int
    n_articles_fetched: int
    n_tp_article: int
    n_fp_article: int
    n_fetch_failed: int
    article_context_precision: float
    runtime_seconds: float
    results: list[ArticleContextResult] = field(default_factory=list)


def _fetch_article_text(pmc_id: str, client: httpx.Client) -> str | None:
    """Fetch article full-text XML from Europe PMC, return extracted plain text."""
    cache_path = ARTICLE_TEXT_DIR / f"{pmc_id}.txt"
    if cache_path.exists():
        return cache_path.read_text(encoding="utf-8")

    url = EPMC_FULLTEXT_URL.format(pmc_id=pmc_id)
    try:
        resp = client.get(url, follow_redirects=True, timeout=60)
        if resp.status_code != 200:
            return None
        if len(resp.content) < 500:
            return None
    except httpx.HTTPError:
        return None

    # Parse JATS XML and pull plain text from <body>, <abstract>, <title>
    try:
        root = ET.fromstring(resp.content)
    except ET.ParseError:
        return None

    pieces: list[str] = []
    for elem_path in (".//title", ".//abstract", ".//body", ".//article-title"):
        for elem in root.iter():
            if elem.tag.endswith(elem_path[3:]):
                if elem.text:
                    pieces.append(elem.text)
                for child in elem.iter():
                    if child.text:
                        pieces.append(child.text)
                    if child.tail:
                        pieces.append(child.tail)
    text = " ".join(pieces)

    ARTICLE_TEXT_DIR.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(text, encoding="utf-8")
    return text


def _expand_gene_to_canonical_set(suggestion: str, hgnc) -> list[str]:
    """Expand a suggestion (e.g. 'SEPT2 | SEP2') to all HGNC-equivalent forms."""
    genes = {g.strip() for g in suggestion.split("|") if g.strip()}
    expanded: set[str] = set()
    for g in genes:
        expanded.add(g)
        if g in hgnc.prev_symbol_to_current:
            expanded.add(hgnc.prev_symbol_to_current[g])
        if g in hgnc.alias_to_current:
            expanded.add(hgnc.alias_to_current[g])
        for old, new in hgnc.prev_symbol_to_current.items():
            if new == g:
                expanded.add(old)
        for alias, new in hgnc.alias_to_current.items():
            if new == g:
                expanded.add(alias)
    return sorted(expanded)


def _gene_in_text(gene: str, text: str) -> bool:
    """Word-boundary, case-sensitive match. Gene symbols are case-sensitive
    (MARCH1 ≠ March1) so we don't lowercase."""
    if not gene or not text:
        return False
    # \b at start; gene name may end in digit so word-boundary at end is fine
    pattern = r"\b" + re.escape(gene) + r"\b"
    return bool(re.search(pattern, text))


def validate_inconclusive_flags_against_articles(
    inconclusive_flags: list[dict],
    rate_limit_per_second: float = 2.0,
) -> ArticleContextRun:
    """For each inconclusive flag, check if its suggested gene appears in the
    article's full text.

    Each input flag dict should have keys: pmc_id, suggestion.
    """
    hgnc = load_hgnc()
    results: list[ArticleContextResult] = []
    delay = 1.0 / rate_limit_per_second
    t0 = time.monotonic()

    # Group by pmc_id to avoid re-fetching
    by_pmc: dict[str, list[dict]] = {}
    for f in inconclusive_flags:
        by_pmc.setdefault(f["pmc_id"], []).append(f)

    headers = {"User-Agent": USER_AGENT}
    with httpx.Client(headers=headers) as client:
        n_pmc = len(by_pmc)
        for i, (pmc_id, flags) in enumerate(by_pmc.items(), 1):
            text = _fetch_article_text(pmc_id, client)
            time.sleep(delay)
            for f in flags:
                sugg = f.get("suggestion") or ""
                expanded = _expand_gene_to_canonical_set(sugg, hgnc)
                if not expanded:
                    results.append(ArticleContextResult(
                        pmc_id=pmc_id, suggestion=sugg, expanded_genes=[],
                        article_text_chars=0, matched_genes=[],
                        outcome="no-genes-in-suggestion",
                    ))
                    continue
                if text is None:
                    results.append(ArticleContextResult(
                        pmc_id=pmc_id, suggestion=sugg, expanded_genes=expanded,
                        article_text_chars=0, matched_genes=[],
                        outcome="fetch-failed",
                    ))
                    continue
                matched = [g for g in expanded if _gene_in_text(g, text)]
                outcome = "TP-article" if matched else "FP-article"
                results.append(ArticleContextResult(
                    pmc_id=pmc_id, suggestion=sugg, expanded_genes=expanded,
                    article_text_chars=len(text), matched_genes=matched,
                    outcome=outcome,
                ))
            if i % 20 == 0:
                import sys
                tp = sum(1 for r in results if r.outcome == "TP-article")
                fp = sum(1 for r in results if r.outcome == "FP-article")
                fail = sum(1 for r in results if r.outcome == "fetch-failed")
                print(f"[article-ctx] {i}/{n_pmc} pmcids, flags={len(results)}, "
                      f"TP={tp} FP={fp} fail={fail}", file=sys.stderr, flush=True)

    n_fetched = len({r.pmc_id for r in results if r.article_text_chars > 0})
    n_tp = sum(1 for r in results if r.outcome == "TP-article")
    n_fp = sum(1 for r in results if r.outcome == "FP-article")
    n_fail = sum(1 for r in results if r.outcome == "fetch-failed")
    grounded = n_tp + n_fp
    precision = n_tp / max(1, grounded)

    return ArticleContextRun(
        timestamp=datetime.now(UTC).isoformat(),
        n_inconclusive_input=len(inconclusive_flags),
        n_articles_fetched=n_fetched,
        n_tp_article=n_tp, n_fp_article=n_fp,
        n_fetch_failed=n_fail,
        article_context_precision=precision,
        runtime_seconds=time.monotonic() - t0,
        results=results,
    )


def append_to_ledger(run: ArticleContextRun, path: Path = LEDGER_PATH) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(json.dumps(asdict(run)) + "\n")
