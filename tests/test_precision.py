"""Unit tests for the precision-measurement module.

Network-dependent calls (esearch, EPMC fetch) are exercised by the live run,
not unit tests. We test the structural pieces here so CI doesn't depend on
NCBI uptime.
"""
from __future__ import annotations

from uncorrupt.corpus import load_ziemann_2021_corpus
from uncorrupt.precision import _exclude_ziemann_set, GENETICS_QUERY_TERM


def test_exclude_ziemann_set_returns_numeric_ids() -> None:
    """All excluded IDs should be numeric (no 'PMC' prefix)."""
    excluded = _exclude_ziemann_set()
    assert len(excluded) > 3000  # ~5086 in S2
    for pmc_id in list(excluded)[:50]:
        assert pmc_id.isdigit(), f"non-numeric id leaked through: {pmc_id!r}"


def test_exclude_ziemann_set_matches_corpus_size() -> None:
    excluded = _exclude_ziemann_set()
    entries = load_ziemann_2021_corpus()
    # Must match unique PMC ID count (some duplicates possible but small)
    unique_pmcs = {e.pmc_id for e in entries}
    assert len(excluded) == len(unique_pmcs)


def test_genetics_query_term_filters_to_open_access() -> None:
    """The query intentionally restricts to open-access content for fetchability."""
    assert "open access[filter]" in GENETICS_QUERY_TERM
    assert "[Title/Abstract]" in GENETICS_QUERY_TERM
