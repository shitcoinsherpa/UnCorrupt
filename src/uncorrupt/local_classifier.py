"""Local LLM classifier for review packets : replaces external sub-agent calls.

Uses llama-cpp-python with a quantized Qwen2.5-3B-Instruct model. Runs
in-process (no separate daemon), pip-installable, CPU/GPU agnostic.

Model is downloaded on first use (~2.1 GB GGUF, cached in HF Hub cache).
"""
from __future__ import annotations

import json
import re
import time
from dataclasses import dataclass
from functools import cache
from typing import Literal

Verdict = Literal["TP", "FP", "INCONCLUSIVE"]

DEFAULT_MODEL_REPO = "Qwen/Qwen2.5-3B-Instruct-GGUF"
DEFAULT_MODEL_FILE = "*q4_k_m.gguf"

SYSTEM_PROMPT = """You are a strict reviewer for an Excel-style gene-name-corruption \
detector. You receive Cell Review Packets : structured evidence about one \
flagged cell in one supplementary spreadsheet. You return one of three \
verdicts and a one-sentence rationale.

Verdicts:
- TP: the cell IS a real Excel-style gene-name corruption (a gene symbol \
like SEPT3, MARCH10, OCT2, or DEC1 was auto-converted by Excel to a \
date/number/string). The detector's suggestion is the original gene.
- FP: the cell is legitimate data : what is stored is what authors intended. \
The detector was wrong to flag it.
- INCONCLUSIVE: insufficient evidence either way.

Decision rules:
- Strong TP signal: the column is dominated by gene symbols AND the cell \
matches Excel's auto-convert pattern for a gene from the same family as \
the suggestion.
- Row-context corroboration: RefSeq IDs (NM_/NR_/NP_), Ensembl IDs \
(ENSG/ENST/ENSMUSG), UniProt accessions, or explicit gene-symbol aliases \
in the same row that match the suggestion confirm TP.
- Multi-species: mouse (TitleCase), zebrafish/yeast (lowercase) gene \
symbols are equally vulnerable to Excel corruption. Suggestions are in \
human-uppercase form.
- Integer cells decoding to uncorrupt dates in gene-symbol columns \
ARE strong TP signals (Excel stores dates as serial integers).
- Excel serial 41897 = 2014-09-15 (Windows epoch). Suggested gene must \
match the *day-of-month* of the decoded date.

Respond in this exact format:

VERDICT: <TP|FP|INCONCLUSIVE>
REASON: <one sentence>
"""


@dataclass
class ClassificationResult:
    verdict: Verdict
    reason: str
    raw_response: str
    inference_seconds: float


@cache
def _load_model(repo_id: str = DEFAULT_MODEL_REPO,
                filename: str = DEFAULT_MODEL_FILE,
                n_ctx: int = 4096,
                n_threads: int | None = None):
    """Lazy-load the GGUF model. Cached for the lifetime of the process."""
    try:
        from llama_cpp import Llama
    except ImportError as exc:
        raise RuntimeError(
            "llama-cpp-python not installed. "
            "Install with: pip install llama-cpp-python"
        ) from exc
    return Llama.from_pretrained(
        repo_id=repo_id,
        filename=filename,
        n_ctx=n_ctx,
        n_threads=n_threads,
        verbose=False,
    )


def classify_packet(packet_md: str,
                    model_repo: str = DEFAULT_MODEL_REPO,
                    model_file: str = DEFAULT_MODEL_FILE,
                    temperature: float = 0.0,
                    max_tokens: int = 256) -> ClassificationResult:
    """Classify a single review packet using the local LLM."""
    llm = _load_model(model_repo, model_file)
    t0 = time.monotonic()
    resp = llm.create_chat_completion(
        messages=[
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": packet_md},
        ],
        temperature=temperature,
        max_tokens=max_tokens,
    )
    elapsed = time.monotonic() - t0
    raw = resp["choices"][0]["message"]["content"]
    verdict, reason = _parse_response(raw)
    return ClassificationResult(
        verdict=verdict, reason=reason,
        raw_response=raw, inference_seconds=elapsed,
    )


def _parse_response(raw: str) -> tuple[Verdict, str]:
    """Parse model output into (verdict, reason). Accepts multiple formats:

    - Strict: "VERDICT: TP\nREASON: ..."
    - Loose: "TP: ..." or "TP\nReason..."
    - Embedded: any line starting with TP/FP/INCONCLUSIVE word
    """
    verdict: Verdict = "INCONCLUSIVE"
    reason = ""

    # 1. Strict VERDICT: format
    m = re.search(r"VERDICT\s*:\s*(TP|FP|INCONCLUSIVE)", raw, re.IGNORECASE)
    if m:
        verdict = m.group(1).upper()  # type: ignore
        rm = re.search(r"REASON\s*:\s*(.+?)(?:\n\n|$)", raw, re.IGNORECASE | re.DOTALL)
        if rm:
            reason = rm.group(1).strip().split("\n")[0].strip()
        return verdict, reason

    # 2. Loose: line starting with the verdict word
    for line in raw.splitlines():
        s = line.strip()
        m2 = re.match(r"^\**\s*(TP|FP|INCONCLUSIVE)\b\**\s*[:\-]?\s*(.*)$",
                      s, re.IGNORECASE)
        if m2:
            verdict = m2.group(1).upper()  # type: ignore
            rest = m2.group(2).strip()
            if rest:
                reason = rest
            return verdict, reason

    # 3. Last resort: keyword anywhere in response
    m3 = re.search(r"\b(TP|FP|INCONCLUSIVE)\b", raw, re.IGNORECASE)
    if m3:
        verdict = m3.group(1).upper()  # type: ignore
    return verdict, raw[:200].strip()


def classify_packets_batch(packets: list[str],
                           progress_callback=None) -> list[ClassificationResult]:
    """Classify multiple packets sequentially. Reuses loaded model."""
    results: list[ClassificationResult] = []
    for i, pkt in enumerate(packets, 1):
        results.append(classify_packet(pkt))
        if progress_callback:
            progress_callback(i, len(packets), results[-1])
    return results
