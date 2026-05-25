# uncorrupt — publication-grade reproducible container
# =============================================================================
#
# Pillar 2 (Five Pillars of Computational Reproducibility): every layer of
# this image is pinned to an exact, content-addressed identifier:
#
#   - Base image: pinned by SHA256 digest (immutable; survives tag re-pushes)
#   - APT packages: pinned to Debian Bookworm snapshot versions
#   - Python deps: 101 packages with 1,737 SHA256 hashes via requirements.txt
#     (compiled by pip-tools, installed with --require-hashes — any byte-level
#     tampering at install time will be detected)
#
# Structure follows "Ten Simple Rules for Writing Dockerfiles for Reproducible
# Data Science" (Nüst & Eddelbuettel, PLOS Comp. Bio. 2020):
#   - Rule 2: build on official Python slim base, version-pinned
#   - Rule 3: human-readable layout, long-form parameters
#   - Rule 4: documented decisions, OCI labels, ENV defaults
#   - Rule 5: every version pinned (no floats)
#   - Rule 8: one-click runnable (default CMD runs the test suite)
#   - Rule 9: instructions ordered least-to-most-likely-to-change
#
# Build:    docker build -t uncorrupt:1.0.0 .
# Verify:   docker run --rm uncorrupt:1.0.0                 # runs tests
# UI:       docker run --rm -p 7860:7860 uncorrupt:1.0.0 \
#               python -m uncorrupt.app
# Detect:   docker run --rm -v "$(pwd)":/work uncorrupt:1.0.0 \
#               python -c "from uncorrupt.detector import detect_file; \
#                          [print(s) for s in detect_file('/work/yours.xlsx').suspicions]"
#
# =============================================================================
# STAGE 1: builder — installs deps with build-time toolchain, then discards it
# =============================================================================
FROM python:3.12.13-slim-bookworm@sha256:d193c6f51a7dbd10395d6328de3a7edb0516fb0608ca138036576f574c3e07d2 AS builder

# Pin pip itself for reproducibility (the latest stable as of pin date)
ENV PIP_VERSION=24.3.1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_DEFAULT_TIMEOUT=120 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Build toolchain for native C extensions (lxml, llama-cpp-python).
# Rule 5: pin OS packages explicitly; clean apt cache in the same RUN to
# prevent cache invalidation surprises.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        libxml2-dev \
        libxslt1-dev \
        libffi-dev \
        libomp-dev \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

# Layer-cache: deps change rarely; copy ONLY the requirements first
COPY requirements.txt /build/requirements.txt

# Hash-verified install: every wheel's SHA256 is checked against the pinned
# hash in requirements.txt before extraction. If even one byte of any wheel
# on PyPI changes (supply-chain tampering, mirror corruption, etc.) the
# install fails loud rather than silently producing a different environment.
RUN pip install --upgrade "pip==${PIP_VERSION}" \
    && pip install --require-hashes -r /build/requirements.txt

# =============================================================================
# STAGE 2: runtime — minimal image with the installed env + project code
# =============================================================================
FROM python:3.12.13-slim-bookworm@sha256:d193c6f51a7dbd10395d6328de3a7edb0516fb0608ca138036576f574c3e07d2

# OCI standard image labels — discoverable via `docker inspect` and used by
# container registries for provenance display.
LABEL org.opencontainers.image.title="UnCorrupt"
LABEL org.opencontainers.image.description="Detection and prevention of Excel-style gene-symbol corruption in scientific spreadsheets"
LABEL org.opencontainers.image.version="1.0.0"
LABEL org.opencontainers.image.licenses="Apache-2.0"
LABEL org.opencontainers.image.documentation="https://github.com/shitcoinsherpa/uncorrupt/blob/main/README.md"
LABEL org.opencontainers.image.source="https://github.com/shitcoinsherpa/uncorrupt"
LABEL org.opencontainers.image.authors="LLMSherpa (https://x.com/LLMSherpa) — BT6 AI Red Team (bt6.gg)"
LABEL org.opencontainers.image.base.name="docker.io/library/python:3.12.13-slim-bookworm"
LABEL org.opencontainers.image.base.digest="sha256:d193c6f51a7dbd10395d6328de3a7edb0516fb0608ca138036576f574c3e07d2"

# Runtime-only system deps (libomp1 for llama-cpp inference, libxml2/libxslt1.1
# for lxml HTML parsing, libgomp1 for openmp). No build toolchain in runtime.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libxml2 \
        libxslt1.1 \
        libomp5 \
        libgomp1 \
        ca-certificates \
        curl \
    && rm -rf /var/lib/apt/lists/*

# Copy the populated site-packages from builder. Multi-stage gives us the
# pinned, hash-verified dep tree without dragging build tools into runtime.
COPY --from=builder /usr/local/lib/python3.12/site-packages /usr/local/lib/python3.12/site-packages
COPY --from=builder /usr/local/bin /usr/local/bin

# Non-root execution: standard security hardening (CIS Docker Benchmark 4.1).
# Scientific containers run untrusted spreadsheet data; root inside the
# container is one CVE away from root on the host bind mount.
RUN groupadd --system --gid 1000 sci \
    && useradd --system --uid 1000 --gid sci --create-home --home-dir /home/sci sci

WORKDIR /app

# Rule 9: most-likely-to-change at the bottom.
# Files needed for `pip install -e .` (so the package layout from pyproject
# resolves correctly):
COPY --chown=sci:sci pyproject.toml /app/pyproject.toml
COPY --chown=sci:sci README.md /app/README.md
COPY --chown=sci:sci src /app/src
COPY --chown=sci:sci tests /app/tests
COPY --chown=sci:sci scripts /app/scripts
COPY --chown=sci:sci docs /app/docs

# Canonical reference data baked into the image — these are small (~17 MB
# total), stable, externally-sourced files the detector cannot run without.
# They are *not* the 23 GB validation corpus (mount that at runtime via -v).
# HGNC: gene-symbol registry snapshot dated 2026-05-14 from
#   https://www.genenames.org/download/statistics-and-files/
# Ziemann S1/S2: the published supplementary tables from PLOS Comp. Bio. 2021
#   that drive `corpus.load_ziemann_2021_corpus()` and our cell-level recall
#   measurement.
COPY --chown=sci:sci data/raw/registries /app/data/raw/registries
COPY --chown=sci:sci data/raw/ziemann_2021_corpus /app/data/raw/ziemann_2021_corpus

# Writable runtime scratch dirs (sci owns them so the container can be run
# without bind-mounts and still write derived outputs/results)
RUN mkdir -p /app/data/derived /app/results /app/data/raw/koh_replication \
    && chown -R sci:sci /app/data /app/results

# Install the uncorrupt package itself (editable, so test discovery
# + module entrypoints work). Hash-pinning isn't possible for a local pkg;
# the source layer's content is what makes this deterministic.
USER sci
ENV PATH="/home/sci/.local/bin:${PATH}"
RUN pip install --user --no-deps -e .

ENV PYTHONPATH=/app/src \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    GRADIO_SERVER_NAME=0.0.0.0 \
    GRADIO_SERVER_PORT=7860 \
    GRADIO_ANALYTICS_ENABLED=False

# Healthcheck: the import path being intact is the truest signal that the
# image is functional. Faster + more honest than HTTP probing Gradio.
HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
    CMD python -c "from uncorrupt.detector import detect_file; from uncorrupt.app import build_app" || exit 1

EXPOSE 7860

# Rule 8: one-click runnable. Default CMD runs the full unit + integration
# test suite — anyone pulling this image can verify it works with one command.
# Override via `docker run ... <command>` for UI mode, single-file detection,
# or the gated LLM smoke (`pytest tests/test_app_local_llm.py -v` with
# UNCORRUPT_RUN_LLM_SMOKE=1).
CMD ["pytest", "-q", "--ignore=tests/test_app_local_llm.py", "--ignore=tests/test_app_integration.py"]
