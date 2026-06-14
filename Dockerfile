# UnCorrupt - reproducible container image
#
# This image ships the UnCorrupt Python package, the Gradio UI, and the
# bundled registry data so it runs offline. Default command launches the
# UI on port 7860; override the CMD to run the CLI:
#
#   # UI (drop-in for the HuggingFace Space)
#   docker run --rm -p 7860:7860 ghcr.io/shitcoinsherpa/uncorrupt:<version>
#
#   # CLI on a local file
#   docker run --rm -v "$(pwd)":/work \
#       ghcr.io/shitcoinsherpa/uncorrupt:<version> \
#       uncorrupt detect /work/yours.xlsx
#
#   # Write the schema sidecar
#   docker run --rm -v "$(pwd)":/work \
#       ghcr.io/shitcoinsherpa/uncorrupt:<version> \
#       uncorrupt schema /work/yours.xlsx
#
# Base image is pinned by SHA256 digest. The Python deps come from
# pyproject.toml (resolved at build time); the registry data is bundled
# inside the package via hatch's shared-data, so no large external corpus
# is required for normal use.

FROM python:3.14.0-slim-bookworm@sha256:d13fa0424035d290decef3d575cea23d1b7d5952cdf429df8f5542c71e961576

LABEL org.opencontainers.image.title="UnCorrupt"
LABEL org.opencontainers.image.description="Detection and repair of Excel-corrupted gene symbols in genomics spreadsheets."
LABEL org.opencontainers.image.licenses="Apache-2.0"
LABEL org.opencontainers.image.documentation="https://github.com/shitcoinsherpa/UnCorrupt/blob/main/README.md"
LABEL org.opencontainers.image.source="https://github.com/shitcoinsherpa/UnCorrupt"
LABEL org.opencontainers.image.authors="LLMSherpa (https://x.com/LLMSherpa)"
LABEL org.opencontainers.image.base.name="docker.io/library/python:3.12.13-slim-bookworm"
LABEL org.opencontainers.image.base.digest="sha256:d193c6f51a7dbd10395d6328de3a7edb0516fb0608ca138036576f574c3e07d2"

ENV PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    GRADIO_SERVER_NAME=0.0.0.0 \
    GRADIO_SERVER_PORT=7860 \
    GRADIO_ANALYTICS_ENABLED=False \
    HF_HOME=/tmp/huggingface

# Runtime libs only. No build-essential / cmake: nothing here compiles from
# source (the optional [llm] extra is research-only and is not installed).
RUN apt-get update && apt-get install -y --no-install-recommends \
        libxml2 libxslt1.1 libomp5 libgomp1 \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Non-root execution: standard hardening. UnCorrupt reads user-supplied
# spreadsheets; root inside the container is one CVE away from root on
# the host bind mount.
RUN groupadd --system --gid 1000 sci \
    && useradd --system --uid 1000 --gid sci --create-home --home-dir /home/sci sci

WORKDIR /app

# Build context, ordered least-to-most-likely-to-change.
COPY --chown=sci:sci pyproject.toml README.md LICENSE NOTICE /app/
COPY --chown=sci:sci src /app/src

# Install the package itself. Registry data is bundled inside the wheel
# via hatch's shared-data / force-include (see pyproject.toml), so the
# detector runs offline with no external corpus.
RUN pip install --upgrade pip \
    && pip install /app

USER sci

EXPOSE 7860

HEALTHCHECK --interval=30s --timeout=5s --start-period=15s --retries=3 \
    CMD python -c "from uncorrupt.detector import detect_file; from uncorrupt.app import build_app" || exit 1

# Default: launch the Gradio UI on port 7860. Override to run the CLI.
CMD ["uncorrupt-app"]
