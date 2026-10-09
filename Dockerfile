FROM python:3.11-slim AS basis

# set environment varibles
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1


FROM basis AS builder
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /usr/local/bin/

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       build-essential \
       fasttree \
       libxml2 \
       mafft \
    && rm -rf /var/lib/apt/lists/*

ENV UV_PROJECT_ENVIRONMENT=/app \
    UV_PYTHON_DOWNLOADS=never \
    UV_LOCKED=1

WORKDIR /src

# Install dependencies only, before copying the source, so this (slow) layer
# is cached across source-only changes.
COPY pyproject.toml uv.lock ./
RUN uv sync --no-install-project --all-extras

# Install the project itself. --no-editable is required: /app is copied into
# the runtime image while /src is discarded, so the venv must be self-contained.
COPY . .
RUN uv sync --no-editable --all-extras


FROM basis AS runtime
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       fasttree \
       libxml2 \
       mafft \
       mmseqs2 \
       procps \
       time \
    && apt-get -y autoremove \
    && apt-get -y autoclean \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /app /app
ENV PATH="/app/bin:$PATH"

RUN python -c "import FastOMA; print(FastOMA.__version__)"