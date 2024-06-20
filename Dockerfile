FROM python:3.11-slim as basis

# set environment varibles
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1


FROM basis as builder
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       build-essential \
       fasttree \
       libxml2 \
       mafft \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
RUN pip install --upgrade hatch pip wheel
COPY pyproject.toml .
RUN python -m venv /app \
    && hatch dep show requirements --all > requirements.txt \
    && /app/bin/pip install -r requirements.txt

COPY . .
RUN ls -la \
    && hatch build \
    && ls -la dist/ \
    && /app/bin/pip install dist/*.whl


FROM basis as runtime
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       fasttree \
       libxml2 \
       mafft \
       mmseqs2 \
       procps \
    && apt-get -y autoremove \
    && apt-get -y autoclean \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /app /app
ENV PATH="/app/bin:$PATH"
