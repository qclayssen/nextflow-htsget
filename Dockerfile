FROM python:3.12-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends curl procps \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
        requests \
        htsget
