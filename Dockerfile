FROM python:3.12-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bash \
        curl \
        procps \
    fastqc \
    samtools \
    bcftools \
    openjdk-21-jre-headless \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
        requests \
        htsget \
        multiqc

# Copy and setup htsget_fetch.py script
COPY bin/htsget_fetch.py /usr/local/bin/htsget_fetch.py
RUN chmod +x /usr/local/bin/htsget_fetch.py