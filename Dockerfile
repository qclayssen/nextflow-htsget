FROM python:3.12-slim

# Install system dependencies
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        curl \
        procps \
        wget \
        fastqc \
        samtools \
        bcftools \
        default-jre \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
# Install Python packages
RUN pip install --no-cache-dir \
        htsget \
        requests \
        multiqc

# Copy and setup htsget_fetch.py script
COPY bin/htsget_fetch.py /usr/local/bin/htsget_fetch.py
RUN chmod +x /usr/local/bin/htsget_fetch.py
