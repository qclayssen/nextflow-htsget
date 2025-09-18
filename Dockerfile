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
RUN pip install --no-cache-dir \
        htsget==0.1.0a2 \
        requests \
        multiqc

# Copy and setup fetch_htsget.py script
COPY bin/fetch_htsget.py /bin/fetch_htsget.py
RUN chmod +x /bin/fetch_htsget.py
