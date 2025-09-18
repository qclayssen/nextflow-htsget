FROM python:3.12-slim

# Install system dependencies
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        curl \
        procps \
        wget \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip install --no-cache-dir \
        requests \
        htsget

# Install htsget CLI tool
RUN wget -O /usr/local/bin/htsget https://github.com/ga4gh/htsget-refserver/releases/download/v1.6.2/htsget-linux-amd64 \
    && chmod +x /usr/local/bin/htsget

# Add local bin directory to PATH
ENV PATH="/opt/bin:${PATH}"

# Copy local scripts
COPY bin/ /opt/bin/

# Set working directory
WORKDIR /work
