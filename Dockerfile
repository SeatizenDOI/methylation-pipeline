# Use a lightweight Python image
FROM python:3.10-slim

# Set environment variables to avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    tar \
    perl \
    bowtie2 \
    samtools \
    default-jre \
    cutadapt \
    fastqc \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Create a non-root user (e.g., 'defaultuser')
RUN useradd -m -s /bin/bash defaultuser && \
    mkdir -p /home/methylation/tools /data && \
    chown -R defaultuser:defaultuser /home/methylation /data

# Set working directory
WORKDIR /home/methylation/tools

# Install TrimGalore!
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -o trim_galore.tar.gz && \
    tar xvzf trim_galore.tar.gz && \
    mv TrimGalore-* TrimGalore && \
    chmod -R u+rwx TrimGalore && \
    chown -R defaultuser:defaultuser TrimGalore && \
    ln -s /home/methylation/tools/TrimGalore/trim_galore /usr/local/bin/trim_galore

# Install Bismark
RUN curl -fsSL https://github.com/FelixKrueger/Bismark/archive/refs/tags/v0.24.2.tar.gz -o Bismark.tar.gz && \
    tar xvzf Bismark.tar.gz && \
    mv Bismark-* Bismark && \
    chmod -R u+rwx Bismark && \
    chown -R defaultuser:defaultuser Bismark && \
    ln -s /home/methylation/tools/Bismark/bismark /usr/local/bin/bismark

# Copy the pipeline script
COPY run_pipeline.sh /home/methylation/run_pipeline.sh
RUN chmod +x /home/methylation/run_pipeline.sh && chown defaultuser:defaultuser /home/methylation/run_pipeline.sh

# Switch to the non-root user
USER defaultuser

# Set working directory
WORKDIR /home/methylation/

# Set default command
# CMD ["bash", "run_pipeline.sh", "/home/methylation/data/...", "-o", "/home/methylation/output", "-g", "/home/methylation/data/NCBI_genome"]
ENTRYPOINT ["bash", "run_pipeline.sh"]
