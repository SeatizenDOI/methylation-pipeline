# Use a lightweight Python image
FROM python:3.10-slim

# Set environment variables to avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    tar \
    unzip \
    perl \
    bowtie2 \
    samtools \
    default-jre \
    && apt-get clean

# Create a non-root user (e.g., 'bioinfo')
RUN useradd -m -s /bin/bash bioinfo && \
    mkdir -p /home/bioinfo/tools /data && \
    chown -R bioinfo:bioinfo /home/bioinfo /data

# Set working directory
WORKDIR /home/bioinfo/tools

# Install TrimGalore!
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -o trim_galore.tar.gz && \
    tar xvzf trim_galore.tar.gz && \
    mv TrimGalore-* TrimGalore && \
    chmod -R u+rwx TrimGalore && \
    chown -R bioinfo:bioinfo TrimGalore && \
    ln -s /home/bioinfo/tools/TrimGalore/trim_galore /usr/local/bin/trim_galore

# Install Bismark
RUN curl -fsSL https://github.com/FelixKrueger/Bismark/archive/refs/tags/v0.24.2.tar.gz -o Bismark.tar.gz && \
    tar xvzf Bismark.tar.gz && \
    mv Bismark-* Bismark && \
    chmod -R u+rwx Bismark && \
    chown -R bioinfo:bioinfo Bismark && \
    ln -s /home/bioinfo/tools/Bismark/bismark /usr/local/bin/bismark

# Install Python dependencies (if any)
COPY requirements.txt /home/bioinfo/
RUN pip install --no-cache-dir -r /home/bioinfo/requirements.txt

# Switch to the non-root user
USER bioinfo

# Set default command
CMD ["bash"]
