# Use a lightweight Debian-based image
FROM ubuntu:24.04

# Set environment variables to avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Copy the run_pipeline.sh script to the container
COPY run_pipeline.sh /usr/local/bin/run_pipeline.sh
RUN chmod +x /usr/local/bin/run_pipeline.sh

# Install dependencies
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    perl \
    python3 \
    bowtie2 \
    samtools \
    default-jre \
    && apt-get clean

# Install TrimGalore!
RUN wget https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.10.zip -O TrimGalore.zip && \
    unzip TrimGalore.zip && \
    mv TrimGalore-* /opt/TrimGalore && \
    chmod +x /opt/TrimGalore/trim_galore && \
    ln -s /opt/TrimGalore/trim_galore /usr/local/bin/trim_galore

# Install Bismark
RUN wget https://github.com/FelixKrueger/Bismark/archive/refs/tags/0.24.2.zip -O Bismark.zip && \
    unzip Bismark.zip && \
    mv Bismark-* /opt/Bismark && \
    chmod +x /opt/Bismark/bismark && \
    ln -s /opt/Bismark/bismark /usr/local/bin/bismark

# Set working directory
WORKDIR /data

# Command to keep container alive (modify as needed)
CMD ["bash"]