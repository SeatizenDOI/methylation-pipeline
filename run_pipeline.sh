#!/bin/bash

# Default directories
OUTPUT_DIR="/data/output"
GENOME_DIR="/data/genome"

# Define the usage function
usage() {
    echo "Usage: $0 <input_fastq> [-o <output_dir>] [-g <genome_dir>]"
    exit 1
}

# Check if input file is provided
if [ $# -lt 1 ]; then
    usage
fi

INPUT_FILE=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -g|--genome)
            GENOME_DIR="$2"
            shift 2
            ;;
        *)
            if [[ -z "$INPUT_FILE" ]]; then
                INPUT_FILE="$1"
                shift
            else
                echo "Unknown option: $1"
                usage
            fi
            ;;
    esac
done

if [[ -z "$INPUT_FILE" ]]; then
    echo "Error: No input FASTQ file provided."
    usage
fi

TRIM_DIR="$OUTPUT_DIR/trim_galore"
BISMARK_DIR="$OUTPUT_DIR/bismark"

# Create output directories if they don't exist
mkdir -p "$TRIM_DIR/reports" "$TRIM_DIR/trimmed_datasets" "$BISMARK_DIR/bams" "$BISMARK_DIR/report"

if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file $INPUT_FILE not found!"
    exit 1
fi

echo "Processing file: $INPUT_FILE"

# Run TrimGalore!
trim_galore --fastqc -o "$TRIM_DIR" "$INPUT_FILE"

# Get the trimmed filename
TRIMMED_FILE=$(basename "$INPUT_FILE" .fastq.gz)_trimmed.fq.gz

# Check if trimming was successful
if [ ! -f "$TRIM_DIR/$TRIMMED_FILE" ]; then
    echo "TrimGalore! failed for $INPUT_FILE. Skipping."
    exit 1
fi

# Run Bismark
echo "Running Bismark alignment on $TRIM_DIR/$TRIMMED_FILE using genome from $GENOME_DIR..."
bismark "$GENOME_DIR" -o "$BISMARK_DIR" --fastq "$TRIM_DIR/$TRIMMED_FILE"

# Find the BAM file produced by Bismark
BAM_FILE=$(find "$BISMARK_DIR" -maxdepth 1 -name "*.bam" | head -n 1)

if [[ -z "$BAM_FILE" ]]; then
    echo "Error: No BAM file found after Bismark alignment!"
    exit 1
fi

# Run Bismark Methylation Extractor
bismark_methylation_extractor -s --comprehensive "$BAM_FILE"

# Organize results
mv "$TRIM_DIR"/*.txt "$TRIM_DIR/reports/"
mv "$TRIM_DIR"/*.fq.gz "$TRIM_DIR/trimmed_datasets/"
mv "$BISMARK_DIR"/*.txt "$BISMARK_DIR/report/"
mv "$BISMARK_DIR"/*.bam "$BISMARK_DIR/bams/"

echo "Finished processing $INPUT_FILE. Results are in $OUTPUT_DIR."
