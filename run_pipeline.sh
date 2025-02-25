#!/bin/bash

# Default directories
OUTPUT_DIR="/data/output"
GENOME_DIR="/data/genome"

# Check if input file is provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 <input_fastq> [-o <output_dir>] [-g <genome_dir>]"
    exit 1
fi

INPUT_FILE=$1
shift  # Move to next argument

# Parse optional arguments
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
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Define output directories
TRIM_DIR="$OUTPUT_DIR/trim_galore"
BISMARK_DIR="$OUTPUT_DIR/bismark"

# Create necessary directories
mkdir -p "$TRIM_DIR/reports" "$TRIM_DIR/trimmed_datasets" "$BISMARK_DIR/bams" "$BISMARK_DIR/report"

# Check if the input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file $INPUT_FILE not found!"
    exit 1
fi

echo "Processing file: $INPUT_FILE"

# Run TrimGalore! and save results in trim_galore folder
trim_galore --fastqc -o "$TRIM_DIR" "$INPUT_FILE"

# Get the trimmed filename (assumes default TrimGalore! naming convention)
TRIMMED_FILE=$(basename "$INPUT_FILE" .fastq.gz)_trimmed.fq.gz

# Ensure trimming was successful
if [ ! -f "$TRIM_DIR/$TRIMMED_FILE" ]; then
    echo "Error: TrimGalore! failed for $INPUT_FILE."
    exit 1
fi

echo "Running Bismark alignment on $TRIM_DIR/$TRIMMED_FILE using genome from $GENOME_DIR..."
bismark "$GENOME_DIR" -o "$BISMARK_DIR" --fastq "$TRIM_DIR/$TRIMMED_FILE"

# Find the latest BAM file in Bismark output
BAM_FILE=$(find "$BISMARK_DIR" -name "*.bam" | sort -r | head -n 1)

# Ensure Bismark produced a BAM file before proceeding
if [ -z "$BAM_FILE" ]; then
    echo "Error: Bismark failed to produce a BAM file."
    exit 1
fi

echo "Running Bismark methylation extractor on $BAM_FILE..."
bismark_methylation_extractor -s --comprehensive -o "$BISMARK_DIR" "$BAM_FILE"

# Organize TrimGalore! results
mv "$TRIM_DIR"/*.txt "$TRIM_DIR/reports/"
mv "$TRIM_DIR"/*.fq.gz "$TRIM_DIR/trimmed_datasets/"

# Organize Bismark results
mv "$BISMARK_DIR"/*.txt "$BISMARK_DIR/report/"
mv "$BISMARK_DIR"/*.bam "$BISMARK_DIR/bams/"

echo "Finished processing $INPUT_FILE."
echo "TrimGalore! results are in $TRIM_DIR, Bismark results are in $BISMARK_DIR."
