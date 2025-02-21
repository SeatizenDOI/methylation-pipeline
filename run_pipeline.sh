#!/bin/bash

# Default directories
OUTPUT_DIR="/data/output"
GENOME_DIR="/data/genome"

# Check if input file is provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 <input_fastq>"
    exit 1
fi

# Check if input directory is provided
if [ $# -lt 1 ]; then
    usage
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
            usage
            ;;
    esac
done

TRIM_DIR="$OUTPUT_DIR/trim_galore"
BISMARK_DIR="$OUTPUT_DIR/bismark"

# Create output directories if they don't exist
mkdir -p "$TRIM_DIR" "$BISMARK_DIR"

if [ ! -f "$INPUT_FILE" ]; then
    echo "No FASTQ files found in $INPUT_DIR. Exiting."
    exit 1
fi

echo "Processing file: $INPUT_FILE"

# Run TrimGalore! and save results in trim_galore folder
trim_galore --fastqc -o "$TRIM_DIR" "$INPUT_FILE"

# Get the trimmed filename (assumes default TrimGalore! naming convention)
TRIMMED_FILE=$(basename "$INPUT_FILE" .fastq.gz)_trimmed.fq.gz

# Check if trimming was successful
if [ ! -f "$TRIM_DIR/$TRIMMED_FILE" ]; then
    echo "TrimGalore! failed for $INPUT_FILE. Skipping."
    continue
fi

# Run Bismark and save results in bismark folder
echo "Running Bismark alignment on $TRIM_DIR/$TRIMMED_FILE using genome from $GENOME_DIR..."
bismark "$GENOME_DIR" -o "$BISMARK_DIR" --fastq "$TRIM_DIR/$TRIMMED_FILE"

echo "Finished processing $INPUT_FILE. TrimGalore! results are in $TRIM_DIR, Bismark results are in $BISMARK_DIR."
