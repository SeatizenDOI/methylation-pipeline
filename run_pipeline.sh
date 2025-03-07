#!/bin/bash

# Default directories
OUTPUT_DIR="/home/methylation/output"
GENOME_DIR="/home/methylation/data/genome"

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

# Create output directories if they don't exist
mkdir -p "$OUTPUT_DIR/tg_reports" "$OUTPUT_DIR/bismark_reports" "$OUTPUT_DIR/bismark_met_reports" "$OUTPUT_DIR/trimmed_datasets" "$OUTPUT_DIR/bams"

if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file $INPUT_FILE not found!"
    exit 1
fi

echo "Processing file: $INPUT_FILE"

# Run TrimGalore!
trim_galore -o "$OUTPUT_DIR" "$INPUT_FILE"

# Get the trimmed filename
TRIMMED_FILE=$(basename "$INPUT_FILE" .fastq.gz)_trimmed.fq.gz

# Check if trimming was successful
if [ ! -f "$OUTPUT_DIR/$TRIMMED_FILE" ]; then
    echo "TrimGalore! failed for $INPUT_FILE. Skipping."
    exit 1
fi

mv "$OUTPUT_DIR"/*.txt "$OUTPUT_DIR/tg_reports/"

# Run Bismark
echo "Running Bismark alignment on $OUTPUT_DIR/$TRIMMED_FILE using genome from $GENOME_DIR..."
bismark "$GENOME_DIR" -o "$OUTPUT_DIR" --temp_dir "$OUTPUT_DIR" --fastq "$OUTPUT_DIR/$TRIMMED_FILE"

mv "$OUTPUT_DIR"/*.txt "$OUTPUT_DIR/bismark_reports/"

# Extract the unique identifier from the input file
BASE_NAME=$(basename "$INPUT_FILE" | sed -E 's/(_QCfiltered)?\.fastq\.gz//')
echo "Base name: $BASE_NAME"

# Find the corresponding BAM file
BAM_FILE=$(find "$OUTPUT_DIR" -maxdepth 1 -name "${BASE_NAME}_trimmed_bismark_bt2.bam" | head -n 1)
echo "BAM file: $BAM_FILE"

if [[ -z "$BAM_FILE" ]]; then
    echo "Error: No BAM file found after Bismark alignment!"
    exit 1
fi

# Run Bismark Methylation Extractor
# bismark_methylation_extractor -s --comprehensive "$BAM_FILE"

# Organize results
mv "$OUTPUT_DIR"/*.fq.gz "$OUTPUT_DIR/trimmed_datasets/"
mv "$OUTPUT_DIR"/*.bam "$OUTPUT_DIR/bams/"
mv "$OUTPUT_DIR"/*.txt "$OUTPUT_DIR/bismark_met_reports/"

echo "Finished processing $INPUT_FILE. Results are in $OUTPUT_DIR."
