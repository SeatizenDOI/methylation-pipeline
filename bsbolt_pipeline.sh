#!/bin/bash

# Default directories
OUTPUT_DIR="/home/methylation/output"
GENOME_DIR="/home/methylation/genomes"
DEBUG_BAM=""

# Define the usage function
usage() {
    echo "Usage: $0 <input_fastq> [-o <output_dir>] [-g <genome_dir>] [--debug-bam <bam_file>]"
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
        --debug-bam)
            DEBUG_BAM="$2"
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

# Create output directories if they don't exist
mkdir -p "$OUTPUT_DIR/reports" "$OUTPUT_DIR/trimmed_datasets" "$OUTPUT_DIR/bams"

if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file $INPUT_FILE not found!"
    exit 1
fi

echo "Processing file: $INPUT_FILE"

# Run TrimGalore!
trim_galore --rrbs -o "$OUTPUT_DIR" "$INPUT_FILE" 
# Get the trimmed filename
TRIMMED_FILE=$(basename "$INPUT_FILE" .fastq.gz)_trimmed.fq.gz

# Check if trimming was successful
if [ ! -f "$OUTPUT_DIR/$TRIMMED_FILE" ]; then
    echo "TrimGalore! failed for $INPUT_FILE. Skipping."
    exit 1
fi

# Extract the unique identifier from the input file
BASE_NAME=$(basename "$INPUT_FILE" | sed -E 's/(_QCfiltered)?\.fastq\.gz//')
echo "Base name: $BASE_NAME"

# Find fasta file in genome directory
GENOME_FASTA=$(find "$GENOME_DIR" -name "*.fa" | head -n 1)
echo "Genome fasta: $GENOME_FASTA"

# Name DB for indexing
BSBOLT_DB="$GENOME_DIR/BSBOLT_DB"
echo "DB name: $BSBOLT_DB"

# Run BSBolt
echo "Running BSBolt alignment on $OUTPUT_DIR/$TRIMMED_FILE using genome from $GENOME_DIR..."
python -m bsbolt Align -F1 "$OUTPUT_DIR/$TRIMMED_FILE" -DB $BSBOLT_DB -O "$OUTPUT_DIR/$BASE_NAME"

# fixmates to prepare for duplicate removal, use -p to disable proper pair check
samtools fixmate -p -m "$OUTPUT_DIR/${BASE_NAME}.bam" "$OUTPUT_DIR/${BASE_NAME}.fixmates.bam"
# sort bam by coordinates for duplicate calling
samtools sort -@ 2 -o "$OUTPUT_DIR/${BASE_NAME}.sorted.bam" "$OUTPUT_DIR/${BASE_NAME}.fixmates.bam"
# remove duplicate reads
samtools markdup "$OUTPUT_DIR/${BASE_NAME}.sorted.bam" "$OUTPUT_DIR/${BASE_NAME}.dup.bam"
# index bam file for methylation calling
samtools index "$OUTPUT_DIR/${BASE_NAME}.dup.bam"

# Run Bismark Methylation Extractor
python -m bsbolt CallMethylation -I "$OUTPUT_DIR/${BASE_NAME}.sorted.bam" -O "$OUTPUT_DIR/${BASE_NAME}" -DB ${GENOME_DIR} -t 2 -verbose > "$OUTPUT_DIR/${BASE_NAME}_stats.txt"

echo "Finished processing. Results are in $OUTPUT_DIR."