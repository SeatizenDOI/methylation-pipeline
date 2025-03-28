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
mkdir -p "$OUTPUT_DIR/tg_reports" "$OUTPUT_DIR/trimmed_datasets" "$OUTPUT_DIR/bams" "$OUTPUT_DIR/bedGraphs"

if [[ -n "$DEBUG_BAM" ]]; then
    # Debug mode: Use provided BAM file
    if [[ ! -f "$DEBUG_BAM" ]]; then
        echo "Error: Debug BAM file $DEBUG_BAM not found!"
        exit 1
    fi
    BAM_FILE="$DEBUG_BAM"
    echo "Debug mode: Using provided BAM file $BAM_FILE"
else
    # Normal mode: Process input FASTQ
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

    # Run BSBolt
    echo "Running BSBolt alignment on $OUTPUT_DIR/$TRIMMED_FILE using genome from $GENOME_DIR..."
    python -m bsbolt Align -F1 "$OUTPUT_DIR/$TRIMMED_FILE" -DB "$GENOME_DIR" -O "$BASE_NAME"

    

    # Find the corresponding BAM file
    BAM_FILE=$(find "$OUTPUT_DIR" -name "${BASE_NAME}.bam" | head -n 1)
    echo "BAM file: $BAM_FILE"

    if [[ -z "$BAM_FILE" ]]; then
        echo "Error: No BAM file found after Bismark alignment!"
        exit 1
    fi
fi

# fixmates to prepare for duplicate removal, use -p to disable proper pair check
samtools fixmate -p -m "$BAM_FILE" "${BASE_NAME}.fixmates.bam"
# sort bam by coordinates for duplicate calling
samtools sort -@ 2 -o "${BASE_NAME}.sorted.bam" "${BASE_NAME}.fixmates.bam"
# remove duplicate reads
samtools markdup "${BASE_NAME}.sorted.bam" "${BASE_NAME}.dup.bam"
# index bam file for methylation calling
samtools index "${BASE_NAME}.dup.bam"

# Run Bismark Methylation Extractor
python -m bsbolt CallMethylation -I "${BASE_NAME}.sorted.bam" -O ${BASE_NAME} -DB ${GENOME_DIR} -t 2 -verbose > "${BASE_NAME}_stats.txt"

echo "Finished processing. Results are in $OUTPUT_DIR."
