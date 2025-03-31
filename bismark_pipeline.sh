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

GENOME_FASTA=$(find "$GENOME_DIR" -name "*.fa" | head -n 1)
echo "Genome fasta: $GENOME_FASTA"

# Check if genome is prepped; if not, prepare it
if [[ ! -d "$GENOME_DIR/Bisulfite_Genome" ]]; then
    echo "Bisulfite_Genome directory not found in $GENOME_DIR. Running Bismark genome preparation..."
    bismark_genome_preparation "$GENOME_DIR"
    
    if [[ ! -d "$GENOME_DIR/Bisulfite_Genome" ]]; then
        echo "Error: Bismark genome preparation failed!"
        exit 1
    fi
    echo "Genome preparation completed."
fi

# Create output directories if they don't exist
mkdir -p "$OUTPUT_DIR/reports" "$OUTPUT_DIR/trimmed_datasets" "$OUTPUT_DIR/bams" "$OUTPUT_DIR/bedGraphs"

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
    trim_galore -o "$OUTPUT_DIR" "$INPUT_FILE" 
    # Get the trimmed filename
    TRIMMED_FILE=$(basename "$INPUT_FILE" .fastq.gz)_trimmed.fq.gz

    # Check if trimming was successful
    if [ ! -f "$OUTPUT_DIR/$TRIMMED_FILE" ]; then
        echo "TrimGalore! failed for $INPUT_FILE. Skipping."
        exit 1
    fi

    # Run Bismark
    echo "Running Bismark alignment on $OUTPUT_DIR/$TRIMMED_FILE using genome from $GENOME_DIR..."
    bismark "$GENOME_DIR" -o "$OUTPUT_DIR" --temp_dir "$OUTPUT_DIR" --fastq "$OUTPUT_DIR/$TRIMMED_FILE"

    # Extract the unique identifier from the input file
    BASE_NAME=$(basename "$INPUT_FILE" | sed -E 's/(_QCfiltered)?\.fastq\.gz//')
    
    # Find the corresponding BAM file
    BAM_FILE=$(find "$OUTPUT_DIR" -name "${BASE_NAME}_QCfiltered_trimmed_bismark_bt2.bam" | head -n 1)

    if [[ -z "$BAM_FILE" ]]; then
        echo "Error: No BAM file found after Bismark alignment!"
        exit 1
    fi
fi

# Run Bismark Methylation Extractor
bismark_methylation_extractor -o "$OUTPUT_DIR" --bedGraph "$BAM_FILE" 

# Convert bismark output to CGmap format
cgmaptools convert bam2cgmap -b "$BAM_FILE" -g "$GENOME_FASTA" -o "$OUTPUT_DIR"

echo "Finished processing. Results are in $OUTPUT_DIR."
