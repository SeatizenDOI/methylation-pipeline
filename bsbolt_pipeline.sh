#!/bin/bash

# Default directories
OUTPUT_DIR="/home/methylation/output"
DB_DIR="/home/methylation/DB"
DEBUG_BAM=""

# Define the usage function
usage() {
    echo "Usage: $0 <input_fastq> [-f2 <input_fastq2>] [-o <output_dir>] [-d <db_dir>]"
    exit 1
}

# Check if input file is provided
if [ $# -lt 1 ]; then
    usage
fi

INPUT_FILE=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        -f2|--fasta2)
            INPUT_FILE2="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -d|--db)
            DB_DIR="$2"
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

# Extract the unique identifier from the input file
BASE_NAME=$(basename "$INPUT_FILE" | sed -E 's/(_QCfiltered)?\.fastq\.gz//')
echo "Base name: $BASE_NAME"

# Run BSBolt
if [ ! -f "$INPUT_FILE2" ]; then
    echo "Running BSBolt alignment on $INPUT_FILE using genome from $GENOME_DIR..."
    python -m bsbolt Align -t 8 -F1 "$INPUT_FILE" -DB $DB_DIR -O "$OUTPUT_DIR/$BASE_NAME" -T 110 >& ${BASE_NAME}.log
else 
    echo "Running BSBolt alignment on $INPUT_FILE using genome from $GENOME_DIR..."
    python -m bsbolt Align -t 8 -F1 "$INPUT_FILE" -F2 "$INPUT_FILE2" -DB $DB_DIR -O "$OUTPUT_DIR/$BASE_NAME" -T 110 >& ${BASE_NAME}.log
fi

# fixmates to prepare for duplicate removal, use -p to disable proper pair check
echo "samtools fixmate"
samtools fixmate -p -m "$OUTPUT_DIR/${BASE_NAME}.bam" "$OUTPUT_DIR/${BASE_NAME}.fixmates.bam"
# sort bam by coordinates for duplicate calling
echo "samtools sort"
samtools sort -@ 2 -o "$OUTPUT_DIR/${BASE_NAME}.sorted.bam" "$OUTPUT_DIR/${BASE_NAME}.fixmates.bam"
# remove duplicate reads
echo "samtools markdup"
samtools markdup "$OUTPUT_DIR/${BASE_NAME}.sorted.bam" "$OUTPUT_DIR/${BASE_NAME}.dup.bam"
# index bam file for methylation calling
echo "samtools index"
samtools index "$OUTPUT_DIR/${BASE_NAME}.dup.bam"

# Run BSBolt Methylation Extractor
echo "Running BSBolt Methylation Extractor on $OUTPUT_DIR/${BASE_NAME}.sorted.bam"
python -m bsbolt CallMethylation -t 8 -BQ 10 -MQ 10 -IO -min 5 -CG -I "$OUTPUT_DIR/${BASE_NAME}.sorted.bam" -O "$OUTPUT_DIR/${BASE_NAME}" -DB $DB_DIR -t 2 -verbose > "$OUTPUT_DIR/${BASE_NAME}_stats.txt"

echo "Finished processing. Results are in $OUTPUT_DIR."
