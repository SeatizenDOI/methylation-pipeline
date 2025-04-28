BAM_DIR="/home/bdl/Documents/project-hub/STAGE/methylation-pipeline/genomes/bams"

for file in $BAM_DIR/*.bam
do
    BASE_NAME=$(basename "$file" | sed -E 's/(_R1_QCfiltered_trimmed_bismark_bt2)?\.bam//')
    echo "Base name: $BASE_NAME"
    # fixmates to prepare for duplicate removal, use -p to disable proper pair check
    echo "samtools fixmate"
    samtools fixmate -p -m $file "$BAM_DIR/${BASE_NAME}.fixmates.bam"
    # sort bam by coordinates for duplicate calling
    echo "samtools sort"
    samtools sort -@ 2 -o "$BAM_DIR/${BASE_NAME}.sorted.bam" "$BAM_DIR/${BASE_NAME}.fixmates.bam"
    # remove duplicate reads
    echo "samtools markdup"
    samtools markdup "$BAM_DIR/${BASE_NAME}.sorted.bam" "$BAM_DIR/${BASE_NAME}.dup.bam"
    # index bam file for methylation calling
    echo "samtools index"
    samtools index "$BAM_DIR/${BASE_NAME}.dup.bam"
done

# cgmaptools convert bam2cgmap -b "$BAM_DIR/${BASE_NAME}.dup.bam" -g $GENOME -o "$OUTPUT_DIR/${BASE_NAME}" 