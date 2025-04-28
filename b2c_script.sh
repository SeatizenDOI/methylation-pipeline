BAM_DIR="/home/bdl/Documents/project-hub/STAGE/methylation-pipeline/genomes/bams"
OUTPUT_DIR="/home/bdl/Documents/project-hub/STAGE/methylation-pipeline/genomes/cgmapstest"
GENOME="/home/bdl/Documents/project-hub/STAGE/methylation-pipeline/genomes/NCBI_Albacares/GCF_914725855.1_fThuAlb1.1_genomic.fa"

for file in $BAM_DIR/*.dup.bam
do
    BASE_NAME=$(basename "$file" | sed -n 's/.dup.bam//p')
    echo "Base name: $BASE_NAME"
    cgmaptools convert bam2cgmap -b $file -g $GENOME -o "$OUTPUT_DIR/${BASE_NAME}" 
done