# Methylation data analysis via Bisulfite Sequencing
TrimGalore and Bismark
BSBolt

## TrimGalore + Bismark pipeline

Run the bismark_pipeline.sh script to run the pipeline.

## BSBolt pipeline

Run the bsbolt_pipeline.sh script to run the pipeline.

## Optional conversion of Bismark output to CGmap

For some analysis tools, CGmap is required. The CGmap conversion script consists of two parts:

1. The first part is the dup_bam_script.sh script, which is used to duplicate and sort the Bismark BAM files with samtools.
2. The second part is the bam2cgmap_script.sh script, which is used to convert the Bismark sorted BAM files to CGmap format.