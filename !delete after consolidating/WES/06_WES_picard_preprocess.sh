#!/bin/bash -l

#$ -N picard                    # Set Job Name
#$ -l mem_free=40G              # Request memory
#$ -j y                         # Merge standard output and standard error
#$ -cwd                         # Set current working directory
#$ -t 1-3                       # Submit an array job with n tasks; n=number of samples = number of fastq files/2 if paired end

species="Mouse"
proj="WES_Hany"
INPUT_DIR=$HOME/scratch/$proj/alignment_results  
OUTPUT_DIR=$INPUT_DIR
REF_FILE=$HOME/NGSTools/Reference_Genomes/$species/*.fa

# Create an array of input files
input=($INPUT_DIR/*.sam)

# Use the SGE_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value
index=$(($SGE_TASK_ID-1))
taskinput=${input[$index]}
sample=${taskinput//$INPUT_DIR/$OUTPUT_DIR}
sample=${sample//.sam/}

unsorted_sam=$sample.sam
query_sorted_bam=$sample.query.sorted.bam
marked_dup_bam=$sample.marked.dup.bam
marked_dup_metrics=$sample.marked.dup.metrics.txt
coord_sorted_bam=$sample.coord.sorted.bam
align_metrics=$sample.align.metrics.txt

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date    : $(date)"
echo "Job name      : $JOB_NAME"
echo "Job ID        : $JOB_ID"  
echo "SGE TASK ID   : $SGE_TASK_ID"
echo "Project       : $proj"
echo "Input Folder	: $INPUT_DIR"
echo "Output Folder : $OUTPUT_DIR"
echo "Index         : $index"
echo "Sample		: $sample"
echo "Task output   : $unsorted_sam"
echo "Task output   : $query_sorted_bam"
echo "Task output   : $coord_sorted_bam"
echo "Task output   : $marked_dup_bam"
echo "=========================================================="

# You need to add path for Java
PATH=$HOME/NGSTools/jdk-23.0.1/bin/:$PATH

# Sort reads by queryname
# NOTE: Picard can mark and remove duplicates in either coordinate-sorted or query-sorted BAM/SAM files,
# however, if the alignments are query-sorted it can test secondary alignments for duplicates. A brief
# discussion of this nuance is discussed in the MarkDuplicates manual of Picard.
java -jar ~/NGSTools/picard.jar SortSam \
--INPUT $unsorted_sam \
--OUTPUT $query_sorted_bam \
--SORT_ORDER queryname
	 	 
# Mark and remove duplicates	 
java -jar ~/NGSTools/picard.jar MarkDuplicates \
--INPUT $query_sorted_bam \
--OUTPUT $marked_dup_bam \
--METRICS_FILE $marked_dup_metrics \
--REMOVE_DUPLICATES true
	  
# Sort reads by co-ordinates
java -jar ~/NGSTools/picard.jar SortSam \
--INPUT $marked_dup_bam \
--OUTPUT $coord_sorted_bam \
--SORT_ORDER coordinate \
--CREATE_INDEX true	 

# Calculate alignment metrics
java -jar ~/NGSTools/picard.jar CollectAlignmentSummaryMetrics \
--INPUT $coord_sorted_bam \
--REFERENCE_SEQUENCE $REF_FILE \
--OUTPUT $align_metrics