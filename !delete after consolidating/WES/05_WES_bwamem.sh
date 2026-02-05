#!/bin/bash -l

#$ -N bwa-align                 # Set Job Name
#$ -l mem_free=40G              # Request memory
#$ -j y                         # Merge standard output and standard error
#$ -cwd                         # Set current working directory
#$ -t 1-3                       # Submit an array job with n tasks; n=number of samples = number of fastq files/2 if paired end

proj="WES_Hany"
species="Mouse"
INPUT_DIR=$HOME/scratch/$proj/cutadapt_results
OUTPUT_DIR=$HOME/scratch/$proj/alignment_results  
REF_FILE=$HOME/NGSTools/Reference_Genomes/$species/*.fa

# Create an array of input files
input1=($INPUT_DIR/*_1*.gz)
input2=($INPUT_DIR/*_2*.gz)

# Use the SGE_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value
index=$(($SGE_TASK_ID-1))
taskinput1=${input1[$index]}
taskinput2=${input2[$index]}
taskoutput=${taskinput1//$INPUT_DIR/$OUTPUT_DIR}
taskoutput=${taskoutput//_1.*/}
sample=$(basename $taskoutput)
taskoutput_sam=$taskoutput.sam

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date    : $(date)"
echo "Job name      : $JOB_NAME"
echo "Job ID        : $JOB_ID"  
echo "SGE TASK ID   : $SGE_TASK_ID"
echo "Index         : $index"
echo "Task input 1  : $taskinput1"
echo "Task input 2  : $taskinput2"
echo "Task output   : $taskoutput"
echo "Project       : $proj"
echo "Sample        : $sample"
echo "Output Folder : $OUTPUT_DIR"
echo "=========================================================="

# You need to add path for Conda
PATH=$HOME/miniconda3/bin/:$PATH
source activate MAGECK

# Read group tags:
# ID : Read group identifier
# SM : The name of the biological sample the reads originated from
# LB : Identifier for the library preparation used to generate the reads
# PL : The sequencing technology used (e.g., Illumina, Ion Torrent)
# PU : Unique identifier for the specific sequencing machine used
# CN : Name of the sequencing facility
# DT : Date of the sequencing run
# PG : Software used for alignment or processing

bwa mem \
-M \
-R "@RG\tID:$proj\tSM:$sample\tPL:ILLUMINA\tPU:$sample" \
$REF_FILE \
$taskinput1 \
$taskinput2 \
-o $taskoutput_sam

# -M: mark shorter split hits as secondary (for Picard compatibility)
# The BWA-MEM algorithm performs local alignment. It may produce multiple primary 
# alignments for different part of a query sequence. This is a crucial feature for
# long sequences. However, some tools such as Picard’s markDuplicates does not work
# with split alignments. One may consider to use option -M to flag shorter split 
# hits as secondary.
# -R: add read group info for GATK compatibility. GATK Mutect2 module uses the SM 
# field value wihtin the read group to match sample names