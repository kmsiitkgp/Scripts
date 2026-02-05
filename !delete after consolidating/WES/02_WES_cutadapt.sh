#!/bin/bash -l

#$ -N cutadapt                  # Set Job Name
#$ -l mem_free=8G               # Request memory
#$ -j y                         # Merge standard output and standard error
#$ -cwd                         # Set current working directory
#$ -t 1-2                     # Submit an array job with n tasks; n=number of samples = number of fastq files/2 if paired end

proj="RNASeq_GSE75192"
OUTPUT_DIR=$HOME/scratch/$proj/fastqc_results/   #output directory to store fastqc results

# Create an array of input files
input=($HOME/common/$proj/*.gz)

# Use the SGE_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value
index=$(($SGE_TASK_ID-1))
taskinput=${input[$index]}

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date    : $(date)"
echo "Job name      : $JOB_NAME"
echo "Job ID        : $JOB_ID"  
echo "SGE TASK ID   : $SGE_TASK_ID"
echo "Index         : $index"
echo "Task input    : $taskinput"
echo "Project       : $proj"
echo "Output Folder : $OUTPUT_DIR"
echo "=========================================================="

# You need to add path for Conda
PATH=$HOME/miniconda3/bin/:$PATH
source activate MAGECK

cutadapt \
-m 25 \
-a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
-A ADAPT2 \
-g GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG \
-G \
-o $OUTPUT1 \
-p $OUTPUT2 \
$INPUT1 $INPUT2