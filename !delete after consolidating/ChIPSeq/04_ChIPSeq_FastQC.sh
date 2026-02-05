#!/bin/bash -l

# Set Project Name
# Set Job Name
#$ -N ChIPSeq_FastQC
# Request memory. You can find maxvmem from qstat -j <jobid> & adjust mem_free accordingly in final run.
#$ -l mem_free=2G
# Merge standard output and standard error
#$ -j y
# Set current working directory
#$ -cwd
# Submit an array job with n tasks where n = number of fastq files
#$ -t 1-12

# You need to add path for java as FastQC needs Java to run
PATH=$HOME/NGSTools/jre1.8.0_271/bin:$PATH

input=($HOME/scratch/ChIPSeq/raw_reads/*.fastq*)
OUTPUT_DIR=$HOME/scratch/ChIPSeq/fastqc_results/		#output directory to store fastqc results

# Use the SGE_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value
index=$(($SGE_TASK_ID-1))
taskinput=${input[$index]}

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID  $taskinput"
echo "=========================================================="

# --outdir=~/Hany/fastqc_results/ wont work. ~ isnt recognized by fastqc
$HOME/NGSTools/FastQC/fastqc $taskinput --outdir=$OUTPUT_DIR