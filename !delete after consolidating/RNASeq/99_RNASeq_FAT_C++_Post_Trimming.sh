#!/bin/bash

# Set Project Name
# Set Job Name
#$ -N RNA_Seq_FAT_Trimming
# Request memory. You can find maxvmem from qstat -j <jobid> & adjust mem_free accordingly in final run.
#$ -l mem_free=1G
# Merge standard output and standard error
#$ -j y
# Set current working directory
#$ -cwd
# Submit an array job with n tasks. Find n using: ls *.fastq.gz | wc -l /2
#$ -t 1-5

# Use the SGE_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value
index=$(($SGE_TASK_ID-1))
OUTPUT_DIR=$HOME/scratch/Neeraj/trimmed_reads         #DO NOT END with /. Script will add /

# Create an array of read1 and read2 filenames with path generated after trimming for compression
output1=($OUTPUT_DIR/trimmed*_1.fq)
output2=($OUTPUT_DIR/UP*_1.fq)
output3=($OUTPUT_DIR/trimmed*_2.fq)
output4=($OUTPUT_DIR/UP*_2.fq)
taskoutput1=${output1[$index]}
taskoutput2=${output2[$index]}
taskoutput3=${output3[$index]}
taskoutput4=${output4[$index]}

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID  $taskoutput1  $taskoutput2  $taskoutput3  $taskoutput4"
echo "=========================================================="

gzip $taskoutput1
gzip $taskoutput2
gzip $taskoutput3
gzip $taskoutput4
