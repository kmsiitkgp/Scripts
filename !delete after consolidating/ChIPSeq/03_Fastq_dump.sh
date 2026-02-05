#!/bin/bash -l

# Set Project Name
# Set Job Name
#$ -N ChIPSeq_get_Fastq
# Request memory. You can find maxvmem from qstat -j <jobid> & adjust mem_free accordingly in final run.
#$ -l mem_free=8G
# Merge standard output and standard error
#$ -j y
# Set current working directory
#$ -cwd
# Submit an array job with n tasks. Find n using: ls *.fastq.gz | wc -l 
#$ -t 1-12

# You need to add path for fastq-dump to run
PATH=$HOME/NGSTools/sratoolkit.2.11.3-centos_linux64/bin:$PATH

# Read a file line by line and store each line as an array
readarray -t input < ~/projects/ChIPSeq/SRR_Acc_List.txt

# Verify if the array stores correctly
# echo ${input[0]}

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

# Using \ enables us to split the long command into easy to read format.
fastq-dump \
--split-3 \
--origfmt \
--clip \
--readids \
--dumpbase \
--skip-technical \
--read-filter pass \
--gzip \
--outdir ~/scratch/ChIPSeq/raw_reads/ $taskinput
