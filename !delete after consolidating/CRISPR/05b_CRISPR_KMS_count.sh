#!/bin/bash -l

#$ -N CRISPR_Map_Count          # Set Job Name
#$ -l mem_free=128G             # Request memory
#$ -j y                         # Merge standard output and standard error
#$ -cwd                         # Set current working directory
#$ -t 1-110                       # Submit an array job with 1 task for each forward read fastq.gz files 

# Compile CRISPR scripts by copying and running the below 2 lines in shell. 
# module load gcc/11.1.0
# g++ $HOME/projects/CRISPR/05b_CRISPR_KMS_count_stable/*.cpp -std=c++2b -O3 -o $HOME/projects/CRISPR/CRISPR
# g++ -I $HOME/NGSTools/boost_1_85_0/ $HOME/projects/CRISPR/05b_CRISPR_KMS_count_dev/*.cpp -std=c++2b -O3 -o $HOME/projects/CRISPR/CRISPR_new
# -O1: Enables basic optimization. This includes optimizations such as common subexpression elimination and instruction scheduling. It's a good balance between optimization and compilation time.
# -O2: Enables more aggressive optimization, including inlining functions, loop optimizations, and better code scheduling. It provides a significant performance boost.
# -O3: Enables even more aggressive optimizations. It can lead to faster code but may increase compilation time and the size of the executable.
# Once you have compiled, you will see an executable file "CRISPR".
# Next, qsub this script.

# The guides.csv file MUST have gene name in 1st column and guide sequence in 2nd column

# Follow the backslash / at end of variables declared below to avoid errors in array.
proj="CRISPR_Lena"
libraries=(CDH12KO CDH12ACT ImmuneKO ImmuneACT)

for library in ${libraries[*]}
do

CORRECTED_GUIDE=$HOME/projects/CRISPR/$proj.$library.corrected.kms.csv
FASTQ_DIR=/common/theodorescudlab/$proj/							# directory with fastq files
#FASTQ_DIR=$HOME/scratch/$proj/cutadapt_results/					
OUTPUT_DIR=$HOME/scratch/$proj/count_results/$library/			# directory to store count results

# Create an array of fq_files filenames with path before decompression
fq_gz_files=($FASTQ_DIR*q.gz)

# Create an array of fq_files filenames with path after decompression 
fq_files=(${fq_gz_files[*]//.gz/})

# Use the SGE_TASK_ID environment variable to select the appropriate fq_files file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value
index=$(($SGE_TASK_ID-1))

# Declare fq_gz_files and fq_files with path for gunzip
fq_gz_file=${fq_gz_files[$index]}
fq_file=${fq_files[$index]}
#gunzip -c $fq_gz_file > $fq_file

# for i in ${!fq_gz_files[*]}
# do
# gunzip -c ${fq_gz_files[i]} > ${fq_files[i]}
# done

# Declare the variables
filenames=(${fq_files[*]//$FASTQ_DIR/})
#filenames=(${filenames[*]//.fq/.csv})  # for old code
filenames=(${filenames[*]//.fq/})
filenames=(${filenames[*]//.fastq/})
filenames=(${filenames[*]%%_*})
filenames=(${filenames[*]/#/$OUTPUT_DIR})
filename=${filenames[$index]}
#filename=/home/kailasamms/scratch/CRISPR_Lena/count_results/CDH12KO/PCR-01DO

head_trim=0
tail_trim=25
max_errors=0  #3 for old code
troubleshoot_mode="F"
#orientation="F"   # for old code. new code can detect orientation
# NOTE: We already have the reverse complement in 2nd column which matches with read. So, we use "F".

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date    : $(date)"
echo "Job name      : $JOB_NAME"
echo "Job ID        : $JOB_ID"  
echo "SGE TASK ID   : $SGE_TASK_ID"
echo "Project       : $proj"
echo "library		: $library"
echo "Output Folder : $OUTPUT_DIR"
echo "Fasta file	: $fq_file"
echo "Filename		: $filename"
echo "Guide file	: $CORRECTED_GUIDE_FILE"
echo "=========================================================="

$HOME/projects/CRISPR/CRISPR_new $fq_file $filename $CORRECTED_GUIDE $head_trim $tail_trim $max_errors $troubleshoot_mode
###$HOME/projects/CRISPR/CRISPR $taskfq_files $taskoutput $CORRECTED_GUIDE $head_trim $tail_trim $max_errors $orientation
# After counts for all samples have been calculated, merge into one txt file manually.

done