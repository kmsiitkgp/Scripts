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

# Follow the backslash / at end of variables declared below to avoid errors in array.
INPUT_DIR=$HOME/scratch/Neeraj/raw_reads/

# Create an array of input read1 and read2 filenames with path before decompression
preinput1=($INPUT_DIR*_1.fq.gz)
preinput2=($INPUT_DIR*_2.fq.gz)

# Create an array of input read1 and read2 filenames with path after decompression 
input1=(${preinput1[*]//fq.gz/fq})
input2=(${preinput2[*]//fq.gz/fq})

# Use the SGE_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value
index=$(($SGE_TASK_ID-1))

# Declare input and output filenames with path for gunzip for each job
taskpreinput1=${preinput1[$index]}
taskpreinput2=${preinput2[$index]}
taskinput1=${input1[$index]}
taskinput2=${input2[$index]}

gunzip -c $taskpreinput1 > $taskinput1
gunzip -c $taskpreinput2 > $taskinput2

# Declare the adapter sequence for read1 and read2, minimum read length after trimming
# Declare minimum quality score for trimming and base
sequence1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
sequence2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
type=PE
length=36
score=20
base=33
OUTPUT_DIR=$HOME/scratch/Neeraj/trimmed_reads         #DO NOT END with /. Script will add /

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID  $taskinput1  $taskinput2"
echo "=========================================================="

$HOME/scratch/Neeraj/scripts/FAT $type $length $score $base $OUTPUT_DIR $taskinput1 $sequence1 $taskinput2 $sequence2

rm $taskinput1
rm $taskinput2

