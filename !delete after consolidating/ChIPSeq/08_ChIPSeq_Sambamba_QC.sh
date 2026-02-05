#!/bin/bash -l

# Set Project Name
# Set Job Name
#$ -N ChIPSeq_Sambamba_Filtering_QC
# Request memory
# You can find maxvmem from qstat -j <jobid> & adjust mem_free accordingly in final run. It needs ~2G.
#$ -l mem_free=8G
# Merge standard output and standard error
#$ -j y
# Set current working directory
#$ -cwd
# Submit an array job with n tasks. Find n using: ls *.fastq.gz | wc -l divided by 2
#$ -t 1-12

# Follow / at end of variables declareed below to avoid errors in array.
INPUT_DIR=$HOME/scratch/ChIPSeq/STAR_alignment_results		#input directory with STAR co-ordinate sorted BAM files
OUTPUT_DIR=$HOME/scratch/ChIPSeq/Sambamba_results			#output directory to store uniquely mapped reads

# Create an array of input BAM filenames with path
input=($INPUT_DIR/*.bam)

# Create an array of input BAM filenames with path
output=(${input[*]//$INPUT_DIR/$OUTPUT_DIR})
output1=(${output[*]//fastqAligned.sortedByCoord.out/filtered})

# Use the SGE_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value
index=$(($SGE_TASK_ID-1))

# Declare input and output filenames with path for each job
taskinput=${input[$index]}
taskoutput=${output[$index]}

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID  $taskinput  $taskoutput"
echo "=========================================================="

$HOME/NGSTools/sambamba-0.8.2-linux-amd64-static flagstat $taskoutput 
