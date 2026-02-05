#!/bin/bash -l

# Set Project Name
# Set Job Name
#$ -N ChIPSeq_Peak_Calling
# Request memory
# You can find maxvmem from qstat -j <jobid> & adjust mem_free accordingly in final run. It needs ~2G.
#$ -l mem_free=16G
# Merge standard output and standard error
#$ -j y
# Set current working directory
#$ -cwd
# Submit an array job with n tasks. Find n using: ls *.fastq.gz | wc -l divided by 2
#$ -t 1-5

# You need to add path for macs2
PATH=$HOME/Python3/bin/:$PATH

# Follow / at end of variables declareed below to avoid errors in array.
INPUT_DIR=$HOME/scratch/ChIPSeq/Sambamba_results 		#input directory with STAR co-ordinate sorted BAM files
OUTPUT_DIR=$HOME/scratch/ChIPSeq/MACS2_results			#output directory to store MACS2 results

# Create an array of input BAM filenames with path
chip=($INPUT_DIR/ChIP_AR*.bam)
input=$INPUT_DIR/Input_AR*.bam
#chip=($INPUT_DIR/ChIP_IgG*.bam)
#input=$INPUT_DIR/Input_IgG*.bam

# Use the SGE_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value
index=$(($SGE_TASK_ID-1))

# Declare input and output filenames with path for each job
taskchip=${chip[$index]}
sample=$(basename $taskchip)
taskoutput=${sample//.filtered.bam/}

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID  $taskchip"  $input
echo "=========================================================="

macs2 callpeak \
--treatment $taskchip \
--control $input \
--format AUTO \
--gsize hs \
--keep-dup all \
--name $taskoutput \
--outdir $OUTPUT_DIR