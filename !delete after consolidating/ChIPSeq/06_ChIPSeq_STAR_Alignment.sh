#!/bin/bash -l

# Set Project Name
# Set Job Name
#$ -N ChIPSeq_STAR_Align_Reads
# Request memory. You can find maxvmem from qstat -j <jobid> & adjust mem_free accordingly in final run.
#$ -l mem_free=30G
# Merge standard output and standard error
#$ -j y
# Set current working directory
#$ -cwd
# Submit an array job with n tasks. Find n using: ls *.fastq.gz | wc -l divided by 2
#$ -t 1-12

# You need to add path for STAR
# Follow / at end of variables declareed below to avoid errors in array.
PATH=$HOME/NGSTools/STAR-2.7.10a/bin/Linux_x86_64_static:$PATH
GENOME_INDICES=$HOME/scratch/ChIPSeq/indexed_ref_genome/	#directory containing genome indices from previous step
INPUT_DIR=$HOME/scratch/ChIPSeq/raw_reads					#input directory with trimmed reads
OUTPUT_DIR=$HOME/scratch/ChIPSeq/STAR_alignment_results		#output directory to store read alignment results

# Create an array of read1 and read2 filenames with path
input1=($INPUT_DIR/*.fastq)
#input1=($INPUT_DIR/*_1.fq.gz)
#input2=($INPUT_DIR/*_2.fq.gz)

# Create an array of filenames with output path
# These are prefix with output path. STAR will append "SortedByCoordinate.bam" to filename 
output1=(${input1[*]//$INPUT_DIR/$OUTPUT_DIR})

# Use the SGE_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value
index=$(($SGE_TASK_ID-1))

# Declare input and output filenames with path for each job
taskinput1=${input1[$index]}
taskoutput=${output1[$index]}

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID  $taskinput1  $taskinput2"
echo "=========================================================="

# Generate Genome Indices. You can write all lines below in a single line without \ also. 
# Using \ enables us to split the long command into easy to read format.
# It is key to set alignIntronMax to 1 to align reads to contiguos sequence in the reference genome. 
# Else, reads will be considered to made of 2 or more exons and split between 2 exons.
STAR \
--runMode alignReads \
--genomeDir $GENOME_INDICES \
--readFilesIn $taskinput1 \
--outFileNamePrefix $taskoutput \
--outFilterMismatchNoverReadLmax 0.05 \
--alignIntronMax 1 \
--outSAMtype BAM SortedByCoordinate