#!/bin/bash -l

# Set Project Name
# Set Job Name
#$ -N ChIPSeq_Genome_Indexing
# Request memory
# You can find maxvmem from qstat -j <jobid> & adjust mem_free accordingly in final run.
#$ -l mem_free=40G
# Merge standard output and standard error
#$ -j y
# Set current working directory
#$ -cwd
# Submit an array job with 1 task. 
#$ -t 1

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID "
echo "=========================================================="

OUTPUT_DIR=$HOME/scratch/ChIPSeq/indexed_ref_genome/		#output directory to store indexed genome. This will be created by STAR
INPUT_DIR=$HOME/scratch/ChIPSeq/ref_genome/			 		#input directory with fa and gtf files of reference genome
MAX_READ_LENGTH_MINUS_1=50									#read length is 51bp. So, we set it to 50
# Calculate max read length using zcat ~/scratch/ChIPSeq/raw_reads/*.fq.gz | head -2 | awk 'END {print}'| wc -c

# You need to add path for STAR
PATH=$HOME/NGSTools/STAR-2.7.10a/bin/Linux_x86_64_static:$PATH

# Generate Genome Indices. You can write all lines below in a single line without \ also. 
# Using \ enables us to split the long command into easy to read format.
STAR \
--runMode genomeGenerate \
--genomeDir $OUTPUT_DIR \
--genomeFastaFiles $INPUT_DIR*.fa \
--sjdbGTFfile $INPUT_DIR*.gtf \
--sjdbOverhang $MAX_READ_LENGTH_MINUS_1