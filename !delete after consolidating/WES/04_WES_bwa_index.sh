#!/bin/bash -l

# RUN THIS ONCE FOR MOUSE AND HUMAN

#$ -N bwa-index                 # Set Job Name
#$ -l mem_free=96G              # Request memory
#$ -j y                         # Merge standard output and standard error
#$ -cwd                         # Set current working directory
#$ -t 1-1                       # Submit an array job with 1 task.

# NOTE: For BWA and most softwares like Picard, Mutect2 to work, the ref.fa 
# and indexed files MUST BE in same location
species="Mouse"
INPUT_DIR=$HOME/NGSTools/Reference_Genomes/$species			#input directory with fa and gtf files of reference genome
REF_FILE=$INPUT_DIR/*.fa

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date    : $(date)"
echo "Job name      : $JOB_NAME"
echo "Job ID        : $JOB_ID"  
echo "SGE TASK ID   : $SGE_TASK_ID"
echo "Project       : $proj"
echo "Input Folder  : $INPUT_DIR"
echo "Output Folder : $OUTPUT_DIR"
echo "=========================================================="

# You need to add path for Conda
PATH=$HOME/miniconda3/bin/:$PATH
source activate MAGECK

cd $INPUT_DIR
# Generate Genome Indices. You can write all lines below in a single line without \ also. 
# Using \ enables us to split the long command into easy to read format.
bwa index $REF_FILE

# Also, create a fasta index for the reference genome which is required by Mutect2
# NOTE: The fai index file MUST be stored in same directory as the reference fasta file
samtools faidx $REF_FILE

# Also, create a fasta dict for the reference genome which is required by Mutect2
java -jar ~/NGSTools/picard.jar CreateSequenceDictionary \
--REFERENCE $REF_FILE
