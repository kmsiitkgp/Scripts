#!/bin/bash -l

# RUN THIS ONCE FOR MOUSE AND HUMAN

#$ -N bowtie2-build             # Set Job Name
#$ -l mem_free=96G              # Request memory
#$ -j y                         # Merge standard output and standard error
#$ -cwd                         # Set current working directory
#$ -t 1-1                       # Submit an array job with 1 task.

species="Mouse"
INPUT_DIR=$HOME/NGSTools/Reference_Genomes/$species	                  #input directory with fa and gtf files of reference genome
OUTPUT_DIR=$HOME/NGSTools/Reference_Genomes_Indexed_bowtie2/$species  #output directory to store indexed genome
base_name=$species.reference

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

# To see the sequences alone for first 100 lines
# zcat $HOME/common/$proj/*.fastq.gz | head -100 | awk 'NR % 4 == 2'
# zcat $HOME/common/$proj/*.fastq.gz | grep -n TTTGTTGTCTCGCTCA | head -15

# You need to add path for Conda
PATH=$HOME/miniconda3/bin/:$PATH
source activate MAGECK

cd $OUTPUT_DIR
# Generate Genome Indices. You can write all lines below in a single line without \ also. 
# Using \ enables us to split the long command into easy to read format.
bowtie2-build -f \
$INPUT_DIR/*.fa \
$base_name