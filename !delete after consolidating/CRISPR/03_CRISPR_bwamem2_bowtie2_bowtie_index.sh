#!/bin/bash -l

#$ -N CRISPR_index              # Set Job Name
#$ -l mem_free=96G              # Request memory
#$ -j y                         # Merge standard output and standard error
#$ -cwd                         # Set current working directory
#$ -t 1-1                       # Submit an array job with 1 task.

proj="CRISPR_Lena"
libraries=(CDH12KO CDH12ACT ImmuneKO ImmuneACT)

proj="CRISPR_Prince"
libraries=(DTKPA1 DTKPA2)
INPUT_DIR=$HOME/projects/CRISPR                            # directory containing guide.csv
REF_DIR_BWAMEM2=$HOME/scratch/$proj/index_results/bwamem2  # directory to store indexed genome
REF_DIR_BOWTIE2=$HOME/scratch/$proj/index_results/bowtie2  # directory to store indexed genome
REF_DIR_BOWTIE=$HOME/scratch/$proj/index_results/bowtie    # directory to store indexed genome

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date    : $(date)"
echo "Job name      : $JOB_NAME"
echo "Job ID        : $JOB_ID"  
echo "SGE TASK ID   : $SGE_TASK_ID"
echo "Project       : $proj"
echo "Input Folder  : $INPUT_DIR"
echo "Output Folder : $REF_DIR_BWAMEM2"
echo "Output Folder : $REF_DIR_BOWTIE2"
echo "Output Folder : $REF_DIR_BOWTIE"
echo "=========================================================="

for library in ${libraries[*]}
do

GUIDE_FILE=$INPUT_DIR/$proj.$library.csv
CORRECTED_GUIDE_FILE=$INPUT_DIR/$proj.$library.corrected.csv
#CORRECTED_GUIDE_FILE=$INPUT_DIR/$proj.$library.corrected.kms.csv

# STOP: MOST IMPORTANT ISSUE
# Your csv file containing guides looks fine if you open in excel or notepad.
# Try cat in linux to view your csv file. You will notice every line is quoted.
# The reason for the quotation marks is that your values contain commas 
# Eg: sgRNA_safe_739,GTATAGATGCTGCATTA,safe and therefore, the CSV writer must 
# put them in quotes to distinguish those commas from the record separators.
# So, first step is to remove these quotes.
tr -d '"' <$GUIDE_FILE >$CORRECTED_GUIDE_FILE

# Create fa file for indexing
awk -F ',' '{print ">"$1"\n"$2}' $CORRECTED_GUIDE_FILE > $INPUT_DIR/$proj.$library.fa

# You need to add path for Python3 for MAGeCK to run. Not needed if installed using Conda
# PATH=$HOME/Python3/bin/:$PATH

# You need to add path for Conda
PATH=$HOME/miniconda3/bin/:$PATH
source activate MAGECK

# Generate Genome Indices. You can write all lines below in a single line without \ also. 
# Using \ enables us to split the long command into easy to read format.
printf "\n\n***********Creating BWAMEM2 Index******************\n\n"
cd $REF_DIR_BWAMEM2
bwa-mem2 index \
-p $proj.$library \
$INPUT_DIR/$proj.$library.fa

printf "\n\n***********Creating BOWTIE2 Index******************\n\n"
cd $REF_DIR_BOWTIE2
bowtie2-build -f \
$INPUT_DIR/$proj.$library.fa \
$proj.$library

printf "\n\n***********Creating BOWTIE Index******************\n\n"
cd $REF_DIR_BOWTIE
bowtie-build -f \
$INPUT_DIR/$proj.$library.fa \
$proj.$library

done