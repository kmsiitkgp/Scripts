#!/bin/bash -l

# RUN THIS ONCE FOR MOUSE AND HUMAN

#$ -N STAR_Genome_Indexing      # Set Job Name
#$ -l mem_free=64G              # Request memory
#$ -j y                         # Merge standard output and standard error
#$ -cwd                         # Set current working directory
#$ -t 1                         # Submit an array job with 1 task.

proj="RNASeq_Hany_Male_ImmuneEditing"
proj="RNASeq_Hany_CRISPR_LOY"
proj="RNASeq_Sandrine.Supriya"
species="Mouse" 

#proj="RNASeq_Vera"
#proj="RNASeq_Manish"
species="Human"

FASTQ_DIR=/common/theodorescudlab/Sequencing_Data/$proj/01.RawData	# directory with trimmed/raw reads
INPUT_DIR=$HOME/NGSTools/Reference_Genomes/$species					# directory with fa and gtf files of reference genome
OUTPUT_DIR=$HOME/NGSTools/Reference_Genomes_STAR/$species	        # directory to store indexed genome. This will be created by STAR
RSEM_DIR=$HOME/NGSTools/Reference_Genomes_RSEM/$species
RSEQC_DIR=$HOME/NGSTools/Reference_Genomes_RSEQC/$species
GTF=$(basename $(ls $INPUT_DIR/*gtf))
BED=${GTF/gtf/bed}

# Calculate max read length using zcat $HOME/common/$proj/*q.gz  | head -2 | awk 'END {print}'| wc -c
# Even if number of bases is 150, wc -c will give  result as 151. So, we subtract 2 to get 149 from output of wc -c
# If read length is 150bp, set it to 149.

# 1. All the conditional expressions should be placed inside square braces with spaces around them. So, 
# ensure whitespaces between the brackets and the actual check/comparison statement.
# For example, if [$x==0] will not work.
# Bash will report an error about a missing ].
# 2. Always end the line before adding a new keyword, such as “then.”
# If, then, else, elif, and fi are all shell keywords, which means they can’t be used on the same line. 
# Put a “;” between the previous statement and the keyword, or start a new line with the keyword.
# 3. To use many conditions in one statement, use logical operators.
# We can use logical AND(&&) or logical OR(||) operators to use multiple conditions. For example,
# if [[ $x -ge $y ]] && [[ $x -ge $z ]]; then
# echo "x is greatest"
# fi
# 4. You cannot use <= or >= within if condition. 
# Use-ge ("greater-than or equal"), -lt ("less-than"), -le ("less-than or equal"), -eq ("equal") and -ne ("not equal").

MAX_READ_LENGTH_MINUS_1=$(($(zcat $FASTQ_DIR/*q.gz | head -2 | awk 'END {print}'| wc -c)-2))

if [ $MAX_READ_LENGTH_MINUS_1 -le 100 ]
then
    MAX_READ_LENGTH_MINUS_1=100
elif [ $MAX_READ_LENGTH_MINUS_1 -le 150 ]
then
    MAX_READ_LENGTH_MINUS_1=150
else
 echo $MAX_READ_LENGTH_MINUS_1
fi

echo $MAX_READ_LENGTH_MINUS_1

MAX_READ_LENGTH_MINUS_1=150

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date    : $(date)"
echo "Job name      : $JOB_NAME"
echo "Job ID        : $JOB_ID"  
echo "SGE TASK ID   : $SGE_TASK_ID"
echo "Project       : $proj"
echo "Output Folder : $OUTPUT_DIR"
echo "Read length-1 : $MAX_READ_LENGTH_MINUS_1"
echo "=========================================================="

# To see the sequences alone for first 100 lines
# zcat $HOME/common/$proj/*.fastq.gz | head -100 | awk 'NR % 4 == 2'
# zcat $HOME/common/$proj/*.fastq.gz | grep -n TTTGTTGTCTCGCTCA | head -15

# You need to add path for STAR
PATH=$HOME/NGSTools/STAR-2.7.11b/bin/Linux_x86_64_static/:$PATH

# Generate Genome Indices. You can write all lines below in a single line without \ also. 
# Using \ enables us to split the long command into easy to read format.
STAR \
--runMode genomeGenerate \
--genomeDir $OUTPUT_DIR \
--genomeFastaFiles $INPUT_DIR/*.fa \
--sjdbGTFfile $INPUT_DIR/*.gtf \
--sjdbOverhang $MAX_READ_LENGTH_MINUS_1

#--sjdbGTFtagExonParentTranscript Parent # use this ONLY if you use GFF instead of GTF 
#--genomeSAindexNbases 11                # use this ONLY if you using small genomes as reference

# You need to add path for Conda
PATH=$HOME/miniconda3/bin/:$PATH
source activate MAGECK

# Create RSEM reference
cd $RSEM_DIR

rsem-prepare-reference \
--gtf $INPUT_DIR/*.gtf \
$INPUT_DIR/*.fa \
$RSEM_DIR/RSEM_

# Convert GTF to BED for use in RSEQC read_distribution.py module
source activate R

gxf2bed \
--input $INPUT_DIR/*.gtf \
--output $RSEQC_DIR/$BED