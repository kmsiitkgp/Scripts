#!/bin/bash -l

#$ -N RSEM_count              	# Set Job Name
#$ -l mem_free=40G              # Request memory
#$ -j y                         # Merge standard output and standard error
#$ -cwd                         # Set current working directory
#$ -t 1-12                      # Submit an array job with n tasks where n = number of samples = number of fastq files/2 if paired end

#STOP. Check strandedness from STAR output (ReadsPerGene.out.tab) before running this script

# Follow / at end of variables declareed below to avoid errors in array.
species="Mouse"
proj="RNASeq_Hany_Male_ImmuneEditing"
proj="RNASeq_Hany_CRISPR_LOY"
proj="RNASeq_Bhowmick"

#species="Human"
#proj="RNASeq_Vera"
INPUT_DIR=$HOME/scratch/$proj/alignment_results 		# directory with STAR co-ordinate sorted BAM files
OUTPUT_DIR=$HOME/scratch/$proj/count_results			# directory to store TPM counts
GTF=$HOME/NGSTools/Reference_Genomes/$species/*.gtf		# path to GTF file we downloaded
RSEM_REF=$HOME/NGSTools/Reference_Genomes_RSEM/$species/RSEM_

# Create an array of input BAM filenames with path
input=($INPUT_DIR/*Aligned.toTranscriptome.out.bam)

# Create an array of filenames with output path
# Create 3 files: unstranded, stranded +, stranded -
read=(${input[*]//$INPUT_DIR/$OUTPUT_DIR})
output=(${read[*]//Aligned.toTranscriptome.out.bam/})

# Use the SGE_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value
index=$(($SGE_TASK_ID-1))

# Declare input and output filenames with path for each job
taskinput=${input[$index]}
taskoutput=${output[$index]}

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date              : $(date)"
echo "Job name                : $JOB_NAME"
echo "Job ID                  : $JOB_ID"  
echo "SGE TASK ID             : $SGE_TASK_ID"
echo "Input BAM	              : $taskinput"
echo "Output TPM	          : $taskoutput"
echo "=========================================================="

# You need to add path for Conda and Sambamba
PATH=$HOME/miniconda3/bin/:$PATH
source activate MAGECK

# Check the ReadsPerGene.out.tab file to determine strandedness
# column 1: gene ID
# column 2: counts for unstranded RNA-seq
# column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
# column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)
# If unstranded, column 2 > column 3 AND column 2 > column 4 AND column 3 ~ column 4
# If stranded forward, column 3 > column 4
# If stranded reverse, column 4 > column 3

STRAND=reverse
rsem-calculate-expression \
--no-bam-output \
--alignments \
--paired-end \
--strandedness $STRAND \
$taskinput $RSEM_REF $taskoutput