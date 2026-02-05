#!/bin/bash -l

#$ -N HTSeq_count               # Set Job Name
#$ -l mem_free=40G              # Request memory
#$ -j y                         # Merge standard output and standard error
#$ -cwd                         # Set current working directory
#$ -t 1-18                      # Submit an array job with n tasks where n = number of samples = number of fastq files/2 if paired end

# Follow / at end of variables declareed below to avoid errors in array.
proj="RNASeq_Hany_Antigen"
species="Mouse"
OUTPUT_DIR=$HOME/scratch/$proj/count_results		#output directory to store read counts
INPUT_DIR=$HOME/scratch/$proj/alignment_results 	#input directory with STAR co-ordinate sorted BAM files
GTF=$HOME/NGSTools/Reference_Genomes/$species/*.gtf		#path to GTF file we downloaded

# Create an array of input BAM filenames with path
input=($INPUT_DIR/*Aligned.sortedByCoord.out.bam)

# Create an array of filenames with output path
# Create 3 files: unstranded, stranded +, stranded -
read=(${input[*]//$INPUT_DIR/$OUTPUT_DIR})
output1=(${read[*]//Aligned.sortedByCoord.out.bam/.txt})
output2=(${read[*]//Aligned.sortedByCoord.out.bam/_pos.txt})
output3=(${read[*]//Aligned.sortedByCoord.out.bam/_neg.txt})

# Use the SGE_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value
index=$(($SGE_TASK_ID-1))

# Declare input and output filenames with path for each job
taskinput=${input[$index]}
taskoutput_un=${output1[$index]}
taskoutput_pos=${output2[$index]}
taskoutput_neg=${output3[$index]}


# Keep track of information related to the current job
echo "=========================================================="
echo "Start date              : $(date)"
echo "Job name                : $JOB_NAME"
echo "Job ID                  : $JOB_ID"  
echo "SGE TASK ID             : $SGE_TASK_ID"
echo "Task input              : $taskinput"
echo "Unstranded file         : $taskoutput_un"
echo "Stranded positive file  : $taskoutput_pos"
echo "Stranded negative file  : $taskoutput_neg"
echo "=========================================================="

# You need to add path for Conda and Sambamba
PATH=$HOME/miniconda3/bin/:$PATH
PATH=$HOME/NGSTools/sambamba-0.8.2-linux-amd64-static/:$PATH
# You have to use source activate instead of conda activate
source activate R
# Check which version of htseq-count is being used
# htseq-count is a python package. We can use conda to invoke htseq-count
which htseq-count
htseq-count --version
which sambamba

# Create a bai file for bam file being processed. Having no index file WONT affect HTSeq results
# It will just give a warning that index file is absent
cd $INPUT_DIR
sambamba index $taskinput

# HTSeq-count is counting the reads, which align to the given exons. 
# If stranded kit was used, choose --stranded=yes or --stranded=reverse
# If unstranded kit was used, choose --stranded=no
# If --stranded=yes, HTSeq-count checks whether the reads are in the same orientation as the transcript.
# If --stranded=reverse, HTSeq-count checks whether the reads are in the reverse orientation as the transcript.
# The vast majority of stranded libraries (e.g. Illumina's TruSeq Stranded protocol) produces libraries, 
# which are in reverse orientation to the transcripts' one. So, use --stranded=reverse.
# However, most bulk RNA seq which are used solely for gene expression are unstranded.

# You can easily determine if which type of protocol was used for library preparation by comparing the 
# count file generated for each of 3 options: --stranded=no, --stranded=yes, --stranded=reverse
# If unstranded kit was used, --stranded=yes (AND) --stranded=reverse counts will be almost equal
# If stranded kit was used, either --stranded=yes (OR) --stranded=reverse will have all the reads.

htseq-count \
--format=bam \
--order=pos \
--stranded=no \
--type=exon \
--idattr=gene_id \
--additional-attr=gene_name\
--mode=union \
--nonunique=none \
--samout-format=sam \
$taskinput \
$GTF > $taskoutput_un

htseq-count \
--format=bam \
--order=pos \
--stranded=yes \
--type=exon \
--idattr=gene_id \
--additional-attr=gene_name\
--mode=union \
--nonunique=none \
--samout-format=sam \
$taskinput \
$GTF > $taskoutput_pos

htseq-count \
--format=bam \
--order=pos \
--stranded=reverse \
--type=exon \
--idattr=gene_id \
--additional-attr=gene_name\
--mode=union \
--nonunique=none \
--samout-format=sam \
$taskinput \
$GTF > $taskoutput_neg