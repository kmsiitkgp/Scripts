#!/bin/bash -l

#$ -N FastQC                    # Set Job Name
#$ -l mem_free=8G               # Request memory
#$ -j y                         # Merge standard output and standard error
#$ -cwd                         # Set current working directory
#$ -t 1-6                       # Submit an array job with 1 task for each fq.gz file. 
# Find n using: ls $HOME/scratch/RNASeq_$proj/raw_reads/*.gz | wc -l
# Note: Adjust number of input files properly

proj="WES_Hany"
OUTPUT_DIR=$HOME/scratch/$proj/fastqc_results/post_trimming/   #output directory to store fastqc results
mkdir $OUTPUT_DIR

# Create an array of input files
input=($HOME/scratch/$proj/cutadapt_results/*.gz)

# Use the SGE_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value
index=$(($SGE_TASK_ID-1))
taskinput=${input[$index]}

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date    : $(date)"
echo "Job name      : $JOB_NAME"
echo "Job ID        : $JOB_ID"  
echo "SGE TASK ID   : $SGE_TASK_ID"
echo "Index         : $index"
echo "Task input    : $taskinput"
echo "Project       : $proj"
echo "Output Folder : $OUTPUT_DIR"
echo "=========================================================="

# You need to add path for FastQC, java, Conda
PATH=$HOME/NGSTools/FastQC/:$PATH
PATH=$HOME/NGSTools/jdk-23.0.1/bin/:$PATH
PATH=$HOME/miniconda3/bin/:$PATH

# You have to use source activate instead of conda activate
source activate MAGECK

# --outdir=~/Hany/fastqc_results/ wont work. ~ isnt recognized by fastqc
fastqc $taskinput --outdir=$OUTPUT_DIR

# Run multiqc once after fastqc is complete on all samples
# multiqc --force \
# --filename MULTIQC_Report \
# --outdir=$OUTPUT_DIR \
# --fullnames \
# $HOME/scratch/$proj/

# Interpreting fastqc output
# https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/qc_fastqc_assessment.html

#“Per base sequence quality” plot:
# This plot provides the distribution of quality scores at each position in the read across all reads. The y-axis gives the quality scores, while the x-axis represents the position in the read. The color coding of the plot denotes what are considered high, medium and low quality scores. A drop in quality towards the ends of the reads, which could be explained by signal decay or phasing is expected for Illumina sequencing.  Any sudden drop in quality or a large percentage of low quality reads across the read could indicate a problem at the facility. NOTE: The quality of the nucleotide base calls are related to the signal intensity and purity of the fluorescent signal. Low intensity fluorescence or the presence of multiple different fluorescent signals (due to overclustering in the flow cell) can lead to a drop in the quality score assigned to the nucleotide. 

#“Per sequence quality scores” plot:
# This plot gives you the average quality score on the x-axis and the number of sequences with that average on the y-axis. We hope the majority of our reads have a high average quality score with no large bumps at the lower quality values.

# “Per base sequence content” plot:
# This plot always gives a FAIL for RNA-seq data because the first 10-12 bases result from the ‘random’ hexamer priming that occurs during RNA-seq library preparation. This priming is not as random as we might hope, giving an enrichment in particular bases for these intial nucleotides.

# “Per sequence GC content” plot :
# This plot gives the GC distribution over all sequences. Generally is a good idea to note whether the GC content of the central peak corresponds to the expected % GC for the organism. Also, the distribution should be normal unless over-represented sequences (sharp peaks on a normal distribution) or contamination with another organism (broad peak).

# “Sequence Duplication Levels” plot:
# This plot explores numbers of duplicated sequences in the library. This plot can help identify a low complexity library, which could result from too many cycles of PCR amplification or too little starting material. For RNA-seq we don’t normally do anything to address this in the analysis, but if this were a pilot experiment, we might adjust the number of PCR cycles, amount of input, or amount of sequencing for future libraries.

# “Overrepresented sequences” table :
#  This table displays the sequences (at least 20 bp) that occur in more than 0.1% of the total number of sequences. This table aids in identifying contamination, such as vector or adapter sequences. If the %GC content was off in the above module, this table can help identify the source. If not listed as a known adapter or vector, it can help to BLAST the sequence to determine the identity.