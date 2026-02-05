#!/bin/bash -l

#$ -N bowtie2                   # Set Job Name
#$ -l mem_free=40G              # Request memory
#$ -j y                         # Merge standard output and standard error
#$ -cwd                         # Set current working directory
#$ -t 1-2                      # Submit an array job with n tasks; n=number of samples = number of fastq files/2 if paired end

proj="WES_Hany"
species="Mouse"
INPUT_DIR=$HOME/scratch/$proj/cutadapt_results
INPUT_DIR=$HOME/scratch/$proj/raw_reads
OUTPUT_DIR=$HOME/scratch/$proj/alignment_results  
INDEXED_REF_DIR=$HOME/NGSTools/Reference_Genomes_Indexed_bowtie2/$species
base_name=$species.reference

# Create an array of input files
input1=($INPUT_DIR/*_1*.gz)
input2=($INPUT_DIR/*_2*.gz)

# Use the SGE_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value
index=$(($SGE_TASK_ID-1))
taskinput1=${input1[$index]}
taskinput2=${input2[$index]}
taskoutput=${taskinput1//$INPUT_DIR/$OUTPUT_DIR}
taskoutput=${taskoutput//.*/}
taskoutput_sam=$taskoutput.sam

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date    : $(date)"
echo "Job name      : $JOB_NAME"
echo "Job ID        : $JOB_ID"  
echo "SGE TASK ID   : $SGE_TASK_ID"
echo "Index         : $index"
echo "Task input 1  : $taskinput1"
echo "Task input 2  : $taskinput2"
echo "Task output   : $taskoutput"
echo "Project       : $proj"
echo "Output Folder : $OUTPUT_DIR"
echo "=========================================================="

# You need to add path for Conda
PATH=$HOME/miniconda3/bin/:$PATH
source activate MAGECK
PATH=$HOME/NGSTools/sambamba-0.8.2-linux-amd64-static/:$PATH

# if reads are paired, use -1 and -2. Else, use -U
# Align reads
bowtie2 \
--local \
-q \
-x $INDEXED_REF_DIR/$base_name \
-1 $taskinput1 \
-2 $taskinput2 \
-S $taskoutput_sam 
#2> $HOME/scratch/$proj/bowtie2.log

# NOTE: 
# Bowtie2 supports (i) gapped, ungapped (ii) local, global and (iii) single-end,paired-end alignment modes.
# Bowtie1 supports only ungapped alignment.
# Bowtie2 works best for reads that are greater than 50 bp.
# Bowtie1 works best for reads that are lesser than 50bp.

# Read about bowtie2
# https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#how-is-bowtie-2-different-from-bowtie-1
# https://hbctraining.github.io/Intro-to-ChIPseq-flipped/lessons/04_alignment_using_bowtie2.html
# By default, Bowtie2 will perform a global end-to-end read alignment, which aligns from the first to the last base of the read. This alignment is best for reads that have already been trimmed for quality and adapters (e.g. reads where nucleotide bases of poor quality or matching adapter sequences have been removed from the ends of the reads prior to alignment). However, Bowtie2 also has a local alignment mode, which, in contrast to end-to-end alignment, ignores portions at the ends of the reads that do match well to the reference genome. This is referred to as soft-clipping and allows for a more accurate alignment. The procedure can carry a small penalty for each soft-clipped base, but amounts to a significantly smaller penalty than mismatching bases. In contrast to trimming, which removes the unwanted sequence (hard-clipping), soft-clipping retains the soft-clipped base in the sequence and simply marks it. We will use this option since we did not trim our reads.

# To perform the Bowtie2 alignment, a genome index is required. The index is analagous to the index in a book. By indexing the genome, we have organized it in a manner that now allows for efficient search and retrieval of matches of the query (sequence read) to the genome. Bowtie2 indexes the genome with an FM Index based on the Burrows-Wheeler Transform method to keep memory requirements low for the alignment process.


# If you did not trim reads, use --local when runing Bowtie2. This will perform “soft-clipping” to ignore parts of the reads that are of low quality.

