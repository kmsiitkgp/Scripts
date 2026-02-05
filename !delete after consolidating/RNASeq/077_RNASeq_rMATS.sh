#!/bin/bash -l

#$ -N rMATS                     # Set Job Name
#$ -l mem_free=40G              # Request memory
#$ -j y                         # Merge standard output and standard error
#$ -cwd                         # Set current working directory
#$ -t 1-10                      # Submit an array job with n tasks where n = number of samples = number of fastq files/2 if paired end

# Follow / at end of variables declareed below to avoid errors in array.
proj="RNASeq_GSE75192"
species="Mouse"
OUTPUT_DIR=$HOME/scratch/$proj/rMATS_results		    #output directory to store rMATS results
INPUT_DIR=$HOME/scratch/$proj/STAR_alignment_results 	#input directory with STAR co-ordinate sorted BAM files
GTF=$HOME/NGSTools/RNASeq_Reference_Genomes/$species/*.gtf		#path to GTF file we downloaded

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

MAX_READ_LENGTH=$(($(zcat $HOME/common/$proj/*q.gz | head -2 | awk 'END {print}'| wc -c)-1))

#-t {paired,single}
#--libType {fr-unstranded,fr-firststrand,fr-secondstrand}
python rmats.py \
--gtf \
--b1 $! \
--b2 $! \
--od $OUTPUT_DIR\
--tmp $OUTPUT_DIR\
--task both \
--readLength $MAX_READ_LENGTH \
--variable-read-length \
--anchorLength 1 \
-t paired \
--libType fr-unstranded \
#--paired-stats  # use ONLY if samples in --b1 and --b2 are paired
#--statoff       # use ONLY if using --b1 without --b2
