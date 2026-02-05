#!/bin/bash -l

#$ -N cutadapt                  # Set Job Name
#$ -l mem_free=40G              # Request memory
#$ -j y                         # Merge standard output and standard error
#$ -cwd                         # Set current working directory
#$ -t 1-12                      # Submit an array job with n tasks; n=number of samples = number of fastq files/2 if paired end


# NOTE: You SHOULD make directories with appropriate library names within "cutadapt_results" folder
#proj="CRISPR_Jinfen"
proj="CRISPR_Prince"
libraries=(DTKPA1 DTKPA2 DTKPB1)

for library in ${libraries[*]}
do

INPUT_DIR=/common/theodorescudlab/Sequencing_Data/$proj/$library
OUTPUT_DIR=$HOME/scratch/$proj/cutadapt_results/$library 

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

# You need to add path for Conda
PATH=$HOME/miniconda3/bin/:$PATH
source activate MAGECK

# Read about cutadapt below
# https://cutadapt.readthedocs.io/en/stable/guide.html

# --quality-cutoff 15,10 removes low quality bases BEFORE adapter trimming from 5'end (score < 15) and 3' end (score < 10)
# --cut removes specified number of bases BEFORE adapter trimming from 3' end if negative and 5' end if positive
# --length removes specified number of bases AFTER adapter trimming from 3' end if negative and 5' end if positive

# -a removes adpater from 3' end of read1
# -A removes adpater from 3' end of read2
# -g removes adpater from 5' end of read1
# -G removes adpater from 5' end of read2

# --overlap  : the minimum length of overlap between read and adapter to be achieved for trimming to take place
# --times    : the number of occurances of adapter to be removed. Default=1
# --action   : action to be performed when matching adapter is found in read. Can be trim,retain,mask,lowercase,none. 
# --poly-a   : removes polyA tails from read1 and polyT heads from read2 AFTER adapter trimming.
# --trim-n   : removes flanking N bases from each read AFTER adapter trimming. 
# Use --quality-cutoff to remove them before adapter trimming since N bases typically also have a low quality value associated with them.

# --revcomp  : checks for adapter in read and also it's rev comp. 
# DO NOT USE --revcomp and MAKE SURE adapter is in read and not its reverse complement.

# --minimum-length : discards reads shorter than minimum length
# --maximum-length : discards reads longer than maximum length

MIN_READ_LENGTH=17
#MAX_READ_LENGTH=25
REMOVE_BEFORE_ADAPTER_TRIMMING=0
READ1_5_PRIME_QUALITY_CUTOFF=10
READ1_3_PRIME_QUALITY_CUTOFF=10
READ2_5_PRIME_QUALITY_CUTOFF=10
READ2_3_PRIME_QUALITY_CUTOFF=10
MIN_READ_ADAPTER_OVERLAP=10
TIMES=1

# For DTKPA1, DTKPA2, Guides are in R1.fastq.gz
# 5' and 3' adapters are different from DTKPB1
if [ $library == DTKPA1 ] || [ $library == DTKPA2 ]
then

# Create an array of input files
# Usually either R1 or R2 has the guide sequences. So, ignore the other fastq.gz file.
input1=($INPUT_DIR/*_R1.*.gz)
#input2=($INPUT_DIR/*_1.*.gz)

# Use the SGE_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value
index=$(($SGE_TASK_ID-1))
taskinput1=${input1[$index]}
taskoutput1=${taskinput1//$INPUT_DIR/$OUTPUT_DIR}
taskoutput3=${taskinput1//.fq.gz/.tsv}
taskoutput3=${taskoutput3//.fastq.gz/.tsv}

#NOTE: Make sure there is no space after the "comma" in --quality-cutoff
cutadapt \
--cut $REMOVE_BEFORE_ADAPTER_TRIMMING \
--quality-cutoff $READ1_5_PRIME_QUALITY_CUTOFF,$READ1_3_PRIME_QUALITY_CUTOFF \
--overlap $MIN_READ_ADAPTER_OVERLAP \
--times $TIMES \
--quality-base 33 \
-g GTATCCCTTGGAGAACCACCTTGTTG \
-a GTTTAAGAGCTAAGCTGGAAA \
--poly-a \
--trim-n \
--minimum-length $MIN_READ_LENGTH \
--output $taskoutput1 \
$taskinput1

# For DTKPB1, Guides are in R2.fastq.gz
# 5' and 3' adapters are different from DTKPA1 and DTKPA2
elif [ $library == DTKPB1 ]
then

input1=($INPUT_DIR/*_R2.*.gz)

index=$(($SGE_TASK_ID-1))
taskinput1=${input1[$index]}
taskoutput1=${taskinput1//$INPUT_DIR/$OUTPUT_DIR}
taskoutput3=${taskinput1//.fq.gz/.tsv}
taskoutput3=${taskoutput3//.fastq.gz/.tsv}

cutadapt \
--cut $REMOVE_BEFORE_ADAPTER_TRIMMING \
--quality-cutoff $READ1_5_PRIME_QUALITY_CUTOFF,$READ1_3_PRIME_QUALITY_CUTOFF \
--overlap $MIN_READ_ADAPTER_OVERLAP \
--times $TIMES \
--quality-base 33 \
-g TGCTGTTTCCAGCTTAGCTCTTAAAC \
-a CAACAAGGTGGTTCTCCAAGG \
--poly-a \
--trim-n \
--minimum-length $MIN_READ_LENGTH \
--output $taskoutput1 \
$taskinput1

else
echo $library
fi

done




# --info-file=info1.tsv \
# $taskinput2
#--maximum-length $MAX_READ_LENGTH \
#-g 
#-A \
#-G \
#-p $taskoutput2 \
# -Q $READ2_5_PRIME_QUALITY_CUTOFF, $READ2_3_PRIME_QUALITY_CUTOFF \
#--too-short-output $taskoutput3 \

# NOTE: First time, run the script on single file with --info-file=info.tsv
# https://cutadapt.readthedocs.io/en/stable/reference.html#info-file-format
# Extract certain columns to text file and analyze in R.
# cat ~/scratch/CRISPR_Jinfen/index_results/bowtie/info.tsv | cut -f 3,4,5,6 > t3.txt
# data <- read.table("...t3.txt", sep="\t", fill=TRUE)
# colnames(data) <- c("START", "END", "GUIDE", "ADAPTER")
# data <- data %>% 
# dplyr::mutate(START=as.numeric(START), END=as.numeric(END), LENGTH=END-START) %>% 
# dplyr::select(START, END, LENGTH, GUIDE, ADAPTER)
# Check length of adapter identified data %>% dplyr::count(LENGTH)
# Check starting position of adapters data %>% dplyr::count(START)
# Choose --minimum-length based on this result