#!/bin/bash -l

#$ -N STAR_Align_Reads          # Set Job Name
#$ -l mem_free=40G              # Request memory
#$ -j y                         # Merge standard output and standard error
#$ -cwd                         # Set current working directory
#$ -t 1-24                     	# Submit an array job with n tasks; n=number of samples = number of fastq files/2 if paired end

# Follow / at end of variables declareed below to avoid errors in array.
proj="RNASeq_Hany_Male_ImmuneEditing"
proj="RNASeq_Hany_CRISPR_LOY"
proj="RNASeq_Sandrine.Supriya"
species="Mouse" 

#proj="RNASeq_Vera"
#proj="RNASeq_Manish_22RV1_ARCaPM"
#proj="RNASeq_Manish_22RV1_Xenograft"
proj="RNASeq_Manish_"
species="Human"

FASTQ_DIR=$HOME/scratch/$proj/01.RawData
OUTPUT_DIR=$HOME/scratch/$proj/alignment_results					    # directory to store read alignment results
GENOME_INDICES=$HOME/NGSTools/Reference_Genomes_STAR/$species/			# directory containing genome indices from previous step
GENOME_BED=$HOME/NGSTools/Reference_Genomes_RSEQC/$species/*bed

# Create an array of read1 and read2 filenames with path. If single end, then there will be no read2.
input1=($FASTQ_DIR/*_1.*q.gz)
input2=($FASTQ_DIR/*_2.*q.gz)

# Sometimes fastqs are labelled _R1 and _R2 instead of _1 and _2
if [ ${#input1[*]} == 1 ] && [ ${#input2[*]} == 1 ]
then
echo howdy
  input1=($FASTQ_DIR/*_R1.*q.gz)
  input2=($FASTQ_DIR/*_R2.*q.gz)
fi

# Create an array of filenames with output path
# These are prefix with output path. STAR will append "SortedByCoordinate.bam" to filename 
read1=(${input1[*]//$FASTQ_DIR/$OUTPUT_DIR})
output=(${read1[*]//_R1.*q.gz/})
output=(${output[*]//_1.*q.gz/.})

# Verify if the array stores correctly
echo ${input1[*]}

# Number of items in array
echo ${#input1[*]}

# Single end filenames wont have _2 or _R2. Hence, use below to create array.
if [ ${#input1[*]} == 1 ] || [ ${#input2[*]} == 1 ]
then
  input1=($FASTQ_DIR/*q.gz)
  read1=(${input1[*]//$FASTQ_DIR/$OUTPUT_DIR})
  output=(${read1[*]//.*q.gz/})
fi
  
# Use the SGE_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value
index=$(($SGE_TASK_ID-1))

# #For Imvigor010, there were 728 samples. HPC couldnt handle 728 jobs.
# # So, we run 15 samples in for loop in each job and submit 49 such jobs (15*49=735)
# START=$(((($index+1)*15)-15))
# END=$(((($index+1)*15)-1))
# STEP=1

# for index in $(seq $START $STEP $END)
# do
# echo $index
# #done

# Declare input and output filenames with path for each job
taskinput1=${input1[$index]}
taskinput2=${input2[$index]}
taskoutput=${output[$index]}

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date   	: $(date)"
echo "Job name     	: $JOB_NAME"
echo "Job ID       	: $JOB_ID"  
echo "SGE TASK ID	: $SGE_TASK_ID"
echo "Read 1		: $taskinput1"
echo "Read 2  		: $taskinput2"
echo "output	  	: $taskoutput"
echo "=========================================================="

# You need to add path for STAR and sambamba
PATH=$HOME/NGSTools/STAR-2.7.11b/bin/Linux_x86_64_static/:$PATH
PATH=$HOME/NGSTools/sambamba-1.0.1-linux-amd64-static/:$PATH

# Generate Genome Indices. You can write all lines below in a single line without \ also. 
# Using \ enables us to split the long command into easy to read format.

if((${#input2[*]}>1)) 
 then 
  STAR \
 --runMode alignReads \
 --genomeDir $GENOME_INDICES \
 --readFilesIn $taskinput1 $taskinput2 \
 --readFilesCommand zcat \
 --outFileNamePrefix $taskoutput \
 --outFilterMultimapNmax 10 \
 --outFilterMismatchNoverReadLmax 0.04 \
 --alignIntronMin 20 \
 --alignIntronMax 1000000 \
 --alignMatesGapMax 1000000 \
 --quantMode TranscriptomeSAM GeneCounts \
 --outSAMunmapped Within \
 --outSAMtype BAM SortedByCoordinate 
elif((${#input2[*]}==1)) 
 then 
  STAR \
 --runMode alignReads \
 --genomeDir $GENOME_INDICES \
 --readFilesIn $taskinput1 \
 --readFilesCommand zcat \
 --outFileNamePrefix $taskoutput \
 --outFilterMultimapNmax 10 \
 --outFilterMismatchNoverReadLmax 0.04 \
 --alignIntronMin 20 \
 --alignIntronMax 1000000 \
 --alignMatesGapMax 1000000 \
 --quantMode TranscriptomeSAM GeneCounts \
 --outSAMunmapped Within \
 --outSAMtype BAM SortedByCoordinate
fi

# Generate index files
sambamba index ${taskoutput}Aligned.sortedByCoord.out.bam

# You need to add path for Conda to use RSeQC
PATH=$HOME/miniconda3/bin/:$PATH

# You have to use source activate instead of conda activate
source activate MAGECK

# Calculate read distribution using RSeQC modules
read_distribution.py \
--input-file ${taskoutput}Aligned.sortedByCoord.out.bam \
--refgene $GENOME_BED > ${taskoutput}RSeQC.txt

# --outFilterMultimapNmax: default: 10
# maximum number of loci the read is allowed to map to. 
# Alignments (all of them) will be output only if the read maps to no more loci than this value. 
# Otherwise no alignments will be output, and the read will be counted as "mapped to too many loci" in the Log.final.out.

# --outFilterMismatchNoverReadLmax: default: 1.0
# alignment will be output only if its ratio of mismatches to *read* length is less than or equal to this value
# for 1x100b, max number of mismatches is 0.04*100=4 (mismatch/readlength = 0.04 i.e. 4 mismatch per 100 base)
# for 2x100b, max number of mismatches is 0.04*200=8 for the paired read

# --alignIntronMin: default: 21
# minimum intron size
# genomic gap is considered intron if its length>=alignIntronMin, otherwise it is considered Deletion

# --alignIntronMax: default: 0
# maximum intron size
# if 0, max intron size will be determined by (2^winBinNbits)*winAnchorDistNbins

# --alignMatesGapMax: default: 0
# maximum gap between two mates
# if 0, max intron gap will be determined by (2^winBinNbits)*winAnchorDistNbins

# Q: I found in ENCODE options: the max number of multiple alignments allowed for a read is 20 (i.e. --outFilterMultimapNmax 20). 
# Why it was not set as 1 to include only unique mapped reads for downstream analysis?
# A: Because counting tools such as featureCounts and HTSeq can exclude multi-mapping reads, 
# and other tools such as RSEM and eXpress can use the information from multi-mapped reads to provide more accurate count estimates.