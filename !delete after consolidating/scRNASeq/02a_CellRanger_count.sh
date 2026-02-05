#!/bin/bash -l

#$ -N CellRanger_count      # Set Job Name
#$ -l mem_free=120G         # Request memory
#$ -j y                     # Merge standard output and standard error
#$ -cwd                     # Set current working directory
#$ -t 1-25                  # Submit an array job with n tasks
# NOTE: n = number of samples; not number of fastqs. Each sample may 
# have 2 fastqs L001 and L002 if each sample was run on 2 lanes

# NOTE: You need to run cellranger count from your preferred directory as 
# there is no option to store output to a specific directory.

# NOTE: cellranger count automatically detects the chemistry and stores 
# this info in output folder of cellranger count
# Eg: FT1/SC_RNA_COUNTER_CS/SC_MULTI_CORE/DETECT_COUNT_CHEMISTRY/fork0/chnk0/_stdout.txt
# IF there is error in auto-detection, refer to the txt file above and re-run cellranger count with one of
# the chemistries as defined in 
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count
# Setting --chemistry=threeprime will auto-detect if it is 3'v1 or 3'v2 or 3'v3.

# Mouse datasets:
#proj="scRNASeq_BBN_C57B6"      
#proj="scRNASeq_BBN_Rag"
#proj="scRNASeq_GSE164557"
#proj="scRNASeq_GSE217093"
#proj="scRNASeq_Jinfen"
#proj="scRNASeq_NA13_CDH12_C57B6"
#proj="scRNASeq_Jyoti"
proj="scRNASeq_Krizia"

##proj="scRNASeq_GSE132042"
##proj="scRNASeq_GSE129845"

GENOME_DIR=$HOME/NGSTools/refdata-gex-GRCm39-2024-A    # use for mouse

# Human datasets:
#proj="scRNASeq_Chen"
#proj="scRNASeq_GSE222315"
#proj="scRNASeq_HRA003620"

##proj="scRNASeq_Simon"
##proj="scRNASeq_GSE267718"

#GENOME_DIR=$HOME/NGSTools/refdata-gex-GRCh38-2024-A   # use for human

cd $HOME/scratch/$proj/cellranger
FASTQ_DIR=/common/theodorescudlab/Sequencing_Data/$proj/01.RawData/
FASTQ_DIR=$HOME/scratch/$proj/01.RawData/

# Create an array of sample names using R2 fastq files. 
# Use ${array_name[index]} to verify that sample names are correct
# NOTE: L001 and L002 means same library was run on 2 lanes. So, both these fastqs
# must be processed simultaneously for each sample. Although we specify only 1 of 
# these fastqs as inputs, we use the basename() to get filename and %%_* to remove
# all characters following _ in the filename from end of string. This way we feed
# the sample name that is common to both L001 and L002 fastqs to cellranger.
# Also, R2 file is the one that has our sequence not R1 or I1 or I2 fastqs.
# https://www.linuxjournal.com/content/bash-arrays

# https://tldp.org/LDP/abs/html/string-manipulation.html
# ${string##substring} Deletes longest match of substring from front of string.
# ${string%%substring} Deletes longest match of substring from back of string.
# ${string#substring} Deletes shortest match of substring from front of string.
# ${string%substring} Deletes shortest match of substring from back of string.
# ${string/substring/replacement} replaces 1st first match
# ${string//substring/replacement} replaces all matches
# ${string/#substring/replacement} replaces 1st first from front of string.
# ${string/%substring/replacement} replaces 1st first from back of string.

# Example: string = "N5_S5_L001_R2_001.fastq.gz"  ; substring = "_*";
# Possible matches are:
# "_S5_L001_R2_001.fastq.gz" ; "_L001_R2_001.fastq.gz" ; "_R2_001.fastq.gz" ; "_001.fastq.gz"
# So, longest substring "_S5_L001_R2_001.fastq.gz" is removed and we are left with "N5" 

#inputs=($FASTQ_DIR*P3T_*L001_R2*.gz $FASTQ_DIR*P4T*L001_R2*.gz)
inputs=($FASTQ_DIR*L001_R2*.gz)
# Check if inputs are correct using "echo ${inputs[*]}"

# Use the SGE_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value.
# Mathematical calculations MUST be enclosed in double round brackets ((5-4)).
# Use $ for assigning the computed value to a variable 
index=$(($SGE_TASK_ID-1))

# Remove path and extract only the sample ID using basename(). 
# basename() cannot be run on multiple files. Output of basename() is NOT an array.
#filename=$(basename "${inputs[$index]}")
#filename=$(basename ${inputs[$index]}) # SRR14037759_B1246-GEX_S15_L001_R2_001.fastq.gz
#taskinput=${filename%%_*}             	# SRR14037759

# Alternatively, use // to remove path from mutliple files and save as array.
# NOTE: Your array is seen by shell as elements on a single line separated by spaces
# BUT sort() expects input to be on separate lines.
# tr ' ' '\n' : Convert all spaces to newlines. 
# sort -u     : Sort and retain only unique elements
# tr '\n' ' ' : Convert the newlines we added in earlier back to spaces.
# $(...)      : Command Substitution

filenames=(${inputs[*]//$FASTQ_DIR/})    # B1246-GEX25_S15_L001_R2_001.fastq.gz 
taskinputs=(${filenames[*]%%_*})         # B1246-GEX25 
taskinputs=($(echo "${taskinputs[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
taskinput=${taskinputs[$index]}
  
# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name   : $JOB_NAME"
echo "Job ID     : $JOB_ID"  
echo "TaskID     : $SGE_TASK_ID"
echo "index      : $index"
echo "filename   : $filename"
echo "taskinput  : $taskinput"
echo "=========================================================="

# You need to add path for Cellranger
PATH=$HOME/NGSTools/cellranger-9.0.1/bin/:$PATH

# Run cellranger count
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count
# All analysis before CellRanger 7 had --include-introns=false by default. RECOMMENDED to set true.
# Set create-bam to true. velocyto needs bam files created by cell ranger count.
# DO NOT use --nosecondary so that we can see clusters generated by cellranger using html file.

#SAMPLES=(FC1 FT1 FT2 MC1 MC2 MT1 MT2 MT3 MT4 MT5 MT6)
#for taskinput in ${SAMPLES[*]}
#do
#echo $taskinput

cellranger count \
--id=$taskinput \
--transcriptome=$GENOME_DIR \
--fastqs=$FASTQ_DIR \
--sample=$taskinput \
--create-bam=true \
--include-introns=true \
--chemistry=auto \
--check-library-compatibility=true \
--cell-annotation-model=auto
#done


# # Once cellranger count is complete, run the code below once to 
# # copy the raw_feature_bc_matrix folders for Seurat analysis
# # (${array[*]/%/suffix}) to add suffix to add array elements
# # (${array[*]/#/prefix}) to add prefix to add array elements
# LOCATION=$HOME/scratch/$proj
# SAMPLES=($(ls $LOCATION/cellranger/))
# SAMPLES_WITH_PATH=($(ls $LOCATION/cellranger/* -d))

# CURRENT_SAMPLES_WITH_PATH=(${SAMPLES_WITH_PATH[*]/%//outs/raw_feature_bc_matrix}) 
# NEW_SAMPLES_WITH_PATH=(${SAMPLES[*]/#/$LOCATION/raw_feature_bc_matrix/})
# for index in ${!SAMPLES[*]}
# do
  # cp -r ${CURRENT_SAMPLES_WITH_PATH[$index]} ${NEW_SAMPLES_WITH_PATH[$index]}
# done

# CURRENT_SAMPLES_WITH_PATH=(${SAMPLES_WITH_PATH[*]/%//outs/filtered_feature_bc_matrix}) 
# NEW_SAMPLES_WITH_PATH=(${SAMPLES[*]/#/$LOCATION/filt_feature_bc_matrix/})
# for index in ${!SAMPLES[*]}
# do
  # cp -r ${CURRENT_SAMPLES_WITH_PATH[$index]} ${NEW_SAMPLES_WITH_PATH[$index]}
# done

# rm -r $HOME/scratch/$proj/cellranger/*/SC_RNA_COUNTER_CS