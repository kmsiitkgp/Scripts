#!/bin/bash -l

#$ -N CITE_Seq_count       # Set Job Name
#$ -l mem_free=64G         # Request memory.
#$ -j y                    # Merge standard output and standard error
#$ -cwd                    # Set current working directory
#$ -t 1-13                 # Submit an array job with n tasks where n = number of unique SRR ids.
# This MUST also be equal to the number of *whitelist.csv files that have barcode info for each batch

# You need to add path for Conda and use source activate instead of conda activate
PATH=$HOME/miniconda3/bin/:$PATH
source activate citeseq

# Check which version of htseq-count is being used
# htseq-count is a python package. We can use conda to invoke htseq-count
which CITE-seq-Count
CITE-seq-Count --version

# NOTE: First, run CellRanger count on GEX fastq.gz 
# Next, use barcodes from filtered cells as whitelist input for CITE-seq-Count
# Then, run CITE-seq-Count on HTO fastq.gz
# Finally, demux the cells using R

# NOTE: If you get error below, follow recommended steps below:
# conda create -n citeseq python=3.7.16
# conda activate citeseq
# pip install CITE-seq-Count==1.4.5
# conda install pandas=1.3.5

# Traceback (most recent call last):
  # File "/home/kailasamms/miniconda3/envs/NGS/bin/CITE-seq-Count", line 8, in <module>
    # sys.exit(main())
  # File "/home/kailasamms/miniconda3/envs/NGS/lib/python3.10/site-packages/cite_seq_count/__main__.py", line 603, in main
    # io.write_dense(
  # File "/home/kailasamms/miniconda3/envs/NGS/lib/python3.10/site-packages/cite_seq_count/io.py", line 48, in write_dense
    # pandas_dense = pd.DataFrame(sparse_matrix.todense(), columns=columns, index=index)
  # File "/home/kailasamms/miniconda3/envs/NGS/lib/python3.10/site-packages/pandas/core/frame.py", line 639, in __init__
    # raise ValueError("columns cannot be a set")
# ValueError: columns cannot be a set

proj="scRNASeq_Simon"
FASTQ_DIR=$HOME/common/$proj/HTO
WHITELIST_DIR=$HOME/projects/scRNASeq
OUTPUT_DIR=$HOME/scratch/$proj/CITESeq_results
HTO_TAGS=$HOME/projects/scRNASeq/$proj.HTO_tags.csv

# Create an array of R1 and R2 fastq.gz files with path.
input1=($FASTQ_DIR/*R1*.gz)
input2=($FASTQ_DIR/*R2*.gz)
input3=($WHITELIST_DIR/$proj*whitelist.csv)

# Verify if the array stores correctly
# echo ${input1[0]}
# echo ${input2[0]}
# echo ${input3[0]}

# Use the SGE_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value
index=$(($SGE_TASK_ID-1))
taskinput1=${input1[$index]}
taskinput2=${input2[$index]}
whitelist=${input3[$index]}
barcode_count=$(cat $whitelist | wc -l)

# Create an array of whitelist filenames with path
files=(${input3[*]//_whitelist.csv/})
folder=$(basename ${files[$index]})
taskoutput=$OUTPUT_DIR/$folder
mkdir $taskoutput

# Keep track of information related to the current job
echo "======================================================================================================="
echo "Start date   : $(date)"
echo "Job name     : $JOB_NAME"
echo "Job ID       : $JOB_ID"  
echo "SGE TASK ID  : $SGE_TASK_ID"
echo "Task input 1 : $taskinput1"
echo "Task input 2 : $taskinput2"
echo "Task output  : $taskoutput"
echo "Whitelist    : $whitelist"
echo "Barcode      : $barcode_count" 
echo "======================================================================================================="

# Using \ enables us to split the long command into easy to read format.
# NOTE: Cellranger count uses hamming distance of 1 for UMI as well as barcode correction.
# So, set UMI_ERRORS=1 (instead of default value of 2) and set BC_ERRORS=1 

# NOTE: Use ONLY 1 thread. Else, it will split file into multiple parts and map quickly but use 800GB+
# memory while merging the parts back together and cluster will kill the job.

# NOTE: Use "HTO-A" instead of "HTO_A" in HTO_tags.csv file.
# If you use "HTO_A", it automatically gets changed to "HTO-A" in next step of analysis.
# This will create problems later. So, best use "HTO-A" from the beginning.

# NOTE: First test with 2 million reads. If you get low % of unmapped reads, proceed with full data set.
# --first_n 1000000. Check the yaml file to see % mapped vs unmapped
# Check with and without sliding window to see if unmapped reads is reduced significantly

CB_FIRST=1   # First nucleotide of cell barcode in read R1
CB_LAST=16   # Last nucleotide of the cell barcode in read R1.
UMI_FIRST=17 # First nucleotide of the UMI in read R1
UMI_LAST=28  # Last nucleotide of the UMI in read R1
BC_ERRORS=1  # How many errors are allowed between two cell barcodes to collapse them onto one cell.
UMI_ERRORS=1 # How many errors are allowed between two umi within the same cell and TAG to collapse.
FIRST_N=1000000

CITE-seq-Count \
--read1 $taskinput1 \
--read2 $taskinput2 \
--tags $HTO_TAGS \
--cell_barcode_first_base $CB_FIRST \
--cell_barcode_last_base $CB_LAST \
--umi_first_base $UMI_FIRST \
--umi_last_base $UMI_LAST \
--bc_collapsing_dist $BC_ERRORS \
--umi_collapsing_dist $UMI_ERRORS \
--whitelist $whitelist \
--expected_cells $barcode_count \
--threads 1 \
--unmapped-tags unmapped.csv \
--unknown-top-tags 100 \
--sliding-window \
--output $taskoutput



#--first_n $FIRST_N \

# # Once citeseq count is complete, run the code below once to 
# # copy the umi_count folders for HTO_Demux analysis
# # (${array[*]/%/suffix}) to add suffix to add array elements
# # (${array[*]/#/prefix}) to add prefix to add array elements 
# LOCATION=$HOME/scratch/$proj
# SAMPLES=($(ls $LOCATION/CITESeq_results/))
# SAMPLES_WITH_PATH=($(ls $LOCATION/CITESeq_results/* -d))
# CURRENT_SAMPLES_WITH_PATH=(${SAMPLES_WITH_PATH[*]/%//umi_count})
# NEW_SAMPLES_WITH_PATH=(${SAMPLES[*]/#/$LOCATION/raw_hto_bc_matrix/})

# for index in ${!SAMPLES[*]}
# do
  # cp -r ${CURRENT_SAMPLES_WITH_PATH[$index]} ${NEW_SAMPLES_WITH_PATH[$index]}
# done  

# Print first 1000 HTO tags and export them
# zcat ~/common/scRNASeq_Simon/HTO/B1-HTO_S00_L000_R2_001.fastq.gz | head -12000 | awk 'NR % 4 == 2' > output.txt

# Open as excel file and make sure majority of R2 reads first 15 bases match with HTO tags. Color code these R2 reads.
# Copy unmatched R2 reads and paste to new sheet. Sort by cell values. 
# Write value 1 on cell B1. On cell B2 type this formula "=IF(EXACT(A1,A2),SUM(B1,1),1)". This calculates frequency of reads
# Copy values to new column, sort these values from highest to lowest, remove duplicate reads. 
# This gives an idea of the frequency of unmapped reads. You will notice top unmapped reads differ from the HTO tags by 1 base.

# UMI count vs Read count:
# BC1-UMI1-...HTO-A..after library preparation, there could be 10 copies of this read which were generated from
# a single molecule i.e. PCR copies..So read count will be 10 after sequencing, but UMI count will be 1 as there
# is only 1 unique molecule from which all PCR copies were generated.
