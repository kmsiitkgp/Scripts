#!/bin/bash -l

#$ -N Get_Fastq       # Set Job Name
#$ -l mem_free=64G    # Request memory
#$ -j y               # Merge standard output and standard error
#$ -cwd               # Set current working directory
#$ -t 1-10            # Submit an array job with n tasks. Set n to number of SRR ids in the SRR_Acc_List.txt

# Cell Ranger requires FASTQ file names to follow the bcl2fastq file naming convention.
# [Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz
# Where Read Type is one of:
# I1: Sample index read (optional)
# I2: Sample index read (optional)
# R1: Read 1 (required)
# R2: Read 2 (required)

# incompatible file name: SRR9291388_1.fastq.gz
# compatible file name: SRR9291388_S1_L001_R1_001.fastq.gz
# Changing the file names will allow Cell Ranger (version >=2.1.1) to accept this data as inputs. 
# Note that only R1 and R2 FASTQ files are required for Cell Ranger. I1, I2 FASTQ files are optional. 

# You need to add path for Conda and use source activate instead of conda activate
PATH=$HOME/miniconda3/bin/:$PATH
source activate R

proj=scRNASeq_Vera
OUTPUT_DIR=$HOME/scratch/$proj

# Find SRR numbers (not SRX, SAMN, GSM numbers) for each sample from
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169379 and 
# https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA716349&o=bytes_l%3Aa
# Save SRR numbers to a text file SRR_Acc_List.txt
  
# Read a file line by line and store each line as an array
readarray -t input < $OUTPUT_DIR/SRR_Acc_List.txt
#mapfile -t input < $HOME/scratch/$proj/SRR_Acc_List.txt
# Verify if the array stores correctly
# echo ${input[*]}

# Use the SGE_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value
index=$(($SGE_TASK_ID-1))
taskinput=${input[$index]}

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name   : $JOB_NAME"
echo "Job ID     : $JOB_ID"  
echo "TaskID     : $SGE_TASK_ID"
echo "index      : $index"
echo "taskinput  : $taskinput"
echo "=========================================================="

# https://edwards.flinders.edu.au/fastq-dump/
# https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump
# https://www.biostars.org/p/12047/ 
# What is a spot in SRA format?
# Spot is location on flowcell which was imaged during sequencing.
# If you did paired end sequencing with 6bp barcode and 12 bp primer,
# there would be 4 reads coming from the spot:  
# 1 read for barcode,
# 1 read for primer, 
# 1 forward read and 
# 1 reverse read. 
# All 4 reads are merged as spot in SRA format

# NOTE: If "prefetch SRR3316476" works but "prefetch $taskinput" doesnt work. 
# make sure the text file has Unix (LF) and not Windows (CR LF) by opening in
# Notepad++

# STEP 1: PREFETCH
# NCBI recommends running prefetch, then faster-dump or fastq-dump.
# NOTE: prefetch download sra files which are then converted to fastq files using fasterq-dump
# Using \ enables us to split the long command into easy to read format.
prefetch $taskinput \
--max-size 420000000000 \
--progress \
--output-directory $OUTPUT_DIR

# STEP 2: CONVERT TO FASTQ
# You can check the composition of SRA files using
# https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&display=metadata
# Search SRR number in above link to see read length of each fastq file stored in the SRR

# # If downloading single/paired end bulk RNA seq data where each SRA file has upto 2 reads, 
# # you can use --split-3 and --skip-technical
# fasterq-dump $taskinput \
# --split-3 \
# --skip-technical \
# --progress \
# --outdir $OUTPUT_DIR

# If downloading single cell RNA Seq data where each SRA file can have upto 4 reads,
# you MUST use --split-files instead of --split-3 and --include-technical
# https://kb.10xgenomics.com/hc/en-us/articles/115003802691-How-do-I-prepare-Sequence-Read-Archive-SRA-data-from-NCBI-for-Cell-Ranger-
# NOTE: If you do not use --include-technical, then only R2 files will be saved.
fasterq-dump $taskinput \
--split-files \
--include-technical \
--progress \
--outdir $OUTPUT_DIR

# STEP 3: Convert to fastq.gz
gzip $OUTPUT_DIR/$taskinput*.fastq

# STEP 4: Remove unwanted files
# Remove prefetch folders and unpaired reads. 
# Single end reads will have _1.fastq.gz. Sometimes, they dont have _1 suffix.
# Paired end reads will have _1.fastq.gz and _2.fastq.gz
#rm -r $OUTPUT_DIR/$taskinput/

# # Downloading single cell fastq.gz files from Chinese databases
# wget -P $OUTPUT_DIR/ $taskinput

# # Create excel file with original filenames in column 1 and new files names in column 2
# # Copy and paste these 2 columns into a t.txt file in HPC cluster
# # Make sure the t.txt file has Unix(LF) ending shown in bottom right corner.
# cut -f 1 t.txt > t1.txt
# cut -f 2 t.txt > t2.txt
# readarray -t input1 < t1.txt
# readarray -t input2 < t2.txt
# for n in ${!input1[*]}
# do 
# mv ${input1[$n]} ${input2[$n]}
# done

# OLD COMMAND (OUTDATED)
# prefetch combined with faster-dump+gzip results are SAME as fastq-dump results below. 
# Only file sizes are larger with faster-dump since they are not gzipped.
# fastq-dump \
# --split-3 \
# --skip-technical \
# --origfmt \
# --clip \
# --readids \
# --dumpbase \
# --read-filter pass \
# --gzip \
# --outdir $OUTPUT_DIR/raw_reads $taskinput

# # OPTIONAL FOR LOOP
# # Sometimes a single fastq file from step 2 may take upto 300GB or more space.
# # In such cases, an array job will fail due to lack of space. 
# # If this happens, run Step 1 as job array but perform Step 2 through 4 using for loop below.
# # !SRR[*] gives array indices while SRR[*] gives the array elements
# proj=scRNASeq_Chen
# OUTPUT_DIR=$HOME/common/$proj
# SRR=(SRR12603781 SRR12603785 SRR12603786 SRR12603787 SRR12603788)
# for index in ${!SRR[*]}
# do
  # taskinput=${SRR[$index]}
  # echo $taskinput
  
  # fasterq-dump $taskinput \
  # --split-files \
  # --skip-technical \
  # --progress \
  # --outdir $OUTPUT_DIR
  
  # gzip $OUTPUT_DIR/$taskinput*.fastq  
# done

# # Checking length of reads in each fastq file
# proj=scRNASeq_Simon
# files=($HOME/common/$proj/SRR*I1*)
# #files=($HOME/common/$proj/SRR*R1*)
# #files=($HOME/common/$proj/SRR*R2*)
# for n in ${!files[*]}
# do
  # echo ${files[$n]}
  # zcat ${files[$n]} | head -2 | awk 'END {print}'| wc -c
# done

# echo $HOME/common/$proj/SRR*I1* | wc -w
# echo $HOME/common/$proj/SRR*R1* | wc -w
# echo $HOME/common/$proj/SRR*R2* | wc -w
