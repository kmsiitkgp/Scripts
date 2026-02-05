#!/usr/bin/env bash

#$ -N GDC_Download    # Set Job Name
#$ -l mem_free=64G    # Request memory
#$ -j y               # Merge standard output and standard error
#$ -cwd               # Set current working directory
#$ -t 1		          # Submit an array job with n tasks. Set n to number of SRR ids in the SRR_Acc_List.txt

# Add path for GDC Client
PATH=$HOME/NGSTools/gdc-client_v1.6.1_Ubuntu_x64/:$PATH
gdc-client --version

DOWNLOAD_DIR=$HOME/scratch/ICB_TCGA_datasets/original/downloads
OUTPUT_DIR=$HOME/scratch/ICB_TCGA_datasets/original/raw_counts
MANIFEST_FILE=$HOME/projects/TCGA_GDC/gdc_manifest.2025-03-06.txt

# Download the TCGA data into data folder
gdc-client download \
--dir $DOWNLOAD_DIR \
--manifest $MANIFEST_FILE \
--no-related-files \
--no-annotations \
--latest

# If the downloaded files are already in .tsv format, just copy them to raw_counts folder 
mv $DOWNLOAD_DIR/*/*_star_*.tsv $OUTPUT_DIR

# IMPORTANT NOTE: Make sure ls $OUTPUT_DIR | wc -l gives same number as number of files.
# I have noticed  once instead of 11499 tsv files, only 11442 were downloaded.

# If the downloaded files are in .gz format, use code below to extract them and store in counts folder 
# gunzip -c $HOME/projects/TCGA-GDC/data/*/*.htseq.counts.gz > $HOME/projects/TCGA-GDC/counts/

# Run either above line (i.e. gunzip -c ... > ... ) (OR) below 2 lines
# gunzip $HOME/projects/TCGA_GDC/data/*/*.htseq.counts.gz
# mv $HOME/projects/TCGA_GDC/data/*/*_star_* $HOME/projects/TCGA_GDC/counts/

# Download the counts folder and analyze in R

# # You can run this script using qsub . If not, remove the qsub paramters (lines 3-7) and run code below
# #sh ~/projects/TCGA-GDC/TCGA_GDC.sh




