#!/bin/bash -l

#$ -N SpaceRanger_count     # Set Job Name
#$ -l mem_free=120G         # Request memory
#$ -j y                     # Merge standard output and standard error
#$ -cwd                     # Set current working directory
#$ -t 1-2                   # Submit an array job with n tasks
# NOTE: n = number of samples; not number of fastqs. Each sample may 
# have 2 fastqs L001 and L002 if each sample was run on 2 lanes

# NOTE: You need to run spaceranger count from your preferred directory as 
# there is no option to store output to a specific directory.

# Refer link below to understand image input for your data
# https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/inputs/image-image-recommendation

#Mouse datasets:
#proj="Spatial_VisiumHD"
#GENOME=$HOME/NGSTools/refdata-gex-mm10-2020-A    # use for mouse
#PROBE_SET=$HOME/NGSTools/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv
#PROBE_SET=$HOME/NGSTools/Visium_Mouse_Transcriptome_Probe_Set_v2.0_mm10-2020-A.csv

#Human datasets:
#proj="visium_GSE171351" 
proj="Spatial_Bhowmick" 
GENOME=$HOME/NGSTools/refdata-gex-GRCh38-2020-A   # use for human
#PROBE_SET=$HOME/NGSTools/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv
PROBE_SET=$HOME/NGSTools/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv

cd $HOME/scratch/$proj/spaceranger_results
FASTQ_DIR=$HOME/scratch/$proj/01.RawData/

inputs=($FASTQ_DIR*L005_R2*.gz)
cyt_images=($FASTQ_DIR*CytAssist*.tif)
he_images=($FASTQ_DIR*HE*.tif)
# Check if inputs are correct using "echo ${inputs[*]}"

# Use the SGE_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value.
# Mathematical calculations MUST be enclosed in double round brackets ((5-4)).
# Use $ for assigning the computed value to a variable 
index=$(($SGE_TASK_ID-1))

# Create an array of sample names using R2 fastq files. 
# Use ${array_name[index]} to verify that sample names are correct
# NOTE: L001 and L002 means same library was run on 2 lanes. So, both these fastqs
# must be processed simultaneously for each sample. Although we specify only 1 of 
# these fastqs as inputs, we use the basename() to get filename and %%_* to remove
# all characters following _ in the filename from end of string. This way we feed
# the sample name that is common to both L001 and L002 fastqs to spaceranger.
# Also, R2 file is the one that has our sequence not R1 or I1 fastqs.
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

filenames=(${inputs[*]//$FASTQ_DIR/})    # B1246-GEX_S15_L001_R2_001.fastq.gz 
#taskinputs=(${filenames[*]%%_*})         # B1246-GEX 
#taskinputs=($(echo "${taskinputs[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
taskinputs=(${filenames[*]%%_S*[0-9]_L[0-9][0-9][0-9]*.gz})
taskinput=${taskinputs[$index]}
CYTASSIST_IMAGE=${cyt_images[$index]}
HE_IMAGE=${he_images[$index]}

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name   : $JOB_NAME"
echo "Job ID     : $JOB_ID"  
echo "TaskID     : $SGE_TASK_ID"
echo "index      : $index"
echo "filename   : $filename"
echo "taskinput  : $taskinput"
echo "image      : $CYTASSIST_IMAGE"
echo "image      : $HE_IMAGE"
echo "=========================================================="

# You need to add path for Cellranger
PATH=$HOME/NGSTools/spaceranger-3.1.3:$PATH

# Run spaceranger count
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count

# DO NOT use --nosecondary. We can see clusters generated by spaceranger using html file.

#SAMPLES=(FC1 FT1 FT2 MC1 MC2 MT1 MT2 MT3 MT4 MT5 MT6)
#for taskinput in ${SAMPLES[*]}
#do
#echo $taskinput
if [ $index == 0 ]
then
spaceranger count \
--id=$taskinput \
--transcriptome=$GENOME \
--probe-set=$PROBE_SET \
--filter-probes=true \
--fastqs=$FASTQ_DIR \
--slide=H1-M2CQPKR \
--area=A1 \
--cytaimage=$CYTASSIST_IMAGE \
--sample=$taskinput \
--reorient-images=true \
--create-bam=false

elif [ $index == 1 ]
then
spaceranger count \
--id=$taskinput \
--transcriptome=$GENOME \
--probe-set=$PROBE_SET \
--filter-probes=true \
--fastqs=$FASTQ_DIR \
--slide=H1-M2CQPKR \
--area=D1 \
--cytaimage=$CYTASSIST_IMAGE \
--sample=$taskinput \
--reorient-images=true \
--create-bam=false

else
echo "Done"
fi

#--image=$HE_IMAGE \
#--unknown-slide=visium-1 \
#done

# Once spaceranger count is complete, run the code below once to 
# copy the raw_feature_bc_matrix folders for Seurat analysis
#(${array[*]/%/suffix}) to add suffix to add array elements
# (${array[*]/#/prefix}) to add prefix to add array elements

# # For visiumHD, the images in 
# # outs/binned_outputs/square_002um/spatial, 
# # outs/binned_outputs/square_008um/spatial, 
# # outs/binned_outputs/square_016um/spatial are usually empty.
# # Delete them and replace with images from outs/spatial.

# cp $HOME/scratch/$proj/spaceranger_results/*/outs/spatial/* $HOME/scratch/$proj/spaceranger_results/*/outs/binned_outputs/square_002um/spatial/
# cp $HOME/scratch/$proj/spaceranger_results/*/outs/spatial/* $HOME/scratch/$proj/spaceranger_results/*/outs/binned_outputs/square_008um/spatial/
# cp $HOME/scratch/$proj/spaceranger_results/*/outs/spatial/* $HOME/scratch/$proj/spaceranger_results/*/outs/binned_outputs/square_016um/spatial/

# LOCATION=$HOME/scratch/$proj
# SAMPLES=($(ls $LOCATION/spaceranger_results/))
# SAMPLES_WITH_PATH=($(ls $LOCATION/spaceranger_results/* -d))
# #SPATIAL_FOLDER_WITH_PATH=(${SAMPLES_WITH_PATH[*]/%//outs/spatial}) 
# #H5_FILES_WITH_PATH=(${SAMPLES_WITH_PATH[*]/%//outs/raw_feature_bc_matrix.h5}) 
# #SPATIAL_FOLDER_WITH_PATH=(${SAMPLES_WITH_PATH[*]/%//outs/binned_outputs}) 
# NEW_LOCATION=(${SAMPLES[*]/#/$LOCATION/raw_feature_bc_matrix/})

# rm $HOME/scratch/$proj/spaceranger_results/*/outs/binned_outputs/*/spatial/*jpg
# rm $HOME/scratch/$proj/spaceranger_results/*/outs/binned_outputs/*/spatial/*png
# rm $HOME/scratch/$proj/spaceranger_results/*/outs/binned_outputs/*/spatial/*tiff

# for index in ${!SAMPLES[*]}
# do
  
 # cp -r ${SAMPLES_WITH_PATH[$index]}/outs/spatial/* ${SAMPLES_WITH_PATH[$index]}/outs/binned_outputs/square_002um/spatial/
 # cp -r ${SAMPLES_WITH_PATH[$index]}/outs/spatial/* ${SAMPLES_WITH_PATH[$index]}/outs/binned_outputs/square_008um/spatial/
 # cp -r ${SAMPLES_WITH_PATH[$index]}/outs/spatial/* ${SAMPLES_WITH_PATH[$index]}/outs/binned_outputs/square_016um/spatial/
 # cp -r ${SAMPLES_WITH_PATH[$index]}/outs/binned_outputs/* ${NEW_LOCATION[$index]}
 # #cp -r ${H5_FILES_WITH_PATH[$index]} ${NEW_LOCATION[$index]}
 
# done