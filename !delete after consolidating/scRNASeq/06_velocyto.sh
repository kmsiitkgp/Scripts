#!/bin/bash -l

#$ -N velocyto         # Set Job Name
#$ -l mem_free=64G     # Request memory.
#$ -j y                # Merge standard output and standard error
#$ -cwd                # Set current working directory
#$ -t 1-12             # Submit an array job with n tasks
# NOTE: n = number of folders; each folder corresponds to one sample

# NOTE: Just like cellranger count which creates an individual folder for each sample,
# velocyto also stores its results as loom file within these individual folders.
proj=scRNASeq_BBN_C57B6

GTF_FILE=$HOME/NGSTools/refdata-gex-mm10-2020-A/genes/genes.gtf
MASK_FILE=$HOME/NGSTools/mm10_rmsk.gtf
SAMPLE_DIR=$HOME/scratch/$proj/cellranger_results

# Create an array of cellranger output folders
input=($SAMPLE_DIR/*)

# Verify if the array stores correctly
# echo ${input[0]}

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
echo "GTF file   : $GTF_FILE"
echo "MASK file  : $MASK_FILE"
echo "taskinput  : $taskinput"
echo "=========================================================="

# You need to add path for Conda
PATH=$HOME/miniconda3/bin/:$PATH
# You have to use source activate instead of conda activate
source activate scVelo
# Check which version of velocyto is being used
# velocyto is a python package. We can use conda to invoke velocyto
which velocyto
velocyto --version

# Using \ enables us to split the long command into easy to read format.
# Masking gtf is needed to mask expressed repetitive elements, since those 
# count could constitute a confounding factor in the downstream analysis. 
# Download the masking file from https://genome.ucsc.edu/cgi-bin/hgTables
# clade:Mammal; genome:Mouse; assembly:mm10; Track:RepeatMasker; Table:rmsk; output format:GTF
# (Although mm39 is latest, the assembly used as input for cellranger count is mm10. 
# So, we use mm10 for generating the masking file)
# Press “get output” to download and upload to HPC cluster using cyberduck and gunzip

velocyto run10x \
--mask $MASK_FILE \
$taskinput \
$GTF_FILE

# # Once velocyto run10x is complete, run the code below once to 
# # copy the loom files for scvelo analysis
# # (${array[*]/%/suffix}) to add suffix to add array elements
# # (${array[*]/#/prefix}) to add prefix to add array elements
# LOCATION=$HOME/scratch/$proj
# SAMPLES=($(ls $LOCATION/cellranger_results/))
# SAMPLES_WITH_PATH=($(ls $LOCATION/cellranger_results/* -d))
# CURRENT_SAMPLES_WITH_PATH=(${SAMPLES_WITH_PATH[*]/%//velocyto/*})

# t=($LOCATION/cellranger_results/${SAMPLES[$index]}/velocyto/*)
# for index in ${!SAMPLES[*]}
# do
  # NEW_SAMPLES=$(basename ${CURRENT_SAMPLES_WITH_PATH[$index]})
  # NEW_SAMPLES_WITH_PATH=${NEW_SAMPLES[*]/#/$LOCATION/velocyto_results/}
  # cp ${CURRENT_SAMPLES_WITH_PATH[$index]} $NEW_SAMPLES_WITH_PATH
# done