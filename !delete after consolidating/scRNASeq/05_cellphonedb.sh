#!/bin/bash -l

#$ -N cellphonedb          # Set Job Name
#$ -l mem_free=120G        # Request memory [>100G recommended]
#$ -pe orte 18             # Request n cores [>12 cores recommended]
#$ -j y                    # Merge standard output and standard error 
#$ -cwd                    # Set current working directory
#$ -t 1                    # Submit an array job with n tasks. 

proj=scRNASeq_BBN_C57B6_paper
META_FILE=$HOME/scratch/$proj/cellphonedb_results/Male_meta.tsv
MTX_FOLDER=$HOME/scratch/$proj/cellphonedb_results/Male
MTX_FILE=$HOME/scratch/$proj/cellphonedb_results/Male/matrix.mtx
MICRO_FILE=$HOME/scratch/$proj/cellphonedb_results/Male_microenv.tsv
DEG_FILE=$HOME/scratch/$proj/cellphonedb_results/Male_DEGs.tsv
OUTPUT_DIR=$HOME/scratch/$proj/cellphonedb_results
PROJ_NAME=Male_results   # cellphonedb will create a folder in output directory

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date         : $(date)"
echo "Job name           : $JOB_NAME"
echo "Job ID             : $JOB_ID"  
echo "TaskID             : $SGE_TASK_ID"
echo "Meta file          : $META_FILE"
echo "MTX folder         : $MTX_FOLDER"
echo "DEG file           : $DEG_FILE"
echo "=========================================================="

# You need to add path for Conda
PATH=$HOME/miniconda3/bin/:$PATH
# You have to use source activate instead of conda activate
source activate cpdb
# Check which version of Python is being used
python --version
echo "PATH ="; echo $PATH
echo "PYTHONPATH ="; echo $PYTHONPATH

# Latest database is automatically downloaded
# This is version 3 command. Works with python 3.7 only
# cellphonedb method statistical_analysis  \
# $META_FILE \
# $MTX_FOLDER \
# --counts-data gene_name \
# --project-name $PROJ_NAME \
# --output-path $OUTPUT_DIR \
# --verbose

# This is version 3 command. Works with python 3.7 only
# cellphonedb method degs_analysis  \
# $META_FILE \
# $MTX_FOLDER \
# $DEG_FILE \
# --counts-data gene_name \
# --project-name $PROJ_NAME \
# --output-path $OUTPUT_DIR \
# --verbose

date
echo "Interactions Inferred"

# If using statistical method:
MEANS_FILE=$OUTPUT_DIR/$PROJ_NAME/means.txt
PVAL_FILE=$OUTPUT_DIR/$PROJ_NAME/pvalues.txt

# If using degs method:
MEANS_FILE=$OUTPUT_DIR/$PROJ_NAME/significant_means.txt
PVAL_FILE=$OUTPUT_DIR/$PROJ_NAME/relevant_interactions.txt

# Keep track of information related to the current job
echo "=========================================================="
echo "Mean file           : $MEANS_FILE"
echo "pval file           : $PVAL_FILE"
echo "=========================================================="

# cellphonedb plot dot_plot \
# --means-path $MEANS_FILE \
# --pvalues-path $PVAL_FILE \
# --output-path $OUTPUT_DIR \
# --verbose



