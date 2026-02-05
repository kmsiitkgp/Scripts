#!/bin/bash -l

# Set Job Name
#$ -N Rscript       	# Project name 
#$ -l mem_free=120G     # Request memory
#$ -j y                 # Merge standard output and standard error
#$ -cwd                 # Set current working directory
#$ -t 1                 # Submit an array job with 1 task 

# Check R verrion available in cluster modules using "module avail".
# It is NOT RECOMMENDED to use R version available in cluster modules as:
# (i) it may not be latest version of R
# (ii) you need to declare location of R packages
#module load R/4.3.1
#export R_LIBS="/hpc/home/kailasamms/NGSTools/R_packages/4.3.1/"

# It is RECOMMENDED to install R using conda.
# To use R from conda, you need to add path for Conda and 
# use "source activate <env_name>" instead of "conda activate <env_name>"
PATH=$HOME/miniconda3/bin/:$PATH
source activate R

# Know the location from where Rscript is running
which Rscript
#SCRIPT=~/projects/TCGA_GDC/access_TCGA.R
SCRIPT=~/saswat.R

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name   : $JOB_NAME"
echo "Job ID     : $JOB_ID"
echo "TaskID     : $SGE_TASK_ID"
echo "Project	 : $PROJ"
echo "=========================================================="

chmod u+x $SCRIPT
Rscript $SCRIPT
