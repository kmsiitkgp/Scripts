#!/bin/bash -l

#$ -N DESeq2                    # Set Job Name
#$ -l mem_free=120G             # Request memory
#$ -j y                         # Merge standard output and standard error
#$ -cwd                         # Set current working directory
#$ -t 1-1                       # Submit an array job with 1 task.

# Check R verrsion available in cluster modules using "module avail".
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
module load hdf5/1.8.18

# Know the location from where Rscript is running
which Rscript

# (RELEVANT ONLY IF USING R FROM CLUSTER MODULE)
# If the cluster doesnt have a package you need to use, you have to install
# it in your local directory as you dont have root access to install in the clsuter
# Some packages that can be installed from within R scipt that has been submitted as a job 
# Eg: install.packages("viridis") works fine
# Others cannot be installed from within R scipt that has been submitted as a job
# Eg: install.packages("multtest") will fail because these packages are written in C/C++/Fortran
# and appropriate modules need to be loaded at the time of loading R.
# Such packages need to be manually installed before submitting a job. 
# First, create a variable in UNIX to store path to install package i.e. library location 
# export R_LIBS="/hpc/home/kailasamms/NGSTools/R_packages" 
# module load R/4.3.1
# R
# >install.packages("hdf5r", configure.args="--with-hdf5=/hpc/apps/hdf5/1.8.18/bin/h5cc")  #locate where h5cc manually
# >install.packages("harmony")

# (GENERAL NOTE)
# The cluster has an R library where Biobase, etc are already available. 
# So, we have to use "force=TRUE", else BiocManager will not install to our local library. 
# Also, set "lib=<library location>", else BiocManager will not update dependencies of Biobase etc.
# >BiocManager::install("DESeq2", lib = "/hpc/home/kailasamms/NGSTools/R_packages", force = TRUE)
# >install.packages("metap")		# You need to install multtest before metap
# >q()								# Quit R and qsub your bash script to submit the job 


SCRIPT="access_TCGA.R"

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name   : $JOB_NAME"
echo "Job ID     : $JOB_ID"
echo "TaskID     : $SGE_TASK_ID"
echo "Script     : $SCRIPT"
echo "Project	 : $PROJ"
echo "=========================================================="

chmod u+x $SCRIPT
Rscript $SCRIPT
