#!/bin/bash -l

# Set Job Name
#$ -N Seurat
# Request memory. You can find maxvmem from qstat -j <jobid> & adjust mem_free accordingly in final run.
#$ -l mem_free=120G
# Merge standard output and standard error
#$ -j y
# Set current working directory
#$ -cwd 
# Submit an array job with n tasks. Find n using: ls *.fastq.gz | wc -l 
#$ -t 1

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

# SCRIPT="scRNASeq_Seurat_Analysis_Import_Data_new.R"
# SCRIPT="scRNASeq_Seurat_Analysis_HTO_Demux_PhaseI.R"
# SCRIPT="scRNASeq_Seurat_Analysis_HTO_Demux_PhaseII.R"
# SCRIPT="scRNASeq_Seurat_Analysis_Initial_Annotation.R"
# SCRIPT="scRNASeq_Seurat_Analysis_Primary_Cleanup.R"
# SCRIPT="scRNASeq_Seurat_Analysis_Secondary_Cleanup.R"
# SCRIPT="scRNASeq_Seurat_Analysis_Subtype.R"
SCRIPT="Spatial_Seurat_Analysis_Import_Data_new.R"

#PROJ="scRNASeq_BBN_C57BL6"      
#PROJ="scRNASeq_BBN_Rag"
#PROJ="scRNASeq_Chen"
#PROJ="scRNASeq_GSE164557"
#PROJ="scRNASeq_GSE217093"
#PROJ="scRNASeq_GSE222315"
#PROJ="scRNASeq_HRA003620"
#PROJ="scRNASeq_Jinfen"
#PROJ="scRNASeq_NA13_CDH12_C57BL6"
#PROJ="scRNASeq_Jyoti"

# PROJ="scRNASeq_Koltsova"
# PROJ="scRNASeq_Koltsova_sn"
# PROJ="scRNASeq_GSE145216"
# PROJ="visium_GSE171351"
PROJ="Spatial_Bhowmick"
SPECIES="Homo sapiens"  #"Mus musculus"
ASSAY="RNA"   #"Spatial"

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name   : $JOB_NAME"
echo "Job ID     : $JOB_ID"
echo "TaskID     : $SGE_TASK_ID"
echo "Script     : $SCRIPT"
echo "Project	 : $PROJ"
echo "=========================================================="

chmod u+x ~/projects/scRNASeq/$SCRIPT
Rscript ~/projects/scRNASeq/$SCRIPT --args proj=$PROJ species=$SPECIES assay=$ASSAY
