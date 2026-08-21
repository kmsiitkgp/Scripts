#!/bin/bash
#SBATCH --job-name=tcga_dds
#SBATCH --time=48:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --output=tcga_dds_%j.log

#===============================================================================
# Environment Initialization
#===============================================================================

# 1. Source the absolute setup path from your miniconda3 folder
source /home/kailasamms/miniconda3/etc/profile.d/conda.sh

# 2. Activate your specific environment
conda activate omics

# 3. Run your R script
Rscript /home/kailasamms/scripts/nextflow/projects/TCGA/03_create_tcga_dds.R