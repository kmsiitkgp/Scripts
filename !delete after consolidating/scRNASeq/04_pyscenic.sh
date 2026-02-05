#!/bin/bash -l

#$ -N pySCENIC             # Set Job Name
#$ -l mem_free=120G        # Request memory [>100G recommended]
#$ -pe orte 18             # Request n cores [>12 cores recommended]
#$ -j y                    # Merge standard output and standard error 
#$ -cwd                    # Set current working directory
#$ -t 1-4                  # Submit an array job with n tasks. 

# # Create results directory
# res_dir=outs_$(date +"%Y-%m-%d")
# mkdir -p $res_dir
# date

# NOTE: If the loom file is named "Epithelial Cells.loom", then arboreto_with_multiprocessing
# will throw error due to the blank space in loom filename. So, rename loom files such
# that there are no blank spaces in the filenames

SPECIES=Mouse
#SPECIES=Human
TF_FILE=$HOME/projects/scRNASeq/SCENIC/TFs/$SPECIES/*.txt
ANNO_FILE=$HOME/projects/scRNASeq/SCENIC/Annotations/$SPECIES/*.tbl
MOTIF_FILE=($HOME/projects/scRNASeq/SCENIC/Motifs/$SPECIES/*.feather)

proj=scRNASeq_BBN_C57B6
OUTPUT_DIR=$HOME/scratch/$proj/results_pyscenic
LOOM_FILES=($HOME/scratch/$proj/results_pyscenic/*.loom)

# Use the SGE_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value
index=$(($SGE_TASK_ID-1))

SAMPLE=$(basename ${LOOM_FILES[$index]})
SAMPLE=(${SAMPLE[*]//pySCENIC_input_/})
SAMPLE=(${SAMPLE[*]//.loom/})

cd $OUTPUT_DIR
LOOM_FILE=${LOOM_FILES[$index]}
GRN_OUTPUT_FILE=pySCENIC_GRN_adjacencies_$SAMPLE.csv
CTX_OUTPUT_FILE=pySCENIC_CTX_regulons_$SAMPLE.csv
AUCELL_OUTPUT_FILE=pySCENIC_AUCELL_filtered_$SAMPLE.loom
LOOM_VIZ_FILE=pySCENIC_AUCELL_filtered_viz_$SAMPLE.loom

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date         : $(date)"
echo "Job name           : $JOB_NAME"
echo "Job ID             : $JOB_ID"  
echo "TaskID             : $SGE_TASK_ID"
echo "TF file            : " $TF_FILE
echo "Annotation file    : " $ANNO_FILE
echo "Motif file#1       : ${MOTIF_FILE[0]}"
echo "Motif file#2       : ${MOTIF_FILE[1]}"
echo "Output directory   : $OUTPUT_DIR"
echo "Loom file          : $LOOM_FILE"
echo "GRN output file    : $GRN_OUTPUT_FILE"
echo "CTX output file    : $CTX_OUTPUT_FILE"
echo "AUCELL output file : $AUCELL_OUTPUT_FILE"
echo "VIZ file           : $LOOM_VIZ_FILE"
echo "=========================================================="

# You need to add path for Conda
PATH=$HOME/miniconda3/bin/:$PATH
# You have to use source activate instead of conda activate
source activate pySCENIC
# Check which version of Python is being used
python --version
echo "PATH ="; echo $PATH
echo "PYTHONPATH ="; echo $PYTHONPATH

# # Step 1: Get regulatory networks
# # NOTE: If you type "pyscenic grn --help", you will see pyscenic grn uses "expression_mtx_fname"  
# # and "tfs_fname" as positional arguments.  So, first specify "expression_mtx_fname" which is 
# # either csv or loom file and then specify "tfs_fname" which is list of TFs in txt file. 
# # For small datasets (less than 2k cells x 10k genes), pyscenic grn will complete without errors.
# # For large datasets (more than 10k cells x 10k genes), dask scheduler will give errors.
# pyscenic grn \
# $LOOM_FILE \
# $TF_FILE \
# --method grnboost2 \
# --seed 787878 \
# --num_workers 18 \
# --output $GRN_OUTPUT_FILE 

# pyscenic grn uses dask scheduler which seems to cause trouble with large datasets.
# To avoid this use "arboreto_with_multiprocessing.py" which circumvents the use of dask. 
# https://pyscenic.readthedocs.io/en/latest/faq.html
# The "arboreto_with_multiprocessing.py" scrip is available on the path once pySCENIC is installed. 
arboreto_with_multiprocessing.py \
$LOOM_FILE \
$TF_FILE \
--method grnboost2 \
--seed 787878 \
--num_workers 18 \
--output $GRN_OUTPUT_FILE

date
echo "Regulatory Networks Inferred"

# Step 2: Define modules of TF-regulons
# Directly continuing from the same slurm script as above, run pyscenic ctx. 
# This step cleans up the identified co-expression networks by drawing information from the
# feather databases and motif annotation files to identify target genes with significantly enriched motifs.
# NOTE: If you type "pyscenic ctx --help", you will see pyscenic ctx uses "module_fname"  and "database_fname"
# as positional arguments.  So, first specify "module_fname" which is output of pyscenic grn (STEP 1) and
# then specify "database_fname" which are feather databases
pyscenic ctx \
$GRN_OUTPUT_FILE \
${MOTIF_FILE[0]} \
${MOTIF_FILE[1]} \
--annotations_fname $ANNO_FILE \
--expression_mtx_fname $LOOM_FILE \
--output $CTX_OUTPUT_FILE \
--mask_dropouts \
--num_workers 18

date
echo "Regulon Modules Defined"

# Step 3:  Calculate regulon modules score for each cell
# NOTE: If you type "pyscenic aucell --help", you will see pyscenic aucell uses "expression_mtx_fname" 
# and "signatures_fname" as positional arguments.  So, first specify "expression_mtx_fname" which
# is loom file and then specify "signatures_fname" which is output of pyscenic ctx (gmt, yaml or dat files)
# pyscenic aucell also seems to accept csv file for "signatures_fname"
pyscenic aucell \
$LOOM_FILE \
$CTX_OUTPUT_FILE \
--num_workers 18 \
--seed 787878 \
--output $AUCELL_OUTPUT_FILE

date   
echo "Cell Scoring Completed"

# Add dimensionality reduction (derived from pySCENIC AUCell matrix)
# I always skipped this step, as it can also be performed in R (with proper sample integration)
python $HOME/projects/scRNASeq/pySCENIC_add_visualization.py \
   --loom_input $AUCELL_OUTPUT_FILE \
   --loom_output $LOOM_VIZ_FILE \
   --num_workers 20

date
echo "Dimensionality Reduction Added"

# Binarize the AUC scores (do in R)
# python $HOME/projects/scRNASeq/pySCENIC_binarize.py