#!/bin/bash -l

# Set Project Name
# Set Job Name
#$ -N CRISPR_MAGeCK_mle
# Request memory. You can find maxvmem from qstat -j <jobid> & adjust mem_free accordingly in final run.
#$ -l mem_free=8G
# Merge standard output and standard error
#$ -j y
# Set current working directory
#$ -cwd
# Submit an array job with n tasks where n = 1 as we are running one by one
#$ -t 1

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID "
echo "=========================================================="

# You need to add path for Python3 for MAGeCK to run
PATH=$HOME/Python3/bin/:$PATH

cd $HOME/scratch/CRISPR/MAGeCK_mle_results

#--count-table: 	output of mageck count
#--design-matrix:	Provide a design matrix. The row of the design matrix must match the order of the samples in the count table.
#--norm-method: 	Use "total" or "control". If using "control", use SAFE sgRNAs as "control", not Non-targeting sgRNAs
#--output-prefix:	The prefix to be added o output of mageck test
#--day0-label:		Use this ONLY if you dont specify --control-id. DO NOT use both parameters simultaneously.
# 					If you want to compare more than 2 sets of treated samples (say Day 7 and Day 14) to the same set of control samples (Day 0), then 
# 					use this parameter. If you want to compare only 2 sets of samples: "Day 7 vs Day 0"/"Day 14 vs Day 0"/"Day 7 vs Day 14", use --control id. 
# 					Turns on negative selection QC. Usually day 0 or plasmid sample is given this label. For every other sample label,
#					the negative selection QC will compare it with day0 sample, and estimate the degree of negative selections in essential genes

#--control-sgrna: 	Only use if --norm-method is control. Use SAFE sgRNAs
#--control-gene: 	Only use if --norm-method is control. Use genes corresponding to SAFE sgRNAs
#--cnv-norm: 		A matrix of copy number variation data across cell lines to normalize CNV-biased sgRNA scores prior to gene ranking
#--cell-line:		Cell line to be used for cnv normalization. Must match one of the column names in the file provided by --cnv-norm


#--count-table $HOME/scratch/CRISPR/MAGeCK_count_results/MAGeCK_count.count.txt 
mageck mle \
--count-table $HOME/scratch/CRISPR/count_results/MAGeCK_kms.count.txt \
--design-matrix $HOME/projects/CRISPR/Jinfen_design_matrix_1.txt \
--control-sgrna $HOME/projects/CRISPR/control_sgRNAs.txt \
--norm-method control \
--output-prefix MAGeCK_mle_1
#--day0-label \
#--control-gene ?? \
#--cnv-norm ?? \
#--cell-line ?? \
#--cnv-est ?? \
#--genes-varmodeling ?? \
#--permutation-round ?? \
#--no-permutation-by-group ?? \
#--max-sgrnapergene-permutation ?? \
#--remove-outliers \
#--threads ?? \
#--adjust-method {fdr,holm,pounds} \
#--sgrna-efficiency ?? \
#--sgrna-eff-name-column ?? \
#--sgrna-eff-score-column ?? \
#--update-efficiency \
#--bayes \
#--PPI-prior \
#--PPI-weighting ?? \
#--negative-control ??

#--count-table $HOME/scratch/CRISPR/MAGeCK_count_results/MAGeCK_count.count.txt
mageck mle \
--count-table $HOME/scratch/CRISPR/count_results/MAGeCK_kms.count.txt \
--design-matrix $HOME/projects/CRISPR/Jinfen_design_matrix_2.txt \
--control-sgrna $HOME/projects/CRISPR/control_sgRNAs.txt \
--norm-method control \
--output-prefix MAGeCK_mle_2