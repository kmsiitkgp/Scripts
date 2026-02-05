#!/usr/bin/env python3

# from pyscenic.binarization import binarize
# import pandas as pd
# import numpy as np

# # Read csv file
# auc_mat = pd.read_csv(filepath_or_buffer = "/hpc/home/kailasamms/scratch/proj_Fibroblast/pySCENIC_results/pySCENIC_AUC.csv")

# # Check that this matches dim(AUC_mat) in R before creating csv file
# print(auc_mat.shape)              
# print(auc_mat.head())
# print(type(auc_mat))

# # Calculate binary AUC scores
# auc_mat_binary = binarize(auc_mat)

# # Save the results as csv file
# auc_mat_binary.to_csv(path_or_buf = "/hpc/home/kailasamms/scratch/proj_Fibroblast/pySCENIC_results/pySCENIC_AUC_binary1.csv", index=False)

# #test = pd.DataFrame(auc_mtx_binary[0])

import pandas as pd
import sys
import os

# Declare variables
cpdb_version = "v4.1.0"                                        # Version of CPDB databse
target_dir = "/hpc/home/kailasamms/projects/scRNASeq/CPDB"     # Directory to store all versions of CPDB database
cpdb_target_dir = os.path.join(target_dir, cpdb_version)           # Directory to store current version of CPDB database

# Download the database
from cellphonedb.utils import db_utils
db_utils.download_database(cpdb_target_dir, cpdb_version)

# Declare variables
home_dir = "/hpc/home/kailasamms/scratch/"                        # Directory to store all results
proj = "scRNASeq_Chen" #"scRNASeq_BBN_C57B6"                      # Folder containing input files
cpdb_file_path = os.path.join(cpdb_target_dir,"cellphonedb.zip")  # CPDB database files
counts_file_path = os.path.join(home_dir, proj, "results_cellphonedb", "Male")
meta_file_path = os.path.join(home_dir, proj, "results_cellphonedb", "Male_meta.tsv")
microenvs_file_path = os.path.join(home_dir, proj, "results_cellphonedb", "Male_micro.tsv")
out_path = os.path.join(home_dir, proj, "results_cellphonedb")

def create_argument_parser():
    parser = argparse.ArgumentParser(
        description="Run Arboreto using a multiprocessing pool"
    )

    parser.add_argument(
        "expression_mtx_fname",
        type=str,
        help="The name of the file that contains the expression matrix for the single cell experiment."
        " Two file formats are supported: csv (rows=cells x columns=genes) or loom (rows=genes x columns=cells).",
    )
    parser.add_argument(
        "tfs_fname",
        type=str,
        help="The name of the file that contains the list of transcription factors (TXT; one TF per line).",
    )
    parser.add_argument(
        "-m",
        "--method",
        choices=["genie3", "grnboost2"],
        default="grnboost2",
        help="The algorithm for gene regulatory network reconstruction (default: grnboost2).",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=sys.stdout,
        help="Output file/stream, i.e. a table of TF-target genes (TSV).",
    )
    
    return parser




def main():
    parser = create_argument_parser()
    args = parser.parse_args()



# Run Statistical method
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

deconvoluted, means, pvalues, significant_means = cpdb_statistical_analysis_method.call(
    cpdb_file_path = cpdb_file_path,                 # mandatory: CellPhoneDB database zip file.
    meta_file_path = meta_file_path,                 # mandatory: tsv file defining barcodes to cell label.
    counts_file_path = counts_file_path,             # mandatory: normalized count matrix.
    counts_data = 'gene_name',                       # defines the gene annotation in counts matrix.
    microenvs_file_path = microenvs_file_path,       # optional (default: None): defines cells per microenvironment.
    iterations = 1000,                               # denotes the number of shufflings performed in the analysis.
    threshold = 0.1,                                 # defines the min % of cells expressing a gene for this to be employed in the analysis.
    threads = 4,                                     # number of threads to use in the analysis.
    debug_seed = 42,                                 # debug randome seed. To disable >=0.
    result_precision = 3,                            # Sets the rounding for the mean values in significan_means.
    pvalue = 0.05,                                   # P-value threshold to employ for significance.
    subsampling = False,                             # To enable subsampling the data (geometri sketching).
    subsampling_log = False,                         # (mandatory) enable subsampling log1p for non log-transformed data inputs.
    subsampling_num_pc = 100,                        # Number of componets to subsample via geometric skectching (dafault: 100).
    subsampling_num_cells = 1000,                    # Number of cells to subsample (integer) (default: 1/3 of the dataset).
    separator = '|',                                 # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    debug = False,                                   # Saves all intermediate tables employed during the analysis in pkl format.
    output_path = out_path,                          # Path to save results.
    output_suffix = None                             # Replaces the timestamp in the output files by a user defined string in the  (default: None).
    )
