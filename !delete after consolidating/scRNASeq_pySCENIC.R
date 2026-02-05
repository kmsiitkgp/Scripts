#******************************************************************************#
# There are several ways of performing SCENIC analysis
# (i) Using R for pre-SCENIC, SCENIC and post-SCENIC analysis (NOT RECOMMENDED as GENIE3 takes weeks)
# (ii) Using python for pre-SCENIC, SCENIC and post-SCENIC analysis (NOT RECOMMENDED as it is too complex)
# (iii) Using R for pre- and post-SCENIC analysis but pySCENIC for SCENIC analysis (RECOMMENDED)

# For (iii), we do the following
# We perform all initial filtering as needed by SCENIC and create a loom file.
# Then, we use pySCENIC to perform SCENIC analysis in python.
# Finally, we import the output into R and complete the post-SCENIC analysis.

# # Read these links to learn how to perform SCENIC analysis
# # https://github.com/aertslab/SCENIC
# # https://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Setup.html
# # https://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html
# # https://github.com/aertslab/SCENIC/issues/364
# # https://raw.githubusercontent.com/aertslab/SCENIC/master/vignettes/SCENIC_Running.Rmd
# # https://pyscenic.readthedocs.io/en/latest/faq.html

#*************STEP 1: DOWNLOAD RELEVANT FILES FOR pySCENIC ANALYSIS************#

# Download species-specific databases for RcisTarget (the motif rankings) from
# https://resources.aertslab.org/cistarget/
# NOTE: v2 version is NEW and v1 version is OLD
# NOTE: USE v2 version for Python (RECOMMENDED) and v1 version for R (NOT RECOMMENDED).
# NOTE: Use gene based for SCENIC and region based for SCENIC+
# Run the commands below outside of R in terminal

# Download the list of transcription factors (TFs)
tf_fly_dir="$HOME/projects/scRNASeq/SCENIC/TFs/Fly/"
tf_fly_url="https://resources.aertslab.org/cistarget/tf_lists/allTFs_dmel.txt"
wget -P $tf_fly_dir $tf_fly_url --no-check-certificate

tf_human_dir="$HOME/projects/scRNASeq/SCENIC/TFs/Human/"
tf_human_url="https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt"
wget -P $tf_human_dir $tf_human_url --no-check-certificate

tf_mouse_dir="$HOME/projects/scRNASeq/SCENIC/TFs/Mouse/"
tf_mouse_url="https://resources.aertslab.org/cistarget/tf_lists/allTFs_mm.txt"
wget -P $tf_mouse_dir $tf_mouse_url --no-check-certificate

# Download the motif annotations for TFs
anno_fly_dir="$HOME/projects/scRNASeq/SCENIC/Annotations/Fly/"
anno_fly_url="https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl"
wget -P $anno_fly_dir $anno_fly_url --no-check-certificate

anno_human_dir="$HOME/projects/scRNASeq/SCENIC/Annotations/Human/"
anno_human_url="https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
wget -P $anno_human_dir $anno_human_url --no-check-certificate

anno_mouse_dir="$HOME/projects/scRNASeq/SCENIC/Annotations/Mouse/"
anno_mouse_url="https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"
wget -P $anno_mouse_dir $anno_mouse_url --no-check-certificate

# Download motif enrichment databases
motif_fly_dir="$HOME/projects/scRNASeq/SCENIC/Motifs/Fly/"
motif_fly_url="https://resources.aertslab.org/cistarget/databases/drosophila_melanogaster/dm6/flybase_r6.02/mc_v10_clust/gene_based/dm6_v10_clust.genes_vs_motifs.rankings.feather"
motif_fly_checksum="https://resources.aertslab.org/cistarget/databases/drosophila_melanogaster/dm6/flybase_r6.02/mc_v10_clust/gene_based/dm6_v10_clust.genes_vs_motifs.rankings.feather.sha1sum.txt"
wget -P $motif_fly_dir $motif_fly_url --no-check-certificate
wget -P $motif_fly_dir $motif_fly_checksum --no-check-certificate
cd $motif_fly_dir
sha1sum -c *sha1sum.txt

motif_human_dir="$HOME/projects/scRNASeq/SCENIC/Motifs/Human/"
motif_human_url1="https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
motif_human_url2="https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
motif_human_checksum1="https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather.sha1sum.txt"
motif_human_checksum2="https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather.sha1sum.txt"
wget -P $motif_human_dir $motif_human_url1 --no-check-certificate
wget -P $motif_human_dir $motif_human_url2 --no-check-certificate
wget -P $motif_human_dir $motif_human_checksum1 --no-check-certificate
wget -P $motif_human_dir $motif_human_checksum2 --no-check-certificate
cd $motif_human_dir
sha1sum -c *sha1sum.txt

motif_mouse_dir="$HOME/projects/scRNASeq/SCENIC/Motifs/Mouse/"
motif_mouse_url1="https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
motif_mouse_url2="https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
motif_mouse_checksum1="https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather.sha1sum.txt"
motif_mouse_checksum2="https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather.sha1sum.txt"
wget -P $motif_mouse_dir $motif_mouse_url1 --no-check-certificate
wget -P $motif_mouse_dir $motif_mouse_url2 --no-check-certificate
wget -P $motif_mouse_dir $motif_mouse_checksum1 --no-check-certificate
wget -P $motif_mouse_dir $motif_mouse_checksum2 --no-check-certificate
cd $motif_mouse_dir
sha1sum -c *sha1sum.txt

#********************STEP 2: PREPARE LOOM FILE FOR pySCENIC********************#

for (celltype in c("Epithelial", "Fibroblasts", "Myeloid", "Lymphoid")){
  
  # Load the integrated seurat object
  integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn", 
                                      dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  integrated_seurat <- subset(x = integrated_seurat,
                              cell_class %in% c("Mixed", "Unclassified"),
                              invert = TRUE)
  
  # Extract raw count matrix from RNA assay [can also use normalized]
  expr_mat <- integrated_seurat@assays$RNA@counts
  print(dim(expr_mat))
  
  # Extract cell-level metadata
  cell_info <- integrated_seurat@meta.data
  print(head(cell_info))
  
  # Filter out genes expressed at low levels. 
  # These genes can be identified by adding their UMIs across all cells.
  umi_cutoff <- 3*0.01*ncol(expr_mat)
  genes_kept1 <- rownames(expr_mat[rowSums(expr_mat) > umi_cutoff,])
  
  # Filter out genes expressed in less than 1% of cells
  min_cells_cutoff <- 0.01*ncol(expr_mat)
  genes_kept2 <- data.frame("n_cells" = rowSums(expr_mat > 0)) %>%
    dplyr::filter(n_cells > min_cells_cutoff) %>%
    rownames(.)
  
  # Filter out genes absent in RcisTarget database
  if (species == "Mus musculus"){
    db1 <- arrow::read_feather("/hpc/home/kailasamms/projects/scRNASeq/SCENIC/Motifs/Mouse/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
    db2 <- arrow::read_feather("/hpc/home/kailasamms/projects/scRNASeq/SCENIC/Motifs/Mouse/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
  } else if (species == "Homo sapiens"){
    db1 <- arrow::read_feather("/hpc/home/kailasamms/projects/scRNASeq/SCENIC/Motifs/Human/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
    db2 <- arrow::read_feather("/hpc/home/kailasamms/projects/scRNASeq/SCENIC/Motifs/Human/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
  } else{
    print("Check species parameter. Only Mus musculus or Homo sapiens are valid values")
  }
  
  rcistarget_genes <- base::union(colnames(db1), colnames(db2))
  genes_kept3 <- base::intersect(rownames(expr_mat), rcistarget_genes)
  
  # # Filter out mitochondrial genes, ribosomal genes
  # # NOTE: Some mitochondrial and ribosomal genes are present in Rcistarget 
  # # database. So, it is NOT RECOMMENDED to remove these genes.
  # genes_kept4 <- rownames(expr_mat)[!grepl(pattern = dplyr::if_else(species=="Homo sapiens", 
  #                                                                   "^RP[SL]|^MT-", 
  #                                                                   "^Rp[sl]|^mt-"),
  #                                          x = rownames(expr_mat))]
  
  # Prepare filtered expression matrix using only common genes
  common_genes <- base::intersect(genes_kept1, genes_kept2)
  common_genes <- base::intersect(common_genes, genes_kept3)
  #common_genes <- base::intersect(common_genes, genes_kept4)
  filt_expr_mat <- expr_mat[common_genes,]
  print(dim(filt_expr_mat))
  
  # Extract UMAP embeddings
  umap.reduc <- Seurat::Embeddings(object = integrated_seurat, reduction = "umap")
  
  # Save the data as loom file
  SCopeLoomR::build_loom(file.name = paste0(pyscenic_results, "pySCENIC_input_", celltype, ".loom"),
                         dgem = filt_expr_mat,
                         default.embedding = umap.reduc,
                         default.embedding.name = "UMAP",
                         title = "Filtered RNA assay for pySCENIC",
                         genome = species)
  
  SCopeLoomR::close_loom(loom = loom)
  
  # ### OPTIONAL (necessary if using pySCENIC for dimensional analyses)
  # # Refer https://github.com/hbc/knowledgebase/blob/master/scrnaseq/pySCENIC.md
  # 
  # # Add further information to loom file (e.g. other embeddings)
  # loom <- SCopeLoomR::open_loom(file.path = paste0(pyscenic_results, "pySCENIC_input_", celltype, ".loom"),
  #                               mode = "r+")
  # # Add PCA embedding
  # pca.reduc <- Seurat::Embeddings(integrated_seurat, reduction = "pca")
  # SCopeLoomR::add_embedding(loom = loom, embedding = pca.reduc, name = "PCA")
  # 
  # # Add Seurat cluster definition (DOESNT WORK CURRENTLY)
  # SCopeLoomR::add_seurat_clustering(loom = loom, seurat = integrated_seurat,
  #                                   seurat.clustering.prefix = "integrated_snn_res.",
  #                                   default.clustering.resolution = "1.4")
  # 
  # # Save and close the given .loom file handler.
  # SCopeLoomR::finalize(loom = loom) 
}

#***************STEP 3: PREPARE ENVIRONMENT FOR pySCENIC ANALYSIS**************#

# Create new conda env with latest python and install necessary packages
# NOTE: Sometimes pyscenic may not work with latest python version.
# Run the commands below outside of R in terminal

# conda search python
conda create --name pySCENIC python=3.9
conda activate pySCENIC
pip install pyscenic
conda install -c conda-forge multicore-tsne

# Download export_to_loom.py and add.visualization.py
wget https://raw.githubusercontent.com/vib-singlecell-nf/vsn-pipelines/master/src/scenic/bin/add_visualization.py
wget https://raw.githubusercontent.com/vib-singlecell-nf/vsn-pipelines/master/src/scenic/bin/export_to_loom.py

# Rename export_to_loom.py as pySCENIC_export_to_loom.py and 
# replace "from pyscenic.genesig import Regulon, GeneSignature" with
from ctxcore.genesig import Regulon, GeneSignature
# Similarly, replace "from sklearn.manifold.t_sne import TSNE" with
from sklearn.manifold._t_sne import TSNE

# Rename add.visualization.py as pySCENIC_add.visualization.py and
# replace "import export_to_loom" with
import sys 
import os
sys.path.append(os.path.abspath("/hpc/home/kailasamms/projects/scRNASeq"))
import pySCENIC_export_to_loom
# Also, change line 82
scope_loom = pySCENIC_export_to_loom.SCopeLoom.read_loom(filename=args.loom_input)

#*****************STEP 4: RUN pySCENIC ANALYSIS IN HPC CLUSTER*****************#

# pySCENIC analysis consists of running 
# (i) pyscenic grn (alternative is arboreto_with_multiprocessing.py)
# (ii) pyscenic ctx 
# (iii) pyscenic aucell
# (iv) binarize the AUC scores (do in R as it is simple)

# Since pyscenic grn uses dask scheduler by default which gives error, we 
# use "arboreto_with_multiprocessing.py" instead of pyscenic grn.
# Assign "SPECIES" & "proj" in the sh file and run
# "qsub $HOME/projects/scRNASeq/02d_pyscenic.sh" for steps (i), (ii), (iii).

#******************************************************************************#

# Read below to understand how to access loom file info using SCopeLoomR package

# Each loom file has 7 parts to it: 
# (i) attrs       : contains global attributes
# (ii) matrix     : contains raw expression data
# (iii) layers    : contains alternative representations of the data in matrix
# (iv) row_attrs  : contains gene metadata
# (v) col_attrs   : contains cell metadata
# (vi) row_graphs : contains gene-based cluster graphs
# (vii) col_graphs: contains cell-based cluster graphs 

#*********View the 7 parts of loom file simply using loom file handler*********#

# > loom
# Class: H5File
# Filename: C:\Users\KailasammS\Desktop\pySCENIC_AUCELL_filtered.loom
# Access type: H5F_ACC_RDONLY
# Attributes: last_modified
# Listing:
#         name    obj_type  dataset.dims dataset.type_class
#        attrs   H5I_GROUP          <NA>               <NA>
#    col_attrs   H5I_GROUP          <NA>               <NA>
#   col_graphs   H5I_GROUP          <NA>               <NA>
#       layers   H5I_GROUP          <NA>               <NA>
#       matrix H5I_DATASET 32380 x 11364          H5T_FLOAT
#    row_attrs   H5I_GROUP          <NA>               <NA>
#   row_graphs   H5I_GROUP          <NA>               <NA>

#*************View the 7 parts of loom along with their attributes*************#

# > lookup_loom(loom)
#                       name     link.type    obj_type num_attrs group.nlinks group.mounted dataset.rank  dataset.dims dataset.maxdims dataset.type_class
# 1                    attrs H5L_TYPE_HARD   H5I_GROUP         0            6             0           NA          <NA>            <NA>               <NA>
# 2       attrs/CreationDate H5L_TYPE_HARD H5I_DATASET         0           NA            NA            0             0               0         H5T_STRING
# 3             attrs/Genome H5L_TYPE_HARD H5I_DATASET         0           NA            NA            0             0               0         H5T_STRING
# 4  attrs/LOOM_SPEC_VERSION H5L_TYPE_HARD H5I_DATASET         0           NA            NA            0             0               0         H5T_STRING
# 5           attrs/MetaData H5L_TYPE_HARD H5I_DATASET         0           NA            NA            0             0               0         H5T_STRING
# 6           attrs/RVersion H5L_TYPE_HARD H5I_DATASET         0           NA            NA            0             0               0         H5T_STRING
# 7              attrs/title H5L_TYPE_HARD H5I_DATASET         0           NA            NA            0             0               0         H5T_STRING
# 8                col_attrs H5L_TYPE_HARD   H5I_GROUP         1            4             0           NA          <NA>            <NA>               <NA>
# 9         col_attrs/CellID H5L_TYPE_HARD H5I_DATASET         0           NA            NA            1         32380             Inf         H5T_STRING
# 10   col_attrs/RegulonsAUC H5L_TYPE_HARD H5I_DATASET         1           NA            NA            1         32380           32380       H5T_COMPOUND
# 11         col_attrs/nGene H5L_TYPE_HARD H5I_DATASET         0           NA            NA            1         32380             Inf        H5T_INTEGER
# 12          col_attrs/nUMI H5L_TYPE_HARD H5I_DATASET         0           NA            NA            1         32380             Inf          H5T_FLOAT
# 13              col_graphs H5L_TYPE_HARD   H5I_GROUP         0            0             0           NA          <NA>            <NA>               <NA>
# 14                  layers H5L_TYPE_HARD   H5I_GROUP         0            0             0           NA          <NA>            <NA>               <NA>
# 15                  matrix H5L_TYPE_HARD H5I_DATASET         0           NA            NA            2 32380 x 11364       Inf x Inf          H5T_FLOAT
# 16               row_attrs H5L_TYPE_HARD   H5I_GROUP         1            2             0           NA          <NA>            <NA>               <NA>
# 17          row_attrs/Gene H5L_TYPE_HARD H5I_DATASET         0           NA            NA            1         11364             Inf         H5T_STRING
# 18      row_attrs/Regulons H5L_TYPE_HARD H5I_DATASET         1           NA            NA            1         11364           11364       H5T_COMPOUND
# 19              row_graphs H5L_TYPE_HARD   H5I_GROUP         0            0             0           NA          <NA>            <NA>               <NA>

#************View attributes associated with each part individually************#

# > loom[["col_attrs"]] 
# Class: H5Group
# Filename: C:\Users\KailasammS\Desktop\pySCENIC_AUCELL_filtered.loom
# Group: /col_attrs
# Attributes: last_modified
# Listing:
#        name    obj_type dataset.dims dataset.type_class
#      CellID H5I_DATASET        32380         H5T_STRING
# RegulonsAUC H5I_DATASET        32380       H5T_COMPOUND
#       nGene H5I_DATASET        32380        H5T_INTEGER
#        nUMI H5I_DATASET        32380          H5T_FLOAT

# > loom[["row_attrs"]]
# Class: H5Group
# Filename: C:\Users\KailasammS\Desktop\pySCENIC_AUCELL_filtered.loom
# Group: /row_attrs
# Attributes: last_modified
# Listing:
#     name    obj_type dataset.dims dataset.type_class
#     Gene H5I_DATASET        11364         H5T_STRING
# Regulons H5I_DATASET        11364       H5T_COMPOUND

# Similarly, you can check loom[["attrs"]], loom[["col_graphs"]], 
# loom[["layers"]], loom[["row_graphs"]]

#*******************Retrieve info associated with attributes*******************#

# g <- get_global_attr(loom, key="CreationDate")
# a <- get_col_attr_by_key(loom, "CellID")
# b <- get_col_attr_by_key(loom, "nGene")
# c <- get_col_attr_by_key(loom, "nUMI")
# d <- get_col_attr_by_key(loom, "RegulonsAUC")
# e <- get_row_attr_by_key(loom, "Gene")
# f <- get_row_attr_by_key(loom, "Regulons")

## Attributes ("CellID", "nGene", "nUMI", "Gene") have simple dataset.type_class
## like "H5T_STRING" or "H5T_FLOAT" and are properly retrieved using the above 
## functions. However, attributes ("RegulonsAUC", "Regulons") have complex 
## dataset.type_class like "H5T_COMPOUND" and are not properly retrieved. If you
## check the corresponding R object they lack rownames.

## Use specialized functions to retrieve these attributes like get_regulons() or
## get_regulons_AUC().

## The regulon thresholds are stored as global metadata
# h <- get_global_meta_data(loom)
## You can see the thresholds for each TF inside the above object

## To easily retrieve all thresholds as a named list, use
# i <- get_regulon_thresholds(loom)

#**********************STEP 5: EXTRACT INFO FROM LOOM FILE*********************#

# Refer https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/scenic.html
# Refer https://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/Tutorials_JupyterNotebooks/SCENIC_tutorial_2-ExploringOutput.html

# Finding regulators for known cell types

# NOTE: get_regulons() returns a matrix [TFs as rows, gene as columns]
# Adjust column.attr.name based on loom[["row_attrs"]]
# All genes positively regulated by a TF will have value 1, else 0.

# NOTE: get_regulons_AUC() returns a matrix [Cells as rows, TFs as columns]
# Adjust column.attr.name based on loom[["col_attrs"]]
# This can be accomplished directly with SCopeLoomR itself as follows:
# regulons_AUC <- SCopeLoomR::get_col_attr_by_key(loom = loom, key = "RegulonsAUC")
# regulons_AUC <- t(regulons_AUC)
# colnames(regulons_AUC) <- SCopeLoomR::get_col_attr_by_key(loom = loom,
#                                                           key = "CellID")

# NOTE: get_regulon_thresholds() stores the thresholds as names of the list and 
# TFs/regulons as the list elements. We can reverse this for binarization.

analyze_pySCENIC <- function(celltype){ 
  
  #***************************IMPORT PYSCENIC RESULTS**************************#
  
  # Open loom file and return a .loom file handler
  loom <- SCopeLoomR::open_loom(file.path = paste0(pyscenic_results, "pySCENIC_AUCELL_filtered_viz_", celltype, ".loom"), 
                                mode = "r")
  
  # Read expression data
  exprMat <- SCopeLoomR::get_dgem(loom = loom)
  
  # Normalize expression data (log normalization RECOMMENDED)
  exprMat_log <- base::log2(exprMat+1)
  
  # Get incidence matrix of regulons
  regulons_incidence_mat <- SCopeLoomR::get_regulons(loom = loom,  
                                                     column.attr.name = "Regulons",
                                                     tf.as.name = FALSE,
                                                     tf.sep = "_")
  
  # Convert the incidence matrix of regulons to gene sets
  regulons <- SCENIC::regulonsToGeneLists(incidMat = regulons_incidence_mat)
  
  # Get AUCell scores of each regulon for each cell
  regulons_AUC <- SCopeLoomR::get_regulons_AUC(loom = loom,
                                               column.attr.name = "RegulonsAUC",
                                               rows = "regulons",
                                               columns = "cells")
  
  # Extract underlying matrix from regulons_AUC object
  AUC_mat <- AUCell::getAUC(object = regulons_AUC)
  
  # Extract thresholds for each TF. 
  # NOTE: output of get_regulon_thresholds() is weird. The thresholds are the
  # names while the TFs are values. We reverse this weird format below.
  AUC_thresholds <- SCopeLoomR::get_regulon_thresholds(loom = loom)
  l1 <- names(AUC_thresholds)
  l2 <- as.vector(AUC_thresholds)
  AUC_thresholds <- l1
  names(AUC_thresholds) <- l2
  
  # Get UMAP embeddings
  embeddings <- SCopeLoomR::get_embeddings(loom)
  
  # # Get cluster identities
  # cellClusters <- SCopeLoomR::get_clusterings(loom)
  
  # Close the loom file
  SCopeLoomR::close_loom(loom = loom)
  
  # Load motif enrichment results ??which file??
  #motifEnrichment <- data.table::fread(motifEnrichmentFile, header=T, skip=1)[-3,]
  #colnames(motifEnrichment)[1:2] <- c("TF", "MotifID")
  
  # Identify cells that have AUC score greater than threshold for each regulon.
  # NOTE: The loop takes 10 times longer than SCENIC authors code below
  # AUC_mat_binary <- AUC_mat
  # for (row in 1:nrow(AUC_mat_binary)){
  #   for (col in 1:ncol(AUC_mat_binary)){
  #     if(AUC_mat_binary[row, col] >= as.numeric(AUC_thresholds[row])){
  #       AUC_mat_binary[row, col] <- 1
  #     } else {
  #       AUC_mat_binary[row, col] <- 0
  #     }
  #   }
  # }
  
  #****************************CALCULATE BINARY AUC****************************#
  
  # Calculate binary AUC scores. Output is a named list with barcodes.
  regulonsCells <- base::lapply(X = names(AUC_thresholds), 
                                FUN = function(x) {
                                  trh <- as.numeric(AUC_thresholds[x])
                                  names(which(AUC_mat[x,]>trh))})
  names(regulonsCells) <- names(AUC_thresholds)
  
  # Convert to matrix (regulons with zero assigned cells are lost)
  regulonActivity <- reshape2::melt(regulonsCells)
  AUC_mat_binary <- table(regulonActivity[,2], regulonActivity[,1])
  class(AUC_mat_binary) <- "matrix"
  
  #***********************FIND KEY DIFFERENTIAL REGULONS***********************#
  
  # Read seurat file
  integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn", 
                                      dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  integrated_seurat <- subset(x = integrated_seurat,
                              cell_class %in% c("Mixed", "Unclassified"),
                              invert = TRUE)
  
  # Add pyscenic results as assay to Seurat object
  integrated_seurat[['AUC']] <- SeuratObject::CreateAssayObject(data = AUC_mat)
  
  # Remove unwanted samples
  integrated_seurat <- subset(x = integrated_seurat,
                              Condition == "Tumor")
  
  # Find Differential Regulons using FindMarkers()
  DefaultAssay(integrated_seurat) <- "AUC"
  Idents(integrated_seurat) <- "Sex"
  deg <- FindMarkers(object = integrated_seurat,
                     ident.1 = "Male",
                     ident.2 = "Female",
                     assay = "AUC",
                     slot = "data",
                     logfc.threshold = 0.005,
                     only.pos = FALSE)
  
  deg <- deg %>% 
    dplyr::filter(avg_log2FC != 0 & p_val_adj < 0.05)
  
  DEG_regulons <- rownames(deg) %>% gsub(pattern="-", replacement="_", x=.)
  
  # Find key regulons i.e. regulons with atleast 10 genes
  a <- c()
  for (i in 1:length(regulons)){
    a <- c(a, length(regulons[[i]]))
  }
  key_regulons <- names(regulons)[a>=10]
  
  # Find key DEG regulons
  key_DEG_regulons <- intersect(DEG_regulons, key_regulons)
  AUC_mat <- AUC_mat[key_DEG_regulons,]
  AUC_mat_binary <- AUC_mat_binary[key_DEG_regulons,]
  
  print(length(DEG_regulons))
  print(length(key_regulons))
  print(length(key_DEG_regulons))
  
  #***************************VISUALIZE SCENIC UMAP****************************#
  
  # # Recalcualte UMAP based on pySCENIC values
  # DefaultAssay(object = integrated_seurat) <- "AUC"
  # 
  # # # Scale the data using all features i.e. all regulons
  # # integrated_seurat <- Seurat::ScaleData(object = integrated_seurat,
  # #                                        features = rownames(integrated_seurat@assays$AUC@data),
  # #                                        assay = "AUC",
  # #                                        vars.to.regress = NULL)
  # 
  # # # Run PCA using all features i.e all regulons
  # # integrated_seurat <- Seurat::RunPCA(object = integrated_seurat,
  # #                                     assay = "AUC",
  # #                                     features = rownames(integrated_seurat@assays$AUC@data),
  # #                                     npcs = 50, 
  # #                                     ndims.print = 1,
  # #                                     nfeatures.print = 1,
  # #                                     verbose = TRUE)
  # 
  # # Run UMAP
  # # Since there are only 315 regulons, we can RunUMAP without ScaleData()+RunPCA()
  # integrated_seurat <- Seurat::RunUMAP(object = integrated_seurat,
  #                                      assay = "AUC",
  #                                      #dims = 1:40,
  #                                      #reduction = "pca",
  #                                      features = rownames(integrated_seurat@assays$AUC@data),
  #                                      reduction.name = "AUC_umap",
  #                                      reduction.key = "AUC_UMAP_")
  # 
  # # Also, save UMAP determined by pySCENIC
  # # List of embeddings available:
  # cat(names(embeddings), sep="\n")
  # selectedEmbedding <- embeddings[["SCENIC AUC UMAP"]]
  # colnames(selectedEmbedding) <- c("SCENIC_UMAP_1", "SCENIC_UMAP_2")
  # 
  # # Create DimReducObject and add to object
  # integrated_seurat[['SCENIC_umap']] <- CreateDimReducObject(embeddings = as.matrix(selectedEmbedding), 
  #                                                            key = "SCENIC_UMAP_", 
  #                                                            global = F,
  #                                                            assay = "RNA")
  # 
  # #saveRDS(integrated_seurat, file = paste0(seurat_results, "pySCENIC_seurat_", celltype, ".rds"))
  # 
  # # # Plot UMAP
  # # # NOTE: SCENIC_UMAP (calculated using AUC scores by pyscenic) SHOULD BE SIMILAR
  # # # to AUC_umap (calculated using AUC scores by Seurat).
  # # # umap is default umap calculated using integrated assay in Seurat
  # # split <- "Sex"
  # # Seurat::DimPlot(object = integrated_seurat,
  # #                 reduction = "SCENIC_umap", #c("SCENIC_umap", "umap", "AUC_umap")
  # #                 #group.by = "Treatment",
  # #                 split.by = split,
  # #                 shape.by = NULL,
  # #                 cols = NULL,
  # #                 pt.size = NULL,
  # #                 label = FALSE,
  # #                 label.size = 4,
  # #                 label.color = "black",
  # #                 label.box = FALSE,
  # #                 repel = FALSE,
  # #                 raster = FALSE,
  # #                 ncol = NULL,
  # #                 combine = TRUE)
  # # 
  # # ggplot2::ggsave(filename = paste0("pySCENIC_UMAP", celltype, ".tiff"),
  # #                 plot = last_plot(),
  # #                 device = "jpeg",
  # #                 path = seurat_results,
  # #                 scale = 1,
  # #                 width = dplyr::if_else(is.null(split), 9, 14),
  # #                 height = 11,
  # #                 units = c("in"),
  # #                 dpi = 600,
  # #                 limitsize = TRUE,
  # #                 bg = "white")
  
  #****************************PLOT REGULON HEATMAP****************************#
  
  # NOTE: Since we didnt include clustering info in original loom file use for
  # SCENIC input, we import these from seurat object
  
  # Plot average regulon activity for each cluster/sample
  grouping_var <- "Sample"  #"integrated_snn_res.1.4"
  
  # Calculate average AUC of each regulon in each cluster/sample
  regulons_AUC_cluster <- dplyr::inner_join(x = AUC_mat %>%
                                              t() %>%
                                              as.data.frame() %>%
                                              tibble::rownames_to_column(var = "Cell"),
                                            y = integrated_seurat@meta.data %>%
                                              dplyr::select(grouping_var) %>%
                                              tibble::rownames_to_column(var = "Cell"),
                                            by=c("Cell"="Cell")) %>%
    dplyr::select(everything(), -c("Cell")) %>%
    dplyr::group_by(!!rlang::sym(grouping_var)) %>%
    dplyr::summarise_all(mean, na.rm = TRUE) %>%
    tibble::column_to_rownames(var = grouping_var) %>%
    scale(center = T, scale=T) %>%
    t()
  
  # Calculate average binary AUC of each regulon in each cluster/sample
  # NOTE: If you plot binary AUC without performing scale(center = T, scale=T),
  # you get a heatmap that shows % of cells in each sample which hight regulon 
  # activity.
  regulons_AUC_cluster_binary <- dplyr::inner_join(x = AUC_mat_binary %>%
                                                     t() %>%
                                                     as.data.frame() %>%
                                                     tibble::rownames_to_column(var = "Cell"),
                                                   y = integrated_seurat@meta.data %>%
                                                     dplyr::select(grouping_var) %>%
                                                     tibble::rownames_to_column(var = "Cell"),
                                                   by=c("Cell"="Cell")) %>%
    dplyr::select(everything(), -c("Cell")) %>%
    dplyr::group_by(!!rlang::sym(grouping_var)) %>%
    dplyr::summarise_all(mean, na.rm = TRUE) %>%
    tibble::column_to_rownames(var = grouping_var) %>%
    scale(center = T, scale=T) %>%
    t()
  
  # Save the results
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "AUC")
  openxlsx::writeData(wb, sheet = "AUC", x = regulons_AUC_cluster, rowNames = TRUE)
  openxlsx::addWorksheet(wb, sheetName = "AUC_binary")
  openxlsx::writeData(wb, sheet = "AUC_binary", x = regulons_AUC_cluster_binary, rowNames = TRUE)
  openxlsx::saveWorkbook(wb, 
                         file = paste0(pyscenic_results, "pySCENIC_AUC_matrices.xlsx"), 
                         overwrite = TRUE)
  
  # Plot heatmaps
  mat <- regulons_AUC_cluster
  mat_binary <- regulons_AUC_cluster_binary
  
  # Define column annotation
  col_annotation  <- integrated_seurat@meta.data %>% 
    dplyr::select(grouping_var, Sex) %>% 
    dplyr::distinct_at(grouping_var, .keep_all = TRUE) %>%
    dplyr::arrange(Sample) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = grouping_var)
  
  # Define row annotation
  row_annotation <- data.frame("TFs"= matrix(data="", nrow=nrow(mat), ncol=1))
  rownames(row_annotation) <- rownames(mat)
  
  # Define colors for heatmap
  my_palette <- colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100)
  
  # Define colors for row and column annotation
  groups <- sort(unique(col_annotation[[1]]))
  colors <- c("#BF812D", "#35978F", "#C51B7D", "#7FBC41", "#762A83",
              "#E08214", "#542788", "#D6604D", "#4393C3", "#878787",
              "#1A1A1A", "#FFFFBF", "#9E0142", "#E41A1C", "#377EB8",
              "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628",
              "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#000000")
  colors <- c("#F6D2E0", "#C8E7F5")    # Female: Pink,   Male:Blue (BBN paper)
  colors <- colors[1:length(groups)]
  names(colors) <- groups
  ann_colors = list(colors)
  names(ann_colors) <- colnames(col_annotation)
  
  # # Find TFs that are statistically different between 2 groups
  # # NOTE: We assume normal distribution
  # # Use F.test() to determine if variances are equal or not.
  # sig_genes <- c()
  # sig_pval <- c()
  # for (j in 1:nrow(mat)){
  #   
  #   # DO NOT USE BIND_COLS UNLESS ROW ORDERS IS SAME
  #   data <- dplyr::left_join(x=as.data.frame(mat[j,]) %>% tibble::rownames_to_column(var = "id"),
  #                            y = col_annotation %>% tibble::rownames_to_column(var = "id"),
  #                            by=c("id"="id")) %>%
  #     tibble::column_to_rownames(var = "id")
  #   colnames(data) <- c("values", "Condition")
  #   
  #   if (!is.na(sum(data$values))){
  #     # Check if variances are equal (p <0.05 => unequal variance)
  #     f_test <- var.test(formula = values ~ Condition, 
  #                        data = data,
  #                        alternative = "two.sided")
  #     
  #     if (f_test$p.value < 0.05){
  #       t_test <- t.test(formula = values ~ Condition, 
  #                        data = data,
  #                        alternative = "two.sided",
  #                        var.equal= FALSE)
  #     } else {
  #       t_test <- t.test(formula = values ~ Condition, 
  #                        data = data,
  #                        alternative = "two.sided",
  #                        var.equal= TRUE)
  #     }
  #     if (t_test$p.value < 0.05){
  #       sig_genes <- c(sig_genes, rownames(mat)[j])
  #       sig_pval <- c(sig_pval, t_test$p.value)
  #     }
  #   }
  # }
  # print(data.frame(sig_genes, sig_pval))
  # mat <- mat[sig_genes, sort(rownames(col_annotation))]
  
  # Define how samples will be arranged in the heatmap.
  # Set cluster_cols=FALSE, if you want to arrange samples in specific order
  # Set cluster_cols=TRUE, if you want to arrange samples based on clustering
  mat <- mat[, rownames(col_annotation)]
  mat_binary <- mat_binary[, rownames(col_annotation)]
  
  row <- c("EGR3_(+)", "ELF2_(+)", "ETV4_(+)", "FOS_(+)", "FOSB_(+)", 
           "FOSL2_(+)", "FOXO3_(+)", "GATA3_(+)", "MAFF_(+)", "MAFK_(+)",
           "MEIS1_(+)", "NFATC1_(+)", "NFATC3_(+)", "NFKB1_(+)", 
           "NR2F6_(+)", "PBX1_(+)", "PDLIM5_(+)", "POU2F1_(+)", 
           "SOX15_(+)", "TFAP2A_(+)")
  row <- c("Egr3_(+)", "Elf2_(+)", "Etv4_(+)", "Fos_(+)", "Fosb_(+)",
           "Fosl2_(+)", "Foxo3_(+)", "Gata3_(+)", "Maff_(+)", "Mafk_(+)",
           "Meis1_(+)", "Nfatc1_(+)", "Nfatc3_(+)", "Nfkb1_(+)", 
           "Nr2f6_(+)", "Pbx1_(+)", "Pdlim5_(+)", "Pou2f1_(+)", 
           "Sox15_(+)", "Tfap2a_(+)")
  mat <- mat[row,]
  
  # List genes and samples you want to display in the plot
  display_row <- rownames(mat)
  display_col <- colnames(mat)
  
  # List where you want to have gaps in the heatmap
  gaps_row <- NULL
  gaps_col <- NULL
  # gaps_col <- col_annotation %>% 
  #   dplyr::count(get(colnames(.)[1])) %>% 
  #   dplyr::mutate(n = cumsum(n)) %>%
  #   dplyr::select(n) %>% 
  #   unlist(use.names=FALSE)
  
  # Determine breaks for heatmap color scale.
  # breaks correspond to numerical ranges for the color palette's bins .i.e. 0 to length(my_palette)
  if(max(mat) == 0){
    breaks <- c(seq(from=min(mat), to=0, length.out=ceiling(100/2) + 1), seq(from=1/100, to=1, length.out=floor(100/2)))
  } else if (min(mat) == 0){
    breaks <- c(seq(from=-1, to=0, length.out=ceiling(100/2) + 1), seq(from=max(mat)/100, to=max(mat), length.out=floor(100/2)))
  } else if(min(mat) < -5 & max(mat) > 5){
    breaks <- c(seq(-3, 0, length.out=ceiling(100/2) + 1), seq(max(mat)/100, 3, length.out=floor(100/2)))
  } else{
    breaks <- c(seq(from=min(mat), to=0, length.out=ceiling(100/2) + 1), seq(from=max(mat)/100, to=max(mat), length.out=floor(100/2)))
  }
  
  pheatmap::pheatmap(mat = as.matrix(mat),
                     color = my_palette,
                     kmeans_k = NA,
                     breaks = breaks, 
                     border_color = "white", #"grey90"
                     cellwidth = 5, 
                     cellheight = 5, 
                     scale = "none",        #scaling already done
                     cluster_rows = TRUE,   #cluster the rows
                     cluster_cols = FALSE,  #cluster the columns
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",
                     clustering_method = "complete",
                     legend = TRUE, 
                     legend_breaks = NA,
                     legend_labels = NA, 
                     annotation_row = row_annotation,
                     annotation_col = col_annotation,
                     annotation_colors = ann_colors,
                     annotation_legend = TRUE,
                     annotation_names_row = FALSE,
                     annotation_names_col = FALSE,
                     show_rownames = dplyr::if_else(nrow(mat)<200, TRUE, FALSE, missing = NULL), 
                     show_colnames = dplyr::if_else(ncol(mat)<200, TRUE, FALSE, missing = NULL),
                     fontsize = 8, 
                     fontsize_row = 8, 
                     fontsize_col = 8,
                     gaps_row = NULL,
                     gaps_col = gaps_col,
                     angle_col = "45",
                     fontsize_number = 0.8*fontsize, 
                     labels_row = c(display_row, rep(x="", times=nrow(mat) - length(display_row))),
                     labels_col = c(display_col, rep(x="", times=ncol(mat) - length(display_col))),
                     #width = 15,
                     #height = 10,
                     filename = paste0(pyscenic_results, "pySCENIC_", grouping_var, "_Heatmap1.pdf"))
  
  # cluster and re-order rows
  rowclust = hclust(dist(mat))
  reordered = mat[rowclust$order,]
  
  # # cluster and re-order columns
  # colclust = hclust(dist(t(mat)))
  # reordered = reordered[, colclust$order]
  
  # Save the clustered scores in xlsx
  openxlsx::addWorksheet(wb, sheetName = paste0(celltype, "_AUC"))
  openxlsx::writeData(wb, sheet = paste0(celltype, "_AUC"), x = reordered, rowNames = TRUE)
  
  pheatmap::pheatmap(mat = as.matrix(mat_binary),
                     color = c("white", "black"),
                     kmeans_k = NA,
                     breaks = breaks, 
                     border_color = "grey90", #"white"
                     cellwidth = 5, 
                     cellheight = 5, 
                     scale = "none",        #scaling already done
                     cluster_rows = TRUE,   #cluster the rows
                     cluster_cols = FALSE,  #cluster the columns
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",
                     clustering_method = "complete",
                     legend = TRUE, 
                     legend_breaks = NA,
                     legend_labels = NA, 
                     annotation_row = row_annotation,
                     annotation_col = col_annotation,
                     annotation_colors = ann_colors,
                     annotation_legend = TRUE,
                     annotation_names_row = FALSE,
                     annotation_names_col = FALSE,
                     show_rownames = dplyr::if_else(nrow(mat)<200, TRUE, FALSE, missing = NULL), 
                     show_colnames = dplyr::if_else(ncol(mat)<200, TRUE, FALSE, missing = NULL),
                     fontsize = 8, 
                     fontsize_row = 8, 
                     fontsize_col = 8,
                     gaps_row = NULL,
                     gaps_col = gaps_col,
                     angle_col = "45",
                     fontsize_number = 0.8*fontsize, 
                     labels_row = c(display_row, rep(x="", times=nrow(mat) - length(display_row))),
                     labels_col = c(display_col, rep(x="", times=ncol(mat) - length(display_col))),
                     width = 15,
                     height = 10,
                     filename = paste0(pyscenic_results, "pySCENIC_", grouping_var, "_Heatmap_binary.pdf"))
  
  # cluster and re-order rows
  rowclust = hclust(dist(mat_binary))
  reordered = mat_binary[rowclust$order,]
  
  # # cluster and re-order columns
  # colclust = hclust(dist(t(mat_binary)))
  # reordered = reordered[, colclust$order]
  
  # Save the clustered scores in xlsx
  openxlsx::addWorksheet(wb, sheetName = paste0(celltype, "_AUC_binary"))
  openxlsx::writeData(wb, sheet = paste0(celltype, "_AUC_binary"), x = reordered, rowNames = TRUE)
  
  #****************************************************************************# 
  
  # # Identify regulators for each cluster
  # topRegulators <- reshape2::melt(regulons_AUC_cluster) %>%
  #   dplyr::rename("Regulon" = identity(1), "CellType" = identity(2), "RelativeActivity" = identity(3)) %>%
  #   dplyr::filter(RelativeActivity > 0)
  # 
  # # Save all positive regulators per cluster
  
  #*************************Cell-type specific regulators************************#
  
  # To identify cluster-specific regulons (especially for analyses with many cell
  # types, where some regulons are common to multiple of them) we find specially
  # useful the Regulon Specificity Score (RSS)
  
  # Find RSS across clusters to identify top regulons in each cluster or each sex
  rss <- SCENIC::calcRSS(AUC = AUC_mat,
                         cellAnnotation = integrated_seurat@meta.data$integrated_snn_res.1.4,
                         #cellAnnotation = integrated_seurat@meta.data$Sex,
                         cellTypes = NULL)
  
  # Get list of top 5 regulons for each clsuter
  top_regulons <- rss %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "Regulons") %>%
    tidyr::pivot_longer(!Regulons, names_to = "Cluster", values_to ="AUC_Score") %>%
    dplyr::arrange(Cluster, desc(AUC_Score)) %>%
    dplyr::group_by(Cluster) %>%
    dplyr::slice_head(n=5) %>%
    #dplyr::slice_head(n=30) %>%
    as.data.frame()
  
  # Save the top 5 regulons per cluster
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "Top Regulons")
  openxlsx::writeData(wb, sheet = "Top Regulons", x = top_regulons)
  openxlsx::saveWorkbook(wb, 
                         file = paste0(pyscenic_results, "pySCENIC_Top_Regulons_cluster_", celltype, ".xlsx"), 
                         overwrite = TRUE)
  
  # Plot a dot plot of top regulons of each cluster
  find_top_regulons <- function(col_id){
    
    data <- rss %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "Regulons") %>%
      dplyr::select(Regulons, all_of(col_id)) %>%
      dplyr::rename(AUC_Score = col_id) %>%
      dplyr::arrange(desc(AUC_Score))
    
    ggplot2::ggplot(data = data, 
                    mapping = aes(x=reorder(Regulons, desc(AUC_Score)), y = AUC_Score)) +
      geom_point() +
      geom_point(data = data[1:5,], pch = 21, fill = "red", colour = "red") +
      geom_text_repel(data = data[1:5,], 
                      mapping = aes(label = Regulons), 
                      size = 3.5,                    
                      box.padding = 0.5,
                      max.overlaps = Inf,
                      xlim = c(NA, NA),
                      ylim = c(-Inf,NA),
                      min.segment.length = 1000,
                      position = position_quasirandom()) +
      theme_classic() + 
      labs(x = "Regulons", y = "AUC_Score",	 title = paste0("Cluster", col_id)) +
      my_theme + 
      theme(axis.text.x = element_blank()) +
      scale_color_viridis(option = "D", limits = c(0, 1))
  }
  
  purrr::map(.x = colnames(rss),
             .f = find_top_regulons) %>% 
    cowplot::plot_grid(plotlist = .,
                       align = "hv",
                       axis = "tblr",
                       nrow = 5,  
                       ncol = 5, 
                       rel_widths = 1,
                       rel_heights = 1,
                       labels = NULL,
                       label_size = 14,
                       label_fontfamily = NULL,
                       label_fontface = "bold",
                       label_colour = NULL,
                       label_x = 0,
                       label_y = 1,
                       hjust = -0.5,
                       vjust = 1.5,
                       scale = 1,
                       greedy = TRUE,
                       byrow = TRUE,
                       cols = NULL,
                       rows = NULL)  
  
  ### Save the plot
  ggplot2::ggsave(filename = paste0("pySCENIC_RSS_cluster_", celltype, ".pdf"),
                  plot = last_plot(),
                  device = "pdf",
                  path = pyscenic_results,
                  scale = 1,
                  width = 17,
                  height = 22,
                  units = c("in"),	 
                  dpi = 600,
                  limitsize = TRUE,
                  bg = NULL)
  
  #*****************Finding cell states based on the GRN activity****************#
  
  # To explore whether there are groups of cells that tend to have the same 
  # regulons active, and reveal the network states that are recurrent across 
  # multiple cells, we can cluster the cells based on this regulon activity 
  # (either the continuous or binary AUC matrix). These states would be 
  # equivalent to the attractor states of the network, and they can match the 
  # known cell types, or provide a new/complementary view of the data.
  
  # Remove duplicated regulons
  top_regulons <- unique(top_regulons$Regulons)
  
  # # List of embeddings available:
  # cat(names(embeddings), sep="\n")
  # 
  # # Plot top TFs/Regulons using UMAP
  # selectedEmbedding <- embeddings[["SCENIC AUC UMAP"]]
  # tfsToPlot <- base::gsub("[_()+]", "", top_regulons_list)
  # regulonsToPlot <- top_regulons_list
  # 
  # for (i in 1:ceiling(length(top_regulons)/5)){
  #   
  #   pdf(paste0(seurat_results, "pySCENIC_RegulonActivity_Expression_", i, ".pdf"), width=20, height=15)
  #   par(mfrow=c(2, 5))  # 2 rows, 5 columns
  #   
  #   j <- i*5-4
  #   k <- i*5
  #   # Plot expression:
  #   AUCell::AUCell_plotTSNE(tSNE = selectedEmbedding,
  #                           exprMat = exprMat_log[tfsToPlot[j:k],], 
  #                           plots = c("Expression"))
  #   
  #   # Plot regulon activity:
  #   AUCell::AUCell_plotTSNE(tSNE = selectedEmbedding,
  #                           exprMat = exprMat_log, 
  #                           cellsAUC = regulons_AUC[regulonsToPlot[j:k],], 
  #                           plots  = c("AUC"))
  #   dev.off()
  #   
  # }
  
  # Plot density plot to detect most likely stable states
  library(KernSmooth)
  library(RColorBrewer)
  
  dens2d <- bkde2D(integrated_seurat@reductions$umap@cell.embeddings, 0.5)$fhat
  pdf(paste0(pyscenic_results, "pySCENIC_Density_plot_", celltype, ".pdf"), width=20, height=15)
  image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
  contour(dens2d, add=TRUE, nlevels=6, drawlabels=FALSE)
  
  dev.off()
  
  # In order to find GRN-based cell states, sometimes it is also useful to 
  # binarize the Regulon activity score into "on"/"off". The resulting binarized
  # activity matrix can be used for performing a new 2D projection or clustering,
  # which might reveal more specific GRN-based cell states. We found this 
  # specially helpful for the analysis of cancer datasets, which allowed to 
  # overcome the patient-of origin "batch" effect, and find meaningful cancer 
  # states.
  
  # # Plot all regulons and adjust the threshold if needed
  # # use http://scope.aertslab.org/ for interactive
  # regulonsToPlot <- "ETS1_(+)-motif"
  # options(repr.plot.width=10, repr.plot.height=4) # To set the figure size in Jupyter
  # par(mfrow=c(1,3))
  # AUCell::AUCell_plotTSNE(selectedEmbedding, exprMat_log, 
  #                         regulonAUC[regulonsToPlot,], thresholds = regulonAucThresholds[regulonsToPlot],
  #                         plots=c("AUC", "histogram", "binary"), cex = .5)
  
  
  # Plot binary heatmap
  # # Get clusters enriched for BBN and normal cells
  # c <- base::seq(from = 0, to = 24, by = 1)
  # normal <- c(0, 1, 2, 4, 5, 6, 7, 18, 19, 22, 23)
  # bbn <- c[!(c %in% normal)]
  
  # Define column annotation
  col_annotation <- integrated_seurat@meta.data %>% 
    #dplyr::filter(Treatment == "BBN") %>%
    dplyr::select(integrated_snn_res.1.4) %>%
    dplyr::rename(Cluster = identity(1)) %>%
    dplyr::arrange(Cluster)
  #dplyr::arrange(factor(Cluster, levels=c(normal, bbn)))
  #dplyr::select(Sex) %>%
  #dplyr::arrange(Sex)
  
  # Define row annotation
  row_annotation <- NULL
  
  # Define how samples will be arranged in the heatmap.
  # Set cluster_cols=FALSE, if you want to arrange samples in specific order
  # Set cluster_cols=TRUE, if you want to arrange samples based on clustering
  mat <- AUC_mat_binary
  mat <- mat[top_regulons,]
  mat <- mat[, rownames(col_annotation)]
  
  my_palette <- c("white", "black")
  Cluster = c(`0` = "#BF812D", `1` = "#35978F", `2` = "#C51B7D", `3` = "#7FBC41", `4` = "#762A83", 
              `5` = "#E08214", `6` = "#542788", `7` = "#D6604D", `8` = "#4393C3", `9` = "#878787", 
              `10` = "#1A1A1A", `11` = "#FFFFBF", `12` = "#9E0142", `13` = "#E41A1C", `14` = "#377EB8",
              `15` = "#4DAF4A", `16` = "#984EA3", `17` = "#FF7F00", `18` = "#FFFF33", `19` = "#A65628",
              `20` = "#F781BF", `21` = "#999999", `22` = "#66C2A5", `23` = "#FC8D62", `24` = "#000000")
  ann_colors <- list(Cluster = Cluster[1:length(levels(as.factor(integrated_seurat@meta.data$integrated_snn_res.1.4)))],
                     Sex = c("Female" = "#ee7ae9", "Male" = "#1c86ee"))
  # Plot heatmap
  pheatmap::pheatmap(mat = mat,
                     color = my_palette,
                     border_color = "white", #"grey60",
                     cellwidth = NA, 
                     cellheight = NA, 
                     scale = "none",   
                     cluster_rows = TRUE,   #cluster the rows
                     cluster_cols = TRUE,   #cluster the columns
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",
                     clustering_method = "complete",
                     cutree_rows = NA,
                     cutree_cols = NA,
                     legend = TRUE, 
                     legend_breaks = NA,
                     legend_labels = NA, 
                     annotation_row = row_annotation,
                     annotation_col = col_annotation,
                     annotation_colors = ann_colors,
                     annotation_legend = TRUE,
                     annotation_names_row = TRUE,
                     annotation_names_col = TRUE,
                     show_rownames = dplyr::if_else(nrow(mat)<80, TRUE, FALSE, missing = NULL), 
                     show_colnames = dplyr::if_else(ncol(mat)<30, TRUE, FALSE, missing = NULL),
                     fontsize = 8, 
                     fontsize_row = 8, 
                     fontsize_col = 8,
                     angle_col = c("270", "0", "45", "90", "315"),
                     fontsize_number = 0.8*fontsize, 
                     #labels_row = display_row,
                     #labels_col = display_col, 
                     filename = paste0(pyscenic_results, "pySCENIC_binary_Heatmap_cluster_", celltype, ".pdf"))
}

# Save the clustered scores in xlsx
wb <- openxlsx::createWorkbook()

purrr::map(.x = c("Epithelial", "Fibroblasts", "Myeloid", "Lymphoid"),
           .f = analyze_pySCENIC)

openxlsx::saveWorkbook(wb, file = paste0(pyscenic_results, "pyscenic.xlsx"), 
                       overwrite = TRUE)

# Plot UMAP of regulon expression as well as TF expression identified from 
# binary heatmap
# Exploring the networks in detail (i.e. TFs, target genes and motifs)
# The regulons generated with VSN (and pySCENIC) follow this naming scheme:
#   TF name
# _(+) or _(-): Whether the target genes are positively or negatively correlated with the TF
# -motif or -track sufix: Whether the regulon is based on TF-motifs or ChIP-seq tracks
# Number of regulons:
length(regulons)

# Note than only regulons with 10 genes or more are scored with AUCell
sum(lengths(regulons)>=10)

# Number of target genes in each regulon:
viewTable(cbind(nGenes=lengths(regulons)), options=list(pageLength=10))

# To check whether a specific TF has a regulon
grep("EBF1", names(regulons), value=T) # paste0("^","EBF1","_")

# potential target genes for EBF1
regulons[["EBF1_(+)-motif"]]

# Potential regulators for a given gene
gene <- "CD3E"
names(regulons)[which(sapply(regulons, function(x) "CD3E" %in% x))]

# It might be easier to subset the incidence matrix (TF by target):
genes <- c("CD3E", "MS4A1", "CD14") 
incidMat_subset <- regulons_incidence_mat[,genes]
incidMat_subset <- incidMat_subset[rowSums(incidMat_subset)>0,]
incidMat_subset

# #{
#   # To save the regulons as dataframe, we need to do some data wrangling
#   list_len <- c()
#   for (i in 1:length(regulons)){
#     list_len <- c(list_len, length(regulons[[i]]))
#   }
#   
#   regulons_df <- as.data.frame(matrix(data = NA, nrow = max(list_len), ncol = length(regulons)))
#   for (col in 1:length(regulons)){
#     for(row in 1:length(regulons[[col]])){
#       regulons_df[row, col] <- regulons[[col]][[row]]
#       colnames(regulons_df)[col] <- names(regulons[col])
#     }
#   }
# }


# # Extract AUC scores and AUC binary scores
# AUC_mat <- integrated_seurat@assays$AUC@data
# AUC_mat_binary <- integrated_seurat@assays$AUC_binary@data