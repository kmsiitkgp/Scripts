# Refer https://github.com/ventolab/CellphoneDB/blob/master/Docs/RESULTS-DOCUMENTATION.md
# Refer https://github.com/zktuong/ktplots/
# Refer  https://github.com/Teichlab/cellphonedb/issues/15

conda remove -n cpdb --all
conda create -n cpdb python=3.9
conda activate cpdb
conda install -c r r-essentials  # you need R for rpy2 installation, 
pip install -U cellphonedb

#********************PREPARE FILES FOR CELLPHONEDB ANALYSIS********************#

# Your gene ids must be HUMAN. If you are working with another species, you need
# to convert the gene ids to their corresponding orthologs. It is BEST to do 
# this in the seurat object so all subsequent files generated from the seurat 
# object automatically will have HUMAN names

celltype <- NULL
# Load the integrated seurat object
integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn",
                                    dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))

# Export mouse features and rename mouse features
write(x = rownames(integrated_seurat@assays$RNA@data), file = paste0(cellphonedb_results, "features.tsv"))

# Use gprofiler to rename mouse features
# 1. Export gprofiler results. Keep ONLY initial alias and name columns. Remove
# duplicates using ONLY "initial alias". Then, paste mouse genes from features.tsv
# initial alias                     name              Mouse features
# (gprofiler makes them UPPERcasse) (Human ortholog)  (copied from features.tsv)
# XKR4	                            XKR4	            Xkr4
# GM1992	                          None	            Gm1992
# GM19938	                          None	            Gm19938
# GM37381	                          None	            Gm37381
# RP1	                              RP1	              Rp1
# 2. First, do =A2=C2, find FALSE values & correct genes in "initial alias" 
# 3. Then, do =IF(B2="None", TRUE, EXACT(A2,B2)), find FALSE values & correct
# genes in "initial alias"
# 4. Copy the genes from initial alias excluding header & paste into features.tsv

# Import features.tsv after replacing mouose gene names with human gene names
features <- read.table(file = paste0(cellphonedb_results, "features.tsv"), sep = '\t', header = F)
features <- unlist(features, use.names=FALSE)

# Rename genes in Seurat object
DefaultAssay(integrated_seurat) <- "RNA"
integrated_seurat[["SCT"]] <- NULL
integrated_seurat[["integrated"]] <- NULL
rownames(integrated_seurat@assays$RNA@data) <- features
rownames(integrated_seurat@assays$RNA@counts) <- features

# Identify cell types to study. If analyzing all cell-cell interactions, ignore
cell_types <- unique(integrated_seurat@meta.data$cell_type)
cell_types <- cell_types[grepl(pattern = "Epithelial|Fibroblasts|Myeloid|Lymphoid", x = cell_types)]

# Subset the object to contain only cell types of interest
integrated_seurat_m <- subset(x = integrated_seurat,
                              cell_type == "Mixed" | Treatment == "Normal" | Sex == "Female",
                              invert = TRUE)
integrated_seurat_f <- subset(x = integrated_seurat,
                              cell_type  == "Mixed" |Treatment == "Normal" | Sex == "Male",
                              invert = TRUE)
integrated_seurat_m <- subset(x = integrated_seurat_m,
                              subset = cell_type %in% cell_types)
integrated_seurat_f <- subset(integrated_seurat_f,
                              subset = cell_type %in% cell_types)

integrated_seurat_m@meta.data %>% count(Treatment, Sex, cell_type, sub_type)
integrated_seurat_f@meta.data %>% count(Treatment, Sex, cell_type, sub_type)

# Save normalised counts (NOT raw or scaled counts) in mtx format
write(x = rownames(integrated_seurat_m@assays$RNA@data), file = paste0(cellphonedb_results, "Male/features.tsv"))
write(x = colnames(integrated_seurat_m@assays$RNA@data), file = paste0(cellphonedb_results, "Male/barcodes.tsv"))
Matrix::writeMM(obj = integrated_seurat_m@assays$RNA@data, file = paste0(cellphonedb_results, "Male/matrix.mtx"))

write(x = rownames(integrated_seurat_f@assays$RNA@data), file = paste0(cellphonedb_results, "Female/features.tsv"))
write(x = colnames(integrated_seurat_f@assays$RNA@data), file = paste0(cellphonedb_results, "Female/barcodes.tsv"))
Matrix::writeMM(obj = integrated_seurat_f@assays$RNA@data, file = paste0(cellphonedb_results, "Female/matrix.mtx"))

#***************************Generate your meta files***************************#

# IMPORTANT: If you are comparing 2 conditions, YOU MUST specify this in the 
# metafile so that ktplots can make appropriate plots
df <- integrated_seurat_m@meta.data[, c("Cell", "cell_type", "Sample")]
df <- df %>% dplyr::mutate(cell_type = paste0("Male_", cell_type))
write.table(df, file = paste0(cellphonedb_results, "Male_meta.tsv"), sep = '\t', quote = F, row.names = F)

df <- integrated_seurat_f@meta.data[, c("Cell", "cell_type", "Sample")]
df <- df %>% dplyr::mutate(cell_type = paste0("Female_", cell_type))
write.table(df, file = paste0(cellphonedb_results, "Female_meta.tsv"), sep = '\t', quote = F, row.names = F)

#*************************Generate your microenv files*************************#

# ALWAYS generate this file. When analyzing multiple samples, we want to limit
# cell-cell interactions within each sample. So, create the file with cell type
# on column 1 and sample on column 2

df <- integrated_seurat_m@meta.data %>% 
  distinct(cell_type, Sample) %>% 
  tibble::remove_rownames(.) %>%
  dplyr::mutate(cell_type = paste0("Male_", cell_type))
write.table(df, file = paste0(cellphonedb_results, "Male_micro.tsv"), sep = '\t', quote = F, row.names = F)

df <- integrated_seurat_f@meta.data %>% 
  distinct(cell_type, Sample) %>% 
  tibble::remove_rownames(.) %>%
  dplyr::mutate(cell_type = paste0("Female_", cell_type))
write.table(df, file = paste0(cellphonedb_results, "Female_micro.tsv"), sep = '\t', quote = F, row.names = F)

#***************************Generate your DEG files****************************#
# # Compute DEGs if you want to use DEG method
# # If you want to study how epithelial subtypes interact with fibroblast subtypes,
# # then you need to run FindAllMarkers() on (i) epithelial cells and then 
# # (ii) Fibroblasts and then merge the 2 outputs and feed it to cellphonedb.
# # IF you want to see how all cell subtypes interact with each other, then run
# # FindMarkers() on Epi, Fibro, Mye, T cell etc and merge all output.
# 
# # Set "RNA" assay to be the default assay.
# DefaultAssay(integrated_seurat_m) <- "RNA"
# DefaultAssay(integrated_seurat_f) <- "RNA"
# 
# # Change active.ident to major clusters we annotated
# Idents(object = integrated_seurat_m) <- "cell_type"
# Idents(object = integrated_seurat_f) <- "cell_type"
# 
# # Get all markers for any given cluster
# all_markers <- Seurat::FindAllMarkers(object=integrated_seurat_m,
#                                       assay = "RNA",
#                                       features = NULL,
#                                       logfc.threshold = 0.25,
#                                       test.use = "wilcox",
#                                       slot = "data",
#                                       min.pct = 0.1,
#                                       min.diff.pct = 0.1,
#                                       node = NULL,
#                                       verbose = TRUE,
#                                       only.pos = TRUE,
#                                       max.cells.per.ident = Inf,
#                                       random.seed = 1,
#                                       latent.vars = NULL,
#                                       min.cells.feature = 3,
#                                       min.cells.group = 1,
#                                       pseudocount.use = 1,
#                                       mean.fxn = NULL,
#                                       fc.name = NULL,
#                                       base = 2,
#                                       return.thresh = 0.05,
#                                       densify = FALSE)
# 
# # Subset the DEGs and re-arrange columns in specific format below
# all_markers <- all_markers %>% 
#   dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0.1) %>%
#   dplyr::select(cluster, gene, p_val_adj, p_val, avg_log2FC, pct.1, pct.2) 
# 
# # Save the DEGs
# write.table(all_markers, file = paste0(cellphonedb_results, "Male_DEGs.tsv"), sep = '\t', quote = F, row.names = F)
# write.table(all_markers, file = paste0(cellphonedb_results, "Female_DEGs.tsv"), sep = '\t', quote = F, row.names = F)

#***************************RUN CELLPHONEDB ANALYSIS***************************#

# Run 02e_cellphonedb.py in HPC cluster. Then, continue below.

#*****************************PERFORM POST ANALYSIS****************************#

# If you used degs_method, use "relevant_interactions.txt", not "pvalues.txt"
# NOTE: means.txt has mean expression of molecules involved in the interaction
# NOTE: pvalues.txt has pvalues for the interaction
# NOTE: significant_means.txt has ONLY means that are significant (merge of 
# means.txt & pvalues.txt) along with an additional column "rank". Interactions 
# that are rare across all cell-cell pairs have lowest rank value. Interactions
# that have no significant interactions in any cell-cell pair have rank > 1. 
# NOTE: ktplots online tutorials use "means.txt" so that scaling can be done
# properly. All versions of plot_cpdb() plot the non-significant interactions &
# highlight the significant interactions in red cirlce DESPITE setting 
# "keep_significant_only=TRUE" & "p.adjust.method="BH"". While it seems better
# to use significant_means.txt instead of means.txt, doing so will affect 
# scaling and eventually size of dots in the resulting plots. So, use sig_means
# for filtering purpose but plot using means.txt

sig_m1 <- read.table(file = paste0(cellphonedb_results, "statistical_analysis_significant_means_Male.txt"), sep = '\t', header = T)
p1 <- read.table(file = paste0(cellphonedb_results, "statistical_analysis_pvalues_Male.txt"), sep = '\t', header = T)
m1 <- read.table(file = paste0(cellphonedb_results, "statistical_analysis_means_Male.txt"), sep = '\t', header = T)
d1 <- read.table(file = paste0(cellphonedb_results, "statistical_analysis_deconvoluted_Male.txt"), sep = '\t', header = T)

sig_m2 <- read.table(file = paste0(cellphonedb_results, "statistical_analysis_significant_means_Female.txt"), sep = '\t', header = T)
p2 <- read.table(file = paste0(cellphonedb_results, "statistical_analysis_pvalues_Female.txt"), sep = '\t', header = T)
m2 <- read.table(file = paste0(cellphonedb_results, "statistical_analysis_means_Female.txt"), sep = '\t', header = T)
d2 <- read.table(file = paste0(cellphonedb_results, "statistical_analysis_deconvoluted_Female.txt"), sep = '\t', header = T)

sig_means <- combine_cpdb(sig_m1, sig_m2)
means <- combine_cpdb(m1,m2)
pvals <- combine_cpdb(p1,p2)
decon <- combine_cpdb(d1,d2)

# # Remove non-significant interactions and dplyr::arrange() so that order of 
# # interacting pairs is same between p1, sig_m1 and m1
# # NOTE: cellphonedb significance is affected by cell number in each cell type.
# # So, it might be best to use means.txt and focus on interacting pairs with 
# # differences in expression.
# sig_means <- sig_means %>% 
#   dplyr::filter(rank.x < 1 | rank.y < 1) %>% 
#   dplyr::select(everything(), -c(rank.x, rank.y)) %>% 
#   replace(is.na(.), 0.00001) %>%
#   dplyr::arrange(id_cp_interaction)
# 
# sig_m1 <- sig_m1 %>% 
#   dplyr::filter(rank < 1) %>% 
#   dplyr::select(everything(), -c(rank)) %>% 
#   replace(is.na(.), 0.00001) %>%
#   dplyr::arrange(id_cp_interaction)
# 
# sig_m2 <- sig_m2 %>% 
#   dplyr::filter(rank < 1) %>% 
#   dplyr::select(everything(), -c(rank)) %>% 
#   replace(is.na(.), 0.00001) %>%
#   dplyr::arrange(id_cp_interaction)
# 
# for (df in c("pvals", "means", "decon", "p1", "m1", "d1", "p2", "m2", "d2")){
#   test <- get(df)
#   
#   if (df %in% c("pvals", "means", "decon")){
#     test <- test %>% 
#       dplyr::filter(id_cp_interaction %in% sig_means$id_cp_interaction) %>% 
#       dplyr::arrange(id_cp_interaction)
#   } 
#   
#   if(df %in% c("p1", "m1", "d1")){
#     test <- test %>% 
#       dplyr::filter(id_cp_interaction %in% sig_m1$id_cp_interaction) %>% 
#       dplyr::arrange(id_cp_interaction)
#   }
#   
#   if(df %in% c("p2", "m2", "d2")){
#     test <- test %>% 
#       dplyr::filter(id_cp_interaction %in% sig_m2$id_cp_interaction) %>% 
#       dplyr::arrange(id_cp_interaction)
#   }
#   
#   assign(x=df, value=test)
# }
# 
# # Check if all dfs are proper
# if (all(sig_means$interacting_pair == pvals$interacting_pair) & 
#     all(sig_means$interacting_pair == means$interacting_pair) &
#     all(colnames(sig_means) == colnames(pvals)) &
#     all(colnames(sig_means) == colnames(means)) &
#     all(sig_m1$interacting_pair == p1$interacting_pair) & 
#     all(sig_m1$interacting_pair == m1$interacting_pair) & 
#     all(colnames(sig_m1) == colnames(p1)) &
#     all(colnames(sig_m1) == colnames(m1)) & 
#     all(sig_m2$interacting_pair == p2$interacting_pair) &
#     all(sig_m2$interacting_pair == m2$interacting_pair) &
#     all(colnames(sig_m2) == colnames(p2)) & 
#     all(colnames(sig_m2) == colnames(m2))) {
#   print("All OK")
# }

# " " and "|" in column names get converted to "." when importing files into R.
# "Myeloid Cells" becomes "Myeloid.Cells" while "Fibroblasts|Epithelial Cells" 
# becomes "Fibroblasts.Epithelial.Cells". So, first, we correct these errors.
# This apparently happens only with pvalues.txt & significant_means.txt but not 
# with means.txt
for (df in c("sig_means", "pvals", "means", "decon", "sig_m1", "p1", "m1", "d1", "sig_m2", "p2", "m2", "d2")){
  test <- get(df)
  
  colnames(test) <- stringi::stri_replace_all_regex(colnames(test),
                                                    pattern=c("\\.", " Male", " Female", "   ", "  DCs"),
                                                    replacement=c(" ", "|Male", "|Female", " - ", ", DCs"),
                                                    vectorize_all = FALSE)
  assign(x=df, value=test)
}


# Identify all "Differential Interactions" (DI) between 2 conditions
# Mean1/Mean2 > 1.5 or Mean1/Mean2 < 0.66 are considered biologically significant
#celltypes <- setdiff(unique(integrated_seurat@meta.data$cell_type), c("Batch specific", "Unknown"))
celltypes <- c("Epithelial Cells", "Fibroblasts", "Lymphoid - B Cells", 
               "Lymphoid - NK Cells", "Lymphoid - T Cells", "Myeloid - Mast Cells", 
               "Myeloid - Macrophages, DCs")
count <- 0
for (celltype1 in c("Epithelial Cells", "Fibroblasts", "Lymphoid - B Cells", 
                    "Lymphoid - NK Cells", "Lymphoid - T Cells", "Myeloid - Mast Cells", 
                    "Myeloid - Macrophages, DCs")){
  
  celltypes <- setdiff(celltypes, celltype1)
  for (celltype2 in celltypes){ 
    
    count <- count+1
    cat("\n", count, "\t", celltype1, "\t", celltype2)
    
    pattern1 <- paste0("Male_", celltype1, "|", "Male_", celltype2)
    pattern2 <- paste0("Female_", celltype1, "|", "Female_", celltype2)
    cat("\n", count, "\t", pattern1, "\t", pattern2)
    
    pattern3 <- paste0("Male_", celltype2, "|", "Male_", celltype1)
    pattern4 <- paste0("Female_", celltype2, "|", "Female_", celltype1)
    cat("\n", count, "\t", pattern3, "\t", pattern4)
    
    test1 <- means %>%   #sig_means %>%
      dplyr::select("interacting_pair", all_of(pattern1), all_of(pattern2)) %>%
      dplyr::distinct_at("interacting_pair", .keep_all = TRUE) %>%
      tibble::column_to_rownames(colnames(.)[1]) %>%
      replace(.==0, 0.001) %>%
      dplyr::mutate(fc = get(colnames(.)[1])/get(colnames(.)[2]), diff = abs(get(colnames(.)[1]) - get(colnames(.)[2]))) %>% 
      dplyr::filter((fc >= 1.5 | fc <= 0.66) & (diff >= 0.25))
    
    test2 <- means %>%  #sig_means %>%
      dplyr::select("interacting_pair", all_of(pattern3), all_of(pattern4)) %>%
      dplyr::distinct_at("interacting_pair", .keep_all = TRUE) %>%
      tibble::column_to_rownames(colnames(.)[1]) %>%
      replace(.==0, 0.001) %>%
      dplyr::mutate(fc = get(colnames(.)[1])/get(colnames(.)[2]), diff = abs(get(colnames(.)[1]) - get(colnames(.)[2]))) %>% 
      dplyr::filter((fc >= 1.5 | fc <= 0.66) & (diff >= 0.25))
    
    # Initialize DI for each cell-cell pair
    DI <- c(rownames(test1), rownames(test2))
    DI <- unique(DI)
    
    if (length(DI) > 0){
      means_subset <- means %>% dplyr::filter(interacting_pair %in% DI)
      pvals_subset <- pvals %>% dplyr::filter(interacting_pair %in% DI)
      
      # To see interactions between 2 specific cell types
      # NOTE: scale = TRUE is  (xi - mean)/sd i.e Mean normalization
      # standard_scale = TRUE is (xi - min(x))/(max(x)-min(x)) i.e. Min-Max normalization
      plot_cpdb(cell_type1 = celltype1, 
                cell_type2 = celltype2,
                scdata = integrated_seurat,
                idents = "cell_type",     # column in metadata defining celltypes
                split.by = "Sex",         # column in metadata defining groups 
                means = means_subset, 
                pvals = pvals_subset,
                degs_analysis = FALSE,    # adjust this based on analysis you did
                #p.adjust.method = "BH",
                #keep_significant_only = TRUE,
                scale = TRUE,
                #standard_scale = TRUE,
                #gene.family = "chemokines",  # c("chemokines", "Th1", "Th2", "Th17", "Treg", "costimulatory", "coinhibitory", "niche")
                genes = NULL,
                cluster_rows = TRUE,
                max_size = 5,
                col_option = viridis::viridis(50),
                highlight = "red",      # color for highlighting p <0.05
                highlight_size = 1)     # stroke size for highlight if p < 0.05
      
      ggplot2::ggsave(filename = paste0(count, "_", celltype1, "_", celltype2, ".jpeg"),
                      plot = last_plot(),
                      device = "jpeg",
                      path = cellphonedb_results,
                      scale = 1,
                      width = 8.5,
                      #height = dplyr::if_else(length(DI) < 15, 7, length(DI)*0.2+3),
                      height = length(DI)*0.2+3,
                      units = c("in"),
                      dpi = 300,
                      limitsize = TRUE,
                      bg = "white")
    }
  }
}

# To see which cell types interact the most
plot_cpdb_heatmap(scdata = integrated_seurat,
                  idents = "cell_type",
                  pvals = pvals,
                  degs_analysis = FALSE,
                  log1p_transform = FALSE,
                  show_rownames = TRUE,
                  show_colnames = TRUE,
                  scale = "none",
                  cluster_cols = TRUE,
                  cluster_rows = TRUE,
                  border_color = "white",
                  fontsize_row = 11,
                  fontsize_col = 11,
                  family = "Arial",
                  low_col = "dodgerblue4",
                  mid_col = "peachpuff",
                  high_col = "deeppink4",
                  alpha = 0.05,
                  return_tables = FALSE,
                  verbose = FALSE,
                  file = "2.pdf")

# celltypes <- setdiff(unique(integrated_seurat@meta.data$cell_type), c("Batch specific", "Unknown"))
# count <- 0
# for (celltype1 in setdiff(unique(integrated_seurat@meta.data$cell_type), c("Batch specific", "Unknown"))){
#   celltypes <- setdiff(celltypes, celltype1)
#   
#   for (celltype2 in celltypes){ 
#     count <- count+1
#     cat("\n", count, "\t", celltype1, "\t", celltype2)
#     
#     # To see interactions between 2 specific cell types
#     # NOTE: scale = TRUE is  (xi - mean)/sd i.e Mean normalization
#     # standard_scale = TRUE is (xi - min(x))/(max(x)-min(x)) i.e. Min-Max normalization
#     plot_cpdb(cell_type1 = celltype1, 
#               cell_type2 = celltype2,
#               scdata = integrated_seurat,
#               idents = "cell_type",     # column in metadata defining celltypes
#               split.by = "Sex",         # column in metadata defining groups 
#               means = means, 
#               pvals = pvals,
#               degs_analysis = FALSE,    # adjust this based on analysis you did
#               p.adjust.method = "BH",
#               keep_significant_only = TRUE,
#               scale = TRUE,
#               #standard_scale = TRUE,
#               #gene.family = "chemokines",       # c("chemokines", "Th1", "Th2", "Th17", "Treg", "costimulatory", "coinhibitory", "niche")
#               genes = NULL,
#               cluster_rows = TRUE,
#               max_size = 5,
#               col_option = viridis::viridis(50),
#               highlight = "red",      # color for highlighting p <0.05
#               highlight_size = 1)     # stroke size for highlight if p < 0.05
#     #+ small_guide() + small_axis()
#     
#     ggplot2::ggsave(filename = paste0(count, "_", celltype1, "_", celltype2, ".jpeg"),
#                     plot = last_plot(),
#                     device = "jpeg",
#                     #path = diagnostics_path,
#                     scale = 1,
#                     width = 8.5,
#                     height = 11,
#                     units = c("in"),
#                     dpi = 300,
#                     limitsize = TRUE,
#                     bg = "white")
#   }
# }
# 
# # To plot in circos format
# plot_cpdb2(cell_type1 = "CD4_Tem|CD4_Tcm|CD4_Treg", # same usage style as plot_cpdb
#            cell_type2 = "cDC",
#            idents = 'fine_clustering',
#            split.by = 'treatment_group_1',
#            scdata = sce,
#            means = m1,
#            pvals = p1,
#            deconvoluted = deconvoluted, # new options from here on specific to plot_cpdb2
#            gene_symbol_mapping = 'index', # column name in rowData holding the actual gene symbols if the row names is ENSG Ids. Might be a bit buggy
#            desiredInteractions = list(c('CD4_Tcm', 'cDC1'), c('CD4_Tcm', 'cDC2'), c('CD4_Tem', 'cDC1'), c('CD4_Tem', 'cDC2	'), c('CD4_Treg', 'cDC1'), c('CD4_Treg', 'cDC2')),
#            interaction_grouping = interaction_grouping,
#            edge_group_colors = c("Activating" = "#e15759", "Chemotaxis" = "#59a14f", "Inhibitory" = "#4e79a7", "   Intracellular trafficking" = "#9c755f", "DC_development" = "#B07aa1"),
#            node_group_colors = c("CD4_Tcm" = "#86bc86", "CD4_Tem" = "#79706e", "CD4_Treg" = "#ff7f0e", "cDC1" = "#bcbd22"  ,"cDC2" = "#17becf"),
#            keep_significant_only = TRUE,
#            standard_scale = TRUE,
#            remove_self = TRUE)
# 
# 
