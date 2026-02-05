#!/usr/bin/env Rscript

.libPaths("/hpc/home/kailasamms/NGSTools/R_packages")
.libPaths()
chooseCRANmirror(ind=1)

#******************************************************************************#
#                           LOAD NECESSARY PACKAGES                            #
#******************************************************************************#

# Data analysis packages
library("dorothea")
library("viper")
library("Seurat")

# Data wrangling packages     
library("openxlsx")
library("dplyr")
library("tibble")
library("stringr")
library("matrixStats")

# Specialized Graph plotting packages
library("pheatmap")
library("RColorBrewer")
library("ggplot2")

#https://github.com/saezlab/transcriptutorial/blob/master/scripts/04_TranscriptionFactor_activity_with_Dorothea.md
#https://rdrr.io/github/christianholland/dorothea/f/vignettes/single_cell_vignette.Rmd

# Choose the database for Dorothea analysis
#dbs <- "dorothea_hs"
#dbs <- "dorothea_hs_pancancer"
dbs <- "dorothea_mm"
#dbs <- "dorothea_mm_pancancer"

# Read Dorothea Regulons:
#dorothea_regulons <- get(data("dorothea_mm", package = "dorothea"))
dorothea_regulons <- get(dbs)

# Obtain the regulons based on interactions with confidence level A, B and C
regulon <- dorothea_regulons %>%
  dplyr::filter(confidence %in% c("A","B","C"))

###############SINGLE CELL

# Read seurat file
celltype <- NULL
integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn", 
                                    dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))

integrated_seurat <- subset(x = integrated_seurat,
                            subset = (cell_class %in% c("Mixed", "Unclassified")),
                            invert = TRUE)

integrated_seurat <- subset(x = integrated_seurat,
                            subset = (Condition == "Tumor"))

celltypes <- unique(integrated_seurat@meta.data$cell_type)
celltypes <- celltypes[grepl(pattern = "Epithelial|Fibroblasts|Myeloid|Lymphoid", x=celltypes)]

plot_dorothea <- function(i){
  
  obj <- subset(x = integrated_seurat,
                subset = (cell_type == i))
  
  # Compute Viper Scores 
  viper_object <- run_viper(input = obj, 
                            regulons = regulon,
                            tidy = FALSE,
                            assay_key = "RNA",
                            options = list(method = "scale", minsize = 4, 
                                           eset.filter = FALSE, cores = 1, 
                                           verbose = TRUE))
  
  # Find Differential Regulons using FindMarkers()
  DefaultAssay(viper_object) <- "dorothea"
  Idents(viper_object) <- "Sex"
  deg <- FindMarkers(object = viper_object,
                     ident.1 = "Male",
                     ident.2 = "Female",
                     assay = "dorothea",
                     slot = "data",
                     logfc.threshold = 0.005,
                     only.pos = FALSE)
  
  deg <- deg %>% 
    dplyr::filter(avg_log2FC != 0 & p_val_adj < 0.05)
  
  # Get Viper Scores
  viper_scores <- viper_object@assays$dorothea@data
  viper_scores <- viper_scores[rownames(deg),]
  
  #*********************************
  
  # Center and scale the values for each TF across cells and rearrange the data
  #viper_scores_df <- base::scale(t(viper_scores), center = TRUE, scale = TRUE) %>%
  viper_scores_df <- viper_scores %>% t() %>%
    data.frame() %>%
    tibble::rownames_to_column("Cell") %>%
    tidyr::pivot_longer(!Cell, names_to = "TF", values_to = "Activity") %>%
    data.frame()
  
  # Get cluster-cell info
  # Cell column      : FB1_AAACGAAAGATTGATG-1 (for merging using inner_join)
  # sample, sex      : has Sex and sample names (for heatmap)
  cell_cluster_info <- obj@meta.data %>% 
    dplyr::select(Cell, Sample, Sex)
  
  # Merge the data
  viper_scores_df <- viper_scores_df %>% 
    dplyr::inner_join(cell_cluster_info, by=c("Cell"="Cell"))
  
  # Calculate average scores of each TF for each cluster/sample
  summarized_viper_scores <- viper_scores_df %>% 
    #dplyr::group_by(integrated_snn_res.1.4, TF) %>% 
    dplyr::group_by(Sample, TF) %>% 
    dplyr::summarise(avg = mean(Activity))
  
  # # Select the 50 most variable TFs across all clusters.
  # highly_variable_tfs <- summarized_viper_scores %>%
  #   dplyr::group_by(TF) %>%
  #   dplyr::mutate(variance = var(avg))  %>%
  #   dplyr::ungroup() %>% 
  #   dplyr::arrange(desc(variance)) %>%
  #   dplyr::distinct(TF) %>%
  #   dplyr::slice_head(n = 50)
  
  # # Prepare the data for the plot
  # summarized_viper_scores <- summarized_viper_scores %>%
  #   dplyr::semi_join(highly_variable_tfs, by =c("TF"="TF")) %>%
  #   tidyr::pivot_wider(names_from = "integrated_snn_res.1.4" , values_from = "avg") %>%
  #   tibble::column_to_rownames("TF") %>%
  #   as.matrix()
  
  # Plot heatmap
  mat <- tidyr::pivot_wider(data = summarized_viper_scores,
                            names_from = "TF",
                            values_from = "avg") %>%
    tibble::column_to_rownames(var = colnames(.)[1]) %>%
    t()
  
  # Define column annotation
  col_annotation <- obj@meta.data %>% 
    dplyr::select(Sample, Sex) %>% 
    dplyr::distinct_at("Sample", .keep_all=TRUE) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "Sample")
  
  # Define row annotation
  row_annotation <- data.frame("Genes"= matrix(data="", nrow=nrow(mat), ncol=1))
  rownames(row_annotation) <- rownames(mat)
  
  # Define colors for heatmap
  my_palette <- c(colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100)[1:49], 
                  "#FFFFFF",
                  colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100)[50:99])
  
  # Define colors for row and column annotation
  groups <- sort(unique(col_annotation[[1]]))
  colors <- c("#BF812D", "#35978F", "#C51B7D", "#7FBC41", "#762A83",
              "#E08214", "#542788", "#D6604D", "#4393C3", "#878787",
              "#1A1A1A", "#FFFFBF", "#9E0142", "#E41A1C", "#377EB8",
              "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628",
              "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#000000")
  #colors <- c("#74006F", "#FFC606")    # Female:purple, Male:Golds
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
  # 
  # if (length(sig_genes) != 0){
  
  # Define how samples will be arranged in the heatmap.
  # Set cluster_cols=FALSE, if you want to arrange samples in specific order
  # Set cluster_cols=TRUE, if you want to arrange samples based on clustering
  #mat <- mat[sig_genes, sort(rownames(col_annotation))]
  mat <- mat[, sort(rownames(col_annotation))]
  
  # List genes and samples you want to display in the plot
  display_row <- rownames(mat)      
  display_col <- colnames(mat)
  
  # BBN paper specific Fig 4E
  # if (i == "Myeloid - Macrophages, DCs" & species == "Homo sapiens"){
  #   display_row <- c("MYC","THAP1","MAX","E2F4","KLF6","AR","NFYB","TBP",
  #                    "ZBTB33", "SP1", "EGR1","SMAD5","CREB3L1", "CUX1",
  #                    "TFAP2A", "HIF1A", "NR5A1", "TGIF2", "HBP1", "MXI1")
  #   
  # } else if (i == "Myeloid - Macrophages, DCs" & species == "Mus musculus"){
  #   display_row <- c("Myc", "Thap1", "Max", "E2f4", "Klf6", "Ar", "Nfyb", 
  #                    "Tbp", "Zbtb33", "Sp1", "Egr1", "Smad5", "Creb3l1",
  #                    "Cux1", "Tfap2a", "Hif1a", "Nr5a1", "Tgif2", "Hbp1", "Mxi1")
  #   
  # }
  # mat <- mat[display_row, sort(rownames(col_annotation))]
  
  # List where you want to have gaps in the heatmap
  gaps_row <- NULL
  gaps_col <- NULL
  # gaps_col <- col_annotation %>% 
  #   dplyr::count(get(colnames(.)[1])) %>% 
  #   dplyr::mutate(n = cumsum(n)) %>%
  #   dplyr::select(n) %>% 
  #   unlist(use.names=FALSE)
  
  mat <- t(scale(t(mat)))
  
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
  
  # Plot heatmap
  # Error in check.length("fill") :  'gpar' element 'fill' must not be length 0
  # This error means sample_annotation dataframe doesnt match with columns of mat
  # Try using mat = t(mat) to see if it fixes the error.
  # NOTE: If you set scale = none, then you are plotting exact progeny scores
  # NOTE: If you set scale = row, then colors do not represent progeny scores 
  # (i.e. actual activity) but relative activity
  pheatmap::pheatmap(mat = as.matrix(mat),
                     color = my_palette,
                     breaks = breaks, 
                     border_color = "white", #"grey90",
                     cellwidth = 10, 
                     cellheight = 10, 
                     scale = "none",   
                     cluster_rows = TRUE,   #cluster the rows
                     cluster_cols = FALSE,   #cluster the columns
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
                     show_rownames = TRUE, #dplyr::if_else(nrow(mat)<100, TRUE, FALSE, missing = NULL), 
                     show_colnames = TRUE, #dplyr::if_else(ncol(mat)<50, TRUE, FALSE, missing = NULL),
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
                     height = 30,
                     filename = paste0(dorothea_results, "Dorothea_", i, ".pdf"))
  
  # cluster and re-order rows
  rowclust = hclust(dist(mat))
  reordered = mat[rowclust$order,]
  
  # # cluster and re-order columns
  # colclust = hclust(dist(t(mat)))
  # reordered = reordered[, colclust$order]
  
  # Save the clustered scores in xlsx
  openxlsx::addWorksheet(wb, sheetName = i)
  openxlsx::writeData(wb, sheet = i, x = reordered, rowNames = TRUE)
}

# Save the clustered scores in xlsx
wb <- openxlsx::createWorkbook()

purrr::map(.x = celltypes, 
           .f = plot_dorothea)

openxlsx::saveWorkbook(wb, file = paste0(dorothea_results, "dorothea.xlsx"), 
                       overwrite = TRUE)

#######################EXTRA ANALYSIS FOR BBN PAPER

# Find common TF in human and mouse and plot their heatmaps
celltypes <- c("Lymphoid - T", "Myeloid - Macrophages, DCs", "Epithelial", "Fibroblasts", "Myeloid - MDSC","Lymphoid - B","Lymphoid - NK")

wb <- openxlsx::createWorkbook()
for (i in celltypes){
  
  bbn <- openxlsx::read.xlsx("/hpc/home/kailasamms/scratch/scRNASeq_BBN_C57B6/results_dorothea/dorothea.xlsx", sheet = i, rowNames = TRUE)
  bbn_genes <- rownames(bbn)
  chen <- openxlsx::read.xlsx("/hpc/home/kailasamms/scratch/scRNASeq_Chen/results_dorothea/dorothea.xlsx", sheet = i, rowNames = TRUE)
  chen_genes <- rownames(chen)
  
  genes <- intersect(toupper(bbn_genes), toupper(chen_genes))
  index <- which(toupper(bbn_genes) %in% genes)
  bbn_genes <- sort(bbn_genes[index])
  index <- which(toupper(chen_genes) %in% genes)
  chen_genes <- sort(chen_genes[index])
  
  # Plot these common TFs in heatmap
  mat <- bbn[bbn_genes,]
  
  # Define column annotation
  col_annotation <- data.frame("Sample" = c("FB1", "FB2", "FB3", "FB4", "FB5", "MB1", "MB2", "MB3", "MB4", "MB5"),
                               "Sex" = c(rep("Female", 5), c(rep("Male",5)))) %>%
    tibble::column_to_rownames(var = "Sample")
  
  # Define row annotation
  row_annotation <- data.frame("Genes"= matrix(data="", nrow=nrow(mat), ncol=1))
  rownames(row_annotation) <- rownames(mat)
  
  # Define colors for heatmap
  my_palette <- c(colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100)[1:49], 
                  "#FFFFFF",
                  colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100)[50:99])
  
  # Define colors for row and column annotation
  groups <- sort(unique(col_annotation[[1]]))
  colors <- c("#BF812D", "#35978F", "#C51B7D", "#7FBC41", "#762A83",
              "#E08214", "#542788", "#D6604D", "#4393C3", "#878787",
              "#1A1A1A", "#FFFFBF", "#9E0142", "#E41A1C", "#377EB8",
              "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628",
              "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#000000")
  #colors <- c("#74006F", "#FFC606")    # Female:purple, Male:Golds
  colors <- c("#F6D2E0", "#C8E7F5")    # Female: Pink,   Male:Blue (BBN paper)
  colors <- colors[1:length(groups)]
  names(colors) <- groups
  ann_colors = list(colors)
  names(ann_colors) <- colnames(col_annotation)
  
  mat <- mat[, sort(rownames(col_annotation))]
  
  # List genes and samples you want to display in the plot
  display_row <- rownames(mat)      #c("SMAD6","EFNB2","CDX2)
  display_col <- colnames(mat)
  
  # List where you want to have gaps in the heatmap
  gaps_row <- NULL
  gaps_col <- NULL
  
  mat <- t(scale(t(mat)))
  
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
                     breaks = breaks, 
                     border_color = "white", #"grey90",
                     cellwidth = 10, 
                     cellheight = 10, 
                     scale = "none",   
                     cluster_rows = TRUE,   #cluster the rows
                     cluster_cols = FALSE,   #cluster the columns
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
                     show_rownames = TRUE, #dplyr::if_else(nrow(mat)<100, TRUE, FALSE, missing = NULL), 
                     show_colnames = TRUE, #dplyr::if_else(ncol(mat)<50, TRUE, FALSE, missing = NULL),
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
                     height = 30,
                     filename = paste0(dorothea_results, "Dorothea_BBN_common", i, ".pdf"))
  
  # Plot these common TFs in heatmap
  mat <- chen[chen_genes,]
  
  # Define column annotation
  col_annotation <- data.frame("Sample" = c("TF1", "TF2", "TM1", "TM2", "TM3", "TM4", "TM5","TM6"),
                               "Sex" = c(rep("Female", 2), c(rep("Male",6)))) %>%
    tibble::column_to_rownames(var = "Sample")
  
  if (i %in% c("Myeloid - MDSC","Lymphoid - B","Lymphoid - NK")){
    col_annotation <- data.frame("Sample" = c("TF1", "TF2", "TM1", "TM2", "TM3", "TM4", "TM6"),
                                 "Sex" = c(rep("Female", 2), c(rep("Male",5)))) %>%
      tibble::column_to_rownames(var = "Sample")
  }
  
  # Define row annotation
  row_annotation <- data.frame("Genes"= matrix(data="", nrow=nrow(mat), ncol=1))
  rownames(row_annotation) <- rownames(mat)
  
  # Define colors for heatmap
  my_palette <- c(colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100)[1:49], 
                  "#FFFFFF",
                  colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100)[50:99])
  
  # Define colors for row and column annotation
  groups <- sort(unique(col_annotation[[1]]))
  colors <- c("#BF812D", "#35978F", "#C51B7D", "#7FBC41", "#762A83",
              "#E08214", "#542788", "#D6604D", "#4393C3", "#878787",
              "#1A1A1A", "#FFFFBF", "#9E0142", "#E41A1C", "#377EB8",
              "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628",
              "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#000000")
  #colors <- c("#74006F", "#FFC606")    # Female:purple, Male:Golds
  colors <- c("#F6D2E0", "#C8E7F5")    # Female: Pink,   Male:Blue (BBN paper)
  colors <- colors[1:length(groups)]
  names(colors) <- groups
  ann_colors = list(colors)
  names(ann_colors) <- colnames(col_annotation)
  
  mat <- mat[, sort(rownames(col_annotation))]
  
  # List genes and samples you want to display in the plot
  display_row <- rownames(mat)      #c("SMAD6","EFNB2","CDX2)
  display_col <- colnames(mat)
  
  # List where you want to have gaps in the heatmap
  gaps_row <- NULL
  gaps_col <- NULL
  
  mat <- t(scale(t(mat)))
  
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
                     breaks = breaks, 
                     border_color = "white", #"grey90",
                     cellwidth = 10, 
                     cellheight = 10, 
                     scale = "none",   
                     cluster_rows = TRUE,   #cluster the rows
                     cluster_cols = FALSE,   #cluster the columns
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
                     show_rownames = TRUE, #dplyr::if_else(nrow(mat)<100, TRUE, FALSE, missing = NULL), 
                     show_colnames = TRUE, #dplyr::if_else(ncol(mat)<50, TRUE, FALSE, missing = NULL),
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
                     height = 30,
                     filename = paste0(dorothea_results, "Dorothea_Chen_common", i, ".pdf"))
}

################################################################################

#******************************DOROTHEA FOR DESEQ******************************#

# For using dorothea with DESeq2 results, first find top TFs using NES. 
# Then, plot heatmap for these TFs in all samples.

proj <- "Hany_Y"
#proj <- "Hany_YKO" 

# parent directory : directory where input files, results, etc are stored
results_path <- "C:/Users/KailasammS/Box/Saravana@cedars/10. Ongoing Projects/Prince project/"

parent_path <- paste0("C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/", proj, "/")
results_path <- paste0("C:/Users/KailasammS/Desktop/")

# Define if data is log transformed already
already_log <- TRUE

# Define if data is scaled already
already_scaled <- TRUE

# Define if you want to perform unsupervised row and column clustering
# NOTE: If set to FALSE, samples (columns) and genes(rows) will be arranged 
# in alphabetical order (default) in heatmap. If you want to arranged in 
# specific order, define below.
row_clustering <- TRUE    # Usually TRUE
col_clustering <- TRUE   # Usually FALSE

# Define if you want genes or samples to be arranged in alphabetical order
# NOTE: If set to FALSE, write the plot_genes in order you want in heatmap
# NOTE: If row_clustering==TRUE, then row_clustering_alphabetical is irrelevant
row_clustering_alphabetical <- FALSE
col_clustering_alphabetical <- FALSE

# List annotations you want on heatmap
# NOTE: anno_columns MUST match one of the column names in metadata
anno_columns <- c("Condition")

# Define colors for heatmap
my_palette <- colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(200)
#my_palette <- viridis(200)
# my_palette <- c(colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100)[1:49], "#FFFFFF",
#                 colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100)[50:99])

#**********************************INPUT DATA**********************************#

# (I) Read expr data
normalized_counts <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "Results__Normalized_Counts.xlsx"))
normalized_counts <- normalized_counts[, -1]
colnames(normalized_counts)[1] <- "SYMBOL"                                           

# (II) Read metadata
metadata <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "Metadata.xlsx")) %>%
  dplyr::mutate(Sample = make.names(Sample))

# (V) Define any filename you want added to final file
label <- NULL

# Format the input matrices properly
normalized_counts <- normalized_counts %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  tibble::remove_rownames(.) %>%
  tibble::column_to_rownames(var = "SYMBOL") %>%
  dplyr::mutate(across(.cols = everything(), .fns = as.numeric)) %>%
  base::replace(is.na(.), 0)

#**********************************RUN VIPER***********************************#

# Compute Viper Scores 
viper_object <- run_viper(input = normalized_counts, 
                          regulons = regulon,
                          tidy = FALSE,
                          assay_key = "RNA",
                          options = list(method = "scale", minsize = 4, 
                                         eset.filter = FALSE, cores = 1, 
                                         verbose = TRUE))

# Get Viper Scores
viper_scores <- viper_object

# Center and scale the values for each TF across cells and rearrange the data
viper_scores <- base::scale(t(viper_scores), center = TRUE, scale = TRUE) %>%
  t() 

# Plot heatmap
normalized_counts <- viper_scores %>% 
  data.frame() %>% 
  tibble::rownames_to_column("SYMBOL")
colnames(normalized_counts)[1] <- "SYMBOL"
metadata_column <- data.frame("Sample" = colnames(normalized_counts[,-1]), 
                              "Condition" = c(rep(x="scr",times=3), rep(x="yko", times=3)))

metadata_row <- NULL
file_suffix <- ""
#plot_genes <-  normalized_counts$SYMBOL
#disp_genes <- plot_genes
already_log <- FALSE
already_scaled <- TRUE
row_clustering <- TRUE    # Usually TRUE
col_clustering <- TRUE    # Usually FALSE
row_clustering_alphabetical <- FALSE
col_clustering_alphabetical <- FALSE
gaps_in_col <- FALSE
gap_columns <- "Sample"   # Irrelevant if gaps_in_col is FALSE
gaps_in_row <- FALSE
gap_rows <- "Pathway"    # Irrelevant if gaps_in_row is FALSE
columns <- "Sample"
anno_columns <- c("Condition")
anno_rows <- NA
color_by_cols <- TRUE
color_by_rows <- FALSE
my_palette <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100)

# Identify which TFs are significantly different
# Define the target and reference groups
# NOTE: These values MUST be present in "Condition" column of metadata
Target <- "yko" #"YKO_centromere"
Reference <- "scr" #"scr_KO_centromere"
sig_df <- calc_stats(normalized_counts %>% tibble::column_to_rownames("SYMBOL"), metadata_column, file_suffix)

# (III) Define genes to plot
plot_genes <- sig_df %>% 
  #dplyr::filter(padj < 0.05) %>%
  dplyr::filter(abs(log2FC) >= 1 & pval <= 0.05) %>%
  dplyr::select(Gene) %>%
  unlist(use.names = FALSE)
disp_genes <- plot_genes


plot_heatmap(normalized_counts, metadata_column, metadata_row, plot_genes, disp_genes, file_suffix)







# if (length(disp_genes) > 0){
#   
#   # Run function
#   cell_width <- 10
#   cell_height <- 5
#   plot_heatmap(viper_scores, metadata, plot_genes, disp_genes, label)
# }

# # Is data normally distributed? If yes, use F-test.
# # Else use Fligner test
# # Note that, normality test is sensitive to sample size. Small samples most
# # often pass normality tests. Therefore, it’s important to combine visual 
# # inspection and significance test in order to take the right decision.

#   my_data <- as.data.frame(mat[i,])
#   my_data <-dplyr::bind_cols(my_data, col_annotation)
#   colnames(my_data) <- c("weight", "group")

#   # # Shapiro test for normality
#   # shapiro <- shapiro.test(my_data$weight)
#   # if(shapiro$p.value < 0.05){
#   #   print("Data is not normally distributed")
#   #   cat(i, shapiro$p.value)
#   #   #break
#   # } else if(i == nrow(mat)){
#   #   print("Data is normally distributed")
#   # }
#   
#   # Fligner-Killeen test for differences
#   fligner <- fligner.test(weight ~ group, data = my_data)
#   if(fligner$p.value < 0.05){
#      print("Variance is different")
#     cat(i, fligner$p.value)
#   }

###############################################################################

# Interpreting dorothea results

## We read Dorothea Regulons:
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))

## We obtain the regulons based on interactions with confidence level A, B and C
regulon <- dorothea_regulon_human %>%        #dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A","B","C", "D", "E"))  #dplyr::filter(confidence %in% c("A","B","C"))

targets_interest <- regulon %>% 
  dplyr::filter(target == "STAT3") %>% 
  dplyr::select(target)

tfs_interest <- regulon %>% 
  dplyr::filter(target == "CDH12")

volcano_nice(as.data.frame(ttop_KOvsWT[ttop_KOvsWT$ID %in% targets_STAT3,]), 
             FCIndex = 2, pValIndex = 5, IDIndex = 1,nlabels = 20, label = TRUE, 
             straight = FALSE) 

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "TFs")
openxlsx::writeData(wb, sheet = "TFs", x = tfs_interest)
openxlsx::saveWorkbook(wb, file = "C:/Users/KailasammS/Desktop/CDH12.xlsx", overwrite=TRUE)

# Here is how to retrieve all regulons from human:
net <- dorothea::dorothea_hs
head(net)

# We can observe the total number of genes per TF:
n_genes <- net %>%
  group_by(tf) %>%
  summarize(n = n())

ggplot(data=n_genes, aes(x=n)) +
  geom_density() +
  theme(text = element_text(size=12)) +
  xlab('Number of target genes') +
  ylab('densities') +
  theme_bw() +
  theme(legend.position = "none")

# we can visualize how many edges each confidence level adds:
n_edges <- net %>%
  group_by(confidence) %>%
  summarize(n = n())

ggplot(data=n_edges, aes(x=confidence, y=log10(n), color=confidence, fill=confidence)) +
  geom_bar(stat="identity") +
  theme(text = element_text(size=12)) +
  xlab('log10(Number of edges)') +
  ylab('densities') +
  theme_bw() +
  theme(legend.position = "none")

# We can also check how many TFs are repressors, TFs with most of their edges
# with negative mode of regulation (mor), and how many are activators, TFs with
# most of their edges with positive mor:
prop <- net %>%
  group_by(tf, mor) %>%
  summarize(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  filter(mor == 1)
#> `summarise()` has grouped output by 'tf'. You can override using the `.groups`
#> argument.

ggplot(data=prop, aes(x=freq)) +
  geom_density() +
  theme(text = element_text(size=12)) +
  xlab('% of positive edges') +
  ylab('densities') +
  theme_bw() +
  theme(legend.position = "none")


