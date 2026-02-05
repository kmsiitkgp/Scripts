#!/usr/bin/env Rscript

.libPaths("/hpc/home/kailasamms/NGSTools/R_packages")
.libPaths()
chooseCRANmirror(ind=1)

#******************************************************************************#
#                           LOAD NECESSARY PACKAGES                            #
#******************************************************************************#

# Data analysis packages
library("progeny")

# Data wrangling packages     
library("openxlsx")
library("dplyr")
library("tibble")
library("stringr")

# Graph plotting packages
library("RColorBrewer")
library("cowplot")
library("viridis")
library("ggplot2")

# Specialized Graph plotting packages
library("pheatmap")

# Single cell analysis packages
library("Seurat")

# 

#******************************************************************************#
#                           DECLARE GLOBAL VARIABLES                           #
#******************************************************************************#

# parent directory : directory where input files, results, etc are stored
results_path <- "C:/Users/KailasammS/Box/Saravana@cedars/10. Ongoing Projects/Prince project/"

# Choose if data is from human or mice. We will adjust gene names accordingly.
species <- "Homo sapiens"
species <- "Mus musculus"

#******************************************************************************#

# Calculate Progeny pathway scores for each cell
# NOTE: https://saezlab.github.io/progeny/articles/ProgenySingleCell.html
# uses seurat to scale. However, I observed that mean and std dev are not 0 and 
# 1 after scaling using Seurat. So, stick to progeny() for scaling by setting
# scale = TRUE or manually do scaling by setting scale=FALSE.
# Bulk RNA Seq has higher coverage than single cell RNA seq. So, you can model
# using just 100 top genes using top = 100. For single cell, use top=200 to 500.

analyze_progeny <- function(celltype){
  
  # Read seurat file
  integrated_seurat <- readRDS(paste0(results_path, "integrated_snn_", celltype, ".rds"))
  
  progeny_object <- progeny::progeny(expr = integrated_seurat,
                                     scale = FALSE, 
                                     organism = dplyr::if_else(species == "Homo sapiens", "Human", "Mouse"),
                                     top = 500, # if you use more genes to model, activity may be inaccurate
                                     perm = 1,
                                     verbose = FALSE,
                                     z_scores = FALSE,
                                     get_nulldist = FALSE,
                                     assay_name = "RNA",
                                     return_assay = TRUE)
  # Get Progeny Scores
  progeny_scores <- progeny_object@assays$progeny@data
  
  # Center and scale the values for each pathway across cells and rearrange the data
  progeny_scores_df <- scale(t(progeny_scores), center = TRUE, scale = TRUE) %>%
    data.frame() %>%
    tibble::rownames_to_column("Cell") %>%
    tidyr::pivot_longer(!Cell, names_to = "Pathway", values_to = "Activity")
  
  # Get cell info
  # Cell column has FB1_AAACGAAAGATTGATG-1 for inner_join
  # cell_type has BBN, Normal, Both etc for annotation of clusters in heatmap
  # integrated_snn_res.0.8 has cluster numbers
  cell_info <- integrated_seurat@meta.data %>% 
    dplyr::select(Cell, cell_type, Sex, integrated_snn_res.0.8)
  
  # Merge the data
  progeny_scores_df <- dplyr::inner_join(progeny_scores_df, cell_info, by=c("Cell"="Cell"))
  
  # Calculate average scores of each pathway for each cluster (sex)
  summarized_progeny_scores <- progeny_scores_df %>%
    dplyr::group_by(Sex, Pathway) %>%
    dplyr::summarise(avg = mean(Activity)) %>%
    tidyr::pivot_wider(names_from = Sex, values_from = avg) %>%
    tibble::column_to_rownames("Pathway")
  
  # # Plot all cells
  # summarized_progeny_scores <- progeny_scores_df %>%
  #   tidyr::pivot_wider(id_cols = Pathway, names_from = Cell, values_from = Activity) %>%
  #   tibble::column_to_rownames("Pathway")
  
  # Define column annotation
  col_annotation <- integrated_seurat@meta.data %>% 
    dplyr::select(Sex) %>%
    dplyr::arrange(Sex)
  
  # sample_annotation <- cell_info %>%
  #   dplyr::count(integrated_snn_res.0.8, cell_type) %>%
  #   tibble::column_to_rownames("integrated_snn_res.0.8") %>%
  #   dplyr::select(cell_type)
  # #rownames(sample_annotation) <- paste0("c",rownames(sample_annotation))
  
  # Define row annotation
  row_annotation <- NULL
  
  # Define how samples will be arranged in the heatmap.
  # Set cluster_cols=FALSE, if you want to arrange samples in specific order
  # Set cluster_cols=TRUE, if you want to arrange samples based on clustering
  mat <- summarized_progeny_scores
  #mat <- mat[, rownames(col_annotation)]
  
  # Define colors for heatmap
  my_palette <- colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100)
  
  # Define colors for row and column annotation
  Cluster = c(`0` = "#BF812D", `1` = "#35978F", `2` = "#C51B7D", `3` = "#7FBC41", `4` = "#762A83", 
              `5` = "#E08214", `6` = "#542788", `7` = "#D6604D", `8` = "#4393C3", `9` = "#878787", 
              `10` = "#1A1A1A", `11` = "#FFFFBF", `12` = "#9E0142", `13` = "#E41A1C", `14` = "#377EB8",
              `15` = "#4DAF4A", `16` = "#984EA3", `17` = "#FF7F00", `18` = "#FFFF33", `19` = "#A65628",
              `20` = "#F781BF", `21` = "#999999", `22` = "#66C2A5", `23` = "#FC8D62", `24` = "#000000")
  my_ann_colors <- list(Cluster = Cluster[1:length(levels(as.factor(integrated_seurat@meta.data$integrated_snn_res.0.8)))],
                        Sex = c("Female" = "#ee7ae9", "Male" = "#1c86ee"))
  
  # Plot heatmap
  # NOTE: Error in check.length("fill"):'gpar' element 'fill' must not be length 0
  # This error means sample_annotation dataframe doesnt match with columns of mat
  # Try using mat = t(mat) to see if it fixes the error.
  # NOTE: If you set scale = none, then you are plotting exact progeny scores
  # NOTE: If you set scale = row, then colors do not represent progeny scores 
  # (i.e. actual activity) but relative activity
  pheatmap::pheatmap(mat = as.matrix(mat),
                     color = my_palette,
                     # breaks correspond to numerical ranges for the color palette's bins .i.e. 0 to length(my_palette)
                     breaks = c(seq(from=min(mat), to=0, length.out=ceiling(100/2) + 1), seq(from=max(mat)/100, to=max(mat), length.out=floor(100/2))),
                     # breaks = c(seq(-3, 0, length.out=ceiling(100/2) + 1), seq(max(mat)/100, 3, length.out=floor(100/2))), 
                     border_color = "white", #"grey60",
                     cellwidth = NA, 
                     cellheight = NA, 
                     scale = "none",   
                     cluster_rows = TRUE,   #cluster the rows
                     cluster_cols = TRUE,   #cluster the columns
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",
                     clustering_method = "complete",
                     legend = TRUE, 
                     legend_breaks = NA,
                     legend_labels = NA, 
                     annotation_row = row_annotation,
                     annotation_col = col_annotation,
                     annotation_colors = dplyr::if_else(nrow(col_annotation)+nrow(row_annotation) > 0, ann_colors, c()),
                     annotation_legend = TRUE,
                     annotation_names_row = FALSE,
                     annotation_names_col = FALSE,
                     show_rownames = dplyr::if_else(nrow(mat)<80, TRUE, FALSE, missing = NULL), 
                     show_colnames = dplyr::if_else(ncol(mat)<50, TRUE, FALSE, missing = NULL),
                     fontsize = 8, 
                     fontsize_row = 8, 
                     fontsize_col = 8,
                     angle_col = c("270", "0", "45", "90", "315"),
                     fontsize_number = 0.8*fontsize, 
                     labels_row = display_row,
                     labels_col = display_col,
                     width = ncol(mat)*0.075+5,
                     filename = paste0(results_path, "progeny_", celltype, ".pdf"))
  
  # Save the scores in xlsx
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "progeny-scores")
  openxlsx::writeData(wb, sheet = "progeny-scores", x = summarized_progeny_scores, rowNames = TRUE)
  
  openxlsx::saveWorkbook(wb, file = paste0(results_path, "progeny_", celltype, ".xlsx"), 
                         overwrite = TRUE)
}

purrr::map(.x = c("Epithelial Cells", "Fibroblasts", "Myeloid Cells", "T Cells"),
           .f = analyze_progeny)

################################################################################

# Performing progeny of DESEq2 results are a bit different
# normalized_counts should have genes as rownames; all samples as column
# DEGs should have genes as rownames; logFC (or t in case of limma) as the only column
# We keep all DEGs irrespective of padj value

files1 <- c("C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/Hany_LOY/Results__id_Kdm5d_OE_vs_scr_OE_DEGs.xlsx",
            "C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/Hany_LOY/Results__id_Kdm5d_KO_vs_scr_KO_DEGs.xlsx",
            "C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/Hany_LOY/Results__id_Uty_OE_vs_scr_OE_DEGs.xlsx",
            "C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/Hany_LOY/Results__id_Uty_KO_vs_scr_KO_DEGs.xlsx",
            "C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/Hany_LOY/Results__id_Kdm5d_OE_Tumor_vs_scr_OE_Tumor_DEGs.xlsx",
            "C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/Hany_LOY/Results__id_Uty_OE_Tumor_vs_scr_OE_Tumor_DEGs.xlsx",
            "C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/Hany_LOY/Results__id_Y_pos_vs_Y_neg_DEGs.xlsx",
            "C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/Hany_LOY/Results__id_YKO_centromere_vs_scr_KO_centromere_DEGs.xlsx")

files2 <- c("C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/Hany_LOY/Results__id_Kdm5d_OE_vs_scr_OE_Normalized_Counts.xlsx",
            "C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/Hany_LOY/Results__id_Kdm5d_KO_vs_scr_KO_Normalized_Counts.xlsx",
            "C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/Hany_LOY/Results__id_Uty_OE_vs_scr_OE_Normalized_Counts.xlsx",
            "C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/Hany_LOY/Results__id_Uty_KO_vs_scr_KO_Normalized_Counts.xlsx",
            "C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/Hany_LOY/Results__id_Kdm5d_OE_Tumor_vs_scr_OE_Tumor_Normalized_Counts.xlsx",
            "C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/Hany_LOY/Results__id_Uty_OE_Tumor_vs_scr_OE_Tumor_Normalized_Counts.xlsx",
            "C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/Hany_LOY/Results__id_Y_pos_vs_Y_neg_Normalized_Counts.xlsx",
            "C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/Hany_LOY/Results__id_YKO_centromere_vs_scr_KO_centromere_Normalized_Counts.xlsx")

names <- c("Kdm5d_OE",
           "Kdm5d_KO",
           "Uty_OE",
           "Uty_KO",
           "Kdm5d_OE_Tumor",
           "Uty_OE_Tumor",
           "Y_pos",
           "YKO")

meta_data <- openxlsx::read.xlsx(xlsxFile = "C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/Hany_LOY/Metadata.xlsx",
                                 rowNames=TRUE)

#*******************MERGE ALL DATA TOGETHER (NOT RECOMMENDED)******************#
# # Find total number of unique genes
# normalized_genes <-c()
# deg_genes <- c()
# for (i in 1:length(files1)){
#   
#   n <- openxlsx::read.xlsx(xlsxFile = files2[i], rowNames=TRUE)
#   d <- openxlsx::read.xlsx(xlsxFile = files1[i], rowNames=TRUE)
#   normalized_genes <- c(normalized_genes, n$SYMBOL)
#   deg_genes <- c(deg_genes, d$SYMBOL)
# }
# normalized_genes <-unique(normalized_genes)
# deg_genes <- unique(deg_genes)
# 
# # Create empty dataframes to store the data
# normalized_counts <- data.frame("SYMBOL" = normalized_genes)
# DEGs <- data.frame("SYMBOL" = normalized_genes)

# 
# # Merge normalized counts of all samples and log2FoldChange of all comparisons 
# for (i in 1:length(files1)){
#   
#   n <- openxlsx::read.xlsx(xlsxFile = files2[i], rowNames=TRUE) %>%
#     dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
#     tibble::remove_rownames(.)
#   
#   d <- openxlsx::read.xlsx(xlsxFile = files1[i], rowNames=TRUE) %>% 
#     dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
#     tibble::remove_rownames(.) %>%
#     dplyr::select(SYMBOL, log2FoldChange) %>%
#     dplyr::rename(!!names[i] := log2FoldChange)
#   
#   normalized_counts <- dplyr::left_join(normalized_counts, n, by=c("SYMBOL"="SYMBOL"))
#   DEGs <- dplyr::left_join(DEGs, d, by=c("SYMBOL"="SYMBOL"))
# }
# 
# # Remove duplicated samples
# normalized_counts <- normalized_counts[,!grepl(pattern = ".y$", x = colnames(normalized_counts))]
# colnames(normalized_counts) <- gsub(pattern = ".x$", replacement = "", x = colnames(normalized_counts))

normalized_counts <- read.xlsx(paste0(results_path, "pg_matrix.xlsx")) %>%
  dplyr::mutate(SYMBOL = make.names(SYMBOL, unique=TRUE)) %>%
  tibble::column_to_rownames("SYMBOL") %>%
  dplyr::mutate(across(.cols = everything(), .fns = as.numeric)) %>%
  base::replace(is.na(.), 0)

normalized_counts <- normalized_counts[,1:6]
normalized_counts <- normalized_counts[,7:12]
  
#*********************IMPORT INDIVIDUAL FILES (RECOMMENDED)********************#

for (i in 1:length(files2)){
  
  normalized_counts <- openxlsx::read.xlsx(xlsxFile = files2[i], rowNames=TRUE)
  
  # Format the input matrices properly
  normalized_counts <- normalized_counts %>%
    # dplyr::mutate(SYMBOL = base::gsub(pattern = "[- ]", replacement = ".", x = SYMBOL)) %>%
    # dplyr::rename_with(.fn = ~base::gsub(pattern = "-| ", replacement = ".", x = .x), .cols = everything()) %>%
    dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
    tibble::remove_rownames(.) %>%
    tibble::column_to_rownames(var = "SYMBOL") %>%
    dplyr::mutate(across(.cols = everything(), .fns = as.numeric)) %>%
    base::replace(is.na(.), 0)
  
  # Calculate progeny scores using normalized counts and scale the values
  progeny_scores <- progeny::progeny(expr = as.matrix(normalized_counts), 
                                     scale = TRUE, 
                                     organism = dplyr::if_else(species=="Homo sapiens", "Human", "Mouse"), 
                                     top = 100, 
                                     perm = 1,
                                     verbose = FALSE,
                                     z_scores = FALSE,
                                     get_nulldist = FALSE,
                                     assay_name = "RNA",
                                     return_assay = FALSE) %>%
    t()
  
  # Plot heatmap
  normalized_counts <- progeny_scores %>% 
    data.frame() %>% 
    tibble::rownames_to_column("SYMBOL")
  colnames(normalized_counts)[1] <- "SYMBOL"
  metadata_column <- data.frame("Sample" = colnames(normalized_counts[,-1]), 
                                "Condition" = c(rep(x="scr",times=3), rep(x="yko", times=3)))
  metadata_row <- NULL
  file_suffix <- ""
  plot_genes <-  normalized_counts$SYMBOL
  disp_genes <- plot_genes
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
  plot_heatmap(normalized_counts, metadata_column, metadata_row, plot_genes, disp_genes, file_suffix)

  # Define column annotation
  col_annotation <- meta_data[colnames(mat),] %>% dplyr::select(Condition)
  rownames(col_annotation) <- colnames(mat)
  
  # Define row annotation
  row_annotation <- data.frame("Genes"= matrix(data="", nrow=nrow(mat), ncol=1))
  rownames(row_annotation) <- rownames(mat)
  
  # Define colors for heatmap
  my_palette <- colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100)
  
  # Define colors for row and column annotation
  #ann_colors = list(id = c(Female = "#74006F", Male = "#FFC606")) # purple/gold
  #ann_colors = list(Sex = c(Female = "#F6D2E0", Male = "#C8E7F5"))  # pink/blue
  groups <- sort(unique(col_annotation[[1]]))
  colors <- c("#BF812D", "#35978F", "#C51B7D", "#7FBC41", "#762A83",
              "#E08214", "#542788", "#D6604D", "#4393C3", "#878787",
              "#1A1A1A", "#FFFFBF", "#9E0142", "#E41A1C", "#377EB8",
              "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628",
              "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#000000")
  colors <- colors[1:length(groups)]
  names(colors) <- groups
  ann_colors = list(colors)
  names(ann_colors) <- colnames(col_annotation)
  
  # Define how samples will be arranged in the heatmap.
  # Set cluster_cols=FALSE, if you want to arrange samples in specific order
  # Set cluster_cols=TRUE, if you want to arrange samples based on clustering
  mat <- mat[, rownames(col_annotation)]
  
  # List genes and samples you want to display in the plot
  display_row <- rownames(mat)      #c("SMAD6","EFNB2","CDX2)
  display_col <- colnames(mat)
  
  # Plot heatmap
  # NOTE: Error in check.length("fill"):'gpar' element 'fill' must not be length 0
  # This error means sample_annotation dataframe doesnt match with columns of mat
  # Try using mat = t(mat) to see if it fixes the error.
  # NOTE: If you set scale = none, then you are plotting exact progeny scores
  # NOTE: If you set scale = row, then colors do not represent progeny scores 
  # (i.e. actual activity) but relative activity
  pheatmap::pheatmap(mat = as.matrix(mat),
                     color = my_palette,
                     # breaks correspond to numerical ranges for the color palette's bins .i.e. 0 to length(my_palette)
                     breaks = c(seq(from=min(mat), to=0, length.out=ceiling(100/2) + 1), seq(from=max(mat)/100, to=max(mat), length.out=floor(100/2))),
                     # breaks = c(seq(-3, 0, length.out=ceiling(100/2) + 1), seq(max(mat)/100, 3, length.out=floor(100/2))), 
                     border_color = "white", #"grey60",
                     cellwidth = 15, 
                     cellheight = 5, 
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
                     annotation_colors = ann_colors, #dplyr::if_else(nrow(col_annotation)+nrow(row_annotation) > 0, ann_colors, c()),
                     annotation_legend = TRUE,
                     annotation_names_row = FALSE,
                     annotation_names_col = FALSE,
                     show_rownames = dplyr::if_else(nrow(mat)<80, TRUE, FALSE, missing = NULL), 
                     show_colnames = dplyr::if_else(ncol(mat)<50, TRUE, FALSE, missing = NULL),
                     fontsize = 8, 
                     fontsize_row = 8, 
                     fontsize_col = 8,
                     angle_col = c("270", "0", "45", "90", "315"),
                     fontsize_number = 0.8*fontsize, 
                     labels_row = c(display_row, rep(x="", times=nrow(mat) - length(display_row))),
                     labels_col = c(display_col, rep(x="", times=ncol(mat) - length(display_col))),
                     #width = 8.5, #ncol(mat)*0.075+5,
                     #height = 11,
                     filename = paste0(results_path, "Progeny_", names[i], ".pdf"))
}

for (i in 1:length(files1)){
  
  DEGs <- openxlsx::read.xlsx(xlsxFile = files1[i], rowNames=TRUE)
  
  DEGs <- DEGs %>%
    # dplyr::mutate(SYMBOL = base::gsub(pattern = "[- ]", replacement = ".", x = SYMBOL)) %>%
    # dplyr::rename_with(.fn = ~base::gsub(pattern = "-| ", replacement = ".", x = .x), .cols = everything()) %>%
    dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
    tibble::remove_rownames(.) %>%
    tibble::column_to_rownames(var = "SYMBOL") %>%
    dplyr::mutate(across(.cols = everything(), .fns = as.numeric)) %>%
    base::replace(is.na(.), 0) %>%
    dplyr::select(log2FoldChange)
  
  # Calculate progeny scores using DEGs
  progeny_zscores <- progeny::progeny(expr = as.matrix(DEGs), 
                                      scale = TRUE, 
                                      organism = dplyr::if_else(species=="Homo sapiens", "Human", "Mouse"), 
                                      top = 100, 
                                      perm = 10000,
                                      verbose = FALSE,
                                      z_scores = TRUE,
                                      get_nulldist = FALSE,
                                      assay_name = "RNA",
                                      return_assay = FALSE) %>%
    t() %>%
    data.frame() %>%
    dplyr::rename(NES = identity(1)) %>%
    tibble::rownames_to_column(var = "Pathway") %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(Pathway = factor(Pathway))
  
  # Plot NES plot
  ggplot2::ggplot(data = progeny_zscores,
                  mapping = aes(x = reorder(x = Pathway, X = NES), y = NES)) + 
    geom_bar(aes(fill = NES), stat = "identity") +
    scale_fill_gradient2(low = "darkblue", high = "indianred",
                         mid = "whitesmoke", midpoint = 0) +
    # viridis::scale_fill_viridis(name = "NES", alpha = 1, begin = 0, end = 1,
    #                             direction = 1, discrete = FALSE, option = "D") +
    theme_classic() +
    labs(x = "Pathways", y = "NES", title = stringr::str_wrap(paste0(""), 30)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10), #adjust size of x-axis text
          axis.text.y = element_text(size = 10),                                 #adjust size of y-axis text
          axis.title = element_text(size = 14),                                  #adjust size of axis label
          plot.title = element_text(hjust=0.5, size = 16, face="bold"))          #adjust size of graph title
  
  # Save the plot
  ggplot2::ggsave(filename = paste0("NES_", names[i], ".pdf"),
                  plot = last_plot(),
                  device = "pdf",
                  path = results_path,
                  scale = 1,
                  #width = 8.5,
                  #height = 11,
                  units = c("in"),	 
                  dpi = 600,
                  limitsize = TRUE,
                  bg = NULL)
} 

# # Based on NES plot, select top pathway and see what genes are being perturbed
# prog_matrix <- getModel(dplyr::if_else(species=="Homo sapiens", "Human", "Mouse"), top=100) %>% 
#   data.frame()  %>%
#   tibble::rownames_to_column("GeneID")
# 
# DEGs <- DEGs %>% 
#   data.frame() %>% 
#   tibble::rownames_to_column("GeneID")
# 
# scat_plots <- progeny::progenyScatter(df = DEGs, 
#                                       weight_matrix = prog_matrix, 
#                                       statName = "log2FC", 
#                                       verbose = FALSE)
# 
# plot(scat_plots[[1]]$`MAPK`)

# for (i in 1:14){
#  
#   #class(plot(scat_plots[[1]]$`MAPK`))
# 
#   ggplot2::ggsave(filename = paste0(names(scat_plots[[1]][i]), ".pdf"),
#                   plot = scat_plots[[1]][i],
#                   device = "pdf",
#                   path = results_path,
#                   scale = 1,
#                   #width = 8.5,
#                   #height = 11,
#                   units = c("in"),
#                   dpi = 600,
#                   limitsize = TRUE,
#                   bg = NULL)
# }

################################################################################

# Get top 100 significant genes per pathway
model_100 <- model %>%
  group_by(pathway) %>%
  slice_min(order_by = p.value, n = 100)

# Plot
ggplot2::ggplot(data=model_100, aes(x=weight, color=pathway, fill=pathway)) +
  geom_density() +
  theme(text = element_text(size=12)) +
  facet_wrap(~ pathway, scales='free') +
  xlab('scores') +
  ylab('densities') +
  theme_bw() +
  theme(legend.position = "none")

################################################################################

