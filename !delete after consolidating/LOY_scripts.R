
# #******************************************************************************#
# #                             STEP 1: TCGA HEATMAP                             #
# #******************************************************************************#
# 
# g1 <- c("AMELY","BPY2","BPY2B","BPY2C","CDY1","CDY1B","CDY2A","CDY2B","DAZ1",
#         "DAZ2","DAZ3","DAZ4","DDX3Y","EIF1AY","HSFY1","HSFY2","KDM5D","NLGN4Y",
#         "PCDH11Y","PRY","PRY2","PRYP3","RBMY1A1","RBMY1B","RBMY1D","RBMY1E",
#         "RBMY1F","RBMY1J","RPS4Y1","RPS4Y2","SRY","TBL1Y","TGIF2LY","TMSB4Y",
#         "TSPY1","TSPY10","TSPY2","TSPY3","TSPY4","TSPY8","TSPY9P","USP9Y","UTY",
#         "VCY","VCY1B","ZFY")
# 
# norm_counts <- read.xlsx("C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/TCGA/Results__id_Male_vs_Female_Normalized_Counts.xlsx")
# norm_counts <- norm_counts[,-1]
# norm_counts <- norm_counts %>%
#   dplyr::filter(SYMBOL %in% g$SYMBOL) %>%
#   dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
#   #dplyr::group_by(SYMBOL) %>%
#   #dplyr::summarise(across(.cols = everything(), .fns = sum)) %>%
#   #dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
#   tibble::column_to_rownames("SYMBOL") %>%
#   dplyr::mutate(across(.cols = everything(), .fns = as.numeric)) 
# 
# # Y_RNA and SNORA70 were repeated 756 and 27 times. So,we remove them.
# norm_counts <- norm_counts[!(rownames(norm_counts) %in% c("Y_RNA", "SNORA70")),]
# 
# meta_data <- read.xlsx("C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/TCGA/Metadata.xlsx")
# 
# mat <- norm_counts
# mat <- log(1+mat, base = 2)
# mat <- data.frame(t(scale(t(mat))))
# mat[is.na(mat)] <- 0
# 
# # Define column annotation
# col_annotation <- meta_data %>% dplyr::select(Sex)
# rownames(col_annotation) <- colnames(mat)
# 
# # Define row annotation
# row_annotation <- data.frame("Genes"= matrix(data="", nrow=nrow(mat), ncol=1))
# rownames(row_annotation) <- rownames(mat)
# 
# # Define colors for heatmap
# my_palette <- colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100)
# 
# # Define colors for row and column annotation
# groups <- sort(unique(col_annotation[[1]]))
# colors <- c("#BF812D", "#35978F", "#C51B7D", "#7FBC41", "#762A83",
#             "#E08214", "#542788", "#D6604D", "#4393C3", "#878787",
#             "#1A1A1A", "#FFFFBF", "#9E0142", "#E41A1C", "#377EB8",
#             "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628",
#             "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#000000")
# colors <- c("#74006F", "#FFC606")    # Female:purple, Male:Golds
# colors <- c("#F6D2E0", "#C8E7F5")    # Female:pink,   Male:Blue
# colors <- colors[1:length(groups)]
# names(colors) <- groups
# ann_colors = list(colors)
# names(ann_colors) <- colnames(col_annotation)
# 
# # # Define how samples will be arranged in the heatmap.
# # # Set cluster_cols=FALSE, if you want to arrange samples in specific order
# # # Set cluster_cols=TRUE, if you want to arrange samples based on clustering
# # # meta_data <- meta_data %>% dplyr::arrange(desc(gender))
# # # mat <- normalized_counts[, !!!!]
# 
# # List genes and samples you want to display in the plot
# display_row <- rownames(mat)      #c("SMAD6","EFNB2","CDX2)
# display_col <- colnames(mat)
# 
# pheatmap::pheatmap(mat = as.matrix(mat),
#                    color = my_palette, #c("white", "black"), #,
#                    # breaks correspond to numerical ranges for the color palette's bins .i.e. 0 to length(my_palette)
#                    #breaks = c(seq(from=min(mat), to=0, length.out=ceiling(100/2) + 1), seq(from=max(mat)/100, to=max(mat), length.out=floor(100/2))),
#                    #breaks = c(seq(-3, 0, length.out=ceiling(100/2) + 1), seq(max(mat)/100, 3, length.out=floor(100/2))),
#                    breaks = c(seq(0, 0.1, length.out=ceiling(100/2) + 1), seq(max(mat)/100, 3, length.out=floor(100/2))),
#                    #breaks = c(0, max(mat), length.out=2),
#                    border_color = "white", #"grey60",
#                    cellwidth = NA, 
#                    cellheight = NA, 
#                    scale = "none",   
#                    cluster_rows = TRUE,   #cluster the rows
#                    cluster_cols = TRUE,   #cluster the columns
#                    clustering_distance_rows = "euclidean",
#                    clustering_distance_cols = "euclidean",
#                    clustering_method = "complete",
#                    legend = TRUE, 
#                    legend_breaks = NA,
#                    legend_labels = NA, 
#                    annotation_row = row_annotation,
#                    annotation_col = col_annotation,
#                    annotation_colors = ann_colors,
#                    annotation_legend = TRUE,
#                    annotation_names_row = FALSE,
#                    annotation_names_col = FALSE,
#                    #show_rownames = dplyr::if_else(nrow(mat)<80, TRUE, FALSE, missing = NULL), 
#                    #show_colnames = dplyr::if_else(ncol(mat)<30, TRUE, FALSE, missing = NULL),
#                    fontsize = 8, 
#                    fontsize_row = 8, 
#                    fontsize_col = 8,
#                    angle_col = c("270", "0", "45", "90", "315"),
#                    fontsize_number = 0.8*fontsize, 
#                    labels_row = c(display_row, rep(x="", times=nrow(mat) - length(display_row))),
#                    labels_col = c(display_col, rep(x="", times=ncol(mat) - length(display_col))),
#                    height = 15,
#                    filename = paste0("C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/TCGA/Heatmap_Ygenes.pdf"))
# 
# # cluster and re-order rows
# rowclust = hclust(dist(mat))
# reordered = mat[rowclust$order,]
# 
# # cluster and re-order columns
# colclust = hclust(dist(t(mat)))
# reordered = reordered[, colclust$order]
# 
# # Save the clustered matrix
# wb <- openxlsx::createWorkbook()
# openxlsx::addWorksheet(wb, sheetName = "Heatmap_matrix")
# openxlsx::writeData(wb, sheet = "Heatmap_matrix", x = reordered, rowNames = TRUE)
# openxlsx::saveWorkbook(wb, 
#                        file = "C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/TCGA/Heatmap_matrix.xlsx", 
#                        overwrite = TRUE)

#******************************************************************************#
#                  STEP 3: HEATMAP Y GENES FROM SINGLE CELL DATA                  #
#******************************************************************************#
#!/usr/bin/env Rscript

for (proj in c("scRNASeq_BBN_C57B6", "scRNASeq_BBN_Rag", "scRNASeq_Chen",
               "scRNASeq_HRA003620", "scRNASeq_Jinfen", "scRNASeq_GSE222315")){
  
  source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")
  
  # Load the integrated seurat object
  integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn.rds"))
  
  # Get Y genes from ensembl
  annotations <- get_annotations(species)
  y_all <- annotations %>% 
    dplyr::filter(CHR == "Y", 
                  nchar(SYMBOL) > 0, 
                  !grepl("predicted|pseudo|rRNA", DESCRIPTION), 
                  !grepl("-ps|Rik", SYMBOL)) %>%
    dplyr::select(SYMBOL) %>%
    unlist(use.names=FALSE) %>%
    unique()
  
  y_18 <- c("DDX3Y", "EIF1AY", "HSFY2", "KDM5D", "UTY", "NLGN4Y", "PCDH11Y", 
            "RPS4Y1", "TBL1Y", "TMSB4Y", "USP9Y", "ZFY", "DAZ1", "DAZ2", "DAZ3",
            "DAZ4", "PRY2", "RBMY1A1")
  
  # Keep y genes common between seurat and ensembl
  y_seurat <- rownames(integrated_seurat@assays$RNA$counts)
  y_all <- y_seurat[tolower(y_seurat) %in% tolower(y_all)]
  y_18 <- y_seurat[tolower(y_seurat) %in% tolower(y_18)]
  
  # # Add Y score
  # DefaultAssay(integrated_seurat) <- "RNA"
  # x <- UCell::AddModuleScore_UCell(obj = integrated_seurat, 
  #                                  features = list(Yscore = y_all),
  #                                  assay = "RNA",
  #                                  slot = "data",
  #                                  ties.method = "average",
  #                                  name = "Yall_UCell")
  # x <- UCell::AddModuleScore_UCell(obj = x, 
  #                                  features = list(Yscore = y_18),
  #                                  assay = "RNA",
  #                                  slot = "data",
  #                                  ties.method = "average",
  #                                  name = "Y18_UCell")
  # x <- Seurat::AddModuleScore(obj = x, 
  #                                  features = list(y_all),
  #                                  assay = "RNA",
  #                                  slot = "data",
  #                                  name = "Yall")
  # x <- Seurat::AddModuleScore(obj = x, 
  #                             features = list(y_18),
  #                             assay = "RNA",
  #                             slot = "data",
  #                             name = "Y18")
  
  #y_genes <- setdiff(y_genes, c("Gm47283","TTTY14","PCDH11Y"))
  
  celltypes <- unique(integrated_seurat@meta.data$cell_type)
  for (c in celltypes){
    
    seurat_obj <- subset(x=integrated_seurat,
                         cell_type == c)
    
    # Extract expression inf, keep only Y genes that have expression
    df <- seurat_obj@assays$RNA$data
    df_all <- df[y_all,]
    df_all <- df_all[rowSums(df_all) != 0,]
    
    df_all <- df_all %>% 
      t() %>% 
      data.frame() %>% 
      tibble::rownames_to_column("Sample") %>%
      dplyr::mutate(Sample = gsub("_.*", "", Sample)) %>%
      dplyr::group_by(Sample) %>%
      dplyr::summarize(across(.cols=everything(), .fns=mean)) %>%
      tibble::column_to_rownames("Sample") %>%
      t() %>%
      data.frame() %>%
      tibble::rownames_to_column("SYMBOL")
    
    df_18 <- df_all %>% 
      dplyr::filter(SYMBOL %in% y_18)
    
    for (mat in c("df_all", "df_18")){
      
      normalized_counts <- get(mat)
      plot_genes <- normalized_counts$SYMBOL
      disp_genes <- plot_genes
      file_suffix <- paste0(proj, "_", c)
      file_format <- "jpeg" 
      results_path <- ""
      bar_width <- 5
      bar_height <- 5
      expr_legend <- TRUE
      metadata_column <- integrated_seurat@meta.data %>% 
        dplyr::select(Sample, Sex) %>% 
        dplyr::distinct_at("Sample", .keep_all = TRUE)
      metadata_row <- data.frame(SYMBOL = normalized_counts$SYMBOL)
      
      perform_log_transform <- FALSE
      perform_scaling <- TRUE
      anno_columns <- "Sex"
      anno_rows <- NA
      columns <- "Sample"
      color_by_cols <- TRUE
      color_by_rows <- FALSE
      gaps_in_row <- FALSE
      gaps_in_col <- FALSE
      row_clustering_alphabetical <- FALSE
      col_clustering_alphabetical <- TRUE
      row_clustering <- TRUE
      col_clustering <- FALSE
      col_clustering_within_group <- FALSE
      row_clustering_within_group <- FALSE
      
      plot_heatmap(normalized_counts, metadata_column, metadata_row, 
                   plot_genes, disp_genes, file_suffix, file_format, 
                   results_path, bar_width, bar_height, expr_legend)
      
      
    }
  }
}

# Convert count matrix to binary format [1=Expressed, 0=Not expressed]
#df <- as.data.frame(df) %>% replace(.> 0, 1)
df[df>0] <- 1

gene_count_per_cell <- data.frame(counts = colSums(df)) %>%
  tibble::rownames_to_column("Cell") %>%
  dplyr::mutate(Cell = gsub(pattern= "_.*" , replacement="", x=Cell)) %>%
  dplyr::group_by(Cell, counts) %>%
  dplyr::summarise(n = n()) %>%
  #tidyr::pivot_wider(id_cols=Cell, names_from=counts, values_from=n) %>%
  dplyr::rename(Sample=Cell, n_gene=counts) %>%
  dplyr::mutate(Percent = 100*n/sum(n))

ggplot(data = gene_count_per_cell, aes(x=n_gene, y=Percent, group=Sample, fill=Sample)) +
  geom_col() +
  my_theme +
  facet_wrap(~Sample, nrow=1) +
  scale_x_continuous(breaks=seq(0,10,1))
theme(legend.position="right",
      panel.spacing = unit(0.1, "lines"),
      axis.ticks.x=element_blank())

ggsave(paste0(celltype, "_ngene_bar.tiff"))

gene_count_per_cell <- data.frame(counts = colSums(df)) %>%
  tibble::rownames_to_column("Cell") %>%
  dplyr::mutate(Cell = gsub(pattern= "_.*" , replacement="", x=Cell)) %>%
  dplyr::group_by(Cell, counts) %>%
  dplyr::rename(Sample=Cell, n_gene=counts)

ggplot(data = gene_count_per_cell, aes(x=n_gene, group=Sample, fill=Sample)) +
  geom_density(adjust=1.5) +
  my_theme +
  facet_wrap(~Sample, nrow=1) +
  scale_x_continuous(breaks=seq(0,10,1))
theme(legend.position="right",
      panel.spacing = unit(0.1, "lines"),
      axis.ticks.x=element_blank())

ggsave(paste0(celltype, "_ngene_density.tiff"))

x <- Seurat::AddModuleScore(object = integrated_seurat,
                            features = list(rownames(df)),
                            name = "Y_score",
                            pool = NULL,
                            nbin = 24,
                            ctrl = 5, #100,
                            k = FALSE,
                            assay = "RNA",
                            seed = 1,
                            search = FALSE)
y_score <- x@meta.data %>% 
  dplyr::select(Sample, Y_score1) %>%
  tibble::remove_rownames()

ggplot(data = y_score, aes(x=Y_score1, group=Sample, fill=Sample)) +
  geom_density() +
  my_theme +
  facet_wrap(~Sample, nrow=1) +
  scale_x_continuous(breaks=seq(0,10,1))
theme(legend.position="right",
      panel.spacing = unit(0.1, "lines"),
      axis.ticks.x=element_blank())

ggsave(paste0(celltype, "_yscore.tiff"))

#BiocManager::install("UCell")   
DefaultAssay(integrated_seurat) <- "RNA"
x <- UCell::AddModuleScore_UCell(obj = integrated_seurat, 
                                 features = list(Yscore = rownames(df)),
                                 maxRank = 1500,
                                 chunk.size = 1000,
                                 BPPARAM = NULL,
                                 ncores = 1,
                                 storeRanks = FALSE,
                                 w_neg = 1,
                                 assay = NULL,
                                 slot = "data",
                                 ties.method = "average",
                                 force.gc = FALSE,
                                 name = "_UCell")

y_score_ucell <- x@meta.data %>% 
  dplyr::select(Sample, Yscore_UCell) %>%
  tibble::remove_rownames()

ggplot(data = y_score_ucell, aes(x=Yscore_UCell, group=Sample, fill=Sample)) +
  geom_density() +
  my_theme +
  facet_wrap(~Sample, nrow=1) +
  scale_x_continuous(breaks=seq(0,10,1)) + 
  ylim(0,40)
theme(legend.position="right",
      panel.spacing = unit(0.1, "lines"),
      axis.ticks.x=element_blank())

ggsave(paste0(celltype, "_yscore_ucell.tiff"))

}

# # Prepare dataframe for plotting
# gene_count_per_cell <- data.frame(counts = colSums(df)) %>%
#   tibble::rownames_to_column("Cell") %>%
#   #dplyr::left_join(integrated_seurat@meta.data %>% dplyr::select(barcode, Patient), by=c("Cell"="barcode")) %>%
#   dplyr::select(Patient, counts) %>%
#   dplyr::group_by(Patient, counts) %>%
#   dplyr::summarise(n = n()) %>%
#   tidyr::pivot_wider(id_cols=Patient, names_from=counts, values_from=n) %>%
#   dplyr::rename(Sample=Patient)

#gsub(pattern= "^.*?-" , replacement="", x="AAAGATGGTGTCCTCT-1-24-0-0")

data <- gene_count_per_cell
file_suffix <- paste0(celltype, "_gene_number")
label_percent <- "TRUE"
already_percent <- "FALSE"
#results_path <- seurat_results

plot_param <- list(
  
  # column to be plotted on x axis
  "data_x" = c("nGene"),
  
  # Define title of x axis
  "title_x" = c(),
  
  # column to be plotted on y axis
  "data_y" = c("counts"),
  
  # Define title of y axis
  "title_y" = c("Number of Cells"),
  
  # column to be used for filling bar colors
  "data_fill" = c(""),
  
  # column to be used for coloring the dots
  "data_color" = c(""),
  
  # column to be used for determining the size of dots
  "data_size" = c(),
  
  # Define title of legend
  "title_legend_fill" = c("Number of Genes Expressed"),
  "title_legend_color" = c(),
  "title_legend_size" = c(),
  
  # Define title of plot
  "title_plot" = c("Gene Expression Frequency")
)

plot_param <- list(
  
  # column to be plotted on x axis
  "data_x" = c("Sample"),
  
  # Define title of x axis
  "title_x" = c(),
  
  # column to be plotted on y axis
  "data_y" = c(),
  
  # Define title of y axis
  "title_y" = c("Percent of Cells"),
  
  # column to be used for filling bar colors
  "data_fill" = c("n_gene"),
  
  # column to be used for coloring the dots
  "data_color" = c(),
  
  # column to be used for determining the size of dots
  "data_size" = c(),
  
  # Define title of legend
  "title_legend_fill" = c("Number of Genes Expressed"),
  "title_legend_color" = c(),
  "title_legend_size" = c(),
  
  # Define title of plot
  "title_plot" = c("Gene Expression Frequency")
)
plot_stackedbar(data, plot_param, label_percent)

cell_count_per_gene <- df %>%
  t() %>%
  data.frame() %>%
  tibble::rownames_to_column("Cell") %>%
  dplyr::left_join(integrated_seurat@meta.data %>% dplyr::select(barcode, Patient), by=c("Cell"="barcode")) %>%
  dplyr::select(Patient, everything(), -Cell) %>%
  dplyr::group_by(Patient) %>%
  dplyr::summarize(across(.cols = everything(), .fns= ~ mean(.x*100))) %>%
  tibble::column_to_rownames("Patient")

#dplyr::mutate(Cell = gsub(pattern= "_.*" , replacement="", x=Cell)) %>%
dplyr::mutate(Cell = gsub(pattern= "[^_]*$" , replacement="", x=Cell)) %>%
  dplyr::group_by(Cell) %>%
  dplyr::summarize(across(.cols = everything(), .fns= ~ mean(.x*100))) %>%
  tibble::column_to_rownames("Cell")

cell_count_per_gene <- cell_count_per_gene[,colSums(cell_count_per_gene) > 1]

cell_count_per_gene <- cell_count_per_gene %>%
  tibble::rownames_to_column("Sample") %>%
  tidyr::pivot_longer(cols=!Sample, names_to="n_gene", values_to="Percent")

data <- cell_count_per_gene
file_suffix <- paste0(celltype, "_gene_composition")
label_percent <- "TRUE"
already_percent <- "TRUE"
results_path <- seurat_results
plot_param <- list(
  
  # column to be plotted on x axis
  "data_x" = c("Sample"),
  
  # Define title of x axis
  "title_x" = c(),  
  
  # column to be plotted on y axis
  "data_y" = c(), 
  
  # Define title of y axis
  "title_y" = c("Percent of Cells"),
  
  # column to be used for filling bar colors
  "data_fill" = c("n_gene"),
  
  # column to be used for coloring the dots
  "data_color" = c(), 
  
  # column to be used for determining the size of dots
  "data_size" = c(), 
  
  # Define title of legend
  "title_legend_fill" = c("Genes Expressed"),
  "title_legend_color" = c(),
  "title_legend_size" = c(), 
  
  # Define title of plot
  "title_plot" = c("Relative Gene Composition")
)
plot_stackedbar(data, plot_param, label_percent)
}
}

#******************************************************************************#
#                        Find DEGs between Y- vs Y+                       #                       
#******************************************************************************#

# Find Yneg cells
celltype <- "Epithelial"
celltype <- "Fibroblast"
celltype <- "Lymphoid"
celltype <- "Myeloid"

integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn",
                                    dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))

integrated_seurat <- subset(x = integrated_seurat,
                            subset = (cell_class %in% c("Mixed")),
                            invert = TRUE)

annotations <- get_annotations(species)

# Get Y genes from ensembl
y_ensembl <- annotations %>% 
  dplyr::filter(CHR == "Y") %>%
  dplyr::filter(nchar(SYMBOL) > 0)

y_genes <- base::intersect(unique(y_ensembl$SYMBOL), rownames(integrated_seurat@assays$RNA@data))
y_genes <- setdiff(y_genes, "Gm47283")

# Extract expression info
df <- integrated_seurat@assays$RNA@data

# Keep only expression of Y genes
df <- df[y_genes,]

# Remove genes which have no expression
df <- df[rowSums(df) != 0,]

# # Convert count matrix to binary format [1=Expressed, 0=Not expressed]
# df <- as.data.frame(df) %>% replace(.> 0, 1)
# 
# # Prepare dataframe for plotting
# gene_count_per_cell <- data.frame(counts = colSums(df)) %>%
#   tibble::rownames_to_column("Cell") %>%
#   dplyr::filter(counts == 0)

# Identify cells with zero Y expression
gene_count_per_cell <- as.data.frame(colSums(df[, colSums(df) == 0])) %>%
  tibble::rownames_to_column("Cell") 

meta_data <-integrated_seurat@meta.data %>% 
  dplyr::mutate(Ystatus = dplyr::case_when(Cell %in% gene_count_per_cell$Cell & Sex == "Male" ~ "Yneg", 
                                           Sex == "Male" ~ "Ypos", 
                                           TRUE ~ "Female"))

integrated_seurat@meta.data <- meta_data

integrated_seurat <- subset(x=integrated_seurat,
                            Sex == "Male" &
                              Condition == "Tumor")
integrated_seurat <- subset(x=integrated_seurat,
                            cell_type == "Lymphoid - T")
integrated_seurat <- subset(x=integrated_seurat,
                            sub_type == "Myeloid - Macrophage")

# Option 1:
Idents(integrated_seurat) <- "Ystatus"
DefaultAssay(object = integrated_seurat) <- "RNA"
DE_genes <- FindMarkers(object = integrated_seurat,
                        slot = "data",
                        ident.1 = "Yneg",
                        ident.2 = "Ypos",
                        logfc.threshold = 0.05)

# integrated_seurat@assays$RNA@counts <- as.matrix(integrated_seurat@assays$RNA@counts)+1
# DE_genes <- FindMarkers(object = integrated_seurat,
#                         slot = "counts",
#                         ident.1 = "Ypos",
#                         ident.2 = "Yneg", 
#                         test.use = "DESeq2")

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb = wb, sheetName = "DEGs")
openxlsx::writeData(wb = wb, sheet = "DEGs", x = DE_genes, rowNames = TRUE)
openxlsx::saveWorkbook(wb = wb, file = paste0("FindMarker_Results_", celltype, "_.xlsx"), overwrite = TRUE)

# option 2:
# Subset metadata
meta_data_yneg <-integrated_seurat@meta.data %>%
  dplyr::filter(Ystatus == "Yneg") %>%
  dplyr::distinct(Sample, .keep_all = TRUE) %>%
  tibble::remove_rownames() %>%
  dplyr::select(Sample, Sex, Condition, Ystatus) %>%
  dplyr::mutate(Batch = 1,
                Sample = paste0(Sample, "_", Ystatus))

meta_data_ypos <-integrated_seurat@meta.data %>%
  dplyr::filter(Ystatus == "Ypos") %>%
  dplyr::distinct(Sample, .keep_all = TRUE) %>%
  tibble::remove_rownames() %>%
  dplyr::select(Sample, Sex, Condition, Ystatus) %>%
  dplyr::mutate(Batch = 1,
                Sample = paste0(Sample, "_", Ystatus))

subset_seurat_yneg <- subset(x = integrated_seurat,
                             subset = Cell %in% gene_count_per_cell$Cell)

subset_seurat_ypos <- subset(x = integrated_seurat,
                             subset = Cell %in% gene_count_per_cell$Cell,
                             invert=TRUE)

# The read data will have "the reads of all cells belonging to a single 
# sample" merged together in each column. First, create a list of samples
samples_yneg <- subset_seurat_yneg@meta.data %>% 
  dplyr::select(Sample) %>% 
  unlist(., use.names=FALSE) %>% 
  unique()

samples_ypos <- subset_seurat_ypos@meta.data %>% 
  dplyr::select(Sample) %>% 
  unlist(., use.names=FALSE) %>% 
  unique()

Variable <- "Ystatus"
Comparisons <- list(Target = c("Yneg"),
                    Reference = c("Ypos"))    
Variable2 <- "Condition"
Variable2_value <- "Tumor"

# Second, create an empty dataframe with rows=genes and columns=samples
read_data_yneg <- data.frame(matrix(NA, nrow = nrow(subset_seurat_yneg@assays$RNA@counts), ncol = length(samples_yneg)))
rownames(read_data_yneg) <- rownames(subset_seurat_yneg@assays$RNA@counts)
colnames(read_data_yneg) <- samples_yneg

read_data_ypos <- data.frame(matrix(NA, nrow = nrow(subset_seurat_ypos@assays$RNA@counts), ncol = length(samples_ypos)))
rownames(read_data_ypos) <- rownames(subset_seurat_ypos@assays$RNA@counts)
colnames(read_data_ypos) <- samples_ypos

# Thirdly, we will add row-wise, the counts of each gene for each sample
for(i in samples_yneg){
  
  # Create a list of cells for each sample
  cells_subset <- rownames(subset_seurat_yneg@meta.data %>% dplyr::filter(Sample == i))
  
  # Use data.frame to convert "." in sparse matrix to "0"
  subset <- data.frame(subset_seurat_yneg@assays$RNA@counts[,cells_subset])
  read_data_yneg[,i]  <- rowSums(subset)
}

for(i in samples_ypos){
  
  # Create a list of cells for each sample
  cells_subset <- rownames(subset_seurat_ypos@meta.data %>% dplyr::filter(Sample == i))
  
  # Use data.frame to convert "." in sparse matrix to "0"
  subset <- data.frame(subset_seurat_ypos@assays$RNA@counts[,cells_subset])
  read_data_ypos[,i]  <- rowSums(subset)
}

colnames(read_data_yneg) <- paste0(colnames(read_data_yneg), "_Yneg")
colnames(read_data_ypos) <- paste0(colnames(read_data_ypos), "_Ypos")

read_data_yneg <- read_data_yneg %>% 
  tibble::rownames_to_column("SYMBOL")

read_data_ypos <- read_data_ypos %>% 
  tibble::rownames_to_column("SYMBOL")

meta_data <- dplyr::bind_rows(meta_data_yneg, meta_data_ypos)
read_data <- dplyr::left_join(read_data_yneg, read_data_ypos, by=c("SYMBOL"="SYMBOL"))

results_path <- "/hpc/home/kailasamms/"
analyze_DESeq2(meta_data, read_data, celltype)


#******************************************************************************#
#Fibroblast subtype and LOY

integrated_seurat <- subset(x= integrated_seurat,
                            Condition == "Tumor")

# Classify fibroblast clusters in C57B6 mice
feature_subset <- c("Pdgfrb", "Pde1a", "Gjc1", "Mcam", "Angpt2", "Gucy1a2")
feature_subset <- c("Fap","Postn","Col1a1","Col3a1","Actb")
feature_subset <- c("Pdpn", "Robo2", "Bmp5", "Pdgfc", "Vcan", "Thbs1","Fn1")
feature_subset <- c("Tgfbr1","Tagln","Acta2","Smtn","Des","Actg2")

x <- Seurat::AddModuleScore(object = integrated_seurat,
                            features = list(feature_subset),
                            name = make.names(i),
                            pool = NULL,
                            nbin = 24,
                            ctrl = 5, #100,
                            k = FALSE,
                            assay = "RNA",
                            seed = 1,
                            search = FALSE)

Seurat::FeaturePlot(object = x,
                    slot = "data",
                    features = paste0(make.names(i),1),
                    cols =  c("grey", viridis(n = 10, option = "C", direction = -1)),
                    pt.size = 0.4,
                    order = TRUE,
                    min.cutoff = 'q10',
                    reduction = "umap",
                    label = TRUE,
                    combine = TRUE,
                    raster = FALSE)

ggsave("pdgfrb.jpg")
ggsave("fap.jpg")
ggsave("pdpn.jpg")
ggsave("sma.jpg")

#*******************
pdgfrb <- c("8","10","14","15")
fap <- c("5","6","9","13","16","24","26","28")
pdpn <- c("11","18","19","23","25")
sma <- c("4","27")

celltype <- "Fibroblasts"

integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn",
                                    dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))

integrated_seurat <- subset(x = integrated_seurat,
                            subset = (cell_class %in% c("Mixed", "Unclassified")),
                            invert = TRUE)

#annotations <- get_annotations(species)

# Get Y genes from ensembl
y_ensembl <- annotations %>% 
  dplyr::filter(CHR == "Y") %>%
  dplyr::filter(nchar(SYMBOL) > 0)

y_genes <- base::intersect(unique(y_ensembl$SYMBOL), rownames(integrated_seurat@assays$RNA@data))
y_genes <- setdiff(y_genes, "Gm47283")

# Extract expression info
df <- integrated_seurat@assays$RNA@data

# Keep only expression of Y genes
df <- df[y_genes,]

# Remove genes which have no expression
df <- df[rowSums(df) != 0,]

# Convert count matrix to binary format [1=Expressed, 0=Not expressed]
df <- as.data.frame(df) %>% replace(.> 0, 1)

# Prepare dataframe for plotting
gene_count_per_cell <- data.frame(counts = colSums(df)) %>%
  tibble::rownames_to_column("Cell") %>%
  dplyr::filter(counts == 0)

meta_data <-integrated_seurat@meta.data %>% 
  dplyr::mutate(Condition = dplyr::case_when(Condition == "BBN" ~ "Tumor",
                                             Condition == "Tumor" ~ "Tumor",
                                             Condition == "Normal" ~ "Normal"), 
                Ystatus = dplyr::case_when(Cell %in% gene_count_per_cell$Cell & Sex == "Male" ~ "Yneg", 
                                           Sex == "Male" ~ "Ypos", 
                                           TRUE ~ "Female"),
                Kenny_group = dplyr::case_when(integrated_snn_res.1.4 %in% pdgfrb ~ "Pdgfrb",
                                               integrated_snn_res.1.4 %in% fap ~ "Fap",
                                               integrated_snn_res.1.4 %in% pdpn ~ "Pdpn",
                                               integrated_snn_res.1.4 %in% sma ~ "Sma",
                                               TRUE ~ "Others"))

integrated_seurat@meta.data <- meta_data
base_seurat <- integrated_seurat
goi <- c("Itga5","Cd44","Tgfbr1","Cntn1","Nrp1","Fap","Itga1","Itga3","Sdc2",
         "Flt1","Pdgfrb","Plxnd1","Itga2","Itgb5","Lrp6","Egfr","Fgfr1","Erap1",
         "Robo2","Axl","ABca1","Smad3","Tgfbr2","Itga8")

# Perform DEseq for each subtype of fibroblasts
for (t in c("Pdgfrb", "Fap", "Pdpn", "Sma")){
  
  integrated_seurat <- subset(x=base_seurat,
                              Kenny_group == t &
                                Sex == "Male" & 
                                Condition == "Tumor")
  
  # Option 1:
  Idents(integrated_seurat) <- "Ystatus"
  DefaultAssay(object = integrated_seurat) <- "RNA"
  DE_genes <- FindMarkers(object = integrated_seurat,
                          slot = "data",
                          ident.1 = "Yneg",
                          ident.2 = "Ypos")
  
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb = wb, sheetName = "DEGs")
  openxlsx::writeData(wb = wb, sheet = "DEGs", x = DE_genes, rowNames = TRUE)
  openxlsx::saveWorkbook(wb = wb, file = paste0(results_path, "FindMarker_Results_", t, "_.xlsx"), overwrite = TRUE)
  
  # Option 2:
  meta_data_yneg <-integrated_seurat@meta.data %>%
    dplyr::filter(Ystatus == "Yneg") %>%
    dplyr::distinct(Sample, .keep_all = TRUE) %>%
    tibble::remove_rownames() %>%
    dplyr::select(Sample, Sex, Condition, Ystatus) %>%
    dplyr::mutate(Batch = 1,
                  Sample = paste0(Sample, "_", Ystatus))
  
  meta_data_ypos <-integrated_seurat@meta.data %>%
    dplyr::filter(Ystatus == "Ypos") %>%
    dplyr::distinct(Sample, .keep_all = TRUE) %>%
    tibble::remove_rownames() %>%
    dplyr::select(Sample, Sex, Condition, Ystatus) %>%
    dplyr::mutate(Batch = 1,
                  Sample = paste0(Sample, "_", Ystatus))
  
  subset_seurat_yneg <- subset(x = integrated_seurat,
                               subset = Cell %in% gene_count_per_cell$Cell)
  
  subset_seurat_ypos <- subset(x = integrated_seurat,
                               subset = Cell %in% gene_count_per_cell$Cell,
                               invert=TRUE)
  
  # The read data will have "the reads of all cells belonging to a single 
  # sample" merged together in each column. First, create a list of samples
  samples_yneg <- subset_seurat_yneg@meta.data %>% 
    dplyr::select(Sample) %>% 
    unlist(., use.names=FALSE) %>% 
    unique()
  
  samples_ypos <- subset_seurat_ypos@meta.data %>% 
    dplyr::select(Sample) %>% 
    unlist(., use.names=FALSE) %>% 
    unique()
  
  # Second, create an empty dataframe with rows=genes and columns=samples
  read_data_yneg <- data.frame(matrix(NA, nrow = nrow(subset_seurat_yneg@assays$RNA@counts), ncol = nrow(meta_data_yneg)))
  rownames(read_data_yneg) <- rownames(subset_seurat_yneg@assays$RNA@counts)
  colnames(read_data_yneg) <- samples_yneg
  
  read_data_ypos <- data.frame(matrix(NA, nrow = nrow(subset_seurat_ypos@assays$RNA@counts), ncol = nrow(meta_data_ypos)))
  rownames(read_data_ypos) <- rownames(subset_seurat_ypos@assays$RNA@counts)
  colnames(read_data_ypos) <- samples_ypos
  
  # Thirdly, we will add row-wise, the counts of each gene for each sample
  for(i in samples_yneg){
    
    # Create a list of cells for each sample
    cells_subset <- rownames(subset_seurat_yneg@meta.data %>% dplyr::filter(Sample == i))
    
    # Use data.frame to convert "." in sparse matrix to "0"
    subset <- data.frame(subset_seurat_yneg@assays$RNA@counts[,cells_subset])
    read_data_yneg[,i]  <- rowSums(subset)
  }
  
  for(i in samples_ypos){
    
    # Create a list of cells for each sample
    cells_subset <- rownames(subset_seurat_ypos@meta.data %>% dplyr::filter(Sample == i))
    
    # Use data.frame to convert "." in sparse matrix to "0"
    subset <- data.frame(subset_seurat_ypos@assays$RNA@counts[,cells_subset])
    read_data_ypos[,i]  <- rowSums(subset)
  }
  
  colnames(read_data_yneg) <- paste0(colnames(read_data_yneg), "_Yneg")
  colnames(read_data_ypos) <- paste0(colnames(read_data_ypos), "_Ypos")
  
  read_data_yneg <- read_data_yneg %>% 
    tibble::rownames_to_column("SYMBOL")
  
  read_data_ypos <- read_data_ypos %>% 
    tibble::rownames_to_column("SYMBOL")
  
  meta_data <- dplyr::bind_rows(meta_data_yneg, meta_data_ypos)
  read_data <- dplyr::left_join(read_data_yneg, read_data_ypos, by=c("SYMBOL"="SYMBOL"))
  
  results_path <- "/hpc/home/kailasamms/"
  analyze_DESeq2(meta_data, read_data, t)
}

#******************************************************************************#
# Plot UMAP of Yneg vs Ypos cells

for (celltype in c("Epithelial", "Fibroblasts", "Myeloid", "Lymphoid")){
  
  integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn",
                                      dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  integrated_seurat <- subset(x = integrated_seurat,
                              subset = (cell_class %in% c("Mixed", "Unclassified")),
                              invert = TRUE)
  
  #annotations <- get_annotations(species)
  
  # Get Y genes from ensembl
  y_ensembl <- annotations %>% 
    dplyr::filter(CHR == "Y") %>%
    dplyr::filter(nchar(SYMBOL) > 0)
  
  y_genes <- base::intersect(unique(y_ensembl$SYMBOL), rownames(integrated_seurat@assays$RNA@data))
  y_genes <- setdiff(y_genes, "Gm47283")
  
  # Extract expression info
  df <- integrated_seurat@assays$RNA@data
  
  # Keep only expression of Y genes
  df <- df[y_genes,]
  
  # Remove genes which have no expression
  df <- df[rowSums(df) != 0,]
  
  # Convert count matrix to binary format [1=Expressed, 0=Not expressed]
  df <- as.data.frame(df) %>% replace(.> 0, 1)
  
  # Prepare dataframe for plotting
  gene_count_per_cell <- data.frame(counts = colSums(df)) %>%
    tibble::rownames_to_column("Cell") %>%
    dplyr::filter(counts == 0)
  
  meta_data <-integrated_seurat@meta.data %>% 
    dplyr::mutate(Condition = dplyr::case_when(Condition == "BBN" ~ "Tumor",
                                               Condition == "Tumor" ~ "Tumor",
                                               Condition == "Normal" ~ "Normal"), 
                  Ystatus = dplyr::case_when(Cell %in% gene_count_per_cell$Cell & Sex == "Male" ~ "Yneg", 
                                             Sex == "Male" ~ "Ypos", 
                                             TRUE ~ "Female"),
                  id = paste0(Condition,"_", Sex, "_", Ystatus))
  
  integrated_seurat@meta.data <- meta_data
  
  # Set identity to an existing column in meta data
  Idents(object = integrated_seurat) <- "id"
  split_by <- NULL
  
  Seurat::DimPlot(object = subset(integrated_seurat, Sex=="Male" & Condition == "Tumor"),
                  reduction = "umap",
                  split.by = split_by,
                  label = FALSE,
                  raster = FALSE,
                  ncol = 1,
                  combine = TRUE)
  ggsave(paste0("Ystatus_", celltype, ".jpg"))
}

#******************************************************************************#

# Correlation of y score of subtypes
# main_df <- data.frame("sub_type"="", "mean_LOY_score"=0, "n"=0)
main_df <- data.frame(matrix(ncol = 1+length(y_genes), nrow = 1))
colnames(main_df) <- c("sub_type", y_genes)
main_df[1,2:(1+length(y_genes))] <- c(rep(x=0,times=length(y_genes)))
main_df[1,1] <- "None"


for (celltype in c("Epithelial", "Fibroblasts", "Myeloid", "Lymphoid")){
  
  integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn",
                                      dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  integrated_seurat <- subset(x = integrated_seurat,
                              subset = (cell_class %in% c("Mixed", "Unclassified") | Sex == "Female" | Condition == "Normal"),
                              invert = TRUE)
  
  #annotations <- get_annotations(species)
  
  # Get Y genes from ensembl
  y_ensembl <- annotations %>% 
    dplyr::filter(CHR == "Y") %>%
    dplyr::filter(nchar(SYMBOL) > 0)
  
  y_genes <- base::intersect(unique(y_ensembl$SYMBOL), rownames(integrated_seurat@assays$RNA@data))
  y_genes <- setdiff(y_genes, "Gm47283")
  
  # Extract expression info
  df <- integrated_seurat@assays$RNA@data
  
  # Keep only expression of Y genes
  df <- df[y_genes,]
  
  # Remove genes which have no expression
  #df <- df[rowSums(df) != 0,]
  
  # Convert count matrix to binary format [1=Expressed, 0=Not expressed]
  #df <- as.data.frame(df) %>% replace(.> 0, 1)
  
  #format df
  df <- t(df) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("Cell") %>% 
    left_join(integrated_seurat@meta.data %>% dplyr::select(Cell, sub_type), by=c("Cell"="Cell")) %>%
    dplyr::select(sub_type, everything(), -Cell) %>%
    dplyr::group_by(sub_type) %>%
    dplyr::mutate(across(.cols = everything(), .fns = ~mean(.x, na.rm=TRUE))) %>%
    dplyr::distinct()
  
  
  main_df <- dplyr::bind_rows(main_df, df)
}

main_df <- main_df[-1,] %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("sub_type")

# Remove genes with no expression
main_df <- main_df[, colSums(main_df) != 0]

# Remove genes with Gm
main_df <- main_df[, !grepl(pattern ="Gm", x=colnames(main_df))]

# Remove genes which have expresssion in only few subtypes
remove_cols <- c()
for (i in 1:ncol(main_df)){
  if(sum(unique(main_df[,i]) > 0) <= 2){
    remove_cols <- c(remove_cols, colnames(main_df)[i])
  }
}

main_df <- main_df %>% dplyr::select(everything(), -remove_cols)

# Transpose and calculate correlation
correlation <- t(main_df) %>% cor(.)
pheatmap(correlation)


#******************************************************************************#
#                    Plot heatmap of X and Y genes in epithelial cells         #
#******************************************************************************#

x_ensembl <- annotations%>% 
  dplyr::filter(CHR == "X")

# Get Y genes from ensembl
y_ensembl <- annotations %>% 
  dplyr::filter(CHR == "Y")

celltype <- "Epithelial"
integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn", 
                                    dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))

file_suffix <- NULL
columns <- "Cell"  # which column of metadata to plot as heatmap columns
already_log <- TRUE
already_scaled <- FALSE
row_clustering <- TRUE    # Usually TRUE
col_clustering <- FALSE   # Usually FALSE
row_clustering_alphabetical <- FALSE
col_clustering_alphabetical <- FALSE
gaps_in_col <- FALSE
gaps_column <- "Sample"
anno_columns <- c("Sample", "Sex", "Condition")
my_palette <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(200)
my_palette <- viridis(200)
my_palette <- c(colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100)[1:49], "#FFFFFF",
                colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100)[50:99])
results_path <- "/hpc/home/kailasamms/"

metadata <- integrated_seurat@meta.data %>% 
  dplyr::select(Cell, Sample, Condition, Sex) %>%
  #dplyr::distinct_at("Sample", .keep_all = TRUE) %>%
  tibble::remove_rownames()

normalized_counts <- integrated_seurat@assays$RNA@data %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Cell") %>%
  dplyr::mutate(Sample = gsub(pattern="_.*", replacement="", x=Cell),
                Expr = Uty + Kdm5d + Ddx3y + Eif2s3y) %>%
  dplyr::arrange(Sample, Expr) %>%
  dplyr::select(everything(), -c(Expr, Sample)) %>%
  tibble::column_to_rownames("Cell") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("SYMBOL")

# Keep x and y genes common between seurat and ensembl
y_genes <- base::intersect(unique(y_ensembl$SYMBOL), normalized_counts$SYMBOL)
x_genes <- base::intersect(unique(x_ensembl$SYMBOL), normalized_counts$SYMBOL)
plot_genes <- c(y_genes, x_genes)
plot_genes <- plot_genes[!grepl(pattern="^Gm|Rik$", x=plot_genes)]

# # Remove genes that are expressed in less than 10% of cells
# # Convert count matrix to binary format [1=Expressed, 0=Not expressed]
# binary <- normalized_counts %>% 
#   tibble::remove_rownames() %>%
#   tibble::column_to_rownames("SYMBOL") %>% 
#   replace(.> 0, 1) %>%
#   dplyr::mutate(percent = rowSums(.)/ncol(.)*100) %>%
#   dplyr::arrange(desc(percent)) %>%
#   dplyr::filter(percent > 20) %>%
#   tibble::rownames_to_column("SYMBOL") %>%
#   dplyr::filter(SYMBOL %in% x_genes)
normalized_counts <- normalized_counts[rowSums(normalized_counts[,-1]) != 0,]

plot_genes <- c("Ddx3y","Eif2s3y", "Kdm5d", "Uty", "Rbmy", "Sly", "Sry", "Uba1y",
                "Usp9y", "Zfy1", "Zfy2", "H2al2b", "H2al2c", "Orly", "Rbm31y",
                "Srsy", "Ssty1", "Ssty2", "Ddx3x", "Eif2s3x", "Kdm5c", "Kdm6a",
                "Rbmx", "Slx", "Sox3", "Uba1", "Usp9x", "Zfx")
disp_genes <- plot_genes

plot_heatmap(normalized_counts, metadata, plot_genes, disp_genes)


#******************************************************************************#
#                    Plot UMAP of TMB and LOY signature                        #
#******************************************************************************#

celltype <- "Epithelial"
integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn",
                                    dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds")))) 
integrated_seurat <- subset(x = integrated_seurat,
                            subset = Sex == "Male")

y_genes <- c("Ddx3y", "Eif2s3y", "Kdm5d","Uty","Rbmy","Sly","Sry","Uba1y",
             "Usp9y","Zfy1","Zfy2","H2al2b","H2al2c","Orly","Rbm31y","Srsy",
             "Ssty1","Ssty2")
features <- intersect(y_genes, rownames(integrated_seurat@assays$RNA$data))
label <- "y_genes"

x <- Seurat::AddModuleScore(object = integrated_seurat,
                            features = features,
                            name = label,
                            pool = NULL,
                            nbin = 24,
                            ctrl = 5, #100,
                            k = FALSE,
                            assay = "RNA",
                            seed = 1,
                            search = FALSE)

Seurat::FeaturePlot(object = x,
                    slot = "data",
                    features = paste0(label, "1"),
                    cols =  c("grey", viridis(n = 10, option = "C", direction = -1)),
                    pt.size = 0.4,
                    order = TRUE,
                    min.cutoff = 'q10',
                    reduction = "umap.harmony",
                    label = TRUE,
                    combine = TRUE,
                    raster = FALSE)


ggplot2::ggsave(filename = paste0("Module_plot(", label, ").pdf"),
                plot = last_plot(),
                device = "pdf",
                path = seurat_results,
                scale = 1,
                width = 8.5*4,
                height = 11*2,
                units = c("in"),
                dpi = 300,
                limitsize = FALSE,
                bg = "white")


#******************************************************************************#
#               CORRELATION OF Y EXPRESSION vs EACH GENE IN SINGLE CELL        #  
#******************************************************************************#

celltype <- "Epithelial"
integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn",
                                    dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))

expr_mat <- integrated_seurat@assays$RNA@data
expr_mat <- expr_mat[rowSums(expr_mat) != 0,]
expr_mat <- t(expr_mat)
expr_mat <- as.data.frame(expr_mat) %>%
  dplyr::mutate(Y = rowSums(select(., intersect(y_ensembl, colnames(.)))))

columnnames <- setdiff(colnames(expr_mat), y_ensembl)
columnnames <- setdiff(columnnames, "Y")

sig_genes <- c()
sig_pval <- c()
sig_cor <- c()
for (i in columnnames){
  
  subset <- expr_mat %>% select(all_of(i), "Y")
  subset <- subset[rowSums(subset) != 0, ]
  cor_results <- cor.test(x=subset[[1]], y=subset[[2]], method=c("pearson"))
  sig_genes <- c(sig_genes, i)
  sig_pval <- c(sig_pval, cor_results$p.value)
  sig_cor <- c(sig_cor, cor_results$estimate[[1]])
  
}

#******************************************************************************#
#                  BOX PLOT Y GENES FROM TRAMP C1 BULK RNA SEQ                #
#******************************************************************************#

read_data <- openxlsx::read.xlsx("C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/TRAMP_GSE79756/Results__Normalized_Counts.xlsx")
colnames(read_data)[1] <- "ID"

# Get Y genes from ensembl
y_ensembl <- annotations%>% 
  dplyr::filter(chr == "Y") %>%
  dplyr::filter(!grepl(pattern = "RIKEN|predicted|pseudogene|rRNA", x = DESCRIPTION)) %>%
  dplyr::filter(nchar(SYMBOL) > 0) %>%
  dplyr::arrange(SYMBOL) %>%
  dplyr::select(SYMBOL) %>% unlist(., use.names = FALSE)

# Get genes from bulk RNA Seq
y_seurat <- read_data$SYMBOL

# Keep y genes common between seurat and ensembl
y_genes <- base::intersect(y_ensembl, y_seurat)

# Keep expr data for only Y genes
expr_data <- read_data %>% dplyr::filter(SYMBOL %in% y_genes) %>%
  dplyr::select(everything(), -ID) %>%
  tibble::column_to_rownames("SYMBOL") %>%
  dplyr::filter(rowSums(.) > 0)

#******************************************************************************#
#                        VENN DIAGRAM ANALYSIS OF DEGS                         #
#******************************************************************************#

results_path <- "C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/"

# Merge count data
data1 <- openxlsx::read.xlsx(paste0(results_path, "Hany_YKO/Results__id_YKO_centromere_vs_scr_KO_centromere_DEGs.xlsx"), rowNames=TRUE)
data2 <- openxlsx::read.xlsx(paste0(results_path, "Hany_Y/Results__id_Y_pos_vs_Y_neg_DEGs.xlsx"), rowNames=TRUE)
data3 <- openxlsx::read.xlsx(paste0(results_path, "Hany_OE_Tumor/Results__id_Kdm5d_OE_Tumor_vs_scr_OE_Tumor_DEGs.xlsx"), rowNames=TRUE)
data4 <- openxlsx::read.xlsx(paste0(results_path, "Hany_OE_Tumor/Results__id_Uty_OE_Tumor_vs_scr_OE_Tumor_DEGs.xlsx"), rowNames=TRUE)
data5 <- openxlsx::read.xlsx(paste0(results_path, "Hany_OE/Results__id_Kdm5d_OE_vs_scr_OE_DEGs.xlsx"), rowNames=TRUE)
data6 <- openxlsx::read.xlsx(paste0(results_path, "Hany_OE/Results__id_Uty_OE_vs_scr_OE_DEGs.xlsx"), rowNames=TRUE)
data7 <- openxlsx::read.xlsx(paste0(results_path, "Hany_KO/Results__id_Kdm5d_KO_vs_scr_KO_DEGs.xlsx"), rowNames=TRUE)
data8 <- openxlsx::read.xlsx(paste0(results_path, "Hany_KO/Results__id_Uty_KO_vs_scr_KO_DEGs.xlsx"), rowNames=TRUE)
data9 <- openxlsx::read.xlsx(paste0(results_path, "Hany_Y_Tumor/Results__id_Y_pos_vs_Y_neg_DEGs.xlsx"), rowNames=TRUE)

# Save to excel
wb <- openxlsx::createWorkbook()

YKO_high <- unique(data1 %>% 
                     dplyr::filter(padj < 0.05 & log2FoldChange >= 0.58) %>% 
                     dplyr::select(SYMBOL) %>% unlist(., use.names = FALSE))
YKO_low <- unique(data1 %>% 
                    dplyr::filter(padj < 0.05 & log2FoldChange <= -0.58) %>% 
                    dplyr::select(SYMBOL) %>% unlist(., use.names = FALSE))
Yneg_low <- unique(data2 %>% 
                     dplyr::filter(padj < 0.05 & log2FoldChange >= 0.58) %>% 
                     dplyr::select(SYMBOL) %>% unlist(., use.names = FALSE))
Yneg_high <- unique(data2 %>% 
                      dplyr::filter(padj < 0.05 & log2FoldChange <= -0.58) %>% 
                      dplyr::select(SYMBOL) %>% unlist(., use.names = FALSE))

l <- max(length(YKO_high), length(YKO_low), length(Yneg_low), length(Yneg_high))
YKO_high <- c(YKO_high, rep(NA, times = l-length(YKO_high)))
YKO_low <- c(YKO_low, rep(NA, times = l-length(YKO_low)))
Yneg_low <- c(Yneg_low, rep(NA, times = l-length(Yneg_low)))
Yneg_high <- c(Yneg_high, rep(NA, times = l-length(Yneg_high)))
df <- data.frame(YKO_high, YKO_low, Yneg_high, Yneg_low)
openxlsx::addWorksheet(wb, sheetName = "YKO")
openxlsx::writeData(wb, sheet = "YKO", x = df, rowNames = FALSE)

Kdm5dOET_high <- unique(data3 %>% 
                          dplyr::filter(padj < 0.05 & log2FoldChange >= 0.58) %>% 
                          dplyr::select(SYMBOL) %>% unlist(., use.names = FALSE))
Kdm5dOET_low <- unique(data3 %>% 
                         dplyr::filter(padj < 0.05 & log2FoldChange <= -0.58) %>% 
                         dplyr::select(SYMBOL) %>% unlist(., use.names = FALSE))
UtyOET_high <- unique(data4 %>% 
                        dplyr::filter(padj < 0.05 & log2FoldChange >= 0.58) %>% 
                        dplyr::select(SYMBOL) %>% unlist(., use.names = FALSE))
UtyOET_low <- unique(data4 %>% 
                       dplyr::filter(padj < 0.05 & log2FoldChange <= -0.58) %>% 
                       dplyr::select(SYMBOL) %>% unlist(., use.names = FALSE))

l <- max(length(Kdm5dOET_high), length(Kdm5dOET_low), length(UtyOET_high), length(UtyOET_low))
Kdm5dOET_high <- c(Kdm5dOET_high, rep(NA, times = l-length(Kdm5dOET_high)))
Kdm5dOET_low <- c(Kdm5dOET_low, rep(NA, times = l-length(Kdm5dOET_low)))
UtyOET_high <- c(UtyOET_high, rep(NA, times = l-length(UtyOET_high)))
UtyOET_low <- c(UtyOET_low, rep(NA, times = l-length(UtyOET_low)))
df <- data.frame(Kdm5dOET_high, Kdm5dOET_low, UtyOET_high, UtyOET_low)
openxlsx::addWorksheet(wb, sheetName = "Tumor")
openxlsx::writeData(wb, sheet = "Tumor", x = df, rowNames = FALSE)


Kdm5dOE_high <- unique(data5 %>% 
                         dplyr::filter(padj < 0.05 & log2FoldChange >= 0.58) %>% 
                         dplyr::select(SYMBOL) %>% unlist(., use.names = FALSE))
Kdm5dOE_low <- unique(data5 %>% 
                        dplyr::filter(padj < 0.05 & log2FoldChange <= -0.58) %>% 
                        dplyr::select(SYMBOL) %>% unlist(., use.names = FALSE))
Kdm5dKO_high <- unique(data7 %>% 
                         dplyr::filter(padj < 0.05 & log2FoldChange >= 0.58) %>% 
                         dplyr::select(SYMBOL) %>% unlist(., use.names = FALSE))
Kdm5dKO_low <- unique(data7 %>% 
                        dplyr::filter(padj < 0.05 & log2FoldChange <= -0.58) %>% 
                        dplyr::select(SYMBOL) %>% unlist(., use.names = FALSE))

l <- max(length(Kdm5dOE_high), length(Kdm5dOE_low), length(Kdm5dKO_high), length(Kdm5dKO_low))
Kdm5dOE_high <- c(Kdm5dOE_high, rep(NA, times = l-length(Kdm5dOE_high)))
Kdm5dOE_low <- c(Kdm5dOE_low, rep(NA, times = l-length(Kdm5dOE_low)))
Kdm5dKO_high <- c(Kdm5dKO_high, rep(NA, times = l-length(Kdm5dKO_high)))
Kdm5dKO_low <- c(Kdm5dKO_low, rep(NA, times = l-length(Kdm5dKO_low)))
df <- data.frame(Kdm5dOE_high, Kdm5dOE_low, Kdm5dKO_low, Kdm5dKO_high)
openxlsx::addWorksheet(wb, sheetName = "Kdm5d")
openxlsx::writeData(wb, sheet = "Kdm5d", x = df, rowNames = FALSE)

UtyOE_high <- unique(data6 %>% 
                       dplyr::filter(padj < 0.05 & log2FoldChange >= 0.58) %>% 
                       dplyr::select(SYMBOL) %>% unlist(., use.names = FALSE))
UtyOE_low <- unique(data6 %>% 
                      dplyr::filter(padj < 0.05 & log2FoldChange <= -0.58) %>% 
                      dplyr::select(SYMBOL) %>% unlist(., use.names = FALSE))
UtyKO_high <- unique(data8 %>% 
                       dplyr::filter(padj < 0.05 & log2FoldChange >= 0.58) %>% 
                       dplyr::select(SYMBOL) %>% unlist(., use.names = FALSE))
UtyKO_low <- unique(data8 %>% 
                      dplyr::filter(padj < 0.05 & log2FoldChange <= -0.58) %>% 
                      dplyr::select(SYMBOL) %>% unlist(., use.names = FALSE))

l <- max(length(UtyOE_high), length(UtyOE_low), length(UtyKO_high), length(UtyKO_low))
UtyOE_high <- c(UtyOE_high, rep(NA, times = l-length(UtyOE_high)))
UtyOE_low <- c(UtyOE_low, rep(NA, times = l-length(UtyOE_low)))
UtyKO_high <- c(UtyKO_high, rep(NA, times = l-length(UtyKO_high)))
UtyKO_low <- c(UtyKO_low, rep(NA, times = l-length(UtyKO_low)))
df <- data.frame(UtyOE_high, UtyOE_low, UtyKO_low, UtyKO_high)
openxlsx::addWorksheet(wb, sheetName = "Uty")
openxlsx::writeData(wb, sheet = "Uty", x = df, rowNames = FALSE)

Yneg_low <- unique(data2 %>% 
                     dplyr::filter(padj < 0.05 & log2FoldChange >= 0.58) %>% 
                     dplyr::select(SYMBOL) %>% unlist(., use.names = FALSE))
Yneg_high <- unique(data2 %>% 
                      dplyr::filter(padj < 0.05 & log2FoldChange <= -0.58) %>% 
                      dplyr::select(SYMBOL) %>% unlist(., use.names = FALSE))
Yneg_T_low <- unique(data9 %>% 
                       dplyr::filter(padj < 0.05 & log2FoldChange >= 0.58) %>% 
                       dplyr::select(SYMBOL) %>% unlist(., use.names = FALSE))
Yneg_T_high <- unique(data9 %>% 
                        dplyr::filter(padj < 0.05 & log2FoldChange <= -0.58) %>% 
                        dplyr::select(SYMBOL) %>% unlist(., use.names = FALSE))

l <- max(length(Yneg_low), length(Yneg_high), length(Yneg_T_low), length(Yneg_T_high))
Yneg_low <- c(Yneg_low, rep(NA, times = l-length(Yneg_low)))
Yneg_high <- c(Yneg_high, rep(NA, times = l-length(Yneg_high)))
Yneg_T_low <- c(Yneg_T_low, rep(NA, times = l-length(Yneg_T_low)))
Yneg_T_high <- c(Yneg_T_high, rep(NA, times = l-length(Yneg_T_high)))
df <- data.frame(Yneg_T_high, Yneg_T_low, Yneg_high, Yneg_low)
openxlsx::addWorksheet(wb, sheetName = "Y_Tumor")
openxlsx::writeData(wb, sheet = "Y_Tumor", x = df, rowNames = FALSE)

openxlsx::saveWorkbook(wb, 
                       file = paste0(results_path, "Venn_input1.xlsx"), 
                       overwrite = TRUE)

#******************************************************************************#
#                 MERGE NORMALIZED & METADATA OF BULK RNA SEQ                  #
#******************************************************************************#

results_path <- "C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/"

# Merge count data
data1 <- openxlsx::read.xlsx(paste0(results_path, "Hany_YKO/Results__Normalized_Counts.xlsx"))
data2 <- openxlsx::read.xlsx(paste0(results_path, "Hany_Y/Results__Normalized_Counts.xlsx"))
data3 <- openxlsx::read.xlsx(paste0(results_path, "Hany_OE_Tumor/Results__Normalized_Counts.xlsx"))
data4 <- openxlsx::read.xlsx(paste0(results_path, "Hany_OE/Results__Normalized_Counts.xlsx"))
data5 <- openxlsx::read.xlsx(paste0(results_path, "Hany_KO/Results__Normalized_Counts.xlsx"))

colnames(data1)[1] <- "ID"
colnames(data2)[1] <- "ID"
colnames(data3)[1] <- "ID"
colnames(data4)[1] <- "ID"
colnames(data5)[1] <- "ID"

data <- left_join(data1,data2 %>% dplyr::select(everything(), -SYMBOL), by=c("ID"="ID"))
data <- left_join(data,data3 %>% dplyr::select(everything(), -SYMBOL), by=c("ID"="ID"))
data <- left_join(data,data4 %>% dplyr::select(everything(), -SYMBOL), by=c("ID"="ID"))
data <- left_join(data,data5 %>% dplyr::select(everything(), -SYMBOL), by=c("ID"="ID"))

data <- data[rowSums(data[,3:ncol(data)]) != 0,]
data <- data %>% 
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>% 
  select(everything(), -ID) %>%
  tibble::column_to_rownames("SYMBOL")

# Merge metadata
data1 <- openxlsx::read.xlsx(paste0(results_path, "Hany_YKO/Metadata.xlsx"))
data2 <- openxlsx::read.xlsx(paste0(results_path, "Hany_Y/Metadata.xlsx"))
data3 <- openxlsx::read.xlsx(paste0(results_path, "Hany_OE_Tumor/Metadata.xlsx"))
data4 <- openxlsx::read.xlsx(paste0(results_path, "Hany_OE/Metadata.xlsx"))
data5 <- openxlsx::read.xlsx(paste0(results_path, "Hany_KO/Metadata.xlsx"))
metadata <- dplyr::bind_rows(data1,data2,data3,data4,data5)
metadata <- metadata %>% filter(Sample %in% intersect(metadata$Sample, colnames(data)))

# Correct for batch effects
batch <- metadata$Batch
mod <- stats::model.matrix(object = ~Condition, data=metadata)
modcombat = stats::model.matrix(object = ~1, data=metadata)
combat_edata <- sva::ComBat(dat = data, batch = batch, mod = modcombat, par.prior=TRUE, prior.plots=FALSE)

# Save results for progeny and dorothea analysis
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Normalized")
openxlsx::writeData(wb, sheet = "Normalized", x = combat_edata, rowNames = TRUE)
openxlsx::saveWorkbook(wb, 
                       file = paste0(results_path, "Normalized.xlsx"), 
                       overwrite = TRUE)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Metadata")
openxlsx::writeData(wb, sheet = "Metadata", x = metadata, rowNames = TRUE)
openxlsx::saveWorkbook(wb, 
                       file = paste0(results_path, "Metadata.xlsx"), 
                       overwrite = TRUE)

#******************************************************************************#
#                          DOROTHEA & PROGENY ANLAYSIS                         #
#******************************************************************************#

meta_data <- openxlsx::read.xlsx(xlsxFile = paste0(results_path, "Hany_bulk_RNASeq_Metadata.xlsx"))
normalized_counts <- openxlsx::read.xlsx(xlsxFile = paste0(results_path, "Hany_bulk_RNASeq_Normalized.xlsx"), rowNames=TRUE)

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
mat <- as.matrix(progeny_scores)

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
                   filename = paste0(results_path, "Progeny.pdf"))

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
mat <- as.matrix(viper_scores)

# Define column annotation
col_annotation <- meta_data %>% dplyr::select(Condition)
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

# Find TFs that are statistically different between 2 groups
# NOTE: We assume normal distribution
# Use F.test() to determine if variances are equal or not.
sig_genes <- c()
sig_pval <- c()
for (j in 1:nrow(mat)){
  
  # DO NOT USE BIND_COLS UNLESS ROW ORDERS IS SAME
  data <- dplyr::left_join(x=as.data.frame(viper_object[j,]) %>% tibble::rownames_to_column(var = "id"),
                           y = col_annotation %>% tibble::rownames_to_column(var = "id"),
                           by=c("id"="id")) %>%
    tibble::column_to_rownames(var = "id")
  colnames(data) <- c("values", "Condition")
  
  # Check if variances are equal (p <0.05 => unequal variance)
  f_test <- var.test(formula = values ~ Condition, 
                     data = data,
                     alternative = "two.sided")
  
  if (f_test$p.value < 0.05){
    t_test <- t.test(formula = values ~ Condition, 
                     data = data,
                     alternative = "two.sided",
                     var.equal= FALSE)
  } else {
    t_test <- t.test(formula = values ~ Condition, 
                     data = data,
                     alternative = "two.sided",
                     var.equal= TRUE)
  }
  if (t_test$p.value < 0.05){
    sig_genes <- c(sig_genes, rownames(mat)[j])
    sig_pval <- c(sig_pval, t_test$p.value)
  }
}
print(data.frame(sig_genes, sig_pval))

# Define how samples will be arranged in the heatmap.
# Set cluster_cols=FALSE, if you want to arrange samples in specific order
# Set cluster_cols=TRUE, if you want to arrange samples based on clustering
mat <- mat[sig_genes, rownames(col_annotation)]

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
                   cellwidth = 10, 
                   #cellheight = 10, 
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
                   annotation_colors = dplyr::if_else(nrow(col_annotation)+nrow(row_annotation) > 0, ann_colors, c()),
                   annotation_legend = TRUE,
                   annotation_names_row = FALSE,
                   annotation_names_col = FALSE,
                   show_rownames = dplyr::if_else(nrow(mat)<100, TRUE, FALSE, missing = NULL), 
                   show_colnames = dplyr::if_else(ncol(mat)<50, TRUE, FALSE, missing = NULL),
                   fontsize = 8, 
                   fontsize_row = 5, 
                   fontsize_col = 8,
                   angle_col = c("270", "0", "45", "90", "315"),
                   fontsize_number = 0.8*fontsize, 
                   labels_row = c(display_row, rep(x="", times=nrow(mat) - length(display_row))),
                   labels_col = c(display_col, rep(x="", times=ncol(mat) - length(display_col))),
                   #width = 8.5, #ncol(mat)*0.075+5,
                   #height = 30,
                   filename = paste0(results_path, "Dorothea1.pdf"))

#******************************************************************************#

integrated_seurat <- subset(x = integrated_seurat,
                            subset = Sex == "Male")

y_genes <- c("Ddx3y", "Eif2s3y", "Kdm5d","Uty","Rbmy","Sly","Sry","Uba1y",
             "Usp9y","Zfy1","Zfy2","H2al2b","H2al2c","Orly","Rbm31y","Srsy",
             "Ssty1","Ssty2")
features <- intersect(y_genes, rownames(integrated_seurat@assays$RNA$data))

i <- 1
celltypes <- unique(integrated_seurat@meta.data$seurat_class)
for (c in celltypes){
  
  x <- paste0("p",i)
  i <- i+1;
  
  seurat_obj <- subset(x=integrated_seurat,
                       seurat_class == c)
  
  # Extract expression inf, keep only Y genes that have expression
  df <- seurat_obj@assays$RNA$data
  df <- df[features,]
  df <- df[rowSums(df) != 0,]
  df[df>0] <- 1
  
  
  gene_count_per_cell <- data.frame(counts = colSums(df)) %>%
    tibble::rownames_to_column("Cell") %>%
    dplyr::mutate(Cell = gsub(pattern= "_.*" , replacement="", x=Cell)) %>%
    dplyr::group_by(Cell, counts) %>%
    dplyr::summarise(n = n()) %>%
    #tidyr::pivot_wider(id_cols=Cell, names_from=counts, values_from=n) %>%
    dplyr::rename(Sample=Cell, n_gene=counts) %>%
    dplyr::mutate(Percent = 100*n/sum(n)) %>%
    dplyr::filter(n_gene==0)
  
  print(c)
  print(gene_count_per_cell)
  
  p <- ggplot(data = gene_count_per_cell, aes(x=n_gene, y=Percent, group=Sample, fill=Sample)) +
    geom_col() +
    my_theme +
    facet_wrap(~Sample, nrow=1) +
    scale_x_continuous(breaks=seq(0,6,1)) + 
    labs(x="Number of Y genes expressed in a single cell", y="% of cells expressing indicated number of Y genes",
         title = c) 
    # theme(legend.position="right",
    #       panel.spacing = unit(0.1, "lines"),
    #       axis.ticks.x=element_blank())
  
  assign(x,p)
}

# Save all plots
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,
          align = c("hv"),
          axis = c("tblr"),
          nrow = 2,
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
          byrow = TRUE)

ggsave(filename = "Y Expression Frequency.pdf",
       plot = last_plot(),
       device = "pdf",
       width = 33,
       height = 11,
       units = c("in"),	 
       dpi = 300,
       limitsize = TRUE,
       bg = "white")

i <- 1
celltypes <- unique(integrated_seurat@meta.data$seurat_class)
for (c in celltypes){
  
  x <- paste0("p",i)
  i <- i+1;
  
  seurat_obj <- subset(x=integrated_seurat,
                       seurat_class == c)
  seurat_obj <- Seurat::AddModuleScore(obj = seurat_obj,
                                       features = list(features),
                                       assay = "RNA",
                                       slot = "data",
                                       name = "Yscore")
  
  gene_count_per_cell <- seurat_obj@meta.data %>%
    dplyr::select(Yscore1) %>%
    tibble::rownames_to_column("Cell") %>%
    dplyr::mutate(Cell = gsub(pattern= "_.*" , replacement="", x=Cell)) %>%
    dplyr::group_by(Cell) %>%
    dplyr::summarise(n = mean(Yscore1)) %>%
    #tidyr::pivot_wider(id_cols=Cell, names_from=counts, values_from=n) %>%
    dplyr::rename(Sample=Cell, Yscore=n)
  
  p <- ggplot(data = gene_count_per_cell, aes(y=Sample, x=Yscore, group=Sample, fill=Sample)) +
    geom_col(width=0.25) +
    my_theme + 
    coord_cartesian(xlim = c(-0.125,0.125), clip = "off") +
    labs(x="Y Score", y="Sample", title = c) + 
    geom_text(aes(label=round(Yscore,2)), position = position_stack(vjust = .5))
  
  assign(x,p)
  
}

# Save all plots
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,
          align = c("hv"),
          axis = c("tblr"),
          nrow = 2,
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
          byrow = TRUE)

ggsave(filename = "Y Score Distribution.pdf",
       plot = last_plot(),
       device = "pdf",
       width = 33,
       height = 11,
       units = c("in"),	 
       dpi = 300,
       limitsize = TRUE,
       bg = "white")
  

######## LOY % in each subtype

proj <- "scRNASeq_Koltsova_sn"
proj <- "scRNASeq_Koltsova"
source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")

integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn.rds"))
integrated_seurat <- subset(x = integrated_seurat,
                            subset = Sex == "Male")

y_genes <- c("Ddx3y", "Eif2s3y", "Kdm5d","Uty","Rbmy","Sly","Sry","Uba1y",
             "Usp9y","Zfy1","Zfy2","H2al2b","H2al2c","Orly","Rbm31y","Srsy",
             "Ssty1","Ssty2")

features <- intersect(y_genes, rownames(integrated_seurat@assays$RNA$data))
celltypes <- unique(integrated_seurat@meta.data$cell_type)
celltypes <- celltypes[!is.na(celltypes)]
final_df <- data.frame(Sample="",n_gene=0,n=0, Percent=0, celltype="")

for (c in celltypes){

  seurat_obj <- subset(x=integrated_seurat,
                       cell_type == c)
  
  # Extract expression inf, keep only Y genes that have expression
  df <- seurat_obj@assays$RNA$data
  df <- df[features,]
  df <- df[rowSums(df) != 0,]
  
  # Convert count matrix to binary format [1=Expressed, 0=Not expressed]
  #df <- as.data.frame(df) %>% replace(.> 0, 1)
  df[df>0] <- 1
  
  gene_count_per_cell <- data.frame(counts = colSums(df)) %>%
    tibble::rownames_to_column("Cell") %>%
    dplyr::mutate(Cell = gsub(pattern= "_.*" , replacement="", x=Cell)) %>%
    dplyr::group_by(Cell, counts) %>%
    dplyr::summarise(n = n()) %>%
    #tidyr::pivot_wider(id_cols=Cell, names_from=counts, values_from=n) %>%
    dplyr::rename(Sample=Cell, n_gene=counts) %>%
    dplyr::mutate(Percent = 100*n/sum(n)) %>%
    dplyr::filter(n_gene == 0) %>%
    dplyr::mutate(celltype = c)
  
  final_df <- dplyr::bind_rows(final_df, gene_count_per_cell)
}

final_df <- final_df[-1,]
final_df <- final_df %>% 
  dplyr::select(Sample, Percent, celltype) %>% 
  tidyr::pivot_wider(id_cols=celltype, names_from=Sample, values_from=Percent, values_fill = NA)

# Save the clustered matrix
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "LOY %")
openxlsx::writeData(wb, sheet = "LOY %", x = final_df, rowNames = FALSE)
openxlsx::saveWorkbook(wb, file = "LOY Percent1.xlsx", overwrite = TRUE)

library(openxlsx)

annotations <- get_annotations(species)

x_GOI <- annotations %>% 
  dplyr::filter(CHR=="X" ) %>%
  dplyr::mutate(SYMBOL = dplyr::case_when(is.na(SYMBOL) ~ ENSEMBL_ID, !is.na(SYMBOL) ~ SYMBOL)) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  dplyr::mutate(SYMBOL = make.names(SYMBOL))

y_GOI <- annotations %>% 
  dplyr::filter(CHR=="Y") %>%
  dplyr::mutate(SYMBOL = dplyr::case_when(is.na(SYMBOL) ~ ENSEMBL_ID, !is.na(SYMBOL) ~ SYMBOL)) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  dplyr::mutate(SYMBOL = make.names(SYMBOL))

GOI <- annotations %>% 
  #dplyr::filter(CHR=="X" | CHR == "Y") %>%
  dplyr::mutate(SYMBOL = dplyr::case_when(is.na(SYMBOL) ~ ENSEMBL_ID, !is.na(SYMBOL) ~ SYMBOL)) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  dplyr::mutate(SYMBOL = make.names(SYMBOL))

#TCGA
tcga <- read.xlsx("/hpc/home/kailasamms/scratch/biomarker/TCGA_Normalized_Counts.xlsx") %>%
  dplyr::mutate(SYMBOL = dplyr::case_when(is.na(SYMBOL) ~ ENSEMBL_ID, !is.na(SYMBOL) ~ SYMBOL)) %>%
  dplyr::select(everything(), -c("ENSEMBL_ID", "ENTREZ_ID", "SYMBOL_ENTREZ")) %>%
  # If there are duplicated genes, keep only data for highest expressing copy
  dplyr::mutate(n = rowSums(.[,-1])) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::slice_max(n) %>%
  dplyr::ungroup() %>%
  # Duplicated genes with 0 expression in all samples still remain, remove them
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  dplyr::select(everything(), -n) %>%
  dplyr::mutate(SYMBOL = make.names(SYMBOL)) %>%
  dplyr::filter(SYMBOL %in% GOI$SYMBOL) %>%
  data.frame()

tcga_metadata <- read.xlsx("/hpc/home/kailasamms/scratch/biomarker/TCGA_Metadata.xlsx")

# Nat Comms
celltype <- "Epithelial"
seurat_results <- "/hpc/home/kailasamms/scratch/scRNASeq_Simon_old/results_seurat/"
integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn",
                                    dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))

male_integrated_seurat <- subset(x=integrated_seurat,
                                 Sex == "Male")
female_integrated_seurat <- subset(x=integrated_seurat,
                                   Sex == "Female")

male_nat_comm <- male_integrated_seurat@assays$RNA@data %>%
  data.frame() %>%
  tibble::rownames_to_column("SYMBOL") %>%
  dplyr::mutate(SYMBOL = make.names(SYMBOL)) %>%
  dplyr::filter(SYMBOL %in% GOI$SYMBOL)

female_nat_comm <- female_integrated_seurat@assays$RNA@data %>%
  data.frame() %>%
  tibble::rownames_to_column("SYMBOL") %>%
  dplyr::mutate(SYMBOL = make.names(SYMBOL)) %>%
  dplyr::filter(SYMBOL %in% GOI$SYMBOL)

male_nat_comm_metadata <- male_integrated_seurat@meta.data %>%
  dplyr::rename(Sample1 = Sample) %>%
  tibble::rownames_to_column("Sample")

female_nat_comm_metadata <- female_integrated_seurat@meta.data %>%
  dplyr::rename(Sample1 = Sample) %>%
  tibble::rownames_to_column("Sample")

# DepMap
model <- read.csv("/hpc/home/kailasamms/scratch/biomarker/Model.csv")
model <- model %>% dplyr::filter(DepmapModelType=="BLCA")

profile <- read.csv("/hpc/home/kailasamms/scratch/biomarker/OmicsProfiles.csv")
profile <- dplyr::inner_join(profile, model, by=c("ModelID"="ModelID")) %>%
  dplyr::select(StrippedCellLineName, ModelID, ProfileID, Sex, PrimaryOrMetastasis)

ccle <- read.csv("/hpc/home/kailasamms/scratch/biomarker/OmicsExpressionAllGenesTPMLogp1Profile.csv")
colnames(ccle) <- gsub("..ENSG.*$","",colnames(ccle))
colnames(ccle)[1] <- "ProfileID"
ccle <- ccle %>% dplyr::select(ProfileID, intersect(GOI$SYMBOL, make.names(colnames(ccle))))
ccle <- dplyr::inner_join(profile, ccle, by=c("ProfileID"="ProfileID")) 

ccle_metadata <- ccle[,1:5] %>%
  dplyr::rename(Sample = StrippedCellLineName)

ccle <- ccle[,-c(2,3,4,5)] %>% 
  tibble::column_to_rownames("StrippedCellLineName") %>%
  t() %>%
  data.frame() %>%
  tibble::rownames_to_column("SYMBOL") %>%
  dplyr::mutate(SYMBOL = make.names(SYMBOL))

# Find good Y genes
for (file_suffix in c("tcga", "ccle", "nat_comm")){
  
  results_path <- "/hpc/home/kailasamms/"
  # Save the clustered scores in xlsx
  wb <- openxlsx::createWorkbook()
  
  if(file_suffix != "nat_comm"){
    
    metadata_column <- get(paste0(file_suffix, "_metadata"))
    normalized_counts <- get(file_suffix)
    colnames(normalized_counts)[1] <- "SYMBOL"
    colnames(normalized_counts) <- make.names(colnames(normalized_counts))
    
    for (sex in c("Male", "Female")){
      
      samples <- metadata_column %>%
        dplyr::filter(Sex == sex)
      
      avg_exp <- rowMeans(as.matrix(normalized_counts[,make.names(samples$Sample)])) %>%
        data.frame()
      colnames(avg_exp)[1] <- "avg_exp"
      rownames(avg_exp) <- normalized_counts$SYMBOL
      avg_exp <- avg_exp %>% 
        tibble::rownames_to_column("SYMBOL")
      
      # median_exp <- rowMedians(x = as.matrix(normalized_counts[,make.names(samples$Sample)])) %>%
      #   data.frame() 
      # colnames(median_exp)[1] <- "median_exp"
      # rownames(median_exp) <- normalized_counts$SYMBOL
      # median_exp <- median_exp %>% 
      #   tibble::rownames_to_column("SYMBOL")
      
      quantile_exp <- rowQuantiles(x = as.matrix(normalized_counts[,make.names(samples$Sample)])) %>%
        data.frame()
      colnames(quantile_exp) <- c("Q0", "Q25", "Q50", "Q75", "Q100")
      rownames(quantile_exp) <- normalized_counts$SYMBOL
      quantile_exp <- quantile_exp %>% 
        tibble::rownames_to_column("SYMBOL")
      
      percent_exp <- normalized_counts %>% 
        tibble::column_to_rownames("SYMBOL") %>% 
        t() %>% 
        scale() %>%
        t() %>%
        data.frame() %>%
        tibble::rownames_to_column("SYMBOL") %>%
        dplyr::mutate_if(is.numeric, ~ dplyr::if_else(. > 0, 1, 0)) %>%
        dplyr::select(SYMBOL, make.names(samples$Sample)) %>%
        dplyr::mutate(percent_exp = rowMeans(select_if(., is.numeric),na.rm=TRUE)*100) %>%
        dplyr::select(SYMBOL, percent_exp)
      
      avg <- normalized_counts %>% 
        dplyr::select(SYMBOL, make.names(samples$Sample)) %>%
        dplyr::left_join(GOI %>% dplyr::select(SYMBOL, CHR, BIOTYPE), by=c("SYMBOL"="SYMBOL")) %>%
        dplyr::left_join(avg_exp,   by=c("SYMBOL"="SYMBOL")) %>%
        dplyr::left_join(percent_exp, by=c("SYMBOL"="SYMBOL")) %>%
        dplyr::left_join(quantile_exp,   by=c("SYMBOL"="SYMBOL")) %>%
        #dplyr::left_join(median_exp,  by=c("SYMBOL"="SYMBOL")) %>%
        dplyr::select(SYMBOL, CHR, BIOTYPE, avg_exp, percent_exp, Q0, Q25, Q50, Q75, Q100) #,everything())
      
      openxlsx::addWorksheet(wb, sheetName = sex)
      openxlsx::writeData(wb, sheet = sex, x = avg, rowNames = FALSE)
    }
    openxlsx::saveWorkbook(wb, file = paste0(results_path, file_suffix, "_Summary.xlsx"), 
                           overwrite = TRUE)
  }
  
  if (file_suffix == "nat_comm"){ 
    
    for (sex in c("Male", "Female")){
      
      normalized_counts <- get(paste0(tolower(sex), "_", file_suffix))
      
      avg_exp <- rowMeans(x = as.matrix(normalized_counts[,-1])) %>%
        data.frame()
      colnames(avg_exp)[1] <- "avg_exp"
      rownames(avg_exp) <- normalized_counts$SYMBOL
      avg_exp <- avg_exp %>% 
        tibble::rownames_to_column("SYMBOL")
      
      # median_exp <- rowMedians(x = as.matrix(normalized_counts[,-1])) %>%
      #   data.frame() 
      # colnames(median_exp)[1] <- "median_exp"
      # rownames(median_exp) <- normalized_counts$SYMBOL
      # median_exp <- median_exp %>% 
      #   tibble::rownames_to_column("SYMBOL")
      
      quantile_exp <- rowQuantiles(x = as.matrix(normalized_counts[,-1])) %>%
        data.frame()
      colnames(quantile_exp) <- c("Q0", "Q25", "Q50", "Q75", "Q100")
      rownames(quantile_exp) <- normalized_counts$SYMBOL
      quantile_exp <- quantile_exp %>% 
        tibble::rownames_to_column("SYMBOL")
      
      percent_exp <- normalized_counts %>% 
        tibble::column_to_rownames("SYMBOL") %>% 
        t() %>% 
        scale() %>%
        t()
      
      percent_exp[percent_exp > 0] <- 1
      percent_exp[percent_exp <= 0] <- 0
      
      percent_exp <- percent_exp %>%
        data.frame() %>%
        dplyr::mutate(percent_exp = rowMeans(., na.rm=TRUE)*100) %>%
        tibble::rownames_to_column("SYMBOL") %>%
        dplyr::select(SYMBOL, percent_exp)
      
      avg <- normalized_counts %>% 
        dplyr::left_join(GOI %>% dplyr::select(SYMBOL, CHR, BIOTYPE), by=c("SYMBOL"="SYMBOL")) %>%
        dplyr::left_join(avg_exp,   by=c("SYMBOL"="SYMBOL")) %>%
        dplyr::left_join(percent_exp, by=c("SYMBOL"="SYMBOL")) %>%
        dplyr::left_join(quantile_exp,   by=c("SYMBOL"="SYMBOL")) %>%
        #dplyr::left_join(median_exp,  by=c("SYMBOL"="SYMBOL")) %>%
        dplyr::select(SYMBOL, CHR, BIOTYPE, avg_exp, percent_exp, Q0, Q25, Q50, Q75, Q100)
      
      openxlsx::addWorksheet(wb, sheetName = sex)
      openxlsx::writeData(wb, sheet = sex, x = avg, rowNames = FALSE)
    }
    openxlsx::saveWorkbook(wb, file = paste0(results_path, file_suffix, "_Summary.xlsx"), 
                           overwrite = TRUE)
  }
}

# Find good reference gene
m <- read.xlsx("/hpc/home/kailasamms/scratch/biomarker/Mutated,CNV_Genes.xlsx", sheet="Mutated_Genes")
c <- read.xlsx("/hpc/home/kailasamms/scratch/biomarker/Mutated,CNV_Genes.xlsx", sheet="CNA_Genes")

goi <- setdiff(annotations$SYMBOL, c$Gene)
goi <- setdiff(goi, m$Gene)
goi <- as.data.frame(goi)
colnames(goi) <- "SYMBOL"

GOI <- annotations %>% 
  dplyr::filter(SYMBOL %in% goi$SYMBOL) %>%
  dplyr::mutate(SYMBOL = dplyr::case_when(is.na(SYMBOL) ~ ENSEMBL_ID, !is.na(SYMBOL) ~ SYMBOL)) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  dplyr::mutate(SYMBOL = make.names(SYMBOL))

# Run lines 23 through 94





wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Unmutated")
openxlsx::writeData(wb, sheet = "Unmutated", x = goi, rowNames = FALSE)
openxlsx::saveWorkbook(wb, file = "C:/Users/KailasammS/Desktop/Unmutated_Genes.xlsx", 
                       overwrite = TRUE)

#..good reference genes
gse <- "TCGA_BLCA"
parent_path <- "C:/Users/KailasammS/Box/Saravana@cedars/10. Ongoing Projects/BBN project/BLCA_Cohorts_correct/"

# Import read_data
read_data <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, gse, "_Normalized.xlsx"))
colnames(read_data)[1] <- "SYMBOL"
read_data <- read_data %>% 
  dplyr::mutate(SYMBOL = make.names(SYMBOL, unique=TRUE)) %>%
  tibble::column_to_rownames("SYMBOL")

mean <- rowMeans(read_data)
sd <- apply(X=read_data,MARGIN=1,FUN=sd)
df <- data.frame(mean,sd)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Unmutated")
openxlsx::writeData(wb, sheet = "Unmutated", x = df, rowNames = TRUE)
openxlsx::saveWorkbook(wb, file = "C:/Users/KailasammS/Desktop/Unmutated_Genes.xlsx", 
                       overwrite = TRUE)

ccle <- ccle %>%
  tibble::column_to_rownames("SYMBOL")
mean <- rowMeans(ccle)
sd <- apply(X=ccle,MARGIN=1,FUN=sd)
df <- data.frame(mean,sd)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Unmutated")
openxlsx::writeData(wb, sheet = "Unmutated", x = df, rowNames = TRUE)
openxlsx::saveWorkbook(wb, file = "/hpc/home/kailasamms/Unmutated_Genes.xlsx", 
                       overwrite = TRUE)


# # Remove genes with zero expression in all samples
# read_data <- read_data[rowSums(read_data) != 0,]
#
# # Remove genes not expressed in all samples
# binary1 <- read_data
# binary1[binary1 > 0] <- 1
# sum(rowSums(binary1) >= 0.9*ncol(binary1))
# g1 <- rownames(binary1[rowSums(binary1) >= 0.9*ncol(binary1),])
# 
# # Remove low expressed genes
# binary2 <- read_data
# q <- quantile(as.vector(as.matrix(binary2)))
# binary2[binary2 < 1000] <- 0
# binary2[binary2 > 1000] <- 1
# sum(rowSums(binary2) >= 0.9*ncol(binary2))
# g2 <- rownames(binary2[rowSums(binary2) >= 0.9*ncol(binary2),])
# 
# g <- intersect(g1,g2)
# t <- read_data[g,]
# 
# 
# 
# t <- apply(X=t,MARGIN=1,FUN=sd)
# sort(t)[1:100]
# 
# wb <- openxlsx::createWorkbook()
# openxlsx::addWorksheet(wb, sheetName = "Unmutated")
# openxlsx::writeData(wb, sheet = "Unmutated", x = as.data.frame(t), rowNames = TRUE)
# openxlsx::saveWorkbook(wb, file = "C:/Users/KailasammS/Desktop/Unmutated_Genes.xlsx", 
#                        overwrite = TRUE)
# 



