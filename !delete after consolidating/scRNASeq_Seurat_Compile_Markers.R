# This script is used to compare cluster markers from multiple NGS datasets
# from time to time to generate a master list of markers that can be used to 
# annotate single cell datasets

# We will use the following approach to generate our own markers:

# (1) FindAllMarkers() using harmony at resolution of 0.4
# NOTE: Even a low resolution, you will often find that epithelial cluster has 
# multiple subclusters. We can manually group these subclusters of epithelial
# cells as well as other celltypes and re-run FindAllMarkers() on the new 
# clusters. However, this takes a lot of time when you keep adding datasets over
# time and it becomes troublesome to manually group subclusters into megaclusters.

# (2) We will next identify similar clusters across datasets
# Keep ONLY markers that are 
# (i) highly expressed log2FC >= 0.58 and 
# (ii) padj <=0.05 and 
# (iii) expressed in most of cells in cluster pct.1 >=0.5
# We are not using any cutoffs on pct.2 because like pointed out earlier
# there are multiple epithelial subclusters. So, epithelial markers will be 
# expressed in most of cells in other epithelial subclusters and pct.2 will be
# high. Setting cutoffs like pct.2 < 0.2 will lead to loss of good markers.

# NOTE: If you use all markers you get from output of Seurat::FindAllMarkers() 
# directly to identify similar clusters across datasets, results are completely 
# wrong. YOU MUST FILTER OUT POOR/WEAK markers and use ONLY top 100 markers for 
# each cluster to identify similar clusters. I have observed that the more 
# markers you use, the poorer the ability to identify similar clusters 
# across datasets.

# NOTE: I tried an alternative way to identify similar clusters across datasets.
# I got normalized counts for cells from each cluster in every dataset
# and calculated the mean normalized counts for each gene for every cluster. 
# Then, I removed the genes with zero counts in all samples and scaled the values
# such that each column i.e. cluster has mean 0 and std 1. Then, I calculated
# the euclidean distance between the clusters and identified similar clusters.
# The results were completely wrong. Hence, it is better to find identical 
# clusters across datasets using Seurat::FindAllMarkers(). 

#******************************************************************************#
#                                STEP 1: FIND MARKERS                          #
#******************************************************************************#

# This plots UMAP for each dataset at harmony resolution 0.4 and outputs an 
# excel file with markers for all datasets
find_markers_standard <- function(projects){
  
  # Create dataframe to store output of FindMarkers()
  df_standard <- data.frame(cluster = "", 
                          gene = "",
                          avg_log2FC = 0,
                          p_val_adj = 0,
                          pct.1 = 0.0, 
                          pct.2 = 0.0,
                          ratio = 0.0)
  
  # Iterate through all datasets and find markers
  for (proj in projects){
    
    seurat_results <- paste0("/hpc/home/kailasamms/scratch/", proj, "/results_seurat/")
    res <- 0.4
    reduc <- "harmony"
    celltype <- NULL
    
    # Load the integrated seurat object generated using "scRNASeq_Seurat_Analysis_Import_Data"
    integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn", 
                                        dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
    
    # Plot UMAP
    Seurat::DimPlot(object=integrated_seurat,
                    reduction="umap.harmony",
                    cols=my_palette,
                    label=TRUE,
                    group.by="cluster.0.4.harmony",
                    split.by=NULL,
                    shape.by=NULL,
                    pt.size=0.2,
                    label.size=5,
                    repel=FALSE,
                    raster=FALSE) +
      ggplot2::labs(fill="CLUSTERS",
                    title = proj,
                    x="UMAP_1",
                    y="UMAP_2") +
      my_theme
    ggsave(paste0("UMAP_0.4_harmony_", proj, ".tiff"))
    
    # Change active.ident to major clusters we annotated
    DefaultAssay(integrated_seurat) <- "RNA"
    idents <- "cluster.0.4.harmony"
    Idents(object = integrated_seurat) <- idents
    
    # Find ALL markers
    all_markers_standard <- Seurat::FindAllMarkers(object=integrated_seurat,
                                                   assay = "RNA",
                                                   features = NULL,
                                                   logfc.threshold = 0.25,
                                                   test.use = "wilcox",
                                                   slot = "data",
                                                   min.pct = 0.1,
                                                   min.diff.pct = 0.1,
                                                   node = NULL,
                                                   verbose = TRUE,
                                                   only.pos = TRUE,
                                                   max.cells.per.ident = Inf,
                                                   random.seed = 1,
                                                   latent.vars = NULL,
                                                   min.cells.feature = 3,
                                                   min.cells.group = 1,
                                                   pseudocount.use = 1,
                                                   mean.fxn = NULL,
                                                   fc.name = NULL,
                                                   base = 2,
                                                   return.thresh = 0.01,
                                                   densify = FALSE)
    
    all_markers_standard <- all_markers_standard %>% 
      dplyr::mutate(pct.1 = dplyr::if_else(pct.1 == 0, 0.001, pct.1),
                    pct.2 = dplyr::if_else(pct.2 == 0, 0.001, pct.2),
                    ratio = pct.1/pct.2,
                    cluster = paste0(proj, "_", cluster)) %>%
      dplyr::relocate(cluster, gene, avg_log2FC, p_val_adj, pct.1, pct.2, ratio)
    
    df_standard <- dplyr::bind_rows(df_standard, all_markers_standard)
  }
  
  # Save results
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "Standard markers")
  openxlsx::writeData(wb, sheet = "Standard markers", x = df_standard, rowNames = FALSE)
  openxlsx::saveWorkbook(wb = wb, file = paste0("Markers_standard_ ", length(projects), "_datasets.xlsx"), overwrite = TRUE)
}

find_markers_standard(projects)

#******************************************************************************#
#                     STEP 2: FILTER OUT BAD MARKERS                           #
#******************************************************************************#

# Filter out weak/non-specific markers & identify top 100 markers for every cluster
df_good <- read.xlsx(paste0("Markers_standard_ ", length(projects), "_datasets.xlsx"),
                     sheet = "Standard markers") %>%
  dplyr::filter(p_val_adj <= 0.05 & 
                  pct.1 >= 0.4 & 
                  avg_log2FC >= 0.58) %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(desc(avg_log2FC)) %>%  # log2FC better than ratio
  dplyr::slice_head(n=50) %>%     
  ungroup()

#******************************************************************************#
#            STEP 3: IDENTIFY SIMILAR CLUSTERS ACROSS DATASETS                 #
#******************************************************************************#

# This will give an xlsx file with 30 groups of similar clusters
identify_similar_clusters <- function(df_markers, k, suffix){
  
  # Create similarity matrix
  # Rows and columns MUST have dataset_cluster_number
  # Values indicate number of common genes between any 2 clusters across datasets
  sim_matrix <- matrix(data=0, 
                       nrow=length(unique(df_markers$cluster)), 
                       ncol=length(unique(df_markers$cluster)))
  colnames(sim_matrix) <- unique(df_markers$cluster)
  rownames(sim_matrix) <- unique(df_markers$cluster)
  
  for (i in 1:ncol(sim_matrix)){
    for (j in 1:nrow(sim_matrix)){
      sim_matrix[i,j] <- length(intersect(tolower(df_markers %>% 
                                                    dplyr::filter(cluster == colnames(sim_matrix)[i]) %>%
                                                    dplyr::select(gene) %>%
                                                    unlist(use.names = FALSE)),
                                          tolower(df_markers %>% 
                                                    dplyr::filter(cluster == rownames(sim_matrix)[j]) %>%
                                                    dplyr::select(gene) %>%
                                                    unlist(use.names = FALSE))))
    }
  }
  
  # Group similar clusters together
  mat <- sim_matrix
  #rowclust <- stats::hclust(d = dist(mat))
  #reordered <- mat[rowclust$order,]
  colclust <- stats::hclust(d = dist(t(mat)))
  reordered <- mat[colclust$order, colclust$order]
  
  # Assuming we have 30 different cell types, set k=30
  groups<-cutree(colclust, k=k)
  mat1 <-cbind(mat,groups)
  
  l <- 0
  for (i in 1:k){
    x1 <- subset(mat1, groups==i)
    l1 <- length(rownames(x1))
    l <- max(l, l1)
  }
  g <- list()
  for (i in 1:k){
    
    x1 <- subset(mat1, groups==i)
    l1 <- length(rownames(x1))
    r1 <- c(rownames(x1), rep(x="", times=l-l1))
    g <- c(g, list(r1))
  }
  names(g) <- paste0("Group", seq(1:k))
  
  # Save the clustered similarity matrix
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "Heatmap_groups")
  openxlsx::writeData(wb, sheet = "Heatmap_groups", x = data.frame(g), rowNames = FALSE)
  openxlsx::addWorksheet(wb, sheetName = "Heatmap_standard")
  openxlsx::writeData(wb, sheet = "Heatmap_standard", x = reordered, rowNames = TRUE)
  openxlsx::saveWorkbook(wb = wb, file = paste0("Markers_Compiled_", suffix, ".xlsx"), overwrite = TRUE)
}

df <- df_good
n_celltypes <- 30
file_suffix <- length(projects)
identify_similar_clusters(df, n_celltypes, file_suffix)

####OPTIONAL SHORT CUT METHOD (UNDER TESTING)

df_final <- data.frame(celltype = "",
                       gene = "",
                       avg_log2FC = 0,
                       p_val_adj = 0,
                       pct.1 = 0, 
                       pct.2 = 0,
                       ratio = 0,
                       p_val = 0, 
                       Freq = 0, 
                       Total = 0,
                       Percent = 0)

class <-  read.xlsx(paste0("Markers_Compiled_", length(projects), ".xlsx"),
          sheet = "Heatmap_groups")

for (i in 1:ncol(class)){
  
  clusters <- class[,i] %>% unique()
  df <- df_good %>% 
    dplyr::filter(cluster %in% clusters) %>%
    dplyr::mutate(cluster = gsub(pattern="scRNASeq_|BBN_", replacement="", x=cluster)) %>%
    tidyr::separate_wider_delim(cols=cluster, delim="_", names=c("dataset","celltype")) %>%
    dplyr::mutate(gene = tolower(gene), celltype = colnames(class)[i]) %>%  #make human and mouse genes identical before adding counts
    dplyr::add_count(celltype,gene) %>%
    dplyr::rename(Freq=n)
  
  df <- df %>%
    dplyr::mutate(Total = length(unique(df$dataset)))
  
  df <- df %>%
    dplyr::mutate(Percent = 100*Freq/Total) %>%
    dplyr::filter(Freq > 1, pct.2 < 0.1) %>%  # reduce non-specific markers
    #dplyr::filter(Freq > 1) %>%  # reduce non-specific markers
    dplyr::select(everything(), -dataset) %>%
    dplyr::group_by(celltype) %>%
    dplyr::distinct_at("gene", .keep_all = TRUE) %>%
    dplyr::ungroup()
  
  df_final <- dplyr::bind_rows(df_final, df)
}
df_final <- df_final[-1,]
# df_final %>% 
#   dplyr::filter(Percent >= 75) %>% 
#   dplyr::select(everything(), -p_val) %>%
#   dplyr::count(celltype) %>%
#   dplyr::arrange(celltype)

marker_list <- list()
for (c in unique(df_final$celltype)){
  
  g <- df_final %>% 
    dplyr::filter(celltype == c) %>%
    dplyr::select(gene) %>%
    unlist(use.names=FALSE) %>%
    list()
  names(g) <- c
  
  marker_list <- c(marker_list, g)
}

max_l <- max(lengths(marker_list))
marker_list <- lapply(X=marker_list, FUN=function(x) { c(x, rep(x="", times=max_l-length(x)))})
df_final_marker <- data.frame(marker_list)

# Save the markers
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Markers")
openxlsx::writeData(wb, sheet = "Markers", x = df_final_marker, rowNames = FALSE)
openxlsx::saveWorkbook(wb = wb, file = "Evaluated_Markers.xlsx", overwrite = TRUE)







#******************************************************************************#
#   STEP 4: MERGE CLUSTERS INTO MACROCLUSTERS & FIND MARKER GENES AGAIN        #
#******************************************************************************#

# "Compiled_Markers.xlsx" has identified similar clusters automatically.
# However, you can notice from the heatmap matrix contained in the xlsx file that
# epithelial and fibroblasts are split across multiple groups. 
# Using the UMAP, heatmap matrix and FindAllMarkers() output, identify
# celltypes and merge the clusters into macroclusters. 

# Once this is complete we will re-run FindAllMarkers using the macro-clusters
# we have identified

# Label the macroclusters with proper celltype names in a new sheet 
# "Heatmap_groups_merged". Delete clusters that you are unable to identify from
# the excel sheet. We will later remove these unclassified clusters before 
# runnign FindAllMarkers() again.

# This outputs excel file with markers of macroclusters from all datasets
find_markers_custom <- function(class){
  
  # Create dataframe to store output of FindMarkers()
  df_custom <- data.frame(cluster = "", 
                          gene = "",
                          avg_log2FC = 0,
                          p_val_adj = 0,
                          pct.1 = 0.0, 
                          pct.2 = 0.0,
                          ratio = 0.0)
  
  # Iterate through all datasets and find markers
  for (proj in projects){
    
    seurat_results <- paste0("/hpc/home/kailasamms/scratch/", proj, "/results_seurat/")
    res <- 0.4
    reduc <- "harmony"
    celltype <- NULL
    
    # Load the integrated seurat object generated using "scRNASeq_Seurat_Analysis_Import_Data"
    integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn", 
                                        dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
    
    integrated_seurat@meta.data <- integrated_seurat@meta.data %>%
      dplyr::mutate(New_class = cluster.0.4.harmony)
    
    dataset <- gsub("scRNASeq_", "", proj)
    for (col in 1:ncol(class)){
      
      clusters <- class[,col][grepl(dataset, class[,col])]
      clusters <- as.numeric(gsub(paste0(dataset,"_"), "", clusters))
      cat(colnames(class)[col], ":", clusters, "\n")
      
      integrated_seurat@meta.data <- integrated_seurat@meta.data %>%
        dplyr::mutate(New_class = dplyr::case_when(New_class %in% clusters ~ colnames(class)[col],
                                                   TRUE ~ New_class))
    }
    
    # Remove unclassified cells before finding markers
    integrated_seurat <- subset(integrated_seurat, 
                                New_class %in% colnames(class))
    
    # Change active.ident to major clusters we annotated
    DefaultAssay(integrated_seurat) <- "RNA"
    idents <- "New_class"
    Idents(object = integrated_seurat) <- idents
    
    # Find ALL markers
    all_markers_standard <- Seurat::FindAllMarkers(object=integrated_seurat,
                                                   assay = "RNA",
                                                   features = NULL,
                                                   logfc.threshold = 0.25,
                                                   test.use = "wilcox",
                                                   slot = "data",
                                                   min.pct = 0.1,
                                                   min.diff.pct = 0.1,
                                                   node = NULL,
                                                   verbose = TRUE,
                                                   only.pos = TRUE,
                                                   max.cells.per.ident = Inf,
                                                   random.seed = 1,
                                                   latent.vars = NULL,
                                                   min.cells.feature = 3,
                                                   min.cells.group = 1,
                                                   pseudocount.use = 1,
                                                   mean.fxn = NULL,
                                                   fc.name = NULL,
                                                   base = 2,
                                                   return.thresh = 0.01,
                                                   densify = FALSE)
    
    all_markers_standard <- all_markers_standard %>% 
      dplyr::mutate(pct.1 = dplyr::if_else(pct.1 == 0, 0.001, pct.1),
                    pct.2 = dplyr::if_else(pct.2 == 0, 0.001, pct.2),
                    ratio = pct.1/pct.2,
                    cluster = paste0(proj, "_", cluster)) %>%
      dplyr::relocate(cluster, gene, avg_log2FC, p_val_adj, pct.1, pct.2, ratio)
    
    df_custom <- dplyr::bind_rows(df_custom, all_markers_standard)
  }
  
  # Save results
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "Custom markers")
  openxlsx::writeData(wb, sheet = "Custom markers", x = df_custom, rowNames = FALSE)
  openxlsx::saveWorkbook(wb = wb, file = paste0("Markers_merged_ ", length(projects), "_datasets.xlsx"), overwrite = TRUE)
}

# "scRNASeq_Simon" is not a good dataset. So, we ignore it.
projects <- c("scRNASeq_BBN_C57BL6", 
              "scRNASeq_BBN_Rag", 
              "scRNASeq_Chen", 
              "scRNASeq_GSE164557",
              "scRNASeq_GSE217093", 
              "scRNASeq_GSE222315",
              "scRNASeq_HRA003620",
              "scRNASeq_Jinfen")

# Read the xlsx file with macroclusters and proper celltype names. 
class <- read.xlsx(paste0("Markers_Compiled_", length(projects), ".xlsx"), 
                   sheet = "Heatmap_groups_merged")
find_markers_custom(class)

#******************************************************************************#
#                     STEP 5: FILTER OUT BAD MARKERS                           #
#******************************************************************************#

# Filter out weak/non-specific markers & identify top 50 markers for every cluster
df_good <- read.xlsx(paste0("Markers_merged_ ", length(projects), "_datasets.xlsx"),
                     sheet = "Custom markers") %>%
  dplyr::filter(p_val_adj <= 0.05 & 
                  pct.1 > pct.2 &  
                  avg_log2FC >= 0.58) %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(desc(avg_log2FC)) %>%  # log2FC better than ratio
  dplyr::slice_head(n=50) %>%     
  ungroup()

#******************************************************************************#
#            STEP 6: IDENTIFY SIMILAR CLUSTERS ACROSS DATASETS                 #
#******************************************************************************#

df_good <- df_good %>%
  dplyr::mutate(cluster = gsub(pattern="scRNASeq_|BBN_", replacement="", x=cluster)) %>%
  tidyr::separate_wider_delim(cols=cluster, delim="_", names=c("dataset","celltype")) %>%
  dplyr::mutate(gene = tolower(gene)) %>%  #make human and mouse genes identical before adding counts
  dplyr::add_count(celltype,gene) %>%
  dplyr::rename(Freq=n)
dim(df_good)

# Calculate total datasets for each celltype. 
# Eg: Myofibroblast cells were defined from 3 datasets while Mast cells were
# defined from 7 datasets
total <- c()
for (i in 1:nrow(df_good)){
  total <- c(total, length(unique(df_good %>% 
                                    dplyr::filter(celltype == df_good$celltype[i]) %>% 
                                    dplyr::select(dataset) %>% 
                                    unlist(use.names=FALSE))))
}

df_good$Total <- total
dim(df_good)

# Add percent and remove duplicate markers within each celltype
df_good <- df_good %>%
  dplyr::mutate(Percent = 100*Freq/Total) %>%
  dplyr::filter(Freq > 1, pct.2 < 0.1) %>%  # reduce non-specific markers
  #dplyr::filter(Freq > 1) %>%  # reduce non-specific markers
  dplyr::select(everything(), -dataset) %>%
  dplyr::group_by(celltype) %>%
  dplyr::distinct_at("gene", .keep_all = TRUE) %>%
  dplyr::ungroup()
dim(df_good)

marker_list <- list()
for (c in unique(df_good$celltype)){
  
  g <- df_good %>% 
    dplyr::filter(celltype == c) %>%
    dplyr::select(gene) %>%
    unlist(use.names=FALSE) %>%
    list()
  names(g) <- c
  
  marker_list <- c(marker_list, g)
}

max_l <- max(lengths(marker_list))
marker_list <- lapply(X=marker_list, FUN=function(x) { c(x, rep(x="", times=max_l-length(x)))})
df_final <- data.frame(marker_list)

# Save the clustered similarity matrix
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Markers")
openxlsx::writeData(wb, sheet = "Markers", x = df_final, rowNames = FALSE)
openxlsx::saveWorkbook(wb = wb, file = "Evaluated_Markers.xlsx", overwrite = TRUE)














































# 
# # If a gene is detected as marker for 2 cell types, use it as marker for the 
# # celltype where it has highest percent
# marker_list <- list()
# percent_cols <- colnames(final_df)[grepl("Percent", colnames(final_df))] 
# freq_cols <- colnames(final_df)[grepl("Freq", colnames(final_df))]
# celltype_cols <- setdiff(colnames(final_df), c("Gene", percent_cols, freq_cols))
# 
# for (i in 1:length(celltype_cols)){
#   
#   f <- final_df %>% 
#     dplyr::filter(get(percent_cols[i]) > 50) %>%
#     group_by(Gene) %>%
#     dplyr::mutate(max_percent = max(!!!rlang::syms(percent_cols), na.rm=TRUE),
#                   max_count = max(!!!rlang::syms(freq_cols), na.rm=TRUE)) %>%
#     dplyr::filter(get(percent_cols[i]) == max_percent) %>%
#     dplyr::select(Gene, all_of(celltype), all_of(paste0(celltype, "_count")), max_percent) %>%
#     dplyr::arrange(desc(max_percent)) %>%
#     data.frame() %>%
#     #dplyr::slice_head(n=50) %>%
#     dplyr::select(Gene) %>%
#     unlist(use.names=FALSE)
#   
#   final_g <- list(final_g)
#   names(final_g) <- celltype
#   marker_list <- c(marker_list, final_g)
# }
# 
# max_l <- max(lengths(marker_list))
# marker_list <- lapply(X=marker_list, FUN=function(x) { c(x, rep(x="", times=max_l-length(x)))})
# df_final_standard <- data.frame(marker_list)
# colnames(df_final_standard) <- names(marker_list)
# 
# wb <- openxlsx::createWorkbook()
# openxlsx::addWorksheet(wb = wb, sheetName = "Standard_Markers")
# openxlsx::writeData(wb = wb, sheet = "Standard_Markers", x = df_final_standard)
# openxlsx::addWorksheet(wb = wb, sheetName = "Intermediate step")
# openxlsx::writeData(wb = wb, sheet = "Intermediate step", x = final_df)
# openxlsx::saveWorkbook(wb = wb, file = "Evaluated_Markers.xlsx", overwrite = TRUE)


for (list in c("Epithelial_1", "Epithelial_2", "Fibroblast_1", "Fibroblast_2", 
               "Lymphatic.Endothelial", "Endothelial", "Granulocytes", "Plasma",
               "B", "TCells", "Macrophages", "DCs", "Mast", "Myofibroblasts")){
  
  # Create a list to store markers of each cell type
  final_g <- c()
  
  # Iterate through each dataset
  for (i in 1:length(get(list))){
    
    # Get unique markers for a cell type from each dataset
    genes <- df_custom %>% 
      dplyr::filter(cluster %in% paste0(names(get(list)[i]),"_",  get(list)[[i]])) %>%
      dplyr::select(gene) %>%
      unlist(use.names = FALSE) %>%
      tolower() %>%
      unique()
    
    # Merge all unique markers from each dataset
    final_g <- c(final_g, genes)
  }
  # Calculate number of datasets each marker was identified
  n_data <- sum(lengths(get(list)) > 0)  # number of datasets used for each cell type
  final_g <- data.frame(table(final_g))
  final_g <- final_g %>%
    dplyr::mutate(Percent = Freq*100/n_data)
  colnames(final_g) <- c(list, paste0(list, "_count"), paste0(list, "_percent"))
  
  # Remove markers identified in only a single dataset
  final_g <- final_g %>% 
    dplyr::filter(get(paste0(list, "_count")) > 1) %>%
    dplyr::mutate(Gene = get(list))
  
  # Merge all markers for each cell type present in multiple daasets to final dataframe
  final_df <- final_df %>%
    left_join(final_g,by=c("Gene"="Gene"))
}

final_df[is.na(final_df)] <- 0

# If a gene is detected as marker for 2 cell types, use it as marker for the 
# celltype where it has highest count
marker_list <- list()

for (celltype in c("Epithelial_1", "Epithelial_2", "Fibroblast_1", "Fibroblast_2", 
                   "Lymphatic.Endothelial", "Endothelial", "Granulocytes", "Plasma",
                   "B", "TCells", "Macrophages", "DCs", "Mast", "Myofibroblasts")){
  
  final_g <- final_df %>% 
    group_by(Gene) %>%
    dplyr::mutate(max_percent = max(Epithelial_1_percent, Epithelial_2_percent, 
                                    Fibroblast_1_percent, Fibroblast_2_percent, 
                                    Lymphatic.Endothelial_percent, Endothelial_percent,
                                    Granulocytes_percent, Plasma_percent,
                                    B_percent, TCells_percent, Macrophages_percent,
                                    DCs_percent, Mast_percent, Myofibroblasts_percent, na.rm=TRUE),
                  max_count = max(Epithelial_1_count, Epithelial_2_count, 
                                  Fibroblast_1_count, Fibroblast_2_count, 
                                  Lymphatic.Endothelial_count, Endothelial_count,
                                  Granulocytes_count, Plasma_count,
                                  B_count, TCells_count, Macrophages_count,
                                  DCs_count, Mast_count, Myofibroblasts_count, na.rm=TRUE)) %>%
    dplyr::filter(get(paste0(celltype, "_percent")) == max_percent) %>%
    dplyr::filter(get(paste0(celltype, "_percent")) >= 50) %>%
    dplyr::select(Gene, all_of(celltype), all_of(paste0(celltype, "_count")), max_percent) %>%
    dplyr::arrange(desc(max_percent)) %>%
    data.frame() %>%
    #dplyr::slice_head(n=50) %>%
    dplyr::select(Gene) %>%
    unlist(use.names=FALSE)
  
  final_g <- list(final_g)
  names(final_g) <- celltype
  marker_list <- c(marker_list, final_g)
}

max_l <- max(lengths(marker_list))
marker_list <- lapply(X=marker_list, FUN=function(x) { c(x, rep(x="", times=max_l-length(x)))})
df_final_standard <- data.frame(marker_list)
colnames(df_final_standard) <- names(marker_list)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb = wb, sheetName = "Standard_Markers")
openxlsx::writeData(wb = wb, sheet = "Standard_Markers", x = df_final_standard)
openxlsx::addWorksheet(wb = wb, sheetName = "Intermediate step")
openxlsx::writeData(wb = wb, sheet = "Intermediate step", x = final_df)
openxlsx::saveWorkbook(wb = wb, file = "Evaluated_Markers.xlsx", overwrite = TRUE)

#******************************************************************************#
#      STEP 5: VALIDATE THE SPECIFICITY FOR EACH MARKER GENE USING UMAPS       #
#******************************************************************************#

for (proj in c("scRNASeq_BBN_C57B6", "scRNASeq_BBN_Rag", "scRNASeq_Jinfen",
               "scRNASeq_Chen", "scRNASeq_GSE217093", "scRNASeq_GSE222315",
               "scRNASeq_HRA003620")){
  
  
  seurat_results <- paste0("/hpc/home/kailasamms/scratch/", proj, "/results_seurat/")
  res <- 0.4
  reduc <- "harmony"
  celltype <- NULL
  idents <- paste0("cluster.", res, ".", base::tolower(reduc))
  
  # Load the integrated seurat object generated using "scRNASeq_Seurat_Analysis_Import_Data"
  integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn", 
                                      dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  # Set default assay
  DefaultAssay(integrated_seurat) <- "RNA"
  
  # Set identity to an existing column in meta data
  Idents(object=integrated_seurat) <- idents
  
  df_final <- read.xlsx("Evaluated_Markers.xlsx")
  
  for (col in 1:ncol(df_final)){
    
    features <- unlist(df_final[col], use.names=FALSE)
    features <- rownames(integrated_seurat@assays$RNA$data)[tolower(rownames(integrated_seurat@assays$RNA$data)) %in% 
                                                              tolower(features)]
    features <- sort(features)
    
    feature_plot <- function(i){
      split <- NULL
      Seurat::FeaturePlot(object=integrated_seurat,
                          slot="data",
                          features=i,
                          split.by=split,
                          #cols= c("grey", viridis(n=10, option="C", direction=-1)),
                          pt.size=0.4,
                          order=TRUE,
                          min.cutoff='q10',
                          reduction=paste0("umap.", base::tolower(reduc)),
                          label=TRUE,
                          combine=TRUE,
                          raster=FALSE) +
        scale_colour_gradientn(colours=rev(brewer.pal(n=11, name="RdBu"))[5:11])
    }
    
    # Plot UMAPs
    for (n in 1:ceiling((length(features)/50))){
      j <- 50*n-49
      k <- 50*n
      vec <- features[j:k]
      purrr::map(.x=vec[!is.na(vec)], 
                 .f=feature_plot) %>%
        cowplot::plot_grid(plotlist=.,
                           align="hv",
                           axis="tblr",
                           nrow=5,
                           ncol=10,
                           rel_widths=1,
                           rel_heights=1,
                           greedy=TRUE,
                           byrow=TRUE)
      
      ggplot2::ggsave(filename=paste0(proj, "_", colnames(df_final[col]), "_", n, ".jpg"),
                      plot=last_plot(),
                      device="jpg",
                      path="/hpc/home/kailasamms/plots/",
                      width=8.5*5,
                      height=11*2,
                      units=c("in"),
                      dpi=300,
                      limitsize=FALSE,
                      bg="white")
    }
  }
}

#####END#######################
































Epithelial_1 <- list("scRNASeq_BBN_C57B6" = c(2),
                     "scRNASeq_BBN_Rag" = c(4,9,12,15,17,18,28,7,10),
                     "scRNASeq_Chen" = c(4),
                     "scRNASeq_GSE217093" = c(0,1,3,16),
                     "scRNASeq_GSE222315" = c(14,15,16),
                     "scRNASeq_HRA003620" = c(4),
                     "scRNASeq_Jinfen" = c())
Epithelial_2 <- list("scRNASeq_BBN_C57B6" = c(),
                     "scRNASeq_BBN_Rag" = c(),
                     "scRNASeq_Chen" = c(0,1),
                     "scRNASeq_GSE217093" = c(),
                     "scRNASeq_GSE222315" = c(4,11),
                     "scRNASeq_HRA003620" = c(9,11,12,25),
                     "scRNASeq_Jinfen" = c())

Fibroblast_1 <- list("scRNASeq_BBN_C57B6" = c(1,9),
                     "scRNASeq_BBN_Rag" = c(0,25,27),
                     "scRNASeq_Chen" = c(8),
                     "scRNASeq_GSE217093" = c(14,25,27,33),
                     "scRNASeq_GSE222315" = c(7),
                     "scRNASeq_HRA003620" = c(3),
                     "scRNASeq_Jinfen" = c(13))
Fibroblast_2 <- list("scRNASeq_BBN_C57B6" = c(6),
                     "scRNASeq_BBN_Rag" = c(14),
                     "scRNASeq_Chen" = c(),
                     "scRNASeq_GSE217093" = c(29),
                     "scRNASeq_GSE222315" = c(),
                     "scRNASeq_HRA003620" = c(),
                     "scRNASeq_Jinfen" = c())
Lymphatic.Endothelial <- list("scRNASeq_BBN_C57B6" = c(18),
                              "scRNASeq_BBN_Rag" = c(30),
                              "scRNASeq_Chen" = c(),
                              "scRNASeq_GSE217093" = c(),
                              "scRNASeq_GSE222315" = c(),
                              "scRNASeq_HRA003620" = c(19),
                              "scRNASeq_Jinfen" = c())
Endothelial <- list("scRNASeq_BBN_C57B6" = c(7),
                    "scRNASeq_BBN_Rag" = c(5,24),
                    "scRNASeq_Chen" = c(5,6,20),
                    "scRNASeq_GSE217093" = c(8),
                    "scRNASeq_GSE222315" = c(5),
                    "scRNASeq_HRA003620" = c(15),
                    "scRNASeq_Jinfen" = c())
Granulocytes <- list("scRNASeq_BBN_C57B6" = c(),
                     "scRNASeq_BBN_Rag" = c(6),
                     "scRNASeq_Chen" = c(11),
                     "scRNASeq_GSE217093" = c(),
                     "scRNASeq_GSE222315" = c(8),
                     "scRNASeq_HRA003620" = c(13),
                     "scRNASeq_Jinfen" = c(3))
Plasma <- list("scRNASeq_BBN_C57B6" = c(23),
               "scRNASeq_BBN_Rag" = c(),
               "scRNASeq_Chen" = c(18,25,16,23),
               "scRNASeq_GSE217093" = c(),
               "scRNASeq_GSE222315" = c(19),
               "scRNASeq_HRA003620" = c(16),
               "scRNASeq_Jinfen" = c())
B <- list("scRNASeq_BBN_C57B6" = c(8),
          "scRNASeq_BBN_Rag" = c(),
          "scRNASeq_Chen" = c(12),
          "scRNASeq_GSE217093" = c(17),
          "scRNASeq_GSE222315" = c(3),
          "scRNASeq_HRA003620" = c(6),
          "scRNASeq_Jinfen" = c(14))
TCells <- list("scRNASeq_BBN_C57B6" = c(3,16),
               "scRNASeq_BBN_Rag" = c(),
               "scRNASeq_Chen" = c(2,3,17),
               "scRNASeq_GSE217093" = c(7),
               "scRNASeq_GSE222315" = c(0,1,9,13,17),
               "scRNASeq_HRA003620" = c(0,7,10,14,8,20),
               "scRNASeq_Jinfen" = c(1,5,9,4,15))
Macrophages <- list("scRNASeq_BBN_C57B6" = c(5,12),
                    "scRNASeq_BBN_Rag" = c(2,8,32),
                    "scRNASeq_Chen" = c(9),
                    "scRNASeq_GSE217093" = c(9),
                    "scRNASeq_GSE222315" = c(6),
                    "scRNASeq_HRA003620" = c(5),
                    "scRNASeq_Jinfen" = c(2,12))
DCs <- list("scRNASeq_BBN_C57B6" = c(17),
            "scRNASeq_BBN_Rag" = c(),
            "scRNASeq_Chen" = c(),
            "scRNASeq_GSE217093" = c(30),
            "scRNASeq_GSE222315" = c(),
            "scRNASeq_HRA003620" = c(18),
            "scRNASeq_Jinfen" = c(16))
Mast <- list("scRNASeq_BBN_C57B6" = c(0),
             "scRNASeq_BBN_Rag" = c(3),
             "scRNASeq_Chen" = c(15),
             "scRNASeq_GSE217093" = c(6,31),
             "scRNASeq_GSE222315" = c(10),
             "scRNASeq_HRA003620" = c(2),
             "scRNASeq_Jinfen" = c(0,6))
Myofibroblasts <- list("scRNASeq_BBN_C57B6" = c(14),
                       "scRNASeq_BBN_Rag" = c(19),
                       "scRNASeq_Chen" = c(7,27),
                       "scRNASeq_GSE217093" = c(),
                       "scRNASeq_GSE222315" = c(2),
                       "scRNASeq_HRA003620" = c(22),
                       "scRNASeq_Jinfen" = c())

# Create a dummy dataframe to store markers from standard approach
final_df <- data.frame("Gene"=unique(tolower(df_good$gene)))

for (list in c("Epithelial_1", "Epithelial_2", "Fibroblast_1", "Fibroblast_2", 
               "Lymphatic.Endothelial", "Endothelial", "Granulocytes", "Plasma",
               "B", "TCells", "Macrophages", "DCs", "Mast", "Myofibroblasts")){
  
  # Create a list to store markers of each cell type
  final_g <- c()
  
  # Iterate through each dataset
  for (i in 1:length(get(list))){
    
    # Get unique markers for a cell type from each dataset
    genes <- df_good %>% 
      dplyr::filter(cluster %in% paste0(names(get(list)[i]),"_",  get(list)[[i]])) %>%
      dplyr::select(gene) %>%
      unlist(use.names = FALSE) %>%
      tolower() %>%
      unique()
    
    # Merge all unique markers from each dataset
    final_g <- c(final_g, genes)
  }
  # Calculate number of datasets each marker was identified
  n_data <- sum(lengths(get(list)) > 0)  # number of datasets used for each cell type
  final_g <- data.frame(table(final_g))
  final_g <- final_g %>%
    dplyr::mutate(Percent = Freq*100/n_data)
  colnames(final_g) <- c(list, paste0(list, "_count"), paste0(list, "_percent"))
  
  # Remove markers identified in only a single dataset
  final_g <- final_g %>% 
    dplyr::filter(get(paste0(list, "_count")) > 1) %>%
    dplyr::mutate(Gene = get(list))
  
  # Merge all markers for each cell type present in multiple daasets to final dataframe
  final_df <- final_df %>%
    left_join(final_g,by=c("Gene"="Gene"))
}

final_df[is.na(final_df)] <- 0

# If a gene is detected as marker for 2 cell types, use it as marker for the 
# celltype where it has highest count
marker_list <- list()

for (celltype in c("Epithelial_1", "Epithelial_2", "Fibroblast_1", "Fibroblast_2", 
                   "Lymphatic.Endothelial", "Endothelial", "Granulocytes", "Plasma",
                   "B", "TCells", "Macrophages", "DCs", "Mast", "Myofibroblasts")){
  
  final_g <- final_df %>% 
    group_by(Gene) %>%
    dplyr::mutate(max_percent = max(Epithelial_1_percent, Epithelial_2_percent, 
                                    Fibroblast_1_percent, Fibroblast_2_percent, 
                                    Lymphatic.Endothelial_percent, Endothelial_percent,
                                    Granulocytes_percent, Plasma_percent,
                                    B_percent, TCells_percent, Macrophages_percent,
                                    DCs_percent, Mast_percent, Myofibroblasts_percent, na.rm=TRUE),
                  max_count = max(Epithelial_1_count, Epithelial_2_count, 
                                  Fibroblast_1_count, Fibroblast_2_count, 
                                  Lymphatic.Endothelial_count, Endothelial_count,
                                  Granulocytes_count, Plasma_count,
                                  B_count, TCells_count, Macrophages_count,
                                  DCs_count, Mast_count, Myofibroblasts_count, na.rm=TRUE)) %>%
    dplyr::filter(get(paste0(celltype, "_percent")) == max_percent) %>%
    dplyr::filter(get(paste0(celltype, "_percent")) >= 50) %>%
    dplyr::select(Gene, all_of(celltype), all_of(paste0(celltype, "_count")), max_percent) %>%
    dplyr::arrange(desc(max_percent)) %>%
    data.frame() %>%
    #dplyr::slice_head(n=50) %>%
    dplyr::select(Gene) %>%
    unlist(use.names=FALSE)
  
  final_g <- list(final_g)
  names(final_g) <- celltype
  marker_list <- c(marker_list, final_g)
}

max_l <- max(lengths(marker_list))
marker_list <- lapply(X=marker_list, FUN=function(x) { c(x, rep(x="", times=max_l-length(x)))})
df_final_standard <- data.frame(marker_list)
colnames(df_final_standard) <- names(marker_list)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb = wb, sheetName = "Standard_Markers")
openxlsx::writeData(wb = wb, sheet = "Standard_Markers", x = df_final_standard)
openxlsx::addWorksheet(wb = wb, sheetName = "Intermediate step")
openxlsx::writeData(wb = wb, sheet = "Intermediate step", x = final_df)
openxlsx::saveWorkbook(wb = wb, file = "Evaluated_Markers.xlsx", overwrite = TRUE)


####### We can further refine these markers

# We re-run Seurat::FindAllMarkers() using new groups we have listed above
# to get accurate pct1 and pct2 values and perform a final filtering of these
# markers

for (list in c("Epithelial_1", "Epithelial_2", "Fibroblast_1", "Fibroblast_2", 
               "Lymphatic.Endothelial", "Endothelial", "Granulocytes", "Plasma",
               "B", "TCells", "Macrophages", "DCs", "Mast", "Myofibroblasts")){
  
  
  
}

# df_custom <- df_custom %>%
#   dplyr::filter(avg_log2FC > 0.58, # custom expression cutoff 
#                 p_val_adj < 0.05,  # standard cutoff
#                 pct.1 > pct.2,     # standard cutoff
#                 pct.1 > 0.4) %>%   # custom percent cutoff
#   dplyr::group_by(cluster) %>%
#   dplyr::add_count(cluster) %>%
#   dplyr::filter(n >= 10) %>%
#   dplyr::arrange(desc(avg_log2FC)) # desc(ratio)
#   dplyr::slice_head(n=100) %>%
#   ungroup()

# ############ Create custom grouping of major clusters
# if (proj == "scRNASeq_BBN_C57B6"){
#   clusters <- list("cluster_1" = c(0,22),
#                    "cluster_2" = c(3,16),
#                    "cluster_3" = c(5,12),
#                    "cluster_4" = c(1,6,9),
#                    "cluster_5" = c(2,4,10,11,13),
#                    "cluster_6" = c(7),
#                    "cluster_7" = c(8),
#                    "cluster_8" = c(14),
#                    "cluster_9" = c(15),
#                    "cluster_10" = c(17),
#                    "cluster_11" = c(18),
#                    "cluster_12" = c(19),
#                    "cluster_13" = c(20),
#                    "cluster_14" = c(21),
#                    "cluster_15" = c(23),
#                    "cluster_16" = c(24))
# }
# 
# if (proj == "scRNASeq_BBN_Rag"){
#   clusters <- list("cluster_1" = c(0,4,7,8,9,10,12,13,15,17,18,20,26,27,31),
#                    "cluster_2" = c(3,16,22),
#                    "cluster_3" = c(5,24),
#                    "cluster_4" = c(1,11,23,25,28),
#                    "cluster_5" = c(2,6,21,29,32,33),
#                    "cluster_6" = c(14),
#                    "cluster_7" = c(19),
#                    "cluster_8" = c(30))
# }
# 
# # if (proj == "scRNASeq_Simon"){
# #   clusters <- list("cluster_1" = c(0,2,4,10,11,17,19),
# #                    "cluster_2" = c(3,5,27),
# #                    "cluster_3" = c(1,16,25,28),
# #                    "cluster_4" = c(7,13,15),
# #                    "cluster_5" = c(8,14,20),
# #                    "cluster_6" = c(6),
# #                    "cluster_7" = c(9),
# #                    "cluster_8" = c(12),
# #                    "cluster_9" = c(18),
# #                    "cluster_10" = c(21),
# #                    "cluster_11" = c(22),
# #                    "cluster_12" = c(23),
# #                    "cluster_13" = c(24),
# #                    "cluster_14" = c(26),
# #                    "cluster_15" = c(29))
# # }
# 
# if (proj == "scRNASeq_Jinfen"){
#   clusters <- list("cluster_1" = c(0,6),
#                    "cluster_2" = c(1,3,8),
#                    "cluster_3" = c(4),
#                    "cluster_4" = c(5),
#                    "cluster_5" = c(9),
#                    "cluster_6" = c(2),
#                    "cluster_7" = c(7),
#                    "cluster_8" = c(10),
#                    "cluster_9" = c(11),
#                    "cluster_10" = c(12),
#                    "cluster_11" = c(13),
#                    "cluster_12" = c(14),
#                    "cluster_13" = c(15),
#                    "cluster_14" = c(16))
# }
# 
# if (proj == "scRNASeq_Chen"){
#   clusters <- list("cluster_1" = c(0,1,4,12,13,22),
#                    "cluster_2" = c(2,3,17),
#                    "cluster_3" = c(9,15),
#                    "cluster_4" = c(16,28),
#                    "cluster_5" = c(5),
#                    "cluster_6" = c(6),                     
#                    "cluster_7" = c(7),
#                    "cluster_8" = c(8),
#                    "cluster_9" = c(10),
#                    "cluster_10" = c(11),
#                    "cluster_11" = c(14),
#                    "cluster_12" = c(18),
#                    "cluster_13" = c(19),
#                    "cluster_14" = c(20),
#                    "cluster_15" = c(21),
#                    "cluster_16" = c(23),
#                    "cluster_17" = c(24),
#                    "cluster_18" = c(25),
#                    "cluster_19" = c(26),
#                    "cluster_20" = c(27),
#                    "cluster_21" = c(29))
# }
# 
# integ <- annotate_data(res, reduc, celltype, clusters)
# 
# DefaultAssay(integ) <- "RNA"
# Idents(integ) <- "cell_type"
# 
# all_markers_custom <- Seurat::FindAllMarkers(object=integ,
#                                              assay = "RNA",
#                                              features = NULL,
#                                              logfc.threshold = 0.25,
#                                              test.use = "wilcox",
#                                              slot = "data",
#                                              min.pct = 0.1,
#                                              min.diff.pct = 0.1,
#                                              node = NULL,
#                                              verbose = TRUE,
#                                              only.pos = TRUE,
#                                              max.cells.per.ident = Inf,
#                                              random.seed = 1,
#                                              latent.vars = NULL,
#                                              min.cells.feature = 3,
#                                              min.cells.group = 1,
#                                              pseudocount.use = 1,
#                                              mean.fxn = NULL,
#                                              fc.name = NULL,
#                                              base = 2,
#                                              return.thresh = 0.01,
#                                              densify = FALSE)
# 
# all_markers_custom <- all_markers_custom %>% 
#   dplyr::mutate(pct.1 = dplyr::if_else(pct.1 == 0, 0.001, pct.1),
#                 pct.2 = dplyr::if_else(pct.2 == 0, 0.001, pct.2),
#                 ratio = pct.1/pct.2,
#                 cluster = paste0(proj, "_", cluster)) %>%
#   dplyr::relocate(cluster, gene, avg_log2FC, p_val_adj, pct.1, pct.2, ratio)






l1 <- list("scRNASeq_Jinfen" = c(),
           "scRNASeq_BBN_C57B6" = c(2),
           "scRNASeq_BBN_Rag" = c(8,9,12,15,17,18,26),
           "scRNASeq_Chen" = c(0,13,23,25,1,4,22),
           "scRNASeq_Simon" = c(0,13))

l2 <- list("scRNASeq_Jinfen" = c(10),
           "scRNASeq_BBN_C57B6" = c(4),
           "scRNASeq_BBN_Rag" = c(0),
           "scRNASeq_Chen" = c(),
           "scRNASeq_Simon" = c())

l3 <- list("scRNASeq_Jinfen" = c(),
           "scRNASeq_BBN_C57B6" = c(18),
           "scRNASeq_BBN_Rag" = c(30),
           "scRNASeq_Chen" = c(),
           "scRNASeq_Simon" = c())

l4 <- list("scRNASeq_Jinfen" = c(),
           "scRNASeq_BBN_C57B6" = c(7),
           "scRNASeq_BBN_Rag" = c(5,24),
           "scRNASeq_Chen" = c(5,6,21),
           "scRNASeq_Simon" = c(12))

l5 <- list("scRNASeq_Jinfen" = c(),
           "scRNASeq_BBN_C57B6" = c(),
           "scRNASeq_BBN_Rag" = c(),
           "scRNASeq_Chen" = c(14),
           "scRNASeq_Simon" = c(21))

l6 <- list("scRNASeq_Jinfen" = c(),
           "scRNASeq_BBN_C57B6" = c(23),
           "scRNASeq_BBN_Rag" = c(),
           "scRNASeq_Chen" = c(28),
           "scRNASeq_Simon" = c())

l7 <- list("scRNASeq_Jinfen" = c(3),
           "scRNASeq_BBN_C57B6" = c(),
           "scRNASeq_BBN_Rag" = c(),
           "scRNASeq_Chen" = c(12),
           "scRNASeq_Simon" = c())

l8 <- list("scRNASeq_Jinfen" = c(14),
           "scRNASeq_BBN_C57B6" = c(8),
           "scRNASeq_BBN_Rag" = c(),
           "scRNASeq_Chen" = c(11),
           "scRNASeq_Simon" = c())

l9 <- list("scRNASeq_Jinfen" = c(4,15),
           "scRNASeq_BBN_C57B6" = c(3),
           "scRNASeq_BBN_Rag" = c(),
           "scRNASeq_Chen" = c(2,3,17),
           "scRNASeq_Simon" = c(6))

l10 <- list("scRNASeq_Jinfen" = c(1,8,9),
            "scRNASeq_BBN_C57B6" = c(16),
            "scRNASeq_BBN_Rag" = c(),
            "scRNASeq_Chen" = c(),
            "scRNASeq_Simon" = c())

l11 <- list("scRNASeq_Jinfen" = c(2,12),
            "scRNASeq_BBN_C57B6" = c(5,12),
            "scRNASeq_BBN_Rag" = c(2,6),
            "scRNASeq_Chen" = c(9),
            "scRNASeq_Simon" = c())

l12 <- list("scRNASeq_Jinfen" = c(11,16),
            "scRNASeq_BBN_C57B6" = c(17),
            "scRNASeq_BBN_Rag" = c(),
            "scRNASeq_Chen" = c(),
            "scRNASeq_Simon" = c())

l13 <- list("scRNASeq_Jinfen" = c(0,6),
            "scRNASeq_BBN_C57B6" = c(0,22),
            "scRNASeq_BBN_Rag" = c(3),
            "scRNASeq_Chen" = c(15),
            "scRNASeq_Simon" = c())

l14 <- list("scRNASeq_Jinfen" = c(),
            "scRNASeq_BBN_C57B6" = c(14),
            "scRNASeq_BBN_Rag" = c(19),
            "scRNASeq_Chen" = c(7,26),
            "scRNASeq_Simon" = c(3,5,22))

l15 <- list("scRNASeq_Jinfen" = c(),
            "scRNASeq_BBN_C57B6" = c(6,9),
            "scRNASeq_BBN_Rag" = c(14),
            "scRNASeq_Chen" = c(),
            "scRNASeq_Simon" = c())

l16 <- list("scRNASeq_Jinfen" = c(13),
            "scRNASeq_BBN_C57B6" = c(1),
            "scRNASeq_BBN_Rag" = c(1),
            "scRNASeq_Chen" = c(8),
            "scRNASeq_Simon" = c())

l17 <- list("scRNASeq_Jinfen" = c(),
            "scRNASeq_BBN_C57B6" = c(),
            "scRNASeq_BBN_Rag" = c(28),
            "scRNASeq_Chen" = c(),
            "scRNASeq_Simon" = c(8,14))

# Create a dummy list to store markers from standard approach
marker_list <- list()
final_g <- tolower(unique(df_good$gene))

for (list in c("l1", "l2", "l3", "l4", "l5", "l6", "l7", "l8", "l9", "l10", 
               "l11", "l12", "l13", "l14", "l15", "l16", "l17")){
  
  final_g <- tolower(unique(df_good$gene))
  for (i in 1:length(get(list))){
    
    genes <- df_good %>% 
      dplyr::filter(cluster %in% paste0(names(get(list)[i]),"_",  get(list)[[i]])) %>%
      dplyr::select(gene) %>%
      unlist(use.names = FALSE) %>%
      tolower()
    
    if(length(genes) > 0){
      final_g <- intersect(final_g, genes)
    }
  }
  marker_list <- c(marker_list, list(final_g))
}

max_l <- max(lengths(marker_list))
marker_list <- lapply(X=marker_list, FUN=function(x) { c(x, rep(x="", times=max_l-length(x)))})
df_final_standard <- data.frame(marker_list)

colnames(df_final_standard) <- c("Epithelial_1", "Epithelial_2", "Lymphatic.Endothelial", 
                                 "Endothelial", "Granulocytes", "Plasma Cells", "Mitotic Cells", "B.Cells", "T.Cells",
                                 "NK cells", "Macrophages", "Dendritic.Cells", "Mast Cells", 
                                 "Myofibroblasts", "Fibroblasts_1", "Fibroblasts_2", 
                                 "Fibroblasts_3")

# # For custom approach
# l1 <- list("scRNASeq_Jinfen" = c(),
#            "scRNASeq_BBN_C57B6" = c(8),
#            "scRNASeq_BBN_Rag" = c(7),
#            "scRNASeq_Chen" = c(7,19),
#            "scRNASeq_Simon" = c(2,5,11))
# 
# l2 <- list("scRNASeq_Jinfen" = c(11),
#            "scRNASeq_BBN_C57B6" = c(4,12),
#            "scRNASeq_BBN_Rag" = c(4,6),
#            "scRNASeq_Chen" = c(8),
#            "scRNASeq_Simon" = c(9))
# 
# l3 <- list("scRNASeq_Jinfen" = c(),
#            "scRNASeq_BBN_C57B6" = c(11),
#            "scRNASeq_BBN_Rag" = c(8),
#            "scRNASeq_Chen" = c(),
#            "scRNASeq_Simon" = c())
# 
# l4 <- list("scRNASeq_Jinfen" = c(),
#            "scRNASeq_BBN_C57B6" = c(6),
#            "scRNASeq_BBN_Rag" = c(3),
#            "scRNASeq_Chen" = c(5,6,15),
#            "scRNASeq_Simon" = c(8))
# 
# l5 <- list("scRNASeq_Jinfen" = c(6,10),
#            "scRNASeq_BBN_C57B6" = c(3),
#            "scRNASeq_BBN_Rag" = c(5),
#            "scRNASeq_Chen" = c(3),
#            "scRNASeq_Simon" = c(7))
# 
# l6 <- list("scRNASeq_Jinfen" = c(),
#            "scRNASeq_BBN_C57B6" = c(5),
#            "scRNASeq_BBN_Rag" = c(1),
#            "scRNASeq_Chen" = c(1,16,18),
#            "scRNASeq_Simon" = c(1,4))
# 
# l7 <- list("scRNASeq_Jinfen" = c(1),
#            "scRNASeq_BBN_C57B6" = c(1),
#            "scRNASeq_BBN_Rag" = c(2),
#            "scRNASeq_Chen" = c(),
#            "scRNASeq_Simon" = c())
# 
# l8 <- list("scRNASeq_Jinfen" = c(),
#            "scRNASeq_BBN_C57B6" = c(),
#            "scRNASeq_BBN_Rag" = c(),
#            "scRNASeq_Chen" = c(11),
#            "scRNASeq_Simon" = c(10))
# 
# l9 <- list("scRNASeq_Jinfen" = c(12),
#             "scRNASeq_BBN_C57B6" = c(7),
#             "scRNASeq_BBN_Rag" = c(),
#             "scRNASeq_Chen" = c(10),
#             "scRNASeq_Simon" = c())
# 
# l10 <- list("scRNASeq_Jinfen" = c(3,13),
#             "scRNASeq_BBN_C57B6" = c(2),
#             "scRNASeq_BBN_Rag" = c(2,17),
#             "scRNASeq_Chen" = c(6),
#             "scRNASeq_Simon" = c())
# 
# l11 <- list("scRNASeq_Jinfen" = c(9,14),
#             "scRNASeq_BBN_C57B6" = c(10),
#             "scRNASeq_BBN_Rag" = c(),
#             "scRNASeq_Chen" = c(),
#             "scRNASeq_Simon" = c())
# 
# # Create a dummy list to store markers from custom approach
# marker_list <- list()
# final_g <- tolower(unique(df_custom$gene))
# 
# for (list in c("l1", "l2", "l3", "l4", "l5", "l6", "l7", "l8", "l9", "l10", 
#                "l11")){
#   
#   final_g <- tolower(unique(df_custom$gene))
#   for (i in 1:length(get(list))){
#     
#     genes <- df_custom %>% 
#       dplyr::filter(cluster %in% paste0(names(get(list)[i]),"_cluster_",  get(list)[[i]])) %>%
#       dplyr::select(gene) %>%
#       unlist(use.names = FALSE) %>%
#       tolower()
#     
#     if(length(genes) > 0){
#       final_g <- intersect(final_g, genes)
#     }
#   }
#   marker_list <- c(marker_list, list(final_g))
# }
# 
# max_l <- max(lengths(marker_list))
# marker_list <- lapply(X=marker_list, FUN=function(x) { c(x, rep(x="", times=max_l-length(x)))})
# df_final_custom <- data.frame(marker_list)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb = wb, sheetName = "Standard_Markers")
openxlsx::writeData(wb = wb, sheet = "Standard_Markers", x = df_final_standard)
# openxlsx::addWorksheet(wb = wb, sheetName = "Custom_Markers")
# openxlsx::writeData(wb = wb, sheet = "Custom_Markers", x = df_final_custom)
openxlsx::saveWorkbook(wb = wb, file = "Evaluated_Markers.xlsx", overwrite = TRUE)

#******************************************************************************#
#      STEP 5: VALIDATE THE SPECIFICITY FOR EACH MARKER GENE USING UMAPS       #
#******************************************************************************#

for (proj in c("scRNASeq_BBN_C57B6", "scRNASeq_BBN_Rag", "scRNASeq_Jinfen", 
               "scRNASeq_Chen", "scRNASeq_Simon")){
  
  seurat_results <- paste0("/hpc/home/kailasamms/scratch/", proj, "/results_seurat/")
  res <- 0.4
  reduc <- "rpca" #"harmony"
  celltype <- NULL
  #genes <- read.xlsx("/hpc/home/kailasamms/Markers.xlsx")
  genes <- read.xlsx("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Markers.xlsx")
  
  # Load the integrated seurat object
  integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn", 
                                      dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  # Change active.ident to major clusters we annotated
  DefaultAssay(integrated_seurat) <- "RNA"
  idents <- paste0("cluster.", res, ".", base::tolower(reduc))
  Idents(object = integrated_seurat) <- idents
  
  feature_plot <- function(i){
    
    split <- "Condition" #NULL
    split <- NULL
    Seurat::FeaturePlot(object = integrated_seurat,
                        slot = "data",
                        features = i,
                        split.by = split,
                        cols =  c("#D1E5F0", "#F7F7F7", "#67001F"),
                        pt.size = 0.4,
                        order = TRUE,
                        min.cutoff = 'q10',
                        reduction = paste0("umap.", base::tolower(reduc)),
                        label = TRUE,
                        combine = TRUE,
                        raster = FALSE) +
      scale_colour_gradientn(colours = c("#D1E5F0", "#F7F7F7", "#67001F"), limits=c(0,2))
    # rev(brewer.pal(n = 11, name = "RdBu"))[5:11]
  } 
  
  for (i in 1:ncol(genes)){
    features <- genes[,i] %>% unlist(use.names=FALSE)
    features <- sort(rownames(integrated_seurat@assays$RNA$data)[tolower(rownames(integrated_seurat@assays$RNA$data)) %in% tolower(features)])
    
    # Plot UMAPs
    for (n in 1:ceiling((length(features)/12))){
      j <- 12*n-11
      k <- 12*n
      vec <- features[j:k]
      purrr::map(.x = vec[!is.na(vec)], 
                 .f = feature_plot) %>%
        cowplot::plot_grid(plotlist = .,
                           align = "hv",
                           axis = "tblr",
                           nrow = 3,
                           ncol = 4,
                           rel_widths = 1,
                           rel_heights = 1,
                           greedy = TRUE,
                           byrow = TRUE)
      
      ggplot2::ggsave(filename = paste0(proj,"_", colnames(genes)[i], "_", n, ".jpg"),
                      plot = last_plot(),
                      device = "jpeg",
                      scale = 1,
                      width = dplyr::if_else(is.null(split), 9*3, 9*6),
                      height = 11,
                      units = c("in"),
                      dpi = 300,
                      limitsize = FALSE,
                      bg = "white")
    }
  }
}


#******************************************************************************#

# # Find out clusters that are similar and genes that are common
# combinations_standard <- utils::combn(x=unique(df_standard$cluster), m=2)
# combinations_custom <- utils::combn(x=unique(df_custom$cluster), m=2)
# final_standard <- list()
# final_custom <- list()
# 
# for (i in 1:ncol(combinations_standard)){
#   
#   l <- intersect(df_standard %>% 
#                    dplyr::filter(cluster == combinations_standard[1,i]) %>%
#                    ungroup() %>%
#                    dplyr::select(gene) %>%
#                    unlist(use.names=FALSE) %>%
#                    tolower(),
#                  df_standard %>% 
#                    dplyr::filter(cluster == combinations_standard[2,i]) %>%
#                    ungroup() %>%
#                    dplyr::select(gene) %>%
#                    unlist(use.names=FALSE) %>%
#                    tolower())
#   
#   # Custom cutoff of 20 genes to be similar between 2 clusters
#   if (length(l) > 20){
#     l <-  list(l)
#     names(l) <- paste0(combinations_standard[1,i],".", combinations_standard[2,i])
#     final_standard <- c(final_standard, l)
#   }
# }
# 
# for (i in 1:ncol(combinations_custom)){
#   
#   l <- intersect(df_custom %>% 
#                    dplyr::filter(cluster == combinations_custom[1,i]) %>%
#                    ungroup() %>%
#                    dplyr::select(gene) %>%
#                    unlist(use.names=FALSE) %>%
#                    tolower(),
#                  df_custom %>% 
#                    dplyr::filter(cluster == combinations_custom[2,i]) %>%
#                    ungroup() %>%
#                    dplyr::select(gene) %>%
#                    unlist(use.names=FALSE) %>%
#                    tolower())
#   
#   # Custom cutoff of 20 genes to be similar between 2 clusters
#   if (length(l) > 20){
#     l <-  list(l)
#     names(l) <- paste0(combinations_custom[1,i],".", combinations_custom[2,i])
#     final_custom <- c(final_custom, l)
#   }
# }
# 
# # Convert list to dataframe
# for (i in 1:length(final_standard)){
#   n <- names(final_standard[i])
#   final_standard[i] <- list(c(unlist(final_standard[i], use.names=FALSE), rep(x="", times=max(lengths(final_standard))-length(final_standard[[i]]))))
#   names(final_standard[i]) <- names(final_standard[i])
# }
# 
# for (i in 1:length(final_custom)){
#   n <- names(final_custom[i])
#   final_custom[i] <- list(c(unlist(final_custom[i], use.names=FALSE), rep(x="", times=max(lengths(final_custom))-length(final_custom[[i]]))))
#   names(final_custom[i]) <- names(final_custom[i])
# }
# 
# common_markers_standard <- data.frame(final_standard)
# common_markers_custom <- data.frame(final_custom)
# 
# # Save all the markers and correlation
# wb <- openxlsx::createWorkbook()
# openxlsx::addWorksheet(wb = wb, sheetName = "Standard_Markers")
# openxlsx::writeData(wb = wb, sheet = "Standard_Markers", x = common_markers_standard )
# openxlsx::addWorksheet(wb = wb, sheetName = "Custom_Markers")
# openxlsx::writeData(wb = wb, sheet = "Custom_Markers", x = common_markers_custom)
# openxlsx::saveWorkbook(wb = wb, file = "Compiled_Markers.xlsx", overwrite = TRUE)

# This script is used to compare cluster markers from multiple NGS datasets
# from time to time to generate a master list of markers that can be used to 
# annotate single cell datasets

# We will use the following approach to generate our own markers:

# (1) FindAllMarkers() using harmony at resolution of 0.4
# NOTE: Even a low resolution, you will often find that epithelial cluster has 
# multiple subclusters. We can manually group these subclusters of epithelial
# cells as well as other celltypes and re-run FindAllMarkers() on the new 
# clusters. However, this takes a lot of time when you keep adding datasets over
# time and it becomes troublesome to manually group subclusters into megaclusters.


# (2) We will next identify similar clusters across datasets
# Keep ONLY markers that are 
# (i) highly expressed log2FC >= 0.58 and 
# (ii) padj <=0.05 and 
# (iii) expressed in most of cells in cluster pct.1 >=0.5
# We are not using any cutoffs on pct.2 because like pointed out earlier
# there are multiple epithelial subclusters. So, epithelial markers will be 
# expressed in most of cells in other epithelial subclusters and pct.2 will be
# high. Setting cutoffs like pct.2 < 0.2 will lead to loss of good markers.

# NOTE: If you use all markers you get from output of Seurat::FindAllMarkers() 
# directly to identify similar clusters across datasets, results are completely 
# wrong. YOU MUST FILTER OUT POOR/WEAK markers and use ONLY top 100 markers for 
# each cluster to identify similar clusters. I have observed that the more 
# markers you use, the poorer the ability to identify similar clusters 
# across datasets.
# NOTE: I tried an alternative way to identify similar clusters across datasets.
# I got normalized counts for cells from each cluster in every dataset
# and calculated the mean normalized counts for each gene for every cluster. 
# Then, I removed the genes with zero counts in all samples and scaled the values
# such that each column i.e. cluster has mean 0 and std 1. Then, I calculated
# the euclidean distance between the clusters and identified similar clusters.
# The results were completely wrong. Hence, it is better to find identical 
# clusters across datasets using Seurat::FindAllMarkers(). 

#******************************************************************************#
#                                STEP 1: FIND MARKERS                          #
#******************************************************************************#

proj <- "scRNASeq_BBN_C57B6"
source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")

# "scRNASeq_Simon" is not a good dataset. So, we ignore it.
projects <- c("scRNASeq_BBN_C57B6", 
              "scRNASeq_BBN_Rag", 
              "scRNASeq_Chen", 
              "scRNASeq_GSE164557",
              "scRNASeq_GSE217093", 
              "scRNASeq_GSE222315",
              "scRNASeq_HRA003620",
              "scRNASeq_Jinfen",
              "scRNASeq_NA13_CDH12_C57B6")

# Create dataframe to store output of FindMarkers()
df_standard <- data.frame(cluster = "", 
                          gene = "",
                          avg_log2FC = 0,
                          p_val_adj = 0,
                          pct.1 = 0.0, 
                          pct.2 = 0.0,
                          ratio = 0.0)

# Iterate through all datasets and find markers
for (proj in projects){
  
  seurat_results <- paste0("/hpc/home/kailasamms/scratch/", proj, "/results_seurat/")
  res <- 0.4
  reduc <- "harmony"
  celltype <- NULL
  
  # Load the integrated seurat object generated using "scRNASeq_Seurat_Analysis_Import_Data"
  integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn", 
                                      dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  # Change active.ident to major clusters we annotated
  DefaultAssay(integrated_seurat) <- "RNA"
  idents <- paste0("cluster.", res, ".", base::tolower(reduc))
  Idents(object = integrated_seurat) <- idents
  
  # Find ALL markers
  all_markers_standard <- Seurat::FindAllMarkers(object=integrated_seurat,
                                                 assay = "RNA",
                                                 features = NULL,
                                                 logfc.threshold = 0.25,
                                                 test.use = "wilcox",
                                                 slot = "data",
                                                 min.pct = 0.1,
                                                 min.diff.pct = 0.1,
                                                 node = NULL,
                                                 verbose = TRUE,
                                                 only.pos = TRUE,
                                                 max.cells.per.ident = Inf,
                                                 random.seed = 1,
                                                 latent.vars = NULL,
                                                 min.cells.feature = 3,
                                                 min.cells.group = 1,
                                                 pseudocount.use = 1,
                                                 mean.fxn = NULL,
                                                 fc.name = NULL,
                                                 base = 2,
                                                 return.thresh = 0.01,
                                                 densify = FALSE)
  
  all_markers_standard <- all_markers_standard %>% 
    dplyr::mutate(pct.1 = dplyr::if_else(pct.1 == 0, 0.001, pct.1),
                  pct.2 = dplyr::if_else(pct.2 == 0, 0.001, pct.2),
                  ratio = pct.1/pct.2,
                  cluster = paste0(proj, "_", cluster)) %>%
    dplyr::relocate(cluster, gene, avg_log2FC, p_val_adj, pct.1, pct.2, ratio)
  
  df_standard <- dplyr::bind_rows(df_standard, all_markers_standard)
  
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "All markers")
  openxlsx::writeData(wb, sheet = "All markers", x = df_standard, rowNames = FALSE)
  openxlsx::saveWorkbook(wb = wb, file = paste0("All Markers ", length(projects), " datasets.xlsx"), overwrite = TRUE)
}

#******************************************************************************#
#                     STEP 2: FILTER OUT BAD MARKERS                           #
#******************************************************************************#

df_standard <- read.xlsx("All Markers for 7 datasets.xlsx")

# Filter out weak/non-specific markers & identify top 100 markers for every cluster
df_standard <- df_standard %>%
  dplyr::filter(p_val_adj <= 0.05 & 
                  pct.1 >= 0.5 & 
                  avg_log2FC >= 0.58) %>%
  #(pct.1 >= 0.7 | ratio >= 2 | avg_log2FC >= 0.58 )) %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(desc(avg_log2FC)) %>%  # log2FC better than ratio
  dplyr::slice_head(n=50) %>%     
  ungroup()

#******************************************************************************#
#            STEP 3: IDENTIFY SIMILAR CLUSTERS ACROSS DATASETS                 #
#******************************************************************************#

# Create similarity matrix.
# Rows and columns MUST have dataset_cluster_number
# Values indicate number of common genes between any 2 clusters across datasets
sim_matrix <- matrix(data=0, 
                     nrow=length(unique(df_standard$cluster)), 
                     ncol=length(unique(df_standard$cluster)))
colnames(sim_matrix) <- unique(df_standard$cluster)
rownames(sim_matrix) <- unique(df_standard$cluster)

for (i in 1:ncol(sim_matrix)){
  for (j in 1:nrow(sim_matrix)){
    sim_matrix[i,j] <- length(intersect(tolower(df_standard %>% 
                                                  dplyr::filter(cluster == colnames(sim_matrix)[i]) %>%
                                                  dplyr::select(gene) %>%
                                                  unlist(use.names = FALSE)),
                                        tolower(df_standard %>% 
                                                  dplyr::filter(cluster == rownames(sim_matrix)[j]) %>%
                                                  dplyr::select(gene) %>%
                                                  unlist(use.names = FALSE))))
  }
}

# Group similar clusters together
mat <- sim_matrix
#rowclust <- stats::hclust(d = dist(mat))
#reordered <- mat[rowclust$order,]
colclust <- stats::hclust(d = dist(t(mat)))
reordered <- mat[colclust$order, colclust$order]

# Assuming we have 22 different cell types, set k=22
groups<-cutree(colclust, k=22)
mat1 <-cbind(mat,groups)

l <- 0
for (i in 1:22){
  
  x1 <- subset(mat1, groups==i)
  l1 <- length(rownames(x1))
  l <- max(l, l1)
}
g <- list()
for (i in 1:22){
  
  x1 <- subset(mat1, groups==i)
  l1 <- length(rownames(x1))
  r1 <- c(rownames(x1), rep(x="", times=l-l1))
  g <- c(g, list(r1))
}
names(g) <- paste0("Group", seq(1:22))

# Save the clustered similarity matrix
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Heatmap_groups")
openxlsx::writeData(wb, sheet = "Heatmap_groups", x = data.frame(g), rowNames = FALSE)
openxlsx::addWorksheet(wb, sheetName = "Heatmap_standard")
openxlsx::writeData(wb, sheet = "Heatmap_standard", x = reordered, rowNames = TRUE)
openxlsx::saveWorkbook(wb = wb, file = "Compiled_Markers.xlsx", overwrite = TRUE)

#******************************************************************************#
#   STEP 4: IDENTIFY COMMOM GENES AMONG SIMILAR CLUSTERS ACROSS DATASETS       #
#******************************************************************************#

# Look at heatmap matrix excel file as well as UMAP at resolution 0.4 for each
# dataset. Identify similar clusters across all experiments as demonstrated in
# the pptx.

# Group the similar clusters in new sheet "Heatmap_groups_final" and give them 
# proper celltype names.

# Create new macro clusters for each dataset based on the above classification
# and again redo step 1 to 3.

# "scRNASeq_Simon" is not a good dataset. So, we ignore it.
projects <- c("scRNASeq_BBN_C57B6", 
              "scRNASeq_BBN_Rag", 
              "scRNASeq_Jinfen", 
              "scRNASeq_Chen", 
              "scRNASeq_GSE217093", 
              "scRNASeq_GSE222315",
              "scRNASeq_HRA003620")

class <- read.xlsx("Compiled_Markers.xlsx", 
                   sheet = "Heatmap_groups_final")

# Create dataframe to store output of FindMarkers()
df_custom <- data.frame(cluster = "", 
                        gene = "",
                        avg_log2FC = 0,
                        p_val_adj = 0,
                        pct.1 = 0.0, 
                        pct.2 = 0.0,
                        ratio = 0.0)

# Iterate through all datasets and find markers
for (proj in projects){
  
  seurat_results <- paste0("/hpc/home/kailasamms/scratch/", proj, "/results_seurat/")
  res <- 0.4
  reduc <- "harmony"
  celltype <- NULL
  
  # Load the integrated seurat object generated using "scRNASeq_Seurat_Analysis_Import_Data"
  integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn", 
                                      dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  integrated_seurat@meta.data <- integrated_seurat@meta.data %>%
    dplyr::mutate(New_class = cluster.0.4.harmony)
  
  dataset <- gsub("scRNASeq_|BBN_", "", proj)
  for (col in 1:ncol(class)){
    
    clusters <- class[,col][grepl(dataset, class[,col])]
    clusters <- as.numeric(gsub(paste0(dataset,"_"), "", clusters))
    cat(colnames(class)[col], ":", clusters, "\n")
    
    integrated_seurat@meta.data <- integrated_seurat@meta.data %>%
      dplyr::mutate(New_class = dplyr::case_when(New_class %in% clusters ~ colnames(class)[col],
                                                 TRUE ~ New_class))
  }
  
  # Change active.ident to major clusters we annotated
  DefaultAssay(integrated_seurat) <- "RNA"
  idents <- "New_class"
  Idents(object = integrated_seurat) <- idents
  
  # Find ALL markers
  all_markers_standard <- Seurat::FindAllMarkers(object=integrated_seurat,
                                                 assay = "RNA",
                                                 features = NULL,
                                                 logfc.threshold = 0.25,
                                                 test.use = "wilcox",
                                                 slot = "data",
                                                 min.pct = 0.1,
                                                 min.diff.pct = 0.1,
                                                 node = NULL,
                                                 verbose = TRUE,
                                                 only.pos = TRUE,
                                                 max.cells.per.ident = Inf,
                                                 random.seed = 1,
                                                 latent.vars = NULL,
                                                 min.cells.feature = 3,
                                                 min.cells.group = 1,
                                                 pseudocount.use = 1,
                                                 mean.fxn = NULL,
                                                 fc.name = NULL,
                                                 base = 2,
                                                 return.thresh = 0.01,
                                                 densify = FALSE)
  
  all_markers_standard <- all_markers_standard %>% 
    dplyr::mutate(pct.1 = dplyr::if_else(pct.1 == 0, 0.001, pct.1),
                  pct.2 = dplyr::if_else(pct.2 == 0, 0.001, pct.2),
                  ratio = pct.1/pct.2,
                  cluster = paste0(proj, "_", cluster)) %>%
    dplyr::relocate(cluster, gene, avg_log2FC, p_val_adj, pct.1, pct.2, ratio)
  
  df_custom <- dplyr::bind_rows(df_custom, all_markers_standard)
  
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "Custom markers")
  openxlsx::writeData(wb, sheet = "Custom markers", x = df_custom, rowNames = FALSE)
  openxlsx::saveWorkbook(wb = wb, file = "Custom Markers for 7 datasets.xlsx", overwrite = TRUE)
}

df_custom <- read.xlsx("Custom Markers for 7 datasets.xlsx")
dim(df_custom)

# Filter out weak/non-specific markers & identify top 100 markers for every cluster
df_custom <- df_custom %>%
  dplyr::filter(p_val_adj <= 0.05 & 
                  pct.1 > pct.2 & 
                  avg_log2FC >= 0.58) %>%
  #(pct.1 >= 0.7 | ratio >= 2 | avg_log2FC >= 0.58 )) %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(desc(avg_log2FC)) %>%  # lof2FC better than ratio
  dplyr::slice_head(n=100) %>%     
  ungroup()
dim(df_custom)

# Keep ONLY clusters which contain Epithelial, Endothelial etc
df_custom <- df_custom %>%
  dplyr::mutate(cluster = gsub(pattern="scRNASeq_|BBN_", replacement="", x=cluster)) %>%
  dplyr::mutate(cluster = gsub(pattern="_1", replacement="1", x=cluster)) %>%
  dplyr::mutate(cluster = gsub(pattern="_2", replacement="2", x=cluster)) %>%
  dplyr::mutate(cluster = gsub(pattern="_type", replacement="type", x=cluster)) %>%
  dplyr::filter(grepl(pattern="Mast|Macrophages|Dendritic|_B|Plasma|T|Endothelial|Myofibroblasts|Epithelial|Fibroblasts|Novel",
                      x=cluster)) %>%
  tidyr::separate_wider_delim(cols=cluster, delim="_", names=c("dataset","celltype")) %>%
  dplyr::mutate(gene = tolower(gene)) %>%  #make human and mouse genes identical before adding counts
  dplyr::add_count(celltype,gene) %>%
  dplyr::rename(Freq=n)
dim(df_custom)

# Calculate total datasets for each celltype. 
# Eg: Myofibroblast cells were defined from 3 datasets while Mast cells were
# defined from 7 datasets
total <- c()
for (i in 1:nrow(df_custom)){
  total <- c(total, length(unique(df_custom %>% dplyr::filter(celltype == df_custom$celltype[i]) %>% dplyr::select(dataset) %>% unlist(use.names=FALSE))))
}

df_custom$Total <- total
dim(df_custom)

# Add percent and remove duplicate markers within each celltype
df_custom <- df_custom %>%
  dplyr::mutate(Percent = 100*Freq/Total) %>%
  #dplyr::filter(Freq > 1) %>%  
  dplyr::filter(Freq > 1, pct.2 < 0.1) %>%  # reduce non-specific markers
  dplyr::select(everything(), -dataset) %>%
  dplyr::group_by(celltype) %>%
  dplyr::distinct_at("gene", .keep_all = TRUE) %>%
  dplyr::ungroup()
dim(df_custom)

marker_list <- list()
for (c in unique(df_custom$celltype)){
  
  g <- df_custom %>% 
    dplyr::filter(celltype == c) %>%
    dplyr::select(gene) %>%
    unlist(use.names=FALSE) %>%
    list()
  names(g) <- c
  
  marker_list <- c(marker_list, g)
  
  
}
max_l <- max(lengths(marker_list))
marker_list <- lapply(X=marker_list, FUN=function(x) { c(x, rep(x="", times=max_l-length(x)))})
df_final <- data.frame(marker_list)

# Save the clustered similarity matrix
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Markers")
openxlsx::writeData(wb, sheet = "Markers", x = df_final, rowNames = FALSE)
openxlsx::saveWorkbook(wb = wb, file = "Evaluated_Markers.xlsx", overwrite = TRUE)













































# # Create similarity matrix with number of common genes between any 2 clusters across datasets
# sim_matrix <- matrix(data=0, 
#                      nrow=length(unique(df_custom$cluster)), 
#                      ncol=length(unique(df_custom$cluster)))
# colnames(sim_matrix) <- unique(df_custom$cluster)
# rownames(sim_matrix) <- unique(df_custom$cluster)
# 
# for (i in 1:ncol(sim_matrix)){
#   for (j in 1:nrow(sim_matrix)){
#     sim_matrix[i,j] <- length(intersect(tolower(df_custom %>% 
#                                                   dplyr::filter(cluster == colnames(sim_matrix)[i]) %>%
#                                                   dplyr::select(gene) %>%
#                                                   unlist(use.names = FALSE)),
#                                         tolower(df_custom %>% 
#                                                   dplyr::filter(cluster == rownames(sim_matrix)[j]) %>%
#                                                   dplyr::select(gene) %>%
#                                                   unlist(use.names = FALSE))))
#   }
# }
# 
# # Group similar clusters together
# mat <- sim_matrix
# #rowclust <- stats::hclust(d = dist(mat))
# #reordered <- mat[rowclust$order,]
# colclust <- stats::hclust(d = dist(t(mat)))
# reordered <- mat[colclust$order, colclust$order]
# 
# # Assuming we have 22 different cell types, set k=22
# groups<-cutree(colclust, k=22)
# mat1 <-cbind(mat,groups)
# 
# l <- 0
# for (i in 1:22){
#   
#   x1 <- subset(mat1, groups==i)
#   l1 <- length(rownames(x1))
#   l <- max(l, l1)
# }
# g <- list()
# for (i in 1:22){
#   
#   x1 <- subset(mat1, groups==i)
#   l1 <- length(rownames(x1))
#   r1 <- c(rownames(x1), rep(x="", times=l-l1))
#   g <- c(g, list(r1))
# }
# names(g) <- paste0("Group", seq(1:22))
# 
# # Save the clustered similarity matrix
# wb <- openxlsx::createWorkbook()
# openxlsx::addWorksheet(wb, sheetName = "Heatmap_groups")
# openxlsx::writeData(wb, sheet = "Heatmap_groups", x = data.frame(g), rowNames = FALSE)
# openxlsx::addWorksheet(wb, sheetName = "Heatmap_standard")
# openxlsx::writeData(wb, sheet = "Heatmap_standard", x = reordered, rowNames = TRUE)
# openxlsx::saveWorkbook(wb = wb, file = "Compiled_Markers_custom.xlsx", overwrite = TRUE)


##############
# 
# df_custom <- read.xlsx("Custom Markers for 7 datasets.xlsx")
# df_custom_groups <- read.xlsx("Compiled_Markers_custom.xlsx")
# 
# # Create a dummy dataframe to store markers
# final_df <- data.frame("Gene"=unique(tolower(df_custom$gene)))
# 
# # Iterate through each celltype
# for (i in 1:ncol(df_custom_groups)){
#   
#   # Create a list to store markers of each cell type
#   final_g <- c()
#   
#   clusters <- df_custom_groups[,i]
#   datasets <- unique(gsub(pattern="_[0-9]*$",replacement="",x=clusters))
#   
#   for (d in datasets){
#     # Get unique markers for a cell type from each dataset
#     genes <- df_custom %>% 
#       dplyr::filter(grepl(d, cluster)) %>%
#       dplyr::filter(cluster %in% clusters) %>%
#       dplyr::select(gene) %>%
#       unlist(use.names = FALSE) %>%
#       tolower() %>%
#       unique()
#     
#     # Merge all unique markers from each dataset
#     final_g <- c(final_g, genes)
#   }
#   
#   # Calculate number of datasets each marker was identified
#   n_data <- length(datasets)  # number of datasets used for each cell type
#   final_g <-final_g %>% 
#     table() %>% 
#     data.frame() %>%
#     dplyr::mutate(Percent = Freq*100/n_data) %>%
#     dplyr::filter(Freq > 1, Percent >= 50) %>%
#     dplyr::arrange(desc(Percent))
#   colnames(final_g) <- paste0(colnames(df_custom_groups)[i],"_", colnames(final_g))
#   final_g$Gene <- final_g[,1]
#   
#   # Merge all markers from different cell types
#   final_df <- final_df %>%
#     dplyr::left_join(final_g,by=c("Gene"="Gene"))
# }
# 
# final_df[is.na(final_df)] <- 0
# 
# # If a gene is detected as marker for 2 cell types, use it as marker for the 
# # celltype where it has highest percent
# marker_list <- list()
# percent_cols <- colnames(final_df)[grepl("Percent", colnames(final_df))] 
# freq_cols <- colnames(final_df)[grepl("Freq", colnames(final_df))]
# celltype_cols <- setdiff(colnames(final_df), c("Gene", percent_cols, freq_cols))
# 
# for (i in 1:length(celltype_cols)){
#   
#   f <- final_df %>% 
#     dplyr::filter(get(percent_cols[i]) > 50) %>%
#     group_by(Gene) %>%
#     dplyr::mutate(max_percent = max(!!!rlang::syms(percent_cols), na.rm=TRUE),
#                   max_count = max(!!!rlang::syms(freq_cols), na.rm=TRUE)) %>%
#     dplyr::filter(get(percent_cols[i]) == max_percent) %>%
#     dplyr::select(Gene, all_of(celltype), all_of(paste0(celltype, "_count")), max_percent) %>%
#     dplyr::arrange(desc(max_percent)) %>%
#     data.frame() %>%
#     #dplyr::slice_head(n=50) %>%
#     dplyr::select(Gene) %>%
#     unlist(use.names=FALSE)
#   
#   final_g <- list(final_g)
#   names(final_g) <- celltype
#   marker_list <- c(marker_list, final_g)
# }
# 
# max_l <- max(lengths(marker_list))
# marker_list <- lapply(X=marker_list, FUN=function(x) { c(x, rep(x="", times=max_l-length(x)))})
# df_final_standard <- data.frame(marker_list)
# colnames(df_final_standard) <- names(marker_list)
# 
# wb <- openxlsx::createWorkbook()
# openxlsx::addWorksheet(wb = wb, sheetName = "Standard_Markers")
# openxlsx::writeData(wb = wb, sheet = "Standard_Markers", x = df_final_standard)
# openxlsx::addWorksheet(wb = wb, sheetName = "Intermediate step")
# openxlsx::writeData(wb = wb, sheet = "Intermediate step", x = final_df)
# openxlsx::saveWorkbook(wb = wb, file = "Evaluated_Markers.xlsx", overwrite = TRUE)


for (list in c("Epithelial_1", "Epithelial_2", "Fibroblast_1", "Fibroblast_2", 
               "Lymphatic.Endothelial", "Endothelial", "Granulocytes", "Plasma",
               "B", "TCells", "Macrophages", "DCs", "Mast", "Myofibroblasts")){
  
  # Create a list to store markers of each cell type
  final_g <- c()
  
  # Iterate through each dataset
  for (i in 1:length(get(list))){
    
    # Get unique markers for a cell type from each dataset
    genes <- df_custom %>% 
      dplyr::filter(cluster %in% paste0(names(get(list)[i]),"_",  get(list)[[i]])) %>%
      dplyr::select(gene) %>%
      unlist(use.names = FALSE) %>%
      tolower() %>%
      unique()
    
    # Merge all unique markers from each dataset
    final_g <- c(final_g, genes)
  }
  # Calculate number of datasets each marker was identified
  n_data <- sum(lengths(get(list)) > 0)  # number of datasets used for each cell type
  final_g <- data.frame(table(final_g))
  final_g <- final_g %>%
    dplyr::mutate(Percent = Freq*100/n_data)
  colnames(final_g) <- c(list, paste0(list, "_count"), paste0(list, "_percent"))
  
  # Remove markers identified in only a single dataset
  final_g <- final_g %>% 
    dplyr::filter(get(paste0(list, "_count")) > 1) %>%
    dplyr::mutate(Gene = get(list))
  
  # Merge all markers for each cell type present in multiple daasets to final dataframe
  final_df <- final_df %>%
    left_join(final_g,by=c("Gene"="Gene"))
}

final_df[is.na(final_df)] <- 0

# If a gene is detected as marker for 2 cell types, use it as marker for the 
# celltype where it has highest count
marker_list <- list()

for (celltype in c("Epithelial_1", "Epithelial_2", "Fibroblast_1", "Fibroblast_2", 
                   "Lymphatic.Endothelial", "Endothelial", "Granulocytes", "Plasma",
                   "B", "TCells", "Macrophages", "DCs", "Mast", "Myofibroblasts")){
  
  final_g <- final_df %>% 
    group_by(Gene) %>%
    dplyr::mutate(max_percent = max(Epithelial_1_percent, Epithelial_2_percent, 
                                    Fibroblast_1_percent, Fibroblast_2_percent, 
                                    Lymphatic.Endothelial_percent, Endothelial_percent,
                                    Granulocytes_percent, Plasma_percent,
                                    B_percent, TCells_percent, Macrophages_percent,
                                    DCs_percent, Mast_percent, Myofibroblasts_percent, na.rm=TRUE),
                  max_count = max(Epithelial_1_count, Epithelial_2_count, 
                                  Fibroblast_1_count, Fibroblast_2_count, 
                                  Lymphatic.Endothelial_count, Endothelial_count,
                                  Granulocytes_count, Plasma_count,
                                  B_count, TCells_count, Macrophages_count,
                                  DCs_count, Mast_count, Myofibroblasts_count, na.rm=TRUE)) %>%
    dplyr::filter(get(paste0(celltype, "_percent")) == max_percent) %>%
    dplyr::filter(get(paste0(celltype, "_percent")) >= 50) %>%
    dplyr::select(Gene, all_of(celltype), all_of(paste0(celltype, "_count")), max_percent) %>%
    dplyr::arrange(desc(max_percent)) %>%
    data.frame() %>%
    #dplyr::slice_head(n=50) %>%
    dplyr::select(Gene) %>%
    unlist(use.names=FALSE)
  
  final_g <- list(final_g)
  names(final_g) <- celltype
  marker_list <- c(marker_list, final_g)
}

max_l <- max(lengths(marker_list))
marker_list <- lapply(X=marker_list, FUN=function(x) { c(x, rep(x="", times=max_l-length(x)))})
df_final_standard <- data.frame(marker_list)
colnames(df_final_standard) <- names(marker_list)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb = wb, sheetName = "Standard_Markers")
openxlsx::writeData(wb = wb, sheet = "Standard_Markers", x = df_final_standard)
openxlsx::addWorksheet(wb = wb, sheetName = "Intermediate step")
openxlsx::writeData(wb = wb, sheet = "Intermediate step", x = final_df)
openxlsx::saveWorkbook(wb = wb, file = "Evaluated_Markers.xlsx", overwrite = TRUE)

#******************************************************************************#
#      STEP 5: VALIDATE THE SPECIFICITY FOR EACH MARKER GENE USING UMAPS       #
#******************************************************************************#

for (proj in c("scRNASeq_BBN_C57B6", "scRNASeq_BBN_Rag", "scRNASeq_Jinfen",
               "scRNASeq_Chen", "scRNASeq_GSE217093", "scRNASeq_GSE222315",
               "scRNASeq_HRA003620")){
  
  
  seurat_results <- paste0("/hpc/home/kailasamms/scratch/", proj, "/results_seurat/")
  res <- 0.4
  reduc <- "harmony"
  celltype <- NULL
  idents <- paste0("cluster.", res, ".", base::tolower(reduc))
  
  # Load the integrated seurat object generated using "scRNASeq_Seurat_Analysis_Import_Data"
  integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn", 
                                      dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  # Set default assay
  DefaultAssay(integrated_seurat) <- "RNA"
  
  # Set identity to an existing column in meta data
  Idents(object=integrated_seurat) <- idents
  
  df_final <- read.xlsx("Evaluated_Markers.xlsx")
  
  for (col in 1:ncol(df_final)){
    
    features <- unlist(df_final[col], use.names=FALSE)
    features <- rownames(integrated_seurat@assays$RNA$data)[tolower(rownames(integrated_seurat@assays$RNA$data)) %in% 
                                                              tolower(features)]
    features <- sort(features)
    
    feature_plot <- function(i){
      split <- NULL
      Seurat::FeaturePlot(object=integrated_seurat,
                          slot="data",
                          features=i,
                          split.by=split,
                          #cols= c("grey", viridis(n=10, option="C", direction=-1)),
                          pt.size=0.4,
                          order=TRUE,
                          min.cutoff='q10',
                          reduction=paste0("umap.", base::tolower(reduc)),
                          label=TRUE,
                          combine=TRUE,
                          raster=FALSE) +
        scale_colour_gradientn(colours=rev(brewer.pal(n=11, name="RdBu"))[5:11])
    }
    
    # Plot UMAPs
    for (n in 1:ceiling((length(features)/50))){
      j <- 50*n-49
      k <- 50*n
      vec <- features[j:k]
      purrr::map(.x=vec[!is.na(vec)], 
                 .f=feature_plot) %>%
        cowplot::plot_grid(plotlist=.,
                           align="hv",
                           axis="tblr",
                           nrow=5,
                           ncol=10,
                           rel_widths=1,
                           rel_heights=1,
                           greedy=TRUE,
                           byrow=TRUE)
      
      ggplot2::ggsave(filename=paste0(proj, "_", colnames(df_final[col]), "_", n, ".jpg"),
                      plot=last_plot(),
                      device="jpg",
                      path="/hpc/home/kailasamms/plots/",
                      width=8.5*5,
                      height=11*2,
                      units=c("in"),
                      dpi=300,
                      limitsize=FALSE,
                      bg="white")
    }
  }
}

#####END#######################
































Epithelial_1 <- list("scRNASeq_BBN_C57B6" = c(2),
                     "scRNASeq_BBN_Rag" = c(4,9,12,15,17,18,28,7,10),
                     "scRNASeq_Chen" = c(4),
                     "scRNASeq_GSE217093" = c(0,1,3,16),
                     "scRNASeq_GSE222315" = c(14,15,16),
                     "scRNASeq_HRA003620" = c(4),
                     "scRNASeq_Jinfen" = c())
Epithelial_2 <- list("scRNASeq_BBN_C57B6" = c(),
                     "scRNASeq_BBN_Rag" = c(),
                     "scRNASeq_Chen" = c(0,1),
                     "scRNASeq_GSE217093" = c(),
                     "scRNASeq_GSE222315" = c(4,11),
                     "scRNASeq_HRA003620" = c(9,11,12,25),
                     "scRNASeq_Jinfen" = c())

Fibroblast_1 <- list("scRNASeq_BBN_C57B6" = c(1,9),
                     "scRNASeq_BBN_Rag" = c(0,25,27),
                     "scRNASeq_Chen" = c(8),
                     "scRNASeq_GSE217093" = c(14,25,27,33),
                     "scRNASeq_GSE222315" = c(7),
                     "scRNASeq_HRA003620" = c(3),
                     "scRNASeq_Jinfen" = c(13))
Fibroblast_2 <- list("scRNASeq_BBN_C57B6" = c(6),
                     "scRNASeq_BBN_Rag" = c(14),
                     "scRNASeq_Chen" = c(),
                     "scRNASeq_GSE217093" = c(29),
                     "scRNASeq_GSE222315" = c(),
                     "scRNASeq_HRA003620" = c(),
                     "scRNASeq_Jinfen" = c())
Lymphatic.Endothelial <- list("scRNASeq_BBN_C57B6" = c(18),
                              "scRNASeq_BBN_Rag" = c(30),
                              "scRNASeq_Chen" = c(),
                              "scRNASeq_GSE217093" = c(),
                              "scRNASeq_GSE222315" = c(),
                              "scRNASeq_HRA003620" = c(19),
                              "scRNASeq_Jinfen" = c())
Endothelial <- list("scRNASeq_BBN_C57B6" = c(7),
                    "scRNASeq_BBN_Rag" = c(5,24),
                    "scRNASeq_Chen" = c(5,6,20),
                    "scRNASeq_GSE217093" = c(8),
                    "scRNASeq_GSE222315" = c(5),
                    "scRNASeq_HRA003620" = c(15),
                    "scRNASeq_Jinfen" = c())
Granulocytes <- list("scRNASeq_BBN_C57B6" = c(),
                     "scRNASeq_BBN_Rag" = c(6),
                     "scRNASeq_Chen" = c(11),
                     "scRNASeq_GSE217093" = c(),
                     "scRNASeq_GSE222315" = c(8),
                     "scRNASeq_HRA003620" = c(13),
                     "scRNASeq_Jinfen" = c(3))
Plasma <- list("scRNASeq_BBN_C57B6" = c(23),
               "scRNASeq_BBN_Rag" = c(),
               "scRNASeq_Chen" = c(18,25,16,23),
               "scRNASeq_GSE217093" = c(),
               "scRNASeq_GSE222315" = c(19),
               "scRNASeq_HRA003620" = c(16),
               "scRNASeq_Jinfen" = c())
B <- list("scRNASeq_BBN_C57B6" = c(8),
          "scRNASeq_BBN_Rag" = c(),
          "scRNASeq_Chen" = c(12),
          "scRNASeq_GSE217093" = c(17),
          "scRNASeq_GSE222315" = c(3),
          "scRNASeq_HRA003620" = c(6),
          "scRNASeq_Jinfen" = c(14))
TCells <- list("scRNASeq_BBN_C57B6" = c(3,16),
               "scRNASeq_BBN_Rag" = c(),
               "scRNASeq_Chen" = c(2,3,17),
               "scRNASeq_GSE217093" = c(7),
               "scRNASeq_GSE222315" = c(0,1,9,13,17),
               "scRNASeq_HRA003620" = c(0,7,10,14,8,20),
               "scRNASeq_Jinfen" = c(1,5,9,4,15))
Macrophages <- list("scRNASeq_BBN_C57B6" = c(5,12),
                    "scRNASeq_BBN_Rag" = c(2,8,32),
                    "scRNASeq_Chen" = c(9),
                    "scRNASeq_GSE217093" = c(9),
                    "scRNASeq_GSE222315" = c(6),
                    "scRNASeq_HRA003620" = c(5),
                    "scRNASeq_Jinfen" = c(2,12))
DCs <- list("scRNASeq_BBN_C57B6" = c(17),
            "scRNASeq_BBN_Rag" = c(),
            "scRNASeq_Chen" = c(),
            "scRNASeq_GSE217093" = c(30),
            "scRNASeq_GSE222315" = c(),
            "scRNASeq_HRA003620" = c(18),
            "scRNASeq_Jinfen" = c(16))
Mast <- list("scRNASeq_BBN_C57B6" = c(0),
             "scRNASeq_BBN_Rag" = c(3),
             "scRNASeq_Chen" = c(15),
             "scRNASeq_GSE217093" = c(6,31),
             "scRNASeq_GSE222315" = c(10),
             "scRNASeq_HRA003620" = c(2),
             "scRNASeq_Jinfen" = c(0,6))
Myofibroblasts <- list("scRNASeq_BBN_C57B6" = c(14),
                       "scRNASeq_BBN_Rag" = c(19),
                       "scRNASeq_Chen" = c(7,27),
                       "scRNASeq_GSE217093" = c(),
                       "scRNASeq_GSE222315" = c(2),
                       "scRNASeq_HRA003620" = c(22),
                       "scRNASeq_Jinfen" = c())

# Create a dummy dataframe to store markers from standard approach
final_df <- data.frame("Gene"=unique(tolower(df_standard$gene)))

for (list in c("Epithelial_1", "Epithelial_2", "Fibroblast_1", "Fibroblast_2", 
               "Lymphatic.Endothelial", "Endothelial", "Granulocytes", "Plasma",
               "B", "TCells", "Macrophages", "DCs", "Mast", "Myofibroblasts")){
  
  # Create a list to store markers of each cell type
  final_g <- c()
  
  # Iterate through each dataset
  for (i in 1:length(get(list))){
    
    # Get unique markers for a cell type from each dataset
    genes <- df_standard %>% 
      dplyr::filter(cluster %in% paste0(names(get(list)[i]),"_",  get(list)[[i]])) %>%
      dplyr::select(gene) %>%
      unlist(use.names = FALSE) %>%
      tolower() %>%
      unique()
    
    # Merge all unique markers from each dataset
    final_g <- c(final_g, genes)
  }
  # Calculate number of datasets each marker was identified
  n_data <- sum(lengths(get(list)) > 0)  # number of datasets used for each cell type
  final_g <- data.frame(table(final_g))
  final_g <- final_g %>%
    dplyr::mutate(Percent = Freq*100/n_data)
  colnames(final_g) <- c(list, paste0(list, "_count"), paste0(list, "_percent"))
  
  # Remove markers identified in only a single dataset
  final_g <- final_g %>% 
    dplyr::filter(get(paste0(list, "_count")) > 1) %>%
    dplyr::mutate(Gene = get(list))
  
  # Merge all markers for each cell type present in multiple daasets to final dataframe
  final_df <- final_df %>%
    left_join(final_g,by=c("Gene"="Gene"))
}

final_df[is.na(final_df)] <- 0

# If a gene is detected as marker for 2 cell types, use it as marker for the 
# celltype where it has highest count
marker_list <- list()

for (celltype in c("Epithelial_1", "Epithelial_2", "Fibroblast_1", "Fibroblast_2", 
                   "Lymphatic.Endothelial", "Endothelial", "Granulocytes", "Plasma",
                   "B", "TCells", "Macrophages", "DCs", "Mast", "Myofibroblasts")){
  
  final_g <- final_df %>% 
    group_by(Gene) %>%
    dplyr::mutate(max_percent = max(Epithelial_1_percent, Epithelial_2_percent, 
                                    Fibroblast_1_percent, Fibroblast_2_percent, 
                                    Lymphatic.Endothelial_percent, Endothelial_percent,
                                    Granulocytes_percent, Plasma_percent,
                                    B_percent, TCells_percent, Macrophages_percent,
                                    DCs_percent, Mast_percent, Myofibroblasts_percent, na.rm=TRUE),
                  max_count = max(Epithelial_1_count, Epithelial_2_count, 
                                  Fibroblast_1_count, Fibroblast_2_count, 
                                  Lymphatic.Endothelial_count, Endothelial_count,
                                  Granulocytes_count, Plasma_count,
                                  B_count, TCells_count, Macrophages_count,
                                  DCs_count, Mast_count, Myofibroblasts_count, na.rm=TRUE)) %>%
    dplyr::filter(get(paste0(celltype, "_percent")) == max_percent) %>%
    dplyr::filter(get(paste0(celltype, "_percent")) >= 50) %>%
    dplyr::select(Gene, all_of(celltype), all_of(paste0(celltype, "_count")), max_percent) %>%
    dplyr::arrange(desc(max_percent)) %>%
    data.frame() %>%
    #dplyr::slice_head(n=50) %>%
    dplyr::select(Gene) %>%
    unlist(use.names=FALSE)
  
  final_g <- list(final_g)
  names(final_g) <- celltype
  marker_list <- c(marker_list, final_g)
}

max_l <- max(lengths(marker_list))
marker_list <- lapply(X=marker_list, FUN=function(x) { c(x, rep(x="", times=max_l-length(x)))})
df_final_standard <- data.frame(marker_list)
colnames(df_final_standard) <- names(marker_list)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb = wb, sheetName = "Standard_Markers")
openxlsx::writeData(wb = wb, sheet = "Standard_Markers", x = df_final_standard)
openxlsx::addWorksheet(wb = wb, sheetName = "Intermediate step")
openxlsx::writeData(wb = wb, sheet = "Intermediate step", x = final_df)
openxlsx::saveWorkbook(wb = wb, file = "Evaluated_Markers.xlsx", overwrite = TRUE)


####### We can further refine these markers

# We re-run Seurat::FindAllMarkers() using new groups we have listed above
# to get accurate pct1 and pct2 values and perform a final filtering of these
# markers

for (list in c("Epithelial_1", "Epithelial_2", "Fibroblast_1", "Fibroblast_2", 
               "Lymphatic.Endothelial", "Endothelial", "Granulocytes", "Plasma",
               "B", "TCells", "Macrophages", "DCs", "Mast", "Myofibroblasts")){
  
  
  
}

# df_custom <- df_custom %>%
#   dplyr::filter(avg_log2FC > 0.58, # custom expression cutoff 
#                 p_val_adj < 0.05,  # standard cutoff
#                 pct.1 > pct.2,     # standard cutoff
#                 pct.1 > 0.4) %>%   # custom percent cutoff
#   dplyr::group_by(cluster) %>%
#   dplyr::add_count(cluster) %>%
#   dplyr::filter(n >= 10) %>%
#   dplyr::arrange(desc(avg_log2FC)) # desc(ratio)
#   dplyr::slice_head(n=100) %>%
#   ungroup()

# ############ Create custom grouping of major clusters
# if (proj == "scRNASeq_BBN_C57B6"){
#   clusters <- list("cluster_1" = c(0,22),
#                    "cluster_2" = c(3,16),
#                    "cluster_3" = c(5,12),
#                    "cluster_4" = c(1,6,9),
#                    "cluster_5" = c(2,4,10,11,13),
#                    "cluster_6" = c(7),
#                    "cluster_7" = c(8),
#                    "cluster_8" = c(14),
#                    "cluster_9" = c(15),
#                    "cluster_10" = c(17),
#                    "cluster_11" = c(18),
#                    "cluster_12" = c(19),
#                    "cluster_13" = c(20),
#                    "cluster_14" = c(21),
#                    "cluster_15" = c(23),
#                    "cluster_16" = c(24))
# }
# 
# if (proj == "scRNASeq_BBN_Rag"){
#   clusters <- list("cluster_1" = c(0,4,7,8,9,10,12,13,15,17,18,20,26,27,31),
#                    "cluster_2" = c(3,16,22),
#                    "cluster_3" = c(5,24),
#                    "cluster_4" = c(1,11,23,25,28),
#                    "cluster_5" = c(2,6,21,29,32,33),
#                    "cluster_6" = c(14),
#                    "cluster_7" = c(19),
#                    "cluster_8" = c(30))
# }
# 
# # if (proj == "scRNASeq_Simon"){
# #   clusters <- list("cluster_1" = c(0,2,4,10,11,17,19),
# #                    "cluster_2" = c(3,5,27),
# #                    "cluster_3" = c(1,16,25,28),
# #                    "cluster_4" = c(7,13,15),
# #                    "cluster_5" = c(8,14,20),
# #                    "cluster_6" = c(6),
# #                    "cluster_7" = c(9),
# #                    "cluster_8" = c(12),
# #                    "cluster_9" = c(18),
# #                    "cluster_10" = c(21),
# #                    "cluster_11" = c(22),
# #                    "cluster_12" = c(23),
# #                    "cluster_13" = c(24),
# #                    "cluster_14" = c(26),
# #                    "cluster_15" = c(29))
# # }
# 
# if (proj == "scRNASeq_Jinfen"){
#   clusters <- list("cluster_1" = c(0,6),
#                    "cluster_2" = c(1,3,8),
#                    "cluster_3" = c(4),
#                    "cluster_4" = c(5),
#                    "cluster_5" = c(9),
#                    "cluster_6" = c(2),
#                    "cluster_7" = c(7),
#                    "cluster_8" = c(10),
#                    "cluster_9" = c(11),
#                    "cluster_10" = c(12),
#                    "cluster_11" = c(13),
#                    "cluster_12" = c(14),
#                    "cluster_13" = c(15),
#                    "cluster_14" = c(16))
# }
# 
# if (proj == "scRNASeq_Chen"){
#   clusters <- list("cluster_1" = c(0,1,4,12,13,22),
#                    "cluster_2" = c(2,3,17),
#                    "cluster_3" = c(9,15),
#                    "cluster_4" = c(16,28),
#                    "cluster_5" = c(5),
#                    "cluster_6" = c(6),                     
#                    "cluster_7" = c(7),
#                    "cluster_8" = c(8),
#                    "cluster_9" = c(10),
#                    "cluster_10" = c(11),
#                    "cluster_11" = c(14),
#                    "cluster_12" = c(18),
#                    "cluster_13" = c(19),
#                    "cluster_14" = c(20),
#                    "cluster_15" = c(21),
#                    "cluster_16" = c(23),
#                    "cluster_17" = c(24),
#                    "cluster_18" = c(25),
#                    "cluster_19" = c(26),
#                    "cluster_20" = c(27),
#                    "cluster_21" = c(29))
# }
# 
# integ <- annotate_data(res, reduc, celltype, clusters)
# 
# DefaultAssay(integ) <- "RNA"
# Idents(integ) <- "cell_type"
# 
# all_markers_custom <- Seurat::FindAllMarkers(object=integ,
#                                              assay = "RNA",
#                                              features = NULL,
#                                              logfc.threshold = 0.25,
#                                              test.use = "wilcox",
#                                              slot = "data",
#                                              min.pct = 0.1,
#                                              min.diff.pct = 0.1,
#                                              node = NULL,
#                                              verbose = TRUE,
#                                              only.pos = TRUE,
#                                              max.cells.per.ident = Inf,
#                                              random.seed = 1,
#                                              latent.vars = NULL,
#                                              min.cells.feature = 3,
#                                              min.cells.group = 1,
#                                              pseudocount.use = 1,
#                                              mean.fxn = NULL,
#                                              fc.name = NULL,
#                                              base = 2,
#                                              return.thresh = 0.01,
#                                              densify = FALSE)
# 
# all_markers_custom <- all_markers_custom %>% 
#   dplyr::mutate(pct.1 = dplyr::if_else(pct.1 == 0, 0.001, pct.1),
#                 pct.2 = dplyr::if_else(pct.2 == 0, 0.001, pct.2),
#                 ratio = pct.1/pct.2,
#                 cluster = paste0(proj, "_", cluster)) %>%
#   dplyr::relocate(cluster, gene, avg_log2FC, p_val_adj, pct.1, pct.2, ratio)






l1 <- list("scRNASeq_Jinfen" = c(),
           "scRNASeq_BBN_C57B6" = c(2),
           "scRNASeq_BBN_Rag" = c(8,9,12,15,17,18,26),
           "scRNASeq_Chen" = c(0,13,23,25,1,4,22),
           "scRNASeq_Simon" = c(0,13))

l2 <- list("scRNASeq_Jinfen" = c(10),
           "scRNASeq_BBN_C57B6" = c(4),
           "scRNASeq_BBN_Rag" = c(0),
           "scRNASeq_Chen" = c(),
           "scRNASeq_Simon" = c())

l3 <- list("scRNASeq_Jinfen" = c(),
           "scRNASeq_BBN_C57B6" = c(18),
           "scRNASeq_BBN_Rag" = c(30),
           "scRNASeq_Chen" = c(),
           "scRNASeq_Simon" = c())

l4 <- list("scRNASeq_Jinfen" = c(),
           "scRNASeq_BBN_C57B6" = c(7),
           "scRNASeq_BBN_Rag" = c(5,24),
           "scRNASeq_Chen" = c(5,6,21),
           "scRNASeq_Simon" = c(12))

l5 <- list("scRNASeq_Jinfen" = c(),
           "scRNASeq_BBN_C57B6" = c(),
           "scRNASeq_BBN_Rag" = c(),
           "scRNASeq_Chen" = c(14),
           "scRNASeq_Simon" = c(21))

l6 <- list("scRNASeq_Jinfen" = c(),
           "scRNASeq_BBN_C57B6" = c(23),
           "scRNASeq_BBN_Rag" = c(),
           "scRNASeq_Chen" = c(28),
           "scRNASeq_Simon" = c())

l7 <- list("scRNASeq_Jinfen" = c(3),
           "scRNASeq_BBN_C57B6" = c(),
           "scRNASeq_BBN_Rag" = c(),
           "scRNASeq_Chen" = c(12),
           "scRNASeq_Simon" = c())

l8 <- list("scRNASeq_Jinfen" = c(14),
           "scRNASeq_BBN_C57B6" = c(8),
           "scRNASeq_BBN_Rag" = c(),
           "scRNASeq_Chen" = c(11),
           "scRNASeq_Simon" = c())

l9 <- list("scRNASeq_Jinfen" = c(4,15),
           "scRNASeq_BBN_C57B6" = c(3),
           "scRNASeq_BBN_Rag" = c(),
           "scRNASeq_Chen" = c(2,3,17),
           "scRNASeq_Simon" = c(6))

l10 <- list("scRNASeq_Jinfen" = c(1,8,9),
            "scRNASeq_BBN_C57B6" = c(16),
            "scRNASeq_BBN_Rag" = c(),
            "scRNASeq_Chen" = c(),
            "scRNASeq_Simon" = c())

l11 <- list("scRNASeq_Jinfen" = c(2,12),
            "scRNASeq_BBN_C57B6" = c(5,12),
            "scRNASeq_BBN_Rag" = c(2,6),
            "scRNASeq_Chen" = c(9),
            "scRNASeq_Simon" = c())

l12 <- list("scRNASeq_Jinfen" = c(11,16),
            "scRNASeq_BBN_C57B6" = c(17),
            "scRNASeq_BBN_Rag" = c(),
            "scRNASeq_Chen" = c(),
            "scRNASeq_Simon" = c())

l13 <- list("scRNASeq_Jinfen" = c(0,6),
            "scRNASeq_BBN_C57B6" = c(0,22),
            "scRNASeq_BBN_Rag" = c(3),
            "scRNASeq_Chen" = c(15),
            "scRNASeq_Simon" = c())

l14 <- list("scRNASeq_Jinfen" = c(),
            "scRNASeq_BBN_C57B6" = c(14),
            "scRNASeq_BBN_Rag" = c(19),
            "scRNASeq_Chen" = c(7,26),
            "scRNASeq_Simon" = c(3,5,22))

l15 <- list("scRNASeq_Jinfen" = c(),
            "scRNASeq_BBN_C57B6" = c(6,9),
            "scRNASeq_BBN_Rag" = c(14),
            "scRNASeq_Chen" = c(),
            "scRNASeq_Simon" = c())

l16 <- list("scRNASeq_Jinfen" = c(13),
            "scRNASeq_BBN_C57B6" = c(1),
            "scRNASeq_BBN_Rag" = c(1),
            "scRNASeq_Chen" = c(8),
            "scRNASeq_Simon" = c())

l17 <- list("scRNASeq_Jinfen" = c(),
            "scRNASeq_BBN_C57B6" = c(),
            "scRNASeq_BBN_Rag" = c(28),
            "scRNASeq_Chen" = c(),
            "scRNASeq_Simon" = c(8,14))

# Create a dummy list to store markers from standard approach
marker_list <- list()
final_g <- tolower(unique(df_standard$gene))

for (list in c("l1", "l2", "l3", "l4", "l5", "l6", "l7", "l8", "l9", "l10", 
               "l11", "l12", "l13", "l14", "l15", "l16", "l17")){
  
  final_g <- tolower(unique(df_standard$gene))
  for (i in 1:length(get(list))){
    
    genes <- df_standard %>% 
      dplyr::filter(cluster %in% paste0(names(get(list)[i]),"_",  get(list)[[i]])) %>%
      dplyr::select(gene) %>%
      unlist(use.names = FALSE) %>%
      tolower()
    
    if(length(genes) > 0){
      final_g <- intersect(final_g, genes)
    }
  }
  marker_list <- c(marker_list, list(final_g))
}

max_l <- max(lengths(marker_list))
marker_list <- lapply(X=marker_list, FUN=function(x) { c(x, rep(x="", times=max_l-length(x)))})
df_final_standard <- data.frame(marker_list)

colnames(df_final_standard) <- c("Epithelial_1", "Epithelial_2", "Lymphatic.Endothelial", 
                                 "Endothelial", "Granulocytes", "Plasma Cells", "Mitotic Cells", "B.Cells", "T.Cells",
                                 "NK cells", "Macrophages", "Dendritic.Cells", "Mast Cells", 
                                 "Myofibroblasts", "Fibroblasts_1", "Fibroblasts_2", 
                                 "Fibroblasts_3")

# # For custom approach
# l1 <- list("scRNASeq_Jinfen" = c(),
#            "scRNASeq_BBN_C57B6" = c(8),
#            "scRNASeq_BBN_Rag" = c(7),
#            "scRNASeq_Chen" = c(7,19),
#            "scRNASeq_Simon" = c(2,5,11))
# 
# l2 <- list("scRNASeq_Jinfen" = c(11),
#            "scRNASeq_BBN_C57B6" = c(4,12),
#            "scRNASeq_BBN_Rag" = c(4,6),
#            "scRNASeq_Chen" = c(8),
#            "scRNASeq_Simon" = c(9))
# 
# l3 <- list("scRNASeq_Jinfen" = c(),
#            "scRNASeq_BBN_C57B6" = c(11),
#            "scRNASeq_BBN_Rag" = c(8),
#            "scRNASeq_Chen" = c(),
#            "scRNASeq_Simon" = c())
# 
# l4 <- list("scRNASeq_Jinfen" = c(),
#            "scRNASeq_BBN_C57B6" = c(6),
#            "scRNASeq_BBN_Rag" = c(3),
#            "scRNASeq_Chen" = c(5,6,15),
#            "scRNASeq_Simon" = c(8))
# 
# l5 <- list("scRNASeq_Jinfen" = c(6,10),
#            "scRNASeq_BBN_C57B6" = c(3),
#            "scRNASeq_BBN_Rag" = c(5),
#            "scRNASeq_Chen" = c(3),
#            "scRNASeq_Simon" = c(7))
# 
# l6 <- list("scRNASeq_Jinfen" = c(),
#            "scRNASeq_BBN_C57B6" = c(5),
#            "scRNASeq_BBN_Rag" = c(1),
#            "scRNASeq_Chen" = c(1,16,18),
#            "scRNASeq_Simon" = c(1,4))
# 
# l7 <- list("scRNASeq_Jinfen" = c(1),
#            "scRNASeq_BBN_C57B6" = c(1),
#            "scRNASeq_BBN_Rag" = c(2),
#            "scRNASeq_Chen" = c(),
#            "scRNASeq_Simon" = c())
# 
# l8 <- list("scRNASeq_Jinfen" = c(),
#            "scRNASeq_BBN_C57B6" = c(),
#            "scRNASeq_BBN_Rag" = c(),
#            "scRNASeq_Chen" = c(11),
#            "scRNASeq_Simon" = c(10))
# 
# l9 <- list("scRNASeq_Jinfen" = c(12),
#             "scRNASeq_BBN_C57B6" = c(7),
#             "scRNASeq_BBN_Rag" = c(),
#             "scRNASeq_Chen" = c(10),
#             "scRNASeq_Simon" = c())
# 
# l10 <- list("scRNASeq_Jinfen" = c(3,13),
#             "scRNASeq_BBN_C57B6" = c(2),
#             "scRNASeq_BBN_Rag" = c(2,17),
#             "scRNASeq_Chen" = c(6),
#             "scRNASeq_Simon" = c())
# 
# l11 <- list("scRNASeq_Jinfen" = c(9,14),
#             "scRNASeq_BBN_C57B6" = c(10),
#             "scRNASeq_BBN_Rag" = c(),
#             "scRNASeq_Chen" = c(),
#             "scRNASeq_Simon" = c())
# 
# # Create a dummy list to store markers from custom approach
# marker_list <- list()
# final_g <- tolower(unique(df_custom$gene))
# 
# for (list in c("l1", "l2", "l3", "l4", "l5", "l6", "l7", "l8", "l9", "l10", 
#                "l11")){
#   
#   final_g <- tolower(unique(df_custom$gene))
#   for (i in 1:length(get(list))){
#     
#     genes <- df_custom %>% 
#       dplyr::filter(cluster %in% paste0(names(get(list)[i]),"_cluster_",  get(list)[[i]])) %>%
#       dplyr::select(gene) %>%
#       unlist(use.names = FALSE) %>%
#       tolower()
#     
#     if(length(genes) > 0){
#       final_g <- intersect(final_g, genes)
#     }
#   }
#   marker_list <- c(marker_list, list(final_g))
# }
# 
# max_l <- max(lengths(marker_list))
# marker_list <- lapply(X=marker_list, FUN=function(x) { c(x, rep(x="", times=max_l-length(x)))})
# df_final_custom <- data.frame(marker_list)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb = wb, sheetName = "Standard_Markers")
openxlsx::writeData(wb = wb, sheet = "Standard_Markers", x = df_final_standard)
# openxlsx::addWorksheet(wb = wb, sheetName = "Custom_Markers")
# openxlsx::writeData(wb = wb, sheet = "Custom_Markers", x = df_final_custom)
openxlsx::saveWorkbook(wb = wb, file = "Evaluated_Markers.xlsx", overwrite = TRUE)

#******************************************************************************#
#      STEP 5: VALIDATE THE SPECIFICITY FOR EACH MARKER GENE USING UMAPS       #
#******************************************************************************#

for (proj in c("scRNASeq_BBN_C57B6", "scRNASeq_BBN_Rag", "scRNASeq_Jinfen", 
               "scRNASeq_Chen", "scRNASeq_Simon")){
  
  seurat_results <- paste0("/hpc/home/kailasamms/scratch/", proj, "/results_seurat/")
  res <- 0.4
  reduc <- "rpca" #"harmony"
  celltype <- NULL
  #genes <- read.xlsx("/hpc/home/kailasamms/Markers.xlsx")
  genes <- read.xlsx("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Markers.xlsx")
  
  # Load the integrated seurat object
  integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn", 
                                      dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  # Change active.ident to major clusters we annotated
  DefaultAssay(integrated_seurat) <- "RNA"
  idents <- paste0("cluster.", res, ".", base::tolower(reduc))
  Idents(object = integrated_seurat) <- idents
  
  feature_plot <- function(i){
    
    split <- "Condition" #NULL
    split <- NULL
    Seurat::FeaturePlot(object = integrated_seurat,
                        slot = "data",
                        features = i,
                        split.by = split,
                        cols =  c("#D1E5F0", "#F7F7F7", "#67001F"),
                        pt.size = 0.4,
                        order = TRUE,
                        min.cutoff = 'q10',
                        reduction = paste0("umap.", base::tolower(reduc)),
                        label = TRUE,
                        combine = TRUE,
                        raster = FALSE) +
      scale_colour_gradientn(colours = c("#D1E5F0", "#F7F7F7", "#67001F"), limits=c(0,2))
    # rev(brewer.pal(n = 11, name = "RdBu"))[5:11]
  } 
  
  for (i in 1:ncol(genes)){
    features <- genes[,i] %>% unlist(use.names=FALSE)
    features <- sort(rownames(integrated_seurat@assays$RNA$data)[tolower(rownames(integrated_seurat@assays$RNA$data)) %in% tolower(features)])
    
    # Plot UMAPs
    for (n in 1:ceiling((length(features)/12))){
      j <- 12*n-11
      k <- 12*n
      vec <- features[j:k]
      purrr::map(.x = vec[!is.na(vec)], 
                 .f = feature_plot) %>%
        cowplot::plot_grid(plotlist = .,
                           align = "hv",
                           axis = "tblr",
                           nrow = 3,
                           ncol = 4,
                           rel_widths = 1,
                           rel_heights = 1,
                           greedy = TRUE,
                           byrow = TRUE)
      
      ggplot2::ggsave(filename = paste0(proj,"_", colnames(genes)[i], "_", n, ".jpg"),
                      plot = last_plot(),
                      device = "jpeg",
                      scale = 1,
                      width = dplyr::if_else(is.null(split), 9*3, 9*6),
                      height = 11,
                      units = c("in"),
                      dpi = 300,
                      limitsize = FALSE,
                      bg = "white")
    }
  }
}


#******************************************************************************#

# # Find out clusters that are similar and genes that are common
# combinations_standard <- utils::combn(x=unique(df_standard$cluster), m=2)
# combinations_custom <- utils::combn(x=unique(df_custom$cluster), m=2)
# final_standard <- list()
# final_custom <- list()
# 
# for (i in 1:ncol(combinations_standard)){
#   
#   l <- intersect(df_standard %>% 
#                    dplyr::filter(cluster == combinations_standard[1,i]) %>%
#                    ungroup() %>%
#                    dplyr::select(gene) %>%
#                    unlist(use.names=FALSE) %>%
#                    tolower(),
#                  df_standard %>% 
#                    dplyr::filter(cluster == combinations_standard[2,i]) %>%
#                    ungroup() %>%
#                    dplyr::select(gene) %>%
#                    unlist(use.names=FALSE) %>%
#                    tolower())
#   
#   # Custom cutoff of 20 genes to be similar between 2 clusters
#   if (length(l) > 20){
#     l <-  list(l)
#     names(l) <- paste0(combinations_standard[1,i],".", combinations_standard[2,i])
#     final_standard <- c(final_standard, l)
#   }
# }
# 
# for (i in 1:ncol(combinations_custom)){
#   
#   l <- intersect(df_custom %>% 
#                    dplyr::filter(cluster == combinations_custom[1,i]) %>%
#                    ungroup() %>%
#                    dplyr::select(gene) %>%
#                    unlist(use.names=FALSE) %>%
#                    tolower(),
#                  df_custom %>% 
#                    dplyr::filter(cluster == combinations_custom[2,i]) %>%
#                    ungroup() %>%
#                    dplyr::select(gene) %>%
#                    unlist(use.names=FALSE) %>%
#                    tolower())
#   
#   # Custom cutoff of 20 genes to be similar between 2 clusters
#   if (length(l) > 20){
#     l <-  list(l)
#     names(l) <- paste0(combinations_custom[1,i],".", combinations_custom[2,i])
#     final_custom <- c(final_custom, l)
#   }
# }
# 
# # Convert list to dataframe
# for (i in 1:length(final_standard)){
#   n <- names(final_standard[i])
#   final_standard[i] <- list(c(unlist(final_standard[i], use.names=FALSE), rep(x="", times=max(lengths(final_standard))-length(final_standard[[i]]))))
#   names(final_standard[i]) <- names(final_standard[i])
# }
# 
# for (i in 1:length(final_custom)){
#   n <- names(final_custom[i])
#   final_custom[i] <- list(c(unlist(final_custom[i], use.names=FALSE), rep(x="", times=max(lengths(final_custom))-length(final_custom[[i]]))))
#   names(final_custom[i]) <- names(final_custom[i])
# }
# 
# common_markers_standard <- data.frame(final_standard)
# common_markers_custom <- data.frame(final_custom)
# 
# # Save all the markers and correlation
# wb <- openxlsx::createWorkbook()
# openxlsx::addWorksheet(wb = wb, sheetName = "Standard_Markers")
# openxlsx::writeData(wb = wb, sheet = "Standard_Markers", x = common_markers_standard )
# openxlsx::addWorksheet(wb = wb, sheetName = "Custom_Markers")
# openxlsx::writeData(wb = wb, sheet = "Custom_Markers", x = common_markers_custom)
# openxlsx::saveWorkbook(wb = wb, file = "Compiled_Markers.xlsx", overwrite = TRUE)