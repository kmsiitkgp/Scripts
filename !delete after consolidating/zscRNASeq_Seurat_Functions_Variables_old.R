# NOTE: Perform cluster identification using conserved markers. If you cannot
# identify any cluster using conserved markers, then use the markers from
# FindAllMarkers() to identify the unidentified clusters.
# FindConservedMarkers() is more accurate as it gives genes conserved across
# multiple samples.
plot_conserved_modules <- function(res, reduc, celltype, sheetname){
  
  # Load the integrated seurat object
  integrated_seurat <- base::readRDS(paste0(seurat_results, "integrated_seurat_snn", 
                                            dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  # Set default assay
  DefaultAssay(integrated_seurat) <- "RNA"
  
  # Set identity to an existing column in meta data
  idents <- paste0("cluster.", res, ".", base::tolower(reduc))
  Idents(object=integrated_seurat) <- idents
  
  # Read markers
  marker_df <- openxlsx::read.xlsx(xlsxFile=paste0(scripts_path, "scRNASeq_Markers.xlsx"),
                                   sheet=sheetname)
  
  module_score_seurat <- function(celltype){
    
    Seurat::FeaturePlot(object=integrated_seurat,
                        slot="data",
                        features=paste0(make.names(celltype),1),
                        #cols= c("grey", viridis(n=10, option="C", direction=1)),
                        pt.size=0.4,
                        order=TRUE,
                        min.cutoff='q10',
                        reduction=paste0("umap.", base::tolower(reduc)),
                        label=TRUE,
                        combine=TRUE,
                        raster=FALSE) +  
      #BUG: if raster=TRUE, order=TRUE is ignored. So, set raster=FALSE
      scale_colour_gradientn(colours=rev(brewer.pal(n=11, name="RdBu"))[5:11])
  }
  
  module_score_ucell <- function(celltype){ 
    
    Seurat::FeaturePlot(object=integrated_seurat,
                        slot="data",
                        features=paste0(make.names(celltype), "_UCell"),
                        #cols= c("grey", viridis(n=10, option="C", direction=1)),
                        pt.size=0.4,
                        order=TRUE,
                        min.cutoff='q10',
                        reduction=paste0("umap.", base::tolower(reduc)),
                        label=TRUE,
                        combine=TRUE,
                        raster=FALSE) +  
      #BUG: if raster=TRUE, order=TRUE is ignored. So, set raster=FALSE
      scale_colour_gradientn(colours=rev(brewer.pal(n=11, name="RdBu"))[5:11])
  }
  
  module_score_ucell_violin <- function(celltype){
    
    Idents(integrated_seurat) <- idents
    Seurat::VlnPlot(object=integrated_seurat,
                    features=paste0(make.names(celltype), "_UCell"),
                    assay=NULL,
                    pt.size=0,
                    sort=TRUE,
                    combine=TRUE,
                    raster=FALSE)
  }
  
  funcs <- c("module_score_seurat", "module_score_ucell", "module_score_ucell_violin") 
  
  for (j in 1:length(funcs)){
    
    purrr::map(.x=colnames(marker_df), 
               .f=get(funcs[j])) %>%
      cowplot::plot_grid(plotlist=.,
                         align="hv",
                         axis="tblr",
                         nrow=NULL,
                         ncol=dplyr::if_else(ncol(marker_df) > 10, 5, ceiling(ncol(marker_df)/3)),
                         rel_widths=1,
                         rel_heights=1,
                         greedy=TRUE,
                         byrow=TRUE)
    
    filename <- dplyr::case_when(j == 1 ~ paste0("Module_plot(", sheetname, ")_", celltype, "_", reduc, ".jpg"),
                                 j == 2 ~ paste0("Module_plot(", sheetname, ")_", celltype, "_", reduc, "_UCell.jpg"),
                                 j == 3 ~ paste0("Module_plot(", sheetname, ")_", celltype, "_", reduc, "_UCell_violin.jpg"))
    
    ggplot2::ggsave(filename=filename,
                    plot=last_plot(),
                    device="jpeg",
                    path=diagnostics_path,
                    width=8.5*4,
                    height=11*2,
                    units=c("in"),
                    dpi=300,
                    limitsize=FALSE,
                    bg="white")
  }
}

#******************************************************************************#
#                   ANNOTATE CLUSTERS AND SAVE SEURAT OBJECT                   #
#******************************************************************************#

# Annotate cells based on UMAP of module scores from plot_conserved_modules()
# Here we annotate all cells belonging to a cluster to a 'single' cell type.
# while in reality there could be some contaminating cells (eg: Myeloid cells in
# epithelial cluster etc). We next have to remove these contaminating cells
# (before performing analysis on subtypes) resulting in loss of cells. 
annotate_data_umap <- function(res, reduc, celltype, clusters){
  
  # Load the integrated seurat object
  integrated_seurat <- base::readRDS(paste0(seurat_results, "integrated_seurat_snn", 
                                            dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  idents <- paste0("cluster.", res, ".", base::tolower(reduc))
  
  # Make sure you have assigned all clusters to one of the cell types
  # NOTE: "integrated_snn_res.1.4" etc are stored as factors. 
  # So, use as.character() and then as.numeric() to get accurate cluster values
  list_1 <- integrated_seurat@meta.data %>% 
    dplyr::count(get(idents)) %>% 
    dplyr::select(identity(1)) %>% 
    unlist(use.names=FALSE) %>% 
    as.character() %>% 
    as.numeric() %>% 
    sort()
  
  list_2 <- clusters %>% 
    unlist(., use.names=FALSE) %>% 
    sort()
  
  # Proceed with annotation ONLY if all clusters have been renamed
  if (identical(list_1, list_2)){
    print("All Clusters have been annotated")
    
    # Extract metadata from Seurat object, assign appropriate resolution to
    # seurat_clusters column and add cell_type, sub_type, cell_class columns
    data <- integrated_seurat@meta.data %>% 
      dplyr::mutate(seurat_clusters=get(idents),
                    cell_type=NA, sub_type=NA)
    
    # Assign cell type based on cluster numbers within seurat_clusters column
    for (j in 1:nrow(data)){
      for (i in 1:length(clusters)){
        if (as.numeric(as.character(data$seurat_clusters[j])) %in% clusters[[i]]){
          data$cell_type[j] <- names(clusters[i])
        }
      }
    }
    
    # Check summary of cell counts
    print(data %>% dplyr::count(cell_type) %>% dplyr::arrange(n))
    cat("\n")
    
    # Import the metadata into Seurat object
    integrated_seurat@meta.data <- data
  } else {
    cat("\nYou missed annotating these clusters:\t", setdiff(list_1, list_2))
    cat("\nThese clusters are not present in data:\t", setdiff(list_2, list_1))
    cat("\nThese clusters have duplicate annotation:\t", list_2[duplicated(list_2)])
  }
  
  return(integrated_seurat)
}

# Annotate cells based on scores calculated by 
# (i) UCell::AddModuleScore_UCell() : scores lie between [0,1]
# (ii) Seurat::AddModuleScore()     : scores can be positive or negative
# Cells in a cluster may belong to 'multiple' cell types. We can retain these 
# contaminating cells while performing subtype analysis but remove them from 
# the final UMAP plot of all cell types. This way we retain most cells for 
# subtype analysis but also identify these contaminants and remove them before
# final visualization.

# NOTE: DO NOT CHANGE column names in scRNASeq_Markers.xlsx as the varaibles
# defined within the function are based on column names in scRNASeq_Markers.xlsx
annotate_data_score <- function(integrated_seurat, celltype){
  
  integrated_seurat@meta.data <- integrated_seurat@meta.data %>%
    # you can also use rowwise() instead of using group_by() and ungroup()
    group_by(Cell) %>%  # only 1 row per cell after grouping
    dplyr::mutate(ucell_class = max(B_UCell, Dendritic_UCell, 
                                    Endothelial_UCell, Epithelial_UCell,
                                    Fibroblasts_UCell, Granulocytes_UCell,
                                    Lymphatic.Endothelial_UCell, Macrophages_UCell,
                                    Mast_UCell, Myofibroblasts_UCell, NK_UCell,
                                    Plasma_UCell, T_UCell, Neurons_UCell, 
                                    Erythrocytes_UCell)) %>%
    dplyr::mutate(seurat_class = max(B1, Dendritic1, Endothelial1, Epithelial1,
                                     Fibroblasts1, Granulocytes1, 
                                     Lymphatic.Endothelial1, Macrophages1,
                                     Mast1, Myofibroblasts1, NK1, Plasma1, T1, 
                                     Neurons1, Erythrocytes1)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(ucell_class = dplyr::case_when(Epithelial_UCell           == ucell_class & ucell_class > 0 ~ "Epithelial",
                                                 Fibroblasts_UCell           == ucell_class & ucell_class > 0 ~ "Fibroblasts",
                                                 Myofibroblasts_UCell        == ucell_class & ucell_class > 0 ~ "Myofibroblasts",
                                                 Dendritic_UCell             == ucell_class & ucell_class > 0 ~ "Myeloid - DCs",
                                                 Granulocytes_UCell          == ucell_class & ucell_class > 0 ~ "Myeloid - Granulocytes",
                                                 Macrophages_UCell           == ucell_class & ucell_class > 0 ~ "Myeloid - Macrophages", 
                                                 Mast_UCell                  == ucell_class & ucell_class > 0 ~ "Myeloid - Mast", 
                                                 B_UCell                     == ucell_class & ucell_class > 0 ~ "Lymphoid - B", 
                                                 Plasma_UCell                == ucell_class & ucell_class > 0 ~ "Lymphoid - Plasma",
                                                 T_UCell                     == ucell_class & ucell_class > 0 ~ "Lymphoid - T", 
                                                 NK_UCell                    == ucell_class & ucell_class > 0 ~ "Lymphoid - NK",
                                                 Endothelial_UCell           == ucell_class & ucell_class > 0 ~ "Endothelial",
                                                 Lymphatic.Endothelial_UCell == ucell_class & ucell_class > 0 ~ "Endothelial - Lymphatic",
                                                 Neurons_UCell               == ucell_class & ucell_class > 0 ~ "Neurons",
                                                 Erythrocytes_UCell          == ucell_class & ucell_class > 0 ~ "Erythrocytes",
                                                 TRUE ~ "Unclassified")) %>%
    dplyr::mutate(seurat_class = dplyr::case_when(Epithelial1           == seurat_class & seurat_class > 0 ~ "Epithelial",
                                                  Fibroblasts1           == seurat_class & seurat_class > 0 ~ "Fibroblasts",
                                                  Myofibroblasts1        == seurat_class & seurat_class > 0 ~ "Myofibroblasts",
                                                  Dendritic1             == seurat_class & seurat_class > 0 ~ "Myeloid - DCs",
                                                  Granulocytes1          == seurat_class & seurat_class > 0 ~ "Myeloid - Granulocytes",
                                                  Macrophages1           == seurat_class & seurat_class > 0 ~ "Myeloid - Macrophages", 
                                                  Mast1                  == seurat_class & seurat_class > 0 ~ "Myeloid - Mast", 
                                                  B1                     == seurat_class & seurat_class > 0 ~ "Lymphoid - B", 
                                                  Plasma1                == seurat_class & seurat_class > 0 ~ "Lymphoid - Plasma",
                                                  T1                     == seurat_class & seurat_class > 0 ~ "Lymphoid - T", 
                                                  NK1                    == seurat_class & seurat_class > 0 ~ "Lymphoid - NK",
                                                  Endothelial1           == seurat_class & seurat_class > 0 ~ "Endothelial",
                                                  Lymphatic.Endothelial1 == seurat_class & seurat_class > 0 ~ "Endothelial - Lymphatic",
                                                  Neurons1               == seurat_class & seurat_class > 0 ~ "Neurons",
                                                  Erythrocytes1          == seurat_class & seurat_class > 0 ~ "Erythrocytes",
                                                  TRUE ~ "Unclassified")) %>%
    dplyr::mutate(rnames = Cell) %>%
    tibble::column_to_rownames("rnames")
  
  return(integrated_seurat)
}

#******************************************************************************#
#                   REMOVE CLUSTERS WITH MULTIPLE CELLTYPES                    #
#******************************************************************************#

# NOTE: Mixed is a list of cluster numbers of mixed clusters for each cell type.
# Mixed <- list("Epithelial"=c(25,26), "Fibroblasts"=c(13,14)}

remove_mixed_clusters <- function(res, reduc, celltype, Mixed){
  
  # Load the integrated seurat object
  integrated_seurat <- base::readRDS(paste0(seurat_results, "integrated_seurat_snn", 
                                            dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  # Make a copy of original seurat object before subsetting
  saveRDS(integrated_seurat, paste0(seurat_results, "integrated_seurat_snn", 
                                    dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, "_original.rds"))))
  
  # Set identity to an existing column in meta data
  idents <- paste0("cluster.", res, ".", base::tolower(reduc))
  Idents(integrated_seurat) <- idents
  
  cat("\n", celltype, " present initially:", nrow(integrated_seurat@meta.data))
  
  # Subset out the Mixed clusters
  integrated_seurat <- subset(x=integrated_seurat,
                              !!rlang::sym(idents) %in% Mixed[[celltype]],
                              invert=TRUE)
  
  cat("\n", celltype, " present finally:", nrow(integrated_seurat@meta.data))
  
  saveRDS(integrated_seurat, paste0(seurat_results, "integrated_seurat_snn", 
                                    dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
}

#******************************************************************************#
#                SUBSET SPECIFIC CELLTYPE FOR SUBTYPE ANALYSIS                 #
#******************************************************************************#

# prep_data() filters celltypes based on cell_class column of metadata which 
# is based on UCell scoring. To filter celltypes based on cell_type column
# of metadata which is based on UMAP classification, change the function code

prep_data <- function(integrated_seurat, celltype){
  
  # Get major cell types present in seurat object
  # major_celltypes <- unique(integrated_seurat@meta.data$cell_type)
  #major_celltypes <- unique(integrated_seurat@meta.data$ucell_class)
  major_celltypes <- unique(integrated_seurat@meta.data$seurat_class)
  
  # Identify all sub types of relevance. 
  # Since, we want to subset all Myeloid subtypes like "Myeloid-MDSC", 
  # "Myeloid-Macrophages", set celltype == "Myeloid"
  celltypes_of_interest <- major_celltypes[grepl(pattern=celltype, x=major_celltypes, ignore.case=TRUE)]
  
  # Keep ONLY necessary celltype being analyzed
  filtered_seurat <- subset(x=integrated_seurat, 
                            seurat_class %in% celltypes_of_interest)
  
  # Print cell numbers to double check
  print(filtered_seurat@meta.data %>% dplyr::count(Condition, seurat_class, cell_type, sub_type, Sample))
  
  # Remove unwanted assays after changing default assay
  DefaultAssay(filtered_seurat) <- "RNA"
  filtered_seurat[["SCT"]] <- NULL
  filtered_seurat@graphs <- list()
  filtered_seurat@reductions <- list()
  
  # Remove samples with less than 50 cells so PCA/integration dont throw errors
  samples_with_few_cells <- filtered_seurat@meta.data %>% 
    dplyr::count(Sample) %>% 
    dplyr::filter(n<50) %>% 
    dplyr::select(Sample) %>% 
    unlist(use.names=FALSE)
  
  filtered_seurat <- subset(x=filtered_seurat,
                            Sample %in% samples_with_few_cells,
                            invert=TRUE)
  
  return(filtered_seurat)
}

#******************************************************************************#
#                      STEP 12: IDENTIFY SIMILAR CLUSTERS                      # 
#******************************************************************************#

identify_lineage <- function(celltype){
  
  cat("\n**********************", celltype, "Analysis**********************\n")
  
  # Load the integrated seurat object
  integrated_seurat <- base::readRDS(paste0(seurat_results, "integrated_seurat_snn", 
                                            dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  # Create matrix similar to read data for each cluster to identify closely related clusters
  cor_meta_data <- integrated_seurat@meta.data %>%
    dplyr::distinct(integrated_snn_res.1.4, .keep_all=TRUE) %>%
    tibble::remove_rownames(.) %>%
    tibble::column_to_rownames("integrated_snn_res.1.4")
  
  clusters <- as.vector(unlist(integrated_seurat@meta.data %>% dplyr::distinct(integrated_snn_res.1.4)))
  cor_read_data <- matrix(NA, nrow=nrow(integrated_seurat@assays$RNA@counts), ncol=length(clusters))
  rownames(cor_read_data) <- rownames(integrated_seurat@assays$RNA@counts)
  colnames(cor_read_data) <- clusters
  
  for(i in 1:length(clusters)){
    
    # Create a list of cells for each cluster
    cor_cells_subset <- rownames(integrated_seurat@meta.data %>% dplyr::filter(integrated_snn_res.1.4 == as.numeric(clusters[i])))
    
    # Get the raw counts of these cells. counts are stored in sparse matrix format
    # So, use data.frame to convert "." in sparse matrix to "0"
    cor_subset <- data.frame(integrated_seurat@assays$RNA@counts[,cor_cells_subset])
    cor_read_data[,i]  <- rowSums(cor_subset)
  }
  
  # Remove unwanted genes
  unwanted_genes <- rownames(cor_read_data)[grep(pattern="^RP[SL]|RIK$|^MT-|^GM[0-9.]+$", x= rownames(cor_read_data), ignore.case=TRUE)]
  keep_genes <- setdiff(rownames(cor_read_data), unwanted_genes)  
  cor_read_data <- cor_read_data[keep_genes,]
  
  # Create DESeq2 object so we can use DESeq2::rlog() or DESeq2::vst() on it
  cor_dds <- DESeqDataSetFromMatrix(countData=cor_read_data,
                                    colData=cor_meta_data,
                                    design=~1)
  
  #Determine number of cuts in heatmap
  cuts <- ceiling(length(clusters)/3)
  
  # Use vst as well as rld
  vst <- DESeq2::vst(cor_dds)
  head(assay(vst))
  sampleDists <- dist(t(assay(vst)))
  sampleDistMatrix <- as.matrix(sampleDists)
  
  pheatmap::pheatmap(mat=sampleDistMatrix,
                     scale="none",
                     cutree_rows=cuts,
                     cutree_cols=cuts,
                     cluster_rows=TRUE,   #cluster the rows
                     cluster_cols=TRUE,   #cluster the columns
                     fontsize=8, 
                     fontsize_row=8, 
                     fontsize_col=8,
                     angle_col=c("270", "0", "45", "90", "315"),
                     fontsize_number=0.8*fontsize,
                     filename=paste0(diagnostics_path, "Cluster_Correlation_vst_", celltype, ".jpg"))
  
  rlog <- DESeq2::rlog(cor_dds)
  head(assay(rlog))
  sampleDists <- dist(t(assay(rlog)))
  sampleDistMatrix <- as.matrix(sampleDists)
  
  pheatmap::pheatmap(mat=sampleDistMatrix,
                     scale="none",
                     cutree_rows=cuts,
                     cutree_cols=cuts,
                     cluster_rows=TRUE,   #cluster the rows
                     cluster_cols=TRUE,   #cluster the columns
                     fontsize=8, 
                     fontsize_row=8, 
                     fontsize_col=8,
                     angle_col=c("270", "0", "45", "90", "315"),
                     fontsize_number=0.8*fontsize,
                     filename=paste0(diagnostics_path, "Cluster_Correlation_rlog_", celltype, ".jpg"))
}

#******************************************************************************#
#       RE-ANNOTATE ORIGINAL SEURAT OBJECT WITH SUBTYPE, CELLCLASS INFO        #
#******************************************************************************#

re_annotate <- function(celltypes){
  
  celltype <- celltypes[1]
  # Load the integrated seurat object
  integrated_seurat <- base::readRDS(paste0(seurat_results, "integrated_seurat_snn",
                                            dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  merged_metadata <- integrated_seurat@meta.data
  
  
  for (celltype in celltypes[-1]){
    
    # Load the integrated seurat object
    integrated_seurat <- base::readRDS(paste0(seurat_results, "integrated_seurat_snn",
                                              dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
    
    merged_metadata <- dplyr::bind_rows(merged_metadata, integrated_seurat@meta.data)
  }
  
  # Load the full integrated seurat object
  integrated_seurat <- base::readRDS(paste0(seurat_results, "integrated_seurat_snn.rds"))
  
  metadata <- integrated_seurat@meta.data
  
  # You can see there are more cells in metadata than merged_metadata.   
  # Some of these cells were likely removed during subtype analysis while others
  # cells belong to cell_type like Muscle, Neurons, Unclassified. So, we need to
  # generate sub_type and cell_class info for these cells.
  missing_metadata <- dplyr::anti_join(x=metadata,
                                       y=merged_metadata, 
                                       by=c("Cell"="Cell")) %>%
    dplyr::mutate(cell_class=dplyr::case_when(cell_type %in% base::setdiff(x=unique(metadata$cell_type), y=unique(merged_metadata$cell_type)) ~ cell_type,
                                              cell_type == "Unclassified" ~ "Unclassified",
                                              TRUE ~ "Mixed"),
                  sub_type=dplyr::case_when(cell_type %in% base::setdiff(x=unique(metadata$cell_type), y=unique(merged_metadata$cell_type)) ~ cell_type,
                                            cell_type == "Unclassified" ~ "Unclassified",
                                            TRUE ~ "Mixed"))
  
  # Import cell_class, sub_type and cell_type from merged_metadata
  subtyped_metadata <-  dplyr::inner_join(x=metadata %>% dplyr::select(everything(), -cell_type, -sub_type, -cell_class),
                                          y=merged_metadata %>% dplyr::select(Cell, cell_type, sub_type, cell_class), 
                                          by=c("Cell"="Cell"))
  
  final_metadata <- dplyr::bind_rows(subtyped_metadata, missing_metadata)
  rownames(final_metadata) <- final_metadata$Cell   #this is most important step, else UMAP labels will be wrong
  
  # Print counts before annotation
  integrated_seurat@meta.data %>% count(cell_type, sub_type, cell_class)
  
  # Re-annotate
  integrated_seurat@meta.data <- final_metadata
  
  # Print counts after annotation
  integrated_seurat@meta.data %>% count(cell_type,sub_type, cell_class)
  
  # # Remove cross labelled cells and cell_class to metadata
  # # Some cells may have cell_type "Myeloid - MDSC" but sub_type as "Myeloid- cDcs".
  # integrated_seurat@meta.data <- integrated_seurat@meta.data %>%
  #   dplyr::mutate(cell_class=dplyr::if_else(cell_type == sub_type |
  #                                               cell_type == "Lymphoid - B" & grepl(pattern="B|Plasma", x=sub_type) |
  #                                               cell_type == "Lymphoid - T" & grepl(pattern="CD4|CD8|Gamma|NKT", x=sub_type) |
  #                                               cell_type == "Myeloid - Macrophages, DCs" & grepl(pattern="Macrophage|cDCs|pDCs", x=sub_type) |
  #                                               cell_type == "Epithelial" & grepl(pattern="Epithelial", x=sub_type) |
  #                                               cell_type == "Fibroblasts" & grepl(pattern="Fibroblasts", x=sub_type),
  #                                             gsub(pattern="\ -.*",replacement="",x=cell_type), "Mixed"))
  
  # Save the object after checking UMAP
  saveRDS(integrated_seurat, paste0(seurat_results, "integrated_seurat_snn.rds"))
}

#******************************************************************************#
#              VISUALIZE CELLTYPES, SUBTYPES, CELLCLASSES IN UMAP              #
#******************************************************************************#

# You must define celltype and split as global variables
visualize_UMAP <- function(){
  
  # Load the integrated seurat object
  integrated_seurat <- base::readRDS(paste0(seurat_results, "integrated_seurat_snn",
                                            dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  # Remove unwanted cells and samples
  integrated_seurat <- subset(x=integrated_seurat,
                              subset=(seurat_class %in% c("Mixed", "Unclassified")),
                              invert=TRUE)
  
  gg1 <- Seurat::DimPlot(object=integrated_seurat,
                         reduction="umap.harmony",
                         cols=my_palette,
                         label=FALSE,
                         group.by="seurat_class",
                         split.by=NULL, #split,
                         shape.by=NULL,
                         pt.size=0.2,
                         label.size=5,
                         repel=FALSE,
                         raster=FALSE) +
    ggplot2::labs(fill="CLUSTERS",
                  x="UMAP_1",
                  y="UMAP_2") +
    my_theme
  
  gg2 <- Seurat::DimPlot(object=integrated_seurat,
                         reduction="umap.harmony",
                         cols=my_palette,
                         label=FALSE,
                         group.by="cell_type",
                         split.by=NULL,#split,
                         shape.by=NULL,
                         pt.size=0.2,
                         label.size=5,
                         repel=FALSE,
                         raster=FALSE) +
    ggplot2::labs(fill="CLUSTERS",
                  x="UMAP_1",
                  y="UMAP_2") +
    my_theme
  
  gg3 <- Seurat::DimPlot(object=integrated_seurat,
                         reduction="umap.harmony",
                         cols=my_palette,
                         label=FALSE,
                         group.by="sub_type",
                         split.by=split,
                         shape.by=NULL,
                         pt.size=0.2,
                         label.size=5,
                         repel=FALSE,
                         raster=FALSE) +
    ggplot2::labs(fill="CLUSTERS",
                  x="UMAP_1",
                  y="UMAP_2") +
    my_theme
  
  ggplot2::ggsave(filename=paste0("UMAP_", celltype, "_final_", split, ".pdf"),
                  plot=gg1+gg2+gg3,
                  device="pdf",
                  path=seurat_results,
                  scale=1,
                  width=dplyr::if_else(is.null(split), 9*3, 9*2),
                  height=dplyr::if_else(is.null(split), 11, 11*3),
                  units=c("in"),
                  dpi=300,
                  limitsize=TRUE,
                  bg=NULL)
}

#******************************************************************************#
#                  VERIFY CLUSTER ANNOTATION USING DOT PLOTS                   #
#******************************************************************************#

# Excel file MUST contain column (i) CELL_TYPE (ii) HUMAN_GENE (iii) MOUSE_GENE
# Excel file MUST be named KMS_Markers.xlsx and located in scripts_path
# Excel file MUST have sheetnames identical to celltype variable (global)
visualize_dotplot <- function(){
  
  # Read file containing marker genes
  kms_markers <- openxlsx::read.xlsx(xlsxFile=paste0(scripts_path, "KMS_Markers.xlsx"),
                                     sheet=celltype)
  
  # Remove all duplicated genes
  feature <- kms_markers %>%
    dplyr::select(dplyr::if_else(species == "Homo sapiens", "HUMAN_GENE", "MOUSE_GENE")) %>%
    dplyr::distinct() %>%
    unlist(use.names=FALSE)
  
  feature <- feature[!is.na(feature)]
  
  if (celltype == "All Markers"){
    celltype <- NULL
  }
  
  # Load the integrated seurat object
  integrated_seurat <- base::readRDS(paste0(seurat_results, "integrated_seurat_snn",
                                            dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  # Remove unwanted cells and samples
  integrated_seurat <- subset(x=integrated_seurat,
                              subset=(cell_class %in% c("Mixed", "Unclassified")),
                              invert=TRUE)
  
  # Set Idents to cell_type so UMAPs show the Idents instead of cluster number
  if (is.null(celltype)){
    Idents(integrated_seurat) <- "cell_type"
  } else {
    Idents(integrated_seurat) <- "sub_type"
  }
  
  # We re-order the active ident alphabetically
  Idents(integrated_seurat) <- base::factor(x=integrated_seurat@active.ident, 
                                            levels=sort(levels(integrated_seurat@active.ident)))
  
  plot_dot <- function(i){
    feature_subset <- feature[((i*20)-19):(i*20)]
    feature_subset <- feature_subset[!is.na(feature_subset)]
    
    Seurat::DotPlot(object=integrated_seurat,
                    assay="RNA",
                    features=feature_subset,
                    #cols=c("blue", "red"),
                    #col.min=-2,
                    #col.max=2,
                    dot.min=0,
                    dot.scale=4,
                    idents=NULL,
                    group.by=NULL,
                    split.by=NULL,
                    cluster.idents=FALSE,
                    scale=TRUE,
                    scale.by="size",
                    scale.min=0,
                    scale.max=100) +
      ggplot2::geom_point(aes(size=pct.exp), shape=21, colour="black", stroke=0.25) + #stroke is width of circle
      ggplot2::scale_colour_distiller(type="div", palette="RdYlGn", direction=-1) +
      ggplot2::guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white", stroke=0.75))) +
      my_theme
  }
  
  purrr::map(.x=1:ceiling(length(feature)/20), 
             .f=plot_dot) %>%
    cowplot::plot_grid(plotlist=.,
                       align="hv",
                       axis="tblr",
                       nrow=ceiling(length(feature)/20),
                       ncol=1,
                       rel_widths=1,
                       rel_heights=1,
                       greedy=TRUE,
                       byrow=TRUE)
  
  ggplot2::ggsave(filename=paste0("Cluster_Verification_Dot_plot_", celltype, ".pdf"),
                  plot=last_plot(),
                  device="pdf",
                  path=seurat_results,
                  scale=1,
                  width=11, #5+0.5*length(feature),
                  height=4*ceiling(length(feature)/20), #1+0.5*length(levels(Idents(integrated_seurat))),
                  units=c("in"),
                  dpi=300,
                  limitsize=FALSE,
                  bg="white")
}

#******************************************************************************#
#                  STEP 13: CALCULATE STATS FOR EACH CLUSTER                   #               
#******************************************************************************#

# NOTE: We identify sparse as well as dominant cells at this step but we 
# remove them after checking all plots in prior steps (10A-C)

calc_stats <- function(){
  wb <- openxlsx::createWorkbook()
  dominant_cells <- data.frame(Cell=NA)
  isolated_cells <- data.frame(Cell=NA)
  
  for (x in 1:length(col_id)){  
    
    # Calculate cells per cluster
    cells_per_cluster <- integrated_seurat@meta.data %>% 
      dplyr::group_by(!!rlang::sym(col_id[x]), Sample) %>%
      dplyr::count() %>%
      dplyr::rename(cluster=col_id[x], nCells=n) %>%
      tidyr::pivot_wider(id_cols=Sample, names_from=cluster, values_from=nCells) %>%
      base::replace(is.na(.), 0)
    
    # Calculate percent of cells in each cluster for each sample
    # Divide by rowSums to adjust for difference in cell number between samples
    cells_per_cluster_percent <- cells_per_cluster %>%
      dplyr::mutate(across(.cols=everything())*100/rowSums(across()))
    
    # Calculate total cells in each cluster across all samples
    cells_per_cluster_total <- c(list(Sample="Total Cells per cluster"), colSums(cells_per_cluster[-1]))
    
    # Merge all data
    cells_per_cluster <- dplyr::bind_rows(cells_per_cluster, data.frame(data=" "), 
                                          cells_per_cluster_percent, data.frame(data=" "),
                                          cells_per_cluster_total) %>%
      dplyr::select(-data)
    
    # Save data
    openxlsx::addWorksheet(wb=wb, sheetName=paste0("Resolution_",all_res[x]))
    openxlsx::writeData(wb=wb, sheet=paste0("Resolution_",all_res[x]), x=cells_per_cluster)
    
    #****************************************************************************#
    # Set a value for outlier classification (outlier > 75% + lenient*IQR)
    # High lenient value means higher upper cutoff. So, filtering is lenient
    lenient <- 25
    
    # Calculate quartiles for each cluster excluding Sample column
    # Column names of transposed quartile  are:  0%, 25%,  50%,  75% and  100%
    # If value > 75% + 100*IQR, then we denote these clusters as outlier 
    quartiles <- apply(X=cells_per_cluster_percent[,-1], MARGIN=2, FUN=stats::quantile) %>%
      t() %>%
      data.frame() %>%
      dplyr::mutate(upper_cutoff=.[[4]] + lenient*(.[[4]] - .[[2]])) %>%
      dplyr::select(upper_cutoff) %>%
      t()
    rownames(quartiles) <- NULL
    
    # Identify dominant cells
    for (i in 1:nrow(cells_per_cluster_percent)) {
      for (j in 1:ncol(quartiles)) {
        if (cells_per_cluster_percent[i, (j+1)] > quartiles[1, j]) {
          
          cluster_id <- colnames(cells_per_cluster_percent)[j+1]
          sample_id <- cells_per_cluster_percent$Sample[i]
          
          cells <- integrated_seurat@meta.data %>%
            dplyr::filter(!!rlang::sym(col_id[x]) == cluster_id & Sample == sample_id) %>%
            dplyr::select(Cell)
          
          dominant_cells <- rbind(dominant_cells, cells)
          
          cat("\nSample:", sample_id, 
              "\tCluster:", cluster_id,
              "\tCluster:", colnames(quartiles)[j],
              "\tcutoff:", quartiles[1,j])
        }
      }
    }
    
    # Identify sparse cells
    cells <- integrated_seurat@meta.data %>% 
      dplyr::add_count(get(col_id[x])) %>%
      dplyr::filter(n < 5) %>%
      dplyr::select(Cell)
    
    isolated_cells <- rbind(isolated_cells, cells)
    
    print(dim(isolated_cells))
  }
  
  # Identify cells consistently dominant in 3 of 6 resolutions)
  dominant_cells <- dominant_cells %>% 
    dplyr::count(Cell) %>%
    dplyr::filter(n >= 3)
  
  # Identify cells consistently sparsely clustered in all 6 resolutions
  isolated_cells <- isolated_cells %>% 
    dplyr::count(Cell) %>% 
    dplyr::filter(n == 6)
  
  # Merge
  if (nrow(dominant_cells) > 0 & nrow(isolated_cells) > 0){
    bad_cells <- dplyr::bind_rows(dominant_cells, isolated_cells) %>%
      dplyr::distinct_at("Cell")
  } else if (nrow(dominant_cells) > 0 & nrow(isolated_cells) == 0){
    bad_cells <- dominant_cells
  } else{
    bad_cells <- data.frame(Cell="NA") 
  }
  
  # Save data
  openxlsx::addWorksheet(wb=wb, sheetName=paste0("Bad cells"))
  openxlsx::writeData(wb=wb, sheet=paste0("Bad cells"), x=bad_cells)
  openxlsx::saveWorkbook(wb=wb, file=paste0(seurat_results, "Cells_per_cluster_", celltype, ".xlsx"), overwrite=TRUE)
  
  # # Next, explore how well our clusters separate by the different PCs.
  # # To visualize this information, we need to extract the UMAP coordinate 
  # # information for the cells along with their corresponding scores for each of 
  # # the PCs to view by UMAP
  # 
  # for(i in 1:5){
  #   j <- 9*i-8
  #   k <- 9*i
  #   Seurat::FeaturePlot(object=integrated_seurat,
  #                       features=c(paste0("PC_",j:k)),
  #                       pt.size=0.4,
  #                       order=TRUE,
  #                       min.cutoff='q10',
  #                       reduction="umap",
  #                       label=TRUE,
  #                       combine=TRUE,
  #                       raster=FALSE)
  #   
  #   ggplot2::ggsave(filename=paste0("UMAP_for_PC_", j, "_through_", k, ".jpg"),
  #                   plot=last_plot(),
  #                   device="jpeg",
  #                   path=diagnostics_path,
  #                   scale=1,
  #                   width=8.5,
  #                   height=11,
  #                   units=c("in"),
  #                   dpi=600,
  #                   limitsize=TRUE,
  #                   bg="white")
  # }
  
  #*******STEP 10E: CHECK IF CELL CYCLE IS A SOURCE OF VARIATION VIA PCs*******#
  
  # # Save the top 30 (positive) and bottom 30 (negative) genes for each of 50 PCs
  # # Rows will be genes and columns will be PCs
  # pc_genes <- data.frame()[1:30, ]
  # rownames(pc_genes) <- paste0("Gene_",1:30)
  # for (i in 1:50){
  #   pc_genes[,2*i-1] <- names(sort(integrated_seurat@reductions$pca@feature.loadings[,i], decreasing=FALSE)[1:30])
  #   colnames(pc_genes)[2*i-1] <- paste0("Negative_PC_",i)
  #   pc_genes[,2*i] <- names(sort(integrated_seurat@reductions$pca@feature.loadings[,i], decreasing=TRUE)[1:30])
  #   colnames(pc_genes)[2*i] <- paste0("Positive_PC_",i)
  # }
  # 
  # # Find which PCs are affected by cell cycle genes
  # PCs_affected <- pc_genes
  # for(i in 1:ncol(PCs_affected)){
  #   PCs_affected[,i] <- pc_genes[,i] %in% cell_cycle_genes
  # }
  # 
  # # Populate the TRUE values with genes names
  # for(i in 1:ncol(PCs_affected)){
  #   for(j in 1:nrow(PCs_affected)){
  #     if(PCs_affected[j,i] == TRUE){
  #       PCs_affected[j,i] <- pc_genes[j,i]
  #     } else{
  #       PCs_affected[j,i] <- ""
  #     }
  #   }
  # }
  # 
  # # Remove PCs with no genes
  # PCs_affected <- PCs_affected[, colSums(PCs_affected != "") != 0]
  # 
  # # Save the dataframe in xlsx format
  # wb <- openxlsx::createWorkbook()
  # openxlsx::addWorksheet(wb=wb, sheetName="PCA_genes")
  # openxlsx::writeData(wb=wb, sheet="PCA_genes", x=pc_genes)
  # openxlsx::addWorksheet(wb=wb, sheetName="Cell_cycle_affected_PCs")
  # openxlsx::writeData(wb=wb, sheet="Cell_cycle_affected_PCs", x=PCs_affected)
  # openxlsx::saveWorkbook(wb=wb, file=paste0(diagnostics_path, "PCA_Genes_After_Integration.xlsx"), overwrite=TRUE)
  
  #*****STEP 10F: CHECK IF CELL CYCLE IS A SOURCE OF VARIATION VIA GRAPHS******#
  
  # integrated_seurat_cc <- RunPCA(object=integrated_seurat,
  #                                assay=NULL,
  #                                features=c(s_genes, g2m_genes),
  #                                npcs=50,
  #                                rev.pca=FALSE,
  #                                weight.by.var=TRUE,
  #                                verbose=TRUE,
  #                                ndims.print=1:5,
  #                                nfeatures.print=30,
  #                                reduction.name="pca",
  #                                reduction.key="PC_",
  #                                seed.use=42)
  # 
  # # Plot the PCA colored by cell cycle phase
  # Seurat::DimPlot(object=integrated_seurat_cc, 
  #                 reduction="pca", 
  #                 group.by="Phase",
  #                 split.by="Sample",
  #                 raster=FALSE,
  #                 ncol=dplyr::if_else(length(unique(integrated_seurat_cc@meta.data$Sample)) < 3, 2, 3),
  #                 combine=TRUE) +
  #   #facet_wrap(.~Sample, nrow=4) +				
  #   NoLegend()  
  # 
  # # Save the plot
  # ggplot2::ggsave(filename="PCA_using_cell_cycle_genes.jpg",
  #                 plot=last_plot(),
  #                 device="jpeg",
  #                 path=diagnostics_path,
  #                 scale=1,
  #                 width=8.5,
  #                 height=11,
  #                 units=c("in"),
  #                 dpi=600,
  #                 limitsize=TRUE,
  #                 bg="white")
  
  #********STEP 10G: DETERMINE THE 'DIMENSIONALITY' FROM HEATMAP OF PCS********#
  
  # # cells parameter in DimHeatmap() specifies the number of cells with most 
  # # negative or postive PCA scores to be used for plotting. The idea is that we 
  # # are looking for a PC where the heatmap starts to look more “fuzzy”, i.e. where
  # # the distinctions between the groups of genes is not so distinct.
  # for(i in 1:5){
  #   j <- 9*i-8
  #   k <- 9*i
  #   Seurat::DimHeatmap(object=integrated_seurat,
  #                      dims=j:k,
  #                      nfeatures=30,
  #                      cells=1000,
  #                      reduction="pca",
  #                      disp.min=-2.5,
  #                      disp.max=NULL,
  #                      balanced=TRUE,
  #                      projected=FALSE,
  #                      ncol=3,
  #                      fast=FALSE,
  #                      raster=FALSE,
  #                      slot="scale.data",
  #                      assays=NULL,
  #                      combine=TRUE)
  #   
  #   # Save the plot
  #   ggplot2::ggsave(filename=paste0("PCA_", j, "_through_", k, ".jpg"),
  #                   plot=last_plot(),
  #                   device="jpeg",
  #                   path=diagnostics_path,
  #                   scale=1,
  #                   width=8.5,
  #                   height=11,
  #                   units=c("in"),
  #                   dpi=600,
  #                   limitsize=TRUE,
  #                   bg="white")
  # }
  
  #**********STEP 10H: DETERMINE THE 'DIMENSIONALITY' FROM ELBOW PLOTS*********#
  
  # # We can quantitatively determine the dimensionality using Elbow plot
  # # We can calculate where the principal components start to elbow by taking the
  # # larger value of:
  # # 1. The point where the principal components only contribute 5% of standard 
  # # deviation and the principal components cumulatively contribute 90% of the 
  # # standard deviation.
  # # 2. The point where the percent change in variation between the consecutive 
  # # PCs is less than 0.1%
  # 
  # # Determine percent of variation associated with each PC
  # pct <- integrated_seurat@reductions$pca@stdev / sum(integrated_seurat@reductions$pca@stdev) * 100
  # # Calculate cumulative percents for each PC
  # cumu <- cumsum(pct)
  # # Determine at which PC the cumulative percent is greater than 90% and 
  # # % variation associated with the PC is less than 5%
  # co1 <- which(cumu > 90 & pct < 5)[1]
  # 
  # # Determine the difference between variation of PC and subsequent PC
  # # We add a 0 to end of diff so that there are equal number of rows in diff, pct & cumu
  # diff <- c(pct[1:length(pct)-1]-pct[2:length(pct)],0) 
  # # Find the last PC which varies from its next PC by more than 0.1%
  # co2 <- max(which(diff > 0.1))
  # 
  # # Minimum of the two calculation
  # pcs <- min(co1, co2+1)
  # cat("\nThe minimum number of PCs to used is ", pcs, "\n")
  # 
  # # Elbow plot to visualize
  # # Create a dataframe with values
  # plot_df <- data.frame(pct=pct,
  #                       cumu=cumu,
  #                       diff=diff,
  #                       rank=1:length(pct))
  # 
  # ggplot(data=plot_df, aes(x=cumu, y=pct, label=rank, color=rank>pcs)) +
  #   geom_text() +
  #   theme_classic() +       
  #   labs(x="Cumulative Sum of Std Dev",
  #        y=" Percent Difference in Std Dev",
  #        title=paste0("Elbow plot of PCs")) +
  #   theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1, size=10), 
  #         axis.text.y=element_text(size=10),                                 
  #         axis.title=element_text(size=14),                                  
  #         plot.title=element_text(hjust=0.5, size=16, face="bold")) +       
  #   geom_vline(xintercept=90, linetype=2, color="red") +
  #   geom_hline(yintercept=min(pct[pct > 5]), linetype=2, color="red")
  # 
  # ggplot2::ggsave(filename="Elbow_plot_of_PCA_1_through_50.jpg",
  #                 plot=last_plot(),
  #                 device="jpeg",
  #                 path=diagnostics_path,
  #                 scale=1,
  #                 width=8.5,
  #                 height=11,
  #                 units=c("in"),
  #                 dpi=600,
  #                 limitsize=TRUE,
  #                 bg="white")
  # 
  # # NOTE: We have to use a minimum dimension as determined above. However, using
  # # more PCs is not a bad thing. Clustering will just take more time with more PCs.
  
  return(bad_cells)
}

plot_genes <- function(res, reduc, celltype, sheetname){
  
  idents <- paste0("cluster.", res, ".", base::tolower(reduc))
  
  # Load the integrated seurat object
  integrated_seurat <- base::readRDS(paste0(seurat_results, "integrated_seurat_snn", 
                                            dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  # Set default assay
  DefaultAssay(integrated_seurat) <- "RNA"
  
  # Set identity to an existing column in meta data
  Idents(object=integrated_seurat) <- idents
  
  # Read markers
  features <- openxlsx::read.xlsx(xlsxFile=paste0(scripts_path, "KMS_Markers.xlsx"),
                                  sheet=sheetname)
  
  features <- as.list(features %>%
                        dplyr::select(dplyr::if_else(species == "Homo sapiens", "HUMAN_GENE", "MOUSE_GENE")))[[1]]
  
  features <- rownames(integrated_seurat@assays$RNA$data)[tolower(rownames(integrated_seurat@assays$RNA$data)) %in% 
                                                            tolower(features)]
  
  feature_plot <- function(i){
    
    split <- "Condition" #NULL
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
  for (n in 1:ceiling((length(features)/10))){
    j <- 10*n-9
    k <- 10*n
    vec <- features[j:k]
    purrr::map(.x=vec[!is.na(vec)], 
               .f=feature_plot) %>%
      cowplot::plot_grid(plotlist=.,
                         align="hv",
                         axis="tblr",
                         nrow=2,
                         ncol=5,
                         rel_widths=1,
                         rel_heights=1,
                         greedy=TRUE,
                         byrow=TRUE)
    
    ggplot2::ggsave(filename=paste0("Feature_plot_",n, ".jpg"),
                    plot=last_plot(),
                    device="jpeg",
                    #path=diagnostics_path,
                    width=8.5*4,
                    height=11*2,
                    units=c("in"),
                    dpi=300,
                    limitsize=FALSE,
                    bg="white")
  }
}

#******************************************************************************#
#              STEP 14: IDENTIFY CONSERVED MARKERS FOR EACH CLUSTER            #
#******************************************************************************#

# # Set default assay
# DefaultAssay(integrated_seurat) <- "RNA"
# 
# # Set identity to an existing column in meta data
# Idents(object=integrated_seurat) <- res
# 
# # Create function to get conserved markers for any given cluster
# get_conserved <- function(cluster){
#   # min.cells.group=1 doesnt work in FindConservedMarkers.
#   # So, we manually skip clusters that have less than 3 cells in any sample
#   if (nrow(integrated_seurat@meta.data %>%
#            dplyr::select(Sample, !!rlang::sym(res)) %>%
#            dplyr::filter(!!rlang::sym(res) == cluster) %>%
#            count(Sample) %>%
#            dplyr::filter(n<3)) == 0){
#     df <- Seurat::FindConservedMarkers(object=integrated_seurat,
#                                        ident.1=cluster,
#                                        ident.2=NULL,
#                                        grouping.var="Sample",
#                                        assay="RNA",
#                                        slot="data",
#                                        min.cells.group=1,
#                                        meta.method=metap::minimump,
#                                        verbose=TRUE,
#                                        features=NULL,
#                                        logfc.threshold=0.25,
#                                        test.use="wilcox",
#                                        min.pct=0.1,
#                                        min.diff.pct=0.25,
#                                        node=NULL,
#                                        only.pos=TRUE,
#                                        max.cells.per.ident=Inf,
#                                        random.seed=1,
#                                        latent.vars=NULL,
#                                        min.cells.feature=3,
#                                        pseudocount.use=1,
#                                        mean.fxn=NULL,
#                                        fc.name=NULL,
#                                        base=2,
#                                        return.thresh=0.01,
#                                        densify=FALSE)
#     
#     if(nrow(df) != 0){
#       df %>%  tibble::rownames_to_column(var="gene") %>%    #add rownames as 1st column named "gene"
#         dplyr::left_join(y=annotations,
#                          by=c("gene"="SYMBOL")) %>%
#         # merge the dataframe generated by FindConservedMarkers() with unique
#         #genes from "annotations" dataframe columns - gene_name and description
#         cbind(cluster=cluster, .)     # add the cluster variable as 1st column named cluster
#     }
#   }
# }
# 
# # Iterate function across all clusters. We want the output of the map family
# # of functions to be a dataframe with each cluster output bound together by
# # rows, we will use the map_dfr() function
# # map_dfr(inputs_to_function, name_of_function)
# conserved_markers <-
#   map_dfr(rev(levels(integrated_seurat@active.ident)), get_conserved)
# 
# # Save conserved markers
# wb <- openxlsx::createWorkbook()
# openxlsx::addWorksheet(wb=wb1, sheetName="Conserved_Markers")
# openxlsx::writeData(wb=wb1, sheet="Conserved_Markers", x=conserved_markers)
# openxlsx::saveWorkbook(wb=wb, file=paste0(seurat_results,"Markers_Conserved.xlsx"), overwrite=TRUE)

#***************(RECOMMENDED): REMOVING BATCH SPECIFIC CLUSTERS**************#

# Sometimes, you will see sample/batch specific clusters where cells of a 
# particular sample dominate. Rather than removing all cells from such a cluster
# we can remove ONLY the cells belonging to dominating sample at all 
# resolutions and re-run entire pipeline for better clustering.

# In example below, we can just remove cells from N6-0, N5-1, N5-11 and N4-12

# Sample  0	            1	            10	          11  	       12
# N1	    40	          38	          240	          0	            130
# N2	    20	          222	          230	          4   	        78
# N3	    632	          350	          2639	        8	            325
# N4	    13	          34	          872	          0	            "5190"
# N5	    399	          "28556"       487	          "6241"	      96
# N6	    "52043"	       562	        1900	        47	          300
# 
# N1	    0.536408743	  0.509588306	  3.218452461	  0	            1.743328416
# N2	    0.109110747	  1.211129296	  1.254773595	  0.021822149	  0.425531915
# N3	    0.854920528	  0.473452824	  3.569834292	  0.010821779	  0.439634765
# N4	    0.098776689	  0.258339032 	6.62563635	  0	            39.43469341
# N5	    0.56627874	  40.52795913	  0.691172296	  8.857507806	  0.136247516
# N6	    71.15045458	  0.768336865	  2.597580149	  0.06425593	  0.410144234
# 
# Total   53147	        29762	        6368	        6300	        6119 

# # Skip round 2 if no bad cells are present
# if (nrow(bad_cells) == 1){
#   break
# } else{
#   # Remove dominant and isolated cells
#   filtered_seurat <- subset(x=integ_data,
#                             Cell %in% bad_cells$Cell, 
#                             invert=TRUE)
#   
#   # Change default assay to RNA and remove other assays.
#   # You need to perform SCTransform() again as removing cells will alter 
#   # Pearson residuals -> UMI corrected value for some genes in some cells that
#   # were earlier 0.
#   DefaultAssay(object=filtered_seurat) <- "RNA"
#   filtered_seurat[["SCT"]] <- NULL   
#   filtered_seurat[["integrated"]] <- NULL
# }

#******************************************************************************#
#               OPTIONAL: IDENTIFYING GOOD MARKERS FOR CELL TYPES              # 
#******************************************************************************#

ident_good_markers <- function(){
  # NOTE: You can run this section ONCE to generate a list of good markers and
  # use these markers for ALL SUBSEQUENT data sets. If these markers DO NOT work
  # well for a particular dataset, you can re-run this section using that dataset
  # to generate good markers specific for that data set.
  
  # NOTE: Epithelial Cells and Fibroblasts show vast changes upon BBN condition.
  # So, it might be better to subset the seurat object before finding good markers
  
  # Load the integrated seurat object
  integrated_seurat <- base::readRDS(paste0(seurat_results, "integrated_seurat_snn", 
                                            dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  # Remove unwanted cells and samples before finding good markers
  integrated_seurat1 <- subset(x=integrated_seurat,
                               Condition == "Tumor")
  
  # Create matrix similar to read data for each cluster to identify closely related clusters
  cor_meta_data <- integrated_seurat1@meta.data %>%
    dplyr::distinct(integrated_snn_res.1.4, .keep_all=TRUE) %>%
    tibble::remove_rownames(.) %>%
    tibble::column_to_rownames("integrated_snn_res.1.4")
  
  clusters <- as.vector(unlist(integrated_seurat1@meta.data %>% dplyr::distinct(integrated_snn_res.1.4)))
  cor_read_data <- matrix(NA, nrow=nrow(integrated_seurat1@assays$RNA@counts), ncol=length(clusters))
  rownames(cor_read_data) <- rownames(integrated_seurat1@assays$RNA@counts)
  colnames(cor_read_data) <- clusters
  for(i in 1:length(clusters)){
    
    # Create a list of cells for each sample
    cor_cells_subset <- rownames(integrated_seurat1@meta.data %>% dplyr::filter(integrated_snn_res.1.4 == as.numeric(clusters[i])))
    
    # Get the raw counts of these cells. counts are stored in sparse matrix format
    # So, use data.frame to convert "." in sparse matrix to "0"
    cor_subset <- data.frame(integrated_seurat1@assays$RNA@counts[,cor_cells_subset])
    cor_read_data[,i]  <- rowSums(cor_subset)
  }
  
  # Remove unwanted genes
  unwanted_genes <- rownames(cor_read_data)[grep(pattern="^RP[SL]|RIK$|^MT-|^GM[0-9.]+$", x= rownames(cor_read_data), ignore.case=TRUE)]
  keep_genes <- setdiff(rownames(cor_read_data), unwanted_genes)  
  cor_read_data <- cor_read_data[keep_genes,]
  
  # Create DESeq2 object so we can use DESeq2::rlog() or DESeq2::vst() on it
  cor_dds <- DESeqDataSetFromMatrix(countData=cor_read_data,
                                    colData=cor_meta_data,
                                    design=~1)
  # Use vst or rld
  vst <- DESeq2::vst(cor_dds)
  head(assay(vst))
  sampleDists <- dist(t(assay(vst)))
  sampleDistMatrix <- as.matrix(sampleDists)
  pheatmap::pheatmap(mat=sampleDistMatrix,
                     scale="none",
                     cutree_rows=20,
                     cutree_cols=20,
                     cluster_rows=TRUE,   #cluster the rows
                     cluster_cols=TRUE,   #cluster the columns
                     fontsize=8, 
                     fontsize_row=8, 
                     fontsize_col=8,
                     angle_col=c("270", "0", "45", "90", "315"),
                     fontsize_number=0.8*fontsize,
                     filename=paste0(seurat_results, "Cluster_Correlation", "_", celltype, ".jpg"))
  
  # Use the UMAP at res 1.4 as well as above correlation heatmap to group clusters
  res <- "integrated_snn_res.1.4"
  
  # BBN_C57B6 scRNASeq res 1.4 (only BBN samples)
  # Refer ppt to understand how the clusters were grouped
  clusters <- list("cluster1"=c(13,16,17,23,29,35),     #Fibro_1
                   "cluster2"=c(5,6,26,30,33,42),       #Fibro_2
                   "cluster3"=c(36),                    #Myocytes
                   "cluster4"=c(0,8,11,14,15,18,25,34,39), #Epithelial
                   "cluster5"=c(4,7,12,21),             #MDSCs
                   "cluster6"=c(9,10,20),               #Macrophages
                   "cluster7"=c(32,37),                 #DCs
                   "cluster8"=c(19,22,24,27,28),        #T cells
                   "cluster9"=c(31),                    #NK cells
                   "cluster10"=c(2,48),                 #B cells
                   "cluster11"=c(1,46),                 #Endo
                   "cluster12"=c(38),                   #Lymp. Endo
                   "cluster13"=c(43),                   #Neurons??
                   "mixed"=c(3,40,41,50,51,52,53,54,44,45,47,49)) #contaminants
  
  # BBN_C57B6 scRNASeq Epithelial res 1.4 (only BBN samples)
  # Refer ppt to understand how the clusters were grouped
  clusters <- list("cluster1"=c(4,18),     
                   "cluster2"=c(1,8,9,20,21),      
                   "cluster3"=c(3,12,19,22), 
                   "cluster4"=c(7,23),             
                   "cluster5"=c(2,10),               
                   "cluster6"=c(0,15),                 
                   "cluster7"=c(5,11,17),        
                   "cluster8"=c(6,16),                    
                   "cluster9"=c(13),                 
                   "cluster10"=c(14),                 
                   "mixed"=c(24,25,26,27,28,30))
  
  # BBN_C57B6 scRNASeq Fibroblast res 1.4 (only BBN samples)
  # Refer ppt to understand how the clusters were grouped
  clusters <- list("cluster1"=c(0,4,5),     
                   "cluster2"=c(3,6,12,14,18),      
                   "cluster3"=c(7,8,9,13,15), 
                   "cluster4"=c(2,10,17),             
                   "cluster5"=c(11),               
                   "cluster6"=c(1,16),                 
                   "cluster7"=c(20),
                   "mixed"=c(19))
  
  # BBN_C57B6 scRNASeq Myeloid res 1.4 (only BBN samples)
  # Refer ppt to understand how the clusters were grouped
  clusters <- list("cluster1"=c(2,4),     
                   "cluster2"=c(0,1,8,11),      
                   "cluster3"=c(5,6,9,10,19,20), 
                   "cluster4"=c(3,7,12,21,24),             
                   "cluster5"=c(13,14,15,16,22),               
                   "cluster6"=c(17), 
                   "cluster7"=c(18), 
                   "cluster8"=c(25),
                   "cluster9"=c(26),
                   "mixed"=c(23,27,28,29))
  
  
  # Make sure you have assigned all clusters to one of the cell types
  # We subtract 1 because (...find explanation..)
  list_1 <- sort(unique(as.numeric(unlist(integrated_seurat1@meta.data %>% 
                                            dplyr::select(all_of(res)), use.names=FALSE))))-1
  list_2 <- sort(unlist(clusters, use.names=FALSE))
  
  # Proceed with annotation ONLY if all clusters have been renamed
  if (identical(list_1, list_2)){
    "All Clusters have been annotated"
    
    # Extract metadata from Seurat object, assign appropriate resolution to
    # seurat_clusters column and add cell_type column
    data <- integrated_seurat1@meta.data %>% 
      dplyr::mutate(seurat_clusters=!!rlang::sym(res),
                    cell_type=NA)
    
    # Assign cell type based on cluster numbers within seurat_clusters column
    for (j in 1:nrow(data)){
      for (i in 1:length(clusters)){
        if (as.numeric(as.character(data$seurat_clusters[j])) %in% clusters[[i]]){
          data$cell_type[j] <- names(clusters[i])
        }
      }
    }
    
    # Check summary of cell counts
    print(data %>% dplyr::count(cell_type) %>% dplyr::arrange(n))
    cat("\n")
    
    # Import the metadata into Seurat object and save it
    integrated_seurat1@meta.data <- data
  } else {
    "Some clusters have not been annotated yet"
  }
  
  # STEP 2: Run FindAllMarkers() on these major clusters
  
  # Set default assay
  DefaultAssay(integrated_seurat1) <- "RNA"
  
  # Set identity to an existing column in meta data
  Idents(object=integrated_seurat1) <- "cell_type"
  
  #Remove unwanted cells and samples before finding good markers
  integrated_seurat1 <- subset(x=integrated_seurat1,
                               subset=(cell_type == "mixed"), 
                               invert=TRUE)
  
  # Get all markers for any given cluster
  all_markers <- Seurat::FindAllMarkers(object=integrated_seurat1,
                                        assay="RNA",
                                        features=NULL,
                                        logfc.threshold=0.25,
                                        test.use="wilcox",
                                        slot="data",
                                        min.pct=0.1,
                                        min.diff.pct=0.1,
                                        node=NULL,
                                        verbose=TRUE,
                                        only.pos=TRUE,
                                        max.cells.per.ident=Inf,
                                        random.seed=1,
                                        latent.vars=NULL,
                                        min.cells.feature=3,
                                        min.cells.group=1,
                                        pseudocount.use=1,
                                        mean.fxn=NULL,
                                        fc.name=NULL,
                                        base=2,
                                        return.thresh=0.01,
                                        densify=FALSE)
  
  all_markers <- all_markers %>% 
    dplyr::mutate(pct.1=dplyr::if_else(pct.1 == 0, 0.001, pct.1),
                  pct.2=dplyr::if_else(pct.2 == 0, 0.001, pct.2),
                  ratio=pct.1/pct.2) %>%
    dplyr::left_join(y=unique(annotations[, c("SYMBOL", "chr", "description")]), by=c("gene"="SYMBOL")) %>%
    dplyr::relocate(cluster, gene, chr, avg_log2FC, p_val, p_val_adj, pct.1, pct.2, ratio, description)
  
  # STEP 3: Find top markers for each major cluster
  # (i) statistically significant (p_val_adj < 0.05)
  # (ii) low expression in other cell types (high ratio)
  # (iii) expressed at decent levels (avg_log2FC >= 0.58)
  # NOTE: We are not using pct.1 as a good marker may not be expressed in most of
  # cells in the major cluster as we have merged many smaller clusters together
  top_markers <- all_markers %>%
    dplyr::filter(avg_log2FC >= 0.58 & p_val_adj < 0.05) %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(desc(ratio), desc(avg_log2FC)) %>%
    dplyr::slice_head(n=20) %>%
    ungroup()
  
  # Save all the markers
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb=wb, sheetName="All_Markers")
  openxlsx::writeData(wb=wb, sheet="All_Markers", x=all_markers)
  openxlsx::addWorksheet(wb=wb, sheetName="Top_Markers")
  openxlsx::writeData(wb=wb, sheet="Top_Markers", x=top_markers)
  openxlsx::saveWorkbook(wb=wb, file=paste0(seurat_results, "Markers_for_major_clusters_", celltype, ".xlsx"), overwrite=TRUE)
  
  # Step 4: Manually rename clusters with appropriate cell types using info 
  # below in the xlsx file and use it for future annotations. Save this info
  # to "KMS_Markers.xlsx" in sheet named "All Markers"
  
  # Fibroblasts express several types of collagens (Col*)
  # Epithelial cells express several types of keratins (Krt*)
  # T cells express CD3 genes (Cd3d, Cd3e, Cd3g)
  # B cells express CD79 genes (Cd79a, Cd79b)
  # Macrophages express C1q genes (C1qa, C1qb, C1qc)
  # Mast Cells express S100A8, S100A9, Retnlg at high levels
  # Myocytes express myosin genes (Myh11, Myl9)
  # Endothelial cells express endothelial proteins (Esam, Ecsr, Emcn, Nos3)
  # Neurons express several potassium gated channels (Kcna1, Kcna2, Kcna6)
  
  # Look at UMAPs plotted from FeaturePlot() below and remove bad markers.
  
  features <- read.xlsx(paste0(seurat_results, "Markers_for_major_clusters_", celltype, ".xlsx"), sheet="Top_Markers")
  features <- feature$gene
  features <- intersect(features, rownames(integrated_seurat@assays$RNA$data))
  length(features)
  
  for(i in 1:ceiling(length(features)/10)){
    j <- 10*i-9
    k <- 10*i
    Seurat::FeaturePlot(object=plot_object,  #integrated_seurat,
                        features=features[j:k],
                        pt.size=0.4,
                        order=TRUE,
                        min.cutoff='q10',
                        reduction="umap.harmony",
                        # cols can handle only upto 3 colors.
                        #cols = rev(brewer.pal(n=11, name="RdBu"))[5:11],
                        cols= c("grey", viridis(n=10, option="C", direction=-1)),
                        ncol =5,
                        label=FALSE,
                        combine=TRUE,
                        raster=FALSE) +
      scale_colour_gradientn(colours=rev(brewer.pal(n=11, name="RdBu"))[5:11])
    
    ggplot2::ggsave(filename=paste0(i,".jpg"),
                    plot=last_plot(),
                    device="jpeg",
                    #path=diagnostics_path,
                    scale=1,
                    width=22,
                    height=11,
                    units=c("in"),
                    dpi=600,
                    limitsize=TRUE,
                    bg="white")
  }
}

#******************************************************************************#
#         PREPARE READ AND META DATA FROM SEURAT OBJECT AND RUN DESEQ2         #
#******************************************************************************#

# This function reads a seurat object and generates read_data, meta_data, file_suffix
# for use by analyze_DESeq2(). It automatically calls analyze_DESeq2() too.

# prep_DESeq2 <- function(celltype){ 
#   
#   cat("\n**********************", celltype, "Analysis**********************\n")
#   
#   # Load the integrated seurat object
#   integrated_seurat <- base::readRDS(paste0(seurat_results, "integrated_seurat_snn", 
#                                             dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
#   
#   #******************************IMPORT META DATA******************************#
#   
#   integrated_seurat <- subset(x=integrated_seurat,
#                               subset=cell_class %in% c("Mixed", "Unclassified"),
#                               invert=TRUE)
#   
#   integrated_seurat <- subset(x=integrated_seurat,
#                               subset=cell_type %in% celltype)
#   
#   
#   # Perform analysis on cell types like T cell, B cell rather than Myeloid, Lymphoid
#   
#   subtypes <- integrated_seurat@meta.data %>% 
#     dplyr::select(cell_type) %>% 
#     unlist(use.names=FALSE) %>%
#     unique()
#   
#   for (subtype in subtypes){
#     
#     subset_seurat <- subset(x=integrated_seurat,
#                             subset=cell_type == subtype)
#     
#     # Subset metadata
#     meta_data <- subset_seurat@meta.data %>%
#       dplyr::distinct(Sample, .keep_all=TRUE) %>%
#       dplyr::mutate(Batch=1)
#     
#     #******************************IMPORT READ DATA******************************#                                           
#     
#     # The read data will have "the reads of all cells belonging to a single 
#     # sample" merged together in each column. First, create a list of samples
#     samples <- subset_seurat@meta.data %>% 
#       dplyr::select(Sample) %>% 
#       unlist(., use.names=FALSE) %>% 
#       unique()
#     
#     # Second, create an empty dataframe with rows=genes and columns=samples
#     read_data <- data.frame(matrix(NA, nrow=nrow(subset_seurat@assays$RNA$counts), ncol=nrow(meta_data)))
#     rownames(read_data) <- rownames(subset_seurat@assays$RNA$counts)
#     colnames(read_data) <- samples
#     
#     # Thirdly, we will add row-wise, the counts of each gene for each sample
#     for(i in samples){
#       
#       # Create a list of cells for each sample
#       cells_subset <- rownames(subset_seurat@meta.data %>% dplyr::filter(Sample == i))
#       
#       # Use data.frame to convert "." in sparse matrix to "0"
#       subset <- data.frame(subset_seurat@assays$RNA$counts[,cells_subset])
#       read_data[,i]  <- rowSums(subset)
#     }
#     
#     read_data <- read_data %>% 
#       tibble::rownames_to_column("SYMBOL")
#     
#     file_suffix <- subtype
#     #analyze_DESeq2(meta_data, read_data, file_suffix)
#     
#     annotations <- get_annotations(species)
#     meta_data <- prep_metadata(meta_data, Variable)
#     read_data <- prep_readdata(read_data, meta_data)
#     l <- check_data(read_data, meta_data)
#     meta_data <- l[[2]]
#     read_data <- l[[1]]
#   }
# }

#*********************************CLEAR MEMORY*********************************#



#******************************************************************************#

#********STEP 7A: DECLARE REFERENCE SAMPLES FOR INTEGRATING THE DATA*********#

# # NOTE: We assign index id and then select the index ids for reference samples
# 
# if (is.na(ref1) & is.na(ref2)){
#   ref_samples <- NULL
# } else if(!is.na(ref1) & is.na(ref2)){
#   ref_samples <- as.vector(sct@meta.data %>% 
#                              dplyr::distinct_at("Sample", .keep_all=TRUE) %>%
#                              dplyr::filter(!(Sample %in% samples_with_few_cells)) %>%
#                              dplyr::mutate(index=row_number()) %>%
#                              dplyr::filter(!!rlang::sym(ref1) == ref1_value) %>%
#                              dplyr::select(index))[[1]]
# } else{
#   ref_samples <- as.vector(sct@meta.data %>% 
#                              dplyr::distinct_at("Sample", .keep_all=TRUE) %>%
#                              dplyr::filter(!(Sample %in% samples_with_few_cells)) %>%
#                              dplyr::mutate(index=row_number()) %>%
#                              dplyr::filter(!!rlang::sym(ref1) == ref1_value, !!rlang::sym(ref2) == ref2_value) %>%
#                              dplyr::select(index))[[1]]
# }  
# 
# cat("\nReference samples are:")
# for (i in ref_samples){
#   cat("\n", i, ":", unique(sct@meta.data$Sample)[i])  
# }

#******************************************************************************#
#                   REMAP CLUSTERS WITH CELL TYPES                             #
#******************************************************************************#

# # We annotated cells using AddModuleScore() and saved them in seurat_class
# df <- table(integrated_seurat@meta.data$seurat_class,
#             integrated_seurat@meta.data$cluster.0.4.harmony)

# # When you convert the table to data.frame(), it automatically calculates frequencies
# df <- df %>%
#   data.frame() %>%
#   dplyr::group_by(Var2) %>%
#   dplyr::filter(Freq == max(Freq)) %>%
#   dplyr::ungroup() %>%
#   dplyr::select(Var1, Var2) %>%
#   data.frame() %>%
#   dplyr::rename(cell_type=Var1, cluster.0.4.harmony=Var2)

# # Remap celltypes to cluster.1.4.harmony
# metadata <- integrated_seurat@meta.data %>%
#   dplyr::left_join(df,by=("cluster.1.4.harmony"="cluster.1.4.harmony"))

# rownames(metadata) <- metadata$Cell
# integrated_seurat@meta.data <- metadata
# 
# # Plot UMAP
# Seurat::DimPlot(object=integrated_seurat,
#                 reduction="umap.harmony",
#                 cols=my_palette,
#                 label=FALSE,
#                 group.by="cell_type",
#                 split.by=NULL, #split,
#                 shape.by=NULL,
#                 pt.size=0.2,
#                 label.size=5,
#                 repel=FALSE,
#                 raster=FALSE) +
#   ggplot2::labs(fill="CLUSTERS",
#                 x="UMAP_1",
#                 y="UMAP_2") +
#   my_theme