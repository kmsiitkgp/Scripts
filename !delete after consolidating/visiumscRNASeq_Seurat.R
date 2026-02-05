#!/usr/bin/env Rscript

# NOTE: All variables and functions are defined within the file below

source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")

# NOTE: In visium v1 spatial transcriptomics, there are 4992 barcoded spots 
# which capture RNA from tissue placed above these spots. 
# So, unlike scRNASeq, each barcode corresponds to a spot, not a cell.
# https://www.youtube.com/watch?v=VwNk4d-0RJc
# https://kb.10xgenomics.com/hc/en-us/articles/360035848191-How-many-spots-are-
# within-a-single-capture-area-on-the-Visium-v1-Spatial-Gene-Expression-Slide

#******************************************************************************#
#                       STEP 1: SETUP THE SEURAT OBJECT                        #
#******************************************************************************#

#************************IMPORTING DATA FROM h5AD FILE*************************#

# data.dir MUST have a folder named "spatial" with the image and H5 file 
# specified by filename parameter.

# Create a list of samples which will be added to each barcode.
# Since folder names correspond to sample name, we just use list.files()
samples <- list.files(path = feature_matrix_path)

# Loop through each of the individual folders in parent directory & import data
for(i in samples){
  
  sample.seurat <- Load10X_Spatial(data.dir = paste0(feature_matrix_path, i),
                                   filename = "raw_feature_bc_matrix.h5",
                                   assay = "Spatial",
                                   slice = i,
                                   filter.matrix = TRUE,
                                   to.upper = FALSE,
                                   image = NULL)
  # if filter.matrix is set to FALSE, all 4992 spots will be recorded in 
  #sample.seurat@images$B8@coordinates. Else, only spots that are over tissue,
  # will be recorded in sample.seurat@images$B8@coordinates.
  
  # Unlike Seurat::Read10X(), Seurat::Load10X_Spatial doesnt have ability to 
  # specify project parameter. So, we manually do it.
  sample.seurat@meta.data <- sample.seurat@meta.data %>% 
    dplyr::mutate(orig.ident = i)
  
  #print(dim(sample.seurat@meta.data))
  
  # Assign the seurat object to its corresponding variable
  assign(i, sample.seurat)
  
  # Explore the meta.data slot
  cat("\nFirst few rows of ", i, "\n")
  print(head(sample.seurat@meta.data))
  cat("\nLast few rows of ", i, "\n")
  print(tail(sample.seurat@meta.data))
}

#******************************************************************************#
#                           STEP 2: QUALITY CONTROL                            #
#******************************************************************************#

#**********************STEP 2A: CALCULATE ALL QC METRICS***********************#

# NOTE:  We are going to do the same QC on each sample. So, we can merge the 
# individual seurat objects into a single seurat object and calculate metrics on
# the merged seurat object. However, if the number of cells exceeds 1.5 million,
# then this step will fail as there are simply too many cells. So, it is better
# to calculate QC metrics on individual samples, perform filtering on individual
# samples and then merge the filtered seurat objects.

# Initialize an empty dataframe where the class of each column resembles those 
# of raw_metadata. We will use this dataframe for making plots later.
raw_metadata <- data.frame(Cell = c(""), Sample = as.factor(1), 
                           nUMIs = c(0), nGenes = c(0), 
                           MitoRatio = c(0), RiboRatio = c(0), Novelty = c(0))

# Calculate QC metrics for each sample individually  
for(i in samples){
  
  sample.seurat <- get(i)
  
  # Compute percent mito percent
  sample.seurat <- Seurat::PercentageFeatureSet(object = sample.seurat,
                                                pattern = dplyr::if_else(species=="Homo sapiens", "MT-", "mt-"),
                                                features = NULL,
                                                col.name = "MitoPercent",
                                                assay = NULL)
  
  # Compute percent ribo percent
  sample.seurat <- Seurat::PercentageFeatureSet(object = sample.seurat,
                                                pattern = dplyr::if_else(species=="Homo sapiens", "^RP[SL]", "^Rp[sl]"),
                                                features = NULL,
                                                col.name = "RiboPercent",
                                                assay = NULL)
  
  # Extract metadata
  sample_metadata <- sample.seurat@meta.data
  
  # Rename columns to be more intuitive and add the additional QC metrics:
  # (i)     Cell      : Unique identifiers corresponding to each cell = barcodes
  # (ii)    Sample    : sample name
  # (iii)   nUMIs     : number of transcripts per cell
  # (iv)    nGenes    : number of genes per cell
  # (v)     nHTO_UMIs : number of HTO reads per cell
  # (vi)    nHTOs     : number of HTOs types per cell
  # (vii)   MitoRatio : MitoPercent/100
  # (viii)	RiboRatio : RiboPercent/100  
  # (ix)    Novelty   : log ratio of genes per UMI
  sample_metadata <- sample_metadata %>% 
    dplyr::transmute(Cell = paste0(orig.ident, "_", rownames(sample_metadata)),
                     Sample = orig.ident,
                     nUMIs = nCount_Spatial,
                     nGenes = nFeature_Spatial,
                     MitoRatio = MitoPercent/100,
                     RiboRatio = RiboPercent/100,
                     Novelty = log10(nGenes)/log10(nUMIs))
  
  # Replace the metadata in raw Seurat object
  sample.seurat@meta.data <- sample_metadata
  
  # Save raw metadata
  raw_metadata <- dplyr::bind_rows(raw_metadata, sample_metadata)
  
  # Assign the seurat object to its corresponding variable
  assign(i, sample.seurat)
}

#******************************STEP 2B: PERFORM QC*****************************#

# Perform QC for each sample individually
for(i in samples){
  
  sample.seurat <- get(i)
  
  gene_cutoff <- 50              # reduced to 50 from 250 used for scRNASeq
  umi_cutoff <- 500
  mito_cutoff <- 0.2
  ribo_cutoff <- 0.05
  novelty_cutoff <- 0.8  	        # use 0.8 as starting point. Maximum 0.9
  
  sample.seurat <- base::subset(x = sample.seurat,
                                subset = ((nGenes >= gene_cutoff) &
                                            (nUMIs >= umi_cutoff) &
                                            (MitoRatio <= mito_cutoff) &
                                            # (RiboRatio >= ribo_cutoff) &
                                            # (!is.na(filtered_seurat@meta.data$Patient)) &
                                            (Novelty >= novelty_cutoff)))
  
  # Assign the seurat object to its corresponding variable
  assign(i, sample.seurat)
}

# Create a merged Seurat object.
# NOTE: Samples will have same barcodes. To keep track of cell identities 
# (i.e.barcodes) coming from each sample after merging, we add a prefix 
# (i.e. sample name) to each barcode using add.cell.ids.
filtered_seurat <- base::merge(x = get(samples[1]),
                               y = lapply(samples[2:length(samples)], get),
                               add.cell.ids = samples,
                               merge.data = FALSE)

#****************************STEP 2C: SAVE THE DATA****************************#

# Create .rds object for filtered seurat object to load at any time
saveRDS(filtered_seurat, file=paste0(seurat_results, "filtered_seurat.rds"))

for(i in samples){
  saveRDS(get(i), file=paste0(seurat_results, i,".rds"))
}

#******************************************************************************#

#*****************STEP 2D: VISUALIZE DATA BEFORE AND AFTER QC******************#

# Remove dummy first row from raw_metadata
raw_metadata <- raw_metadata[-1,]
rownames(raw_metadata) <- raw_metadata$Cell

filtered_metadata <- filtered_seurat@meta.data

# Visualize the number of cell counts per sample
cell_qc <- function(metadata, tag){
  
  ggplot(data = get(metadata), aes(x=Sample, fill=Sample)) + 
    geom_bar() +              
    theme_classic() +         #display with x and y axis lines and no gridlines
    labs(x = "Sample", y = "Number of Cells", title = stringr::str_wrap(paste0("Number of Cells ", tag), 30)) +
    #coord_cartesian(ylim = c(0, 20000)) +
    geom_text(stat="count", aes(label=after_stat(count)), vjust = -1) +
    my_theme
}

# Visualize the number of UMIs per cell
umi_qc <- function(metadata, tag){
  
  ggplot(data = get(metadata), aes(x=Sample, y=nUMIs, fill=Sample)) +
    geom_violin() +        
    theme_classic() +       
    labs(x = "Sample", y = "Number of UMIs", title = stringr::str_wrap(paste0("Distribution of UMIs ", tag),30)) +
    coord_cartesian(ylim = c(100, 1000000)) +
    my_theme +        
    scale_y_log10(breaks = c(100, 1000, 10000, 100000, 1000000)) +  		     #display y axis in log scale
    geom_boxplot(width=0.1) + 
    geom_hline(yintercept = gene_cutoff, linetype = 2)
}

# Visualize the number of genes per cell
gene_qc <- function(metadata, tag){
  
  ggplot(data = get(metadata), aes(x=Sample, y=nGenes, fill = Sample)) +
    geom_violin() +         
    theme_classic() +       
    labs(x = "Sample", y = "Number of Genes", title = stringr::str_wrap(paste0("Distribution of Genes ", tag),30)) +
    coord_cartesian(ylim = c(1, 30000)) +
    my_theme +        
    scale_y_log10() +    	
    geom_boxplot(width=0.1) + 
    geom_hline(yintercept = gene_cutoff, linetype = 2)
}

# Visualize the MitoRatio of each cell
mito_qc <- function(metadata, tag){
  
  ggplot(data = get(metadata), aes(x=Sample, y=MitoRatio, fill = Sample)) +
    geom_violin() +         
    theme_classic() +       
    labs(x = "Sample", y = "MitoRatio", title = stringr::str_wrap(paste0("Distribution of MitoRatio ", tag),30)) +
    coord_cartesian(ylim = c(0, 1)) +
    my_theme +        
    geom_boxplot(width=0.1) + 
    geom_hline(yintercept = mito_cutoff, linetype = 2)
}

# Visualize the RiboRatio of each cell
ribo_qc <- function(metadata, tag){
  
  ggplot(data = get(metadata), aes(x=Sample, y=RiboRatio, fill = Sample)) +
    geom_violin() +         
    theme_classic() +       
    labs(x = "Sample", y = "RiboRatio", title = stringr::str_wrap(paste0("Distribution of RiboRatio ", tag),30)) +
    coord_cartesian(ylim = c(0, 1)) +
    my_theme +        
    geom_boxplot(width=0.1) + 
    geom_hline(yintercept = ribo_cutoff, linetype = 2)
}

# Visualize the novelty or complexity of each cell
novelty_qc <- function(metadata, tag){
  
  ggplot(data = get(metadata), aes(x=Sample, y=Novelty, fill = Sample)) +
    geom_violin() +     
    theme_classic() + 
    labs(x = "Sample", y = "Novelty Score", title = stringr::str_wrap(paste0("Distribution of Novelty Score ", tag),30)) +
    coord_cartesian(ylim = c(0, 1)) +
    my_theme +        
    geom_boxplot(width=0.1) + 
    geom_hline(yintercept = novelty_cutoff, linetype = 2)
}

# Visualize number of genes/cell, number of UMIs/cell & MitoRatio together.
# Bottom left quadrant : Poor quality cells with low genes & UMIs per cell 
# Top right quadrant   : Good quality cells with high genes & UMIs per cell
# Bottom right quadrant: Cells with low genes but high UMIs per cell. These 
# could be dying cells or population of low complexity cells (i.e erythrocytes)
gene_umi_mito_qc <- function(metadata, tag){
  
  ggplot(data = get(metadata), aes(x=nUMIs, y=nGenes, color = MitoRatio)) +
    geom_point() +
    theme_classic() + 
    labs(x = "Number of UMIs", y = "Number of Genes",	 title = paste0("Distribution of UMIs, Genes & MitoRatio ", tag)) +
    my_theme + 
    coord_cartesian(xlim = c(100, 1000000), ylim = c(100, 20000)) +
    scale_x_log10(breaks = c(100, 1000, 10000, 100000, 1000000)) + 
    scale_y_log10() + 
    facet_wrap(.~Sample, nrow = 4) +   #split the plot by X-axis label
    stat_smooth(method=lm, color="yellow") +
    geom_vline(xintercept = umi_cutoff) +    	#draw a vertical line at x=500 i.e.UMIs cutoff
    geom_hline(yintercept = gene_cutoff) +    #draw a horizontal line at y =250 i.e. Genes cutoff
    scale_color_viridis(option = "D", limits = c(0, 1)) 		# limits sets max and min values of gradient 
}

# Plot all QC metrics before and after QC
funcs <- c("cell_qc", "umi_qc", "gene_qc", "mito_qc", "ribo_qc", "novelty_qc",
           "gene_umi_mito_qc")

filenames <- c("Cell_Counts", "UMI_Distribution", "Gene_Distribution",
               "MitoRatio_Distribution", "RiboRatio_Distribution", 
               "Novelty_Score_Distribution", "Genes_UMI_MitoRatio_Distribution")

for (i in 1:length(funcs)){
  
  # Plot QC metrics
  purrr::map2(.x = c("raw_metadata", "filtered_metadata"),
              .y = c("Pre QC", "Post QC"),
              .f = get(funcs[i])) %>% 
    cowplot::plot_grid(plotlist = .,
                       align = "hv",
                       axis = "tblr",
                       nrow = 2,  
                       ncol = 1, 
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
                       #scale = 1,
                       greedy = TRUE,
                       byrow = TRUE,
                       cols = NULL,
                       rows = NULL)  
  
  # Save the plot
  ggplot2::ggsave(filename = paste0("QC_", filenames[i], ".pdf"),
                  plot = last_plot(),
                  device = "pdf",
                  path = seurat_results,
                  #scale = 1,
                  #width = dplyr::if_else(i==7, 17,length(samples)/2),
                  height = dplyr::if_else(i==7, 22, 11),
                  units = c("in"),	 
                  dpi = 600,
                  limitsize = TRUE,
                  bg = NULL)
}

#******************************************************************************#
#                       STEP 3: RUN THE STANDARD PIPELINE                      #
#******************************************************************************#

# Load rds file of seurat objects
for (i in samples){
  sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  sample.seurat <- sctransform_spatial_data(sample.seurat)
  sample.seurat <- cluster_spatial_data(sample.seurat)
  saveRDS(sample.seurat, file=paste0(seurat_results, i,".rds"))
  
  DefaultAssay(sample.seurat) <- "SCT"
  p1 <- DimPlot(sample.seurat, reduction = "umap", label = TRUE)
  p2 <- SpatialDimPlot(sample.seurat, label = TRUE, label.size = 3)
  p <- p1 + p2
  ggsave(filename = paste0(seurat_results, "Clusters_on_slide_", i, ".jpg"), 
         plot = p)
}

# # There is no integration like scRNA Seq. We analyse each slide individually.
# # Unfortunately, all images are stored within each sample. So, we remove 
# # unwanted images from each sample. If more than 1 image is present in each 
# # sample, SpatialDimPlot() will give error.
# integ_data <- sct_data
# for (i in 1:length(sct_data)){
#   integ_data[[i]] <- cluster_data(sct_data[[i]])
#   integ_data[[i]]@images <- integ_data[[i]]@images[names(integ_data[[i]]@images) == names(sct_data)[[i]]]
# }

# Color spots were GOI are present based on expression
for (i in samples){
  sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  DefaultAssay(sample.seurat) <- "SCT"
  Seurat::SpatialFeaturePlot(object = sample.seurat, 
                             features = c("CD8A", "CD8B", "NPEPPS", "CDH12"),
                             ncol = 4,
                             slot = "data")
  ggsave(filename = paste0(seurat_results, "Feature_plot_", i, ".jpg"),
         plot = last_plot(),
         units = c("in"),
         width = 11,
         height = 8)
  
  # This can ONLY plot 2 genes at a time
  SpatialFeaturePlotBlend(sample.seurat, "CD8A", "NPEPPS")
}

# Color spots were CD8 and NPEPPS are expressed in same plot. 
# NOTE: This is NOT expression based. We just color the cells that express our 
# GOI in different colors.
# Use SCT assay to get counts, not Spatial assay as it has too much background.
for (i in samples){
  
  sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  
  # Make sure the genes you want are present in the assay
  GOI <- intersect(c("CD8A", "CD8B"), rownames(sample.seurat@assays$SCT@data))
  cd8_df <- sample.seurat@assays$SCT$counts[GOI, ]
  
  if (length(GOI) > 1){
    cd8_cells <- colnames(cd8_df[,colSums(cd8_df) > 0])
  } else {
    cd8_cells <- names(cd8_df[(cd8_df>0)])
  }
  
  GOI <- intersect(c("NPEPPS"), rownames(sample.seurat@assays$SCT@data))
  npepps_df <- sample.seurat@assays$SCT$counts[c("NPEPPS"), ]
  npepps_cells <- names(npepps_df[(npepps_df>0)])
  
  SpatialPlot(object = sample.seurat, 
              cells.highlight = list(CD8 = cd8_cells, NPEPPS = npepps_cells), 
              cols.highlight =  c("green", "blue", "grey"),
              pt.size.factor = 2)
  
  ggsave(filename = paste0(seurat_results, "Location_plot_", i, ".jpg"),
         plot = last_plot(),
         units = c("in"),
         width = 11,
         height = 8)
}
  
  




# # Classify each spot as CD8+ or CD8-
# #for (i in samples){
# i <- "B8"
#   
#   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
#   #rownames(sample.seurat@meta.data) <- sample.seurat@meta.data$Cell
#   
#   df <- sample.seurat@assays$Spatial$counts[c("CD8A", "CD8B"), ]
#   pos_cells  <- colnames(df[,colSums(df) > 0])
#   neg_cells <- setdiff(colnames(df), pos_cells)
#   pos_cells <- paste0(i,"_", pos_cells)
#   neg_cells <- paste0(i,"_", neg_cells)
#   sample.seurat@meta.data <- sample.seurat@meta.data %>% 
#     dplyr::mutate(CD8_status = dplyr::case_when(Cell %in% pos_cells ~ "CD8_pos",
#                                                 TRUE ~ "CD8_neg"))
# }


i <- "B8"
sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))

a <- GetTissueCoordinates(sample.seurat)
df <- data.frame(sample.seurat@images$B8@coordinates)

### CHECK THAT WE ARE SUCCESSFULLY ABLE TO IDENTIFY NEIGHBORS ###
# Find a cell in the center of the tissue
row_cell <- median(df$row)
col_cell <- median(df$col)
cell <- rownames(df %>% dplyr::filter(row==row_cell, col==col_cell))

# Find adjacent cells
cell1 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell-2)))
cell2 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell+2)))
cell3 <- rownames(df %>% dplyr::filter(row==(row_cell-2), col==col_cell))
cell4 <- rownames(df %>% dplyr::filter(row==(row_cell+2), col==col_cell))
cell5 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell-1)))
cell6 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell+1)))
cell7 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell-1)))
cell8 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell+1)))

SpatialPlot(object = sample.seurat, 
            cells.highlight = list(cell, unlist(lapply(paste0("cell", seq(1,8)), get), use.names=FALSE)), 
            #facet.highlight = TRUE,
            cols.highlight =  c("yellow", "red", "black"), 
            pt.size.factor = 2)

#******************************************************************************#

i <- "B8"
sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
df <- data.frame(sample.seurat@images$B8@coordinates)
expr_df <- sample.seurat@assays$SCT$data[rownames(sample.seurat@assays$SCT$data) %in% c("NPEPPS", "CD8A", "CD8B"),]

# Identify cells that express NPEPPS
npepps <- expr_df["NPEPPS",]
npepps <- npepps[npepps > 0]
npepps <- names(npepps)
# Not all npepps cells have co-ordinates. Remove those that lack co-ordinates.
npepps <- intersect(npepps, rownames(df))

expr_df <- data.frame(expr_df)
# For every cell that expresses NPEPPS, calculate total expression of CD8A and 
# CD8B in 1st neighbor
t_expr <- c()
npepps_expr <- c()
for (i in npepps){
  row_cell <- df[rownames(df) == i,]$row
  col_cell <- df[rownames(df) == i,]$col
  
  # Find adjacent cells
  cell1 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell-2)))
  cell2 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell+2)))
  cell3 <- rownames(df %>% dplyr::filter(row==(row_cell-2), col==col_cell))
  cell4 <- rownames(df %>% dplyr::filter(row==(row_cell+2), col==col_cell))
  cell5 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell-1)))
  cell6 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell+1)))
  cell7 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell-1)))
  cell8 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell+1)))
  
  # SpatialPlot(object = sample.seurat, 
  #             cells.highlight = list(i, unlist(lapply(paste0("cell", seq(1,8)), get), use.names=FALSE)), 
  #             #facet.highlight = TRUE,
  #             cols.highlight =  c("yellow", "red", "black"), 
  #             pt.size.factor = 2)
  
  cells <- intersect(make.names(unlist(lapply(paste0("cell", seq(1,8)), get))), colnames(expr_df))
  print(cells)
  
  t <- expr_df %>% 
    dplyr::select(all_of(cells)) %>%
    tibble::rownames_to_column("Gene") %>%
    dplyr::filter(Gene != "NPEPPS") %>%
    tibble::column_to_rownames("Gene")
  
  n <- expr_df %>% 
    dplyr::select(all_of(cells)) %>%
    tibble::rownames_to_column("Gene") %>%
    dplyr::filter(Gene == "NPEPPS") %>%
    tibble::column_to_rownames("Gene")
  
  t_expr <- c(t_expr, sum(t))
  npepps_expr <- c(npepps_expr, sum(n))
}

df <- data.frame(npepps_expr, t_expr)
# Save batch corrected normalized counts for entire dataset
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "spatial")
openxlsx::writeData(wb, sheet = "spatial", x = df, rowNames = FALSE)
openxlsx::saveWorkbook(wb,
                       file = paste0(seurat_results, "Correlation.xlsx"),
                       overwrite = TRUE)

















# pseudobulk the counts based on donor-condition-celltype
pseudo_ifnb <- AggregateExpression(ifnb, assays = "RNA", return.seurat = T, group.by = c("stim", "donor_id", "seurat_annotations"))

# each 'cell' is a donor-condition-celltype pseudobulk profile
tail(Cells(pseudo_ifnb))




































# Identification of Spatially Variable Features

# Seurat offers two workflows to identify molecular features that correlate 
# with spatial location within a tissue. The first is to perform differential 
# expression based on pre-annotated anatomical regions within the tissue, which
# may be determined either from unsupervised clustering or prior knowledge. 
# This strategy works will in this case, as the clusters above exhibit clear 
# spatial restriction.

de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = brain, 
                   features = rownames(de_markers)[1:3], 
                   alpha = c(0.1, 1), 
                   ncol = 3)


# An alternative approach, implemented in FindSpatiallyVariables(), is to 
# search for features exhibiting spatial patterning in the absence of 
# pre-annotation. The default method (method = 'markvariogram), is inspired by 
# the Trendsceek, which models spatial transcriptomics data as a mark point 
# process and computes a ‘variogram’, which identifies genes whose expression 
# level is dependent on their spatial location. More specifically, this process
# calculates gamma(r) values measuring the dependence between two spots a 
# certain “r” distance apart. By default, we use an r-value of ‘5’ in these 
# analyses, and only compute these values for variable genes (where variation 
# is calculated independently of spatial location) to save time.


brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000],
                                       selection.method = "moransi")
top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "moransi"), 6)
SpatialFeaturePlot(brain, 
                   features = top.features, 
                   alpha = c(0.1, 1),
                   ncol = 3)



# Integration with single-cell data i.e. label transfer from single-cell data
# NOTE: At ~50um size, spots from the visium assay will encompass the expression 
# profiles of multiple cells. 

# Users may be interested to ‘deconvolute’ each of the spatial voxels to predict
# the underlying composition of cell types. We tested a wide variety of 
# deconvolution and integration methods, using a reference scRNA-seq dataset of
# ~14,000 adult mouse cortical cell taxonomy from the Allen Institute, 
# generated with the SMART-Seq2 protocol. 

# We consistently found superior performance using integration methods (as 
# opposed to deconvolution methods), likely because of substantially different 
# noise models that characterize spatial and single-cell datasets, and 
# integration methods are specifically designed to be robust to 
# these differences. 

# We therefore apply the ‘anchor’-based integration workflow introduced in 
# Seurat v3, that enables the probabilistic transfer of annotations from a 
# reference to a query set. 

# While many of the methods are conserved (both procedures begin by identifying 
# anchors), there are two important distinctions between data transfer and 
# integration:
# (i) In data transfer, Seurat does not correct or modify query expression data.
# (ii) In data transfer, Seurat has an option (set by default) to project the 
# PCA structure of a reference onto the query, instead of learning a joint 
# structure with CCA. We generally suggest using this option when projecting 
# data between scRNA-seq datasets.

# After finding anchors, we use the TransferData() function to classify the 
# query cells based on reference data (a vector of reference cell type labels). 
# TransferData() returns a matrix with predicted IDs and prediction scores, 
# which we can add to the query metadata.

# NOTE: Make sure reference seurat object is SCTransformed.
integrated_seurat <- readRDS("/hpc/home/kailasamms/scratch/scRNASeq_Chen/results_seurat/integrated_seurat_snn.rds")

# Seurat v3 vs Seurat v5 issues: Check if cells in graph are in same order
all(Cells(integrated_seurat@graphs$integrated_snn) == colnames(integrated_seurat))
# if FALSE, reorder the cells in all existing graphs
integrated_seurat@graphs$integrated_snn <- 
  integrated_seurat@graphs$integrated_snn[colnames(integrated_seurat), colnames(integrated_seurat)]
integrated_seurat@graphs$integrated_nn <- 
  integrated_seurat@graphs$integrated_nn[colnames(integrated_seurat), colnames(integrated_seurat)]

# Remove ambiguous cells before label transfer
integrated_seurat <- subset(integrated_seurat, 
                            sub_type == "Unclassified", 
                            invert=TRUE)

# NOTE: For FindTransferAnchors() to work, both reference and query MUST be
# SCtransformed()  https://github.com/satijalab/seurat/issues/3937
# So, run SCTransform() on RNA assay.
DefaultAssay(integrated_seurat) <- "integrated"

# Make sure all assays in Seuratv3 object has same order of cells
acells <- colnames(x = integrated_seurat[["integrated"]])
ocells <- colnames(x = integrated_seurat)

integrated_seurat@assays$RNA@counts <- integrated_seurat@assays$RNA@counts[,ocells]
integrated_seurat@assays$RNA@data <- integrated_seurat@assays$RNA@data[,ocells]
integrated_seurat@assays$integrated@data <- integrated_seurat@assays$integrated@data[,ocells]
integrated_seurat@assays$integrated@scale.data <- integrated_seurat@assays$integrated@scale.data[,ocells]
integrated_seurat@assays$SCT@counts <- integrated_seurat@assays$SCT@counts[,ocells]
integrated_seurat@assays$SCT@data <- integrated_seurat@assays$SCT@data[,ocells]
integrated_seurat@assays$SCT@scale.data <- integrated_seurat@assays$SCT@scale.data[,ocells]

# NOTE: If you subset the spatial object, re-run SCTranform() on it.
for (i in samples){
  sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  DefaultAssay(sample.seurat) <- "SCT"
  
  # Run ONLY if you had subset the spatial data
  # sample.seurat <- SCTransform(sample.seurat, 
  #                              assay = "Spatial",
  #                              verbose = FALSE)
  # sample.seurat <- RunPCA(sample.seurat,
  #                         verbose = FALSE)
  
  # Find anchors between reference and query
  anchors <- FindTransferAnchors(reference = integrated_seurat,
                                 query = sample.seurat,
                                 normalization.method = "SCT",
                                 reference.assay = "integrated",
                                 query.assay = "SCT")
  
  predictions.assay <- TransferData(anchorset = anchors, 
                                    refdata = integrated_seurat$cell_type,
                                    prediction.assay = TRUE,
                                    weight.reduction = sample.seurat[["pca"]],
                                    dims = 1:30)
  sample.seurat[["predictions"]] <- predictions.assay
  
  DefaultAssay(sample.seurat) <- "predictions"
  
  SpatialFeaturePlot(sample.seurat, 
                     features = rownames(predictions.assay@data), 
                     pt.size.factor = 1.6,
                     #ncol = 4,
                     crop = FALSE)
  
  ggsave(filename = paste0(seurat_results, "Labeltransfer_plot_", i, ".jpg"),
         plot = last_plot(),
         units = c("in"),
         width = 11,
         height = 8)
  
  
  
  # sample.seurat <- FindSpatiallyVariableFeatures(sample.seurat,
  #                                                assay = "predictions",
  #                                                selection.method = "moransi",
  #                                                features = rownames(sample.seurat),
  #                                                r.metric = 5,
  #                                                slot = "data")
  
  # top.clusters <- head(SpatiallyVariableFeatures(sample.seurat, 
  #                                                selection.method = "moransi"), 4)
  # SpatialPlot(object = sample.seurat, 
  #             features = top.clusters,
  #             ncol = 2)
}

# NOTE: In single cell experiment, we combine multiple samples into a single 
# seurat object and do analysis. In spatial experiment, we analyze each sample
# individually.

# NOTE: Reads are obtained from both background (non-tissue) and tissue specific
# areas and stored in raw_feature_bc_matrix and filtered_fearure_bc_matrix
# respectively. Since, we are ONLY interested in reads obtained from tissues,
# use ONLY filtered_feature_bc_matrix.

# NOTE: Size of smallest mammalian cell is ~4um. 

# NOTE: In visium v1 spatial transcriptomics, there are 4992 barcoded spots 
# which capture RNA from tissue placed above these spots. 
# So, unlike scRNASeq, each barcode corresponds to a spot, not a cell.
# https://www.youtube.com/watch?v=VwNk4d-0RJc
# https://kb.10xgenomics.com/hc/en-us/articles/360035848191-How-many-spots-are-within-a-single-capture-area-on-the-Visium-v1-Spatial-Gene-Expression-Slide

# DIRECTORY STRUCTURE IS IMPORTANT:
# -> filt_feature_bc_matrix  (dir)
#   -> TMA1-A1 (dir)
#     -> binned_outputs (dir)
#       -> square_002um (dir)      
#         -> filtered_feature_bc_matrix.h5
#         -> spatial (dir)
#           -> aligned_fiducials.jpg
#           -> cytassist_image.tiff
#           -> detected_tissue_image.jpg
#           -> scalefactors_json.json
#           -> tissue_hires_image.png
#           -> tissue_lowres_image.png
#           -> tissue_positions.parquet
#       -> square_008um (dir)
#         -> filtered_feature_bc_matrix.h5
#         -> spatial (dir)
#       -> square_016um (dir)
#         -> filtered_feature_bc_matrix.h5
#         -> spatial (dir)

# data.dir MUST have 
# (i) a H5 file specified by filename parameter
# (ii) a folder named "spatial" containing the image

#**********************STEP 2A: CALCULATE ALL QC METRICS***********************#

# NOTE: PercentageFeatureSet() calculates metrics for cells ONLY within the 
# default assay, if no assay is indicated. So, we loop through each assay and
# calculate metrics

# In single cell analysis, we have to perform QC to remove bad cells in order to
# accurately annotate each cell. In spatial analysis, RNA from tissue is 
# captured within the spot below the tissue. Since sequencing depth is usually 
# way lower than single cell, we use almost all reads (lenient cutoffs) captured
# in the spot to annotate the spots.

# NOTE: Unlike PercentageFeatureSet(), subset() uses all the cells irrespective 
# of the default assay. 2um bins are too small to contain 100 UMIs or 50 genes 
# and the entire Spatial.002um assay will be removed. So, we use lenient cutoffs

# STEP 3: RUN THE STANDARD PIPELINE                      #

# Create workbook to save markers
wb <- openxlsx::createWorkbook()
for (i in samples){
  
  object <- readRDS(paste0(seurat_results, i, ".filt.rds"))
  
  # Run workflow on all assays
  for (assay in c("Spatial.008um", "Spatial.016um")){
    
    DefaultAssay(object) <- assay
    object <- Seurat::NormalizeData(       object = object, assay = assay, normalization.method = "LogNormalize")
    object <- Seurat::FindVariableFeatures(object = object, assay = assay, nfeatures = 2000)
    object <- Seurat::ScaleData(           object = object, assay = assay, features = VariableFeatures(object))
    
    # Create a new 'sketch' assay using 50k cells
    sample_n <- 50000
    e <- "Error"
    while (class(e) == "character" & sample_n > 0){
      e <- tryCatch(Seurat::SketchData(object = object, 
                                       assay = assay,
                                       ncells = sample_n,
                                       method = "LeverageScore",
                                       features = VariableFeatures(object),
                                       sketched.assay = paste0(assay, ".sketch")), error = function(msg){
                                         print("Reducing number of cells sampled")
                                         return("Error")})
      sample_n <- sample_n-5000
      cat(assay, ":", sample_n, "\n")
    }
    
    object <- e
    
    # Switch analysis to sketched cells
    DefaultAssay(object) <- paste0(assay, ".sketch")
    object <- Seurat::FindVariableFeatures(object = object, 
                                           assay = paste0(assay, ".sketch"), 
                                           nfeatures = 2000)
    
    object <- Seurat::ScaleData(object = object, 
                                assay = paste0(assay, ".sketch"), 
                                features = VariableFeatures(object))
    
    object <- Seurat::RunPCA(object = object, 
                             assay = paste0(assay, ".sketch"), 
                             reduction.name = paste0(assay, ".pca.sketch"))
    
    object <- Seurat::FindNeighbors(object = object, 
                                    assay = paste0(assay, ".sketch"), 
                                    reduction = paste0(assay, ".pca.sketch"), 
                                    dims = 1:50)
    
    object <- Seurat::FindClusters(object = object, 
                                   algorithm = 4,
                                   cluster.name = paste0(assay, ".cluster.sketch"), 
                                   resolution = 0.6)
    
    object <- Seurat::RunUMAP(object = object, 
                              reduction.name = paste0(assay, ".umap.sketch"), 
                              reduction = paste0(assay, ".pca.sketch"), 
                              dims = 1:50)
    
    # Project the cluster labels, dimensional reductions (PCA and UMAP) that we
    # learned from the 50,000 sketched cells to the entire dataset
    object <- Seurat::ProjectData(object = object, 
                                  assay = assay,
                                  sketched.assay = paste0(assay, ".sketch"),
                                  sketched.reduction = paste0(assay, ".pca.sketch"),
                                  full.reduction = paste0(assay, ".pca.full"),
                                  dims = 1:50,
                                  refdata = list(seurat_cluster.projected = paste0(assay, ".cluster.sketch")),
                                  umap.model = paste0(assay, ".umap.sketch"))
    
    object <- Seurat::RunUMAP(object = object, 
                              reduction.name = paste0(assay, ".umap.full"), 
                              reduction = paste0(assay, ".pca.full"), 
                              dims = 1:50)
    
    # Append assay name to the newly generated metadata columns so they dont get overwritten
    object@meta.data <- object@meta.data %>%
      dplyr::rename(!!rlang::sym(paste0(assay, ".cluster.full")) := "seurat_cluster.projected",
                    !!rlang::sym(paste0(assay, ".cluster.full.score")) := "seurat_cluster.projected.score")
    
    # We can visualize the clustering results for the sketched cells, as well as the
    # projected clustering results for the full dataset:
    DefaultAssay(object) <- paste0(assay, ".sketch")
    Idents(object) <- paste0(assay, ".cluster.sketch")
    p1 <- DimPlot(object, reduction = paste0(assay, ".umap.sketch"), label = T) + 
      ggtitle("Sketched clustering") + 
      theme(legend.position = "none")
    
    # switch to full dataset
    DefaultAssay(object) <- assay
    Idents(object) <- paste0(assay, ".cluster.full")
    p2 <- DimPlot(object, reduction = paste0(assay, ".umap.full"), label = T) + 
      ggtitle("Projected clustering (full dataset)") + 
      theme(legend.position = "none")
    
    # Save the plot
    ggplot2::ggsave(filename=paste0(i, assay, ".umap.tiff"),
                    plot=p1+p2,
                    device="jpeg",
                    path=diagnostics_path,
                    width=11,
                    height=8.5,
                    units=c("in"),
                    dpi=600,
                    limitsize=TRUE,
                    bg="white")
    
    # Find markers
    markers <- Seurat::FindAllMarkers(object = object, 
                                      assay = assay,
                                      only.pos = TRUE)
    
    top10 <- markers %>%
      dplyr::mutate(cluster = as.numeric(as.character(cluster))) %>%
      dplyr::group_by(cluster) %>%
      dplyr::filter(avg_log2FC > 1, p_val_adj < 0.05) %>%
      slice_max(order_by=avg_log2FC*pct.1, n = 10) %>%
      ungroup() %>%
      data.frame() %>%
      dplyr::arrange(cluster, desc(avg_log2FC*pct.1))
    
    # Save markers to worksheet
    openxlsx::addWorksheet(wb, sheetName = paste0(i, assay))
    openxlsx::writeData(wb, sheet = paste0(i, assay), x = markers, rowNames = FALSE)
    openxlsx::addWorksheet(wb, sheetName = paste0(i, assay, "top10"))
    openxlsx::writeData(wb, sheet = paste0(i, assay, "top10"), x = top10, rowNames = FALSE)
  }
  # Save xlsx file with markers
  openxlsx::saveWorkbook(wb, file = paste0(seurat_results, i, "_Markers.xlsx"), overwrite = TRUE)
  
  # Save the object
  saveRDS(object, file=paste0(seurat_results, i,".rds"))
}
  


#******************REMOVE BAD SAMPLES

for(i in samples){
  
  sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
 
  if (i == "TMA1-A1"){
    sample.seurat <- subset(sample.seurat, 
                            Sample %in% c("TMA1.D1.A10", "TMA1.D1.A11","TMA1.D1.A12",
                                          "TMA1.D1.B10", "TMA1.D1.C09", 
                                          "TMA1.D1.D09", "TMA1.D1.D12"))
  }
  
  if (i == "TMA1-D1"){
    sample.seurat <- subset(sample.seurat, 
                            Sample %in% c("TMA1.D1.A10", "TMA1.D1.A11","TMA1.D1.A12",
                                          "TMA1.D1.B10", "TMA1.D1.C09", 
                                          "TMA1.D1.D09", "TMA1.D1.D12"))

  }
}

#******************FIND CELL TYPES BASED ON MARKERS

df <- openxlsx::read.xlsx("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Markers.All.TMA1-A1.Spatial.016um.cluster.0.4.xlsx") %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::distinct(cluster, gene, avg_log2FC, pct.1, pct.2, .keep_all = TRUE) %>%
  tidyr::pivot_wider(id_cols=cluster, names_from=gene, values_from=avg_log2FC, values_fill=0.0) %>%
  tibble::column_to_rownames("cluster") %>%
  scale()

row_dist <- stats::dist(x=df, method = "euclidean", diag = TRUE, upper = TRUE)
row_clust <- hclust(row_dist)
row_order <- rownames(df[row_clust$order,])
col_dist <- stats::dist(x=t(df), method = "euclidean", diag = TRUE, upper = TRUE)
col_clust <- hclust(col_dist)
col_order <- colnames(df[,col_clust$order])
mat_A1 <- df[row_order, col_order]

df <- openxlsx::read.xlsx("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Markers.All.TMA1-A1.Spatial.008um.cluster.0.4.xlsx") %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::distinct(cluster, gene, avg_log2FC, pct.1, pct.2, .keep_all = TRUE) %>%
  tidyr::pivot_wider(id_cols=cluster, names_from=gene, values_from=avg_log2FC, values_fill=0) %>%
  tibble::column_to_rownames("cluster") %>%
  scale()

row_dist <- stats::dist(x=df, method = "euclidean", diag = TRUE, upper = TRUE)
row_clust <- hclust(row_dist)
row_order <- rownames(df[row_clust$order,])
col_dist <- stats::dist(x=t(df), method = "euclidean", diag = TRUE, upper = TRUE)
col_clust <- hclust(col_dist)
col_order <- colnames(df[,col_clust$order])
mat_D1 <- df[row_order, col_order]

mat_A1[mat_A1 < 0] <- 0 
mat_D1[mat_D1 < 0] <- 0 

mat_A1 <- mat_A1 %>% 
  data.frame() %>% 
  mutate(across(where(is.numeric), ~round(., 2)))

mat_D1 <- mat_D1 %>% 
  data.frame() %>% 
  mutate(across(where(is.numeric), ~round(., 2)))

# Create workbook to save markers
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "A1")
openxlsx::writeData(wb, sheet = "A1", x = t(mat_A1), rowNames = TRUE)
openxlsx::addWorksheet(wb, sheetName = "D1")
openxlsx::writeData(wb, sheet = "D1", x = t(mat_D1), rowNames = TRUE)
openxlsx::saveWorkbook(wb, file ="C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/TMA-A1.marker.new.xlsx")
openxlsx::saveWorkbook(wb, file = paste0(seurat_results, i, "_Markers.xlsx"), overwrite = TRUE)

# my_palette <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)
# # Define breaks
# if(max(mat) == 0){
#   breaks <- c(seq(from = floor(min(mat)), to = 0, length.out = 100))
#   my_palette <- my_palette[1:50]
# } else if (min(mat) == 0){
#   breaks <- c(seq(from = 0, to = ceiling(max(mat)), length.out = 100))
#   my_palette <- my_palette[50:100]
# } else if(min(mat) < -3 | max(mat) > 3){
#   breaks <- c(seq(-3, 0, length.out = 50), seq(3/100, 3, length.out = 50))
# } else{
#   breaks <- c(seq(from = floor(min(mat)), to = 0, length.out = 50), seq(from = max(mat)/100, to = ceiling(max(mat)), length.out = 50))
# }
# pheatmap::pheatmap(df, breaks=breaks)

#******************ANNOTATE 8um BINS

Pancreatic.Exocrine <- c("AMY2A", "CEL", "CELA2A", "CELA2B", "CELA3A", "CELA3B",
                         "CLPS", "CPA1", "CPA2", "CPB1", "CTRB1", "CTRB2", 
                         "CTRC", "CTRL", "CUZD1", "ERP27", "GP2", "KLK1", 
                         "PDIA2", "PLA2G1B", "PNLIP", "PNLIPRP1", "PRSS1", 
                         "PRSS3", "SERPINI2", "SYCN")
Hepatocytes <- c("AHSG", "ALB", "ALDOB", "AMBP", "ANGPTL3", "AOX1", "APCS", 
                 "APOA1", "APOA2", "APOA5", "APOC2", "APOC3", "APOH", "ARG1", 
                 "BAAT", "C3", "C4BPB", "C5", "C9", "CP", "CPB2", "CPS1", "F2",
                 "FGA", "FGB", "FGG", "FGL1", "GC", "HAMP", "HAO1", "HP", "HPR",
                 "HPX", "HRG", "ITIH3", "MAT1A", "ORM1", "ORM2", "PAH", 
                 "SERPINC1", "SLC25A47", "SLC2A2", "SULT2A1", "VTN")

for (f in c("Pancreatic.Exocrine", "Hepatocytes")){
  features <- get(f)
  features <- rownames(sample.seurat@assays$Spatial.008um$data)[tolower(rownames(sample.seurat@assays$Spatial.008um$data)) %in% 
                                                                  tolower(features)]
  features <- list(sort(features))
  
  # Calculate module scores
  sample.seurat <- Seurat::AddModuleScore(object=sample.seurat,
                                          features=features,
                                          assay="Spatial.008um",
                                          slot="data",
                                          name=f)
}
DefaultAssay(sample.seurat) <- "Spatial.008um"
SpatialFeaturePlot(object = sample.seurat,
                   features="Pancreatic.Exocrine1",
                   image.scale="lowres",   # "hires" hides show H&E image
                   pt.size.factor = 4)     # default value 1.6 gives weird shapes. Values greater than 4 makes dot large

SpatialFeaturePlot(object = sample.seurat,
                   features="Hepatocytes1",
                   pt.size.factor = 4)  

  # DefaultAssay(object) <- "Spatial.008um"
  # 
  # vln.plot <- VlnPlot(object, features = "nCount_Spatial.008um", pt.size = 0) + 
  #   theme(axis.text = element_text(size = 4)) + 
  #   NoLegend()
  # count.plot <- SpatialFeaturePlot(object, features = "nCount_Spatial.008um") + 
  #   theme(legend.position = "right")
  
  # # switch back to 8um
  # #object@assays$Spatial.008um@features@dimnames[[1]][25:30]
  # DefaultAssay(object) <- "Spatial.008um"
  # p2 <- SpatialFeaturePlot(object, features = "GAPDH") + ggtitle("Hpca expression (8um)")
  
  
  #   # Visualize the unsupervised clusters based on their spatial location. 
  #   SpatialDimPlot(object, label = T, repel = T, label.size = 4)
  #   
  #   # Plot the spatial location of different clusters individually.
  #   Idents(object) <- "seurat_cluster.projected"
  #   cells <- CellsByIdentities(object, idents = c(0, 4, 32, 34, 35))
  #   p <- SpatialDimPlot(object,
  #                       cells.highlight = cells[setdiff(names(cells), "NA")],
  #                       cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
  #   ) + NoLegend()
  #   p
  #   
  #   
  #   # visualize the top gene expression markers for each cluster:
  #   # Create downsampled object to make visualization either
  #   DefaultAssay(object) <- "Spatial.008um"
  #   Idents(object) <- "seurat_cluster.projected"
  #   object_subset <- subset(object, cells = Cells(object[["Spatial.008um"]]), downsample = 1000)
  #   
  #   
  #   
  #   object_subset <- ScaleData(object_subset, assay = "Spatial.008um", features = top5$gene)
  #   p <- DoHeatmap(object_subset, assay = "Spatial.008um", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
  #   p
  #   
  #   # Identifying spatially-defined tissue domains
  #   object <- Banksy::RunBanksy(object,
  #                               lambda = 0.8, verbose = TRUE,
  #                               assay = "Spatial.008um", slot = "data", features = "variable",
  #                               k_geom = 50)
  #   
  #   DefaultAssay(object) <- "BANKSY"
  #   object <- RunPCA(object, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(object), npcs = 30)
  #   object <- FindNeighbors(object, reduction = "pca.banksy", dims = 1:30)
  #   object <- FindClusters(object, cluster.name = "banksy_cluster", resolution = 0.5)
  #   
  #   Idents(object) <- "banksy_cluster"
  #   p <- SpatialDimPlot(object, group.by = "banksy_cluster", label = T, repel = T, label.size = 4)
  #   p
  #   
  #   # highlight the spatial location of each tissue domain individually
  #   banksy_cells <- CellsByIdentities(object)
  #   p <- SpatialDimPlot(object, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T) + NoLegend()
  #   p
  #   
  # }
  # 

  # 
  # # Load rds file of seurat objects
  # for (i in samples){
  #   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  #   sample.seurat <- sctransform_spatial_data(sample.seurat)
  #   sample.seurat <- cluster_spatial_data(sample.seurat)
  #   saveRDS(sample.seurat, file=paste0(seurat_results, i,".rds"))
  #   
  #   DefaultAssay(sample.seurat) <- "SCT"
  #   p1 <- DimPlot(sample.seurat, reduction = "umap", label = TRUE)
  #   p2 <- SpatialDimPlot(sample.seurat, label = TRUE, label.size = 3)
  #   p <- p1 + p2
  #   ggsave(filename = paste0(seurat_results, "Clusters_on_slide_", i, ".jpg"), 
  #          plot = p)
  # }
  # 
  # # # There is no integration like scRNA Seq. We analyse each slide individually.
  # # # Unfortunately, all images are stored within each sample. So, we remove 
  # # # unwanted images from each sample. If more than 1 image is present in each 
  # # # sample, SpatialDimPlot() will give error.
  # # integ_data <- sct_data
  # # for (i in 1:length(sct_data)){
  # #   integ_data[[i]] <- cluster_data(sct_data[[i]])
  # #   integ_data[[i]]@images <- integ_data[[i]]@images[names(integ_data[[i]]@images) == names(sct_data)[[i]]]
  # # }
  # 
  # # Color spots were GOI are present based on expression
  # for (i in samples){
  #   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  #   DefaultAssay(sample.seurat) <- "SCT"
  #   Seurat::SpatialFeaturePlot(object = sample.seurat, 
  #                              features = c("CD8A", "CD8B", "NPEPPS", "CDH12"),
  #                              ncol = 4,
  #                              slot = "data")
  #   ggsave(filename = paste0(seurat_results, "Feature_plot_", i, ".jpg"),
  #          plot = last_plot(),
  #          units = c("in"),
  #          width = 11,
  #          height = 8)
  #   
  #   # This can ONLY plot 2 genes at a time
  #   SpatialFeaturePlotBlend(sample.seurat, "CD8A", "NPEPPS")
  # }
  # 
  # # Color spots were CD8 and NPEPPS are expressed in same plot. 
  # # NOTE: This is NOT expression based. We just color the cells that express our 
  # # GOI in different colors.
  # # Use SCT assay to get counts, not Spatial assay as it has too much background.
  # for (i in samples){
  #   
  #   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  #   
  #   # Make sure the genes you want are present in the assay
  #   GOI <- intersect(c("CD8A", "CD8B"), rownames(sample.seurat@assays$SCT@data))
  #   cd8_df <- sample.seurat@assays$SCT$counts[GOI, ]
  #   
  #   if (length(GOI) > 1){
  #     cd8_cells <- colnames(cd8_df[,colSums(cd8_df) > 0])
  #   } else {
  #     cd8_cells <- names(cd8_df[(cd8_df>0)])
  #   }
  #   
  #   GOI <- intersect(c("NPEPPS"), rownames(sample.seurat@assays$SCT@data))
  #   npepps_df <- sample.seurat@assays$SCT$counts[c("NPEPPS"), ]
  #   npepps_cells <- names(npepps_df[(npepps_df>0)])
  #   
  #   SpatialPlot(object = sample.seurat, 
  #               cells.highlight = list(CD8 = cd8_cells, NPEPPS = npepps_cells), 
  #               cols.highlight =  c("green", "blue", "grey"),
  #               pt.size.factor = 2)
  #   
  #   ggsave(filename = paste0(seurat_results, "Location_plot_", i, ".jpg"),
  #          plot = last_plot(),
  #          units = c("in"),
  #          width = 11,
  #          height = 8)
  # }
  # 
  # # # Classify each spot as CD8+ or CD8-
  # # #for (i in samples){
  # # i <- "B8"
  # #   
  # #   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  # #   #rownames(sample.seurat@meta.data) <- sample.seurat@meta.data$Cell
  # #   
  # #   df <- sample.seurat@assays$Spatial$counts[c("CD8A", "CD8B"), ]
  # #   pos_cells  <- colnames(df[,colSums(df) > 0])
  # #   neg_cells <- setdiff(colnames(df), pos_cells)
  # #   pos_cells <- paste0(i,"_", pos_cells)
  # #   neg_cells <- paste0(i,"_", neg_cells)
  # #   sample.seurat@meta.data <- sample.seurat@meta.data %>% 
  # #     dplyr::mutate(CD8_status = dplyr::case_when(Cell %in% pos_cells ~ "CD8_pos",
  # #                                                 TRUE ~ "CD8_neg"))
  # # }
  # 
  # 
  # i <- "B8"
  # sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  # 
  a <- GetTissueCoordinates(sample.seurat)
  df <- data.frame(sample.seurat@images$B8@coordinates)

  ### CHECK THAT WE ARE SUCCESSFULLY ABLE TO IDENTIFY NEIGHBORS ###
  # Find a cell in the center of the tissue
  row_cell <- median(df$x)
  col_cell <- median(df$y)
  cell <- rownames(df %>% dplyr::filter(x==row_cell, y==col_cell))

  # # Find adjacent cells
  # cell1 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell-2)))
  # cell2 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell+2)))
  # cell3 <- rownames(df %>% dplyr::filter(row==(row_cell-2), col==col_cell))
  # cell4 <- rownames(df %>% dplyr::filter(row==(row_cell+2), col==col_cell))
  # cell5 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell-1)))
  # cell6 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell+1)))
  # cell7 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell-1)))
  # cell8 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell+1)))
  # 
  SpatialPlot(object = sample.seurat,
              cells.highlight = list(cell, unlist(lapply(paste0("cell", seq(1,8)), get), use.names=FALSE)),
              #facet.highlight = TRUE,
              cols.highlight =  c("yellow", "red", "black"),
              pt.size.factor = 2)
  
  SingleImagePlot(data = sample.seurat,
  # 
  # #******************************************************************************#
  # 
  # i <- "B8"
  # sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  # df <- data.frame(sample.seurat@images$B8@coordinates)
  # expr_df <- sample.seurat@assays$SCT$data[rownames(sample.seurat@assays$SCT$data) %in% c("NPEPPS", "CD8A", "CD8B"),]
  # 
  # # Identify cells that express NPEPPS
  # npepps <- expr_df["NPEPPS",]
  # npepps <- npepps[npepps > 0]
  # npepps <- names(npepps)
  # # Not all npepps cells have co-ordinates. Remove those that lack co-ordinates.
  # npepps <- intersect(npepps, rownames(df))
  # 
  # expr_df <- data.frame(expr_df)
  # # For every cell that expresses NPEPPS, calculate total expression of CD8A and 
  # # CD8B in 1st neighbor
  # t_expr <- c()
  # npepps_expr <- c()
  # for (i in npepps){
  #   row_cell <- df[rownames(df) == i,]$row
  #   col_cell <- df[rownames(df) == i,]$col
  #   
  #   # Find adjacent cells
  #   cell1 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell-2)))
  #   cell2 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell+2)))
  #   cell3 <- rownames(df %>% dplyr::filter(row==(row_cell-2), col==col_cell))
  #   cell4 <- rownames(df %>% dplyr::filter(row==(row_cell+2), col==col_cell))
  #   cell5 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell-1)))
  #   cell6 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell+1)))
  #   cell7 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell-1)))
  #   cell8 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell+1)))
  #   
  #   # SpatialPlot(object = sample.seurat, 
  #   #             cells.highlight = list(i, unlist(lapply(paste0("cell", seq(1,8)), get), use.names=FALSE)), 
  #   #             #facet.highlight = TRUE,
  #   #             cols.highlight =  c("yellow", "red", "black"), 
  #   #             pt.size.factor = 2)
  #   
  #   cells <- intersect(make.names(unlist(lapply(paste0("cell", seq(1,8)), get))), colnames(expr_df))
  #   print(cells)
  #   
  #   t <- expr_df %>% 
  #     dplyr::select(all_of(cells)) %>%
  #     tibble::rownames_to_column("Gene") %>%
  #     dplyr::filter(Gene != "NPEPPS") %>%
  #     tibble::column_to_rownames("Gene")
  #   
  #   n <- expr_df %>% 
  #     dplyr::select(all_of(cells)) %>%
  #     tibble::rownames_to_column("Gene") %>%
  #     dplyr::filter(Gene == "NPEPPS") %>%
  #     tibble::column_to_rownames("Gene")
  #   
  #   t_expr <- c(t_expr, sum(t))
  #   npepps_expr <- c(npepps_expr, sum(n))
  # }
  # 
  # df <- data.frame(npepps_expr, t_expr)
  # # Save batch corrected normalized counts for entire dataset
  # wb <- openxlsx::createWorkbook()
  # openxlsx::addWorksheet(wb, sheetName = "spatial")
  # openxlsx::writeData(wb, sheet = "spatial", x = df, rowNames = FALSE)
  # openxlsx::saveWorkbook(wb,
  #                        file = paste0(seurat_results, "Correlation.xlsx"),
  #                        overwrite = TRUE)
  # 

  # 
  # # pseudobulk the counts based on donor-condition-celltype
  # pseudo_ifnb <- AggregateExpression(ifnb, assays = "RNA", return.seurat = T, group.by = c("stim", "donor_id", "seurat_annotations"))
  # 
  # # each 'cell' is a donor-condition-celltype pseudobulk profile
  # tail(Cells(pseudo_ifnb))
  # 
  # # Identification of Spatially Variable Features
  # 
  # # Seurat offers two workflows to identify molecular features that correlate 
  # # with spatial location within a tissue. The first is to perform differential 
  # # expression based on pre-annotated anatomical regions within the tissue, which
  # # may be determined either from unsupervised clustering or prior knowledge. 
  # # This strategy works will in this case, as the clusters above exhibit clear 
  # # spatial restriction.
  # 
  # de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
  # SpatialFeaturePlot(object = brain, 
  #                    features = rownames(de_markers)[1:3], 
  #                    alpha = c(0.1, 1), 
  #                    ncol = 3)
  # 
  # 
  # # An alternative approach, implemented in FindSpatiallyVariables(), is to 
  # # search for features exhibiting spatial patterning in the absence of 
  # # pre-annotation. The default method (method = 'markvariogram), is inspired by 
  # # the Trendsceek, which models spatial transcriptomics data as a mark point 
  # # process and computes a ‘variogram’, which identifies genes whose expression 
  # # level is dependent on their spatial location. More specifically, this process
  # # calculates gamma(r) values measuring the dependence between two spots a 
  # # certain “r” distance apart. By default, we use an r-value of ‘5’ in these 
  # # analyses, and only compute these values for variable genes (where variation 
  # # is calculated independently of spatial location) to save time.
  # 
  # 
  # brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000],
  #                                        selection.method = "moransi")
  # top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "moransi"), 6)
  # SpatialFeaturePlot(brain, 
  #                    features = top.features, 
  #                    alpha = c(0.1, 1),
  #                    ncol = 3)
  # 
  # 
  # 
  # # Integration with single-cell data i.e. label transfer from single-cell data
  # # NOTE: At ~50um size, spots from the visium assay will encompass the expression 
  # # profiles of multiple cells. 
  # 
  # # Users may be interested to ‘deconvolute’ each of the spatial voxels to predict
  # # the underlying composition of cell types. We tested a wide variety of 
  # # deconvolution and integration methods, using a reference scRNA-seq dataset of
  # # ~14,000 adult mouse cortical cell taxonomy from the Allen Institute, 
  # # generated with the SMART-Seq2 protocol. 
  # 
  # # We consistently found superior performance using integration methods (as 
  # # opposed to deconvolution methods), likely because of substantially different 
  # # noise models that characterize spatial and single-cell datasets, and 
  # # integration methods are specifically designed to be robust to 
  # # these differences. 
  # 
  # # We therefore apply the ‘anchor’-based integration workflow introduced in 
  # # Seurat v3, that enables the probabilistic transfer of annotations from a 
  # # reference to a query set. 
  # 
  # # While many of the methods are conserved (both procedures begin by identifying 
  # # anchors), there are two important distinctions between data transfer and 
  # # integration:
  # # (i) In data transfer, Seurat does not correct or modify query expression data.
  # # (ii) In data transfer, Seurat has an option (set by default) to project the 
  # # PCA structure of a reference onto the query, instead of learning a joint 
  # # structure with CCA. We generally suggest using this option when projecting 
  # # data between scRNA-seq datasets.
  # 
  # # After finding anchors, we use the TransferData() function to classify the 
  # # query cells based on reference data (a vector of reference cell type labels). 
  # # TransferData() returns a matrix with predicted IDs and prediction scores, 
  # # which we can add to the query metadata.
  # 
  # # NOTE: Make sure reference seurat object is SCTransformed.
  # integrated_seurat <- readRDS("/hpc/home/kailasamms/scratch/scRNASeq_Chen/results_seurat/integrated_seurat_snn.rds")
  # 
  # # Seurat v3 vs Seurat v5 issues: Check if cells in graph are in same order
  # all(Cells(integrated_seurat@graphs$integrated_snn) == colnames(integrated_seurat))
  # # if FALSE, reorder the cells in all existing graphs
  # integrated_seurat@graphs$integrated_snn <- 
  #   integrated_seurat@graphs$integrated_snn[colnames(integrated_seurat), colnames(integrated_seurat)]
  # integrated_seurat@graphs$integrated_nn <- 
  #   integrated_seurat@graphs$integrated_nn[colnames(integrated_seurat), colnames(integrated_seurat)]
  # 
  # # Remove ambiguous cells before label transfer
  # integrated_seurat <- subset(integrated_seurat, 
  #                             sub_type == "Unclassified", 
  #                             invert=TRUE)
  # 
  # # NOTE: For FindTransferAnchors() to work, both reference and query MUST be
  # # SCtransformed()  https://github.com/satijalab/seurat/issues/3937
  # # So, run SCTransform() on RNA assay.
  # DefaultAssay(integrated_seurat) <- "integrated"
  # 
  # # Make sure all assays in Seuratv3 object has same order of cells
  # acells <- colnames(x = integrated_seurat[["integrated"]])
  # ocells <- colnames(x = integrated_seurat)
  # 
  # integrated_seurat@assays$RNA@counts <- integrated_seurat@assays$RNA@counts[,ocells]
  # integrated_seurat@assays$RNA@data <- integrated_seurat@assays$RNA@data[,ocells]
  # integrated_seurat@assays$integrated@data <- integrated_seurat@assays$integrated@data[,ocells]
  # integrated_seurat@assays$integrated@scale.data <- integrated_seurat@assays$integrated@scale.data[,ocells]
  # integrated_seurat@assays$SCT@counts <- integrated_seurat@assays$SCT@counts[,ocells]
  # integrated_seurat@assays$SCT@data <- integrated_seurat@assays$SCT@data[,ocells]
  # integrated_seurat@assays$SCT@scale.data <- integrated_seurat@assays$SCT@scale.data[,ocells]
  # 
  # # NOTE: If you subset the spatial object, re-run SCTranform() on it.
  # for (i in samples){
  #   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  #   DefaultAssay(sample.seurat) <- "SCT"
  #   
  #   # Run ONLY if you had subset the spatial data
  #   # sample.seurat <- SCTransform(sample.seurat, 
  #   #                              assay = "Spatial",
  #   #                              verbose = FALSE)
  #   # sample.seurat <- RunPCA(sample.seurat,
  #   #                         verbose = FALSE)
  #   
  #   # Find anchors between reference and query
  #   anchors <- FindTransferAnchors(reference = integrated_seurat,
  #                                  query = sample.seurat,
  #                                  normalization.method = "SCT",
  #                                  reference.assay = "integrated",
  #                                  query.assay = "SCT")
  #   
  #   predictions.assay <- TransferData(anchorset = anchors, 
  #                                     refdata = integrated_seurat$cell_type,
  #                                     prediction.assay = TRUE,
  #                                     weight.reduction = sample.seurat[["pca"]],
  #                                     dims = 1:30)
  #   sample.seurat[["predictions"]] <- predictions.assay
  #   
  #   DefaultAssay(sample.seurat) <- "predictions"
  #   
  #   SpatialFeaturePlot(sample.seurat, 
  #                      features = rownames(predictions.assay@data), 
  #                      pt.size.factor = 1.6,
  #                      #ncol = 4,
  #                      crop = FALSE)
  #   
  #   ggsave(filename = paste0(seurat_results, "Labeltransfer_plot_", i, ".jpg"),
  #          plot = last_plot(),
  #          units = c("in"),
  #          width = 11,
  #          height = 8)
  #   
  #   
  #   
  #   # sample.seurat <- FindSpatiallyVariableFeatures(sample.seurat,
  #   #                                                assay = "predictions",
  #   #                                                selection.method = "moransi",
  #   #                                                features = rownames(sample.seurat),
  #   #                                                r.metric = 5,
  #   #                                                slot = "data")
  #   
  #   # top.clusters <- head(SpatiallyVariableFeatures(sample.seurat, 
  #   #                                                selection.method = "moransi"), 4)
  #   # SpatialPlot(object = sample.seurat, 
  #   #             features = top.clusters,
  #   #             ncol = 2)
  # }
  
  #!/usr/bin/env Rscript
  
  # Read and store variables from command line interface (CLI)
  cli <- base::commandArgs(trailingOnly = TRUE) 
  args <- base::strsplit(x = cli, split = "=", fixed = TRUE)
  
  for (e in args){
    argname <- e[1]
    argval <- e[2]
    assign(argname, argval)
  }
  
  # NOTE: All variables and functions are defined within the file below
  source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")
  
  # NOTE: In single cell experiment, we combine multiple samples into a single 
  # seurat object and do analysis. In spatial experiment, we analyze each sample
  # individually.
  
  # NOTE: In visium v1 spatial transcriptomics, there are 4992 barcoded spots 
  # which capture RNA from tissue placed above these spots. 
  # So, unlike scRNASeq, each barcode corresponds to a spot, not a cell.
  # https://www.youtube.com/watch?v=VwNk4d-0RJc
  # https://kb.10xgenomics.com/hc/en-us/articles/360035848191-How-many-spots-are-within-a-single-capture-area-on-the-Visium-v1-Spatial-Gene-Expression-Slide
  
  # NOTE: Reads are obtained from both background (non-tissue) and tissue specific
  # areas and stored in raw_feature_bc_matrix and filtered_fearure_bc_matrix
  # respectively. Since, we are ONLY interested in reads obtained from tissues,
  # use ONLY filtered_feature_bc_matrix.
  
  #******************************************************************************#
  #                       STEP 1: SETUP THE SEURAT OBJECT                        #
  #******************************************************************************#
  
  #*********************IMPORTING GEX (GENE EXPRESSION) DATA*********************#
  
  # raw_matrix_path
  # ->Sample1
  #   ->binned_outputs
  #     ->square_002um
  #       ->raw_feature_bc_matrix.h5
  #       ->filtered_feature_bc_matrix.h5
  #       ->spatial
  #     ->square_008um
  #       ->raw_feature_bc_matrix.h5
  #       ->filtered_feature_bc_matrix.h5
  #       ->spatial
  #     ->square_016um
  #       ->raw_feature_bc_matrix.h5
  #       ->filtered_feature_bc_matrix.h5
  #       ->spatial
  
  # DIRECTORY STRUCTURE IS IMPORTANT:
  # data.dir MUST have a H5 file specified by filename parameter as well as folder
  # named "spatial" containing the image
  
  # Create a list of sample names which will be added to each barcode.
  # Since folder names correspond to sample name, we just use list.files()
  samples <- list.files(path = filt_matrix_path)
  
  # Loop through each of the individual folders in parent directory & import data
  for(i in samples){
    
    sample.seurat <- Load10X_Spatial(data.dir = paste0(filt_matrix_path, i),
                                     filename = "filtered_feature_bc_matrix.h5",
                                     assay = "Spatial",
                                     slice = i,
                                     bin.size = c(2,8,16),
                                     filter.matrix = TRUE,
                                     to.upper = FALSE,
                                     image = NULL)
    
    # Since we use filtered_matrix which ONLY has reads from tissues, 
    # filter.matrix = TRUE or FALSE didnt show any difference in sample.seurat
    
    # [i think wrong] If filter.matrix is set to FALSE, all 4992 spots will be recorded in 
    # sample.seurat@images$B8@coordinates. Else, only spots that are over tissue,
    # will be recorded in sample.seurat@images$B8@coordinates.
    
    # Unlike Seurat::Read10X(), Seurat::Load10X_Spatial doesnt have ability to 
    # specify project parameter. So, we manually do it.
    sample.seurat@meta.data <- sample.seurat@meta.data %>% 
      dplyr::mutate(orig.ident = i)
    
    # Assign the seurat object to its corresponding variable
    assign(paste0(i, ".filt"), sample.seurat)
    cat("DATA IMPORTED FOR ", i, ".filt dataset\n")
  }
  
  #******************************************************************************#
  #                           STEP 2: QUALITY CONTROL                            #
  #******************************************************************************#
  
  #**********************STEP 2A: CALCULATE ALL QC METRICS***********************#
  
  # NOTE:  We are going to do the same QC on each sample. So, we can merge the 
  # individual seurat objects into a single seurat object and calculate metrics on
  # the merged seurat object. However, if the number of cells exceeds 1.5 million,
  # then this step will fail as there are simply too many cells. So, it is better
  # to calculate QC metrics on individual samples, perform filtering on individual
  # samples and then merge the filtered seurat objects which have fewer cells.
  
  # Initialize an empty dataframe where the class of each column resembles those 
  # of raw_metadata. We will use this dataframe for making QC plots later.
  raw_metadata <- data.frame(Cell = c(""), 
                             Sample = as.factor(1), 
                             nUMIs = c(0), 
                             nGenes = c(0), 
                             MitoRatio = c(0), 
                             RiboRatio = c(0), 
                             Novelty = c(0))
  
  # Calculate QC metrics for each sample individually using raw matrices
  for(i in paste0(samples, ".filt")){
    
    sample.seurat <- get(i)
    
    # NOTE: PercentageFeatureSet() ONLY calculates metrics for cells within the 
    # indicated assay. If no assay is indicated, it will use default assay. So, we
    # loop through each assay and calculate metrics
    
    # Get list of assays
    assay.list <- Assays(sample.seurat)
    
    for (assay in assay.list){
      # Compute percent mito percent
      sample.seurat <- Seurat::PercentageFeatureSet(object = sample.seurat,
                                                    pattern = "^[Mm][Tt]-",
                                                    features = NULL,
                                                    col.name = "MitoPercent",
                                                    assay = assay)
      
      # Compute percent ribo percent
      sample.seurat <- Seurat::PercentageFeatureSet(object = sample.seurat,
                                                    pattern = "^[Rr][Pp][SsLl]", 
                                                    features = NULL,
                                                    col.name = "RiboPercent",
                                                    assay = assay)
    }
    
    # Extract metadata
    sample_metadata <- sample.seurat@meta.data
    
    # Rename columns to be more intuitive and add the additional QC metrics:
    # (i)     Cell      : unique identifiers corresponding to each cell i.e. barcodes
    # (ii)    Sample    : sample names
    # (iii)   nUMIs     : number of transcripts per cell
    # (iv)    nGenes    : number of genes per cell
    # (v)     nHTO_UMIs : number of HTO reads per cell
    # (vi)    nHTOs     : number of HTO types per cell
    # (vii)   MitoRatio : MitoPercent/100
    # (viii)	RiboRatio : RiboPercent/100  
    # (ix)    Novelty   : log ratio of genes per UMI
    sample_metadata <- sample_metadata %>% 
      dplyr::mutate(Cell = paste0(orig.ident, "_", rownames(sample_metadata)),
                    Sample = orig.ident,
                    nUMIs = dplyr::case_when(!is.na(nFeature_Spatial.002um) & !is.na(nCount_Spatial.002um) ~ nCount_Spatial.002um,
                                             !is.na(nFeature_Spatial.008um) & !is.na(nCount_Spatial.008um) ~ nCount_Spatial.008um,
                                             !is.na(nFeature_Spatial.016um) & !is.na(nCount_Spatial.016um) ~ nCount_Spatial.016um,
                                             TRUE ~ NA),
                    nGenes = dplyr::case_when(!is.na(nFeature_Spatial.002um) & !is.na(nCount_Spatial.002um) ~ nFeature_Spatial.002um,
                                              !is.na(nFeature_Spatial.008um) & !is.na(nCount_Spatial.008um) ~ nFeature_Spatial.008um,
                                              !is.na(nFeature_Spatial.016um) & !is.na(nCount_Spatial.016um) ~ nFeature_Spatial.016um,
                                              TRUE ~ NA),
                    MitoRatio = MitoPercent/100,
                    RiboRatio = RiboPercent/100,
                    Novelty = log10(nGenes)/log10(nUMIs)) %>%
      dplyr::select(Cell, Sample, contains(c("nFeature", "nCount")), nUMIs, nGenes, MitoRatio, RiboRatio, Novelty)
    
    # Replace the metadata in raw Seurat object
    sample.seurat@meta.data <- sample_metadata
    
    # Append raw metadata of each seurat object which will be used for QC plots later
    raw_metadata <- dplyr::bind_rows(raw_metadata, sample_metadata)
    
    # Assign the seurat object to its corresponding variable
    assign(i, sample.seurat)
  }
  
  #******************************STEP 2B: PERFORM QC*****************************#
  
  # In single cell analysis, we have to perform QC to remove bad cells in order to
  # accurately annotate each cell. In spatial analysis, RNA from tissue is 
  # captured within the spot below the tissue. Since sequencing depth is usually 
  # way lower than single cell, we use almost all reads (lenient cutoffs) captured
  # in the spot to annotate the spots.
  
  # NOTE: Unlike PercentageFeatureSet(), subset() uses all the cells irrespective 
  # of the default assay. 2um bins are too small to contain 100 UMIs or 50 genes 
  # and the entire Spatial.002um assay will be removed. So, we use lenient cutoffs
  
  # Perform QC for each sample individually
  for(i in paste0(samples, ".filt")){
    
    sample.seurat <- get(i)
    
    # If default assay is 2um and no cells of default assay pass QC, then Seurat 
    # will give error. So, change default assay to 8um
    DefaultAssay(sample.seurat) <- "Spatial.008um"
    
    gene_cutoff <- 1           # reduced to 1 from 250 used for spatial
    umi_cutoff <- 1            # reduced to 1 from 500 used for spatial
    #mito_cutoff <- 0.2
    #ribo_cutoff <- 0.05
    #novelty_cutoff <- 0.8  	    # use 0.8 as starting point. Maximum 0.9
    
    sample.seurat <- base::subset(x = sample.seurat,
                                  subset = (nGenes >= gene_cutoff) &
                                    (nUMIs >= umi_cutoff))
    #(MitoRatio <= mito_cutoff) &
    # (RiboRatio >= ribo_cutoff) &
    # (Novelty >= novelty_cutoff))
    
    # Assign the seurat object to its corresponding variable
    assign(i, sample.seurat)
  }
  
  # Create a merged Seurat object
  # NOTE: Samples will have same barcodes. To keep track of cell identities 
  # (i.e. barcodes) coming from each sample after merging, we add a prefix 
  # (i.e. sample name) to each barcode using "add.cell.ids"
  
  # IMPORTANT: We DO NOT merge samples in spatial experiments. We analyze each
  # sample individually.
  
  # NOTE: spatial objects have multiple assays. add.cell.ids modifies rownames
  # ONLY for default assay. So, sample name is not added to cells from other
  # assays. This creates problem in identifying common_bc. So, we exclude
  # add.cell.ids
  
  #****************************STEP 2C: SAVE THE DATA****************************#
  
  for(i in paste0(samples, ".filt")){
    
    sample.seurat <- get(i)
    saveRDS(sample.seurat, file=paste0(seurat_results, i,".rds"))
  }
  
  #******************************************************************************#
  #                       STEP 3: RUN THE STANDARD PIPELINE                      #
  #******************************************************************************#
  
  # Create workbook to save markers
  wb <- openxlsx::createWorkbook()
  for (i in samples){
    
    object <- get(paste0(i, ".filt"))
    
    # Run workflow on all assays
    for (assay in c("Spatial.002um", "Spatial.008um", "Spatial.016um")){
      
      DefaultAssay(object) <- assay
      object <- Seurat::NormalizeData(       object = object, assay = assay, normalization.method = "LogNormalize")
      object <- Seurat::FindVariableFeatures(object = object, assay = assay, nfeatures = 2000)
      object <- Seurat::ScaleData(           object = object, assay = assay, features = VariableFeatures(object))
      
      # Create a new 'sketch' assay using 50k cells
      sample_n <- 50000
      e <- "Error"
      while (class(e) == "character"){
        e <- tryCatch(Seurat::SketchData(object = object, assay = assay,
                                         ncells = sample_n,
                                         method = "LeverageScore",
                                         features = VariableFeatures(object),
                                         sketched.assay = paste0(assay, ".sketch")), error = function(msg){
                                           print("Reducing number of cells sampled")
                                           return("Error")})
        sample_n <- sample_n-5000
      }
      
      object <- e
      
      # Switch analysis to sketched cells
      DefaultAssay(object) <- paste0(assay, ".sketch")
      object <- Seurat::FindVariableFeatures(object = object, 
                                             assay = paste0(assay, ".sketch"), 
                                             nfeatures = 2000)
      
      object <- Seurat::ScaleData(object = object, 
                                  assay = paste0(assay, ".sketch"), 
                                  features = VariableFeatures(object))
      
      object <- Seurat::RunPCA(object = object, 
                               assay = paste0(assay, ".sketch"), 
                               reduction.name = paste0(assay, ".pca.sketch"))
      
      object <- Seurat::FindNeighbors(object = object, 
                                      assay = paste0(assay, ".sketch"), 
                                      reduction = paste0(assay, ".pca.sketch"), 
                                      dims = 1:50)
      
      object <- Seurat::FindClusters(object = object, 
                                     algorithm = 4,
                                     cluster.name = paste0(assay, ".cluster.sketch"), 
                                     resolution = 0.6)
      
      object <- Seurat::RunUMAP(object = object, 
                                reduction.name = paste0(assay, ".umap.sketch"), 
                                reduction = paste0(assay, ".pca.sketch"), 
                                dims = 1:50)
      
      # Project the cluster labels, dimensional reductions (PCA and UMAP) that we
      # learned from the 50,000 sketched cells to the entire dataset
      object <- Seurat::ProjectData(object = object, 
                                    assay = assay,
                                    sketched.assay = paste0(assay, ".sketch"),
                                    sketched.reduction = paste0(assay, ".pca.sketch"),
                                    full.reduction = paste0(assay, ".pca.full"),
                                    dims = 1:50,
                                    refdata = list(seurat_cluster.projected = paste0(assay, ".cluster.sketch")),
                                    umap.model = paste0(assay, ".umap.sketch"))
      
      object <- Seurat::RunUMAP(object = object, 
                                reduction.name = paste0(assay, ".umap.full"), 
                                reduction = paste0(assay, ".pca.full"), 
                                dims = 1:50)
      
      # Append assay name to the newly generated metadata columns so they dont get overwritten
      object@meta.data <- object@meta.data %>%
        dplyr::rename(!!rlang::sym(paste0(assay, ".cluster.full")) := "seurat_cluster.projected",
                      !!rlang::sym(paste0(assay, ".cluster.full.score")) := "seurat_cluster.projected.score")
      
      # We can visualize the clustering results for the sketched cells, as well as the
      # projected clustering results for the full dataset:
      DefaultAssay(object) <- paste0(assay, ".sketch")
      Idents(object) <- paste0(assay, ".cluster.sketch")
      p1 <- DimPlot(object, reduction = paste0(assay, ".umap.sketch"), label = T) + 
        ggtitle("Sketched clustering") + 
        theme(legend.position = "none")
      
      # switch to full dataset
      DefaultAssay(object) <- assay
      Idents(object) <- paste0(assay, ".cluster.full")
      p2 <- DimPlot(object, reduction = paste0(assay, ".umap.full"), label = T) + 
        ggtitle("Projected clustering (full dataset)") + 
        theme(legend.position = "none")
      
      # Save the plot
      ggplot2::ggsave(filename=paste0(assay, ".umap.tiff"),
                      plot=p1+p2,
                      device="jpeg",
                      path=diagnostics_path,
                      width=11,
                      height=8.5,
                      units=c("in"),
                      dpi=600,
                      limitsize=TRUE,
                      bg="white")
      
      # Find markers
      markers <- Seurat::FindAllMarkers(object = object, 
                                        assay = assay,
                                        only.pos = TRUE)
      
      top10 <- markers %>%
        dplyr::mutate(cluster = as.numeric(as.character(cluster))) %>%
        dplyr::group_by(cluster) %>%
        dplyr::filter(avg_log2FC > 1, p_val_adj < 0.05) %>%
        slice_max(order_by=avg_log2FC*pct.1, n = 10) %>%
        ungroup() %>%
        data.frame() %>%
        dplyr::arrange(cluster, desc(avg_log2FC*pct.1))
      
      # Save markers to worksheet
      openxlsx::addWorksheet(wb, sheetName = paste0(i, assay))
      openxlsx::writeData(wb, sheet = paste0(i, assay), x = markers, rowNames = FALSE)
      openxlsx::addWorksheet(wb, sheetName = paste0(i, assay, "top10"))
      openxlsx::writeData(wb, sheet = paste0(i, assay, "top10"), x = top10, rowNames = FALSE)
      
    }
    
    # Save xlsx file with markers
    openxlsx::saveWorkbook(wb, file = paste0(seurat_results, i, "_Markers.xlsx"), overwrite = TRUE)
    
    # Save the object
    saveRDS(object, file=paste0(seurat_results, i,".rds"))
  }
  
  # DefaultAssay(object) <- "Spatial.008um"
  # 
  # vln.plot <- VlnPlot(object, features = "nCount_Spatial.008um", pt.size = 0) + 
  #   theme(axis.text = element_text(size = 4)) + 
  #   NoLegend()
  # count.plot <- SpatialFeaturePlot(object, features = "nCount_Spatial.008um") + 
  #   theme(legend.position = "right")
  
  # # switch back to 8um
  # #object@assays$Spatial.008um@features@dimnames[[1]][25:30]
  # DefaultAssay(object) <- "Spatial.008um"
  # p2 <- SpatialFeaturePlot(object, features = "GAPDH") + ggtitle("Hpca expression (8um)")
  
  
  #   # Visualize the unsupervised clusters based on their spatial location. 
  #   SpatialDimPlot(object, label = T, repel = T, label.size = 4)
  #   
  #   # Plot the spatial location of different clusters individually.
  #   Idents(object) <- "seurat_cluster.projected"
  #   cells <- CellsByIdentities(object, idents = c(0, 4, 32, 34, 35))
  #   p <- SpatialDimPlot(object,
  #                       cells.highlight = cells[setdiff(names(cells), "NA")],
  #                       cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
  #   ) + NoLegend()
  #   p
  #   
  #   
  #   # visualize the top gene expression markers for each cluster:
  #   # Create downsampled object to make visualization either
  #   DefaultAssay(object) <- "Spatial.008um"
  #   Idents(object) <- "seurat_cluster.projected"
  #   object_subset <- subset(object, cells = Cells(object[["Spatial.008um"]]), downsample = 1000)
  #   
  #   
  #   
  #   object_subset <- ScaleData(object_subset, assay = "Spatial.008um", features = top5$gene)
  #   p <- DoHeatmap(object_subset, assay = "Spatial.008um", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
  #   p
  #   
  #   # Identifying spatially-defined tissue domains
  #   object <- Banksy::RunBanksy(object,
  #                               lambda = 0.8, verbose = TRUE,
  #                               assay = "Spatial.008um", slot = "data", features = "variable",
  #                               k_geom = 50)
  #   
  #   DefaultAssay(object) <- "BANKSY"
  #   object <- RunPCA(object, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(object), npcs = 30)
  #   object <- FindNeighbors(object, reduction = "pca.banksy", dims = 1:30)
  #   object <- FindClusters(object, cluster.name = "banksy_cluster", resolution = 0.5)
  #   
  #   Idents(object) <- "banksy_cluster"
  #   p <- SpatialDimPlot(object, group.by = "banksy_cluster", label = T, repel = T, label.size = 4)
  #   p
  #   
  #   # highlight the spatial location of each tissue domain individually
  #   banksy_cells <- CellsByIdentities(object)
  #   p <- SpatialDimPlot(object, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T) + NoLegend()
  #   p
  #   
  # }
  # 
  
  # 
  # # Load rds file of seurat objects
  # for (i in samples){
  #   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  #   sample.seurat <- sctransform_spatial_data(sample.seurat)
  #   sample.seurat <- cluster_spatial_data(sample.seurat)
  #   saveRDS(sample.seurat, file=paste0(seurat_results, i,".rds"))
  #   
  #   DefaultAssay(sample.seurat) <- "SCT"
  #   p1 <- DimPlot(sample.seurat, reduction = "umap", label = TRUE)
  #   p2 <- SpatialDimPlot(sample.seurat, label = TRUE, label.size = 3)
  #   p <- p1 + p2
  #   ggsave(filename = paste0(seurat_results, "Clusters_on_slide_", i, ".jpg"), 
  #          plot = p)
  # }
  # 
  # # # There is no integration like scRNA Seq. We analyse each slide individually.
  # # # Unfortunately, all images are stored within each sample. So, we remove 
  # # # unwanted images from each sample. If more than 1 image is present in each 
  # # # sample, SpatialDimPlot() will give error.
  # # integ_data <- sct_data
  # # for (i in 1:length(sct_data)){
  # #   integ_data[[i]] <- cluster_data(sct_data[[i]])
  # #   integ_data[[i]]@images <- integ_data[[i]]@images[names(integ_data[[i]]@images) == names(sct_data)[[i]]]
  # # }
  # 
  # # Color spots were GOI are present based on expression
  # for (i in samples){
  #   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  #   DefaultAssay(sample.seurat) <- "SCT"
  #   Seurat::SpatialFeaturePlot(object = sample.seurat, 
  #                              features = c("CD8A", "CD8B", "NPEPPS", "CDH12"),
  #                              ncol = 4,
  #                              slot = "data")
  #   ggsave(filename = paste0(seurat_results, "Feature_plot_", i, ".jpg"),
  #          plot = last_plot(),
  #          units = c("in"),
  #          width = 11,
  #          height = 8)
  #   
  #   # This can ONLY plot 2 genes at a time
  #   SpatialFeaturePlotBlend(sample.seurat, "CD8A", "NPEPPS")
  # }
  # 
  # # Color spots were CD8 and NPEPPS are expressed in same plot. 
  # # NOTE: This is NOT expression based. We just color the cells that express our 
  # # GOI in different colors.
  # # Use SCT assay to get counts, not Spatial assay as it has too much background.
  # for (i in samples){
  #   
  #   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  #   
  #   # Make sure the genes you want are present in the assay
  #   GOI <- intersect(c("CD8A", "CD8B"), rownames(sample.seurat@assays$SCT@data))
  #   cd8_df <- sample.seurat@assays$SCT$counts[GOI, ]
  #   
  #   if (length(GOI) > 1){
  #     cd8_cells <- colnames(cd8_df[,colSums(cd8_df) > 0])
  #   } else {
  #     cd8_cells <- names(cd8_df[(cd8_df>0)])
  #   }
  #   
  #   GOI <- intersect(c("NPEPPS"), rownames(sample.seurat@assays$SCT@data))
  #   npepps_df <- sample.seurat@assays$SCT$counts[c("NPEPPS"), ]
  #   npepps_cells <- names(npepps_df[(npepps_df>0)])
  #   
  #   SpatialPlot(object = sample.seurat, 
  #               cells.highlight = list(CD8 = cd8_cells, NPEPPS = npepps_cells), 
  #               cols.highlight =  c("green", "blue", "grey"),
  #               pt.size.factor = 2)
  #   
  #   ggsave(filename = paste0(seurat_results, "Location_plot_", i, ".jpg"),
  #          plot = last_plot(),
  #          units = c("in"),
  #          width = 11,
  #          height = 8)
  # }
  # 
  # # # Classify each spot as CD8+ or CD8-
  # # #for (i in samples){
  # # i <- "B8"
  # #   
  # #   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  # #   #rownames(sample.seurat@meta.data) <- sample.seurat@meta.data$Cell
  # #   
  # #   df <- sample.seurat@assays$Spatial$counts[c("CD8A", "CD8B"), ]
  # #   pos_cells  <- colnames(df[,colSums(df) > 0])
  # #   neg_cells <- setdiff(colnames(df), pos_cells)
  # #   pos_cells <- paste0(i,"_", pos_cells)
  # #   neg_cells <- paste0(i,"_", neg_cells)
  # #   sample.seurat@meta.data <- sample.seurat@meta.data %>% 
  # #     dplyr::mutate(CD8_status = dplyr::case_when(Cell %in% pos_cells ~ "CD8_pos",
  # #                                                 TRUE ~ "CD8_neg"))
  # # }
  # 
  # 
  # i <- "B8"
  # sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  # 
  # a <- GetTissueCoordinates(sample.seurat)
  # df <- data.frame(sample.seurat@images$B8@coordinates)
  # 
  # ### CHECK THAT WE ARE SUCCESSFULLY ABLE TO IDENTIFY NEIGHBORS ###
  # # Find a cell in the center of the tissue
  # row_cell <- median(df$row)
  # col_cell <- median(df$col)
  # cell <- rownames(df %>% dplyr::filter(row==row_cell, col==col_cell))
  # 
  # # Find adjacent cells
  # cell1 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell-2)))
  # cell2 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell+2)))
  # cell3 <- rownames(df %>% dplyr::filter(row==(row_cell-2), col==col_cell))
  # cell4 <- rownames(df %>% dplyr::filter(row==(row_cell+2), col==col_cell))
  # cell5 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell-1)))
  # cell6 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell+1)))
  # cell7 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell-1)))
  # cell8 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell+1)))
  # 
  # SpatialPlot(object = sample.seurat, 
  #             cells.highlight = list(cell, unlist(lapply(paste0("cell", seq(1,8)), get), use.names=FALSE)), 
  #             #facet.highlight = TRUE,
  #             cols.highlight =  c("yellow", "red", "black"), 
  #             pt.size.factor = 2)
  # 
  # #******************************************************************************#
  # 
  # i <- "B8"
  # sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  # df <- data.frame(sample.seurat@images$B8@coordinates)
  # expr_df <- sample.seurat@assays$SCT$data[rownames(sample.seurat@assays$SCT$data) %in% c("NPEPPS", "CD8A", "CD8B"),]
  # 
  # # Identify cells that express NPEPPS
  # npepps <- expr_df["NPEPPS",]
  # npepps <- npepps[npepps > 0]
  # npepps <- names(npepps)
  # # Not all npepps cells have co-ordinates. Remove those that lack co-ordinates.
  # npepps <- intersect(npepps, rownames(df))
  # 
  # expr_df <- data.frame(expr_df)
  # # For every cell that expresses NPEPPS, calculate total expression of CD8A and 
  # # CD8B in 1st neighbor
  # t_expr <- c()
  # npepps_expr <- c()
  # for (i in npepps){
  #   row_cell <- df[rownames(df) == i,]$row
  #   col_cell <- df[rownames(df) == i,]$col
  #   
  #   # Find adjacent cells
  #   cell1 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell-2)))
  #   cell2 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell+2)))
  #   cell3 <- rownames(df %>% dplyr::filter(row==(row_cell-2), col==col_cell))
  #   cell4 <- rownames(df %>% dplyr::filter(row==(row_cell+2), col==col_cell))
  #   cell5 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell-1)))
  #   cell6 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell+1)))
  #   cell7 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell-1)))
  #   cell8 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell+1)))
  #   
  #   # SpatialPlot(object = sample.seurat, 
  #   #             cells.highlight = list(i, unlist(lapply(paste0("cell", seq(1,8)), get), use.names=FALSE)), 
  #   #             #facet.highlight = TRUE,
  #   #             cols.highlight =  c("yellow", "red", "black"), 
  #   #             pt.size.factor = 2)
  #   
  #   cells <- intersect(make.names(unlist(lapply(paste0("cell", seq(1,8)), get))), colnames(expr_df))
  #   print(cells)
  #   
  #   t <- expr_df %>% 
  #     dplyr::select(all_of(cells)) %>%
  #     tibble::rownames_to_column("Gene") %>%
  #     dplyr::filter(Gene != "NPEPPS") %>%
  #     tibble::column_to_rownames("Gene")
  #   
  #   n <- expr_df %>% 
  #     dplyr::select(all_of(cells)) %>%
  #     tibble::rownames_to_column("Gene") %>%
  #     dplyr::filter(Gene == "NPEPPS") %>%
  #     tibble::column_to_rownames("Gene")
  #   
  #   t_expr <- c(t_expr, sum(t))
  #   npepps_expr <- c(npepps_expr, sum(n))
  # }
  # 
  # df <- data.frame(npepps_expr, t_expr)
  # # Save batch corrected normalized counts for entire dataset
  # wb <- openxlsx::createWorkbook()
  # openxlsx::addWorksheet(wb, sheetName = "spatial")
  # openxlsx::writeData(wb, sheet = "spatial", x = df, rowNames = FALSE)
  # openxlsx::saveWorkbook(wb,
  #                        file = paste0(seurat_results, "Correlation.xlsx"),
  #                        overwrite = TRUE)
  # 
  
  # 
  # # pseudobulk the counts based on donor-condition-celltype
  # pseudo_ifnb <- AggregateExpression(ifnb, assays = "RNA", return.seurat = T, group.by = c("stim", "donor_id", "seurat_annotations"))
  # 
  # # each 'cell' is a donor-condition-celltype pseudobulk profile
  # tail(Cells(pseudo_ifnb))
  # 
  # # Identification of Spatially Variable Features
  # 
  # # Seurat offers two workflows to identify molecular features that correlate 
  # # with spatial location within a tissue. The first is to perform differential 
  # # expression based on pre-annotated anatomical regions within the tissue, which
  # # may be determined either from unsupervised clustering or prior knowledge. 
  # # This strategy works will in this case, as the clusters above exhibit clear 
  # # spatial restriction.
  # 
  # de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
  # SpatialFeaturePlot(object = brain, 
  #                    features = rownames(de_markers)[1:3], 
  #                    alpha = c(0.1, 1), 
  #                    ncol = 3)
  # 
  # 
  # # An alternative approach, implemented in FindSpatiallyVariables(), is to 
  # # search for features exhibiting spatial patterning in the absence of 
  # # pre-annotation. The default method (method = 'markvariogram), is inspired by 
  # # the Trendsceek, which models spatial transcriptomics data as a mark point 
  # # process and computes a ‘variogram’, which identifies genes whose expression 
  # # level is dependent on their spatial location. More specifically, this process
  # # calculates gamma(r) values measuring the dependence between two spots a 
  # # certain “r” distance apart. By default, we use an r-value of ‘5’ in these 
  # # analyses, and only compute these values for variable genes (where variation 
  # # is calculated independently of spatial location) to save time.
  # 
  # 
  # brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000],
  #                                        selection.method = "moransi")
  # top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "moransi"), 6)
  # SpatialFeaturePlot(brain, 
  #                    features = top.features, 
  #                    alpha = c(0.1, 1),
  #                    ncol = 3)
  # 
  # 
  # 
  # # Integration with single-cell data i.e. label transfer from single-cell data
  # # NOTE: At ~50um size, spots from the visium assay will encompass the expression 
  # # profiles of multiple cells. 
  # 
  # # Users may be interested to ‘deconvolute’ each of the spatial voxels to predict
  # # the underlying composition of cell types. We tested a wide variety of 
  # # deconvolution and integration methods, using a reference scRNA-seq dataset of
  # # ~14,000 adult mouse cortical cell taxonomy from the Allen Institute, 
  # # generated with the SMART-Seq2 protocol. 
  # 
  # # We consistently found superior performance using integration methods (as 
  # # opposed to deconvolution methods), likely because of substantially different 
  # # noise models that characterize spatial and single-cell datasets, and 
  # # integration methods are specifically designed to be robust to 
  # # these differences. 
  # 
  # # We therefore apply the ‘anchor’-based integration workflow introduced in 
  # # Seurat v3, that enables the probabilistic transfer of annotations from a 
  # # reference to a query set. 
  # 
  # # While many of the methods are conserved (both procedures begin by identifying 
  # # anchors), there are two important distinctions between data transfer and 
  # # integration:
  # # (i) In data transfer, Seurat does not correct or modify query expression data.
  # # (ii) In data transfer, Seurat has an option (set by default) to project the 
  # # PCA structure of a reference onto the query, instead of learning a joint 
  # # structure with CCA. We generally suggest using this option when projecting 
  # # data between scRNA-seq datasets.
  # 
  # # After finding anchors, we use the TransferData() function to classify the 
  # # query cells based on reference data (a vector of reference cell type labels). 
  # # TransferData() returns a matrix with predicted IDs and prediction scores, 
  # # which we can add to the query metadata.
  # 
  # # NOTE: Make sure reference seurat object is SCTransformed.
  # integrated_seurat <- readRDS("/hpc/home/kailasamms/scratch/scRNASeq_Chen/results_seurat/integrated_seurat_snn.rds")
  # 
  # # Seurat v3 vs Seurat v5 issues: Check if cells in graph are in same order
  # all(Cells(integrated_seurat@graphs$integrated_snn) == colnames(integrated_seurat))
  # # if FALSE, reorder the cells in all existing graphs
  # integrated_seurat@graphs$integrated_snn <- 
  #   integrated_seurat@graphs$integrated_snn[colnames(integrated_seurat), colnames(integrated_seurat)]
  # integrated_seurat@graphs$integrated_nn <- 
  #   integrated_seurat@graphs$integrated_nn[colnames(integrated_seurat), colnames(integrated_seurat)]
  # 
  # # Remove ambiguous cells before label transfer
  # integrated_seurat <- subset(integrated_seurat, 
  #                             sub_type == "Unclassified", 
  #                             invert=TRUE)
  # 
  # # NOTE: For FindTransferAnchors() to work, both reference and query MUST be
  # # SCtransformed()  https://github.com/satijalab/seurat/issues/3937
  # # So, run SCTransform() on RNA assay.
  # DefaultAssay(integrated_seurat) <- "integrated"
  # 
  # # Make sure all assays in Seuratv3 object has same order of cells
  # acells <- colnames(x = integrated_seurat[["integrated"]])
  # ocells <- colnames(x = integrated_seurat)
  # 
  # integrated_seurat@assays$RNA@counts <- integrated_seurat@assays$RNA@counts[,ocells]
  # integrated_seurat@assays$RNA@data <- integrated_seurat@assays$RNA@data[,ocells]
  # integrated_seurat@assays$integrated@data <- integrated_seurat@assays$integrated@data[,ocells]
  # integrated_seurat@assays$integrated@scale.data <- integrated_seurat@assays$integrated@scale.data[,ocells]
  # integrated_seurat@assays$SCT@counts <- integrated_seurat@assays$SCT@counts[,ocells]
  # integrated_seurat@assays$SCT@data <- integrated_seurat@assays$SCT@data[,ocells]
  # integrated_seurat@assays$SCT@scale.data <- integrated_seurat@assays$SCT@scale.data[,ocells]
  # 
  # # NOTE: If you subset the spatial object, re-run SCTranform() on it.
  # for (i in samples){
  #   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  #   DefaultAssay(sample.seurat) <- "SCT"
  #   
  #   # Run ONLY if you had subset the spatial data
  #   # sample.seurat <- SCTransform(sample.seurat, 
  #   #                              assay = "Spatial",
  #   #                              verbose = FALSE)
  #   # sample.seurat <- RunPCA(sample.seurat,
  #   #                         verbose = FALSE)
  #   
  #   # Find anchors between reference and query
  #   anchors <- FindTransferAnchors(reference = integrated_seurat,
  #                                  query = sample.seurat,
  #                                  normalization.method = "SCT",
  #                                  reference.assay = "integrated",
  #                                  query.assay = "SCT")
  #   
  #   predictions.assay <- TransferData(anchorset = anchors, 
  #                                     refdata = integrated_seurat$cell_type,
  #                                     prediction.assay = TRUE,
  #                                     weight.reduction = sample.seurat[["pca"]],
  #                                     dims = 1:30)
  #   sample.seurat[["predictions"]] <- predictions.assay
  #   
  #   DefaultAssay(sample.seurat) <- "predictions"
  #   
  #   SpatialFeaturePlot(sample.seurat, 
  #                      features = rownames(predictions.assay@data), 
  #                      pt.size.factor = 1.6,
  #                      #ncol = 4,
  #                      crop = FALSE)
  #   
  #   ggsave(filename = paste0(seurat_results, "Labeltransfer_plot_", i, ".jpg"),
  #          plot = last_plot(),
  #          units = c("in"),
  #          width = 11,
  #          height = 8)
  #   
  #   
  #   
  #   # sample.seurat <- FindSpatiallyVariableFeatures(sample.seurat,
  #   #                                                assay = "predictions",
  #   #                                                selection.method = "moransi",
  #   #                                                features = rownames(sample.seurat),
  #   #                                                r.metric = 5,
  #   #                                                slot = "data")
  #   
  #   # top.clusters <- head(SpatiallyVariableFeatures(sample.seurat, 
  #   #                                                selection.method = "moransi"), 4)
  #   # SpatialPlot(object = sample.seurat, 
  #   #             features = top.clusters,
  #   #             ncol = 2)
  # }

  
  # NOTE: In single cell experiment, we combine multiple samples into a single 
  # seurat object and do analysis. In spatial experiment, we analyze each sample
  # individually.
  
  # NOTE: Reads are obtained from both background (non-tissue) and tissue specific
  # areas and stored in raw_feature_bc_matrix and filtered_fearure_bc_matrix
  # respectively. Since, we are ONLY interested in reads obtained from tissues,
  # use ONLY filtered_feature_bc_matrix.
  
  # NOTE: Size of smallest mammalian cell is ~4um. 
  
  # NOTE: In visium v1 spatial transcriptomics, there are 4992 barcoded spots 
  # which capture RNA from tissue placed above these spots. 
  # So, unlike scRNASeq, each barcode corresponds to a spot, not a cell.
  # https://www.youtube.com/watch?v=VwNk4d-0RJc
  # https://kb.10xgenomics.com/hc/en-us/articles/360035848191-How-many-spots-are-within-a-single-capture-area-on-the-Visium-v1-Spatial-Gene-Expression-Slide
  
  # DIRECTORY STRUCTURE IS IMPORTANT:
  # -> filt_feature_bc_matrix  (dir)
  #   -> TMA1-A1 (dir)
  #     -> binned_outputs (dir)
  #       -> square_002um (dir)      
  #         -> filtered_feature_bc_matrix.h5
  #         -> spatial (dir)
  #           -> aligned_fiducials.jpg
  #           -> cytassist_image.tiff
  #           -> detected_tissue_image.jpg
  #           -> scalefactors_json.json
  #           -> tissue_hires_image.png
  #           -> tissue_lowres_image.png
  #           -> tissue_positions.parquet
  #       -> square_008um (dir)
  #         -> filtered_feature_bc_matrix.h5
  #         -> spatial (dir)
  #       -> square_016um (dir)
  #         -> filtered_feature_bc_matrix.h5
  #         -> spatial (dir)
  
  # data.dir MUST have 
  # (i) a H5 file specified by filename parameter
  # (ii) a folder named "spatial" containing the image
  
  #**********************STEP 2A: CALCULATE ALL QC METRICS***********************#
  
  # NOTE: PercentageFeatureSet() calculates metrics for cells ONLY within the 
  # default assay, if no assay is indicated. So, we loop through each assay and
  # calculate metrics
  
  # In single cell analysis, we have to perform QC to remove bad cells in order to
  # accurately annotate each cell. In spatial analysis, RNA from tissue is 
  # captured within the spot below the tissue. Since sequencing depth is usually 
  # way lower than single cell, we use almost all reads (lenient cutoffs) captured
  # in the spot to annotate the spots.
  
  # NOTE: Unlike PercentageFeatureSet(), subset() uses all the cells irrespective 
  # of the default assay. 2um bins are too small to contain 100 UMIs or 50 genes 
  # and the entire Spatial.002um assay will be removed. So, we use lenient cutoffs
  
  # STEP 3: RUN THE STANDARD PIPELINE                      #
  
  # Create workbook to save markers
  wb <- openxlsx::createWorkbook()
  for (i in samples){
    
    object <- readRDS(paste0(seurat_results, i, ".filt.rds"))
    
    # Run workflow on all assays
    for (assay in c("Spatial.008um", "Spatial.016um")){
      
      DefaultAssay(object) <- assay
      object <- Seurat::NormalizeData(       object = object, assay = assay, normalization.method = "LogNormalize")
      object <- Seurat::FindVariableFeatures(object = object, assay = assay, nfeatures = 2000)
      object <- Seurat::ScaleData(           object = object, assay = assay, features = VariableFeatures(object))
      
      # Create a new 'sketch' assay using 50k cells
      sample_n <- 50000
      e <- "Error"
      while (class(e) == "character" & sample_n > 0){
        e <- tryCatch(Seurat::SketchData(object = object, 
                                         assay = assay,
                                         ncells = sample_n,
                                         method = "LeverageScore",
                                         features = VariableFeatures(object),
                                         sketched.assay = paste0(assay, ".sketch")), error = function(msg){
                                           print("Reducing number of cells sampled")
                                           return("Error")})
        sample_n <- sample_n-5000
        cat(assay, ":", sample_n, "\n")
      }
      
      object <- e
      
      # Switch analysis to sketched cells
      DefaultAssay(object) <- paste0(assay, ".sketch")
      object <- Seurat::FindVariableFeatures(object = object, 
                                             assay = paste0(assay, ".sketch"), 
                                             nfeatures = 2000)
      
      object <- Seurat::ScaleData(object = object, 
                                  assay = paste0(assay, ".sketch"), 
                                  features = VariableFeatures(object))
      
      object <- Seurat::RunPCA(object = object, 
                               assay = paste0(assay, ".sketch"), 
                               reduction.name = paste0(assay, ".pca.sketch"))
      
      object <- Seurat::FindNeighbors(object = object, 
                                      assay = paste0(assay, ".sketch"), 
                                      reduction = paste0(assay, ".pca.sketch"), 
                                      dims = 1:50)
      
      object <- Seurat::FindClusters(object = object, 
                                     algorithm = 4,
                                     cluster.name = paste0(assay, ".cluster.sketch"), 
                                     resolution = 0.6)
      
      object <- Seurat::RunUMAP(object = object, 
                                reduction.name = paste0(assay, ".umap.sketch"), 
                                reduction = paste0(assay, ".pca.sketch"), 
                                dims = 1:50)
      
      # Project the cluster labels, dimensional reductions (PCA and UMAP) that we
      # learned from the 50,000 sketched cells to the entire dataset
      object <- Seurat::ProjectData(object = object, 
                                    assay = assay,
                                    sketched.assay = paste0(assay, ".sketch"),
                                    sketched.reduction = paste0(assay, ".pca.sketch"),
                                    full.reduction = paste0(assay, ".pca.full"),
                                    dims = 1:50,
                                    refdata = list(seurat_cluster.projected = paste0(assay, ".cluster.sketch")),
                                    umap.model = paste0(assay, ".umap.sketch"))
      
      object <- Seurat::RunUMAP(object = object, 
                                reduction.name = paste0(assay, ".umap.full"), 
                                reduction = paste0(assay, ".pca.full"), 
                                dims = 1:50)
      
      # Append assay name to the newly generated metadata columns so they dont get overwritten
      object@meta.data <- object@meta.data %>%
        dplyr::rename(!!rlang::sym(paste0(assay, ".cluster.full")) := "seurat_cluster.projected",
                      !!rlang::sym(paste0(assay, ".cluster.full.score")) := "seurat_cluster.projected.score")
      
      # We can visualize the clustering results for the sketched cells, as well as the
      # projected clustering results for the full dataset:
      DefaultAssay(object) <- paste0(assay, ".sketch")
      Idents(object) <- paste0(assay, ".cluster.sketch")
      p1 <- DimPlot(object, reduction = paste0(assay, ".umap.sketch"), label = T) + 
        ggtitle("Sketched clustering") + 
        theme(legend.position = "none")
      
      # switch to full dataset
      DefaultAssay(object) <- assay
      Idents(object) <- paste0(assay, ".cluster.full")
      p2 <- DimPlot(object, reduction = paste0(assay, ".umap.full"), label = T) + 
        ggtitle("Projected clustering (full dataset)") + 
        theme(legend.position = "none")
      
      # Save the plot
      ggplot2::ggsave(filename=paste0(i, assay, ".umap.tiff"),
                      plot=p1+p2,
                      device="jpeg",
                      path=diagnostics_path,
                      width=11,
                      height=8.5,
                      units=c("in"),
                      dpi=600,
                      limitsize=TRUE,
                      bg="white")
      
      # Find markers
      markers <- Seurat::FindAllMarkers(object = object, 
                                        assay = assay,
                                        only.pos = TRUE)
      
      top10 <- markers %>%
        dplyr::mutate(cluster = as.numeric(as.character(cluster))) %>%
        dplyr::group_by(cluster) %>%
        dplyr::filter(avg_log2FC > 1, p_val_adj < 0.05) %>%
        slice_max(order_by=avg_log2FC*pct.1, n = 10) %>%
        ungroup() %>%
        data.frame() %>%
        dplyr::arrange(cluster, desc(avg_log2FC*pct.1))
      
      # Save markers to worksheet
      openxlsx::addWorksheet(wb, sheetName = paste0(i, assay))
      openxlsx::writeData(wb, sheet = paste0(i, assay), x = markers, rowNames = FALSE)
      openxlsx::addWorksheet(wb, sheetName = paste0(i, assay, "top10"))
      openxlsx::writeData(wb, sheet = paste0(i, assay, "top10"), x = top10, rowNames = FALSE)
    }
    # Save xlsx file with markers
    openxlsx::saveWorkbook(wb, file = paste0(seurat_results, i, "_Markers.xlsx"), overwrite = TRUE)
    
    # Save the object
    saveRDS(object, file=paste0(seurat_results, i,".rds"))
  }
  
  
  
  #******************REMOVE BAD SAMPLES
  
  for(i in samples){
    
    sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
    
    if (i == "TMA1-A1"){
      sample.seurat <- subset(sample.seurat, 
                              Sample %in% c("TMA1.D1.A10", "TMA1.D1.A11","TMA1.D1.A12",
                                            "TMA1.D1.B10", "TMA1.D1.C09", 
                                            "TMA1.D1.D09", "TMA1.D1.D12"))
    }
    
    if (i == "TMA1-D1"){
      sample.seurat <- subset(sample.seurat, 
                              Sample %in% c("TMA1.D1.A10", "TMA1.D1.A11","TMA1.D1.A12",
                                            "TMA1.D1.B10", "TMA1.D1.C09", 
                                            "TMA1.D1.D09", "TMA1.D1.D12"))
      
    }
  }
  
  #******************FIND CELL TYPES BASED ON MARKERS
  
  df <- openxlsx::read.xlsx("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Markers.All.TMA1-A1.Spatial.016um.cluster.0.4.xlsx") %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::distinct(cluster, gene, avg_log2FC, pct.1, pct.2, .keep_all = TRUE) %>%
    tidyr::pivot_wider(id_cols=cluster, names_from=gene, values_from=avg_log2FC, values_fill=0.0) %>%
    tibble::column_to_rownames("cluster") %>%
    scale()
  
  row_dist <- stats::dist(x=df, method = "euclidean", diag = TRUE, upper = TRUE)
  row_clust <- hclust(row_dist)
  row_order <- rownames(df[row_clust$order,])
  col_dist <- stats::dist(x=t(df), method = "euclidean", diag = TRUE, upper = TRUE)
  col_clust <- hclust(col_dist)
  col_order <- colnames(df[,col_clust$order])
  mat_A1 <- df[row_order, col_order]
  
  df <- openxlsx::read.xlsx("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Markers.All.TMA1-A1.Spatial.008um.cluster.0.4.xlsx") %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::distinct(cluster, gene, avg_log2FC, pct.1, pct.2, .keep_all = TRUE) %>%
    tidyr::pivot_wider(id_cols=cluster, names_from=gene, values_from=avg_log2FC, values_fill=0) %>%
    tibble::column_to_rownames("cluster") %>%
    scale()
  
  row_dist <- stats::dist(x=df, method = "euclidean", diag = TRUE, upper = TRUE)
  row_clust <- hclust(row_dist)
  row_order <- rownames(df[row_clust$order,])
  col_dist <- stats::dist(x=t(df), method = "euclidean", diag = TRUE, upper = TRUE)
  col_clust <- hclust(col_dist)
  col_order <- colnames(df[,col_clust$order])
  mat_D1 <- df[row_order, col_order]
  
  mat_A1[mat_A1 < 0] <- 0 
  mat_D1[mat_D1 < 0] <- 0 
  
  mat_A1 <- mat_A1 %>% 
    data.frame() %>% 
    mutate(across(where(is.numeric), ~round(., 2)))
  
  mat_D1 <- mat_D1 %>% 
    data.frame() %>% 
    mutate(across(where(is.numeric), ~round(., 2)))
  
  # Create workbook to save markers
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "A1")
  openxlsx::writeData(wb, sheet = "A1", x = t(mat_A1), rowNames = TRUE)
  openxlsx::addWorksheet(wb, sheetName = "D1")
  openxlsx::writeData(wb, sheet = "D1", x = t(mat_D1), rowNames = TRUE)
  openxlsx::saveWorkbook(wb, file ="C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/TMA-A1.marker.new.xlsx")
  openxlsx::saveWorkbook(wb, file = paste0(seurat_results, i, "_Markers.xlsx"), overwrite = TRUE)
  
  # my_palette <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)
  # # Define breaks
  # if(max(mat) == 0){
  #   breaks <- c(seq(from = floor(min(mat)), to = 0, length.out = 100))
  #   my_palette <- my_palette[1:50]
  # } else if (min(mat) == 0){
  #   breaks <- c(seq(from = 0, to = ceiling(max(mat)), length.out = 100))
  #   my_palette <- my_palette[50:100]
  # } else if(min(mat) < -3 | max(mat) > 3){
  #   breaks <- c(seq(-3, 0, length.out = 50), seq(3/100, 3, length.out = 50))
  # } else{
  #   breaks <- c(seq(from = floor(min(mat)), to = 0, length.out = 50), seq(from = max(mat)/100, to = ceiling(max(mat)), length.out = 50))
  # }
  # pheatmap::pheatmap(df, breaks=breaks)
  
  #******************ANNOTATE 8um BINS
  
  Pancreatic.Exocrine <- c("AMY2A", "CEL", "CELA2A", "CELA2B", "CELA3A", "CELA3B",
                           "CLPS", "CPA1", "CPA2", "CPB1", "CTRB1", "CTRB2", 
                           "CTRC", "CTRL", "CUZD1", "ERP27", "GP2", "KLK1", 
                           "PDIA2", "PLA2G1B", "PNLIP", "PNLIPRP1", "PRSS1", 
                           "PRSS3", "SERPINI2", "SYCN")
  Hepatocytes <- c("AHSG", "ALB", "ALDOB", "AMBP", "ANGPTL3", "AOX1", "APCS", 
                   "APOA1", "APOA2", "APOA5", "APOC2", "APOC3", "APOH", "ARG1", 
                   "BAAT", "C3", "C4BPB", "C5", "C9", "CP", "CPB2", "CPS1", "F2",
                   "FGA", "FGB", "FGG", "FGL1", "GC", "HAMP", "HAO1", "HP", "HPR",
                   "HPX", "HRG", "ITIH3", "MAT1A", "ORM1", "ORM2", "PAH", 
                   "SERPINC1", "SLC25A47", "SLC2A2", "SULT2A1", "VTN")
  
  for (f in c("Pancreatic.Exocrine", "Hepatocytes")){
    features <- get(f)
    features <- rownames(sample.seurat@assays$Spatial.008um$data)[tolower(rownames(sample.seurat@assays$Spatial.008um$data)) %in% 
                                                                    tolower(features)]
    features <- list(sort(features))
    
    # Calculate module scores
    sample.seurat <- Seurat::AddModuleScore(object=sample.seurat,
                                            features=features,
                                            assay="Spatial.008um",
                                            slot="data",
                                            name=f)
  }
  DefaultAssay(sample.seurat) <- "Spatial.008um"
  SpatialFeaturePlot(object = sample.seurat,
                     features="Pancreatic.Exocrine1",
                     image.scale="lowres",   # "hires" hides show H&E image
                     pt.size.factor = 4)     # default value 1.6 gives weird shapes. Values greater than 4 makes dot large
  
  SpatialFeaturePlot(object = sample.seurat,
                     features="Hepatocytes1",
                     pt.size.factor = 4)  
  
  # DefaultAssay(object) <- "Spatial.008um"
  # 
  # vln.plot <- VlnPlot(object, features = "nCount_Spatial.008um", pt.size = 0) + 
  #   theme(axis.text = element_text(size = 4)) + 
  #   NoLegend()
  # count.plot <- SpatialFeaturePlot(object, features = "nCount_Spatial.008um") + 
  #   theme(legend.position = "right")
  
  # # switch back to 8um
  # #object@assays$Spatial.008um@features@dimnames[[1]][25:30]
  # DefaultAssay(object) <- "Spatial.008um"
  # p2 <- SpatialFeaturePlot(object, features = "GAPDH") + ggtitle("Hpca expression (8um)")
  
  
  #   # Visualize the unsupervised clusters based on their spatial location. 
  #   SpatialDimPlot(object, label = T, repel = T, label.size = 4)
  #   
  #   # Plot the spatial location of different clusters individually.
  #   Idents(object) <- "seurat_cluster.projected"
  #   cells <- CellsByIdentities(object, idents = c(0, 4, 32, 34, 35))
  #   p <- SpatialDimPlot(object,
  #                       cells.highlight = cells[setdiff(names(cells), "NA")],
  #                       cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
  #   ) + NoLegend()
  #   p
  #   
  #   
  #   # visualize the top gene expression markers for each cluster:
  #   # Create downsampled object to make visualization either
  #   DefaultAssay(object) <- "Spatial.008um"
  #   Idents(object) <- "seurat_cluster.projected"
  #   object_subset <- subset(object, cells = Cells(object[["Spatial.008um"]]), downsample = 1000)
  #   
  #   
  #   
  #   object_subset <- ScaleData(object_subset, assay = "Spatial.008um", features = top5$gene)
  #   p <- DoHeatmap(object_subset, assay = "Spatial.008um", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
  #   p
  #   
  #   # Identifying spatially-defined tissue domains
  #   object <- Banksy::RunBanksy(object,
  #                               lambda = 0.8, verbose = TRUE,
  #                               assay = "Spatial.008um", slot = "data", features = "variable",
  #                               k_geom = 50)
  #   
  #   DefaultAssay(object) <- "BANKSY"
  #   object <- RunPCA(object, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(object), npcs = 30)
  #   object <- FindNeighbors(object, reduction = "pca.banksy", dims = 1:30)
  #   object <- FindClusters(object, cluster.name = "banksy_cluster", resolution = 0.5)
  #   
  #   Idents(object) <- "banksy_cluster"
  #   p <- SpatialDimPlot(object, group.by = "banksy_cluster", label = T, repel = T, label.size = 4)
  #   p
  #   
  #   # highlight the spatial location of each tissue domain individually
  #   banksy_cells <- CellsByIdentities(object)
  #   p <- SpatialDimPlot(object, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T) + NoLegend()
  #   p
  #   
  # }
  # 
  
  # 
  # # Load rds file of seurat objects
  # for (i in samples){
  #   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  #   sample.seurat <- sctransform_spatial_data(sample.seurat)
  #   sample.seurat <- cluster_spatial_data(sample.seurat)
  #   saveRDS(sample.seurat, file=paste0(seurat_results, i,".rds"))
  #   
  #   DefaultAssay(sample.seurat) <- "SCT"
  #   p1 <- DimPlot(sample.seurat, reduction = "umap", label = TRUE)
  #   p2 <- SpatialDimPlot(sample.seurat, label = TRUE, label.size = 3)
  #   p <- p1 + p2
  #   ggsave(filename = paste0(seurat_results, "Clusters_on_slide_", i, ".jpg"), 
  #          plot = p)
  # }
  # 
  # # # There is no integration like scRNA Seq. We analyse each slide individually.
  # # # Unfortunately, all images are stored within each sample. So, we remove 
  # # # unwanted images from each sample. If more than 1 image is present in each 
  # # # sample, SpatialDimPlot() will give error.
  # # integ_data <- sct_data
  # # for (i in 1:length(sct_data)){
  # #   integ_data[[i]] <- cluster_data(sct_data[[i]])
  # #   integ_data[[i]]@images <- integ_data[[i]]@images[names(integ_data[[i]]@images) == names(sct_data)[[i]]]
  # # }
  # 
  # # Color spots were GOI are present based on expression
  # for (i in samples){
  #   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  #   DefaultAssay(sample.seurat) <- "SCT"
  #   Seurat::SpatialFeaturePlot(object = sample.seurat, 
  #                              features = c("CD8A", "CD8B", "NPEPPS", "CDH12"),
  #                              ncol = 4,
  #                              slot = "data")
  #   ggsave(filename = paste0(seurat_results, "Feature_plot_", i, ".jpg"),
  #          plot = last_plot(),
  #          units = c("in"),
  #          width = 11,
  #          height = 8)
  #   
  #   # This can ONLY plot 2 genes at a time
  #   SpatialFeaturePlotBlend(sample.seurat, "CD8A", "NPEPPS")
  # }
  # 
  # # Color spots were CD8 and NPEPPS are expressed in same plot. 
  # # NOTE: This is NOT expression based. We just color the cells that express our 
  # # GOI in different colors.
  # # Use SCT assay to get counts, not Spatial assay as it has too much background.
  # for (i in samples){
  #   
  #   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  #   
  #   # Make sure the genes you want are present in the assay
  #   GOI <- intersect(c("CD8A", "CD8B"), rownames(sample.seurat@assays$SCT@data))
  #   cd8_df <- sample.seurat@assays$SCT$counts[GOI, ]
  #   
  #   if (length(GOI) > 1){
  #     cd8_cells <- colnames(cd8_df[,colSums(cd8_df) > 0])
  #   } else {
  #     cd8_cells <- names(cd8_df[(cd8_df>0)])
  #   }
  #   
  #   GOI <- intersect(c("NPEPPS"), rownames(sample.seurat@assays$SCT@data))
  #   npepps_df <- sample.seurat@assays$SCT$counts[c("NPEPPS"), ]
  #   npepps_cells <- names(npepps_df[(npepps_df>0)])
  #   
  #   SpatialPlot(object = sample.seurat, 
  #               cells.highlight = list(CD8 = cd8_cells, NPEPPS = npepps_cells), 
  #               cols.highlight =  c("green", "blue", "grey"),
  #               pt.size.factor = 2)
  #   
  #   ggsave(filename = paste0(seurat_results, "Location_plot_", i, ".jpg"),
  #          plot = last_plot(),
  #          units = c("in"),
  #          width = 11,
  #          height = 8)
  # }
  # 
  # # # Classify each spot as CD8+ or CD8-
  # # #for (i in samples){
  # # i <- "B8"
  # #   
  # #   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  # #   #rownames(sample.seurat@meta.data) <- sample.seurat@meta.data$Cell
  # #   
  # #   df <- sample.seurat@assays$Spatial$counts[c("CD8A", "CD8B"), ]
  # #   pos_cells  <- colnames(df[,colSums(df) > 0])
  # #   neg_cells <- setdiff(colnames(df), pos_cells)
  # #   pos_cells <- paste0(i,"_", pos_cells)
  # #   neg_cells <- paste0(i,"_", neg_cells)
  # #   sample.seurat@meta.data <- sample.seurat@meta.data %>% 
  # #     dplyr::mutate(CD8_status = dplyr::case_when(Cell %in% pos_cells ~ "CD8_pos",
  # #                                                 TRUE ~ "CD8_neg"))
  # # }
  # 
  # 
  # i <- "B8"
  # sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  # 
  a <- GetTissueCoordinates(sample.seurat)
  df <- data.frame(sample.seurat@images$B8@coordinates)
  
  ### CHECK THAT WE ARE SUCCESSFULLY ABLE TO IDENTIFY NEIGHBORS ###
  # Find a cell in the center of the tissue
  row_cell <- median(df$x)
  col_cell <- median(df$y)
  cell <- rownames(df %>% dplyr::filter(x==row_cell, y==col_cell))
  
  # # Find adjacent cells
  # cell1 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell-2)))
  # cell2 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell+2)))
  # cell3 <- rownames(df %>% dplyr::filter(row==(row_cell-2), col==col_cell))
  # cell4 <- rownames(df %>% dplyr::filter(row==(row_cell+2), col==col_cell))
  # cell5 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell-1)))
  # cell6 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell+1)))
  # cell7 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell-1)))
  # cell8 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell+1)))
  # 
  SpatialPlot(object = sample.seurat,
              cells.highlight = list(cell, unlist(lapply(paste0("cell", seq(1,8)), get), use.names=FALSE)),
              #facet.highlight = TRUE,
              cols.highlight =  c("yellow", "red", "black"),
              pt.size.factor = 2)
  
  SingleImagePlot(data = sample.seurat,
                  # 
                  # #******************************************************************************#
                  # 
                  # i <- "B8"
                  # sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
                  # df <- data.frame(sample.seurat@images$B8@coordinates)
                  # expr_df <- sample.seurat@assays$SCT$data[rownames(sample.seurat@assays$SCT$data) %in% c("NPEPPS", "CD8A", "CD8B"),]
                  # 
                  # # Identify cells that express NPEPPS
                  # npepps <- expr_df["NPEPPS",]
                  # npepps <- npepps[npepps > 0]
                  # npepps <- names(npepps)
                  # # Not all npepps cells have co-ordinates. Remove those that lack co-ordinates.
                  # npepps <- intersect(npepps, rownames(df))
                  # 
                  # expr_df <- data.frame(expr_df)
                  # # For every cell that expresses NPEPPS, calculate total expression of CD8A and 
                  # # CD8B in 1st neighbor
                  # t_expr <- c()
                  # npepps_expr <- c()
                  # for (i in npepps){
                  #   row_cell <- df[rownames(df) == i,]$row
                  #   col_cell <- df[rownames(df) == i,]$col
                  #   
                  #   # Find adjacent cells
                  #   cell1 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell-2)))
                  #   cell2 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell+2)))
                  #   cell3 <- rownames(df %>% dplyr::filter(row==(row_cell-2), col==col_cell))
                  #   cell4 <- rownames(df %>% dplyr::filter(row==(row_cell+2), col==col_cell))
                  #   cell5 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell-1)))
                  #   cell6 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell+1)))
                  #   cell7 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell-1)))
                  #   cell8 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell+1)))
                  #   
                  #   # SpatialPlot(object = sample.seurat, 
                  #   #             cells.highlight = list(i, unlist(lapply(paste0("cell", seq(1,8)), get), use.names=FALSE)), 
                  #   #             #facet.highlight = TRUE,
                  #   #             cols.highlight =  c("yellow", "red", "black"), 
                  #   #             pt.size.factor = 2)
                  #   
                  #   cells <- intersect(make.names(unlist(lapply(paste0("cell", seq(1,8)), get))), colnames(expr_df))
                  #   print(cells)
                  #   
                  #   t <- expr_df %>% 
                  #     dplyr::select(all_of(cells)) %>%
                  #     tibble::rownames_to_column("Gene") %>%
                  #     dplyr::filter(Gene != "NPEPPS") %>%
                  #     tibble::column_to_rownames("Gene")
                  #   
                  #   n <- expr_df %>% 
                  #     dplyr::select(all_of(cells)) %>%
                  #     tibble::rownames_to_column("Gene") %>%
                  #     dplyr::filter(Gene == "NPEPPS") %>%
                  #     tibble::column_to_rownames("Gene")
                  #   
                  #   t_expr <- c(t_expr, sum(t))
                  #   npepps_expr <- c(npepps_expr, sum(n))
                  # }
                  # 
                  # df <- data.frame(npepps_expr, t_expr)
                  # # Save batch corrected normalized counts for entire dataset
                  # wb <- openxlsx::createWorkbook()
                  # openxlsx::addWorksheet(wb, sheetName = "spatial")
                  # openxlsx::writeData(wb, sheet = "spatial", x = df, rowNames = FALSE)
                  # openxlsx::saveWorkbook(wb,
                  #                        file = paste0(seurat_results, "Correlation.xlsx"),
                  #                        overwrite = TRUE)
                  # 
                  
                  # 
                  # # pseudobulk the counts based on donor-condition-celltype
                  # pseudo_ifnb <- AggregateExpression(ifnb, assays = "RNA", return.seurat = T, group.by = c("stim", "donor_id", "seurat_annotations"))
                  # 
                  # # each 'cell' is a donor-condition-celltype pseudobulk profile
                  # tail(Cells(pseudo_ifnb))
                  # 
                  # # Identification of Spatially Variable Features
                  # 
                  # # Seurat offers two workflows to identify molecular features that correlate 
                  # # with spatial location within a tissue. The first is to perform differential 
                  # # expression based on pre-annotated anatomical regions within the tissue, which
                  # # may be determined either from unsupervised clustering or prior knowledge. 
                  # # This strategy works will in this case, as the clusters above exhibit clear 
                  # # spatial restriction.
                  # 
                  # de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
                  # SpatialFeaturePlot(object = brain, 
                  #                    features = rownames(de_markers)[1:3], 
                  #                    alpha = c(0.1, 1), 
                  #                    ncol = 3)
                  # 
                  # 
                  # # An alternative approach, implemented in FindSpatiallyVariables(), is to 
                  # # search for features exhibiting spatial patterning in the absence of 
                  # # pre-annotation. The default method (method = 'markvariogram), is inspired by 
                  # # the Trendsceek, which models spatial transcriptomics data as a mark point 
                  # # process and computes a ‘variogram’, which identifies genes whose expression 
                  # # level is dependent on their spatial location. More specifically, this process
                  # # calculates gamma(r) values measuring the dependence between two spots a 
                  # # certain “r” distance apart. By default, we use an r-value of ‘5’ in these 
                  # # analyses, and only compute these values for variable genes (where variation 
                  # # is calculated independently of spatial location) to save time.
                  # 
                  # 
                  # brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000],
                  #                                        selection.method = "moransi")
                  # top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "moransi"), 6)
                  # SpatialFeaturePlot(brain, 
                  #                    features = top.features, 
                  #                    alpha = c(0.1, 1),
                  #                    ncol = 3)
                  # 
                  # 
                  # 
                  # # Integration with single-cell data i.e. label transfer from single-cell data
                  # # NOTE: At ~50um size, spots from the visium assay will encompass the expression 
                  # # profiles of multiple cells. 
                  # 
                  # # Users may be interested to ‘deconvolute’ each of the spatial voxels to predict
                  # # the underlying composition of cell types. We tested a wide variety of 
                  # # deconvolution and integration methods, using a reference scRNA-seq dataset of
                  # # ~14,000 adult mouse cortical cell taxonomy from the Allen Institute, 
                  # # generated with the SMART-Seq2 protocol. 
                  # 
                  # # We consistently found superior performance using integration methods (as 
                  # # opposed to deconvolution methods), likely because of substantially different 
                  # # noise models that characterize spatial and single-cell datasets, and 
                  # # integration methods are specifically designed to be robust to 
                  # # these differences. 
                  # 
                  # # We therefore apply the ‘anchor’-based integration workflow introduced in 
                  # # Seurat v3, that enables the probabilistic transfer of annotations from a 
                  # # reference to a query set. 
                  # 
                  # # While many of the methods are conserved (both procedures begin by identifying 
                  # # anchors), there are two important distinctions between data transfer and 
                  # # integration:
                  # # (i) In data transfer, Seurat does not correct or modify query expression data.
                  # # (ii) In data transfer, Seurat has an option (set by default) to project the 
                  # # PCA structure of a reference onto the query, instead of learning a joint 
                  # # structure with CCA. We generally suggest using this option when projecting 
                  # # data between scRNA-seq datasets.
                  # 
                  # # After finding anchors, we use the TransferData() function to classify the 
                  # # query cells based on reference data (a vector of reference cell type labels). 
                  # # TransferData() returns a matrix with predicted IDs and prediction scores, 
                  # # which we can add to the query metadata.
                  # 
                  # # NOTE: Make sure reference seurat object is SCTransformed.
                  # integrated_seurat <- readRDS("/hpc/home/kailasamms/scratch/scRNASeq_Chen/results_seurat/integrated_seurat_snn.rds")
                  # 
                  # # Seurat v3 vs Seurat v5 issues: Check if cells in graph are in same order
                  # all(Cells(integrated_seurat@graphs$integrated_snn) == colnames(integrated_seurat))
                  # # if FALSE, reorder the cells in all existing graphs
                  # integrated_seurat@graphs$integrated_snn <- 
                  #   integrated_seurat@graphs$integrated_snn[colnames(integrated_seurat), colnames(integrated_seurat)]
                  # integrated_seurat@graphs$integrated_nn <- 
                  #   integrated_seurat@graphs$integrated_nn[colnames(integrated_seurat), colnames(integrated_seurat)]
                  # 
                  # # Remove ambiguous cells before label transfer
                  # integrated_seurat <- subset(integrated_seurat, 
                  #                             sub_type == "Unclassified", 
                  #                             invert=TRUE)
                  # 
                  # # NOTE: For FindTransferAnchors() to work, both reference and query MUST be
                  # # SCtransformed()  https://github.com/satijalab/seurat/issues/3937
                  # # So, run SCTransform() on RNA assay.
                  # DefaultAssay(integrated_seurat) <- "integrated"
                  # 
                  # # Make sure all assays in Seuratv3 object has same order of cells
                  # acells <- colnames(x = integrated_seurat[["integrated"]])
                  # ocells <- colnames(x = integrated_seurat)
                  # 
                  # integrated_seurat@assays$RNA@counts <- integrated_seurat@assays$RNA@counts[,ocells]
                  # integrated_seurat@assays$RNA@data <- integrated_seurat@assays$RNA@data[,ocells]
                  # integrated_seurat@assays$integrated@data <- integrated_seurat@assays$integrated@data[,ocells]
                  # integrated_seurat@assays$integrated@scale.data <- integrated_seurat@assays$integrated@scale.data[,ocells]
                  # integrated_seurat@assays$SCT@counts <- integrated_seurat@assays$SCT@counts[,ocells]
                  # integrated_seurat@assays$SCT@data <- integrated_seurat@assays$SCT@data[,ocells]
                  # integrated_seurat@assays$SCT@scale.data <- integrated_seurat@assays$SCT@scale.data[,ocells]
                  # 
                  # # NOTE: If you subset the spatial object, re-run SCTranform() on it.
                  # for (i in samples){
                  #   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
                  #   DefaultAssay(sample.seurat) <- "SCT"
                  #   
                  #   # Run ONLY if you had subset the spatial data
                  #   # sample.seurat <- SCTransform(sample.seurat, 
                  #   #                              assay = "Spatial",
                  #   #                              verbose = FALSE)
                  #   # sample.seurat <- RunPCA(sample.seurat,
                  #   #                         verbose = FALSE)
                  #   
                  #   # Find anchors between reference and query
                  #   anchors <- FindTransferAnchors(reference = integrated_seurat,
                  #                                  query = sample.seurat,
                  #                                  normalization.method = "SCT",
                  #                                  reference.assay = "integrated",
                  #                                  query.assay = "SCT")
                  #   
                  #   predictions.assay <- TransferData(anchorset = anchors, 
                  #                                     refdata = integrated_seurat$cell_type,
                  #                                     prediction.assay = TRUE,
                  #                                     weight.reduction = sample.seurat[["pca"]],
                  #                                     dims = 1:30)
                  #   sample.seurat[["predictions"]] <- predictions.assay
                  #   
                  #   DefaultAssay(sample.seurat) <- "predictions"
                  #   
                  #   SpatialFeaturePlot(sample.seurat, 
                  #                      features = rownames(predictions.assay@data), 
                  #                      pt.size.factor = 1.6,
                  #                      #ncol = 4,
                  #                      crop = FALSE)
                  #   
                  #   ggsave(filename = paste0(seurat_results, "Labeltransfer_plot_", i, ".jpg"),
                  #          plot = last_plot(),
                  #          units = c("in"),
                  #          width = 11,
                  #          height = 8)
                  #   
                  #   
                  #   
                  #   # sample.seurat <- FindSpatiallyVariableFeatures(sample.seurat,
                  #   #                                                assay = "predictions",
                  #   #                                                selection.method = "moransi",
                  #   #                                                features = rownames(sample.seurat),
                  #   #                                                r.metric = 5,
                  #   #                                                slot = "data")
                  #   
                  #   # top.clusters <- head(SpatiallyVariableFeatures(sample.seurat, 
                  #   #                                                selection.method = "moransi"), 4)
                  #   # SpatialPlot(object = sample.seurat, 
                  #   #             features = top.clusters,
                  #   #             ncol = 2)
                  # }
                  
                  #!/usr/bin/env Rscript
                  
                  # Read and store variables from command line interface (CLI)
                  cli <- base::commandArgs(trailingOnly = TRUE) 
                  args <- base::strsplit(x = cli, split = "=", fixed = TRUE)
                  
                  for (e in args){
                    argname <- e[1]
                    argval <- e[2]
                    assign(argname, argval)
                  }
                  
                  # NOTE: All variables and functions are defined within the file below
                  source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")
                  
                  # NOTE: In single cell experiment, we combine multiple samples into a single 
                  # seurat object and do analysis. In spatial experiment, we analyze each sample
                  # individually.
                  
                  # NOTE: In visium v1 spatial transcriptomics, there are 4992 barcoded spots 
                  # which capture RNA from tissue placed above these spots. 
                  # So, unlike scRNASeq, each barcode corresponds to a spot, not a cell.
                  # https://www.youtube.com/watch?v=VwNk4d-0RJc
                  # https://kb.10xgenomics.com/hc/en-us/articles/360035848191-How-many-spots-are-within-a-single-capture-area-on-the-Visium-v1-Spatial-Gene-Expression-Slide
                  
                  # NOTE: Reads are obtained from both background (non-tissue) and tissue specific
                  # areas and stored in raw_feature_bc_matrix and filtered_fearure_bc_matrix
                  # respectively. Since, we are ONLY interested in reads obtained from tissues,
                  # use ONLY filtered_feature_bc_matrix.
                  
                  #******************************************************************************#
                  #                       STEP 1: SETUP THE SEURAT OBJECT                        #
                  #******************************************************************************#
                  
                  #*********************IMPORTING GEX (GENE EXPRESSION) DATA*********************#
                  
                  # raw_matrix_path
                  # ->Sample1
                  #   ->binned_outputs
                  #     ->square_002um
                  #       ->raw_feature_bc_matrix.h5
                  #       ->filtered_feature_bc_matrix.h5
                  #       ->spatial
                  #     ->square_008um
                  #       ->raw_feature_bc_matrix.h5
                  #       ->filtered_feature_bc_matrix.h5
                  #       ->spatial
                  #     ->square_016um
                  #       ->raw_feature_bc_matrix.h5
                  #       ->filtered_feature_bc_matrix.h5
                  #       ->spatial
                  
                  # DIRECTORY STRUCTURE IS IMPORTANT:
                  # data.dir MUST have a H5 file specified by filename parameter as well as folder
                  # named "spatial" containing the image
                  
                  # Create a list of sample names which will be added to each barcode.
                  # Since folder names correspond to sample name, we just use list.files()
                  samples <- list.files(path = filt_matrix_path)
                  
                  # Loop through each of the individual folders in parent directory & import data
                  for(i in samples){
                    
                    sample.seurat <- Load10X_Spatial(data.dir = paste0(filt_matrix_path, i),
                                                     filename = "filtered_feature_bc_matrix.h5",
                                                     assay = "Spatial",
                                                     slice = i,
                                                     bin.size = c(2,8,16),
                                                     filter.matrix = TRUE,
                                                     to.upper = FALSE,
                                                     image = NULL)
                    
                    # Since we use filtered_matrix which ONLY has reads from tissues, 
                    # filter.matrix = TRUE or FALSE didnt show any difference in sample.seurat
                    
                    # [i think wrong] If filter.matrix is set to FALSE, all 4992 spots will be recorded in 
                    # sample.seurat@images$B8@coordinates. Else, only spots that are over tissue,
                    # will be recorded in sample.seurat@images$B8@coordinates.
                    
                    # Unlike Seurat::Read10X(), Seurat::Load10X_Spatial doesnt have ability to 
                    # specify project parameter. So, we manually do it.
                    sample.seurat@meta.data <- sample.seurat@meta.data %>% 
                      dplyr::mutate(orig.ident = i)
                    
                    # Assign the seurat object to its corresponding variable
                    assign(paste0(i, ".filt"), sample.seurat)
                    cat("DATA IMPORTED FOR ", i, ".filt dataset\n")
                  }
                  
                  #******************************************************************************#
                  #                           STEP 2: QUALITY CONTROL                            #
                  #******************************************************************************#
                  
                  #**********************STEP 2A: CALCULATE ALL QC METRICS***********************#
                  
                  # NOTE:  We are going to do the same QC on each sample. So, we can merge the 
                  # individual seurat objects into a single seurat object and calculate metrics on
                  # the merged seurat object. However, if the number of cells exceeds 1.5 million,
                  # then this step will fail as there are simply too many cells. So, it is better
                  # to calculate QC metrics on individual samples, perform filtering on individual
                  # samples and then merge the filtered seurat objects which have fewer cells.
                  
                  # Initialize an empty dataframe where the class of each column resembles those 
                  # of raw_metadata. We will use this dataframe for making QC plots later.
                  raw_metadata <- data.frame(Cell = c(""), 
                                             Sample = as.factor(1), 
                                             nUMIs = c(0), 
                                             nGenes = c(0), 
                                             MitoRatio = c(0), 
                                             RiboRatio = c(0), 
                                             Novelty = c(0))
                  
                  # Calculate QC metrics for each sample individually using raw matrices
                  for(i in paste0(samples, ".filt")){
                    
                    sample.seurat <- get(i)
                    
                    # NOTE: PercentageFeatureSet() ONLY calculates metrics for cells within the 
                    # indicated assay. If no assay is indicated, it will use default assay. So, we
                    # loop through each assay and calculate metrics
                    
                    # Get list of assays
                    assay.list <- Assays(sample.seurat)
                    
                    for (assay in assay.list){
                      # Compute percent mito percent
                      sample.seurat <- Seurat::PercentageFeatureSet(object = sample.seurat,
                                                                    pattern = "^[Mm][Tt]-",
                                                                    features = NULL,
                                                                    col.name = "MitoPercent",
                                                                    assay = assay)
                      
                      # Compute percent ribo percent
                      sample.seurat <- Seurat::PercentageFeatureSet(object = sample.seurat,
                                                                    pattern = "^[Rr][Pp][SsLl]", 
                                                                    features = NULL,
                                                                    col.name = "RiboPercent",
                                                                    assay = assay)
                    }
                    
                    # Extract metadata
                    sample_metadata <- sample.seurat@meta.data
                    
                    # Rename columns to be more intuitive and add the additional QC metrics:
                    # (i)     Cell      : unique identifiers corresponding to each cell i.e. barcodes
                    # (ii)    Sample    : sample names
                    # (iii)   nUMIs     : number of transcripts per cell
                    # (iv)    nGenes    : number of genes per cell
                    # (v)     nHTO_UMIs : number of HTO reads per cell
                    # (vi)    nHTOs     : number of HTO types per cell
                    # (vii)   MitoRatio : MitoPercent/100
                    # (viii)	RiboRatio : RiboPercent/100  
                    # (ix)    Novelty   : log ratio of genes per UMI
                    sample_metadata <- sample_metadata %>% 
                      dplyr::mutate(Cell = paste0(orig.ident, "_", rownames(sample_metadata)),
                                    Sample = orig.ident,
                                    nUMIs = dplyr::case_when(!is.na(nFeature_Spatial.002um) & !is.na(nCount_Spatial.002um) ~ nCount_Spatial.002um,
                                                             !is.na(nFeature_Spatial.008um) & !is.na(nCount_Spatial.008um) ~ nCount_Spatial.008um,
                                                             !is.na(nFeature_Spatial.016um) & !is.na(nCount_Spatial.016um) ~ nCount_Spatial.016um,
                                                             TRUE ~ NA),
                                    nGenes = dplyr::case_when(!is.na(nFeature_Spatial.002um) & !is.na(nCount_Spatial.002um) ~ nFeature_Spatial.002um,
                                                              !is.na(nFeature_Spatial.008um) & !is.na(nCount_Spatial.008um) ~ nFeature_Spatial.008um,
                                                              !is.na(nFeature_Spatial.016um) & !is.na(nCount_Spatial.016um) ~ nFeature_Spatial.016um,
                                                              TRUE ~ NA),
                                    MitoRatio = MitoPercent/100,
                                    RiboRatio = RiboPercent/100,
                                    Novelty = log10(nGenes)/log10(nUMIs)) %>%
                      dplyr::select(Cell, Sample, contains(c("nFeature", "nCount")), nUMIs, nGenes, MitoRatio, RiboRatio, Novelty)
                    
                    # Replace the metadata in raw Seurat object
                    sample.seurat@meta.data <- sample_metadata
                    
                    # Append raw metadata of each seurat object which will be used for QC plots later
                    raw_metadata <- dplyr::bind_rows(raw_metadata, sample_metadata)
                    
                    # Assign the seurat object to its corresponding variable
                    assign(i, sample.seurat)
                  }
                  
                  #******************************STEP 2B: PERFORM QC*****************************#
                  
                  # In single cell analysis, we have to perform QC to remove bad cells in order to
                  # accurately annotate each cell. In spatial analysis, RNA from tissue is 
                  # captured within the spot below the tissue. Since sequencing depth is usually 
                  # way lower than single cell, we use almost all reads (lenient cutoffs) captured
                  # in the spot to annotate the spots.
                  
                  # NOTE: Unlike PercentageFeatureSet(), subset() uses all the cells irrespective 
                  # of the default assay. 2um bins are too small to contain 100 UMIs or 50 genes 
                  # and the entire Spatial.002um assay will be removed. So, we use lenient cutoffs
                  
                  # Perform QC for each sample individually
                  for(i in paste0(samples, ".filt")){
                    
                    sample.seurat <- get(i)
                    
                    # If default assay is 2um and no cells of default assay pass QC, then Seurat 
                    # will give error. So, change default assay to 8um
                    DefaultAssay(sample.seurat) <- "Spatial.008um"
                    
                    gene_cutoff <- 1           # reduced to 1 from 250 used for spatial
                    umi_cutoff <- 1            # reduced to 1 from 500 used for spatial
                    #mito_cutoff <- 0.2
                    #ribo_cutoff <- 0.05
                    #novelty_cutoff <- 0.8  	    # use 0.8 as starting point. Maximum 0.9
                    
                    sample.seurat <- base::subset(x = sample.seurat,
                                                  subset = (nGenes >= gene_cutoff) &
                                                    (nUMIs >= umi_cutoff))
                    #(MitoRatio <= mito_cutoff) &
                    # (RiboRatio >= ribo_cutoff) &
                    # (Novelty >= novelty_cutoff))
                    
                    # Assign the seurat object to its corresponding variable
                    assign(i, sample.seurat)
                  }
                  
                  # Create a merged Seurat object
                  # NOTE: Samples will have same barcodes. To keep track of cell identities 
                  # (i.e. barcodes) coming from each sample after merging, we add a prefix 
                  # (i.e. sample name) to each barcode using "add.cell.ids"
                  
                  # IMPORTANT: We DO NOT merge samples in spatial experiments. We analyze each
                  # sample individually.
                  
                  # NOTE: spatial objects have multiple assays. add.cell.ids modifies rownames
                  # ONLY for default assay. So, sample name is not added to cells from other
                  # assays. This creates problem in identifying common_bc. So, we exclude
                  # add.cell.ids
                  
                  #****************************STEP 2C: SAVE THE DATA****************************#
                  
                  for(i in paste0(samples, ".filt")){
                    
                    sample.seurat <- get(i)
                    saveRDS(sample.seurat, file=paste0(seurat_results, i,".rds"))
                  }
                  
                  #******************************************************************************#
                  #                       STEP 3: RUN THE STANDARD PIPELINE                      #
                  #******************************************************************************#
                  
                  # Create workbook to save markers
                  wb <- openxlsx::createWorkbook()
                  for (i in samples){
                    
                    object <- get(paste0(i, ".filt"))
                    
                    # Run workflow on all assays
                    for (assay in c("Spatial.002um", "Spatial.008um", "Spatial.016um")){
                      
                      DefaultAssay(object) <- assay
                      object <- Seurat::NormalizeData(       object = object, assay = assay, normalization.method = "LogNormalize")
                      object <- Seurat::FindVariableFeatures(object = object, assay = assay, nfeatures = 2000)
                      object <- Seurat::ScaleData(           object = object, assay = assay, features = VariableFeatures(object))
                      
                      # Create a new 'sketch' assay using 50k cells
                      sample_n <- 50000
                      e <- "Error"
                      while (class(e) == "character"){
                        e <- tryCatch(Seurat::SketchData(object = object, assay = assay,
                                                         ncells = sample_n,
                                                         method = "LeverageScore",
                                                         features = VariableFeatures(object),
                                                         sketched.assay = paste0(assay, ".sketch")), error = function(msg){
                                                           print("Reducing number of cells sampled")
                                                           return("Error")})
                        sample_n <- sample_n-5000
                      }
                      
                      object <- e
                      
                      # Switch analysis to sketched cells
                      DefaultAssay(object) <- paste0(assay, ".sketch")
                      object <- Seurat::FindVariableFeatures(object = object, 
                                                             assay = paste0(assay, ".sketch"), 
                                                             nfeatures = 2000)
                      
                      object <- Seurat::ScaleData(object = object, 
                                                  assay = paste0(assay, ".sketch"), 
                                                  features = VariableFeatures(object))
                      
                      object <- Seurat::RunPCA(object = object, 
                                               assay = paste0(assay, ".sketch"), 
                                               reduction.name = paste0(assay, ".pca.sketch"))
                      
                      object <- Seurat::FindNeighbors(object = object, 
                                                      assay = paste0(assay, ".sketch"), 
                                                      reduction = paste0(assay, ".pca.sketch"), 
                                                      dims = 1:50)
                      
                      object <- Seurat::FindClusters(object = object, 
                                                     algorithm = 4,
                                                     cluster.name = paste0(assay, ".cluster.sketch"), 
                                                     resolution = 0.6)
                      
                      object <- Seurat::RunUMAP(object = object, 
                                                reduction.name = paste0(assay, ".umap.sketch"), 
                                                reduction = paste0(assay, ".pca.sketch"), 
                                                dims = 1:50)
                      
                      # Project the cluster labels, dimensional reductions (PCA and UMAP) that we
                      # learned from the 50,000 sketched cells to the entire dataset
                      object <- Seurat::ProjectData(object = object, 
                                                    assay = assay,
                                                    sketched.assay = paste0(assay, ".sketch"),
                                                    sketched.reduction = paste0(assay, ".pca.sketch"),
                                                    full.reduction = paste0(assay, ".pca.full"),
                                                    dims = 1:50,
                                                    refdata = list(seurat_cluster.projected = paste0(assay, ".cluster.sketch")),
                                                    umap.model = paste0(assay, ".umap.sketch"))
                      
                      object <- Seurat::RunUMAP(object = object, 
                                                reduction.name = paste0(assay, ".umap.full"), 
                                                reduction = paste0(assay, ".pca.full"), 
                                                dims = 1:50)
                      
                      # Append assay name to the newly generated metadata columns so they dont get overwritten
                      object@meta.data <- object@meta.data %>%
                        dplyr::rename(!!rlang::sym(paste0(assay, ".cluster.full")) := "seurat_cluster.projected",
                                      !!rlang::sym(paste0(assay, ".cluster.full.score")) := "seurat_cluster.projected.score")
                      
                      # We can visualize the clustering results for the sketched cells, as well as the
                      # projected clustering results for the full dataset:
                      DefaultAssay(object) <- paste0(assay, ".sketch")
                      Idents(object) <- paste0(assay, ".cluster.sketch")
                      p1 <- DimPlot(object, reduction = paste0(assay, ".umap.sketch"), label = T) + 
                        ggtitle("Sketched clustering") + 
                        theme(legend.position = "none")
                      
                      # switch to full dataset
                      DefaultAssay(object) <- assay
                      Idents(object) <- paste0(assay, ".cluster.full")
                      p2 <- DimPlot(object, reduction = paste0(assay, ".umap.full"), label = T) + 
                        ggtitle("Projected clustering (full dataset)") + 
                        theme(legend.position = "none")
                      
                      # Save the plot
                      ggplot2::ggsave(filename=paste0(assay, ".umap.tiff"),
                                      plot=p1+p2,
                                      device="jpeg",
                                      path=diagnostics_path,
                                      width=11,
                                      height=8.5,
                                      units=c("in"),
                                      dpi=600,
                                      limitsize=TRUE,
                                      bg="white")
                      
                      # Find markers
                      markers <- Seurat::FindAllMarkers(object = object, 
                                                        assay = assay,
                                                        only.pos = TRUE)
                      
                      top10 <- markers %>%
                        dplyr::mutate(cluster = as.numeric(as.character(cluster))) %>%
                        dplyr::group_by(cluster) %>%
                        dplyr::filter(avg_log2FC > 1, p_val_adj < 0.05) %>%
                        slice_max(order_by=avg_log2FC*pct.1, n = 10) %>%
                        ungroup() %>%
                        data.frame() %>%
                        dplyr::arrange(cluster, desc(avg_log2FC*pct.1))
                      
                      # Save markers to worksheet
                      openxlsx::addWorksheet(wb, sheetName = paste0(i, assay))
                      openxlsx::writeData(wb, sheet = paste0(i, assay), x = markers, rowNames = FALSE)
                      openxlsx::addWorksheet(wb, sheetName = paste0(i, assay, "top10"))
                      openxlsx::writeData(wb, sheet = paste0(i, assay, "top10"), x = top10, rowNames = FALSE)
                      
                    }
                    
                    # Save xlsx file with markers
                    openxlsx::saveWorkbook(wb, file = paste0(seurat_results, i, "_Markers.xlsx"), overwrite = TRUE)
                    
                    # Save the object
                    saveRDS(object, file=paste0(seurat_results, i,".rds"))
                  }
                  
                  # DefaultAssay(object) <- "Spatial.008um"
                  # 
                  # vln.plot <- VlnPlot(object, features = "nCount_Spatial.008um", pt.size = 0) + 
                  #   theme(axis.text = element_text(size = 4)) + 
                  #   NoLegend()
                  # count.plot <- SpatialFeaturePlot(object, features = "nCount_Spatial.008um") + 
                  #   theme(legend.position = "right")
                  
                  # # switch back to 8um
                  # #object@assays$Spatial.008um@features@dimnames[[1]][25:30]
                  # DefaultAssay(object) <- "Spatial.008um"
                  # p2 <- SpatialFeaturePlot(object, features = "GAPDH") + ggtitle("Hpca expression (8um)")
                  
                  
                  #   # Visualize the unsupervised clusters based on their spatial location. 
                  #   SpatialDimPlot(object, label = T, repel = T, label.size = 4)
                  #   
                  #   # Plot the spatial location of different clusters individually.
                  #   Idents(object) <- "seurat_cluster.projected"
                  #   cells <- CellsByIdentities(object, idents = c(0, 4, 32, 34, 35))
                  #   p <- SpatialDimPlot(object,
                  #                       cells.highlight = cells[setdiff(names(cells), "NA")],
                  #                       cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
                  #   ) + NoLegend()
                  #   p
                  #   
                  #   
                  #   # visualize the top gene expression markers for each cluster:
                  #   # Create downsampled object to make visualization either
                  #   DefaultAssay(object) <- "Spatial.008um"
                  #   Idents(object) <- "seurat_cluster.projected"
                  #   object_subset <- subset(object, cells = Cells(object[["Spatial.008um"]]), downsample = 1000)
                  #   
                  #   
                  #   
                  #   object_subset <- ScaleData(object_subset, assay = "Spatial.008um", features = top5$gene)
                  #   p <- DoHeatmap(object_subset, assay = "Spatial.008um", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
                  #   p
                  #   
                  #   # Identifying spatially-defined tissue domains
                  #   object <- Banksy::RunBanksy(object,
                  #                               lambda = 0.8, verbose = TRUE,
                  #                               assay = "Spatial.008um", slot = "data", features = "variable",
                  #                               k_geom = 50)
                  #   
                  #   DefaultAssay(object) <- "BANKSY"
                  #   object <- RunPCA(object, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(object), npcs = 30)
                  #   object <- FindNeighbors(object, reduction = "pca.banksy", dims = 1:30)
                  #   object <- FindClusters(object, cluster.name = "banksy_cluster", resolution = 0.5)
                  #   
                  #   Idents(object) <- "banksy_cluster"
                  #   p <- SpatialDimPlot(object, group.by = "banksy_cluster", label = T, repel = T, label.size = 4)
                  #   p
                  #   
                  #   # highlight the spatial location of each tissue domain individually
                  #   banksy_cells <- CellsByIdentities(object)
                  #   p <- SpatialDimPlot(object, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T) + NoLegend()
                  #   p
                  #   
                  # }
                  # 
                  
                  # 
                  # # Load rds file of seurat objects
                  # for (i in samples){
                  #   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
                  #   sample.seurat <- sctransform_spatial_data(sample.seurat)
                  #   sample.seurat <- cluster_spatial_data(sample.seurat)
                  #   saveRDS(sample.seurat, file=paste0(seurat_results, i,".rds"))
                  #   
                  #   DefaultAssay(sample.seurat) <- "SCT"
                  #   p1 <- DimPlot(sample.seurat, reduction = "umap", label = TRUE)
                  #   p2 <- SpatialDimPlot(sample.seurat, label = TRUE, label.size = 3)
                  #   p <- p1 + p2
                  #   ggsave(filename = paste0(seurat_results, "Clusters_on_slide_", i, ".jpg"), 
                  #          plot = p)
                  # }
                  # 
                  # # # There is no integration like scRNA Seq. We analyse each slide individually.
                  # # # Unfortunately, all images are stored within each sample. So, we remove 
                  # # # unwanted images from each sample. If more than 1 image is present in each 
                  # # # sample, SpatialDimPlot() will give error.
                  # # integ_data <- sct_data
                  # # for (i in 1:length(sct_data)){
                  # #   integ_data[[i]] <- cluster_data(sct_data[[i]])
                  # #   integ_data[[i]]@images <- integ_data[[i]]@images[names(integ_data[[i]]@images) == names(sct_data)[[i]]]
                  # # }
                  # 
                  # # Color spots were GOI are present based on expression
                  # for (i in samples){
                  #   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
                  #   DefaultAssay(sample.seurat) <- "SCT"
                  #   Seurat::SpatialFeaturePlot(object = sample.seurat, 
                  #                              features = c("CD8A", "CD8B", "NPEPPS", "CDH12"),
                  #                              ncol = 4,
                  #                              slot = "data")
                  #   ggsave(filename = paste0(seurat_results, "Feature_plot_", i, ".jpg"),
                  #          plot = last_plot(),
                  #          units = c("in"),
                  #          width = 11,
                  #          height = 8)
                  #   
                  #   # This can ONLY plot 2 genes at a time
                  #   SpatialFeaturePlotBlend(sample.seurat, "CD8A", "NPEPPS")
                  # }
                  # 
                  # # Color spots were CD8 and NPEPPS are expressed in same plot. 
                  # # NOTE: This is NOT expression based. We just color the cells that express our 
                  # # GOI in different colors.
                  # # Use SCT assay to get counts, not Spatial assay as it has too much background.
                  # for (i in samples){
                  #   
                  #   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
                  #   
                  #   # Make sure the genes you want are present in the assay
                  #   GOI <- intersect(c("CD8A", "CD8B"), rownames(sample.seurat@assays$SCT@data))
                  #   cd8_df <- sample.seurat@assays$SCT$counts[GOI, ]
                  #   
                  #   if (length(GOI) > 1){
                  #     cd8_cells <- colnames(cd8_df[,colSums(cd8_df) > 0])
                  #   } else {
                  #     cd8_cells <- names(cd8_df[(cd8_df>0)])
                  #   }
                  #   
                  #   GOI <- intersect(c("NPEPPS"), rownames(sample.seurat@assays$SCT@data))
                  #   npepps_df <- sample.seurat@assays$SCT$counts[c("NPEPPS"), ]
                  #   npepps_cells <- names(npepps_df[(npepps_df>0)])
                  #   
                  #   SpatialPlot(object = sample.seurat, 
                  #               cells.highlight = list(CD8 = cd8_cells, NPEPPS = npepps_cells), 
                  #               cols.highlight =  c("green", "blue", "grey"),
                  #               pt.size.factor = 2)
                  #   
                  #   ggsave(filename = paste0(seurat_results, "Location_plot_", i, ".jpg"),
                  #          plot = last_plot(),
                  #          units = c("in"),
                  #          width = 11,
                  #          height = 8)
                  # }
                  # 
                  # # # Classify each spot as CD8+ or CD8-
                  # # #for (i in samples){
                  # # i <- "B8"
                  # #   
                  # #   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
                  # #   #rownames(sample.seurat@meta.data) <- sample.seurat@meta.data$Cell
                  # #   
                  # #   df <- sample.seurat@assays$Spatial$counts[c("CD8A", "CD8B"), ]
                  # #   pos_cells  <- colnames(df[,colSums(df) > 0])
                  # #   neg_cells <- setdiff(colnames(df), pos_cells)
                  # #   pos_cells <- paste0(i,"_", pos_cells)
                  # #   neg_cells <- paste0(i,"_", neg_cells)
                  # #   sample.seurat@meta.data <- sample.seurat@meta.data %>% 
                  # #     dplyr::mutate(CD8_status = dplyr::case_when(Cell %in% pos_cells ~ "CD8_pos",
                  # #                                                 TRUE ~ "CD8_neg"))
                  # # }
                  # 
                  # 
                  # i <- "B8"
                  # sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
                  # 
                  # a <- GetTissueCoordinates(sample.seurat)
                  # df <- data.frame(sample.seurat@images$B8@coordinates)
                  # 
                  # ### CHECK THAT WE ARE SUCCESSFULLY ABLE TO IDENTIFY NEIGHBORS ###
                  # # Find a cell in the center of the tissue
                  # row_cell <- median(df$row)
                  # col_cell <- median(df$col)
                  # cell <- rownames(df %>% dplyr::filter(row==row_cell, col==col_cell))
                  # 
                  # # Find adjacent cells
                  # cell1 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell-2)))
                  # cell2 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell+2)))
                  # cell3 <- rownames(df %>% dplyr::filter(row==(row_cell-2), col==col_cell))
                  # cell4 <- rownames(df %>% dplyr::filter(row==(row_cell+2), col==col_cell))
                  # cell5 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell-1)))
                  # cell6 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell+1)))
                  # cell7 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell-1)))
                  # cell8 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell+1)))
                  # 
                  # SpatialPlot(object = sample.seurat, 
                  #             cells.highlight = list(cell, unlist(lapply(paste0("cell", seq(1,8)), get), use.names=FALSE)), 
                  #             #facet.highlight = TRUE,
                  #             cols.highlight =  c("yellow", "red", "black"), 
                  #             pt.size.factor = 2)
                  # 
                  # #******************************************************************************#
                  # 
                  # i <- "B8"
                  # sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
                  # df <- data.frame(sample.seurat@images$B8@coordinates)
                  # expr_df <- sample.seurat@assays$SCT$data[rownames(sample.seurat@assays$SCT$data) %in% c("NPEPPS", "CD8A", "CD8B"),]
                  # 
                  # # Identify cells that express NPEPPS
                  # npepps <- expr_df["NPEPPS",]
                  # npepps <- npepps[npepps > 0]
                  # npepps <- names(npepps)
                  # # Not all npepps cells have co-ordinates. Remove those that lack co-ordinates.
                  # npepps <- intersect(npepps, rownames(df))
                  # 
                  # expr_df <- data.frame(expr_df)
                  # # For every cell that expresses NPEPPS, calculate total expression of CD8A and 
                  # # CD8B in 1st neighbor
                  # t_expr <- c()
                  # npepps_expr <- c()
                  # for (i in npepps){
                  #   row_cell <- df[rownames(df) == i,]$row
                  #   col_cell <- df[rownames(df) == i,]$col
                  #   
                  #   # Find adjacent cells
                  #   cell1 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell-2)))
                  #   cell2 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell+2)))
                  #   cell3 <- rownames(df %>% dplyr::filter(row==(row_cell-2), col==col_cell))
                  #   cell4 <- rownames(df %>% dplyr::filter(row==(row_cell+2), col==col_cell))
                  #   cell5 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell-1)))
                  #   cell6 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell+1)))
                  #   cell7 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell-1)))
                  #   cell8 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell+1)))
                  #   
                  #   # SpatialPlot(object = sample.seurat, 
                  #   #             cells.highlight = list(i, unlist(lapply(paste0("cell", seq(1,8)), get), use.names=FALSE)), 
                  #   #             #facet.highlight = TRUE,
                  #   #             cols.highlight =  c("yellow", "red", "black"), 
                  #   #             pt.size.factor = 2)
                  #   
                  #   cells <- intersect(make.names(unlist(lapply(paste0("cell", seq(1,8)), get))), colnames(expr_df))
                  #   print(cells)
                  #   
                  #   t <- expr_df %>% 
                  #     dplyr::select(all_of(cells)) %>%
                  #     tibble::rownames_to_column("Gene") %>%
                  #     dplyr::filter(Gene != "NPEPPS") %>%
                  #     tibble::column_to_rownames("Gene")
                  #   
                  #   n <- expr_df %>% 
                  #     dplyr::select(all_of(cells)) %>%
                  #     tibble::rownames_to_column("Gene") %>%
                  #     dplyr::filter(Gene == "NPEPPS") %>%
                  #     tibble::column_to_rownames("Gene")
                  #   
                  #   t_expr <- c(t_expr, sum(t))
                  #   npepps_expr <- c(npepps_expr, sum(n))
                  # }
                  # 
                  # df <- data.frame(npepps_expr, t_expr)
                  # # Save batch corrected normalized counts for entire dataset
                  # wb <- openxlsx::createWorkbook()
                  # openxlsx::addWorksheet(wb, sheetName = "spatial")
                  # openxlsx::writeData(wb, sheet = "spatial", x = df, rowNames = FALSE)
                  # openxlsx::saveWorkbook(wb,
                  #                        file = paste0(seurat_results, "Correlation.xlsx"),
                  #                        overwrite = TRUE)
                  # 
                  
                  # 
                  # # pseudobulk the counts based on donor-condition-celltype
                  # pseudo_ifnb <- AggregateExpression(ifnb, assays = "RNA", return.seurat = T, group.by = c("stim", "donor_id", "seurat_annotations"))
                  # 
                  # # each 'cell' is a donor-condition-celltype pseudobulk profile
                  # tail(Cells(pseudo_ifnb))
                  # 
                  # # Identification of Spatially Variable Features
                  # 
                  # # Seurat offers two workflows to identify molecular features that correlate 
                  # # with spatial location within a tissue. The first is to perform differential 
                  # # expression based on pre-annotated anatomical regions within the tissue, which
                  # # may be determined either from unsupervised clustering or prior knowledge. 
                  # # This strategy works will in this case, as the clusters above exhibit clear 
                  # # spatial restriction.
                  # 
                  # de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
                  # SpatialFeaturePlot(object = brain, 
                  #                    features = rownames(de_markers)[1:3], 
                  #                    alpha = c(0.1, 1), 
                  #                    ncol = 3)
                  # 
                  # 
                  # # An alternative approach, implemented in FindSpatiallyVariables(), is to 
                  # # search for features exhibiting spatial patterning in the absence of 
                  # # pre-annotation. The default method (method = 'markvariogram), is inspired by 
                  # # the Trendsceek, which models spatial transcriptomics data as a mark point 
                  # # process and computes a ‘variogram’, which identifies genes whose expression 
                  # # level is dependent on their spatial location. More specifically, this process
                  # # calculates gamma(r) values measuring the dependence between two spots a 
                  # # certain “r” distance apart. By default, we use an r-value of ‘5’ in these 
                  # # analyses, and only compute these values for variable genes (where variation 
                  # # is calculated independently of spatial location) to save time.
                  # 
                  # 
                  # brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000],
                  #                                        selection.method = "moransi")
                  # top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "moransi"), 6)
                  # SpatialFeaturePlot(brain, 
                  #                    features = top.features, 
                  #                    alpha = c(0.1, 1),
                  #                    ncol = 3)
                  # 
                  # 
                  # 
                  # # Integration with single-cell data i.e. label transfer from single-cell data
                  # # NOTE: At ~50um size, spots from the visium assay will encompass the expression 
                  # # profiles of multiple cells. 
                  # 
                  # # Users may be interested to ‘deconvolute’ each of the spatial voxels to predict
                  # # the underlying composition of cell types. We tested a wide variety of 
                  # # deconvolution and integration methods, using a reference scRNA-seq dataset of
                  # # ~14,000 adult mouse cortical cell taxonomy from the Allen Institute, 
                  # # generated with the SMART-Seq2 protocol. 
                  # 
                  # # We consistently found superior performance using integration methods (as 
                  # # opposed to deconvolution methods), likely because of substantially different 
                  # # noise models that characterize spatial and single-cell datasets, and 
                  # # integration methods are specifically designed to be robust to 
                  # # these differences. 
                  # 
                  # # We therefore apply the ‘anchor’-based integration workflow introduced in 
                  # # Seurat v3, that enables the probabilistic transfer of annotations from a 
                  # # reference to a query set. 
                  # 
                  # # While many of the methods are conserved (both procedures begin by identifying 
                  # # anchors), there are two important distinctions between data transfer and 
                  # # integration:
                  # # (i) In data transfer, Seurat does not correct or modify query expression data.
                  # # (ii) In data transfer, Seurat has an option (set by default) to project the 
                  # # PCA structure of a reference onto the query, instead of learning a joint 
                  # # structure with CCA. We generally suggest using this option when projecting 
                  # # data between scRNA-seq datasets.
                  # 
                  # # After finding anchors, we use the TransferData() function to classify the 
                  # # query cells based on reference data (a vector of reference cell type labels). 
                  # # TransferData() returns a matrix with predicted IDs and prediction scores, 
                  # # which we can add to the query metadata.
                  # 
                  # # NOTE: Make sure reference seurat object is SCTransformed.
                  # integrated_seurat <- readRDS("/hpc/home/kailasamms/scratch/scRNASeq_Chen/results_seurat/integrated_seurat_snn.rds")
                  # 
                  # # Seurat v3 vs Seurat v5 issues: Check if cells in graph are in same order
                  # all(Cells(integrated_seurat@graphs$integrated_snn) == colnames(integrated_seurat))
                  # # if FALSE, reorder the cells in all existing graphs
                  # integrated_seurat@graphs$integrated_snn <- 
                  #   integrated_seurat@graphs$integrated_snn[colnames(integrated_seurat), colnames(integrated_seurat)]
                  # integrated_seurat@graphs$integrated_nn <- 
                  #   integrated_seurat@graphs$integrated_nn[colnames(integrated_seurat), colnames(integrated_seurat)]
                  # 
                  # # Remove ambiguous cells before label transfer
                  # integrated_seurat <- subset(integrated_seurat, 
                  #                             sub_type == "Unclassified", 
                  #                             invert=TRUE)
                  # 
                  # # NOTE: For FindTransferAnchors() to work, both reference and query MUST be
                  # # SCtransformed()  https://github.com/satijalab/seurat/issues/3937
                  # # So, run SCTransform() on RNA assay.
                  # DefaultAssay(integrated_seurat) <- "integrated"
                  # 
                  # # Make sure all assays in Seuratv3 object has same order of cells
                  # acells <- colnames(x = integrated_seurat[["integrated"]])
                  # ocells <- colnames(x = integrated_seurat)
                  # 
                  # integrated_seurat@assays$RNA@counts <- integrated_seurat@assays$RNA@counts[,ocells]
                  # integrated_seurat@assays$RNA@data <- integrated_seurat@assays$RNA@data[,ocells]
                  # integrated_seurat@assays$integrated@data <- integrated_seurat@assays$integrated@data[,ocells]
                  # integrated_seurat@assays$integrated@scale.data <- integrated_seurat@assays$integrated@scale.data[,ocells]
                  # integrated_seurat@assays$SCT@counts <- integrated_seurat@assays$SCT@counts[,ocells]
                  # integrated_seurat@assays$SCT@data <- integrated_seurat@assays$SCT@data[,ocells]
                  # integrated_seurat@assays$SCT@scale.data <- integrated_seurat@assays$SCT@scale.data[,ocells]
                  # 
                  # # NOTE: If you subset the spatial object, re-run SCTranform() on it.
                  # for (i in samples){
                  #   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
                  #   DefaultAssay(sample.seurat) <- "SCT"
                  #   
                  #   # Run ONLY if you had subset the spatial data
                  #   # sample.seurat <- SCTransform(sample.seurat, 
                  #   #                              assay = "Spatial",
                  #   #                              verbose = FALSE)
                  #   # sample.seurat <- RunPCA(sample.seurat,
                  #   #                         verbose = FALSE)
                  #   
                  #   # Find anchors between reference and query
                  #   anchors <- FindTransferAnchors(reference = integrated_seurat,
                  #                                  query = sample.seurat,
                  #                                  normalization.method = "SCT",
                  #                                  reference.assay = "integrated",
                  #                                  query.assay = "SCT")
                  #   
                  #   predictions.assay <- TransferData(anchorset = anchors, 
                  #                                     refdata = integrated_seurat$cell_type,
                  #                                     prediction.assay = TRUE,
                  #                                     weight.reduction = sample.seurat[["pca"]],
                  #                                     dims = 1:30)
                  #   sample.seurat[["predictions"]] <- predictions.assay
                  #   
                  #   DefaultAssay(sample.seurat) <- "predictions"
                  #   
                  #   SpatialFeaturePlot(sample.seurat, 
                  #                      features = rownames(predictions.assay@data), 
                  #                      pt.size.factor = 1.6,
                  #                      #ncol = 4,
                  #                      crop = FALSE)
                  #   
                  #   ggsave(filename = paste0(seurat_results, "Labeltransfer_plot_", i, ".jpg"),
                  #          plot = last_plot(),
                  #          units = c("in"),
                  #          width = 11,
                  #          height = 8)
                  #   
                  #   
                  #   
                  #   # sample.seurat <- FindSpatiallyVariableFeatures(sample.seurat,
                  #   #                                                assay = "predictions",
                  #   #                                                selection.method = "moransi",
                  #   #                                                features = rownames(sample.seurat),
                  #   #                                                r.metric = 5,
                  #   #                                                slot = "data")
                  #   
                  #   # top.clusters <- head(SpatiallyVariableFeatures(sample.seurat, 
                  #   #                                                selection.method = "moransi"), 4)
                  #   # SpatialPlot(object = sample.seurat, 
                  #   #             features = top.clusters,
                  #   #             ncol = 2)
                  # }

