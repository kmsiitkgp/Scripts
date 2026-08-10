# 1. Load Libraries
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(RANN) # For fast distance math
library(patchwork)

plan("multisession", workers = 24) 

# Set a massive RAM limit for data transfer (50GB per worker)
options(future.globals.maxSize = 50000 * 1024^2)

# 1. SETUP PATHS
# ------------------------------------------------------------------------------
data_dir <- "/hpc/home/kailasamms/scratch/Xenium_Krizia/output-XETG00402__0088455__TMA__20260212__220101/"
sample_name <- "SampleA"

# 2. LOAD DATA
# ------------------------------------------------------------------------------
# Note: segmentations = NULL keeps the object light by excluding cell boundaries
xenium_obj <- Seurat::LoadXenium(data.dir = data_dir,
                                 fov = "fov",
                                 assay = "Xenium",
                                 mols.qv.threshold = 20,
                                 cell.centroids = TRUE,
                                 molecule.coordinates = TRUE,
                                 segmentations = NULL,  # Change to "segmentation" if you need cell polygons
                                 flip.xy = FALSE)

# 3. PRE-PROCESSING & GLOBAL FILTERING
# ------------------------------------------------------------------------------
# Add X/Y coordinates to metadata as a permanent reference
coords <- GetTissueCoordinates(xenium_obj)
xenium_obj$x_centroid <- coords$x
xenium_obj$y_centroid <- coords$y

# Global Filter: Remove "empty" cells (fewer than 5 transcripts)
# This ensures both Expression and QC files share the same high-quality cell list
xenium_obj <- subset(xenium_obj, subset = nCount_Xenium > 5)

# 4. CREATE EXPRESSION OBJECT ("Science" File)
# ------------------------------------------------------------------------------
spatial_expr_obj <- subset(xenium_obj, assays = "Xenium")

# Hard delete technical assays to free up RAM and clarify the object structure
spatial_expr_obj[["BlankCodeword"]]   <- NULL
spatial_expr_obj[["ControlCodeword"]] <- NULL
spatial_expr_obj[["ControlProbe"]]    <- NULL
spatial_expr_obj[["GenomicControl"]]  <- NULL

# 5. CREATE QUALITY CONTROL OBJECT ("Audit" File)
# ------------------------------------------------------------------------------
noise_assays <- c("ControlProbe", "BlankCodeword", 
                  "ControlCodeword", "GenomicControl")

qc_noise_obj <- subset(xenium_obj, assays = noise_assays)

# Remove the real gene assay and heavy image data to make the QC file tiny
DefaultAssay(qc_noise_obj) <- "ControlProbe"
qc_noise_obj[["Xenium"]] <- NULL
qc_noise_obj@images      <- list() 

# 6. EXPORT
# ------------------------------------------------------------------------------
saveRDS(spatial_expr_obj, paste0(sample_name, "_spatial_expression.rds"))
saveRDS(qc_noise_obj,     paste0(sample_name, "_quality_control.rds"))

message("Success: Processed ", ncol(spatial_expr_obj), " cells.")
message("Files saved: ", sample_name, "_spatial_expression.rds & ", sample_name, "_quality_control.rds")

options(future.globals.maxSize = 8000 * 1024^2) # Increase to 8GB per worker

# 4. Normalization and Dimensional Reduction
# Use vst.flavor = "v2" for the modern, faster version of SCT
xenium_obj <- spatial_expr_obj
xenium_obj <- SCTransform(xenium_obj, 
                          assay = "Xenium", 
                          vst.flavor = "v2", 
                          verbose = TRUE)

# PCA
xenium_obj <- RunPCA(xenium_obj, npcs = 30, verbose = FALSE)

# UMAP
# For massive datasets, 'uwot' is the engine; 'n.neighbors' at 30 helps with 
# global structure for 400k cells.
xenium_obj <- RunUMAP(xenium_obj, dims = 1:30, verbose = FALSE)

# 5. Clustering
xenium_obj <- FindNeighbors(xenium_obj, dims = 1:30, verbose = FALSE)
xenium_obj <- FindClusters(xenium_obj, resolution = 0.5, verbose = FALSE)

# 6. Save the Analyzed Object
saveRDS(xenium_obj, paste0(sample_name, "_analyzed.rds"))


# 1. Standard UMAP (Gene Space)
# This shows you how similar cells are based on their RNA
p1 <- DimPlot(xenium_obj, reduction = "umap", label = TRUE) + 
  ggtitle("UMAP: Transcriptional Clusters")

# 2. Spatial Plot (Physical Space)
# This shows you where those same clusters live on the slide
# 'mols.size = 0.01' or similar can help if the plot feels too "crowded"
p2 <- ImageDimPlot(xenium_obj, 
                   fov = "fov", 
                   cols = "polychrome", 
                   axes = TRUE, 
                   size = 0.75) + 
  ggtitle("Spatial: Cluster Distribution")

# Combine them side-by-side to compare
# 2. Save as a High-Res PNG (Best for presentations/Slack)
ggsave(filename = "SampleA_UMAP_Spatial_Clusters.png", 
       plot = p1 + p2, 
       width = 16,     # Wide enough for two plots side-by-side
       height = 8, 
       dpi = 300)     # High resolution for 373k dots

# 3. Save as a PDF (Best for zooming in on 373k dots without pixelation)
ggsave(filename = "SampleA_UMAP_Spatial_Clusters.pdf", 
       plot = p1 + p2, 
       width = 16, 
       height = 8)

# --- STOP HERE: You must annotate your clusters now ---
# Example: xenium_obj$cell_type <- "Tumor" (based on cluster markers)
# ------------------------------------------------------



### 6. Niche Analysis Example (Neighborhood based)
# This creates a NEW assay where "genes" are actually neighboring cell types.
# This helps identify which tumor cells are in a "Vascular Niche."

xenium_obj <- BuildNicheAssay(
  object = xenium_obj, 
  fov = "fov", 
  group.by = "cell_type", # Use the annotations you made
  neighbors.k = 20,       # Look at 20 closest cells
  niches.k = 5            # Group those environments into 5 types
)

# You can now see which "Niche" every cell belongs to
# Idents(xenium_obj) <- "niches"
# ImageDimPlot(xenium_obj, group.by = "niches")

### 7. Custom Distance Analysis (Tumor vs Blood Vessel)
# Let's say you want to manually compare tumor cells near vs far from vessels.

# A. Extract Cores
tumor_coords <- xenium_obj@meta.data %>% filter(cell_type == "Tumor") %>% select(x_centroid, y_centroid)
vessel_coords <- xenium_obj@meta.data %>% filter(cell_type == "Endothelial") %>% select(x_centroid, y_centroid)

# B. Calculate Distance
nn <- nn2(data = vessel_coords, query = tumor_coords, k = 1)
xenium_obj$dist_to_vessel <- NA
xenium_obj@meta.data[rownames(tumor_coords), "dist_to_vessel"] <- nn$nn.dists

# C. Differential Expression (DE) by Location
# We subset to ONLY tumor cells to compare them to each other
tumor_subset <- subset(xenium_obj, subset = cell_type == "Tumor")
tumor_subset$spatial_group <- ifelse(tumor_subset$dist_to_vessel < 20, "Near", "Far")

Idents(tumor_subset) <- "spatial_group"
spatial_markers <- FindMarkers(tumor_subset, ident.1 = "Near", ident.2 = "Far")

# 8. Visualization
# Show the expression of a top marker spatially
ImageFeaturePlot(xenium_obj, features = "VEGFA")

# ------

# =============================================================================

# Xenium Distance-Based Spatial Analysis

# =============================================================================

library(Seurat)
library(dplyr)
library(RANN)

# -----------------------------------------------------------------------------

# USER SETTINGS

# -----------------------------------------------------------------------------

tumor_label <- "Tumor"

target_cell_types <- c(
  "Endothelial",
  "CD8_T",
  "Treg",
  "Macrophage",
  "Fibroblast"
)

output_dir <- "Spatial_Distance_Results"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------

# CHECK CELL TYPES

# -----------------------------------------------------------------------------

cat("\nAvailable cell types:\n")
print(sort(unique(xenium_obj$cell_type)))

# -----------------------------------------------------------------------------

# MAIN LOOP

# -----------------------------------------------------------------------------

for(target_type in target_cell_types){
  
  cat("\n=====================================================\n")
  cat("Analyzing:", target_type, "\n")
  cat("=====================================================\n")
  
  meta <- xenium_obj@meta.data
           
  tumor_meta <- meta %>%
    filter(cell_type == tumor_label)
  
  target_meta <- meta %>%
    filter(cell_type == target_type)
  
  if(nrow(tumor_meta) < 10){
    cat("Skipping:", target_type,
        "- too few tumor cells\n")
    next
  }
  
  if(nrow(target_meta) < 10){
    cat("Skipping:", target_type,
        "- too few target cells\n")
    next
  }
  
  tumor_coords <- tumor_meta %>%
    select(x_centroid, y_centroid)
  
  target_coords <- target_meta %>%
    select(x_centroid, y_centroid)
  
  # ---------------------------------------------------------------------------
  # Calculate nearest-neighbor distance
  # ---------------------------------------------------------------------------
  
  nn <- nn2(
    data = as.matrix(target_coords),
    query = as.matrix(tumor_coords),
    k = 1
  )
  
  dist_col <- paste0(
    "dist_to_",
    gsub("[^A-Za-z0-9]", "_", target_type)
  )
  
  xenium_obj[[dist_col]] <- NA_real_
  
  [xenium_obj@meta.data](mailto:xenium_obj@meta.data)[
    rownames(tumor_meta),
    dist_col
  ] <- nn$nn.dists[,1]
  
  # ---------------------------------------------------------------------------
  # Create tumor-only object
  # ---------------------------------------------------------------------------
  
  tumor_subset <- subset(
    xenium_obj,
    subset = cell_type == tumor_label
  )
  
  distances <- tumor_subset[[dist_col]][,1]
  
  q25 <- quantile(
    distances,
    0.25,
    na.rm = TRUE
  )
  
  q75 <- quantile(
    distances,
    0.75,
    na.rm = TRUE
  )
  
  keep_cells <- rownames(
    tumor_subset@meta.data[distances <= q25 | distances >= q75,]
  )
  
  tumor_subset <- subset(
    tumor_subset,
    cells = keep_cells
  )
  
  tumor_subset$spatial_group <- ifelse(
    tumor_subset[[dist_col]][,1] <= q25,
    "Near",
    "Far"
  )
  
  Idents(tumor_subset) <- "spatial_group"
  
  # ---------------------------------------------------------------------------
  # Differential Expression
  # ---------------------------------------------------------------------------
  
  DefaultAssay(tumor_subset) <- "SCT"
  
  markers <- FindMarkers(
    tumor_subset,
    ident.1 = "Near",
    ident.2 = "Far",
    min.pct = 0.10,
    logfc.threshold = 0.25
  )
  
  markers$gene <- rownames(markers)
  
  write.csv(
    markers,
    file = file.path(
      output_dir,
      paste0(
        "Tumor_vs_",
        target_type,
        "_Distance_DE.csv"
      )
    ),
    row.names = FALSE
  )
  
  cat(
    "Completed:",
    target_type,
    "\n"
  )
}

# -----------------------------------------------------------------------------
# SAVE UPDATED OBJECT
# -----------------------------------------------------------------------------

saveRDS(
  xenium_obj,
  file.path(
    output_dir,
    "xenium_with_distance_metrics.rds"
  )
)

cat("\nAll analyses completed.\n")
cat("Results written to:", output_dir, "\n")
cat("Updated Seurat object saved.\n")
