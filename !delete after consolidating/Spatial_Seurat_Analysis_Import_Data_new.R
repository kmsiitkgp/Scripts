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
source("/hpc/home/kailasamms/projects/scRNASeq/Custom.Functions.R.Analysis.R")

# If integrating multiple samples, use below codes
# NOTE: Run workflow for each bin size separately
for (assay in c("Spatial.016um")){
  
  filt.obj <- readRDS(paste0(seurat_results, "filtered.seurat.", assay, ".rds"))
  
  # Perform SCTransformation at sample level, not slide level
  sct.obj <- sctransform_sc.sp(filt.obj, assay, seurat_results)
  
  # Integrate all samples belonging to specific bin size
  reference.samples <- NULL
  kweight <- min(sct.obj@meta.data %>% dplyr::count(Sample) %>% dplyr::select(n) %>% min()/2, 100) 
  integ <- integrate_sc.sp(sct.obj, assay, reference.samples, kweight, seurat_results)
  
  # Perform clustering at various resolution using different reductions
  integ.clust <- cluster_sc.sp(integ, assay, seurat_results)
  
  # Remove sparse clusters and save the results
  integ.final <- remove_sparse_clusters_sc.sp(integ.clust, assay, seurat_results)
  
  # Plot metrics post integration (DO this ONLY after adding sample info)
  plot_metrics_post_integration_sc.sp(integ.final, assay, diagnostics_path)
  
  # Plot spatial map of samples and groups
  plot.seurat <- Seurat::SplitObject(object = integ.final,  
                                     split.by = "Sample")
  for(x in 1:length(plot.seurat)){
    plot_spatial_map(plot.seurat[[x]], assay, diagnostics_path)
  }
  
  # Find markers
  for (res in c(0.4, 0.8)){  
    resolution <- res
    reduction <- "harmony"
    identify_markers_sc.sp(integ.final, assay, resolution, reduction, assay, seurat_results)
  }
}

# Decide which reduction gives (i) clear separation of clusters 
# (ii) clusters cells together i.e. B cells are not split between clusters.
# In all cases so far, RPCA and Harmony seemed to work best.

# Visualize by UMAP the contribution of each sample based on groups
for (assay in c("Spatial.008um", "Spatial.016um")){
  
  integ.final <- readRDS(paste0(seurat_results, "integrated.seurat.", assay, ".rds"))
  
  # File names, reductions, splits  for each of the figures
  groups <- unique(integ.final@meta.data$Group)
  filenames <- paste0("UMAP.Group.", groups, ".", assay)
  reductions <- rep(x="umap.harmony", times=length(filenames))
  splits <- rep(x="Sample", times=length(filenames))
  
  for (i in 1:length(groups)){ 
    
    plot.seurat <- subset(integ.final, Group == groups[i])
    
    plot.seurat <- Seurat::SplitObject(object = plot.seurat,
                                       split.by = splits[i])
    
    purrr::map(.x = c(1:length(plot.seurat)),
               .f = function(x){  
                 Idents(plot.seurat[[x]]) <- "cluster.0.4.harmony"
                 Seurat::DimPlot(object = plot.seurat[[x]],
                                 reduction = reductions[i],
                                 group.by = "cluster.0.4.harmony",
                                 pt.size = 0.1,
                                 order = TRUE,  # plot positive cells above negative cells
                                 label = TRUE,
                                 raster = FALSE,
                                 combine = TRUE) +
                   NoLegend() +
                   my_theme + 
                   ggplot2::labs(title = names(plot.seurat)[x]) 
               }) %>% cowplot::plot_grid(plotlist=.,
                                         align="hv",
                                         axis="tblr",
                                         nrow=ceiling(sqrt(length(plot.seurat))),
                                         ncol=floor(sqrt(length(plot.seurat))),
                                         rel_widths=1,
                                         rel_heights=1,
                                         greedy=TRUE,
                                         byrow=TRUE)
    # Save the plot
    ggplot2::ggsave(filename = paste0(filenames[i], ".tiff"),
                    plot = last_plot(),
                    device = "jpeg",
                    path = diagnostics_path,
                    scale = 1,
                    width = 4*floor(sqrt(length(plot.seurat))),
                    height = 4*ceiling(sqrt(length(plot.seurat))),
                    units = c("in"),
                    dpi = 600,
                    limitsize = TRUE,
                    bg = "white")
  }
}

#*****************************************************************************#



