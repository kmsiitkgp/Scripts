gmt.dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/GSEA_genesets"

# ---- LOAD PACKAGES ----

pkgs <- c(
  "BiocManager", "remotes", "AnnotationHub", "ensembldb", "org.Hs.eg.db",
  "org.Mm.eg.db", "fgsea", "clusterProfiler", "progeny", "dorothea", 
  "viper", "DESeq2", "sva", "GSVA", "RcisTarget", "glmGamPoi", "Seurat",
  "harmony", "hdf5r", "arrow", "leidenAlg", "scCustomize", "reticulate", 
  "ashr", "infercnv", "UCell", "scDblFinder", "DropletUtils", "batchelor", 
  "ClusterFoldSimilarity", "CellChat", "Banksy", "SeuratWrappers", "presto", 
  "DoubletFinder", "SeuratData", "oligo", "oligoData", "illuminaHumanv4.db", 
  "hgu133plus2.db", "GEOquery", "affy", "lumi", "openxlsx", "dplyr", "tibble", 
  "stringr", "purrr", "ggplot2", "ggplotify", "ggrepel", "ggpubr", "ggfortify", 
  "ggridges", "ggbeeswarm", "pheatmap", "VennDiagram", "survival", "survminer", 
  "UpSetR", "umap", "plot3D", "cowplot", "viridis", "RColorBrewer", "colorspace", 
  "enrichplot", "ComplexHeatmap", "NanoStringNCTools", "GeomxTools", 
  "GeoMxWorkflows", "networkD3"
)

for (pkg in pkgs) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    message(paste("Loaded", pkg))
  } else {
    message(paste("Package", pkg, "is not installed — skipping"))
  }
}

# NOTE: survminer handles %++% while dplyr handles %>%

# ---- Custom ggplot Theme ----

custom_theme <- ggplot2::theme(
  plot.title    = element_text(family = "sans", face = "bold",  colour = "black", size = 15, hjust = 0.5),
  legend.title  = element_text(family = "sans", face = "bold",  colour = "black", size = 12, hjust = 0,   vjust = 1, angle = 0),
  axis.title.x  = element_text(family = "sans", face = "bold",  colour = "black", size = 12, hjust = 0.5, vjust = 0, angle = 0),
  axis.title.y  = element_text(family = "sans", face = "bold",  colour = "black", size = 12, hjust = 0.5, vjust = 1, angle = 90),
  legend.text   = element_text(family = "sans", face = "plain", colour = "black", size = 10, hjust = 0.5),
  axis.text.x   = element_text(family = "sans", face = "plain", colour = "black", size = 10, hjust = 0.5, vjust = 0.5, angle = 45),
  axis.text.y   = element_text(family = "sans", face = "plain", colour = "black", size = 10, hjust = 0.5, vjust = 0.5, angle = 0)
)

# ---- Custom Palette ----

custom_palette <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173", "#5254a3", "#6b6ecf", "#9c9ede", "#cedb9c", "#8ca252",
  "#a55194", "#e5e56f", "#66a61e", "#e6ab02", "#a6761d", "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ffff33",
  "#f781bf", "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3", "#1e90ff",
  "#ff4500", "#32cd32", "#ff0000", "#8a2be2", "#a0522d", "#ff66cc", "#333333", "#adff2f", "#00ced1", "#ffd700",
  "#6699cc", "#cc6644", "#66aa66", "#cc6666", "#9966cc", "#996633", "#cc99cc", "#99cc44", "#66cccc", "#cccc66",
  "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
  "#ffffcc", "#e7ba52", "#ce6dbd", "#d6616b", "#b5cf6b", "#dbdb5c", "#e7cb94", "#ad494a", "#bd9e39", "#de9ed6",
  "#e7969c", "#cedb9c", "#33a02c", "#b2df8a", "#fdbf6f", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#8dd3c7",
  "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5",
  "#ffed6f", "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#999999"
)

tab_palettes <- c(
  # tab20
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#ffff33",
  "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#ffffcc",
  # tab20 (saturated)
  "#1e90ff", "#ff4500", "#32cd32", "#ff0000", "#8a2be2", "#a0522d", "#ff33ff", "#000000", "#adff2f", "#00ced1", "#ffd700",
  "#63b3ff", "#ff8052", "#66e066", "#ff6666", "#b199e8", "#cd8866", "#ff99ff", "#999999", "#d4ff7f", "#66f0f0", "#ffee66",
  # tab20 (muted)
  "#6699cc", "#cc6644", "#66aa66", "#cc6666", "#9966cc", "#996633", "#cc66cc", "#666666", "#99cc44", "#66cccc","#cccc66",
  "#99bbee", "#ddbbaa", "#aacc99", "#ee9999", "#bbaadd", "#ccaa88", "#ee99ee", "#bbbbbb", "#ddffaa", "#66dddd","#ffffaa",
  # tab20c (modified)
  "#393b79", "#e7ba52", "#637939", "#843c39", "#6b6ecf", "#8c6d31", "#ce6dbd", "#d6616b", "#b5cf6b", "#7b4173", "#dbdb5c",
  "#5254a3", "#e7cb94", "#8ca252", "#ad494a", "#9c9ede", "#bd9e39", "#de9ed6", "#e7969c", "#cedb9c", "#a55194", "#e5e56f"
)

# ---- BULK RNA ANALYSIS RELATED FUNCTIONS ----

# # User Override Project Directories & Parameters
# 
# proj <- ""
# species <- ""
# contrasts <- c("Treatment1-Control1",
#                "Treatment2-Control2")
# 
# parent.dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data"
# gmt.dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/GSEA_genesets"
# 
# # DESeq2 overrides
# deseq2.override <- list(
#   contrasts     = contrasts
#   #design        = "Comparisons",            # DESeq2 design formula or column name
#   #lfc.cutoff    = 0,                        # Log fold change cutoff for significance
#   #padj.cutoff   = 0.1,                      # Adjusted p-value cutoff for significance
#   #batch.correct = FALSE                     # Boolean, whether to apply batch correction
# )
# 
# # Heatmap overrides
# heatmap.override <- list(
#   #force.log        = TRUE,                  # Force log transformation
#   col.ann          = NULL,                  # Column annotation
#   #row.ann          = NULL,                  # Row annotation
#   col.gaps         = NULL,                  # Column gaps
#   #row.gaps         = NULL,                  # Row gaps
#   col.cluster      = "all",                 # Column clustering
#   #row.cluster      = "all",                 # Row clustering
#   #palette         = "rdbu",                # Heatmap palette
#   #ann.palette     = "discrete",            # Annotation palette
#   #border.color    = NA,                    # Cell border color
#   #show.expr.legend = TRUE,                  # Show expression legend
#   #title           = "",                    # Heatmap title
#   #format           = "tiff"                 # Output file format
# )
# 
# # Volcano plot overrides
# volcano.override <- list(
#   #lfc.cutoff  = 0.58,                         # Minimum log2 fold-change to highlight
#   #padj.cutoff = 0.05,                      # Adjusted p-value cutoff
#   #color       = "vrds",                    # Color palette
#   label.genes = c()                         # Genes to label on the plot
# )
# # Setup project
# proj.params <- setup_project(
#   proj = proj,
#   species = species,
#   contrasts = contrasts,
#   parent.dir = parent.dir,
#   gmt.dir = gmt.dir,
#   deseq2.override = deseq2.override,
#   heatmap.override = heatmap.override,
#   volcano.override = volcano.override
# )

# Default Project Directories & Parameters
setup_project <- function(proj, species, contrasts, parent.dir, gmt.dir,
                          deseq2.override  = list(),
                          heatmap.override = list(),
                          volcano.override = list(),
                          pathway.override = list()) {
  
  # Default DESeq2 Parameters
  default.deseq2 <- list(
    contrasts     = c("Treatment-Control"),   # Vector of contrasts for DE analysis
    design        = "Comparisons",            # DESeq2 design formula or column name
    lfc.cutoff    = 0,                        # Log fold change cutoff for significance
    padj.cutoff   = 0.1,                      # Adjusted p-value cutoff for significance
    batch.correct = FALSE                     # Boolean, whether to apply batch correction
  )
  
  # Default Heatmap Parameters
  default.heatmap <- list(
    force.log        = FALSE,       # Force log transform (default FALSE i.e. auto detect)
    col.ann          = NULL,        # Single/multiple columns from metadata_col for column annotation
    row.ann          = NULL,        # Single/multiple columns from metadata_row for row annotation
    col.gaps         = NULL,        # Single column from metadata_col to define column gaps in heatmap
    row.gaps         = NULL,        # Single column from metadata_row to define row gaps in heatmap
    col.cluster      = "all",       # Single column from metadata_col, "all", "alphabetical" for clustering columns
    row.cluster      = "all",       # Single column from metadata_row, "all", "alphabetical" for clustering columns
    palette          = "rdbu",      # Color palette for heatmap matrix ("rdbu" or "vrds")
    ann.palette      = "discrete",  # Color palette for heatmap annotation ("discrete" or "sequential")
    border.color     = NA,          # Color of heatmap cell borders (default NA i.e. no border)
    show.expr.legend = TRUE,        # Show expression legend (set FALSE if annotations overlap)
    title            = NA,          # Title for heatmap (default NA i.e. no title)
    format           = "tiff"       # Output format for heatmap ("tiff", "jpeg", "pdf")
  )
  
  # Default Volcano Parameters
  default.volcano <- list(
    lfc.cutoff   = 0.58,            # Log fold change cutoff for volcano plot
    padj.cutoff  = 0.05,            # Adjusted p-value cutoff for volcano plot
    color        = "vrds",          # Color palette for volcano plot ("vrds", etc.)
    label.genes  = NULL             # Optional vector of genes to label on volcano plot
  )
  
  # Merge Overrides
  deseq2.params  <- modifyList(default.deseq2,  deseq2.override)
  heatmap.params <- modifyList(default.heatmap, heatmap.override)
  volcano.params <- modifyList(default.volcano, volcano.override)
  
  # Derived Directories
  proj.dir      <- file.path(parent.dir, proj)
  counts.dir    <- file.path(proj.dir, "counts")                # Directory containing count files
  contrast.dir  <- file.path(proj.dir, contrasts)               # Directory to store results for each contrast
  deseq2.dir    <- file.path(contrast.dir, "DEG_Analysis")     # Directory to store DESeq2 results
  pathway.dir   <- file.path(contrast.dir, "Pathway_Analysis") # Directory to store Pathway analysis results
  tf.dir        <- file.path(contrast.dir, "TF_Analysis")      # Directory to store TF analysis results
  
  # Final Project Parameters
  proj.params <- list(
    proj.dir    = normalizePath(proj.dir, mustWork = FALSE),
    counts.dir  = normalizePath(counts.dir, mustWork = FALSE),
    gmt.dir     = normalizePath(gmt.dir, mustWork = FALSE),
    contrast.dir= normalizePath(contrast.dir, mustWork = FALSE),
    deseq2.dir  = normalizePath(file.path(contrast.dir, "DEG_Analysis"), mustWork = FALSE),
    pathway.dir = normalizePath(file.path(contrast.dir, "Pathway_Analysis"), mustWork = FALSE),
    tf.dir      = normalizePath(file.path(contrast.dir, "TF_Analysis"), mustWork = FALSE),
    
    proj        = proj,
    species     = species,
    
    deseq2      = deseq2.params,
    heatmap     = heatmap.params,
    volcano     = volcano.params
  )
  
  if (proj == "") warning("⚠️ Project name is empty")
  if (species == "") warning("⚠️ Species is not set")
  if (length(contrasts) == 0) warning("⚠️ No contrasts specified")
  
  return(proj.params)
}

# Helper functions #1
save_xlsx <- function(data, file, sheet_name, row_names) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = sheet_name)
  openxlsx::writeData(wb, sheet = sheet_name, x = data, rowNames = row_names)
  if(row_names){
    openxlsx::writeData(wb, sheet = sheet_name, x = "SYMBOL", startCol = 1, startRow = 1)
  }
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
}

# Helper functions #2
filter_samples_by_contrast <- function(metadata, contrast) {
  contrast <- base::gsub(pattern = "\\(|\\)", replacement = "", x = contrast)
  metadata %>%
    dplyr::filter(Comparisons %in% stringr::str_split(contrast, "-")[[1]]) %>%
    dplyr::pull(Sample_ID) %>%
    as.character()
}

merge_counts <- function(proj.params) {
  
  # Check required proj.params attributes ----
  required_attrs <- c("counts.dir", "proj.dir", "proj")
  for (attr in required_attrs) {
    if (is.null(proj.params[[attr]])) {
      stop("⚠️ Missing required proj.params attribute: ", attr)
    } else {
      message("✔ Found proj.params attribute: ", attr)
    }
  }
  
  set.seed(1234)
  
  # Get count files
  count_files <- list.files(path = proj.params$counts.dir, pattern = "\\.txt$|ReadsPerGene\\.out\\.tab$", full.names = TRUE)
  if (length( count_files) == 0) {
    stop("No count files found in the directory.")
  }
  
  # ---- Initialize ----
  all_counts <- list()
  gene_lists <- list()
  sample_ids <- character()
  
  # ---- Define special counters generated by HTSeq and STAR outputs ----
  special_counters <- c("__no_feature", "__ambiguous", "__too_low_aQual", 
                        "__not_aligned", "__alignment_not_unique", "__assignment_counts",
                        "N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
  
  # ---- Parse Files ----
  for (count_file in count_files) {
    
    # Read count file
    df <- tryCatch({
      read.table(file = count_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    }, error = function(e) {
      stop("Error reading file: ", count_file, " — ", e$message)
    })
    
    if (ncol(df) < 4) stop("count file does not have expected 4 columns: ", count_files[i])
    
    # Remove special counters
    df <- df %>% dplyr::filter(!(.[[1]] %in% special_counters))
    
    sample_id <- gsub("\\..*$|ReadsPerGene\\.out\\.tab", "", basename(count_file))
    gene_ids <- df[[1]]
    strand_sums <- colSums(df[2:4], na.rm = TRUE)
    
    # Determine strandedness
    if (abs((strand_sums[1]/strand_sums[2]) - (strand_sums[1]/strand_sums[3])) < 2) {
      message("Detected unstranded library for: ", sample_id)
      counts <- df[[2]]
    } else if (strand_sums[2] > 3 * strand_sums[3]) {
      message("Detected positively stranded library for: ", sample_id)
      counts <- df[[3]]
    } else if (strand_sums[3] > 3 * strand_sums[2]) {
      message("Detected negatively stranded library for: ", sample_id)
      counts <- df[[4]]
    } else {
      stop("Could not determine strandedness for: ", f)
    }
    
    all_counts[[sample_id]] <- counts
    gene_lists[[sample_id]] <- gene_ids
    sample_ids <- c(sample_ids, sample_id)
  }
  
  # ---- Check Gene Consistency ----
  ref_genes <- gene_lists[[1]]
  for (i in seq_along(gene_lists)) {
    if (!identical(ref_genes, gene_lists[[i]])) {
      stop("Gene mismatch detected in file: ", names(gene_lists)[i])
    }
  }
  
  # ---- Build Count Matrix ----
  count_matrix <- do.call(cbind, all_counts)
  colnames(count_matrix) <- sample_ids
  count_matrix <- data.frame(SYMBOL = ref_genes, count_matrix, stringsAsFactors = FALSE)
  
  # ---- Filter Rows and Columns with All Zeros ----
  count_matrix <- count_matrix[rowSums(count_matrix[,-1]) > 0, , drop = FALSE]
  count_matrix <- count_matrix[, c(TRUE, colSums(count_matrix[,-1]) > 0), drop = FALSE]
  
  # ---- Export ----
  filename <- file.path(proj.params$proj.dir, paste0(proj.params$proj, "_Raw_counts.xlsx"))
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb = wb, sheetName = "Raw_counts")
  openxlsx::writeData(wb = wb, sheet = "Raw_counts", x = count_matrix)
  openxlsx::saveWorkbook(wb = wb, file = filename, overwrite = TRUE)
  message("✅ Saved counts to: ", file.path(proj.params$proj.dir, proj.params$proj))
  
  return(count_matrix)
}

prepare_deseq2_input <- function(meta_data, read_data, project_params) {
  
  set.seed(1234)
  
  # ---- Input checks ----
  if (!"Sample_ID" %in% colnames(meta_data)) {
    stop("`meta_data` must contain a 'Sample_ID' column.")
  }
  
  if (!"SYMBOL" %in% colnames(read_data)) {
    stop("`read_data` must contain a 'SYMBOL' column.")
  }
  if (anyDuplicated(read_data$SYMBOL)) {
    stop("The 'SYMBOL' column in `read_data` must not contain duplicate values.")
  }
  
  if (is.null(project_params$deseq2$design)) {
    stop("project_params$deseq2$design is missing. Please define the DESeq2 design formula.")
  }
  
  # ---- Initial summary ----
  message("Samples before filtering: ", nrow(meta_data))
  message("Genes before filtering: ", nrow(read_data))
  
  # ---- Clean and align meta_data ----
  meta_data <- meta_data %>%
    dplyr::filter(!is.na(Sample_ID)) %>%
    dplyr::mutate(Sample_ID = make.names(Sample_ID, unique = TRUE)) %>%
    dplyr::filter(Sample_ID %in% make.names(colnames(read_data)))
  rownames(meta_data) <- make.names(meta_data$Sample_ID)
  
  # Add Batch column if missing
  if (!"Batch" %in% colnames(meta_data)) {
    meta_data$Batch <- 1
    warning("No 'Batch' column found. Assigning all samples to Batch 1.")
  }
  
  # ---- Clean and align read_data ----
  colnames(read_data) <- make.names(colnames(read_data))
  valid_samples <- intersect(colnames(read_data), rownames(meta_data))
  
  read_data <- read_data %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "SYMBOL") %>%
    dplyr::select(all_of(valid_samples)) %>%
    replace(is.na(.), 0)
  
  read_data <- read_data[rowSums(read_data) != 0, ]
  
  # ---- Remove zero-count samples ----
  zero_samples <- which(colSums(read_data) == 0)
  if (length(zero_samples) > 0) {
    read_data <- read_data[, -zero_samples, drop = FALSE]
    meta_data <- meta_data[-zero_samples, , drop = FALSE]
  }
  
  # ---- Remove samples with NA in design variables ----
  design_vars <- stringr::str_split(string = project_params$deseq2$design, pattern = "[+*:]")[[1]]
  na_samples <- unique(unlist(lapply(design_vars, function(var) {
    if (var %in% colnames(meta_data)) {
      which(is.na(meta_data[[var]]))
    } else {
      warning(sprintf("Variable '%s' not found in meta_data.", var))
      integer(0)
    }
  })))
  
  if (length(na_samples) > 0) {
    read_data <- read_data[, -na_samples, drop = FALSE]
    meta_data <- meta_data[-na_samples, , drop = FALSE]
  }
  
  # ---- Remove sizeFactor column if present ----
  meta_data <- meta_data[, colnames(meta_data) != "sizeFactor", drop = FALSE]
  
  # ---- Convert meta_data columns to factors ----
  message("Structure of meta_data before conversion:")
  str(meta_data)
  meta_data[] <- lapply(meta_data, as.factor)
  message("Structure of meta_data after conversion:")
  str(meta_data)
  
  # ---- Reorder read_data to match meta_data ----
  read_data <- read_data[, rownames(meta_data), drop = FALSE]
  
  # ---- Sanity checks ----
  if (!is.data.frame(read_data) || !is.data.frame(meta_data)) {
    stop("`read_data` and `meta_data` must be data.frames")
  }
  if (!all(colnames(read_data) %in% rownames(meta_data))) {
    stop("Some samples in `read_data` are missing in `meta_data`")
  }
  if (!all(colnames(read_data) == rownames(meta_data))) {
    stop("Sample order mismatch between `read_data` and `meta_data`")
  }
  
  # ---- Return cleaned data ----
  return(invisible(list(meta_data = meta_data, read_data = read_data)))
}

plot_pca <- function(meta_data, read_data, output_path, top_n_genes = 500, file_name = "PCA_Plots.pdf") {
  
  set.seed(1234)
  
  # ---- Input Checks ----
  if (!is.data.frame(meta_data)) {
    stop("`meta_data` must be a data.frame")
  }
  if (!is.data.frame(read_data)) {
    stop("`read_data` must be a data.frame")
  }
  if (!"Sample_ID" %in% colnames(meta_data)) {
    stop("`meta_data` must contain a 'Sample_ID' column.")
  }
  
  # Warn if output path does not exist
  if (!dir.exists(output_path)) {
    warning("Output path does not exist. Attempting to create: ", output_path)
    dir.create(output_path, recursive = TRUE)
  }
  
  # ---- DESeq2 Object Preparation ----
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = read_data,
    colData = meta_data,
    design = ~1
  )
  dds <- DESeq2::estimateSizeFactors(dds)
  
  # ---- VST Transformation ----
  vsd <- DESeq2::vst(dds, blind = TRUE)
  
  # ---- PCA Data Preparation ----
  vst_mat <- SummarizedExperiment::assay(vsd) %>%
    as.data.frame() %>%
    dplyr::mutate(row_variance = matrixStats::rowVars(as.matrix(.))) %>%
    dplyr::slice_max(order_by = row_variance, n = top_n_genes) %>%
    dplyr::select(-row_variance)
  
  # ---- PCA Calculation ----
  pca <- stats::prcomp(x = t(vst_mat), center = TRUE, scale. = FALSE)
  
  # ---- Merge PCA Output with Metadata ----
  pca_df <- pca$x %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Sample_ID")
  
  df <- dplyr::inner_join(meta_data, pca_df, by=c("Sample_ID"="Sample_ID"))
  
  # ---- Variance explained ----
  percentVar <- round(100 * summary(pca)$importance[2, 1:2])
  comp_variables <- setdiff(colnames(meta_data), "Sample_ID")
  
  # ---- Save all plots to a single PDF ----
  output_file <- file.path(output_path, file_name)
  message("Saving PCA plots to: ", output_file)
  
  grDevices::pdf(file = output_file, width = 8.5, height = 11)
  
  for (var in comp_variables) {
    
    n_unique <- length(unique(meta_data[[var]]))
    if (n_unique < 2 || n_unique == nrow(meta_data)) {
      warning("Variable ", var, " is either constant or has one unique value 
      per sample; skipping PCA plot.")
      next
    }
    
    # Define color palette
    pca_palette <- custom_palette[1:length(unique(meta_data[[var]]))]
    names(pca_palette) <- as.character(unique(meta_data[[var]]))
    
    # Create PCA plot
    p <- ggplot2::ggplot(data = df, mapping = aes(x = PC1, y = PC2, color = get(var))) +
      ggplot2::geom_point(size = 3, shape = 16) +
      ggrepel::geom_text_repel(ggplot2::aes(label = Sample_ID), show.legend = FALSE) +
      ggplot2::theme_light() +
      ggplot2::labs(
        color = var,
        x = paste0("PC1: ", percentVar[1], "% variance"),
        y = paste0("PC2: ", percentVar[2], "% variance")
      ) +
      custom_theme +
      ggplot2::scale_color_manual(values = pca_palette)
    
    print(p)  # print to PDF (one page per plot)
  }
  
  dev.off()
}

plot_ma <- function(dds, output_path, file_name = "MA_Plot.pdf") {
  
  set.seed(1234)
  
  # ---- Input Checks ----
  if (!inherits(dds, "DESeqDataSet")) {
    stop("`dds` must be a DESeqDataSet object.")
  }
  
  if (!dir.exists(output_path)) {
    warning("Output path does not exist. Creating: ", output_path)
    dir.create(output_path, recursive = TRUE)
  }
  
  # ---- MA Plot ----
  output_file <- file.path(output_path, file_name)
  message("Saving MA plot to: ", output_file)
  
  grDevices::pdf(file = output_file, width = 8.5, height = 11)
  
  DESeq2::plotMA(
    object = dds,
    alpha  = 0.1, # FDR threshold (blue dots for significant genes)
    main   = "MA Plot",
    xlab   = "Mean of Normalized Counts",
    MLE    = FALSE
  )
  
  grDevices::dev.off()
  message("MA plot completed successfully.")
}

plot_dispersion <- function(dds, output_path, file_name = "Dispersion_Plot.pdf") {
  
  set.seed(1234)
  
  # ---- Input Checks ----
  if (!inherits(dds, "DESeqDataSet")) {
    stop("`dds` must be a DESeqDataSet object.")
  }
  
  if (!dir.exists(output_path)) {
    warning("Output path does not exist. Creating: ", output_path)
    dir.create(output_path, recursive = TRUE)
  }
  
  # ---- Dispersion Plot ----
  output_file <- file.path(output_path, file_name)
  message("Saving dispersion plot to: ", output_file)
  
  grDevices::pdf(file = output_file, width = 8.5, height = 11)
  
  DESeq2::plotDispEsts(
    object  = dds,
    genecol = "black",
    fitcol  = "red",
    finalcol= "dodgerblue",
    legend  = TRUE,
    xlab    = "Mean of Normalized Counts",
    ylab    = "Dispersion",
    log     = "xy",
    cex     = 0.45
  )
  
  grDevices::dev.off()
  
  message("Dispersion plot completed successfully.")
  # Expected results: Higher the mean, lower the dispersion
}

plot_volcano <- function(DEGs_df, proj.params, contrast = "Target-Reference", 
                         output_path, top_n = 5, file_prefix = "Volcano_Plot") {
  
  set.seed(1234)
  
  # ---- Input Checks ----
  required_cols <- c("log2FoldChange", "padj", "SYMBOL")
  if (!all(required_cols %in% colnames(DEGs_df))) {
    stop("DEGs_df must contain columns: log2FoldChange, padj, SYMBOL")
  }
  
  if (!dir.exists(output_path)) {
    warning("Output path does not exist. Creating: ", output_path)
    dir.create(output_path, recursive = TRUE)
  }
  
  # ---- Check required proj.params attributes ----
  required_attrs <- c("lfc.cutoff", "padj.cutoff")
  if (is.null(proj.params$volcano)) {
    stop("⚠️ proj.params must contain a 'volcano' list")
  }
  
  for (attr in required_attrs) {
    if (is.null(proj.params$volcano[[attr]])) {
      stop("⚠️ Missing required proj.params$volcano attribute: ", attr)
    } else {
      val <- proj.params$volcano[[attr]]
      message("✔ Found proj.params$volcano$", attr, " = ", paste(val, collapse = ", "))
    }
  }
  
  
  # ---- Volcano Parameters ----
  padj_cutoff <- proj.params$volcano$padj.cutoff
  lfc_cutoff  <- proj.params$volcano$lfc.cutoff
  label_genes <- if(!is.null(proj.params$volcano$label.genes)) proj.params$volcano$label.genes else NULL
  
  target <- stringr::str_split(string = contrast, pattern = "-")[[1]][1]
  reference <- stringr::str_split(string = contrast, pattern = "-")[[1]][2]
  
  # ---- Format DEGs ----
  DEGs_df <- DEGs_df %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::mutate(
      Direction = dplyr::case_when(
        padj < padj_cutoff & log2FoldChange > lfc_cutoff  ~ paste0("Up in ", target),
        padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ paste0("Up in ", reference),
        TRUE ~ "Not Significant"
      ),
      padj = dplyr::case_when(padj == 0 ~ min(padj[padj > 0], na.rm = TRUE), 
                              TRUE ~ padj),
      Significance = dplyr::case_when(
        abs(log2FoldChange) >= lfc_cutoff & padj <= 0.001 ~ "FDR < 0.001",
        abs(log2FoldChange) >= lfc_cutoff & padj <= 0.01  ~ "FDR < 0.01",
        abs(log2FoldChange) >= lfc_cutoff & padj <= 0.05  ~ "FDR < 0.05",
        TRUE ~ "Not Significant"
      ),
      Relevance = abs(log2FoldChange) * -log10(padj)
    )
  
  # ---- Color Palettes ----
  volcano_palette <- c(
    viridis::viridis(100)[50],  # Up in reference
    viridis::viridis(100)[1],   # Up in target
    viridis::viridis(100)[100]  # Not Significant
  )
  names(volcano_palette) <- c(
    paste0("Up in ", reference),
    paste0("Up in ", target),
    "Not Significant"
  )
  
  alpha_palette <- c("FDR < 0.001" = 1, "FDR < 0.01" = 0.8, 
                     "FDR < 0.05" = 0.6, "Not Significant" = 0.4)
  
  # ---- Axis limits and breaks ----
  x_vals <- DEGs_df$log2FoldChange
  y_vals <- -log10(DEGs_df$padj)
  
  # Keep only finite values
  x_vals <- x_vals[is.finite(x_vals)]
  y_vals <- y_vals[is.finite(y_vals)]
  
  # Round to nearest integer for nice axis limits
  x_min <- floor(min(x_vals, na.rm = TRUE))
  x_max <- ceiling(max(x_vals, na.rm = TRUE))
  y_min <- 0  # Start y-axis at 0 for volcano
  y_max <- ceiling(max(y_vals, na.rm = TRUE))
  
  # Set x-axis breaks in reasonable bins (approx 5 units each)
  x_bin <- max(abs(floor(x_min / 5)), abs(ceiling(x_max / 5)))
  x_breaks <- seq(from = -max(x_max, abs(x_min)), to = max(x_max, abs(x_min)), by = x_bin)
  x_breaks <- x_breaks[!x_breaks <= x_min-x_bin]
  x_breaks <- x_breaks[!x_breaks >= x_max+x_bin] 
  
  # Set y-axis breaks (dynamic, based on magnitude)
  y_bin <- if (y_max > 100) 100 else if (y_max > 10) 10 else 1
  y_breaks <- seq(from = y_min, to = ceiling(y_max/y_bin)*y_bin, by = y_bin)
  
  # ---- Build Base Plot ----
  p <- ggplot2::ggplot(data = DEGs_df, 
                       mapping = aes(x = log2FoldChange, y = -log10(padj),
                                     color = Direction, alpha = Significance,
                                     size = Relevance)) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.05, height = 0.05)) +
    ggplot2::theme_classic() +
    custom_theme +
    ggplot2::labs(x = expression("log"[2]*"FC"),
                  y = expression("-log"[10]*"padj"),
                  fill = "Direction",
                  title = contrast) +
    ggplot2::geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), 
                        linetype = "dotted", color = "black", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = -log10(padj_cutoff), 
                        linetype = "dotted", color = "black", linewidth = 0.5) +
    ggplot2::scale_color_manual(values = volcano_palette) +
    ggplot2::scale_alpha_manual(values = alpha_palette) +
    ggplot2::scale_size_continuous(range = c(0, 3)) +
    ggplot2::coord_cartesian(xlim = c(min(x_breaks), max(x_breaks)),
                             ylim = c(min(y_breaks), max(y_breaks))) +
    ggplot2::scale_x_continuous(breaks = x_breaks, 
                                labels = function(x) { base::ifelse(x %% 1 == 0, as.integer(x), format(x, digits = 2)) }) +
    ggplot2::scale_y_continuous(breaks = y_breaks) +
    ggplot2::guides(size = "none",
                    shape = guide_legend(override.aes = list(size = 3)),
                    fill = guide_colourbar(theme = theme(legend.key.width = unit(0.75, "lines"),
                                                         legend.key.height = unit(10, "lines"),
                                                         legend.ticks = element_blank(),
                                                         legend.frame = element_rect(colour = "Black", linewidth = 1)))) 
  
  # ---- Identify Top 5 Up and Down regulated Genes ----
  predicted_gene_pattern <- "^(Gm[0-9]+|ENSMUSG[0-9]+|ENSG[0-9]+|LOC[0-9]+|C[0-9]+orf[0-9]+|RP[0-9]+-)|Rik$"
  
  # Filter out predicted/placeholder genes and non-significant genes
  filtered_DEGs <- DEGs_df %>%
    dplyr::filter(!stringr::str_detect(string = SYMBOL, pattern = predicted_gene_pattern)) %>%
    dplyr::filter(padj < padj_cutoff)
  
  # Top 5 up-regulated (log2FC > lfc_cutoff)
  top_up <- filtered_DEGs %>%
    dplyr::filter(log2FoldChange > lfc_cutoff) %>%
    #dplyr::slice_min(padj, n = 5, with_ties = FALSE) %>%
    dplyr::slice_max(Relevance, n = 5, with_ties = FALSE) %>%
    dplyr::pull(SYMBOL)
  
  # Top 5 down-regulated (log2FC < -lfc_cutoff)
  top_down <- filtered_DEGs %>%
    dplyr::filter(log2FoldChange < -lfc_cutoff) %>%
    #dplyr::slice_min(padj, n = 5, with_ties = FALSE) %>%
    dplyr::slice_max(Relevance, n = 5, with_ties = FALSE) %>%
    dplyr::pull(SYMBOL)
  
  # Combine
  top_genes <- unique(c(top_up, top_down))
  
  # ---- Label Top Genes ----
  if (!is.null(label_genes)) {
    genes_to_label <- intersect(label_genes, DEGs_df$SYMBOL)
  } else {
    genes_to_label <- top_genes
  }
  
  q <- p + ggrepel::geom_text_repel(
    data = DEGs_df %>% dplyr::filter(SYMBOL %in% genes_to_label),
    aes(label = SYMBOL),
    direction = "both",
    box.padding = 0.8,                # ↓ smaller padding around label
    point.padding = 0.1,              # minimal space between point and line start
    max.overlaps = nrow(DEGs_df),
    show.legend = FALSE,
    min.segment.length = 0,           # Only draw segments longer than this
    segment.curvature = -0.5,         # Negative = curve upward, positive = downward
    segment.ncp = 50,                 # More control points = smoother curves
    segment.angle = 20,               # Affects entry/exit angles
    segment.size = 0.5,               # Optional: line thickness
    size = 4,                         # text size in mm (1 mm = 2.83 points)
    position = ggbeeswarm::position_quasirandom(width = 0.1, varwidth = TRUE)
  )
  
  # ---- Save Plots ----
  ggplot2::ggsave(file.path(output_path, paste0(file_prefix, "_", contrast, ".pdf")),
                  plot = p, width = 7, height = 7, device = "pdf")
  
  ggplot2::ggsave(file.path(output_path, paste0(file_prefix, "_top_", contrast, ".pdf")),
                  plot = q, width = 7, height = 7, device = "pdf")
  
  message("Volcano plots saved successfully.")
  
  return(invisible(p))
}

# metadata has column Sample_ID and columns defined in proj.params$heatmap.col.ann
# metadata_row has column SYMBOL and columns defined in proj.params$heatmap.row.ann
# metadata_row <- NULL if no row annotations needed
# disp_genes either empty vector c() or vector of genes
plot_heatmap <- function(norm_counts, proj.params, metadata_col = NULL, metadata_row = NULL, disp_genes = c()) {
  
  set.seed(1234)
  
  # ---- Check required proj.params attributes ----
  if (is.null(proj.params$heatmap)) {
    stop("proj.params must contain a 'heatmap' list")
  }
  required_attrs <- c("force.log","col.ann","ann.palette","palette",
                      "col.cluster","row.cluster",
                      "title","border.color","show.expr.legend")
  
  for (attr in required_attrs){
    if(is.null(proj.params$heatmap[[attr]])){
      stop("Missing proj.params$heatmap$", attr)
    } else {
      message("✔ Found proj.params$heatmap$", attr, " = ", paste(proj.params$heatmap[[attr]], collapse=", "))
    }
  }
  
  # ---- Input check ----
  if (nrow(norm_counts) < 2){
    message("Input data frame has less than 2 genes. Skipping plotting.")
    return(NULL)
  }
  
  # ---- Prepare matrix ----
  mat <- norm_counts %>%
    as.data.frame() %>%
    tibble::rownames_to_column("SYMBOL") %>%
    base::replace(is.na(.), 0) %>%
    # This section retains gene copy with highest expression in case of duplicates
    # dplyr::mutate(n = rowSums(.[, -1])) %>%
    # dplyr::group_by(SYMBOL) %>%
    # dplyr::slice_max(n) %>%
    # dplyr::ungroup() %>%
    # dplyr::filter(n != 0) %>%
    # dplyr::select(-n) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("SYMBOL") %>%
    dplyr::mutate(across(.cols = everything(), .fns = as.numeric))
  
  rownames(mat) <- make.names(rownames(mat), unique = TRUE)
  colnames(mat) <- make.names(colnames(mat), unique = TRUE)
  
  # Log transform if necessary
  quantiles <- stats::quantile(x = as.vector(as.matrix(mat)), probs = c(0, 0.01, 0.99, 1), na.rm = TRUE)
  huge_range <- (quantiles[4] - quantiles[1]) > 100   # Range of values greater than 100
  only_pos <- quantiles[1] >= 0                       # Min value greater than 0
  if ((huge_range & only_pos) | proj.params$heatmap$force.log){
    mat <- log2(1 + mat)
  }
  
  # Scale every feature across samples 
  mat_scaled <- t(scale(t(mat)))
  mat_scaled[is.na(mat_scaled)] <- 0
  
  # ---- Annotations ----
  # Define column annotations 
  if(is.null(metadata_col)) {
    col_annotation <- NULL
  } else{
    col_annotation <- metadata_col %>%
      dplyr::select(Sample_ID, all_of(proj.params$heatmap$col.ann)) %>%
      dplyr::mutate(Sample_ID = make.names(Sample_ID)) %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames("Sample_ID") %>%
      as.data.frame() %>%
      mutate(across(where(is.factor), as.character))
  }
  
  # Define row annotations
  if(is.null(metadata_row)) {
    row_annotation <- NULL
  } else{
    row_annotation <- metadata_row %>%
      dplyr::select(SYMBOL, all_of(proj.params$heatmap$row.ann)) %>%
      dplyr::mutate(SYMBOL = make.names(SYMBOL, unique = TRUE)) %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames("SYMBOL") %>%
      as.data.frame() %>%
      mutate(across(where(is.factor), as.character))
  }
  
  # ---- Annotation Palettes ----
  # Define Color Palette for Annotation 
  # This is an example of how ann_colors should be specified
  ann_colors <- list(CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
                     GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E"))
  
  ann_colors <- list()
  base_colors <- c("#E08214", "#762A83", "#C51B7D", "#7FBC41", "#35978F", "#BF812D", "#542788",
                   "#D6604D", "#4393C3", "#878787", "#E41A1C", "#F781BF", "#4DAF4A", "#FFFFBF",
                   "#377EB8", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#999999", "#66C2A5",
                   "#FC8D62", "#000000", "#9E0142", "#1A1A1A", "#74006F", "#FFC606", "#F6D2E0",
                   "#C8E7F5")
  base_colors <- custom_palette
  
  col_list <- base::lapply(X = as.list(col_annotation), FUN = function(x) { as.character(x) %>% unique})
  row_list <- base::lapply(X = as.list(row_annotation), FUN = function(x) { as.character(x) %>% unique})
  ann_list <- c(col_list, row_list)
  
  color_index <- 1
  for (i in seq_along(ann_list)) {  # Iterate through each annotation variable (Eg: CellType) 
    levels <- sort(ann_list[[i]])  # Get levels within each annotation variable (Eg: CT1, CT2)
    n_levels <- length(levels)    # Get number of levels within each annotation variable
    
    palette_colors <- if (proj.params$heatmap$ann.palette == "discrete" | n_levels == 1){
      base_colors[color_index:(color_index + n_levels - 1)]
    } else{
      alphas <- seq(1 / n_levels, 1, length.out = n_levels)
      base::sapply(X = alphas, 
                   FUN = function(x) { colorspace::adjust_transparency(col = base_colors[color_index], alpha = x) })
    }
    
    names(palette_colors) <- levels          # Name each color with levels
    ann_colors <- c(ann_colors, list(palette_colors)) # Append named color palette
    names(ann_colors)[i] <- names(ann_list)[i]     # Name the color palette with corresponding annotation variable name
    color_index <- color_index + n_levels       # Move to next color
  }
  
  # ---- Heatmap Palette ----
  # Define Color Palette for Heatmap 
  valid_palettes <- c("vrds", "rdbu")
  
  if (!proj.params$heatmap$palette %in% valid_palettes) {
    stop("Invalid heatmap palette. proj.params$heatmap$palette must be either 'vrds' or 'rdbu'.")
  }
  
  heatmap_palette <- switch(proj.params$heatmap$palette,
                            vrds = viridis::viridis(100),
                            rdbu = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100))
  
  # ---- Breaks ----
  # Define Color Breaks 
  n_breaks <- 100
  heatmap_palette <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(n_breaks)
  
  # Handle min and max thresholds with soft clamping
  mat_min <- min(mat_scaled, na.rm = TRUE)
  mat_max <- max(mat_scaled, na.rm = TRUE)
  mat_min <- dplyr::case_when(mat_min >= 0 ~ 0,
                              mat_min <= -3 ~ -3,
                              TRUE ~ mat_min)
  mat_max <- dplyr::case_when(mat_max <= 0 ~ 0,
                              mat_max >= 3 ~ 3,
                              TRUE ~ mat_max)
  
  if (mat_max == 0){
    breaks <- seq(from = floor(mat_min), to = 0, length.out = n_breaks)
  } else if (mat_min == 0){
    breaks <- seq(from = 0, to = ceiling(mat_max), length.out = n_breaks)
  } else{
    breaks <- c(seq(from = floor(mat_min),   to = 0,        length.out = n_breaks / 2),
                seq(from = mat_max / n_breaks, to = ceiling(mat_max), length.out = n_breaks / 2))
  }
  
  # ---- Gaps ----
  # Define gaps in heatmap 
  gaps_col <- if (!gtools::invalid(proj.params$heatmap$col.gaps) & proj.params$heatmap$col.cluster %in% colnames(col_annotation)) {
    if (all(proj.params$heatmap$col.gaps %in% colnames(col_annotation))) {
      col_annotation %>%
        dplyr::count(get(proj.params$heatmap$col.gaps)) %>%
        dplyr::mutate(n = cumsum(n)) %>%
        dplyr::pull(n) %>%
        .[. < ncol(mat_scaled)]
    } else{
      message("Gaps not introduced as gap identifiers missing in column annotation")
    }
  }
  
  gaps_row <- if (!gtools::invalid(proj.params$heatmap$row.gaps) & proj.params$heatmap$row.cluster %in% colnames(row_annotation)) {
    if (all(proj.params$heatmap$row.gaps %in% colnames(row_annotation))) {
      row_annotation %>%
        dplyr::count(get(proj.params$heatmap$row.gaps)) %>%
        dplyr::mutate(n = cumsum(n)) %>%
        dplyr::pull(n) %>%
        .[. < nrow(mat_scaled)]
    } else{
      message("Gaps not introduced as gap identifiers missing in row annotation")
    }
  }
  
  # ---- Clustering ----
  # Determine Ordering 
  if (proj.params$heatmap$col.cluster == "all"){
    colclust <- hclust(dist(t(mat_scaled)))
    col_order <- colnames(mat_scaled)[colclust$order]
  }
  if (proj.params$heatmap$col.cluster == "alphabetical"){
    col_order <- sort(colnames(mat_scaled))
  }
  if (proj.params$heatmap$col.cluster %in% colnames(col_annotation)){
    
    # NOTE: While calculating gaps_col, we use count(). It sorts alphabetically.
    # So, WE MUST sort col_elements to match gaps_col
    col_order <- c()
    col_elements <- col_annotation %>% 
      dplyr::pull(proj.params$heatmap$col.cluster) %>%
      unique() %>% sort()
    
    for (g in col_elements){
      
      samples <- rownames(col_annotation)[col_annotation %>% dplyr::pull(proj.params$heatmap$col.cluster) == g]
      if (length(samples) == 0) next
      temp_mat <- mat_scaled[, samples]
      
      if (length(samples) > 1){
        colclust <- hclust(dist(t(temp_mat)))
        col_order <- c(col_order, colnames(temp_mat)[colclust$order])
      } else if(length(samples) == 1){
        col_order <- c(col_order, samples)
      }
    }
  }
  
  if (proj.params$heatmap$row.cluster == "all"){
    rowclust <- hclust(dist(mat_scaled))
    row_order <- rownames(mat_scaled)[rowclust$order]
  }
  if (proj.params$heatmap$row.cluster == "alphabetical"){
    row_order <- sort(rownames(mat_scaled))
  }
  if (proj.params$heatmap$row.cluster %in% colnames(row_annotation)){
    
    # NOTE: While calculating gaps_row, we use count(). It sorts alphabetically.
    # So, WE MUST sort row_elements to match gaps_row
    row_order <- c()
    row_elements <- row_annotation %>% 
      dplyr::pull(proj.params$heatmap$row.cluster) %>%
      unique() %>% sort()
    
    for (g in row_elements){
      
      genes <- rownames(row_annotation)[row_annotation %>% dplyr::pull(proj.params$heatmap$row.cluster) == g]
      if (length(genes) == 0) next
      temp_mat <- mat_scaled[rownames(mat_scaled) %in% genes,]
      
      if (length(genes) > 1){
        rowclust <- hclust(dist(temp_mat))
        row_order <- c(row_order, rownames(temp_mat)[rowclust$order])
      } else if(length(genes) == 1){
        row_order <- c(row_order, genes)
      }
    }
  }
  
  # ---- Prepare Matrix for Plotting ---- 
  reordered <- mat_scaled[row_order, col_order]
  
  # ---- Set font sizes ----
  fontsize <- 10
  fontsize.row <- fontsize*1
  fontsize.col <- fontsize*1
  fontsize.number <- fontsize*0.8
  angle.col <- 45       # column label angle
  
  # ---- Set cell width, height (in points) dynamically ----
  cell.width <- dplyr::if_else(ncol(reordered) <= 60, fontsize+2, NA)
  cell.height <- dplyr::if_else(nrow(reordered) <= 48, fontsize+2, NA)
  
  # ---- Truncate long row and column labels ----
  main.title <- stringr::str_wrap(string = proj.params$heatmap$title, width = 20)
  labels.col <- stringr::str_trunc(string = colnames(reordered), width = 15, side = "right", ellipsis = "…")
  labels.row <- stringr::str_trunc(string = rownames(reordered), width = 15, side = "right", ellipsis = "…")
  if(length(disp_genes) > 0){
    labels.row <- dplyr::if_else(rownames(reordered) %in% make.names(disp_genes), rownames(reordered), " ")
  }
  
  # ---- Heatmap plotting ----
  ph <- pheatmap::pheatmap(mat               = reordered,
                           color             = heatmap_palette,
                           breaks            = breaks,
                           annotation_row    = row_annotation,
                           annotation_col    = col_annotation,
                           annotation_colors = ann_colors,
                           gaps_row          = gaps_row,
                           gaps_col          = gaps_col,
                           
                           cellwidth         = cell.width,     
                           cellheight        = cell.height,  
                           show_rownames     = !is.na(cell.height),
                           show_colnames     = !is.na(cell.width),
                           labels_row        = labels.row,
                           labels_col        = labels.col,
                           angle_col         = angle.col,        # column label angle
                           fontsize          = fontsize,         # points; 72 points = 1 inch
                           fontsize_row      = fontsize.row,     # points
                           fontsize_col      = fontsize.col,     # points
                           fontsize_number   = fontsize.number,  # points
                           silent            = TRUE, 
                           
                           main              = main.title,
                           border_color      = proj.params$heatmap$border.color,
                           legend            = proj.params$heatmap$show.expr.legend,
                           
                           scale                    = "none",
                           cluster_rows             = FALSE,
                           cluster_cols             = FALSE,
                           clustering_distance_rows = "correlation", #"euclidean",
                           clustering_distance_cols = "correlation", #"euclidean",
                           clustering_method        = "average",     #complete",
                           annotation_legend        = TRUE,
                           annotation_names_row     = FALSE,
                           annotation_names_col     = FALSE,
                           width                    = NA,               # inches
                           height                   = NA,               # inches
                           filename                 = NA)
  
  # ---- Return matrix for Excel ----
  ph_mat <- if(ncol(reordered) > nrow(reordered)) t(reordered) else reordered
  
  return(invisible(list(ph = ph, mat = ph_mat)))
}

run_deseq2 <- function(meta_data, read_data, proj.params, n = 1) {
  
  set.seed(1234)
  
  # ---- Input Checks ----
  if (!is.list(proj.params) || !"deseq2" %in% names(proj.params)) {
    stop("⚠️ proj.params must contain a 'deseq2' list")
  }
  
  required_attrs <- c("design", "lfc.cutoff", "padj.cutoff", "contrasts")
  for (attr in required_attrs) {
    if (is.null(proj.params$deseq2[[attr]])) {
      stop("⚠️ Missing required proj.params$deseq2 attribute: ", attr)
    } else {
      val <- proj.params$deseq2[[attr]]
      message("✔ Found proj.params$deseq2$", attr, " = ", paste(val, collapse = ", "))
    }
  }
  
  contrast <- proj.params$deseq2$contrasts[n]
  contrast.dir <- proj.params$contrast.dir[n]
  
  if (!dir.exists(contrast.dir)) {
    warning("Contrast directory does not exist. Creating: ", contrast.dir)
    dir.create(contrast.dir, recursive = TRUE)
  }
  
  # ---- DESeq2 Object Preparation ----
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = read_data,
                                        colData   = meta_data,
                                        design    = ~1)
  design(dds) <- as.formula(paste0("~", proj.params$deseq2$design))
  
  # ---- Size Factor Estimation ----
  if (all(rowSums(read_data == 0) > 0)) {
    message("Cannot compute log geometric means as every gene contains at least
            one zero. Using 'poscounts' for size factor estimation.")
    dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
  }
  
  # ---- Pre-filter Lowly Expressed Genes ----
  # NOTE: This improves sizefactor estimation
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  
  # ---- DESeq2 Fit: Parametric and Local ----
  dds_para <- DESeq2::DESeq(object = dds, test = "Wald", fitType = "parametric",
                            betaPrior = FALSE, minReplicatesForReplace = 7)
  dds_local <- DESeq2::DESeq(object = dds, test = "Wald", fitType = "local",
                             betaPrior = FALSE, minReplicatesForReplace = 7)
  
  # ---- Select Best Fit Based on Residuals ----
  residual_para  <- mcols(dds_para)$dispGeneEst - mcols(dds_para)$dispFit
  residual_local <- mcols(dds_local)$dispGeneEst - mcols(dds_local)$dispFit
  dds <- if (median(residual_para^2, na.rm = TRUE) <= median(residual_local^2, na.rm = TRUE)) {
    dds_para
  } else {
    dds_local
  }
  
  # ---- Prepare Contrast Vector ----
  mod_mat <- model.matrix(design(dds), colData(dds))
  design_factors <- stringr::str_split(string = proj.params$deseq2$design, 
                                       pattern = "[+*:]")[[1]] %>% unique()
  
  # Create a replicate of meta_data with all possible groups that could be compared based on the design
  df <- colData(dds) %>%
    as.data.frame() %>%
    tidyr::unite(col = "Groups", all_of(design_factors), sep = ".")
  
  # Define all possible groups that could be compared based on the design
  groups <- unique(df$Groups)
  
  # Get all possible coefficient vectors
  group_coef_list <- lapply(groups, function(i) colMeans(as.matrix(mod_mat[df$Groups == i, , drop = FALSE])))
  names(group_coef_list) <- groups
  
  target <- stringr::str_split(string = contrast, pattern = "-")[[1]][1]
  ref    <- stringr::str_split(string = contrast, pattern = "-")[[1]][2]
  
  # Below line works for simple contrasts like (A - B). Fails if contrast is (A-B)-(C-D)
  # contrast_vec <- group_coef_list[[target]] - group_coef_list[[ref]] 
  
  replace_symbols <- function(node) {
    if (is.symbol(node)) {
      nm <- as.character(node)
      if (nm %in% names(group_coef_list)) {
        return(group_coef_list[[nm]])
      } else {
        # it's an operator like "-" or "+" → return unchanged
        return(node)
      }
    } else if (is.call(node)) {
      return(as.call(lapply(node, replace_symbols)))
    } else {
      return(node)
    }
  }
  
  parsed <- base::parse(text = contrast)[[1]]
  expr_sub <- replace_symbols(parsed)
  contrast_vec <- base::eval(expr_sub)
  
  # # [ALTERNATIVE METHOD] Get all possible coefficient vectors
  # for (i in groups) {
  #   val <- colMeans(as.matrix(mod_mat[df$Groups == i, , drop =FALSE]))
  #   assign(x = i, value = val)
  # }
  # # Define the contrast vector
  # contrast_vec <- base::eval(expr = base::parse(text = contrast))
  
  # ---- DESeq2 Results with LFC Threshold & Shrinkage ----
  res <- DESeq2::results(object = dds, 
                         contrast = contrast_vec,
                         lfcThreshold = proj.params$deseq2$lfc.cutoff,
                         altHypothesis = "greaterAbs",
                         cooksCutoff = TRUE,
                         independentFiltering = TRUE,
                         alpha = proj.params$deseq2$padj.cutoff,
                         pAdjustMethod = "BH")
  
  # Safe LFC shrinkage
  set.seed(1234)
  res <- DESeq2::lfcShrink(dds = dds, res = res, type = "ashr")
  summary(res)
  
  # ---- Differential Expressed Genes (DEGs) ----
  DEGs_df <- res %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID") %>%
    add_annotation() %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarize(across(.cols = where(is.numeric), .fns = mean, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(padj = ifelse(padj == 0, min(padj[padj > 0], na.rm = TRUE), padj))
  
  save_xlsx(DEGs_df, file.path(contrast.dir, "DEGs.xlsx"), "DEGs", row_names = FALSE)
  
  # ---- VST Counts (Non-blind) ----
  vsd <- DESeq2::vst(dds, blind = FALSE)
  vst_counts <- SummarizedExperiment::assay(vsd) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID") %>%
    add_annotation() %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarize(across(.cols = where(is.numeric), .fns = mean, na.rm = TRUE), .groups = "drop") %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("SYMBOL") %>%
    dplyr::select(-dplyr::starts_with("ENSEMBL"), -dplyr::starts_with("ENTREZ")) %>%
    as.matrix()
  
  save_xlsx(vst_counts, file.path(contrast.dir, "VST_counts.xlsx"), "VST_Nonblind", row_names = TRUE)
  
  # ---- Return Results ----
  invisible(list(degs = DEGs_df, vst = vst_counts, dds = dds))
}

# DEGs_df with column SYMBOL, padj, log2FoldChange
# k <- # overlapping genes between input and pathway
# n <- # overlapping genes between input and collection
# K <- # genes in pathway
# N <- # genes in collection
pathway_analysis <- function(DEGs_df, proj.params, output_path) {
  
  # Check required proj.params attributes
  required_attrs <- c("gmt.dir", "species")
  for (attr in required_attrs) {
    if (is.null(proj.params[[attr]])) {
      stop("⚠️ Missing required proj.params attribute: ", attr)
    } else {
      val <- proj.params[[attr]]
      message("✔ Found proj.params attribute: ", attr, " = ", val)
    }
  }
  
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  
  set.seed(1234)
  
  gmt_files <- list.files(file.path(proj.params$gmt.dir, proj.params$species), full.names = TRUE)
  
  # Initialize result dataframes 
  fgsea_df <- data.frame()
  gsea_df <- data.frame()
  ora_df_up <- data.frame()
  ora_df_down <- data.frame()
  concise_fgsea_df <- data.frame()
  
  # Define input genes for GSEA (Ranked list of all genes) 
  # IMPORTANT: Rank genes from high LFC to low lFC, so +NES ~ up-regulated
  ranked_df <- DEGs_df %>%
    dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
    dplyr::filter(!is.na(padj), !is.na(SYMBOL)) %>%
    dplyr::mutate(log2FoldChange = as.numeric(log2FoldChange),
                  padj = as.numeric(padj)) %>%
    dplyr::arrange(desc(log2FoldChange))  
  
  ranked_list <- ranked_df$log2FoldChange
  names(ranked_list) <- ranked_df$SYMBOL
  
  # Define input and universe genes for ORA (significant genes only) 
  sig_genes_up <- DEGs_df %>%
    dplyr::filter(padj <= 0.05, log2FoldChange > 0) %>%
    dplyr::pull(SYMBOL)
  
  sig_genes_down <- DEGs_df %>%
    dplyr::filter(padj <= 0.05, log2FoldChange < 0) %>%
    dplyr::pull(SYMBOL)
  
  universe_genes <- unique(DEGs_df$SYMBOL)
  
  # Iterate over GMT files and perform enrichment analyses
  for (gmt_file in gmt_files) {
    
    # Extract gene set name
    # gmt_name <- gsub(pattern = "^.*/|.v[0-9].*$", replacement = "", x = gmt_file)
    gmt_name <- gsub(pattern = "^.*/|", replacement = "", x = gmt_file)
    
    # Format gene sets for fgsea and keep only genes present in ranked_list
    gmt <- fgsea::gmtPathways(gmt_file)
    gmt <- lapply(X = gmt, FUN = function(x){x[x %in% names(ranked_list)]})
    
    # Format gene sets for clusterProfiler and keep only genes present in ranked_list
    pathway_gene_df <- data.frame(pathways = base::rep(x = names(gmt), times = base::unname(obj = lengths(gmt))),
                                  genes = unlist(gmt, use.names = FALSE))
    
    # Run fgseaMultilevel (GSEA)
    fgsea_res <- fgsea::fgseaMultilevel(pathways = gmt,
                                        stats = ranked_list,
                                        scoreType = dplyr::case_when(min(ranked_list) > 0 ~ "pos",
                                                                     max(ranked_list) < 0 ~ "neg",
                                                                     TRUE ~ "std"),
                                        sampleSize = 101,
                                        minSize = 1,
                                        maxSize = 500, # recommended 500 genes max
                                        eps = 1e-50,
                                        nproc = 0,
                                        gseaParam = 1,
                                        BPPARAM = NULL,
                                        nPermSimple = 10000)
    
    # Identify overlapping pathways and collapse into major pathways
    concise_fgsea_res <- fgsea::collapsePathways(fgseaRes = fgsea_res,
                                                 pathways = gmt,
                                                 stats = ranked_list)
    concise_fgsea_res <- fgsea_res %>%
      dplyr::filter(pathway %in% concise_fgsea_res$mainPathways)
    
    # Run clusterProfiler GSEA
    gsea_res <- clusterProfiler::GSEA(geneList = ranked_list,
                                      exponent = 1,
                                      minGSSize = 10,
                                      maxGSSize = 500,
                                      eps = 1e-10,
                                      pvalueCutoff = 0.05,
                                      pAdjustMethod = "BH",
                                      TERM2GENE = pathway_gene_df,
                                      TERM2NAME = NA,
                                      verbose = TRUE,
                                      seed = FALSE,
                                      by = "fgsea")
    
    # Run clusterProfiler ORA (enricher)
    # NOTE: Avoid using clusterProfiler::enrichGO() as it doesnt use proper 
    # background in universe parameter and includes GO terms outside the 
    # intended gene set collection.
    ora_res_up <- clusterProfiler::enricher(gene = sig_genes_up,
                                            pvalueCutoff = 0.05,
                                            pAdjustMethod = "BH",
                                            universe = universe_genes,
                                            minGSSize = 10,
                                            maxGSSize = 500,
                                            qvalueCutoff = 0.2,
                                            TERM2GENE = pathway_gene_df,
                                            TERM2NAME = NA)
    
    ora_res_down <- clusterProfiler::enricher(gene = sig_genes_down,
                                              pvalueCutoff = 0.05,
                                              pAdjustMethod = "BH",
                                              universe = universe_genes,
                                              minGSSize = 10,
                                              maxGSSize = 500,
                                              qvalueCutoff = 0.2,
                                              TERM2GENE = pathway_gene_df,
                                              TERM2NAME = NA)
    
    # Bind significant results to respective dataframes
    fgsea_df <- dplyr::bind_rows(fgsea_df, fgsea_res)
    concise_fgsea_df <- dplyr::bind_rows(concise_fgsea_df, concise_fgsea_res)
    
    if (!is.null(gsea_res)) {
      gsea_df <- dplyr::bind_rows(gsea_df, gsea_res@result)
    }
    
    if (!is.null(ora_res_up)) {
      ora_df_up <- dplyr::bind_rows(ora_df_up, ora_res_up@result)
    }
    
    if (!is.null(ora_res_down)) {
      ora_df_down <- dplyr::bind_rows(ora_df_down, ora_res_down@result)
    }
  } 
  
  # Format the results
  # Rename columns consistently across different methods
  lookup <- c(pathway = "ID",
              geneID = "leadingEdge", geneID = "core_enrichment",
              K = "size", K = "setSize",
              padj = "p.adjust", 
              pval = "pvalue")
  
  ora_df <- dplyr::bind_rows(ora_df_up %>% dplyr::mutate(Direction = "Upregulated"), 
                             ora_df_down %>% dplyr::mutate(Direction = "Downregulated"))
  if (nrow(ora_df) > 0){
    ora_df <- ora_df %>%
      tidyr::separate(col = GeneRatio, into = c("k", "n")) %>%
      tidyr::separate(col = BgRatio, into = c("K", "N")) %>%
      dplyr::mutate_at(c("k", "n", "K", "N"), as.numeric) %>%
      dplyr::mutate(GeneRatio = k / n,
                    BackgroundRatio = K / N,
                    EnrichmentRatio = GeneRatio / BackgroundRatio,
                    combined_score = GeneRatio * -log10(p.adjust),
                    NES = NA_integer_) %>%
      dplyr::rename(any_of(lookup))
  }
  
  for (i in c("fgsea_df", "gsea_df")){
    df <- get(i) %>% 
      dplyr::mutate(Direction = dplyr::case_when(NES > 0 ~ "Upregulated",
                                                 NES < 0 ~ "Downregulated",
                                                 TRUE ~ "No change")) %>%
      dplyr::rename(any_of(lookup))
    assign(x = i, value = df)
  }
  
  for (i in c("fgsea_df", "gsea_df", "ora_df")){
    
    df <- get(i)
    if (nrow(df) > 0){
      df <- df %>%
        dplyr::filter(padj <= 0.05) %>%
        tibble::remove_rownames() %>%
        tidyr::separate(col = pathway, into = c("Collection", "Description"), sep = "_", extra = "merge") %>%
        dplyr::mutate(Description = base::gsub(pattern = "_", replacement = " ", x= Description),
                      geneID = base::sapply(X = geneID, FUN = paste, collapse = "/")) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(leading_edge_size = length(unlist(stringr::str_split(geneID, "/")))) %>%
        dplyr::ungroup() %>%
        as.data.frame() %>%
        dplyr::select(Collection, Description, leading_edge_size, K, pval, padj, NES, Direction, everything(), -geneID, geneID)
      
      max_len <-  max(df$leading_edge_size, na.rm = TRUE)
      sep_char <- ifelse(grepl(",", df$geneID[1], fixed = TRUE), ",", "/")
      
      # across() allows you to apply a function(s) to multiple columns at once. 
      # map_chr() iterates over a list/vector and applies a function to each 
      # element (in this context, each cell within the selected columns),
      # returning a character vector.
      if (is.finite(max_len) & max_len > 0){
        df <- df %>% 
          tidyr::separate(col = geneID, into = paste0("gene", 1:max_len), sep = sep_char, remove = TRUE, fill = "right") %>%
          dplyr::mutate(across(.cols = starts_with("gene", ignore.case = FALSE), 
                               .fns = function (col) { col %>% 
                                   purrr::map_chr(.f = function (x){gsub('c\\(', "", x) }) %>%
                                   purrr::map_chr(.f = function (x){gsub('\\)', "", x) }) %>%
                                   purrr::map_chr(.f = function (x){gsub('"', "", x) }) %>%
                                   trimws() })) %>%
          dplyr::select(Collection, Description, leading_edge_size, K, padj, NES, Direction, everything())
      }
      assign(x = i, value = df)
    }
  }
  
  # Summarize results from fgsea, gsea and ora 
  consensus_df <- dplyr::bind_rows(fgsea_df %>% dplyr::mutate(method = "FGSEA"), 
                                   gsea_df %>% dplyr::mutate(method = "GSEA"), 
                                   ora_df %>% dplyr::mutate(method = "ORA")) %>%
    dplyr::add_count(Collection, Description, Direction, name = "n_methods") %>%
    dplyr::filter(n_methods > 1) %>%
    dplyr::mutate(Consensus = Direction) %>%
    dplyr::arrange(Collection, Description, desc(NES)) %>%
    dplyr::select(n_methods, method, Consensus, Collection, Description, 
                  leading_edge_size, K, padj, NES, Direction, everything(), 
                  -starts_with("gene",ignore.case = FALSE),
                  starts_with("gene", ignore.case = FALSE)) 
  
  gsea_res_list <- list(consensus = consensus_df,
                        fgsea = fgsea_df, 
                        gsea = gsea_df, 
                        ora = ora_df)
  # Save results
  wb <- openxlsx::createWorkbook()
  for (i in seq_along(gsea_res_list)) {
    openxlsx::addWorksheet(wb, sheetName = names(gsea_res_list)[i])
    openxlsx::writeData(wb, sheet = names(gsea_res_list)[i], x = gsea_res_list[[i]], rowNames = FALSE)
  }
  openxlsx::saveWorkbook(wb, file.path(output_path, "Pathway_results.xlsx"), overwrite = TRUE)
  
  # Return results
  return(invisible(list(fgsea = fgsea_df, 
                        gsea = gsea_df, 
                        ora = ora_df,
                        consensus = consensus_df)))
}

### decoupleR Analysis [RECOMMENDED as it can run multiple algorithms & give consensus]
# input can be vst_counts or DEGs with t-statistic (-log10padj*log2FC)
tf_analysis <- function(input, species = "Homo sapiens", top = 500) {
  
  set.seed(1234)
  
  # --- Input checks ---
  if (!is.matrix(input)) stop("`input` must be a matrix.")
  if (!species %in% c("Homo sapiens", "Mus musculus")) {
    stop("`species` must be either 'Homo sapiens' or 'Mus musculus'.")
  }
  
  # --- Map species to decoupleR format ---
  organism <- dplyr::case_when(species == "Homo sapiens" ~ "human",
                               species == "Mus musculus" ~ "mouse",
                               TRUE ~ "rat")
  
  # Load network models 
  pathway_net <- decoupleR::get_progeny(organism = organism, top = top)
  tf_net <- decoupleR::get_collectri(organism = organism, split_complexes = FALSE)
  
  # Run decoupleR
  # stats <- c("aucell", fgsea", "gsva", "mdt", "mlm", "ora", "udt", "ulm", "viper", "wmean", "wsum")
  stats <- c("ulm", "mlm", "viper") # 3 best methods to get consensus
  pathway_df <- decoupleR::decouple(mat = input, network = pathway_net,
                                    statistics = stats, minsize = 5)
  
  tf_df <- decoupleR::decouple(mat = input, network = tf_net,
                               statistics = stats, minsize = 5)
  
  # Remove insignificant entries
  pathway_sig <- pathway_df %>%
    dplyr::group_by(statistic) %>%
    dplyr::mutate(padj = stats::p.adjust(p_value, method = "BH")) %>%
    dplyr::ungroup() %>% 
    dplyr::filter(padj <= 0.05)
  
  tf_sig <- tf_df %>%
    dplyr::group_by(statistic) %>%
    dplyr::mutate(padj = stats::p.adjust(p_value, method = "BH")) %>%
    dplyr::ungroup() %>% 
    dplyr::filter(padj <= 0.05)
  
  # Return results 
  result <- list(all_pathways = pathway_df,
                 all_tfs = tf_df,
                 sig_pathways = pathway_sig,
                 sig_tf = tf_sig)
  
  return(invisible(result))
} 

# No "_" in Description column so that str_wrap works
# Collection column needed
# (combined_score, GeneRatio, k) OR (NES, leading_edge_size) columns needed
plot_pathway <- function(pathway_df, vst_counts, meta_data, samples, method, output_path){
  
  set.seed(1234)
  
  if (nrow(pathway_df) == 0){
    message("Input data frame is empty. Skipping plotting.")
    return(NULL)
  }
  
  if (is.null(method) || method == "") {
    stop("⚠️ 'method' must be defined before creating output directories.")
  }
  
  # ---- Create output directories ----
  bar_output_path <- file.path(output_path, method, "Bar_Plots")
  dot_output_path <- file.path(output_path, method, "Dot_Plots")
  heatmap_output_path <- file.path(output_path, method, "Heatmaps")
  dirs <- c(output_path,
            bar_output_path,
            dot_output_path,
            heatmap_output_path)
  
  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE)
      message("📂 Created directory: ", d)
    }
  }
  
  # ----   ----
  dot_plot_list <- list()
  bar_plot_list <- list()
  plot_colors <- c("Upregulated" = "#E69F00", "Downregulated" = "#56B4E9")
  
  # str_wrap() wraps only between words, and it defines "words" based on spaces (" ").
  # If all words are connected by "_", it wont split.
  pathway_df <- pathway_df %>%
    dplyr::mutate(Description = gsub(pattern = "_", replacement = " ", x = Description),
                  Description = stringr::str_wrap(string = Description, width = 30))
  
  # Decide if there are multiple collections or single collection.
  # If multiple, plot each collection separately,
  collections <- unique(pathway_df$Collection)
  n_collections <- length(unique(pathway_df$Collection))
  
  for (i in seq_len(n_collections)){
    
    # Save heatmap for pathways in each collection
    heatmap_plot_list <- list()
    
    plot_df <- pathway_df %>% dplyr::filter(Collection == collections[i])
    
    # Choose x axis column, size column and labels
    if (method == "ORA"){
      x_col <- sym("GeneRatio")
      x_label <- "GeneRatio"
      size_col <- sym("leading_edge_size")
      color_col <- sym("Direction")
      score_var <- dplyr::case_when("GeneRatio" %in% colnames(pathway_df) ~ "GeneRatio",
                                    "combined_score" %in% colnames(pathway_df) ~ "combined_score", 
                                    TRUE ~ NA_character_)
      x_limits <- c(0, NA)  # start at 0, auto end
      
      pathway_df <- pathway_df %>%
        dplyr::filter(!is.na(.data[[score_var]])) %>%
        dplyr::arrange(dplyr::desc(.data[[score_var]]))
      
    } else if (method == "GSEA"){
      x_col <- sym("NES")
      x_label <- "Normalized Enrichment Score (NES)"
      size_col <- sym("leading_edge_size")
      color_col <- sym("Direction")
      score_var <- "NES"
      x_min <- ifelse(floor(min(plot_df$NES, na.rm = TRUE)) > 0, 0, floor(min(plot_df$NES, na.rm = TRUE)))
      x_limits <- c(x_min, NA)
      
    }
    
    # Pad with empty rows if fewer than 15 pathways 
    n_missing <- 20 - nrow(plot_df)
    if(n_missing > 0){
      empty_df <- matrix(data = "", 
                         nrow = n_missing,
                         ncol = ncol(plot_df)) %>%
        as.data.frame() %>%
        dplyr::mutate(Description = paste0("", seq_len(n_missing)))
      plot_df <- dplyr::bind_rows(plot_df, empty_df)
    }
    
    max_label_len <- max(nchar(plot_df$Description), na.rm = TRUE)
    y_text_size <- dplyr::case_when(max_label_len > 50 ~ 6,
                                    max_label_len > 35 ~ 7,
                                    max_label_len > 25 ~ 8,
                                    TRUE ~ 10)
    
    # ---- Plot bar plot ---- 
    p1 <- ggplot2::ggplot(data = plot_df,
                          aes(x = !!x_col,
                              y = reorder(Description, !!x_col),
                              fill = !!color_col,
                              alpha = -log10(padj))) +
      ggplot2::geom_col(width = 0.75, na.rm = TRUE) +
      ggplot2::theme_classic() +
      ggplot2::labs(x = x_label,
                    y = "",
                    title = paste("Top", collections[i], "Pathways"),
                    fill = "Direction") +
      custom_theme +
      theme(axis.text.y = element_text(size = y_text_size)) +
      coord_cartesian(clip = "off") +
      scale_x_continuous(limits = x_limits, expand = expansion(mult = c(0, 0.05))) +
      scale_alpha_continuous(range = c(0.5, 1)) +
      scale_fill_manual(values = plot_colors) +
      guides(fill = guide_legend(override.aes = list(shape = 22, size = 6)),
             color = guide_legend(override.aes = list(shape = 22, size = 6)),
             alpha = guide_legend(override.aes = list(shape = 22, size = 6))) +
      ggplot2::geom_text(aes(label =  !!size_col), x = 0, hjust = -0.1, size = 3, show.legend = FALSE)
    
    bar_plot_list <- c(bar_plot_list, list(p1))
    
    # ---- Plot dot plot ---- 
    vals <- c(min(plot_df[[size_col]], na.rm = TRUE), max(plot_df[[size_col]], na.rm = TRUE))
    breaks <- as.vector(floor(quantile(vals) / 10) * 10)
    
    p2 <- ggplot2::ggplot(data = plot_df,
                          aes(x = !!x_col,
                              y = reorder(Description, !!x_col),
                              fill = !!color_col,
                              alpha = -log10(padj),
                              color = !!color_col,
                              size = !!size_col)) +
      ggplot2::geom_point() +
      ggplot2::theme_classic() +
      ggplot2::labs(x = x_label ,
                    y = "",
                    title = paste("Top", collections[i], "Pathways"),
                    color = "Direction",
                    size = "Counts") +
      custom_theme +
      theme(axis.text.y = element_text(size = y_text_size)) +
      coord_cartesian(clip = "off") + 
      scale_x_continuous(limits = x_limits, expand = expansion(mult = c(0, 0.05))) +
      scale_alpha_continuous(range = c(0.5, 1)) +
      scale_color_manual(values = plot_colors) +
      scale_fill_manual(values = plot_colors) + # need for coloring the legend
      guides(fill = guide_legend(override.aes = list(shape = 22, size = 6)),
             color = guide_legend(override.aes = list(shape = 22, size = 6)),
             alpha = guide_legend(override.aes = list(shape = 15, size = 6))) +
      ggplot2::scale_size(breaks = breaks) 
    
    dot_plot_list <- c(dot_plot_list, list(p2))
    
    # ---- Plot heatmap ----
    
    pathways <- pathway_df %>% 
      dplyr::filter(Collection == collections[i]) %>%
      dplyr::pull(Description) %>%
      unique()
    
    for (pathway in pathways) {
      plot_genes <- pathway_df %>% 
        dplyr::filter(Collection == collections[i], Description == pathway) %>%
        dplyr::select(dplyr::starts_with("gene", ignore.case = FALSE)) %>%
        unlist(use.names = FALSE) %>%
        trimws() %>%
        na.omit() %>% 
        unique()
      
      if (length(plot_genes) >= 2){
        norm_counts <- vst_counts[plot_genes, samples, drop = FALSE]
        metadata_col <- meta_data %>% dplyr::filter(Sample_ID %in% samples)
        metadata_row <- NULL
        disp_genes <- c()
        proj.params$heatmap$title <- stringr::str_wrap(string = pathway, width = 30)
        
        ph <- plot_heatmap(norm_counts, proj.params, metadata_col, metadata_row, disp_genes)
        heatmap_plot_list[[length(heatmap_plot_list) + 1]] <- ph$ph$gtable
      }
    }
    
    if (length(heatmap_plot_list) > 0) {
      pdf(file.path(heatmap_output_path, paste0("Heatmaps_", collections[i], "_", method, ".pdf")),
          width = 10, height = 8)
      for (ht in heatmap_plot_list) {
        grid::grid.newpage()
        grid::grid.draw(ht)
      }
      dev.off()
    }
  }
  
  bar_plots <- cowplot::plot_grid(plotlist = bar_plot_list, ncol = 3, nrow = 3)
  ggplot2::ggsave(filename = file.path(bar_output_path, paste0("Bar_plot_pathways_", method, ".tiff")),
                  plot = bar_plots,
                  device = "jpeg",
                  width = 3*7,
                  height = 3*7,
                  units = "in",
                  dpi = 300,
                  bg = "white")
  
  dot_plots <- cowplot::plot_grid(plotlist = dot_plot_list, ncol = 3, nrow = 3)
  ggplot2::ggsave(filename = file.path(dot_output_path, paste0("Dot_plot_pathways_", method, ".tiff")),
                  plot = dot_plots,
                  device = "jpeg",
                  width = 3*7,
                  height = 3*7,
                  units = "in",
                  dpi = 300,
                  bg = "white")
}

plot_tf <- function(input, contrast, meta_data, samples, output_path, n_tfs = 20){
  
  set.seed(1234)
  
  # ---- Create output directories ----
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
    message("📂 Created directory: ", output_path)
  }
  
  # NOTE: wsum returns wsum, norm_wsum and corr_wsum.
  # wsum (DONT USE): Biased toward larger gene sets (more genes → bigger sum)
  # norm_wsum (DONT USE): Adjusts for pathway length so small and large gene sets are comparable
  # corr_sum (DONT USE): corrects for high correlation as it can make enrichment appear stronger
  
  target <- stringr::str_split(string = contrast, pattern = "-")[[1]][1]
  ref <- stringr::str_split(string = contrast, pattern = "-")[[1]][2]
  stats <- c("consensus", "ulm", "mlm", "viper") #unique(input$all_tfs$statistic)
  bar_plot_list <- list()
  heatmap_plot_list <- list()
  
  for (stat in stats) {
    
    if (length(unique(input$all_tfs[["condition"]])) == 1){
      top_tf <- input$all_tfs %>%
        dplyr::mutate(Direction = dplyr::case_when(score < 0 ~ "Downregulated",
                                                   score > 0 ~ "Upregulated",
                                                   TRUE ~ "No change")) %>%
        dplyr::group_by(statistic, Direction) %>%
        dplyr::slice_max(order_by = abs(score), n = n_tfs, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::filter(statistic == stat)
      
      p <- ggplot(data = top_tf, aes(x = reorder(source, score), y = score, fill = score)) +
        geom_col(width = 0.75, na.rm = TRUE) +
        scale_fill_gradient2(low = "darkblue", high = "indianred", mid = "whitesmoke", midpoint = 0) +
        labs(x = "", y = "Score", title = paste0("Top TFs (", stat, " method)"), fill = "Score") +
        custom_theme +
        coord_cartesian(clip = "off") +
        geom_text(label = paste0("Activated in ", target),
                  x = top_tf$source[which.max(top_tf$score)],
                  y = ceiling(max(top_tf$score)) + 1, hjust = 1, color = "indianred", fontface = "bold") +
        geom_text(label = paste0("Activated in ", ref),
                  x = top_tf$source[which.min(top_tf$score)],
                  y = ceiling(max(top_tf$score)) + 1, hjust = 0, color = "darkblue", fontface = "bold")
      
      bar_plot_list <- c(bar_plot_list, list(p))
    }
    
    else{
      top_tf <- input$all_tfs %>%
        dplyr::filter(statistic == stat) %>%
        dplyr::group_by(source) %>%
        dplyr::summarise(std = sd(score, na.rm = TRUE), .groups = "drop") %>%
        dplyr::slice_max(order_by = abs(std), n = n_tfs, with_ties = FALSE) %>%
        dplyr::pull(source)
      
      tf_mat <- input$all_tfs %>%
        dplyr::filter(statistic == stat) %>%
        tidyr::pivot_wider(id_cols = "condition", names_from = "source", values_from = "score") %>%
        tibble::column_to_rownames("condition") %>%
        dplyr::select(all_of(top_tf)) %>%
        as.matrix() %>%
        t()
      
      norm_counts <- tf_mat[top_tf, samples, drop = FALSE]
      metadata_col <- meta_data %>% dplyr::filter(Sample_ID %in% samples)
      metadata_row <- NULL
      disp_genes <- c()
      proj.params$heatmap$title <- paste0("Top TFs (", stat, ") method")
      
      ph <- plot_heatmap(norm_counts, proj.params, metadata_col, metadata_row, disp_genes)
      heatmap_plot_list[[length(heatmap_plot_list) + 1]] <- ph$ph$gtable
    }
  }
  
  # Save combined heatmaps
  if (length(heatmap_plot_list) > 0) {
    pdf(file.path(output_path, "Heatmap_TFs.pdf"),
        width = 10, height = 8)
    for (ht in heatmap_plot_list) {
      grid::grid.newpage()
      grid::grid.draw(ht)
    }
    dev.off()
  }
  
  # Save combined barplots
  if (length(bar_plot_list) > 0) {
    bar_plots <- cowplot::plot_grid(plotlist = bar_plot_list, ncol = 1, nrow = length(bar_plot_list))
    ggplot2::ggsave(filename = file.path(output_path, "Bar_plot_TFs.tiff"),
                    plot = bar_plots,
                    device = "jpeg",
                    width = 10,
                    height = 4 * 3,
                    units = "in",
                    dpi = 300,
                    bg = "white")
  }
  
  return(invisible(list(barplots = bar_plot_list, heatmaps = heatmap_plot_list)))
}

main_analysis <- function(meta_data, read_data, proj.params, trial) {
  
  # Compile raw counts if read_data not available
  if (is.null(read_data)) {
    if (dir.exists(proj.params$counts.dir)) {
      read_data <- merge_counts(proj.params)
    } else {
      stop("The specified count directory does not exist: ", proj.params$counts.dir)
    }
  }
  
  # Prepare input
  deseq2_inputs <- prepare_deseq2_input(meta_data, read_data, proj.params)
  
  contrasts <- proj.params$deseq2$contrast
  
  # ---- Visualization: PCA plot ----
  meta_data <- deseq2_inputs$meta_data
  read_data <- deseq2_inputs$read_data
  output_path <- proj.params$proj.dir
  plot_pca(meta_data, read_data, output_path)
  
  
  if(!trial){
    
    # Loop through contrasts
    for (n in seq_along(contrasts)) {
      
      # ---- Differential Expression (DESeq2) ----
      meta_data <- deseq2_inputs$meta_data
      read_data <- deseq2_inputs$read_data
      deseq2_results <- run_deseq2(meta_data, read_data, proj.params, n)
      
      # ---- Visualization: MA Plot ----
      dds <- deseq2_results$dds
      output_path <- proj.params$deseq2.dir[n]
      plot_ma(dds, output_path)
      
      # ---- Visualization: Volcano Plot ----
      DEGs_df <- deseq2_results$degs
      contrast <- proj.params$deseq2$contrasts[n]
      output_path <- proj.params$deseq2.dir[n]
      plot_volcano(DEGs_df, proj.params, contrast, output_path)
      
      
      # ---- Visualization: Heatmap ----
      contrast <- proj.params$deseq2$contrasts[n]
      DEGs_df <- deseq2_results$degs
      vst_counts <- deseq2_results$vst
      samples <- filter_samples_by_contrast(meta_data, contrast)
      samples <- intersect(colnames(vst_counts), samples)
      sig_genes <- DEGs_df %>% dplyr::filter(padj <= 0.05) %>% dplyr::pull(SYMBOL)
      
      if (length(sig_genes) != 0){
        
        norm_counts <- vst_counts[sig_genes, samples, drop = FALSE]
        metadata_col <- meta_data %>% dplyr::filter(Sample_ID %in% samples)
        metadata_row <- NULL
        disp_genes <- c()
        output_path <- proj.params$deseq2.dir[n]
        
        ph <- plot_heatmap(norm_counts, proj.params, metadata_col, metadata_row, disp_genes)
        
        jpeg(file.path(output_path, "Heatmap.jpeg"), width = 5, height = 7, units = "in", res = 300)
        gridExtra::grid.arrange(grobs = list(ph$ph$gtable), ncol = 1)
        dev.off()
        
        save_xlsx(ph$mat, file.path(output_path, "Heatmap_Matrix.xlsx"), "Heatmap_matrix", row_names = TRUE)
      } 
      
      # ---- Pathway Analysis ----
      
      # Perform pathway analysis
      output_path <- proj.params$pathway.dir[n]
      gsea_res_list <- pathway_analysis(DEGs_df, proj.params, output_path)
      
      # Plot Bar, Dot, Heatmap for top 10 Up & Down Pathways
      top_pathways_gsea <- gsea_res_list$consensus %>%
        dplyr::group_by(Collection, Consensus, Description) %>%
        dplyr::slice_min(order_by = match(method, c("FGSEA", "GSEA", "ORA")), n = 1) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(Collection, Consensus) %>%
        dplyr::slice_max(order_by = abs(NES), n = 10, with_ties = FALSE) %>%
        dplyr::ungroup()
      
      top_pathways_ora <- gsea_res_list$consensus %>%
        dplyr::filter(method == "ORA") %>%
        dplyr::group_by(Collection, Consensus) %>%
        dplyr::slice_min(order_by = padj, n = 10, with_ties = FALSE) %>%
        dplyr::ungroup()
      
      output_path <- proj.params$pathway.dir[n]
      samples <- filter_samples_by_contrast(meta_data, contrast)
      samples <- intersect(colnames(vst_counts), samples)
      plot_pathway(top_pathways_gsea, vst_counts, meta_data, samples, "GSEA", output_path)
      plot_pathway(top_pathways_ora, vst_counts, meta_data, samples, "ORA", output_path)
      
      # ---- Transcription Factor (TF) Analysis ----
      
      # Perform TF analysis on DEGs using t-statistics
      t_stats_mat <- DEGs_df %>% 
        as.data.frame() %>%
        dplyr::mutate(t = -log10(padj) * log2FoldChange) %>%
        dplyr::filter(!is.na(t)) %>%
        dplyr::select(SYMBOL, t) %>%
        tibble::column_to_rownames("SYMBOL") %>%
        as.matrix()
      
      tf_res_degs <- tf_analysis(t_stats_mat, proj.params$species)
      
      # Perform TF analysis on vst counts 
      samples <- filter_samples_by_contrast(meta_data, contrast)
      samples <- intersect(colnames(vst_counts), samples)
      norm_counts_sub <- vst_counts[, samples]
      
      tf_res_counts <- tf_analysis(norm_counts_sub, proj.params$species)
      
      # Plot Bar, Heatmap for top 20 Up and Down TFs
      contrast <- contrasts[n]
      output_path <- proj.params$tf.dir[n]
      samples <- filter_samples_by_contrast(meta_data, contrast)
      samples <- intersect(colnames(vst_counts), samples)
      plot_tf(tf_res_degs, contrast, meta_data, samples, output_path, n_tfs = 20)
      plot_tf(tf_res_counts, contrast, meta_data, samples, output_path, n_tfs = 20)
    }
    
    # ---- Visualization: Dispersion Estimates ----
    dds <- deseq2_results$dds
    output_path <- proj.params$proj.dir
    plot_dispersion(dds, output_path)
    
    # ---- VST Counts (Blind) ----
    # NOTE: vst counts are affected by design ONLY when blind=FALSE
    dds <- deseq2_results$dds
    output_path <- proj.params$proj.dir
    vsd_blind <- DESeq2::vst(dds, blind = TRUE)
    vst_counts_blind <- SummarizedExperiment::assay(vsd_blind) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("ID") %>%
      add_annotation() %>%
      dplyr::filter(!is.na(SYMBOL)) %>%
      dplyr::group_by(SYMBOL) %>%
      dplyr::summarize(across(.cols = where(is.numeric), .fns = mean, na.rm = TRUE), .groups = "drop") %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames("SYMBOL") %>%
      dplyr::select(-dplyr::starts_with("ENSEMBL"), -dplyr::starts_with("ENTREZ")) %>%
      as.matrix()
    
    save_xlsx(vst_counts_blind, file.path(output_path, "VST_counts_blind.xlsx"), "VST_blind", row_names = TRUE)
    
    # ---- Perform LRT Test ----
    dds <- deseq2_results$dds
    dds_LRT <- DESeq2::DESeq(dds, test = "LRT", reduced = ~1)
    res_LRT <- DESeq2::results(dds_LRT) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("ID") %>%
      add_annotation()
    
    # ---- Visualization: Heatmap ----
    output_path <- proj.params$proj.dir
    samples <- colnames(vst_counts_blind)
    sig_genes <- res_LRT %>% dplyr::filter(padj <= 0.05) %>% dplyr::pull(SYMBOL)
    norm_counts <- vst_counts_blind[sig_genes, samples, drop = FALSE]
    metadata_col <- meta_data %>% dplyr::filter(Sample_ID %in% samples)
    metadata_row <- NULL
    disp_genes <- c()
    proj.params$heatmap.title <- ""
    if (length(sig_genes) != 0){
      
      ph <- plot_heatmap(norm_counts, proj.params, metadata_col, metadata_row, disp_genes)
      
      jpeg(file.path(output_path, "Heatmap.jpeg"), width = 5, height = 7, units = "in", res = 300)
      gridExtra::grid.arrange(grobs = list(ph$ph$gtable), ncol = 1)
      dev.off()
      
      save_xlsx(ph$mat, file.path(output_path, "Heatmap_Matrix.xlsx"), "Heatmap_matrix", row_names = TRUE)
    }
  }
}

get_annotations <- function() {
  
  set.seed(1234)
  
  # Initialize 
  species_list <- c("Homo sapiens", "Mus musculus")
  annotations_list <- list()
  
  for (species in species_list) {
    
    # Connect to AnnotationHub and Fetch Ensembl DB 
    ah <- AnnotationHub::AnnotationHub()
    ah_db <- AnnotationHub::query(x = ah, 
                                  pattern = c(species, "EnsDb"), 
                                  ignore.case = TRUE)
    
    # Acquire the latest annotation files
    latest <- ah_db %>%
      mcols() %>%
      rownames() %>%
      tail(n=1)
    
    # Download the appropriate Ensembldb database
    edb <- ah[[latest]]
    
    # Extract ENSEMBL Annotations 
    ensembl <- ensembldb::genes(x = edb, 
                                return.type = "data.frame") %>%
      dplyr::rename(ENSEMBL_ID    = gene_id,
                    ENSEMBL_SYMBOL  = gene_name,
                    ENSEMBL_BIOTYPE  = gene_biotype,
                    START       = gene_seq_start,
                    END        = gene_seq_end,
                    CHR        = seq_name,
                    STRAND      = seq_strand,
                    DESCRIPTION    = description,
                    ENSEMBL_TRANSCRIPT = canonical_transcript) %>%
      dplyr::mutate(ENSEMBL_SYMBOL = dplyr::case_when(nchar(ENSEMBL_SYMBOL) == 0 ~ NA,
                                                      TRUE ~ ENSEMBL_SYMBOL)) %>%
      dplyr::select(ENSEMBL_ID, ENSEMBL_TRANSCRIPT, ENSEMBL_SYMBOL, ENSEMBL_BIOTYPE,
                    START, END, CHR, STRAND, DESCRIPTION)
    
    # Extract ENTREZ Annotations 
    org_db <- if (species == "Homo sapiens") org.Hs.eg.db else org.Mm.eg.db
    entrez <- AnnotationDbi::select(x = org_db,
                                    keys = AnnotationDbi::keys(org_db),
                                    columns = c("ENSEMBL", "SYMBOL", "GENETYPE")) %>%
      dplyr::rename(ENTREZ_ID    = ENTREZID,
                    ENSEMBL_ID   = ENSEMBL,
                    ENTREZ_SYMBOL  = SYMBOL,
                    ENTREZ_BIOTYPE = GENETYPE)
    
    # Merge Annotations 
    annotations <- dplyr::full_join(ensembl, entrez, by=c("ENSEMBL_ID"="ENSEMBL_ID")) %>%
      dplyr::select(ENSEMBL_ID, ENSEMBL_TRANSCRIPT, ENTREZ_ID, ENSEMBL_SYMBOL, 
                    ENTREZ_SYMBOL, ENSEMBL_BIOTYPE, ENTREZ_BIOTYPE,START, END,
                    CHR, STRAND, DESCRIPTION)
    
    # Store Output 
    annotations_list[[species]] <- annotations
  }
  
  # Return: list(human, mouse) 
  return(invisible(annotations_list))
}

add_annotation <- function(normalized_counts) {
  
  set.seed(1234)
  
  # Retrieve Annotations 
  annotations <- get_annotations()
  
  # Flatten all annotation data into one named list 
  named_lists <- list(ensembl_id_human   = annotations$`Homo sapiens`$ENSEMBL_ID,
                      entrez_id_human   = annotations$`Homo sapiens`$ENTREZ_ID,
                      ensembl_symbol_human = annotations$`Homo sapiens`$ENSEMBL_SYMBOL,
                      entrez_symbol_human = annotations$`Homo sapiens`$ENTREZ_SYMBOL,
                      ensembl_id_mouse   = annotations$`Mus musculus`$ENSEMBL_ID,
                      entrez_id_mouse   = annotations$`Mus musculus`$ENTREZ_ID,
                      ensembl_symbol_mouse = annotations$`Mus musculus`$ENSEMBL_SYMBOL,
                      entrez_symbol_mouse = annotations$`Mus musculus`$ENTREZ_SYMBOL)
  
  # Compute intersection counts 
  overlap_counts <- sapply(X = named_lists, 
                           FUN = function(x) {length(intersect(x, normalized_counts$ID))})
  
  # Identify the best matching column 
  best_match <- names(which.max(overlap_counts))
  message("Best match: ", best_match)
  message("Overlap counts:\n", paste(names(overlap_counts), overlap_counts, sep = " \t: ", collapse = "\n"))
  
  # Select Species-Specific Annotation 
  if (grepl(pattern = "human", x = best_match)) {
    ann <- annotations$`Homo sapiens`
  } else if (grepl(pattern = "mouse", x = best_match)) {
    ann <- annotations$`Mus musculus`
  } else {
    stop("Unable to determine organism (human/mouse) from ID column.")
  }
  
  # Select ID and SYMBOL Columns 
  if (grepl(pattern = "ensembl", x = best_match)) {
    id_col <- "ENSEMBL_ID"
    symbol_col <- "ENSEMBL_SYMBOL"
  } else if (grepl(pattern ="entrez", x = best_match)) {
    id_col <- "ENTREZ_ID"
    symbol_col <- "ENTREZ_SYMBOL"
  } else {
    stop("Unable to determine ID type (Ensembl/Entrez) from ID column.")
  }
  
  # Finalize Columns 
  keep_ids <- c(id_col, symbol_col)
  drop_ids <- setdiff(colnames(ann), keep_ids)
  
  # Join Annotations and Create SYMBOL Column 
  normalized_counts <- ann %>%
    dplyr::right_join(normalized_counts, by = stats::setNames("ID", id_col), multiple = "all") %>%
    dplyr::mutate(SYMBOL = dplyr::case_when(is.na(ENSEMBL_SYMBOL) & !is.na(ENSEMBL_ID) ~ ENSEMBL_ID,
                                            is.na(ENSEMBL_SYMBOL) & !is.na(ENTREZ_ID) ~ ENTREZ_ID,
                                            TRUE ~ ENSEMBL_SYMBOL)) %>%
    dplyr::select(SYMBOL, all_of(keep_ids), dplyr::everything(), -all_of(drop_ids)) %>%
    dplyr::distinct_at(id_col, .keep_all = TRUE)
  
  # Return norm counts with annotation added 
  return(normalized_counts)
}

norm_counts_DESeq2 <- function(meta_data, read_data, proj.params) {
  
  set.seed(1234)
  
  # Create DESeq2 Dataset 
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = read_data,
                                        colData = meta_data,
                                        design = ~1)
  dds <- DESeq2::estimateSizeFactors(object = dds)
  
  # Normalized Counts 
  normalized_counts <- DESeq2::counts(dds, normalized = TRUE)
  
  # VST-Transformed Counts 
  vsd <- DESeq2::vst(dds, blind = TRUE)
  vst_counts <- SummarizedExperiment::assay(vsd)
  
  # Batch Correction (if applicable) 
  if ("Batch" %in% colnames(meta_data) && length(unique(meta_data$Batch)) > 1) {
    normalized_counts_batch <- limma::removeBatchEffect(x = log2(normalized_counts + 1),
                                                        batch = dds$Batch)
  } else {
    normalized_counts_batch <- normalized_counts
    message("No batch correction performed")
  }
  
  # Convert to Data Frames 
  normalized_counts_df <- normalized_counts %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID") %>%
    add_annotation()
  
  normalized_counts_batch_df <- normalized_counts_batch %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID") %>%
    add_annotation()
  
  vst_counts_df <- vst_counts %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID") %>%
    add_annotation()
  
  # Write to Excel 
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "VST_counts(blind=TRUE)")
  openxlsx::writeData(wb, sheet = "VST_counts(blind=TRUE)", x = vst_counts_df, rowNames = FALSE)
  
  openxlsx::addWorksheet(wb, sheetName = "Norm_counts")
  openxlsx::writeData(wb, sheet = "Norm_counts", x = normalized_counts_df, rowNames = FALSE)
  
  openxlsx::addWorksheet(wb, sheetName = "Norm_counts_batch_corrected")
  openxlsx::writeData(wb, sheet = "Norm_counts_batch_corrected", x = normalized_counts_batch_df, rowNames = FALSE)
  
  openxlsx::saveWorkbook(wb, file = file.path(proj.params$dir, proj.params$proj, "Normalized.counts.DESeq2.xlsx"), overwrite = TRUE)
  
  
  return(invisible(NULL))
}

add_major_pathway <- function(df){
  
  set.seed(1234)
  
  df <- df %>% 
    dplyr::mutate(MAJOR_PATHWAY = dplyr::case_when(
      grepl(pattern = "hallmark", x = base::tolower(Description)) ~ "HALLMARK",
      grepl(pattern = "interferon|interleukin|cytokine|chemokine|immune|toll|antigen|leukocyte|lymphocyte|macrophage", x = base::tolower(Description)) ~ "IMMUNE RELATED",
      grepl(pattern = "metaboli|purine|carbohydrate", x = base::tolower(Description)) ~ "METABOLISM",
      grepl(pattern = "translation", x = base::tolower(Description)) ~ "PROTEIN REGULATION",
      grepl(pattern = "transcription", x = base::tolower(Description)) ~ "GENE REGULATION",
      grepl(pattern = "mitotic|cell cycle", x = base::tolower(Description)) ~ "CELL CYCLE",
      grepl(pattern = "muscle", x = base::tolower(Description)) ~ "MUSCLE",
      grepl(pattern = "cardiac", x = base::tolower(Description)) ~ "HEART",
      grepl(pattern = "angiogenesis|blood_vessel", x = base::tolower(Description)) ~ "ANGIOGENESIS",
      grepl(pattern = "actin_", x = base::tolower(Description)) ~ "CYTOSKELETAN",
      grepl(pattern = "glycosyl", x = base::tolower(Description)) ~ "GLYCOSYLATION",
      grepl(pattern = "dna_", x = base::tolower(Description)) ~ "DNA DAMAGE/REPAIR",
      grepl(pattern = "rna_", x = base::tolower(Description)) ~ "RNA REGULATION",
      grepl(pattern = "ubiquitin|proteasome", x = base::tolower(Description)) ~ "PROTEIN DEGRATION",
      grepl(pattern = "transport", x = base::tolower(Description)) ~ "LOCALIZATION",
      grepl(pattern = "phagy|apopto", x = base::tolower(Description)) ~ "CELL DEATH",
      grepl(pattern = "ribosom", x = base::tolower(Description)) ~ "RIBOSOME",
      grepl(pattern = "gtpase", x = base::tolower(Description)) ~ "GTPASE",
      TRUE ~ "UNCLASSIFIED")) %>%
    dplyr::select(MAJOR_PATHWAY, everything())
  
  return(invisible(df))
}

### PROGENy analysis [NOT RECOMMENDED as decoupleR is better]
progeny_analysis <- function(norm_counts, assay = "RNA", species = "Homo sapiens") {
  
  set.seed(1234)
  
  # --- Error checking ---
  if (!is.matrix(norm_counts)) stop("`norm_counts` must be a numeric matrix (genes x samples).")
  if (!species %in% c("Homo sapiens", "Mus musculus")) stop("`species` must be 'Homo sapiens' or 'Mus musculus'.")
  
  organism <- dplyr::if_else(species == "Homo sapiens", "Human", "Mouse")
  
  # --- Raw PROGENy scores (no permutation) ---
  progeny_scores <- progeny(expr = norm_counts,
                            scale = FALSE,  # Refer section below to understand how it works
                            organism = organism,
                            top = 500,
                            perm = 1,
                            verbose = FALSE,
                            z_scores = FALSE,
                            get_nulldist = FALSE,
                            assay_name = assay,
                            return_assay = FALSE)
  
  # --- Significance scores from permutations ---
  progeny_sig_scores <- progeny(expr = norm_counts,
                                scale = FALSE,  # Refer section below to understand how it works
                                organism = organism,
                                top = 100,
                                perm = 10000,  # Returns list where [[1]] has p values
                                verbose = FALSE,
                                z_scores = FALSE,
                                get_nulldist = FALSE, # TRUE swaps rows and columns
                                assay_name = assay,
                                return_assay = FALSE)
  
  # --- Empirical two-sided p-values ---
  progeny_pvals <- (1-abs(progeny_sig_scores))/2
  
  # --- Adjusted p-values (FDR) ---
  pval_vec <- as.vector(progeny_pvals)
  padj_vec <- stats::p.adjust(pval_vec, method = "BH")
  
  progeny_padj <- matrix(data = padj_vec, 
                         nrow = nrow(progeny_pvals), 
                         ncol = ncol(progeny_pvals),
                         dimnames = dimnames(progeny_pvals))
  
  # --- Return as named list ---
  progeny_list <- list(scores = progeny_scores,
                       padj = progeny_padj,
                       significance = progeny_sig_scores,
                       pval = progeny_pvals)
  
  return(invisible(progeny_list))
  
  # NOTE: Below section is for understanding how scale parameter works
  if (FALSE) {
    
    # Sample expr data
    human_input <- as.matrix(read.csv(system.file("extdata", "human_input.csv", 
                                                  package = "progeny"), 
                                      row.names = 1))
    mouse_input <- as.matrix(read.csv(system.file("extdata", "mouse_input.csv", 
                                                  package = "progeny"), 
                                      row.names = 1))
    
    # Expected pathway scores with scale=TRUE and top = 10
    human_def_expected <- read.csv(system.file("extdata", "human_def_expected.csv", 
                                               package = "progeny"),
                                   row.names = 1,
                                   check.names = FALSE)
    mouse_def_expected <- read.csv(system.file("extdata", "mouse_def_expected.csv", 
                                               package = "progeny"),
                                   row.names = 1, 
                                   check.names = FALSE)
    
    # With scaling
    scores_scaled <- progeny(human_input, scale = TRUE, organism = "Human", top = 10)
    scores_scaled_subset <- progeny(human_input[,1:5], scale = TRUE, organism = "Human", top = 10)
    
    # Without scaling
    scores_raw <- progeny(human_input, scale = FALSE, organism = "Human", top = 10)
    scores_raw_subset <- progeny(human_input[,1:5], scale = FALSE, organism = "Human", top = 10)
    scores_raw_perm <- progeny(human_input, scale = FALSE, organism = "Human", top = 10, perm=10)
    
    # Check mean and sd
    apply(X=scores_scaled, MARGIN=2, FUN=function(x){mean(x) %>% round(.,2)})
    apply(X=scores_scaled, MARGIN=2, FUN=function(x){sd(x) %>% round(.,2)})
    apply(X=scores_raw, MARGIN=2, FUN=function(x){mean(x) %>% round(.,2)})
    apply(X=scores_raw, MARGIN=2, FUN=function(x){sd(x) %>% round(.,2)})
    
    # Scale the raw scores
    scores_raw_scaled <- scale(x=scores_raw, center=TRUE, scale=TRUE)
    
    # Verify
    glimpse(scores_raw_scaled)
    glimpse(scores_scaled)
    identical(scores_raw_scaled %>% as.data.frame(), scores_scaled %>% as.data.frame())
    
    # [1] This shows scale parameter ONLY scales the calculated pathway scores.
    # and doesnt scale the input gene expression data.
    # [2] raw scores DONT CHANGE based on number of samples but scaled scores
    # vary if number of samples change.
    # RECOMMENDATION: Set scale=FALSE and scale the scores later when plotting
    
  }
}

# ---- SINGLE CELL & SPATIAL ANALYSIS RELATED FUNCTIONS ----

read_cellranger <- function(sample, path){
  
  set.seed(1234)
  
  # Construct the full data directory path
  data_dir <- base::file.path(path, sample)
  
  # Check if path exists
  if (!dir.exists(data_dir)) {
    stop("The directory '", data_dir, "' does not exist!")
  }
  
  # Check for HDF5 (.h5) files
  h5_files <- list.files(data_dir, pattern = "\\.h5$", full.names = TRUE, ignore.case = TRUE)
  
  # Read the feature-barcode matrix with gene symbols
  # gene.column = 1 → Ensembl IDs
  # gene.column = 2 → Gene symbols (required for mito/ribo/heme ratios)
  counts <- NULL
  if (length(h5_files) > 0) {
    # Use the first .h5 file found
    message("HDF5 file detected. Reading '", basename(h5_files[1]), "' for sample '", sample, "'.")
    counts <- tryCatch({
      Seurat::Read10X_h5(filename = h5_files[1], 
                         use.names = TRUE, 
                         unique.features = TRUE)
    }, error = function(e) {
      stop("Failed to read HDF5 matrix from '", h5_files[1], "': ", e$message)
    })
  } else {
    # Fall back to standard 10X directory
    message("No HDF5 file found. Reading standard 10X matrix from directory.")
    counts <- tryCatch({
      Seurat::Read10X(data.dir = data_dir,
                      gene.column = 2,  # gene symbols
                      cell.column = 1,
                      unique.features = TRUE,
                      strip.suffix = FALSE
      )
    }, error = function(e) {
      stop("Failed to read 10X matrix from '", data_dir, "': ", e$message)
    })
  }
  
  # Create Seurat object
  sample_seurat <- SeuratObject::CreateSeuratObject(counts = counts,
                                                    project = sample,
                                                    assay = "RNA",
                                                    names.field = 1,
                                                    names.delim = "_",
                                                    meta.data = NULL,
                                                    min.cells = 0,
                                                    min.features = 0)
  
  # Log matrix type
  if (length(h5_files) > 0) {
    message("HDF5 Feature Barcode Matrix imported for '", sample, "'.")
  } else if (grepl("raw", path, ignore.case = TRUE)) {
    message("Raw Feature Barcode Matrix imported for '", sample, "'.")
  } else if (grepl("filt", path, ignore.case = TRUE)) {
    message("Filtered Feature Barcode Matrix imported for '", sample, "'.")
  } else {
    message("Feature-barcode matrix (raw/filtered) imported for '", sample, "'.")
  }
  
  return(invisible(sample_seurat))
}

read_spaceranger <- function(sample, bin, path){
  
  set.seed(1234)
  
  # Construct the full data directory path
  data_dir <- base::file.path(path, sample)
  
  # Check if path exists
  if (!dir.exists(data_dir)) {
    stop("The directory '", data_dir, "' does not exist!")
  }
  
  # Define expected matrix file path
  bin_dir <- dplyr::case_when(b == 2 ~ "binned_outputs/square_002um",
                              b == 8 ~ "binned_outputs/square_008um",
                              b == 16 ~ "binned_outputs/square_016um")
  matrix_file <- base::file.path(data_dir, bin_dir, "filtered_feature_bc_matrix.h5")
  
  # Check that matrix file exists
  if (!file.exists(matrix_file)) {
    stop("Matrix file 'filtered_feature_bc_matrix.h5' not found in: ", file.path(data_dir, bin_dir))
  }
  # Create Seurat object with filtered matrix
  sample_seurat <- Seurat::Load10X_Spatial(data.dir = data_dir,
                                           filename = "filtered_feature_bc_matrix.h5",
                                           assay = "Spatial",
                                           slice = sample,
                                           bin.size = bin,
                                           filter.matrix = TRUE,
                                           to.upper = FALSE,
                                           image = NULL)
  
  # Add orig.ident manually as Load10X_Spatial() doesn't set project/sample name
  sample_seurat@meta.data <- sample_seurat@meta.data %>% 
    dplyr::mutate(orig.ident = sample)
  
  message("Filtered Barcode Matrix imported for ", sample, " at bin size ", bin, "µm.")
  
  return(invisible(sample_seurat))
}

mark_empty_droplets_dropletutils <- function(sample_seurat){ 
  
  set.seed(1234)
  
  # Convert to SingleCellExperiment object
  sce <- Seurat::as.SingleCellExperiment(x = sample_seurat)
  
  # Set reproducible seed for emptyDrops
  set.seed(100)
  
  # NOTE: If FDR > 0.05 for some droplets AND Limited == TRUE, it indicates that
  # with more iterations, the FDR of these droplets can be reduced.
  
  # Iteratively run emptyDrops until limited droplets with high FDR are resolved
  niters <- 10000
  n_improve <- 1  # trigger loop
  
  while (n_improve > 0){
    e.out <- DropletUtils::emptyDrops(m = SingleCellExperiment::counts(sce),
                                      niters = niters)
    
    df <- as.data.frame(e.out)
    n_improve <- nrow(df %>%
                        dplyr::filter(Limited == TRUE, FDR > 0.05))
    message("emptyDrops check: ", n_improve, " droplets need more iterations (FDR > 0.05 & Limited=TRUE); niters = ", niters)
    niters <- niters + 10000
  }
  
  # Identify true (non-empty) cells with FDR ≤ 0.05
  true_cells <- df %>%
    dplyr::filter(FDR <= 0.05) %>% 
    rownames()
  
  # Annotate metadata with droplet classification
  sample_seurat@meta.data <- sample_seurat@meta.data %>%
    dplyr::mutate(Cell = rownames(.)) %>%
    dplyr::mutate(DropletUtils = dplyr::case_when(Cell %in% true_cells ~ "Non-Empty Droplet",
                                                  TRUE ~ "Empty Droplet"))
  
  # Log output
  ident <- as.character(unique(sample_seurat@meta.data$orig.ident))
  message("DropletUtils classification complete for sample: '", ident, "'")
  
  return(invisible(sample_seurat))
}

mark_empty_droplets_cellranger <- function(sample_seurat, filt_matrix_path){
  
  set.seed(1234)
  
  # Get sample name
  sample <- sample_seurat@meta.data$orig.ident %>% unique() %>% as.character()
  
  # Read the filtered barcode-feature matrix output of CellRanger
  sample_seurat_filt <- read_cellranger(sample, filt_matrix_path)
  
  # Mark cells absent in filtered barcode-feature matrix as empty droplets
  sample_seurat@meta.data <- sample_seurat@meta.data %>%
    dplyr::mutate(Cell = rownames(.)) %>%
    dplyr::mutate(CellRanger = dplyr::case_when(Cell %in% colnames(sample_seurat_filt) ~ "Non-Empty Droplet",
                                                TRUE ~ "Empty Droplet"))
  
  # Log output
  message("CellRanger empty droplets identified for sample: '", sample, "'")
  
  return(sample_seurat)
}

doublet_finder <- function(sample_seurat){
  
  set.seed(1234)
  
  # Populate missing columns with NA
  required_cols <- c("DropletUtils", "CellRanger")
  
  for (col in required_cols) {
    if (!col %in% colnames(sample_seurat@meta.data)) {
      sample_seurat@meta.data[[col]] <- NA  # assign NA if missing
    }
  }
  
  # Filter out empty droplets before doublet identification
  subset_seurat <- subset(x = sample_seurat,
                          subset = (DropletUtils == "Empty Droplet" & CellRanger == "Empty Droplet"),
                          invert = TRUE)
  
  # Preprocess
  subset_seurat <- subset_seurat %>% 
    Seurat::NormalizeData() %>% 
    Seurat::FindVariableFeatures() %>%
    Seurat::ScaleData() %>%
    Seurat::RunPCA()
  
  # Find significant PCs
  stdev_pc <- subset_seurat@reductions$pca@stdev
  percent_stdev_pc <- (stdev_pc / sum(stdev_pc)) * 100
  cumulative_stdev_pc <- cumsum(percent_stdev_pc)
  
  pc1 <- which(cumulative_stdev_pc > 90 & percent_stdev_pc < 5)[1]
  pc2 <- sort(which((percent_stdev_pc[1:(length(percent_stdev_pc)-1)] - 
                       percent_stdev_pc[2:length(percent_stdev_pc)]) > 0.1),
              decreasing = TRUE)[1] + 1
  n_pcs <- min(pc1, pc2)
  
  # pK Identification (no prior info)
  # Introduces artificial doublets in varying proportions into real dataset,
  # preprocesses the data and calculates proportion of artificial nearest 
  # neighbors. Output is a list of proportions of artificial nearest neighbors
  # for varying combinations of pK and pN. Optimal pK is the max of bimodality
  # coefficient (BCmvn) distribution
  
  # Run UMAP, neighbors, clustering for DoubletFinder paramSweep
  subset_seurat <- subset_seurat %>%
    Seurat::RunUMAP(dims = 1:n_pcs) %>%
    Seurat::FindNeighbors(dims = 1:n_pcs) %>%
    Seurat::FindClusters(resolution = 0.1)
  
  # Find optimal pK for DoubletFinder
  sweep.res <- DoubletFinder::paramSweep(subset_seurat, PCs = 1:n_pcs, sct = FALSE)
  sweep.stats <- DoubletFinder::summarizeSweep(sweep.res, GT=FALSE)
  bcmvn <- DoubletFinder::find.pK(sweep.stats)
  
  optimal_pK <- bcmvn %>% 
    dplyr::slice_max(order_by = BCmetric, n=1) %>%
    dplyr::pull(pK) %>%
    as.numeric()
  #optimal_pK <- as.numeric(as.character(optimal_pK[[1]]))
  
  # Estimate expected homotypic doublets based on 10X multiplet rate
  # https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled
  # From the above link, we can see that the multiplet rate is 8*10^-6 per cell
  multiplet_rate <- 8e-6
  n_cells <- nrow(subset_seurat@meta.data)
  n_exp <- round(multiplet_rate * n_cells)
  
  # Adjust for homotypic doublets
  homotypic.prop <- DoubletFinder::modelHomotypic(subset_seurat@meta.data$seurat_clusters)
  n_exp_adj <- round(n_exp * (1 - homotypic.prop))
  
  # Run DoubletFinder
  subset_seurat <- DoubletFinder::doubletFinder(seu = subset_seurat,
                                                PCs = 1:n_pcs,
                                                pN = 0.25,  #default
                                                pK = optimal_pK,
                                                nExp = n_exp_adj)
  
  # Rename classification column to 'DoubletFinder'
  colnames(subset_seurat@meta.data)[grepl(pattern = "DF.classifications", x = colnames(subset_seurat@meta.data))] <- "DoubletFinder"
  
  # Prepare classification dataframe for merging
  doublet_df <- subset_seurat@meta.data %>% 
    dplyr::select(DoubletFinder) %>%
    tibble::rownames_to_column(var = "Cell") %>%
    dplyr::mutate(DoubletFinder = stringr::str_to_title(DoubletFinder))
  
  # Merge classification back to original sample metadata
  sample_metadata <- sample_seurat@meta.data %>%
    dplyr::mutate(Cell = rownames(.)) %>%
    dplyr::left_join(doublet_df, by=c("Cell"="Cell")) %>%
    dplyr::mutate(DoubletFinder = dplyr::case_when(is.na(DoubletFinder) ~ "Empty Droplet",
                                                   TRUE ~ DoubletFinder)) %>%
    dplyr::mutate(index = Cell) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "index")
  
  # Assign updated metadata
  sample_seurat@meta.data <- sample_metadata
  
  # Log output
  ident <- as.character(unique(sample_seurat@meta.data$orig.ident))
  message("DoubletFinder doublets identified for sample: '", ident, "'")
  
  return(invisible(sample_seurat))
}

scdbl_finder <- function(sample_seurat){
  
  set.seed(1234)
  
  # Populate missing columns with NA
  required_cols <- c("DropletUtils", "CellRanger")
  
  for (col in required_cols) {
    if (!col %in% colnames(sample_seurat@meta.data)) {
      sample_seurat@meta.data[[col]] <- NA  # assign NA if missing
    }
  }
  
  # Filter out empty droplets before doublet identification
  subset_seurat <- subset(x = sample_seurat,
                          subset = (DropletUtils == "Empty Droplet" & CellRanger == "Empty Droplet"),
                          invert = TRUE)
  
  # Convert to SingleCellExperiment
  sce <- Seurat::as.SingleCellExperiment(x = subset_seurat)
  
  # Run scDblFinder
  scDbl <- scDblFinder::scDblFinder(sce = sce, 
                                    clusters = NULL,
                                    samples = NULL,
                                    dbr = NULL)
  
  # Extract classifications
  dbl_df <- scDbl@colData@listData %>% 
    data.frame() %>% 
    dplyr::rename(scDblFinder = scDblFinder.class) %>% 
    dplyr::mutate(Cell = scDbl@colData@rownames) %>%
    dplyr::mutate(scDblFinder = stringr::str_to_title(scDblFinder)) %>%
    dplyr::select(Cell, scDblFinder)
  
  # Merge into original metadata
  sample_metadata <- sample_seurat@meta.data %>%
    dplyr::mutate(Cell = rownames(.)) %>%
    dplyr::left_join(dbl_df, by=c("Cell"="Cell")) %>%
    dplyr::mutate(scDblFinder = dplyr::case_when(is.na(scDblFinder) ~ "Empty Droplet",
                                                 TRUE ~ scDblFinder)) %>%
    dplyr::mutate(index = Cell) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "index")
  
  # Assign updated metadata
  sample_seurat@meta.data <-  sample_metadata
  
  # Log output
  ident <- as.character(unique(sample_seurat@meta.data$orig.ident))
  message("scDblFinder doublets identified for sample: '", ident, "'")
  
  return(invisible(sample_seurat))
}

calc_qc_metrics_sc_sp <- function(sample_seurat, assay){
  
  set.seed(1234)
  
  # Compute mitochondrial percent
  sample_seurat <- Seurat::PercentageFeatureSet(object = sample_seurat,
                                                pattern = "^[Mm][Tt]-",
                                                col.name = "MitoPercent",
                                                assay = assay)
  
  # Compute ribosomal percent
  sample_seurat <- Seurat::PercentageFeatureSet(object = sample_seurat,
                                                pattern = "^[Rr][Pp][SsLl]", 
                                                col.name = "RiboPercent",
                                                assay = assay)
  
  # Compute hemoglobin percent
  sample_seurat <- Seurat::PercentageFeatureSet(object = sample_seurat,
                                                pattern = "^[Hh][Bb][AaBb]-",                                                 
                                                col.name = "HemePercent",
                                                assay = assay)
  
  # Rename columns to be more intuitive and add the QC metrics:
  # (i)    Cell      : unique identifiers corresponding to each cell i.e. barcodes
  # (ii)   Sample    : sample names
  # (iii)  nUMIs     : number of transcripts per cell
  # (iv)   nGenes    : number of genes per cell
  # (v)    nHTO_UMIs : number of HTO reads per cell
  # (vi)   nHTOs     : number of HTO types per cell
  # (vii)  MitoRatio : MitoPercent/100
  # (viii) RiboRatio : RiboPercent/100 
  # (ix)   HemeRatio : HemePercent/100
  # (x)    Novelty   : log ratio of genes per UMI
  
  sample_metadata <- sample_seurat@meta.data %>% 
    tibble::rownames_to_column(var = "barcode") %>%
    dplyr::mutate(Cell = paste0(orig.ident, "_", barcode),
                  Sample = orig.ident,
                  nUMIs = get(paste0("nCount_", assay)),
                  nGenes = get(paste0("nFeature_", assay)),
                  MitoRatio = MitoPercent / 100,
                  RiboRatio = RiboPercent / 100,
                  HemeRatio = HemePercent / 100,
                  Novelty = log10(nGenes) / log10(nUMIs))
  
  # Handle HTO metadata if present
  hto_cols <- c("nCount_HTO", "nFeature_HTO", "HTO_Final")
  if (any(hto_cols %in% colnames(sample_metadata))){
    sample_metadata <- sample_metadata %>% 
      dplyr::rename(nHTO_UMIs = nCount_HTO,
                    nHTOs = nFeature_HTO)
  }
  
  # Add spatial coordinates if available
  images <- names(sample_seurat@images)
  if (length(sample_seurat@images) > 0){
    df_coords <- data.frame()
    for (n in 1:length(images)){
      image <-  gsub(pattern=".002|0.008|.016|um", replacement="", x=images[n])
      df <- data.frame("barcode" = sample_seurat@images[[n]]@boundaries$centroids@cells, 
                       X = sample_seurat@images[[n]]@boundaries$centroids@coords[,1], 
                       Y = sample_seurat@images[[n]]@boundaries$centroids@coords[,2])
      
      df_coords <- dplyr::bind_rows(df_coords, df)
    }
    
    sample_metadata <- sample_metadata %>% 
      dplyr::left_join(df_coords, by=c("barcode"="barcode")) 
  }
  
  # Select columns to keep
  keep_cols <- c("Cell", "Sample", "nUMIs", "nGenes", "nHTO_UMIs", "nHTOs", 
                 "HTO_Final", "MitoRatio", "RiboRatio", "HemeRatio", "Novelty",
                 "DropletUtils", "CellRanger", "DoubletFinder", "scDblFinder", 
                 "X", "Y", "barcode")
  keep_cols <- base::intersect(keep_cols, colnames(sample_metadata))
  
  # Restore rownames and subset columns
  sample_metadata <- sample_metadata %>%
    dplyr::select(all_of(keep_cols)) %>%
    tibble::column_to_rownames(var = "barcode")
  
  # Assign updated metadata
  sample_seurat@meta.data <- sample_metadata
  
  # Log output
  ident <- as.character(unique(sample_seurat@meta.data$Sample))
  message("Cell-level QC metrics calculated for sample: '", ident, "'")
  
  return(invisible(sample_seurat))
}

mark_low_quality_sc_sp <- function(sample_seurat){
  
  set.seed(1234)
  
  # Select columns to keep
  required_cols <- c("Cell", "Sample", "nUMIs", "nGenes", "nHTO_UMIs", "nHTOs", 
                     "HTO_Final", "MitoRatio", "RiboRatio", "HemeRatio", "Novelty",
                     "DropletUtils", "CellRanger", "DoubletFinder", "scDblFinder", 
                     "X", "Y")
  
  for (col in required_cols) {
    if (!col %in% colnames(sample_seurat@meta.data)) {
      sample_seurat@meta.data[[col]] <- NA  # assign NA if missing
    }
  }
  
  # Define lenient cutoff thresholds for QC metrics
  if(length(names(sample_seurat@images)) == 0){
    
    # Single-cell cutoffs
    gene_cutoff <- 250
    umi_cutoff <- 500
    mito_cutoff <- 0.2
    ribo_cutoff <- 0.05
    novelty_cutoff <- 0.8
    
    sample_seurat@meta.data <- sample_seurat@meta.data %>% 
      dplyr::mutate(QC = dplyr::case_when((DropletUtils == "Empty Droplet" & CellRanger == "Empty Droplet") ~ "Empty Droplet",
                                          (DoubletFinder == "Doublet" & scDblFinder == "Doublet") ~ "Doublet",
                                          (nGenes >= gene_cutoff & 
                                             nUMIs >= umi_cutoff & 
                                             MitoRatio <= mito_cutoff & 
                                             Novelty >= novelty_cutoff) ~ "Singlet", #RiboRatio >= ribo_cutoff &
                                          TRUE ~ "Low Quality"))
  } else {
    
    # Spatial cutoffs (more lenient)
    gene_cutoff <- 25   # reduced from 250 used for single cell
    umi_cutoff <- 50    # reduced from 500 used for single cell
    mito_cutoff <- 0.2
    ribo_cutoff <- 0.05
    novelty_cutoff <- 0.8
    
    sample_seurat@meta.data <- sample_seurat@meta.data %>% 
      dplyr::mutate(QC = dplyr::case_when((nGenes >= gene_cutoff & 
                                             nUMIs >= umi_cutoff & 
                                             MitoRatio <= mito_cutoff & 
                                             Novelty >= novelty_cutoff) ~ "Singlet", #RiboRatio >= ribo_cutoff &
                                          TRUE ~ "Low Quality"))
  }
  
  # Log output
  ident <- as.character(unique(sample_seurat@meta.data$Sample))
  message("Good quality singlets identified for sample: '", ident, "'")
  
  return(invisible(sample_seurat))
}

generate_plotdata <- function(sample_seurat, raw_metadata){
  
  set.seed(1234)
  
  # Append metadata from current sample and filter out any rows with missing 
  # Sample info. This will be used for QC plots later
  raw_metadata <- dplyr::bind_rows(raw_metadata, sample_seurat@meta.data) %>%
    dplyr::filter(!is.na(Sample))
  
  # Log output
  ident <- as.character(unique(sample_seurat@meta.data$Sample))
  message("Raw metadata (for plotting) appended for sample: '", ident, "'")
  
  return(invisible(raw_metadata))
}

filter_singlets_sc_sp <- function(sample_seurat){
  
  set.seed(1234)
  
  # Check for required column
  if (!"QC" %in% colnames(sample_seurat@meta.data)) {
    stop("Missing 'QC' column in metadata. Please run `mark_low_quality_sc_sp()` first.")
  }
  
  # Subset to retain only singlets
  sample_seurat <- base::subset(x = sample_seurat, subset = QC == "Singlet")
  
  # Log output
  ident <- as.character(unique(sample_seurat@meta.data$Sample))
  message("Retained high-quality singlets for sample: ", ident)
  
  return(invisible(sample_seurat))
}

merge_filtered_sc_sp <- function(samples, assay, proj.params, output_path){
  
  set.seed(1234)
  
  # Create a merged seurat object after all the above filtering
  # NOTE: Samples will have the same barcodes. To keep track of cell identities
  # (i.e. barcodes) coming from each sample after merging, we add a prefix
  # (i.e. sample name) to each barcode using "add.cell.ids"
  samples.seurat <- lapply(samples, get)
  
  # Merge seurat objects with unique prefixes
  cell_prefixes <- gsub(pattern=".Spatial.*", replacement="", x=samples)
  merged_seurat <- base::merge(x = samples.seurat[[1]],   #get(paste0(samples[1])
                               y = samples.seurat[-1],    #lapply(paste0(samples[2:length(samples)]), get)
                               add.cell.ids = cell_prefixes,
                               merge.data = FALSE)
  
  # Remove HTO assay if present to avoid complications during integration
  if ("HTO" %in% Assays(merged_seurat)){
    merged_seurat[["HTO"]] <- NULL
  }
  
  # Add extra metadata if any to seurat object
  if (file.exists(proj.params$meta.file)) {
    extra_metadata <-  openxlsx::read.xlsx(xlsxFile = proj.params$meta.file) %>%
      dplyr::select(-dplyr::any_of("Comments"))
  } else{
    warning("Metadata file not found: ", proj.params$meta.file)
    extra_metadata <- data.frame()
  }
  
  if (nrow(extra_metadata) > 0){
    meta_data <- merged_seurat@meta.data
    
    if (assay == "RNA"){
      meta_data <- meta_data %>%
        dplyr::mutate(Unique_ID = dplyr::case_when(!is.na(HTO_Final) ~ paste0(Sample, "_", HTO_Final), 
                                                   is.na(HTO_Final) ~ paste0(Sample))) %>%
        dplyr::left_join(extra_metadata, by = ("Unique_ID" = "Unique_ID"))
    } else {
      meta_data <- meta_data %>%
        dplyr::left_join(extra_metadata, by = c("Sample" = "Slide")) %>% 
        dplyr::filter(dplyr::between(X, Xmin, Xmax) & dplyr::between(Y, Ymin, Ymax))
    }
    
    # Add row names before replacing metadata in Seurat object as dplyr::left_join() will remove row names.
    if (!"Cell" %in% colnames(meta_data)) stop("Column 'Cell' containing rownames is not found in metadata after dplyr join.")
    merged_seurat@meta.data <- meta_data %>%
      dplyr::mutate(index = Cell) %>%
      tibble::column_to_rownames(var = "index")  # VERY VERY IMPORTANT
  }
  
  # Save merged seurat object
  if (assay == "RNA"){
    filename <- file.path(output_path, "filtered.seurat.rds")
  } else{  # For Spatial.008um and Spatial.016um assays
    filename <- file.path(output_path, paste0("filtered.seurat.", assay, ".rds"))
  }
  saveRDS(merged_seurat, file = filename)
  
  # Log output
  message("Filtered Seurat object saved to: ", filename)
  
  return(merged_seurat)
}

plot_qc <- function(raw_metadata, output_path){
  
  set.seed(1234)
  
  ### Helper function definitions
  qc_levels <- c("Empty Droplet", "Doublet", "Low Quality", "Singlet")
  fill_colors <- c("Empty Droplet" = "#FFC61E", "Singlet" = "#2CA02C",
                   "Doublet" = "#1F77B4", "Low Quality" = "#D62728")
  
  ### Visualize cell counts per sample
  cell_qc <- function(meta){
    
    df <- meta %>%
      dplyr::count(Sample, QC) %>%
      data.frame() %>%
      dplyr::mutate(QC = factor(QC, levels = qc_levels))
    
    p <- ggplot(data = df, aes(x = Sample, y = n, fill = QC)) + 
      # position = "dodge" for grouped; "stack" for stacked
      # stat = "identity" if y axis defined; "count" if y axis determined based on X axis frequency
      geom_bar(stat = "identity", position = position_dodge(0.9), , drop = FALSE) +             
      theme_classic() +               # display with x and y axis lines and no gridlines
      custom_theme +
      labs(x = "Sample", y = "Cell Counts", title = "Number of Cells") +
      coord_cartesian(ylim = c(1,10000000), clip = "off", expand = FALSE) +
      scale_y_log10(breaks = c(10, 100, 1000, 10000, 100000, 1000000)) +
      scale_fill_manual(values = fill_colors) +
      geom_text(aes(label = n, ymin = 0.1, ymax = 1), 
                position = position_dodge(width = 0.9), y = 0.1, hjust = 0, angle = 90)
    #geom_text(stat ="count", aes(label = after_stat(count)), y = 0, hjust = 0, angle = 90)
    
    return(p)
  }
  
  ### Visualize nUMIs, nGenes, MitoRatio, RiboRatio, Novelty per sample
  violin_qc <- function(meta, yvar, ylab, title, cutoff = NULL, ylog = TRUE, ylim = NULL) {
    df <- meta %>%
      dplyr::mutate(QC = factor(QC, levels = qc_levels))
    
    p <- ggplot(data = df, aes(x = Sample, y = .data[[yvar]], fill = QC)) +
      geom_violin(position = position_dodge(0.9), scale = "width", drop = FALSE) +
      geom_boxplot(position = position_dodge(0.9), width = 0.15, outlier.size = 0.5, drop = FALSE) +
      theme_classic() + custom_theme +
      labs(x = "Sample", y = ylab, title = title) +
      scale_fill_manual(values = fill_colors)
    
    if (!is.null(cutoff)) p <- p + geom_hline(yintercept = cutoff, linetype = 2)
    if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim, clip = "off")
    if (ylog) p <- p + scale_y_log10()
    
    return(p)
  }
  
  ### Visualize number of genes/cell, number of UMIs/cell & MitoRatio together.
  # Bottom left quadrant : Poor quality cells with low genes & UMIs per cell 
  # Top right quadrant   : Good quality cells with high genes & UMIs per cell
  # Bottom right quadrant: Cells with low genes but high UMIs per cell. These 
  # could be dying cells or population of low complexity cells (i.e erythrocytes)
  gene_umi_mito_qc <- function(meta){
    
    umi_cutoff <- 500
    gene_cutoff <- 250
    df <- meta %>%
      dplyr::mutate(QC = factor(QC, levels = qc_levels))
    
    ggplot(data = df, aes(x = nUMIs, y = nGenes, color = MitoRatio)) +
      geom_point(alpha = 0.5) +
      theme_classic() + 
      custom_theme + 
      labs(x = "Number of UMIs", y = "Number of Genes",	 title = "UMIs vs Genes (Colored by MitoRatio)") +
      coord_cartesian(xlim = c(1, 1000000), ylim = c(1, 20000), clip = "off") +
      scale_x_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000)) + 
      scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000)) + 
      scale_color_viridis(option = "D", limits = c(0, 1)) + 		# limits sets max and min values of gradient
      #facet_wrap(.~Sample, nrow = 4) +   #split the plot by X-axis label
      facet_wrap(Sample ~ QC, ncol = 4) +
      geom_vline(xintercept = umi_cutoff, linetype = "dashed") +    	#draw a vertical line at x=500 i.e.UMIs cutoff
      geom_hline(yintercept = gene_cutoff, linetype = "dashed") +    #draw a horizontal line at y =250 i.e. Genes cutoff
      stat_smooth(method=lm, color="yellow", se = FALSE)
  }
  
  # === Plot generators
  plot_list <- list(
    Cell_Counts                = cell_qc,
    UMI_Distribution           = function(x) violin_qc(x, "nUMIs", "Number of UMIs", "UMI Distribution", 500, TRUE, c(1, 1e6)),
    Gene_Distribution          = function(x) violin_qc(x, "nGenes", "Number of Genes", "Gene Distribution", 250, TRUE, c(1, 30000)),
    MitoRatio_Distribution     = function(x) violin_qc(x, "MitoRatio", "MitoRatio", "MitoRatio Distribution", 0.2, TRUE, c(1e-5, 1)),
    RiboRatio_Distribution     = function(x) violin_qc(x, "RiboRatio", "RiboRatio", "RiboRatio Distribution", 0.05, TRUE, c(1e-4, 1)),
    Novelty_Score_Distribution = function(x) violin_qc(x, "Novelty", "Novelty", "Novelty Score Distribution", 0.8, FALSE, c(0.3, 1)),
    Genes_UMI_MitoRatio_Distribution = gene_umi_mito_qc)
  
  # === Generate and save plots
  for (plot_name in names(plot_list)) {
    p <- plot_list[[plot_name]](raw_metadata)   #  p <- get(funcs[i])(raw_metadata)
    ggplot2::ggsave(filename = paste0("QC_", plot_name, ".pdf"),
                    plot = p,
                    device = "pdf",
                    path = output_path,
                    width = 11,
                    height = 8,
                    dpi = 600,
                    units = "in")
  }
  
  # Log output
  message("QC plots generated and saved to: ", output_path)
}

sctransform_sc_sp <- function(filtered.seurat, assay, output_path){
  
  set.seed(1234)
  
  stopifnot(is(filtered.seurat, "Seurat"))
  if (!assay %in% names(filtered.seurat@assays)) {
    stop("Assay '", assay, "' not found in Seurat object.")
  }
  
  message("Normalizing data before cell cycle scoring...")
  filtered.seurat <- Seurat::NormalizeData(object = filtered.seurat,
                                           assay = assay,
                                           normalization.method = "LogNormalize",
                                           scale.factor = 10000,
                                           margin = 1,
                                           verbose = TRUE)
  
  message("Joining layers for cell cycle scoring...")
  # NOTE: CellCycleScoring uses a single data layer (i.e. log norm counts) but
  # currently, data layer for each sample is stored separately. Join them first.
  filtered.seurat@assays[[assay]] <- SeuratObject::JoinLayers(filtered.seurat@assays[[assay]])
  
  message("Scoring cell cycle...")
  filtered.seurat <- Seurat::CellCycleScoring(object = filtered.seurat,
                                              s.features = intersect(s_genes,rownames(filtered.seurat@assays[[assay]]@features)),
                                              g2m.features = intersect(g2m_genes, rownames(filtered.seurat@assays[[assay]]@features)),
                                              ctrl = NULL)
  
  # Calculate CC.Score to regress out the difference between G2M & S scores
  # https://satijalab.org/seurat/archive/v3.1/cell_cycle_vignette
  filtered.seurat$CC.Score <- filtered.seurat$G2M.Score-filtered.seurat$S.Score
  
  message("Splitting assay by sample before SCTransform...")
  # NOTE: All cells within same batch MUST be analyzed together
  filtered.seurat@assays[[assay]] <- base::split(x = filtered.seurat@assays[[assay]],
                                                 f = filtered.seurat@meta.data[["Sample"]])
  
  message("Running SCTransform with CC.Score and MitoRatio regression...") 
  # https://github.com/satijalab/seurat/issues/7342
  sct.seurat <- Seurat::SCTransform(object = filtered.seurat,
                                    assay = assay,
                                    new.assay.name = "SCT",
                                    do.correct.umi = TRUE,
                                    ncells = 5000,
                                    variable.features.n = 3000,
                                    vars.to.regress = c("CC.Score","MitoRatio"),
                                    do.scale = FALSE,
                                    do.center = TRUE,
                                    vst.flavor = "v2",
                                    return.only.var.genes = TRUE,
                                    verbose = TRUE)
  
  message("Filtering out ribosomal, mitochondrial, RIKEN, predicted genes from variable features...")
  # NOTE: Do this to so PCA, UMAP and clustering are not influenced by these genes.
  var_f <- sct.seurat@assays[["SCT"]]@var.features
  var_f <- var_f[!grepl(pattern = "^[Rr][Pp][SsLl]|R[Ii][Kk]$|^[Mm][Tt]-|^G[Mm][0-9.]+$", 
                        x = var_f)]
  sct.seurat@assays[["SCT"]]@var.features <- var_f
  cat("\nFinal number of variable features:", length(var_f), "\n")
  
  message("Scaling and running PCA on original assay...")
  # NOTE: Needed for scVI integration
  sct.seurat <- Seurat::ScaleData(object = sct.seurat, 
                                  assay = assay,
                                  features = VariableFeatures(sct.seurat))
  sct.seurat <- Seurat::RunPCA(object = sct.seurat,
                               assay = assay,
                               features = Seurat::VariableFeatures(sct.seurat),
                               reduction.name = paste0(assay, ".pca"),
                               reduction.key = "PC_")
  
  message("Running PCA on SCT assay...")
  sct.seurat <- Seurat::RunPCA(object = sct.seurat,
                               assay = "SCT",
                               features = Seurat::VariableFeatures(sct.seurat),
                               reduction.name = "sct.pca",
                               reduction.key = "PC_")
  
  # # Create .rds object for sct seurat object
  # if (assay == "RNA"){
  #   saveRDS(sct.seurat, file=paste0(output_path, "sct.seurat.rds"))
  # } else{
  #   # For Spatial.008um and Spatial.016um assays
  #   saveRDS(sct.seurat, file=paste0(output_path, "sct.seurat.", assay, ".rds"))
  # }
  
  message("SCTransform workflow completed successfully.")
  return(invisible(sct.seurat))
}

integrate_sc_sp <- function(sct.seurat, assay, reference.samples, output_path){
  
  set.seed(1234)
  
  stopifnot(is(sct.seurat, "Seurat"))
  if (!assay %in% names(sct.seurat@assays)) {
    stop("Assay '", assay, "' not found in Seurat object.")
  }
  if (!"sct.pca" %in% names(sct.seurat@reductions)) {
    stop("Reduction 'sct.pca' is missing. Please run PCA on SCT assay before integration.")
  }
  
  kweight <- min(sct.seurat@meta.data %>% 
                   dplyr::count(Sample) %>% 
                   dplyr::select(n) %>% 
                   min()/2, 100) 
  
  message("Reference samples: ", paste(reference.samples, collapse = ", "))
  message("k.weight: ", kweight)
  
  integrated.seurat <- sct.seurat
  
  # NOTE: While Seurat::CCAIntegration(), Seurat::RPCAIntegration() &
  # Seurat::JointPCAIntegration() all have dims=1:30 as default, only 
  # Seurat::JointPCAIntegration() gives error when integrated.seurat has fewer 
  # than 30 dims. So, we explicitly specify this.
  # Determine the maximum number of dimensions available for integration
  max_dims <- min(30, ncol(integrated.seurat@reductions$sct.pca@cell.embeddings))
  
  integration.methods <- c("CCA", "RPCA", "Harmony", "JointPCA")
  
  for (method in integration.methods){
    
    message("Running ", method, " integration...")
    
    reduction.name <- paste0("integrated.", tolower(method))
    
    integrated.seurat <- Seurat::IntegrateLayers(object = integrated.seurat,
                                                 method = paste0(method, "Integration"),
                                                 normalization.method = "SCT",
                                                 orig.reduction = "sct.pca", 
                                                 new.reduction =  reduction.name,
                                                 reference = reference.samples,
                                                 k.weight = kweight,    # for RPCA
                                                 dims = 1:max_dims,
                                                 verbose = TRUE)
  }
  
  # NOTE: scVI needs raw counts. Vignette also uses it on RNA assay
  # NOTE: We use variable features of SCT assay for integration.
  # NOTE: We use pca reduction from RNA assay (derived using variable features of SCT assay)
  # FastMNN throws error "Error in checkBatchConsistency(batches, cells.in.columns = TRUE)"
  # for (method in c("scVI", "FastMNN")){
  #   DefaultAssay(integrated.seurat) <- assay
  #   integrated.seurat <- Seurat::IntegrateLayers(object = integrated.seurat,
  #                                                method = paste0(r, "Integration"),
  #                                                normalization.method = "LogNormalize",
  #                                                orig.reduction = paste0(assay, ".pca"),
  #                                                features = integrated.seurat@assays$SCT@var.features,
  #                                                new.reduction = paste0("integrated.", base::tolower(method)),
  #                                                reference = reference.samples,
  #                                                k.weight = kweight,                                    # for RPCA
  #                                                conda_env = "/hpc/home/kailasamms/miniconda3/envs/R",  # for scVI
  #                                                verbose = FALSE)
  # }
  
  # Merge layers after integration
  integrated.seurat@assays[[assay]] <- SeuratObject::JoinLayers(integrated.seurat@assays[[assay]])
  
  # Optional RDS saving block (currently commented out)
  rds_path <- if (assay == "RNA") {
    file.path(output_path, "integrated.seurat.rds")
  } else {
    file.path(output_path, paste0("integrated.seurat.", assay, ".rds"))
  }
  saveRDS(integrated.seurat, file = rds_path)
  
  message("Integration completed.")
  return(invisible(integrated.seurat))
}

cluster_sc_sp <- function(integrated.seurat, assay, output_path){
  
  set.seed(1234)
  
  #***************STEP 8A: FIND NEAREST NEIGHBORS FOR EVERY CELL***************#
  
  # Determine the K-nearest neighbor graph
  for (r in c("CCA", "RPCA", "Harmony", "JointPCA")){
    max_dims <- min(40, ncol(integrated.seurat@reductions[[paste0("integrated.", base::tolower(r))]]@cell.embeddings))
    integrated.seurat <- Seurat::FindNeighbors(object = integrated.seurat,
                                               reduction = paste0("integrated.", base::tolower(r)),
                                               dims = 1:max_dims,
                                               k.param = 30,
                                               graph.name = c(paste0("graph_nn.", base::tolower(r)),
                                                              paste0("graph_snn.", base::tolower(r))))
  }
  
  #**********STEP 8B: SEPARATE CELLS INTO CLUSTERS BASED ON SNN GRAPH**********#
  
  # Determine the clusters for various resolutions
  for (r in c("CCA", "RPCA", "Harmony", "JointPCA")){
    for (res in c(0.4, 0.6, 0.8, 1, 1.2, 1.4)){
      integrated.seurat <- Seurat::FindClusters(object = integrated.seurat,
                                                resolution = res,
                                                graph.name = paste0("graph_snn.", base::tolower(r)),
                                                cluster.name = paste0("cluster.", res, ".", base::tolower(r)),
                                                modularity.fxn = 1,
                                                algorithm = 4)     #4=Leiden is best
    }
  }
  
  #**********STEP 8C: PERFORM DIMENSIONAL REDUCTION FOR VISUALIZATION**********#
  
  # Run UMAP
  for (r in c("CCA", "RPCA", "Harmony", "JointPCA")){ 
    max_dims <- min(40, ncol(integrated.seurat@reductions[[paste0("integrated.", base::tolower(r))]]@cell.embeddings))
    integrated.seurat <- Seurat::RunUMAP(object = integrated.seurat,
                                         dims = 1:max_dims,
                                         n.neighbors = 30L,
                                         reduction = paste0("integrated.", base::tolower(r)),
                                         reduction.name = paste0("umap.", base::tolower(r)))
  }
  
  # # Create .rds object for integrated seurat object
  # if (assay == "RNA"){
  #   saveRDS(integrated.seurat, file=paste0(output_path, "integrated.seurat.rds"))
  # } else{
  #   # For Spatial.008um and Spatial.016um assays
  #   saveRDS(integrated.seurat, file=paste0(output_path, "integrated.seurat.", assay, ".rds"))
  # }
  
  cat("Clustering completed", "\n")
  return(integrated.seurat)
}

remove_sparse_clusters_sc_sp <- function(integrated.seurat, assay, output_path){
  
  set.seed(1234)
  
  # Get all available resolutions at different reductions
  col_id <- colnames(integrated.seurat@meta.data %>% 
                       dplyr::select(starts_with("cluster.")))
  
  sparse_cells <- c()
  for (id in col_id){
    # Identify clusters which have less than 5 cells
    sparse_clusters <- integrated.seurat@meta.data %>%
      dplyr::count(get(id)) %>%
      dplyr::filter(n <=5) %>%
      dplyr::select(identity(1)) %>%
      unlist(.,use.names=FALSE) %>%
      as.character() %>%
      as.numeric()
    
    print(sparse_clusters)
    
    # Identify the cells in these clusters
    cells <- integrated.seurat@meta.data %>%
      dplyr::filter(get(id) %in% sparse_clusters) %>%
      dplyr::select(Cell) %>%
      unlist(.,use.names=FALSE)
    
    # Create a list of cells identified in sparse clusters at all resolutions 
    # and reductions
    sparse_cells <- c(sparse_cells, cells)
  }
  
  # Remove sparse_cells
  integrated.seurat <- subset(x=integrated.seurat,
                              subset = (Cell %in% unique(sparse_cells)),
                              invert=TRUE)
  
  cat("\nCells removed:", length(unique(sparse_cells)), "\n")
  
  # Create .rds object for integrated seurat object
  if (assay == "RNA"){
    saveRDS(integrated.seurat, file=file.path(output_path, "integrated.seurat.rds"))
  } else{
    # For Spatial.008um and Spatial.016um assays
    saveRDS(integrated.seurat, file=file.path(output_path, paste0("integrated.seurat.", assay, ".rds")))
  }
  
  cat("Integrated seurat object saved after removing sparse clusters (below 5 cells)", "\n")
  return(integrated.seurat) 
}

plot_metrics_post_integration_sc_sp <- function(integrated.seurat, assay, output_path){
  
  set.seed(1234)
  
  #                   Filename                       = c(reduction,             split.by,idents)  
  plot.params <- list(Pre.Integration.PCA.           = c("sct.pca",            "Sample", "cluster.0.8.harmony"),
                      Post.Integration.PCA.          = c("integrated.harmony", "Sample", "cluster.0.8.harmony"),
                      UMAP.Sample.                   = c("umap.harmony",       "Sample", "cluster.0.8.harmony"),
                      UMAP.Phase.                    = c("umap.harmony",       "Phase",  "cluster.0.8.harmony"),
                      UMAP.All.Resolutions.CCA       = c("umap.cca",            NA,      "All"),
                      UMAP.All.Resolutions.RPCA.     = c("umap.rpca",           NA,      "All"),
                      UMAP.All.Resolutions.JointPCA. = c("umap.jointpca",       NA,      "All"),
                      UMAP.All.Resolutions.Harmony.  = c("umap.harmony",        NA,      "All"),
                      UMAP.Singlets.Doublets.        = c("umap.harmony",        NA,      "cluster.0.8.harmony"),
                      UMAP.Numerical.Metrics.        = c("umap.harmony",        NA,      "cluster.0.8.harmony"))
  
  for (i in 1:length(plot.params)){
    
    if(names(plot.params)[i] %in% c("Pre.Integration.PCA.", 
                                    "Post.Integration.PCA.", "UMAP.Sample.", "UMAP.Phase.")){
      
      plot.seurat <- Seurat::SplitObject(object = integrated.seurat,
                                         split.by = plot.params[[i]][2])
      
      purrr::map(.x = c(1:length(plot.seurat)),
                 .f = function(x){  
                   Idents(plot.seurat[[x]]) <- plot.params[[i]][3]   # define name of cluster
                   Seurat::DimPlot(object = plot.seurat[[x]],
                                   reduction = plot.params[[i]][1],  # define reduction to use
                                   group.by = plot.params[[i]][3],   # define color of cluster
                                   pt.size = 0.1,
                                   order = TRUE,  # plot positive cells above negative cells
                                   label = TRUE,
                                   raster = FALSE,
                                   combine = TRUE) +
                     NoLegend() +
                     custom_theme + 
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
    }
    
    else if(grepl(pattern="All.Resolutions", x=names(plot.params)[i])){
      
      purrr::map(.x = c(0.4, 0.6, 0.8, 1, 1.2, 1.4),
                 .f = function(x){  
                   idents <- paste0("cluster.", x, gsub(pattern="umap", replacement="", x=plot.params[[i]][1]))
                   Idents(integrated.seurat) <- idents
                   Seurat::DimPlot(object = integrated.seurat,
                                   reduction = plot.params[[i]][1],
                                   group.by = idents,
                                   pt.size = 0.1,
                                   order = TRUE,  # plot positive cells above negative cells
                                   label = TRUE,
                                   raster = FALSE,
                                   combine = TRUE) +
                     NoLegend() +
                     custom_theme + 
                     ggplot2::labs(title = base::gsub(pattern="cluster.", replacement="", x=idents))
                 }) %>% cowplot::plot_grid(plotlist=.,
                                           align="hv",
                                           axis="tblr",
                                           nrow=2,
                                           ncol=3,
                                           rel_widths=1,
                                           rel_heights=1,
                                           greedy=TRUE,
                                           byrow=TRUE)
    }
    
    else if (names(plot.params)[i] == "UMAP.Singlets.Doublets."){
      
      purrr::map(.x = c("DropletUtils", "CellRanger", "DoubletFinder", "scDblFinder", "QC"),
                 .f = function(x){ 
                   Idents(integrated.seurat) <- plot.params[[i]][3]
                   Seurat::DimPlot(object = integrated.seurat,
                                   reduction =  plot.params[[i]][1],
                                   group.by = x,
                                   pt.size = 0.1,
                                   order = c("Doublet"),  # plot doublets on above rest of cells
                                   label = FALSE,
                                   raster = FALSE,
                                   combine = TRUE) +
                     #NoLegend() +
                     custom_theme + 
                     ggplot2::labs(title = x)
                 }) %>% cowplot::plot_grid(plotlist=.,
                                           align="hv",
                                           axis="tblr",
                                           nrow=2,
                                           ncol=3,
                                           rel_widths=1,
                                           rel_heights=1,
                                           greedy=TRUE,
                                           byrow=TRUE)
    }
    
    else if (names(plot.params)[i] == "UMAP.Numerical.Metrics."){
      
      purrr::map(.x = c("nUMIs", "nGenes", "S.Score", "G2M.Score", "CC.Score", "MitoRatio"),
                 .f = function(x){ 
                   Idents(integrated.seurat) <- plot.params[[i]][3]
                   Seurat::FeaturePlot(object = integrated.seurat,
                                       slot = "data",
                                       features = x,
                                       reduction = plot.params[[i]][1],
                                       pt.size = 0.1,
                                       order = TRUE, 
                                       label = FALSE,
                                       raster = FALSE,
                                       combine = TRUE) +
                     ggplot2::labs(title = x, x="UMAP_1", y="UMAP_2") +
                     custom_theme }) %>% cowplot::plot_grid(plotlist=.,
                                                            align="hv",
                                                            axis="tblr",
                                                            nrow=2,
                                                            ncol=3,
                                                            rel_widths=1,
                                                            rel_heights=1,
                                                            greedy=TRUE,
                                                            byrow=TRUE)
    }
    else {
      ggplot() + geom_blank()
      cat("Invalid params:", names(plot.params)[i], ":", plot.params[[i]][1], ":",  plot.params[[i]][2], ":", plot.params[[i]][3], "\n")
    }
    
    # Save the plot
    ggplot2::ggsave(filename = paste0(names(plot.params)[i], assay, ".tiff"),
                    plot = last_plot(),
                    device = "jpeg",
                    path = output_path,
                    scale = 1,
                    width = dplyr::case_when(!is.na(plot.params[[i]][2]) ~ 4*floor(sqrt(length(plot.seurat))),
                                             TRUE ~ 4*3),
                    height = dplyr::case_when(!is.na(plot.params[[i]][2]) ~ 4*ceiling(sqrt(length(plot.seurat))),
                                              TRUE ~ 4*2),
                    units = c("in"),
                    dpi = 600,
                    limitsize = TRUE,
                    bg = "white")
  }
}

calc_scores <- function(integrated.seurat, assay, proj.params, output_path){
  
  set.seed(1234)
  
  # Read marker file
  marker_df <- read.xlsx(xlsxFile = proj.params$cell.type.marker.file)
  
  # Set default assay
  DefaultAssay(integrated.seurat) <- assay
  
  # Iterate through each celltype and plot its module scores
  for (i in 1:ncol(marker_df)){
    
    features <- marker_df[[i]] %>% na.omit()
    features <- rownames(integrated.seurat@assays[[assay]]$data)[tolower(rownames(integrated.seurat@assays[[assay]]$data)) %in% 
                                                                   tolower(features)]
    features <- list(sort(features))
    
    # Calculate module scores
    if (lengths(features) > 1){
      integrated.seurat <- Seurat::AddModuleScore(object=integrated.seurat,
                                                  features=features,
                                                  assay=assay,
                                                  slot="data",
                                                  name=make.names(colnames(marker_df)[i]))
      
      # names(features) <- make.names(colnames(marker_df)[i])
      # integrated.seurat <- UCell::AddModuleScore_UCell(obj=integrated.seurat,
      #                                                  features=features,
      #                                                  assay="RNA",
      #                                                  slot="data",
      #                                                  name="_UCell")
    }
  }
  
  # Create .rds object for integrated seurat object
  if (assay == "RNA"){
    saveRDS(integrated.seurat, file=file.path(output_path, "integrated.seurat.rds"))
  } else{
    # For Spatial.008um and Spatial.016um assays
    saveRDS(integrated.seurat, file=file.path(output_path, paste0("integrated.seurat.", assay, ".rds")))
  }
  
  return(integrated.seurat)
}

plot_scores <- function(integrated.seurat, assay, resolution, reduction, proj.params, output_path){
  
  set.seed(1234)
  
  # Read marker file
  marker_df <- read.xlsx(xlsxFile = proj.params$cell.type.marker.file)
  
  modules <- intersect(colnames(integrated.seurat@meta.data), paste0(make.names(colnames(marker_df)),1))
  
  # Change active.ident
  idents <- paste0("cluster.", resolution, ".", base::tolower(reduction))
  Idents(object = integrated.seurat) <- idents
  
  # Set default assay
  DefaultAssay(integrated.seurat) <- assay
  
  module_score_seurat <- function(module){
    
    Seurat::FeaturePlot(object = integrated.seurat,
                        slot = "data",
                        features = module,
                        reduction = paste0("umap.", base::tolower(reduction)),
                        #cols= c("grey", viridis(n = 10, option="C", direction = 1)),
                        pt.size = 0.2,
                        order = TRUE,
                        label = TRUE,
                        raster = FALSE,
                        combine = TRUE) +  
      #BUG: if raster = TRUE, order = TRUE is ignored. So, set raster = FALSE
      scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name="RdBu"))[5:11])
    # Dont set limits like in plot_features as module scores can be negative too
  }
  
  combined_plot <- purrr::map(.x = modules, 
                              .f = module_score_seurat) %>%         #get(funcs[j]))
    cowplot::plot_grid(plotlist=.,
                       align="hv",
                       axis="tblr",
                       nrow = NULL,
                       ncol = dplyr::if_else(ncol(marker_df) > 10, 5, ceiling(ncol(marker_df)/3)),
                       rel_widths = 1,
                       rel_heights = 1,
                       greedy = TRUE,
                       byrow = TRUE)
  
  ggplot2::ggsave(filename = "Module_plot.jpg",
                  plot = combined_plot,
                  device = "jpeg",
                  path = output_path,
                  width = 8.5*4,
                  height = 11*2,
                  units = c("in"),
                  dpi = 300,
                  limitsize = FALSE,
                  bg = "white")
}

identify_markers_sc_sp <- function(integrated.seurat, assay, resolution, reduction, output_path){
  
  set.seed(1234)
  
  # Set default assay
  DefaultAssay(integrated.seurat) <- assay
  
  # Change active.ident
  idents <- paste0("cluster.", resolution, ".", base::tolower(reduction))
  Idents(object=integrated.seurat) <- idents
  
  # Find ALL markers
  all_markers <- Seurat::FindAllMarkers(object=integrated.seurat,
                                        assay=assay,
                                        logfc.threshold=0.25,
                                        test.use="wilcox",
                                        slot="data",
                                        min.pct=0.1,
                                        min.diff.pct= -Inf,
                                        only.pos=TRUE)
  
  # Get annotations from ENSEMBL
  annotations_list <- get_annotations()
  if (length(intersect(annotations_list[[1]]$SYMBOL, all_markers$gene)) > 
      length(intersect(annotations_list[[2]]$SYMBOL, all_markers$gene))){
    annotations <- annotations_list[[1]]
  } else {
    annotations <- annotations_list[[2]]
  }
  
  # Add gene descriptions
  all_markers <- all_markers %>% 
    dplyr::mutate(pct.1=dplyr::if_else(pct.1 == 0, 0.001, pct.1),
                  pct.2=dplyr::if_else(pct.2 == 0, 0.001, pct.2),
                  ratio=pct.1/pct.2) %>%
    dplyr::left_join(y=annotations, by=c("gene"="ENSEMBL_SYMBOL")) %>%
    dplyr::relocate(cluster, gene, CHR, avg_log2FC, p_val, p_val_adj, pct.1, pct.2, ratio, DESCRIPTION) %>%
    dplyr::distinct(cluster, gene, avg_log2FC, pct.1, pct.2, .keep_all = TRUE)
  
  # Find top markers for each major cluster
  top_markers <- all_markers %>%
    dplyr::filter(avg_log2FC >= 0.58 & p_val_adj < 0.05) %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(desc(avg_log2FC)) %>%  #desc(ratio)
    dplyr::slice_head(n=30) %>%
    ungroup()
  
  # Create a matrix of markers in heatmap format
  mat <- all_markers %>% 
    dplyr::filter(p_val_adj < 0.05) %>%
    tidyr::pivot_wider(id_cols=cluster, names_from = gene, values_from = avg_log2FC, values_fill = 0.0) %>%
    tibble::column_to_rownames("cluster") %>%
    scale()   # column wise scaling so each gene has mean 0, stdev = 1
  
  row_dist <- stats::dist(x=mat, method = "euclidean", diag = TRUE, upper = TRUE)  # distance between rows based on columns
  row_clust <- hclust(row_dist)                                                    # clustering based on distance calculated
  row_order <- rownames(mat[row_clust$order,])
  col_dist <- stats::dist(x=t(mat), method = "euclidean", diag = TRUE, upper = TRUE)
  col_clust <- hclust(col_dist)
  col_order <- colnames(mat[,col_clust$order])
  mat <- mat[row_order, col_order]
  
  mat[mat < 0] <- 0 
  mat <- mat %>% 
    t() %>%
    data.frame() %>% 
    mutate(across(where(is.numeric), ~round(., 2)))
  
  # Save all the markers
  filename <- paste0(proj, ".Markers.All.", idents, ".xlsx")
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb=wb, sheetName="Matrix")
  openxlsx::writeData(wb=wb, sheet="Matrix", x=mat, rowNames = TRUE)
  openxlsx::addWorksheet(wb=wb, sheetName="All_Markers")
  openxlsx::writeData(wb=wb, sheet="All_Markers", x=all_markers)
  openxlsx::addWorksheet(wb=wb, sheetName="Top_Markers")
  openxlsx::writeData(wb=wb, sheet="Top_Markers", x=top_markers)
  openxlsx::saveWorkbook(wb=wb, file=file.path(output_path, filename), overwrite=TRUE)
} 

### Annotate based on clusters variable defined by user
# clusters <- list("Hepatocytes"         = c(),
#                  "Pancreatic.Acinar"   = c(),
#                  "Pancreatic.Islet"    = c(),
#                  "B.Plasma"            = c(),
#                  "T.NK"                = c(),
#                  "Fibroblasts"         = c(),
#                  "Macrophages"         = c(),
#                  "Dendritic"           = c(),
#                  "Endothelial"         = c(),
#                  "Lymph.Endothelial"   = c(),
#                  "Myocytes"            = c(),
#                  "CAFs"                = c(),
#                  "Epithelial"          = c(),
#                  "Neurons"             = c(),
#                  "Epithelial.I"        = c(),
#                  "Epithelial.II"       = c(),
#                  "Epithelial.III"      = c(),
#                  "Epithelial.IV"       = c(),
#                  "Unclassified"        = c())

annotate_manual_sc_sp <- function(integrated.seurat, assay, clusters, resolution, reduction, output_path){
  
  set.seed(1234)
  
  idents <- paste0("cluster.", resolution, ".", base::tolower(reduction))
  
  # Make sure you have assigned all clusters to one of the cell types
  # NOTE: "integrated_snn_res.1.4" etc are stored as factors. 
  # So, use as.character() and then as.numeric() to get accurate cluster values
  list_1 <- integrated.seurat@meta.data %>% 
    dplyr::count(get(idents)) %>% 
    dplyr::select(identity(1)) %>% 
    unlist(use.names=FALSE) %>% 
    as.character() %>% 
    as.numeric() %>% 
    sort()
  
  list_2 <- clusters %>% 
    unlist(., use.names=FALSE) %>% 
    sort()
  
  # Proceed with annotation ONLY if all clusters have been annotated
  if (identical(list_1, list_2)){
    print("All Clusters have been annotated")
    
    # Extract metadata from Seurat object, assign appropriate resolution to
    # seurat_clusters column and add Cell.Type, Cell.Subtype columns
    data <- integrated.seurat@meta.data %>% 
      dplyr::mutate(seurat_clusters=get(idents),
                    Cell.Type=NA, Cell.Subtype=NA)
    
    # Assign cell type based on cluster numbers within seurat_clusters column
    for (j in 1:nrow(data)){
      for (i in 1:length(clusters)){
        if (as.numeric(as.character(data$seurat_clusters[j])) %in% clusters[[i]]){
          data$Cell.Type[j] <- names(clusters[i])
        }
      }
    }
    
    # Check summary of cell counts
    print(data %>% dplyr::count(Cell.Type) %>% dplyr::arrange(n))
    cat("\n")
    
    # Import the metadata into Seurat object
    integrated.seurat@meta.data <- data
    
    # Create .rds object for integrated seurat object
    if (assay == "RNA"){
      saveRDS(integrated.seurat, file=file.path(output_path, "integrated.seurat.ann.rds"))
    } else{
      # For Spatial.008um and Spatial.016um assays
      saveRDS(integrated.seurat, file=file.path(output_path, paste0("integrated.seurat.", assay, ".ann.rds")))
    }
    
  } else {
    cat("\nYou missed annotating these clusters:\t", setdiff(list_1, list_2))
    cat("\nThese clusters are not present in data:\t", setdiff(list_2, list_1))
    cat("\nThese clusters have duplicate annotation:\t", list_2[duplicated(list_2)])
  }
  
  return(integrated.seurat)
} 

estimate_cell_type <- function(integrated.seurat){
  
  set.seed(1234)
  
  # Get columns containing module scores fro each cell type
  score_cols <- integrated.seurat@meta.data %>% 
    dplyr::select(ends_with("1")) %>% 
    colnames()
  
  # Extract numeric score matrix
  score_mat <- as.matrix(integrated.seurat@meta.data[, score_cols])
  
  # Get index of max per row
  max_idx <- max.col(score_mat, ties.method = "first")
  
  # Map to column names
  integrated.seurat@meta.data$Est.Cell.Type <- colnames(score_mat)[max_idx]
  
  cleaned_score_cols <- gsub("1$", "", score_cols)
  cell.types <- unique(integrated.seurat@meta.data$Cell.Type)
  
  # Key:value dataframe
  map_df <- list("Epithelial1" = c("Epithelial.I", "Epithelial.II"),
                 "Macrophages1" = c("Macrophages"),
                 "Myofibroblasts1" = c("Myocytes.Myofibroblasts"),
                 "Endothelial1" = c("Endothelial"),
                 "Lymphatic.Endothelial1" = c("Endothelial"),
                 "T1" = c("T.NK"), 
                 "NK1" = c("T.NK"),
                 "Fibroblasts1" = c("Fibroblasts"),
                 "B1" = c("B.Plasma"),
                 "Plasma1" = c("B.Plasma"),
                 "Granulocytes1" = c("Granulocytes", "Mast"),
                 "Erythrocytes1" = c("Erythrocytes"),
                 "Neurons1" = c("Neurons"),
                 "Dendritic1" = c("Dendritic"),
                 "Mast1" = c("Mast"))
  
  
  integrated.seurat@meta.data <- integrated.seurat@meta.data %>%
    dplyr::mutate(Confidence = ifelse(mapply(function(est, true){
      true %in% map_df[[est]]}, Est.Cell.Type, Cell.Type), "HIGH", "LOW"))
  
  
  return(integrated.seurat) 
}

plot_dot_plot <- function(integrated.seurat, idents, features, filename, output_path, gene_y_axis = FALSE, split.col = NULL){
  
  set.seed(1234)
  
  # Re-order the active ident alphabetically
  Idents(integrated.seurat) <- idents
  Idents(integrated.seurat) <- base::factor(x = integrated.seurat@active.ident, 
                                            levels= sort(levels(integrated.seurat@active.ident)))
  
  # Determine grouping variable
  if (is.null(split.col)) {
    groups <- "All"
  } else {
    groups <- unique(integrated.seurat@meta.data[[split.col]])
  }
  
  # Loop over groups
  plot_list <- list()
  for (g in groups) {
    if (g == "All") {
      subset_obj <- integrated.seurat
      plot_title <- "All"
    } else {
      subset_obj <- integrated.seurat[, integrated.seurat@meta.data[[split.col]] == g]
      plot_title <- g
    }
    
    p <- Seurat::DotPlot(object = subset_obj,
                         assay = "RNA",
                         features = features,
                         dot.min = 0,
                         dot.scale = 4,
                         scale = TRUE,         # if plotting 1+ gene, always SCALE
                         scale.by = "size",    # scale percent cells expressed i.e. size(area) of dots
                         scale.min = 0,        # set this to 0, so pct.size = 0 has smallest size
                         scale.max = NA) +     # if you set this to 100, points become too small
      ggplot2::geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.25) +
      ggplot2::labs(fill = "Avg. Expression",         # rename color legend
                    size  = "% Cells Expressed",
                    title = plot_title) +    # rename size legend   
      ggplot2::scale_colour_distiller(type = "div", palette = "RdYlGn", direction = -1) +
      ggplot2::guides(size = ggplot2::guide_legend(
        override.aes = list(shape = 21, colour = "black", fill = "white", stroke = 0.75))) +
      custom_theme +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
    
    if (gene_y_axis){
      p <- p + ggplot2::coord_flip()
    }
    
    plot_list[[plot_title]] <- p
  }
  
  # Per panel dimensions
  n_y <- length(unique(integrated.seurat@active.ident))  
  height_per_element <- 0.5
  width_per_element <- 0.2 
  if (gene_y_axis){
    height_panel <- length(features) * width_per_element + 2
    width_panel <- n_y * height_per_element + 6
  } else{
    height_panel <- n_y * height_per_element + 2              # extra for x-axis labels
    width_panel <- length(features) * width_per_element + 6 #extra for y-axis labels, legend
  }
  
  # Automatically calculate ncol and nrow for cowplot
  n_plots <- length(plot_list)
  ncol_plots <- ceiling(sqrt(n_plots))
  nrow_plots <- ceiling(n_plots / ncol_plots)
  
  # Total height
  total_height <- height_panel * nrow_plots
  total_width  <- width_panel * ncol_plots
  
  
  
  # Combine all plots
  combined_plot <- cowplot::plot_grid(plotlist = plot_list, 
                                      ncol = ncol_plots, 
                                      nrow = nrow_plots)
  
  # Save combined plot
  ggsave(filename = file.path(output_path, paste0(filename,".pdf")),
         plot = combined_plot,
         width = total_width,
         height = total_height,
         limitsize = FALSE,
         units = "in",
         bg        = "white")
}

find_top_features <- function(integrated.seurat, assay, proj.params){
  
  set.seed(1234)
  
  marker_df <- openxlsx::read.xlsx(file.path(proj.params$cell.type.marker.file))
  expr_data <- Seurat::GetAssayData(object = integrated.seurat, assay = assay, layer = "data")
  
  top_features <- lapply(marker_df, FUN = function(f){
    
    # Subset features and expression data
    f <- f[base::tolower(f) %in% base::tolower(SeuratObject::Features(integrated.seurat))] %>% 
      na.omit() %>%
      as.vector()
    
    if(length(f) == 0) return(NULL)
    
    expr_sub <- expr_data[f, , drop = FALSE]
    
    # Find top 3 expressing markers
    rowMeans(expr_sub, na.rm = TRUE) %>% 
      sort(., decreasing = TRUE) %>% 
      head(3) %>% 
      names()
  })
  
  # Keep only cell types present in your dataset
  celltypes_in_seurat <- unique(integrated.seurat$Cell.Type)
  top_features <- top_features[intersect(names(top_features), celltypes_in_seurat)]
  
  return(unlist(top_features, use.names = FALSE))
}

plot_dot_custom <- function(integrated.seurat, assay, proj.params, ident.1, ident.2, features, filename, output_path, split.col=NULL){
  
  set.seed(1234)
  
  # Keep only cell types present in your dataset and marker file
  marker_df <- openxlsx::read.xlsx(file.path(proj.params$cell.type.marker.file))
  celltypes_in_seurat <- unique(integrated.seurat$Cell.Type)
  celltypes_to_keep <- intersect(colnames(marker_df), celltypes_in_seurat)
  
  # Determine groups based on split.col
  groups <- if (!is.null(split.col) && split.col %in% colnames(integrated.seurat@meta.data)) {
    as.character(unique(integrated.seurat@meta.data[[split.col]]))
  } else {
    "All"
  }
  
  # Determine features present in data set
  features <- features[base::tolower(features) %in% base::tolower(SeuratObject::Features(integrated.seurat))] %>% 
    na.omit() %>%
    as.vector()
  
  if(length(features) == 0) stop("No matching features found in the assay.")
  
  # Loop over groups
  plot_list <- list()
  for (g in groups) {
    
    # Subset object
    if (g == "All") {
      subset_obj <- integrated.seurat
      plot_title <- ""
    } else {
      subset_obj <- integrated.seurat[, integrated.seurat@meta.data[[split.col]] == g]
      plot_title <- g
    }
    
    # Determine levels based on ident.1
    levels <- if (!is.null(ident.1) && ident.1 %in% colnames(integrated.seurat@meta.data)) {
      subset_obj@meta.data[[ident.1]] %>% unique() %>% as.character()
    } else {
      stop(ident.1, " missing in metadata")
    }
    
    # Determine levels.extra based on ident.2
    levels.extra <- if (!is.null(ident.2) && ident.2 %in% colnames(integrated.seurat@meta.data)) {
      subset_obj@meta.data[[ident.2]] %>% unique() %>% as.character()
    } else {
      levels.extra <- NA # set NA so loop runs once later below
    }
    
    # Limit to one feature if ident.2 is used
    if (length(features) > 1 & !is.null(ident.2)){
      warning("You cannot plot more than 1 feature using variables ", ident.1,  " and ", ident.2, ".\n")
      warning("Proceeding to generate dot plot using only one feature ", features[[1]], ".")
      features <- features[1]
    }
    
    # Prepare dot plot data
    pct.exp <- list()
    avg.exp <- list()
    features.plot <- list()
    id <- list()
    id.extra <- list()
    
    for (le in levels.extra){
      for (l in levels){
        for (f in features){
          
          # Cells to retain
          cells <- subset_obj@meta.data %>% 
            dplyr::filter(.data[[ident.1]] == l) 
          if (!is.null(ident.2)){
            cells <- cells %>% 
              dplyr::filter(.data[[ident.2]] == le) 
          }
          cells <- rownames(cells)
          
          # Subset features and expression data
          expr_data <- Seurat::GetAssayData(object = subset_obj, assay = assay, layer = "data")
          subset_data <- expr_data[f, cells, drop = FALSE]
          pct <- sum(subset_data > 0)/length(subset_data) * 100
          avg <- mean(expm1(subset_data), na.rm = TRUE)   # back-transform from log1p
          
          pct.exp <- c(pct.exp, pct)
          avg.exp <- c(avg.exp, avg)
          features.plot <- c(features.plot, f)
          id <- c(id, l)
          id.extra <- c(id.extra, le)
        }
      }
    }
    
    # Combine and scale data
    dotplot_data <- data.frame(avg.exp = unlist(avg.exp),
                               pct.exp = unlist(pct.exp),
                               features.plot = unlist(features.plot),
                               id = unlist(id),
                               id.extra = unlist(id.extra)) %>%
      dplyr::group_by(features.plot, id.extra) %>%
      dplyr::mutate(avg.exp.scaled = scale(avg.exp)) %>%
      dplyr::mutate(avg.exp.scaled = dplyr::case_when(avg.exp.scaled > 2.5 ~ 2.5,
                                                      avg.exp.scaled < -2.5 ~ -2.5,
                                                      TRUE ~ avg.exp.scaled)) %>%
      as.data.frame() %>%
      dplyr::mutate(features.plot = factor(features.plot, levels = features))
    
    if (ident.1 == "Cell.Type"){
      dotplot_data <- dotplot_data %>%
        filter(id %in% celltypes_to_keep)
    }
    
    # Determine X and Y variables
    if (is.null(ident.2) && length(features) > length(levels) & length(features) > length(levels.extra)){
      y_var <- "features.plot"
      x_var <- "id"
    } else if (!is.null(ident.2) && length(levels.extra) > length(levels)){
      y_var <- "id.extra"
      x_var <- "id"
    } else if (!is.null(ident.2) && length(levels.extra) < length(levels)){
      y_var <- "id"
      x_var <- "id.extra"
    } else{
      y_var <- "id"
      x_var <- "features.plot"
    }
    
    # Create ggplot
    p <- ggplot(data = dotplot_data, 
                mapping = aes(x = !!sym(x_var), y = !!sym(y_var), 
                              size = pct.exp, fill = avg.exp.scaled)) +
      geom_point(shape = 21, colour = "black", stroke = 0.25) +
      labs(x = "", y = "", title = plot_title,
           fill = "Avg. Expression",
           size = "% Cells Expressed") +
      scale_fill_distiller(type = "div", palette = "RdYlGn", direction = -1) +
      guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white", stroke = 0.75))) +
      scale_size_area(max_size = 15) + #min_size = 1) +
      custom_theme +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    
    plot_list[[plot_title]] <- p
  }
  
  # Calculate panel dimensions
  height_per_element <- 0.5
  width_per_element <- 0.2
  height_panel <- length(unique(dotplot_data[[y_var]])) * height_per_element + 2
  width_panel <- length(unique(dotplot_data[[x_var]])) * width_per_element + 6
  
  # Determine layout
  n_plots <- length(plot_list)
  ncol_plots <- ceiling(sqrt(n_plots))
  nrow_plots <- ceiling(n_plots / ncol_plots)
  
  # Total dimensions
  total_height <- height_panel * nrow_plots
  total_width  <- width_panel * ncol_plots
  
  # Combine plots
  combined_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = ncol_plots, nrow = nrow_plots)  
  
  # Save plot
  filename <- gsub("[/\\?<>\\:*|\"]", "_", filename)
  ggsave(filename = file.path(output_path, paste0(filename,".pdf")),
         plot = combined_plot,
         width = total_width,
         height = total_height,
         limitsize = FALSE,
         units = "in",
         bg = "white")
}

plot_umap <- function(integrated.seurat, reduction, color.col, filename, output_path, split.col = NULL){
  
  set.seed(1234)
  
  # Determine levels of color.col
  all_levels <- unique(integrated.seurat@meta.data[[color.col]])
  
  # Determine groups
  if (is.null(split.col)) {
    groups <- "All"
  } else {
    groups <- unique(integrated.seurat@meta.data[[split.col]])
  }
  
  # Loop over groups
  all_plots <- list()
  for (g in groups) {
    if (g == "All") {
      subset_obj <- integrated.seurat
      group_label <- ""
    } else {
      subset_obj <- integrated.seurat[, integrated.seurat@meta.data[[split.col]] == g]
      group_label <- g
    }
    
    # Fix factor levels for consistent coloring, in case a cell type is absent in one of the splits 
    subset_obj@meta.data[[color.col]] <- factor(subset_obj@meta.data[[color.col]], levels = all_levels)
    
    p <- Seurat::DimPlot(object = subset_obj,
                         reduction = reduction, #paste0("umap.", tolower(reduction)),
                         cols      = custom_palette,
                         label     = FALSE,
                         group.by  = color.col,
                         pt.size   = 0.2,
                         label.size= 5,
                         repel     = FALSE,
                         raster    = FALSE) +
      ggplot2::labs(color = "CLUSTERS", x = "UMAP_1", y = "UMAP_2", title = group_label) +
      custom_theme +
      ggplot2::coord_fixed(ratio = 1)  # 1 unit on x-axis = 1 unit on y-axis, ensuring a square aspect
    
    all_plots[[group_label]] <- p
  }
  
  # Automatically calculate ncol and nrow for cowplot
  n_plots <- length(all_plots)
  ncol_plots <- ceiling(sqrt(n_plots))
  nrow_plots <- ceiling(n_plots / ncol_plots)
  
  # Combine all plots
  combined_plot <- cowplot::plot_grid(plotlist = all_plots, ncol = ncol_plots, nrow = nrow_plots)
  
  # Save
  ggplot2::ggsave(filename  = file.path(output_path, paste0(filename, ".pdf")),
                  plot      = combined_plot,
                  width     = ncol_plots * 10, # 2 inch for legend
                  height    = nrow_plots * 8, 
                  units     = "in",
                  limitsize = FALSE,
                  bg        = "white")
}

plot_features <- function(integrated.seurat, features, reduction, filename, output_path, split.col=NULL){
  
  set.seed(1234)
  
  # Input checks
  existing_features <- intersect(features, rownames(integrated.seurat@assays$RNA))
  missing_features <- setdiff(features, existing_features)
  
  if (length(existing_features) == 0) {
    stop("None of the provided features are found in the RNA assay. Nothing to plot.")
  }
  
  if (length(missing_features) > 0) {
    warning("The following feature(s) were not found and will be skipped: ", paste(missing_features, collapse = ", "))
  }
  
  if (!is.null(split.col) && !split.col %in% colnames(integrated.seurat@meta.data)) {
    stop(paste("Column", split.col, "not found in metadata."))
  }
  
  # Determine groups
  if (is.null(split.col)) {
    groups <- "All"
  } else {
    groups <- unique(integrated.seurat@meta.data[[split.col]])
  }
  
  # Loop over groups and features
  all_plots <- list()
  for (feature in existing_features) { 
    for (g in groups) {
      if (g == "All") {
        subset_obj <- integrated.seurat
        group_label <- ""
      } else {
        subset_obj <- integrated.seurat[, integrated.seurat@meta.data[[split.col]] == g]
        group_label <- paste0(" (", g, ")")
      }
      
      max_expr <- max(subset_obj@assays$RNA@data[feature, ], na.rm = TRUE)
      p <- Seurat::FeaturePlot(object = subset_obj,
                               slot = "data",
                               features = feature,
                               reduction = reduction, #paste0("umap.", tolower(reduction)),
                               #cols = c("grey", viridis(n = 10, option = "C", direction = -1)),
                               pt.size = 0.2,
                               label = FALSE,
                               order = TRUE,
                               raster = FALSE,
                               combine = TRUE) +
        ggplot2::labs(title = paste0(feature, group_label), x="UMAP_1", y="UMAP_2") +
        custom_theme +
        ggplot2::scale_colour_gradientn(
          colours = rev(RColorBrewer::brewer.pal(11, "RdBu"))[5:11],
          limits = c(0, max(max_expr, 1e-2)))
      
      all_plots[[paste0(feature, group_label)]] <- p
    }
  }
  
  # Automatically calculate ncol and nrow for cowplot
  n_plots <- length(all_plots)
  ncol_plots <- ceiling(sqrt(n_plots))
  nrow_plots <- ceiling(n_plots / ncol_plots)
  
  # Restrict unsupported combination
  if (ncol_plots > 10 && nrow_plots > 10) {
    stop("Image size too large. More than 100 plots cannot be viewed in single figure")
  }
  
  # Combine all plots
  combined_plot <- cowplot::plot_grid(plotlist = all_plots, ncol = ncol_plots, nrow = nrow_plots)
  
  # Save
  ggplot2::ggsave(filename  = file.path(output_path, paste0(filename, ".jpg")),
                  plot      = combined_plot,
                  width     = ncol_plots * 6,
                  height    = nrow_plots * 6,
                  units     = "in",
                  limitsize = FALSE,
                  bg        = "white")
}

# Input is seurat object of a single slide with columns X, Y, Sample, Group
plot_spatial_map <- function(plot.seurat, x1, y1, x2, y2, suffix, output_path){
  
  set.seed(1234)
  
  Sample <- unique(plot.seurat@meta.data$Sample)
  
  # Get an idea of co-ordinates
  xmin <- (min(plot.seurat@meta.data$X)%/%100-1)*100
  xmax <- (max(plot.seurat@meta.data$X)%/%100+1)*100
  ymin <- (min(plot.seurat@meta.data$Y)%/%100-1)*100
  ymax <- (max(plot.seurat@meta.data$Y)%/%100+1)*100
  
  if (Sample == "TMA1-A1"){
    # Flip on X axis to get correct orientation for TMA1-A1
    p1 <- ggplot(plot.seurat@meta.data, mapping=aes(x=X, y=Y, color = Sample)) +
      geom_point() +
      scale_x_reverse(limits =c(xmax, xmin), breaks=seq(xmax, xmin, -100), position="top") +
      scale_y_continuous(limits =c(ymin, ymax), breaks=seq(ymin, ymax, 100), position="left") +
      geom_vline(xintercept = x1) +
      geom_hline(yintercept = y1) +
      scale_color_manual(values=custom_palette)
    
    p2 <- ggplot(plot.seurat@meta.data, mapping=aes(x=X, y=Y, color = Treatment)) +
      geom_point() +
      scale_x_reverse(limits =c(xmax, xmin), breaks=seq(xmax, xmin, -100), position="top") +
      scale_y_continuous(limits =c(ymin, ymax), breaks=seq(ymin, ymax, 100), position="left") +
      geom_vline(xintercept = x1) +
      geom_hline(yintercept = y1) +
      scale_color_manual(values=custom_palette)
  } else if (Sample == "TMA1-D1"){
    # Flip on Y axis to get correct orientation for TMA1-D1
    p1 <- ggplot(plot.seurat@meta.data, mapping=aes(x=X, y=Y, color = Sample)) +
      geom_point() +
      scale_x_continuous(limits =c(xmin, xmax), breaks=seq(xmin, xmax, 100), position="top") +
      scale_y_reverse(limits =c(ymax, ymin), breaks=seq(ymax, ymin, -100), position="left") +
      geom_vline(xintercept = x2) +
      geom_hline(yintercept = y2) +
      scale_color_manual(values=custom_palette)
    
    p2 <- ggplot(plot.seurat@meta.data, mapping=aes(x=X, y=Y, color = Treatment)) +
      geom_point() +
      scale_x_continuous(limits =c(xmin, xmax), breaks=seq(xmin, xmax, 100), position="top") +
      scale_y_reverse(limits =c(ymax, ymin), breaks=seq(ymax, ymin, -100), position="left") +
      geom_vline(xintercept = x2) +
      geom_hline(yintercept = y2) +
      scale_color_manual(values=custom_palette)
  }
  
  # Save the plot
  ggplot2::ggsave(filename=paste0("Sample.Map.", Sample, ".", suffix, ".tiff"),
                  plot=p1+p2,
                  device="jpeg",
                  path=output_path,
                  width=25,
                  height=8.5,
                  units=c("in"),
                  dpi=600,
                  limitsize=TRUE,
                  bg="white")
}  

tabulate_frequency <- function(integrated.seurat, split.cols, output_path){
  
  set.seed(1234)
  
  # Create a new workbook
  wb <- createWorkbook()
  
  for (split.col in split.cols){
    if(!is.null(integrated.seurat@meta.data[[split.col]] %>% unique())){
      
      counts <- integrated.seurat@meta.data %>% 
        dplyr::count(.data[[split.col]], Cell.Type) %>% 
        dplyr::group_by(.data[[split.col]]) %>% 
        dplyr::mutate(Percent = round(100 * n / sum(n), 2)) %>%
        tidyr::pivot_wider(names_from = .data[[split.col]], values_from = c(n, Percent)) %>%
        as.data.frame()
      
      
      # Make sheet name safe
      sheet_name <- substr(gsub("[/\\?<>\\:*|\"]", "_", split.col), 1, 31)  # Excel sheet name limit = 31
      addWorksheet(wb, sheetName = sheet_name)
      writeData(wb, sheet = sheet_name, counts)
    }
  }
  
  # Save workbook
  saveWorkbook(wb, file = file.path(output_path, "Population_Frequencies.xlsx"), overwrite = TRUE)
}

prep_pseudobulk <- function(integrated.seurat, comparison.col, assay = "RNA"){
  
  set.seed(1234)
  
  # Input checks
  if (!inherits(integrated.seurat, "Seurat")) {
    stop("'integrated.seurat' must be a Seurat object.")
  }
  if (!comparison.col %in% colnames(integrated.seurat@meta.data)) {
    stop(paste0("'", comparison.col, "' not found in metadata columns."))
  }
  if (!"Sample" %in% colnames(integrated.seurat@meta.data)) {
    stop("'Sample' column is required in metadata.")
  }
  if (!assay %in% names(integrated.seurat@assays)) {
    stop(paste0("Assay '", assay, "' not found in Seurat object."))
  }
  
  # Initialize storage
  counts_matrix <- integrated.seurat@assays[[assay]]$counts
  meta_data_full <- data.frame()
  read_data_full <- data.frame(SYMBOL = rownames(counts_matrix))
  
  # Extract unique groups from comparison.col that will be compared in DE analysis
  groups <- unique(integrated.seurat@meta.data[[comparison.col]])
  groups <- groups[!is.na(groups)]  # remove NA groups if present
  
  # Loop through each group
  for (g in groups) {
    subset.seurat <- integrated.seurat[, integrated.seurat@meta.data[[comparison.col]] == g]
    
    # Generate metadata
    meta_data <- subset.seurat@meta.data %>%
      tibble::rownames_to_column("barcode") %>%
      dplyr::filter(!is.na(Sample)) %>%
      dplyr::mutate(Sample_ID = paste0(Sample, "_", g),
                    Comparisons = .data[[comparison.col]]) %>%
      dplyr::add_count(Sample_ID) %>%
      dplyr::filter(n >= 100)
    
    ### Generate read data
    # read data will have "the reads of all cells belonging to a single sample" 
    # merged together in each column. 
    
    # First, create a list of samples
    samples <-  meta_data[["Sample_ID"]] %>%
      unique()
    
    # Second, initialize read counts matrix
    read_data <- matrix(0,
                        nrow = nrow(counts_matrix),
                        ncol = length(samples),
                        dimnames = list(rownames(counts_matrix), samples))
    
    # Third, add row-wise, the counts of each gene for each sample
    for(i in samples){
      
      # Create a list of cells for each sample
      cells_subset <- meta_data %>% 
        dplyr::filter(Sample_ID == i) %>% 
        dplyr::pull(barcode)
      
      # Subset counts
      subset_counts <- counts_matrix[, cells_subset, drop=FALSE]
      
      # Sum counts
      read_data[, i] <- Matrix::rowSums(subset_counts)
    }
    
    # Fourth, format metadata and readdata
    cols_to_keep <- intersect(c("Sample_ID", "Comparisons", "Patient", "Condition", "Treatment", "Disease", 
                                "Tissue", "Strain", "Cell_line", "Sex", "n", "barcode"),
                              colnames(meta_data))
    meta_data <- meta_data %>%
      dplyr::distinct(Sample_ID, .keep_all = TRUE) %>%
      dplyr::select(all_of(cols_to_keep))
    
    read_data <- as.data.frame(read_data) %>% 
      tibble::rownames_to_column("SYMBOL")
    
    # Finally, merge results
    meta_data_full <- dplyr::bind_rows(meta_data_full, meta_data)
    read_data_full <- dplyr::left_join(read_data_full, read_data, by=c("SYMBOL"="SYMBOL"))
  }
  
  return(list(meta_data = meta_data_full, read_data = read_data_full))
}

find_degs_seurat <- function(integrated.seurat, comparison.col, output_path, celltype.col = "Cell.Type", assay = "RNA"){
  
  set.seed(1234)
  
  # Setup log file
  log_file <- file.path(output_path, "find_degs_seurat.log")
  log_con <- file(log_file, open = "wt")
  log_msg <- function(...) {
    msg <- paste0(...)
    message(msg)                 # print to console
    writeLines(msg, log_con)     # also save to log file
  }
  on.exit(close(log_con), add = TRUE)  # close file safely when function ends
  
  # Initialize storage
  deg_MAST_list <- list()
  deg_Wilcox_list <- list()
  deg_DESeq2_list <- list()
  counts_list <- list()
  vst_counts_list <- list()
  skipped <- c()
  
  # Define cell types and comparisons
  celltype_levels <- unique(integrated.seurat@meta.data[[celltype.col]])
  comparison_levels <- unique(integrated.seurat@meta.data[[comparison.col]])
  comparisons <- utils::combn(x = comparison_levels, m = 2, simplify = FALSE)
  
  for(celltype in celltype_levels){
    for(comparison in comparisons){
      
      target <- comparison[1]
      reference <- comparison[2]
      
      # Subset Seurat object
      subset_obj <- integrated.seurat[, integrated.seurat@meta.data[[celltype.col]] == celltype]
      subset_obj <- subset_obj[, subset_obj@meta.data[[comparison.col]] %in% c(target, reference)]
      
      # Count cells
      n_cells_ref <- sum(subset_obj@meta.data[[comparison.col]] == reference)
      n_cells_target <- sum(subset_obj@meta.data[[comparison.col]] == target)
      
      # Count samples with >= 100 cells
      n_samples_ref <- subset_obj@meta.data %>%
        dplyr::filter(.data[[comparison.col]] == reference) %>%
        dplyr::count(Sample) %>% dplyr::filter(n >= 100) %>% nrow()
      n_samples_target <- subset_obj@meta.data %>%
        dplyr::filter(.data[[comparison.col]] == target) %>%
        dplyr::count(Sample) %>% dplyr::filter(n >= 100) %>% nrow()
      
      # Skip if too few cells
      if(n_cells_ref < 100 | n_cells_target < 100){
        log_msg("Skipping subset ", celltype, " | ",
                target, " (", n_cells_target, " cells) vs ",
                reference, " (", n_cells_ref, " cells) due to <100 cells")
        skipped <- c(skipped, paste(celltype, target, reference, sep="_"))
        next
      }
      
      # ---- MAST (always run) ----
      deg_obj <- subset_obj
      log_msg("Running MAST for ", celltype, " | ",
              target, " (", n_cells_target, " cells) vs ",
              reference, " (", n_cells_ref, " cells)")
      
      mast_degs <- FindMarkers(
        object = deg_obj,
        ident.1 = target,
        ident.2 = reference,
        group.by = comparison.col,
        assay = assay,
        test.use = "MAST",
        slot = "data",         # MAST must be run on log transformed data
        min.cells.group = 2,
        min.pct = 0.1
      )
      mast_degs <- mast_degs %>%
        tibble::rownames_to_column("SYMBOL") %>%
        dplyr::mutate(Comparison = paste0(celltype, ".", target, ".vs.", reference))
      deg_MAST_list[[paste0(celltype, ".", target, ".vs.", reference)]] <- mast_degs
      
      # ---- Wilcoxon (always run) ----
      deg_obj <- subset_obj
      log_msg("Running Wilcoxon for ", celltype, " | ",
              target, " (", n_cells_target, " cells) vs ",
              reference, " (", n_cells_ref, " cells)")
      wilcox_degs <- FindMarkers(
        object = deg_obj,
        ident.1 = target,
        ident.2 = reference,
        group.by = comparison.col,
        assay = assay,
        test.use = "wilcox",
        slot = "data",
        min.cells.group = 2,
        min.pct = 0.1
      )
      wilcox_degs <- wilcox_degs %>%
        tibble::rownames_to_column("SYMBOL") %>%
        dplyr::mutate(Comparison = paste0(celltype, ".", target, ".vs.", reference))
      deg_Wilcox_list[[paste0(celltype, ".", target, ".vs.", reference)]] <- wilcox_degs
      
      
      # ---- Counts (always run) ----
      deg_obj <- Seurat::AggregateExpression(
        object = subset_obj,
        group.by = c("Sample", comparison.col),
        assays = assay,
        slot = "counts",
        return.seurat = TRUE
      )
      
      raw_counts <- as.data.frame(as.matrix(GetAssayData(deg_obj, assay = assay, slot = "counts"))) 
      counts_list[[paste0(celltype, ".", target, ".vs.", reference)]] <- raw_counts %>%
        tibble::rownames_to_column(var = "SYMBOL")
      
      # VST normalized counts
      meta_data <- deg_obj@meta.data[colnames(raw_counts), , drop = FALSE]
      rownames(meta_data) <- colnames(raw_counts)
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = raw_counts,
                                            colData = meta_data,
                                            design = ~1)
      vst_counts <- tryCatch(
        {
          vst <- DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
          as.data.frame(SummarizedExperiment::assay(vst)) %>%
            rownames_to_column(var = "SYMBOL")
        },
        error = function(e) {
          message("VST failed due to zeros in counts. Returning empty data frame instead.")
          data.frame(SYMBOL = character(0))  # ensures at least SYMBOL column exists
        }
      )
      vst_counts_list[[paste0(celltype, ".", target, ".vs.", reference)]] <- vst_counts
      
      # ---- DESeq2 pseudobulk ----
      if(n_samples_ref >= 2 & n_samples_target >= 2){
        log_msg("Running DESeq2 for ", celltype, " | ",
                target, " (", n_cells_target, " cells) vs ",
                reference, " (", n_cells_ref, " cells)")
        
        deseq2_degs <- FindMarkers(
          object = deg_obj,
          ident.1 = target,
          ident.2 = reference,
          group.by = comparison.col,
          assay = assay,
          test.use = "DESeq2",
          slot = "counts",      # DESeq2 must be run on raw counts
          min.cells.group = 2   # DESeq2 uses this argument to determine number of groups
        )
        deseq2_degs <- deseq2_degs %>%
          tibble::rownames_to_column("SYMBOL") %>%
          dplyr::mutate(Comparison = paste0(celltype, ".", target, ".vs.", reference))
        deg_DESeq2_list[[paste0(celltype, ".", target, ".vs.", reference)]] <- deseq2_degs
      }
    }
  }
  
  # ---- Save to separate Excel files ----
  save_list_to_xlsx <- function(lst, file_path) {
    wb <- openxlsx::createWorkbook()
    
    # Create index sheet first
    index_df <- data.frame(
      Sheet = paste0("Sheet", seq_along(lst)),
      Comparison = names(lst),
      stringsAsFactors = FALSE
    )
    openxlsx::addWorksheet(wb, "Index")
    openxlsx::writeData(wb, "Index", index_df)
    
    # Add each comparison as Sheet1, Sheet2, ...
    for (i in seq_along(lst)) {
      sheet_nm <- paste0("Sheet", i)
      openxlsx::addWorksheet(wb, sheetName = sheet_nm)
      openxlsx::writeData(wb, sheet = sheet_nm, lst[[i]])
    }
    openxlsx::saveWorkbook(wb, file_path, overwrite = TRUE)
  }
  
  all_DESeq2 <- dplyr::bind_rows(deg_DESeq2_list)
  all_MAST <- dplyr::bind_rows(deg_MAST_list)
  all_Wilcox <- dplyr::bind_rows(deg_Wilcox_list)
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "DEGs_DESeq2")
  openxlsx::writeData(wb, "DEGs_DESeq2", all_DESeq2)
  openxlsx::addWorksheet(wb, "DEGs_MAST")
  openxlsx::writeData(wb, "DEGs_MAST", all_MAST)
  openxlsx::addWorksheet(wb, "DEGs_Wilcox")
  openxlsx::writeData(wb, "DEGs_Wilcox", all_Wilcox)
  openxlsx::saveWorkbook(wb, file.path(output_path, "Seurat_DEGs.xlsx"), overwrite = TRUE)
  
  
  if (length(counts_list) > 0) save_list_to_xlsx(counts_list, file.path(output_path, "Seurat_Raw_Counts.xlsx"))
  if (length(vst_counts_list) > 0) save_list_to_xlsx(vst_counts_list, file.path(output_path, "Seurat_VST_Counts.xlsx"))
  
  # ---- Report skipped ----
  if (length(skipped) > 0) {
    log_msg("Skipped comparisons (low cell counts):")
    log_msg(paste(" -", skipped, collapse = "\n"))
  }
  
  # ---- Return ----
  return(list(
    DESeq2 = deg_DESeq2_list,
    MAST   = deg_MAST_list,
    Wilcox = deg_Wilcox_list
  ))
}

# We need to run SCT on each assay. But for cell cycle scoring, NA is present
# in barcodes from other assays like 002um and 016um. So, it throws error.
# Therefore, create a separate seurat object for each bin size so that each
# seurat object has only one assay equivalent to "RNA" assay of single cell data


# ---- SURVIVAL FUNCTIONS ----

# Display pre-defined survival analysis scenarios
show_survival_scenarios <- function() {
  scenarios <- tibble::tribble(
    ~Scenario, ~Description, ~stratify_var, ~substratify_var, ~facet_var, ~multiple_cutoff, ~Curves_per_plot, ~Plots_facets, ~Curve_labels,
    "(i)", "Survival based on gene A", "Gene A", "–", "–", "–", 2, 1, "HIGH vs LOW",
    "(ii)", "Survival based on gene A + Sex", "Gene A", "Sex", "–", "FALSE / TRUE", 4, 1, "HIGH/LOW × Male/Female",
    "(iii)", "Survival based on gene A, faceted by Sex", "Gene A", "–", "Sex", "FALSE / TRUE", 2, 2, "HIGH vs LOW",
    "(iv)", "Survival based on gene A + Smoking, faceted by Sex", "Gene A", "Smoking", "Sex", "FALSE / TRUE", 4, 2, "HIGH/LOW × Yes/No",
    "(v)", "Survival based on Sex", "Sex", "–", "–", "–", 2, 1, "Male vs Female",
    "(vi)", "Survival based on Sex + Race", "Sex", "Race", "–", "–", 4, 1, "Male/Female × White/Black",
    "(vii)", "Survival based on Sex, faceted by Race", "Sex", "–", "Race", "–", 2, 2, "Male vs Female",
    "(viii)", "Survival based on Sex + Smoking, faceted by Race", "Sex", "Smoking", "Race", "FALSE / TRUE", 4, 2, "Male/Female × Yes/No"
  )
  return(scenarios)
}

# Calculate multi-gene signature scores
advanced_z <- function(gene_set, expr_matrix) {
  ix <- toupper(rownames(expr_matrix)) %in% toupper(gene_set)
  cat("Genes found:", sum(ix), "\n")
  
  if (sum(ix) > 1) {
    avg_gene_set <- base::apply(X = expr_matrix[ix, , drop = FALSE], MARGIN = 2, FUN = mean, na.rm = TRUE)
    avg_all      <- base::apply(X = expr_matrix,                     MARGIN = 2, FUN = mean, na.rm = TRUE)
    sd_all       <- base::apply(X = expr_matrix,                     MARGIN = 2, FUN = sd,   na.rm = TRUE)
    
    z <- (avg_gene_set - avg_all) * sqrt(sum(ix)) / sd_all
  } else {
    stop("Cannot use sum(ix) genes for signature score calculation")
  }
  return(z)
}

# Calculate cutoffs (expression based survival) for each facet
calc_cutoffs <- function(cutoff_df, survival_params){
  
  stratify_var  <- survival_params$stratify_var
  expr_values <- cutoff_df[[stratify_var]]
  
  # Compute quantiles
  qs <- stats::quantile(expr_values,
                        probs = c(0, 0.25, 0.33, 0.5, 0.66, 0.75),
                        na.rm = TRUE)
  iqr <- stats::IQR(expr_values, na.rm = TRUE)
  
  # Determine cutoffs
  cutoffs <- switch(
    survival_params$cutoff_method,
    "median"   = list(lower = qs["50%"], upper = qs["50%"], middle = NA),
    "tertile"  = list(lower = qs["33%"], upper = qs["66%"], middle = NA),
    "quartile" = list(lower = qs["25%"], upper = qs["75%"], middle = qs["50%"]),
    "thirds"   = {
      lower <- max(qs["0%"] - 1.5*iqr, 0)
      upper <- qs["75%"] + 1.5*iqr
      list(lower = lower + (upper-lower)/3, upper = lower + (upper-lower)*2/3, middle = NA)
    },
    "optimal" = tryCatch({
      res <- survminer::surv_cutpoint(data = cutoff_df,
                                      time = survival_params$time_col,
                                      event = survival_params$status_col,
                                      variables = survival_params$stratify_var)
      list(lower = res$cutpoint$cutpoint, upper = res$cutpoint$cutpoint, middle = NA)
    }, error = function(e) {
      list(lower = NA, upper = NA, middle = NA)
    }))
  
  # Categorize expression into bins
  cutoff_df <- cutoff_df %>%
    dplyr::filter(!is.na(.data[[stratify_var]])) %>%
    dplyr::mutate(model = dplyr::case_when(.data[[stratify_var]] > cutoffs$upper ~ "HIGH",
                                           .data[[stratify_var]] <= cutoffs$lower ~ "LOW",
                                           !is.na(cutoffs$middle) && .data[[stratify_var]] > cutoffs$middle ~ "MED_HIGH",
                                           !is.na(cutoffs$middle) && .data[[stratify_var]] <= cutoffs$middle ~ "MED_LOW",
                                           TRUE ~ "MID"),
                  model = factor(model)) %>%
    dplyr::select(Sample_ID, model)
  
  # Optionally plot only HIGH and LOW
  if (!survival_params$show_all_bins) {
    cutoff_df <- cutoff_df %>%
      dplyr::filter(model %in% c("HIGH", "LOW"))
  }
  return(cutoff_df)
}

# Calculate cox model stats (HR, CI, p vals) for each facet
calc_cox_stats <- function(facet_df, surv_formula, survival_params){
  
  # ---- Cox model ----
  
  # Ensure model column is a factor
  facet_df$model <- factor(facet_df$model)
  
  # Fit Cox model
  cox_model <- survival::coxph(formula = surv_formula, data = facet_df)
  
  # Cox model coefficients
  cox_coef_df <- summary(cox_model)$coefficients
  cox_ci_df <- as.data.frame(confint(cox_model))
  baseline <- levels(facet_df$model)[1]  # the factor baseline
  cox_df <- data.frame(Target = gsub("^model", "", rownames(cox_coef_df)), # remove "model" prefix if present
                       Reference = baseline,
                       HR = exp(cox_coef_df[, "coef"]),
                       CI_lower = exp(cox_ci_df[, 1]),
                       CI_upper = exp(cox_ci_df[, 2]),
                       pval = cox_coef_df[, "Pr(>|z|)"],
                       stringsAsFactors = FALSE) %>%
    dplyr::mutate(contrast = paste0(Target, " / ", Reference)) %>%
    dplyr::select(contrast, HR, pval, CI_lower, CI_upper, Target, Reference) %>%
    tibble::remove_rownames() 
  
  # ---- Optional: emmeans for pairwise contrasts ----
  # Step 1: estimated marginal means on log-hazard scale
  emm <- emmeans::emmeans(object = cox_model, specs = "model")
  
  # Step 2: pairwise contrasts on HR scale
  pairwise <- emmeans::contrast(object = emm, method = "pairwise", type = "response")
  
  # Step 3: add confidence intervals
  pairwise_ci <- stats::confint(pairwise)
  
  # Step 4: combine into clean data.frame
  pairwise_df <- dplyr::left_join(x = as.data.frame(pairwise) %>%
                                    dplyr::select(contrast, ratio, p.value),
                                  y = as.data.frame(pairwise_ci) %>%
                                    dplyr::select(contrast, asymp.LCL, asymp.UCL),
                                  by = c("contrast"="contrast")) %>%
    dplyr::rename(HR = ratio,
                  CI_lower = asymp.LCL,
                  CI_upper = asymp.UCL,
                  pval = p.value) %>%
    dplyr::mutate(Target = sub(" / .*", "", contrast),
                  Reference = sub(".* / ", "", contrast))
  
  # Step 5: Calculate HR, CI for reverse contrasts
  reversed <- pairwise_df %>%
    dplyr::mutate(Reference = sub(" / .*", "", contrast),
                  Target = sub(".* / ", "", contrast),
                  contrast = paste(Target, "/", Reference),
                  # reciprocal HR
                  HR = 1 / HR,
                  # swap CI bounds
                  CI_lower = 1 / CI_upper,
                  CI_upper = 1 / CI_lower)
  
  # Step 6: Merge all stats
  emmeans_df <- dplyr::bind_rows(pairwise_df, reversed)
  
  # Step 7: Compute non-parametric p-values
  emmeans_df <- sapply(X = emmeans_df$contrast, FUN = calc_pvals,
                       data = facet_df, survival_params = survival_params) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(contrast = emmeans_df$contrast) %>%
    tibble::remove_rownames() %>%
    dplyr::right_join(emmeans_df, by=c("contrast" = "contrast"))
  
  cox_df <- sapply(X = cox_df$contrast, FUN = calc_pvals,
                   data = facet_df, survival_params = survival_params) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(contrast = cox_df$contrast) %>%
    tibble::remove_rownames() %>%
    dplyr::right_join(cox_df, by=c("contrast" = "contrast"))
  
  # ---- Return both ----
  list(cox_model_df = cox_df,
       emmeans_df = emmeans_df)
}

# Calculate all 7 non-parametric p-values for each contrast in cox_df or emmeans_df
calc_pvals <- function(contrast, facet_df, survival_params) {
  
  g1 <- sub(".* / ", "", contrast)
  g2 <- sub(" / .*", "", contrast)
  df_pair <- subset(facet_df, model %in% c(g1, g2))
  surv_obj <- survival::Surv(
    time   = df_pair[[survival_params$time_col]],
    event  = df_pair[[survival_params$status_col]],
    type   = "right",
    origin = 0
  )
  surv_form <- surv_obj ~ model
  fit <- survminer::surv_fit(formula = surv_form , data = df_pair)
  
  test_methods <- c("survdiff", "1", "n", "sqrtN", "S1", "S2", "FH_p=1_q=1")
  sapply(X = test_methods, FUN = function(test.method) {
    res <- tryCatch(
      survminer::surv_pvalue(fit = fit,
                             method = test.method,
                             test.for.trend = FALSE,
                             combine = FALSE)[[2]],
      error = function(e) NA
    )})
}

# Plot survival for each facet
plot_facets <- function(facet_df, surv_curve, cox_df, survival_params, surv_type, facet_group){
  
  if (survival_params$plot_curve) {
    
    # Legend labels
    legend_label <- facet_df %>%
      dplyr::pull(model) %>% 
      unique()
    
    # Legend title
    if (surv_type %in% c("single_gene", "multi_gene")){
      legend_title <- paste0(c("Expression", survival_params$substratify_var, survival_params$facet_var), collapse = "_")
    } else if (surv_type == "meta"){
      legend_title <- paste0(c(survival_params$stratify_var, survival_params$substratify_var, survival_params$facet_var), collapse = "_")
    }
    
    # X-axis breaks
    # We want a maximum of 10 timepoint intervals that are multiples of 12
    max_time <- max(facet_df[[survival_params$time_col]], na.rm = TRUE)
    n <- max(floor(max_time / 10 / 12) * 12, 1)
    breaks <- if (max_time %/% n <= 10) n else n + 12
    
    # Plot KM curve
    # NOTE: Use survminer::ggsurvplot() instead of base::plot()
    # ggsurvplot() produces a list of ggplot objects: survival curve and risk table
    surv_plot <- survminer::ggsurvplot(
      fit = surv_curve,
      pval = FALSE,
      palette = survival_params$color_palette,
      linetype = "solid",
      size = 1.5,                  # thickness of line
      
      # Format the legend
      legend = "top",
      legend.title = legend_title,
      legend.labs = legend_label,
      
      # Format the axes
      break.time.by = breaks,       # break X axis in time intervals of 12 months
      xlab = "Time",
      ylab = "Survival Probability",
      title = survival_params$stratify_var,
      
      # Format confidence intervals
      conf.int = survival_params$conf_interval,
      #conf.int.fill = ?,               # color to fill confidence interval
      conf.int.style = "ribbon",        # confidence interval style
      conf.int.alpha = 0.3,             # confidence fill color transparency
      
      # Format the risk table
      risk.table = survival_params$plot_risk_table,
      risk.table.title = "Number at risk",
      risk.table.y.text.col = TRUE,     # color of risk table text annotations
      risk.table.pos = "out",           # draw risk table outside survival plot
      
      # Format the censor points
      censor = TRUE,
      censor.shape = '|',
      censor.size = 4
    )
    
    # Adjust x-axis limits for survival plot and risk table
    surv_plot$table <- surv_plot$table +
      coord_cartesian(x = c(0, ceiling(max_time / breaks) * breaks), clip = "off")
    
    surv_plot$plot <- surv_plot$plot +
      coord_cartesian(x = c(0, ceiling(max_time / breaks) * breaks), clip = "off")
    
    # Create a text annotation showing p-value, HR with CI, and method
    # [Accurate ONLY if 2 curves are present. If more, refer excel file and manually add to image]
    method_plot <- "log-rank"
    p_plot <- formatC(cox_df$pval[1], format = "e", digits = 1)
    hr_plot <- round(cox_df$HR[1], 1)
    ci_lower_plot <- round(cox_df$CI_lower[1], 1)
    ci_upper_plot <- round(cox_df$CI_upper[1], 1)
    
    survplot_stats_grob <- grobTree(textGrob(
      label = paste0("p = ", p_plot,
                     "\nHR = ", hr_plot, " [", ci_lower_plot, ", ", ci_upper_plot, "]",
                     "\nMethod = ", method_plot),
      x = 0.50, y = 0.90, hjust = 0,
      gp = grid::gpar(fontfamily = "Times", fontface = "bold", col = "black", fontsize = 10)
    ))
    
    # Merge text annotation with survival plot
    surv_plot$plot <- surv_plot$plot %++%
      ggplot2::annotation_custom(survplot_stats_grob)
    
    
    # NOTE: Using cowplot() and then ggsave() works nicely as compared to 
    # saving directly using ggsave()
    cowplot::plot_grid(plotlist = surv_plot,
                       align = "hv",
                       axis = "tblr",
                       nrow = 2,
                       ncol = 1,
                       rel_widths = 1,
                       rel_heights = c(1, 0.45),
                       labels = NULL,
                       label_size = 14,
                       label_fontface = "bold")
    
    # File name for plot 
    facet_group <- if(is.na(facet_group)) ""
    method <- survival_params$stratify_criteria
    stratify_var <- survival_params$stratify_var
    file_name <- paste0(c("KM_curve", stratify_var, facet_group, method, ".pdf"), collapse = "_")
    file_name <- base::gsub(pattern = "/", replacement = "-", x = file_name)
    
    # Save the plot
    ggplot2::ggsave(filename = file.path(survival_params$output_path, file_name),
                    plot = last_plot(),
                    device = "pdf",
                    width = 7, height = 7, units = "in")
  }
}

survival_analysis <- function(meta_data, expr_data = NULL, survival_params) {
  
  # ---- Input checks & parameter extraction ----
  
  # Extract parameters
  stratify_var    <- survival_params$stratify_var
  substratify_var <- survival_params$substratify_var
  facet_var       <- survival_params$facet_var
  time_col        <- survival_params$time_col
  status_col      <- survival_params$status_col
  cutoff_method   <- survival_params$cutoff_method
  multiple_cutoff <- survival_params$multiple_cutoff
  show_all_bins   <- survival_params$show_all_bins
  
  # Must provide stratify_var
  if (is.null(stratify_var) || length(stratify_var) == 0) {
    stop("Must provide a non-empty stratify_var")
  }
  
  # substratify_var, if provided, must exist in meta_data
  if (!is.null(substratify_var) && !substratify_var %in% colnames(meta_data)) {
    stop("substratify_var not found in meta_data columns")
  }
  
  # facet_var, if provided, must exist in meta_data
  if (!is.null(facet_var) && !facet_var %in% colnames(meta_data)) {
    stop("facet_var not found in meta_data columns")
  }
  
  # Must define substratify_var if multiple_cutoff = TRUE
  if (isTRUE(multiple_cutoff) && is.null(substratify_var)) {
    stop("multiple_cutoff = TRUE requires substratify_var")
  }
  
  # Check stratify_var against both metadata and expression data and determine
  # type of survival analysis
  missing_genes <- setdiff(stratify_var, rownames(expr_data))
  valid_genes   <- intersect(stratify_var, rownames(expr_data))
  in_meta       <- stratify_var %in% colnames(meta_data)
  
  if (!in_meta) {
    # Expression-based stratification
    if (length(valid_genes) == 0) {
      stop("stratify_var must match either a column in meta_data OR genes in expr_data.")
    } else if (length(valid_genes) == 1) {
      if (length(missing_genes) > 0) {
        warning("Some requested genes not found in expr_data: ", paste(missing_genes, collapse = ", "))
      }
      message("Proceeding with single-gene expression-based survival analysis.")
      surv_type <- "single_gene"
    } else if (length(valid_genes) > 1) {
      if (length(missing_genes) > 0) {
        warning("Some requested genes not found in expr_data: ", paste(missing_genes, collapse = ", "))
      }
      message("Proceeding with multi-gene expression-based survival analysis (", length(valid_genes), " valid genes).")
      surv_type <- "multi_gene"
    }
  } else {
    # Metadata-based stratification
    if (length(valid_genes) == 0) {
      message("Proceeding with metadata-based survival analysis.")
      surv_type <- "meta"
    } else {
      stop("stratify_var matches BOTH a column in meta_data and genes in expr_data — ambiguous.")
    }
  }
  
  
  # ---- Prepare expression data ----
  if (surv_type == "multi_gene") {
    # Multi-gene signature score survival (whole dataset, even if substratify_var defined)
    sig_scores <- advanced_z(gene_set = stratify_var, expr_matrix = expr_data)
    expr_df    <- as.data.frame(sig_scores, check.names = FALSE) %>%
      dplyr::rename(sig_score = identity(1)) %>%
      tibble::rownames_to_column(var = "Sample_ID")
    
    # Update local variable and survival_params
    stratify_var <- "sig_score"
    survival_params$stratify_var <- stratify_var
    
  } else if (surv_type == "single_gene") {
    # Single-gene survival (expression-based)
    expr_df <- expr_data[stratify_var, , drop = FALSE] %>%
      t() %>%
      as.data.frame(check.names = FALSE) %>%
      tibble::rownames_to_column(var = "Sample_ID")
    
  } else if (surv_type == "meta") {
    # Metadata-based survival
    expr_df <- meta_data %>%
      dplyr::select(Sample_ID, dplyr::all_of(stratify_var))
    
  } else {
    stop("Invalid stratify_var: must be gene(s) in expr_data or a column in meta_data.")
  }
  
  
  # ---- Merge expression data with metadata ----
  
  keep_cols <- c(time_col, status_col, stratify_var, substratify_var, facet_var)
  
  surv_df <- expr_df %>%
    dplyr::inner_join(meta_data, by = c("Sample_ID"="Sample_ID")) %>%
    dplyr::select(Sample_ID, dplyr::all_of(keep_cols))
  
  
  # ---- Define model column for metadata based survival ----
  
  if (surv_type == "meta") {
    surv_df <- surv_df %>%
      dplyr::mutate(model = paste(.data[[stratify_var]],
                                  if (!is.null(substratify_var)) .data[[substratify_var]] else NULL,
                                  if (!is.null(facet_var)) .data[[facet_var]] else NULL,
                                  sep = "_"))
  }

  
  # ---- Define model column for expression based survival ----
  
  if (surv_type %in% c("multi_gene", "single_gene")) {
    
    # Define facet groups
    if (!is.null(facet_var) && facet_var %in% colnames(surv_df)) {
      facet_groups  <- unique(surv_df[[facet_var]])
    } else {
      facet_groups  <- NA_character_ # placeholder for whole dataset
    }
    
    merged_df <- tibble::tibble()
    
    # Loop over each facet group to calculate cutoffs
    for (facet_group in facet_groups){
      
      # Subset by facet (or use full data if NA)
      if (is.na(facet_group)){
        facet_df <- surv_df
      } else{
        facet_df <- surv_df %>%
          dplyr::filter(.data[[facet_var]] == facet_group)
      }
      
      # Determine cutoff groups
      if (!is.null(substratify_var) && substratify_var %in% colnames(surv_df) && isTRUE(multiple_cutoff)) {
        cutoff_groups <- unique(facet_df[[substratify_var]])
      } else {
        cutoff_groups <- NA_character_
      }
      
      # Calculate cutoffs for each cutoff_group
      for (cutoff_group in cutoff_groups){
        
        # Subset by cutoff_group (or use facet data if NA)
        if (is.na(cutoff_group)) {
          cutoff_df <- facet_df
        } else {
          cutoff_df <- facet_df %>% 
            dplyr::filter(.data[[substratify_var]] == cutoff_group)
        }
        
        # Calculate cutoffs for this subset
        df <- calc_cutoffs(cutoff_df, survival_params)
        
        merged_df <- dplyr::bind_rows(merged_df, df)
      }
    }
    
    # Merge classifications into surv_df
    surv_df <- surv_df %>%
      dplyr::left_join(merged_df, by=c("Sample_ID"="Sample_ID"))
    
    # Append substratify_var to model if defined
    if (!is.null(substratify_var)) {
      surv_df <- surv_df %>%
        dplyr::mutate(model = paste0(model, "_", .data[[substratify_var]]))
    }
  }
  
  
  # ---- Main loop: facets, compute stats, save, plot ----
  
  # Define facet groups
  if (!is.null(facet_var) && facet_var %in% colnames(surv_df)) {
    facet_groups  <- unique(surv_df[[facet_var]])
  } else {
    facet_groups  <- NA  # placeholder for whole dataset
  }
  
  # Plot a survival curve for each facet
  for (facet_group in facet_groups){
    
    # Subset by facet (or use full data if NA)
    if (is.na(facet_group)){
      facet_df <- surv_df
    } else{
      facet_df <- surv_df %>%
        dplyr::filter(.data[[facet_var]] == facet_group)
    }
    
    # Check if each facet has least 2 groups in model column
    if (dplyr::n_distinct(facet_df$model) < 2) {
      message("Skipping facet (single group): ", facet_group)
      next
    } else if (dplyr::n_distinct(facet_df$model) == 2){
      survival_params$color_palette <- c("#d73027","#0c2c84")
    }
    
    # Create a survival object (Alive = 0, Dead = 1)
    surv_object <- survival::Surv(
      time   = facet_df[[time_col]],
      event  = facet_df[[status_col]],
      type   = "right",
      origin = 0
    )
    
    # Create a formula for survival analysis
    surv_formula <- surv_object ~ model
    
    # Create a fit for Kaplan-Meier curve
    # NOTE: survival::survfit() gives error in ggsurvplot()
    surv_curve <- survminer::surv_fit(
      formula   = surv_formula,
      data      = facet_df
      )
    
    # Compute Cox & emmeans stats for each facet
    stats_list <- calc_cox_stats(facet_df, surv_formula, survival_params)
    
    # Save stats for each facet
    wb <- createWorkbook()
    openxlsx::addWorksheet(wb, sheetName = "cox_stats")
    openxlsx::writeData(wb, sheet = "cox_stats", x = stats_list$cox_model_df)
    openxlsx::addWorksheet(wb, sheetName = "emmeans_stats")
    openxlsx::writeData(wb, sheet = "emmeans_stats", x = stats_list$emmeans_df)
    openxlsx::addWorksheet(wb, sheetName = "facet_df")
    openxlsx::writeData(wb, sheet = "facet_df", x = facet_df)
    file_name <- paste0("Survival_data_", ifelse(is.na(facet_group), "all", facet_group), ".xlsx")
    openxlsx::saveWorkbook(wb, file = file.path(survival_params$output_path, file_name), overwrite = TRUE)
    
    # Plot survival curve for each facet
    cox_df <- stats_list$cox_model_df
    plot_facets(facet_df, surv_curve, cox_df, survival_params, surv_type, facet_group)
    }
}

survival_params <- list(
  
  # Expression + Metadata-based survival parameters
  stratify_var    = NULL,                 # one or more genes or metadata column for stratifying survival curves
  substratify_var = NULL,                 # optional metadata column for additional stratification
  facet_var       = NULL,                 # optional metadata column for faceting plots
  
  # Expression-based survival parameters
  cutoff_method   = "optimal",            # method for defining cutoffs: mean, median, quartile, tertile, thirds, optimal
  show_all_bins   = FALSE,                # TRUE = plot all bins (LOW, HIGH, MID/MED_HIGH/MED_LOW)
  multiple_cutoff = FALSE,                # TRUE = compute cutoffs separately for substratify_var
  
  # Common parameters
  conf_interval   = FALSE,                # TRUE = show confidence interval in survival curve
  plot_curve      = TRUE,                 # TRUE = plot the survival curve
  plot_risk_table = TRUE,                 # TRUE = plot the risk table below the curve
  time_col        = "Time",               # metadata column containing Time values
  status_col      = "Status",             # metadata column containing Status values
  color_palette   = custom_palette,       # vector of colors for groups c("#d73027","#0c2c84")
  prefix          = "",                   # file prefix
  output_path     = "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop"
)

# Run to understand how to define parameters for the survival function
df <- show_survival_scenarios()

# Run your survival analysis
#survival_analysis(meta_data, expr_data, survival_params)

#******************************************************************************#
#                       SURVIVAL CURVE RELATED FUNCTIONS                       #
#******************************************************************************#

# Read this paper for survival analysis
# https://doi.org/10.1093/jncimonographs/lgu024

# NOTE:  Output of prep_expr_df is df
#log_norm_counts is matrix with SYMBOLS as rownames
prep_expr_df <- function(log_norm_counts, meta_data, plot_genes, survival_params){
  
  # Merge expression data with survival data
  if (survival_params$gene_sig_score == TRUE){
    
    # Calculate gene signature score
    expr_df <- as.data.frame(advanced_Z(plot_genes, log_norm_counts))
    
    expr_df <- expr_df %>%
      data.frame() %>%
      dplyr::rename(combined.exp = identity(1)) %>%
      tibble::rownames_to_column("Sample_ID") %>%
      dplyr::inner_join(meta_data, by=c("Sample_ID"="Sample_ID")) %>%
      dplyr::select(Sample_ID, combined.exp, Time, Status)
  } else {
    expr_df <- log_norm_counts %>%
      t() %>%
      data.frame() %>%
      tibble::rownames_to_column("Sample_ID") %>%
      dplyr::inner_join(meta_data, by=c("Sample_ID"="Sample_ID")) %>%
      dplyr::select(Sample_ID, all_of(plot_genes), Time, Status)
  }
  
  # Add split_by column to expr_df to define groups in order to calculate multiple_cutoff
  if (!is.na(survival_params$split_by)){
    expr_df <- expr_df %>% 
      dplyr::left_join(meta_data %>% dplyr::select(Sample_ID, survival_params$split_by),
                       by=c("Sample_ID"="Sample_ID"))
  }
  
  return(expr_df)
}

# NOTE:  Output of calc_cutoffs is list(df,ls)
# If plotting by Sex, make sure to create column "model" based on which lines will be split
calc_cutoffs <- function(df, gene, group, survival_params){
  
  # Identify upper & lower cutoffs based on stratify_criteria
  #*************************Split samples by median**************************#
  if(survival_params$stratify_criteria == "m"){
    quartiles <- stats::quantile(x = df[[gene]],
                                 probs = c(0, 0.25, 0.50, 0.75, 1),
                                 na.rm=TRUE)
    
    cutoff_lower_end <- quartiles[[3]]
    cutoff_upper_end <- quartiles[[3]]
    cutoff_middle <- "NA"
  }
  
  #****************Split samples into top and bottom tertiles****************#
  else if(survival_params$stratify_criteria == "t"){
    tertiles <- stats::quantile(x = df[[gene]],
                                probs = c(0, 0.33, 0.66, 1),
                                na.rm=TRUE)
    
    cutoff_lower_end <- tertiles[[2]]
    cutoff_upper_end <- tertiles[[3]]
    cutoff_middle <- "NA"
  }
  
  #***************Split samples into top and bottom quartiles****************#
  else if(survival_params$stratify_criteria == "q"){
    quartiles <- stats::quantile(x = df[[gene]], 
                                 probs = c(0, 0.25, 0.50, 0.75, 1),
                                 na.rm=TRUE)
    
    cutoff_lower_end <- quartiles[[2]]
    cutoff_upper_end <- quartiles[[4]]
    cutoff_middle <- quartiles[[3]]
  }
  
  #*********************Split expression range by thirds*********************#
  else if(survival_params$stratify_criteria == "th"){
    quartiles <- stats::quantile(x = df[[gene]], 
                                 probs = c(0, 0.25, 0.50, 0.75, 1),
                                 na.rm=TRUE)
    iqr <- stats::IQR(x = df[[gene]],
                      na.rm=TRUE)
    
    # Normal range of expression values lie between cutoff_lower & cutoff_upper
    cutoff_upper <- quartiles[[4]]+1.5*iqr
    cutoff_lower <- dplyr::if_else(quartiles[[1]]-1.5*iqr > 0, quartiles[[1]]-1.5*iqr, 0)
    
    # Based on normal range of expression, identify onethird & twothird cutoff
    cutoff_lower_end <- cutoff_lower + (cutoff_upper-cutoff_lower)/3
    cutoff_upper_end <- cutoff_lower + (cutoff_upper-cutoff_lower)*2/3
    cutoff_middle <- "NA"
  }
  
  #***************Split expression range using optimal cutoff****************#
  else if(survival_params$stratify_criteria == "o"){
    
    # Sometimes quartiles will look like: 
    # 0%       25%      50%      75%     100% 
    # 0.000000 0.000000 0.000000 0.000000 3.495493 
    # In such cases, surv_cutpoint() will fail. So, we add extra if() here.
    quartiles <- stats::quantile(x = df[[gene]], 
                                 probs = c(0, 0.25, 0.50, 0.75, 1),
                                 na.rm=TRUE)
    
    if (quartiles[[4]] > quartiles[[2]]){
      res.cut <- survminer::surv_cutpoint(data = df,
                                          time = "Time",
                                          event = "Status",
                                          variables = gene)
      
      cutoff_lower_end <- res.cut$cutpoint$cutpoint
      cutoff_upper_end <- res.cut$cutpoint$cutpoint
      cutoff_middle <- "NA"
    } else{
      #cat("Surv cutpoint unable to detect optimum cutoff")
      cutoff_lower_end <- "NA"
      cutoff_upper_end <- "NA"
      cutoff_middle <- "NA"
    }
  }
  
  # Categorize the sample based on above cutoffs
  if (survival_params$plot_all_bins == TRUE & survival_params$stratify_criteria == "q"){
    df <- df %>% 
      dplyr::mutate(Expression = dplyr::case_when(get(gene) > cutoff_upper_end ~ "HIGH",
                                                  get(gene) <= cutoff_lower_end ~ "LOW",
                                                  get(gene) <= cutoff_middle ~ "MED_LOW",
                                                  TRUE ~ "MED_HIGH"))
    
  } else if (survival_params$plot_all_bins == TRUE) {
    df <- df %>% 
      dplyr::mutate(Expression = dplyr::case_when(get(gene) > cutoff_upper_end ~ "HIGH",
                                                  get(gene) <= cutoff_lower_end ~ "LOW",
                                                  TRUE ~ "MID"))
    
  } else if (survival_params$stratify_criteria == "none") {
    #When plotting by Sex, Treatment response, we dont use expression data.
    df <- df %>% 
      dplyr::mutate(Expression = model)
    cutoff_lower_end <- NA
    cutoff_upper_end <- NA
    cutoff_middle <- NA
    
  } else {
    df <- df %>% 
      dplyr::mutate(Expression = dplyr::case_when(get(gene) > cutoff_upper_end ~ "HIGH",
                                                  get(gene) <= cutoff_lower_end ~ "LOW",
                                                  TRUE ~ "MID")) %>%
      dplyr::filter(Expression != "MID")
  }
  
  # # Print the cutoffs
  # cat("\nGene         :", gene)
  # cat("\nGroup        :", group)
  # cat("\nLower cutoff :", cutoff_lower_end)
  # cat("\nUpper cutoff :", cutoff_upper_end)
  # cat("\nMiddle cutoff:", cutoff_middle)
  
  # Create a list to store cutoff values
  ls <- list("group" = c(), 
             "gene" = c(), 
             "lower" = c(), 
             "upper" = c(), 
             "middle" = c())
  
  ls$group <- c(group)
  ls$gene <- c(gene)
  ls$lower <- c(cutoff_lower_end)
  ls$upper <- c(cutoff_upper_end)
  ls$middle <- c(cutoff_middle)
  
  # Return the df and the cutoffs
  return(list(df, ls))
}

# NOTE:  Output of calc_surv_stats is list. 
# It also generate survival plot with risk table
calc_surv_stats <- function(df, group, prefix, output_path){
  
  # If all samples belong to one group (like LOW or HIGH or males or female),
  # then quit the function as comparison cannot be done
  if (nrow(df %>% dplyr::count(model)) > 1){
    
    # Create a survival object where Alive = 0, Dead = 1
    surv_object <- survival::Surv(time = df$Time,
                                  event = df$Status,
                                  type = "right",
                                  origin = 0)
    
    # Create a formula for plotting survival curve
    surv_formula <- surv_object ~ model
    
    # Create a fit for survival curve.
    # NOTE: survival::survfit() gives error in ggsurvplot(). Use survminer::surv_fit()
    surv_curve <- survminer::surv_fit(formula = surv_formula,
                                      data = df,
                                      type = "kaplan-meier",
                                      group.by = NULL,
                                      match.fd = FALSE)
    
    # # Check summary of the survival curve with time duration of our interest
    # cat("\nRange of survival (months):", range(df$Time, na.rm=TRUE), "\n")
    # base::summary(surv_curve, times = base::seq(from = floor(range(df$Time, na.rm=TRUE)[[1]]),
    #                                                 to = ceiling(range(df$Time, na.rm=TRUE)[[2]]),
    #                                                 by = 3))
    
    # Create a Cox model for the survival curve and calculate stats
    cox_model <- survival::coxph(formula = surv_formula,
                                 data = df)
    #print(summary(cox_model))
    cat("\n")
    
    # Calculate HR, 95% CI for HR, p-val
    # NOTE: Variable mentioned in summary(cox_model) is numerator in h1(t)/h0(t).
    # The reference variable h0(t) will not be mentioned in co-efficients.
    # Make sure this is not the reference level i.e. low expression. If this is
    # the reference, then we need to reverse the HR ratio, legend labels
    #print(names(cox_model$coefficients))  
    
    # If samples belong to more than 2 groups (like LOW, MID, HIGH), then we 
    # cannot have survival stats. So, we set them to 0.
    if (nrow(df %>% dplyr::count(model)) == 2){
      # Store HR and CI
      if (stringr::str_detect(names(cox_model$coefficients), survival_params$reference)){
        HR <- round(exp(-cox_model$coefficients[[1]]), 2)
        CI <- round(exp(-confint(cox_model)), 2)
        CI_1 <- CI[1]
        CI[1] <- CI[2]
        CI[2] <- CI_1
      } else {
        HR <- round(exp(cox_model$coefficients[[1]]),2)
        CI <- round(exp(confint(cox_model)),2)
      }
      
      # Store pvalues by different methods
      pvals <- c()
      for (test.method in c("survdiff", "1", "n", "sqrtN", "S1","S2", "FH_p=1_q=1")){
        
        # Some of the methods give error. So, we catch them and skip
        if (sum(str_detect(string = class(tryCatch(survminer::surv_pvalue(fit = surv_curve,
                                                                          method = test.method,
                                                                          test.for.trend = FALSE,
                                                                          combine = FALSE), error = function(e) e)),
                           pattern = "error")) == 0){
          p_val <- survminer::surv_pvalue(fit = surv_curve,
                                          method = test.method,
                                          test.for.trend = FALSE,
                                          combine = FALSE)
          pvals <- c(pvals, p_val[[2]])
        } else{
          pvals <- c(pvals, 1)
          print(test.method)
        } 
      } 
    } else {
      HR <- 0
      CI <- c(0, 0)
      pvals <- c(0,0,0,0,0,0,0) # 7 pvalues for each method
    }
    
    # Plot the survival curve using survminer::ggsurvplot() instead of base::plot()
    # ggsurvplot() produces a list of ggplot objects: a survival curve and a risk table
    # Saving it using cowplot() first and then using ggsave() works nicely as
    # compared to saving directly using ggsave()
    if(survival_params$plot_curve == TRUE){
      
      # Plot the survival curve
      legend_label <- df %>% 
        dplyr::count(model) %>% 
        dplyr::select(model) %>% 
        unlist(.,use.names=FALSE)
      
      # We identify proper breaks based on max duration of the dataset
      # We want a maximum of 10 timepoint intervals that are multiples of 12
      max_time <- max(df$Time,na.rm=TRUE)
      n <- floor(max_time/10/12)*12
      if(max_time %/% n <= 10){
        breaks <- n
      } else{
        breaks <- n+12
      }
      
      surv_plot <- survminer::ggsurvplot(fit = surv_curve,
                                         pval = FALSE,
                                         palette = survival_params$color_palette,
                                         linetype = "solid",
                                         size = 1.5,                       # thickness of line
                                         
                                         # Format the legend
                                         legend  = "top",                  # position of legend
                                         legend.title = survival_params$legend_title,
                                         legend.labs = survival_params$legend_label,
                                         
                                         # Format the axes
                                         break.time.by = breaks,           # break X axis in time intervals of 12 months
                                         xlab = "Time (Months)",           # customize X axis label
                                         ylab = "Survival Probability",    # customize Y axis label
                                         title = dplyr::if_else(gene == "combined.exp", "", gene),
                                         
                                         # Format confidence intervals
                                         conf.int = survival_params$conf_interval,
                                         #conf.int.fill = ?,               # color to fill confidence interval
                                         conf.int.style = "ribbon",        # confidence interval style
                                         conf.int.alpha = 0.3,             # confidence fill color transparency
                                         
                                         # Format the risk table
                                         risk.table = survival_params$plot_risk_table,
                                         risk.table.title = "Number at risk",
                                         risk.table.y.text.col = TRUE,     # color of risk table text annotations
                                         risk.table.pos = "out",           # draw risk table outside survival plot
                                         
                                         # Format the censor points
                                         censor = TRUE,
                                         censor.shape = '|',
                                         censor.size = 5)
      
      surv_plot$table <- surv_plot$table + 
        coord_cartesian(x=c(0,ceiling(max_time/breaks)*breaks), clip = "off")
      
      surv_plot$plot <- surv_plot$plot + 
        coord_cartesian(x=c(0,ceiling(max_time/breaks)*breaks), clip = "off")
      
      # Plot p and HR value
      method_plot <- "log-rank"
      p_plot <- pvals[1]  
      
      survplot_stats_grob <- grobTree(textGrob(label = paste0("p = ", formatC(p_plot, format = "e", digits = 1),
                                                "\nHR = ", round(HR,1), " [", round(CI[1],1), ", ", round(CI[2],1), "]",
                                                "\nMethod = ", method_plot),
                                 x = 0.50,
                                 y = 0.90,
                                 hjust = 0,
                                 gp = grid::gpar(fontfamily="Times", fontface="bold", col="black", fontsize=10)))
      
      # Add p values and HR values to plot
      surv_plot$plot <- surv_plot$plot %++%
        ggplot2::annotation_custom(survplot_stats_grob)
      
      cowplot::plot_grid(plotlist = surv_plot,
                         align = "hv",
                         axis = "tblr",
                         nrow = 2,  
                         ncol = 1, 
                         rel_widths = 1,
                         rel_heights = c(1,0.45),
                         labels = NULL,
                         label_size = 14,
                         label_fontfamily = NULL,
                         label_fontface = "bold",
                         label_colour = NULL)
      
      f_name <- paste0(prefix, "_", group, "_", survival_params$stratify_criteria, ".pdf")
      f_name <- gsub("/", "-", x=f_name)
      
      # Save the plot
      ggplot2::ggsave(filename = f_name,
                      plot = last_plot(),
                      device = "pdf",
                      path = output_path,
                      width = 7,
                      height = 7,
                      units = c("in"),
                      dpi = 300,
                      limitsize = TRUE,
                      bg = NULL)
    }
  } else {
    HR <- 0
    CI <- c(0, 0)
    pvals <- c(0,0,0,0,0,0,0) # 7 pvalues for each method
  }
  
  # Create a list to store survival stats
  ls <- list("group" = c(), 
             "HR" = c(), 
             "CI_lower" = c(), 
             "CI_upper" = c(), 
             "pvalue" =c(), 
             "logrank" = c(), 
             "reg_logrank.late" = c(), 
             "Gehan_Breslow.early" = c(),
             "Tarone_Ware.early" = c(), 
             "Peto_Peto.early" = c(),  
             "modified_Peto_Peto" = c(), 
             "Fleming_Harrington" = c())
  
  ls$group               <- c(group)
  ls$HR                  <- c(HR)
  ls$CI_lower            <- c(CI[1])
  ls$CI_upper            <- c(CI[2])
  ls$logrank             <- c(pvals[1])
  ls$reg_logrank.late    <- c(pvals[2])
  ls$Gehan_Breslow.early <- c(pvals[3])
  ls$Tarone_Ware.early   <- c(pvals[4])
  ls$Peto_Peto.early     <- c(pvals[5])
  ls$modified_Peto_Peto  <- c(pvals[6])
  ls$Fleming_Harrington  <-c(pvals[7])
  
  return(ls)
}

# NOTE: Output of plot_survival is list(df,ls)
plot_survival <- function(expr_df, gene, survival_params, prefix, output_path){
  
  # Create an empty dataframe to store expr_df and classification from calculate_cutoffs()
  survival_data <- data.frame(model = " ")
  
  # Create a list to store results of calculate_cutoffs, surv_plot et
  stats <- list("gene" = c(), 
                "group" = c(),
                "lower_cutoff" = c(),
                "middle_cutoff" = c(),
                "upper_cutoff" = c(),
                "HR" = c(), 
                "CI_lower" = c(), 
                "CI_upper" = c(), 
                "logrank" = c(), 
                "reg_logrank.late" = c(), 
                "Gehan_Breslow.early" = c(),
                "Tarone_Ware.early" = c(), 
                "Peto_Peto.early" = c(),  
                "modified_Peto_Peto" = c(), 
                "Fleming_Harrington" = c())
  
  # Create a list of groups for multiple_cutoff calculation 
  if (survival_params$multiple_cutoff == TRUE & !is.na(survival_params$split_by)){
    cutoff_groups <- expr_df %>% 
      dplyr::add_count(get(survival_params$split_by)) %>%
      dplyr::filter(n>2) %>%
      dplyr::select(all_of(survival_params$split_by)) %>% 
      unlist(use.names=FALSE) %>% 
      unique()
  } else if (survival_params$multiple_cutoff == TRUE & is.na(survival_params$split_by)){
    print("Please define survival_params$split_by to calculate multiple cutoffs")
  } else {
    cutoff_groups <- c(NA)
  }
  
  # STEP 1: Calculate cutoffs
  # If cutoffs need to be calculated for each group, subset the expr_df and pass
  # it to calculate_cutoffs(). Else, pass entire expr_df to calculate_cutoffs()
  for (group in cutoff_groups){
    
    # Subset the expr_df for each group to calculate cutoffs
    if (survival_params$multiple_cutoff == TRUE & !is.na(survival_params$split_by)){
      df <- expr_df %>% dplyr::filter(get(survival_params$split_by) == group)
    } else if (survival_params$multiple_cutoff == TRUE & is.na(survival_params$split_by)){
      cat("\n 'split_by' variable is undefined but 'multiple_cutoff' is set to TRUE")
    } else{
      df <- expr_df
    }
    
    # Calculate cutoffs for each group
    mat <- calc_cutoffs(df, gene, group, survival_params)
    
    ##### Save the data from output of calculate_cutoffs()
    survival_data       <- dplyr::bind_rows(survival_data, mat[[1]])
    stats$gene          <- c(stats$gene,          mat[[2]]$gene)
    #stats$group         <- c(stats$group,         mat[[2]]$group)
    stats$lower_cutoff  <- c(stats$lower_cutoff,  mat[[2]]$lower)
    stats$middle_cutoff <- c(stats$middle_cutoff, mat[[2]]$middle)
    stats$upper_cutoff  <- c(stats$upper_cutoff,  mat[[2]]$upper)
  }
  
  # Populate the model variable by concatenating "Expression" and "split_by"
  if (!is.na(survival_params$split_by)){
    survival_data <- survival_data %>%
      dplyr::mutate(model = paste0(Expression, "_", get(survival_params$split_by))) %>%
      dplyr::filter(!is.na(Sample_ID))
  } else {
    survival_data <- survival_data %>%
      dplyr::mutate(model = Expression) %>%
      dplyr::filter(!is.na(Sample_ID))
  }
  
  # Create a list of groups for plotting survival curves 
  if (!is.na(survival_params$split_by)){
    plot_groups <- expr_df %>% 
      dplyr::add_count(get(survival_params$split_by)) %>%
      dplyr::filter(n>2) %>%
      dplyr::select(all_of(survival_params$split_by)) %>% 
      unlist(use.names=FALSE) %>% 
      unique()
  } else {
    plot_groups <- c(NA)
  }
  
  # STEP 2: Calculate survival stats
  # If each group has to be plotted in separate plots, subset the survival_data
  # and pass it to calc_surv_stats(). Else, pass entire survival_data to 
  # calc_surv_stats().
  for (group in plot_groups){
    
    # Subset the survival_data for each group to generate separate plots
    if (survival_params$split_plot == TRUE & !is.na(survival_params$split_by)){
      df <- survival_data %>% dplyr::filter(get(survival_params$split_by) == group)
    } else if (survival_params$split_plot == TRUE & is.na(survival_params$split_by)){
      cat("\n 'split_by' variable is undefined but 'split_plot' is set to TRUE")
    } else{
      df <- survival_data
    }
    
    # Calculate survival stats for each group
    cox_stats <- calc_surv_stats(df, group, prefix, output_path)
    
    ##### Save the data from output of calc_surv_stats()
    stats$group               <- c(stats$group,               cox_stats$group)
    stats$HR                  <- c(stats$HR,                  cox_stats$HR)
    stats$CI_lower            <- c(stats$CI_lower,            cox_stats$CI_lower)
    stats$CI_upper            <- c(stats$CI_upper,            cox_stats$CI_upper)
    stats$logrank             <- c(stats$logrank,             cox_stats$logrank)
    stats$reg_logrank.late    <- c(stats$reg_logrank.late,    cox_stats$reg_logrank.late)
    stats$Gehan_Breslow.early <- c(stats$Gehan_Breslow.early, cox_stats$Gehan_Breslow.early)
    stats$Tarone_Ware.early   <- c(stats$Tarone_Ware.early,   cox_stats$Tarone_Ware.early)
    stats$Peto_Peto.early     <- c(stats$Peto_Peto.early,     cox_stats$Peto_Peto.early)
    stats$modified_Peto_Peto  <- c(stats$modified_Peto_Peto,  cox_stats$modified_Peto_Peto)
    stats$Fleming_Harrington  <- c(stats$Fleming_Harrington,  cox_stats$Fleming_Harrington)
  }
  return(list(survival_data, stats))
}

# We need to calculate a combined expression value of all genes in the gene
# signature for each sample. Next, we use surv_cutpoint() to classify samples 
# into high and low groups and plot survival curves.
# There are several approaches to calculate the combined expression value. While
# normal z-scaling is logical, https://doi.org/10.1186/gb-2006-7-10-r93 
# recommends a more advanced z-scaling which is implemented below. Refer
# "Measuring gene set expression" in "Materials & stratify_criteria" section.


normal_Z <- function(gset, eset) {
  
  # Compute z-score for gene set of interest
  eset <- eset[gset,]
  a <- t(scale(t(eset)))
  z <- colSums(a, na.rm=TRUE) 
  
  return(z)
}


# ---- DEPRECATED FUNCTIONS ----

# DEPRECATED (used during Seurat v3)
v3_sctransform_singlecell <- function(filtered.seurat){
  
  # Seurat v5 stores counts of each sample in separate layers. Merge them.
  filtered.seurat@assays$RNA <- SeuratObject::JoinLayers(filtered.seurat@assays$RNA)
  
  # Split each sample into a seurat object to get a list of seurat object
  split.seurat <- Seurat::SplitObject(object = filtered.seurat,
                                      split.by = "Sample")
  
  # Remove samples with less than 50 cells so that RunPCA() doesnt give error
  split.seurat <- split.seurat[names(split.seurat)[sapply(split.seurat, ncol) > 50]]
  
  for (i in 1:length(split.seurat)){
    
    # Normalize the data
    split.seurat[[i]] <- Seurat::NormalizeData(object = split.seurat[[i]],
                                               assay = "RNA",
                                               normalization.method = "LogNormalize",
                                               scale.factor = 10000,
                                               margin = 1,
                                               verbose = TRUE)
    
    # Perform cell cycle scoring
    split.seurat[[i]]  <- Seurat::CellCycleScoring(object = split.seurat[[i]],
                                                   s.features = intersect(s_genes,rownames(split.seurat[[i]]@assays$RNA@features)),
                                                   g2m.features = intersect(g2m_genes, rownames(split.seurat[[i]]@assays$RNA@features)),
                                                   ctrl = NULL)
    
    split.seurat[[i]]$CC.Score <- split.seurat[[i]]$G2M.Score-split.seurat[[i]]$S.Score
    
    # SCTransform() is better than FindVariableFeatures() & ScaleData()
    # split.seurat[[i]] <- Seurat::FindVariableFeatures(object = split.seurat[[i]],
    #                                                   assay = "RNA",
    #                                                   selection.method = "vst",
    #                                                   nfeatures = 2000)
    # split.seurat[[i]] <- Seurat::ScaleData(object = split.seurat[[i]],
    #                                        features = NULL,
    #                                        assay = "RNA",
    #                                        vars.to.regress = NULL)
    
    # Perform scaling & variable feature identification usign SCTransform()
    split.seurat[[i]] <- Seurat::SCTransform(object =  split.seurat[[i]],
                                             assay = "RNA",
                                             new.assay.name = "SCT",
                                             do.correct.umi = TRUE,
                                             ncells = 5000,
                                             variable.features.n = 3000,
                                             vars.to.regress = c("CC.Score","MitoRatio"),
                                             do.scale = FALSE,
                                             do.center = TRUE,
                                             vst.flavor = "v2",
                                             return.only.var.genes = TRUE)
    
    
    # Remove ribosomal, Riken, predicted and mitochondrial genes from
    # VariableFeatures so that PCA, UMAP and hence clustering are not affected
    var_f <- split.seurat[[i]]@assays$SCT@var.features
    var_f <- var_f[!grepl(pattern = "^[Rr][Pp][SsLl]|R[Ii][Kk]$|^[Mm][Tt]-|^G[Mm][0-9.]+$", 
                          x=var_f)]
    
    split.seurat[[i]]@assays$SCT@var.features <- var_f
    cat("\nFinal number of variable features:", length(var_f), "\n")
    
    # Perform dimensional reduction using PCA on SCT assay variable features
    split.seurat[[i]] <- Seurat::RunPCA(object = split.seurat[[i]],
                                        assay = "SCT",
                                        features = NULL,
                                        ndims.print = 1,
                                        nfeatures.print = 1,
                                        reduction.name = "pca",
                                        reduction.key = "PC_")
    
    # Perform dimensional reduction using UMAP on PCA dimensions
    split.seurat[[i]] <- Seurat::RunUMAP(object = split.seurat[[i]],
                                         dims = 1:40,
                                         reduction = "pca",
                                         reduction.name = "umap",
                                         reduction.key = "UMAP_")
    
  }
  
  return(split.seurat)
}

# DEPRECATED (used during Seurat v3)
v3_integrate_singlecell <- function(sct, ref_samples){
  
  # NOTE: In v3, Harmony integration was not possible as FindIntegrationAnchors()
  # doesnt support reduction="harmony". Moreover, integration using rpca, cca, 
  # jpca couldnt be stored in same object since FindIntegrationAnchors() output
  # varies for each method
  
  #***STEP 7B: SELECT 3000 MOST VARIABLE GENES TO USE FOR INTEGRATING THE DATA***#
  integ_features <- Seurat::SelectIntegrationFeatures(object.list=split.seurat,
                                                      nfeatures=3000,
                                                      assay=NULL) #c("SCT", "SCT"),
  
  
  #******************STEP 7C: FIND RESIDUALS FOR MISSING GENES*******************#
  split.seurat <- Seurat::PrepSCTIntegration(object.list=split.seurat,
                                             assay="SCT",
                                             anchor.features=integ_features)
  
  #******STEP 7D: FIND COMMON ANCHORS BETWEEN SAMPLES TO INTEGRATE THE DATA******#
  integ_anchors <- Seurat::FindIntegrationAnchors(object.list=split.seurat,
                                                  reference=ref_samples,
                                                  anchor.features=integ_features,
                                                  scale=TRUE,
                                                  normalization.method="SCT",
                                                  sct.clip.range=NULL,
                                                  reduction="rpca", #"cca", "jpca", "rlsi"
                                                  l2.norm=TRUE,
                                                  dims=1:30)
  
  #******STEP 7E: FIND OPTIMUM k.weight FOR USE IN Seurat::IntegrateData()*****#
  # Find minimum anchors between 2 datasets
  kweight1 <- as.data.frame(integ_anchors@anchors) %>%
    dplyr::group_by(dataset1, dataset2) %>%
    distinct_at("cell1", .keep_all=TRUE) %>%
    dplyr::summarize(n=n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(n) %>%
    unlist(use.names=FALSE) %>%
    min()
  
  # Find half of number of cells in sample with least cell count
  kweight2 <- filtered.seurat@meta.data %>%
    dplyr::count(Sample) %>%
    dplyr::filter(n >=50) %>%
    dplyr::select(n) %>%
    unlist(use.names=FALSE) %>%
    min()
  
  kweight2 <- floor(kweight2/2)
  
  kweight <- base::min(kweight1, kweight2, 100)
  dplyr::if_else(kweight >= 100, 100, kweight)
  cat("\n", celltype, "\tkweight1:", kweight1, "\tkweight2:", kweight2, "\tkweight:",  kweight, "\n")
  
  # NOTE: Integration will not fail anymore. If it fails, identify the 2
  # datasets that are involved in the error and use kweight=number of anchors
  # for these 2 datasets.
  cat("\nNumber of unique anchors between datasets\n")
  print(as.data.frame(integ_anchors@anchors) %>%
          dplyr::group_by(dataset1, dataset2) %>%
          distinct_at("cell1", .keep_all=TRUE) %>%
          dplyr::summarize(n=n()), n=1000)
  
  #************************STEP 7F: INTEGRATE THE DATA*************************#
  # NOTE: weight.reduction=NULL means new PCA will be calculated & used to
  # calculate anchor weights
  integrated.seurat.rpca <- Seurat::IntegrateData(anchorset=integ_anchors.rpca,
                                                  new.assay.name="integrated",
                                                  normalization.method="SCT",
                                                  features=NULL,
                                                  features.to.integrate=NULL,
                                                  dims=1:30,
                                                  k.weight=kweight, #default is 100
                                                  weight.reduction=NULL,
                                                  sd.weight=1)
  
  #**STEP 7G: RUN PCA USING 3000 INTEGRATION FEATURES & UMAP USING FIRST 40 PCs**#
  integrated.seurat <- Seurat::RunPCA(object=integrated.seurat,
                                      assay="integrated",
                                      features=NULL)
  
  integrated.seurat <- Seurat::RunUMAP(object=integrated.seurat,
                                       dims=1:40,
                                       reduction="pca")
  return(integrated.seurat)
}

# DEPRECATED (used during Seurat v3)
### Generate whitelist for CITESeq
# Input is filtered seurat object
# Output is a list of csv files - one per batch containing valid barcodes
v3_generate_whitelist <- function(filtered.seurat, output_path){
  
  # Extract barcodes and split by "_"
  bc <- filtered.seurat@meta.data$Cell
  
  # Adjust this based on how your samples are named
  # NOTE: There will be multiple samples within each batch
  barcodes <- data.frame(stringr::str_split_fixed(bc, "_", 2)) %>%
    dplyr::rename(Batch = identity(1), Barcodes = identity(2)) %>%
    dplyr::mutate(Barcodes = stringr::str_replace(Barcodes, "-1", ""),
                  Batch = gsub(pattern="-GEX.*", replacement="", x=Batch))
  
  # Remove duplicate barcodes within each batch
  barcodes <- barcodes %>% 
    dplyr::group_by(Batch) %>% 
    dplyr::distinct_at("Barcodes", .keep_all=TRUE) %>% 
    as.data.frame()
  
  # Check how many barcodes are present in each batch
  barcodes %>% dplyr::group_by(Batch) %>% dplyr::count()
  
  # Save barcodes from each batch to individual csv files
  for (i in unique(barcodes$Batch)){
    whitelist <- barcodes %>%
      dplyr::filter(Batch == i) %>%
      dplyr::select(Barcodes)
    
    write.table(x = whitelist,
                file = paste0(scripts_path, proj, "_", i, "_whitelist.csv"),
                row.names = FALSE,
                col.names = FALSE)
  }
}

# DEPRECATED (used during Seurat v3)
# Input is path to folder containing h5ad
# Output is a raw seurat object
v3_read_h5ad <- function(input_path){
  
  # Load h5ad (useful if analyzing collaborator data in h5ad format)
  SeuratDisk::Convert(source = paste0(input_path, proj, ".h5ad"),
                      dest = "h5seurat",
                      assay="RNA",
                      overwrite = FALSE)
  
  raw_seurat <- SeuratDisk::LoadH5Seurat(file = paste0(input_path, proj, ".h5seurat"))
  
  return(raw_seurat)
}

### Import data from output of citeseq
# Input is path to demux_results folder
# Ouput is seurat object of each sample
v3_read_citeseq <- function(input_path){
  
  # Create a list of samples that have been demultiplexed already
  files <- list.files(path = paste0(input_path, "singlets/"),
                      full.names = FALSE)
  samples <- gsub(pattern="\\..*", replacement="", x=files)
  
  # Loop through each of the individual object in demux directory & import data
  for (i in 1:length(files)){
    
    # Read the seurat object containing demultiplexed singlets
    sample.seurat <- readRDS(file = paste0(demux_results, "singlets/", files[i]))
    
    # Assign the seurat object to its corresponding variable
    assign(samples[i], sample.seurat)
  }
}

SpatialFeaturePlotBlend <- function(cells_obj, column_1, column_2, combine=TRUE){
  
  # Convert decimal number to hexadecimal. Pad with 0s if only a single
  # character following conversion.
  as_hex <- function(num) {
    hex_str <- as.character(as.hexmode(num))
    if (nchar(hex_str) == 1) {
      hex_str <- paste0("0", hex_str)
    }
    
    return(hex_str)
  }
  
  metadata_to_hexadecimal <- function(in_dat) {
    apply(in_dat, 2,
          function(x) {
            # Make minimum 0
            x - min(x)
          }) %>%
      apply(2,
            function(x) {
              # Constrain to range [0, 255]
              round(255 * (x / max(x)))
            }) %>%
      apply(1,
            function(x) {
              # Convert to hexadecimal codes
              toupper(paste0("#", as_hex(x[1]), as_hex(x[2]), "00"))
            })
  }
  
  blend_plot_theme <- theme(legend.position="none",
                            plot.title=element_text(hjust=0.5))
  
  plot_list <- lapply(c(column_1, column_2),
                      function(column) {
                        max_color <- if_else(column == column_1,
                                             "#FF0000", "#00FF00")
                        SpatialFeaturePlot(cells_obj, column) +
                          scale_fill_gradient(low="#000000",
                                              high=max_color) +
                          ggtitle(column) +
                          blend_plot_theme
                      })
  
  dat <- FetchData(cells_obj, c(column_1, column_2))
  colors <- as.matrix(dat) %>% metadata_to_hexadecimal()
  
  new_md_column <- paste0(column_1, "_vs_", column_2)
  cells_obj[[new_md_column]] <- colors
  names(colors) <- as.character(colors)
  
  plot_list[[3]] <- SpatialDimPlot(cells_obj, new_md_column, cols=colors) +
    ggtitle(paste0(column_1, "_", column_2)) +
    blend_plot_theme
  
  side_length <- 100
  legend_grid <- expand.grid(seq(from=min(dat[, column_1]),
                                 to=max(dat[, column_1]),
                                 length.out=side_length),
                             seq(from=min(dat[, column_2]),
                                 to=max(dat[, column_2]),
                                 length.out=side_length))
  colnames(legend_grid) <- c(column_1, column_2)
  legend_colors <- metadata_to_hexadecimal(legend_grid)
  legend_grid$color <- legend_colors
  names(legend_colors) <- legend_colors
  
  legend <- ggplot(legend_grid,
                   aes(x=.data[[column_1]], y=.data[[column_2]],
                       color=color)) +
    geom_point(shape=15, size=1.9) +
    scale_color_manual(values=legend_colors) +
    coord_cartesian(expand=FALSE) +
    theme(legend.position="none", aspect.ratio=1,
          panel.background=element_blank())
  
  plot_list[[4]] <- wrap_plots(ggplot() + theme_void(), legend,
                               ggplot() + theme_void(), ncol=1,
                               heights=c(0.2, 0.6, 0.2))
  
  if (combine == FALSE) {
    return(plot_list)
  } else {
    p <- wrap_plots(plot_list, nrow=1,
                    widths=c(0.28, 0.28, 0.28, 0.16))
    return(p)
  }
}
