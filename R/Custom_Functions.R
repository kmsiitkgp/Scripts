# ---- 📦 LOAD PACKAGES ----

# NOTE: survminer handles %++% while dplyr handles %>%

pkgs <- c(
  "BiocManager", "remotes", "AnnotationHub", "ensembldb", "org.Hs.eg.db",
  "org.Mm.eg.db", "fgsea", "clusterProfiler", "DESeq2", "sva", "GSVA", 
  "RcisTarget", "glmGamPoi", "tximport", "Seurat", "harmony", "hdf5r", "scCustomize", 
  "reticulate", "ashr", "infercnv", "UCell", "scDblFinder", "DropletUtils", 
  "CellChat", "SeuratWrappers", "presto", "DoubletFinder", "SeuratData", 
  "oligo", "oligoData", "illuminaHumanv4.db", "hgu133plus2.db", "GEOquery", 
  "affy", "lumi", "openxlsx", "dplyr", "tibble", "stringr", "purrr", "ggplot2",
  "ggplotify", "ggrepel", "ggpubr", "ggfortify", "ggridges", "ggbeeswarm",
  "pheatmap", "VennDiagram", "survival", "survminer", "UpSetR", "umap", 
  "plot3D", "cowplot", "viridis", "RColorBrewer", "colorspace", 
  "enrichplot", "ComplexHeatmap", "NanoStringNCTools", "GeomxTools", 
  "GeoMxWorkflows", "networkD3", "httr", "decoupleR", "OmnipathR", "SeuratDisk",
  "clustree", "crayon", "uwot"
)

for (pkg in pkgs) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    message(paste("Loaded", pkg))
  } else {
    message(paste("Package", pkg, "is not installed — skipping"))
  }
}



# ---- 🎨 CUSTOM PALETTE &  THEME ----

# custom_theme <- ggplot2::theme(
#   plot.title    = element_text(hjust = 0.5),
#   legend.title  = element_text(hjust = 0,   vjust = 1, angle = 0),
#   axis.title.x  = element_text(hjust = 0.5, vjust = 0, angle = 0),
#   axis.title.y  = element_text(hjust = 0.5, vjust = 1, angle = 90),
#   legend.text   = element_text(hjust = 0.5),
#   axis.text.x   = element_text(hjust = 0.5, vjust = 0.5, angle = 45),
#   axis.text.y   = element_text(hjust = 0.5, vjust = 0.5, angle = 0)
# )

# NOTE: Text Justification and Spacing:
# hjust: Horizontal alignment (0=left, 0.5=center, 1=right).
# vjust: Vertical alignment (0=bottom, 0.5=middle, 1=top).
# margin: Controls the external padding (empty space) around the text element.
#         It takes arguments in the order T(op), R(ight), B(ottom), L(eft).

custom_theme <- ggplot2::theme_classic() + 
  ggplot2::theme(
    
    # Plot title
    plot.title    = ggplot2::element_text(family = "sans", face = "bold",  colour = "black", size = 16, hjust = 0.5),
    plot.subtitle = ggplot2::element_text(family = "sans", face = "plain", colour = "black", size = 12, hjust = 0.5, margin = ggplot2::margin(b = 5)),
    
    # Axis titles
    axis.title.x  = ggplot2::element_text(family = "sans", face = "bold",  colour = "black", size = 12, margin = ggplot2::margin(t = 10)),
    axis.title.y  = ggplot2::element_text(family = "sans", face = "bold",  colour = "black", size = 12, margin = ggplot2::margin(r = 10)), 
    
    # Axis text
    axis.text.x   = ggplot2::element_text(family = "sans", face = "plain", colour = "black", size = 10, angle = 45, hjust = 1),
    axis.text.y   = ggplot2::element_text(family = "sans", face = "plain", colour = "black", size = 10),
    
    # Legend styling
    legend.title  = ggplot2::element_text(family = "sans", face = "bold",  colour = "black", size = 12),
    legend.text   = ggplot2::element_text(family = "sans", face = "plain", colour = "black", size = 10),
    legend.position = "right", 
    legend.key = ggplot2::element_rect(fill = "white", colour = NA),
    
    # Additional spacing (optional but improves readability)
    plot.margin = ggplot2::margin(3, 3, 3, 3, unit = "mm"))

custom_palette <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173", "#5254a3", "#6b6ecf", "#333333", "#cedb9c", "#8ca252",
  "#a55194", "#e5e56f", "#66a61e", "#e6ab02", "#a6761d", "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ffff33",
  "#f781bf", "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3", "#1e90ff",
  "#ff4500", "#32cd32", "#ff0000", "#8a2be2", "#a0522d", "#ff66cc", "#9c9ede", "#adff2f", "#00ced1", "#ffd700",
  "#6699cc", "#cc6644", "#66aa66", "#cc6666", "#9966cc", "#996633", "#cc99cc", "#99cc44", "#66cccc", "#cccc66",
  "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
  "#ffffcc", "#e7ba52", "#ce6dbd", "#d6616b", "#b5cf6b", "#dbdb5c", "#e7cb94", "#ad494a", "#bd9e39", "#de9ed6",
  "#e7969c", "#33a02c", "#b2df8a", "#fdbf6f", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#8dd3c7", "#ffffb3",
  "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f",
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#999999"
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

# ---- 📝 LOGGING RELATED FUNCTIONS ----

# Suppress warnings and messages
quiet_msg <- function(expr) {
  # create a temp file and open a connection
  tmp <- tempfile()
  con <- file(tmp, open = "wt")
  
  # sink both output and message streams to the same connection
  sink(con, type = "output")
  sink(con, type = "message")
  
  result <- NULL
  tryCatch({
    result <- expr
  }, finally = {
    # restore normal output in the correct order
    sink(type = "message")
    sink(type = "output")
    close(con)
  })
  
  return(result)
}

# Log info messages (green)
log_info <- function(sample, step, msg) {
  sample <- sample %||% ""  # fallback to empty string if NULL
  prefix <- green(formatC("[INFO]", width = 7, flag = " "))
  message(glue::glue("{prefix} [{sample} | {toupper(step)}] {msg}"))
}

# Log warning messages (yellow)
log_warn <- function(sample, step, msg) {
  sample <- sample %||% ""  # fallback to empty string if NULL
  prefix <- yellow(formatC("[WARN]", width = 7, flag = " "))
  message(glue::glue("{prefix} [{sample} | {toupper(step)}] {msg}"))
}

# Log error messages (red)
log_error <- function(sample, step, msg) {
  sample <- sample %||% ""  # fallback to empty string if NULL
  prefix <- red(formatC("[ERROR]", width = 7, flag = " "))
  message(glue::glue("{prefix} [{sample} | {toupper(step)}] {msg}"))
  stop("Workflow stopped.", call. = FALSE)
}

# Optional: header for sample processing
log_sample_header <- function(sample) {
  cat(blue$bold(glue::glue("\n--- Processing Sample: {sample} ---\n\n")))
  
}

# Optional: section header
log_section <- function(section_name) {
  cat(magenta$bold(glue::glue("\n[{toupper(section_name)}]\n")))
}

# ---- BULK RNA SEQ ANALYSIS RELATED FUNCTIONS ----

# --- DESeq2 Design Notes ---
# The design formula specifies the variables to model in differential expression.
# 
# 1. When to merge columns (e.g., Batch + Condition → "Batch_Condition"):
#    - Use when some combinations of variables are missing (unbalanced design).
#    - Use when you want to directly compare specific groups (pairwise contrasts).
#    - Simplifies the model and avoids rank-deficiency errors.
#
# 2. When not to merge (keep separate columns, possibly with interactions):
#    - Use when all combinations exist (balanced design) and you want to estimate
#      the overall effect of a variable while adjusting for others (e.g., Age, Sex).
#    - Use when you want to explicitly model interactions (e.g., Sex:Condition).
#
# Example:
#   design = ~ Age + Sex + Condition        # keeps variables separate
#   design = ~ Batch_Condition              # merged groups for direct pairwise comparisons
#
# Notes:
# - Merging is safer for unbalanced datasets.
# - Interactions allow testing variable-specific effects, but contrasts become more complex.

# --- CONTRAST NAMING RESTRICTIONS FOR all.vars() ---
# To avoid errors in as.formula(), group names in 'contrast' MUST NOT:
# 1. Start with a number (e.g., '10_Control' is invalid)
# 2. Start with a dot followed by a number (e.g., '.2Group' is invalid)
# 3. Contain spaces (e.g., 'Wild Type' is invalid)
# 4. Contain special characters: ! @ # $ % ^ & * ( ) + = [ ] { } | \ ; : ' " < > , / ?
# 5. Use R Reserved Words: if, else, repeat, while, function, for, in, next, break, 
#    TRUE, FALSE, NULL, Inf, NaN, NA, NA_integer_, etc.
#
# VALID CHARACTERS: Only letters, numbers, dots (.), and underscores (_).
# MUST START WITH: A letter or a dot (if not followed by a number).

# # Setup project
# proj.params <- setup_project(
#   proj = proj,
#   species = species,
#   contrasts = contrasts,
#   parent_dir = parent_dir,
#   gmt_dir = gmt_dir,
#   deseq2.override = deseq2.override,
#   heatmap.override = heatmap.override,
#   volcano.override = volcano.override
# )

# Default Project Directories & Parameters
setup_project <- function(proj, species, contrasts,
                          parent_dir, gmt_dir, scripts_dir = NULL,
                          deseq2.override  = list(), heatmap.override = list(),
                          volcano.override = list(), pathway.override = list()) {
  
  # ---- 🛠️ Global Environment Configuration ----
  
  options(future.globals.maxSize = 1e15)            # Increase future globals size for parallelization
  options(Seurat.object.assay.version = "v5")       # Ensure Seurat uses v5 assay format
  set.seed(1234)                                    # Set seed for reproducibility
  
  # ---- 📊 Default Parameter Definition ----
  
  # Default DESeq2 Parameters
  default.deseq2 <- list(
    contrasts     = c("Treatment-Control"),   # Vector of contrasts for DE analysis
    design        = "Comparisons",            # DESeq2 design formula or column name
    lfc_cutoff    = 0,                        # Log fold change cutoff for significance
    padj_cutoff   = 0.1,                      # Adjusted p-value cutoff for significance
    batch_correct = FALSE                     # Boolean, whether to apply batch correction
  )
  
  # Default Heatmap Parameters
  default.heatmap <- list(
    label_genes        = NULL,        # NULL, 1 or more genes to be labelled
    filename           = NULL,        # Suffix added to filename
    output_dir         = NULL,        # Location to save plot and matrix
    metadata_col       = NULL,        # Dataframe with samples as rownames and columns for annotation
    metadata_row       = NULL,        # Dataframe with genes as rownames and columns for annotation
    col_annotations    = NULL,        # NULL, 1 or more columns from metadata_col for column annotation
    row_annotations    = NULL,        # NULL, 1 or more columns from metadata_row for row annotation
    col_gap_by         = NULL,        # NULL, 1 column from metadata_col to define column gaps in heatmap
    row_gap_by         = NULL,        # NULL, 1 column from metadata_row to define row gaps in heatmap
    col_cluster_by     = "all",       # NULL, 1 column from metadata_col, "all", "alphabetical" for clustering columns
    row_cluster_by     = "all",       # NULL, 1 column from metadata_row, "all", "alphabetical" for clustering columns
    plot_title         = NULL,        # NULL, Title for heatmap (default NULL i.e. no title)
    heatmap_palette    = "rdbu",      # Color palette for heatmap matrix ("rdbu" or "vrds")
    annotation_palette = "discrete",  # Color palette for heatmap annotation ("discrete" or "continuous")
    border_color       = NA,          # Color of heatmap cell borders (default NA i.e. no border)
    force_log          = FALSE,       # Force log transform (default FALSE i.e. auto detect)
    show_expr_legend   = TRUE,        # Show expression legend (set FALSE if annotations overlap)
    save_plot          = FALSE,        # Save the heatmap plot as pdf (default FALSE i.e. no save)
    save_matrix        = FALSE         # Save the heatmap matrix as xlsx (default FALSE i.e. no save)
  )
  
  # Default Volcano Parameters
  default.volcano <- list(
    filename     = NULL,               # Suffix added to filename
    contrast     = "Target-Reference", 
    label_genes  = NULL,               # NULL, 1 or more genes to be labelled
    top_n        = 5,                  # If label_genes = NULL, label top_n genes (default 5 genes)
    lfc_cutoff   = 0.58,               # Log fold change cutoff (default 0.58)
    padj_cutoff  = 0.05                # Adjusted p-value cutoff (default 0.05)
  )
  
  # Apply Overrides
  deseq2.params  <- modifyList(default.deseq2,  deseq2.override)
  heatmap.params <- modifyList(default.heatmap, heatmap.override)
  volcano.params <- modifyList(default.volcano, volcano.override)
  
  # ---- 📂 Directory Structure Setup ----
  
  proj_dir      <- file.path(parent_dir, proj)
  
  # Bulk RNA-seq Directories
  counts_dir    <- file.path(proj_dir, "counts")               # Directory containing count files
  salmon_dir    <- file.path(proj_dir, "salmon")               # Directory containing salmon quant.sf files
  contrast_dir  <- file.path(proj_dir, contrasts)              # Directory to store results for each contrast
  pathway_dir   <- file.path(contrast_dir, "Pathway_Analysis") # Directory to store Pathway analysis results
  tf_dir        <- file.path(contrast_dir, "TF_Analysis")      # Directory to store TF analysis results
  
  # Single-Cell RNA-seq Directories
  diagnostics_dir        <- file.path(proj_dir, "04.Diagnostics")
  demux_dir              <- file.path(proj_dir, "05.Demux")
  raw_matrix_dir         <- file.path(proj_dir, "06.Matrix", "raw_feature_bc_matrix")
  filt_matrix_dir        <- file.path(proj_dir, "06.Matrix", "filt_feature_bc_matrix")
  hto_matrix_dir         <- file.path(proj_dir, "06.Matrix", "HTO_bc_matrix")
  seurat_dir             <- file.path(proj_dir, "Seurat")
  sc_deseq2_dir          <- file.path(proj_dir, "DESeq2")
  pyscenic_dir           <- file.path(proj_dir, "pySCENIC")
  scvelo_dir             <- file.path(proj_dir, "scVelo")
  velocyto_dir           <- file.path(proj_dir, "velocyto")
  cellphonedb_dir        <- file.path(proj_dir, "CellphoneDB")
  cellchat_dir           <- file.path(proj_dir, "CellChat")
  
  # ---- 🧬 Cell Marker and Gene Set Extraction ----
  
  s_genes <- c()
  g2m_genes <- c()
  
  if (!is.null(scripts_dir)){ 
    metafile               <- file.path(scripts_dir, "scRNASeq_Metadata", paste0(proj, "_Metadata.xlsx"))
    markerfile             <- file.path(scripts_dir, "Cell_Type_Markers.xlsx")
    cell.cycle.marker.file <- file.path(scripts_dir, "Cell_Cycle_Markers.xlsx")
    
    # Extract Cell Cycle Genes (Human and Mouse)
    if (file.exists(cell.cycle.marker.file)) {
      cc_markers <- openxlsx::read.xlsx(cell.cycle.marker.file)
      
      # Combine Human and Mouse genes for each phase
      s_genes <- c(cc_markers$Human_Gene[cc_markers$Phase == "S"],
                   cc_markers$Mouse_Gene[cc_markers$Phase == "S"])
      
      g2m_genes <- c(cc_markers$Human_Gene[cc_markers$Phase == "G2/M"],
                     cc_markers$Mouse_Gene[cc_markers$Phase == "G2/M"])
      
      # Clean up potential NA/empty values resulting from unlist/subsetting
      s_genes <- s_genes[!is.na(s_genes) & s_genes != ""]
      g2m_genes <- g2m_genes[!is.na(g2m_genes) & g2m_genes != ""]
    }
  }
  
  # ---- 📦 Final Project Parameters Construction ----
  
  proj.params <- c(
    list(
      # Project info
      proj        = proj,
      species     = species,
      
      # Bulk RNA-seq directories
      proj_dir    = normalizePath(proj_dir,     mustWork = FALSE),
      counts_dir  = normalizePath(counts_dir,   mustWork = FALSE),
      salmon_dir  = normalizePath(salmon_dir,   mustWork = FALSE),
      gmt_dir     = normalizePath(gmt_dir,      mustWork = FALSE),
      contrast_dir= normalizePath(contrast_dir, mustWork = FALSE),
      pathway_dir = normalizePath(pathway_dir,  mustWork = FALSE),
      tf_dir      = normalizePath(tf_dir,       mustWork = FALSE),
      
      # Bulk RNA-seq parameters
      contrasts   = contrasts,
      deseq2      = deseq2.params,
      heatmap     = heatmap.params,
      volcano     = volcano.params,
      
      # Single cell RNA-seq directories
      diagnostics_dir  = normalizePath(diagnostics_dir, mustWork = FALSE),
      demux_dir        = normalizePath(demux_dir,       mustWork = FALSE),
      raw_matrix_dir   = normalizePath(raw_matrix_dir,  mustWork = FALSE), 
      filt_matrix_dir  = normalizePath(filt_matrix_dir, mustWork = FALSE), 
      hto_matrix_dir   = normalizePath(hto_matrix_dir,  mustWork = FALSE), 
      seurat_dir       = normalizePath(seurat_dir,      mustWork = FALSE), 
      sc_deseq2_dir    = normalizePath(sc_deseq2_dir,   mustWork = FALSE), 
      pyscenic_dir     = normalizePath(pyscenic_dir,    mustWork = FALSE), 
      scvelo_dir       = normalizePath(scvelo_dir,      mustWork = FALSE), 
      velocyto_dir     = normalizePath(velocyto_dir,    mustWork = FALSE), 
      cellphonedb_dir  = normalizePath(cellphonedb_dir, mustWork = FALSE), 
      cellchat_dir     = normalizePath(cellchat_dir,    mustWork = FALSE)
    )
  )
  
  # Conditionally add single cell data only if scripts_dir is not NULL
  if (!is.null(scripts_dir)) {
    proj.params$metafile   <- normalizePath(metafile,   mustWork = FALSE)
    proj.params$markerfile <- normalizePath(markerfile, mustWork = FALSE)
    proj.params$cell_cycle <- list(S = s_genes, G2M = g2m_genes)
  }
  
  # ---- ⚠️ Final Checks and Return ----
  
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

merge_counts <- function(counts_dir, filename = NULL, output_dir) {
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(counts_dir = counts_dir, filename = filename, output_dir = output_dir)
  
  # ---- 🔎 Identify Count Files ----
  
  # Pattern matches STAR (ReadsPerGene.out.tab) or standard HTSeq .txt files
  count_files <- list.files(path = counts_dir, 
                            pattern = "\\.txt$|ReadsPerGene\\.out\\.tab$", 
                            full.names = TRUE)
  
  if (length(count_files) == 0) {
    log_warn(sample = "", 
             step    = "merge_counts", 
             msg     = glue::glue("No count files found in: '{counts_dir}'.
                                 Provide raw counts as excel for analysis."))
    return(NULL)
  }
  
  # ---- 🧪 Initialize Containers ----
  
  all_counts <- list()
  gene_lists <- list()
  sample_ids <- character()
  
  # Define metadata rows to exclude (HTSeq/STAR stats)
  special_counters <- c("__no_feature", "__ambiguous", "__too_low_aQual", 
                        "__not_aligned", "__alignment_not_unique", "__assignment_counts",
                        "N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
  
  # ---- 🔄 Parse Files and Detect Strandedness ----
  
  for (count_file in count_files) {
    
    sample_id <- gsub("\\..*$|ReadsPerGene\\.out\\.tab", "", basename(count_file))
    
    # Read count file
    df <- tryCatch({
      read.table(file = count_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    }, error = function(e) {
      log_error(sample = sample_id, 
                step   = "merge_counts", 
                msg    = glue::glue("Error reading file: {e$message}"))
    })
    
    if (ncol(df) < 4) {
      log_error(sample = sample_id, 
                step = "merge_counts", 
                msg ="STAR/HTSeq count file does not have expected 4 columns.")
    }
    
    # Remove special counters/metadata rows
    df <- df %>% dplyr::filter(!(.[[1]] %in% special_counters))
    
    gene_ids <- df[[1]]
    strand_sums <- colSums(df[2:4], na.rm = TRUE)
    
    # Logic to determine strandedness based on column sums
    # Col 2: Unstranded 
    # Col 3: Forward 
    # Col 4: Reverse
    if (abs((strand_sums[1]/strand_sums[2]) - (strand_sums[1]/strand_sums[3])) < 2) {
      log_info(sample = sample_id, step = "merge_counts", msg ="Detected unstranded library.")
      counts <- df[[2]]
    } else if (strand_sums[2] > 3 * strand_sums[3]) {
      log_info(sample = sample_id, step = "merge_counts", msg ="Detected positively stranded library.")
      counts <- df[[3]]
    } else if (strand_sums[3] > 3 * strand_sums[2]) {
      log_info(sample = sample_id, step = "merge_counts", msg ="Detected negatively stranded library.")
      counts <- df[[4]]
    } else {
      log_error(sample = sample_id, 
                step = "merge_counts", 
                msg ="Could not confidently determine strandedness.")
    }
    
    all_counts[[sample_id]] <- counts
    gene_lists[[sample_id]] <- gene_ids
    sample_ids <- c(sample_ids, sample_id)
  }
  
  # ---- ⚖️ Check Gene Consistency ----
  
  ref_genes <- gene_lists[[1]]
  for (i in seq_along(gene_lists)) {
    if (!identical(ref_genes, gene_lists[[i]])) {
      log_error(sample = names(gene_lists)[i], 
                step = "merge_counts", 
                msg ="Gene ID mismatch. Ensure all samples were mapped to the same reference.")
    }
  }
  
  # ---- 📊 Build and Filter Count Matrix ----
  
  count_matrix <- do.call(cbind, all_counts)
  colnames(count_matrix) <- sample_ids
  count_matrix <- data.frame(SYMBOL = ref_genes, count_matrix, stringsAsFactors = FALSE)
  
  # Remove genes with 0 counts across all samples
  count_matrix <- count_matrix[rowSums(count_matrix[,-1]) > 0, , drop = FALSE]
  rownames(count_matrix) <- NULL    # Reset row names after subsetting so they are sequential (1, 2, 3, ...)
  
  # Remove samples with 0 counts across all genes
  count_matrix <- count_matrix[, c(TRUE, colSums(count_matrix[,-1]) > 0), drop = FALSE]
  
  # ---- 💾 Export to Excel ----
  
  file_name <- file.path(output_dir, paste0(filename, "_Raw_counts.xlsx"))
  
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb = wb, sheetName = "Raw_counts")
  openxlsx::writeData(wb = wb, sheet = "Raw_counts", x = count_matrix)
  openxlsx::saveWorkbook(wb = wb, file = file_name, overwrite = TRUE)
  
  # ---- 🪵 Log Output and Return Count Matrix ----
  
  log_info(sample = "", 
           step   = "merge_counts", 
           msg    = glue::glue("Merged counts successfully. Saved to: '{file_name}'"))
  
  return(invisible(count_matrix))
}

prep_txi <- function(salmon_dir, species, output_dir, 
                     db_version = NULL, filename = NULL){
  
  # Connect to AnnotationHub 
  hub <- AnnotationHub::AnnotationHub()
  
  # ---- 🔍 Query Database ----
  
  log_info(sample = species, 
           step   = "prep_txi", 
           msg    = "Fetching Ensembl Database...")
  
  hub_db <- AnnotationHub::query(x           = hub, 
                                 pattern     = c("EnsDb", species), 
                                 ignore.case = TRUE)
  
  # Glimpse of hub_db
  print(hub_db %>%
          mcols() %>%
          as.data.frame() %>%
          dplyr::select(title, species, genome, rdatadateadded, sourcetype))
  
  # Acquire the latest version available in the hub
  latest_id <- hub_db %>%
    mcols() %>%
    as.data.frame() %>%
    { 
      if (!is.null(db_version)) {
        filter(., grepl(db_version, .data$title))
      } else {
        .
      }
    } %>%
    dplyr::arrange(desc(rdatadateadded)) %>%
    head(n = 1) %>%
    rownames()
  
  if (length(latest_id) == 0) {
    log_error(sample = species, 
              step   = "prep_txi", 
              msg    = "Could not find a valid EnsDb in AnnotationHub.")
  }
  
  # Download the appropriate Ensembldb database
  ensdb <- hub_db[[latest_id]]
  
  # Extract transcript and gene info
  tx2gene <- GenomicFeatures::transcripts(x = ensdb) %>%
    as.data.frame() %>%
    dplyr::select("tx_id", "gene_id") %>%
    data.frame()
  
  
  # Get the salmon files
  quant_files <- list.files(path       = salmon_dir, 
                            pattern    = "quant.sf$", 
                            recursive  = TRUE,
                            full.names = TRUE)
  
  # Name the files using sample IDs
  sample_names <- list.files(path       = salmon_dir, 
                             pattern    = "quant.sf$", 
                             recursive  = TRUE,
                             full.names = FALSE)
  sample_names <- gsub(pattern = "\\.quant.sf$", replacement = "", x = sample_names)
  names(quant_files) <- make.names(sample_names)
  
  txi <- tximport::tximport(files           = quant_files, 
                            type            = "salmon", 
                            tx2gene         = tx2gene,  
                            ignoreTxVersion = TRUE)
  
  # txi is your tximport object
  file_name <- paste("txi", filename, sep = "_")
  file_extension <- ".rds"
  saveRDS(txi, file = file.path(output_dir, paste0(file_name, file_extension)))
  
  # Save the Raw Counts (The "Must-Have" for GEO)
  txi$counts %>%
    as.data.frame() %>%
    rownames_to_column(var = "GeneID") %>%
    write.table(file.path(output_dir, "Processed_Raw_Gene_Counts.txt"), 
                sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Save the TPMs (The "Nice-to-Have" for Heatmaps/Reviewers)
  txi$abundance %>%
    as.data.frame() %>%
    rownames_to_column(var = "GeneID") %>%
    write.table(file.path(output_dir, "Processed_TPM_Values.txt"), 
                sep = "\t", quote = FALSE, row.names = FALSE)
  
  return(txi)
}

prepare_deseq2_input <- function(expr_mat, txi, metadata, design) {
  
  if (is.null(expr_mat) & !is.null(txi)){
    expr_mat <- txi$counts
  }
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(expr_mat = expr_mat, metadata = metadata)
  
  if (!"Sample_ID" %in% colnames(metadata)) {
    log_error(sample = "",
              step   = "prepare_deseq2_input",
              msg    = glue::glue("`metadata` MUST contain 'Sample_ID' column."))
  }
  
  if (base::anyDuplicated(metadata$Sample_ID)) {
    log_error(sample = "", 
              step   = "prepare_deseq2_input", 
              msg    = "Duplicate samples detected in 'Sample_ID' column in metadata.")
  }
  
  if (base::anyDuplicated(rownames(expr_mat))) {
    log_error(sample = "", 
              step   = "prepare_deseq2_input", 
              msg    = "Duplicate gene symbols detected in row names of count matrix.")
  }
  
  if (is.null(design)) {
    log_error(sample = "", 
              step   = "prepare_deseq2_input", 
              msg    = "DESeq2 design formula is NULL. Please fix it.")
  }
  
  # Alert the user if samples in expr_mat and metadata differ
  clean_meta_names <- make.names(metadata$Sample_ID)
  clean_mat_names  <- make.names(colnames(expr_mat))
  common <- intersect(clean_meta_names, clean_mat_names)
  
  if (length(common) < length(clean_meta_names)) {
    lost <- length(clean_meta_names) - length(common)
    log_warn(sample = "",
             step   = "prepare_deseq2_input",
             msg    = glue::glue("Alignment lost {lost} samples. Check for naming mismatches (e.g. '-' vs '_')."))
  }
  
  
  # ⚠️ Detect potential gene ID column in expr_mat
  # Compute column sums
  col_sums <- colSums(expr_mat, na.rm = TRUE)
  first_col_name <- colnames(expr_mat)[1]
  first_col <- expr_mat[,1]
  
  # Heuristic checks
  is_numeric_outlier <- is.numeric(first_col) && col_sums[1] > 10 * median(col_sums[-1])
  is_name_not_in_meta <- !(first_col_name %in% make.names(metadata$Sample_ID))
  
  if (is_numeric_outlier || is_name_not_in_meta) {
    log_error(sample = "",
              step    = "prepare_deseq2_input",
              msg     = glue::glue("First column '{first_col_name}' looks like gene IDs (Entrez or SYMBOL), not sample counts."))
  }
  
  # ---- 📝️ Metadata Preparation & Design Validation ----
  
  log_info(sample = "", 
           step   = "prepare_deseq2_input", 
           msg    = glue::glue("Initial counts: {nrow(expr_mat)} genes x {nrow(metadata)} samples."))
  
  # Remove "sizeFactor" column if present and rows with missing sample names
  metadata <- metadata %>%
    dplyr::select(-any_of("sizeFactor")) %>%
    dplyr::filter(!is.na(Sample_ID))
  
  # Assign "Sample_ID" column as row names and convert to valid R names
  rownames(metadata) <- make.names(metadata$Sample_ID)
  
  # Add Batch column if absent
  if (!"Batch" %in% colnames(metadata)) {
    metadata$Batch <- 1
    log_warn(sample = "", 
             step   = "prepare_deseq2_input", 
             msg    = "No 'Batch' column found. Assigning default value '1'.")
  }
  
  # Extract variables from formula (e.g., ~Condition + Batch -> c("Condition", "Batch"))
  design_formula <- if (grepl("^~", design)) {
    as.formula(design)
  } else {
    as.formula(paste0("~", design))
  }
  design_vars <- all.vars(design_formula)
  missing_vars <- setdiff(design_vars, colnames(metadata))
  
  if (length(missing_vars) > 0) {
    log_error(sample = "", 
              step   = "prepare_deseq2_input", 
              msg    = glue::glue("Design variable(s) not found in metadata: '{paste(missing_vars, collapse = ', ')}'"))
  }
  
  # Identify samples with NA in the design columns
  na_idx <- apply(X = metadata[, design_vars, drop = FALSE],
                  MARGIN = 1, 
                  FUN = function(x) any(is.na(x)))
  
  # Only subset metadata if there are any NAs
  if (any(na_idx)) {
    na_samples <- rownames(metadata)[na_idx]
    log_warn(sample = "",
             step   = "prepare_deseq2_input",
             msg    = glue::glue("Removing {length(na_samples)} samples with NA in design columns: {paste(na_samples, collapse = ', ')}"))
    metadata <- metadata[!na_idx, , drop = FALSE]
  }
  
  # ---- 🧮 Expression Matrix Cleaning & Filtering ----
  
  # Convert column names to valid R names
  colnames(expr_mat) <- make.names(colnames(expr_mat))
  
  # Retain ONLY samples in metadata for accurate zero count calculations
  expr_mat <- expr_mat[, intersect(colnames(expr_mat), rownames(metadata)), drop = FALSE]
  
  # Remove rows with NA gene names
  na_rows <- is.na(rownames(expr_mat))
  expr_mat <- expr_mat[!na_rows, , drop = FALSE]
  log_info(sample = "",
           step   = "prepare_deseq2_input",
           msg    = glue::glue("Removed {sum(na_rows)} rows with missing gene names."))
  
  # Replace missing counts with 0
  expr_mat[is.na(expr_mat)] <- 0
  
  # Convert all columns to numeric safely
  expr_mat <- matrix(as.numeric(as.matrix(expr_mat)),
                     nrow = nrow(expr_mat),
                     ncol = ncol(expr_mat),
                     dimnames = dimnames(expr_mat))
  
  if (any(is.na(expr_mat))) {
    log_error(sample = "",
              step   = "plot_heatmap",
              msg    = "`expr_mat` contains non-numeric values that could not be converted.")
  }
  
  # Remove genes with zero counts across all samples
  zero_genes   <- rownames(expr_mat)[which(rowSums(expr_mat) == 0)]
  expr_mat <- expr_mat[rowSums(expr_mat) != 0, , drop = FALSE]
  log_warn(sample = "",
           step   = "prepare_deseq2_input",
           msg    = glue::glue("Removed {length(zero_genes)} genes with zero counts across all samples."))
  
  # Remove samples with zero total reads
  zero_samples <- colnames(expr_mat)[which(colSums(expr_mat) == 0)]
  expr_mat <- expr_mat[, colSums(expr_mat) != 0, drop = FALSE]
  log_warn(sample = "",
           step   = "prepare_deseq2_input",
           msg    = glue::glue("Removed {length(zero_samples)} samples with zero total counts."))
  
  
  # ---- 🧩 Final Data Structuring for DESeq2 ----
  
  # Synchronize samples between expr_mat and metadata
  
  common_samples <- intersect(colnames(expr_mat), rownames(metadata))
  expr_mat <- expr_mat[              , common_samples, drop = FALSE]
  metadata <- metadata[common_samples,               , drop = FALSE]
  
  if (length(common_samples) == 0) {
    log_error(sample = "",
              step   = "prepare_deseq2_input",
              msg    = "No overlapping samples between expr_mat and metadata.")
  }
  
  # Convert all metadata to factors for DESeq2 compatibility
  message("Structure of metadata before conversion:")
  str(metadata)
  
  #metadata <- as.data.frame(unclass(metadata), stringsAsFactors = TRUE)
  metadata[] <- lapply(metadata, as.factor)
  
  message("Structure of metadata after conversion:")
  str(metadata)
  
  # Reset expr_mat to null
  if (!is.null(txi)){
    expr_mat <- NULL
  }
  
  # ---- 🪵 Log Output and Return Cleaned Data ----
  
  if (!is.null(expr_mat)) {
    n_genes   <- nrow(expr_mat)
    n_samples <- ncol(expr_mat)
  } else if (!is.null(txi) && !is.null(txi$counts)) {
    n_genes   <- nrow(txi$counts)
    n_samples <- ncol(txi$counts)
    
  } else {
    log_error(sample = "",
              step   = "prepare_deseq2_input",
              msg    = "Both expr_mat and txi$counts are NULL. Cannot determine genes or samples.")
  }
  
  log_info(sample = "", 
           step   = "prepare_deseq2_input", 
           msg    = glue::glue("Successfully prepared DESeq2 input: {n_genes} genes and {n_samples} samples."))
  
  return(invisible(list(metadata = metadata, expr_mat = expr_mat, txi = txi)))
}

fit_deseq2_model <- function(expr_mat, txi, metadata, design) {
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(expr_mat = expr_mat, metadata = metadata)
  
  # ---- 🧪 DESeq2 Object Preparation & Filtering ----
  
  # Standardize DESeq2 design formula
  if (inherits(design, "formula")) {
    design_formula <- design
    
  } else if (is.numeric(design) && design == 1) {
    # Intercept-only model
    design_formula <- stats::as.formula("~ 1")
    
  } else if (is.character(design) && length(design) == 1 && nzchar(design)) {
    # nzchar = “non-zero number of characters”
    # Remove leading ~ and whitespace, then standardize
    clean_design <- sub("^\\s*~\\s*", "", design)
    design_formula <- stats::as.formula(paste0("~", clean_design))
    
  } else {
    log_error(sample = "",
              step   = "prepare_deseq2_input",
              msg    = "`design` must be a formula, the number 1, or a non-empty character string.")
  }
  
  # Prepare DESeq2 object
  if (!is.null(expr_mat)){
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = expr_mat,
                                          colData   = metadata,
                                          design    = design_formula)
  } else if (!is.null(txi)){
    dds <- DESeq2::DESeqDataSetFromTximport(txi     = txi,
                                            colData = metadata,
                                            design  = design_formula)
  }
  
  # Auto-detect for 'poscounts' if many zeros exist
  if (!is.null(expr_mat)){
    is_scRNA <- all(rowSums(expr_mat == 0) > 0)
  } else if (!is.null(txi$counts)){
    is_scRNA <- all(rowSums(txi$counts == 0) > 0)
  }
  
  if (is_scRNA) {
    log_warn(sample = "",
             step   = "run_deseq2",
             msg    = "scRNA-seq–like sparsity detected. Estimating size factors using 'poscounts'.")
    dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
  } else {
    # Pre-filter lowly expressed genes to improve sizefactor estimation in next step
    keep <- rowSums(DESeq2::counts(dds)) >= 10
    dds <- dds[keep, ]
  }
  
  # ---- 📉 Fit Selection (Parametric vs. Local) ----
  
  # We fit both to determine which dispersion model better suits the data distribution
  dds_para <- DESeq2::DESeq(object = dds, 
                            test = "Wald", 
                            fitType = "parametric", 
                            betaPrior = FALSE, 
                            minReplicatesForReplace = 7, 
                            quiet = TRUE)
  
  dds_local <- DESeq2::DESeq(object = dds, 
                             test = "Wald", 
                             fitType = "local", 
                             betaPrior = FALSE, 
                             minReplicatesForReplace = 7, 
                             quiet = TRUE)
  
  residual_para  <- mcols(dds_para)$dispGeneEst - mcols(dds_para)$dispFit
  residual_local <- mcols(dds_local)$dispGeneEst - mcols(dds_local)$dispFit
  
  if (median(residual_para^2, na.rm = TRUE) <= median(residual_local^2, na.rm = TRUE)) {
    fit <- "Parametric"
    dds <- dds_para
  } else {
    fit <- "Local"
    dds <- dds_local
  }
  
  # ---- 🪵 Log Output and Return DESeq2 object ----
  
  log_info(sample = "", 
           step   = "fit_deseq2_model", 
           msg    = glue::glue("DESeq2 model built using '{fit}' fit."))
  
  return(invisible(dds))
}
  
get_deseq2_results <- function(dds, contrast, output_dir, 
                               lfc_cutoff = 0, padj_cutoff = 0.1){
  
  # For ashr
  set.seed(1234)
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(dds = dds, output_dir = output_dir)
  
  # ---- 🧮 Dynamic Contrast Parsing ----
  
  # Extract design variables (e.g., c("Condition", "Batch"))
  # all.vars is robust; it handles ~Condition + Batch and ~Condition*Batch
  design_vars <- all.vars(DESeq2::design(dds))
  
  # Extract the model matrix 
  # This is useful for downstream checking of rank deficiency
  mod_mat <- stats::model.matrix(DESeq2::design(dds), data = SummarizedExperiment::colData(dds))
  
  # Create "Groups" column for easy subsetting/contrasts
  if (length(design_vars) > 0) {
    df_groups <- SummarizedExperiment::colData(dds) %>% 
      as.data.frame() %>% 
      tidyr::unite(col = "Groups", dplyr::all_of(design_vars), sep = ".", remove = FALSE)
    
    # Update the dds colData to include this combined Groups column
    SummarizedExperiment::colData(dds)$Groups <- as.factor(df_groups$Groups)
  } else {
    # If design is ~1, Groups is just "All"
    SummarizedExperiment::colData(dds)$Groups <- as.factor("All")
  }
  
  # Define all possible groups that could be compared based on the design
  groups <- unique(df_groups$Groups)
  
  # Map groups to mean model coefficients (Get all possible coefficient vectors)
  group_coef_list <- lapply(groups, function(i) colMeans(as.matrix(mod_mat[df_groups$Groups == i, , drop = FALSE])))
  names(group_coef_list) <- groups
  
  # Recursive function to evaluate the contrast string (e.g., "A-B") as a vector operation
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
  
  parsed_expr <- base::parse(text = contrast)[[1]]
  expr_sub <- replace_symbols(parsed_expr)
  contrast_vec <- base::eval(expr_sub)
  
  # ---- ⚠️ Sanity Check: contrast_vec Alignment ----
  
  if (length(contrast_vec) != ncol(mod_mat)) {
    log_error(sample = contrast,
              step   = "run_deseq2",
              msg    = glue::glue("Length of contrast_vec ({length(contrast_vec)}) does 
                                   not match the number of columns in the design matrix ({ncol(mod_mat)})."))
  }
  
  if (all(contrast_vec == 0)) {
    log_error(sample = contrast,
              step   = "run_deseq2",
              msg    = "contrast_vec is all zeros — invalid contrast. Check the user-provided contrast string.")
  }
  
  # ---- 🧬 Results & LFC Shrinkage ----
  
  res <- DESeq2::results(object               = dds, 
                         contrast             = contrast_vec,
                         lfcThreshold         = lfc_cutoff,
                         altHypothesis        = "greaterAbs",
                         cooksCutoff          = TRUE,
                         independentFiltering = TRUE,
                         alpha                = padj_cutoff,
                         pAdjustMethod        = "BH")
  
  # Use 'ashr' for shrinkage as it is robust for varied effect sizes
  res <- DESeq2::lfcShrink(dds = dds, res = res, type = "ashr", quiet = TRUE)
  summary(res)
  
  # ---- 📊 Data Formatting & Export ----
  
  # Consolidate duplicate symbols if any (using lowest padj value)
  DEGs_df <- res %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID") %>%
    add_annotation() %>% 
    dplyr::filter(!is.na(SYMBOL)) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::slice_min(order_by = padj, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(padj = case_when(is.na(padj) ~ 1,
                                   padj == 0   ~ min(padj[padj > 0]),
                                   TRUE        ~ padj))
  
  save_xlsx(DEGs_df, file.path(output_dir, "DEGs.xlsx"), "DEGs", row_names = FALSE)
  
  # ---- 🪵 Log Output and Return DESeq2 Results ----
  
  log_info(sample = contrast, 
           step   = "get_deseq2_results", 
           msg    = glue::glue("DESeq2 complete. Total genes: {nrow(dds)}. Contrast: {contrast}"))
  
  return(invisible(list(degs = DEGs_df, dds = dds)))
}

plot_ma <- function(dds, output_dir, filename = NULL) {
  
  # For ggrepel
  set.seed(1234)
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(dds = dds, output_dir = output_dir, filename = filename)
  
  # ---- 🖼️ Generate Plots ----
  
  # Prepare the data with "squished" values
  res_df <- as.data.frame(DESeq2::results(dds)) %>%
    dplyr::mutate(significant = padj < 0.1 & !is.na(padj),
                  log2FoldChange_capped = dplyr::case_when(log2FoldChange > 5  ~ 5,
                                                           log2FoldChange < -5 ~ -5,
                                                           TRUE                ~ log2FoldChange),
                  # Assign shapes: 16 is a solid circle, 17 is a solid triangle
                  is_outlier = log2FoldChange > 5 | log2FoldChange < -5,
                  point_shape = ifelse(is_outlier, 17, 16))
  
  p <- ggplot(data = res_df, 
              mapping = aes(x = baseMean, y = log2FoldChange, color = significant)) +
    geom_point(alpha = 0.5, size = 1.25) +
    theme_classic() +
    custom_theme +
    theme(legend.position = "none") +
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_color_manual(values = c("grey70", "firebrick3")) + 
    labs(title = "MA Plot",
         subtitle = "Red points indicate FDR < 0.1",
         x = "Mean of Normalized Counts",
         y = "Log2 Fold Change")
  
  capped_p <- ggplot(data = res_df, 
                     mapping = aes(x = baseMean, y = log2FoldChange_capped, color = significant)) +
    geom_point(alpha = 0.5, size = 1.25, aes(shape = point_shape)) +
    scale_shape_identity() + # Tells ggplot to use the actual numeric shape codes
    theme_classic() +
    custom_theme +
    theme(legend.position = "none") +
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_color_manual(values = c("grey70", "firebrick3")) +
    coord_cartesian(ylim = c(-5, 5)) + 
    labs(title = "MA Plot (with Capped Fold Changes)",
         subtitle = "Red points indicate FDR < 0.1",
         x = "Mean of Normalized Counts",
         y = "Log2 Fold Change")
  
  # ---- 💾 Save Plots ----
  
  file_extension <- ".pdf"
  file_name <- file.path(output_dir, paste0("MA_Plot_", filename, file_extension))
  
  # Open multi-page PDF
  grDevices::pdf(file = file_name, width = 8, height = 11.5, onefile = TRUE)  
  
  # PAGE 1: Traditional (DESeq2)
  # This function draws directly to the open PDF device
  DESeq2::plotMA(object = dds,
                 alpha  = 0.1, 
                 main   = "MA Plot (Traditional DESeq2)",
                 xlab   = "Mean of Normalized Counts",
                 MLE    = FALSE)
  
  print(capped_p)    # PAGE 2: Capped View (ggplot2)
  print(p)           # PAGE 3: Full View (ggplot2)
  
  dev.off() 
  
  # ---- 🪵 Log Output and Return ----
  
  log_info(sample = "",
           step   = "plot_ma",
           msg    = glue::glue("MA plot saved successfully to : '{file_name}'."))
  
  return(invisible(NULL))
}

plot_dispersion <- function(dds, output_dir, filename = NULL) {
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(dds = dds, output_dir = output_dir, filename = filename)
  
  # ---- 💾 Save Plots ----
  
  file_extension <- ".pdf"
  file_name <- file.path(output_dir, paste0("Dispersion_Plot_", filename, file_extension))
  
  # Open multi-page PDF
  grDevices::pdf(file = file_name, width = 8, height = 11.5, onefile = TRUE)  
  
  # PAGE 1: Traditional (DESeq2)
  # This function draws directly to the open PDF device
  DESeq2::plotDispEsts(object   = dds,
                       genecol  = "black",
                       fitcol   = "red",
                       finalcol = "dodgerblue",
                       legend   = TRUE,
                       xlab     = "Mean of Normalized Counts",
                       ylab     = "Dispersion",
                       log      = "xy",
                       cex      = 0.45)
  
  dev.off() 
  
  # ---- 🪵 Log Output and Return ----
  
  log_info(sample = "",
           step   = "plot_ma",
           msg    = glue::glue("Dispersion plot saved successfully to : '{file_name}'."))
  log_info(sample = "",
           step   = "plot_ma",
           msg    = "Expected results: Higher the mean, lower the dispersion")
  
  return(invisible(NULL))
}

plot_volcano <- function(res_df, output_dir, filename = NULL, 
                         contrast    = "Target-Reference", 
                         label_genes = NULL,
                         top_n       = 5,
                         lfc_cutoff  = 0.58, 
                         padj_cutoff = 0.05) {
  
  # For ggrepel
  set.seed(1234)
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(res_df = res_df, filename = filename, output_dir = output_dir)
  
  # ---- 🧪 Setup Volcano Parameters ----
  
  target    <- stringr::str_split(string = contrast, pattern = "-")[[1]][1]
  reference <- stringr::str_split(string = contrast, pattern = "-")[[1]][2]
  
  # ---- 🔄 Data Formatting & Relevance Scoring ----
  
  res_df <- res_df %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::mutate(Direction = dplyr::case_when(padj < padj_cutoff & log2FoldChange > lfc_cutoff  ~ paste0("Up in ", target),
                                               padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ paste0("Up in ", reference),
                                               TRUE ~ "Not Significant"),
                  # Handle padj = 0 for log scaling
                  padj = dplyr::case_when(padj == 0 ~ min(padj[padj > 0], na.rm = TRUE), 
                                          TRUE ~ padj),
                  Significance = dplyr::case_when(abs(log2FoldChange) >= lfc_cutoff & padj <= 0.001 ~ "FDR < 0.001",
                                                  abs(log2FoldChange) >= lfc_cutoff & padj <= 0.01  ~ "FDR < 0.01",
                                                  abs(log2FoldChange) >= lfc_cutoff & padj <= 0.05  ~ "FDR < 0.05",
                                                  TRUE ~ "Not Significant"),
                  Relevance = {
                    r <- abs(log2FoldChange) * -log10(padj)
                    pmin(r, quantile(r, 0.99, na.rm = TRUE)) # Cap outliers for visualization
                  })
  
  # ---- 🎨 Aesthetics & Palettes ----
  
  volcano_palette <- c(setNames(viridis::viridis(100)[1], paste0("Up in ", target)),
                       setNames(viridis::viridis(100)[50], paste0("Up in ", reference)),
                       #setNames(viridis::viridis(100)[100], "Not Significant"),
                       "Not Significant" = "grey80")
  
  alpha_palette <- c("FDR < 0.001" = 1, "FDR < 0.01" = 0.8, 
                     "FDR < 0.05" = 0.6, "Not Significant" = 0.4)
  
  # ---- 📐 Axis limits and breaks ----
  
  # Calculate dynamic axis limits
  x_vals <- res_df$log2FoldChange
  y_vals <- -log10(res_df$padj)
  
  # Keep only finite values
  x_vals <- x_vals[is.finite(x_vals)]
  y_vals <- y_vals[is.finite(y_vals)]
  
  x_max <- ceiling(max(abs(res_df$log2FoldChange), na.rm = TRUE))
  y_max <- ceiling(max(-log10(res_df$padj), na.rm = TRUE))
  
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
  
  # ---- 🖼️ Generate Plots ----
  
  p <- ggplot2::ggplot(data = res_df, 
                       mapping = aes(x = log2FoldChange, y = -log10(padj),
                                     color = Direction, 
                                     alpha = Significance,
                                     size = Relevance)) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.05, height = 0.05)) +
    ggplot2::theme_classic() +
    custom_theme + 
    ggplot2::labs(x = expression("log"[2]*"FC"),
                  y = expression("-log"[10]*"padj"),
                  fill = "Direction",
                  title = contrast) +
    ggplot2::geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dotted", color = "black", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = -log10(padj_cutoff), linetype = "dotted", color = "black", linewidth = 0.5) +
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
  
  # ---- 🏷️ Gene Labeling ----
  
  # Exclude non-informative gene symbols
  pred_pattern <- "^(Gm[0-9]+|ENSMUSG[0-9]+|ENSG[0-9]+|LOC[0-9]+|C[0-9]+orf[0-9]+|RP[0-9]+-)|Rik$"
  
  # Identify top genes
  top_genes <- res_df %>%
    dplyr::filter(!stringr::str_detect(string = SYMBOL, pattern = pred_pattern), 
                  padj < padj_cutoff,
                  abs(log2FoldChange) > lfc_cutoff) %>%
    dplyr::group_by(Direction) %>%
    dplyr::filter(Direction != "Not Significant") %>%
    #dplyr::slice_max(order_by = padj, n = top_n) %>%
    dplyr::slice_max(order_by = Relevance, n = top_n) %>%
    dplyr::pull(SYMBOL)
  
  # Add a column identifying genes to label
  genes_to_label <- if (!is.null(label_genes)) {
    intersect(label_genes, res_df$SYMBOL) 
  } else {
    top_genes
  }
  
  q <- p + ggrepel::geom_text_repel(
    data = res_df %>% dplyr::filter(SYMBOL %in% genes_to_label),
    aes(label = SYMBOL), 
    direction = "both",
    box.padding = 0.8,                # ↓ smaller padding around label
    point.padding = 0.1,              # minimal space between point and line start
    max.overlaps = nrow(res_df),
    show.legend = FALSE,
    min.segment.length = 0,           # Only draw segments longer than this
    segment.curvature = -0.5,         # Negative = curve upward, positive = downward
    segment.ncp = 50,                 # More control points = smoother curves
    segment.angle = 20,               # Affects entry/exit angles
    segment.size = 0.5,               # Optional: line thickness
    size = 4,                         # text size in mm (1 mm = 2.83 points)
    position = ggbeeswarm::position_quasirandom(width = 0.1, varwidth = TRUE)
  )
  
  # ---- 💾 Save Plots ----
  
  file_extension <- ".pdf"
  file_name <- file.path(output_dir,
                         paste0("Volcano_Plot_", filename, file_extension))
  
  # Open multi-page PDF
  grDevices::cairo_pdf(filename = file_name, width = 8, height = 11.5, onefile = TRUE) 
  
  print(p)  # Page 1
  print(q)  # Page 2
  
  grDevices::dev.off()
  
  # ---- 🪵 Log Output and Return ----
  
  log_info(sample = contrast, 
           step   = "plot_volcano", 
           msg    = glue::glue("Volcano plots saved to '{file_name}'"))
  
  return(invisible(NULL))
}

plot_heatmap <- function(expr_mat, 
                         label_genes         = NULL,
                         filename            = NULL, 
                         output_dir          = NULL,
                         metadata_col        = NULL, 
                         metadata_row        = NULL,
                         col_annotations     = NULL,
                         row_annotations     = NULL,
                         col_gap_by          = NULL,
                         row_gap_by          = NULL,
                         col_cluster_by      = "all",
                         row_cluster_by      = "all",
                         plot_title          = NULL,
                         heatmap_palette     = "rdbu",
                         annotation_palette  = "discrete",
                         border_color        = NA,
                         force_log           = FALSE,
                         show_expr_legend    = TRUE,
                         save_plot           = FALSE,
                         save_matrix         = FALSE) {
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(expr_mat = expr_mat, 
                  logic_vars = list(force_log        = force_log,
                                    show_expr_legend = show_expr_legend,
                                    save_plot        = save_plot,
                                    save_matrix      = save_matrix),
                  filename = filename,
                  output_dir = output_dir)
  
  # Validate `expr_mat`
  if (nrow(expr_mat) < 2) {
    log_warn(sample = "", 
             step   = "plot_heatmap", 
             msg    = "Input has fewer than 2 genes. Skipping heatmap generation.")
    return(NULL)
  }
  
  if (base::anyDuplicated(rownames(expr_mat))) {
    log_warn(sample = "", 
             step   = "prepare_deseq2_input", 
             msg    = "Duplicate gene symbols detected in row names of count matrix.
             Only highest expressing copy will be retained.")
  }
  
  # Validate `label_genes`
  if (!is.null(label_genes)){
    
    if (!is.character(label_genes)) {
      log_error(sample = "", 
                step   = "plot_heatmap",
                msg    = "`label_genes` must be a character vector of gene names.")
    } else {
      missing_genes <- base::setdiff(label_genes, rownames(expr_mat))
      if (length(missing_genes) > 0) {
        log_error(sample = "", 
                  step   = "plot_heatmap",
                  msg    = glue::glue("Following `label_genes` are missing in the expr_mat: {paste(missing_genes, collapse=', ')}"))
      }
    }
  }
  
  # Validate `metadata_col`
  if (!is.null(metadata_col)) {
    
    # (i) Data Frame
    if (!is.data.frame(metadata_col)) {
      log_error(sample = "", 
                step   = "plot_heatmap",
                msg    = "metadata_col must be a data.frame.")
    }
    
    # (ii) Required columns
    if (!"Sample_ID" %in% colnames(metadata_col)) {
      log_error(sample = "", 
                step   = "plot_heatmap",
                msg    = "`metadata_col` must contain a 'Sample_ID' column.")
    }
    
    # (iii) Duplicates
    if (any(duplicated(metadata_col$Sample_ID))) {
      log_error(sample = "", 
                step   = "plot_heatmap",
                msg    = "Duplicate Sample_ID values detected in `metadata_col`.")
    }
  }
  
  # Validate `metadata_row`
  if (!is.null(metadata_row)) {
    
    # (i) Data Frame
    if (!is.data.frame(metadata_row)) {
      log_error(sample = "",
                step   = "plot_heatmap",
                msg    = "`metadata_row` must be a data.frame.")
    }
    
    # (ii) Required columns
    if (!"SYMBOL" %in% colnames(metadata_row)) {
      log_error(sample = "", 
                step   = "plot_heatmap",
                msg    = "`metadata_row` must contain a 'SYMBOL' column.")
    }
    
    # (iii) Duplicates
    if (any(duplicated(metadata_row$SYMBOL))) {
      log_warn(sample = "",
               step   = "plot_heatmap",
               msg    = "Duplicate SYMBOL values detected in `metadata_row`.")
    }
  }
  
  # Validate `plot_title`
  if (!is.null(plot_title) && (!is.character(plot_title) || length(plot_title) != 1)) {
    log_warn(sample = "",
             step   = "plot_heatmap",
             msg    = "`plot_title` should be a single character string. Using default NULL instead.")
  }
  
  # Validate `heatmap_palette`
  if (!heatmap_palette %in% c("vrds", "rdbu")) {
    log_error(sample = "", 
              step   = "plot_heatmap",
              msg    = "`heatmap_palette` must be either 'vrds' or 'rdbu'.")
  }
  
  # Validate `annotation_palette`
  if (!annotation_palette %in% c("discrete", "continuous")) {
    log_error(sample = "", 
              step   = "plot_heatmap",
              msg    = "`annotation_palette` must be either 'discrete' or 'continuous'.")
  }
  
  # Validate `border_color`
  if (!is.null(border_color) && !is.na(border_color)) {
    
    if (!is.character(border_color) || length(border_color) != 1) {
      log_error(sample = "", 
                step   = "plot_heatmap",
                msg    = "`border_color` must be a single character color or NA.")
    }
    
    if (!border_color %in% colors()) {
      log_warn(sample = "", 
               step   = "plot_heatmap",
               msg    = glue::glue("`border_color` '{border_color}' is not a standard R color."))
    }
  }
  
  # 🛡 Metadata Validation
  
  # Helper function for checking single-character metadata columns
  validate_single_column <- function(col_value, col_name, metadata_df = NULL) {
    if (!gtools::invalid(col_value)) {
      if (!is.character(col_value) || length(col_value) != 1) {
        log_error(sample = "", 
                  step   = "plot_heatmap",
                  msg    = glue::glue("'{col_name}' must be a single character value."))
      }
      if (!is.null(metadata_df) && !(col_value %in% colnames(metadata_df))) {
        log_error(sample = "", 
                  step   = "plot_heatmap",
                  msg    = glue::glue("'{col_name}' '{col_value}' must be a column in the metadata dataframe."))
      }
    }
  }
  
  # Helper function for checking multiple-column metadata annotations
  validate_multi_column <- function(cols, col_name, metadata_df) {
    if (!gtools::invalid(cols)) {
      missing_cols <- setdiff(cols, colnames(metadata_df))
      if (length(missing_cols) > 0) {
        log_error(sample = "", 
                  step  = "plot_heatmap",
                  msg   = glue::glue("{col_name} missing in metadata: {paste(missing_cols, collapse=', ')}"))
      }
    }
  }
  
  # Validate col_annotations, col_gap_by, col_cluster_by
  requires_col_metadata <- (!gtools::invalid(col_annotations) ||
                              !gtools::invalid(col_gap_by) ||
                              (!col_cluster_by %in% c("all", "alphabetical")))
  
  if (requires_col_metadata && is.null(metadata_col)) {
    log_error(sample = "", 
              step   = "plot_heatmap",
              msg    = "col_annotations, col_gap_by, or col_cluster_by require metadata_col, but it is NULL.")
  }
  
  if (!is.null(metadata_col)) {
    validate_single_column(col_cluster_by, "col_cluster_by")
    validate_single_column(col_gap_by, "col_gap_by", metadata_col)
    validate_multi_column(col_annotations, "col_annotations", metadata_col)
  }
  
  # Validate row_annotations, row_gap_by, row_cluster_by
  requires_row_metadata <- (!gtools::invalid(row_annotations) ||
                              !gtools::invalid(row_gap_by) ||
                              (!row_cluster_by %in% c("all", "alphabetical")))
  
  if (requires_row_metadata && is.null(metadata_row)) {
    log_error(sample = "", 
              step   = "plot_heatmap",
              msg    = "row_annotations, row_gap_by, or row_cluster_by require metadata_row, but it is NULL.")
  }
  
  if (!is.null(metadata_row)) {
    validate_single_column(row_cluster_by, "row_cluster_by")
    validate_single_column(row_gap_by, "row_gap_by", metadata_row)
    validate_multi_column(row_annotations, "row_annotations", metadata_row)
  }
  
  # ---- 🧪 Prepare Matrix & Normalization ----
  
  # Convert column names to valid R names
  colnames(expr_mat) <- make.names(colnames(expr_mat))
  
  # # Retain ONLY samples in metadata for accurate zero count calculations
  # if (!is.null(metadata_col)){
  #   expr_mat <- expr_mat[, intersect(colnames(expr_mat), rownames(metadata_col)), drop = FALSE]
  # }
  
  # Remove rows with NA gene names
  na_rows <- is.na(rownames(expr_mat))
  expr_mat <- expr_mat[!na_rows, , drop = FALSE]
  log_info(sample = "",
           step   = "plot_heatmap",
           msg    = glue::glue("Removed {sum(na_rows)} rows with missing gene names."))
  
  # Replace missing counts with 0
  expr_mat[is.na(expr_mat)] <- 0
  
  # Convert all columns to numeric safely
  expr_mat <- matrix(as.numeric(as.matrix(expr_mat)),
                     nrow = nrow(expr_mat),
                     ncol = ncol(expr_mat),
                     dimnames = dimnames(expr_mat))
  
  if (any(is.na(expr_mat))) {
    log_error(sample = "",
              step   = "plot_heatmap",
              msg    = "`expr_mat` contains non-numeric values that could not be converted.")
  }
  
  # Remove genes with zero counts across all samples
  zero_genes   <- rownames(expr_mat)[which(rowSums(expr_mat) == 0)]
  expr_mat <- expr_mat[rowSums(expr_mat) != 0, , drop = FALSE]
  log_info(sample = "",
           step   = "plot_heatmap",
           msg    = glue::glue("Removed {length(zero_genes)} genes with zero counts across all samples."))
  
  # Remove samples with zero total reads
  zero_samples <- colnames(expr_mat)[which(colSums(expr_mat) == 0)]
  expr_mat <- expr_mat[, colSums(expr_mat) != 0, drop = FALSE]
  log_info(sample = "",
           step   = "plot_heatmap",
           msg    = glue::glue("Removed {length(zero_samples)} samples with zero total counts."))
  
  # Handle duplicate gene symbols if any
  if (any(duplicated(rownames(expr_mat)))) {
    expr_mat <- expr_mat %>%
      tibble::rownames_to_column("SYMBOL") %>%  
      dplyr::mutate(total_expr = rowSums(across(-SYMBOL))) %>%  # Sum across all samples
      dplyr::group_by(SYMBOL) %>%
      dplyr::slice_max(order_by = total_expr, n = 1, with_ties = FALSE) %>% # Keep row with highest total expression
      dplyr::ungroup() %>%
      tibble::column_to_rownames("SYMBOL") %>%
      dplyr::select(-total_expr)                                # Drop the temporary column
  }
  
  # Convert rownames and colnames to valid R names
  rownames(expr_mat) <- make.names(rownames(expr_mat))
  colnames(expr_mat) <- make.names(colnames(expr_mat))
  
  # Check for high dynamic range to trigger log transform
  quantiles <- stats::quantile(x = as.vector(as.matrix(expr_mat)), probs = c(0, 0.01, 0.99, 1), na.rm = TRUE)
  huge_range <- (quantiles[4] - quantiles[1]) > 100   # Range of values greater than 100
  only_pos <- quantiles[1] >= 0                       # Min value greater than 0
  if ((huge_range & only_pos) | force_log){
    expr_mat <- log2(1 + expr_mat)
  }
  
  # Perform Z-score scaling (across rows i.e. genes)
  expr_mat_scaled <- expr_mat %>% t() %>% scale() %>% t() 
  expr_mat_scaled[is.na(expr_mat_scaled)] <- 0
  
  # ---- 🏷️ Prepare Annotations & Colors ----
  
  # Column (Sample) Annotation
  col_annotation <- if (!is.null(metadata_col)) {
    metadata_col %>%
      dplyr::select(Sample_ID, all_of(col_annotations)) %>%
      dplyr::mutate(Sample_ID = make.names(Sample_ID)) %>%
      dplyr::filter(Sample_ID %in% colnames(expr_mat)) %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames("Sample_ID") %>%
      as.data.frame() %>%
      mutate(across(where(is.factor), as.character))
  } else NULL
  
  # Row (Gene) Annotation
  row_annotation <- if (!is.null(metadata_row)) {
    metadata_row %>%
      dplyr::select(SYMBOL, all_of(row_annotations)) %>%
      dplyr::mutate(SYMBOL = make.names(SYMBOL)) %>%
      dplyr::filter(SYMBOL %in% rownames(expr_mat)) %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames("SYMBOL") %>%
      as.data.frame() %>%
      mutate(across(where(is.factor), as.character))
  } else NULL
  
  # Generate Annotation Palette
  # This is an example of how ann_colors should be specified
  ann_colors <- list(CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
                     GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E"))
  
  ann_colors <- list()
  base_colors <- if (exists("custom_palette")) {
    custom_palette
  } else { 
   log_error(sample = "",
             step   = "plot_heatmap", 
             msg    = "`custom_palette` not defined. Needed for annotation of heatmaps.")
  }
  
  col_list <- base::lapply(X = as.list(col_annotation), FUN = function(x) { as.character(x) %>% unique})
  row_list <- base::lapply(X = as.list(row_annotation), FUN = function(x) { as.character(x) %>% unique})
  ann_list <- c(row_list, col_list)
  
  color_index <- 1
  for (i in seq_along(ann_list)) {  # Iterate through each annotation variable (Eg: CellType) 
    levels <- sort(ann_list[[i]])   # Get levels within each annotation variable (Eg: CT1, CT2)
    n_levels <- length(levels)      # Get number of levels within each annotation variable
    
    palette_colors <- if (annotation_palette == "discrete" | n_levels == 1){
      base_colors[color_index:(color_index + n_levels - 1)]
    } else if (annotation_palette == "continuous"){
      alphas <- seq(1 / n_levels, 1, length.out = n_levels)
      base::sapply(X = alphas, 
                   FUN = function(x) { colorspace::adjust_transparency(col = base_colors[color_index], alpha = x) })
    }
    
    names(palette_colors) <- levels                   # Name each color with levels
    ann_colors <- c(ann_colors, list(palette_colors)) # Append named color palette
    names(ann_colors)[i] <- names(ann_list)[i]        # Name the color palette with corresponding annotation variable name
    color_index <- color_index + n_levels             # Move to next color
  }
  
  # ---- 🎨 Heatmap Palette & Breaks ----
  
  # Define number of Color Breaks 
  n_breaks <- 100
  
  # Define Color Palette for Heatmap 
  heatmap_palette <- if (heatmap_palette == "vrds") {
    viridis::viridis(n_breaks) 
  } else if (heatmap_palette == "rdbu") {
    colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(n_breaks)
  } else {
    log_error(sample = "", 
              step   = "plot_heatmap", 
              msg    = "`heatmap_palette` '{heatmap_palette}' must be either 'vrds' or 'rdbu'.")
  }
  
  # Handle min and max thresholds with soft clamping
  expr_mat_min <- min(expr_mat_scaled, na.rm = TRUE)
  expr_mat_max <- max(expr_mat_scaled, na.rm = TRUE)
  expr_mat_min <- dplyr::case_when(expr_mat_min >= 0 ~ 0, 
                                   expr_mat_min <= -3 ~ -3, 
                                   TRUE ~ expr_mat_min)
  expr_mat_max <- dplyr::case_when(expr_mat_max <= 0 ~ 0, 
                                   expr_mat_max >= 3 ~ 3, 
                                   TRUE ~ expr_mat_max)
  
  # Custom breaks to ensure zero-centering
  if (expr_mat_max == 0) {
    breaks <- seq(from = floor(expr_mat_min), to = 0, length.out = n_breaks)
  } else if (expr_mat_min == 0) {
    breaks <- seq(from = 0, to = ceiling(expr_mat_max), length.out = n_breaks)
  } else {
    breaks <- c(seq(from = floor(expr_mat_min), to = 0, length.out = n_breaks / 2),
                seq(from = expr_mat_max / n_breaks, to = ceiling(expr_mat_max), length.out = n_breaks / 2))
  }
  
  # ---- 🏁 Gaps Logic ----
  
  # Define gaps in columns
  gaps_col <- NULL
  if (!gtools::invalid(col_gap_by)) {
    if (col_gap_by %in% colnames(col_annotation)) {
      gaps_col <- col_annotation %>%
        dplyr::count(.data[[col_gap_by]]) %>%
        dplyr::mutate(n = cumsum(n)) %>%
        dplyr::pull(n) %>%
        .[. < ncol(expr_mat_scaled)]
    } else {
      log_warn(sample = "", 
               step   = "plot_heatmap", 
               msg    = glue::glue("col_gap_by '{col_gap_by}' is absent in metadata_col."))
    }
  }
  
  # Define gaps in rows
  gaps_row <- NULL
  if (!gtools::invalid(row_gap_by)) {
    if (row_gap_by %in% colnames(row_annotation)) {
      gaps_row <- row_annotation %>%
        dplyr::count(.data[[row_gap_by]]) %>%
        dplyr::mutate(n = cumsum(n)) %>%
        dplyr::pull(n) %>%
        .[. < nrow(expr_mat_scaled)]
    } else {
      log_warn(sample = "",
               step   = "plot_heatmap", 
               msg    = glue::glue("row_gap_by '{row_gap_by}' is absent in metadata_row."))
    }
  }
  
  # ---- 📊 Clustering & Ordering ----
  
  # Clustering guideline:
  #   # - Use Ward.D2 linkage only with Euclidean distances (e.g., after z-scoring rows),
  #   #   as Ward assumes Euclidean geometry to minimize within-cluster variance.
  #   # - Use average or complete linkage when using correlation-based distances (1 - cor(x)),
  #   #   since correlation distances are not strictly Euclidean and Ward may produce unintuitive results.
  
  # Logical ordering for columns (Samples)
  if (col_cluster_by %in% colnames(col_annotation)) {
    
    # Initialize column order
    col_order <- c()
    
    # Identify all groups
    col_vars <- col_annotation %>%
      dplyr::pull(col_cluster_by) %>%
      unique() %>%
      sort()  # While calculating gaps_col, we use count() which sorts 
    # alphabetically. So, WE MUST sort col_vars to match gaps_col
    
    # Iterate over each group
    for (col_var in col_vars){
      
      # Get samples belonging to this group
      samples <- col_annotation %>%
        tibble::rownames_to_column("Sample_ID") %>%
        dplyr::filter(.data[[col_cluster_by]] == col_var) %>%
        dplyr::pull(Sample_ID) %>%
        base::intersect(., colnames(expr_mat_scaled))
      
      if (length(samples) > 1){
        # Hierarchical clustering within cluster
        temp_mat  <- expr_mat_scaled[, samples, drop = FALSE]
        col_dist  <- stats::dist(x = t(temp_mat), method = "euclidean")
        col_clust <- stats::hclust(d = col_dist, method = "ward.D2")  
        col_order <- c(col_order, colnames(temp_mat)[col_clust$order])
      } else if(length(samples) == 1){
        col_order <- c(col_order, samples)
      } else{}
    }
  } else if (col_cluster_by == "all") {
    col_dist  <- stats::dist(x = t(expr_mat_scaled), method = "euclidean")
    col_clust <- stats::hclust(d = col_dist, method = "ward.D2")
    col_order <- colnames(expr_mat_scaled)[col_clust$order]
  } else if (col_cluster_by == "alphabetical"){
    col_order <- sort(colnames(expr_mat_scaled))
  } else {
    col_order <- colnames(expr_mat_scaled)
  }
  
  # Logical ordering for rows (Genes)
  if (row_cluster_by %in% colnames(row_annotation)) {
    
    # Initialize column order
    row_order <- c()
    
    # Identify all groups
    row_vars <- row_annotation %>%
      dplyr::pull(row_cluster_by) %>%
      unique() %>%
      sort()  # While calculating gaps_row, we used count() which sorts
    # alphabetically. So, WE MUST sort row_vars to match gaps_row
    
    # Iterate over each group
    for (row_var in row_vars){
      
      # Get genes belonging to this group
      genes <- row_annotation %>%
        tibble::rownames_to_column("SYMBOL") %>%
        dplyr::filter(.data[[row_cluster_by]] == row_var) %>%
        dplyr::pull(SYMBOL) %>%
        base::intersect(., rownames(expr_mat_scaled))
      
      if (length(genes) > 1){
        # Hierarchical clustering within cluster
        temp_mat  <- expr_mat_scaled[genes, , drop = FALSE]
        row_dist  <- stats::dist(x = temp_mat, method = "euclidean")
        row_clust <- stats::hclust(d = row_dist, method = "ward.D2")
        row_order <- c(row_order, rownames(temp_mat)[row_clust$order])
      } else if(length(genes) == 1){
        row_order <- c(row_order, genes)
      } else{}
    }
  } else if (row_cluster_by == "all") {
    row_dist  <- stats::dist(x = expr_mat_scaled, method = "euclidean")
    row_clust <- stats::hclust(d = row_dist, method = "ward.D2")
    row_order <- rownames(expr_mat_scaled)[row_clust$order]
  } else if (row_cluster_by == "alphabetical"){
    row_order <- sort(rownames(expr_mat_scaled))
  } else {
    row_order <- rownames(expr_mat_scaled) # Default to input order
  }
  
  # Final order of rows and columns
  reordered <- expr_mat_scaled[row_order, col_order]
  
  # ---- 🎨 Heatmap Aesthetics ----
  
  # For text to be readable, ideal fontsize = 10 points. 
  # For plot to be pretty, ideal cell_width and cell_height = fontsize + 5 points
  # In 8.5 x 11 inch page, ideal plot area ~ 6 x 8 inch (excluding figure legend,
  # column/row annotations, plot margins). So, plot area ~ 6*72 points wide and 
  # 8*72 points high. 
  # NOTE: If you use col_gaps, then heatmap can get cutoff. In such cases, 
  # increase page width of pdf from 8inch to 10inch in  ---- 💾 Save Plots 
  
  # Set font sizes
  fontsize <- 10
  fontsize_number <- fontsize * 0.8
  angle_col <- 45                           # column label angle
  
  # Set cell width, height (in points) dynamically
  cell_width <- min(fontsize + 5, (6 * 72) / ncol(reordered))
  cell_height <- min(fontsize + 5, (8 * 72) / nrow(reordered))
  
  # Truncate long plot title
  # NOTE: pheatmap() throws error if plot_title = NULL
  plot_title <- if(!is.null(plot_title)) {
    stringr::str_wrap(string = plot_title, width = 20)
  } else { NA }
  
  # Truncate long row labels and display ONLY is cell_height is sufficient
  labels_row <- if (!is.null(label_genes)){
    dplyr::if_else(condition = rownames(reordered) %in% make.names(label_genes),
                   true = stringr::str_trunc(string = rownames(reordered), width = 15), 
                   false = " ")
  } else if (cell_height == fontsize + 5) {
    stringr::str_trunc(string = rownames(reordered), width = 15)
  } else {
    rep(x = " ", times = nrow(reordered))
  }
  
  # Truncate long column labels and display ONLY is cell_width is sufficient
  labels_col <- if (cell_width == fontsize + 5) {
    stringr::str_trunc(string = colnames(reordered), width = 15)
  } else{
    rep(x = " ", times = ncol(reordered))
  }
  
  # ---- 🖼️ Generate Heatmap Plot ----
  
  ph <- pheatmap::pheatmap(mat               = reordered,
                           color             = heatmap_palette,
                           breaks            = breaks,
                           annotation_row    = row_annotation,
                           annotation_col    = col_annotation,
                           annotation_colors = ann_colors,
                           gaps_row          = gaps_row,
                           gaps_col          = gaps_col,
                           
                           cellwidth         = cell_width,     
                           cellheight        = cell_height,  
                           show_rownames     = cell_height >= fontsize,
                           show_colnames     = cell_width >= fontsize,
                           labels_row        = labels_row,
                           labels_col        = labels_col,
                           angle_col         = angle_col,        # column label angle
                           fontsize          = fontsize,         # points; 72 points = 1 inch
                           fontsize_row      = fontsize,         # points
                           fontsize_col      = fontsize,         # points
                           fontsize_number   = fontsize_number,  # points
                           silent            = TRUE, 
                           
                           main              = plot_title,
                           border_color      = border_color,
                           legend            = show_expr_legend,
                           
                           scale                    = "none",
                           cluster_rows             = FALSE,
                           cluster_cols             = FALSE,
                           # clustering_distance_rows = "correlation", #"euclidean",
                           # clustering_distance_cols = "correlation", #"euclidean",
                           # clustering_method        = "average",     #"ward.D2",
                           annotation_legend        = TRUE,
                           annotation_names_row     = FALSE,
                           annotation_names_col     = FALSE,
                           width                    = NA,               # inches
                           height                   = NA,               # inches
                           filename                 = NA)
  
  # Prepare matrix for saving in excel
  ph_mat <- reordered
  
  # Excel maximum limits
  max_rows <- 1e6
  max_cols <- 16384
  
  # Transpose if needed to fit in Excel
  if(ncol(reordered) > max_cols && ncol(reordered) <= max_rows && nrow(reordered) <= max_cols){
    ph_mat <- t(reordered)
  } else if(nrow(reordered) > max_rows){
    warning("Matrix is too large for Excel even after transpose. Consider saving as CSV or binary format.")
  }
  
  # ---- 💾 Save Plots ----
  
  if (save_plot){
    file_extension <- ".pdf"
    file_name <- file.path(output_dir,
                           paste0("Heatmap_Plot_", filename, file_extension))
    
    # Open multi-page PDF
    grDevices::cairo_pdf(filename = file_name, width = 10, height = 11.5, onefile = TRUE) 
    
    # Page 1
    grid::grid.draw(x = ph$gtable)
    
    grDevices::dev.off()
  }
  
  if (save_matrix){
    file_extension <- ".xlsx"
    file_name <- file.path(output_dir,
                           paste0("Heatmap_Matrix_", filename, file_extension))
    
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheetName = "Heatmap_matrix")
    openxlsx::writeData(wb, sheet = "Heatmap_matrix", x = ph_mat, rowNames = TRUE)
    openxlsx::saveWorkbook(wb, file = file_name, overwrite = TRUE)
  }
  
  # ---- 🪵 Log Output and Return ----
  
  log_info(sample = "", 
           step   = "plot_heatmap", 
           msg    = glue::glue("Generated heatmap with {nrow(reordered)} genes and {ncol(reordered)} samples."))
  
  return(invisible(list(ph = ph, mat = ph_mat)))
}

analyze_pathway <- function(res_df = NULL, gene_list = NULL, gene_direction = "UP",
                            species, gmt_dir, output_dir, minsize = 15, maxsize = 500) {
  
  # For fgsea
  set.seed(1234)
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(res_df = res_df, output_dir = output_dir)
  
  # Validate `gmt_dir`
  if (!is.null(gmt_dir)) {
    if (!dir.exists(gmt_dir)){
      log_error(sample = "",
                step   = "pathway_analysis",
                msg    = "`gmt_dir` location '{gmt_dir}' does not exist.") 
    }
  } else {
    log_error(sample = "",
              step   = "pathway_analysis",
              msg    = "`gmt_dir` is NULL i.e. not defined.")
  }
  
  # ---- 🧪 Data Preparation ----
  
  # Initialize result dataframes 
  fgsea_df         <- data.frame()
  gsea_df          <- data.frame()
  ora_df_up        <- data.frame()
  ora_df_down      <- data.frame()
  concise_fgsea_df <- data.frame()
  
  if (!is.null(res_df)){
    
    # Define input genes for GSEA (Ranked list)
    # IMPORTANT: Rank genes from high LFC to low lFC, so +NES ~ up-regulated
    ranked_df <- res_df %>%
      dplyr::distinct(SYMBOL, .keep_all = TRUE) %>%
      dplyr::filter(!is.na(padj), !is.na(SYMBOL)) %>%
      dplyr::mutate(log2FoldChange = as.numeric(log2FoldChange),
                    padj = as.numeric(padj)) %>%
      dplyr::arrange(desc(log2FoldChange))  
    
    ranked_list <- ranked_df$log2FoldChange
    names(ranked_list) <- ranked_df$SYMBOL
    
    # Define input for ORA (Significant genes only)
    sig_genes_up <-  ranked_df %>% 
      dplyr::filter(padj <= 0.05, log2FoldChange > 0) %>% 
      dplyr::pull(SYMBOL)
    sig_genes_down <-  ranked_df %>% 
      dplyr::filter(padj <= 0.05, log2FoldChange < 0) %>% 
      dplyr::pull(SYMBOL)
    universe_genes <- ranked_df$SYMBOL
  
  } else if (!is.null(gene_list) & gene_direction == "UP"){
    
    ranked_df <- NULL
    sig_genes_up <- gene_list
    sig_genes_down <- NULL
    universe_genes <- NULL
 
  } else if (!is.null(gene_list) & gene_direction == "DOWN"){
    
    ranked_df <- NULL
    sig_genes_up <- NULL
    sig_genes_down <- gene_list
    universe_genes <- NULL
  
  } else {
    
    log_error(sample = "",
              step   = "analyze_pathway",
              msg    = "Provide either `res_df` or `gene_list` and `gene_direction (UP/DOWN)`")
  }

  # ---- 🔄 Enrichment Loop ----
  
  gmt_files <- list.files(file.path(gmt_dir, species), full.names = TRUE)
  
  for (gmt_file in gmt_files) {
    
    # gmt_name <- gsub(pattern = "^.*/|.v[0-9].*$", replacement = "", x = gmt_file)
    gmt_name <- gsub(pattern = "^.*/|", replacement = "", x = gmt_file)
    
    # Format gene sets for fgsea
    gmt <- fgsea::gmtPathways(gmt_file)
    
    # Keep only genes present in ranked_df (i.e., GSEA is possible)
    if (!is.null(ranked_df)) {
      gmt <- lapply(X = gmt, FUN = intersect, y = ranked_df$SYMBOL)
    }
    
    # Format gene sets for clusterProfiler and keep only genes present in ranked_list
    # pathway_gene_df <- data.frame(pathways = base::rep(x = names(gmt), times = base::unname(lengths(gmt))),
    #                                genes = unlist(gmt, use.names = FALSE))
    pathway_gene_df <- utils::stack(x = gmt)
    colnames(pathway_gene_df) <- c("genes", "pathways")
    pathway_gene_df <- pathway_gene_df[, c("pathways", "genes")] # Reorder
    
    if (!is.null(res_df)){
      
    # Run fgseaMultilevel (GSEA)
    fgsea_res <- fgsea::fgseaMultilevel(pathways    = gmt,
                                        stats       = ranked_list,
                                        scoreType   = dplyr::case_when(min(ranked_list) > 0 ~ "pos",
                                                                       max(ranked_list) < 0 ~ "neg",
                                                                       TRUE ~ "std"),
                                        minSize     = minsize,
                                        maxSize     = maxsize, 
                                        nPermSimple = 10000)
    
    # Run clusterProfiler GSEA
    gsea_res <- clusterProfiler::GSEA(geneList      = ranked_list,
                                      TERM2GENE     = pathway_gene_df,
                                      minGSSize     = minsize,
                                      maxGSSize     = maxsize,
                                      pvalueCutoff  = 0.05,
                                      pAdjustMethod = "BH",
                                      verbose       = FALSE,
                                      by            = "fgsea")
    
    # Identify overlapping pathways and collapse into major pathways
    concise_fgsea_res <- fgsea::collapsePathways(fgseaRes = fgsea_res,
                                                 pathways = gmt,
                                                 stats    = ranked_list)
    concise_fgsea_res <- fgsea_res %>%
      dplyr::filter(pathway %in% concise_fgsea_res$mainPathways)
    
    # Accumulate results
    if (!is.null(fgsea_res))         { fgsea_df         <- dplyr::bind_rows(fgsea_df,         fgsea_res) }
    if (!is.null(concise_fgsea_res)) { concise_fgsea_df <- dplyr::bind_rows(concise_fgsea_df, concise_fgsea_res) }
    if (!is.null(gsea_res))          { gsea_df          <- dplyr::bind_rows(gsea_df,          gsea_res@result) }
    
    }
    
    # Run clusterProfiler ORA (enricher)
    # NOTE: Avoid using clusterProfiler::enrichGO() as it doesnt use proper 
    # background in universe parameter and includes GO terms outside of the 
    # intended gene set collection.
    ora_res_up <- clusterProfiler::enricher(gene          = sig_genes_up,
                                            universe      = universe_genes,
                                            TERM2GENE     = pathway_gene_df,
                                            minGSSize     = minsize,
                                            maxGSSize     = maxsize,
                                            pvalueCutoff  = 0.05,
                                            pAdjustMethod = "BH",
                                            qvalueCutoff  = 0.2)
    
    ora_res_down <- clusterProfiler::enricher(gene          = sig_genes_down,
                                              universe      = universe_genes,
                                              TERM2GENE     = pathway_gene_df,
                                              minGSSize     = minsize,
                                              maxGSSize     = maxsize,
                                              pvalueCutoff  = 0.05,
                                              pAdjustMethod = "BH",
                                              qvalueCutoff  = 0.2)
    
    # Accumulate results
    if (!is.null(ora_res_up))        { ora_df_up        <- dplyr::bind_rows(ora_df_up,        ora_res_up@result) }
    if (!is.null(ora_res_down))      { ora_df_down      <- dplyr::bind_rows(ora_df_down,      ora_res_down@result) }
    
  }
  
  # ---- 🧹 Column Standardization & Formatting ----
  
  # Rename columns consistently across different methods
  lookup <- c(pathway = "ID", 
              geneID = "leadingEdge", geneID = "core_enrichment", 
              K = "size", K = "setSize", 
              padj = "p.adjust", 
              pval = "pvalue")
  
  # Put your data frames in a named list
  dfs <- list(fgsea_df      = fgsea_df,
              gsea_df       = gsea_df,
              ora_df_up     = ora_df_up,
              ora_df_down   = ora_df_down)
  df_names <- names(dfs)
  
  dfs <- lapply(X = df_names, FUN = function(df_name) {
    
    # Extract specific df
    df <- dfs[[df_name]]
    if (nrow(df) == 0 || is.null(df)) return(df)  # skip empty data frames
    
    # Add Direction column
    if (df_name == "ora_df_up"){
      df$Direction <- "Upregulated"
    } else if (df_name == "ora_df_down"){
      df$Direction <- "Downregulated"
    } else {
      df <- df %>%
        dplyr::mutate(Direction = dplyr::case_when(NES > 0 ~ "Upregulated",
                                                   NES < 0 ~ "Downregulated",
                                                   TRUE    ~ "No change"))
    }
    
    # ORA df specific formatting
    # k <- # overlap between pathway and input (sig_genes_up/sig_genes_down/gmt/ranked_list)
    # n <- # overlap between collection and input (sig_genes_up/sig_genes_down/gmt/ranked_list)
    # K <- # overlap between pathway and universe
    # N <- # overlap between collection and universe
    if (df_name %in% c("ora_df_up", "ora_df_down")){
      df <- df %>%
        tidyr::separate(col = GeneRatio, into = c("k", "n")) %>%
        tidyr::separate(col = BgRatio,   into = c("K", "N")) %>%
        dplyr::mutate(across(.cols = c(k, n, K, N), .fns = as.numeric)) %>%
        dplyr::mutate(GeneRatio       = k / n,
                      BackgroundRatio = K / N,
                      EnrichmentRatio = GeneRatio / BackgroundRatio,
                      combined_score  = GeneRatio * -log10(p.adjust),
                      NES             = NA_integer_)
    }
    
    # Standardize column names based on look up table
    df <- df %>%
      dplyr::rename(any_of(lookup))
    
    # Format results
    # IMPORTANT: gsea_df and ora_df store geneID as string "CXCL11/CCL2/..."
    # fgsea_df stores geneID as list of vectors which we convert to string "CXCL11/CCL2/..."
    df <- df %>%
      tibble::remove_rownames() %>%
      tidyr::separate(col = pathway, into = c("Collection", "Description"), sep = "_", extra = "merge") %>%
      dplyr::mutate(Description = base::gsub(pattern = "_", replacement = " ", x = Description),
                    geneID = base::sapply(X = geneID, FUN = paste, collapse = "/")) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(leading_edge_size = length(unlist(stringr::str_split(geneID, "/")))) %>%
      dplyr::ungroup() %>%
      as.data.frame() %>%
      dplyr::select(Collection, Description, leading_edge_size, K, padj, NES, Direction, everything(), -geneID, geneID) %>%
      dplyr::filter(padj <= 0.05)
    
    # Generate separate column for each gene in enriched pathways
    max_len <- max(df$leading_edge_size, na.rm = TRUE)
    if (is.finite(max_len) & max_len > 0) {
      df <- df %>%
        tidyr::separate(col = geneID, into = paste0("gene", 1:max_len), sep = "/", remove = TRUE, fill = "right")
    }
    
    return(df)
  })
  
  # Restore names
  dfs <- stats::setNames(dfs, df_names)
  
  # Put them back into separate variables if needed
  fgsea_df <- dfs$fgsea_df
  gsea_df  <- dfs$gsea_df
  ora_df   <- dplyr::bind_rows(dfs$ora_df_up, dfs$ora_df_down)
  
  # ---- 🤝 Consensus Result ----
  
  consensus_df <- dplyr::bind_rows(fgsea_df %>% dplyr::mutate(method = "FGSEA"), 
                                   gsea_df %>% dplyr::mutate(method = "GSEA"), 
                                   ora_df %>% dplyr::mutate(method = "ORA")) %>%
    dplyr::add_count(Collection, Description, Direction, name = "n_methods") %>%
    dplyr::mutate(Consensus = Direction) %>%
    dplyr::arrange(Collection, Description, desc(NES)) %>%
    dplyr::select(n_methods, method, Consensus, Collection, Description, 
                  leading_edge_size, K, padj, NES, Direction, everything(), 
                  -starts_with("gene",ignore.case = FALSE),
                  starts_with("gene", ignore.case = FALSE)) 
  
  # ---- 💾 Save Outputs ----
  
  pathway_results_list <- list(fgsea     = fgsea_df, 
                               gsea      = gsea_df, 
                               ora       = ora_df,
                               consensus = consensus_df)
  # Save as excel
  file_name <- file.path(output_dir, "Pathway_results.xlsx")
  wb <- openxlsx::createWorkbook()
  for (i in seq_along(pathway_results_list)) {
    openxlsx::addWorksheet(wb, sheetName = names(pathway_results_list)[i])
    openxlsx::writeData(wb, sheet = names(pathway_results_list)[i], x = pathway_results_list[[i]], rowNames = FALSE)
  }
  openxlsx::saveWorkbook(wb, file_name, overwrite = TRUE)
  
  # ---- 🪵 Log Output and Return ----
  
  log_info(sample = "", 
           step   = "pathway_analysis", 
           msg    = glue::glue("Pathway analysis results saved to '{file_name}'."))
  
  return(invisible(list(fgsea = fgsea_df, gsea = gsea_df, ora = ora_df, consensus = consensus_df)))
}

analyze_tf <- function(expr_mat, res_df, species, output_dir,
                       stats = c("ulm", "mlm", "viper"), 
                       minsize = 5, top_n = 500) {
  
  # decoupleR analysis RECOMMENDED over progeny/dorothea as it can run 
  # multiple algorithms & give consensus result
  
  # stats can be "aucell", fgsea", "gsva", "mdt", "mlm", "ora", "udt", "ulm", 
  # "viper", "wmean", "wsum". We use ulm, mlm, and viper provide a balance of 
  # sensitivity and specificity
  
  set.seed(1234)
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(expr_mat = expr_mat, res_df = res_df) 
  
  if (!species %in% c("Homo sapiens", "Mus musculus")) {
    log_error(sample = "", 
              step   = "tf_analysis", 
              msg    = "`species` must be either 'Homo sapiens' or 'Mus musculus'.")
  }
  
  if (is.null(expr_mat) && !is.null(res_df)){
    log_info(sample = "", 
             step   = "tf_analysis", 
             msg    = "`expr_mat` is NULL. Using `res_df` for Transcription factor analysis.")
  } else if (!is.null(expr_mat) && is.null(res_df)){
    log_info(sample = "", 
             step   = "tf_analysis", 
             msg    = "`res_df` is NULL. Using `expr_mat` for Transcription factor analysis.")
  } else if (is.null(expr_mat) && is.null(res_df)){
    log_error(sample = "", 
              step   = "tf_analysis", 
              msg    = "Either `expr_mat` or `res_df` is needed for Transcription factor analysis.
              Both cannot be NULL.")
  }
  
  # ---- 🕸️ Load Regulatory Networks ----
  
  # Map species to decoupleR organism format
  organism <- dplyr::case_when(species == "Homo sapiens" ~ "human",
                               species == "Mus musculus" ~ "mouse",
                               TRUE ~ "rat")
  
  # PROGENy: Pathway activity (Top weights per pathway)
  progeny_pathway_net <- decoupleR::get_progeny(organism = organism, top = top_n)
  
  # CollecTRI: Transcription Factor gene regulatory network (GRN)
  # CollecTRI is chosen over DoRothEA for high-confidence, literature-backed interactions
  
  # get_collectri() returns mor column without edge weights from -1 or +1
  collectri_tf_net <- decoupleR::get_collectri(organism = organism, split_complexes = FALSE)
  
  # get_dorothea() returns mor column with edge weights from -1 through +1
  # dorothea_tf_net <- decoupleR::get_dorothea(organism = organism, levels = c('A', 'B', 'C'))
  # viper works better with edge weights but dorothea's edge weights are solely based on confidence
  
  # ---- 📝 Format Input for DecoupleR ----
  
  # Calculate a signed significance score for gene ranking.
  # We combine statistical confidence (-log10 p-value) with directionality
  # using sign(log2FoldChange), rather than the raw fold-change magnitude.
  # This prevents genes with extreme fold changes but weak statistical support
  # from dominating the ranking, and aligns with best practices for
  # transcription factor and pathway activity inference (e.g. VIPER, decoupleR).
  
  if (!is.null(expr_mat)) {
    input_mat <- as.matrix(expr_mat)
  } else {
    # If using res_df, create't-statistic' equivalent using padj & log2FoldChange
    input_mat <- res_df %>%
      as.data.frame() %>%
      dplyr::mutate(t = -log10(pmax(padj, .Machine$double.xmin)) * sign(log2FoldChange)) %>%
      dplyr::filter(!is.na(t), !is.na(SYMBOL)) %>%
      dplyr::select(SYMBOL, t) %>%
      tibble::column_to_rownames("SYMBOL") %>%
      as.matrix()
  }
  
  # ---- 📉 Run DecoupleR (Consensus Statistics) ----
  
  # Calculate Pathway Activity
  pathway_df <- decoupleR::decouple(mat        = input_mat,
                                    network    = progeny_pathway_net,
                                    statistics = stats, 
                                    minsize    = minsize)
  
  # Calculate Transcription Factor Activity
  tf_df <- decoupleR::decouple(mat        = input_mat, 
                               network    = collectri_tf_net,
                               statistics = stats, 
                               minsize    = minsize)
  
  # tf_df <- decoupleR::decouple(mat        = input_mat,
  #                              network    = dorothea_tf_net,
  #                              statistics = "viper",
  #                              minsize    = minsize)
  
  # ---- 🧹 Significance (BH Correction) ----
  
  # Function to calculate padj for each method i.e ulm/mlm/viper
  calc_padj <- function(df) {
    df %>%
      dplyr::group_by(statistic) %>%
      dplyr::mutate(padj = stats::p.adjust(p_value, method = "BH")) %>%
      dplyr::ungroup()
  }
  
  pathway_df <- calc_padj(pathway_df)
  tf_df <- calc_padj(tf_df)
  
  # ---- 💾 Save Outputs ----
  
  tf_results_list <- list(tf      = tf_df, 
                          pathway = pathway_df) 
                              
  # Save as excel
  file_name <- file.path(output_dir, "TF_results.xlsx")
  wb <- openxlsx::createWorkbook()
  for (i in seq_along(tf_results_list)) {
    openxlsx::addWorksheet(wb, sheetName = names(tf_results_list)[i])
    openxlsx::writeData(wb, sheet = names(tf_results_list)[i], x = tf_results_list[[i]], rowNames = FALSE)
  }
  openxlsx::saveWorkbook(wb, file_name, overwrite = TRUE)
  
  # ---- 🪵 Log Output and Return ----
  
  log_info(sample = "", 
           step   = "tf_analysis", 
           msg    = glue::glue("TF analysis results saved to '{file_name}'."))
  
  return(invisible(list(pathway = pathway_df,
                        tf      = tf_df)))
                       
}

plot_pathway <- function(pathway_df, expr_mat, metadata, method, output_dir){
  
 
  # pathway_df$Collection, pathway_df$leading_edge_size,
  # pathway_df$Direction, pathway_df$padj, pathway_df$Description
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(expr_mat = expr_mat, metadata = metadata, output_dir = output_dir)  

  if (is.null(method) || !method %in% c("GSEA", "ORA")) {
    log_error(sample = "",
              step   = "plot_pathway",
              msg    = glue::glue("`method` '{method}' is invalid. Must be 'GSEA' or 'ORA'."))
  }
  
  # ---- 📊 Define Plotting Mapping Logic ----
  
  # Map method -> plotting aesthetics
  x_axis_map    <- c(ORA = "GeneRatio", GSEA = "NES")
  x_label_map   <- c(ORA = "Gene Ratio", GSEA = "Normalized Enrichment Score (NES)")
  
  x_axis  <- x_axis_map[[method]]
  x_label <- x_label_map[[method]]
  
  size_col      <- "leading_edge_size"
  color_col     <- "Direction"
  alpha_col     <- "padj"
  
  plot_colors <- c("Upregulated" = "#E69F00", "Downregulated" = "#56B4E9")
  
  # Pre-process descriptions for better wrapping in plots
  pathway_df <- pathway_df %>%
    dplyr::mutate(Description = base::gsub(pattern = "_", replacement = " ", x = Description),
                  Description = stringr::str_wrap(string = Description, width = 30))
  
  collections   <- unique(pathway_df$Collection)
  dot_plots     <- list()
  bar_plots     <- list()
  
  # ---- 🔄 Iterate Through Pathway Collections ----
  
  for (collection in collections) {
    
    # Subset and rank pathways based on score_var
    plot_df <- pathway_df %>% 
      dplyr::filter(Collection == collection) %>% 
      dplyr::filter(!is.na(.data[[x_axis]])) %>%
      dplyr::arrange(dplyr::desc(.data[[x_axis]]))
    
    if (nrow(plot_df) == 0) next
    
    # Pad with empty rows if fewer than 20 pathways for consistent plot scaling
    n_missing <- 20 - nrow(plot_df)
    if (n_missing > 0) {
      plot_df <- dplyr::bind_rows(plot_df, 
                                  data.frame(Description = as.character(base::seq_len(n_missing))))
    }
    
    # Dynamic Y-axis text sizing
    max_label_len <- base::max(base::nchar(plot_df$Description), na.rm = TRUE)
    y_text_size <- dplyr::case_when(max_label_len > 50 ~ 6,
                                    max_label_len > 35 ~ 7,
                                    max_label_len > 25 ~ 8,
                                    TRUE ~ 10)
    
    # Calculate limits dynamically per collection
    x_min <- if (method == "GSEA") { base::pmin(0, floor(min(plot_df[[x_axis]], na.rm = TRUE))) } else  { 0 }
    x_limits <- c(x_min, NA)
    
    # ---- 📈 Generate Bar Plot ---- 
    
    bar_p <- ggplot2::ggplot(data = plot_df,
                             mapping = aes(x     = .data[[x_axis]],
                                           y     = stats::reorder(Description, .data[[x_axis]]),
                                           fill  = .data[[color_col]],
                                           alpha = -log10(.data[[alpha_col]]))) +
      ggplot2::geom_col(width = 0.75, na.rm = TRUE) +
      ggplot2::labs(x = x_label, 
                    y = "", 
                    title = base::paste("Top", collection, "Pathways"), 
                    fill = "Direction") +
      ggplot2::theme_classic() +
      custom_theme +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = y_text_size)) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::scale_x_continuous(limits = x_limits, expand = ggplot2::expansion(mult = c(0, 0.05))) +
      ggplot2::scale_alpha_continuous(range = c(0.5, 1)) +
      ggplot2::scale_fill_manual(values = plot_colors) +
      guides(fill  = guide_legend(override.aes = list(shape = 22, size = 6)),
             color = guide_legend(override.aes = list(shape = 22, size = 6)),
             alpha = guide_legend(override.aes = list(shape = 22, size = 6))) +
      ggplot2::geom_text(aes(label = .data[[size_col]]), x = 0, hjust = -0.5, size = 3, show.legend = FALSE)
    
    bar_plots[[collection]] <- bar_p
    
    # ---- 🟢 Generate Dot Plot ---- 
    
    #size_vals <- c(min(plot_df[[size_col]], na.rm = TRUE), max(plot_df[[size_col]], na.rm = TRUE))
    size_vals <- plot_df$leading_edge_size[!base::is.na(plot_df$leading_edge_size)]
    breaks <- as.vector(floor(stats::quantile(size_vals, na.rm = TRUE) / 10) * 10)
    
    dot_p <- ggplot2::ggplot(data = plot_df,
                             mapping = aes(x     = .data[[x_axis]],
                                           y     = stats::reorder(Description, .data[[x_axis]]),
                                           fill  = .data[[color_col]],
                                           alpha = -log10(.data[[alpha_col]]),
                                           color = .data[[color_col]],
                                           size  = .data[[size_col]])) +
      ggplot2::geom_point(na.rm = TRUE) +
      ggplot2::labs(x = x_label, 
                    y = "", 
                    title = paste("Top", collection, "Pathways"), 
                    color = "Direction", 
                    size = "Counts") +
      ggplot2::theme_classic() +
      custom_theme +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = y_text_size)) +
      ggplot2::coord_cartesian(clip = "off") + 
      ggplot2::scale_x_continuous(limits = x_limits, expand = ggplot2::expansion(mult = c(0, 0.05))) +
      ggplot2::scale_alpha_continuous(range = c(0.5, 1)) +
      ggplot2::scale_color_manual(values = plot_colors) +
      ggplot2::scale_fill_manual(values = plot_colors) +  # need for coloring the legend
      ggplot2::scale_size(breaks = unique(breaks)) 
    
    dot_plots[[collection]] <- dot_p
    
    # ---- 🔥 ️ Generate Heatmap Multi-page PDF ----
    
    if (!is.null(expr_mat)) {
      
      heatmap_plots <- list()
      pathways <- plot_df %>% 
        dplyr::filter(!is.na(Direction)) %>% 
        dplyr::pull(Description) %>% 
        base::unique()
      
      for (pathway in pathways) {
        
        # Extract genes for this specific pathway from the long-format dataframe
        plot_genes <- pathway_df %>%
          dplyr::filter(Collection == collection, Description == pathway) %>%
          tidyr::pivot_longer(cols = dplyr::matches("^(gene|Gene)[0-9]+$"), 
                              values_to = "gene") %>%
          dplyr::pull(gene) %>%
          base::trimws() %>%
          stats::na.omit() %>%
          .[. != ""] %>%          # remove empty strings
          base::unique() %>%
          base::intersect(rownames(expr_mat))
        
        # Skip plotting if less than 2 genes
        if (base::length(plot_genes) < 2) next
        
        plot_title <- stringr::str_wrap(string = pathway, width = 30)
        
        # Plot heatmap
        ph <- plot_heatmap(expr_mat            = expr_mat[plot_genes, ,drop = FALSE], 
                           label_genes         = NULL,
                           filename            = NULL,
                           output_dir          = NULL,
                           metadata_col        = metadata, 
                           metadata_row        = NULL,
                           col_annotations     = proj.params$heatmap$col_annotations,
                           row_annotations     = proj.params$heatmap$row_annotations,
                           col_gap_by          = proj.params$heatmap$col_gap_by,
                           row_gap_by          = proj.params$heatmap$row_gap_by,
                           col_cluster_by      = proj.params$heatmap$col_cluster_by,
                           row_cluster_by      = proj.params$heatmap$row_cluster_by,
                           plot_title          = plot_title,
                           heatmap_palette     = proj.params$heatmap$heatmap_palette,
                           annotation_palette  = proj.params$heatmap$annotation_palette,
                           border_color        = proj.params$heatmap$border_color,
                           force_log           = proj.params$heatmap$force_log,
                           show_expr_legend    = proj.params$heatmap$show_expr_legend,
                           save_plot           = FALSE,
                           save_matrix         = FALSE)
        
        heatmap_plots[[pathway]] <- ph$ph$gtable
      }
      
      # Save stored heatmaps as pdf
      if (length(heatmap_plots) > 0) {
        
        file_extension <- ".pdf"
        file_name <- file.path(output_dir, paste0("Heatmap_", collection, "_", method, file_extension))
        
        # Open multi-page PDF
        grDevices::cairo_pdf(filename = file_name, width = 8, height = 11.5, onefile = TRUE)  
        
        for (ht in heatmap_plots) {
          grid::grid.newpage()
          grid::grid.draw(ht)
        }
        grDevices::dev.off() 
      }
    }
  } 
  
  # ---- 💾 Save Consolidated Summary Plots ----
  
  summary_plots <- base::list(Bar = bar_plots, Dot = dot_plots)
  
  for (type in names(summary_plots)) {
    
    file_extension <- ".pdf"
    file_name <- file.path(output_dir, paste0(type, "_plot_pathways_", method, file_extension))
    
    ggplot2::ggsave(filename = file_name,
                    plot     = cowplot::plot_grid(plotlist = summary_plots[[type]], ncol = 3, align = "hv"),
                    device   = grDevices::cairo_pdf,
                    width    = 3 * 6, 
                    height   = ceiling(length(summary_plots[[type]]) / 3) * 6, 
                    units    = "in",
                    dpi      = 300,
                    bg       = "white")
  }
  
  # ---- 🪵 Log Output and Return ----
  
  log_info(sample = "", 
           step   = "plot_pathway", 
           msg    = glue::glue("Successfully generated {method} visualizations in {output_dir}"))
  
  return(invisible(NULL))
}

plot_tf <- function(tf_df, metadata, output_dir, 
                    contrast = "Target-Reference", 
                    top_n    = 20) {
  
  set.seed(1234)
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(metadata = metadata, output_dir = output_dir)
  
  # ---- 🧪 Setup Parameters ----
  
  target    <- stringr::str_split(string = contrast, pattern = "-")[[1]][1]
  reference <- stringr::str_split(string = contrast, pattern = "-")[[1]][2]
  
  # Determine plot type based on the 'condition' column
  # Single condition (usually "t") implies DE results -> Bar Plot
  # Multiple conditions imply sample-level scores -> Heatmap
  plot_type <- if (length(unique(tf_df[["condition"]])) == 1) { "bar" } else { "heatmap" }
  
  log_info(sample = "", 
           step   = "plot_tf", 
           msg    = glue::glue("Detected plot mode: '{plot_type}' based on `condition` column in `tf_df`."))
  
  bar_plots     <- list()    
  heatmap_plots <- list()
  
  # NOTE: wsum returns wsum, norm_wsum and corr_wsum.
  # wsum      (DONT USE) : Biased toward larger gene sets (more genes → bigger sum)
  # norm_wsum (DONT USE) : Adjusts for pathway length so small and large gene sets are comparable
  # corr_sum  (DONT USE) : corrects for high correlation as it can make enrichment appear stronger
  
  # ---- 🖼️ Generate Plots ----
  
  # tf_df structure:
  # `statistic` column : has method used like ulm, mlm, viper, etc
  # `source`    column : has transcription factors
  # `condition` column : has samples names or test-statistic used  like "t"
  # `score`     column : has TF activity score
  
  # ---- 🖼️ Generate Plots per Statistic Method ----
  
  for (stat in unique(tf_df$statistic)) {
    
    # CASE A: BAR PLOT (DE-based activity)
    if (plot_type == "bar") {
      
      # Select top_n TFs for both directions (Upregulated and Downregulated)
      top_tf <- tf_df %>%
        dplyr::filter(statistic == stat) %>%
        dplyr::mutate(Direction = dplyr::case_when(score < 0 ~ "Downregulated",
                                                   score > 0 ~ "Upregulated",
                                                   TRUE      ~ "No change")) %>%
        dplyr::group_by(Direction) %>%
        dplyr::slice_max(order_by = abs(score), n = top_n, with_ties = FALSE) %>%
        dplyr::ungroup()
      
      bar_p <- ggplot2::ggplot(data    = top_tf, 
                               mapping = aes(x = stats::reorder(source, score),
                                             y = score, 
                                             fill = score)) +
        ggplot2::geom_col(width = 0.75, na.rm = TRUE) +
        ggplot2::scale_fill_gradient2(low = "darkblue", high = "indianred", mid = "whitesmoke", midpoint = 0) +
        ggplot2::labs(x = "", 
                      y = "Activity Score", 
                      title = paste0("Top TFs (", stat, " method)"), 
                      fill = "Score") +
        custom_theme +
        ggplot2::coord_cartesian(clip = "off") +
        # Add dynamic labels for biological groups
        ggplot2::geom_text(label = paste0("Activated in ", target),
                           x = top_tf$source[which.max(top_tf$score)],
                           y = ceiling(max(top_tf$score)) + 0.5, 
                           hjust = 1, color = "indianred", fontface = "bold", size = 3) +
        ggplot2::geom_text(label = paste0("Activated in ", reference),
                           x = top_tf$source[which.min(top_tf$score)],
                           y = floor(min(top_tf$score)) - 0.5, 
                           hjust = 0, color = "darkblue", fontface = "bold", size = 3)
      
      bar_plots[[stat]] <- bar_p
    }
  
    # CASE B: HEATMAP (Sample-based activity)
    if (plot_type == "heatmap") {
      
      # Select TFs with the highest variance (standard deviation) across samples
      top_tf_names <- tf_df %>%
        dplyr::filter(statistic == stat) %>%
        dplyr::group_by(source) %>%
        dplyr::summarise(std = stats::sd(score, na.rm = TRUE), .groups = "drop") %>%
        dplyr::slice_max(order_by = abs(std), n = top_n, with_ties = FALSE) %>%
        dplyr::pull(source)
      
      # Pivot data to matrix format: TFs (Rows) x Samples (Columns)
      tf_mat <- tf_df %>%
        dplyr::filter(statistic == stat, source %in% top_tf_names) %>%
        tidyr::pivot_wider(id_cols = "condition", names_from = "source", values_from = "score") %>%
        tibble::column_to_rownames("condition") %>%
        as.matrix() %>%
        t()
      
      plot_title <- paste0("Top TFs (", stat, ") method")
      
      # Plot heatmap
      ph <- plot_heatmap(expr_mat            = tf_mat, 
                         label_genes         = NULL,
                         filename            = NULL,
                         output_dir          = NULL,
                         metadata_col        = metadata, 
                         metadata_row        = NULL,
                         col_annotations     = proj.params$heatmap$col_annotations,
                         row_annotations     = proj.params$heatmap$row_annotations,
                         col_gap_by          = proj.params$heatmap$col_gap_by,
                         row_gap_by          = proj.params$heatmap$row_gap_by,
                         col_cluster_by      = proj.params$heatmap$col_cluster_by,
                         row_cluster_by      = proj.params$heatmap$row_cluster_by,
                         plot_title          = plot_title,
                         heatmap_palette     = proj.params$heatmap$heatmap_palette,
                         annotation_palette  = proj.params$heatmap$annotation_palette,
                         border_color        = proj.params$heatmap$border_color,
                         force_log           = proj.params$heatmap$force_log,
                         show_expr_legend    = proj.params$heatmap$show_expr_legend,
                         save_plot           = FALSE,
                         save_matrix         = FALSE)
      
      heatmap_plots[[stat]] <- ph$ph$gtable
    }
  }
  
  # ---- 💾 Save and Export Results ----
  
  # Export Heatmaps
  if (length(heatmap_plots) > 0) {
    
    file_extension <- ".pdf"
    file_name <- file.path(output_dir, paste0("Heatmap_TF_Activity_", contrast, file_extension))
    
    # Open multi-page PDF
    grDevices::cairo_pdf(filename = file_name, width = 8, height = 11.5, onefile = TRUE) 
    
    for (ht in heatmap_plots) {
      grid::grid.newpage()
      grid::grid.draw(ht)
    }
    grDevices::dev.off() 
  }
  
  # Export Bar Plots
  if (length(bar_plots) > 0) {
    
    file_extension <- ".pdf"
    file_name <- file.path(output_dir, paste0("Bar_plot_TF_Activity", contrast, file_extension))
    ggplot2::ggsave(filename = file_name,
                    plot     = cowplot::plot_grid(plotlist = bar_plots, align = "hv", ncol = 1),
                    device   = grDevices::cairo_pdf,
                    width    = 8.5,
                    height   = 4 * length(bar_plots),
                    units    = "in",
                    dpi      = 300,
                    bg       = "white")
  }
  
  # ---- 🪵 Log Output and Return ----
  
  log_info(sample = "", 
           step   = "plot_tf", 
           msg    = glue::glue("TF Activity plotting for {contrast} completed successfully."))
  
  return(invisible(NULL))
}

get_annotations <- function() {
  
  # ---- ⚙️ Initialize & Connection ----
  
  species_list <- c("Homo sapiens", "Mus musculus")
  annotations_list <- list()

  # Connect to AnnotationHub 
  hub <- AnnotationHub::AnnotationHub()
  
  for (species in species_list) {
    
    # ---- 🔍 Query Database ----
    
    log_info(sample = species, 
             step   = "get_annotations", 
             msg    = glue::glue("Fetching Ensembl Database for '{species}'"))
    
    hub_db <- AnnotationHub::query(x           = hub, 
                                   pattern     = c(species, "EnsDb"), 
                                   ignore.case = TRUE)
    
    # Acquire the latest version available in the hub
    latest_id <- hub_db %>%
      mcols() %>%
      as.data.frame() %>%
      dplyr::arrange(desc(rdatadateadded)) %>%
      head(n = 1) %>%
      rownames()
    
    if (length(latest_id) == 0) {
      log_error(sample = species, 
                step   = "get_annotations", 
                msg    = "Could not find a valid EnsDb in AnnotationHub.")
    }
    
    # Download the appropriate Ensembldb database
    ensdb <- hub_db[[latest_id]]
    
    # ---- 🧬 Extract ENSEMBL Annotations ----
    
    ensembl <- ensembldb::genes(x = ensdb, 
                                return.type = "data.frame") %>%
      dplyr::rename(ENSEMBL_ID         = gene_id,
                    ENSEMBL_SYMBOL     = gene_name,
                    ENSEMBL_BIOTYPE    = gene_biotype,
                    START              = gene_seq_start,
                    END                = gene_seq_end,
                    CHR                = seq_name,
                    STRAND             = seq_strand,
                    DESCRIPTION        = description,
                    ENSEMBL_TRANSCRIPT = canonical_transcript) %>%
      dplyr::mutate(ENSEMBL_SYMBOL = dplyr::if_else(nchar(ENSEMBL_SYMBOL) == 0, 
                                                    NA_character_, ENSEMBL_SYMBOL)) %>%
      dplyr::select(ENSEMBL_ID, ENSEMBL_TRANSCRIPT, ENSEMBL_SYMBOL, ENSEMBL_BIOTYPE,
                    START, END, CHR, STRAND, DESCRIPTION)
    
    # ---- 🔢 Extract ENTREZ Annotations ----
    
    log_info(sample = species, 
             step   = "get_annotations", 
             msg    = "Fetching OrgDb mappings...")
    
    org_db <- if (species == "Homo sapiens") {
      requireNamespace("org.Hs.eg.db", quietly = TRUE)
      org.Hs.eg.db::org.Hs.eg.db 
    } else {
      requireNamespace("org.Mm.eg.db", quietly = TRUE)
      org.Mm.eg.db::org.Mm.eg.db
    }
    
    entrez <- AnnotationDbi::select(x = org_db,
                                    keys = AnnotationDbi::keys(org_db),
                                    columns = c("ENSEMBL", "SYMBOL", "GENETYPE")) %>%
      dplyr::rename(ENTREZ_ID      = ENTREZID,
                    ENSEMBL_ID     = ENSEMBL,
                    ENTREZ_SYMBOL  = SYMBOL,
                    ENTREZ_BIOTYPE = GENETYPE)
    
    # ---- 🤝 Merge and Clean ----
    
    annotations <- dplyr::full_join(ensembl, entrez, by = c("ENSEMBL_ID"="ENSEMBL_ID")) %>%
      dplyr::select(ENSEMBL_ID, ENSEMBL_TRANSCRIPT, ENTREZ_ID, ENSEMBL_SYMBOL, 
                    ENTREZ_SYMBOL, ENSEMBL_BIOTYPE, ENTREZ_BIOTYPE, START, END,
                    CHR, STRAND, DESCRIPTION) %>%
      dplyr::distinct(ENSEMBL_ID, .keep_all = TRUE) # Remove potential duplication from 1:M mappings
    
    # Store Output 
    annotations_list[[species]] <- annotations
    
    log_info(sample = species, 
             step   = "get_annotations", 
             msg    = glue::glue("Successfully retrieved {nrow(annotations)} gene annotations."))
  }
  
  # ---- 🪵 Log Output and Return ----
  
  log_info(sample = "", step = "get_annotations", 
           msg =glue::glue("Successfully retrieved gene annotations."))
  
  # Return: list(human, mouse) 
  return(invisible(annotations_list))
}

add_annotation <- function(df, ann_list = NULL, remove_ann_col = TRUE) {
  
  # ---- ⚙️ Retrieve & Flatten Annotations ----
  
  # Retrieve the list of dataframes
  if (is.null(ann_list)) {
    log_info(sample = "", 
             step = "add_annotation", 
             msg = "No annotation list provided. Fetching now...")
    ann_list <- get_annotations()
  }
  
  # Flatten all annotation data into one named list
  named_lists <- list(ensembl_id_human     = ann_list$`Homo sapiens`$ENSEMBL_ID,
                      entrez_id_human      = ann_list$`Homo sapiens`$ENTREZ_ID,
                      ensembl_symbol_human = ann_list$`Homo sapiens`$ENSEMBL_SYMBOL,
                      entrez_symbol_human  = ann_list$`Homo sapiens`$ENTREZ_SYMBOL,
                      ensembl_id_mouse     = ann_list$`Mus musculus`$ENSEMBL_ID,
                      entrez_id_mouse      = ann_list$`Mus musculus`$ENTREZ_ID,
                      ensembl_symbol_mouse = ann_list$`Mus musculus`$ENSEMBL_SYMBOL,
                      entrez_symbol_mouse  = ann_list$`Mus musculus`$ENTREZ_SYMBOL)
  
  # ---- 🔍 Identify ID Type & Species ----
  
  # Compute intersection counts to find the best match for the 'ID' column
  overlap_counts <- sapply(X = named_lists, 
                           FUN = function(x) { length(intersect(x, df$ID)) })
  
  best_match <- names(which.max(overlap_counts))
  
  log_info(sample = "",
           step   = "add_annotation",
           msg    = glue::glue("Auto-detected ID format: {best_match}\n",
                               "Overlap counts:\n",
                               "{paste(names(overlap_counts), overlap_counts, sep = ' \\t: ', collapse = '\\n')}"))
  
  # Determine ID Type (Ensembl vs Entrez)
  if (grepl(pattern = "ensembl", x = best_match)) {
    id_col     <- "ENSEMBL_ID"
    symbol_col <- "ENSEMBL_SYMBOL"
  } else if (grepl(pattern ="entrez", x = best_match)) {
    id_col     <- "ENTREZ_ID"
    symbol_col <- "ENTREZ_SYMBOL"
  }
  
  # Determine Species and Key Columns
  species <- if (grepl(pattern = "human", x = best_match)) "Homo sapiens" else "Mus musculus"
  ann_df      <- ann_list[[species]]
  
  # ---- 🤝 Join & Finalize Columns ----
  
  # Clean the annotation source FIRST
  # NOTE: If id_col is "ENTREZ_ID", then it could map to multiple "ENSEMBL_ID".
  # Since we used multiple = "all" in past to retain most info, duplicates occured after left_join().
  # So, now we first remove duplicates from ann_df based on id_col and then join.
  # This way final annotated df will be free from duplicates.
  ann_df_clean <- ann_df %>%
    # STRICT PRIORITY: symbol_col > Ensembl > Entrez > ID
    # If id_col = ENTREZ_ID, then 647042 and 643707 map to same ENSEMBL_SYMBOL "GOLGA6L10".
    # By prioritizing, symbol_col = ENTREZ_SYMBOL, they map to "GOLGA6L10" and "GOLGA6L4".
    dplyr::mutate(priority = dplyr::case_when(!is.na(.data[[symbol_col]]) & nzchar(.data[[symbol_col]]) ~ 0, # Top Priority
                                              !is.na(ENSEMBL_SYMBOL) & nzchar(ENSEMBL_SYMBOL)           ~ 1,
                                              !is.na(ENTREZ_SYMBOL)  & nzchar(ENTREZ_SYMBOL)            ~ 2,
                                              TRUE                                                      ~ 3)) %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::slice_min(order_by = priority, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(ENSEMBL_SYMBOL, ENTREZ_SYMBOL, all_of(id_col))
  
  # Join with input
  annotated_df <- df %>%
    dplyr::mutate(ID = as.character(ID)) %>%
    dplyr::left_join(ann_df_clean, by = stats::setNames(id_col, "ID"), multiple = "all") %>%
    dplyr::mutate(SYMBOL = dplyr::case_when(!is.na(.data[[symbol_col]]) & nzchar(.data[[symbol_col]]) ~ .data[[symbol_col]],
                                            !is.na(ENSEMBL_SYMBOL) & nzchar(ENSEMBL_SYMBOL)           ~ ENSEMBL_SYMBOL,
                                            !is.na(ENTREZ_SYMBOL)  & nzchar(ENTREZ_SYMBOL)            ~ ENTREZ_SYMBOL,
                                            TRUE                                                      ~ ID)) 
  
  # ---- 🧹 Cleanup & Organization ----
  
  if (remove_ann_col) {
    # Keep SYMBOL and your original numeric sample columns.
    # We drop ID and the helper symbol columns from the join.
    annotated_df <- annotated_df %>% 
      dplyr::select(SYMBOL, dplyr::everything(), -dplyr::any_of(c("ID", colnames(ann_df_clean))))
  } else {
    # Full mode: Put Symbol first, but keep everything else
    annotated_df <- annotated_df %>% 
      dplyr::select(SYMBOL, dplyr::everything())
  }
  
  # ---- 📊 Identify Ambiguous Mappings (Collapse Check) ----
  
  collapse_stats <- annotated_df %>%
    dplyr::count(SYMBOL) %>%
    dplyr::filter(n > 1)
  
  if (nrow(collapse_stats) > 0) {
    num_symbols <- nrow(collapse_stats)
    total_ids   <- sum(collapse_stats$n)
    
    log_warn(sample = "", 
             step   = "add_annotation", 
             msg    = glue::glue("Ambiguity Check: {total_ids} IDs collapsed into {num_symbols} unique Symbols. ",
                                 "Use their mean or sum in downstream analysis."))
  }
  
  # ---- 🪵 Log Output and Return ----
  
  log_info(sample = "", 
           step   = "add_annotation", 
           msg    = glue::glue("Annotations added. ID '{id_col}' converted to SYMBOL."))
  
  return(annotated_df)
}

fit_sva <- function(dds, condition_col) {
  
  # ---- ⚙️ Prepare Data ----
  # SVA requires normalized counts (transformed to stabilize variance)
  dat <- counts(dds, normalized = TRUE)
  
  # Filter out very low power genes for SVA speed (Standard Practice)
  idx  <- rowMeans(dat) > 1
  dat  <- dat[idx, ]
  
  # Create the Model Matrices
  # Full model: includes the biological variable of interest
  # Null model: ignores the biological variable
  mod  <- model.matrix(as.formula(glue::glue("~ {condition_col}")), colData(dds))
  mod0 <- model.matrix(~ 1, colData(dds))
  
  # ---- 🔍 Identify Surrogate Variables ----
  log_info(sample = "", step = "SVA", msg = "Estimating number of surrogate variables...")
  n.sv <- sva::num.sv(dat, mod, method = "leek")
  
  if (n.sv == 0) {
    log_warn(sample = "", step = "SVA", msg = "No surrogate variables found. Proceeding with standard model.")
    return(dds)
  }
  
  log_info(sample = "", step = "SVA", msg = glue::glue("Found {n.sv} surrogate variables. Estimating..."))
  svobj <- sva::sva(dat, mod, mod0, n.sv = n.sv)
  
  # ---- 🤝 Integrate into DDS ----
  # Add the SVs to the metadata (colData)
  for (i in 1:n.sv) {
    colData(dds)[[paste0("SV", i)]] <- svobj$sv[, i]
  }
  
  # Update the DESeq2 Design Formula to account for hidden batches
  # New formula: ~ SV1 + SV2 + ... + condition
  sv_names <- paste0("SV", 1:n.sv)
  new_formula <- as.formula(paste("~", paste(c(sv_names, condition_col), collapse = " + ")))
  design(dds) <- new_formula
  
  log_info(sample = "", step = "SVA", msg = glue::glue("Updated design formula to: {format(new_formula)}"))
  
  return(dds)
}

norm_counts_DESeq2 <- function(metadata, read_data, proj.params) {
  
  # Batch Correction (if applicable) 
  if ("Batch" %in% colnames(metadata) && length(unique(metadata$Batch)) > 1) {
    normalized_counts_batch <- limma::removeBatchEffect(x = log2(normalized_counts + 1),
                                                        batch = dds$Batch)
  }
  
  normalized_counts_batch_df <- normalized_counts_batch %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID") %>%
    add_annotation()

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

# ---- 🧬 SINGLE CELL & SPATIAL ANALYSIS RELATED FUNCTIONS ----

validate_inputs <- function(sample = NULL, 
                            gene.column = NULL, 
                            bin = NULL, 
                            assay = NULL, 
                            n_pcs = NULL, 
                            pN = NULL,
                            reference_samples = NULL,
                            features = NULL, 
                            s_genes = NULL, 
                            g2m_genes = NULL, 
                            res = NULL, 
                            reduction = NULL, 
                            filename = NULL,
                            seurat_object = NULL, 
                            seurat_list = NULL, 
                            metadata = NULL, 
                            matrix_dir = NULL, 
                            output_dir = NULL, 
                            counts_dir = NULL,
                            expr_mat = NULL,
                            dds = NULL,
                            res_df = NULL,
                            xlsx_file = NULL, 
                            logic_vars = NULL,
                            metadata_cols = NULL) {
  
  # Store the calling function as a string
  caller <- tryCatch(expr = as.character(sys.call(-1)[[1]]),
                     error = function(e) { "unknown_caller"})
  step <- paste0(caller, " : validate_inputs")
  
  
  # ---- ⚡ Check if user accidentally passed a undefined variable ----
  
  # Capture user-supplied args safely
  called_args <- as.list(match.call(expand.dots = FALSE))[-1]
  
  # Collect all undefined variables
  bad_vars <- character()
  
  # Iterate through each variable passed to the function
  for (var_name in names(called_args)) {
    
    val <- tryCatch(expr = eval(called_args[[var_name]], envir = parent.frame()),
                    error = function(e) { e })
    
    # Case 1: evaluation failed (object not defined)
    if (inherits(val, "error")) {
      bad_vars <- c(bad_vars, glue::glue("'{var_name}' is not defined in the calling environment"))
      next
    }
    
    # Case 2: user accidentally passed a function (e.g. metadata())
    if (is.function(val)) {
      bad_vars <- c(bad_vars, glue::glue("'{var_name}' refers to a function, not a user-defined object"))
    }
  }
  
  # Print all undefined variables
  if (length(bad_vars) > 0) {
    log_error(sample = sample %||% "",
              step   = step,
              msg    = paste("INPUT ERROR: Invalid arguments detected:\n",
                             paste0("  • ", bad_vars, collapse = "\n")))
  }
  
  # ---- 1️⃣ Basic Parameter Checks ----
  
  if (!is.null(sample) && (!is.character(sample) || length(sample) != 1 || nchar(sample) == 0)) {
    log_error(sample  = sample,
              step    = step, 
              msg =glue::glue("INPUT ERROR: 'sample' must be a single, non-empty character string."))
  }
  
  if (!is.null(gene.column) && (!is.numeric(gene.column) || length(gene.column) != 1 || !(gene.column %in% c(1,2)))) {
    log_error(sample  = sample,
              step    = step, 
              msg =glue::glue("INPUT ERROR: 'gene.column' must be 1 (Ensembl IDs) or 2 (Gene symbols)."))
  }
  
  if (!is.null(bin) && (!is.numeric(bin) || length(bin) != 1 || !(bin %in% c(2,8,16)))) {
    log_error(sample  = sample,
              step    = step, 
              msg =glue::glue("INPUT ERROR: Invalid 'bin' size provided. Must be one of 2, 8, or 16."))
  }
  
  if (!is.null(assay) && (!is.character(assay) || length(assay) != 1 || nchar(assay) == 0)) {
    log_error(sample  = sample,
              step    = step, 
              msg =glue::glue("INPUT ERROR: 'assay' must be a single, non-empty character string."))
  }
  
  if (!is.null(n_pcs) && (!is.numeric(n_pcs) || length(n_pcs) != 1 || n_pcs <= 0 || n_pcs %% 1 != 0 || n_pcs > 50)) {
    log_error(sample  = sample,
              step    = step, 
              msg =glue::glue("INPUT ERROR: 'n_pcs' must be a positive integer less than 50."))
  }
  
  if (!is.null(pN) && (!is.numeric(pN) || length(pN) != 1 || pN <= 0 || pN > 1)) {
    log_error(sample  = sample,
              step    = step, 
              msg =glue::glue("INPUT ERROR: 'pN' must be a numeric value between 0 and 1 (exclusive of 0)."))
  }
  
  if (!is.null(reduction) && (!is.character(reduction) || length(reduction) != 1 || nchar(reduction) == 0)) {
    log_error(sample  = sample,
              step    = step, 
              msg =glue::glue("INPUT ERROR: 'reduction' must be a single, non-empty character string."))
  }
  
  if (!is.null(reference_samples) && (!is.character(reference_samples) || length(reference_samples) == 0)) {
    log_error(sample  = sample,
              step    = step, 
              msg =glue::glue("INPUT ERROR: 'reference_samples' must be a non-empty character vector."))
  }
  
  # ---- 2️⃣ Feature Checks ----
  
  # Helper Function for Validation
  validate_gene_set <- function(gene_set, set_name, available_features) {
    
    if (is.null(gene_set)) return(invisible(NULL))
    
    existing_features <- base::intersect(gene_set, available_features)
    missing_features <- base::setdiff(gene_set, existing_features)
    
    if (length(existing_features) == 0) {
      log_error(sample  = sample,
                step    = step, 
                msg =glue::glue("INPUT ERROR: None of the '{length(gene_set)}' provided features were found in the Seurat object."))
    }
  }
  
  # Check features, s_genes, g2m_genes
  if (!is.null(seurat_object)) {
    
    # Collect ALL available features once (genes from default assay + metadata columns)
    available_features <- c(SeuratObject::Features(seurat_object),
                            colnames(seurat_object@meta.data))
    
    # Run validation for each set using the helper function
    if (!is.null(features)) {
      features <- validate_gene_set(features, "features", available_features)
    }
    
    if (!is.null(s_genes)) {
      s_genes <- validate_gene_set(s_genes, "s_genes", available_features)
    }
    
    if (!is.null(g2m_genes)) {
      g2m_genes <- validate_gene_set(g2m_genes, "g2m_genes", available_features)
    }
  }
  
  # ---- 3️⃣ 'res' Check ----
  
  if (!is.null(res)) {
    allowed_res <- c(0.2,0.4,0.6,0.8,1.0,1.2)
    resolution_num <- as.numeric(res)
    if (is.na(resolution_num) || !(resolution_num %in% allowed_res)) {
      log_error(sample  = sample,
                step    = step, 
                msg =glue::glue("INPUT ERROR: Specified res '{res}' is NOT one of the **valid** res : '{paste(allowed_res, collapse = ", ")}'."))
    }
  }
  
  # ---- 4️⃣ 'filename' Check ----
  
  if (!is.null(filename)) {
    if (!is.character(filename) || length(filename) != 1 || nchar(filename) == 0) {
      log_error(sample  = sample,
                step    = step, 
                msg =glue::glue("'filename' must be a single, non-empty string."))
    }
    if (grepl("[<>:\"/\\\\|?*]", filename)) {
      log_error(sample  = sample,
                step    = step, 
                msg =glue::glue("'filename' contains illegal characters for file paths."))
    }
  }
  
  # ---- 5️⃣ Seurat Object and List Integrity Checks ----
  
  # Check Seurat object class (only if provided)
  if (!is.null(seurat_object)) {
    if (!inherits(seurat_object, "Seurat")) {
      log_error(sample  = sample,
                step    = step, 
                msg =glue::glue("INPUT ERROR: 'seurat_object' must be a Seurat object."))
    }
    
    # Check for Non-Empty Seurat Object
    if (ncol(seurat_object) == 0) {
      log_error(sample  = sample,
                step    = step, 
                msg =glue::glue("DATA ERROR: The input Seurat object is empty (0 cells)."))
    }
    
    if (!is.null(assay)) {
      # First, check if the string is structurally sound (Done in Section 1)
      # Second, check if the assay name exists in the object:
      if (!(assay %in% names(seurat_object@assays))) {
        log_error(sample  = sample,
                  step    = step, 
                  msg =glue::glue("ASSAY ERROR: Specified assay '{assay}' is NOT one of the **valid assays** : '{paste(names(seurat_object@assays), collapse = ', ')}'."))
      }
    }
  }
  
  # Check 'seurat_list' structure and naming (only if provided)
  if (!is.null(seurat_list)) {
    if (!is.list(seurat_list) || length(seurat_list) == 0) {
      log_error(sample  = sample,
                step    = step, 
                msg =glue::glue("INPUT ERROR: 'seurat_list' must be a non-empty list of Seurat objects."))
    }
    if (is.null(names(seurat_list)) || any(nchar(names(seurat_list)) == 0)) {
      log_error(sample  = sample,
                step    = step, 
                msg =glue::glue("INPUT ERROR: 'seurat_list' must be a *named* list with non-empty names."))
    }
    # Check all elements are Seurat objects
    if (!all(vapply(seurat_list, inherits, logical(1), what = "Seurat"))) {
      log_error(sample  = sample,
                step    = step, 
                msg =glue::glue("INPUT ERROR: All elements in 'seurat_list' must be Seurat objects."))
    }
  }
  
  # ---- 6️⃣ 'metadata' Check ----
  
  if (!is.null(metadata)) {
    if (!is.data.frame(metadata)) {
      log_error(sample  = sample,
                step    = step, 
                msg =glue::glue("INPUT ERROR: 'metadata' must be a data.frame."))
    }
    if (nrow(metadata) == 0) {
      log_error(sample  = sample,
                step    = step, 
                msg =glue::glue("DATA ERROR: 'metadata' must be a non-empty data.frame (0 rows detected)."))
    }
  }
  
  # ---- 7️⃣ 'matrix_dir' Check ----
  
  if (!is.null(matrix_dir) && !is.null(sample)) {
    
    data_dir <- base::file.path(matrix_dir, sample)
    
    # Check for existence of the specific data directory
    if (!dir.exists(data_dir)) {
      log_error(sample  = sample,
                step    = step, 
                msg =glue::glue("PATH ERROR: The data directory '{data_dir}' does not exist!"))
    }
    
    # Check for required file structure (HDF5 OR 3-file structure - for CellRanger)
    h5_files <- list.files(data_dir, pattern = "\\.h5$", full.names = TRUE, ignore.case = TRUE)
    expected_files_mtx <- c("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")
    standard_files_found <- all(sapply(file.path(data_dir, expected_files_mtx), file.exists))
    
    if (length(h5_files) == 0 && !standard_files_found) {
      log_error(sample  = sample,
                step    = step, 
                msg =glue::glue("FILE STRUCTURE ERROR: Directory '{data_dir}' does not contain a valid 10X matrix."))
    } else if (length(h5_files) > 1) {
      log_warn(sample  = sample,
               step    = step, 
               msg =glue::glue("⚠️ FILE STRUCTURE WARNING: Multiple HDF5 files found. Only using the first one : '{h5_files[1]}'."))
    }
    
    # # Check for REQUIRED Space Ranger files (for load_spaceranger)
    # required_spatial_files <- c("filtered_feature_bc_matrix.h5", "tissue_positions_list.csv", "tissue_lowres_image.png")
    # required_files_exist <- all(sapply(paste0(data_dir, "/", required_spatial_files), file.exists))
    # 
    # if (!required_files_exist) {
    #   stop("FILE STRUCTURE ERROR: The directory '", data_dir, "' is missing one or more of the required Space Ranger files: ", paste(required_spatial_files, collapse = ", "), ".")
    # }
  }
  
  # ---- 8️⃣ 'output_dir' Check ----
  
  if (!is.null(output_dir)) {
    
    # Check 1: Directory must exist
    if (!dir.exists(output_dir)) {
      log_warn(sample  = sample,
               step    = step, 
               msg =glue::glue("PATH WARN: Output directory '{output_dir}' does not exist.
                                     Attempting to create."))
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Check 2: Directory must be writable
    # file.access(path, mode=2) checks for write permission. Returns 0 on success.
    if (file.access(output_dir, 2) != 0) {
      log_error(sample  = sample,
                step    = step, 
                msg =glue::glue("PATH ERROR: Output directory '{output_dir}' is NOT writable."))
    }
  }
  
  # ---- 9️⃣ 'xlsx file' Check ---- 
  
  if (!is.null(xlsx_file)) {
    
    # Check 1: File must exist
    if (!file.exists(xlsx_file)) {
      log_error(sample  = sample,
                step    = step, 
                msg =glue::glue("FILE ERROR: External metadata file not found at :  '{xlsx_file}'."))
    }
    
    # Check 2: File must be in xlsx format
    if (!grepl("\\.xlsx$", xlsx_file, ignore.case = TRUE)) {
      log_error(sample  = sample,
                step    = step, 
                msg =glue::glue("FILE ERROR: 'metafile' must be an .xlsx file."))
    }
  }
  
  # ---- 🔟 'metadata_cols' Check ----
  
  if (!is.null(seurat_object)) {
    
    required_cols <- unlist(metadata_cols)
    required_cols <- required_cols[!sapply(required_cols, is.null)]
    
    if (length(required_cols) > 0) {
      missing_cols <- base::setdiff(required_cols, colnames(seurat_object@meta.data))
      
      if (length(missing_cols) > 0) {
        log_error(sample  = sample,
                  step    = step, 
                  msg =glue::glue("INPUT ERROR: Required column(s) '{paste(missing_cols, collapse = ", ")}' NOT found in metadata.")) 
      }
    }
  }
  
  # ---- 1️⃣1️⃣ 'counts_dir' Check ----
  
  if (!is.null(counts_dir)) {
    if (!dir.exists(counts_dir)) {
      log_warn(sample  = sample,
               step    = step, 
               msg =glue::glue("Count directory '{counts_dir}' does not exist. 
                                    Provide raw counts as excel for analysis."))
    }
  }
  
  # ---- 1️⃣2️⃣ 'expr_mat' Check ----
  
  if (!is.null(expr_mat)) {
    if (is.data.frame(expr_mat)) {
      log_warn(sample = sample,
               step   = step,
               msg    = "`expr_mat` is a data.frame and will be converted to a numeric matrix.")
      
    } else if (!is.matrix(expr_mat)) {
      log_error(sample = sample,
                step   = step,
                msg    = "`expr_mat` must be a matrix or data.frame.")
    }
    
    expr_mat_numeric <- apply(X = expr_mat, MARGIN = 2, FUN = function(x){as.numeric(as.character(x))})
    if (any(is.na(expr_mat_numeric))) {
      log_error(sample = sample,
                step   = step,
                msg    = "`expr_mat` contains non-numeric columns that cannot be converted safely.")
    }
    
    if (is.null(rownames(expr_mat))) {
      log_error(sample = sample,
                step   = step, 
                msg    = "`expr_mat` must have rownames representing genes.")
    }
    
    if (is.null(colnames(expr_mat))) {
      log_error(sample = sample,
                step   = step, 
                msg    = "`expr_mat` must have colnames representing sample names.")
    }
    
    if (any(duplicated(colnames(expr_mat)))) {
      log_error(sample = sample,
                step   = step, 
                msg    = "Duplicate sample names detected in expr_mat.")
    }
  }
  
  # ---- 1️⃣3️⃣ 'logic_vars' Check ----
  
  if (!is.null(logic_vars)) {
    for (i in seq_along(logic_vars)) {
      val <- logic_vars[[i]]
      logic_var <- names(logic_vars)[i] %||% paste0("logic_var_", i)  # fallback name
      
      if (!is.logical(val) || length(val) != 1) {
        log_error(sample = sample,
                  step   = "validate_inputs",
                  msg    = glue::glue("INPUT ERROR: `{logic_var}` must be a single logical (TRUE/FALSE)."))
      }
    }
  }
  
  # ---- 1️⃣4️⃣ 'DESeq2 object' Check ----
  
  if (!is.null(dds)) {
    if (!inherits(dds, "DESeqDataSet")) {
      log_error(sample = sample,
                step   = "validate_inputs",
                msg    = "`dds` must be a DESeqDataSet object.")
    }
  }
  
  # ---- 1️⃣5️⃣ `res_df` Check ----
  
  if (!is.null(res_df)) {
    
    # (i) Data frame
    if (!is.data.frame(res_df)) {
      log_error(sample = "",
                step   = "plot_volcano",
                msg    = "`res_df` must be a data.frame.")
    }
    
    # (ii) Required columns
    required_cols <- c("log2FoldChange", "padj", "SYMBOL")
    missing_cols  <- base::setdiff(required_cols, colnames(res_df))
    
    if (length(missing_cols) > 0) {
      log_error(sample = "",
                step   = "plot_volcano",
                msg    = glue::glue("`res_df` is missing required column(s): '{paste(missing_cols, collapse = ", ")}'."))
    }
    
    # (iii) Required column types
    if (!is.numeric(res_df$log2FoldChange)) {
      log_error(sample = "",
                step   = "plot_volcano",
                msg    = "`res_df$log2FoldChange` must be numeric.")
    }
    
    if (!is.numeric(res_df$padj)) {
      log_error(sample = "",
                step   = "plot_volcano",
                msg    = "`res_df$padj` must be numeric.")
    }
    
    if (!is.character(res_df$SYMBOL)) {
      log_error(sample = "",
                step   = "plot_volcano",
                msg    = "`res_df$SYMBOL` must be a character vector.")
    }
  }
  
  # ---- 🪵 Log Output and Return ----
  
  log_info(sample = "", 
           step = "validate_inputs",
           msg =glue::glue("Successfully validated input variables for function : '{sample}'."))
  
  return(invisible(NULL))
  
}

load_cellranger <- function(sample, matrix_dir, gene.column = 2){
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(sample = sample, matrix_dir = matrix_dir, gene.column = gene.column)
  
  # ---- 📥 Read the Feature-Barcode Matrix ----
  
  # gene.column = 1 → Ensembl IDs
  # gene.column = 2 → Gene symbols (required for mito/ribo/heme ratios)
  
  data_dir <- base::file.path(matrix_dir, sample)
  
  # Look for HDF5 files in the directory
  h5_files <- list.files(data_dir, pattern = "\\.h5$", full.names = TRUE, ignore.case = TRUE)
  
  counts <- NULL
  
  # If HDF5 files are found, try to read them
  if (length(h5_files) > 0) {
    log_info(sample = sample, 
             step = "load_cellranger",
             msg =glue::glue("Reading HDF5 matrix from directory."))
    
    counts <- tryCatch({
      # Read the gene expression matrix from the HDF5 file
      Seurat::Read10X_h5(filename = h5_files[1],
                         use.names = TRUE,
                         unique.features = TRUE)
    }, error = function(e) {
      log_error(sample = sample, 
                step = "load_cellranger",
                msg =glue::glue("Failed to read HDF5 matrix from '{h5_files[1]}'."))
    })
    
    # Error check: Ensure the read returns a matrix, not a list
    if (is.list(counts)) {
      log_error(sample = sample, 
                step = "load_cellranger",
                msg ="HDF5 file returned a list. Check the file content.")
    }
    
  } else {
    # If no HDF5 files, fall back to the standard 10X directory
    log_info(sample = sample, 
             step = "load_cellranger",
             msg ="Reading standard 10X matrix from directory.")
    
    counts <- tryCatch({
      Seurat::Read10X(data.dir = data_dir,
                      gene.column = gene.column, # Parameterized (Default: 2 for gene symbols)
                      cell.column = 1,
                      unique.features = TRUE,
                      strip.suffix = FALSE
      )
    }, error = function(e) {
      log_error(sample = sample, 
                step = "load_cellranger",
                msg =glue::glue("Failed to read 10X matrix from '{data_dir}'."))
    })
  }
  
  # ---- 🏗️ Create Seurat Object ----
  
  sample_seurat <- SeuratObject::CreateSeuratObject(counts = counts,
                                                    project = sample,
                                                    assay = "RNA",
                                                    names.field = 1,
                                                    names.delim = "_",
                                                    meta.data = NULL,
                                                    min.cells = 0,
                                                    min.features = 0)
  
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  log_info(sample = sample, 
           step = "load_cellranger",
           msg =glue::glue("Successfully created Seurat object for sample : '{step}'."))
  return(invisible(sample_seurat))
}

load_spaceranger <- function(sample, bin, matrix_dir){
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(sample = sample, bin = bin, matrix_dir = matrix_dir)
  
  # ---- 📥 Load Seurat Object with Filtered Matrix ----
  
  data_dir <- base::file.path(matrix_dir, sample)
  
  # Load10X_Spatial automatically finds the necessary files (images, coordinates, matrix)
  sample_seurat <- tryCatch({
    Seurat::Load10X_Spatial(data.dir = data_dir,
                            filename = "filtered_feature_bc_matrix.h5",
                            assay = "Spatial",
                            slice = sample,
                            bin.size = bin,
                            filter.matrix = TRUE,
                            to.upper = FALSE,
                            image = NULL)
  }, error = function(e) {
    pad <- base::strrep(" ", nchar(sample) + 30)
    log_error(sample = sample, 
              step = "load_spaceranger",
              msg =glue::glue("Failed to load binned data from '{data_dir}' for bin size '{bin}'.\n",
                              "{pad}Please verify the Space Ranger output structure."))
  })
  
  # Add sample identifier to Seurat object
  sample_seurat[["orig.ident"]] <- sample
  
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  log_info(sample = sample, 
           step = "load_spaceranger",
           msg =glue::glue("Successfully created Seurat object for sample : '{sample}' ('{bin}' bin)."))
  return(invisible(sample_seurat))
}

classify_dropletutils <- function(sample_seurat){ 
  
  # Set seed for reproducible stochastic processes (emptyDrops)
  set.seed(100)
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = sample_seurat)
  
  sample <- as.character(unique(sample_seurat@meta.data$orig.ident))
  
  # ---- 🧪 Iterative emptyDrops Run ----
  
  # Convert Seurat Object to SingleCellExperiment (MUST contains all barcodes)
  sce <- Seurat::as.SingleCellExperiment(x = sample_seurat)
  
  # Initialize Parameters for Iterative emptyDrops
  niters <- 10000
  
  # Wrapper to safely run emptyDrops
  df <- tryCatch({
    
    # Run emptyDrops to classify droplets as empty or non-empty
    e.out <- DropletUtils::emptyDrops(m = SingleCellExperiment::counts(sce),
                                      niters = niters)
    
    # Convert emptyDrops Output to Data Frame for easier filtering
    df <- as.data.frame(e.out)
    
    # Loop to handle limited droplets if needed
    # NOTE: If FDR > 0.05 AND Limited == TRUE, the FDR of these droplets can be 
    # reduced with more iterations.
    n_improve <- nrow(df %>% filter(Limited == TRUE, FDR > 0.05))
    while (n_improve > 0) {
      niters <- niters + 10000
      e.out <- DropletUtils::emptyDrops(m = SingleCellExperiment::counts(sce),
                                        niters = niters)
      df <- as.data.frame(e.out)
      n_improve <- nrow(df %>% filter(Limited == TRUE, FDR > 0.05))
      pad <- base::strrep(" ", nchar(sample) + 35)
      log_info(sample = sample, 
               step = "classify_dropletutils",
               msg =glue::glue("emptyDrops check: '{n_improve}' droplets need more iterations.\n",
                               "{pad}Current niters = '{niters}'."))
    }
    
    df
  }, error = function(e) {
    
    # Catch the specific ambient profile error
    if (grepl("no counts available to estimate the ambient profile", e$message)) {
      log_warn(sample = sample, 
               step = "classify_dropletutils",
               msg =paste("emptyDrops failed: no counts available for ambient profile.\n",
                          base::strrep(" ", 35), "Setting all barcodes as non-empty (FDR < 0.05)."))
      
      # Create dummy dataframe with FDR < 0.05
      data.frame(Total = colSums(SingleCellExperiment::counts(sce)),
                 LogProb = NA,
                 PValue = NA,
                 FDR = 0,          # Mark all as significant / non-empty
                 Limited = FALSE,
                 row.names = colnames(sce))
    } else {
      # Re-throw other errors
      log_error(sample = sample, 
                step = "classify_dropletutils",
                msg =glue::glue("Caught error ({class(e)[1]}): {e$message}")) 
    }
  })
  
  # Identify true (non-empty) cells with FDR ≤ 0.05
  # NOTE: Filtering out NA ensures only tested barcodes are considered.
  true_barcodes <- df %>%
    dplyr::filter(FDR <= 0.05, !is.na(FDR)) %>% 
    rownames()
  
  # ---- 🔍 Build Classification Vector ----
  
  # All barcodes present in the Seurat object — includes empty droplets
  all_barcodes <- colnames(sample_seurat)
  
  # Vectorized classification
  classification_vector <- dplyr::case_when(all_barcodes %in% true_barcodes ~ "Non-Empty Droplet",
                                            TRUE ~ "Empty Droplet")
  
  # Assign names so Seurat maps the metadata correctly by barcode
  names(classification_vector) <- all_barcodes
  
  # ---- 📥 Add Metadata to Seurat Object ----
  
  # Add the classification vector to the original Seurat object
  sample_seurat <- SeuratObject::AddMetaData(object   = sample_seurat,
                                             metadata = classification_vector,
                                             col.name = "DropletUtils")
  
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  log_info(sample = sample, 
           step = "classify_dropletutils",
           msg = glue::glue("DropletUtils classification complete for sample : '{sample}'."))
  return(invisible(sample_seurat))
}

classify_cellranger <- function(sample_seurat, filt_matrix_dir){
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = sample_seurat, matrix_dir = filt_matrix_dir)
  
  # ---- 📥 Read Filtered Matrix for Barcodes ----
  
  # Get sample name from sample_seurat
  sample <- as.character(unique(sample_seurat@meta.data$orig.ident))
  
  # Read the filtered matrix for the sample
  sample_seurat_filt <- load_cellranger(sample, matrix_dir = filt_matrix_dir)
  
  # ---- 🔍 Build Classification Vector ----
  
  # All barcodes present in the Seurat object — includes empty droplets
  all_barcodes <- colnames(sample_seurat)
  
  # Barcodes identified as 'cells' by Cell Ranger's internal algorithm
  filtered_barcodes <- colnames(sample_seurat_filt)
  
  # Vectorized classification
  classification_vector <- dplyr::case_when(all_barcodes %in% filtered_barcodes ~ "Non-Empty Droplet",
                                            TRUE ~ "Empty Droplet")
  
  # Assign names so Seurat maps the metadata correctly by barcode
  names(classification_vector) <- all_barcodes
  
  # ---- 📥 Add Metadata to Seurat Object ----
  
  # Add the classification vector to the original Seurat object
  sample_seurat <- SeuratObject::AddMetaData(object   = sample_seurat,
                                             metadata = classification_vector,
                                             col.name = "CellRanger")
  
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  log_info(sample = sample, 
           step = "classify_cellranger",
           msg = glue::glue("CellRanger empty droplets identified for sample : '{sample}'."))
  return(invisible(sample_seurat))
}

calc_qc_metrics <- function(sample_seurat, assay){
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = sample_seurat, assay = assay)
  
  sample <- as.character(unique(sample_seurat@meta.data$orig.ident))
  
  # ---- 📥 Populate missing optional metadata columns ----
  
  # NOTE: Use "ND" instead of NA for safer subsetting
  optional_cols <- c("DropletUtils", "CellRanger")
  for (col in optional_cols) {
    if (!col %in% colnames(sample_seurat@meta.data)) {
      sample_seurat <- Seurat::AddMetaData(object   = sample_seurat,
                                           metadata = "ND",
                                           col.name = col)
    }
  }
  
  # ---- 🔍 Calculate Percentage Feature Sets ----
  
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
  
  # ---- 📥 Pull Metadata and Compute QC Metrics ----
  
  # (i)   Cell       : unique identifiers corresponding to each cell i.e. barcodes
  # (ii)  Sample     : sample names
  # (iii) nUMIs      : number of transcripts per cell
  # (iv)  nGenes     : number of genes per cell
  # (v)   nHTO_UMIs  : number of HTO reads per cell
  # (vi)  nHTOs      : number of HTO types per cell
  # (vii) MitoRatio   : MitoPercent / 100
  # (viii) RiboRatio  : RiboPercent / 100
  # (ix)  HemeRatio   : HemePercent / 100
  # (x)   Novelty     : log ratio of genes per UMI
  
  sample_metadata <- sample_seurat@meta.data %>%
    tibble::rownames_to_column(var = "barcode") %>%
    dplyr::mutate(Cell = paste0(orig.ident, "_", barcode),
                  Sample = orig.ident,
                  nUMIs = .data[[paste0("nCount_", assay)]],
                  nGenes = .data[[paste0("nFeature_", assay)]],
                  MitoRatio = MitoPercent / 100,
                  RiboRatio = RiboPercent / 100,
                  HemeRatio = HemePercent / 100,
                  Novelty = log10(nGenes) / log10(nUMIs))
  
  # ---- 📥 Handle HTO Metadata (if present) ----
  
  if ("nCount_HTO" %in% colnames(sample_metadata)){
    sample_metadata <- sample_metadata %>%
      dplyr::rename(nHTO_UMIs = nCount_HTO,
                    nHTOs = nFeature_HTO)
  }
  
  # ---- 📥 Add Spatial Coordinates (if available) ----
  
  if (length(sample_seurat@images) > 0){
    
    # Collect spatial coordinates into a single data frame
    df_coords <- data.frame()
    
    for (image_name in names(sample_seurat@images)){
      df <- data.frame(barcode = sample_seurat@images[[image_name]]@boundaries$centroids@cells,
                       X = sample_seurat@images[[image_name]]@boundaries$centroids@coords[,1],
                       Y = sample_seurat@images[[image_name]]@boundaries$centroids@coords[,2])
      df_coords <- dplyr::bind_rows(df_coords, df)
    }
    
    # Join coordinates based on barcode
    sample_metadata <- sample_metadata %>%
      dplyr::left_join(df_coords, by=c("barcode"="barcode"))
  }
  
  # ---- 🔍 Define QC Cutoffs and Classify Cells ----
  
  is_spatial <- length(names(sample_seurat@images)) > 0
  
  # Define cutoffs based on spatial or single-cell conditions
  if (is_spatial) {
    gene_cutoff     <- 25
    umi_cutoff      <- 50
    mito_cutoff     <- 0.2
    novelty_cutoff  <- 0.8
    ribo_cutoff     <- 0.05
  } else {
    gene_cutoff     <- 250
    umi_cutoff      <- 500
    mito_cutoff     <- 0.2
    novelty_cutoff  <- 0.8
    ribo_cutoff     <- 0.05
  }
  
  # Log the applied cutoffs
  pad <- base::strrep(" ", nchar(sample) + 46)
  log_info(sample = sample, 
           step = "calc_qc_metrics",
           msg = glue::glue("Cutoffs applied: nGenes ≥ {gene_cutoff}\n",
                            "{pad}nUMIs ≥ {umi_cutoff}\n",
                            "{pad}MitoRatio ≤ {mito_cutoff}\n",
                            "{pad}Novelty ≥ {novelty_cutoff}"))
  
  # Create Quality column, drop 'barcode' column and restore it as rownames
  sample_metadata <- sample_metadata %>%
    dplyr::mutate(passes_qc = nGenes  >= gene_cutoff & 
                    nUMIs   >= umi_cutoff &
                    MitoRatio <= mito_cutoff &
                    Novelty   >= novelty_cutoff,
                  Quality = dplyr::case_when(
                    DropletUtils == "Empty Droplet" & CellRanger == "Empty Droplet" ~ "Empty Droplet",
                    !passes_qc                                                      ~ "Low Quality",
                    TRUE                                                            ~ "High Quality")) %>%
    tibble::column_to_rownames(var = "barcode")
  
  # ---- 📥 Add Metadata to Seurat Object ----
  
  # Add the classification vector to the original Seurat object
  sample_seurat <- SeuratObject::AddMetaData(object   = sample_seurat,
                                             metadata = sample_metadata)
  
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  
  log_info(sample = sample, 
           step = "calc_qc_metrics",
           msg = glue::glue("Cell-level QC metrics calculated for sample : '{sample}'."))
  return(invisible(sample_seurat))
}

classify_doubletfinder <- function(sample_seurat, n_pcs = NULL, pN = 0.25){
  
  # Set seed for reproducible stochastic processes (PCA, UMAP, Clustering, DoubletFinder)
  set.seed(100)
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = sample_seurat, n_pcs = n_pcs, pN = pN)
  
  sample <- as.character(unique(sample_seurat@meta.data$Sample))
  
  # Ensure 'Quality' column is present in metadata
  required_cols <- c("Quality")
  if (!all(required_cols %in% colnames(sample_seurat@meta.data))) {
    log_error(sample = sample, 
              step = "classify_doubletfinder",
              msg ="Missing 'Quality' column in metadata. Please run calc_qc_metrics().")
  }
  
  # ---- 🔍 Filter Out Empty Droplets and Low Quality Cells ----
  
  subset_seurat <- subset(x = sample_seurat,
                          subset = (Quality == "High Quality"))
  
  # Ensure enough cells remain
  if (ncol(subset_seurat) < 50) {
    log_warn(sample = sample, 
             step = "classify_doubletfinder",
             msg ="Fewer than 50 cells remain after filtering empty droplets. Skipping DoubletFinder.")
    return(invisible(sample_seurat))
  }
  
  # ---- 🔍 Preprocessing Required for DoubletFinder ----
  
  subset_seurat <- subset_seurat %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    Seurat::FindVariableFeatures(verbose = FALSE) %>%
    Seurat::ScaleData(verbose = FALSE) %>%
    Seurat::RunPCA(npcs = 50, verbose = FALSE)  # Run enough PCs for the sweep/clustering
  
  # ---- 🔍 Find Optimal pK ----
  
  # Artificial doublets are introduced into the real dataset at varying proportions.
  # Data is then preprocessed and the proportion of artificial nearest neighbors 
  # (pANN) is computed for multiple combinations of pK and pN.
  
  # The optimal pK corresponds to the value yielding the maximum bimodality coefficient 
  # (BCmvn), indicating the strongest separation between real and artificial doublets.
  
  # Find significant PCs
  if (is.null(n_pcs)) {
    stdev_pc <- subset_seurat@reductions$pca@stdev
    percent_stdev_pc <- (stdev_pc / sum(stdev_pc)) * 100
    cumulative_stdev_pc <- cumsum(percent_stdev_pc)
    
    pc1 <- which(cumulative_stdev_pc > 90 & percent_stdev_pc < 5)[1]
    pc2 <- sort(which((percent_stdev_pc[1:(length(percent_stdev_pc)-1)] - 
                         percent_stdev_pc[2:length(percent_stdev_pc)]) > 0.1),
                decreasing = TRUE)[1] + 1
    n_pcs <- min(pc1, pc2)
  }
  
  # Run UMAP, FindNeighbors, and FindClusters is only needed for the sweep
  subset_seurat <- subset_seurat %>%
    Seurat::RunUMAP(dims = 1:n_pcs, verbose = FALSE) %>%
    Seurat::FindNeighbors(dims = 1:n_pcs, verbose = FALSE) %>%
    Seurat::FindClusters(res = 0.1, verbose = FALSE)
  
  # Run paramSweep
  sweep.res <- quiet_msg(DoubletFinder::paramSweep(subset_seurat, PCs = 1:n_pcs, sct = FALSE))
  sweep.stats <- quiet_msg(DoubletFinder::summarizeSweep(sweep.res, GT = FALSE))
  bcmvn <- quiet_msg(DoubletFinder::find.pK(sweep.stats))
  
  # Find optimal pK using base R functions (safer than dplyr)
  optimal_pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
  optimal_pK <- as.numeric(as.character(optimal_pK))  # Ensure numeric
  log_info(sample = sample, 
           step = "classify_doubletfinder",
           msg =glue::glue("Optimal pK found : '{optimal_pK}'."))
  
  # ---- 🔍 Calculate Adjusted nExp (Expected Doublets) ----
  
  # https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled
  # From the 10x Genomics documentation, the multiplet rate is ~8*10^-6 per cell
  n_cells <- nrow(subset_seurat@meta.data)
  multiplet_rate_per_cell <- 8e-6               # Increment in doublet fraction per additional recovered cell
  multiplet_fraction <- 8e-6 * n_cells          # Fraction of cells expected to be doublets
  n_exp <- round(multiplet_fraction * n_cells)  # Number of doublets
  
  # Adjust for homotypic doublets
  homotypic.prop <- DoubletFinder::modelHomotypic(subset_seurat@meta.data$seurat_clusters)
  n_exp_adj <- round(n_exp * (1 - homotypic.prop))
  log_info(sample = sample, 
           step = "classify_doubletfinder",
           msg =glue::glue("Expected doublets (nExp_adj) : '{n_exp_adj}'."))
  
  # ---- 🔍 Run DoubletFinder ----
  
  subset_seurat <- quiet_msg(DoubletFinder::doubletFinder(seu = subset_seurat,
                                                          PCs = 1:n_pcs,
                                                          pN = pN,  #0.25 default
                                                          pK = optimal_pK,
                                                          nExp = n_exp_adj))
  
  # Extract DoubletFinder classification column
  df_col <- grep("^DF.classifications", colnames(subset_seurat@meta.data), value = TRUE)
  
  if (length(df_col) != 1) {
    log_error(sample = sample, 
              step = "classify_doubletfinder",
              msg =paste("More than 1 DoubletFinder classification columns present.",
                         "Cannot uniquely identify the DoubletFinder classification column.",
                         "Please check if DoubletFinder was already run.",
                         sep = "\n"))
  }
  
  # ---- 🔍 Build Classification Vector ----
  
  # All barcodes present in the Seurat object — includes empty droplets
  all_barcodes <- colnames(sample_seurat)
  
  # Barcodes identified as 'low quality' as per cutoffs
  quality_calls <- sample_seurat@meta.data[["Quality"]]
  low_quality_barcodes <- rownames(sample_seurat@meta.data)[quality_calls == "Low Quality"]
  
  # Barcodes identified as 'singlets' and 'doublets' by DoubletFinder's algorithm
  df_calls <- subset_seurat@meta.data[[df_col]]
  singlet_barcodes <- rownames(subset_seurat@meta.data)[df_calls == "Singlet"]
  doublet_barcodes <- rownames(subset_seurat@meta.data)[df_calls == "Doublet"]
  
  # Vectorized classification
  classification_vector <- dplyr::case_when(all_barcodes %in% low_quality_barcodes ~ "Low Quality",
                                            all_barcodes %in% doublet_barcodes ~ "Doublet",
                                            all_barcodes %in% singlet_barcodes ~ "Singlet",
                                            TRUE ~ "Empty Droplet")
  
  # Assign names so Seurat maps the metadata correctly by barcode
  names(classification_vector) <- all_barcodes
  
  # ---- 📥 Add Metadata to Seurat Object ----
  
  # Add the classification vector to the original Seurat object
  sample_seurat <- Seurat::AddMetaData(object = sample_seurat,
                                       metadata = classification_vector,
                                       col.name = "DoubletFinder")
  
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  log_info(sample = sample, 
           step = "classify_doubletfinder",
           msg = glue::glue("DoubletFinder doublets identified for sample : '{sample}'."))
  return(invisible(sample_seurat))
}

classify_scdblfinder <- function(sample_seurat){
  
  # Set seed for reproducible stochastic processes 
  set.seed(100)
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = sample_seurat)
  
  sample <- as.character(unique(sample_seurat@meta.data$Sample))
  
  # Ensure 'Quality' column is present in metadata
  required_cols <- c("Quality")
  if (!all(required_cols %in% colnames(sample_seurat@meta.data))) {
    stop("Missing 'Quality' classification columns. Please run calc_qc_metrics first.")
  }
  
  # ---- 🔍 Filter Out Empty Droplets and Low Quality Cells ----
  
  subset_seurat <- subset(x = sample_seurat,
                          subset = (Quality == "High Quality"))
  
  # Ensure enough cells remain
  if (ncol(subset_seurat) < 50) {
    log_warn(sample = sample, 
             step = "classify_scdblfinder",
             msg ="Fewer than 50 cells remain after filtering empty droplets. Skipping scDblFinder.")
    return(invisible(sample_seurat))
  }
  
  # ---- 🔍 Run scDblFinder ----
  
  # Convert Seurat Object to SingleCellExperiment
  sce <- Seurat::as.SingleCellExperiment(x = subset_seurat)
  
  sce <- quiet_msg(scDblFinder::scDblFinder(sce = sce,
                                            clusters = NULL,
                                            samples = NULL,
                                            dbr = NULL))
  
  # ---- 🔍 Build Classification Vector ----
  
  # All barcodes present in the Seurat object — includes empty droplets
  all_barcodes <- colnames(sample_seurat)
  
  # Barcodes identified as 'low quality' as per cutoffs
  quality_calls <- sample_seurat@meta.data[["Quality"]]
  low_quality_barcodes <- rownames(sample_seurat@meta.data)[quality_calls == "Low Quality"]
  
  # Barcodes identified as 'singlets' and 'doublets' by scDblFinder's algorithm
  df_calls <- base::tolower(as.character(sce$scDblFinder.class))   # scDblFinder classification column
  singlet_barcodes <- colnames(sce)[df_calls == "singlet"]
  doublet_barcodes <- colnames(sce)[df_calls == "doublet"]
  
  # Vectorized classification
  classification_vector <- dplyr::case_when(all_barcodes %in% low_quality_barcodes ~ "Low Quality",
                                            all_barcodes %in% doublet_barcodes ~ "Doublet",
                                            all_barcodes %in% singlet_barcodes ~ "Singlet",
                                            TRUE ~ "Empty Droplet")
  
  # Assign names so Seurat maps the metadata correctly by barcode
  names(classification_vector) <- all_barcodes
  
  # ---- 📥 Add Metadata to Seurat Object ----
  
  # Add the classification vector to the original Seurat object
  sample_seurat <- Seurat::AddMetaData(object = sample_seurat,
                                       metadata = classification_vector,
                                       col.name = "scDblFinder")
  
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  log_info(sample = sample, 
           step = "classify_scdblfinder",
           msg =glue::glue("scDblFinder doublets identified for sample : '{sample}'."))
  return(invisible(sample_seurat))
}

filter_singlets <- function(sample_seurat){
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = sample_seurat)
  
  sample <- as.character(unique(sample_seurat@meta.data$Sample))
  
  # ---- 📥 Populate missing optional metadata columns ----
  
  # NOTE: Use "ND" instead of NA for safer subsetting
  optional_cols <- c("DropletUtils", "CellRanger", "HTO_Final")
  for (col in optional_cols) {
    if (!col %in% colnames(sample_seurat@meta.data)) {
      sample_seurat <- Seurat::AddMetaData(object   = sample_seurat,
                                           metadata = "ND",
                                           col.name = col)
    }
  }
  
  # ---- 🧪 Ensure QC columns exist ----
  
  required_cols <- c("Quality", "DoubletFinder", "scDblFinder")
  if (!all(required_cols %in% colnames(sample_seurat@meta.data))) {
    log_error(sample = sample, 
              step = "filter_singlets",
              msg =paste("Missing required columns in metadata: Quality, DoubletFinder, scDblFinder.",
                         "Please run calc_qc_metrics(), classify_doubletfinder() and classify_scdblfinder()",
                         sep = "\n"))
  }
  
  # Columns to retain (if available)
  keep_cols <- c("barcode", "Cell", "Sample", "QC", "nUMIs", "nGenes", 
                 "MitoRatio", "RiboRatio", "HemeRatio", "Novelty", "Quality",
                 "DropletUtils", "CellRanger", "DoubletFinder", "scDblFinder",
                 "nHTO_UMIs", "nHTOs", "HTO_Final", "X", "Y")
  
  # ---- 🔍 Subset to retain only high-quality singlets ----
  
  # QC CONFIDENCE STRATEGY (Final Classification)
  # This strategy defines the final 'QC' status based on a prioritized hierarchy 
  # designed to balance specificity (for retention) and sensitivity (for removal).
  #
  # 1. Empty Droplets: Require BOTH DropletUtils AND CellRanger to agree ('&') 
  #    to maximize **Specificity**. High confidence ensures true empty droplets 
  #    are removed while avoiding accidental exclusion of low-UMI real cells.
  #
  # 2. Low Quality: Flagged by failure of initial quality metrics (e.g., MitoRatio, nUMIs).
  #    These cells are prioritized for removal over doublet calls.
  #
  # 3. Doublets: Remove if EITHER DoubletFinder OR scDblFinder calls doublet ('|') 
  #    to maximize **Sensitivity**. This ensures all suspicious cells are excluded, 
  #    preventing doublets from contaminating downstream analyses.
  #
  # 4. Singlets: Any cell passing the Empty Droplet, Low Quality, and Doublet criteria.
  
  # Assign final classification (QC)
  metadata <- sample_seurat@meta.data %>%
    dplyr::mutate(QC = dplyr::case_when(DropletUtils == "Empty Droplet" & CellRanger == "Empty Droplet" ~ "Empty Droplet",
                                        Quality == "Low Quality" ~ "Low Quality",
                                        DoubletFinder == "Doublet" | scDblFinder == "Doublet" ~ "Doublet",
                                        TRUE ~ "Singlet")) 
  
  
  # Select available columns
  available_cols <- base::intersect(keep_cols, colnames(metadata))
  metadata <- metadata %>%
    dplyr::select(all_of(available_cols)) 
  
  # Assign the cleaned metadata back to the Seurat object
  sample_seurat@meta.data <- metadata
  
  # Subset high quality singlets
  sample_seurat <- base::subset(x = sample_seurat, QC == "Singlet")
  
  # ---- 🪵 Log Output and Return Seurat Object, metadata ----
  
  log_info(sample = sample,
           step = "filter_singlets",
           msg =glue::glue("Retained high-quality singlets for sample : '{sample}'."))
  return(list(sample_seurat = sample_seurat,
              metadata  = metadata))
}

merge_filtered <- function(seurat_list, assay, meta_file, output_dir){
  
  # Set seed for reproducible stochastic processes
  set.seed(100)
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_list = seurat_list, assay = assay, 
                  xlsx_file = meta_file, output_dir = output_dir)
  
  # ---- 🤝 Merge Seurat Objects ----
  
  # Create sample ids from sample names
  samples <- names(seurat_list)
  sample_ids <- gsub(pattern = ".Spatial.*", replacement = "", x = samples)
  
  # NOTE: Since multiple samples can have same barcodes, add sample name to 
  # barcodes to keep track of cell identities (i.e., barcodes) coming from each 
  # sample after merging
  merged_seurat <- tryCatch({
    merge(x = seurat_list[[1]],       # get(paste0(samples[1])
          y = seurat_list[-1],        # lapply(paste0(samples[2:length(samples)]), get)
          add.cell.ids = sample_ids,  # add unique prefix to barcodes
          merge.data = FALSE)
  }, error = function(e){
    log_error(sample = "",
              step = "merge_filtered",
              msg =glue::glue("Merge failed : '{e$message}'."))
  }) 
  
  # ---- 🗑️ Remove HTO assay (if present) ----
  
  if ("HTO" %in% Seurat::Assays(merged_seurat)){
    merged_seurat[["HTO"]] <- NULL
    log_info(sample = "",
             step = "merge_filtered",
             msg ="Removed HTO assay prior to integration.")
  }
  
  # ---- 📥 Load Extra Metadata ----
  
  extra_metadata <- tryCatch({
    openxlsx::read.xlsx(xlsxFile = meta_file) %>%
      dplyr::select(-dplyr::any_of("Comments"))
  }, error = function(e){
    log_warn(sample = "",
             step = "merge_filtered",
             msg =glue::glue("Failed to read metadata file : '{e$message}'."))
    return(data.frame())  # Return empty dataframe to skip merge
  })
  
  # ---- 🔗 Join Extra Metadata ----
  
  if (nrow(extra_metadata) > 0){
    
    metadata <- merged_seurat@meta.data %>%
      tibble::rownames_to_column(var = "temp_barcode_id")
    
    if (assay == "RNA"){
      # NOTE: For droplet-based RNA, join using Sample + HTO if available
      metadata <- metadata %>%
        dplyr::mutate(Unique_ID = dplyr::case_when(HTO_Final != "ND" ~ paste0(Sample, "_", HTO_Final),
                                                   TRUE ~ Sample)) %>%
        dplyr::left_join(extra_metadata, by = c("Unique_ID" = "Unique_ID"))
    } else {
      # NOTE: For Spatial assays, join by Slide/Sample
      metadata <- metadata %>%
        dplyr::left_join(extra_metadata, by = c("Sample" = "Slide"))
      
      log_info(sample = "",
               step = "merge_filtered",
               msg ="Extra metadata joined for Spatial analysis. No cell filtering performed.")
    }
    
    # ---- 🧹 Clean Metadata ----
    
    # Remove columns that are entirely NA
    metadata <- metadata %>%
      dplyr::select(where(~ !all(is.na(.))),         # keep columns that are not entirely NA
                    -dplyr::any_of(c("Initial filename", "Comments")))  # remove unwanted columns
    
    # Re-assign metadata with original barcodes as rownames
    merged_seurat@meta.data <- metadata %>%
      tibble::column_to_rownames(var = "temp_barcode_id")
    
    log_info(sample = "",
             step = "merge_filtered",
             msg ="Extra metadata added successfully.")
  } else {
    log_info(sample = "",
             step = "merge_filtered",
             msg ="No external metadata was loaded or joined.")
  }
  
  # ---- 💾 Save Merged Seurat Object ----
  
  # Determine file name based on assay
  filename <- if (assay == "RNA"){
    "filtered_seurat.rds"
  } else{   # For Spatial.008um and Spatial.016um assays
    paste0("filtered_seurat.", assay, ".rds")
  }
  
  base::saveRDS(object = merged_seurat, file = file.path(output_dir, filename))
  
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  log_info(sample = "",
           step = "merge_filtered",
           msg =glue::glue("Filtered Seurat object saved to : '{filename}'."))
  return(invisible(merged_seurat))
}

plot_qc <- function(metadata, output_dir){
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(metadata = metadata, output_dir = output_dir)
  
  # Ensure required columns exist in metadata
  required_cols <- c("Sample", "QC", "nUMIs", "nGenes", "MitoRatio", "RiboRatio", "Novelty")
  missing_cols <- setdiff(required_cols, colnames(metadata))
  if (length(missing_cols) > 0) {
    log_error(sample = "",
              step = "plot_qc",
              msg =glue::glue("Missing required metadata columns : '{paste(missing_cols, collapse = ", ")}'."))
  }
  
  # ---- 🛠️ Helper Functions ----
  
  # Define classification levels and color palette
  qc_levels <- c("Empty Droplet", "Doublet", "Singlet", "Low Quality")
  fill_colors <- c("Empty Droplet" = "#FFC61E", "Singlet" = "#2CA02C",
                   "Doublet" = "#1F77B4", "Low Quality" = "#D62728")
  
  # Visualize cell counts per Sample
  cell_qc <- function(meta){
    
    # Count cells by Sample and QC category and fill missing levels
    df <- meta %>%
      dplyr::count(Sample, QC) %>%
      data.frame() %>%
      dplyr::mutate(QC = factor(QC, levels = qc_levels)) %>% 
      tidyr::complete(Sample, QC, fill = list(n = 1))
    
    # Numeric positions for vertical lines between samples
    sample_levels <- levels(factor(df$Sample))
    vline_positions <- seq(1.5, length(sample_levels) - 0.5, by = 1)
    
    # NOTE: position = "dodge" for grouped bars; "stack" for stacked bar
    # stat = "identity" because y is precomputed; "count" if y axis determined based on X axis frequency
    p <- ggplot(data = df, aes(x = Sample, y = n, fill = QC)) +
      geom_bar(stat = "identity", position = position_dodge(0.9)) +
      theme_classic() +
      custom_theme +
      labs(x = "Sample", y = "Cell Counts", title = "Number of Cells") +
      coord_cartesian(ylim = c(1,1e7), clip = "off", expand = FALSE) +
      scale_y_log10(breaks = c(10, 100, 1000, 10000, 100000, 1000000)) +
      scale_fill_manual(values = fill_colors) +
      #geom_text(stat ="count", aes(label = after_stat(count)), y = 0, hjust = 0, angle = 90)
      geom_text(mapping = aes(label = n, ymin = 0.1, ymax = 1),
                position = position_dodge(width = 0.9),
                y = 0.1, hjust = 0, angle = 90) +
      geom_vline(xintercept = vline_positions, linetype = "dotted", color = "gray50")
    
    return(p)
  }
  
  # Visualize nUMIs, nGenes, MitoRatio, RiboRatio, Novelty per sample
  violin_qc <- function(meta, yvar, ylab, title, cutoff = NULL, ylog = TRUE, ylim = NULL){
    
    # Fill missing levels
    df <- meta %>% 
      dplyr::mutate(QC = factor(QC, levels = qc_levels))
    
    # Numeric positions for vertical lines between samples
    sample_levels <- levels(factor(df$Sample))
    vline_positions <- seq(1.5, length(sample_levels) - 0.5, by = 1)
    
    p <- ggplot(data = df, aes(x = Sample, y = .data[[yvar]], fill = QC)) +
      geom_violin(position = position_dodge(0.9), scale = "width", drop = FALSE) +
      geom_boxplot(position = position_dodge(0.9), width = 0.05, outlier.size = 0.5) +
      theme_classic() + 
      custom_theme +
      labs(x = "Sample", y = ylab, title = title) +
      scale_fill_manual(values = fill_colors) +
      geom_vline(xintercept = vline_positions, linetype = "dotted", color = "gray50")
    
    if (!is.null(cutoff)) p <- p + geom_hline(yintercept = cutoff, linetype = 2)
    if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim, clip = "off")
    if (ylog) p <- p + scale_y_log10()
    
    return(p)
  }
  
  # Visualize number of genes/cell, number of UMIs/cell & MitoRatio together.
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
      geom_point(alpha = 0.3, size = 0.25) +
      theme_classic() +
      custom_theme +
      labs(x = "Number of UMIs", y = "Number of Genes", title = "UMIs vs Genes (Colored by MitoRatio)") +
      coord_cartesian(xlim = c(1,1e6), ylim = c(1,20000), clip = "off") +
      scale_x_log10(breaks = c(1,10,100,1000,10000,100000,1e6)) + 
      scale_y_log10(breaks = c(1,10,100,1000,10000,100000)) + 
      scale_color_viridis(option = "D", limits = c(0,1)) +
      #facet_wrap(.~Sample, nrow = 4) +   #split the plot by X-axis label
      facet_wrap(Sample ~ QC, ncol = 4, drop = FALSE) +
      geom_vline(xintercept = umi_cutoff, linetype = "dashed") +
      geom_hline(yintercept = gene_cutoff, linetype = "dashed") +
      stat_smooth(method = lm, color = "red", se = FALSE, linewidth = 0.3)
  }
  
  # ---- 💾 Generate and Save Plots ----
  
  # Plot generators
  plot_list <- list(
    Cell_Counts                      = cell_qc,
    UMI_Distribution                 = function(x) violin_qc(x, "nUMIs", "Number of UMIs", "UMI Distribution", 500, TRUE, c(1, 1e6)),
    Gene_Distribution                = function(x) violin_qc(x, "nGenes", "Number of Genes", "Gene Distribution", 250, TRUE, c(1, 30000)),
    MitoRatio_Distribution           = function(x) violin_qc(x, "MitoRatio", "MitoRatio", "MitoRatio Distribution", 0.2, FALSE, c(0, 1)), 
    RiboRatio_Distribution           = function(x) violin_qc(x, "RiboRatio", "RiboRatio", "RiboRatio Distribution", 0.05, FALSE, c(0, 1)),
    Novelty_Score_Distribution       = function(x) violin_qc(x, "Novelty", "Novelty", "Novelty Score Distribution", 0.8, FALSE, c(0, 1)),
    Genes_UMI_MitoRatio_Distribution = gene_umi_mito_qc)
  
  # Generate plots
  for (plot_name in names(plot_list)) {
    p <- plot_list[[plot_name]](metadata)    #  p <- get(funcs[i])(metadata)
    
    # Find number of samples
    n_samples <- length(unique(metadata$Sample))
    
    # Set plot width based on number of samples, e.g., 0.8 inch per sample
    max_pdf_width <- 30
    max_pdf_height <- 30
    
    # Default: 0.8 inch per sample
    desired_width_per_sample <- 0.8
    pdf_width <- min(n_samples * desired_width_per_sample + 2, max_pdf_width)  # 2 inch for margin/legend
    pdf_height <- 8
    
    # Special case for Genes_UMI_MitoRatio_Distribution (4 panels per sample)
    if (plot_name == "Genes_UMI_MitoRatio_Distribution"){
      desired_width_per_panel <- 2
      desired_height_per_sample <- 2
      pdf_width <- min(4 * desired_width_per_panel + 2, max_pdf_width)
      pdf_height <- min(n_samples * desired_height_per_sample + 2, max_pdf_height)  # 2 inch for axis etc
    }
    
    # Optional: enforce minimum dimensions
    pdf_width  <- max(pdf_width, 6)
    pdf_height <- max(pdf_height, 5)
    
    # Save plot
    ggplot2::ggsave(filename = paste0("QC_", plot_name, ".pdf"),
                    plot = p,
                    device = "pdf",
                    path = output_dir,
                    width =  pdf_width,
                    height = pdf_height,
                    dpi = 600,
                    units = "in")
  }
  
  # ---- 🪵 Log Output ----
  
  log_info(sample = "",
           step = "plot_qc",
           msg =glue::glue("QC plots generated and saved to : '{output_dir}'."))
}

run_sctransform <- function(filtered_seurat, assay, s_genes, g2m_genes){
  
  # Set seed for reproducible stochastic processes (PCA, UMAP)
  set.seed(1234)
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = filtered_seurat, assay = assay,
                  s_genes = s_genes, g2m_genes = g2m_genes)
  
  # ---- 🔍 Check for Cell Cycle Genes ----
  
  # Step 1: Remove duplicates ignoring case (currently mix of human and mouse)
  s_genes_unique <- s_genes[!duplicated(toupper(s_genes))]
  g2m_genes_unique <- g2m_genes[!duplicated(toupper(g2m_genes))]
  
  # Step 2: Keep only genes present in the Seurat object
  
  # SeuratObject::Features() returns union of all feature names stored at the 
  # assay level in the assay metadata, not in a specific layer or slot.
  # SeuratObject::GetAssayData() fails if there are multiple layers
  present_features <- SeuratObject::Features(x = filtered_seurat, layer = "counts", simplify = TRUE)
  s_genes <- intersect(s_genes, present_features)
  g2m_genes <- intersect(g2m_genes, present_features)
  
  if (length(s_genes) == 0) {
    pad <- base::strrep(" ", 29)
    log_error(sample = "",
              step = "run_sctransform",
              msg =glue::glue("Zero S-phase genes were found in '{assay}' assay\n",
                              "{pad}Cannot perform cell cycle scoring."))
    
  } else {
    log_info(sample = "",
             step = "run_sctransform",
             msg =glue::glue("{length(s_genes)} of {length(s_genes_unique)} S-phase genes found and will be used for scoring."))
  }
  
  if (length(g2m_genes) == 0) {
    pad <- base::strrep(" ", 29)
    log_error(sample = "",
              step = "run_sctransform",
              msg =glue::glue("Zero G2M-phase genes were found in '{assay}' assay\n",
                              "{pad}Cannot perform cell cycle scoring."))
  } else {
    log_info(sample = "",
             step = "run_sctransform",
             msg =glue::glue("{length(g2m_genes)} of {length(g2m_genes_unique)} G2M-phase genes found and will be used for scoring."))
  }
  
  # ---- SCTransfrom Workflow ----
  
  # 1️⃣ Normalize Data using LogNormalize (for cell cycle scoring)
  filtered_seurat <- Seurat::NormalizeData(object = filtered_seurat,
                                           assay = assay,
                                           normalization.method = "LogNormalize",
                                           scale.factor = 10000,
                                           margin = 1,
                                           verbose = FALSE)  
  
  # 2️⃣ Join Layers (for complete data view before cell cycle scoring)
  # NOTE: CellCycleScoring uses a single data layer (log-normalized counts), 
  # but currently, data layers for each sample may be stored separately. Join them first.
  filtered_seurat@assays[[assay]] <- SeuratObject::JoinLayers(object = filtered_seurat@assays[[assay]])
  
  # 3️⃣ Score Cell Cycle
  filtered_seurat <- Seurat::CellCycleScoring(object = filtered_seurat,
                                              s.features = s_genes,
                                              g2m.features = g2m_genes,
                                              ctrl = NULL)
  
  # 4️⃣ Calculate CC.Score to regress out cell cycle differences
  # https://satijalab.org/seurat/archive/v3.1/cell_cycle_vignette
  filtered_seurat$CC.Score <- filtered_seurat$G2M.Score - filtered_seurat$S.Score
  
  
  # 5️⃣ Split Object by 'Sample' (for batch-aware SCTransform)
  # NOTE: All cells within the same batch MUST be analyzed together
  filtered_seurat@assays[[assay]] <- base::split(x = filtered_seurat@assays[[assay]],
                                                 f = filtered_seurat@meta.data[["Sample"]])
  
  # 6️⃣ Run SCTransform (with CC.Score and MitoRatio regression)
  # https://github.com/satijalab/seurat/issues/7342
  sct_seurat <- quiet_msg(Seurat::SCTransform(object = filtered_seurat,
                                              assay = assay,
                                              new.assay.name = "SCT",
                                              do.correct.umi = TRUE,
                                              ncells = 5000,
                                              variable.features.n = 3000,
                                              vars.to.regress = c("CC.Score", "MitoRatio"),
                                              do.scale = FALSE,
                                              do.center = TRUE,
                                              vst.flavor = "v2",
                                              return.only.var.genes = TRUE,
                                              verbose = FALSE))
  
  # 7️⃣ Prepare SCT Assay for Differential Expression
  sct_seurat <- quiet_msg(Seurat::PrepSCTFindMarkers(object = sct_seurat,
                                                     assay = "SCT",
                                                     verbose = FALSE))
  
  # 8️⃣ Filter Out Unwanted Variable Features (Ribo, Mito, etc.)
  # NOTE: PCA, UMAP, and clustering should not be influenced by these genes.
  var_f <- Seurat::VariableFeatures(sct_seurat, assay = "SCT")
  var_f <- var_f[!grepl(pattern = "^[Rr][Pp][SsLl]|R[Ii][Kk]$|^[Mm][Tt]-|^G[Mm][0-9.]+$", x = var_f)]
  Seurat::VariableFeatures(sct_seurat, assay = "SCT") <- var_f
  #sct_seurat@assays[["SCT"]]@var.features <- var_f
  cat("\nFinal number of variable features:", length(var_f), "\n")
  
  # 9️⃣ Scale and Run PCA on Original Assay (for export compatibility)
  # NOTE: Populates 'scale.data' slot and generates PCA reduction on log-normalized data.
  # This is REQUIRED for seamless export to AnnData/Scanpy/scVI environments.
  sct_seurat <- Seurat::ScaleData(object = sct_seurat,
                                  assay = assay,
                                  features = Seurat::VariableFeatures(sct_seurat, assay = "SCT"), # Use SCT features
                                  verbose = FALSE)
  
  # Use a consistent naming scheme for dimensional reduction keys. See NOTE on 🔟.
  sct_seurat <- Seurat::RunPCA(object = sct_seurat,
                               assay = assay,
                               features = Seurat::VariableFeatures(sct_seurat, assay = "SCT"),
                               reduction.name = paste0(tolower(assay), "_pca"),
                               reduction.key = paste0(toupper(assay), "pca_"),
                               verbose = FALSE)
  
  # 🔟 Run PCA on SCT Assay
  # NOTE: IntegrateLayers() sets reduction.key as "integcca_" when reduction.name = "integ_cca".
  # For a consistent naming, we set reduction.key as "sctpca_" when reduction.name = "sct_pca".
  sct_seurat <- Seurat::RunPCA(object = sct_seurat,
                               assay = "SCT",
                               features = Seurat::VariableFeatures(sct_seurat, assay = "SCT"),
                               reduction.name = "sct_pca",
                               reduction.key = "sctpca_", 
                               verbose = FALSE)
  
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  log_info(sample = "",
           step = "run_sctransform",
           msg ="SCTransform workflow completed successfully.")
  return(invisible(sct_seurat))
}

integrate_sct_data <- function(sct_seurat, assay, reference_samples = NULL){
  
  # Set seed for reproducible stochastic processes (Harmony, PCA)
  set.seed(1234)
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = sct_seurat, assay = assay,
                  reference_samples = reference_samples)
  
  # ---- 🔍 Ensure 'sct_pca' reduction exists ----
  if (!"sct_pca" %in% names(sct_seurat@reductions)) {
    
    log_error(sample = "",
              step = "integrate_sct_data",
              msg ="Reduction 'sct_pca' is missing. Please run PCA on SCT assay before integration.")
  }
  
  # ---- ⚖️ Compute optimal k.weight for integration ----
  
  # NOTE: k.weight is typically half the number of cells in the smallest sample, capped at 100
  kweight <- min(sct_seurat@meta.data %>% 
                   dplyr::count(Sample) %>% 
                   dplyr::pull(n) %>% 
                   min()/2, 100) 
  log_info(sample = "",
           step = "integrate_sct_data",
           msg =glue::glue("Dataset k.weight : '{kweight}'."))
  
  # ---- 📏 Determine maximum number of dimensions for integration ----
  
  # NOTE: Seurat::JointPCAIntegration() gives errors if fewer than 30 dims are provided
  max_dims <- min(30, ncol(Seurat::Embeddings(sct_seurat, reduction = "sct_pca")))
  #max_dims <- min(30, ncol(sct_seurat@reductions$sct_pca@cell.embeddings))
  
  # ---- 🔁 Loop over Integration Methods ----
  
  # Methods supported:
  #   - CCA      : Canonical Correlation Analysis, aligns datasets based on correlated components
  #   - RPCA     : Reciprocal PCA, preserves sample-specific variance while integrating
  #   - Harmony  : Batch correction in PCA space, fast and scalable
  #   - JointPCA : Joint PCA integration, combines PCA embeddings from multiple datasets
  
  # Initialize Seurat object for integration
  integrated_seurat <- sct_seurat
  Seurat::DefaultAssay(integrated_seurat) <- "SCT"
  
  integration_methods <- c("CCA", "RPCA", "Harmony", "JointPCA")
  
  for (method in integration_methods){
    
    # Set the name of the reduction after integration
    reduction_name <- paste0("integ_", base::tolower(method))
    
    # Perform integration using IntegrateLayers
    integrated_seurat <- Seurat::IntegrateLayers(object = integrated_seurat,
                                                 method = paste0(method, "Integration"),
                                                 normalization.method = "SCT",
                                                 orig.reduction = "sct_pca", 
                                                 new.reduction = reduction_name,
                                                 reference = reference_samples,
                                                 k.weight = kweight,    # relevant for RPCA
                                                 dims = 1:max_dims,
                                                 verbose = FALSE)
  }
  
  # ---- Optional integration for scVI and FastMNN ----
  
  # NOTE: scVI needs raw counts. Vignette also uses it on RNA assay
  # NOTE: We use variable features of SCT assay for integration.
  # NOTE: We use pca reduction from RNA assay (derived using variable features of SCT assay)
  # FastMNN throws error "Error in checkBatchConsistency(batches, cells.in.columns = TRUE)"
  
  # for (method in c("scVI", "FastMNN")) {
  #   DefaultAssay(integrated_seurat) <- assay
  #   integrated_seurat <- Seurat::IntegrateLayers(object = integrated_seurat,
  #                                                method = paste0(r, "Integration"),
  #                                                normalization.method = "LogNormalize",
  #                                                orig.reduction = paste0(assay, ".pca"),
  #                                                features = integrated_seurat@assays$SCT@var.features,
  #                                                new.reduction = base::tolower(method),
  #                                                reference = reference_samples,
  #                                                k.weight = kweight,                                    # for RPCA
  #                                                conda_env = "/hpc/home/kailasamms/miniconda3/envs/R",  # for scVI
  #                                                verbose = FALSE)
  # }
  
  # ---- 🔗 Merge layers after integration (for RNA/Spatial assay) ----
  
  # NOTE: Only applicable for RNA or Spatial assays, not for SCT assay itself
  integrated_seurat@assays[[assay]] <- SeuratObject::JoinLayers(integrated_seurat@assays[[assay]])
  
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  log_info(sample = "",
           step = "integrate_sct_data",
           msg ="Integration completed successfully.")
  return(invisible(integrated_seurat))
}

cluster_integrated_data <- function(integrated_seurat, assay){
  
  # Set seed for reproducible stochastic processes (Leiden clustering, UMAP)
  set.seed(1234)
  
  # ---- ⚙️ Validate Input Parameters ----
  validate_inputs(seurat_object = integrated_seurat, assay = assay)
  
  # ---- 🔬 Find Nearest Neighbors (for every cell) ----
  
  # Integration methods used during integration
  integration_methods <- c("CCA", "RPCA", "Harmony", "JointPCA")
  
  # Define clustering resolutions to explore
  resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
  
  # Loop over integration methods to calculate the nearest neighbors
  for (method in integration_methods) {
    
    # Define name of reduction, NN and SNN graphs after integration
    reduction_name <- paste0("integ_", base::tolower(method))
    nn_name <- paste0(base::tolower(method), ".nn")   # kNN graph
    snn_name <- paste0(base::tolower(method), ".snn") # Shared nearest neighbor (SNN) graph
    
    # Determine max dimensions available for this reduction
    max_dims <- min(40, ncol(SeuratObject::Embeddings(integrated_seurat, reduction = reduction_name)))
    #max_dims <- min(40, ncol(integrated_seurat@reductions[[reduction_name]]@cell.embeddings))
    
    # Ensure SCT assay is active
    DefaultAssay(integrated_seurat) <- "SCT"
    
    # Compute nearest neighbors
    integrated_seurat <- Seurat::FindNeighbors(object = integrated_seurat,
                                               reduction = reduction_name,
                                               dims = 1:max_dims,
                                               k.param = 30,   # Default k for kNN
                                               graph.name = c(nn_name, snn_name),
                                               verbose = FALSE)
  }
  
  # ---- 🧩 Find Clusters using the SNN graph ----
  
  # Loop over integration methods and clustering resolutions
  for (method in integration_methods) {
    for (res in resolutions) {
      
      snn_name <- paste0(base::tolower(method), ".snn")    # SNN graph used for clustering
      cluster_name <- paste0(base::tolower(method), res)   # Cluster label name
      
      # Leiden clustering
      integrated_seurat <- Seurat::FindClusters(object = integrated_seurat,
                                                res = res,
                                                graph.name = snn_name,
                                                cluster.name = cluster_name,
                                                modularity.fxn = 1,
                                                algorithm = 4,     # 4 = Leiden algorithm (recommended)
                                                verbose = FALSE)
    }
  }
  
  # ---- 🗺️ Run UMAP for dimensionality reduction (visualization) ----
  
  for (method in integration_methods) {
    
    # Define the name of the reduction after integration
    reduction_name <- paste0("integ_", base::tolower(method))
    
    # Determine max dimensions available for this reduction
    max_dims <- min(40, ncol(SeuratObject::Embeddings(integrated_seurat, reduction = reduction_name)))
    #max_dims <- min(40, ncol(integrated_seurat@reductions[[reduction_name]]@cell.embeddings))
    
    integrated_seurat <- Seurat::RunUMAP(object = integrated_seurat,
                                         dims = 1:max_dims,
                                         n.neighbors = 30L,   # Number of neighbors for UMAP graph
                                         reduction = reduction_name,
                                         reduction.name = paste0("umap_", base::tolower(method)),
                                         verbose = FALSE)
  }
  
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  log_info(sample = "",
           step = "cluster_integrated_data",
           msg ="Clustering completed successfully.")
  return(invisible(integrated_seurat))
}

remove_sparse_clusters <- function(integrated_seurat, assay){
  
  # ---- ⚙️ Validate Input Parameters ----
  validate_inputs(seurat_object = integrated_seurat, assay = assay)
  
  # ---- 🔎 Identify all cluster columns dynamically ----
  
  # Define all possible starting prefixes based on your integration methods
  prefixes <- c("cca", "rpca", "harmony", "jointpca")
  
  # Define pattern of prefix to search in column names
  pattern <- paste0("^(", paste(prefixes, collapse = "|"), ")")
  
  # Select all columns that start with the prefix and end with a numeric cluster ID
  cluster_cols <- colnames(integrated_seurat@meta.data)
  cluster_cols <- cluster_cols[grepl(paste0("^", pattern, "[0-9.]+$"), cluster_cols)]
  
  if (length(cluster_cols) == 0) {
    log_error(sample = "",
              step = "remove_sparse_clusters",
              msg =glue::glue("No clustering columns found with prefix : '{pattern}'."))
  }
  
  # ---- 🗑️ Find and Remove Sparse Cells ----
  
  sparse_cells <- character(0)
  
  for (col in cluster_cols) {
    
    # Count number of cells per cluster in this column 
    # NOTE: Use table() for fast frequency counting
    cluster_counts <- table(integrated_seurat@meta.data[[col]])
    
    # Identify clusters with very few cells (<= 5)
    sparse_clusters_ids <- names(cluster_counts[cluster_counts <= 5])
    
    # If sparse clusters exist, identify the cells belonging to them
    if (length(sparse_clusters_ids) > 0) {
      
      # Get barcodes of cells in sparse clusters
      cells_in_sparse_clusters <- rownames(integrated_seurat@meta.data)[integrated_seurat@meta.data[[col]] %in% sparse_clusters_ids]
      
      # Append these cells to the master sparse cell list
      sparse_cells <- c(sparse_cells, cells_in_sparse_clusters)
      
      # Optional progress message
      log_info(sample = "",
               step = "remove_sparse_clusters",
               msg =glue::glue("Found '{length(cells_in_sparse_clusters)}' cells to remove from column : '{col}'."))
    }
  }
  
  # Identify barcodes in non-sparse clusters
  unique_sparse_cells <- unique(sparse_cells)
  all_cells <- SeuratObject::Cells(x = integrated_seurat)
  keep_cells <- base::setdiff(all_cells, unique_sparse_cells)
  
  # Subset Seurat object safely using barcode name (i.e. rownames)
  if (length(unique_sparse_cells) > 0) {
    integrated_seurat <- subset(x = integrated_seurat,
                                cells = keep_cells)
  }
  
  log_info(sample = "",
           step = "remove_sparse_clusters",
           msg =glue::glue("Total cells removed from sparse clusters : '{length(unique_sparse_cells)}'."))
  
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  log_info(sample = "",
           step = "remove_sparse_clusters",
           msg ="Successfully removed sparse cells.")
  return(invisible(integrated_seurat))
}

calc_optimal_resolution <- function(integrated_seurat, reduction, output_dir){
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = integrated_seurat, reduction = reduction,
                  output_dir = output_dir)
  
  # ---- 🔎 Identify all cluster columns dynamically ----
  
  # Define all possible starting prefixes based on your integration method
  prefixes <- base::tolower(reduction)
  
  # Remove common reduction prefixes
  pattern <- stringr::str_remove(prefixes, "^(umap_|integrated_|sct_|pca_)*")
  
  # Select all columns that start with the prefix and end with a numeric cluster ID
  cluster_cols <- colnames(integrated_seurat@meta.data)
  cluster_cols <- cluster_cols[grepl(paste0("^", pattern, "[0-9.]+$"), cluster_cols)]
  
  if (length(cluster_cols) == 0) {
    log_error(sample = "",
              step = "calc_optimal_resolution",
              msg =glue::glue("No clustering columns found with prefix : '{pattern}'."))
  }
  
  # ---- 🌳 Visualize Clustering Stability (Clustree) ----
  
  # Sort cluster columns by numeric res (REQUIRED for SC3 stability scoring)
  cluster_cols <- cluster_cols[order(as.numeric(stringr::str_extract(cluster_cols, "[0-9.]+")))]
  
  # Run clustree
  clustree_plot <- clustree::clustree(x = integrated_seurat,
                                      prefix = pattern)
  
  # Save clustree plot
  ggplot2::ggsave(filename = paste0("Clustree_", pattern, ".pdf"),
                  plot = clustree_plot,
                  path = output_dir,
                  device = "pdf",
                  width = 8,
                  height = 11,
                  dpi = 300,
                  units = "in")
  
  # ---- 📈 Quantify Clustering Stability (SC3) ----
  
  # Convert factor/character cluster columns to integer IDs
  cluster_df <- integrated_seurat@meta.data[, cluster_cols, drop = FALSE]
  cluster_df[] <- lapply(cluster_df, function(x) as.integer(as.factor(x)))
  cluster_matrix <- as.matrix(cluster_df)
  
  # Calculate SC3 stability
  # NOTE: stability_mat uses res index like 1, 2, 3 instead of actual
  # res names like "harmony0.4" etc...
  stability_mat <- clustree:::calc_sc3_stability(cluster_matrix)
  
  # Compute mean stability per res and map res names
  stability_df <- stability_mat %>%
    as.data.frame() %>%
    dplyr::mutate(across(.cols = everything(), .fns = as.numeric)) %>%
    dplyr::group_by(resolution) %>%
    dplyr::summarise(mean_stability = mean(stability, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(resolution_name = colnames(cluster_matrix)[resolution])
  
  # Pick the res with highest mean stability
  optimal_res <- stability_df %>%
    dplyr::slice_max(mean_stability) %>%
    dplyr::pull(resolution_name)
  
  # Store the optimal res in seurat object
  optimal_res_numeric <- as.numeric(stringr::str_extract(optimal_res, "[0-9.]+"))
  integrated_seurat$optimal_res <- optimal_res_numeric
  
  # Visualize stability across resolutions
  stability_plot <- ggplot(stability_df, aes(x = resolution_name, y = mean_stability, group = 1)) +
    geom_line(color = "blue") +
    geom_point(size = 3, color = "red") +
    theme_classic() +
    labs(title = "Mean Cluster Stability per Resolution",
         x = "Resolution",
         y = "Mean SC3 Stability")
  
  # Save stability scores plot
  ggplot2::ggsave(filename = paste0("Clustree_Stability_scores_", pattern, ".pdf"),
                  plot = stability_plot,
                  path = output_dir,
                  device = "pdf",
                  width = 8,
                  height = 11,
                  dpi = 300,
                  units = "in")
  
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  log_info(sample = "",
           step = "calc_optimal_resolution",
           msg =glue::glue("Optimal res identified : '{optimal_res_numeric}'."))
  return(invisible(integrated_seurat))
}

identify_markers <- function(integrated_seurat, res, reduction, cluster_col = NULL, output_dir){
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = integrated_seurat, res = res,
                  reduction = reduction, metadata_cols = c(cluster_col),
                  output_dir = output_dir)
  
  # ---- 🧪 Detect Appropriate Assay for FindAllMarkers ----
  
  all_assays <- names(integrated_seurat@assays)
  
  if ("SCT" %in% all_assays) {
    active_assay <- "SCT"
  } else if ("RNA" %in% all_assays) {
    active_assay <- "RNA"
  } else if (length(all_assays) > 0) {
    active_assay <- all_assays[1]
    log_info(sample = "",
             step = "identify_markers",
             msg =glue::glue("No SCT/RNA assay found. Using assay : '{active_assay}'."))
  } else {
    log_error(sample = "",
              step = "identify_markers",
              msg ="No assays found in Seurat object for DE analysis.")
  }
  
  # Set active assay
  DefaultAssay(integrated_seurat) <- active_assay
  
  # ---- 🏷️ Set Idents Robustly ----
  
  # Define all possible starting prefixes based on your integration method
  prefixes <- base::tolower(reduction)
  
  # Remove common reduction prefixes
  pattern <- stringr::str_remove(prefixes, "^(umap_|integrated_|sct_|pca_)*")
  
  # Append the res number to get the exact cluster column name
  pattern <- paste0(pattern, res)
  
  # Select column(s) with exact name match
  cluster_cols <- colnames(integrated_seurat@meta.data)
  cluster_cols <- cluster_cols[cluster_cols == pattern]
  
  # Check if identifier column exists
  if (length(cluster_cols) == 0 & is.null(cluster_col)) {
    log_error(sample = "",
              step = "identify_markers",
              msg =glue::glue("Either 'res' + 'reduction' OR 'cluster_col' must be provided!"))
  } else if (length(cluster_cols) > 1) {
    log_error(sample = "",
              step = "identify_markers",
              msg =glue::glue("More than 1 matching column found in metadata : '{paste(cluster_cols, collapse = ", ")}'"))
  }
  
  # Set active ident
  if (!is.null(cluster_col)){
    cluster_cols <- cluster_col
    Idents(object = integrated_seurat) <- cluster_cols
  }
  
  # ---- 🔍 Find ALL Markers ----
  
  all_markers <- quiet_msg(Seurat::FindAllMarkers(object = integrated_seurat,
                                                  assay = active_assay,
                                                  logfc.threshold = 0.25, 
                                                  test.use = "wilcox",
                                                  slot = "data",
                                                  min.pct = 0.1,
                                                  min.diff.pct = -Inf,
                                                  only.pos = TRUE))
  
  if (nrow(all_markers) == 0) {
    
    log_warn(sample = "",
             step = "identify_markers",
             msg =glue::glue("No markers found for res '{res}' using method '{reduction}'."))
    return(invisible(FALSE))
  }
  
  # ---- 🧹 Annotation and Filtering ----
  
  # Get annotations from ENSEMBL
  annotations_list <- get_annotations()
  
  # Compute number of matching genes for each annotation list
  matches <- sapply(annotations_list, function(x) length(intersect(x$SYMBOL, all_markers$gene)))
  
  # Pick the annotation list with the most matches
  annotations <- annotations_list[[which.max(matches)]]
  
  sig_markers <- all_markers %>%
    # Correct division by zero warning before calculating ratio
    dplyr::mutate(pct.1 = dplyr::if_else(pct.1 == 0, 0.001, pct.1),
                  pct.2 = dplyr::if_else(pct.2 == 0, 0.001, pct.2),
                  ratio = pct.1 / pct.2) %>%
    dplyr::filter(p_val_adj <= 0.05) %>%
    dplyr::left_join(y = annotations, by = c("gene" = "ENSEMBL_SYMBOL")) %>%
    dplyr::relocate(cluster, gene, CHR, avg_log2FC, p_val, p_val_adj, pct.1, pct.2, ratio, DESCRIPTION) %>%
    dplyr::distinct(cluster, gene, avg_log2FC, pct.1, pct.2, .keep_all = TRUE)
  
  if (nrow(sig_markers) == 0) {
    log_warn(sample = "",
             step = "identify_markers",
             msg =glue::glue("No markers (p_val_adj < 0.05) identified at res '{cluster_cols}'. Skipping saving."))
    return(invisible(FALSE))
  }
  
  # Find top 30 markers for each major cluster
  top_markers <- sig_markers %>%
    dplyr::filter(avg_log2FC >= 0.58, p_val_adj <= 0.05) %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(desc(avg_log2FC)) %>%
    dplyr::slice_head(n = 30) %>%
    ungroup()
  
  # ---- 📊 Create Ordered Marker Matrix for Heatmap ----
  
  mat <- sig_markers %>%
    tidyr::pivot_wider(id_cols = cluster,
                       names_from = gene,
                       values_from = avg_log2FC,
                       values_fill = 0) %>%
    tibble::column_to_rownames("cluster") %>%
    scale()   # column wise scaling so each gene has mean 0, stdev = 1
  
  # Cluster rows and columns
  row_dist <- stats::dist(x = mat, method = "euclidean")         # distance between rows based on columns
  row_clust <- stats::hclust(d = row_dist, method = "ward.D2")   # clustering based on distance calculated
  col_dist <- stats::dist(x = t(mat), method = "euclidean")
  col_clust <- stats::hclust(d = col_dist, method = "ward.D2")
  
  # Reorder matrix
  row_order <- rownames(mat[row_clust$order,])
  col_order <- colnames(mat[,col_clust$order])
  mat <- mat[row_order, col_order]
  
  # Truncate negative values and round
  mat[mat < 0] <- 0
  mat <- mat %>%
    t() %>%               # transpose so genes are rows
    as.data.frame() %>%   
    dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 2)))
  
  # ---- 💾 Save Markers to Excel ----
  
  file_name <- file.path(output_dir, paste0("Markers.All.", cluster_cols, ".xlsx"))
  wb <- openxlsx::createWorkbook()
  
  openxlsx::addWorksheet(wb, "All_Markers")
  openxlsx::writeData(wb, "All_Markers", sig_markers)
  
  openxlsx::addWorksheet(wb, "Top_Markers")
  openxlsx::writeData(wb, "Top_Markers", top_markers)
  
  openxlsx::addWorksheet(wb, "Matrix")
  openxlsx::writeData(wb, "Matrix", mat, rowNames = TRUE)
  
  openxlsx::saveWorkbook(wb, file = file_name, overwrite = TRUE)
  
  # ---- 🪵 Log Output ----
  
  log_info(sample = "",
           step = "identify_markers",
           msg =glue::glue("Marker analysis complete. Results saved to : '{file_name}'."))
}

plot_seurat <- function(integrated_seurat, reduction, features, filename, output_dir, raster = FALSE, split_col = NULL){
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = integrated_seurat, reduction = reduction, 
                  features = features, metadata_cols = c(split_col), 
                  filename = filename, output_dir = output_dir)
  
  
  # ---- 🔄 Check and Update Reduction Name ----
  
  if (!(reduction %in% names(integrated_seurat@reductions))) {
    log_warn(sample = "",
             step = "plot_umap",
             msg =glue::glue("Reduction '{reduction}' is NOT present in Seurat object."))
    
    # Check alternative reduction name
    alt_reduction <- paste0("umap_", tolower(reduction))
    
    if (alt_reduction %in% names(integrated_seurat@reductions)) {
      log_info(sample = "",
               step = "plot_umap",
               msg =glue::glue("Using alternative reduction : '{alt_reduction}'."))
      reduction <- alt_reduction
    } else {
      log_error(sample = "",
                step = "plot_umap",
                msg =glue::glue("Alternative reduction '{alt_reduction}' is NOT present."))
    }
  }
  
  # ---- 🧪 Detect Appropriate Assay for FeaturePlot ----
  
  all_assays <- names(integrated_seurat@assays)
  
  if ("SCT" %in% all_assays) {
    active_assay <- "SCT"
  } else if ("RNA" %in% all_assays) {
    active_assay <- "RNA"
  } else if (length(all_assays) > 0) {
    active_assay <- all_assays[1]
    log_info(sample = "",
             step = "plot_umap",
             msg =glue::glue("No SCT/RNA assay found. Using assay : '{active_assay}'."))
  } else {
    log_error(sample = "",
              step = "plot_umap",
              msg ="No assays found in Seurat object for DE analysis.")
  }
  
  # Set active assay
  DefaultAssay(integrated_seurat) <- active_assay
  
  # ---- 📊 Feature Property Calculation ----
  
  # Get all possible features in seurat object
  # NOTE: features can be genes in assays, module_score columns or other columns in metadata
  metadata_cols <- colnames(integrated_seurat@meta.data)
  present_features <- rownames(SeuratObject::GetAssayData(integrated_seurat, assay = active_assay, layer = "data"))
  
  bad_prefixes <- c("cca", "rpca", "harmony", "jointpca")
  modules <- metadata_cols[base::grepl(pattern = "[0-9]+$|_UCell", x = metadata_cols)]
  modules <- modules[!grepl(paste0("(", paste(bad_prefixes, collapse = "|"), ")"), modules)]
  
  # Initialize a list to store feature properties
  features_list   <- list()
  values_list    <- list()
  
  for (feature in features) {
    
    # Determine Feature Type and fetch raw values
    if (feature %in% present_features) {
      type <- "gene"
      is_continuous <- TRUE
      values <- SeuratObject::GetAssayData(integrated_seurat, assay = active_assay, layer = "data")[feature, ]
      
    } else if (feature %in% modules) {
      type <- "module"
      is_continuous <- TRUE
      values <- integrated_seurat@meta.data[[feature]]
      
    } else if (feature %in% c("Stability_Score")) {
      type <- "metadata"
      is_continuous <- FALSE
      values <- integrated_seurat@meta.data[[feature]]
      
    } else if (feature %in% metadata_cols) {
      type <- "metadata"
      values <- integrated_seurat@meta.data[[feature]]
      
      # Determine Continuity
      if (is.numeric(values)) {
        has_decimals <- any(abs(values %% 1) > .Machine$double.eps^0.5, na.rm = TRUE)
        has_many_unique_levels <- length(unique(values)) >= 80
        is_continuous <- has_decimals | has_many_unique_levels
      } else {
        is_continuous <- FALSE
      }
      
    } else {
      # Feature not found (Skip missing feature logic)
      next
    }
    
    # Calculate Min/Max (Only necessary for continuous features)
    min_expr <- NA
    max_expr <- NA
    if (is_continuous) {
      if (min(values, na.rm = TRUE) >= 0) {
        # Sequential (non-negative)
        min_expr <- 0
        max_expr <- quantile(values, probs = 0.99, na.rm = TRUE)
      } else {
        # Diverging
        min_expr <- quantile(values, probs = 0.01, na.rm = TRUE)
        max_expr <- quantile(values, probs = 0.99, na.rm = TRUE)
      }
    }
    
    # Save scalars → metadata
    features_list[[feature]] <- list(feature       = feature,
                                     type          = type,
                                     is_continuous = is_continuous,
                                     min_expr      = min_expr,
                                     max_expr      = max_expr)
    # Save long vector separately
    values_list[[feature]] <- values
  }
  
  # Convert to a dataframe for easier global lookup
  features_df <- do.call(base::rbind, lapply(features_list, function(x) data.frame(x, stringsAsFactors = FALSE)))
  
  # Correctly assign the feature names as row names at the end
  rownames(features_df) <- names(features_list)
  
  
  # ---- ⚖️ Global Scale Determination ----
  
  # NOTE: global min and max allows us to set same scale for all plots if 
  # plotting multiple genes or module scores 
  global_min <- NULL
  global_max <- NULL
  
  # If multiple features are to be plotted
  if (length(features) > 1) {
    
    # Check if ALL features belong to the same type
    
    # If all features are genes or module scores
    if (all(features_df$type == "gene") | all(features_df$type == "module")){
      
      # Combine all min/max values for scaling
      all_min_expr <- features_df$min_expr
      all_max_expr <- features_df$max_expr
      
      # Determine if the overall scale should be sequential or diverging
      if (min(all_min_expr, na.rm = TRUE) >= 0) {
        global_min <- 0
      } else {
        # Use the 1st percentile of all MINS to clip outlier minimums
        global_min <- quantile(all_min_expr, probs = 0.01, na.rm = TRUE)
      }
      
      # Use the 99th percentile of all MAXS to clip outlier maximums
      global_max <- quantile(all_max_expr, probs = 0.99, na.rm = TRUE)
      
    } else if (all(features_df$type == "metadata")) {
      
      # If all features are non-module metadata, they must be treated individually
      # because some might be continuous (e.g., QC metrics) and some categorical.
      # We enforce no global scale here.
      global_min <- NULL
      global_max <- NULL
      
    } else {
      # If the features are mixed (e.g., GeneA and MetadataC), or mixed continuous types,
      # the function should likely halt or enforce single-feature plotting.
      
      stop("Cannot plot a mix of feature types (gene, module score, metadata) simultaneously. Please request features of only one type (e.g., all genes, or all module scores).")
    }
  }
  
  # ---- 👥 Determine Groups for Plotting ----
  
  if (is.null(split_col)) {
    # No splitting → single panel
    group_vars <- "All"
  } else if (length(split_col) == 1) {
    # Single-column split → one panel per unique value
    
    vals <- integrated_seurat@meta.data[[split_col]]
    vals_chr <- as.character(vals)
    vals_num <- suppressWarnings(as.numeric(vals_chr))
    
    group_vars <- if (all(!is.na(vals_num))) {
      # Numeric-like → numeric sort, return character
      vals_chr[order(vals_num)] %>%
        unique()
    } else {
      # Non-numeric → character sort
      sort(unique(vals_chr))
    }
  } else {
    stop("Use Only one column for splitting")
    # # Multi-column split → one panel per column name (no subsetting)
    # group_vars <- split_col
  }
  
  # ---- 🖼️ Generate Plots for each group ----
  
  all_plots <- list()
  all_plots_labelled <- list()
  
  for (feature in features) {
    
    # Skip if feature was not processed
    if (!(feature %in% rownames(features_df))) next
    
    is_gene <- features_list[[feature]]$type == "gene"
    is_module <- features_list[[feature]]$type == "module"
    is_metadata <- features_list[[feature]]$type == "metadata"
    is_continuous <- features_list[[feature]]$is_continuous
    expr_values <- values_list[[feature]]
    
    # Get UMAP co-ordinates and add color column/expr values to UMAP co-ordinates
    full_df <- SeuratObject::Embeddings(object = integrated_seurat, reduction = reduction) %>% 
      as.data.frame() %>%
      dplyr::rename(UMAP_1 = 1, UMAP_2 = 2) %>%
      dplyr::mutate(!!feature := expr_values)
    
    for (group_var in group_vars) {
      
      if (length(split_col) == 1) {
        # Single-column split → subset by value
        subset_obj <- integrated_seurat[, integrated_seurat@meta.data[[split_col]] == group_var]
        barcodes <- rownames(subset_obj@meta.data)
        
        df <- full_df %>% 
          dplyr::filter(rownames(.) %in% barcodes)
        
      } else{
        df <- full_df
      }
      
      if (is_metadata & !is_continuous){
        
        # Determine levels of feature
        all_levels <- sort(unique(integrated_seurat@meta.data[[feature]]))
        
        # Convert grouping column to factor to ensure proper ordering for colors
        df[[feature]] <- factor(df[[feature]], levels = all_levels)
        
        # Assign colors, keeping names aligned with all_levels
        umap_palette <- custom_palette[seq_along(all_levels)]
        names(umap_palette) <- all_levels
        
        # [OPTIONAL] Calculate Cluster Centroids for Labeling
        cluster_centroids <- df %>%
          dplyr::group_by(.data[[feature]]) %>%
          dplyr::summarise(UMAP_1 = median(UMAP_1, na.rm = TRUE), 
                           UMAP_2 = median(UMAP_2, na.rm = TRUE),
                           .groups = 'drop') %>%
          dplyr::filter(!is.na(UMAP_1), !is.na(UMAP_2))
        
      } else {
        
        # Define color scale
        cols <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
        #cols = c("grey", viridis(n = 10, option = "C", direction = -1))
        #cols = c("grey", viridis(n = 10, option = "C", direction = 1))
        #cols =  c("#440154FF", viridis(n = 10, option = "C", direction = 1))
        
        # Set min and max using Seurat's default quantiles (1% and 99%)
        if (min(expr_values, na.rm = TRUE) >= 0) {
          # Sequential (non-negative)
          plot_min <- 0
          plot_max <- global_max %||% max(expr_values, na.rm = TRUE)
          
          if (is_gene | is_module){
            # Clip max only
            df[[feature]] <- pmin(df[[feature]], plot_max)
          }
          
          # Gradient params
          mid_frac <- abs(plot_min) / (abs(plot_min) + plot_max)  # position of 0 in the palette
          
          # If min_expr = max_expr = 0, then mid_frac is not finite and gives error
          if (plot_max == 0 && plot_min == 0) {
            plot_max <- 1
            mid_frac <- 0
          }
          
          values <- c(seq(0, mid_frac, length.out = 6), seq(mid_frac, 1, length.out = 6)[-1])
          values <- unique(values)
          cols <- cols[6:11]     # use only red half
          
        } else {
          # Diverging
          plot_min <- global_min %||% min(expr_values, na.rm = TRUE)
          plot_max <- global_max %||% max(expr_values, na.rm = TRUE)
          
          if (is_gene | is_module){
            # Clip min and max
            df[[feature]] <- pmin(df[[feature]], plot_max)
            df[[feature]] <- pmax(df[[feature]], plot_min)
          }
          
          # Map 11 colors to min -> 0 -> max
          mid_frac <- abs(plot_min) / (abs(plot_min) + plot_max)  # position of 0 in the palette
          values <- c(seq(0, mid_frac, length.out = 6), seq(mid_frac, 1, length.out = 6)[-1])
          cols <- cols
        } 
        
        # Order points for plotting
        df <- df[base::order(df[[feature]]), ]   # same ordering Seurat would use
      }
      
      # Plot title
      title <- ifelse(group_var == "All", feature, group_var)
      title <- as.character(title)
      
      # Define and rasterize ONLY the point layer
      point_layer <- geom_point(size = 0.2, stroke = 0)
      if (raster == TRUE) {
        point_layer <- ggrastr::rasterise(point_layer, dpi = 300)
      }
      
      # Plot
      p <- ggplot(data = df, 
                  aes(x = UMAP_1, y = UMAP_2, color = .data[[feature]])) +
        point_layer +
        theme_classic() +
        coord_fixed(ratio = 1) +
        ggplot2::labs(color = feature, x = "UMAP_1", y = "UMAP_2", title = title) +
        custom_theme
      
      if (is_metadata & !is_continuous){
        p <- p + ggplot2::scale_color_manual(values = umap_palette) +
          guides(color = guide_legend(override.aes = list(size = 3)))
      } else {
        p <- p + ggplot2::scale_colour_gradientn(colours = cols,
                                                 values = values,
                                                 name = NULL,
                                                 limits = c(plot_min, plot_max))
      }
      
      # Create another plot with cluster labels
      if (is_metadata & !is_continuous){
        q <- p +
          # Add Labels using ggrepel::geom_text_repel for non-overlapping labels
          ggrepel::geom_text_repel(data = cluster_centroids,
                                   mapping = aes(label = .data[[feature]]),
                                   color = "black",           # Set label color to black
                                   size = 12 / ggplot2::.pt,  # Divide by ggplot2::.pt = 2.845276mm to get 12 point             
                                   point.padding = NA,
                                   segment.colour = 'transparent') # No lines connecting labels to clusters
      } else{
        q <- NULL
      }
      
      all_plots[[title]] <- p
      all_plots_labelled[[title]] <- q
      
      log_info(sample = "",
               step = "plot_features",
               msg =glue::glue("Successfully plotted feature : '{feature}'."))
    }
  }
  
  # ---- 🌐 Combine Plots Using cowplot ----
  
  n_plots <- length(all_plots)
  ncol_plots <- ceiling(sqrt(n_plots))
  nrow_plots <- ceiling(n_plots / ncol_plots)
  
  # Restrict unsupported combination
  if (ncol_plots > 10 && nrow_plots > 10) {
    
    log_error(sample = "",
              step = "plot_features",
              msg ="Image size too large. More than 100 plots cannot be viewed in a single figure")
  }
  
  # Combine all plots
  combined_plot <- cowplot::plot_grid(plotlist = all_plots, 
                                      ncol = ncol_plots, 
                                      nrow = nrow_plots)
  if (length(all_plots_labelled) > 1){
    combined_plot_labelled <- cowplot::plot_grid(plotlist = all_plots_labelled, 
                                                 ncol = ncol_plots, 
                                                 nrow = nrow_plots)
  }
  
  # ---- 💾 Save Combined Plot ----
  
  file_extension <- ".pdf"
  file_name <- file.path(output_dir, paste0(filename, file_extension))
  ggplot2::ggsave(filename  = file_name,
                  plot      = combined_plot,
                  device    = cairo_pdf, 
                  width     = ncol_plots * 8,  # extra 2 inch for legend
                  height    = nrow_plots * 6, 
                  units     = "in",
                  limitsize = FALSE,
                  bg        = "white")
  
  
  if (length(all_plots_labelled) > 1){
    file_name_labelled <- file.path(output_dir, paste0(filename, "_labelled", file_extension))
    ggplot2::ggsave(filename  = file_name_labelled,
                    plot      = combined_plot_labelled,
                    device    = cairo_pdf, 
                    width     = ncol_plots * 8,  # extra 2 inch for legend
                    height    = nrow_plots * 6, 
                    units     = "in",
                    limitsize = FALSE,
                    bg        = "white")
  }
  
  # ---- 🪵 Log Output and Return ----
  
  log_info(sample = "",
           step = "plot_umap",
           msg =glue::glue("UMAP plot saved successfully to : '{file_name}'."))
  
  return(invisible(NULL))
}

plot_metrics_post_integration <- function(integrated_seurat, output_dir){
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = integrated_seurat, output_dir = output_dir)
  
  # ---- 1️⃣ Feature Plots for QC Metrics ----
  
  qc_features <- c("nUMIs", "nGenes", "S.Score", "G2M.Score", "CC.Score", "MitoRatio")
  plot_seurat(integrated_seurat, reduction = "umap_harmony",  features = qc_features,  filename = "UMAP.Numerical.Metrics", output_dir, split_col = NULL)
  
  qc_features <- c("DropletUtils", "CellRanger", "Quality", "DoubletFinder", "scDblFinder", "QC")
  plot_seurat(integrated_seurat, reduction = "umap_harmony",  features = qc_features,  filename = "UMAP.Cell.Metrics",      output_dir, split_col = NULL)
  
  # ---- 2️⃣ UMAP Plots for Cluster Visualization ----
  
  plot_seurat(integrated_seurat, reduction = "sct_pca",       features = "harmony0.8", filename = "Pre.Integ.PCA",          output_dir, split_col = "Sample")
  plot_seurat(integrated_seurat, reduction = "integ_harmony", features = "harmony0.8", filename = "Post.Integ.PCA",         output_dir, split_col = "Sample")
  plot_seurat(integrated_seurat, reduction = "umap_harmony",  features = "harmony0.8", filename = "UMAP.Sample",            output_dir, split_col = "Sample")
  plot_seurat(integrated_seurat, reduction = "umap_harmony",  features = "harmony0.8", filename = "UMAP.Phase",             output_dir, split_col = "Phase")
  
  plot_seurat(integrated_seurat, reduction = "umap_cca",      features = paste0("cca",      seq(0.2, 1.2, by = 0.2)),  filename = "UMAP.CCA",      output_dir, split_col = NULL)
  plot_seurat(integrated_seurat, reduction = "umap_rpca",     features = paste0("rpca",     seq(0.2, 1.2, by = 0.2)),  filename = "UMAP.RPCA",     output_dir, split_col = NULL)
  plot_seurat(integrated_seurat, reduction = "umap_jointpca", features = paste0("jointpca", seq(0.2, 1.2, by = 0.2)),  filename = "UMAP.JointPCA", output_dir, split_col = NULL)
  plot_seurat(integrated_seurat, reduction = "umap_harmony",  features = paste0("harmony",  seq(0.2, 1.2, by = 0.2)),  filename = "UMAP.Harmony",  output_dir, split_col = NULL)
  
  # ---- 🪵 Log Output ----
  
  log_info(sample = "",
           step = "plot_metrics_post_integration",
           msg ="All post-integration QC and UMAP plots generated successfully.")
}

calc_module_scores <- function(integrated_seurat, reduction, marker_file, filename = NULL, output_dir){
  
  # Set seed for reproducible stochastic processes (UCell)
  set.seed(1234)
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = integrated_seurat, reduction = reduction,
                  xlsx_file = marker_file, filename = filename, output_dir = output_dir)
  
  # ---- 🔄 Check and Update Reduction Name ----
  
  if (!(reduction %in% names(integrated_seurat@reductions))) {
    log_warn(sample = "",
             step = "calc_module_scores",
             msg =glue::glue("Reduction '{reduction}' is NOT present in Seurat object."))
    
    # Check alternative reduction name
    alt_reduction <- paste0("umap_", tolower(reduction))
    
    if (alt_reduction %in% names(integrated_seurat@reductions)) {
      log_info(sample = "",
               step = "calc_module_scores",
               msg =glue::glue("Using alternative reduction : '{alt_reduction}'."))
      reduction <- alt_reduction
    } else {
      log_error(sample = "",
                step = "calc_module_scores",
                msg =glue::glue("Alternative reduction '{alt_reduction}' is NOT present."))
    }
  }
  
  # ---- 🧪 Detect Appropriate Assay for FeaturePlot ----
  
  all_assays <- names(integrated_seurat@assays)
  
  if ("SCT" %in% all_assays) {
    active_assay <- "SCT"
  } else if ("RNA" %in% all_assays) {
    active_assay <- "RNA"
  } else if (length(all_assays) > 0) {
    active_assay <- all_assays[1]
    log_info(sample = "",
             step = "calc_module_scores",
             msg =glue::glue("No SCT/RNA assay found. Using assay : '{active_assay}'."))
  } else {
    log_error(sample = "",
              step = "calc_module_scores",
              msg ="No assays found in Seurat object for DE analysis.")
  }
  
  # Set active assay
  DefaultAssay(integrated_seurat) <- active_assay
  
  if (active_assay != "SCT"){
    log_warn(sample = "",
             step = "calc_module_scores",
             msg =glue::glue(" Currently using assay : '{active_assay}'."))
  }
  
  # ---- 📥 Create Feature List from Marker file ----
  
  # Load marker file
  marker_df <- tryCatch({
    openxlsx::read.xlsx(xlsxFile = marker_file)
  }, error = function(e){
    log_error(sample = "",
              step = "calc_module_scores",
              msg =glue::glue("Failed to read marker file : '{e$message}'."))
  })
  
  # Initialize an empty list to store all cell type signatures
  signatures_list <- list()
  
  # Iterate through each column (cell type) in the marker dataframe
  for (i in base::seq_len(ncol(marker_df))){
    
    # Define name of modules
    module_name <- make.names(colnames(marker_df)[i])
    
    # Determine features from the marker_df
    xlsx_features <- marker_df[[i]] %>%
      stats::na.omit() %>%
      as.vector()
    
    # Determine features present in data set
    present_features <- rownames(SeuratObject::GetAssayData(object = integrated_seurat, 
                                                            assay = active_assay, 
                                                            layer = "data"))
    
    # Match features (case-insensitive) and filter to only genes present in the data
    features <- present_features[base::tolower(present_features) %in% base::tolower(xlsx_features)]
    
    # Only add the validated and sorted feature vector to the master list if it contains 2 or more genes
    if (length(features) >= 2) {
      signatures_list[[module_name]] <- sort(features)
    } else {
      log_warn(sample = "",
               step = "calc_module_scores",
               msg =glue::glue("Skipping module '{module_name}'. Fewer than 2 matching markers found."))
    }
  }
  
  # ---- 🧬 Calculate Module Scores per Cell Type ----
  
  if (!"data" %in% SeuratObject::Layers(integrated_seurat[[active_assay]])) {
    log_error(sample = "",
              step = "calc_module_scores",
              msg =glue::glue("'data' layer missing in assay '{active_assay}'. Please normalize data first."))
  }
  
  # Calculate Module Score using Seurat
  integrated_seurat <- Seurat::AddModuleScore(object = integrated_seurat,
                                              features = signatures_list,
                                              assay = active_assay,
                                              layer = "data", # Use normalized expression data for scoring
                                              name = names(signatures_list))
  # Calculate Module Score using uCell
  integrated_seurat <- UCell::AddModuleScore_UCell(obj = integrated_seurat,
                                                   features = signatures_list,
                                                   assay = active_assay,
                                                   slot = "data",
                                                   name = "_UCell")
  log_info(sample = "",
           step = "calc_module_scores",
           msg =glue::glue("Successfully calculated module scores."))
  
  # ---- 🔎 Identify Module Score Columns ----
  
  # Find all columns ending in '.1' (default suffix from Seurat::AddModuleScore)
  score_cols <- colnames(integrated_seurat@meta.data)
  seurat_modules <- score_cols[base::grepl(pattern = "[0-9]+$", x = score_cols)]
  ucell_modules <- score_cols[base::grepl(pattern = "_UCell$", x = score_cols)]
  
  bad_prefixes <- c("cca", "rpca", "harmony", "jointpca")
  seurat_modules <- seurat_modules[!grepl(paste0("(", paste(bad_prefixes, collapse = "|"), ")"), seurat_modules)]
  ucell_modules <- ucell_modules[!grepl(paste0("(", paste(bad_prefixes, collapse = "|"), ")"), ucell_modules)]
  
  if (length(seurat_modules) < 2 & length(ucell_modules) < 2) {
    log_error(sample = "", 
              step = "annotate_cells", 
              msg ="Fewer than 2 module score columns (ending in '1' or '_UCell') found. Aborting.")
    return(integrated_seurat)
  }
  
  # ---- Z-score transformation of module scores ----
  
  # Extract the raw score matrix
  seurat_score_matrix <- integrated_seurat@meta.data[, seurat_modules, drop = FALSE] %>%
    as.matrix()
  ucell_score_matrix <- integrated_seurat@meta.data[, ucell_modules, drop = FALSE] %>%
    as.matrix()
  
  # Set negative module scores to 0 because they reflect lack of enrichment; 
  # keeping them would artificially inflate probabilities during min–max normalization.
  seurat_score_matrix[seurat_score_matrix < 0] <- 0
  ucell_score_matrix[ucell_score_matrix < 0] <- 0   #redundant as all ucell scores always > 0
  
  # Z-score standardization across all modules for each cell (row-wise)
  # The 'scale()' function performs Z-score transformation (mean=0, SD=1).
  # We transpose (t()) the matrix, apply scale(), and then transpose back
  # to ensure the scaling is applied across the genesets (rows in the original matrix).
  seurat_zscore_matrix <- t(scale(t(seurat_score_matrix)))
  ucell_zscore_matrix <- t(scale(t(ucell_score_matrix)))
  
  # ---- 💾 Save integrated Object ----
  
  # Determine file name based on assay
  filename <- if (assay == "RNA") {
    "integrated_seurat.rds"
  } else {
    paste0("integrated_seurat.", assay, ".rds")
  }
  
  base::saveRDS(integrated_seurat, file = file.path(output_dir, filename))
  
  # ---- 🖼️ Generate Plots for each module ----
  
  # Define module names 
  # NOTE: Use make.names() to ensure the score column name is valid and unique in metadata
  module_names <- make.names(colnames(marker_df))
  seurat_modules <- paste0(module_names, 1:length(module_names))        # Seurat::AddModuleScore() adds 1 as suffix
  ucell_modules <- paste0(module_names, "_UCell")                       # We added using name = "_UCell"
  
  # Find available modules in seurat object
  seurat_modules <- seurat_modules[seurat_modules %in% colnames(integrated_seurat@meta.data)]
  ucell_modules <- ucell_modules[ucell_modules %in% colnames(integrated_seurat@meta.data)]
  
  if (length(seurat_modules) > 0){
    plot_seurat(integrated_seurat, reduction = reduction, features = seurat_modules, filename = paste("Module_plot_Seurat", filename, sep = "_"), output_dir = output_dir, split_col = NULL)
  } else{
    log_warn(sample = "",
             step = "calc_module_scores",
             msg =glue::glue("Skipping Seurat module score plotting as no module were found."))
  }
  if (length(ucell_modules) > 0){
    plot_seurat(integrated_seurat, reduction = reduction, features = ucell_modules, filename = paste("Module_plot_UCell", filename, sep = "_"), output_dir = output_dir, split_col = NULL)
  } 
  else{
    log_warn(sample = "",
             step = "calc_module_scores",
             msg =glue::glue("Skipping UCell module score plotting as no module were found."))
  }
  
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  log_info(sample = "",
           step = "calc_module_scores",
           msg ="Module score calculated and Seurat object saved successfully.")
  return(invisible(integrated_seurat))
}

annotate_clusters <- function(integrated_seurat, reduction, assay,
                              AVG_PROB_THRESHOLD = 0.3,
                              output_dir){
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = integrated_seurat, reduction = reduction,
                  assay = assay, output_dir = output_dir)
  
  
  # ---- 🔄 Check and Update Reduction Name ----
  
  if (!(reduction %in% names(integrated_seurat@reductions))) {
    log_warn(sample = "",
             step = "annotate_clusters",
             msg =glue::glue("Reduction '{reduction}' is NOT present in Seurat object."))
    
    # Check alternative reduction name
    alt_reduction <- paste0("umap_", tolower(reduction))
    
    if (alt_reduction %in% names(integrated_seurat@reductions)) {
      log_info(sample = "",
               step = "annotate_clusters",
               msg =glue::glue("Using alternative reduction : '{alt_reduction}'."))
      reduction <- alt_reduction
    } else {
      log_error(sample = "",
                step = "annotate_clusters",
                msg =glue::glue("Alternative reduction '{alt_reduction}' is NOT present."))
    }
  }
  
  # ---- 🔎 Identify all cluster columns dynamically ----
  
  integration_methods <- c("cca", "rpca", "harmony", "jointpca")
  
  # Wrap the methods in parentheses and enforce the start (^) and end ($)
  pattern <- paste0("^(", paste(integration_methods, collapse = "|"), ")[0-9.]+$")
  
  # Select all columns that start with the prefix and end with a numeric cluster ID
  cluster_cols <- colnames(integrated_seurat@meta.data)
  cluster_cols <- cluster_cols[grepl(pattern, cluster_cols, ignore.case = TRUE)]
  
  if (length(cluster_cols) == 0) {
    log_error(sample = "",
              step = "annotate_clusters",
              msg =glue::glue("No clustering columns found with prefix : '{pattern}'."))
  }
  
  # ---- 🔎 Identify Module Score Columns ----
  
  # Find all columns ending in '.1' (default suffix from Seurat::AddModuleScore)
  score_cols <- colnames(integrated_seurat@meta.data)
  seurat_modules <- score_cols[base::grepl(pattern = "[0-9]+$", x = score_cols)]
  ucell_modules <- score_cols[base::grepl(pattern = "_UCell$", x = score_cols)]
  
  bad_prefixes <- c("cca", "rpca", "harmony", "jointpca")
  seurat_modules <- seurat_modules[!grepl(paste0("(", paste(bad_prefixes, collapse = "|"), ")"), seurat_modules)]
  ucell_modules <- ucell_modules[!grepl(paste0("(", paste(bad_prefixes, collapse = "|"), ")"), ucell_modules)]
  
  if (length(seurat_modules) < 2 & length(ucell_modules) < 2) {
    log_error(sample = "", 
              step = "annotate_cells", 
              msg ="Fewer than 2 module score columns (ending in '1' or '_UCell') found. Aborting.")
    return(integrated_seurat)
  }
  
  # ---- 1️⃣ Z-score & Softmax transformation of module scores ----
  
  # Extract the raw score matrix
  seurat_score_matrix <- integrated_seurat@meta.data[, seurat_modules, drop = FALSE] %>%
    as.matrix()
  ucell_score_matrix <- integrated_seurat@meta.data[, ucell_modules, drop = FALSE] %>%
    as.matrix()
  
  # Set negative module scores to 0 because they reflect lack of enrichment; 
  # keeping them would artificially inflate probabilities during min–max normalization.
  seurat_score_matrix[seurat_score_matrix < 0] <- 0
  ucell_score_matrix[ucell_score_matrix < 0] <- 0   #redundant as all ucell scores always > 0
  
  # Z-score standardization across all modules for each cell (row-wise)
  # The 'scale()' function performs Z-score transformation (mean=0, SD=1).
  # We transpose (t()) the matrix, apply scale(), and then transpose back
  # to ensure the scaling is applied across the genesets (rows in the original matrix).
  seurat_zscore_matrix <- t(scale(t(seurat_score_matrix)))
  ucell_zscore_matrix <- t(scale(t(ucell_score_matrix)))
  
  softmax <- function(x) {
    e_x <- exp(x - max(x))  # numerical stability
    e_x / sum(e_x)
  }
  
  # Apply row-wise softmax to convert raw module scores → probabilities
  seurat_prob_matrix <- t(base::apply(X = seurat_zscore_matrix, MARGIN = 1, FUN = softmax))
  ucell_prob_matrix <- t(base::apply(X = ucell_zscore_matrix, MARGIN = 1, FUN = softmax))
  
  # Replace NaN rows (all-zero input scores) with zeros
  seurat_prob_matrix[is.nan(seurat_prob_matrix)] <- 0
  ucell_prob_matrix[is.nan(ucell_prob_matrix)] <- 0
  
  # ---- 2️⃣ Weighted Cluster Consensus ----
  
  # Determine the cluster-assigned type based on mean probabilities of module scores
  # For each cluster, compute the mean probability per cell type (module),
  # then assign the type with the highest mean probability as the cluster's consensus.
  # Max_P stores this highest mean probability.
  
  score_cols_list <- list(Seurat = seurat_modules,
                          UCell = ucell_modules)
  
  prob_matrix_list <- list(Seurat = seurat_prob_matrix,
                           UCell = ucell_prob_matrix)
  
  # Initialize a list to store per-method / per-resolution data
  cluster_summary_list <- list()
  
  for (method in c("Seurat", "UCell")){
    
    score_cols <- score_cols_list[[method]]
    prob_matrix <- prob_matrix_list[[method]]
    
    temp_df <- integrated_seurat@meta.data[, c(score_cols, cluster_cols)]
    temp_df[, score_cols] <- prob_matrix       # Use the calculated probabilities
    
    # Iterate through resolutions
    for (col in cluster_cols) {
      
      # Compute mean probabilities per cluster
      avg_probs <- temp_df %>%
        dplyr::group_by(.data[[col]]) %>%
        dplyr::summarise(across(.cols = all_of(score_cols), 
                                .fns = mean),
                         .groups = "drop")
      
      # Set mean probabilities below 0.3 to 0
      avg_probs[score_cols] <- lapply(X = avg_probs[score_cols],
                                      FUN = function(x) ifelse(x < AVG_PROB_THRESHOLD, 0, x))
      
      # Compute row-wise max
      max_vals <- do.call(pmax, c(avg_probs[, score_cols], na.rm = TRUE))
      avg_probs$Max_P <- base::round(max_vals, digits = 2)
      
      # Assign Assigned_Type
      avg_probs$Assigned_Type <- apply(avg_probs[, score_cols], 1, function(row) {
        if (max(row) == 0) {
          "Unknown"
        } else {
          paste(names(row)[row == max(row)], collapse = ",")
        }
      })
      
      # Keep only relevant columns
      avg_probs <- avg_probs %>%
        dplyr::select(any_of(col), Max_P, Assigned_Type)
      
      # Clean names for assignment (e.g., "T.Cell_UCell" -> "T.Cell")
      cleaned_modes <- base::gsub(pattern = "_UCell$", replacement = "", x = avg_probs$Assigned_Type )
      avg_probs$Assigned_Type <- base::gsub(pattern = "[0-9]+$", replacement = "", x = cleaned_modes)
      
      # Convert to named vectors for easy lookup
      assigned_vec <- stats::setNames(object = avg_probs$Assigned_Type, 
                                      nm = avg_probs[[col]])
      
      meanprob_vec <- stats::setNames(object = avg_probs$Max_P, 
                                      nm = avg_probs[[col]])
      # ---------------------------------------------------------------------------- #
      # > assigned_vec
      # 1              2      3     4               5             6
      # "Epithelial"  "TNK"  "TNK"  "Fibroblasts"   "Epithelial"  "Dendritic"
      # 7                         8               9 
      # "MyocytesMyofibroblasts"  "Endothelial"   "Proliferating"
      # 10             11           12            13          14         15
      # "Epithelial"  "Dendritic"   "Monocytes"   "BPlasma"   "BPlasma"  "Unknown"
      # 16              17
      # "Epithelial"   "Erythrocytes"
      # ---------------------------------------------------------------------------- #
      
      # Assign values to each cell using the cluster column
      
      # NOTE: Lookup vector returns values in the exact order of the keys you provide
      # keys = cell_clusters
      # lookup vectors = assigned_vec and meanprob_vec
      
      cell_clusters <- integrated_seurat@meta.data[[col]]
      cluster_assignment <- assigned_vec[cell_clusters]
      cluster_meanprob   <- meanprob_vec[cell_clusters]
      
      # ---------------------------------------------------------------------------- #
      # > length(cell_clusters)
      # [1] 43727
      
      # > length(cluster_assignment)
      # [1] 43727
      
      # > length(cluster_meanprob)
      # [1] 43727
      
      # > head(cell_clusters)
      # [1] 2  5  12 5  1  4
      # Levels: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17
      
      # > head(cluster_assignment)
      # 2       5             12           5              1             4
      # "TNK"  "Epithelial"   "Monocytes"  "Epithelial"   "Epithelial"  "Fibroblasts"
      
      # > head(cluster_meanprob)
      # 2       5       12      5       1       4
      # 0.48    0.62    0.34    0.62    0.31    0.59
      # ---------------------------------------------------------------------------- #
      
      # Set names as barcodes
      cluster_assignment <- stats::setNames(object = cluster_assignment,
                                            nm = base::rownames(integrated_seurat@meta.data))
      
      cluster_meanprob <- stats::setNames(object = cluster_meanprob,
                                          nm = base::rownames(integrated_seurat@meta.data))
      
      # Store in Seurat metadata
      integrated_seurat[[paste(method, col, sep = "_")]] <- cluster_assignment
      integrated_seurat[[paste(method, col, "MeanProb", sep = "_")]] <- cluster_meanprob
      
      # Store in list
      cluster_summary_list[[paste(method, col, sep = "_")]] <- avg_probs
    }
  }
  
  # Combine all into single dataframe and save as xlsx
  cluster_summary_df <- bind_rows(cluster_summary_list)
  wb <- createWorkbook()
  addWorksheet(wb, "Cluster_Summary")
  writeData(wb, "Cluster_Summary", cluster_summary_df)
  saveWorkbook(wb, file = file.path(output_dir, "Cluster_Summary.xlsx"), overwrite = TRUE)
  
  log_info(sample = "", step = "annotate_clusters", 
           msg =glue::glue("Completed weighted cluster consensus across {length(cluster_cols)} resolutions."))
  
  # ---- 3️⃣ Confidence Scoring per Cell  ----
  
  consensus_cols <- colnames(integrated_seurat@meta.data)
  consensus_cols <- consensus_cols[grepl(pattern = paste0("^(Seurat|UCell)_(", paste(cluster_cols, collapse = "|"), ")$"), 
                                         x = consensus_cols)]
  
  consensus_mat <- integrated_seurat@meta.data[, consensus_cols, drop = FALSE]
  
  # Calculate Final Cell Type: mode across resolutions  and store in Seurat metadata
  # NOTE: use tabulate for extremely fast calculation
  cell_mode_assignment <- base::apply(X = consensus_mat, MARGIN = 1, function(x) {
    
    ux <- stats::na.omit(x)                                  # Remove NAs first
    if(length(ux) == 0) return(NA)                          # If no assignments are available, return NA
    ux_tab <- base::tabulate(base::match(ux, unique(ux)))   # tabulate finds frequency counts,
    # match finds positions of first occurrence
    return(unique(ux)[which.max(ux_tab)])                   # Return value corresponding to max frequency
  })
  
  # Clean names for assignment (e.g., "T.Cell_UCell" -> "T.Cell")
  cleaned_modes <- base::gsub(pattern = "_UCell$", replacement = "", x = cell_mode_assignment)
  cleaned_modes <- base::gsub(pattern = "[0-9]+$", replacement = "", x = cleaned_modes)
  integrated_seurat$CellType <- cleaned_modes
  
  # Create a matrix where each column is the Mode assignment, replicated.
  Mode_Matrix <- matrix(cell_mode_assignment, 
                        nrow = nrow(consensus_mat), 
                        ncol = ncol(consensus_mat), 
                        byrow = FALSE)
  
  # Stability = (Count of matches to the Mode) / (Total resolutions)
  stability_score <- base::rowSums(consensus_mat == Mode_Matrix, na.rm = TRUE) / ncol(consensus_mat)
  stability_score <- base::round(stability_score, digits = 2)
  
  # Store results
  integrated_seurat$Stability_Score <- stability_score
  
  log_info(sample = "", 
           step = "annotate_clusters", 
           msg ="Calculated Stability Scores and Consensus CellType.")
  
  # ---- 4️⃣ UMAP Plots for Cluster Visualization ----
  
  features <- c("CellType", "harmony0.8", "Stability_Score")
  plot_seurat(integrated_seurat, reduction = reduction, features = features, filename = "UMAP.Ann", output_dir = output_dir, split_col = NULL)
  features <- c("Stability_Score")
  plot_seurat(integrated_seurat, reduction = reduction, features = features, filename = "UMAP.Stability", output_dir = output_dir, split_col = features)
  
  # ---- 💾 Save Seurat Object ----
  
  # Determine file name based on assay
  if (assay == "RNA"){
    filename <- "integrated_seurat_ann.rds"
    filename_clean <- "integrated_seurat_clean_ann.rds"
  } else{   # For Spatial.008um and Spatial.016um assays
    filename <- paste0("integrated_seurat_", assay, "_ann.rds")
    filename_clean <- paste0("integrated_seurat_clean_", assay, "_ann.rds")
  }
  
  # Filter out bad data
  subset_seurat <- base::subset(x = integrated_seurat, Stability_Score >= 0.9)
  
  # Save rds
  base::saveRDS(object = integrated_seurat, file = file.path(output_dir, filename))
  base::saveRDS(object = subset_seurat,     file = file.path(output_dir, filename_clean))
  
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  log_info(sample = "", 
           step = "annotate_clusters", 
           msg ="Consensus annotation completed successfully. New metadata columns added: PredictedCellType, Pmax_Seurat, Consensus_*, CellType, Stability_Score, Confidence.")
  
  return(invisible(integrated_seurat))
}

consolidate_gold_standard_markers <- function(data_dir = NULL, K_MODULES = NULL){
  
  # ---- 📁 0. Initial Setup and Data Loading ----
  
  # Define the directory where your scRNA-seq marker gene Excel files are located.
  data_dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/scrna_markers/"
  
  # List all files in the directory.
  marker_files <- list.files(data_dir, full.names = TRUE, pattern = "\\.xlsx$")
  
  # Initialize an empty data frame to store all marker results.
  raw_marker_df <- data.frame()
  
  # Loop through each Excel file to load and process data
  message("✨ Loading and Initial Processing of Marker Files...")
  for (f in marker_files){
    df <- read.xlsx(f) %>%
      # Convert gene names to uppercase for consistency in mouse and human
      dplyr::mutate(gene = toupper(gene),
                    # Extract dataset name from the filename
                    dataset = gsub("Markers.All.harmony0.[0-9]_|\\.xlsx", "", basename(f))) %>%
      # Select and rename essential columns
      dplyr::select(dataset, cluster, gene, p_val_adj, avg_log2FC, pct.1, pct.2, ratio)
    
    # Combine data from all files
    raw_marker_df <- dplyr::bind_rows(raw_marker_df, df)
  }
  message(paste0("   - Total rows loaded: ", nrow(raw_marker_df)))
  
  # ---- 🔎 1. Quality Control and Consensus Gene Filtering ----
  
  consensus_marker_df <- raw_marker_df %>%
    
    # Filter for statistically significant markers
    dplyr::filter(p_val_adj <= 0.05) %>%
    
    # Calculate ranks for three key metrics within each cluster-dataset combination
    dplyr::group_by(dataset, cluster) %>%
    dplyr::mutate(rank_pct1     = rank(desc(pct.1),      ties.method = "first"),
                  rank_log2fc   = rank(desc(avg_log2FC), ties.method = "first"),
                  rank_ratio    = rank(desc(ratio),      ties.method = "first")) %>%
    
    # Keep genes that are top 50 in ANY of the 3 metrics (local filter)
    dplyr::filter(rank_pct1 <= 50 | rank_log2fc <= 50 | rank_ratio <= 50) %>%
    ungroup() %>%
    
    # Count the number of datasets where the gene is a top marker
    dplyr::group_by(gene) %>%
    dplyr::mutate(ndatasets = n_distinct(dataset)) %>%
    
    # Keep only genes that appear in at least 5 different datasets (consensus filter)
    dplyr::filter(ndatasets >= 5)
  
  message(paste0("   - Number of consensus genes retained: ", n_distinct(consensus_marker_df$gene)))
  message(paste0("   - Total gene-condition pairs: ", nrow(consensus_marker_df)))
  
  # ---- 🧬 2. Data Structuring for Clustering ----
  
  # Create the Gene x Condition matrix. Rows are genes, columns are conditions.
  # Values are the log2 Fold Change (avg_log2FC). 
  gene_log2fc_matrix <- consensus_marker_df %>%
    # Create a unique identifier for each condition (Dataset_Cluster)
    dplyr::mutate(dataset_cluster = paste(dataset, cluster, sep = "_")) %>%
    # Reshape the data from long to wide format
    tidyr::pivot_wider(id_cols = gene,
                       names_from = dataset_cluster,
                       values_from = avg_log2FC,
                       values_fill = 0) %>% # Fill non-existent values with 0
    tibble::column_to_rownames('gene') %>%  # Set gene names as row names
    as.matrix()
  
  message(paste0("   - Matrix dimensions (Genes x Conditions): ", paste(dim(gene_log2fc_matrix), collapse = " x ")))
  message("--------------------------------------------------")
  
  # ---- 🎯 3. GENE CLUSTERING APPROACHES ----
  
  # 🚀 Approach 1: Z-Score Scaling + Euclidean Distance (The 'Shape' Approach)
  # This approach focuses on the relative up/down regulation pattern (shape)
  # by normalizing the magnitude across genes.
  
  # 3.1. Z-Score Scaling
  # scale each row so each gene has mean = 0, sd = 1 across datasets+clusters. 
  # - scale() centers and scales columns by default (center = TRUE, scale = TRUE)
  #      -> Subtracts the mean of each column across all rows.
  #      -> Divides each column by its standard deviation.
  #      -> This ensures all columns contribute equally, regardless of their 
  #         magnitude or raw variance.
  scaled_log2fc_matrix <- gene_log2fc_matrix %>%
    t() %>%         # Transpose: Genes are now columns
    scale() %>%     # Z-score each gene's profile (column)
    t()             # Transpose back: Genes are rows again
  
  message("🚀 Approach 1: Z-Score + Euclidean Distance (Ward's Method)")
  message(paste0("   - Scaled matrix dimensions: ", paste(dim(scaled_log2fc_matrix), collapse = " x ")))
  
  # 3.2. Distance Calculation (Euclidean)
  # dist(x) always computes pairwise distances between the rows of x
  gene_dist_euclidean <- stats::dist(scaled_log2fc_matrix, method = "euclidean")
  cond_dist_euclidean <- stats::dist(t(scaled_log2fc_matrix), method = "euclidean")
  
  # 3.3. Hierarchical Clustering using distance object (Ward's Method)
  # Ward's method ('ward.D2') is often preferred for gene data as it aims 
  # to minimize the variance within clusters, leading to tight, compact groups.
  gene_hclust_euclidean <- stats::hclust(gene_dist_euclidean, method = "ward.D2")
  cond_hclust_euclidean <- stats::hclust(cond_dist_euclidean, method = "ward.D2")
  
  message("   - Clustering complete for Approach 1.")
  
  # # 🔗 Approach 2: Pearson Correlation Distance (The 'Trend' Approach) 
  # # This approach directly measures the linear relationship between gene profiles,
  # # ignoring the raw magnitude difference.
  # 
  # message("\n🔗 Approach 2: Pearson Correlation Distance (Average Method)")
  # 
  # # 3.4. Correlation Matrix Calculation (Pearson correlation)
  # # cor(x) always computes correlation between the columns of x
  # # NOTE: why use = "pairwise.complete.obs"?
  # #   - "everything"	           Default. If any NA exists, correlation = NA
  # #   - "complete.obs"	         Only use rows where all values are non-NA. Can discard a lot of data.
  # #   - "pairwise.complete.obs"	 Use all available pairs, ignore missing values for that gene pair
  # gene_cor_matrix <- gene_log2fc_matrix %>%
  #   t() %>%                              # Transpose: Genes are now columns
  #   cor(.,
  #       use = "pairwise.complete.obs",
  #       method = "pearson")
  # 
  # # 3.5. Convert Correlation to Distance (D = 1 - R)
  # # Correlation 1 → distance 0 (very close)
  # # Correlation 0 → distance 1 (medium)
  # # Correlation -1 → distance 2 (farthest)
  # gene_dist_cor_based <- 1 - gene_cor_matrix
  # 
  # # Convert the full symmetric matrix into a 'dist' object for hclust
  # # as.dist() converts a symmetric matrix into a compact format storing only the 
  # # lower triangle, which hclust requires.
  # gene_dist_cor_based <- stats::as.dist(gene_dist_cor_based)
  # 
  # # 3.6. Hierarchical Clustering (Average Method)
  # # 'average' (UPGMA) is standard for correlation-based distance.
  # gene_hclust_cor_based <- stats::hclust(gene_dist_cor_based, method = "average")
  # 
  # message("   - Clustering complete for Approach 2.")
  # message("--------------------------------------------------")
  
  # ---- 📝 4. Final Output Formatting and Export ----
  
  # We will use the results from Approach 1 for the final heatmap reordering.
  
  # 4.1. Reorder Matrix and Assign Modules
  final_gene_order <- gene_hclust_euclidean$order
  final_cond_order <- cond_hclust_euclidean$order
  
  # Reorder the original log2FC matrix based on both gene and condition clustering
  reordered_log2fc_df <- gene_log2fc_matrix[final_gene_order, final_cond_order] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene")
  
  # Cut the gene tree into k=15 modules (arbitrary choice for visualization)
  if (is.null(K_MODULES)){
    K_MODULES <- 15
  }
  gene_modules <- stats::cutree(gene_hclust_euclidean, k = K_MODULES)
  
  # 4.2. Final Data Frame Assembly and Cleanup
  final_output_df <- reordered_log2fc_df %>%
    # Add module assignment
    dplyr::left_join(y = tibble::tibble(gene = names(gene_modules), module = gene_modules),
                     by = c("gene" = "gene")) %>%
    # Move 'module' column to the front for clarity
    dplyr::select(gene, module, dplyr::everything()) %>%
    # Replace all 0s (imputed for NA) with true NA and round numeric values
    dplyr::mutate(dplyr::across(.cols = where(is.numeric) & !module, .fns = ~ dplyr::na_if(., 0))) %>%
    dplyr::mutate(dplyr::across(.cols = where(is.numeric) & !module, .fns = ~ base::round(., 1)))
  
  message(paste0("📌 Final matrix reordered and genes assigned to ", K_MODULES, " modules."))
  
  # 4.3. Export to Excel
  wb <- createWorkbook()
  addWorksheet(wb, "log2fc_Clustered")
  writeData(wb, sheet = "log2fc_Clustered", final_output_df)
  file_name <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/scrna_markers_results_final.xlsx"
  saveWorkbook(wb, file_name, overwrite = TRUE)
  
  message(paste0("\n✅ Success! Clustered log2FC matrix saved to:\n", file_name))
  
}

plot_gold_standard_markers <- function(integrated_seurat, reduction = "Harmony", 
                                       marker_file, output_dir){
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = integrated_seurat, reduction = reduction,
                  xlsx_file = marker_file, output_dir = output_dir)
  
  # ---- 🔄 Check and Update Reduction Name ----
  
  if (!(reduction %in% names(integrated_seurat@reductions))) {
    log_warn(sample = "", 
             step   = "plot_gold_standard_markers",
             msg    = glue::glue("Reduction '{reduction}' is NOT present in Seurat object."))
    
    # Check alternative reduction name
    alt_reduction <- base::paste0("umap_", base::tolower(reduction))
    
    if (alt_reduction %in% names(integrated_seurat@reductions)) {
      log_info(sample = "", 
               step   = "plot_gold_standard_markers",
               msg    = glue::glue("Using alternative reduction: '{alt_reduction}'."))
      reduction <- alt_reduction
    } else {
      log_error(sample = "", 
                step   = "plot_gold_standard_markers",
                msg    = glue::glue("Alternative reduction '{alt_reduction}' is NOT present."))
    }
  }
  
  # ---- 🧪 Detect Appropriate Assay for FeaturePlot ----
  
  all_assays <- names(integrated_seurat@assays)
  
  if ("SCT" %in% all_assays) {
    active_assay <- "SCT"
  } else if ("RNA" %in% all_assays) {
    active_assay <- "RNA"
  } else if (length(all_assays) > 0) {
    active_assay <- all_assays[1]
    log_info(sample = "", 
             step   = "plot_gold_standard_markers",
             msg    = glue::glue("No SCT/RNA assay found. Using assay: '{active_assay}'."))
  } else {
    log_error(sample = "", 
              step   = "plot_gold_standard_markers",
              msg    = "No assays found in Seurat object.")
  }
  
  # Set active assay
  Seurat::DefaultAssay(integrated_seurat) <- active_assay
  
  # ---- 📥 Create Feature List from Marker file ----
  
  # Load marker file
  marker_df <- tryCatch({
    openxlsx::read.xlsx(xlsxFile = marker_file)
  }, error = function(e) {
    log_error(sample = "", 
              step   = "plot_gold_standard_markers",
              msg    = glue::glue("Failed to read marker file: '{e$message}'."))
  })
  
  # Initialize an empty list to store all cell type signatures
  signatures_list <- list()
  
  # Iterate through each column (cell type) in the marker dataframe
  for (i in base::seq_len(base::ncol(marker_df))) {
    
    # Define name of modules
    module_name <- base::make.names(colnames(marker_df)[i])
    
    # Determine features from the marker_df
    xlsx_features <- marker_df[[i]] %>%
      stats::na.omit() %>%
      base::as.vector()
    
    # Determine features present in data set
    present_features <- base::rownames(SeuratObject::GetAssayData(object = integrated_seurat, 
                                                                  assay = active_assay, 
                                                                  layer = "data"))
    
    # Match features (case-insensitive)
    features <- present_features[base::tolower(present_features) %in% base::tolower(xlsx_features)]
    
    # Filter to unique and sorted features
    if (length(features) >= 2) {
      signatures_list[[module_name]] <- sort(unique(features))
    } else {
      log_warn(sample = "", 
               step   = "plot_gold_standard_markers",
               msg    = glue::glue("Skipping module '{module_name}'. Fewer than 2 matching markers found."))
    }
  }
  
  # ---- 🖼️ Plot UMAPs for each module ----
  
  for (module_name in names(signatures_list)) {
    
    # Log progress
    log_info(sample = "", 
             step   = "plot_gold_standard_markers",
             msg    = glue::glue("Plotting genes for module: {module_name}"))
    
    # Call the core plotting utility
    plot_seurat(integrated_seurat, 
                reduction  = reduction, 
                features   = signatures_list[[module_name]], 
                raster     = FALSE,
                filename   = paste("Module_plot_Seurat", proj, module_name, sep = "_"),  
                output_dir = output_dir, 
                split_col  = NULL)
  }
  
  # ---- 🪵 Log Output and Return ----
  
  log_info(sample = "", 
           step   = "plot_gold_standard_markers", 
           msg    = "Gold standard marker plotting completed successfully.")
  
  return(invisible(integrated_seurat))
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

annotate_manual <- function(integrated_seurat, assay, clusters, res, reduction, output_dir){
  
  set.seed(1234)
  
  idents <- paste0("cluster.", res, ".", base::tolower(reduction))
  
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
  
  # Proceed with annotation ONLY if all clusters have been annotated
  if (identical(list_1, list_2)){
    print("All Clusters have been annotated")
    
    # Extract metadata from Seurat object, assign appropriate res to
    # seurat_clusters column and add Cell.Type, Cell.Subtype columns
    data <- integrated_seurat@meta.data %>% 
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
    integrated_seurat@meta.data <- data
    
    # Create .rds object for integrated seurat object
    if (assay == "RNA"){
      saveRDS(integrated_seurat, file=file.path(output_dir, "integrated_seurat.ann.rds"))
    } else{
      # For Spatial.008um and Spatial.016um assays
      saveRDS(integrated_seurat, file=file.path(output_dir, paste0("integrated_seurat.", assay, ".ann.rds")))
    }
    
  } else {
    cat("\nYou missed annotating these clusters:\t", setdiff(list_1, list_2))
    cat("\nThese clusters are not present in data:\t", setdiff(list_2, list_1))
    cat("\nThese clusters have duplicate annotation:\t", list_2[duplicated(list_2)])
  }
  
  return(integrated_seurat)
} 

plot_dot_plot <- function(integrated_seurat, idents, features, filename, output_dir, gene_y_axis = FALSE, split_col = NULL){
  
  set.seed(1234)
  
  # Re-order the active ident alphabetically
  Idents(integrated_seurat) <- idents
  Idents(integrated_seurat) <- base::factor(x = integrated_seurat@active.ident, 
                                            levels= sort(levels(integrated_seurat@active.ident)))
  
  # Determine grouping variable
  if (is.null(split_col)) {
    groups <- "All"
  } else {
    groups <- unique(integrated_seurat@meta.data[[split_col]])
  }
  
  # Loop over groups
  plot_list <- list()
  for (g in groups) {
    if (g == "All") {
      subset_obj <- integrated_seurat
      plot_title <- "All"
    } else {
      subset_obj <- integrated_seurat[, integrated_seurat@meta.data[[split_col]] == g]
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
  n_y <- length(unique(integrated_seurat@active.ident))  
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
  ggsave(filename = file.path(output_dir, paste0(filename,".pdf")),
         plot = combined_plot,
         width = total_width,
         height = total_height,
         limitsize = FALSE,
         units = "in",
         bg        = "white")
}

plot_dot_custom <- function(integrated_seurat, assay, proj.params, ident.1, ident.2, features, filename, output_dir, CellType.col = "Cell.Type", split_col=NULL){
  
  set.seed(1234)
  
  # Keep only cell types present in your dataset and marker file
  marker_df <- openxlsx::read.xlsx(file.path(proj.params$cell.type.marker.file))
  celltypes_in_seurat <- unique(integrated_seurat[[CellType.col]])
  celltypes_to_keep <- intersect(colnames(marker_df), celltypes_in_seurat)
  
  # Determine groups based on split_col
  groups <- if (!is.null(split_col) && split_col %in% colnames(integrated_seurat@meta.data)) {
    as.character(unique(integrated_seurat@meta.data[[split_col]]))
  } else {
    "All"
  }
  
  # Determine features present in data set
  features <- features[base::tolower(features) %in% base::tolower(SeuratObject::Features(integrated_seurat))] %>% 
    na.omit() %>%
    as.vector()
  
  if(length(features) == 0) stop("No matching features found in the assay.")
  
  # Loop over groups
  plot_list <- list()
  for (g in groups) {
    
    # Subset object
    if (g == "All") {
      subset_obj <- integrated_seurat
      plot_title <- ""
    } else {
      subset_obj <- integrated_seurat[, integrated_seurat@meta.data[[split_col]] == g]
      plot_title <- g
    }
    
    # Determine levels based on ident.1
    levels <- if (!is.null(ident.1) && ident.1 %in% colnames(integrated_seurat@meta.data)) {
      subset_obj@meta.data[[ident.1]] %>% unique() %>% as.character()
    } else {
      stop(ident.1, " missing in metadata")
    }
    
    # Determine levels.extra based on ident.2
    levels.extra <- if (!is.null(ident.2) && ident.2 %in% colnames(integrated_seurat@meta.data)) {
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
    
    # Calculate panel dimensions
    height_per_element <- 0.5
    width_per_element <- 0.5
    height_panel <- length(unique(dotplot_data[[y_var]])) * height_per_element + 2
    width_panel <- length(unique(dotplot_data[[x_var]])) * width_per_element + 6
    x_label_size <- 72 * width_per_element * 0.5
    y_label_size <- 72 * height_per_element * 0.5  # 0.5 scaling so it doesnt get too large
    
    # Create ggplot
    p <- ggplot(data = dotplot_data, 
                mapping = aes(x = .data[[x_var]], y = .data[[y_var]], 
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
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = x_label_size ),
            axis.text.y = element_text(size = y_label_size),
            legend.text  = element_text(size = 0.8 * y_label_size),   # slightly smaller than y-axis labels
            legend.title = element_text(size = y_label_size))
    
    plot_list[[plot_title]] <- p
  }
  
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
  ggsave(filename = file.path(output_dir, paste0(filename,".pdf")),
         plot = combined_plot,
         width = total_width,
         height = total_height,
         limitsize = FALSE,
         units = "in",
         bg = "white")
}

# Input is seurat object of a single slide with columns X, Y, Sample, Group
plot_spatial_map <- function(plot.seurat, x1, y1, x2, y2, suffix, output_dir){
  
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
                  path=output_dir,
                  width=25,
                  height=8.5,
                  units=c("in"),
                  dpi=600,
                  limitsize=TRUE,
                  bg="white")
}  

tabulate_frequency <- function(integrated_seurat, split_cols, output_dir){
  
  set.seed(1234)
  
  # Create a new workbook
  wb <- createWorkbook()
  
  for (split_col in split_cols){
    if(!is.null(integrated_seurat@meta.data[[split_col]] %>% unique())){
      
      counts <- integrated_seurat@meta.data %>% 
        dplyr::count(.data[[split_col]], Cell.Type) %>% 
        dplyr::group_by(.data[[split_col]]) %>% 
        dplyr::mutate(Percent = round(100 * n / sum(n), 2)) %>%
        tidyr::pivot_wider(names_from = .data[[split_col]], values_from = c(n, Percent)) %>%
        as.data.frame()
      
      
      # Make sheet name safe
      sheet_name <- substr(gsub("[/\\?<>\\:*|\"]", "_", split_col), 1, 31)  # Excel sheet name limit = 31
      addWorksheet(wb, sheetName = sheet_name)
      writeData(wb, sheet = sheet_name, counts)
    }
  }
  
  # Save workbook
  saveWorkbook(wb, file = file.path(output_dir, "Population_Frequencies.xlsx"), overwrite = TRUE)
}

prep_pseudobulk <- function(integrated_seurat, comparison.col, assay = "RNA"){
  
  set.seed(1234)
  
  # Input checks
  if (!inherits(integrated_seurat, "Seurat")) {
    stop("'integrated_seurat' must be a Seurat object.")
  }
  if (!comparison.col %in% colnames(integrated_seurat@meta.data)) {
    stop(paste0("'", comparison.col, "' not found in metadata columns."))
  }
  if (!"Sample" %in% colnames(integrated_seurat@meta.data)) {
    stop("'Sample' column is required in metadata.")
  }
  if (!assay %in% names(integrated_seurat@assays)) {
    stop(paste0("Assay '", assay, "' not found in Seurat object."))
  }
  
  # Initialize storage
  counts_matrix <- integrated_seurat@assays[[assay]]$counts
  metadata_full <- data.frame()
  read_data_full <- data.frame(SYMBOL = rownames(counts_matrix))
  
  # Extract unique groups from comparison.col that will be compared in DE analysis
  groups <- unique(integrated_seurat@meta.data[[comparison.col]])
  groups <- groups[!is.na(groups)]  # remove NA groups if present
  
  # Loop through each group
  for (g in groups) {
    subset.seurat <- integrated_seurat[, integrated_seurat@meta.data[[comparison.col]] == g]
    
    # Generate metadata
    metadata <- subset.seurat@meta.data %>%
      tibble::rownames_to_column("barcodes") %>%
      dplyr::filter(!is.na(Sample)) %>%
      dplyr::mutate(Sample_ID = paste0(Sample, "_", g),
                    Comparisons = .data[[comparison.col]]) %>%
      dplyr::add_count(Sample_ID) %>%
      dplyr::filter(n >= 100)
    
    ### Generate read data
    # read data will have "the reads of all cells belonging to a single sample" 
    # merged together in each column. 
    
    # First, create a list of samples
    samples <-  metadata[["Sample_ID"]] %>%
      unique()
    
    # Second, initialize read counts matrix
    read_data <- matrix(0,
                        nrow = nrow(counts_matrix),
                        ncol = length(samples),
                        dimnames = list(rownames(counts_matrix), samples))
    
    # Third, add row-wise, the counts of each gene for each sample
    for(i in samples){
      
      # Create a list of cells for each sample
      cells_subset <- metadata %>% 
        dplyr::filter(Sample_ID == i) %>% 
        dplyr::pull(barcodes)
      
      # Subset counts
      subset_counts <- counts_matrix[, cells_subset, drop=FALSE]
      
      # Sum counts
      read_data[, i] <- Matrix::rowSums(subset_counts)
    }
    
    # Fourth, format metadata and readdata
    cols_to_keep <- intersect(c("Sample_ID", "Comparisons", "Patient", "Condition", "Treatment", "Disease", 
                                "Tissue", "Strain", "Cell_line", "Sex", "n", "barcodes"),
                              colnames(metadata))
    metadata <- metadata %>%
      dplyr::distinct(Sample_ID, .keep_all = TRUE) %>%
      dplyr::select(all_of(cols_to_keep))
    
    read_data <- as.data.frame(read_data) %>% 
      tibble::rownames_to_column("SYMBOL")
    
    # Finally, merge results
    metadata_full <- dplyr::bind_rows(metadata_full, metadata)
    read_data_full <- dplyr::left_join(read_data_full, read_data, by=c("SYMBOL"="SYMBOL"))
  }
  
  return(list(metadata = metadata_full, read_data = read_data_full))
}

find_degs_seurat <- function(integrated_seurat, comparison.col, output_dir, celltype.col = "Cell.Type", assay = "RNA"){
  
  set.seed(1234)
  
  # Setup log file
  log_file <- file.path(output_dir, "find_degs_seurat.log")
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
  celltype_levels <- unique(integrated_seurat@meta.data[[celltype.col]])
  comparison_levels <- unique(integrated_seurat@meta.data[[comparison.col]])
  comparisons <- utils::combn(x = comparison_levels, m = 2, simplify = FALSE)
  
  for(celltype in celltype_levels){
    for(comparison in comparisons){
      
      target <- comparison[1]
      reference <- comparison[2]
      
      # Subset Seurat object
      subset_obj <- integrated_seurat[, integrated_seurat@meta.data[[celltype.col]] == celltype]
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
      
      raw_counts <- as.data.frame(as.matrix(SeuratObject::GetAssayData(deg_obj, assay = assay, layer = "counts"))) 
      counts_list[[paste0(celltype, ".", target, ".vs.", reference)]] <- raw_counts %>%
        tibble::rownames_to_column(var = "SYMBOL")
      
      # VST normalized counts
      metadata <- deg_obj@meta.data[colnames(raw_counts), , drop = FALSE]
      rownames(metadata) <- colnames(raw_counts)
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = raw_counts,
                                            colData = metadata,
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
  save_list_to_xlsx <- function(lst, file_name) {
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
    openxlsx::saveWorkbook(wb, file_name, overwrite = TRUE)
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
  openxlsx::saveWorkbook(wb, file.path(output_dir, "Seurat_DEGs.xlsx"), overwrite = TRUE)
  
  
  if (length(counts_list) > 0) save_list_to_xlsx(counts_list, file.path(output_dir, "Seurat_Raw_Counts.xlsx"))
  if (length(vst_counts_list) > 0) save_list_to_xlsx(vst_counts_list, file.path(output_dir, "Seurat_VST_Counts.xlsx"))
  
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

h5ad_to_seurat_batch <- function(path_to_h5ad) {
  
  # Find all .h5ad files (path must be a directory)
  h5ad_files <- list.files(path = path_to_h5ad,
                           pattern = "\\.h5ad$",
                           full.names = TRUE)
  
  if(length(h5ad_files) == 0) {
    stop("No .h5ad files found in the specified path.")
  }
  
  seurat_list <- list()
  
  for(h5ad_file in h5ad_files) {
    
    message("Processing file: ", h5ad_file)
    
    # Prepare .h5seurat filename
    h5seurat_file <- sub("\\.h5ad$", ".h5seurat", h5ad_file)
    
    # Convert .h5ad to .h5seurat
    if(!file.exists(h5seurat_file)) {
      message("Converting .h5ad to .h5seurat...")
      # IMPORTANT CORRECTION: 'dest' must be the full path to the output file
      # or simply the filename if the working directory is correct, but passing 
      # the full path to the source is safer.
      SeuratDisk::Convert(source = h5ad_file,
                          dest = "h5seurat", 
                          assay = "RNA",
                          overwrite = FALSE)
    } else {
      message(".h5seurat file already exists: ", h5seurat_file)
    }
    
    # Load Seurat object
    message("Loading Seurat object...")
    # If you get error, "Error: Missing required datasets 'levels' and 'values'",
    # set meta.data = FALSE, misc = FALSE within SeuratDisk::LoadH5Seurat()
    seurat_obj <- SeuratDisk::LoadH5Seurat(file = h5seurat_file, 
                                           meta.data = FALSE, 
                                           misc = FALSE)
    
    # Extract and Assign Metadata (adata.obs to seurat_obj@meta.data)
    message("Extracting and assigning metadata from AnnData (.obs) to Seurat (@meta.data)...")
    adata <- anndata::read_h5ad(h5ad_file)
    obs_df <- as.data.frame(adata$obs)
    # adata.obs.to_excel("/hpc/home/kailasamms/scratch/scRNASeq_PanCancer/Metadata.xlsx", index=True)
    # seurat_obj@meta.data <- read.xlsx("/hpc/home/kailasamms/scratch/scRNASeq_PanCancer/Metadata.xlsx") %>%
    #   tibble::column_to_rownames("X1")
    
    # Ensure rownames match Seurat object cells
    if (!all(rownames(obs_df) %in% colnames(seurat_obj))) {
      warning("Cell names in metadata do not match Seurat object. Check rownames!")
    }
    seurat_obj@meta.data <- obs_df
    
    # Save Seurat object as RDS
    rds_file <- base::gsub(pattern = "\\.h5ad$", replacement = ".rds", x = h5ad_file)
    message("Saving Seurat object to RDS: ", rds_file)
    saveRDS(object = seurat_obj, file = rds_file)
    
    seurat_list[[basename(rds_file)]] <- seurat_obj
  }
  
  message("All .h5ad files processed. Returning list of Seurat objects.")
  return(seurat_list)
}

# ---- ⏳ SURVIVAL RELATED FUNCTIONS ----

# NOTE: When plotting KM curves for individual genes, non-transformed data,
# log-transformed or median centered log-transformed data give identical results
# for all cutoff methods except thirds. RECOMMEND using vst counts from DESeq2 or
# log transformed counts

# When plotting KM curves for signature scores, RECOMMEND using vst counts from 
# DESeq2 or log transformed counts. advanced_z() automatically does the 
# necessary median centering before signature score calculation.

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
  return(data.frame(scenarios))
}

survival_params <- list(
  
  # Stratification (Expression + Metadata-based survival)
  stratify_var     = NULL,          # one or more genes or metadata columns
  substratify_var  = NULL,          # optional metadata column for sub-stratification
  facet_var        = NULL,          # optional faceting variable
  
  # Cutoff settings (ONLY for Expression-based survival)
  cutoff_method    = "optimal",      # mean, median, quartile, tertile, optimal, thirds
  # median   : splits samples into 2 bins (below 50%, above 50%)
  # tertile  : splits samples into 3 bins (below 33%, 33%-67%, above 67%)
  # quartile : splits samples into 4 bins (below 25%, 25%-50%, 50%-75%, above 75%)
  # optimal  : splits samples into 2 bins (above & below optimum cutoff)
  # thirds   : splits samples into 3 bins (bottom 33%, middle33%, top 33% based on expression range)
  show_all_bins    = FALSE,          # TRUE = plot all bins (LOW, HIGH, MID/MED_HIGH+MED_LOW)
  multiple_cutoff  = FALSE,          # TRUE = compute cutoffs separately for substratify_var
  
  # Plot settings
  sig_score        = FALSE,          # TRUE = combine genes into one signature score
  conf_interval    = FALSE,          # TRUE = show confidence interval in survival curve
  plot_curve       = TRUE,           # TRUE = plot the survival curve
  plot_risk_table  = TRUE,           # TRUE = plot the risk table below the curve
  color_palette    = custom_palette, # vector of colors for groups c("#d73027","#0c2c84")
  
  # Survival data columns
  time_col         = "Time",         # metadata column containing Time values
  status_col       = "Status",       # metadata column containing Status values
  
  # Output
  prefix           = "",
  output_dir      = "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop"
)

# Calculate multi-gene signature scores
# Described in Levine et al https://doi.org/10.1186/gb-2006-7-10-r93
advanced_z <- function(gene_set, expr_matrix) {
  
  ix <- toupper(rownames(expr_matrix)) %in% toupper(gene_set)
  cat("Genes found:", sum(ix), "\n")
  
  # --------------------------------------------------------------------------------------------------------------------------------------------|
  # Method                 | What it does                                 | Goal                                  | Use Case
  # ---------------------- | ------------------------------------------------------------------------------------ | ----------------------------|
  # Row-Wise Centering     | Subtracts the median expression of a gene    | To look at how a gene's expression in | Differential Expression,    |
  # (MARGIN=1)             | across all samples from that gene's          | a specific sample deviates from its   | Clustering genes, Heatmaps  |
  #                        | expression in each sample.                   | typical expression across all samples.| (to compare gene behavior). |
  #
  # Column-Wise Centering  | Subtracts the median expression of all genes | To look at how a gene's expression    | Signature Scoring,          |
  # (MARGIN=2)             | in a sample from each gene's expression in   | average deviates from the overall     | Single-Sample Z-Score       |
  #                        | that sample.                                 | transcriptome of that specific sample.| (ssGSEA-like).              |
  # ---------------------------------------------------------------------------------------------------------------------------------------------
  
  # Median-center each sample across genes
  expr_matrix_centered <- base::sweep(x = expr_matrix, 
                                      MARGIN = 2, 
                                      STATS = apply(expr_matrix, 2, median, na.rm = TRUE), 
                                      FUN = "-")
  
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
calc_cutoffs <- function(cutoff_df, stratify_var, survival_params){
  
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
                                      variables = stratify_var)
      list(lower = res$cutpoint$cutpoint, upper = res$cutpoint$cutpoint, middle = NA)
    }, error = function(e) {
      list(lower = NA, upper = NA, middle = NA)
    }))
  
  # Categorize expression into bins
  model_col <- paste0("model_", stratify_var)
  cutoff_df <- cutoff_df %>%
    dplyr::filter(!is.na(.data[[stratify_var]])) %>%
    dplyr::mutate(!!model_col := dplyr::case_when(.data[[stratify_var]] > cutoffs$upper ~ "HIGH",
                                                  .data[[stratify_var]] <= cutoffs$lower ~ "LOW",
                                                  (!is.na(cutoffs$middle) & .data[[stratify_var]] > cutoffs$middle) ~ "MED_HIGH",
                                                  (!is.na(cutoffs$middle) & .data[[stratify_var]] <= cutoffs$middle) ~ "MED_LOW",
                                                  TRUE ~ "MID"),
                  !!model_col := factor(.data[[model_col]], 
                                        levels = c("LOW", "MED_LOW", "MID", "MED_HIGH", "HIGH"))) %>%
    dplyr::select(Sample_ID, all_of(model_col))
  
  # Optionally plot only HIGH and LOW
  if (!survival_params$show_all_bins) {
    cutoff_df <- cutoff_df %>%
      dplyr::filter(.data[[model_col]] %in% c("HIGH", "LOW"))
  }
  return(cutoff_df)
}

# Calculate cox model stats (HR, CI, p vals) for each facet
calc_cox_stats <- function(facet_df, stratify_var, surv_formula, survival_params){
  
  # ---- Cox model ----
  
  # Ensure model column is a factor
  model_col <- paste0("model_", stratify_var)
  facet_df[[model_col]] <- factor(facet_df[[model_col]])
  
  # Fit Cox model
  cox_model <- survival::coxph(formula = surv_formula, data = facet_df)
  
  # Cox model coefficients
  cox_coef_df <- summary(cox_model)$coefficients
  cox_ci_df <- as.data.frame(confint(cox_model))
  baseline <- levels(facet_df[[model_col]])[1]  # the factor baseline
  cox_df <- data.frame(Gene = stratify_var,
                       Target = gsub(paste0("^model_", stratify_var), "", rownames(cox_coef_df)), # remove "model" prefix if present
                       Reference = baseline,
                       HR = exp(cox_coef_df[, "coef"]),
                       CI_lower = exp(cox_ci_df[, 1]),
                       CI_upper = exp(cox_ci_df[, 2]),
                       pval = cox_coef_df[, "Pr(>|z|)"],
                       stringsAsFactors = FALSE) %>%
    dplyr::mutate(contrast = paste0(Target, " / ", Reference)) %>%
    dplyr::select(Gene, contrast, HR, pval, CI_lower, CI_upper, Target, Reference) %>%
    tibble::remove_rownames() 
  
  # ---- Optional: emmeans for pairwise contrasts ----
  # Step 1: estimated marginal means on log-hazard scale
  emm <- emmeans::emmeans(object = cox_model, specs = as.formula(paste0("~", model_col)))
  
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
                  Reference = sub(".* / ", "", contrast),
                  Gene = stratify_var)
  
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
                       stratify_var =  stratify_var,
                       facet_df = facet_df, survival_params = survival_params) %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("contrast") %>%
    dplyr::right_join(emmeans_df, by=c("contrast" = "contrast")) %>%
    dplyr::select(Gene, contrast, Target, Reference, HR, pval, CI_lower, CI_upper, everything())
  
  cox_df <- sapply(X = cox_df$contrast, FUN = calc_pvals,
                   stratify_var =  stratify_var,
                   facet_df = facet_df, survival_params = survival_params) %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("contrast") %>%
    dplyr::right_join(cox_df, by=c("contrast" = "contrast")) %>%
    dplyr::select(Gene, contrast, Target, Reference, HR, pval, CI_lower, CI_upper, everything())
  
  # ---- Return both ----
  list(cox_model_df = cox_df,
       emmeans_df = emmeans_df)
}

# Calculate all 7 non-parametric p-values for each contrast in cox_df or emmeans_df
calc_pvals <- function(contrast, facet_df, stratify_var, survival_params) {
  
  model_col <- paste0("model_", stratify_var)
  g1 <- sub(".* / ", "", contrast)
  g2 <- sub(" / .*", "", contrast)
  df_pair <- subset(facet_df, facet_df[[model_col]] %in% c(g1, g2))
  surv_obj <- survival::Surv(
    time   = df_pair[[survival_params$time_col]],
    event  = df_pair[[survival_params$status_col]],
    type   = "right",
    origin = 0
  )
  surv_form <- as.formula(paste("surv_obj ~", model_col))
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
plot_facets <- function(facet_df, stratify_var, surv_curve, cox_df, surv_type, facet_group, survival_params){
  
  if (survival_params$plot_curve) {
    
    model_col <- paste0("model_", stratify_var)
    
    # Legend labels
    # IMPORTANT: ggsurvplot() labels the groups in alphabetical order. So,
    # when we want to use custom labels, initialize them in alphabetical order.
    # Eg: c("High", "Low") instead of c("Low, "High")
    legend_label <- sort(unique(facet_df[[model_col]]))
    
    # Legend title
    if (surv_type %in% c("single_gene", "multi_gene", "signature")){
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
      title =  paste(na.omit(c(stratify_var, survival_params$substratify_var, 
                               facet_group, survival_params$cutoff_method)), 
                     collapse = "."),
      
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
    if (dplyr::n_distinct(facet_df[[model_col]]) == 2){
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
    }
    
    # NOTE: Using cowplot() and then ggsave() works nicely as compared to 
    # saving directly using ggsave()
    p <- cowplot::plot_grid(plotlist = surv_plot,
                            align = "hv",
                            axis = "tblr",
                            nrow = 2,
                            ncol = 1,
                            rel_widths = 1,
                            rel_heights = c(1, 0.45),
                            labels = NULL,
                            label_size = 14,
                            label_fontface = "bold")
    
    return(p)
  }
}

survival_analysis <- function(metadata, expr_data = NULL, survival_params) {
  
  # ---- Input checks & parameter extraction ----
  
  # Create output directory if missing
  if (!dir.exists(survival_params$output_dir)) {
    dir.create(survival_params$output_dir, recursive = TRUE)
  }
  
  rownames(expr_data) <- make.names(rownames(expr_data))
  colnames(expr_data) <- make.names(colnames(expr_data))
  metadata$Sample_ID <- make.names(metadata$Sample_ID)
  survival_params$stratify_var <- make.names(survival_params$stratify_var)
  
  # Extract parameters
  stratify_vars   <- survival_params$stratify_var
  substratify_var <- survival_params$substratify_var
  facet_var       <- survival_params$facet_var
  time_col        <- survival_params$time_col
  status_col      <- survival_params$status_col
  cutoff_method   <- survival_params$cutoff_method
  multiple_cutoff <- survival_params$multiple_cutoff
  show_all_bins   <- survival_params$show_all_bins
  sig_score       <- survival_params$sig_score
  
  # Must provide stratify_var
  if (is.null(stratify_vars) || length(stratify_vars) == 0) {
    stop("Must provide a non-empty stratify_var")
  }
  
  # substratify_var, if provided, must exist in metadata
  if (!is.null(substratify_var) && !substratify_var %in% colnames(metadata)) {
    stop("substratify_var not found in metadata columns")
  }
  
  # facet_var, if provided, must exist in metadata
  if (!is.null(facet_var) && !facet_var %in% colnames(metadata)) {
    stop("facet_var not found in metadata columns")
  }
  
  # Must define substratify_var if multiple_cutoff = TRUE
  if (isTRUE(multiple_cutoff) && is.null(substratify_var)) {
    stop("multiple_cutoff = TRUE requires substratify_var")
  }
  
  # Check stratify_var against both metadata and expression data and determine
  # type of survival analysis
  missing_genes <- setdiff(stratify_vars, rownames(expr_data))
  valid_genes   <- intersect(stratify_vars, rownames(expr_data))
  in_meta       <- all(stratify_vars %in% colnames(metadata))
  
  if (!in_meta) {
    # Expression-based stratification
    if (length(missing_genes) > 0) {
      warning("Some requested genes not found in expr_data: ", paste(missing_genes, collapse = ", "))
    }
    
    if (length(valid_genes) == 0) {
      stop("stratify_var must match either a column in metadata OR genes in expr_data.")
      
    } else if (length(valid_genes) > 1 && isTRUE(sig_score)) {
      message("Proceeding with signature expression-based survival analysis (", length(valid_genes), " valid genes).")
      surv_type <- "signature"
      
    } else if (length(valid_genes) > 1 && !isTRUE(sig_score)) {
      message("Proceeding with expression-based survival analysis for ", length(valid_genes), " gene(s).")
      surv_type <- "multi_gene"
      
    } else if (length(valid_genes) == 1) {
      message("Proceeding with expression-based survival analysis for single gene.")
      surv_type <- "single_gene"
    }
    
  } else {
    # Metadata-based stratification
    if (length(valid_genes) == 0) {
      message("Proceeding with metadata-based survival analysis.")
      surv_type <- "meta"
    } else {
      stop("stratify_var matches BOTH a column in metadata and gene(s) in expr_data — ambiguous.")
    }
  }
  
  # Define facet groups
  if (!is.null(facet_var) && facet_var %in% colnames(surv_df)) {
    facet_groups  <- unique(surv_df[[facet_var]])
  } else {
    facet_groups  <- NA_character_ # placeholder for whole dataset
  }
  
  # ---- Format metadata ----
  
  metadata <- metadata %>% 
    dplyr::mutate(Sample_ID = make.names(names = Sample_ID), 
                  !!time_col := as.numeric(.data[[time_col]])) %>%
    dplyr::filter(.data[[time_col]] > 0 & !is.na(.data[[time_col]])) %>%
    dplyr::distinct(Sample_ID, .keep_all = TRUE)
  
  # ---- Prepare expression data ----
  
  if (surv_type == "signature") {
    # Signature score survival (whole dataset, even if substratify_var defined)
    sig_scores <- advanced_z(gene_set = stratify_vars, expr_matrix = expr_data)
    expr_df    <- as.data.frame(sig_scores, check.names = FALSE) %>%
      dplyr::rename(sig_score = identity(1)) %>%
      tibble::rownames_to_column(var = "Sample_ID")
    
    # Update local variable and survival_params
    stratify_vars <- "sig_score"
    survival_params$stratify_var <- stratify_vars
    
  } else if (surv_type %in% c("single_gene", "multi_gene")) {
    # Single/Multi-gene survival (expression-based)
    expr_df <- expr_data[stratify_vars, , drop = FALSE] %>%
      t() %>%
      as.data.frame(check.names = FALSE) %>%
      tibble::rownames_to_column(var = "Sample_ID")
    
  } else if (surv_type == "meta") {
    # Metadata-based survival
    expr_df <- metadata %>%
      dplyr::select(Sample_ID, dplyr::all_of(stratify_vars))
    
  } else {
    stop("Invalid stratify_var: must be gene(s) in expr_data or a column in metadata.")
  }
  
  colnames(expr_df) <- make.names(colnames(expr_df))
  
  # ---- Merge expression data with metadata ----
  
  keep_cols <- unique(c(time_col, status_col, stratify_vars, substratify_var, facet_var))
  
  surv_df <- expr_df %>%
    dplyr::inner_join(metadata, by = c("Sample_ID"="Sample_ID")) %>%
    dplyr::select(Sample_ID, dplyr::all_of(keep_cols))
  
  if (nrow(surv_df) == 0) stop("No overlapping Sample_IDs between expr_data and metadata.")
  
  
  # ---- Define model column for metadata based survival ----
  
  if (surv_type == "meta") {
    model_col <- paste0("model_", stratify_vars)
    surv_df <- surv_df %>%
      dplyr::mutate(!!model_col := paste(.data[[stratify_vars]],
                                         if (!is.null(substratify_var)) .data[[substratify_var]] else NULL,
                                         if (!is.null(facet_var)) .data[[facet_var]] else NULL,
                                         sep = "_"))
  }
  
  # ---- Define model column for expression based survival ----
  
  if (surv_type %in% c("signature", "single_gene", "multi_gene")) {
    
    full_merged_df <- data.frame(Sample_ID = surv_df$Sample_ID)
    # Loop over each gene
    for (stratify_var in stratify_vars){
      
      model_col <- paste0("model_", stratify_var)
      merged_df <- tibble::tibble()
      # Loop over each facet group
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
          df <- calc_cutoffs(cutoff_df = cutoff_df,
                             stratify_var = stratify_var,
                             survival_params = survival_params)
          
          merged_df <- dplyr::bind_rows(merged_df, df)
        }
      }
      
      # Append substratify_var to model if defined
      if (!is.null(substratify_var)) {
        merged_df <- merged_df %>%
          dplyr::mutate(!!model_col := paste0(model_col, "_", .data[[substratify_var]]))
      }
      
      # Add model columns for every gene
      full_merged_df <- dplyr::left_join(x = full_merged_df, 
                                         y = merged_df,
                                         by = c("Sample_ID" = "Sample_ID"))
    }
    
    # Merge classifications into surv_df
    surv_df <- surv_df %>%
      dplyr::left_join(full_merged_df, by=c("Sample_ID"="Sample_ID"))
  }
  
  
  # ---- Main loop: facets, compute stats, save, plot ----
  
  # Initialize empty data frames to store all stats
  all_cox_df     <- tibble::tibble()
  all_emmeans_df <- tibble::tibble()
  
  # Open PDF for all KM plots
  pdf(file = file.path(survival_params$output_dir, "KM_curves.pdf"),
      width = 7, height = 7)  # open PDF device
  
  # Loop over each gene
  for (stratify_var in stratify_vars){
    
    model_col <- paste0("model_", stratify_var)
    
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
      if (dplyr::n_distinct(facet_df[[model_col]]) < 2) {
        message("Skipping facet (single group): ", facet_group)
        next
      } else if (base::setequal(x = unique(facet_df[[model_col]]), y = c("HIGH", "LOW"))) {
        survival_params$color_palette <- c(HIGH = "#d73027", LOW = "#0c2c84")
      }
      
      # Create a survival object (Alive = 0, Dead = 1)
      surv_object <- survival::Surv(time   = facet_df[[time_col]],
                                    event  = facet_df[[status_col]],
                                    type   = "right",
                                    origin = 0)
      
      # Create a formula for survival analysis
      surv_formula <-  as.formula(paste("surv_object ~", model_col))
      
      # Create a fit for Kaplan-Meier curve
      # NOTE: survival::survfit() gives error in ggsurvplot()
      surv_curve <- survminer::surv_fit(formula   = surv_formula,
                                        data      = facet_df)
      
      # Compute Cox & emmeans stats for each facet
      stats_list <- calc_cox_stats(facet_df, stratify_var, surv_formula, survival_params)
      
      # Add identifiers for traceability
      stats_list$cox_model_df$Facet <- facet_group
      stats_list$emmeans_df$Facet   <- facet_group
      
      # Accumulate stats
      all_cox_df     <- dplyr::bind_rows(all_cox_df, stats_list$cox_model_df)
      all_emmeans_df <- dplyr::bind_rows(all_emmeans_df, stats_list$emmeans_df)
      
      # Plot survival curve for each facet
      cox_df <- stats_list$cox_model_df
      p <- plot_facets(facet_df, stratify_var, surv_curve, cox_df, surv_type, facet_group, survival_params)
      
      # Print the plot to the current PDF page
      print(p)
    }
  }
  
  dev.off()  # close PDF device
  
  # Save combined stats to a single Excel workbook
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "cox_stats")
  openxlsx::writeData(wb, sheet = "cox_stats", x = all_cox_df)
  openxlsx::addWorksheet(wb, sheetName = "emmeans_stats")
  openxlsx::writeData(wb, sheet = "emmeans_stats", x = all_emmeans_df)
  openxlsx::addWorksheet(wb, sheetName = "surv_df")
  openxlsx::writeData(wb, sheet = "surv_df", x = surv_df)
  openxlsx::saveWorkbook(wb, file = file.path(survival_params$output_dir, "Survival_Stats.xlsx"),
                         overwrite = TRUE)
}

# Run to understand how to define parameters for the survival function
show_survival_scenarios()

# Run your survival analysis
#survival_analysis(metadata, expr_data, survival_params)

#******************************************************************************#
#                       SURVIVAL CURVE RELATED FUNCTIONS                       #
#******************************************************************************#

# Read this paper for survival analysis
# https://doi.org/10.1093/jncimonographs/lgu024

# NOTE:  Output of prep_expr_df is df
#log_norm_counts is matrix with SYMBOLS as rownames
prep_expr_df <- function(log_norm_counts, metadata, plot_genes, survival_params){
  
  # Merge expression data with survival data
  if (survival_params$gene_sig_score == TRUE){
    
    # Calculate gene signature score
    expr_df <- as.data.frame(advanced_Z(plot_genes, log_norm_counts))
    
    expr_df <- expr_df %>%
      data.frame() %>%
      dplyr::rename(combined.exp = identity(1)) %>%
      tibble::rownames_to_column("Sample_ID") %>%
      dplyr::inner_join(metadata, by=c("Sample_ID"="Sample_ID")) %>%
      dplyr::select(Sample_ID, combined.exp, Time, Status)
  } else {
    expr_df <- log_norm_counts %>%
      t() %>%
      data.frame() %>%
      tibble::rownames_to_column("Sample_ID") %>%
      dplyr::inner_join(metadata, by=c("Sample_ID"="Sample_ID")) %>%
      dplyr::select(Sample_ID, all_of(plot_genes), Time, Status)
  }
  
  # Add split_by column to expr_df to define groups in order to calculate multiple_cutoff
  if (!is.na(survival_params$split_by)){
    expr_df <- expr_df %>% 
      dplyr::left_join(metadata %>% dplyr::select(Sample_ID, survival_params$split_by),
                       by=c("Sample_ID"="Sample_ID"))
  }
  
  return(expr_df)
}

# NOTE:  Output of calc_cutoffs is list(df,ls)
# If plotting by Sex, make sure to create column "model" based on which lines will be split
# calc_cutoffs <- function(df, gene, group, survival_params){
#   
#   # Identify upper & lower cutoffs based on stratify_criteria
#   #*************************Split samples by median**************************#
#   if(survival_params$stratify_criteria == "m"){
#     quartiles <- stats::quantile(x = df[[gene]],
#                                  probs = c(0, 0.25, 0.50, 0.75, 1),
#                                  na.rm=TRUE)
#     
#     cutoff_lower_end <- quartiles[[3]]
#     cutoff_upper_end <- quartiles[[3]]
#     cutoff_middle <- "NA"
#   }
#   
#   #****************Split samples into top and bottom tertiles****************#
#   else if(survival_params$stratify_criteria == "t"){
#     tertiles <- stats::quantile(x = df[[gene]],
#                                 probs = c(0, 0.33, 0.66, 1),
#                                 na.rm=TRUE)
#     
#     cutoff_lower_end <- tertiles[[2]]
#     cutoff_upper_end <- tertiles[[3]]
#     cutoff_middle <- "NA"
#   }
#   
#   #***************Split samples into top and bottom quartiles****************#
#   else if(survival_params$stratify_criteria == "q"){
#     quartiles <- stats::quantile(x = df[[gene]], 
#                                  probs = c(0, 0.25, 0.50, 0.75, 1),
#                                  na.rm=TRUE)
#     
#     cutoff_lower_end <- quartiles[[2]]
#     cutoff_upper_end <- quartiles[[4]]
#     cutoff_middle <- quartiles[[3]]
#   }
#   
#   #*********************Split expression range by thirds*********************#
#   else if(survival_params$stratify_criteria == "th"){
#     quartiles <- stats::quantile(x = df[[gene]], 
#                                  probs = c(0, 0.25, 0.50, 0.75, 1),
#                                  na.rm=TRUE)
#     iqr <- stats::IQR(x = df[[gene]],
#                       na.rm=TRUE)
#     
#     # Normal range of expression values lie between cutoff_lower & cutoff_upper
#     cutoff_upper <- quartiles[[4]]+1.5*iqr
#     cutoff_lower <- dplyr::if_else(quartiles[[1]]-1.5*iqr > 0, quartiles[[1]]-1.5*iqr, 0)
#     
#     # Based on normal range of expression, identify onethird & twothird cutoff
#     cutoff_lower_end <- cutoff_lower + (cutoff_upper-cutoff_lower)/3
#     cutoff_upper_end <- cutoff_lower + (cutoff_upper-cutoff_lower)*2/3
#     cutoff_middle <- "NA"
#   }
#   
#   #***************Split expression range using optimal cutoff****************#
#   else if(survival_params$stratify_criteria == "o"){
#     
#     # Sometimes quartiles will look like: 
#     # 0%       25%      50%      75%     100% 
#     # 0.000000 0.000000 0.000000 0.000000 3.495493 
#     # In such cases, surv_cutpoint() will fail. So, we add extra if() here.
#     quartiles <- stats::quantile(x = df[[gene]], 
#                                  probs = c(0, 0.25, 0.50, 0.75, 1),
#                                  na.rm=TRUE)
#     
#     if (quartiles[[4]] > quartiles[[2]]){
#       res.cut <- survminer::surv_cutpoint(data = df,
#                                           time = "Time",
#                                           event = "Status",
#                                           variables = gene)
#       
#       cutoff_lower_end <- res.cut$cutpoint$cutpoint
#       cutoff_upper_end <- res.cut$cutpoint$cutpoint
#       cutoff_middle <- "NA"
#     } else{
#       #cat("Surv cutpoint unable to detect optimum cutoff")
#       cutoff_lower_end <- "NA"
#       cutoff_upper_end <- "NA"
#       cutoff_middle <- "NA"
#     }
#   }
#   
#   # Categorize the sample based on above cutoffs
#   if (survival_params$plot_all_bins == TRUE & survival_params$stratify_criteria == "q"){
#     df <- df %>% 
#       dplyr::mutate(Expression = dplyr::case_when(get(gene) > cutoff_upper_end ~ "HIGH",
#                                                   get(gene) <= cutoff_lower_end ~ "LOW",
#                                                   get(gene) <= cutoff_middle ~ "MED_LOW",
#                                                   TRUE ~ "MED_HIGH"))
#     
#   } else if (survival_params$plot_all_bins == TRUE) {
#     df <- df %>% 
#       dplyr::mutate(Expression = dplyr::case_when(get(gene) > cutoff_upper_end ~ "HIGH",
#                                                   get(gene) <= cutoff_lower_end ~ "LOW",
#                                                   TRUE ~ "MID"))
#     
#   } else if (survival_params$stratify_criteria == "none") {
#     #When plotting by Sex, Treatment response, we dont use expression data.
#     df <- df %>% 
#       dplyr::mutate(Expression = model)
#     cutoff_lower_end <- NA
#     cutoff_upper_end <- NA
#     cutoff_middle <- NA
#     
#   } else {
#     df <- df %>% 
#       dplyr::mutate(Expression = dplyr::case_when(get(gene) > cutoff_upper_end ~ "HIGH",
#                                                   get(gene) <= cutoff_lower_end ~ "LOW",
#                                                   TRUE ~ "MID")) %>%
#       dplyr::filter(Expression != "MID")
#   }
#   
#   # # Print the cutoffs
#   # cat("\nGene         :", gene)
#   # cat("\nGroup        :", group)
#   # cat("\nLower cutoff :", cutoff_lower_end)
#   # cat("\nUpper cutoff :", cutoff_upper_end)
#   # cat("\nMiddle cutoff:", cutoff_middle)
#   
#   # Create a list to store cutoff values
#   ls <- list("group" = c(), 
#              "gene" = c(), 
#              "lower" = c(), 
#              "upper" = c(), 
#              "middle" = c())
#   
#   ls$group <- c(group)
#   ls$gene <- c(gene)
#   ls$lower <- c(cutoff_lower_end)
#   ls$upper <- c(cutoff_upper_end)
#   ls$middle <- c(cutoff_middle)
#   
#   # Return the df and the cutoffs
#   return(list(df, ls))
# }

# NOTE:  Output of calc_surv_stats is list. 
# It also generate survival plot with risk table
calc_surv_stats <- function(df, group, prefix, output_dir){
  
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
                      path = output_dir,
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
plot_survival <- function(expr_df, gene, survival_params, prefix, output_dir){
  
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
    cox_stats <- calc_surv_stats(df, group, prefix, output_dir)
    
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

# ---- CRISPR FUNCTIONS ----

prep_crispr_guides <- function(data, output_dir){
  
  # Ensure required columns are present
  required_cols <- c("ID", "Seq", "Gene")
  
  if (!all(required_cols %in% colnames(data))) {
    stop("Input data must contain columns: ID, Seq, Gene")
  }
  
  # Define a simple reverse complement function
  reverse_complement <- function(seq) {
    seq_chars <- unlist(strsplit(seq, split = ""))
    comp <- c(A = "T", T = "A", G = "C", C = "G")
    rc <- rev(comp[seq_chars])
    paste(rc, collapse = "")
  }
  
  # Compute reverse complements
  data <- data %>%
    dplyr::mutate(Rev_Seq = vapply(Seq, reverse_complement, character(1)),
                  Mageck_format = paste(ID, Seq, Gene, sep = ","),
                  Mageck_format_rev = paste(ID, Rev_Seq, Gene, sep = ",")
    )
  
  # Write to Excel
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "Guides")
  openxlsx::writeData(wb, sheet = "Guides", x = data, rowNames = FALSE)
  openxlsx::saveWorkbook(wb, 
                         file = file.path(output_dir, "CRISPR_Guides.xlsx"), 
                         overwrite = TRUE)
  
  message("CRISPR guide data saved to ", file.path(output_dir, "CRISPR_Guides.xlsx"))
}

# Dataframe as input
# Column "Gene" MUST be present
# Column Gene MUST have control sgRNAs labelled as "none" and/or "safe"
calc_t_score <- function(data){
  
  # Create a dataframe of control sgRNAs
  data_control <- data %>%
    dplyr::filter(Gene %in% c("none", "safe"))
  
  median_ctrl <- median(data_control$LFC, na.rm=TRUE)
  sd_ctrl <- sd(data_control$LFC, na.rm=TRUE)
  
  # Normalize to control sgRNAs
  data <- data %>%
    dplyr::mutate(pZ = (LFC-median_ctrl)/sd_ctrl)
  
  data_control <- data_control %>%
    dplyr::mutate(pZ = (LFC-median_ctrl)/sd_ctrl)
  
  U_ctrl <- median(data_control$pZ)
  Var_ctrl <- var(data_control$pZ)
  N_ctrl <- mean((data %>% dplyr::count(Gene))$n)
  # Nctrl is the average number of sgRNAs per gene in a given screen
  
  data <- data %>%
    dplyr::group_by(Gene) %>%
    dplyr::mutate(U_gene = median(pZ),
                  Var_gene = var(pZ),
                  N_gene = n(),
                  U_ctrl = U_ctrl,
                  Var_ctrl = Var_ctrl,
                  N_ctrl = N_ctrl,
                  S_gene = (Var_gene*(N_gene-1)) + (Var_ctrl*(N_ctrl-1)),
                  t_score = (U_gene - U_ctrl)/sqrt(S_gene/N_gene + S_gene/N_ctrl),
                  Abs_t_score = abs(t_score)) %>%
    dplyr::select(Gene, U_gene, Var_gene, N_gene, U_ctrl, Var_ctrl, N_ctrl, S_gene, t_score, Abs_t_score) %>%
    dplyr::distinct_at("Gene", .keep_all = TRUE)
  
  return(data)
}

# data is a dataframe output of calc_t_score
plot_t_score <- function(data, disp_genes, suffix, save_path){
  
  y_cutoff <- sort(data$Abs_t_score, decreasing = TRUE)[100]
  xmin <- floor(min(data$U_gene))
  xmax <- ceiling(max(data$U_gene))
  ymin <- 0
  ymax <- max(data$Abs_t_score)
  
  color_breaks <- c(-20,0,20)
  p <- ggplot2::ggplot(data = data,
                       aes(x = U_gene, 
                           y = Abs_t_score,
                           size = Abs_t_score,
                           #color = pz,
                           fill = U_gene)) +
    # Plot dot plot
    ggplot2::geom_point(col="black", 
                        shape=21,
                        stroke=0.5,
                        position=position_jitter(h=0.01,w=0.01)) +
    # Define the theme of plot
    ggplot2::theme_classic() +
    ggplot2::labs(fill = "U_gene") +
    coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax), clip = "off") +
    #scale_x_continuous(breaks = seq(-5, 5, by = 1)) +
    #scale_y_continuous(breaks = seq(0, 5, by = 1)) +
    ggplot2::guides(size = "none",
                    fill = guide_colourbar(theme = theme(legend.key.width  = unit(0.75, "lines"),
                                                         legend.key.height = unit(10, "lines"),
                                                         legend.ticks = element_blank(),
                                                         legend.frame = element_rect(colour = "Black",
                                                                                     linewidth = 0.5)))) +
    # Define the color of the dots
    ggplot2::scale_fill_viridis_c(option="turbo", limits =c(-5,3))
  
  if (length(disp_genes) > 0){
    p <- p + ggrepel::geom_text_repel(data = data %>% dplyr::filter(Gene %in% disp_genes),
                                      mapping = aes(label = Gene),
                                      size = 5,
                                      show.legend = FALSE,
                                      direction = "both",   #"y"
                                      box.padding = 2.5,      # increases line length somehow
                                      point.padding = 0.1,  # distance around point = dist between line and point
                                      max.overlaps = nrow(data),
                                      position = position_quasirandom())
    
    
  }
  #geom_hline(yintercept= y_cutoff, linetype ="dotted")
  
  # scale_fill_gradientn(colors=c("#007ba7", "Black","#FFFF00"), 
  #                      limits=c(-20, 20), 
  #                      values=c(0, scales::rescale(color_breaks, from = range(color_breaks)), 1))
  #scale_fill_gradient2(low="#007ba7", mid="Black", high="Yellow", midpoint = 0, limits=c(-5, 2))
  #scale_fill_continuous_diverging(palette = "Tofino")
  
  ggplot2::ggsave(filename = paste0(suffix, ".jpg"),
                  plot = p,
                  device = "tiff",
                  path = save_path,
                  width = 7,
                  height = 7,
                  units = c("in"),
                  dpi = 300,
                  limitsize = TRUE,
                  bg = NULL)
}

# ---- 🔵🟡 VENN DIAGRAM ----

# ALTERNATIVE: Use http://www.ehbio.com/test/venn/#/  
plot_venn <- function(data, filename, output_dir){
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(metadata = data, output_dir = output_dir, filename = filename)
  
  # Validate column count (upto 4 columns)
  if (ncol(data) < 1 || ncol(data) > 4) {
    stop("`data` must have between 1 and 4 columns.")
  }
  
  # ---- 🛠️ Configure Venn Diagram Settings ----
  
  # Clean column names (replace non-standard characters with space for visualization)
  colnames(data) <- stringr::str_replace_all(colnames(data), c("_" = " ", "\\." = " "))
  
  # Set category position (cat.pos), label distance (cat.dist), font size (cex), and palette
  # Settings are optimized based on the number of sets (columns)
  if (ncol(data) == 4){
    pos <- c(330, 15, 330, 15) # Quadrant positioning
    dist <- c(0.27, 0.25, 0.15, 0.13)
    cex = 2
    palette1 <- c("#C8E7F5", "#00008C", "#F6D2E0", "#E75480")
  } else if (ncol(data) == 3){
    pos <- c(0, 0, 180)
    dist <- c(0.1, 0.1, 0.1)
    cex = 2
    palette1 <- c("#C8E7F5", "#F6D2E0", "#db6d00")
  } else if (ncol(data) == 2){
    pos <- c(0, 0)
    dist <- c(0.05, 0.05)
    cex = 2.75
    palette1 <- c("#C8E7F5", "#db6d00")
  } else if (ncol(data) == 1){
    pos <- c(0)
    dist <- c(0.1)
    cex = 2.75
    palette1 <- c("#F6D2E0")
  }
  
  # Create a dataframe to store the wrapped column names for clean labels
  annotation <- data.frame(Labels = stringr::str_wrap(colnames(data), width = 10))
  
  # Convert the data frame to a named list (VennDiagram input format)
  genes <- base::vector(mode = "list", length = ncol(data))
  names(genes) <- annotation$Labels
  
  for (i in 1:ncol(data)){
    # Remove NA values to generate a clean list of genes for each set (label)
    genes[[i]] <- data[!is.na(data[i]), i]
  }
  
  # ---- 📊 Generate Venn Diagram ----
  
  file_name <- file.path(output_dir, paste0("Venn_Diagram_", filename, ".tiff"))
  VennDiagram::venn.diagram(x = genes,
                            main = filename,
                            category.names = annotation$Labels,
                            filename = file_name,
                            output = TRUE,
                            scaled = FALSE, # Disable proportional scaling (can distort appearance)
                            imagetype = "tiff",
                            height = 11,
                            width = 11,
                            units = "in",
                            resolution = 600,
                            compression = "lzw",
                            margin = 0.3, # Amount of white space around Venn Diagram in grid units
                            
                            # 1️⃣ Line Formatting
                            lwd = 1.5,                 # line thickness
                            lty = 1,                   # line type
                            col = "black",             # line color
                            
                            # 2️⃣ Number Formatting
                            cex = cex,                 # font size (2 or 2.75)
                            fontface = "bold",         # font style  
                            fontfamily = "sans",       # font type
                            
                            # 3️⃣ Main Title Formatting
                            main.cex = 2,              # font size
                            main.fontface = "bold",    # font style 
                            main.fontfamily = "sans",  # font type
                            main.col = "black",        # font color
                            
                            # 4️⃣ Category Label Formatting
                            cat.cex = 2,               # font size
                            cat.fontface = "bold",     # font style
                            cat.fontfamily = "sans",   # font type
                            cat.col = palette1, #"black"
                            cat.pos = pos,
                            cat.dist = dist,
                            
                            # 5️⃣ Fill Colors
                            fill = palette1,
                            alpha = rep(0.5, ncol(data)),    # 50% transparency for fill color
                            ext.text = TRUE,           # Draw external text (labels)
                            #cat.default.pos = "outer",
                            disable.logging = TRUE)
  
  # ---- 💾 Save Overlapping Genes (Excel Output) ----
  
  # 1️⃣ Detect Unique Overlapping Genes Across All Combinations
  
  overlap_list <- list()
  detected_genes <- character(0) # Use character(0) for an empty character vector
  
  # Iterate combinations from the largest (all datasets) down to the smallest (single datasets)
  for (n in seq(from = ncol(data), to = 1, by = -1)){
    
    # Generate all possible n-element combinations of set names
    cmb <- utils::combn(x = names(genes), m = n)
    
    for (col in 1:ncol(cmb)){
      
      # Get the names of the datasets in the current combination
      datasets <- cmb[,col]
      
      # Calculate the intersection (overlapping genes) for the current combination
      overlap <- Reduce(intersect, genes[datasets])
      
      # Remove genes already detected in overlaps of previous (larger) comparisons.
      # This ensures that each gene is counted only in the *largest* combination
      # where it is present (i.e., unique to that specific overlap region).
      overlap_unique <- base::setdiff(overlap, detected_genes)
      
      # Add the newly found unique genes to the master list of detected genes
      detected_genes <- c(detected_genes, overlap_unique)
      
      # Collapse the set names into a clean label (e.g., "Set1.Set2.Set3")
      names <- paste(datasets, collapse = ".")
      
      # Store the vector of unique overlapping genes under the meaningful name
      overlap_list[[names]] <- overlap_unique
    }
  }
  
  # 2️⃣ Format Results into a Data Frame for Saving
  
  # Identify maximum number of genes present in any single unique overlap region
  max_len = max(lengths(overlap_list))
  
  # Create an empty data frame structure to hold all results
  results = data.frame(matrix("", nrow = max_len, ncol = length(overlap_list)))
  
  # Set row names to be generic gene numbers
  rownames(results) <- paste0("Gene#", seq(max_len))
  
  # Set column names to the meaningful set combination labels
  colnames(results) <- names(overlap_list)
  
  # Populate the dataframe (This loop was missing in your original paste)
  for (i in 1:length(overlap_list)){
    if (length(overlap_list[[i]]) > 0){
      results[1:length(overlap_list[[i]]), i] <- overlap_list[[i]]
    }
  }
  
  # 3️⃣ Save Results to Excel Workbook (using openxlsx)
  
  file_name <- file.path(output_dir, paste0("Overlap_", filename, ".xlsx"))
  wb <- openxlsx::createWorkbook()
  
  openxlsx::addWorksheet(wb, sheetName = "Output")
  openxlsx::writeData(wb, sheet = "Output", x = results, rowNames = TRUE, keepNA = FALSE)
  
  openxlsx::addWorksheet(wb, sheetName = "Input")
  openxlsx::writeData(wb, sheet = "Input", x = data, rowNames = FALSE)
  
  openxlsx::saveWorkbook(wb, file = file_name, overwrite = TRUE)
  
  # ---- 🪵 Log Output ----
  
  log_info(sample = "", 
           step   = "plot_venn",
           msg    = glue::glue("Venn plot saved successfully to : '{file_name}'."))
  
  return(invisible(NULL))
}

# ---- 🧩 📊 UPSET PLOT ----

# 💡 Key Features of UpSet Plot Visualization: 
# 1️⃣ Sets/List Names (e.g., Genes/Proteins) are displayed on the Y-axis of the bottom matrix.
# 2️⃣ Intersection size (count of overlapping elements) is based on the Y-values of the top bar graph.
# 3️⃣ All possible intersections (combinations) are displayed along the X-axis of the bottom matrix.

plot_upset <- function(listInput = NULL, selected_sets = NULL,
                       min_intersection_size = NULL, filename = "Upset_Plot",
                       output_dir) {
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(filename = filename, n_pcs = min_intersection_size,
                  output_dir = output_dir)
  
  # Check required inputs after defaulting
  if (!is.list(listInput) || length(listInput) == 0) {
    log_warn(sample = "",
             step = "plot_upset", 
             msg ="Input list is empty or invalid. using example dataset")
  }
  
  # Handle Default Input (Example Data)
  if (is.null(listInput)) {
    message("No listInput provided. Using default example kinase data for demonstration.")
    listInput <- list(
      DMPK = c(1,2,3,4,5,6,7,8,9,11,16,17),
      MAPK1 = c(1,2,3,4,5,6,7,8,9,16,17,18),
      RAF1 = c(1,2,3,4,5,6,7,8,9,16,17,18),
      ROCK1 = c(1,2,3,4,5,6,7,8,9,17,18),
      PIM1 = c(1,2,3,4,5,6,7,8,9,17),
      DYRK2 = c(1,2,3,4,5,6,7,8,9,17),
      STK26 = c(1,2,3,4,5,6,7,8,9,17),
      MAP2K2 = c(1,2,3,4,5,6,7,8,9),
      LIMK1 = c(1,2,3,4,5,6,7,8,9),
      MYO3A = c(1,2,3,4,5,6,7,8,9),
      TSSK1B = c(1,2,3,4,5,6,7,8,9),
      HIPK4 = c(1,2,3,4,5,6,7,8,9),
      DCLK1 = c(1,2,3,4,5,6,7,8,9),
      NEK3 = c(1,2,3,4,5,6,7,8,9),
      CDKL1 = c(1,2,3,4,5,6,7,8,9),
      TRPM7 = c(1,2,3,4,5,6,7,8,10,11,12,13,14,15),
      MAPK7 = c(1,2,3,4,5,6,7,8,17),
      IKBKE = c(1,2,4,5,6,7,8,9,17),
      PRKAG3 = c(1,2,4,5,6,7,8,9),
      TAF1L = c(1,2,4,5),
      JAK2 = c(2,4,5,6,7,8,9,16,17,18)
    )
    
    # Set parameters relevant to the default data
    if (is.null(selected_sets)) {
      selected_sets <- c("DCLK1", "MAPK1", "CDKL1", "ROCK1", "MAPK7",
                         "HIPK4", "MAP2K2", "DYRK2")
    }
    if (is.null(min_intersection_size)) {
      min_intersection_size <- 5
    }
  }
  
  # ---- 🛠 Prepare Data for UpSetR ----
  
  # Convert the named R list of vectors into the binary matrix format required by UpSetR
  upset_data <- UpSetR::fromList(listInput)
  
  # ---- 📊 Create UpSet Plot ----
  
  # Generate the UpSet plot object
  p <- UpSetR::upset(data = upset_data,
                     empty.intersections = "on",
                     cutoff = min_intersection_size,      # minimum intersection size to show
                     mb.ratio = c(0.5, 0.5),              # ratio between main bar and sets bar
                     sets = selected_sets,                # sets to display (order matters)
                     order.by = "freq",                   # order intersections by frequency
                     main.bar.color = "#1f78b4",          # nicer color for bars
                     sets.bar.color = "#33a02c",
                     #nintersects = 5,                    # number of groups on X axis
                     #nsets = 21,                         # number of groups on Y axis
                     text.scale = c(2, 2, 1.5, 1.5, 2, 1.5) # scale axis and text sizes
  )
  
  # Convert to ggplot object for standardized saving
  ggplot_obj <- ggplotify::as.ggplot(p)
  
  # ---- 💾 Save Plot ----
  
  file_name <- file.path(output_dir, paste0("Upset_Plot_", filename, ".pdf"))
  ggplot2::ggsave(filename = file_name,
                  plot = ggplot_obj,
                  height = 11,
                  width = 11)
  
  # ---- 🪵 Log Output ----
  
  log_info(sample = "", 
           step   = "plot_upset",
           msg    = glue::glue("Upset plot saved successfully to : '{file_name}'."))
  
  return(invisible(NULL))
}

# ---- 🧭 PCA & UMAP PLOT ----

# 1️⃣ Data Orientation:
#    - Rows MUST be Observations (i.e. Samples)
#    - Columns MUST be Variables (i.e. Features, Genes, DMRs, etc.)
#    - If your data matrix is Features × Samples, you MUST use t() to transpose
#      data matrix before running prcomp().

# 2️⃣ Centering and Scaling:
#    - prcomp() centers columns by default (center = TRUE)
#      -> Subtracts the mean of each feature across all samples.
#    - Scaling (scale. = TRUE) divides each column by its standard deviation.
#      -> This ensures all features contribute equally, regardless of their 
#         magnitude or raw variance.
#    - Note: prcomp(X, center=T, scale.=T) is equivalent to 
#            prcomp(scale(X), center=F, scale.=F).

# 3️⃣ When to Scale:
#    - Use scale. = TRUE when feature variance is NOT biologically meaningful
#      Eg: raw counts, FPKM/TPM
#    - Use scale. = FALSE when feature variance is biologically meaningful
#      Eg: VST, rlog, fractions, or data that has been variance-stabilized).

# 4️⃣ PCA Interpretation:
#    - PCA computes principal components for Observations (i.e. Samples) in the
#      space of Variables (i.e. features)
#    - The first PC (PC1) captures the largest variance among samples.
#    - '$x' matrix contains the new coordinates for your Observations (i.e. Samples).#    
#    - '$rotation' matrix shows the Loadings, indicating how each original 
#      Feature/Gene contributes to the direction of the Principal Components. 

# 5️⃣ Example Usage (Assuming gene_data is Genes × Samples of VST counts):
#    pca_result <- prcomp(t(gene_data), center = TRUE, scale. = FALSE)
#    # t(gene_data) ensures Samples are rows and Genes are columns.

# IMPORTANT: plot_pca() NEEDS expr_mat to be VST/rlog/logit transformed. It 
# hasn't been coded to work with other types of data. metadata MUST have column
# "Sample_ID"

plot_pca <- function(expr_mat, txi, metadata, filename, output_dir,
                     perform_vst = TRUE, top_n_genes = 500, skip_plot = FALSE){
  # For ggrepel
  set.seed(1234)
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(expr_mat = expr_mat, metadata = metadata, 
                  logic_vars = list(perform_vst = perform_vst,
                                    skip_plot   = skip_plot),
                  filename = filename, output_dir = output_dir)
  
  if (!"Sample_ID" %in% colnames(metadata)) {
    log_error(sample = "",
              step   = "plot_pca",
              msg    = glue::glue("`metadata` must contain 'Sample_ID' column."))
  }
  
  log_warn(sample = "",
           step   = "plot_pca",
           msg    = "'expr_mat' must be genes(rows) x samples(columns) matrix.")
  
  # ---- VST Transformation ----
  
  if (perform_vst){
    
    log_warn(sample = "",
             step   = "plot_pca",
             msg    = glue::glue("Applying VST transformation. VST works on raw counts ONLY.
                                 Please make sure provided expr_mat is raw_counts."))
    
    # Reformat raw_counts_mat and metadata for DESeq2
    deseq2_data <- prepare_deseq2_input(metadata = metadata,
                                        txi      = txi,
                                        expr_mat = expr_mat,
                                        design   = 1)
    # Prepare DESeq2 object
    if (!is.null(deseq2_data$expr_mat)){
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = deseq2_data$expr_mat,
                                            colData   = deseq2_data$metadata,
                                            design    = ~1)
    } else if (!is.null(deseq2_data$txi)){
      dds <- DESeq2::DESeqDataSetFromTximport(txi     = deseq2_data$txi,
                                              colData = deseq2_data$metadata,
                                              design  = ~1)
    }

    dds <- DESeq2::estimateSizeFactors(dds)
    vsd <- DESeq2::vst(dds, blind = TRUE)
    expr_mat <- SummarizedExperiment::assay(vsd)
    
  } else {
    
    log_info(sample = "",
             step   = "plot_pca",
             msg    = glue::glue("Using provided matrix directly (skipping VST)."))
    expr_mat <- as.matrix(expr_mat) 
  }
  
  # ---- Feature Selection ----
  
  if (nrow(expr_mat) > top_n_genes){
    expr_df <- expr_mat %>%
      as.data.frame() %>%
      dplyr::mutate(row_variance = matrixStats::rowVars(as.matrix(.))) %>%
      dplyr::slice_max(order_by = row_variance, n = top_n_genes) %>%
      dplyr::select(-row_variance)
  } else{
    expr_df <- expr_mat %>%
      as.data.frame()
  }
  
  # ---- PCA Calculation ----
  
  pca_results <- stats::prcomp(x = t(expr_df), center = TRUE, scale. = FALSE)
  
  if (skip_plot){
    return(invisible(pca_results))
  }
  
  # ---- Merge PCA co-ordinates with Metadata ----
  
  pca_df <- pca_results$x %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Sample_ID") %>%
    dplyr::inner_join(metadata, by = c("Sample_ID"="Sample_ID"))
  
  # ---- 🖼️ Generate Plots for each group ----
  
  # PC1 and PC2 Variance
  percentVar <- round(100 * summary(pca_results)$importance[2, 1:2])
  group_vars <- base::setdiff(colnames(metadata), "Sample_ID")
  
  all_plots <- list()
  
  for (var in group_vars) {
    
    # Get levels for each group
    n_unique <- length(unique(metadata[[var]]))
    if (n_unique < 2 || n_unique == nrow(metadata)) {
      log_warn(sample = "",
               step = "plot_pca",
               msg =glue::glue("Variable '{var}' is either constant or has one unique value per sample; skipping PCA plot."))
      next
    }
    
    # Define color palette
    pca_palette <- custom_palette[1:length(unique(metadata[[var]]))]
    names(pca_palette) <- as.character(unique(metadata[[var]]))
    
    # PCA Plot
    p <- ggplot2::ggplot(data = pca_df, 
                         mapping = aes(x = PC1, y = PC2, color = factor(.data[[var]]))) +
      ggplot2::geom_point(size = 3, shape = 16) +
      ggrepel::geom_text_repel(ggplot2::aes(label = Sample_ID), show.legend = FALSE) +
      theme_classic() +
      coord_fixed(ratio = 1) +
      ggplot2::labs(color = var, 
                    x = paste0("PC1: ", percentVar[1], "% variance"),
                    y = paste0("PC2: ", percentVar[2], "% variance"), 
                    title = var) +
      custom_theme +
      ggplot2::scale_color_manual(values = pca_palette)
    
    all_plots[[var]] <- p
    
    log_info(sample = "",
             step = "plot_pca",
             msg =glue::glue("Successfully plotted variable : '{var}'."))
  }
  
  # ---- 💾 Save Plots ----
  
  file_extension <- ".pdf"
  file_name <- file.path(output_dir, paste0("PCA_Plot_", filename, file_extension))
  
  # Open multi-page PDF
  grDevices::cairo_pdf(filename = file_name, width = 8, height = 11.5, onefile = TRUE)  
  
  for (var in names(all_plots)) {
    print(all_plots[[var]])  # each plot goes to a new page
  }
  
  dev.off() 
  
  # ---- 🪵 Log Output and Return PCA co-ordinates ----
  
  log_info(sample = "", 
           step = "plot_pca",
           msg =glue::glue("PCA plot saved successfully to : '{file_name}'."))
  
  return(invisible(pca_results))
}

plot_umap <- function(expr_mat, metadata, n_pcs = 50, n_neighbors = NULL,
                      filename, output_dir){
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(expr_mat = expr_mat, metadata = metadata, n_pcs = n_pcs,
                  filename = filename, output_dir = output_dir)
  
  if(!is.null(n_neighbors) && n_neighbors >= ncol(expr_mat)){
    log_error(sample = "",
              step = "plot_umap",
              msg =glue::glue("n_neighbors 'n_neighbors' MUST be lesser than number of samples."))
  }
  
  # Get PC co-ordinates
  pca_results <- plot_pca(expr_mat = expr_mat, metadata = metadata, skip_plot = TRUE,
                          filename = filename, output_dir = output_dir)
  
  # Select top 50 PCs for UMAP
  if(ncol(pca_results$x) > n_pcs){
    X_pcs <- pca_results$x[, 1:n_pcs]
  } else{
    X_pcs <- pca_results$x
  }
  
  # Calculate n_neighbors
  if (is.null(n_neighbors)){
    n_neighbors = base::max(2, base::floor(ncol(expr_mat) * 0.4))
  }
  
  # Get UMAP co-ordinates
  umap_results <- uwot::umap(X = X_pcs, 
                             n_neighbors = n_neighbors)
  
  # Merge UMAP co-ordinates with Metadata ----
  umap_df <- umap_results %>%
    as.data.frame() %>%
    stats::setNames(nm = c("UMAP1", "UMAP2")) %>%
    tibble::rownames_to_column("Sample_ID") %>%
    dplyr::inner_join(metadata, by = c("Sample_ID"="Sample_ID"))
  
  # ---- 🖼️ Generate Plots for each group ----
  
  group_vars <- base::setdiff(colnames(metadata), "Sample_ID")
  
  all_plots <- list()
  
  for (group_var in group_vars) {
    
    # Get levels for each group
    n_unique <- length(unique(metadata[[group_var]]))
    if (n_unique < 2 || n_unique == nrow(metadata)) {
      log_warn(sample = "",
               step = "plot_umap",
               msg =glue::glue("Variable '{group_var}' is either constant or has one unique value per sample; skipping UMAP plot."))
      next
    }
    
    # Define color palette
    umap_palette <- custom_palette[1:length(unique(metadata[[group_var]]))]
    names(umap_palette) <- as.character(unique(metadata[[group_var]]))
    
    # PCA Plot
    p <- ggplot2::ggplot(data = umap_df, 
                         mapping = aes(x = UMAP1, y = UMAP2, color = factor(.data[[group_var]]))) +
      ggplot2::geom_point(size = 3, shape = 16) +
      ggrepel::geom_text_repel(ggplot2::aes(label = Sample_ID), show.legend = FALSE) +
      theme_classic() +
      coord_fixed(ratio = 1) +
      ggplot2::labs(color = group_var, 
                    x = "UMAP1",
                    y = "UMAP2", 
                    title = group_var) +
      custom_theme +
      ggplot2::scale_color_manual(values = umap_palette)
    
    all_plots[[group_var]] <- p
    
    log_info(sample = "",
             step = "plot_umap",
             msg =glue::glue("Successfully plotted variable : '{group_var}'."))
  }
  
  # ---- 💾 Save Plot ----
  
  file_extension <- ".pdf"
  file_name <- file.path(output_dir, paste0("UMAP_Plot_", filename,  file_extension))
  
  # Open multi-page PDF
  grDevices::cairo_pdf(filename = file_name, width = 8, height = 11.5, onefile = TRUE)   
  
  for (var in names(all_plots)) {
    print(all_plots[[var]])  # each plot goes to a new page
  }
  
  dev.off() 
  
  # ---- 🪵 Log Output ----
  
  log_info(sample = "", 
           step = "plot_umap",
           msg =glue::glue("UMAP plot saved successfully to : '{file_name}'."))
  
  return(invisible(NULL))
}

# ---- 🥧 PIE CHART ----

plot_piechart <- function(metadata, segment_col, filename, output_dir, split_col = NULL){
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(metadata = metadata,
                  filename = filename, output_dir = output_dir)
  
  # ---- 👥 Determine Groups for Plotting ----
  
  if (is.null(split_col)) {
    # No splitting → single panel
    group_vars <- "All"
  } else if (length(split_col) == 1) {
    # Single-column split → one panel per unique value
    group_vars <- unique(metadata[[split_col]])
  } else {
    stop("Use Only one column for splitting")
  }
  
  # ---- 🖼️ Generate Plots for each group ----
  
  all_plots <- list()
  
  for (group_var in group_vars) {
    
    if (length(split_col) == 1) {
      # Single-column split → subset by value
      df <- metadata %>% 
        dplyr::filter(.data[[split_col]] == group_var)
    } else{
      df <- metadata
    }
    
    # Determine levels of segment_col
    all_levels <- sort(unique(metadata[[segment_col]]))
    
    # Convert grouping column to factor to ensure proper ordering for colors
    df[[segment_col]] <- factor(df[[segment_col]], levels = all_levels)
    
    # Assign colors, keeping names aligned with all_levels
    pie_palette <- custom_palette[seq_along(all_levels)]
    names(pie_palette) <- all_levels
    
    # Plot title
    title <- ifelse(group_var == "All", "", group_var)
    
    df <- df %>%
      dplyr::count(.data[[segment_col]]) %>%
      dplyr::mutate(Percent = round(100*n/sum(n, na.rm=TRUE), digits = 0), 
                    Percent_label = paste0(Percent,"%")) %>%
      dplyr::arrange(.data[[segment_col]])
    
    # If you run ggplot without themevoid(), you will see 0 and 100% dont overlap.
    p <- ggplot(data = df, 
                mapping = aes(x = "", y = Percent, fill = .data[[segment_col]])) +
      ggplot2::geom_bar(stat = "identity", width = 3, color = "white") +
      ggplot2::coord_polar(theta = "y", start = 0, direction = -1) +
      ggplot2::geom_text(aes(x = 3.2, label = Percent_label), 
                         position = position_stack(vjust = 0.5), 
                         color = "black", size = 3.5, check_overlap = TRUE) +
      scale_fill_manual(values = pie_palette,
                        aesthetics = "fill") +
      ggplot2::labs(title = title,
                    fill = segment_col,
                    x = "",
                    y = "") +
      theme_void() +        #remove background, grid, numeric labels
      custom_theme +
      ggplot2::theme(axis.text.x =  element_blank(),
                     axis.text.y =  element_blank(),
                     strip.text.x = element_text(family = "sans", face = "bold",  colour="black", size=10, hjust = 0.5),
                     legend.position = "none",
                     axis.line = element_blank(),
                     axis.ticks = element_blank())
    
    all_plots[[group_var]] <- p
  }
  
  # ---- Extract shared legend from entire dataset ----
  
  dummy <- ggplot(data = metadata, 
                  mapping = aes(x = 1, fill = .data[[segment_col]])) +
    geom_bar(width = 1) +
    scale_fill_manual(values = pie_palette,
                      aesthetics = "fill") +
    ggplot2::theme(legend.position = "bottom",
                   legend.direction = "horizontal")
  
  shared_legend <- cowplot::get_legend(dummy)
  
  # ---- 🌐 Combine Plots Using cowplot ----
  
  n_plots <- length(all_plots)
  ncol_plots <- ceiling(sqrt(n_plots))
  nrow_plots <- ceiling(n_plots / ncol_plots)
  
  # Restrict unsupported combination
  if (ncol_plots > 10 && nrow_plots > 10) {
    
    log_error(sample = "",
              step = "plot_piechart",
              msg ="Image size too large. More than 100 plots cannot be viewed in a single figure")
  }
  
  # Combine all plots
  combined_plot <- cowplot::plot_grid(plotlist = all_plots, 
                                      ncol = ncol_plots, 
                                      nrow = nrow_plots)
  
  final_plot <- cowplot::plot_grid(plot = combined_plot, 
                                   shared_legend, 
                                   ncol = 1, 
                                   rel_heights = c(1, 0.1))
  
  # ---- 💾 Save Combined Plot ----
  
  file_name <- file.path(output_dir, paste0("Pie_Chart_", filename, ".pdf"))
  ggplot2::ggsave(filename  = file_name,
                  plot      = final_plot,
                  device    = cairo_pdf,
                  width     = ncol_plots * 8, # extra 2 inch for legend
                  height    = nrow_plots * 6, 
                  units     = "in",
                  limitsize = FALSE,
                  bg        = "white")
  
  # ---- 🪵 Log Output ----
  
  log_info(sample = "",
           step = "plot_piechart",
           msg =glue::glue("Pie chart saved successfully to : '{file_name}'."))
}

# ---- DEPRECATED FUNCTIONS ----

# DEPRECATED (used during Seurat v3)
v3_sctransform_singlecell <- function(filtered_seurat){
  
  # Seurat v5 stores counts of each sample in separate layers. Merge them.
  filtered_seurat@assays$RNA <- SeuratObject::JoinLayers(filtered_seurat@assays$RNA)
  
  # Split each sample into a seurat object to get a list of seurat object
  split.seurat <- Seurat::SplitObject(object = filtered_seurat,
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
                                               verbose = FALSE)
    
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
                                        reduction.key = "PC_",
                                        verbose = FALSE)
    
    # Perform dimensional reduction using UMAP on PCA dimensions
    split.seurat[[i]] <- Seurat::RunUMAP(object = split.seurat[[i]],
                                         dims = 1:40,
                                         reduction = "pca",
                                         reduction.name = "umap",
                                         reduction.key = "UMAP_",
                                         verbose = FALSE)
    
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
  kweight2 <- filtered_seurat@meta.data %>%
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
  integrated_seurat.rpca <- Seurat::IntegrateData(anchorset=integ_anchors.rpca,
                                                  new.assay.name="integrated",
                                                  normalization.method="SCT",
                                                  features=NULL,
                                                  features.to.integrate=NULL,
                                                  dims=1:30,
                                                  k.weight=kweight, #default is 100
                                                  weight.reduction=NULL,
                                                  sd.weight=1)
  
  #**STEP 7G: RUN PCA USING 3000 INTEGRATION FEATURES & UMAP USING FIRST 40 PCs**#
  integrated_seurat <- Seurat::RunPCA(object=integrated_seurat,
                                      assay="integrated",
                                      features=NULL)
  
  integrated_seurat <- Seurat::RunUMAP(object=integrated_seurat,
                                       dims=1:40,
                                       reduction="pca")
  return(integrated_seurat)
}

# DEPRECATED (used during Seurat v3)
### Generate whitelist for CITESeq
# Input is filtered seurat object
# Output is a list of csv files - one per batch containing valid barcodes
v3_generate_whitelist <- function(filtered_seurat, output_dir){
  
  # Extract barcodes and split by "_"
  bc <- filtered_seurat@meta.data$Cell
  
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

# DEPRECATED (use decoupleR)
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

# ---- NOTES ----

#******************************************************************************#
#                         TO SUBSET DATA OR NOT BEFORE DESEQ2                          
#******************************************************************************#

# https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
# https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html

# NOTE: If any one group has high within-group variability in PCA plot, we 
# SHOULD exclude those samples by sub-setting before creating dds object and 
# calculating dispersion estimates. Else, use full dataset for modelling.

# NOTE: 
# One factor   : design ~ Condition
# Two factors  : design ~ Condition + Treatment + Condition:Treatment (OR) design ~ Condition*Treatment
# Three factors: design ~ Cell.Line + Condition + Treatment + Condition:Treatment +.. (OR) design ~ Cell.Line*Condition*Treatment
# Using * will include all possible interaction terms and is RECOMMENDED
# The results are the same irrespective of the order i.e 
# Cell.Line*Condition*Treatment gives same results as Cell.Line*Treatment*Condition

# NOTE: If you have Cell.Line, Treatment, Condition variables and want to 
# include them in design, ideally design MUST be ~ Cell.Line*Treatment*Condition.
# However, if one cell line doesnt have a specific treatment, then DESeq2 will
# throw "Full model martix is less than full rank" error. So, ALWAYS create a 
# new column "Comparisons" in metadata and use design ~ Comparisons ALWAYS.
# Populate the "Comparisons" column in metadata by concatenating as below:
# paste0(Cell.Line, Treatment, Condition, sep=".") 

# NOTE: contrast terms CANNOT start with numbers. So, if you have cell line name
# "22RV1" rename it as "RV1_22" etc

# NOTE: You MUST have atleast 2 replicates per group. Else, DESeq2 will throw
# "Error in base::colMeans(x, na.rm = na.rm, dims = dims, ...) : 'x' must be an
# array of at least two dimensions"

# for (n in 1:length(DEG.params$contrast)){
#   
#   if (DEG.params$deseq2.batch.correct == TRUE){
#     
#     # Perform DESeq2() using sva modelled surrogate variables
#     # Create DESeq2 object with surrogate variables in design
#     sva.formula_string <- formula_string
#     sva_dds <- batch_correct_sva(metadata, read_data, sva.formula_string)
#     sva_dds <- run_deseq2(sva_dds, metadata, DEG.params, n, "sva", degs_dir)
#     
#     # Perform DESeq2() using combat corrected counts
#     # Create DESeq2 object with appropiate variables in design
#     # Get combat corrected raw reads
#     read_data_combat <- batch_correct_combat(metadata, read_data, combat.formula_string)
#     combat.formula_string <- formula_string
#     combat_dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data_combat,
#                                                  colData=metadata, 
#                                                  design=~1)
#     design(combat_dds) <- as.formula(combat.formula_string)
#     dds <- run_deseq2(combat_dds, metadata, DEG.params, n, "combat", degs_dir)
#   }
# }

# You have 3 groups A, B & C but want to perform DEG between group A and B only.
# Should you exclude group C samples and perform DEG analysis??
# NOTE: Subsetting data affects (i) sizeFactors (ii) Dispersion estimates and
# hence the final DEGs

# The links below explain when data needs to be subset.
# https://support.bioconductor.org/p/108471/
# https://support.bioconductor.org/p/81038/
# https://support.bioconductor.org/p/9156036/
# https://support.bioconductor.org/p/69341/

# The authors of DESeq2 and EdgeR recommend NOT to subset metadata & 
# read_data but rather use contrasts to get DEGs between groups
# You should subset your data ONLY when one of the groups has large variance 
# (like group C) in link below
# https://master.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#if-i-have-multiple-groups-should-i-run-all-together-or-split-into-pairs-of-groups

# Dispersion of gene X in group A ~ Variance of gene X in group A/mean of gene X group A
# Dispersion of gene X in group B ~ Variance of gene X in group B/mean of gene X group B
# Dispersion of gene X in group C ~ Variance of gene X in group C/mean of gene X group C
# Dispersion of gene X for experiment ~ (Dispersion of gene X in groups A, B, C)
# Samples within Group A and B have low variance in overall expression profile 
# (clustered closely in PCA plot) but samples in group C have different
# expression profile. So, the final dispersion estimates of most genes will be 
# affected a lot by group C. In such a scenario, if we are only doing DEG 
# comparison between groups A & B, then we can subset only group A and B samples
# before doing DESeq2. 
# SO, CHECK PCA PLOT FOR EVERY EXPERIMENT TO DETERMINE IF DATA NEED TO BE SUBSET.

# metadata <- meta %>% dplyr::filter(get(Variable) %in% c(Comparisons$Target[n], Comparisons$Reference[n]))
# corrected_read_data <- read %>% dplyr::select(rownames(metadata))
# NOTE: sizefactors MUST be calculated ONLY using samples being compared for 
# differential expression. So, make sure read_data and metadata ONLY have
# samples being compared for differential expression.

# dds <- DESeq2::DESeqDataSetFromMatrix(countData=corrected_read_data, colData=metadata, design=~ Tissue + Age)
# dds <- DESeq(dds, test="LRT", reduced=~ Tissue)
# res <- DESeq2::results(object=dds)

# Study effect of Treatment ignoring effect of Sex
#dds <- DESeqDataSetFromMatrix(countData=corrected_read_data, colData=metadata, design=~ Treatment)

# Study effect of Treatment after removing effect of Sex (assumes equal effect of Sex on different treatments)
#dds <- DESeqDataSetFromMatrix(countData=corrected_read_data, colData=metadata, design=~ Sex + Treatment)

# Study effect of Sex on treatment?
# Study effect of Treatment after removing effect of Sex (assumes different effect of Sex on different treatments)
#dds <- DESeqDataSetFromMatrix(countData=corrected_read_data, colData=metadata, design=~ Sex + Treatment + Sex:Treatment)

# #*******************************DIAGNOSTIC TESTS DESEQ2 *******************************#

# YET TO IMPLEMENT

# metadata MUST be xlsx file, MUST have "Sample_ID" column with sample names
# without any duplication
# Example: If samples were isolated on different days & prepared using different
# kits, "Batch" column must have values like "1_Ribo", "1_Poly-A", "2_Ribo", etc
# NOTE: Make sure there are no white spaces in the Target and Reference columns
# in excel file. R will change any white space (" ") to dot ("."). So, replace 
# white space (" ") with underscore ("_") before importing into R.

# NOTE: The normalized counts you get from sva_dds and DESeq2 dds are NOT batch 
# corrected and they will be identical. Both DESeq2() and svaseq() ONLY model
# for batch factors/surrogate variables in the design formula, they DO NOT 
# modify counts like combatseq. The normalized counts you get from combatseq
# are batch corrected. 

# for (n in 1:length(Comparisons$Variable)){
#   
#   # This generates a new column "id" that has info on samples being comparared
#   metadata_comp <- metadata %>%
#     dplyr::mutate(id=get(Comparisons$Variable[n]))
#   
#   #combat_corrected_read_data <- combatseq_batch(read_data, metadata_comp)
#   #sva_dds <- svaseq_batch(read_data, metadata_comp)
#   
#   # Perform DESeq2() using in-built batch modelling
#   approach <- "DESeq2_modelled"
#   if (length(unique(metadata_comp$Batch)) > 1){
#     dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
#                                           colData=metadata_comp, 
#                                           design=~ Batch+id)
#   } else {
#     dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
#                                           colData=metadata_comp, 
#                                           design=~ id)
#   }
#   dds <- run_deseq2(dds, metadata_comp, annotations, Comparisons, n, approach, prefix, results_path)
#   deseq2_norm_counts(dds, annotations, approach, suffix) # batch corrected if you more than 1 batch
#   plot_qc(dds, metadata_comp, approach, suffix)
#   
#   # Perform DESeq2() using combat corrected read counts
#   if (!identical(read_data, combat_corrected_read_data)){
#     approach <- "combat_corrected"
#     combat_dds <- DESeq2::DESeqDataSetFromMatrix(countData=combat_corrected_read_data,
#                                                  colData=metadata_comp, 
#                                                  design=~ id)
#     combat_dds <- run_deseq2(combat_dds, metadata_comp, annotations, Comparisons, n, lfc.cutoff, padj.cutoff, suffix)
#     combatseq_norm_counts(combat_dds, annotations, approach, suffix)  #combat batch corrected
#     plot_qc(combat_dds, metadata_comp, approach, suffix)
#   }
#   
#   # Perform DESeq2() using sva modelled surrogate variables SV1 and SV2
#   approach <- "sva_modelled"
#   sva_dds <- run_deseq2(sva_dds, metadata_comp, annotations, Comparisons, n, lfc.cutoff, padj.cutoff, approach, suffix)
#   # calc_norm_counts(sva_dds, annotations, approach, suffix)   # uncorrected
#   svaseq_norm_counts(sva_dds, annotations, approach, suffix)   # sva seq batch corrected
#   plot_qc(sva_dds, metadata_comp, approach, suffix)
# }


# # (i) To view counts of specific gene across samples
# plotCounts(dds, gene=which.min(res$padj), intgroup=Variable)           # gene with lowest padj
# plotCounts(dds, gene=which.min(res$log2FoldChange), intgroup=Variable) # gene with lowest log2FC
# plotCounts(dds, gene=which.max(res$log2FoldChange), intgroup=Variable) # gene with highest log2FC
# plotCounts(dds, gene="ENSMUSG00000030598", intgroup=Variable)          # a specific gene
# 

# # To identify the genes interactively, run the 2 lines below. 
# # Then click on multiple dots and click Finish. A list of genes will be displayed 
# # idx <- identify(res$baseMean, res$log2FoldChange)
# # rownames(res)[idx]

# # (iv) Hierarchical clustering of samples using rld or vst
# rld_mat <- DESeq2::assay(rld)     # extract the matrix from a DESeq2 object
# rld_cor <- cor(x=rld_mat,       # compute pairwise correlation values
#                y=NULL,
#                use="everything",
#                method="pearson") 
# 
# # Check the output of cor(), make note of the rownames and colnames
# head(rld_cor)
# 
# vst_mat <- assay(vst)  
# vst_cor <- cor(vst_mat) 
# 
# # Check the output of cor(), make note of the rownames and colnames
# head(vst_cor)    
# 
# for (process in c("rld", "vst")){
#   
#   pheatmap::pheatmap(mat=get(paste0(process, "_cor")),
#                      color=colorRampPalette(rev(brewer.pal(n=11, name ="RdYlBu")))(100),
#                      breaks=NA, 
#                      border_color="white", #"grey60",
#                      cellwidth=NA, 
#                      cellheight=NA, 
#                      scale="none",   
#                      cluster_rows=TRUE,   #cluster the rows
#                      cluster_cols=TRUE,   #cluster the columns
#                      clustering_distance_rows="euclidean",
#                      clustering_distance_cols="euclidean",
#                      clustering_method="complete",
#                      legend=TRUE, 
#                      legend_breaks=NA,
#                      legend_labels=NA, 
#                      #annotation_row=,  
#                      #annotation_col=, 
#                      annotation_colors=dplyr::if_else(nrow(col_annotation)+nrow(row_annotation) > 0, ann_colors, NA),
#                      annotation_legend=TRUE,
#                      annotation_names_row=TRUE,
#                      annotation_names_col=TRUE,
#                      show_rownames=dplyr::if_else(nrow(mat)<80, TRUE, FALSE, missing=NULL), 
#                      show_colnames=dplyr::if_else(ncol(mat)<30, TRUE, FALSE, missing=NULL),
#                      fontsize=8, 
#                      fontsize_row=8, 
#                      fontsize_col=8,
#                      angle_col=c("270", "0", "45", "90", "315"),
#                      fontsize_number=0.8*fontsize, 
#                      #labels_row=display_row,
#                      #labels_col=display_col,
#                      filename=paste0(diagnostics_path, "Diagnostic_Correlation_Heatmap_using_", process, "_", celltype, ".pdf"))
# }
# 
# # (v) Checking if mean < variance (for NB model) or Mean=Variance (for Poisson model). 
# # Each point is a gene denoted by (x,y) where 
# # x=mean_count of gene across all samples & y=variance_count of gene across all samples
# 
# mean_counts <- apply(read_data[,], 1, mean)
# variance_counts <- apply(read_data[,], 1, var)
# df <- data.frame(mean_counts, variance_counts)
# ggplot(df) +
#   geom_point(aes(x=mean_counts, y=variance_counts)) +
#   geom_line(aes(x=mean_counts, y=mean_counts, color="red")) +
#   scale_y_log10() +
#   scale_x_log10()
# 
# ggplot2::ggsave(filename=paste0("Diagnostic_Scatter_plot_", celltype, ".pdf"),
#                 plot=last_plot(),
#                 device="pdf",
#                 path=diagnostics_path,
#                 scale=1,
#                 #width=8.5,
#                 #height=11,
#                 units=c("in"),
#                 dpi=300,
#                 limitsize=TRUE,
#                 bg="white")
# 
# # # (vii) Plot p-value histogram
# # hist(res$pvalue, col="lightblue", breaks=20)

# # # (ix) Plot a histogram for one of the samples to see how the counts are distributed. Adjust xlim if needed
# # ggplot(read_data, aes(x=results.S10_R1_trimmed.fastq.gz.csv)) +
# #   geom_histogram(stat="bin", bins=200) +
# #   xlim(-5,1000) +
# #   xlab("Raw expression counts") +
# #   ylab("Number of genes") +
# #   scale_color_brewer(palette="Dark2")
# # 
# # # (x) Extract counts, size factors, etc.. from dds
# # dds <- estimateSizeFactors(dds)         #redundant if you already ran dds <- DESeq(dds)
# # counts <- counts(dds)                   #counts[i,j]=raw_count of gene i in sample j
# # sizefactors <- sizeFactors(dds)         #sizefactors[j]=median (counts[,j]/geometric_mean(counts[i,]))
# # colSums(counts(dds))                    #Total number of raw counts per sample
# # colSums(counts(dds, normalized=TRUE))   #Total number of normalized counts per sample

# ---- CODE TESTING ---- 
# 
# DEG.params  <- list(Variable    = c("Comparisons"),
#                     Target      = c("ARCaPM.4Gy.NDRG1_mut"),
#                     Reference   = c("ARCaPM.0Gy.NDRG1_mut"),
#                     contrast    = c("ARCaPM.NDRG1_mut.4Gy-ARCaPM.NDRG1_mut.0Gy"),
#                     lfc.cutoff  = 0,
#                     padj.cutoff = 0.1,
#                     design      = "Cell.Line*Treatment*Condition", #"Comparisons", #"0+Condition:Treatment"
#                     design.ref  = c("Condition:WT", "Treatment:0Gy", "Comparisons:ARCaPM.0Gy.WT"),
#                     deseq2.batch.correct = FALSE,
#                     proj        = "RNASeq_Manish_22RV1_ARCaPM",
#                     species     = "Homo sapiens")
# 
# metadata <- metadata %>% filter(Cell.Line != "22RV1")
# n <- 1
# 
# 
# dds.new <- run_deseq2(dds, metadata, DEG.params, n, approach, data_dir)
# design(dds) <- as.formula(~Comparisons)
# dds.old <- run_deseq2_old(dds, metadata, DEG.params, n, approach, data_dir)
# 
# 
# # The following 3 approaches give identical results
# 
# # Method 1 => Similar way of defining contrasts like method2. Easy to compare 
# # samples but difference of difference not possible
# design <- "0+Condition:Treatment" 
# dds.test <- dds
# contrast1 <- c("Condition", "NDRG1_mut.Treatment4Gy", "NDRG1_mut.Treatment0Gy")
# contrast2 <- c("Condition", "NDRG1_mut.Treatment4Gy", "WT.Treatment0Gy")
# contrast3 <- c("Condition", "NDRG1_mut.Treatment4Gy", "WT.Treatment4Gy")
# contrast4 <- c("Condition", "NDRG1_mut.Treatment0Gy", "WT.Treatment4Gy")
# contrast5 <- c("Condition", "NDRG1_mut.Treatment0Gy", "WT.Treatment0Gy")
# contrast6 <- c("Condition", "WT.Treatment4Gy", "WT.Treatment0Gy")
# res1 <- DESeq2::results(dds, contrast = contrast1)
# res2 <- DESeq2::results(dds, contrast = contrast2)
# res3 <- DESeq2::results(dds, contrast = contrast3)
# res4 <- DESeq2::results(dds, contrast = contrast4)
# res5 <- DESeq2::results(dds, contrast = contrast5)
# res6 <- DESeq2::results(dds, contrast = contrast6)
# 
# # Method 2= > combine COndition and Treatment to new column Comparisons. 
# # Similar way of defining contrasts like method1. Easy to compare 
# # samples but difference of difference not possible
# design <- "Comparisons"
# contrast1A <- c("Comparisons", "ARCaPM.4Gy.NDRG1_mut", "ARCaPM.0Gy.NDRG1_mut")
# contrast2A <- c("Comparisons", "ARCaPM.4Gy.NDRG1_mut", "ARCaPM.0Gy.WT")
# contrast3A <- c("Comparisons", "ARCaPM.4Gy.NDRG1_mut", "ARCaPM.4Gy.WT")
# contrast4A <- c("Comparisons", "ARCaPM.0Gy.NDRG1_mut", "ARCaPM.4Gy.WT")
# contrast5A <- c("Comparisons", "ARCaPM.0Gy.NDRG1_mut", "ARCaPM.0Gy.WT")
# contrast6A <- c("Comparisons", "ARCaPM.4Gy.WT", "ARCaPM.0Gy.WT")
# res1A <- DESeq2::results(dds, contrast = contrast1A)
# res2A <- DESeq2::results(dds, contrast = contrast2A)
# res3A <- DESeq2::results(dds, contrast = contrast3A)
# res4A <- DESeq2::results(dds, contrast = contrast4A)
# res5A <- DESeq2::results(dds, contrast = contrast5A)
# res6A <- DESeq2::results(dds, contrast = contrast6A)
# 
# # Method 3 => conventional design, contrast needs to calculated for each sample
# # difference of difference is easy
# design <- "Condition+Treatment+Condition:Treatment"
# dds.standard <- dds
# mod_mat <- model.matrix(design(dds), colData(dds))
# NDRG1_0Gy <- colMeans(mod_mat[dds$Condition == "NDRG1_mut" & dds$Treatment == "0Gy",])
# NDRG1_4Gy <- colMeans(mod_mat[dds$Condition == "NDRG1_mut" & dds$Treatment == "4Gy",])
# WT_0Gy <- colMeans(mod_mat[dds$Condition == "WT" & dds$Treatment == "0Gy",])
# WT_4Gy <- colMeans(mod_mat[dds$Condition == "WT" & dds$Treatment == "4Gy",])
# 
# contrast1B <- NDRG1_4Gy - NDRG1_0Gy 
# contrast2B <- NDRG1_4Gy - WT_0Gy
# contrast3B <- NDRG1_4Gy - WT_4Gy
# contrast4B <- NDRG1_0Gy - WT_4Gy 
# contrast5B <- NDRG1_0Gy - WT_0Gy
# contrast6B <- WT_4Gy - WT_0Gy
# contrast7B <- (NDRG1_4Gy - NDRG1_0Gy) - (WT_4Gy - WT_0Gy)
# contrast8B <- 
#   
#   res1B. <- DESeq2::results(dds, contrast = contrast1B)
# res2B. <- DESeq2::results(dds, contrast = contrast2B)
# res3B. <- DESeq2::results(dds, contrast = contrast3B)
# res4B. <- DESeq2::results(dds, contrast = contrast4B)
# res5B. <- DESeq2::results(dds, contrast = contrast5B)
# res6B. <- DESeq2::results(dds, contrast = contrast6B)
# res7B. <- DESeq2::results(dds, contrast = contrast7B)
# 
# 
# # Check for differences in results
# df1 <- data.frame(res1) %>% dplyr::mutate_all(~(round(.,2)))
# df1A <- data.frame(res1A) %>% dplyr::mutate_all(~(round(.,2)))
# df1B <- data.frame(res1B) %>% dplyr::mutate_all(~(round(.,2)))
# 
# df2 <- data.frame(res2) %>% dplyr::mutate_all(~(round(.,2)))
# df2A <- data.frame(res2A) %>% dplyr::mutate_all(~(round(.,2)))
# df2B <- data.frame(res2B) %>% dplyr::mutate_all(~(round(.,2)))
# 
# df3 <- data.frame(res3) %>% dplyr::mutate_all(~(round(.,2)))
# df3A <- data.frame(res3A) %>% dplyr::mutate_all(~(round(.,2)))
# df3B <- data.frame(res3B) %>% dplyr::mutate_all(~(round(.,2)))
# 
# df4 <- data.frame(res4) %>% dplyr::mutate_all(~(round(.,2)))
# df4A <- data.frame(res4A) %>% dplyr::mutate_all(~(round(.,2)))
# df4B <- data.frame(res4B) %>% dplyr::mutate_all(~(round(.,2)))
# 
# df5 <- data.frame(res5) %>% dplyr::mutate_all(~(round(.,2)))
# df5A <- data.frame(res5A) %>% dplyr::mutate_all(~(round(.,2)))
# df5B <- data.frame(res5B) %>% dplyr::mutate_all(~(round(.,2)))
# 
# df6 <- data.frame(res6) %>% dplyr::mutate_all(~(round(.,2)))
# df6A <- data.frame(res6A) %>% dplyr::mutate_all(~(round(.,2)))
# df6B <- data.frame(res6B) %>% dplyr::mutate_all(~(round(.,2)))
# 
# # Display only rows and columns that are different
# dplyr::full_join(df1[rowSums(df1-df1A, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
#                  df1A[rowSums(df1-df1A, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
#                  by=c("SYMBOL"="SYMBOL")) %>%
#   dplyr::transmute(SYMBOL              = SYMBOL,
#                    baseMean.diff       = baseMean.x-baseMean.y,
#                    log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
#                    lfcSE.diff          = lfcSE.x-lfcSE.y,
#                    stat.diff           = stat.x-stat.y,
#                    pvalue.diff         = pvalue.x-pvalue.y,
#                    padj.diff           = padj.x-padj.y)
# 
# dplyr::full_join(df2[rowSums(df2-df2A, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
#                  df2A[rowSums(df2-df2A, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
#                  by=c("SYMBOL"="SYMBOL")) %>%
#   dplyr::transmute(SYMBOL              = SYMBOL,
#                    baseMean.diff       = baseMean.x-baseMean.y,
#                    log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
#                    lfcSE.diff          = lfcSE.x-lfcSE.y,
#                    stat.diff           = stat.x-stat.y,
#                    pvalue.diff         = pvalue.x-pvalue.y,
#                    padj.diff           = padj.x-padj.y)
# 
# dplyr::full_join(df3[rowSums(df3-df3A, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
#                  df3A[rowSums(df3-df3A, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
#                  by=c("SYMBOL"="SYMBOL")) %>%
#   dplyr::transmute(SYMBOL              = SYMBOL,
#                    baseMean.diff       = baseMean.x-baseMean.y,
#                    log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
#                    lfcSE.diff          = lfcSE.x-lfcSE.y,
#                    stat.diff           = stat.x-stat.y,
#                    pvalue.diff         = pvalue.x-pvalue.y,
#                    padj.diff           = padj.x-padj.y)
# 
# dplyr::full_join(df4[rowSums(df4-df4A, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
#                  df4A[rowSums(df4-df4A, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
#                  by=c("SYMBOL"="SYMBOL")) %>%
#   dplyr::transmute(SYMBOL              = SYMBOL,
#                    baseMean.diff       = baseMean.x-baseMean.y,
#                    log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
#                    lfcSE.diff          = lfcSE.x-lfcSE.y,
#                    stat.diff           = stat.x-stat.y,
#                    pvalue.diff         = pvalue.x-pvalue.y,
#                    padj.diff           = padj.x-padj.y)
# 
# dplyr::full_join(df5[rowSums(df5-df5A, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
#                  df5A[rowSums(df5-df5A, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
#                  by=c("SYMBOL"="SYMBOL")) %>%
#   dplyr::transmute(SYMBOL              = SYMBOL,
#                    baseMean.diff       = baseMean.x-baseMean.y,
#                    log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
#                    lfcSE.diff          = lfcSE.x-lfcSE.y,
#                    stat.diff           = stat.x-stat.y,
#                    pvalue.diff         = pvalue.x-pvalue.y,
#                    padj.diff           = padj.x-padj.y)
# 
# dplyr::full_join(df6[rowSums(df6-df6A, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
#                  df6A[rowSums(df6-df6A, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
#                  by=c("SYMBOL"="SYMBOL")) %>%
#   dplyr::transmute(SYMBOL              = SYMBOL,
#                    baseMean.diff       = baseMean.x-baseMean.y,
#                    log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
#                    lfcSE.diff          = lfcSE.x-lfcSE.y,
#                    stat.diff           = stat.x-stat.y,
#                    pvalue.diff         = pvalue.x-pvalue.y,
#                    padj.diff           = padj.x-padj.y)
# 
# dplyr::full_join(df1[rowSums(df1-df1B, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
#                  df1B[rowSums(df1-df1B, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
#                  by=c("SYMBOL"="SYMBOL")) %>%
#   dplyr::transmute(SYMBOL              = SYMBOL,
#                    baseMean.diff       = baseMean.x-baseMean.y,
#                    log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
#                    lfcSE.diff          = lfcSE.x-lfcSE.y,
#                    stat.diff           = stat.x-stat.y,
#                    pvalue.diff         = pvalue.x-pvalue.y,
#                    padj.diff           = padj.x-padj.y)
# 
# dplyr::full_join(df2[rowSums(df2-df2B, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
#                  df2B[rowSums(df2-df2B, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
#                  by=c("SYMBOL"="SYMBOL")) %>%
#   dplyr::transmute(SYMBOL              = SYMBOL,
#                    baseMean.diff       = baseMean.x-baseMean.y,
#                    log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
#                    lfcSE.diff          = lfcSE.x-lfcSE.y,
#                    stat.diff           = stat.x-stat.y,
#                    pvalue.diff         = pvalue.x-pvalue.y,
#                    padj.diff           = padj.x-padj.y)
# 
# dplyr::full_join(df3[rowSums(df3-df3B, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
#                  df3B[rowSums(df3-df3B, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
#                  by=c("SYMBOL"="SYMBOL")) %>%
#   dplyr::transmute(SYMBOL              = SYMBOL,
#                    baseMean.diff       = baseMean.x-baseMean.y,
#                    log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
#                    lfcSE.diff          = lfcSE.x-lfcSE.y,
#                    stat.diff           = stat.x-stat.y,
#                    pvalue.diff         = pvalue.x-pvalue.y,
#                    padj.diff           = padj.x-padj.y)
# 
# dplyr::full_join(df4[rowSums(df4-df4B, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
#                  df4B[rowSums(df4-df4B, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
#                  by=c("SYMBOL"="SYMBOL")) %>%
#   dplyr::transmute(SYMBOL              = SYMBOL,
#                    baseMean.diff       = baseMean.x-baseMean.y,
#                    log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
#                    lfcSE.diff          = lfcSE.x-lfcSE.y,
#                    stat.diff           = stat.x-stat.y,
#                    pvalue.diff         = pvalue.x-pvalue.y,
#                    padj.diff           = padj.x-padj.y)
# 
# dplyr::full_join(df5[rowSums(df5-df5B, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
#                  df5B[rowSums(df5-df5B, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
#                  by=c("SYMBOL"="SYMBOL")) %>%
#   dplyr::transmute(SYMBOL              = SYMBOL,
#                    baseMean.diff       = baseMean.x-baseMean.y,
#                    log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
#                    lfcSE.diff          = lfcSE.x-lfcSE.y,
#                    stat.diff           = stat.x-stat.y,
#                    pvalue.diff         = pvalue.x-pvalue.y,
#                    padj.diff           = padj.x-padj.y)
# 
# dplyr::full_join(df6[rowSums(df6-df6B, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
#                  df6B[rowSums(df6-df6B, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
#                  by=c("SYMBOL"="SYMBOL")) %>%
#   dplyr::transmute(SYMBOL              = SYMBOL,
#                    baseMean.diff       = baseMean.x-baseMean.y,
#                    log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
#                    lfcSE.diff          = lfcSE.x-lfcSE.y,
#                    stat.diff           = stat.x-stat.y,
#                    pvalue.diff         = pvalue.x-pvalue.y,
#                    padj.diff           = padj.x-padj.y)
# 
# 
# 
# res1 <- DESeq2::lfcShrink(dds=dds.test, res=res1, type="ashr")
# res1B <- DESeq2::lfcShrink(dds=dds.standard, res=res1B, type="ashr")
# res2 <- DESeq2::lfcShrink(dds=dds.test, res=res2, type="ashr")
# res2B <- DESeq2::lfcShrink(dds=dds.standard, res=res2B, type="ashr")
# res3 <- DESeq2::lfcShrink(dds=dds.test, res=res3, type="ashr")
# res3B <- DESeq2::lfcShrink(dds=dds.standard, res=res3B, type="ashr")
# res4 <- DESeq2::lfcShrink(dds=dds.test, res=res4, type="ashr")
# res4B <- DESeq2::lfcShrink(dds=dds.standard, res=res4B, type="ashr")
# res5 <- DESeq2::lfcShrink(dds=dds.test, res=res5, type="ashr")
# res5B <- DESeq2::lfcShrink(dds=dds.standard, res=res5B, type="ashr")
# res6 <- DESeq2::lfcShrink(dds=dds.test, res=res6, type="ashr")
# res6B <- DESeq2::lfcShrink(dds=dds.standard, res=res6B, type="ashr")
# res7B <- DESeq2::lfcShrink(dds=dds.standard, res=res7B, type="ashr")
# 
# 
# summary(res1)
# summary(res1B)
# summary(res1B.)
# summary(res2)
# summary(res2B)
# summary(res2B.)
# summary(res3)
# summary(res3B)
# summary(res3B.)
# summary(res4)
# summary(res4B)
# summary(res4B.)
# summary(res5)
# summary(res5B)
# summary(res5B.)
# summary(res6)
# summary(res6B)
# summary(res6B.)
# summary(res7B)
# summary(res7B.)
# res1B. <- DESeq2::lfcShrink(dds=dds, res=res1B., type="ashr")
# res2B. <- DESeq2::lfcShrink(dds=dds, res=res2B., type="ashr")
# res3B. <- DESeq2::lfcShrink(dds=dds, res=res3B., type="ashr")
# res4B. <- DESeq2::lfcShrink(dds=dds, res=res4B., type="ashr")
# res5B. <- DESeq2::lfcShrink(dds=dds, res=res5B., type="ashr")
# res6B. <- DESeq2::lfcShrink(dds=dds, res=res6B., type="ashr")
# res7B. <- DESeq2::lfcShrink(dds=dds, res=res7B., type="ashr")
# 
# setdiff(rownames(data.frame(res1B) %>% filter(log2FoldChange > 0.58, padj < 0.1)),
#         rownames(data.frame(res1B.) %>% filter(log2FoldChange > 0.58, padj < 0.1))) 
# 
# setdiff(rownames(data.frame(res1B) %>% filter(log2FoldChange < -0.58, padj < 0.1)),
#         rownames(data.frame(res1B.) %>% filter(log2FoldChange < -0.58, padj < 0.1)))
# 
# setdiff(rownames(data.frame(res2B) %>% filter(log2FoldChange > 0.58, padj < 0.1)),
#         rownames(data.frame(res2B.) %>% filter(log2FoldChange > 0.58, padj < 0.1))) 
# 
# setdiff(rownames(data.frame(res2B) %>% filter(log2FoldChange < -0.58, padj < 0.1)),
#         rownames(data.frame(res2B.) %>% filter(log2FoldChange < -0.58, padj < 0.1))) 
# 
# 
# setdiff(rownames(data.frame(res3B) %>% filter(log2FoldChange > 0.58, padj < 0.1)),
#         rownames(data.frame(res3B.) %>% filter(log2FoldChange > 0.58, padj < 0.1))) 
# 
# setdiff(rownames(data.frame(res3B) %>% filter(log2FoldChange < -0.58, padj < 0.1)),
#         rownames(data.frame(res3B.) %>% filter(log2FoldChange < -0.58, padj < 0.1))) 
# 
# 
# setdiff(rownames(data.frame(res4B) %>% filter(log2FoldChange > 0.58, padj < 0.1)),
#         rownames(data.frame(res4B.) %>% filter(log2FoldChange > 0.58, padj < 0.1))) 
# 
# setdiff(rownames(data.frame(res4B) %>% filter(log2FoldChange < -0.58, padj < 0.1)),
#         rownames(data.frame(res4B.) %>% filter(log2FoldChange < -0.58, padj < 0.1))) 
# 
# setdiff(rownames(data.frame(res5B) %>% filter(log2FoldChange > 0.58, padj < 0.1)),
#         rownames(data.frame(res5B.) %>% filter(log2FoldChange > 0.58, padj < 0.1))) 
# 
# setdiff(rownames(data.frame(res5B) %>% filter(log2FoldChange < -0.58, padj < 0.1)),
#         rownames(data.frame(res5B.) %>% filter(log2FoldChange < -0.58, padj < 0.1))) 
# 
# 
# setdiff(rownames(data.frame(res6B) %>% filter(log2FoldChange > 0.58, padj < 0.1)),
#         rownames(data.frame(res6B.) %>% filter(log2FoldChange > 0.58, padj < 0.1))) 
# 
# setdiff(rownames(data.frame(res6B) %>% filter(log2FoldChange < -0.58, padj < 0.1)),
#         rownames(data.frame(res6B.) %>% filter(log2FoldChange < -0.58, padj < 0.1))) 
# 
# setdiff(rownames(data.frame(res7B) %>% filter(log2FoldChange > 0.58, padj < 0.1)),
#         rownames(data.frame(res7B.) %>% filter(log2FoldChange > 0.58, padj < 0.1))) 
# 
# setdiff(rownames(data.frame(res7B) %>% filter(log2FoldChange < -0.58, padj < 0.1)),
#         rownames(data.frame(res7B.) %>% filter(log2FoldChange < -0.58, padj < 0.1)))