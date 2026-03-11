# ---- 📦 Load Required Libraries ----

suppressPackageStartupMessages({
  library(crayon)  # Provides green, yellow, red, blue, magenta
  library(glue)    # Provides glue
  library(rlang)   # Provides the %||% operator
  library(readr)
  library(readxl)
})

# To make colors work even in Nextflow logs, force crayon to be active
options(crayon.enabled = TRUE)

# ---- 📝 LOGGING RELATED FUNCTIONS ----

# Suppress warnings and messages
quiet_msg <- function(expr) {
  tmp <- tempfile()
  con <- file(tmp, open = "wt")
  
  on.exit({
    sink(type = "message")
    sink(type = "output")
    close(con)
    if (file.exists(tmp)) unlink(tmp)
  }, add = TRUE)
  
  sink(con, type = "output")
  sink(con, type = "message")
  
  # Capture the result explicitly so it's returned out of the sink
  result <- force(expr)
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


# --- SMART FILE LOADER FUNCTION ----

load_smart <- function(input_path, sheet = NULL) {
  
  # 1. If it's already an object (data.frame/tibble), just return it
  if (!is.character(input_path) || length(input_path) != 1) {
    return(input_path)
  }
  
  if (!file.exists(input_path)) {
    log_error(sample = "", step = "load_smart", msg = glue::glue("File does not exist: '{input_path}'"))
  }
  
  ext <- tools::file_ext(tolower(input_path))
  
  # 2. Smart Loading with Sheet Support
  # 2. Smart Loading with Sheet Support
  data <- switch(EXPR = ext,
                 "csv"  = readr::read_csv(input_path, show_col_types = FALSE),
                 "rds"  = readRDS(input_path),
                 "txt"  = readr::read_tsv(input_path, show_col_types = FALSE),
                 "xlsx" = , # Fall through to xls
                 "xls"  = {
                   # If sheet is NULL, read_excel defaults to the first sheet
                   if (!is.null(sheet)) {
                     readxl::read_excel(input_path, sheet = sheet)
                   } else {
                     readxl::read_excel(input_path)
                   }
                 },
                 log_error(sample = "", step = "load_smart", msg = glue::glue("Unsupported file extension: '{ext}'"))
  )
  
  return(data)
}

# ---- 📝 ANNOTATION FUNCTION ----

add_annotation <- function(df, tx2gene) {
  
  # 1. Load and Clean Annotation Map
  tx2gene <- tx2gene %>% 
    dplyr::select(gene_id = 2, gene_name = 3, gene_biotype = 4) %>%
    dplyr::distinct(gene_id, .keep_all = TRUE) %>%
    # Logical fix: fill missing names ONLY for protein_coding
    dplyr::mutate(gene_name = dplyr::case_when((is.na(gene_name) | gene_name == "") & (gene_biotype == "protein_coding") ~ gene_id,
                                               TRUE ~ gene_name))
  
  # 2. Heuristic Column Matching
  # Pull rownames into a column to check them alongside actual columns
  test_df <- df %>% tibble::rownames_to_column("row_id")
  
  overlap_counts <- sapply(test_df, function(x) sum(as.character(x) %in% tx2gene$gene_id))
  
  max_overlap_col <- names(overlap_counts)[which.max(overlap_counts)]
  overlap_pct     <- max(overlap_counts) / nrow(df)
  pct_label       <- round(overlap_pct * 100, 1)
  
  if (overlap_pct >= 0.5) {
    log_warn(sample = "", 
             step = "add_annotation", 
             msg = glue::glue("Success: Matched IDs in '{max_overlap_col}' ({pct_label}% overlap)"))
  } else {
    log_error(sample = "", 
              step = "add_annotation", 
              msg = "CRITICAL ERROR: Annotation failed due to low ID overlap ({pct_label}% overlap)")
  }
  
  # 4. Final Join and Formatting
  annotated_df <- test_df %>%
    dplyr::rename(gene_id = .data[[max_overlap_col]]) %>%
    # Drop "row_id" only if it wasn't the column we just renamed to "gene_id"
    dplyr::select(-any_of("row_id")) %>%
    dplyr::mutate(gene_id = as.character(gene_id)) %>%
    dplyr::left_join(tx2gene, by = c("gene_id" = "gene_id")) %>%
    dplyr::relocate(gene_id, gene_name, gene_biotype) %>%
    dplyr::rename(SYMBOL = gene_name)
  
  return(annotated_df)
}

# ---- FORMAT COUNTS/DEGS AND ADD ANNOTATION FUNCTION ----

process_counts <- function(expr_obj, tx2gene) {
  
  # Calculate row-wise average expression to identify dominant ID
  row_means <- rowMeans(as.matrix(expr_obj), na.rm = TRUE)
  # Optional: handle cases where a gene has 0 counts in all samples
  row_means[is.na(row_means)] <- 0
  
  # Convert to DF and initial cleanup
  df <- as.data.frame(expr_obj) %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::mutate(row_mean = row_means)
  
  ann_df <- df %>%
    add_annotation(tx2gene) %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    dplyr::group_by(SYMBOL) %>%
    # Pick the ID with the highest average expression
    dplyr::slice_max(order_by = row_mean, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(-row_mean)
  
  return(ann_df)
}

process_degs <- function(res_obj, tx2gene) {
  
  # Convert to DF and initial cleanup
  df <- as.data.frame(res_obj) %>%
    tibble::rownames_to_column("gene_id")
  
  # Calculate the minimum non-zero padj safely for this specific dataset
  # We do this outside mutate to avoid the 'padj not found' error
  min_padj <- min(df$padj[df$padj > 0], na.rm = TRUE)
  
  ann_df <- df %>%
    dplyr::mutate(padj = case_when(is.na(padj) ~ 1,
                                   padj == 0   ~ min_padj,
                                   TRUE        ~ padj)) %>%
    add_annotation(tx2gene) %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::slice_min(order_by = padj, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  return(ann_df)
}

# ---- 📝 HEATMAP FUNCTION ----

# Data type & transform control
# data_type controls log transform, z-scoring, and color break behavior.
# Use "auto" to let the function detect the data type from the matrix values.
# Override with explicit type if auto-detection is wrong:
#   "counts"      - Raw or normalized counts (TPM, CPM); wide positive range
#                   → log2(1+x) if not already logged, then z-score
#   "counts_log"  - Already log-transformed counts (VST, rlog, log2TPM)
#                   → skip log, apply z-score only
#   "counts_scaled" - Already log-transformed AND z-scored
#                   → no transforms applied
#   "centered"    - Zero-centered scores (NES, LFC); has negatives, range > ±1
#                   → no log, no z-score; symmetric breaks at ±99th percentile
#   "correlation" - Correlation matrices; values within [-1, 1]
#                   → no log, no z-score; fixed breaks -1 to 1
#   "bounded"     - Values bounded within [0, 1] (beta, p-values, binary)
#                   → no log, no z-score; fixed breaks 0 to 1

# Clustering & ordering
# "all" = cluster all cols together; 
# "none" = preserve input order;
# "alphabetical" = sort by name; 
# a metadata_col column name to cluster within each group separately

# Plot aesthetics
# "rdbu" (diverging, good for centered/correlation data) or
# "viridis" (sequential, good for counts/bounded data)

plot_heatmap <- function(expr_mat,
                         data_type           = "auto",
                         # Row label control
                         label_genes         = NULL,       # Character vector of gene names to label
                         # Sample & gene metadata
                         metadata_col        = NULL,       # data.frame with Sample_ID column for column (sample) annotations
                         metadata_row        = NULL,       # data.frame with SYMBOL column for row (gene) annotations
                         # Annotation column selection
                         col_annotations     = NULL,       # Character vector of metadata_col columns to show as annotations
                         row_annotations     = NULL,       # Character vector of metadata_row columns to show as annotations
                         # Gap control (visual separators between groups)
                         col_gap_by          = NULL,       # Single metadata_col column name to draw column gaps by
                         row_gap_by          = NULL,       # Single metadata_row column name to draw row gaps by
                         # Clustering & ordering
                         col_cluster_by      = "all",      
                         row_cluster_by      = "all",      # Same as col_cluster_by but for rows (genes)
                         # Plot aesthetics
                         plot_title          = NULL,       # Single string; subtitle "(X rows x Y cols)" auto-appended
                         heatmap_palette     = "rdbu",     # "viridis"
                         border_color        = NA,         # Cell border color; NA = no border
                         show_expr_legend    = TRUE){      # Show/hide the heatmap color scale legend
  
  # ---- 🧪 Prepare Matrix ----
  
  # Convert column names to valid R names
  colnames(expr_mat) <- make.names(colnames(expr_mat))
  
  # Remove rows with NA gene names — cannot be used for labeling or annotation matching
  na_rows  <- is.na(rownames(expr_mat))
  expr_mat <- expr_mat[!na_rows, , drop = FALSE]
  if (sum(na_rows) > 0) {
    log_info(sample = "", step = "plot_heatmap",
             msg = glue::glue("Removed {sum(na_rows)} rows with missing gene names."))
  }
  
  # Replace any NA values with 0 to prevent pheatmap errors
  expr_mat[is.na(expr_mat)] <- 0
  
  # Coerce matrix to numeric safely, preserving row/col names
  # This handles cases where input was passed as data.frame with mixed types
  expr_mat <- matrix(as.numeric(as.matrix(expr_mat)),
                     nrow     = nrow(expr_mat),
                     ncol     = ncol(expr_mat),
                     dimnames = dimnames(expr_mat))
  
  # After coercion, flag any values that failed to convert
  if (any(is.na(expr_mat))) {
    log_error(sample = "", step = "plot_heatmap",
              msg = "`expr_mat` contains non-numeric values that could not be converted.")
  }
  
  # Remove genes where ALL values are exactly 0 (no signal, uninformative)
  # abs() used to handle NES/LFC where positive and negative values 
  # could cancel each other out in rowSums, masking non-zero rows
  zero_genes <- rownames(expr_mat)[rowSums(abs(expr_mat)) == 0]
  expr_mat   <- expr_mat[rowSums(abs(expr_mat)) != 0, , drop = FALSE]
  if (length(zero_genes) > 0) {
    log_info(sample = "", step = "plot_heatmap",
             msg = glue::glue("Removed {length(zero_genes)} genes with zero counts across all samples."))
  }
  
  # Remove samples where ALL values are 0 (empty samples, likely QC failures)
  zero_samples <- colnames(expr_mat)[colSums(abs(expr_mat)) == 0]
  expr_mat     <- expr_mat[, colSums(abs(expr_mat)) != 0, drop = FALSE]
  if (length(zero_samples) > 0) {
    log_info(sample = "", step = "plot_heatmap",
             msg = glue::glue("Removed {length(zero_samples)} samples with zero total counts."))
  }
  
  # Guard: need at least 2 genes to draw a meaningful heatmap
  if (nrow(expr_mat) < 2) {
    log_warn(sample = "", step = "plot_heatmap",
             msg = "Input has fewer than 2 genes. Skipping heatmap generation.")
    return(NULL)
  }
  
  # Handle duplicate gene symbols — keep the row with the highest total absolute expression
  # abs() used so NES rows with large negative values are retained over near-zero rows
  if (any(duplicated(rownames(expr_mat)))) {
    log_warn(sample = "", step = "plot_heatmap",
             msg = "Duplicate gene symbols detected. Retaining highest expressing copy per gene.")
    expr_mat <- expr_mat %>%
      tibble::rownames_to_column("SYMBOL") %>%
      dplyr::mutate(total_expr = rowSums(dplyr::across(-SYMBOL, abs))) %>%
      dplyr::group_by(SYMBOL) %>%
      dplyr::slice_max(order_by = total_expr, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      tibble::column_to_rownames("SYMBOL") %>%
      dplyr::select(-total_expr)
  }
  
  # Sanitize row and column names again after deduplication
  rownames(expr_mat) <- make.names(rownames(expr_mat))
  colnames(expr_mat) <- make.names(colnames(expr_mat))
  
  # ---- 🔍 Data Type Detection & Validation ----
  
  # Compute summary statistics used for detection
  mat_min   <- min(expr_mat, na.rm = TRUE)
  mat_max   <- max(expr_mat, na.rm = TRUE)
  mat_range <- mat_max - mat_min
  
  # Auto-detect data type from matrix properties if user did not specify
  # Detection order matters — check most specific conditions first:
  #   1. bounded  : all values in [0,1] — catches beta, p-values, binary
  #   2. correlation: has negatives AND fully within [-1,1]
  #   3. centered : has negatives but range exceeds ±1 (NES, LFC)
  #   4. counts   : all positive, range > 1 (default fallback)
  if (data_type == "auto") {
    
    detected_type <- dplyr::case_when(
      mat_min >= 0  & mat_max <= 1    ~ "bounded",
      mat_min >= -1 & mat_max <= 1    ~ "correlation",
      mat_min <  0  & mat_range > 2   ~ "centered",
      mat_min >= 0  & mat_range > 100 ~ "counts",
      mat_min >= 0  & mat_range > 1   ~ "counts_log",
      TRUE                            ~ "counts")
    
    log_info(sample = "", step = "plot_heatmap",
             msg = glue::glue(
               "data_type = 'auto' → detected '{detected_type}'.
               If incorrect, override with data_type = 'counts', 'counts_log',
               'counts_scaled', 'centered', 'correlation', or 'bounded'."))
    
  } else {
    
    # User-specified data type — validate it is a known type
    valid_types <- c("counts", "counts_log", "counts_scaled", 
                     "centered", "correlation", "bounded")
    
    if (!data_type %in% valid_types) {
      log_error(sample = "", step = "plot_heatmap",
                msg = glue::glue(
                  "`data_type` '{data_type}' is not valid. 
                  Must be one of: {paste(valid_types, collapse = ', ')}."))
    }
    
    # Sanity check: warn if user-specified type conflicts with obvious data properties
    if (data_type == "counts" && mat_min < 0) {
      log_warn(sample = "", step = "plot_heatmap",
               msg = "data_type = 'counts' specified but matrix contains negative values.")
    }
    
    if (data_type == "correlation" && (mat_min < -1 || mat_max > 1)) {
      log_warn(sample = "", step = "plot_heatmap",
               msg = "data_type = 'correlation' specified but values fall outside [-1, 1].")
    }
    
    if (data_type == "bounded" && (mat_min < 0 || mat_max > 1)) {
      log_warn(sample = "", step = "plot_heatmap",
               msg = "data_type = 'bounded' specified but values fall outside [0, 1].")
    }
    
    detected_type <- data_type
    log_info(sample = "", step = "plot_heatmap",
             msg = glue::glue("data_type = '{detected_type}' (user-specified)."))
  }
  
  # ---- ⚖️ Log Transform & Z-score Scaling ----
  
  # Transform decisions by data type:
  #  data_type        | log2(1+x) | z-score
  #  counts           |    YES    |   YES     (raw counts: wide range, needs both)
  #  counts_log       |    no     |   YES     (VST/rlog: already logged, needs scaling)
  #  counts_scaled    |    no     |   no      (already fully processed)
  #  centered         |    no     |   no      (NES/LFC: already meaningful scale)
  #  correlation      |    no     |   no      (bounded [-1,1]: fixed scale)
  #  bounded          |    no     |   no      (bounded [0,1]: fixed scale)
  
  # Step 1: Log transform (counts only)
  if (detected_type == "counts") {
    log_info(sample = "", step = "plot_heatmap", 
             msg = "Applying log2(1+x) transform.")
    expr_mat <- log2(1 + expr_mat)
  } else {
    log_info(sample = "", step = "plot_heatmap", 
             msg = glue::glue("Skipping log transform for data_type = '{detected_type}'."))
  }
  
  # Step 2: Z-score scaling (counts and counts_log only)
  # Transpose → scale (operates on columns) → transpose back to original orientation
  # This z-scores each GENE across samples (row-wise z-score)
  if (detected_type %in% c("counts", "counts_log")) {
    log_info(sample = "", step = "plot_heatmap", 
             msg = "Applying row-wise z-score scaling.")
    expr_mat_scaled <- expr_mat %>% t() %>% scale() %>% t()
  } else {
    log_info(sample = "", step = "plot_heatmap", 
             msg = glue::glue("Skipping z-score for data_type = '{detected_type}'."))
    expr_mat_scaled <- expr_mat
  }
  
  # Replace any NAs introduced by scaling (e.g. zero-variance rows have undefined z-score)
  expr_mat_scaled[is.na(expr_mat_scaled)] <- 0
  
  # ---- 🎨 Color Breaks ----
  
  # Number of color bins in the palette
  n_breaks <- 100
  
  # Compute break range based on data type:
  #   bounded/binary  → fixed [0, 1]
  #   correlation     → fixed [-1, 1]  
  #   counts/centered → data-driven ±99th percentile (after any transforms)
  #                     99th percentile clips extreme outliers that would 
  #                     otherwise compress the color scale for all other values
  
  if (detected_type == "bounded") {
    
    # Fixed scale for all [0,1] data (beta values, p-values, binary)
    scale_min <- 0
    scale_max <- 1
    breaks    <- seq(from = 0, to = 1, length.out = n_breaks + 1)
    
  } else if (detected_type == "correlation") {
    
    # Fixed symmetric scale for correlation matrices
    scale_min <- -1
    scale_max <-  1
    breaks    <- seq(from = -1, to = 1, length.out = n_breaks + 1)
    
  } else {
    
    # Data-driven symmetric breaks for all remaining types:
    #   counts / counts_log → after z-scoring, distribution is centered at 0
    #   counts_scaled       → already centered at 0
    #   centered            → NES/LFC naturally centered at 0
    # All share identical break logic: symmetric around 0 at ±99th percentile
    # of absolute values, clipping outliers without hardcoding thresholds.
    # Two equal-length sequences meeting at 0 ensure zero is always a break
    # point — critical so the color midpoint maps exactly to 0
    p99       <- stats::quantile(abs(as.vector(expr_mat_scaled)), probs = 0.99, na.rm = TRUE)
    scale_min <- -p99
    scale_max <-  p99
    breaks    <- c(seq(from = scale_min, to = 0,         length.out = n_breaks / 2 + 1),
                   seq(from = 0,         to = scale_max, length.out = n_breaks / 2 + 1)[-1])
  }
  
  log_info(sample = "", step = "plot_heatmap",
           msg = glue::glue("Color scale range: [{round(scale_min, 2)}, {round(scale_max, 2)}]"))
  
  # ---- 🖌️ Heatmap Color Palette ----
  
  # Validate heatmap_palette argument
  if (!heatmap_palette %in% c("rdbu", "viridis")) {
    log_error(sample = "", step = "plot_heatmap",
              msg = "`heatmap_palette` must be either 'rdbu' or 'viridis'.")
  }
  
  # Build color vector with n_breaks colors
  # heatmap_palette accepts:
  #   "rdbu"              : Red-Blue diverging — best for centered/correlation data
  #                         where both positive (red) and negative (blue) have meaning
  #   "viridis"           : Viridis sequential — best for counts/bounded data
  #                         where values increase monotonically
  #   Any RColorBrewer    : e.g. "PuOr", "RdYlBu", "Blues" — see RColorBrewer::display.brewer.all()
  #   Any viridis variant : "magma", "plasma", "inferno", "cividis", "mako", "rocket", "turbo"
  #   A character vector  : Custom color vector e.g. c("#FFFFFF", "#FF0000") — interpolated to n_breaks
  
  heatmap_colors <- if (is.character(heatmap_palette) && length(heatmap_palette) > 1) {
    # User passed a custom color vector — interpolate to exactly n_breaks colors
    colorRampPalette(heatmap_palette)(n_breaks)
    
  } else if (heatmap_palette == "rdbu") {
    # Default: reversed RdBu (red = high, blue = low) — most intuitive for biology
    colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(n_breaks)
    
  } else if (heatmap_palette %in% c("viridis", "magma", "plasma", "inferno", "cividis", "mako", "rocket", "turbo")) {
    # Any viridis family palette
    viridis::viridis(n = n_breaks, option = heatmap_palette)
    
  } else if (heatmap_palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    # Any valid RColorBrewer palette name
    n_max <- RColorBrewer::brewer.pal.info[heatmap_palette, "maxcolors"]
    colorRampPalette(RColorBrewer::brewer.pal(n = n_max, name = heatmap_palette))(n_breaks)
    
  } else {
    log_error(sample = "", step = "plot_heatmap",
              msg = glue::glue("`heatmap_palette` '{heatmap_palette}' not recognized. ",
                               "Use 'rdbu', a viridis variant, any RColorBrewer palette name, ",
                               "or a custom character vector of colors."))
  }
  
  # ---- 🏷️ Prepare Annotations ----
  
  # Column (Sample) Annotation
  # Filter to only samples present in the matrix, select requested annotation columns
  # Rownames of col_annotation must match colnames of expr_mat for pheatmap alignment
  col_annotation <- if (!is.null(metadata_col) && !is.null(col_annotations)) {
    
    # Validate that requested annotation columns exist in metadata
    missing_col_anns <- setdiff(col_annotations, colnames(metadata_col))
    if (length(missing_col_anns) > 0) {
      log_error(sample = "", step = "plot_heatmap",
                msg = glue::glue("col_annotations missing in metadata_col: {paste(missing_col_anns, collapse = ', ')}"))
    }
    
    metadata_col %>%
      dplyr::select(Sample_ID, dplyr::all_of(col_annotations)) %>%
      dplyr::mutate(Sample_ID = make.names(Sample_ID)) %>%
      dplyr::filter(Sample_ID %in% colnames(expr_mat_scaled)) %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames("Sample_ID") %>%
      as.data.frame() %>%
      dplyr::mutate(dplyr::across(where(is.factor), as.character))
  } else NULL
  
  # Row (Gene) Annotation
  # Filter to only genes present in the matrix, select requested annotation columns
  # Rownames of row_annotation must match rownames of expr_mat for pheatmap alignment
  row_annotation <- if (!is.null(metadata_row) && !is.null(row_annotations)) {
    
    # Validate that requested annotation columns exist in metadata
    missing_row_anns <- setdiff(row_annotations, colnames(metadata_row))
    if (length(missing_row_anns) > 0) {
      log_error(sample = "", step = "plot_heatmap",
                msg = glue::glue("row_annotations missing in metadata_row: {paste(missing_row_anns, collapse = ', ')}"))
    }
    
    metadata_row %>%
      dplyr::select(SYMBOL, dplyr::all_of(row_annotations)) %>%
      dplyr::mutate(SYMBOL = make.names(SYMBOL)) %>%
      dplyr::filter(SYMBOL %in% rownames(expr_mat_scaled)) %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames("SYMBOL") %>%
      as.data.frame() %>%
      dplyr::mutate(dplyr::across(where(is.factor), as.character))
  } else NULL
  
  # ---- 🎨 Annotation Color Palette ----
  
  # EXAMPLE: Structure of ann_colors for reference
  # ann_colors defines the color for each level of each annotation variable.
  # Each list element is named after an annotation column, and contains
  # a named character vector mapping levels to colors:
  #
  # ann_colors <- list(
  #   CellType  = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
  #   GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E")
  # )
  
  # Use custom_palette if defined globally, otherwise fall back to a built-in
  # 25-color palette with perceptually distinct colors suitable for bioinformatics
  # annotation use cases (cell types, pathways, conditions, etc.)
  base_colors <- if (exists("custom_palette")) {
    custom_palette
  } else {
    log_warn(sample = "", step = "plot_heatmap",
             msg = "`custom_palette` not found in environment. Using built-in default palette.")
    c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
      "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
      "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
      "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
      "#66A61E", "#E6AB02", "#A6761D", "#666666", "#BC80BD")
  }
  
  # Build ann_colors by iterating over all annotation columns (row + col combined)
  # For each annotation variable:
  #   - Character/factor columns → discrete palette (one color per unique level)
  #   - Numeric columns         → continuous palette (gradient from light to saturated)
  # color_index tracks position in base_colors so each variable gets distinct colors
  
  ann_colors  <- list()
  color_index <- 1
  
  # Combine row and col annotation lists — process them uniformly
  col_list <- if (!is.null(col_annotation)) base::lapply(X = as.list(col_annotation), FUN = function(x) { as.character(x) %>% unique() }) else list()
  row_list <- if (!is.null(row_annotation)) base::lapply(X = as.list(row_annotation), FUN = function(x) { as.character(x) %>% unique() }) else list()
  ann_list <- c(col_list, row_list)
  # Since columns are almost always annotated, color columns first then rows, so that
  # color order is maintained across multiple plots in case rows are not annotated. 
  
  # Determine palette type per annotation column based on data type in metadata
  # numeric → continuous gradient; character/factor → discrete distinct colors
  col_types <- if (!is.null(col_annotation)) sapply(col_annotation, class) else c()
  row_types <- if (!is.null(row_annotation)) sapply(row_annotation, class) else c()
  ann_types <- c(row_types, col_types)  # named vector: annotation column → R class
  
  for (i in seq_along(ann_list)) {
    
    levels   <- sort(ann_list[[i]])   # Sorted levels within this annotation variable
    n_levels <- length(levels)        # Number of unique levels
    ann_name <- names(ann_list)[i]    # Annotation variable name (e.g. "CellType")
    is_numeric_ann <- ann_types[ann_name] %in% c("numeric", "integer", "double")
    
    # Assign colors based on annotation type:
    #   Discrete  : assign a unique base color to each level
    #   Continuous: assign alpha-graduated shades of a single base color
    palette_colors <- if (!is_numeric_ann || n_levels == 1) {
      # Discrete: map each level to a distinct color from base_colors
      base_colors[color_index:(color_index + n_levels - 1)]
    } else {
      # Continuous: vary alpha from light (low) to full opacity (high)
      # using a single base color — preserves hue while showing magnitude
      alphas <- seq(1 / n_levels, 1, length.out = n_levels)
      base::sapply(X   = alphas,
                   FUN = function(x) { colorspace::adjust_transparency(col = base_colors[color_index], alpha = x) })
    }
    
    names(palette_colors)  <- levels                    # Name each color with its level
    ann_colors             <- c(ann_colors, list(palette_colors))  # Append to master list
    names(ann_colors)[i]   <- ann_name                  # Name the list element
    color_index            <- color_index + n_levels    # Advance index past used colors
  }
  
  # ---- 🏁 Gap Positions ----
  
  # Gaps visually separate groups of columns/rows in the heatmap.
  # Gap positions are cumulative counts of samples/genes per group,
  # computed from sorted group labels (count() sorts alphabetically),
  # which must match the sort order used in col_order/row_order below.
  # The last gap (at total N) is excluded — pheatmap doesn't need it.
  
  # Column gaps
  gaps_col <- NULL
  if (!gtools::invalid(col_gap_by)) {
    if (col_gap_by %in% colnames(col_annotation)) {
      gaps_col <- col_annotation %>%
        dplyr::count(.data[[col_gap_by]]) %>%
        dplyr::mutate(n = cumsum(n)) %>%
        dplyr::pull(n) %>%
        .[. < ncol(expr_mat_scaled)]
    } else {
      log_warn(sample = "", step = "plot_heatmap",
               msg = glue::glue("col_gap_by '{col_gap_by}' not found in col_annotation columns."))
    }
  }
  
  # Row gaps
  gaps_row <- NULL
  if (!gtools::invalid(row_gap_by)) {
    if (row_gap_by %in% colnames(row_annotation)) {
      gaps_row <- row_annotation %>%
        dplyr::count(.data[[row_gap_by]]) %>%
        dplyr::mutate(n = cumsum(n)) %>%
        dplyr::pull(n) %>%
        .[. < nrow(expr_mat_scaled)]
    } else {
      log_warn(sample = "", step = "plot_heatmap",
               msg = glue::glue("row_gap_by '{row_gap_by}' not found in row_annotation columns."))
    }
  }
  
  
  # ---- 📊 Clustering & Ordering ----
  
  # Clustering strategy:
  #   "all"          → single hierarchical clustering across all columns/rows
  #   "none"         → preserve input order exactly as provided
  #   "alphabetical" → sort by name
  #   <column name>  → cluster within each group defined by that metadata column,
  #                    then concatenate groups in sorted order (must match gap order)
  #
  # Distance metric  : Euclidean (appropriate after z-scoring)
  # Linkage method   : Ward.D2 (minimizes within-cluster variance; requires Euclidean)
  # NOTE: If using correlation-based distances, switch to average/complete linkage
  #       as Ward.D2 assumes Euclidean geometry
  
  # --- Column (Sample) Ordering ---
  
  if (!is.null(col_annotation) && col_cluster_by %in% colnames(col_annotation)) {
    
    # Group-wise clustering: cluster within each group, concatenate groups
    col_order <- c()
    
    # Sort groups alphabetically — MUST match count() sort order used in gaps_col
    col_vars <- col_annotation %>%
      dplyr::pull(col_cluster_by) %>%
      unique() %>%
      sort()
    
    for (col_var in col_vars) {
      
      # Get samples in this group that are also present in the matrix
      samples <- col_annotation %>%
        tibble::rownames_to_column("Sample_ID") %>%
        dplyr::filter(.data[[col_cluster_by]] == col_var) %>%
        dplyr::pull(Sample_ID) %>%
        base::intersect(colnames(expr_mat_scaled))
      
      if (length(samples) > 1) {
        # Hierarchical clustering within this group
        temp_mat  <- expr_mat_scaled[, samples, drop = FALSE]
        col_dist  <- stats::dist(x = t(temp_mat), method = "euclidean")
        col_clust <- stats::hclust(d = col_dist, method = "ward.D2")
        col_order <- c(col_order, colnames(temp_mat)[col_clust$order])
      } else if (length(samples) == 1) {
        # Single sample — no clustering needed
        col_order <- c(col_order, samples)
      }
      # length == 0: group exists in metadata but not in matrix — silently skip
    }
    
  } else if (col_cluster_by == "all") {
    # Global clustering across all samples
    col_dist  <- stats::dist(x = t(expr_mat_scaled), method = "euclidean")
    col_clust <- stats::hclust(d = col_dist, method = "ward.D2")
    col_order <- colnames(expr_mat_scaled)[col_clust$order]
    
  } else if (col_cluster_by == "alphabetical") {
    col_order <- sort(colnames(expr_mat_scaled))
    
  } else if (col_cluster_by == "none") {
    # Preserve exact input order — no clustering or sorting
    col_order <- colnames(expr_mat_scaled)
    
  } else {
    log_warn(sample = "", step = "plot_heatmap",
             msg = glue::glue("col_cluster_by '{col_cluster_by}' not recognized. Preserving input order."))
    col_order <- colnames(expr_mat_scaled)
  }
  
  # --- Row (Gene) Ordering ---
  
  if (!is.null(row_annotation) && row_cluster_by %in% colnames(row_annotation)) {
    
    # Group-wise clustering: cluster within each group, concatenate groups
    row_order <- c()
    
    # Sort groups alphabetically — MUST match count() sort order used in gaps_row
    row_vars <- row_annotation %>%
      dplyr::pull(row_cluster_by) %>%
      unique() %>%
      sort()
    
    for (row_var in row_vars) {
      
      # Get genes in this group that are also present in the matrix
      genes <- row_annotation %>%
        tibble::rownames_to_column("SYMBOL") %>%
        dplyr::filter(.data[[row_cluster_by]] == row_var) %>%
        dplyr::pull(SYMBOL) %>%
        base::intersect(rownames(expr_mat_scaled))
      
      if (length(genes) > 1) {
        # Hierarchical clustering within this group
        temp_mat  <- expr_mat_scaled[genes, , drop = FALSE]
        row_dist  <- stats::dist(x = temp_mat, method = "euclidean")
        row_clust <- stats::hclust(d = row_dist, method = "ward.D2")
        row_order <- c(row_order, rownames(temp_mat)[row_clust$order])
      } else if (length(genes) == 1) {
        # Single gene — no clustering needed
        row_order <- c(row_order, genes)
      }
      # length == 0: group exists in metadata but not in matrix — silently skip
    }
    
  } else if (row_cluster_by == "all") {
    # Global clustering across all genes
    row_dist  <- stats::dist(x = expr_mat_scaled, method = "euclidean")
    row_clust <- stats::hclust(d = row_dist, method = "ward.D2")
    row_order <- rownames(expr_mat_scaled)[row_clust$order]
    
  } else if (row_cluster_by == "alphabetical") {
    row_order <- sort(rownames(expr_mat_scaled))
    
  } else if (row_cluster_by == "none") {
    # Preserve exact input order — no clustering or sorting
    row_order <- rownames(expr_mat_scaled)
    
  } else {
    log_warn(sample = "", step = "plot_heatmap",
             msg = glue::glue("row_cluster_by '{row_cluster_by}' not recognized. Preserving input order."))
    row_order <- rownames(expr_mat_scaled)
  }
  
  # Apply final row and column ordering to the scaled matrix
  reordered <- expr_mat_scaled[row_order, col_order]
  
  # ---- 🎨 Heatmap Aesthetics ----
  
  # Layout philosophy:
  #   Ideal font size for readability : 10pt
  #   Ideal cell size                 : fontsize + 5pt (gives breathing room)
  #   Available plot area on 8.5x11"  : ~6" wide x 8" tall (432pt x 576pt)
  #   Cell sizes are capped at ideal but shrink if many rows/cols don't fit
  #   NOTE: If using col_gaps, heatmap may get cut off — increase PDF width
  #         from 10" to 12"+ in your save_heatmaps() helper function
  
  fontsize        <- 10
  fontsize_number <- fontsize * 0.8   # Smaller font for in-cell numbers if shown
  angle_col       <- 45               # Column label rotation angle
  
  # Dynamically scale cell dimensions to fit available plot area
  # If many genes/samples: cells shrink below ideal and labels are hidden
  cell_width  <- min(fontsize + 5, (6 * 72) / ncol(reordered))
  cell_height <- min(fontsize + 5, (8 * 72) / nrow(reordered))
  
  # ---- Plot Title ----
  # Build a two-line title:
  #   Line 1: user-provided title (word-wrapped at 40 chars)
  #   Line 2: auto-generated subtitle showing matrix dimensions
  # pheatmap requires main != NULL, so use NA when no title provided
  
  n_rows     <- nrow(reordered)
  n_cols     <- ncol(reordered)
  subtitle   <- glue::glue("({n_rows} rows \u00d7 {n_cols} columns)")  # × unicode symbol
  
  plot_title <- if (!is.null(plot_title) && is.character(plot_title) && length(plot_title) == 1) {
    # Combine user title + subtitle on separate lines
    # pheatmap renders \n as a line break in the title
    paste0(stringr::str_wrap(string = plot_title, width = 40), "\n",
           subtitle)
  } else {
    # No user title — show only the dimension subtitle
    subtitle
  }
  
  # ---- Row Labels ----
  # Show row labels only when cells are large enough to be readable
  # If label_genes provided: show only those genes (blank all others)
  # If no label_genes: show all labels only when cell height is sufficient
  
  labels_row <- if (!is.null(label_genes)) {
    # Show only specified genes — blank out all others
    dplyr::if_else(condition = rownames(reordered) %in% make.names(label_genes),
                   true      = stringr::str_trunc(string = rownames(reordered), width = 15),
                   false     = " ")
  } else if (cell_height == fontsize + 5) {
    # Cells are large enough — show all row labels, truncated for length
    stringr::str_trunc(string = rownames(reordered), width = 15)
  } else {
    # Cells too small — hide all row labels to avoid overplotting
    rep(x = " ", times = nrow(reordered))
  }
  
  # ---- Column Labels ----
  # Show column labels only when cells are large enough to be readable
  
  labels_col <- if (cell_width == fontsize + 5) {
    stringr::str_trunc(string = colnames(reordered), width = 15)
  } else {
    rep(x = " ", times = ncol(reordered))
  }
  
  # ---- 🖼️ Generate Heatmap ----
  
  ph <- pheatmap::pheatmap(
    
    # Core matrix & colors
    mat               = reordered,
    color             = heatmap_colors,
    breaks            = breaks,
    
    # Annotations
    annotation_row    = row_annotation,
    annotation_col    = col_annotation,
    annotation_colors = ann_colors,
    
    # Gaps between groups
    gaps_row          = gaps_row,
    gaps_col          = gaps_col,
    
    # Cell sizing & label display
    cellwidth         = cell_width,
    cellheight        = cell_height,
    show_rownames     = cell_height >= fontsize,
    show_colnames     = cell_width  >= fontsize,
    labels_row        = labels_row,
    labels_col        = labels_col,
    angle_col         = angle_col,
    
    # Font sizes
    fontsize          = fontsize,
    fontsize_row      = fontsize,
    fontsize_col      = fontsize,
    fontsize_number   = fontsize_number,
    
    # Plot metadata
    main              = plot_title,
    border_color      = border_color,
    legend            = show_expr_legend,
    
    # Clustering — disabled here because we handle ordering manually above
    # All row/col ordering is done pre-pheatmap via hierarchical clustering
    # or user-specified ordering logic; pheatmap just renders the result
    scale                    = "none",
    cluster_rows             = FALSE,
    cluster_cols             = FALSE,
    annotation_legend        = TRUE,
    annotation_names_row     = FALSE,
    annotation_names_col     = FALSE,
    
    # Dimensions managed externally (in save_heatmaps helper)
    width                    = NA,
    height                   = NA,
    filename                 = NA,  # Never save from within this function
    silent                   = TRUE)
  
  # ---- 📋 Prepare Output Matrix ----
  
  # The output matrix reflects the final reordered, transformed state
  # (same as what is displayed in the heatmap)
  ph_mat <- reordered
  
  # Excel has hard limits: 1,048,576 rows x 16,384 columns
  # Transpose if needed to fit (e.g. many samples, few genes)
  max_rows <- 1048576
  max_cols <- 16384
  
  if (ncol(reordered) > max_cols && ncol(reordered) <= max_rows && nrow(reordered) <= max_cols) {
    log_info(sample = "", step = "plot_heatmap",
             msg = "Matrix transposed to fit within Excel column limit.")
    ph_mat <- t(reordered)
  } else if (nrow(reordered) > max_rows) {
    log_warn(sample = "", step = "plot_heatmap",
             msg = "Matrix exceeds Excel row limit even after transpose. Save as CSV instead.")
  }
  
  # ---- 🪵 Log & Return ----
  
  log_info(sample = "", step = "plot_heatmap",
           msg = glue::glue("Heatmap generated: {nrow(reordered)} rows x {ncol(reordered)} cols | ",
                            "data_type = '{detected_type}' | palette = '{heatmap_palette}'"))
  
  # Return invisibly so assignment is optional but all outputs are accessible:
  #   result$ph        → pheatmap object (pass to grid.draw() for rendering/saving)
  #   result$mat       → final matrix (save to Excel or use for downstream analysis)
  #   result$data_type → detected/used data type (useful for sanity checking upstream)
  return(invisible(list(ph        = ph,
                        mat       = ph_mat,
                        data_type = detected_type)))
}