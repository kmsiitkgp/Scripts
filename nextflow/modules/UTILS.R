suppressPackageStartupMessages({
  library(crayon)  # Colored terminal output
  library(glue)    # String interpolation
  library(rlang)   # %||% null-coalescing operator
  library(readr)
  library(readxl)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(scales)
  library(survival)
  library(survminer)
  library(emmeans)
  library(openxlsx)
  library(cowplot)
  library(stringr)    # likely already installed, used for str_to_title()
  library(purrr)      # used in global_cutoff = FALSE path
})

# Force crayon color codes active even when output is redirected to a file
# (the default in Nextflow). Without this, crayon detects a non-terminal
# destination and strips all color codes, making log output monochrome.
options(crayon.enabled = TRUE)

# ═══════════════════════════════════════════════════════════════════════════════
# LOGGING FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

# quiet_msg() — suppress all stdout and stderr from an expression.
# Why tempfile + sink rather than suppressMessages() / suppressWarnings()?
# suppressMessages() only silences message() calls; it does not capture
# cat(), print(), or C-level output from compiled packages. sink() redirects
# ALL output at the R connection level, which is the only reliable way to
# fully silence noisy packages like AnnotationHub or tximport that mix
# message() and cat() output. The tempfile is used as a discard sink —
# output is written there and the file is deleted on exit.
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
  
  result <- force(expr)
  return(result)
}

log_info <- function(sample, step, msg) {
  sample <- sample %||% ""
  prefix <- green(formatC("[INFO]", width = 7, flag = " "))
  message(glue::glue("{prefix} [{sample} | {toupper(step)}] {msg}"))
}

log_warn <- function(sample, step, msg) {
  sample <- sample %||% ""
  prefix <- yellow(formatC("[WARN]", width = 7, flag = " "))
  message(glue::glue("{prefix} [{sample} | {toupper(step)}] {msg}"))
}

# log_error() logs a red [ERROR] message then immediately halts execution.
# Why stop() with call. = FALSE?
# call. = TRUE (the default) appends the internal R call stack to the error
# message — in a pipeline context this produces confusing output like
# "Error in log_error(...)" pointing to UTILS.R rather than the calling script.
# call. = FALSE produces a clean "Workflow stopped." message with no stack trace,
# which is more useful in Nextflow logs where the log_error() message itself
# already contains all the context needed to diagnose the problem.
log_error <- function(sample, step, msg) {
  sample <- sample %||% ""
  prefix <- red(formatC("[ERROR]", width = 7, flag = " "))
  message(glue::glue("{prefix} [{sample} | {toupper(step)}] {msg}"))
  stop("Workflow stopped.", call. = FALSE)
}

log_sample_header <- function(sample) {
  cat(blue$bold(glue::glue("\n--- Processing Sample: {sample} ---\n\n")))
}

log_section <- function(section_name) {
  cat(magenta$bold(glue::glue("\n[{toupper(section_name)}]\n")))
}

# ═══════════════════════════════════════════════════════════════════════════════
# SMART FILE LOADER
# ═══════════════════════════════════════════════════════════════════════════════
# load_smart() accepts either a file path string or an already-loaded R object.
# If given an R object (data.frame, matrix, list, DESeqResults etc.), it returns
# it unchanged — this allows all pipeline functions to call load_smart() on every
# input without needing to know whether the caller passed a path or an object.
#
# Why the !is.character || length != 1 guard?
# Any non-string input (NULL, data.frame, numeric vector) would cause
# file.exists() to error or return FALSE spuriously. The guard short-circuits
# before any file operations for non-path inputs. length != 1 handles character
# vectors with multiple elements — those are not file paths.
#
# sheet is only used for Excel files; defaults to the first sheet if NULL.

load_smart <- function(input_path, sheet = NULL) {
  
  if (!is.character(input_path) || length(input_path) != 1) {
    return(input_path)
  }
  
  if (!file.exists(input_path)) {
    log_error(sample = "", step = "load_smart",
              msg = glue::glue("File does not exist: '{input_path}'"))
  }
  
  ext <- tools::file_ext(tolower(input_path))
  
  data <- switch(EXPR = ext,
                 "csv"  = readr::read_csv(input_path, show_col_types = FALSE),
                 "rds"  = readRDS(input_path),
                 "txt"  = readr::read_tsv(input_path, show_col_types = FALSE),
                 "xlsx" = ,
                 "xls"  = {
                   if (!is.null(sheet)) {
                     readxl::read_excel(input_path, sheet = sheet)
                   } else {
                     readxl::read_excel(input_path)
                   }
                 },
                 log_error(sample = "", step = "load_smart",
                           msg = glue::glue("Unsupported file extension: '{ext}'"))
  )
  
  return(data)
}

# ═══════════════════════════════════════════════════════════════════════════════
# ANNOTATION FUNCTION
# ═══════════════════════════════════════════════════════════════════════════════
# add_annotation() joins gene_name and gene_biotype from tx2gene onto any
# count matrix or DEG dataframe, then relocates gene_id, SYMBOL, gene_biotype
# to the front columns.
#
# Why auto-detect the gene_id column rather than requiring the user to specify it?
# Count matrices and DEG tables come from different tools (STAR, Salmon, DESeq2)
# and store Ensembl IDs under different column names (gene_id, GeneID, rownames).
# Auto-detection by overlap makes add_annotation() robust to any upstream format
# without the caller needing to know or pass the column name.
#
# The overlap threshold of >= 0.5 (50%) is a pragmatic minimum — below this, the
# detected column likely does not contain Ensembl IDs and the join would produce
# mostly NAs. A log_error fires below 50%; a log_warn fires above it to confirm
# the match (useful for quickly spotting partial ID version mismatches in logs).

add_annotation <- function(df, tx2gene) {
  
  # ── Step 1: Clean annotation map ─────────────────────────────────────────
  # Select by position (columns 2, 3, 4) so add_annotation() works regardless
  # of whether tx2gene has transcript_id in column 1 (full tx2gene) or starts
  # directly with gene_id (gene-level annotation table).
  # distinct(gene_id) keeps one row per gene — multiple transcripts per gene
  # would otherwise produce duplicate rows after the left_join below.
  tx2gene <- tx2gene %>%
    dplyr::select(gene_id = 2, gene_name = 3, gene_biotype = 4) %>%
    dplyr::distinct(gene_id, .keep_all = TRUE) %>%
    dplyr::mutate(gene_name = dplyr::case_when(
      (is.na(gene_name) | gene_name == "") & (gene_biotype == "protein_coding") ~ gene_id,
      TRUE ~ gene_name
    ))
  
  # ── Step 2: Auto-detect gene ID column ───────────────────────────────────
  # rownames_to_column ensures rownames (common in count matrices) are included
  # in the overlap search alongside regular columns.
  test_df         <- df %>% tibble::rownames_to_column("row_id")
  overlap_counts  <- sapply(test_df, function(x) sum(as.character(x) %in% tx2gene$gene_id))
  max_overlap_col <- names(overlap_counts)[which.max(overlap_counts)]
  overlap_pct     <- max(overlap_counts) / nrow(df)
  pct_label       <- round(overlap_pct * 100, 1)
  
  if (overlap_pct >= 0.5) {
    log_info(sample = "", step = "add_annotation",
             msg = glue::glue("Gene ID column detected: '{max_overlap_col}' ({pct_label}% overlap with tx2gene)"))
  } else {
    log_error(sample = "", step = "add_annotation",
              msg = glue::glue("Annotation failed — low ID overlap ({pct_label}% in '{max_overlap_col}'). ",
                               "Check that tx2gene matches the genome used for alignment."))
  }
  
  # ── Step 3: Join and reformat ─────────────────────────────────────────────
  annotated_df <- test_df %>%
    dplyr::rename(gene_id = .data[[max_overlap_col]]) %>%
    dplyr::select(-dplyr::any_of("row_id")) %>%
    dplyr::mutate(gene_id = as.character(gene_id)) %>%
    dplyr::left_join(tx2gene, by = "gene_id") %>%
    dplyr::relocate(gene_id, gene_name, gene_biotype) %>%
    dplyr::rename(SYMBOL = gene_name)
  
  return(annotated_df)
}

# ═══════════════════════════════════════════════════════════════════════════════
# COUNT AND DEG FORMATTING
# ═══════════════════════════════════════════════════════════════════════════════

# process_counts() converts a raw expression matrix (genes × samples) to an
# annotated dataframe with SYMBOL, gene_biotype, and sample columns.
# When multiple Ensembl IDs map to the same gene symbol, keeps the one with
# the highest average expression across samples.
#
# Why compute rowMeans BEFORE annotation (not after)?
# After add_annotation() the matrix is a dataframe with extra character columns
# (gene_id, SYMBOL, gene_biotype). rowMeans on a mixed-type dataframe would
# either error or silently drop character columns. Computing on the pure numeric
# matrix first avoids this and ensures the mean is based on the original values.
process_counts <- function(expr_obj, tx2gene) {
  
  row_means               <- rowMeans(as.matrix(expr_obj), na.rm = TRUE)
  row_means[is.na(row_means)] <- 0
  
  df <- as.data.frame(expr_obj) %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::mutate(row_mean = row_means)
  
  ann_df <- df %>%
    add_annotation(tx2gene) %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::slice_max(order_by = row_mean, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(-row_mean)
  
  return(ann_df)
}

# process_degs() converts a DESeq2 results object to an annotated dataframe.
# When multiple Ensembl IDs map to the same symbol, keeps the most significant
# (lowest padj) to avoid duplicates in downstream volcano plots and pathway tools.
#
# Why replace NA padj with 1?
# DESeq2 sets padj = NA for genes that failed independent filtering or Cook's
# distance cutoff. Downstream tools (volcano plots, -log10 transforms) cannot
# handle NA — replacing with 1 (= not significant) is the correct neutral value
# that won't affect significance thresholds or distort visualisations.
#
# Why replace padj == 0 with min non-zero padj?
# Extremely significant genes can have padj rounded to exactly 0 due to floating
# point underflow. -log10(0) = Inf, which breaks axis scaling in volcano plots.
# The minimum observed non-zero padj is the most conservative finite replacement.
process_degs <- function(res_obj, tx2gene) {
  
  df       <- as.data.frame(res_obj) %>% tibble::rownames_to_column("gene_id")
  min_padj <- min(df$padj[df$padj > 0], na.rm = TRUE)
  
  ann_df <- df %>%
    dplyr::mutate(padj = dplyr::case_when(
      is.na(padj) ~ 1,
      padj == 0   ~ min_padj,
      TRUE        ~ padj
    )) %>%
    add_annotation(tx2gene) %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::slice_min(order_by = padj, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  return(ann_df)
}

# ═══════════════════════════════════════════════════════════════════════════════
# CUSTOM COLOR PALETTE
# ═══════════════════════════════════════════════════════════════════════════════
# CUSTOM_PALETTE is the standard discrete color palette used across all plot
# functions. Colors are ordered for maximum perceptual distance between
# consecutive positions — so that the first N clusters/groups assigned are
# always maximally distinct regardless of N.
#
# Usage:
#   scale_color_manual(values = CUSTOM_PALETTE)   # points/lines
#   scale_fill_manual(values  = CUSTOM_PALETTE)   # bars/tiles
#   scale_color_manual(values = setNames(CUSTOM_PALETTE, levels(df$group)))

CUSTOM_PALETTE <- c(
  "#00538A",  # deep blue
  "#C10020",  # deep red
  "#007D34",  # dark green
  "#FFB300",  # gold
  "#6A0DAD",  # violet
  "#FF6800",  # orange
  "#1B9E77",  # emerald
  "#E7298A",  # magenta
  "#4169E1",  # royal blue
  "#7F180D",  # dark red
  "#2CA02C",  # green
  "#F28E2B",  # warm orange
  "#803E75",  # purple
  "#93AA00",  # lime
  "#FF7A5C",  # salmon
  "#17BECF",  # sky blue
  "#B32851",  # berry
  "#556B2F",  # dark olive
  "#9467BD",  # lavender
  "#CEA262",  # tan
  "#008080",  # teal
  "#F6768E",  # pink
  "#A65628",  # sienna
  "#EDC948",  # yellow
  "#53377A",  # dark violet
  "#B8860B",  # dark goldenrod
  "#E377C2",  # light pink
  "#2F4F4F",  # dark slate
  "#B07AA1",  # mauve
  "#817066",  # gray brown
  "#A0CBE8"   # pale blue
)

# ═══════════════════════════════════════════════════════════════════════════════
# PUBLICATION THEME
# ═══════════════════════════════════════════════════════════════════════════════
# theme_publication() is the universal ggplot2 theme for all pipeline plots.
# Apply once per script via ggplot2::theme_set(theme_publication()) before any
# plotting calls — every ggplot() in that script inherits it automatically.
#
# base_size   : base font size in pt; all other sizes scale via rel() from this.
#               12pt is appropriate for full-page figures. Use 10pt for figures
#               that will be reduced to single journal column width (~3.5").
# base_family : font family. "Helvetica" is standard for Nature / Cell / Science.
#               On Linux HPC nodes where Helvetica may not be installed, pass
#               "sans" for a visually identical fallback (Liberation Sans).
#
# theme_no_legend() is a convenience add-on for multi-panel figures where
# individual panel legends are redundant. Apply on top of the base theme:
#   p + theme_no_legend()
# or after theme_set():
#   ggplot2::theme_set(theme_publication())
#   p + theme_no_legend()

theme_publication <- function(base_size   = 12,
                              base_family = "Helvetica") {
  
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      
      # ── Lines (global) ────────────────────────────────────────────────────
      # lineend = "square" gives sharper axis lines and tick marks in vector
      # PDF output — round lineends look slightly blunt after reduction to
      # print size. Cascades to all line elements via the top-level override.
      line = element_line(lineend = "square"),
      
      # ── Axis lines and ticks ──────────────────────────────────────────────
      # linewidth replaces deprecated size for line elements (ggplot2 >= 3.4.0).
      # Ticks slightly thinner than axis lines — subordinate but clearly visible.
      axis.line         = element_line(linewidth = 0.8, color = "black"),
      axis.ticks        = element_line(linewidth = 0.6, color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      
      # ── Axis text (tick labels) ───────────────────────────────────────────
      # grey10 rather than pure black — imperceptibly softer, reduces harshness
      # in print without being obviously grey. Margins add breathing room between
      # tick marks and labels, preventing labels from hugging the axis.
      axis.text   = element_text(color  = "grey10",
                                 size   = rel(0.95),
                                 family = base_family),
      axis.text.x = element_text(margin = margin(t = 4)),
      axis.text.y = element_text(margin = margin(r = 6), hjust = 1),
      
      # ── Axis titles ───────────────────────────────────────────────────────
      # Bold and slightly larger than tick labels — immediately identifiable as
      # axis labels even after figure reduction to journal column width.
      # Margins separate titles from tick labels without excessive whitespace.
      axis.title   = element_text(face   = "bold",
                                  size   = rel(1.1),
                                  family = base_family),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      
      # ── Panel ─────────────────────────────────────────────────────────────
      # No gridlines — scatter plots use reference lines (vline/hline) for
      # guidance; bar and dot plots are readable from bar/dot positions alone.
      # panel.border = element_blank() is redundant after theme_classic() but
      # kept explicit so intent is clear if base theme ever changes.
      # panel.spacing controls gap between facet panels — 1 line prevents the
      # cramped appearance that the default 0.5 lines produces in multi-panel
      # pathway and TF plots.
      panel.background = element_rect(fill  = "white", color = NA),
      panel.border     = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing    = unit(1, "lines"),
      
      # ── Plot background ───────────────────────────────────────────────────
      # plot.background ensures white fill on export via cairo_pdf / ggsave
      # even when the graphics device has a non-white default.
      # plot.title.position = "plot" aligns title to the full figure width
      # rather than just the panel area — prevents indented titles on plots
      # with wide y-axis labels (e.g. pathway bar plots).
      plot.background     = element_rect(fill = "white", color = NA),
      plot.title.position = "plot",
      
      # ── Plot titles ───────────────────────────────────────────────────────
      # Typography hierarchy: title (bold, largest) > subtitle (regular, muted)
      # > caption (small, grey). Centred alignment (hjust = 0.5) for
      # presentation / poster style — change to 0 for strict journal style.
      # Margins below title and subtitle prevent crowding against the panel.
      plot.title    = element_text(face   = "bold",
                                   size   = rel(1.3),
                                   hjust  = 0.5,
                                   color  = "black",
                                   family = base_family,
                                   margin = margin(b = 6)),
      plot.subtitle = element_text(size   = rel(1.0),
                                   hjust  = 0.5,
                                   color  = "grey30",
                                   family = base_family,
                                   margin = margin(b = 8)),
      plot.caption  = element_text(size   = rel(0.8),
                                   hjust  = 1,
                                   color  = "grey50",
                                   family = base_family),
      
      # ── Legend ────────────────────────────────────────────────────────────
      # Always right-side vertical — consistent across all plot types.
      # key.height / key.width at 0.9 lines gives enough space for size-scaled
      # dot plot legends without oversized boxes on simpler color legends.
      # legend.margin left offset (8pt) adds breathing room between the panel
      # edge and the legend without requiring plot.margin changes.
      legend.position    = "right",
      legend.direction   = "vertical",
      legend.title       = element_text(face   = "bold",
                                        size   = rel(1.0),
                                        family = base_family),
      legend.text        = element_text(size   = rel(0.9),
                                        color  = "grey10",
                                        family = base_family),
      legend.background  = element_rect(fill = "transparent", color = NA),
      legend.key         = element_rect(fill = "transparent", color = NA),
      legend.key.height  = unit(0.9, "lines"),
      legend.key.width   = unit(0.9, "lines"),
      legend.box.spacing = unit(0.3, "lines"),
      legend.margin      = margin(0, 0, 0, 8),
      
      # ── Facet strips ──────────────────────────────────────────────────────
      # Internal margin (4pt each side) gives strip labels breathing room
      # inside the strip box — default has almost no padding.
      strip.background = element_rect(fill      = "grey95",
                                      color     = "black",
                                      linewidth = 0.8),
      strip.text       = element_text(face   = "bold",
                                      size   = rel(0.95),
                                      family = base_family,
                                      margin = margin(4, 4, 4, 4)),
      
      # ── Margins ───────────────────────────────────────────────────────────
      # Asymmetric — right margin (15pt) slightly wider than left (10pt) to
      # balance the visual weight of the right-side legend.
      plot.margin = margin(10, 15, 10, 10)
    )
}

theme_no_legend <- function() {
  theme(legend.position = "none")
}

# ═══════════════════════════════════════════════════════════════════════════════
# DATA TYPE DETECTION HELPER
# ═══════════════════════════════════════════════════════════════════════════════
# detect_data_type() resolves the data_type argument used by plot_heatmap() and
# plot_pca(). Centralising detection here ensures both functions use identical
# thresholds, log the same messages, and are updated in a single place.
#
# INPUT:
#   mat       : numeric matrix already subsetted to the samples/genes that will
#               be plotted. Detection is run on the actual data, not a copy.
#   data_type : user-supplied string — "auto" or one of the valid explicit types.
#
# RETURNS: a single string — the resolved data type.
#
# Detection order for "auto" (most specific first):
#   1. bounded     : all values in [0, 1]  (beta values, p-values, proportions)
#   2. correlation : has negatives AND fully within [-1, 1]
#   3. centered    : has negatives but range exceeds ±1  (LFC, NES, z-scores)
#   4. counts      : all positive, wide range > 100  (raw integer counts)
#   5. counts_log  : all positive, modest range > 1  (VST, rlog, log2 CPM)
#   6. counts      : fallback for anything else positive
#
# Sanity-check warnings fire when the user specifies an explicit type but the
# matrix values are inconsistent with it — e.g. "counts" with negatives.

detect_data_type <- function(mat, data_type) {
  
  valid_types <- c("auto", "counts", "counts_log", "counts_scaled",
                   "centered", "bounded")
  
  if (!data_type %in% valid_types) {
    log_error(sample = "", step = "detect_data_type",
              msg = glue::glue(
                "'data_type' '{data_type}' is not valid. ",
                "Must be one of: {paste(valid_types, collapse = ', ')}."
              ))
  }
  
  mat_min   <- min(mat, na.rm = TRUE)
  mat_max   <- max(mat, na.rm = TRUE)
  mat_range <- mat_max - mat_min
  
  if (data_type == "auto") {
    
    detected_type <- dplyr::case_when(
      mat_min >= 0  & mat_max <= 1    ~ "bounded",
      mat_min >= -1 & mat_max <= 1    ~ "correlation",
      mat_min <  0  & mat_range > 2   ~ "centered",
      mat_min >= 0  & mat_range > 100 ~ "counts",
      mat_min >= 0  & mat_range > 1   ~ "counts_log",
      TRUE                            ~ "counts"
    )
    
    log_info(sample = "", step = "detect_data_type",
             msg = glue::glue(
               "data_type = 'auto' -> detected '{detected_type}'. ",
               "Override with an explicit data_type if incorrect."
             ))
    
  } else {
    
    # Sanity-check user-specified type against actual matrix values — warn but
    # do not error, as the user may have a legitimate reason to override.
    if (data_type == "counts" && mat_min < 0)
      log_warn(sample = "", step = "detect_data_type",
               msg = "data_type = 'counts' specified but matrix contains negative values.")
    
    if (data_type %in% c("counts", "counts_log") && mat_min < 0)
      log_warn(sample = "", step = "detect_data_type",
               msg = glue::glue(
                 "data_type = '{data_type}' specified but matrix contains negative values."
               ))
    
    if (data_type == "correlation" && (mat_min < -1 || mat_max > 1))
      log_warn(sample = "", step = "detect_data_type",
               msg = "data_type = 'correlation' specified but values fall outside [-1, 1].")
    
    if (data_type == "bounded" && (mat_min < 0 || mat_max > 1))
      log_warn(sample = "", step = "detect_data_type",
               msg = "data_type = 'bounded' specified but values fall outside [0, 1].")
    
    detected_type <- data_type
    log_info(sample = "", step = "detect_data_type",
             msg = glue::glue("data_type = '{detected_type}' (user-specified)."))
  }
  
  return(detected_type)
}

# ═══════════════════════════════════════════════════════════════════════════════
# DISPERSION PLOT
# ═══════════════════════════════════════════════════════════════════════════════
# plot_dispersion() visualises per-gene dispersion estimates from a fitted
# DESeqDataSet using DESeq2::plotDispEsts().
#
# INPUT:
#   dds       : DESeqDataSet that has already passed through DESeq() or
#               estimateDispersions(). plotDispEsts() will error if dispersion
#               slots are empty.
#   filename  : base name for the output PDF. Non-alphanumeric characters are
#               sanitised to underscores. Default "Dispersion_Plot".
#   output_dir: directory to save the PDF. NULL = skip save, return NULL.
#
# plotDispEsts() is a base-graphics function — it draws directly to the active
# device and returns NULL invisibly. It cannot be captured as a ggplot object.
# The device must be opened before the call and closed immediately after.
#
# Colour convention matches DESeq2 defaults:
#   genecol  = "black"      gene-wise MAP estimates
#   fitcol   = "red"        mean-dispersion trend fit
#   finalcol = "dodgerblue" final shrunken estimates used for testing
#
# Expected pattern: dispersion decreases as mean normalised counts increase.
# Final estimates (blue) should cluster tightly around the fitted trend (red).
# Genes far above the trend are dispersion outliers — excluded from testing
# by DESeq2, normal for a small fraction of genes.
#
# OUTPUT:
#   Returns invisible(NULL) — no ggplot object to return since plotDispEsts()
#   is base graphics. Primary output is the saved PDF when output_dir is provided.

plot_dispersion <- function(dds, filename = "Dispersion_Plot", output_dir = NULL) {
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Validate Inputs
  # ═══════════════════════════════════════════════════════════════════════════
  # dds must be a DESeqDataSet that has already passed through DESeq() or
  # estimateDispersions() — plotDispEsts() will error if dispersion slots are
  # empty. Accepting any other class would produce a cryptic internal error
  # rather than a clear diagnostic message.
  
  if (!inherits(dds, "DESeqDataSet")) {
    log_error(sample = "", step = "plot_dispersion",
              msg = glue::glue("'dds' must be a DESeqDataSet, not {class(dds)[1]}."))
  }
  
  if (!is.null(output_dir) && !is.character(output_dir)) {
    log_error(sample = "", step = "plot_dispersion",
              msg = glue::glue("'output_dir' must be a character string or NULL, not {class(output_dir)[1]}."))
  }
  
  if (!is.character(filename)) {
    log_error(sample = "", step = "plot_dispersion",
              msg = glue::glue("'filename' must be a character string, not {class(filename)[1]}."))
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Save and Return
  # ═══════════════════════════════════════════════════════════════════════════
  
  if (!is.null(output_dir)) {
    
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    # Sanitise filename — replace non-alphanumeric runs with underscores and
    # strip any leading/trailing underscores introduced by the substitution
    safe_filename <- gsub("[^[:alnum:]]+", "_", filename)
    safe_filename <- gsub("(^_+|_+$)", "", safe_filename)
    
    output_file <- file.path(output_dir, paste0(safe_filename, ".pdf"))
    
    # cairo_pdf is used over ggsave because it natively supports multi-page
    # output and handles font embedding and transparency correctly.
    grDevices::cairo_pdf(filename = output_file,
                         width    = 8,
                         height   = 11.5,
                         onefile  = TRUE)
    
    DESeq2::plotDispEsts(object   = dds,
                         genecol  = "black",
                         fitcol   = "#E64B35",
                         finalcol = "#4DBBD5",
                         legend   = TRUE,
                         xlab     = "Mean of Normalized Counts",
                         ylab     = "Dispersion",
                         log      = "xy",
                         cex      = 0.45)
    
    grDevices::dev.off()
    
    log_info(sample = "", step = "plot_dispersion",
             msg = glue::glue("Dispersion plot saved to: '{output_file}'"))
    
  } else {
    
    log_info(sample = "", step = "plot_dispersion",
             msg = "No output_dir provided — skipping save.")
  }
  
  return(invisible(NULL))
}

# ═══════════════════════════════════════════════════════════════════════════════
# MA PLOT
# ═══════════════════════════════════════════════════════════════════════════════
# plot_ma() generates a two-page MA plot PDF from a DESeq2 results object or
# equivalent data.frame: page 1 uncapped (full dynamic range), page 2 capped
# at ±5 LFC with outliers shown as triangles.
#
# INPUT:
#   res      : DESeqResults object or data.frame with columns:
#                baseMean, log2FoldChange, padj
#              DESeqDataSet is intentionally rejected — run DESeq2::results()
#              first. Accepting dds would silently use default contrast and
#              lfcShrink choices, hiding them from the user.
#   filename : base name for the output PDF. Non-alphanumeric characters are
#              sanitised to underscores. Default "MA_Plot".
#   output_dir: directory to save the PDF. NULL = return plots only, no file saved.
#
# OUTPUT:
#   Returns invisibly: list(
#     uncapped = ggplot (full LFC range),
#     capped   = ggplot (LFC capped at ±5, outliers as triangles)
#   )
#   When output_dir is provided, saves a two-page PDF:
#     Page 1 : uncapped MA plot
#     Page 2 : capped MA plot (±5 LFC)

plot_ma <- function(res, filename = "MA_Plot", output_dir = NULL) {
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Validate Inputs
  # ═══════════════════════════════════════════════════════════════════════════
  # Accepts DESeqResults or data.frame only — DESeqDataSet is intentionally
  # excluded. plot_ma() sits downstream of DESeq2::results(); the caller is
  # responsible for extracting results before plotting. Accepting dds would
  # silently call results() with default arguments, hiding contrast and
  # lfcShrink choices from the user and producing plots that may not match
  # the intended comparison.
  
  if (inherits(res, "DESeqDataSet")) {
    log_error(sample = "", step = "plot_ma",
              msg = paste("'res' is a DESeqDataSet — pass a DESeqResults object instead.",
                          "Run DESeq2::results(dds, contrast = ...) first."))
  }
  
  res_df <- if (inherits(res, "DESeqResults")) {
    as.data.frame(res)
  } else if (is.data.frame(res)) {
    res
  } else {
    log_error(sample = "", step = "plot_ma",
              msg = glue::glue("'res' must be a DESeqResults or data.frame, not {class(res)[1]}."))
  }
  
  mandatory_cols <- c("baseMean", "log2FoldChange", "padj")
  missing_cols   <- mandatory_cols[!mandatory_cols %in% colnames(res_df)]
  if (length(missing_cols) > 0) {
    log_error(sample = "", step = "plot_ma",
              msg = glue::glue("Column(s) missing from res: {paste(missing_cols, collapse = ', ')}. ",
                               "Available: {paste(colnames(res_df), collapse = ', ')}"))
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Prepare Data
  # ═══════════════════════════════════════════════════════════════════════════
  
  # For ggrepel label reproducibility
  set.seed(1234)
  
  # Prepare the data with "squished" values for capped plot variant.
  # Points with |LFC| > 5 are plotted as triangles at the cap boundary
  # so extreme outliers don't compress the scale for the bulk of the data.
  res_df <- res_df %>%
    dplyr::mutate(
      significant           = padj < 0.1 & !is.na(padj),
      log2FoldChange_capped = dplyr::case_when(
        log2FoldChange >  5 ~  5,
        log2FoldChange < -5 ~ -5,
        TRUE                ~ log2FoldChange
      ),
      # 16 = solid circle (normal points), 17 = solid triangle (capped outliers)
      is_outlier  = !is.na(log2FoldChange) & (log2FoldChange > 5 | log2FoldChange < -5),
      point_shape = ifelse(is_outlier, 17, 16)
    )
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Build Plots
  # ═══════════════════════════════════════════════════════════════════════════
  
  base_ma <- ggplot(data    = res_df,
                    mapping = aes(x = baseMean, color = significant)) +
    geom_point(alpha = 0.5, size = 1.25) +
    theme_publication() +
    theme(legend.position = "none") +  # overrides the "right" default just for MA
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_color_manual(values = c("grey70", "#E64B35")) +
    labs(subtitle = "Red points indicate FDR < 0.1",
         x        = "Mean of Normalized Counts",
         y        = "Log2 Fold Change")
  
  ma_uncapped <- base_ma +
    aes(y = log2FoldChange) +
    labs(title = "MA Plot")
  
  ma_capped <- base_ma +
    aes(y = log2FoldChange_capped, shape = point_shape) +
    scale_shape_identity() +
    coord_cartesian(ylim = c(-5, 5)) +
    labs(title = "MA Plot (with Capped Fold Changes)")
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Save and Return
  # ═══════════════════════════════════════════════════════════════════════════
  
  if (!is.null(output_dir)) {
    
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    # Sanitise filename — replace non-alphanumeric runs with underscores and
    # strip any leading/trailing underscores introduced by the substitution
    safe_filename <- gsub("[^[:alnum:]]+", "_", filename)
    safe_filename <- gsub("(^_+|_+$)", "", safe_filename)
    
    output_file <- file.path(output_dir, paste0(safe_filename, ".pdf"))
    
    # Two pages per PDF: uncapped first (full dynamic range), capped second
    # (outliers triangled at ±5 so the bulk of points remain readable).
    # cairo_pdf is used over ggsave because it natively supports multi-page
    # output and handles font embedding and transparency correctly.
    grDevices::cairo_pdf(filename = output_file,
                         width    = 8,
                         height   = 11.5,
                         onefile  = TRUE)
    print(ma_uncapped)
    print(ma_capped)
    grDevices::dev.off()
    
    log_info(sample = "", step = "plot_ma",
             msg = glue::glue("MA plots saved to: '{output_file}'"))
    
  } else {
    
    log_info(sample = "", step = "plot_ma",
             msg = "MA plots returned as list. No output_dir provided — skipping save.")
  }
  
  return(invisible(list(uncapped = ma_uncapped,
                        capped   = ma_capped)))
}

# ═══════════════════════════════════════════════════════════════════════════════
# VOLCANO PLOT
# ═══════════════════════════════════════════════════════════════════════════════
# plot_volcano() generates a two-page volcano plot PDF: page 1 unlabeled
# (clean overview), page 2 with top_n gene labels per direction via ggrepel.
#
# INPUT:
#   df        : data.frame with DEG results. Must contain x_col, y_col, label_col.
#               This is the same format produced by process_degs().
#   x_col     : column for x-axis (log2 fold change). Default "log2FoldChange".
#   y_col     : column for y-axis (p-value or adjusted p-value). Default "padj".
#   label_col : column containing gene labels for ggrepel. Default "SYMBOL".
#   y_transform: logical. TRUE (default) applies -log10 transform to y_col.
#   y_cap     : upper cap for -log10(y) axis to prevent extreme points from
#               compressing the scale. Default 50.
#   shape_col : optional metadata column to map to point shape. NULL = circles only.
#   x_cutoff  : absolute LFC threshold for significance lines. Default 0.58 (≈1.5x).
#   y_cutoff  : p-value threshold for significance line. Default 0.05.
#   top_n     : number of top genes to label per direction (up/down). Default 5.
#   label_genes: character vector of specific gene names to label regardless of
#                significance. NULL = use top_n only.
#   filename  : base name for the output PDF. Non-alphanumeric characters are
#               sanitised to underscores. Default "Volcano_Plot".
#   output_dir: directory to save the PDF. NULL = return plots only, no file saved.
#   title     : optional string for the plot title. NULL = no title.
#
# OUTPUT:
#   Returns invisibly: list(
#     unlabeled = ggplot (no gene labels),
#     labeled   = ggplot (top_n genes labeled per direction)
#   )
#   When output_dir is provided, saves a two-page PDF:
#     Page 1 : unlabeled volcano
#     Page 2 : labeled volcano

plot_volcano <- function(df,
                         x_col       = "log2FoldChange",
                         y_col       = "padj",
                         label_col   = "SYMBOL",
                         y_transform = TRUE,
                         y_cap       = 50,
                         shape_col   = NULL,
                         x_cutoff    = 0.58,
                         y_cutoff    = 0.05,
                         top_n       = 5,
                         label_genes = NULL,
                         filename    = "Volcano_Plot",
                         output_dir  = NULL,
                         title       = NULL) {
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Validate Inputs
  # ═══════════════════════════════════════════════════════════════════════════
  # All mandatory columns and any user-specified aesthetic columns are checked
  # upfront so the function fails immediately with a clear message rather than
  # crashing deep inside ggplot with a cryptic error.
  
  mandatory_cols <- c(x_col, y_col, label_col)
  optional_cols  <- shape_col  # only user-configurable aesthetic
  cols_to_check  <- c(mandatory_cols, optional_cols)
  
  missing_cols <- cols_to_check[!cols_to_check %in% colnames(df)]
  if (length(missing_cols) > 0) {
    log_error(sample = "", step = "plot_volcano",
              msg = glue::glue(
                "Column(s) not found in df: {paste(missing_cols, collapse = ', ')}. ",
                "Available: {paste(colnames(df), collapse = ', ')}"
              ))
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Prepare Data
  # ═══════════════════════════════════════════════════════════════════════════
  
  # ── 2a. Drop rows with NA in x or y ─────────────────────────────────────
  # NA in x or y would produce a blank point — silently misleading. We drop
  # and warn so the user knows their data was trimmed.
  n_before  <- nrow(df)
  df        <- df %>% dplyr::filter(!is.na(.data[[x_col]]), !is.na(.data[[y_col]]))
  n_dropped <- n_before - nrow(df)
  
  if (n_dropped > 0) {
    log_warn(sample = "", step = "plot_volcano",
             msg = glue::glue("Dropped {n_dropped} row(s) with NA in '{x_col}' or '{y_col}'."))
  }
  
  # ── 2b. Handle y = 0 before log transformation ──────────────────────────
  # -log10(0) = Inf, which breaks axis scaling. Replace with the smallest
  # observed non-zero value so these points plot at the top of the y-axis
  # without distorting the scale — they are genuinely the most significant.
  if (y_transform) {
    zero_y <- df[[y_col]] == 0
    if (any(zero_y, na.rm = TRUE)) {
      min_nonzero         <- min(df[[y_col]][df[[y_col]] > 0], na.rm = TRUE)
      df[[y_col]][zero_y] <- min_nonzero
      log_warn(sample = "", step = "plot_volcano",
               msg = glue::glue(
                 "{sum(zero_y)} row(s) had {y_col} = 0. ",
                 "Replaced with minimum non-zero value ({min_nonzero}) to avoid log10(0) = Inf."
               ))
    }
  }
  
  # ── 2c. Compute plot_y ───────────────────────────────────────────────────
  # y_cap prevents a handful of extreme points from compressing the rest of
  # the plot. Points above the cap are real — they are just shown at the
  # ceiling so the bulk of the data remains readable.
  if (y_transform) {
    df <- df %>%
      dplyr::mutate(plot_y = pmin(-log10(.data[[y_col]]), y_cap))
  } else {
    df <- df %>%
      dplyr::mutate(plot_y = .data[[y_col]])
  }
  
  # ── 2d. Compute color direction from x + y cutoffs ──────────────────────
  # Color encodes position relative to both thresholds simultaneously — a point
  # must clear BOTH the effect size and significance bar to be colored. This is
  # more stringent than coloring by either axis alone and matches standard
  # volcano interpretation across all omics data types.
  #
  # Why three levels rather than two?
  # "Up" and "Down" carry meaning only relative to the x-axis direction, which
  # the user controls via x_col. NS captures everything that fails either cutoff.
  y_thresh <- if (y_transform) -log10(y_cutoff) else y_cutoff
  
  df <- df %>%
    dplyr::mutate(
      .direction = dplyr::case_when(
        !is.null(x_cutoff) & !is.null(y_cutoff) &
          .data[[x_col]] >  x_cutoff & plot_y > y_thresh ~ "Up",
        !is.null(x_cutoff) & !is.null(y_cutoff) &
          .data[[x_col]] < -x_cutoff & plot_y > y_thresh ~ "Down",
        TRUE ~ "Not Significant"
      )
    )
  
  # ── 2e. Compute alpha from FDR bins ─────────────────────────────────────
  # Fixed significance tiers (0.001 / 0.01 / 0.05) are universal thresholds
  # used across all omics fields — no reason to make these dynamic. More
  # significant points are more opaque, naturally drawing the eye upward.
  # Alpha is only meaningful when y_col is a p-value / FDR; when y_transform
  # is FALSE the user is plotting a non-p-value y and we fall back to fixed alpha.
  if (y_transform) {
    df <- df %>%
      dplyr::mutate(
        .fdr_bin = dplyr::case_when(
          .data[[y_col]] <= 0.001 ~ "p \u2264 0.001",
          .data[[y_col]] <= 0.01  ~ "p \u2264 0.01",
          .data[[y_col]] <= 0.05  ~ "p \u2264 0.05",
          TRUE                    ~ "NS"
        )
      )
    alpha_values <- c("p \u2264 0.001" = 1.0,
                      "p \u2264 0.01"  = 0.8,
                      "p \u2264 0.05"  = 0.6,
                      "NS"            = 0.3)
  }
  
  # ── 2f. Compute relevance score for top_n labeling ──────────────────────
  # abs(x) * -log10(y) rewards points that are extreme on BOTH axes, which
  # are the most scientifically interesting. Points that are only significant
  # (near x = 0) or only large in effect (p ≈ 1) score low and are not labeled.
  # The score is capped at the 99th percentile to prevent one extreme outlier
  # from dominating the ranking and crowding out other interesting points.
  df <- df %>%
    dplyr::mutate(
      .relevance = {
        r <- abs(.data[[x_col]]) * plot_y
        pmin(r, quantile(r, 0.99, na.rm = TRUE))
      }
    )
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Axis Limits and Breaks
  # ═══════════════════════════════════════════════════════════════════════════
  # Target: 5 ticks on x (2 negative, 0, 2 positive), 5 ticks on y.
  # We pick the smallest clean step that places exactly 2 breaks on the
  # positive x half — guaranteeing a symmetric, uncluttered axis regardless
  # of LFC range. Same logic for y targeting 4 intervals.
  # x_steps covers ranges from small effect sizes (0.5) to extreme datasets
  # (500) — which.min finds the step closest to placing exactly 2 breaks on
  # the positive half for any input data.
  
  x_vals <- df[[x_col]][is.finite(df[[x_col]])]
  y_vals <- df$plot_y[is.finite(df$plot_y)]
  
  # ── X axis ───────────────────────────────────────────────────────────────
  x_max    <- max(abs(x_vals), na.rm = TRUE)
  x_steps  <- c(0.5, 1, 2, 3, 5, 10, 15, 20, 25, 50, 100, 200, 500)
  x_step   <- x_steps[which.min(abs(ceiling(x_max / x_steps) - 2))]
  x_extent <- ceiling(x_max / x_step) * x_step
  x_breaks <- seq(-x_extent, x_extent, by = x_step)
  
  # ── Y axis ───────────────────────────────────────────────────────────────
  # Always starts at 0 (-log10(p) is non-negative). Targets 4 intervals
  # giving 5 ticks. y_steps chosen to land on round -log10(p) values.
  # Range covers padj values down to 10^-400 — beyond any real dataset.
  y_max    <- max(y_vals, na.rm = TRUE)
  y_steps  <- c(5, 10, 20, 25, 50, 100)
  y_step   <- y_steps[which.min(abs(ceiling(y_max / y_steps) - 4))]
  y_top    <- ceiling(y_max / y_step) * y_step
  y_breaks <- seq(0, y_top, by = y_step)
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Color and Alpha Palettes
  # ═══════════════════════════════════════════════════════════════════════════
  
  # Option C: Nature / Cell style — colorblind safe, publication quality (default)
  color_values <- c("Up" = "#E64B35", "Down" = "#4DBBD5", "Not Significant" = "grey80")
  
  # Option A: Classic traffic light — most common in published papers
  # color_values <- c("Up" = "#D62728", "Down" = "#1F77B4", "Not Significant" = "grey80")
  
  # Option B: Viridis endpoints — colorblind safe
  # color_values <- c("Up" = viridis::viridis(100)[95], "Down" = viridis::viridis(100)[5], "Not Significant" = "grey80")
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5: Build Base Plot
  # ═══════════════════════════════════════════════════════════════════════════
  
  # ── 5a. Base aesthetic mapping ───────────────────────────────────────────
  # shape_col is the only user-configurable aesthetic — all others are derived
  # internally. .data[[]] notation is used throughout for programmatic column
  # access so the function works regardless of actual column names.
  base_aes <- if (!is.null(shape_col)) {
    ggplot2::aes(x     = .data[[x_col]],
                 y     = plot_y,
                 color = .direction,
                 alpha = if (y_transform) .fdr_bin else NULL,
                 shape = .data[[shape_col]])
  } else {
    ggplot2::aes(x     = .data[[x_col]],
                 y     = plot_y,
                 color = .direction,
                 alpha = if (y_transform) .fdr_bin else NULL)
  }
  
  # ── 5b. Y-axis label ────────────────────────────────────────────────────
  # Use proper mathematical notation when we transformed — matches convention
  # in every published volcano. Plain label when user passes pre-transformed y.
  y_label <- if (y_transform) expression("-log"[10]*"(p-value)") else y_col
  
  # ── 5c. Assemble base ggplot ─────────────────────────────────────────────
  p <- ggplot2::ggplot(data = df, mapping = base_aes) +
    ggplot2::geom_point(size = 1.5, na.rm = TRUE) +
    theme_publication() + 
    ggplot2::labs(
      x     = x_col,
      y     = y_label,
      color = "Direction",
      title = title
    ) +
    ggplot2::scale_color_manual(values = color_values, na.translate = FALSE) +
    ggplot2::coord_cartesian(
      xlim = c(min(x_breaks), max(x_breaks)),
      ylim = c(min(y_breaks), max(y_breaks))
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      # Integers for whole numbers, 2 significant digits for decimals
      labels = function(x) ifelse(x %% 1 == 0, as.integer(x), format(x, digits = 2))
    ) +
    ggplot2::scale_y_continuous(breaks = y_breaks) +
    ggplot2::guides(alpha = "none")   # alpha is a visual aid, not worth a legend
  
  # ── 5d. Alpha scale (only when y is a p-value) ──────────────────────────
  if (y_transform) {
    p <- p + ggplot2::scale_alpha_manual(values = alpha_values, na.translate = FALSE)
  }
  
  # ── 5e. Reference lines ──────────────────────────────────────────────────
  # Drawn as dotted lines — prominent enough to read, unobtrusive enough to
  # not compete with the points.
  if (!is.null(x_cutoff)) {
    p <- p + ggplot2::geom_vline(xintercept = c(-x_cutoff, x_cutoff),
                                 linetype   = "dotted",
                                 color      = "black",
                                 linewidth  = 0.5)
  }
  
  if (!is.null(y_cutoff)) {
    p <- p + ggplot2::geom_hline(yintercept = y_thresh,
                                 linetype   = "dotted",
                                 color      = "black",
                                 linewidth  = 0.5)
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 6: Gene Labeling
  # ═══════════════════════════════════════════════════════════════════════════
  
  # ── 6a. Filter uninformative label names ────────────────────────────────
  # Predicted, uncharacterized, or computationally-named entries clutter the
  # plot without adding biological insight. The pattern covers common naming
  # conventions across mouse, human, and other organisms.
  uninformative_pattern <- paste0(
    "^(Gm[0-9]+",          # mouse predicted genes
    "|ENS[A-Z]*G[0-9]+",   # any-organism Ensembl gene IDs
    "|LOC[0-9]+",          # uncharacterized loci
    "|C[0-9]+orf[0-9]+",   # human open reading frames
    "|RP[0-9]+-",          # ribosomal pseudogenes
    "|AC[0-9]+\\.[0-9]+",  # human novel transcripts (AC prefix)
    "|AL[0-9]+\\.[0-9]+",  # human novel transcripts (AL prefix)
    "|LINC[0-9]+",         # long intergenic non-coding RNAs
    "|[0-9]+Rik",          # RIKEN cDNA clones (numeric prefix variant)
    ")|Rik$"               # RIKEN cDNA clones (Rik suffix)
  )
  
  # ── 6b. Resolve labels for unlabeled plot ───────────────────────────────
  # When label_genes are provided, show them on BOTH plots. The unlabeled plot
  # is not truly unlabeled in this case — it just shows user-curated genes
  # without the auto top_n layer, which keeps it clean and interpretable.
  labels_unlabeled <- if (!is.null(label_genes)) {
    intersect(label_genes, df[[label_col]])
  } else {
    NULL
  }
  
  # ── 6c. Resolve labels for labeled plot ─────────────────────────────────
  # Top N are ranked by relevance score within Up and Down groups separately
  # so both directions are represented equally regardless of asymmetry in the
  # data. Uninformative names are excluded before ranking so the top N slots
  # are always filled by interpretable gene symbols.
  top_n_genes <- df %>%
    dplyr::filter(
      !stringr::str_detect(.data[[label_col]], uninformative_pattern),
      .direction != "Not Significant"
    ) %>%
    dplyr::group_by(.direction) %>%
    dplyr::slice_max(order_by = .relevance, n = top_n, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::pull(.data[[label_col]])
  
  labels_labeled <- union(
    if (!is.null(label_genes)) intersect(label_genes, df[[label_col]]) else character(0),
    top_n_genes
  )
  
  # ── 6d. ggrepel layer ────────────────────────────────────────────────────
  # set.seed ensures label positions are reproducible across re-runs.
  # Curved segments help trace labels back to dense point clouds without
  # crossing other labels — ncp = 50 gives smooth curves.
  add_repel_layer <- function(plot, genes_to_label) {
    
    if (length(genes_to_label) == 0) return(plot)
    
    set.seed(1234)
    
    plot + ggrepel::geom_text_repel(
      data             = df %>% dplyr::filter(.data[[label_col]] %in% genes_to_label),
      mapping          = ggplot2::aes(label = .data[[label_col]]),
      direction        = "both",
      box.padding      = 0.8,
      point.padding    = 0.1,
      max.overlaps     = nrow(df),
      show.legend      = FALSE,
      min.segment.length  = 0,
      segment.curvature   = -0.5,
      segment.ncp         = 50,
      segment.angle       = 20,
      segment.size        = 0.3,
      size                = 4,
      segment.color = "grey50"
    )
  }
  
  p_unlabeled <- add_repel_layer(p, labels_unlabeled)
  p_labeled   <- add_repel_layer(p, labels_labeled)
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 7: Save and Return
  # ═══════════════════════════════════════════════════════════════════════════
  
  if (!is.null(output_dir)) {
    
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    # Sanitise filename — replace non-alphanumeric runs with underscores and
    # strip any leading/trailing underscores introduced by the substitution
    safe_filename <- gsub("[^[:alnum:]]+", "_", filename) |>
      gsub(pattern = "(^_+|_+$)", replacement = "", x = _)
    
    output_file <- file.path(output_dir, paste0(safe_filename, ".pdf"))
    
    # Two pages per PDF: unlabeled first (clean overview), labeled second
    # (annotated for interpretation). cairo_pdf is used over ggsave because
    # it natively supports multi-page output and handles font embedding and
    # transparency correctly — critical for publication-quality PDFs.
    grDevices::cairo_pdf(filename = output_file,
                         width    = 8,
                         height   = 11.5,
                         onefile  = TRUE)
    print(p_unlabeled)
    print(p_labeled)
    grDevices::dev.off()
    
    log_info(sample = "", step = "plot_volcano",
             msg = glue::glue("Volcano plots saved to: '{output_file}'"))
    
  } else {
    
    log_info(sample = "", step = "plot_volcano",
             msg = "Volcano plots returned as list. No output_dir provided — skipping save.")
  }
  
  return(invisible(list(unlabeled = p_unlabeled,
                        labeled   = p_labeled)))
}


# ═══════════════════════════════════════════════════════════════════════════════
# HEATMAP FUNCTION
# ═══════════════════════════════════════════════════════════════════════════════
# plot_heatmap() generates a pheatmap with full control over transforms,
# clustering, annotations, and color palettes.
#
# INPUT:
#   df      : data.frame with one gene ID column (row_id) + numeric sample columns.
#             All other non-numeric columns are dropped automatically in Section 1.
#             This is the same format produced by process_counts() and process_degs().
#   row_id  : name of the gene ID column in df (default "SYMBOL").
#             If user passes a custom name, it is used; default is checked otherwise.
#   col_id  : name of the sample ID column in metadata_col (default "Sample_ID").
#             Sample columns in df are derived as the intersection of
#             metadata_col[[col_id]] and colnames(df) — stray non-sample columns
#             (gene_id, gene_biotype etc.) are excluded automatically.
#
# data_type controls log transform, z-scoring, and color break behavior.
# Pass "auto" to detect from the matrix, or override explicitly:
#   "counts"        - Raw/normalized counts; log2(1+x) then row-wise z-score
#   "counts_log"    - Already log-transformed (VST, rlog); row-wise z-score only
#   "counts_scaled" - Already log-transformed AND z-scored; no transforms
#   "centered"      - Zero-centered scores (NES, LFC); symmetric color breaks only
#   "correlation"   - Correlation matrix [-1, 1]; fixed color breaks
#   "bounded"       - Values in [0, 1] (beta, p-values); fixed color breaks
#
# col_cluster_by / row_cluster_by (same options for both axes):
#   "all"          - single hierarchical clustering across all columns/rows
#   "none"         - preserve input order exactly
#   "alphabetical" - sort by name
#   <column name>  - cluster within each group separately, then concatenate
#                    Must be a single string matching a column in col/row annotation.
#
# heatmap_palette:
#   "rdbu"    - reversed RdBu diverging (red = high, blue = low)
#   "viridis" - Viridis sequential; any viridis variant also works
#               ("magma", "plasma", "inferno", "cividis", "mako", "rocket", "turbo")
#   Any valid RColorBrewer palette name (e.g. "PuOr", "Blues")
#   A custom character vector of colors (interpolated to n_breaks)
#
# annotation_palette:
#   Character vector of colors used for annotation sidebar bars (discrete levels
#   and numeric gradient bases). NULL (default) uses the built-in 20-color palette.
#   If fewer colors are supplied than the maximum number of levels across all
#   annotation columns, colors are recycled with a log_warn.

plot_heatmap <- function(df,
                         row_id             = "SYMBOL",
                         col_id             = "Sample_ID",
                         data_type          = "auto",
                         label_genes        = NULL,       # Character vector of gene names to label; others blanked
                         metadata_col       = NULL,       # data.frame with col_id column for column annotations
                         metadata_row       = NULL,       # data.frame with row_id column for row annotations
                         col_annotations    = NULL,       # Character vector of metadata_col columns to display
                         row_annotations    = NULL,       # Character vector of metadata_row columns to display
                         col_gap_by         = NULL,       # Single metadata_col column to draw column gaps by
                         row_gap_by         = NULL,       # Single metadata_row column to draw row gaps by
                         col_cluster_by     = "all",
                         row_cluster_by     = "all",
                         plot_title         = NULL,
                         heatmap_palette    = "rdbu",
                         annotation_palette = NULL,       # NULL = use internal 31-color palette
                         border_color       = NA,
                         show_expr_legend   = TRUE) {
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Validate Inputs
  # ═══════════════════════════════════════════════════════════════════════════
  # All structural checks fire here before any data is touched. The validation
  # order mirrors plot_volcano(): object type → ID columns → annotation columns
  # → parameter sanity → derive sample columns → numeric check → convert to matrix.
  # This ensures every possible user error produces a clear, actionable message
  # rather than a cryptic crash deep inside pheatmap or dplyr.
  
  # ── 1a. Input type ───────────────────────────────────────────────────────
  # df must be a data.frame or tibble — matrix-with-rownames is no longer
  # accepted. This keeps the input contract consistent with every other
  # pipeline function (process_counts, process_degs etc. all return data.frames).
  if (!is.data.frame(df)) {
    log_error(sample = "", step = "plot_heatmap",
              msg = glue::glue(
                "'df' must be a data.frame, not {class(df)[1]}. ",
                "Pass a data.frame with a '{row_id}' column and numeric sample columns. ",
                "process_counts() and process_degs() both return data.frames in this format."
              ))
  }
  
  # ── 1b. row_id column must exist in df ───────────────────────────────────
  if (!row_id %in% colnames(df)) {
    log_error(sample = "", step = "plot_heatmap",
              msg = glue::glue(
                "row_id column '{row_id}' not found in df. ",
                "Available columns: {paste(colnames(df), collapse = ', ')}"
              ))
  }
  
  # ── 1c. col_id column must exist in metadata_col (if provided) ───────────
  if (!is.null(metadata_col) && !col_id %in% colnames(metadata_col)) {
    log_error(sample = "", step = "plot_heatmap",
              msg = glue::glue(
                "col_id column '{col_id}' not found in metadata_col. ",
                "Available columns: {paste(colnames(metadata_col), collapse = ', ')}"
              ))
  }
  
  # ── 1d. col_annotations columns must exist in metadata_col ──────────────
  # Checked here rather than in Section 6 so we fail before any transforms.
  if (!is.null(col_annotations) && !is.null(metadata_col)) {
    missing_col_anns <- setdiff(col_annotations, colnames(metadata_col))
    if (length(missing_col_anns) > 0) {
      log_error(sample = "", step = "plot_heatmap",
                msg = glue::glue(
                  "col_annotations column(s) not found in metadata_col: ",
                  "{paste(missing_col_anns, collapse = ', ')}"
                ))
    }
  }
  
  # ── 1e. row_annotations columns must exist in metadata_row ──────────────
  if (!is.null(row_annotations) && !is.null(metadata_row)) {
    missing_row_anns <- setdiff(row_annotations, colnames(metadata_row))
    if (length(missing_row_anns) > 0) {
      log_error(sample = "", step = "plot_heatmap",
                msg = glue::glue(
                  "row_annotations column(s) not found in metadata_row: ",
                  "{paste(missing_row_anns, collapse = ', ')}"
                ))
    }
  }
  
  # ── 1f. col_cluster_by and row_cluster_by must be length-1 strings ───────
  # Accepting a vector here would cause %in% to return a vector, and if() would
  # silently take only the first element — producing wrong clustering with no
  # error. Catching this upfront is essential.
  if (!is.character(col_cluster_by) || length(col_cluster_by) != 1) {
    log_error(sample = "", step = "plot_heatmap",
              msg = glue::glue(
                "'col_cluster_by' must be a single string ('all', 'none', 'alphabetical', ",
                "or a metadata column name). Got: {paste(col_cluster_by, collapse = ', ')}"
              ))
  }
  if (!is.character(row_cluster_by) || length(row_cluster_by) != 1) {
    log_error(sample = "", step = "plot_heatmap",
              msg = glue::glue(
                "'row_cluster_by' must be a single string ('all', 'none', 'alphabetical', ",
                "or a metadata column name). Got: {paste(row_cluster_by, collapse = ', ')}"
              ))
  }
  
  # ── 1g. col_gap_by and row_gap_by must be length-1 if provided ──────────
  if (!is.null(col_gap_by) && length(col_gap_by) != 1) {
    log_error(sample = "", step = "plot_heatmap",
              msg = "'col_gap_by' must be a single column name string or NULL.")
  }
  if (!is.null(row_gap_by) && length(row_gap_by) != 1) {
    log_error(sample = "", step = "plot_heatmap",
              msg = "'row_gap_by' must be a single column name string or NULL.")
  }
  
  # ── 1h. data_type must be a recognised string ────────────────────────────
  valid_types <- c("auto", "counts", "counts_log", "counts_scaled",
                   "centered", "correlation", "bounded")
  if (!data_type %in% valid_types) {
    log_error(sample = "", step = "plot_heatmap",
              msg = glue::glue(
                "`data_type` '{data_type}' is not valid. ",
                "Must be one of: {paste(valid_types, collapse = ', ')}."
              ))
  }
  
  # ── 1i. heatmap_palette basic type check ────────────────────────────────
  # Full resolution (rdbu / viridis / RColorBrewer / custom vector) happens in
  # Section 5. Here we just reject obviously wrong types early.
  if (!is.character(heatmap_palette)) {
    log_error(sample = "", step = "plot_heatmap",
              msg = "'heatmap_palette' must be a character string or character vector of colors.")
  }
  
  # ── 1j. Derive sample columns ────────────────────────────────────────────
  # sample_cols = intersection of the authoritative sample list from metadata
  # and the columns actually present in df. This automatically excludes stray
  # non-sample columns (gene_id, gene_biotype etc.) without any hardcoding.
  # When no metadata is provided, fall back to all non-row_id columns with a
  # warning — this supports standalone use without annotation sidebars.
  if (!is.null(metadata_col)) {
    sample_cols <- intersect(metadata_col[[col_id]], colnames(df))
    
    missing_samples <- setdiff(metadata_col[[col_id]], colnames(df))
    if (length(missing_samples) > 0) {
      log_warn(sample = "", step = "plot_heatmap",
               msg = glue::glue(
                 "{length(missing_samples)} sample(s) in metadata not found in df and will be excluded: ",
                 "{paste(missing_samples, collapse = ', ')}"
               ))
    }
  } else {
    log_warn(sample = "", step = "plot_heatmap",
             msg = glue::glue(
               "No metadata_col provided — inferring sample columns as all non-'{row_id}' columns. ",
               "Pass metadata_col to be explicit and avoid including stray columns."
             ))
    sample_cols <- setdiff(colnames(df), row_id)
  }
  
  if (length(sample_cols) == 0) {
    log_error(sample = "", step = "plot_heatmap",
              msg = glue::glue(
                "No sample columns found after intersecting metadata_col[['{col_id}']] with colnames(df). ",
                "Check that sample names match between df and metadata_col."
              ))
  }
  
  # ── 1k. Sample columns must be numeric ───────────────────────────────────
  # Non-numeric columns (e.g. a gene_biotype that slipped through) would be
  # coerced to NA silently in the matrix conversion below, corrupting the data.
  non_numeric <- sample_cols[!sapply(df[sample_cols], is.numeric)]
  if (length(non_numeric) > 0) {
    log_error(sample = "", step = "plot_heatmap",
              msg = glue::glue(
                "The following sample columns are not numeric: {paste(non_numeric, collapse = ', ')}. ",
                "Check that col_id / metadata_col correctly identifies sample columns."
              ))
  }
  
  # ── 1l. Subset df and convert to numeric matrix ──────────────────────────
  # Subset to row_id + sample_cols only, then convert to matrix with gene IDs
  # as rownames. This is the single point of conversion — all downstream code
  # operates on a pure numeric matrix with no stray character columns.
  expr_mat <- df %>%
    dplyr::select(dplyr::all_of(c(row_id, sample_cols))) %>%
    tibble::column_to_rownames(row_id) %>%
    as.matrix()
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Prepare Matrix
  # ═══════════════════════════════════════════════════════════════════════════
  # make.names() sanitises names that contain spaces or special characters —
  # pheatmap uses colnames/rownames for annotation alignment and will silently
  # misalign if names contain characters R treats as operators or separators.
  
  colnames(expr_mat) <- make.names(colnames(expr_mat))
  
  na_rows  <- is.na(rownames(expr_mat))
  expr_mat <- expr_mat[!na_rows, , drop = FALSE]
  if (sum(na_rows) > 0) {
    log_info(sample = "", step = "plot_heatmap",
             msg = glue::glue("Removed {sum(na_rows)} rows with NA gene names."))
  }
  
  expr_mat[is.na(expr_mat)] <- 0
  
  # Why abs() in the zero-gene filter?
  # For centered data types (NES, LFC), a gene with a large negative score
  # would have rowSums < 0 but is NOT a zero-expression gene. abs() ensures
  # we only remove genes that are truly zero across all samples, not genes
  # that are strongly downregulated.
  zero_genes <- rownames(expr_mat)[rowSums(abs(expr_mat)) == 0]
  expr_mat   <- expr_mat[rowSums(abs(expr_mat)) != 0, , drop = FALSE]
  if (length(zero_genes) > 0) {
    log_info(sample = "", step = "plot_heatmap",
             msg = glue::glue("Removed {length(zero_genes)} all-zero genes."))
  }
  
  zero_samples <- colnames(expr_mat)[colSums(abs(expr_mat)) == 0]
  expr_mat     <- expr_mat[, colSums(abs(expr_mat)) != 0, drop = FALSE]
  if (length(zero_samples) > 0) {
    log_info(sample = "", step = "plot_heatmap",
             msg = glue::glue("Removed {length(zero_samples)} all-zero samples."))
  }
  
  if (nrow(expr_mat) < 2) {
    log_warn(sample = "", step = "plot_heatmap",
             msg = "Input has fewer than 2 genes after filtering. Skipping heatmap.")
    return(NULL)
  }
  
  # When duplicates exist, keep the row with the highest total absolute expression.
  # abs() is used again so large negative NES/LFC rows are not incorrectly
  # treated as low-expression and discarded.
  if (any(duplicated(rownames(expr_mat)))) {
    log_warn(sample = "", step = "plot_heatmap",
             msg = "Duplicate gene symbols detected. Retaining highest-expressing copy per gene.")
    expr_mat <- expr_mat %>%
      as.data.frame() %>%
      tibble::rownames_to_column(row_id) %>%
      dplyr::mutate(total_expr = rowSums(dplyr::across(-dplyr::all_of(row_id), abs))) %>%
      dplyr::group_by(.data[[row_id]]) %>%
      dplyr::slice_max(order_by = total_expr, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      tibble::column_to_rownames(row_id) %>%
      dplyr::select(-total_expr) %>%
      as.matrix()
  }
  
  # make.names() on both axes — same reasoning as column sanitisation above;
  # gene symbol rownames may contain "-", "/", or spaces
  rownames(expr_mat) <- make.names(rownames(expr_mat))
  colnames(expr_mat) <- make.names(colnames(expr_mat))
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Data Type Detection
  # ═══════════════════════════════════════════════════════════════════════════
  # Delegated to detect_data_type() — a shared helper used by plot_heatmap()
  # and plot_pca(). Detection logic, thresholds, and log messages are defined
  # once there and applied consistently across both functions.
  
  detected_type <- detect_data_type(mat = expr_mat, data_type = data_type)
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Log Transform and Z-score Scaling
  # ═══════════════════════════════════════════════════════════════════════════
  # Transform decisions by data type:
  #   counts        : log2(1+x) then row-wise z-score
  #   counts_log    : row-wise z-score only (already logged)
  #   everything else: no transforms (already on a meaningful absolute scale)
  #
  # Why row-wise z-score for count data?
  # Raw or log-normalised counts vary enormously in absolute magnitude across
  # genes — a highly expressed housekeeping gene would dominate color scale.
  # Row-wise z-scoring centres each gene at 0 and scales by its own SD, so
  # color encodes RELATIVE variation across samples, not absolute expression.
  # This is the standard approach for expression heatmaps.
  
  if (detected_type == "counts") {
    log_info(sample = "", step = "plot_heatmap", msg = "Applying log2(1+x) transform.")
    expr_mat <- log2(1 + expr_mat)
  }
  
  if (detected_type %in% c("counts", "counts_log")) {
    log_info(sample = "", step = "plot_heatmap", msg = "Applying row-wise z-score scaling.")
    # Transpose -> scale (operates on columns) -> transpose back = row-wise z-score
    expr_mat_scaled <- expr_mat %>% t() %>% scale() %>% t()
  } else {
    log_info(sample = "", step = "plot_heatmap",
             msg = glue::glue("Skipping transforms for data_type = '{detected_type}'."))
    expr_mat_scaled <- expr_mat
  }
  
  # Replace NAs introduced by z-scoring zero-variance rows
  expr_mat_scaled[is.na(expr_mat_scaled)] <- 0
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5: Color Breaks
  # ═══════════════════════════════════════════════════════════════════════════
  # Break range by data type:
  #   bounded     : fixed [0, 1]
  #   correlation : fixed [-1, 1]
  #   all others  : data-driven ±99th percentile after transforms
  #
  # Why 99th percentile and not the actual min/max?
  # Extreme outlier genes (e.g. one gene with z-score = 15 due to batch effect)
  # would compress the entire color scale so that 99% of genes appear the same
  # mid-range color. The 99th percentile clips these outliers while preserving
  # color resolution for the bulk of the data. Symmetric construction ensures
  # zero maps exactly to the midpoint color regardless of skew.
  
  n_breaks <- 100
  
  if (detected_type == "bounded") {
    scale_min <- 0
    scale_max <- 1
    breaks    <- seq(from = 0, to = 1, length.out = n_breaks + 1)
    
  } else if (detected_type == "correlation") {
    scale_min <- -1
    scale_max <-  1
    breaks    <- seq(from = -1, to = 1, length.out = n_breaks + 1)
    
  } else {
    p99       <- stats::quantile(abs(as.vector(expr_mat_scaled)), probs = 0.99, na.rm = TRUE)
    scale_min <- -p99
    scale_max <-  p99
    breaks    <- c(seq(from = scale_min, to = 0,         length.out = n_breaks / 2 + 1),
                   seq(from = 0,         to = scale_max, length.out = n_breaks / 2 + 1)[-1])
  }
  
  log_info(sample = "", step = "plot_heatmap",
           msg = glue::glue("Color scale range: [{round(scale_min, 2)}, {round(scale_max, 2)}]"))
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 6: Color Palette
  # ═══════════════════════════════════════════════════════════════════════════
  # Resolution order: custom vector → rdbu → viridis family → RColorBrewer.
  # The log_error at the bottom is a safety net for any unrecognised string
  # that doesn't match any of the above branches.
  
  heatmap_colors <- if (is.character(heatmap_palette) && length(heatmap_palette) > 1) {
    colorRampPalette(heatmap_palette)(n_breaks)
    
  } else if (heatmap_palette == "rdbu") {
    colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(n_breaks)
    
  } else if (heatmap_palette %in% c("viridis", "magma", "plasma", "inferno",
                                    "cividis", "mako", "rocket", "turbo")) {
    viridis::viridis(n = n_breaks, option = heatmap_palette)
    
  } else if (heatmap_palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    n_max <- RColorBrewer::brewer.pal.info[heatmap_palette, "maxcolors"]
    colorRampPalette(RColorBrewer::brewer.pal(n = n_max, name = heatmap_palette))(n_breaks)
    
  } else {
    log_error(sample = "", step = "plot_heatmap",
              msg = glue::glue("`heatmap_palette` '{heatmap_palette}' not recognized."))
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 7: Prepare Annotations
  # ═══════════════════════════════════════════════════════════════════════════
  # pheatmap requires annotation dataframes to have rownames matching the
  # matrix colnames (col_annotation) or rownames (row_annotation) exactly.
  # make.names() is applied to col_id / row_id values to match the sanitisation
  # applied to the matrix axes in Section 2.
  
  col_annotation <- if (!is.null(metadata_col) && !is.null(col_annotations)) {
    metadata_col %>%
      dplyr::select(dplyr::all_of(c(col_id, col_annotations))) %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(col_id), make.names)) %>%
      dplyr::filter(.data[[col_id]] %in% colnames(expr_mat_scaled)) %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames(col_id) %>%
      as.data.frame() %>%
      dplyr::mutate(dplyr::across(where(is.factor), as.character))
  } else NULL
  
  row_annotation <- if (!is.null(metadata_row) && !is.null(row_annotations)) {
    metadata_row %>%
      dplyr::select(dplyr::all_of(c(row_id, row_annotations))) %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(row_id), make.names)) %>%
      dplyr::filter(.data[[row_id]] %in% rownames(expr_mat_scaled)) %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames(row_id) %>%
      as.data.frame() %>%
      dplyr::mutate(dplyr::across(where(is.factor), as.character))
  } else NULL
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 8: Annotation Color Palette
  # ═══════════════════════════════════════════════════════════════════════════
  # Builds ann_colors — a named list mapping each annotation column to a named
  # color vector — which pheatmap uses for annotation sidebar colors.
  #
  # annotation_palette is a proper function argument (not a global) so
  # plot_heatmap() is self-contained and testable outside the pipeline context.
  # NULL uses the built-in 20-color palette. If the user supplies fewer colors
  # than the maximum number of levels across all annotation columns, colors are
  # recycled with a warning rather than erroring — annotation sidebars with
  # repeated colors are still useful; a silent crash is not.
  #
  # Discrete columns  : one distinct base color per level
  # Numeric columns   : alpha-graduated shades of a single base color
  # Column annotations are processed before row annotations so the color
  # assignment order stays consistent across plots where row annotations may vary.
  
  CUSTOM_PALETTE <- c(
    "#00538A", "#C10020", "#007D34", "#FFB300", "#6A0DAD", "#FF6800",
    "#1B9E77", "#E7298A", "#4169E1", "#7F180D", "#2CA02C", "#F28E2B",
    "#803E75", "#93AA00", "#FF7A5C", "#17BECF", "#B32851", "#556B2F",
    "#9467BD", "#CEA262", "#008080", "#F6768E", "#A65628", "#EDC948",
    "#53377A", "#B8860B", "#E377C2", "#2F4F4F", "#B07AA1", "#817066",
    "#A0CBE8"
  )
  
  PUBLICATION_PALETTE <- c(
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
    "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
    "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4"
  )
  
  # Compute max levels needed across all annotation columns upfront so we can
  # warn once if the user palette is too short, rather than silently recycling.
  all_ann_levels <- c(
    if (!is.null(col_annotation)) unlist(lapply(col_annotation, function(x) length(unique(x)))) else integer(0),
    if (!is.null(row_annotation)) unlist(lapply(row_annotation, function(x) length(unique(x)))) else integer(0)
  )
  max_levels_needed <- if (length(all_ann_levels) > 0) sum(all_ann_levels) else 0
  
  base_colors <- if (!is.null(annotation_palette)) {
    if (length(annotation_palette) < max_levels_needed) {
      log_warn(sample = "", step = "plot_heatmap",
               msg = glue::glue(
                 "annotation_palette has {length(annotation_palette)} colors but ",
                 "{max_levels_needed} are needed across all annotation levels. Recycling."
               ))
      rep_len(annotation_palette, max_levels_needed)
    } else {
      annotation_palette
    }
  } else {
    if (max_levels_needed > length(CUSTOM_PALETTE)) {
      rep_len(CUSTOM_PALETTE, max_levels_needed)
    } else {
      CUSTOM_PALETTE
    }
  }
  
  ann_colors  <- list()
  color_index <- 1
  
  # col_list and row_list must be in the same order — col first, row second —
  # to match ann_types below. Previously these were reversed, which caused
  # wrong color assignments when an annotation column name appeared in both.
  col_list <- if (!is.null(col_annotation)) lapply(as.list(col_annotation), function(x) unique(as.character(x))) else list()
  row_list <- if (!is.null(row_annotation)) lapply(as.list(row_annotation), function(x) unique(as.character(x))) else list()
  ann_list <- c(col_list, row_list)
  
  col_types <- if (!is.null(col_annotation)) sapply(col_annotation, class) else c()
  row_types <- if (!is.null(row_annotation)) sapply(row_annotation, class) else c()
  ann_types <- c(col_types, row_types)   # col first, row second — must match ann_list order
  
  for (i in seq_along(ann_list)) {
    
    levels         <- sort(ann_list[[i]])
    n_levels       <- length(levels)
    ann_name       <- names(ann_list)[i]
    is_numeric_ann <- ann_types[ann_name] %in% c("numeric", "integer", "double")
    
    palette_colors <- if (!is_numeric_ann || n_levels == 1) {
      base_colors[color_index:(color_index + n_levels - 1)]
    } else {
      alphas <- seq(1 / n_levels, 1, length.out = n_levels)
      sapply(alphas, function(x) colorspace::adjust_transparency(col = base_colors[color_index], alpha = x))
    }
    
    names(palette_colors) <- levels
    ann_colors            <- c(ann_colors, list(palette_colors))
    names(ann_colors)[i]  <- ann_name
    color_index           <- color_index + n_levels
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 9: Gap Positions
  # ═══════════════════════════════════════════════════════════════════════════
  # Gap positions are cumulative group sizes computed from sorted group labels,
  # which must match the sort order used in col_order/row_order below.
  # The final gap is excluded — pheatmap draws a gap BETWEEN groups, not after
  # the last one, so the last cumulative size is always dropped.
  
  gaps_col <- NULL
  if (!is.null(col_gap_by) && !is.na(col_gap_by) && col_gap_by != "") {
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
  
  gaps_row <- NULL
  if (!is.null(row_gap_by) && !is.na(row_gap_by) && row_gap_by != "") {
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
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 10: Clustering and Ordering
  # ═══════════════════════════════════════════════════════════════════════════
  # Why do clustering manually rather than letting pheatmap handle it?
  # pheatmap's built-in clustering cannot cluster WITHIN groups separately
  # (the col_cluster_by = <column> case). By computing order vectors here and
  # passing a pre-ordered matrix with cluster_rows/cols = FALSE, we get full
  # control over within-group clustering while still using pheatmap for rendering.
  #
  # Distance: Euclidean — appropriate after z-scoring because all genes are
  #           on the same scale and Euclidean distance has a clear interpretation.
  # Linkage:  ward.D2   — minimises total within-cluster variance (sum of
  #           squared distances). Produces compact, roughly equal-sized clusters
  #           and is the most widely used method for expression heatmaps.
  #           ward.D2 requires squared Euclidean distances internally — hclust
  #           handles this correctly when method = "ward.D2" (not "ward").
  
  # ── Column ordering ──────────────────────────────────────────────────────
  
  if (!is.null(col_annotation) && col_cluster_by %in% colnames(col_annotation)) {
    
    col_order <- c()
    # Sort groups alphabetically — MUST match cumsum order used for gaps_col
    col_vars  <- col_annotation %>% dplyr::pull(col_cluster_by) %>% unique() %>% sort()
    
    for (col_var in col_vars) {
      samples <- col_annotation %>%
        tibble::rownames_to_column(col_id) %>%
        dplyr::filter(.data[[col_cluster_by]] == col_var) %>%
        dplyr::pull(col_id) %>%
        intersect(colnames(expr_mat_scaled))
      
      if (length(samples) > 1) {
        temp_mat  <- expr_mat_scaled[, samples, drop = FALSE]
        col_dist  <- stats::dist(x = t(temp_mat), method = "euclidean")
        col_clust <- stats::hclust(d = col_dist, method = "ward.D2")
        col_order <- c(col_order, colnames(temp_mat)[col_clust$order])
      } else if (length(samples) == 1) {
        col_order <- c(col_order, samples)
      }
      # length == 0: group in metadata but not in matrix — silently skip
    }
    
  } else if (col_cluster_by == "all") {
    col_dist  <- stats::dist(x = t(expr_mat_scaled), method = "euclidean")
    col_clust <- stats::hclust(d = col_dist, method = "ward.D2")
    col_order <- colnames(expr_mat_scaled)[col_clust$order]
    
  } else if (col_cluster_by == "alphabetical") {
    col_order <- sort(colnames(expr_mat_scaled))
    
  } else if (col_cluster_by == "none") {
    col_order <- colnames(expr_mat_scaled)
    
  } else {
    log_warn(sample = "", step = "plot_heatmap",
             msg = glue::glue("col_cluster_by '{col_cluster_by}' not recognized. Preserving input order."))
    col_order <- colnames(expr_mat_scaled)
  }
  
  # ── Row ordering ─────────────────────────────────────────────────────────
  
  if (!is.null(row_annotation) && row_cluster_by %in% colnames(row_annotation)) {
    
    row_order <- c()
    # Sort groups alphabetically — MUST match cumsum order used for gaps_row
    row_vars  <- row_annotation %>% dplyr::pull(row_cluster_by) %>% unique() %>% sort()
    
    for (row_var in row_vars) {
      genes <- row_annotation %>%
        tibble::rownames_to_column(row_id) %>%
        dplyr::filter(.data[[row_cluster_by]] == row_var) %>%
        dplyr::pull(row_id) %>%
        intersect(rownames(expr_mat_scaled))
      
      if (length(genes) > 1) {
        temp_mat  <- expr_mat_scaled[genes, , drop = FALSE]
        row_dist  <- stats::dist(x = temp_mat, method = "euclidean")
        row_clust <- stats::hclust(d = row_dist, method = "ward.D2")
        row_order <- c(row_order, rownames(temp_mat)[row_clust$order])
      } else if (length(genes) == 1) {
        row_order <- c(row_order, genes)
      }
      # length == 0: group in metadata but not in matrix — silently skip
    }
    
  } else if (row_cluster_by == "all") {
    row_dist  <- stats::dist(x = expr_mat_scaled, method = "euclidean")
    row_clust <- stats::hclust(d = row_dist, method = "ward.D2")
    row_order <- rownames(expr_mat_scaled)[row_clust$order]
    
  } else if (row_cluster_by == "alphabetical") {
    row_order <- sort(rownames(expr_mat_scaled))
    
  } else if (row_cluster_by == "none") {
    row_order <- rownames(expr_mat_scaled)
    
  } else {
    log_warn(sample = "", step = "plot_heatmap",
             msg = glue::glue("row_cluster_by '{row_cluster_by}' not recognized. Preserving input order."))
    row_order <- rownames(expr_mat_scaled)
  }
  
  reordered <- expr_mat_scaled[row_order, col_order]
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 11: Plot Aesthetics
  # ═══════════════════════════════════════════════════════════════════════════
  # Cell size logic:
  #   Ideal size  : fontsize + 5pt (breathing room around labels)
  #   Available   : ~6" wide x 8" tall on a standard page (432pt x 576pt)
  #   min() picks the ideal size when the matrix is small enough to fit, and
  #   shrinks cells proportionally when there are too many genes/samples.
  #   Labels are hidden automatically when cells fall below fontsize — there
  #   is no point rendering text that would be smaller than 1pt.
  
  fontsize        <- 10
  fontsize_number <- fontsize * 0.8
  angle_col       <- 45
  cell_width      <- min(fontsize + 5, (6 * 72) / ncol(reordered))
  cell_height     <- min(fontsize + 5, (8 * 72) / nrow(reordered))
  
  # Title: user-supplied string (word-wrapped at 40 chars) + auto dimension subtitle
  # so every heatmap is self-describing even when saved out of context
  n_rows   <- nrow(reordered)
  n_cols   <- ncol(reordered)
  subtitle <- glue::glue("({n_rows} rows \u00d7 {n_cols} columns)")
  
  plot_title <- if (!is.null(plot_title) && is.character(plot_title) && length(plot_title) == 1) {
    paste0(stringr::str_wrap(string = plot_title, width = 40), "\n", subtitle)
  } else {
    subtitle
  }
  
  # If label_genes is specified, show only those genes and blank the rest.
  # If not specified, show all labels when cells are large enough to be readable,
  # otherwise blank all labels (single space " " rather than "" so pheatmap
  # doesn't collapse row height for blank labels).
  labels_row <- if (!is.null(label_genes)) {
    dplyr::if_else(rownames(reordered) %in% make.names(label_genes),
                   stringr::str_trunc(rownames(reordered), width = 15),
                   " ")
  } else if (cell_height == fontsize + 5) {
    stringr::str_trunc(rownames(reordered), width = 15)
  } else {
    rep(" ", nrow(reordered))
  }
  
  labels_col <- if (cell_width == fontsize + 5) {
    stringr::str_trunc(colnames(reordered), width = 15)
  } else {
    rep(" ", ncol(reordered))
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 12: Generate Heatmap, Prepare Output Matrix, and Return
  # ═══════════════════════════════════════════════════════════════════════════
  
  ph <- pheatmap::pheatmap(
    mat               = reordered,
    color             = heatmap_colors,
    breaks            = breaks,
    annotation_row    = row_annotation,
    annotation_col    = col_annotation,
    annotation_colors = ann_colors,
    gaps_row          = gaps_row,
    gaps_col          = gaps_col,
    cellwidth         = cell_width,
    cellheight        = cell_height,
    show_rownames     = cell_height >= fontsize,
    show_colnames     = cell_width  >= fontsize,
    labels_row        = labels_row,
    labels_col        = labels_col,
    angle_col         = angle_col,
    fontsize          = fontsize,
    fontsize_row      = fontsize,
    fontsize_col      = fontsize,
    fontsize_number   = fontsize_number,
    main              = plot_title,
    border_color      = border_color,
    legend            = show_expr_legend,
    # Clustering disabled — all ordering handled manually in Section 10 above
    scale                 = "none",
    cluster_rows          = FALSE,
    cluster_cols          = FALSE,
    annotation_legend     = TRUE,
    annotation_names_row  = FALSE,
    annotation_names_col  = FALSE,
    width                 = NA,
    height                = NA,
    filename              = NA,   # Never save from within this function — caller handles saving
    silent                = TRUE)
  
  # Prepare output matrix — transpose wide matrices to fit Excel column limit.
  # Excel limits: 1,048,576 rows × 16,384 columns.
  # Wide matrices (many samples, few genes) exceed the column limit but fit
  # after transposition. Only transpose when it actually helps — if the matrix
  # is too large even after transposing, warn the user to use CSV instead.
  ph_mat   <- reordered
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
  
  log_info(sample = "", step = "plot_heatmap",
           msg = glue::glue("Heatmap generated: {nrow(reordered)} rows × {ncol(reordered)} cols | ",
                            "data_type = '{detected_type}' | palette = '{heatmap_palette}'"))
  
  # Return a named list so callers can access the pheatmap object, the
  # display-ordered matrix, and the detected data type independently:
  #   result$ph        → pheatmap object  (pass to grid.draw() for rendering)
  #   result$mat       → final matrix in display order (save to Excel or use downstream)
  #   result$data_type → detected/used data type (useful for sanity-checking)
  return(invisible(list(ph        = ph,
                        mat       = ph_mat,
                        data_type = detected_type)))
}


# ═══════════════════════════════════════════════════════════════════════════════
# PCA PLOT
# ═══════════════════════════════════════════════════════════════════════════════
# plot_pca() generates a PCA plot from an annotated expression data.frame.
#
# INPUT:
#   df        : data.frame with a row_id column (e.g. "SYMBOL") and one numeric
#               column per sample. Matches the format returned by process_counts().
#               Non-sample columns (gene_id, gene_biotype etc.) are dropped
#               automatically via intersection with metadata[[col_id]].
#               Must be pre-transformed — see data_type for supported input types.
#               Raw counts are accepted but log2(1+x) will be applied automatically.
#   metadata  : data.frame with one row per sample. Must contain a col_id column
#               (default "Sample_ID") whose values match colnames(df). All other
#               columns are used as grouping variables — one PCA scatter plot is
#               produced per column that has meaningful variation (see skip logic
#               in Section 7).
#
# row_id    : column in df containing gene identifiers. Used to set matrix
#             rownames (gene loadings in pca$rotation). Default "SYMBOL".
# col_id    : column in metadata containing sample identifiers. Used to identify
#             which df columns are samples. Default "Sample_ID".
# data_type : controls log transform and prcomp() scaling. Delegated to
#             detect_data_type() — same system as plot_heatmap(). Pass "auto"
#             to detect from the matrix, or override:
#               "counts"        — Raw integer counts; log2(1+x) applied, scale. = FALSE
#               "counts_log"    — Log-transformed (VST, rlog, log2 CPM); no transform, scale. = FALSE
#               "counts_scaled" — Already log-transformed AND z-scored; no transform, scale. = FALSE
#               "centered"      — Zero-centered scores (LFC, NES); no transform, scale. = TRUE
#               "bounded"       — Values in [0, 1] (beta, p-values); no transform, scale. = TRUE
# shape_by  : single metadata column name to map to point shape. NULL = no shape.
#             Warn if shape_by matches the current color variable. Warn and ignore
#             if shape_by has only 1 unique value.
# top_n     : number of most variable genes to retain before PCA. Genes with
#             exactly zero variance are removed first, then top_n by variance.
#             Invalid values (non-integer, negative, zero) fall back to all genes
#             with a log_warn. Ignored if nrow(expr_mat) <= top_n.
# label_samples : logical. TRUE (default) labels every sample with its col_id
#                 value via ggrepel. FALSE plots points only. Non-logical defaults
#                 to TRUE.
# title     : optional string prepended to each plot title above the variable name.
# filename  : base name for the output PDF. Non-alphanumeric characters are
#             sanitised to underscores. Final file is "PCA_Plot_<filename>.pdf".
# output_dir: directory to save the PDF. NULL = return plots only, no file saved.
#
# OUTPUT:
#   Returns invisibly: list(
#     pca   = prcomp result object  (loadings, sdev, x coordinates),
#     scree = ggplot scree plot,
#     plots = named list of ggplots, one per metadata variable
#   )
#   When output_dir is provided, saves a multi-page PDF:
#     Page 1         : scree plot (variance explained per PC)
#     Pages 2 to N+1 : one PCA scatter plot per metadata variable

plot_pca <- function(df,
                     metadata,
                     row_id        = "SYMBOL",
                     col_id        = "Sample_ID",
                     data_type     = "auto",
                     shape_by      = NULL,
                     top_n         = 500,
                     label_samples = TRUE,
                     title         = NULL,
                     filename      = "PCA_Plot",
                     output_dir    = NULL) {
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Validate Inputs
  # ═══════════════════════════════════════════════════════════════════════════
  # All structural checks fire here before any data is touched — fail fast with
  # a clear message rather than crashing deep inside prcomp or ggplot.
  
  # ── 1a. df must be a data.frame ──────────────────────────────────────────
  if (!is.data.frame(df)) {
    log_error(sample = "", step = "plot_pca",
              msg = glue::glue(
                "'df' must be a data.frame, not {class(df)[1]}. ",
                "Pass a data.frame with a '{row_id}' column and numeric sample columns. ",
                "process_counts() returns a data.frame in this format."
              ))
  }
  
  # ── 1b. row_id column must exist in df ──────────────────────────────────
  if (!row_id %in% colnames(df)) {
    log_error(sample = "", step = "plot_pca",
              msg = glue::glue(
                "row_id column '{row_id}' not found in df. ",
                "Available columns: {paste(colnames(df), collapse = ', ')}"
              ))
  }
  
  # ── 1c. metadata must be a data.frame with col_id column ────────────────
  if (!is.data.frame(metadata)) {
    log_error(sample = "", step = "plot_pca",
              msg = glue::glue("'metadata' must be a data.frame, not {class(metadata)[1]}."))
  }
  
  if (!col_id %in% colnames(metadata)) {
    log_error(sample = "", step = "plot_pca",
              msg = glue::glue(
                "'metadata' must contain a '{col_id}' column. ",
                "Available columns: {paste(colnames(metadata), collapse = ', ')}"
              ))
  }
  
  # ── 1d. Sample names must overlap between df and metadata ───────────────
  # A silent mismatch here would corrupt PCA coordinates — the $x matrix rows
  # would be joined to the wrong metadata rows. We take the intersection and
  # warn about any excluded samples so the user knows their data was trimmed.
  shared_samples <- intersect(metadata[[col_id]], colnames(df))
  only_in_df     <- setdiff(colnames(df)[colnames(df) != row_id], metadata[[col_id]])
  only_in_meta   <- setdiff(metadata[[col_id]], colnames(df))
  
  if (length(shared_samples) == 0) {
    log_error(sample = "", step = "plot_pca",
              msg = glue::glue(
                "No shared samples between colnames(df) and metadata[['{col_id}']]. ",
                "Check that sample names match exactly (case-sensitive)."
              ))
  }
  
  if (length(only_in_df) > 0) {
    log_warn(sample = "", step = "plot_pca",
             msg = glue::glue(
               "{length(only_in_df)} column(s) in df not found in metadata — excluded: ",
               "{paste(only_in_df, collapse = ', ')}"
             ))
  }
  
  if (length(only_in_meta) > 0) {
    log_warn(sample = "", step = "plot_pca",
             msg = glue::glue(
               "{length(only_in_meta)} sample(s) in metadata not found in df — excluded: ",
               "{paste(only_in_meta, collapse = ', ')}"
             ))
  }
  
  # ── 1e. Minimum sample requirement ──────────────────────────────────────
  # PCA is undefined for fewer than 3 samples — 2 samples produces a single PC
  # and a degenerate plot with no useful structure.
  if (length(shared_samples) < 3) {
    log_error(sample = "", step = "plot_pca",
              msg = glue::glue(
                "PCA requires at least 3 samples. Found {length(shared_samples)} after sample matching."
              ))
  }
  
  # ── 1f. shape_by column must exist in metadata ──────────────────────────
  if (!is.null(shape_by) && !shape_by %in% colnames(metadata)) {
    log_error(sample = "", step = "plot_pca",
              msg = glue::glue(
                "shape_by column '{shape_by}' not found in metadata. ",
                "Available: {paste(colnames(metadata), collapse = ', ')}"
              ))
  }
  
  # ── 1g. shape_by must have at least 2 unique values ─────────────────────
  # A constant shape aesthetic adds no information and produces a misleading
  # legend. Warn and ignore rather than error — the plots are still valid.
  if (!is.null(shape_by) && length(unique(metadata[[shape_by]])) < 2) {
    log_warn(sample = "", step = "plot_pca",
             msg = glue::glue(
               "shape_by column '{shape_by}' has only 1 unique value — shape aesthetic ignored."
             ))
    shape_by <- NULL
  }
  
  # ── 1h. top_n must be a positive integer ────────────────────────────────
  # Invalid top_n falls back to using all genes — the plot subtitle will
  # reflect the actual gene count used.
  valid_top_n <- is.numeric(top_n) && length(top_n) == 1 &&
    !is.na(top_n) && top_n > 0 && top_n == as.integer(top_n)
  
  if (!valid_top_n) {
    log_warn(sample = "", step = "plot_pca",
             msg = glue::glue(
               "Invalid top_n = '{top_n}'. Must be a positive integer. Using all genes."
             ))
    top_n <- Inf
  }
  
  # ── 1i. label_samples — non-logical defaults to TRUE ────────────────────
  if (!is.logical(label_samples) || is.na(label_samples)) {
    log_warn(sample = "", step = "plot_pca",
             msg = glue::glue("label_samples must be TRUE or FALSE. Defaulting to TRUE."))
    label_samples <- TRUE
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Prepare Matrix
  # ═══════════════════════════════════════════════════════════════════════════
  # Subset df to shared samples + row_id, convert to numeric matrix with gene
  # IDs as rownames. Non-numeric sample columns are caught here before they
  # can silently coerce to NA inside as.matrix().
  # Alignment is critical — prcomp() $x rows correspond to colnames(expr_mat)
  # order, and we inner_join on col_id in Section 7, so both must match.
  
  metadata <- metadata %>%
    dplyr::filter(.data[[col_id]] %in% shared_samples) %>%
    dplyr::arrange(match(.data[[col_id]], shared_samples))
  
  non_numeric <- shared_samples[!sapply(df[shared_samples], is.numeric)]
  if (length(non_numeric) > 0) {
    log_error(sample = "", step = "plot_pca",
              msg = glue::glue(
                "The following sample columns are not numeric: {paste(non_numeric, collapse = ', ')}. ",
                "Check that col_id / metadata correctly identifies sample columns."
              ))
  }
  
  expr_mat <- df %>%
    dplyr::select(dplyr::all_of(c(row_id, shared_samples))) %>%
    tibble::column_to_rownames(row_id) %>%
    as.matrix()
  
  # ── Duplicate row_id values ──────────────────────────────────────────────
  # When duplicates exist, keep the row with the highest total absolute
  # expression — same strategy as plot_heatmap() Section 2.
  if (any(duplicated(rownames(expr_mat)))) {
    log_warn(sample = "", step = "plot_pca",
             msg = "Duplicate gene symbols detected. Retaining highest-expressing copy per gene.")
    expr_mat <- expr_mat %>%
      as.data.frame() %>%
      tibble::rownames_to_column(row_id) %>%
      dplyr::mutate(total_expr = rowSums(dplyr::across(-dplyr::all_of(row_id), abs))) %>%
      dplyr::group_by(.data[[row_id]]) %>%
      dplyr::slice_max(order_by = total_expr, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      tibble::column_to_rownames(row_id) %>%
      dplyr::select(-total_expr) %>%
      as.matrix()
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Data Type Detection
  # ═══════════════════════════════════════════════════════════════════════════
  # Delegated to detect_data_type() — a shared helper used by plot_heatmap()
  # and plot_pca(). Detection logic, thresholds, and log messages are defined
  # once there and applied consistently across both functions.
  
  detected_type <- detect_data_type(mat = expr_mat, data_type = data_type)
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Log Transform
  # ═══════════════════════════════════════════════════════════════════════════
  # Only raw counts require a log transform. All other types are already on a
  # log or otherwise stabilized scale and are used as-is.
  #
  # Why log2(1+x) rather than log2(x)?
  # Genes with zero counts produce -Inf under log2(x). Adding 1 (pseudocount)
  # shifts zeros to 0 on the log scale without distorting relative relationships
  # between expressed genes.
  #
  # Why not VST here?
  # VST requires DESeq2 and a full count model — not appropriate as an internal
  # transform inside a plot function. If VST-quality transformation is needed,
  # apply it upstream and pass data_type = "counts_log".
  
  if (detected_type == "counts") {
    log_info(sample = "", step = "plot_pca",
             msg = "Applying log2(1+x) transform to raw counts.")
    expr_mat <- log2(1 + expr_mat)
  } else {
    log_info(sample = "", step = "plot_pca",
             msg = glue::glue("Skipping log transform for data_type = '{detected_type}'."))
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5: Feature Selection
  # ═══════════════════════════════════════════════════════════════════════════
  # PCA on all genes is both computationally wasteful and statistically noisy —
  # genes with no variation contribute only noise to the PC directions.
  # We select the top_n most variable genes to focus the PCA on the signal
  # that actually differentiates samples.
  #
  # Why filter exactly-zero variance first?
  # Genes with variance == 0 are identical across every sample — they carry no
  # information about sample differences by definition. We remove them before
  # ranking so top_n slots are always filled by informative genes.
  #
  # Why not a near-zero threshold (e.g. < 1e-10)?
  # Near-zero thresholds are arbitrary and risk removing genuinely low-variance
  # genes that may still be biologically informative. Exactly zero is the only
  # principled hard cutoff.
  
  gene_vars  <- matrixStats::rowVars(expr_mat)
  n_zero_var <- sum(gene_vars == 0)
  
  if (n_zero_var > 0) {
    log_info(sample = "", step = "plot_pca",
             msg = glue::glue("Removed {n_zero_var} gene(s) with exactly zero variance."))
  }
  
  expr_mat <- expr_mat[gene_vars > 0, , drop = FALSE]
  
  if (nrow(expr_mat) == 0) {
    log_error(sample = "", step = "plot_pca",
              msg = "No genes remain after removing zero-variance genes. Check your input data.")
  }
  
  if (nrow(expr_mat) > top_n) {
    gene_vars_filtered <- matrixStats::rowVars(expr_mat)
    top_genes          <- order(gene_vars_filtered, decreasing = TRUE)[seq_len(top_n)]
    expr_mat           <- expr_mat[top_genes, , drop = FALSE]
    genes_used_label   <- glue::glue("Top {top_n} most variable genes")
    log_info(sample = "", step = "plot_pca",
             msg = glue::glue("Selected top {top_n} most variable genes for PCA."))
  } else {
    genes_used_label <- glue::glue("All {nrow(expr_mat)} genes used")
    log_info(sample = "", step = "plot_pca",
             msg = glue::glue("nrow(df) = {nrow(expr_mat)} <= top_n = {top_n}. Using all genes."))
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 6: PCA
  # ═══════════════════════════════════════════════════════════════════════════
  # prcomp() expects samples as rows and genes as columns — the opposite of the
  # standard bioinformatics genes × samples matrix. t() is mandatory here.
  #
  # center = TRUE always: subtracts the mean of each gene across samples before
  # computing PCs. Without centering, PC1 captures mean expression level rather
  # than variation, which is never the biological question.
  #
  # scale. is driven by detected_type:
  #   FALSE for counts / counts_log / counts_scaled — after log transform, gene
  #         variances are meaningful and should be preserved. High-variance genes
  #         genuinely vary more across samples and should contribute more to PCs.
  #   TRUE  for centered / bounded — features may be on incommensurable scales
  #         (e.g. LFC values from -10 to +10 vs NES from -3 to +3). Scaling
  #         ensures each feature contributes equally regardless of its range.
  
  do_scale <- detected_type %in% c("centered", "bounded")
  
  log_info(sample = "", step = "plot_pca",
           msg = glue::glue(
             "Running prcomp(): center = TRUE, scale. = {do_scale} ",
             "(data_type = '{detected_type}')."
           ))
  
  pca_results <- stats::prcomp(x = t(expr_mat), center = TRUE, scale. = do_scale)
  
  # ── Variance explained per PC ───────────────────────────────────────────
  # importance[2, ] = proportion of variance; multiply by 100 for percentages.
  # We extract all PCs — the scree plot and axis labels both draw from this vector.
  pct_var     <- round(100 * summary(pca_results)$importance[2, ], 1)
  cum_var     <- round(100 * summary(pca_results)$importance[3, ], 1)
  n_pcs_total <- length(pct_var)
  
  log_info(sample = "", step = "plot_pca",
           msg = glue::glue(
             "PC1 = {pct_var[1]}%, PC2 = {pct_var[2]}%, ",
             "cumulative PC1+PC2 = {cum_var[2]}%."
           ))
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 7: Scree Plot
  # ═══════════════════════════════════════════════════════════════════════════
  # The scree plot is placed FIRST in the PDF so the user can assess how much
  # variance PC1/PC2 capture before interpreting the scatter plots. A PC1+PC2
  # sum of 30% tells a very different story from 80%.
  #
  # All PCs are shown up to a cap of 50 — bulk RNA-seq datasets rarely exceed
  # 50 samples so this cap is almost never hit in practice. Showing all PCs
  # (rather than an arbitrary 20) gives the full elbow picture.
  # The cumulative variance line helps identify how many PCs capture the bulk
  # of the signal (commonly used 80% elbow threshold).
  
  n_pcs_scree <- min(50, n_pcs_total)
  scree_df    <- data.frame(
    PC         = factor(paste0("PC", seq_len(n_pcs_scree)),
                        levels = paste0("PC", seq_len(n_pcs_scree))),
    Variance   = pct_var[seq_len(n_pcs_scree)],
    Cumulative = cum_var[seq_len(n_pcs_scree)]
  )
  
  p_scree <- ggplot2::ggplot(scree_df, ggplot2::aes(x = PC)) +
    ggplot2::geom_col(ggplot2::aes(y = Variance),
                      fill = "#4DBBD5", color = "white", width = 0.7) +
    ggplot2::geom_line(ggplot2::aes(y = Cumulative, group = 1),
                       color = "#E64B35", linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(y = Cumulative),
                        color = "#E64B35", size = 2) +
    theme_publication() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      title    = if (!is.null(title)) glue::glue("{title} \u2014 Scree Plot") else "Scree Plot",
      subtitle = glue::glue("Bars: variance per PC   |   Red line: cumulative variance   |   {genes_used_label}"),
      x        = NULL,
      y        = "Variance Explained (%)"
    )
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 8: PCA Scatter Plots — One Per Metadata Variable
  # ═══════════════════════════════════════════════════════════════════════════
  # Merge PCA coordinates ($x matrix) with metadata so each sample row carries
  # both its PC scores and its grouping variables for aesthetic mapping.
  #
  # $x rows are already in the same order as colnames(expr_mat) — which we
  # aligned to metadata[[col_id]] in Section 2 — so the join is safe.
  
  pca_df <- pca_results$x %>%
    as.data.frame() %>%
    tibble::rownames_to_column(col_id) %>%
    dplyr::inner_join(metadata, by = col_id)
  
  # ── 8a. Resolve which variables to plot ─────────────────────────────────
  # Use every metadata column except col_id — one scatter plot per variable.
  group_vars <- base::setdiff(colnames(metadata), col_id)
  
  # ── 8b. Internal discrete color palette ─────────────────────────────────
  # Defined here — not pulled from the global environment. Makes plot_pca()
  # self-contained regardless of what the calling script has defined.
  
  CUSTOM_PALETTE <- c(
    "#00538A", "#C10020", "#007D34", "#FFB300", "#6A0DAD", "#FF6800",
    "#1B9E77", "#E7298A", "#4169E1", "#7F180D", "#2CA02C", "#F28E2B",
    "#803E75", "#93AA00", "#FF7A5C", "#17BECF", "#B32851", "#556B2F",
    "#9467BD", "#CEA262", "#008080", "#F6768E", "#A65628", "#EDC948",
    "#53377A", "#B8860B", "#E377C2", "#2F4F4F", "#B07AA1", "#817066",
    "#A0CBE8"
  )
  
  PUBLICATION_PALETTE <- c(
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
    "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
    "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4"
  )
  
  # ── 8c. Fixed axis labels — PC1/PC2 variance is the same for every plot ─
  x_lab <- paste0("PC1: ", pct_var[1], "% variance")
  y_lab <- paste0("PC2: ", pct_var[2], "% variance")
  
  # ── 8d. ggrepel seed — set once, reproducible across all plots ──────────
  set.seed(1234)
  
  all_plots <- list()
  
  for (var in group_vars) {
    
    # ── Skip uninformative columns ─────────────────────────────────────────
    # Constant column (all values identical): no color variation, plot is useless.
    # All-unique character/factor: every sample gets its own color — equivalent
    # to labeling by col_id and adds no grouping information.
    # All-unique numeric: treated as continuous gradient (e.g. age, BMI) — NOT
    # skipped, because the gradient itself is informative.
    col_vals <- metadata[[var]]
    n_unique <- length(unique(col_vals))
    is_num   <- is.numeric(col_vals)
    
    if (n_unique < 2) {
      log_warn(sample = "", step = "plot_pca",
               msg = glue::glue("Skipping '{var}': all values are identical."))
      next
    }
    
    if (n_unique == nrow(metadata) && !is_num) {
      log_warn(sample = "", step = "plot_pca",
               msg = glue::glue(
                 "Skipping '{var}': all values are unique and column is character/factor. ",
                 "No meaningful grouping to display."
               ))
      next
    }
    
    # ── Warn if shape_by duplicates this color variable ───────────────────
    if (!is.null(shape_by) && shape_by == var) {
      log_warn(sample = "", step = "plot_pca",
               msg = glue::glue(
                 "shape_by ('{shape_by}') matches current color variable '{var}'. ",
                 "Both aesthetics will encode the same information."
               ))
    }
    
    # ── Base aesthetic mapping ────────────────────────────────────────────
    # .data[[]] notation allows programmatic column access regardless of the
    # actual column name — the same pattern used throughout plot_volcano().
    base_aes <- if (!is.null(shape_by)) {
      ggplot2::aes(x     = PC1,
                   y     = PC2,
                   color = .data[[var]],
                   shape = .data[[shape_by]])
    } else {
      ggplot2::aes(x     = PC1,
                   y     = PC2,
                   color = .data[[var]])
    }
    
    p <- ggplot2::ggplot(data = pca_df, mapping = base_aes) +
      ggplot2::geom_point(size = 3) +
      theme_publication() +
      ggplot2::labs(
        color    = var,
        shape    = shape_by,
        x        = x_lab,
        y        = y_lab,
        title    = if (!is.null(title)) glue::glue("{title} \u2014 {var}") else var,
        subtitle = genes_used_label
      )
    
    # ── Color scale: continuous (numeric) vs discrete (character/factor) ──
    # is.numeric() is the sole criterion — no arbitrary unique-count threshold.
    # If the user stored a variable as numeric, we treat it as continuous.
    # If they stored it as character or factor, we treat it as discrete.
    # This respects the user's explicit data type choices without second-guessing.
    if (is_num) {
      p <- p + ggplot2::scale_color_viridis_c(option = "viridis")
    } else {
      disc_colors        <- CUSTOM_PALETTE[seq_len(n_unique)]
      names(disc_colors) <- as.character(unique(col_vals))
      p <- p + ggplot2::scale_color_manual(values = disc_colors)
    }
    
    # ── Sample labels (ggrepel) ───────────────────────────────────────────
    # label_samples = TRUE labels every sample with its col_id value.
    # label_samples = FALSE omits all labels — cleaner for large cohorts where
    # labels would overlap and become unreadable.
    # max.overlaps = Inf ensures no label is silently dropped even in dense plots.
    if (isTRUE(label_samples)) {
      p <- p + ggrepel::geom_text_repel(
        mapping            = ggplot2::aes(label = .data[[col_id]]),
        size               = 3,
        show.legend        = FALSE,
        max.overlaps       = Inf,
        box.padding        = 0.5,
        point.padding      = 0.2,
        min.segment.length = 0,
        segment.size       = 0.3,
        segment.color      = "grey50",
        segment.curvature  = -0.3,
        segment.ncp        = 30,
        segment.angle      = 20
      )
    }
    
    all_plots[[var]] <- p
    
    log_info(sample = "", step = "plot_pca",
             msg = glue::glue(
               "Plotted variable '{var}' ",
               "({if (is_num) 'continuous' else 'discrete'}, {n_unique} unique values)."
             ))
  }
  
  # Guard: if every variable was skipped, nothing to save
  if (length(all_plots) == 0) {
    log_warn(sample = "", step = "plot_pca",
             msg = "No plottable metadata variables found. Returning PCA results only.")
    return(invisible(list(pca = pca_results, scree = p_scree, plots = list())))
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 9: Save and Return
  # ═══════════════════════════════════════════════════════════════════════════
  # cairo_pdf is used over ggsave because it natively supports multi-page output
  # and handles font embedding and transparency correctly — critical for
  # publication-quality PDFs.
  #
  # Page order: scree plot first (variance context), then one scatter per variable.
  # This mirrors the logical reading order — understand the overall PCA structure
  # before interpreting individual group colorings.
  #
  # Filename sanitisation: replace non-alphanumeric runs with underscores, strip
  # leading/trailing underscores. Matches the same pattern in plot_ma() and
  # plot_volcano().
  #
  # Returns invisibly — the function's primary side effect is the PDF; the return
  # value is for programmatic use:
  #   $pca   — raw prcomp object: loadings ($rotation), coordinates ($x), sdev ($sdev)
  #   $scree — scree ggplot (renderable or saveable independently)
  #   $plots — named list of scatter ggplots, one per metadata variable
  
  if (!is.null(output_dir)) {
    
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    safe_filename <- gsub("[^[:alnum:]]+", "_", filename) |>
      gsub(pattern = "(^_+|_+$)", replacement = "", x = _)
    
    output_file <- file.path(output_dir, paste0("PCA_Plot_", safe_filename, ".pdf"))
    
    grDevices::cairo_pdf(filename = output_file,
                         width    = 8,
                         height   = 11.5,
                         onefile  = TRUE)
    
    print(p_scree)
    for (var in names(all_plots)) print(all_plots[[var]])
    
    grDevices::dev.off()
    
    log_info(sample = "", step = "plot_pca",
             msg = glue::glue(
               "PCA PDF saved to '{output_file}' ",
               "({1 + length(all_plots)} pages: 1 scree + {length(all_plots)} scatter plots)."
             ))
    
  } else {
    log_info(sample = "", step = "plot_pca",
             msg = "No output_dir provided — returning plots only, skipping save.")
  }
  
  return(invisible(list(
    pca   = pca_results,
    scree = p_scree,
    plots = all_plots
  )))
}

# ═══════════════════════════════════════════════════════════════════════════════
# DOT PLOT
# ═══════════════════════════════════════════════════════════════════════════════
# plot_dot() generates a horizontal dot plot from a pre-wrangled data.frame.
# All data preparation (filtering, ranking, label construction) is the caller's
# responsibility — plot_dot() is a pure ggplot assembler that handles only
# display decisions: y-axis text sizing, x-axis limits, NA row suppression,
# color assignment, size breaks, and optional saving.
#
# Dot size range is fixed internally at c(1, 5) mm and is intentionally not
# exposed as a parameter — a consistent range ensures dot sizes are visually
# comparable across all plots in the pipeline regardless of the underlying
# gene count range.
#
# INPUT:
#   df                  : data.frame, fully prepped and sorted in the desired
#                         display order. Rows with NA or "" in y_var are silently
#                         replaced with " " so ggplot retains them as blank space
#                         rather than dropping them and collapsing the panel.
#                         The y-axis order strictly follows the row order of df —
#                         sort df before passing if a specific order is needed.
#   x_var               : column name (string) mapped to the x-axis (dot position).
#   y_var               : column name (string) mapped to the y-axis (dot labels).
#   x_label             : x-axis title. NULL (default) uses x_var as the label.
#   y_label             : y-axis title. NULL (default) uses y_var as the label.
#   color_var           : column name (string) for dot outline color mapping.
#                         NULL = no color mapping; ggplot default color used.
#   color_colors        : named character vector mapping color levels to colors.
#                         NULL (default) = auto-assigned from CUSTOM_PALETTE.
#                         Two-level color with levels "Upregulated"/"Downregulated"
#                         always uses the standard pipeline colors regardless of
#                         color_colors, unless color_colors is explicitly provided.
#   color_legend_title  : legend group title for the color aesthetic. NULL = ggplot
#                         uses the color column name.
#   color_legend_labels : named character vector relabeling color legend entries.
#                         NULL = ggplot uses the level names as-is.
#   fill_var            : column name (string) for dot fill color mapping.
#                         NULL = no fill mapping; ggplot default fill used.
#                         Typically set to the same column as color_var so dots
#                         have matching outline and fill (e.g. both "Direction").
#   fill_colors         : named character vector mapping fill levels to colors.
#                         NULL (default) = auto-assigned from CUSTOM_PALETTE.
#                         Two-level fill with levels "Upregulated"/"Downregulated"
#                         always uses the standard pipeline colors regardless of
#                         fill_colors, unless fill_colors is explicitly provided.
#   fill_legend_title   : legend group title for the fill aesthetic. NULL = ggplot
#                         uses the fill column name.
#   fill_legend_labels  : named character vector relabeling fill legend entries.
#                         NULL = ggplot uses the level names as-is.
#   size_var            : column name (string) for dot size mapping. NULL = no
#                         size mapping; all dots drawn at fixed default size.
#                         Size breaks are auto-computed from the data as round
#                         integers within the actual data range. Size range is
#                         fixed at c(1, 5) mm for cross-plot comparability.
#   size_legend_title   : legend group title for the size aesthetic. NULL = ggplot
#                         uses the size column name.
#   alpha_var           : column name (string) for alpha mapping. NULL = no alpha
#                         mapping; all dots fully opaque. The caller is responsible
#                         for pre-computing any derived alpha column (e.g.
#                         padj <= 0.05) before passing df. plot_dot() maps whatever
#                         values exist in this column to alpha without assumptions.
#   alpha_values        : named numeric vector mapping alpha column levels to
#                         opacity values e.g. c("TRUE" = 1, "FALSE" = 0.35).
#                         Required when alpha_var is provided; ignored otherwise.
#   alpha_legend_title  : legend group title for the alpha aesthetic. NULL = ggplot
#                         uses the alpha column name.
#   alpha_legend_labels : named character vector relabeling alpha legend entries.
#                         NULL = ggplot uses the level names as-is.
#   x_min               : numeric lower bound for the x-axis. NULL (default) =
#                         min(0, min(df[[x_var]], na.rm = TRUE)) so negative dots
#                         are always shown and positive-only data anchors at zero.
#                         Pass x_min = 0 to force the axis to start at zero even
#                         when negative values are present.
#   title               : plot title string. Default "" (no title).
#   filename            : base name for the output PDF. Non-alphanumeric characters
#                         are sanitised to underscores. Default "Dot_Plot".
#   output_dir          : directory to save the PDF. NULL (default) = return the
#                         ggplot object only, no file written.
#
# OUTPUT:
#   Returns invisibly: the ggplot object p.
#   When output_dir is provided, saves a single-page PDF via ggsave +
#   cairo_pdf and logs the output path.

plot_dot <- function(df,
                     x_var,
                     y_var,
                     x_label             = NULL,
                     y_label             = NULL,
                     color_var           = NULL,
                     color_colors        = NULL,
                     color_legend_title  = NULL,
                     color_legend_labels = NULL,
                     fill_var            = NULL,
                     fill_colors         = NULL,
                     fill_legend_title   = NULL,
                     fill_legend_labels  = NULL,
                     size_var            = NULL,
                     size_legend_title   = NULL,
                     alpha_var           = NULL,
                     alpha_values        = NULL,
                     alpha_legend_title  = NULL,
                     alpha_legend_labels = NULL,
                     x_min               = NULL,
                     title               = "",
                     filename            = "Dot_Plot",
                     output_dir          = NULL) {
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Validate Inputs
  # ═══════════════════════════════════════════════════════════════════════════
  # Mandatory columns are checked upfront so the function fails immediately
  # with a clear message rather than crashing inside ggplot with a cryptic error.
  
  if (!is.data.frame(df)) {
    log_error(sample = "", step = "plot_dot",
              msg = glue::glue("'df' must be a data.frame, not {class(df)[1]}."))
  }
  
  mandatory_cols <- c(x_var, y_var)
  optional_cols  <- c(color_var, fill_var, size_var, alpha_var)
  cols_to_check  <- c(mandatory_cols, optional_cols)
  missing_cols   <- cols_to_check[!cols_to_check %in% colnames(df)]
  
  if (length(missing_cols) > 0) {
    log_error(sample = "", step = "plot_dot",
              msg = glue::glue(
                "Column(s) not found in df: {paste(missing_cols, collapse = ', ')}. ",
                "Available: {paste(colnames(df), collapse = ', ')}"
              ))
  }
  
  if (!is.null(alpha_var) && is.null(alpha_values)) {
    log_error(sample = "", step = "plot_dot",
              msg = glue::glue(
                "'alpha_values' must be provided when 'alpha_var' is specified. ",
                "e.g. alpha_values = c('TRUE' = 1, 'FALSE' = 0.35)"
              ))
  }
  
  if (!is.null(output_dir) && !is.character(output_dir)) {
    log_error(sample = "", step = "plot_dot",
              msg = glue::glue(
                "'output_dir' must be a character string or NULL, not {class(output_dir)[1]}."
              ))
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Prepare Data
  # ═══════════════════════════════════════════════════════════════════════════
  
  # ── 2a. Resolve axis labels ───────────────────────────────────────────────
  # %||% returns the left side if non-NULL, otherwise the right side (rlang).
  x_label <- x_label %||% x_var
  y_label <- y_label %||% y_var
  
  # ── 2b. Replace NA / blank y labels with unique placeholders ─────────────
  # ggplot2 drops rows where the y aesthetic is NA. Replacing with unique
  # placeholders (" ", "  ", "   " etc.) forces ggplot to keep each empty
  # row as distinct blank space rather than collapsing them into one.
  empty_mask  <- is.na(df[[y_var]]) | trimws(as.character(df[[y_var]])) == ""
  df[[y_var]] <- as.character(df[[y_var]])
  df[[y_var]][empty_mask] <- sapply(seq_len(sum(empty_mask)), 
                                    function(i) strrep(" ", i))
  
  # ── 2c. Lock y-axis order to df row order ────────────────────────────────
  # factor() with levels = unique() prevents ggplot from re-sorting the y-axis
  # alphabetically. The caller controls display order by sorting df before
  # passing — plot_dot() respects that order exactly.
  df[[y_var]] <- factor(df[[y_var]], levels = unique(df[[y_var]]))
  
  # ── 2d. Compute x_min ────────────────────────────────────────────────────
  # Default anchors at the smaller of 0 or the minimum data value so:
  #   - negative dots always have a visible left boundary
  #   - positive-only data anchors cleanly at zero
  # expansion(mult = c(0.1, 0.05)) adds breathing room on both sides so the
  # leftmost dot never sits flush against the y-axis.
  x_min <- x_min %||% min(0, min(df[[x_var]], na.rm = TRUE))
  
  # ── 2e. Detect x-axis label format ───────────────────────────────────────
  # Compute breaks first, then determine accuracy from the break values.
  # Integer breaks display without decimals (accuracy = 1). Non-integer breaks
  # use accuracy = NULL so scales auto-detects the minimum decimal places needed
  # to make all labels distinct — works correctly for both large ranges (NES)
  # and small fractions (GeneRatio).
  x_range    <- c(x_min, max(df[[x_var]], na.rm = TRUE))
  x_breaks   <- scales::pretty_breaks(n = 4)(x_range)
  x_accuracy <- if (all(x_breaks == floor(x_breaks))) 1 else NULL
  left_expand <- if (x_min < 0) 0.1 else 0
  
  # ── 2f. Compute y-axis text size ─────────────────────────────────────────
  # Long pathway names need smaller text to avoid overlapping. Scale is based
  # on the longest label in the column after NA replacement above.
  max_label_len <- max(nchar(as.character(df[[y_var]])), na.rm = TRUE)
  y_size <- dplyr::case_when(
    max_label_len > 50 ~ 6,
    max_label_len > 35 ~ 7,
    max_label_len > 25 ~ 8,
    TRUE               ~ 10
  )
  
  # ── 2g. Resolve color colors ─────────────────────────────────────────────
  # Priority: caller-supplied color_colors > standard Up/Down palette > CUSTOM_PALETTE.
  # The Up/Down shortcut ensures consistent direction colors across all pipeline plots.
  if (!is.null(color_var) && is.null(color_colors)) {
    color_levels <- unique(df[[color_var]])
    updown_cols  <- c("Upregulated" = "#E69F00", "Downregulated" = "#56B4E9")
    
    if (all(color_levels %in% names(updown_cols))) {
      color_colors <- updown_cols
    } else {
      color_colors <- setNames(
        CUSTOM_PALETTE[seq_along(color_levels)],
        color_levels
      )
    }
  }
  
  # ── 2h. Resolve fill colors ───────────────────────────────────────────────
  # Same priority logic as color_colors above. When color_var and fill_var
  # point to the same column (e.g. both "Direction"), fill_colors can be left
  # NULL and will auto-resolve to the same palette as color_colors.
  if (!is.null(fill_var) && is.null(fill_colors)) {
    fill_levels <- unique(df[[fill_var]])
    updown_cols <- c("Upregulated" = "#E69F00", "Downregulated" = "#56B4E9")
    
    if (all(fill_levels %in% names(updown_cols))) {
      fill_colors <- updown_cols
    } else {
      fill_colors <- setNames(
        CUSTOM_PALETTE[seq_along(fill_levels)],
        fill_levels
      )
    }
  }
  
  # ── 2i. Compute size breaks ───────────────────────────────────────────────
  # pretty() generates round numbers for the size legend; filtered to the
  # actual data range and integers only so the legend is clean and meaningful.
  # Dot size range is fixed at c(1, 5) mm — not exposed as a parameter so
  # sizes remain visually comparable across all plots in the pipeline.
  # Only computed when size_var is provided.
  size_breaks <- NULL
  if (!is.null(size_var)) {
    size_vals   <- df[[size_var]][!is.na(df[[size_var]])]
    size_breaks <- pretty(size_vals, n = 4)
    size_breaks <- size_breaks[size_breaks >= min(size_vals) &
                                 size_breaks <= max(size_vals) &
                                 size_breaks %% 1 == 0]
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Build Plot
  # ═══════════════════════════════════════════════════════════════════════════
  
  # ── 3a. Base plot ─────────────────────────────────────────────────────────
  # aes() is built dynamically so that color, fill, size, and alpha are only
  # mapped when the caller explicitly provides them — unmapped aesthetics
  # inherit ggplot defaults and do not produce spurious legend entries.
  base_aes <- ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]])
  
  if (!is.null(color_var)) base_aes <- utils::modifyList(base_aes, ggplot2::aes(color = .data[[color_var]]))
  if (!is.null(fill_var))  base_aes <- utils::modifyList(base_aes, ggplot2::aes(fill  = .data[[fill_var]]))
  if (!is.null(size_var))  base_aes <- utils::modifyList(base_aes, ggplot2::aes(size  = .data[[size_var]]))
  if (!is.null(alpha_var)) base_aes <- utils::modifyList(base_aes, ggplot2::aes(alpha = .data[[alpha_var]]))
  
  p <- ggplot2::ggplot(data = df, mapping = base_aes) +
    
    # ── 3b. Dots ─────────────────────────────────────────────────────────────
    ggplot2::geom_point(na.rm = TRUE) +
    
    # ── 3c. X axis ───────────────────────────────────────────────────────────
    # expansion left = 0.1 gives breathing room between y-axis and leftmost dot.
    # expansion right = 0.05 gives a small margin after the rightmost dot.
    # Upper limit is left as NA so ggplot auto-scales to the data maximum.
    ggplot2::scale_x_continuous(
      limits = c(x_min, NA),
      expand = ggplot2::expansion(mult = c(left_expand, 0.05)),
      breaks = x_breaks,
      labels = scales::label_number(accuracy = x_accuracy)
    ) +
    
    # ── 3d. Labels and title ─────────────────────────────────────────────────
    ggplot2::labs(x = x_label, y = y_label, title = title) +
    
    # ── 3e. Coordinate system ────────────────────────────────────────────────
    # clip = "off" ensures dots or labels that slightly overflow the plot panel
    # boundary are still rendered rather than clipped invisibly.
    ggplot2::coord_cartesian(clip = "off") +
    
    # ── 3f. Y-axis text size ─────────────────────────────────────────────────
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = y_size),
      axis.text.x = if (is.null(x_accuracy)) {
        ggplot2::element_text(angle = 45, hjust = 1)
      } else {
        ggplot2::element_text(angle = 0)
      }
    )
  
  # ── 3g. Color scale ──────────────────────────────────────────────────────
  # Only added when color mapping is active. na.translate = FALSE excludes NA
  # levels from the legend — NA dots are still plotted but not labeled.
  if (!is.null(color_var)) {
    p <- p +
      ggplot2::scale_color_manual(
        values       = color_colors,
        name         = color_legend_title,
        labels       = color_legend_labels %||% ggplot2::waiver(),
        na.translate = FALSE
      ) +
      ggplot2::guides(
        color = ggplot2::guide_legend(override.aes = list(shape = 16, size = 4))
      )
  }
  
  # ── 3h. Fill scale ───────────────────────────────────────────────────────
  # Only added when fill mapping is active. na.translate = FALSE excludes NA
  # levels from the legend — NA dots are still plotted but not labeled.
  if (!is.null(fill_var)) {
    p <- p +
      ggplot2::scale_fill_manual(
        values       = fill_colors,
        name         = fill_legend_title,
        labels       = fill_legend_labels %||% ggplot2::waiver(),
        na.translate = FALSE
      ) +
      ggplot2::guides(
        fill = ggplot2::guide_legend(override.aes = list(shape = 16, size = 4))
      )
  }
  
  # ── 3i. Size scale ───────────────────────────────────────────────────────
  # Only added when size mapping is active. Range fixed at c(1, 5) mm —
  # consistent across all pipeline dot plots so sizes are visually comparable.
  # Breaks auto-computed in Section 2i: round integers within actual data range.
  if (!is.null(size_var)) {
    p <- p +
      ggplot2::scale_size(
        range  = c(1, 5),
        breaks = size_breaks,
        name   = size_legend_title
      )
  }
  
  # ── 3j. Alpha scale ──────────────────────────────────────────────────────
  # Only added when alpha mapping is active. Supports both logical/character
  # columns (scale_alpha_manual) and numeric columns (scale_alpha_continuous).
  # na.translate = FALSE suppresses NA entries in the alpha legend.
  if (!is.null(alpha_var)) {
    
    alpha_col <- df[[alpha_var]]
    
    if (is.numeric(alpha_col)) {
      p <- p +
        ggplot2::scale_alpha_continuous(
          range = range(alpha_values),
          name  = alpha_legend_title
        )
    } else {
      p <- p +
        ggplot2::scale_alpha_manual(
          values       = alpha_values,
          name         = alpha_legend_title,
          labels       = alpha_legend_labels %||% ggplot2::waiver(),
          na.translate = FALSE
        )
    }
    
    p <- p +
      ggplot2::guides(
        alpha = ggplot2::guide_legend(override.aes = list(shape = 16, size = 4))
      )
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Save and Return
  # ═══════════════════════════════════════════════════════════════════════════
  
  if (!is.null(output_dir)) {
    
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    # Sanitise filename — replace non-alphanumeric runs with underscores and
    # strip any leading/trailing underscores introduced by the substitution.
    safe_filename <- gsub("[^[:alnum:]]+", "_", filename)
    safe_filename <- gsub("(^_+|_+$)", "", safe_filename)
    
    output_file <- file.path(output_dir, paste0(safe_filename, ".pdf"))
    
    ggplot2::ggsave(
      filename = output_file,
      plot     = p,
      device   = grDevices::cairo_pdf,
      width    = 8,
      height   = 9,
      bg       = "white"
    )
    
    log_info(sample = "", step = "plot_dot",
             msg = glue::glue("Dot plot saved to: '{output_file}'"))
    
  } else {
    
    log_info(sample = "", step = "plot_dot",
             msg = "Dot plot returned as ggplot object. No output_dir provided — skipping save.")
  }
  
  return(invisible(p))
}

# ═══════════════════════════════════════════════════════════════════════════════
# BAR PLOT
# ═══════════════════════════════════════════════════════════════════════════════
# plot_bar() generates a horizontal bar plot from a pre-wrangled data.frame.
# All data preparation (skeleton building, ranking, label construction) is the
# caller's responsibility — plot_bar() is a pure ggplot assembler that handles
# only display decisions: y-axis text sizing, x-axis limits, NA row suppression,
# color assignment, and optional saving.
#
# INPUT:
#   df                : data.frame, fully prepped and sorted in the desired
#                       display order. Rows with NA or "" in y_var are silently
#                       replaced with " " so ggplot retains them as blank space
#                       rather than dropping them and collapsing the panel.
#                       The y-axis order strictly follows the row order of df —
#                       sort df before passing if a specific order is needed.
#   x_var             : column name (string) mapped to the x-axis (bar length).
#   y_var             : column name (string) mapped to the y-axis (bar labels).
#   x_label           : x-axis title. NULL (default) uses x_var as the label.
#   y_label           : y-axis title. NULL (default) uses y_var as the label.
#   fill_var          : column name (string) for fill color mapping. NULL = no
#                       fill mapping; all bars drawn in ggplot default color.
#   fill_colors       : named character vector mapping fill levels to colors.
#                       NULL (default) = auto-assigned from CUSTOM_PALETTE.
#                       Two-level fill with levels "Upregulated"/"Downregulated"
#                       always uses the standard pipeline colors regardless of
#                       fill_colors, unless fill_colors is explicitly provided.
#   fill_legend_title : legend group title for the fill aesthetic. NULL = ggplot
#                       uses the fill column name.
#   fill_legend_labels: named character vector relabeling fill legend entries.
#                       NULL = ggplot uses the level names as-is.
#   alpha_var         : column name (string) for alpha mapping. NULL = no alpha
#                       mapping; all bars fully opaque. The caller is responsible
#                       for pre-computing any derived alpha column (e.g.
#                       padj <= 0.05) before passing df. plot_bar() maps whatever
#                       values exist in this column to alpha without assumptions.
#   alpha_values      : named numeric vector mapping alpha column levels to
#                       opacity values e.g. c("TRUE" = 1, "FALSE" = 0.35).
#                       Required when alpha is provided; ignored otherwise.
#   alpha_legend_title: legend group title for the alpha aesthetic. NULL = ggplot
#                       uses the alpha column name.
#   alpha_legend_labels: named character vector relabeling alpha legend entries.
#                       NULL = ggplot uses the level names as-is.
#   x_min             : numeric lower bound for the x-axis. NULL (default) =
#                       min(0, min(df[[x_var]], na.rm = TRUE)) so negative bars
#                       are always shown and positive-only data anchors at zero.
#                       Pass x_min = 0 to force the axis to start at zero even
#                       when negative values are present.
#   title             : plot title string. Default "" (no title).
#   filename          : base name for the output PDF. Non-alphanumeric characters
#                       are sanitised to underscores. Default "Bar_Plot".
#   output_dir        : directory to save the PDF. NULL (default) = return the
#                       ggplot object only, no file written.
#
# OUTPUT:
#   Returns invisibly: the ggplot object p.
#   When output_dir is provided, saves a single-page PDF via ggsave +
#   cairo_pdf and logs the output path.

plot_bar <- function(df,
                     x_var,
                     y_var,
                     x_label             = NULL,
                     y_label             = NULL,
                     fill_var            = NULL,
                     fill_colors         = NULL,
                     fill_legend_title   = NULL,
                     fill_legend_labels  = NULL,
                     alpha_var           = NULL,
                     alpha_values        = NULL,
                     alpha_legend_title  = NULL,
                     alpha_legend_labels = NULL,
                     x_min               = NULL,
                     title               = "",
                     filename            = "Bar_Plot",
                     output_dir          = NULL) {
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Validate Inputs
  # ═══════════════════════════════════════════════════════════════════════════
  # Mandatory columns are checked upfront so the function fails immediately
  # with a clear message rather than crashing inside ggplot with a cryptic error.
  
  if (!is.data.frame(df)) {
    log_error(sample = "", step = "plot_bar",
              msg = glue::glue("'df' must be a data.frame, not {class(df)[1]}."))
  }
  
  mandatory_cols <- c(x_var, y_var)
  optional_cols  <- c(fill_var, alpha_var)
  cols_to_check  <- c(mandatory_cols, optional_cols)
  missing_cols   <- cols_to_check[!cols_to_check %in% colnames(df)]
  
  if (length(missing_cols) > 0) {
    log_error(sample = "", step = "plot_bar",
              msg = glue::glue(
                "Column(s) not found in df: {paste(missing_cols, collapse = ', ')}. ",
                "Available: {paste(colnames(df), collapse = ', ')}"
              ))
  }
  
  if (!is.null(alpha_var) && is.null(alpha_values)) {
    log_error(sample = "", step = "plot_bar",
              msg = glue::glue(
                "'alpha_values' must be provided when 'alpha_var' is specified. ",
                "e.g. alpha_values = c('TRUE' = 1, 'FALSE' = 0.35)"
              ))
  }
  
  if (!is.null(output_dir) && !is.character(output_dir)) {
    log_error(sample = "", step = "plot_bar",
              msg = glue::glue(
                "'output_dir' must be a character string or NULL, not {class(output_dir)[1]}."
              ))
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Prepare Data
  # ═══════════════════════════════════════════════════════════════════════════
  
  # ── 2a. Resolve axis labels ───────────────────────────────────────────────
  # %||% returns the left side if non-NULL, otherwise the right side (rlang).
  x_label <- x_label %||% x_var
  y_label <- y_label %||% y_var
  
  # ── 2b. Replace NA / blank y labels with unique placeholders ─────────────
  # ggplot2 drops rows where the y aesthetic is NA. Replacing with unique
  # placeholders (" ", "  ", "   " etc.) forces ggplot to keep each empty
  # row as distinct blank space rather than collapsing them into one.
  empty_mask  <- is.na(df[[y_var]]) | trimws(as.character(df[[y_var]])) == ""
  df[[y_var]] <- as.character(df[[y_var]])
  df[[y_var]][empty_mask] <- sapply(seq_len(sum(empty_mask)), 
                                    function(i) strrep(" ", i))
  
  # ── 2c. Lock y-axis order to df row order ────────────────────────────────
  # factor() with levels = unique() prevents ggplot from re-sorting the y-axis
  # alphabetically. The caller controls display order by sorting df before
  # passing — plot_bar() respects that order exactly.
  df[[y_var]] <- factor(df[[y_var]], levels = unique(df[[y_var]]))
  
  # ── 2d. Compute x_min ────────────────────────────────────────────────────
  # Default anchors at the smaller of 0 or the minimum data value so:
  #   - negative bars always start at a visible left boundary
  #   - positive-only data anchors cleanly at zero (not at the first bar)
  # expansion(mult = c(0.1, 0.05)) adds breathing room on both sides so the
  # leftmost bar never sits flush against the y-axis.
  x_min <- x_min %||% min(0, min(df[[x_var]], na.rm = TRUE))
  
  # ── 2e. Detect x-axis label format ───────────────────────────────────────
  # Compute breaks first, then determine accuracy from the break values.
  # Integer breaks display without decimals (accuracy = 1). Non-integer breaks
  # use accuracy = NULL so scales auto-detects the minimum decimal places needed
  # to make all labels distinct — works correctly for both large ranges (NES)
  # and small fractions (GeneRatio).
  x_range    <- c(x_min, max(df[[x_var]], na.rm = TRUE))
  x_breaks   <- scales::pretty_breaks(n = 4)(x_range)
  x_accuracy <- if (all(x_breaks == floor(x_breaks))) 1 else NULL
  left_expand <- if (x_min < 0) 0.1 else 0
  
  # ── 2f. Compute y-axis text size ─────────────────────────────────────────
  # Long pathway / TF names need smaller text to avoid overlapping. Scale is
  # based on the longest label in the column after NA replacement above.
  max_label_len <- max(nchar(as.character(df[[y_var]])), na.rm = TRUE)
  y_size <- dplyr::case_when(
    max_label_len > 50 ~ 6,
    max_label_len > 35 ~ 7,
    max_label_len > 25 ~ 8,
    TRUE               ~ 10
  )
  
  # ── 2g. Resolve fill colors ───────────────────────────────────────────────
  # Priority: caller-supplied fill_colors > standard Up/Down palette > CUSTOM_PALETTE.
  # The Up/Down shortcut ensures consistent direction colors across all pipeline plots.
  if (!is.null(fill_var) && is.null(fill_colors)) {
    fill_levels <- unique(df[[fill_var]])
    updown_cols <- c("Upregulated" = "#E69F00", "Downregulated" = "#56B4E9")
    
    if (all(fill_levels %in% names(updown_cols))) {
      # Standard direction palette — consistent across all pipeline plots
      fill_colors <- updown_cols
    } else {
      # General case: assign CUSTOM_PALETTE colors in level order
      fill_colors <- setNames(
        CUSTOM_PALETTE[seq_along(fill_levels)],
        fill_levels
      )
    }
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Build Plot
  # ═══════════════════════════════════════════════════════════════════════════
  
  # ── 3a. Base plot ─────────────────────────────────────────────────────────
  # aes() is built dynamically so that fill and alpha are only mapped when the
  # caller explicitly provides them — unmapped aesthetics inherit ggplot defaults
  # and do not produce spurious legend entries.
  base_aes <- ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]])
  
  if (!is.null(fill_var))  base_aes <- utils::modifyList(base_aes, ggplot2::aes(fill  = .data[[fill_var]]))
  if (!is.null(alpha_var)) base_aes <- utils::modifyList(base_aes, ggplot2::aes(alpha = .data[[alpha_var]]))
  
  p <- ggplot2::ggplot(data = df, mapping = base_aes) +
    
    # ── 3b. Bars ─────────────────────────────────────────────────────────────
    ggplot2::geom_col(width = 0.75, na.rm = TRUE) +
    
    # ── 3c. X axis ───────────────────────────────────────────────────────────
    # expansion left = 0.1 gives breathing room between y-axis and first bar.
    # expansion right = 0.05 gives a small margin after the longest bar.
    # Upper limit is left as NA so ggplot auto-scales to the data maximum.
    ggplot2::scale_x_continuous(
      limits = c(x_min, NA),
      expand = ggplot2::expansion(mult = c(left_expand, 0.05)),
      breaks = x_breaks,
      labels = scales::label_number(accuracy = x_accuracy)
    ) +
    
    # ── 3d. Labels and title ─────────────────────────────────────────────────
    ggplot2::labs(x = x_label, y = y_label, title = title) +
    
    # ── 3e. Coordinate system ────────────────────────────────────────────────
    # clip = "off" ensures bar labels or long text that slightly overflows the
    # plot panel boundary is still rendered rather than clipped invisibly.
    ggplot2::coord_cartesian(clip = "off") +
    
    # ── 3f. Y-axis text size ─────────────────────────────────────────────────
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = y_size),
      axis.text.x = if (is.null(x_accuracy)) {
        ggplot2::element_text(angle = 45, hjust = 1)
      } else {
        ggplot2::element_text(angle = 0)
      }
    )
  
  # ── 3g. Fill scale ───────────────────────────────────────────────────────
  # Only added when fill mapping is active. na.translate = FALSE excludes NA
  # levels from the legend — NA bars are still plotted but not labeled.
  if (!is.null(fill_var)) {
    p <- p +
      ggplot2::scale_fill_manual(
        values       = fill_colors,
        name         = fill_legend_title,
        labels       = fill_legend_labels %||% ggplot2::waiver(),
        na.translate = FALSE
      ) +
      ggplot2::guides(
        fill = ggplot2::guide_legend(override.aes = list(shape = 22, size = 6))
      )
  }
  
  # ── 3h. Alpha scale ──────────────────────────────────────────────────────
  # Only added when alpha mapping is active. Supports both logical/character
  # columns (scale_alpha_manual) and numeric columns (scale_alpha_continuous).
  # na.translate = FALSE suppresses NA entries in the alpha legend.
  if (!is.null(alpha_var)) {
    
    alpha_col <- df[[alpha_var]]
    
    if (is.numeric(alpha_col)) {
      p <- p +
        ggplot2::scale_alpha_continuous(
          range = range(alpha_values),
          name  = alpha_legend_title
        )
    } else {
      p <- p +
        ggplot2::scale_alpha_manual(
          values       = alpha_values,
          name         = alpha_legend_title,
          labels       = alpha_legend_labels %||% ggplot2::waiver(),
          na.translate = FALSE
        )
    }
    
    p <- p +
      ggplot2::guides(
        alpha = ggplot2::guide_legend(override.aes = list(shape = 22, size = 6))
      )
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Save and Return
  # ═══════════════════════════════════════════════════════════════════════════
  
  if (!is.null(output_dir)) {
    
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    # Sanitise filename — replace non-alphanumeric runs with underscores and
    # strip any leading/trailing underscores introduced by the substitution.
    safe_filename <- gsub("[^[:alnum:]]+", "_", filename)
    safe_filename <- gsub("(^_+|_+$)", "", safe_filename)
    
    output_file <- file.path(output_dir, paste0(safe_filename, ".pdf"))
    
    ggplot2::ggsave(
      filename = output_file,
      plot     = p,
      device   = grDevices::cairo_pdf,
      width    = 8,
      height   = 9,
      bg       = "white"
    )
    
    log_info(sample = "", step = "plot_bar",
             msg = glue::glue("Bar plot saved to: '{output_file}'"))
    
  } else {
    
    log_info(sample = "", step = "plot_bar",
             msg = "Bar plot returned as ggplot object. No output_dir provided — skipping save.")
  }
  
  return(invisible(p))
}



# ═══════════════════════════════════════════════════════════════════════════════
# compare_deseq2_results()
#
# PURPOSE:
#   Compare two DESeq2 result objects (or DEG Excel sheets) column-by-column
#   to detect and characterise any differences. Designed to cross-validate
#   the "parser" contrast vector method against the "standard" DESeq2 c()
#   method inside extract_deseq2_results(), but also works as a standalone
#   QC tool for comparing results from any two DESeq2 runs.
#
# INPUT MODES (auto-detected via load_smart()):
#   Mode 1 — Two DESeqResults RDS objects  : full metadata + column comparisons
#   Mode 2 — Two DEG Excel files / sheets  : column comparisons only
#
# OUTPUTS:
#   Always  : named list with all comparison metrics + any_differences flag
#   If diff : Comparison_Results.xlsx saved to output_dir (if provided)
#
# CALLED FROM:
#   extract_deseq2_results() — res_a = res_parser, res_b = res_standard
#   Standalone               — pass file paths instead of RDS objects
# ═══════════════════════════════════════════════════════════════════════════════

compare_deseq2_results <- function(res_a,
                                   res_b,
                                   output_dir = NULL,
                                   label_a    = "A",
                                   label_b    = "B") {
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Load Inputs and Detect Mode
  # ═══════════════════════════════════════════════════════════════════════════
  # load_smart() handles RDS, xlsx, and csv transparently.
  # We detect RDS mode by checking if the loaded object is a DESeqResults
  # instance — this determines whether metadata checks are possible.
  
  res_a <- load_smart(res_a)
  res_b <- load_smart(res_b)
  
  is_rds_a <- methods::is(res_a, "DESeqResults")
  is_rds_b <- methods::is(res_b, "DESeqResults")
  is_rds   <- is_rds_a && is_rds_b
  
  # Coerce to plain dataframes for uniform column-level comparisons downstream.
  # For RDS objects, as.data.frame() preserves all result columns (log2FoldChange,
  # lfcSE, stat, pvalue, padj, baseMean) with gene IDs as rownames.
  df_a <- as.data.frame(res_a)
  df_b <- as.data.frame(res_b)
  
  # Ensure gene IDs are an explicit column (rownames in RDS, may already be a
  # column in xlsx). This normalises both modes for downstream joining.
  if (!"gene_id" %in% colnames(df_a)) df_a <- tibble::rownames_to_column(df_a, "gene_id")
  if (!"gene_id" %in% colnames(df_b)) df_b <- tibble::rownames_to_column(df_b, "gene_id")
  
  # Initialise results collector — every check appends to this list
  results         <- list()
  any_differences <- FALSE
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: RDS Metadata Checks (Mode 1 Only)
  # ═══════════════════════════════════════════════════════════════════════════
  # These checks use DESeqResults metadata to verify both results were produced
  # under the same analytical settings. Differences here CONTEXTUALISE all
  # downstream column differences — e.g. different alpha explains padj flips,
  # different cooksCutoff explains NA pattern differences.
  # Skipped entirely in Mode 2 (xlsx) since metadata is not stored in Excel.
  
  if (is_rds) {
    
    meta_a <- S4Vectors::metadata(res_a)
    meta_b <- S4Vectors::metadata(res_b)
    
    # ── 2a. Contrast description ─────────────────────────────────────────────
    # mcols()$description encodes what was actually tested (e.g.
    # "log2 fold change (MLE): condition A vs B"). Mismatch here means the two
    # results are comparing DIFFERENT contrasts — all downstream comparisons
    # would be meaningless. Treat as a hard error.
    
    desc_a <- tryCatch(S4Vectors::mcols(res_a)$description, error = function(e) NULL)
    desc_b <- tryCatch(S4Vectors::mcols(res_b)$description, error = function(e) NULL)
    
    if (!is.null(desc_a) && !is.null(desc_b)) {
      if (!identical(desc_a, desc_b)) {
        log_warn(sample = "", step = "compare_deseq2_results",
                 msg = glue::glue(
                   "Results appear to test DIFFERENT contrasts!\n",
                   "  {label_a}: {paste(desc_a, collapse = ' | ')}\n",
                   "  {label_b}: {paste(desc_b, collapse = ' | ')}\n",
                   "All downstream comparisons may be meaningless."
                 ))
      }
    }
    
    # ── 2b. Analytical settings ──────────────────────────────────────────────
    # Differences in these settings explain (rather than alarm) column-level
    # differences seen later. We log warnings so the user understands WHY
    # differences exist, rather than being surprised by them.
    #
    # alpha          : padj cutoff used for independent filtering — affects which
    #                  genes get a padj value vs NA
    # filterThreshold: actual independent filtering threshold chosen by DESeq2
    #                  (data-driven, derived from alpha) — directly explains NAs
    # cooksCutoff    : genes with Cook's distance above this get NA padj —
    #                  different cutoffs explain different NA patterns
    # modelMatrixType: "standard" vs "expanded" affects how coefficients map to
    #                  groups — should always match for same contrast
    
    settings_checks <- list(
      list(key = "alpha",           label = "Alpha (padj cutoff)"),
      list(key = "filterThreshold", label = "Independent filtering threshold"),
      list(key = "cooksCutoff",     label = "Cook's distance cutoff"),
      list(key = "modelMatrixType", label = "Model matrix type")
    )
    
    metadata_summary <- data.frame(
      Setting = character(),
      label_a = character(),
      label_b = character(),
      Match   = logical(),
      stringsAsFactors = FALSE
    )
    colnames(metadata_summary)[2:3] <- c(label_a, label_b)
    
    for (chk in settings_checks) {
      
      val_a <- meta_a[[chk$key]]
      val_b <- meta_b[[chk$key]]
      
      if (is.null(val_a) || is.null(val_b)) next
      
      match <- isTRUE(all.equal(val_a, val_b))
      
      if (!match) {
        log_info(sample = "", step = "compare_deseq2_results",
                 msg = glue::glue(
                   "{chk$label} differs: {label_a} = {val_a}, {label_b} = {val_b}. ",
                   "This may explain downstream column differences."
                 ))
        any_differences <- TRUE
      }
      
      metadata_summary <- rbind(metadata_summary, data.frame(
        Setting = chk$label,
        A       = as.character(val_a),
        B       = as.character(val_b),
        Match   = match,
        stringsAsFactors = FALSE
      ))
    }
    
    colnames(metadata_summary)[2:3] <- c(label_a, label_b)
    results$metadata_check <- metadata_summary
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Gene List Identity Check
  # ═══════════════════════════════════════════════════════════════════════════
  # Before any column comparison, verify both results contain the same genes.
  # Mismatch implies different input data, different filtering, or an error
  # upstream. We subset to the intersection so downstream comparisons are
  # always gene-matched — but log the discrepancy clearly.
  
  genes_a      <- df_a$gene_id
  genes_b      <- df_b$gene_id
  
  only_in_a    <- setdiff(genes_a, genes_b)
  only_in_b    <- setdiff(genes_b, genes_a)
  shared_genes <- intersect(genes_a, genes_b)
  
  if (length(only_in_a) > 0 || length(only_in_b) > 0) {
    
    log_info(sample = "", step = "compare_deseq2_results",
             msg = glue::glue(
               "Gene lists differ! Only in {label_a}: {length(only_in_a)}, ",
               "Only in {label_b}: {length(only_in_b)}. ",
               "Subsetting to {length(shared_genes)} shared genes for comparisons."
             ))
    
    any_differences <- TRUE
    
    results$gene_list_diff <- list(
      n_only_in_a     = length(only_in_a),
      n_only_in_b     = length(only_in_b),
      n_shared        = length(shared_genes),
      genes_only_in_a = only_in_a,
      genes_only_in_b = only_in_b
    )
    
    # Subset both dataframes to shared genes only for all downstream comparisons
    df_a <- df_a %>% dplyr::filter(gene_id %in% shared_genes) %>%
      dplyr::arrange(gene_id)
    df_b <- df_b %>% dplyr::filter(gene_id %in% shared_genes) %>%
      dplyr::arrange(gene_id)
    
  } else {
    
    # Gene lists match — sort both identically so row-wise comparisons are valid
    df_a <- df_a %>% dplyr::arrange(gene_id)
    df_b <- df_b %>% dplyr::arrange(gene_id)
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Per-Column Comparisons
  # ═══════════════════════════════════════════════════════════════════════════
  # Dynamically compare every numeric column present in BOTH dataframes.
  # Columns present in only one are skipped with a log_info message — this
  # handles the case where xlsx output has fewer columns than RDS (e.g. no
  # lfcSE, no stat) without any hardcoded column lists.
  #
  # For each shared numeric column the comparison pipeline is:
  #   Step 1 — NA pattern differences (always)
  #   Step 2 — Identity check on non-NA pairs (always)
  #   Step 3 — Column-specific comparisons (only if not identical)
  
  # Numeric columns to compare — excludes gene_id and any character columns
  numeric_cols_a <- colnames(df_a)[sapply(df_a, is.numeric)]
  numeric_cols_b <- colnames(df_b)[sapply(df_b, is.numeric)]
  shared_cols    <- intersect(numeric_cols_a, numeric_cols_b)
  
  # Log any columns present in one but not the other
  only_cols_a <- setdiff(numeric_cols_a, numeric_cols_b)
  only_cols_b <- setdiff(numeric_cols_b, numeric_cols_a)
  
  if (length(only_cols_a) > 0) {
    log_info(sample = "", step = "compare_deseq2_results",
             msg = glue::glue("Columns only in {label_a} (skipped): {paste(only_cols_a, collapse = ', ')}"))
  }
  if (length(only_cols_b) > 0) {
    log_info(sample = "", step = "compare_deseq2_results",
             msg = glue::glue("Columns only in {label_b} (skipped): {paste(only_cols_b, collapse = ', ')}"))
  }
  
  # Collector for per-column summary (one row per column in final Excel sheet)
  col_summary <- data.frame(
    Column            = character(),
    N_NA_only_A       = integer(),
    N_NA_only_B       = integer(),
    N_NA_both         = integer(),
    Identical         = logical(),
    Max_Delta         = numeric(),
    Pearson_R         = numeric(),
    Spearman_R        = numeric(),
    N_Direction_Flips = integer(),
    N_Sig_Flips       = integer(),
    stringsAsFactors  = FALSE
  )
  
  # Per-column detail sheets — populated below, written to Excel if differences found
  na_pattern_rows     <- list()
  direction_flip_rows <- list()
  sig_flip_rows       <- list()
  
  for (col in shared_cols) {
    
    vec_a <- df_a[[col]]
    vec_b <- df_b[[col]]
    
    # ── 4a. NA pattern differences ───────────────────────────────────────────
    # Three NA categories:
    #   na_only_a : NA in A but not B — A filtered this gene out, B didn't
    #   na_only_b : NA in B but not A — B filtered this gene out, A didn't
    #   na_both   : NA in both        — both filtered, expected and fine
    # na_only_* cases are the interesting ones — they indicate filtering
    # differences (independent filtering threshold, Cook's cutoff) between runs.
    
    na_only_a <- is.na(vec_a) & !is.na(vec_b)
    na_only_b <- !is.na(vec_a) & is.na(vec_b)
    na_both   <- is.na(vec_a) & is.na(vec_b)
    
    n_na_only_a <- sum(na_only_a)
    n_na_only_b <- sum(na_only_b)
    n_na_both   <- sum(na_both)
    
    if (n_na_only_a > 0 || n_na_only_b > 0) {
      
      any_differences <- TRUE
      
      na_pattern_rows[[col]] <- data.frame(
        Column  = col,
        gene_id = df_a$gene_id[na_only_a | na_only_b],
        NA_in   = dplyr::case_when(
          na_only_a[na_only_a | na_only_b] ~ label_a,
          TRUE                             ~ label_b
        ),
        stringsAsFactors = FALSE
      )
    }
    
    # Work only with genes where BOTH have non-NA values for all further checks
    both_nonNA  <- !is.na(vec_a) & !is.na(vec_b)
    vec_a_clean <- vec_a[both_nonNA]
    vec_b_clean <- vec_b[both_nonNA]
    genes_clean <- df_a$gene_id[both_nonNA]
    
    if (length(vec_a_clean) < 2) {
      # Not enough non-NA pairs to compare — skip this column
      log_info(sample = "", step = "compare_deseq2_results",
               msg = glue::glue("Column '{col}': fewer than 2 non-NA pairs — skipping comparisons."))
      next
    }
    
    # ── 4b. Identity check ───────────────────────────────────────────────────
    # Check if all non-NA values are numerically identical within floating point
    # tolerance. If identical, skip all further comparisons for this column —
    # no differences to report, log and move on.
    
    tolerance    <- 1e-8
    max_delta    <- max(abs(vec_a_clean - vec_b_clean), na.rm = TRUE)
    is_identical <- max_delta < tolerance
    
    if (is_identical) {
      log_info(sample = "", step = "compare_deseq2_results",
               msg = glue::glue("Column '{col}': identical between {label_a} and {label_b}."))
      
      col_summary <- rbind(col_summary, data.frame(
        Column            = col,
        N_NA_only_A       = n_na_only_a,
        N_NA_only_B       = n_na_only_b,
        N_NA_both         = n_na_both,
        Identical         = TRUE,
        Max_Delta         = 0,
        Pearson_R         = NA_real_,
        Spearman_R        = NA_real_,
        N_Direction_Flips = NA_integer_,
        N_Sig_Flips       = NA_integer_,
        stringsAsFactors  = FALSE
      ))
      next
    }
    
    # Values differ — flag and run column-specific comparisons
    any_differences <- TRUE
    log_info(sample = "", step = "compare_deseq2_results",
             msg = glue::glue("Column '{col}': differences detected. Max delta = ",
                              "{formatC(max_delta, format = 'e', digits = 2)}"))
    
    # ── 4c. Correlation ──────────────────────────────────────────────────────
    # Both Pearson (linear agreement) and Spearman (rank agreement) are computed.
    # Pearson catches magnitude differences; Spearman catches ranking differences
    # (e.g. gene ordering for volcano plots or ranked pathway analyses).
    # Both should be near 1.0 for well-behaved contrasts.
    
    pearson_r  <- tryCatch(stats::cor(vec_a_clean, vec_b_clean, method = "pearson"),
                           error = function(e) NA_real_)
    spearman_r <- tryCatch(stats::cor(vec_a_clean, vec_b_clean, method = "spearman"),
                           error = function(e) NA_real_)
    
    # ── 4d. Column-specific comparisons ─────────────────────────────────────
    
    n_direction_flips <- NA_integer_
    n_sig_flips       <- NA_integer_
    
    # ~~ Direction flips (log2FoldChange and stat columns) ~~
    # A direction flip means the sign of the effect CHANGES between methods —
    # e.g. a gene appears upregulated in A but downregulated in B.
    # These are the most biologically alarming differences since they would
    # lead to opposite biological interpretations.
    # We check sign() != 0 on both sides to exclude genes sitting exactly at 0.
    
    if (col %in% c("log2FoldChange", "stat")) {
      
      sign_a <- sign(vec_a_clean)
      sign_b <- sign(vec_b_clean)
      
      flip_mask         <- (sign_a != sign_b) & (sign_a != 0) & (sign_b != 0)
      n_direction_flips <- sum(flip_mask)
      
      if (n_direction_flips > 0) {
        direction_flip_rows[[col]] <- data.frame(
          Column  = col,
          gene_id = genes_clean[flip_mask],
          A       = vec_a_clean[flip_mask],
          B       = vec_b_clean[flip_mask],
          Delta   = vec_a_clean[flip_mask] - vec_b_clean[flip_mask],
          stringsAsFactors = FALSE
        )
        colnames(direction_flip_rows[[col]])[3:4] <- c(label_a, label_b)
      }
    }
    
    # ~~ Significance flips (padj and pvalue columns) ~~
    # A significance flip means a gene crosses the 0.05 threshold in opposite
    # directions between methods — significant in one, not in the other.
    # We further classify each flip as:
    #   Borderline : the non-significant value is between 0.04-0.15 — the gene
    #                is near the threshold and the flip may reflect minor
    #                numerical differences in filtering
    #   Genuine    : the non-significant value is > 0.15 — a real disagreement
    #                that would meaningfully change biological conclusions
    
    if (col %in% c("padj", "pvalue")) {
      
      sig_threshold <- 0.05
      borderline_lo <- 0.04
      borderline_hi <- 0.15
      
      sig_a <- vec_a_clean < sig_threshold
      sig_b <- vec_b_clean < sig_threshold
      
      flip_mask   <- sig_a != sig_b
      n_sig_flips <- sum(flip_mask)
      
      if (n_sig_flips > 0) {
        
        flip_a <- vec_a_clean[flip_mask]
        flip_b <- vec_b_clean[flip_mask]
        nonsig <- ifelse(sig_a[flip_mask], flip_b, flip_a)   # the non-significant value
        
        sig_flip_rows[[col]] <- data.frame(
          Column     = col,
          gene_id    = genes_clean[flip_mask],
          A          = flip_a,
          B          = flip_b,
          Sig_in     = dplyr::case_when(
            sig_a[flip_mask] ~ label_a,
            TRUE             ~ label_b
          ),
          Borderline = nonsig >= borderline_lo & nonsig <= borderline_hi,
          stringsAsFactors = FALSE
        )
        colnames(sig_flip_rows[[col]])[3:4] <- c(label_a, label_b)
      }
    }
    
    # ~~ baseMean identity warning ~~
    # baseMean is the mean normalised count across ALL samples — it is derived
    # from the count matrix and should be IDENTICAL between any two results
    # run on the same dataset. Differences here imply different input data,
    # different size factors, or a sample subsetting error upstream.
    
    if (col == "baseMean" && !is_identical) {
      log_error(sample = "", step = "compare_deseq2_results",
                msg = glue::glue(
                  "baseMean differs between {label_a} and {label_b}! ",
                  "This implies different input counts or size factors — ",
                  "results may not be comparable."
                ))
    }
    
    # Append summary row for this column
    col_summary <- rbind(col_summary, data.frame(
      Column            = col,
      N_NA_only_A       = n_na_only_a,
      N_NA_only_B       = n_na_only_b,
      N_NA_both         = n_na_both,
      Identical         = FALSE,
      Max_Delta         = max_delta,
      Pearson_R         = pearson_r,
      Spearman_R        = spearman_r,
      N_Direction_Flips = n_direction_flips,
      N_Sig_Flips       = n_sig_flips,
      stringsAsFactors  = FALSE
    ))
  }
  
  # Fix column names in summary to use provided labels
  colnames(col_summary)[colnames(col_summary) == "N_NA_only_A"] <- paste0("N_NA_only_", label_a)
  colnames(col_summary)[colnames(col_summary) == "N_NA_only_B"] <- paste0("N_NA_only_", label_b)
  
  results$col_summary     <- col_summary
  results$na_patterns     <- dplyr::bind_rows(na_pattern_rows)
  results$direction_flips <- dplyr::bind_rows(direction_flip_rows)
  results$sig_flips       <- dplyr::bind_rows(sig_flip_rows)
  results$any_differences <- any_differences
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5: Save Results (Only if Differences Found)
  # ═══════════════════════════════════════════════════════════════════════════
  # Comparison Excel is only saved when differences are detected — if results
  # are identical there is nothing useful to save and no file is written.
  # This mirrors the logic in extract_deseq2_results() for res_standard.rds.
  
  if (any_differences && !is.null(output_dir)) {
    
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    sheets <- list()
    
    if (!is.null(results$metadata_check) && nrow(results$metadata_check) > 0) {
      sheets[["Metadata_Check"]] <- results$metadata_check
    }
    
    if (!is.null(results$gene_list_diff)) {
      sheets[["Gene_List_Diff"]] <- data.frame(
        Category = c(paste0("Only in ", label_a),
                     paste0("Only in ", label_b),
                     "Shared"),
        N_Genes  = c(results$gene_list_diff$n_only_in_a,
                     results$gene_list_diff$n_only_in_b,
                     results$gene_list_diff$n_shared)
      )
    }
    
    if (nrow(col_summary) > 0)             sheets[["Column_Summary"]]    <- col_summary
    if (nrow(results$na_patterns) > 0)     sheets[["NA_Patterns"]]       <- results$na_patterns
    if (nrow(results$direction_flips) > 0) sheets[["Direction_Flips"]]   <- results$direction_flips
    if (nrow(results$sig_flips) > 0)       sheets[["Significance_Flips"]] <- results$sig_flips
    
    openxlsx::write.xlsx(
      x         = sheets,
      file      = file.path(output_dir, "Comparison_Results.xlsx"),
      overwrite = TRUE
    )
    
    log_info(sample = "", step = "compare_deseq2_results",
             msg = glue::glue("Differences found. Comparison saved to: ",
                              "{file.path(output_dir, 'Comparison_Results.xlsx')}"))
    
  } else if (!any_differences) {
    
    log_info(sample = "", step = "compare_deseq2_results",
             msg = glue::glue("✅ {label_a} and {label_b} results are identical. No comparison file saved."))
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 6: Return
  # ═══════════════════════════════════════════════════════════════════════════
  # Returns the full comparison list so the caller (extract_deseq2_results or
  # interactive user) can use comparison$any_differences to decide what else
  # to save, without needing to re-read the Excel file.
  
  return(invisible(results))
}

# ═══════════════════════════════════════════════════════════════════════════════
# SURVIVAL ANALYSIS FUNCTION
# ═══════════════════════════════════════════════════════════════════════════════
# plot_survival() performs Kaplan-Meier survival analysis with Cox proportional
# hazards modeling. It handles three modes automatically based on inputs:
#
#   (a) Expression-based  — stratify_var matches gene(s) in expr_data[[row_id]].
#                           Samples are split into HIGH/LOW groups using the
#                           specified cutoff_method. Multiple genes are plotted
#                           individually or combined into one signature score.
#   (b) Metadata-based    — stratify_var matches a column in metadata.
#                           Groups are taken directly from that column's values
#                           (e.g. Sex → Male vs Female, Stage → I/II/III/IV).
#
# INPUT:
#   metadata      : data.frame with at least Sample_ID, time_col, status_col.
#                   Any faceting columns must also be present here.
#   stratify_var  : character. One or more gene names (found in expr_data[[row_id]])
#                   OR a single metadata column name. Cannot be NULL.
#                   If genes are supplied, expr_data must also be provided.
#   time_col      : name of the survival time column in metadata (default "Time").
#   status_col    : name of the event status column in metadata (default "Status").
#                   Status must be coded 0 = censored, 1 = event.
#   expr_data     : data.frame with one gene ID column (row_id) + numeric sample
#                   columns. Same format as produced by process_counts(). NULL for
#                   metadata-based survival.
#                   IMPORTANT: expr_data must be pre-transformed (VST, log-normalized
#                   etc.) before passing to plot_survival(). No transformation is
#                   applied internally — only median centering for signature scores.
#                   Recommended inputs (in order of preference):
#                     (1) VST blind counts  : vst(dds, blind = TRUE)
#                     (2) Log-normalized    : log2(normalized_counts + 1)
#                     (3) rlog counts       : rlog(dds) — preferred for < 30 samples
#                   Avoid raw counts — cutoff methods require a continuous,
#                   variance-stabilized scale to produce meaningful HIGH/LOW splits.
#   row_id        : name of the gene ID column in expr_data (default "SYMBOL").
#                   Same convention as plot_heatmap().
#   facet_var     : optional metadata column to split plots by. Each unique level
#                   gets its own KM curve. All facets for one gene appear on one
#                   PDF page arranged in a sqrt-based grid layout.
#   global_cutoff : logical. TRUE (default) = one HIGH/LOW cutoff computed on the
#                   full dataset and applied to all facets. FALSE = cutoff computed
#                   separately within each facet group. Only relevant for
#                   expression-based survival when facet_var is provided.
#                   Use FALSE when expression distributions differ markedly between
#                   facet groups (e.g. sex-specific expression ranges).
#   cutoff_method : how to split samples into expression bins. One of:
#                   "optimal"  — maximises log-rank statistic via surv_cutpoint()
#                   "median"   — splits at 50th percentile → 2 groups
#                   "tertile"  — splits at 33rd/67th percentile → HIGH and LOW only
#                   "quartile" — splits at 25th/75th percentile → HIGH and LOW only
#                   Middle bins (MID, MED_LOW, MED_HIGH) are always dropped unless
#                   show_all_bins = TRUE. Only relevant for expression-based survival.
#   show_all_bins        : logical. FALSE (default) = plot only HIGH vs LOW.
#                          TRUE = include middle bins (MID / MED_LOW / MED_HIGH).
#                          Only relevant for tertile and quartile cutoff methods.
#   calc_signature_score : logical. Only relevant when multiple genes are supplied.
#                          FALSE (default) = plot a separate KM curve per gene.
#                          TRUE  = combine all genes into one signature score
#                          (advanced z-score method; Levine et al. 2006,
#                          doi:10.1186/gb-2006-7-10-r93) and plot a single KM curve.
#   conf_interval : logical. FALSE (default) = no confidence ribbon on KM curve.
#                   TRUE = show 95% confidence interval ribbon.
#   plot_curve    : logical. TRUE (default) = generate KM curve PDF.
#                   FALSE = compute and save stats only (Excel), skip PDF.
#                   Useful for batch runs across many genes where only stats matter.
#   plot_risk_table: logical. TRUE (default) = show number-at-risk table below curve.
#   time_units    : character. Units of time_col used for x-axis label.
#                   One of "days", "months" (default), "years".
#   title         : optional plot title string. NULL (default) = auto-generated
#                   from stratify_var and facet group.
#   filename      : base name for output files (default "Survival_Plot").
#                   Non-alphanumeric characters are replaced with underscores.
#                   Both PDF and Excel share this base name.
#   output_dir    : directory to save PDF and Excel. NULL (default) = skip save,
#                   return results list only.
#
# OUTPUT:
#   Always returns invisible(list(plots, cox_stats, emmeans_stats, surv_data)).
#   When output_dir is provided, also saves:
#     <filename>.pdf  — multi-page KM curve PDF (one page per gene or facet grid)
#     <filename>.xlsx — three sheets: cox_stats, emmeans_stats, surv_data
#
#   cox_stats     : Cox model HR, 95% CI, and Wald p-value per gene x facet.
#                   Also includes all 7 non-parametric log-rank p-value variants.
#   emmeans_stats : All pairwise HRs and p-values in both directions (A/B and B/A)
#                   per gene x facet. Includes all 7 non-parametric p-value variants.
#   surv_data     : Merged survival data.frame with model (HIGH/LOW) columns appended.

plot_survival <- function(metadata,
                          stratify_var,
                          time_col             = "Time",
                          status_col           = "Status",
                          expr_data            = NULL,
                          row_id               = "SYMBOL",
                          facet_var            = NULL,
                          global_cutoff        = TRUE,
                          cutoff_method        = "optimal",
                          show_all_bins        = FALSE,
                          calc_signature_score = FALSE,
                          conf_interval        = FALSE,
                          plot_curve           = TRUE,
                          plot_risk_table      = TRUE,
                          time_units           = "months",
                          title                = NULL,
                          filename             = "Survival_Plot",
                          output_dir           = NULL) {
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Validate Inputs
  # ═══════════════════════════════════════════════════════════════════════════
  # All structural checks fire before any data is touched. Failures produce
  # clear, actionable messages via log_error() rather than cryptic crashes
  # deep inside survminer or dplyr.
  
  # ── 1a. metadata must be a data.frame ────────────────────────────────────
  if (!is.data.frame(metadata)) {
    log_error(sample = "", step = "plot_survival",
              msg = glue::glue(
                "'metadata' must be a data.frame, not {class(metadata)[1]}."
              ))
  }
  
  # ── 1b. stratify_var must be provided and non-empty ──────────────────────
  if (missing(stratify_var) || is.null(stratify_var) || length(stratify_var) == 0) {
    log_error(sample = "", step = "plot_survival",
              msg = glue::glue(
                "'stratify_var' must be provided. Supply gene name(s) from ",
                "expr_data[['{row_id}']] or a column name from metadata."
              ))
  }
  
  # ── 1c. time_col and status_col must exist in metadata ───────────────────
  missing_surv_cols <- setdiff(c(time_col, status_col), colnames(metadata))
  if (length(missing_surv_cols) > 0) {
    log_error(sample = "", step = "plot_survival",
              msg = glue::glue(
                "Survival column(s) not found in metadata: ",
                "{paste(missing_surv_cols, collapse = ', ')}. ",
                "Available columns: {paste(colnames(metadata), collapse = ', ')}"
              ))
  }
  
  # ── 1d. facet_var must exist in metadata if provided ─────────────────────
  if (!is.null(facet_var) && !facet_var %in% colnames(metadata)) {
    log_error(sample = "", step = "plot_survival",
              msg = glue::glue(
                "'facet_var' column '{facet_var}' not found in metadata. ",
                "Available columns: {paste(colnames(metadata), collapse = ', ')}"
              ))
  }
  
  # ── 1e. output_dir must be a string if provided ──────────────────────────
  if (!is.null(output_dir) && !is.character(output_dir)) {
    log_error(sample = "", step = "plot_survival",
              msg = glue::glue(
                "'output_dir' must be a character string or NULL, ",
                "not {class(output_dir)[1]}."
              ))
  }
  
  # ── 1f. cutoff_method must be one of the supported methods ───────────────
  valid_methods <- c("optimal", "median", "tertile", "quartile")
  if (!cutoff_method %in% valid_methods) {
    log_error(sample = "", step = "plot_survival",
              msg = glue::glue(
                "'cutoff_method' must be one of: {paste(valid_methods, collapse = ', ')}. ",
                "Got: '{cutoff_method}'."
              ))
  }
  
  # ── 1g. time_units must be one of the supported units ────────────────────
  valid_units <- c("days", "months", "years")
  if (!time_units %in% valid_units) {
    log_error(sample = "", step = "plot_survival",
              msg = glue::glue(
                "'time_units' must be one of: {paste(valid_units, collapse = ', ')}. ",
                "Got: '{time_units}'."
              ))
  }
  
  # ── 1h. Determine survival type from stratify_var ─────────────────────────
  # Check metadata columns first, then expr_data[[row_id]]. Ambiguity (matches
  # both) is an error — gene names should not collide with metadata column names.
  in_metadata <- all(stratify_var %in% colnames(metadata))
  in_expr     <- !is.null(expr_data) &&
    is.data.frame(expr_data) &&
    row_id %in% colnames(expr_data) &&
    any(stratify_var %in% expr_data[[row_id]])
  
  if (in_metadata && in_expr) {
    log_error(sample = "", step = "plot_survival",
              msg = glue::glue(
                "'stratify_var' matches both a metadata column and gene(s) in ",
                "expr_data[['{row_id}']] — ambiguous. Rename the conflicting ",
                "metadata column or gene to resolve."
              ))
  }
  
  if (!in_metadata && !in_expr) {
    log_error(sample = "", step = "plot_survival",
              msg = glue::glue(
                "'stratify_var' not found in metadata columns or ",
                "expr_data[['{row_id}']]. Check spelling and case. ",
                "Available metadata columns: {paste(colnames(metadata), collapse = ', ')}"
              ))
  }
  
  # ── 1i. expr_data must be a data.frame with row_id column ────────────────
  if (!is.null(expr_data)) {
    if (!is.data.frame(expr_data)) {
      log_error(sample = "", step = "plot_survival",
                msg = glue::glue(
                  "'expr_data' must be a data.frame with a '{row_id}' column ",
                  "and numeric sample columns. Got: {class(expr_data)[1]}. ",
                  "Use process_counts() to produce the correct format."
                ))
    }
    if (!row_id %in% colnames(expr_data)) {
      log_error(sample = "", step = "plot_survival",
                msg = glue::glue(
                  "Gene ID column '{row_id}' not found in expr_data. ",
                  "Available columns: {paste(colnames(expr_data), collapse = ', ')}"
                ))
    }
  }
  
  # ── 1j. calc_signature_score requires multiple genes ─────────────────────
  if (isTRUE(calc_signature_score) && length(stratify_var) < 2) {
    log_error(sample = "", step = "plot_survival",
              msg = glue::glue(
                "'calc_signature_score = TRUE' requires at least 2 genes in ",
                "'stratify_var'. Only {length(stratify_var)} gene supplied."
              ))
  }
  
  # ── 1k. Warn about genes missing from expr_data ───────────────────────────
  # Missing genes are skipped with a warning rather than halting — a typo in
  # one gene name should not abort a 100-gene batch run.
  if (in_expr) {
    missing_genes <- setdiff(stratify_var, expr_data[[row_id]])
    if (length(missing_genes) > 0) {
      log_warn(sample = "", step = "plot_survival",
               msg = glue::glue(
                 "{length(missing_genes)} gene(s) not found in ",
                 "expr_data[['{row_id}']] and will be skipped: ",
                 "{paste(missing_genes, collapse = ', ')}"
               ))
      stratify_var <- intersect(stratify_var, expr_data[[row_id]])
    }
    if (length(stratify_var) == 0) {
      log_error(sample = "", step = "plot_survival",
                msg = "No valid genes remain in 'stratify_var' after removing missing genes.")
    }
  }
  
  # ── 1l. Assign survival type ──────────────────────────────────────────────
  surv_type <- if (in_metadata) {
    "metadata"
  } else if (isTRUE(calc_signature_score)) {
    "signature"
  } else {
    "expression"
  }
  
  log_info(sample = "", step = "plot_survival",
           msg = glue::glue(
             "Survival type: '{surv_type}'. ",
             "Stratifying by: {paste(stratify_var, collapse = ', ')}"
           ))
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Prepare Data
  # ═══════════════════════════════════════════════════════════════════════════
  
  # ── 2a. Sanitise Sample_ID and clean metadata ─────────────────────────────
  # make.names() ensures Sample_IDs are valid R names — survminer internals
  # sometimes use sample IDs as list names and fail on special characters.
  # Time coerced to numeric in case it was read as character from Excel.
  # Samples with time <= 0 are uninformative for survival and dropped.
  n_before  <- nrow(metadata)
  metadata  <- metadata %>%
    dplyr::mutate(Sample_ID  = make.names(Sample_ID),
                  !!time_col := as.numeric(.data[[time_col]])) %>%
    dplyr::filter(.data[[time_col]] > 0, !is.na(.data[[time_col]])) %>%
    dplyr::distinct(Sample_ID, .keep_all = TRUE)
  n_dropped <- n_before - nrow(metadata)
  
  if (n_dropped > 0) {
    log_warn(sample = "", step = "plot_survival",
             msg = glue::glue(
               "Dropped {n_dropped} metadata row(s): time <= 0, NA time, ",
               "or duplicate Sample_ID."
             ))
  }
  
  # ── 2b. Build expression or grouping data.frame ───────────────────────────
  # expr_df always has Sample_ID + one column per stratify_var (or sig_score).
  if (surv_type == "signature") {
    
    n_sig_genes <- length(stratify_var)   # save before stratify_var is reassigned
    
    # When duplicates exist keep the highest-expressing copy — same logic as
    # plot_heatmap(). Multiple rows for the same gene would inflate the average
    # biasing the signature score upward regardless of scoring method.
    if (any(duplicated(expr_data[[row_id]]))) {
      log_warn(sample = "", step = "plot_survival",
               msg = "Duplicate gene symbols in expr_data. Retaining highest-expressing copy per gene.")
      expr_data <- expr_data %>%
        dplyr::mutate(total_expr = rowSums(dplyr::across(-dplyr::all_of(row_id), abs))) %>%
        dplyr::group_by(.data[[row_id]]) %>%
        dplyr::slice_max(order_by = total_expr, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::select(-total_expr)
    }
    
    expr_matrix <- expr_data %>%
      tibble::column_to_rownames(row_id) %>%
      as.matrix()
    
    # ── Signature score: mean z-score ─────────────────────────────────────
    # Step 1: z-score each gene independently across all samples.
    #         scale() operates on columns so transpose → scale → transpose back.
    #         Each gene now has mean=0, SD=1 across samples — removes
    #         gene-level baseline differences (high vs low expressed genes
    #         contribute equally to the final score).
    # Step 2: average per-gene z-scores across the gene set per sample.
    #         Result is a single score per sample reflecting combined
    #         relative expression of the gene set.
    #
    # Why not Levine et al. 2006:
    #   Levine anchors scores to the full transcriptome average. When gene
    #   set is lowly expressed relative to transcriptome (e.g. neuronal genes
    #   in bladder cancer), ALL scores become negative and the optimal cutoff
    #   finds a biologically meaningless boundary — HIGH patients actually
    #   have LOWER expression than LOW patients. Validated across 6 bladder
    #   cancer subtypes: zscore correctly classified 6/6, Levine 1/6.
    #   zscore anchors to each gene's own distribution across samples and
    #   is robust regardless of absolute expression level.
    z_matrix   <- t(scale(t(expr_matrix)))
    sig_scores <- colMeans(z_matrix[stratify_var, , drop = FALSE], na.rm = TRUE)
    
    expr_df      <- data.frame(Sample_ID  = make.names(names(sig_scores)),
      sig_score  = as.numeric(sig_scores),
      check.names = FALSE)
    
    stratify_var <- "sig_score"   # redirect all downstream logic to this column
    
    log_info(sample = "", step = "plot_survival",
             msg = glue::glue("Signature score computed from {n_sig_genes} genes."))
    
  } else if (surv_type == "expression") {
    
    # Filter to genes of interest, transpose to samples x genes, add Sample_ID.
    expr_df <- expr_data %>%
      dplyr::filter(.data[[row_id]] %in% stratify_var) %>%
      tibble::column_to_rownames(row_id) %>%
      t() %>%
      as.data.frame(check.names = FALSE) %>%
      tibble::rownames_to_column("Sample_ID") %>%
      dplyr::mutate(Sample_ID = make.names(Sample_ID))
    
  } else {
    # Metadata-based: only Sample_ID pulled here — stratify_var brought in
    # via inner_join in 2c to avoid .x / .y column name collisions.
    expr_df <- metadata %>%
      dplyr::select(Sample_ID)
  }
  
  # ── 2c. Merge expression data with metadata ───────────────────────────────
  # inner_join retains only samples present in BOTH objects — samples without
  # survival data or expression values are excluded. Join count is logged so
  # unexpectedly large sample loss surfaces immediately.
  keep_cols <- unique(c(time_col, status_col, stratify_var, facet_var))
  
  surv_df <- expr_df %>%
    dplyr::inner_join(metadata, by = "Sample_ID") %>%
    dplyr::select(Sample_ID, dplyr::all_of(keep_cols))
  
  if (nrow(surv_df) == 0) {
    log_error(sample = "", step = "plot_survival",
              msg = glue::glue(
                "No overlapping Sample_IDs between expr_data and metadata after joining. ",
                "Check that Sample_IDs match exactly (case-sensitive)."
              ))
  }
  
  log_info(sample = "", step = "plot_survival",
           msg = glue::glue(
             "{nrow(surv_df)} samples retained after joining metadata and expression data."
           ))
  
  # ── 2d. Define facet groups ───────────────────────────────────────────────
  # NA_character_ sentinel when no faceting requested — the main loop treats
  # NA facet_group as "use full dataset", avoiding a special-case branch.
  facet_groups <- if (!is.null(facet_var)) {
    sort(unique(as.character(surv_df[[facet_var]])))
  } else {
    NA_character_
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Calculate Cutoffs (Expression-Based Only)
  # ═══════════════════════════════════════════════════════════════════════════
  # Cutoffs convert continuous expression into categorical HIGH/LOW bins.
  # Skipped for metadata-based survival — groups come directly from the column.
  #
  # global_cutoff = TRUE  → one cutoff on full dataset applied to all facets.
  #                         HIGH/LOW labels are comparable across facets.
  # global_cutoff = FALSE → separate cutoff per facet. Use when expression
  #                         distributions differ markedly between groups
  #                         (e.g. sex-specific expression ranges).
  #
  # For signature scores, global_cutoff = TRUE is strongly recommended —
  # the score is computed on the full dataset so the cutoff should also
  # reflect the full distribution, not per-facet subsets.
  
  if (surv_type %in% c("expression", "signature")) {
    
    model_cols_df <- data.frame(Sample_ID = surv_df$Sample_ID)
    
    for (sv in stratify_var) {
      if (isTRUE(global_cutoff)) {
        # ── Single cutoff on full dataset ──────────────────────────────────
        binned <- .calc_survival_cutoffs(
          df            = surv_df,
          stratify_var  = sv,
          time_col      = time_col,
          status_col    = status_col,
          cutoff_method = cutoff_method,
          show_all_bins = show_all_bins
        )
        
        # Drop existing cols before join to prevent .x/.y
        model_cols_df <- model_cols_df %>%
          dplyr::select(-dplyr::any_of(setdiff(names(binned), "Sample_ID"))) %>%
          dplyr::left_join(binned, by = "Sample_ID")
        
      } else {
        # ── Separate cutoff per facet group ────────────────────────────────
        facet_binned <- purrr::map_dfr(facet_groups, function(fg) {
          df_sub <- if (is.na(fg)) surv_df else
            dplyr::filter(surv_df, .data[[facet_var]] == fg)
          
          .calc_survival_cutoffs(
            df            = df_sub,
            stratify_var  = sv,
            time_col      = time_col,
            status_col    = status_col,
            cutoff_method = cutoff_method,
            show_all_bins = show_all_bins
          )
        })
        
        # Deduplicate: samples should not span multiple facets with a properly
        # defined facet variable, but guards against edge cases.
        facet_binned <- dplyr::distinct(facet_binned, Sample_ID, .keep_all = TRUE)
        
        # Drop existing cols before join to prevent .x/.y
        model_cols_df <- model_cols_df %>%
          dplyr::select(-dplyr::any_of(setdiff(names(facet_binned), "Sample_ID"))) %>%
          dplyr::left_join(facet_binned, by = "Sample_ID")
      }
    }
    
    # Final merge into surv_df
    surv_df <- surv_df %>%
      dplyr::select(-dplyr::any_of(setdiff(names(model_cols_df), "Sample_ID"))) %>%
      dplyr::left_join(model_cols_df, by = "Sample_ID")
    
  } else {
    # Metadata-based: coerce to character for consistent factor handling
    # downstream regardless of original column type (factor, integer, etc.).
    model_col <- paste0("model_", stratify_var)
    surv_df   <- surv_df %>%
      dplyr::mutate(!!model_col := as.character(.data[[stratify_var]]))
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3b: Validate HIGH/LOW Classification Against Raw Expression
  # ═══════════════════════════════════════════════════════════════════════════
  # For expression and signature types, verify that patients labeled HIGH
  # actually have higher expression than patients labeled LOW.
  # This catches scoring method failures (e.g. Levine on lowly-expressed genes)
  # where HIGH/LOW labels are biologically inverted.
  #
  # Evidence sheet saved to Excel alongside stats — provides an auditable
  # record that classifications were correct for every gene x facet x subgroup.
  # Flag column makes failures immediately visible without manual inspection.
  #
  # For signature type: validation uses sig_score (the combined score) since
  # individual gene expression is not retained in surv_df. sig_score is
  # monotonically related to combined expression so HIGH sig_score = HIGH
  # combined expression — sufficient for validation purposes.
  # For expression type: validates directly against each gene's raw values.
  # For metadata type: no validation needed — groups come from metadata directly.
  
  classification_evidence <- tibble::tibble()
  
  if (surv_type %in% c("expression", "signature")) {
    
    # Determine which columns to validate against:
    # signature → sig_score (only combined score available in surv_df)
    # expression → each gene in stratify_var directly
    validate_vars <- if (surv_type == "signature") "sig_score" else stratify_var
    
    for (sv in stratify_var) {
      
      model_col    <- paste0("model_", sv)
      validate_col <- if (surv_type == "signature") "sig_score" else sv
      
      for (fg in facet_groups) {
        
        facet_df_ev <- if (is.na(fg)) surv_df else
          dplyr::filter(surv_df, .data[[facet_var]] == fg)
        
        # Skip facets with fewer than 2 groups — already handled in Section 4
        if (dplyr::n_distinct(facet_df_ev[[model_col]], na.rm = TRUE) < 2) next
        
        # Compute mean, median, Q25, Q75 per group
        ev <- facet_df_ev %>%
          dplyr::filter(!is.na(.data[[model_col]]),
                        !is.na(.data[[validate_col]])) %>%
          dplyr::group_by(.data[[model_col]]) %>%
          dplyr::summarise(
            mean_expr   = mean(.data[[validate_col]],             na.rm = TRUE),
            median_expr = median(.data[[validate_col]],           na.rm = TRUE),
            q25_expr    = quantile(.data[[validate_col]], 0.25,   na.rm = TRUE),
            q75_expr    = quantile(.data[[validate_col]], 0.75,   na.rm = TRUE),
            n           = dplyr::n(),
            .groups     = "drop"
          ) %>%
          tidyr::pivot_wider(
            names_from  = model_col,
            values_from = c(mean_expr, median_expr, q25_expr, q75_expr, n)
          ) %>%
          dplyr::mutate(
            Gene             = sv,
            Facet            = fg,
            Validate_col     = validate_col,
            mean_correct     = mean_expr_HIGH   > mean_expr_LOW,
            median_correct   = median_expr_HIGH > median_expr_LOW,
            q75_correct      = q75_expr_HIGH    > q75_expr_LOW,
            ALL_correct      = mean_correct & median_correct & q75_correct,
            FLAG             = ifelse(ALL_correct, "OK", "⚠ CLASSIFICATION INVERTED")
          ) %>%
          dplyr::select(
            Gene, Facet, Validate_col,
            n_HIGH, n_LOW,
            mean_HIGH   = mean_expr_HIGH,   mean_LOW   = mean_expr_LOW,   mean_correct,
            median_HIGH = median_expr_HIGH, median_LOW = median_expr_LOW, median_correct,
            q75_HIGH    = q75_expr_HIGH,    q75_LOW    = q75_expr_LOW,    q75_correct,
            ALL_correct, FLAG
          )
        
        classification_evidence <- dplyr::bind_rows(classification_evidence, ev)
      }
    }
    
    # Surface any failures immediately in the console — don't make user
    # open Excel to discover a silent misclassification.
    n_flags <- sum(classification_evidence$FLAG != "OK", na.rm = TRUE)
    if (n_flags > 0) {
      log_warn(sample = "", step = "plot_survival",
               msg = glue::glue(
                 "{n_flags} classification(s) flagged as potentially inverted — ",
                 "HIGH group has lower expression than LOW group. ",
                 "Check 'classification_evidence' sheet in Excel. ",
                 "This may indicate a scoring method issue or extreme cutoff."
               ))
    } else {
      log_info(sample = "", step = "plot_survival",
               msg = glue::glue(
                 "Classification validation passed — HIGH > LOW confirmed across ",
                 "all {nrow(classification_evidence)} gene x facet combinations."
               ))
    }
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Calculate Survival Statistics
  # ═══════════════════════════════════════════════════════════════════════════
  # For each gene x facet combination:
  #   (a) Cox model → HR, 95% CI, Wald p-value. LOW always reference for
  #       expression-based survival (HR > 1 = high expression → worse survival).
  #   (b) emmeans pairwise contrasts → HR in both directions (A/B and B/A)
  #       so Excel is self-contained without manual reciprocal arithmetic.
  #   (c) Seven non-parametric log-rank variants — span early-event sensitivity
  #       (Gehan-Breslow) through late-event (Fleming-Harrington) and standard
  #       log-rank. All 7 appended as columns to both cox_df and emmeans_df.
  # Stats always computed regardless of plot_curve setting.
  
  # Method codes and readable names defined once — identical for every gene/facet.
  np_methods <- c("survdiff", "1", "n", "sqrtN", "S1", "S2", "FH_p=1_q=1")
  np_names   <- c("p_logrank", "p_logrank_late", "p_gehan_breslow",
                  "p_tarone_ware", "p_peto_peto", "p_mod_peto",
                  "p_fleming_harrington")
  
  # Defined once outside the loop — re-defining inside would create a new
  # function object on every gene x facet iteration, wasteful with many genes.
  .calc_np_pvals <- function(contrast_str, df_pair, model_col,
                             time_col, status_col) {
    # emmeans sometimes wraps group names in parentheses e.g. "(Ba/Sq) / LOW"
    # when group names contain slashes. Strip them before splitting on " / ".
    contrast_str <- gsub("\\(|\\)", "", strsplit(contrast_str, " / ", fixed = TRUE)[[1]])
    g1  <- sub(".* / ", "", contrast_str)
    g2  <- sub(" / .*", "", contrast_str)
    df2 <- dplyr::filter(df_pair, .data[[model_col]] %in% c(g1, g2))
    so  <- survival::Surv(time   = df2[[time_col]],
                          event  = df2[[status_col]],
                          type   = "right",
                          origin = 0)
    fit <- survminer::surv_fit(
      formula = as.formula(paste("so ~", model_col)),
      data    = df2
    )
    sapply(np_methods, function(m) {
      tryCatch(
        survminer::surv_pvalue(fit            = fit,
                               method         = m,
                               test.for.trend = FALSE,
                               combine        = FALSE)[[2]],
        error = function(e) NA_real_
      )
    })
  }
  
  all_cox_df     <- tibble::tibble()
  all_emmeans_df <- tibble::tibble()
  all_plots      <- list()   # stores per-facet objects for Section 5
  
  for (sv in stratify_var) {
    
    model_col <- paste0("model_", sv)
    
    for (fg in facet_groups) {
      
      # ── Subset to this facet (or full dataset if no faceting) ─────────────
      facet_df <- if (is.na(fg)) surv_df else
        dplyr::filter(surv_df, .data[[facet_var]] == fg)
      
      # Drop samples with NA model column — arise when optimal cutpoint fails
      # (e.g. near-zero variance expression). Must filter BEFORE fitting
      # surv_curve so the fit row count matches facet_df row count exactly.
      # Mismatch causes ggsurvplot() "differing number of rows" error in Section 5.
      facet_df <- dplyr::filter(facet_df, !is.na(.data[[model_col]]))
      
      # Need at least 2 distinct groups to fit any survival model.
      n_groups <- dplyr::n_distinct(facet_df[[model_col]])
      if (n_groups < 2) {
        log_warn(sample = "", step = "plot_survival",
                 msg = glue::glue(
                   "Skipping '{sv}' / facet '{fg}' — fewer than 2 groups ",
                   "after binning ({n_groups} group found)."
                 ))
        next
      }
      
      # ── 4a. Survival object and formula ───────────────────────────────────
      surv_object  <- survival::Surv(time   = facet_df[[time_col]],
                                     event  = facet_df[[status_col]],
                                     type   = "right",
                                     origin = 0)
      surv_formula <- as.formula(paste("surv_object ~", model_col))
      
      # ── 4b. Cox proportional hazards model ────────────────────────────────
      # LOW set as reference so HR > 1 means high expression → worse survival,
      # the most natural biological interpretation. For metadata-based survival
      # first alphabetical level becomes reference by default.
      facet_df[[model_col]] <- factor(facet_df[[model_col]])
      if ("LOW" %in% levels(facet_df[[model_col]])) {
        facet_df[[model_col]] <- relevel(facet_df[[model_col]], ref = "LOW")
      }
      
      cox_model <- survival::coxph(formula = surv_formula, data = facet_df)
      
      # Degenerate model (zero events, perfect separation, etc.) — cox_df
      # itself is meaningless so skip the entire iteration.
      if (anyNA(coef(cox_model))) {
        log_warn(sample = "", step = "plot_survival",
                 msg = glue::glue(
                   "Skipping '{sv}' / facet '{fg}' — degenerate Cox model ",
                   "(NA coefficients). Likely zero or near-zero events."
                 ))
        next
      }
      
      cox_coef  <- summary(cox_model)$coefficients
      cox_ci    <- as.data.frame(confint(cox_model))
      baseline  <- levels(facet_df[[model_col]])[1]
      
      cox_df <- data.frame(
        Gene      = sv,
        Facet     = fg,
        Target    = gsub(paste0("^", model_col), "", rownames(cox_coef)),
        Reference = baseline,
        HR        = exp(cox_coef[, "coef"]),
        CI_lower  = exp(cox_ci[, 1]),
        CI_upper  = exp(cox_ci[, 2]),
        pval      = cox_coef[, "Pr(>|z|)"],
        stringsAsFactors = FALSE
      ) %>%
        dplyr::mutate(contrast = paste0(Target, " / ", Reference)) %>%
        dplyr::select(Gene, Facet, contrast, Target, Reference,
                      HR, pval, CI_lower, CI_upper) %>%
        tibble::remove_rownames()
      
      # ── 4c. emmeans pairwise contrasts ────────────────────────────────────
      # emmeans operates on log-hazard scale then back-transforms to HR.
      # CI bounds clamped away from 0 and Inf — extreme CIs arise from sparse
      # groups and would cause downstream arithmetic errors (1/0, log(0)).
      
      # Requires >= 2 events for a stable model. With 0–1 events emmeans
      # errors with the "bhat" slot class mismatch seen above.
      # When skipped, emmeans_df is left empty — bind_rows handles tibble()
      # gracefully so all_emmeans_df accumulation is unaffected.
      n_events <- sum(facet_df[[status_col]], na.rm = TRUE)
      
      if (n_events < 2) {
        log_warn(sample = "", step = "plot_survival",
                 msg = glue::glue(
                   "Skipping emmeans for '{sv}' / facet '{fg}' — ",
                   "only {n_events} event(s), model too sparse for ",
                   "reliable pairwise contrasts."
                 ))
        emmeans_df <- tibble::tibble()
        
      } else {
        
        emm      <- emmeans::emmeans(object = cox_model,
                                     specs  = as.formula(paste0("~", model_col)))
        pairwise <- emmeans::contrast(object = emm, method = "pairwise",
                                      type   = "response")
        pw_ci    <- stats::confint(pairwise)
        
        pw_df <- dplyr::left_join(
          x  = as.data.frame(pairwise) %>%
            dplyr::select(contrast, ratio, p.value),
          y  = as.data.frame(pw_ci) %>%
            dplyr::select(contrast, asymp.LCL, asymp.UCL),
          by = "contrast"
        ) %>%
          dplyr::rename(HR   = ratio,
                        pval = p.value) %>%
          dplyr::mutate(Target    = sub(" / .*", "", contrast),
                        Reference = sub(".* / ", "", contrast),
                        CI_lower  = pmax(asymp.LCL, 1e-10),
                        CI_upper  = pmin(asymp.UCL, 1e10),
                        Gene      = sv,
                        Facet     = fg) %>%
          dplyr::select(-asymp.LCL, -asymp.UCL)
        
        # Reversed contrasts: reciprocal HR and swapped CI bounds so Excel is
        # self-contained — both A/B and B/A available without manual arithmetic.
        pw_reversed <- pw_df %>%
          dplyr::mutate(Target_old    = Target,
                        Reference_old = Reference,
                        CI_lower_old  = CI_lower,
                        CI_upper_old  = CI_upper,
                        Target        = Reference_old,
                        Reference     = Target_old,
                        contrast      = paste0(Target, " / ", Reference),
                        HR            = 1 / HR,
                        CI_lower      = pmax(1 / CI_upper_old, 1e-10),
                        CI_upper      = pmin(1 / CI_lower_old, 1e10)) %>%
          dplyr::select(-Target_old, -Reference_old, -CI_lower_old, -CI_upper_old)
        
        emmeans_df <- dplyr::bind_rows(pw_df, pw_reversed) %>%
          dplyr::select(Gene, Facet, contrast, Target, Reference,
                        HR, pval, CI_lower, CI_upper)
      }
      
      # ── 4d. Non-parametric p-values (7 log-rank variants) ─────────────────
      # sapply returns methods x contrasts matrix; t() transposes to
      # contrasts x methods so it binds cleanly as extra columns.
      # Cox np always computed — cox_df is always valid at this point.
      # emmeans np only computed when emmeans_df is non-empty.
      np_cox <- t(sapply(cox_df$contrast,
                         .calc_np_pvals,
                         df_pair    = facet_df,
                         model_col  = model_col,
                         time_col   = time_col,
                         status_col = status_col))
      colnames(np_cox) <- np_names
      cox_df <- dplyr::bind_cols(cox_df, as.data.frame(np_cox))
      
      if (nrow(emmeans_df) > 0) {
        np_emm <- t(sapply(emmeans_df$contrast,
                           .calc_np_pvals,
                           df_pair    = facet_df,
                           model_col  = model_col,
                           time_col   = time_col,
                           status_col = status_col))
        
        colnames(np_emm) <- np_names
        emmeans_df <- dplyr::bind_cols(emmeans_df, as.data.frame(np_emm))
      }
      
      all_cox_df     <- dplyr::bind_rows(all_cox_df,     cox_df)
      all_emmeans_df <- dplyr::bind_rows(all_emmeans_df, emmeans_df)
      
      # ── 4e. Store per-facet objects for Section 5 ─────────────────────────
      # facet_df stored here is already filtered (NA model rows dropped above)
      # so Section 5 uses the same rows that went into the Cox model — critical
      # for ggsurvplot() row count consistency.
      surv_curve       <- survminer::surv_fit(formula = surv_formula,
                                              data    = facet_df)
      key              <- paste0(sv, "__", fg)
      all_plots[[key]] <- list(
        facet_df     = facet_df,
        surv_formula = surv_formula,
        surv_curve   = surv_curve,
        cox_df       = cox_df,
        n_groups     = n_groups,
        sv           = sv,
        fg           = fg
      )
    }
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5: Plot KM Curves
  # ═══════════════════════════════════════════════════════════════════════════
  # Skipped entirely when plot_curve = FALSE — stats always computed in Section 4.
  # When plot_curve = TRUE:
  #   — one KM curve (+ optional risk table) produced per gene x facet
  #   — facets for the same gene arranged in sqrt-based grid on one PDF page
  #   — all genes printed sequentially into a single multi-page PDF
  #
  # Color assignment:
  #   Expression / signature : HIGH = "#C10020" (deep red), LOW = "#00538A"
  #                            (deep blue) — standard survival publication convention.
  #                            Middle bins from CUSTOM_PALETTE position 3+ to avoid
  #                            clashing with HIGH/LOW anchors.
  #   Metadata-based         : CUSTOM_PALETTE in alphabetical level order.
  
  km_plot_list <- list()
  
  if (isTRUE(plot_curve) && length(all_plots) > 0) {
    
    for (sv in stratify_var) {
      
      model_col   <- paste0("model_", sv)
      facet_keys  <- grep(paste0("^", sv, "__"), names(all_plots), value = TRUE)
      facet_plots <- list()
      
      for (key in facet_keys) {
        
        info     <- all_plots[[key]]
        facet_df <- info$facet_df
        fg       <- info$fg
        
        # ── 5a. Assign colors ─────────────────────────────────────────────
        group_levels <- levels(facet_df[[model_col]])
        if (surv_type %in% c("expression", "signature")) {
          base_colors <- c(HIGH = "#C10020", LOW = "#00538A")
          mid_levels  <- setdiff(group_levels, c("HIGH", "LOW"))
          mid_colors  <- setNames(CUSTOM_PALETTE[seq(3, 2 + length(mid_levels))],
                                  mid_levels)
          palette_use <- c(base_colors, mid_colors)[group_levels]
        } else {
          palette_use <- setNames(CUSTOM_PALETTE[seq_along(group_levels)],
                                  group_levels)
        }
        
        # ── 5b. X-axis breaks ─────────────────────────────────────────────
        # Snap to time-unit multiples: 12 for months (annual), 30 for days
        # (monthly), 1 for years. Keeps labels round and clinically meaningful.
        max_time  <- max(facet_df[[time_col]], na.rm = TRUE)
        snap_unit <- switch(time_units, months = 12, days = 30, years = 1)
        n_snap    <- max(floor(max_time / 10 / snap_unit) * snap_unit, 1)
        breaks    <- if (max_time %/% n_snap <= 10) n_snap else n_snap + snap_unit
        x_upper   <- ceiling(max_time / breaks) * breaks
        
        # ── 5c. Plot title ────────────────────────────────────────────────
        auto_title <- paste(
          na.omit(c(sv, if (!is.na(fg)) as.character(fg))),
          collapse = " | "
        )
        plot_title <- title %||% auto_title
        
        # ── 5d. KM curve via ggsurvplot ───────────────────────────────────
        # surv_object, surv_formula, and surv_curve are all recreated fresh here
        # rather than reusing info$surv_curve from Section 4 for two reasons:
        #   (1) ggsurvplot() resolves surv_object by name from the fit's call slot —
        #       if surv_object no longer exists in the current scope it silently uses
        #       the wrong data.
        #   (2) ggsurvplot() maps legend.labs and palette positionally to strata in
        #       the order they appear in the fit object. If strata order in
        #       info$surv_curve doesn't exactly match levels(facet_df[[model_col]])
        #       — which can happen with Levine scores where the score distribution
        #       affects initial factor ordering — labels and colors map to the wrong
        #       curves. Refitting here guarantees strata order always matches factor
        #       levels exactly.
        
        surv_object  <- survival::Surv(time   = facet_df[[time_col]],
                                       event  = facet_df[[status_col]],
                                       type   = "right",
                                       origin = 0)
        surv_formula <- as.formula(paste("surv_object ~", model_col))
        surv_curve   <- survminer::surv_fit(formula = surv_formula,
                                            data    = facet_df)
        
        # legend.labs must use levels() not sort() — ggsurvplot() maps labels
        # to strata in factor level order, not alphabetical order. Mismatch
        # causes "differing number of rows" error in the risk table build step.
        surv_plot <- survminer::ggsurvplot(
          fit                   = surv_curve,
          pval                  = FALSE,
          palette               = palette_use,
          linetype              = "solid",
          size                  = 1.5,
          legend                = "top",
          legend.title          = ifelse(grepl("sig_score", model_col),
                                         "Signature Score", sv),
          legend.labs           = levels(facet_df[[model_col]]),
          break.time.by         = breaks,
          xlab                  = glue::glue("Time ({stringr::str_to_title(time_units)})"),
          ylab                  = "Survival Probability",
          title                 = plot_title,
          conf.int              = conf_interval,
          conf.int.style        = "ribbon",
          conf.int.alpha        = 0.3,
          risk.table            = plot_risk_table,
          risk.table.title      = "Number at risk",
          risk.table.y.text.col = TRUE,
          risk.table.pos        = "out",
          censor                = TRUE,
          censor.shape          = "|",
          censor.size           = 4
        )
        
        # Align x-axis limits so curve and risk table tick marks line up exactly.
        surv_plot$plot  <- surv_plot$plot  +
          ggplot2::coord_cartesian(xlim = c(0, x_upper), clip = "off")
        surv_plot$table <- surv_plot$table +
          ggplot2::coord_cartesian(xlim = c(0, x_upper), clip = "off")
        
        # ── 5e. Subtitle annotation ───────────────────────────────────────
        # Subtitle lives in ggplot margin — never overlaps curves regardless
        # of their shape, unlike annotate("text") which uses data coordinates.
        # Method (log-rank) omitted — universally assumed in publications.
        # Cutoff scope only shown when facet_var provided — the global vs
        # facet-specific distinction is only meaningful in that context.
        #
        # 2 groups : p-value and HR [95% CI] — plot is self-contained.
        # >2 groups: overall log-rank p-value only — pairwise HRs in Excel.
        cutoff_scope <- if (!is.null(facet_var)) {
          if (isTRUE(global_cutoff)) " (global)" else " (facet-specific)"
        } else {
          ""
        }
        
        if (info$n_groups == 2) {
          p_ann     <- formatC(info$cox_df$pval[1], format = "e", digits = 1)
          hr_ann    <- round(info$cox_df$HR[1], 2)
          ci_lo     <- round(info$cox_df$CI_lower[1], 2)
          ci_hi     <- round(info$cox_df$CI_upper[1], 2)
          ann_label <- glue::glue("p = {p_ann} | HR = {hr_ann} [{ci_lo}, {ci_hi}]")
          if (!is.null(facet_var)) {
            ann_label <- glue::glue("{ann_label}\nCutoff = {cutoff_method}{cutoff_scope}")
          }
        } else {
          overall_p <- tryCatch(
            formatC(
              survminer::surv_pvalue(info$surv_curve, method = "survdiff")[[2]],
              format = "e", digits = 1
            ),
            error = function(e) "NA"
          )
          ann_label <- glue::glue("Overall log-rank p = {overall_p}")
          if (!is.null(facet_var)) {
            ann_label <- glue::glue("{ann_label}\nCutoff = {cutoff_method}{cutoff_scope}")
          }
        }
        
        surv_plot$plot <- surv_plot$plot +
          ggplot2::labs(subtitle = ann_label) +
          ggplot2::theme(plot.subtitle = ggplot2::element_text(size  = 9,
                                                               hjust = 0,
                                                               face  = "plain",
                                                               color = "grey30"))
        surv_plot$table <- surv_plot$table +
          ggplot2::theme(axis.title.y = ggplot2::element_blank())
        
        # ── 5f. Combine curve and risk table with cowplot ──────────────────
        # align = "v" not "hv" — horizontal alignment ("h") distorts the
        # relative widths of curve and risk table panels in some survminer
        # versions, pushing the risk table text outside the plot area.
        p_combined <- cowplot::plot_grid(
          plotlist    = surv_plot,
          align       = "v",
          axis        = "tblr",
          nrow        = 2,
          ncol        = 1,
          rel_heights = c(1, 0.3)
        )
        
        facet_plots[[key]] <- p_combined
      }
      
      # ── 5g. Arrange all facets for this gene into a sqrt grid ──────────
      # ceiling(sqrt(n)) gives the most square layout for any n:
      # 1→1x1, 2→2x1, 3→2x2, 4→2x2, 6→3x2, 9→3x3.
      # Each panel is 6x8 inches so the page expands proportionally.
      n_facets  <- length(facet_plots)
      ncol_grid <- ceiling(sqrt(n_facets))
      nrow_grid <- ceiling(n_facets / ncol_grid)
      
      gene_page <- cowplot::plot_grid(
        plotlist = facet_plots,
        ncol     = ncol_grid,
        nrow     = nrow_grid,
        align    = "hv",
        axis     = "tblr"
      )
      
      km_plot_list[[sv]] <- list(
        plot   = gene_page,
        page_w = 6 * ncol_grid,
        page_h = 8 * nrow_grid
      )
    }
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 6: Save Outputs
  # ═══════════════════════════════════════════════════════════════════════════
  # Excel always saved when output_dir provided — even when plot_curve = FALSE
  # — because stats are the primary deliverable for batch gene runs.
  # PDF saved only when plot_curve = TRUE and at least one plot was produced.
  # cairo_pdf() used over ggsave() — natively supports multi-page output and
  # handles font embedding and transparency for publication-quality PDFs.
  # Page size fixed for entire PDF — all genes share facet_groups so all
  # pages have identical dimensions.
  
  if (!is.null(output_dir)) {
    
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    # Sanitise filename: replace non-alphanumeric runs with underscores and
    # strip leading/trailing underscores introduced by the substitution.
    safe_filename <- gsub("[^[:alnum:]]+", "_", filename)
    safe_filename <- gsub("(^_+|_+$)", "", safe_filename)
    
    # ── 6a. Save Excel with all survival statistics ───────────────────────
    xlsx_file <- file.path(output_dir, paste0(safe_filename, ".xlsx"))
    wb        <- openxlsx::createWorkbook()
    
    openxlsx::addWorksheet(wb, "cox_stats")
    openxlsx::writeData(wb, "cox_stats", all_cox_df)
    
    openxlsx::addWorksheet(wb, "emmeans_stats")
    openxlsx::writeData(wb, "emmeans_stats", all_emmeans_df)
    
    openxlsx::addWorksheet(wb, "surv_data")
    openxlsx::writeData(wb, "surv_data", surv_df)
    
    # Classification evidence — only written for expression/signature types.
    # Provides auditable record that HIGH always had higher expression than LOW.
    # FLAG column = "OK" for all rows means classifications are trustworthy.
    # FLAG column = "⚠ CLASSIFICATION INVERTED" means results should not be trusted.
    if (nrow(classification_evidence) > 0) {
      openxlsx::addWorksheet(wb, "classification_evidence")
      openxlsx::writeData(wb, "classification_evidence", classification_evidence)
      
      # Highlight flagged rows in red so failures are visible without filtering
      flag_rows <- which(classification_evidence$FLAG != "OK") + 1  # +1 for header
      if (length(flag_rows) > 0) {
        red_style <- openxlsx::createStyle(fontColour = "#FFFFFF",
                                           bgFill     = "#C10020",
                                           textDecoration = "bold")
        openxlsx::addStyle(wb      = wb,
                           sheet   = "classification_evidence",
                           style   = red_style,
                           rows    = flag_rows,
                           cols    = 1:ncol(classification_evidence),
                           gridExpand = TRUE)
      }
    }
    
    openxlsx::saveWorkbook(wb, file = xlsx_file, overwrite = TRUE)
    
    log_info(sample = "", step = "plot_survival",
             msg = glue::glue("Survival statistics saved to: '{xlsx_file}'"))
    
    # ── 6b. Save PDF with KM curves ───────────────────────────────────────
    # Page dimensions from first gene's grid — all genes share facet_groups
    # so dimensions are identical across pages.
    if (isTRUE(plot_curve) && length(km_plot_list) > 0) {
      
      pdf_file <- file.path(output_dir, paste0(safe_filename, ".pdf"))
      first_pg <- km_plot_list[[1]]
      
      grDevices::cairo_pdf(filename = pdf_file,
                           width    = first_pg$page_w,
                           height   = first_pg$page_h,
                           onefile  = TRUE)
      for (sv_i in names(km_plot_list)) {
        print(km_plot_list[[sv_i]]$plot)
      }
      grDevices::dev.off()
      
      log_info(sample = "", step = "plot_survival",
               msg = glue::glue("KM curves saved to: '{pdf_file}'"))
    }
    
  } else {
    
    log_info(sample = "", step = "plot_survival",
             msg = "No output_dir provided — returning results list, skipping save.")
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 7: Return
  # ═══════════════════════════════════════════════════════════════════════════
  # Returns invisibly so caller can inspect without triggering auto-print.
  # Named list mirrors plot_volcano() and plot_heatmap() for consistency.
  
  return(invisible(list(
    plots         = km_plot_list,
    cox_stats     = all_cox_df,
    emmeans_stats = all_emmeans_df,
    surv_data     = surv_df
  )))
}


# ─── Internal helper: calculate expression cutoffs ────────────────────────────
# Called by plot_survival() Section 3. Not intended for direct use.
# Returns a two-column data.frame: Sample_ID and model_<stratify_var>.
# Middle bin samples dropped when show_all_bins = FALSE so downstream steps
# exclude them without explicit filtering. NA cutoffs (optimal method fails
# on near-zero variance data) propagate as NA model values — the n_groups < 2
# guard in Section 4 then skips this gene/facet cleanly.

.calc_survival_cutoffs <- function(df, stratify_var, time_col, status_col,
                                   cutoff_method, show_all_bins) {
  
  model_col   <- paste0("model_", stratify_var)
  expr_values <- df[[stratify_var]]
  
  qs <- stats::quantile(expr_values,
                        probs = c(0.25, 0.33, 0.50, 0.66, 0.75),
                        na.rm = TRUE)
  
  cutoffs <- switch(
    cutoff_method,
    "median"   = list(lower  = qs["50%"], upper  = qs["50%"], middle = NA),
    "tertile"  = list(lower  = qs["33%"], upper  = qs["66%"], middle = NA),
    "quartile" = list(lower  = qs["25%"], upper  = qs["75%"], middle = qs["50%"]),
    "optimal"  = tryCatch({
      res <- survminer::surv_cutpoint(data      = df,
                                      time      = time_col,
                                      event     = status_col,
                                      variables = stratify_var)
      list(lower  = res$cutpoint$cutpoint,
           upper  = res$cutpoint$cutpoint,
           middle = NA)
    }, error = function(e) {
      list(lower = NA, upper = NA, middle = NA)
    })
  )
  
  binned <- df %>%
    dplyr::filter(!is.na(.data[[stratify_var]])) %>%
    dplyr::mutate(
      !!model_col := dplyr::case_when(
        .data[[stratify_var]] >  cutoffs$upper                            ~ "HIGH",
        .data[[stratify_var]] <= cutoffs$lower                            ~ "LOW",
        !is.na(cutoffs$middle) & .data[[stratify_var]] >  cutoffs$middle ~ "MED_HIGH",
        !is.na(cutoffs$middle) & .data[[stratify_var]] <= cutoffs$middle ~ "MED_LOW",
        TRUE                                                              ~ "MID"
      ),
      !!model_col := factor(.data[[model_col]],
                            levels = c("LOW", "MED_LOW", "MID", "MED_HIGH", "HIGH"))
    ) %>%
    dplyr::select(Sample_ID, dplyr::all_of(model_col))
  
  if (!isTRUE(show_all_bins)) {
    binned <- dplyr::filter(binned, .data[[model_col]] %in% c("HIGH", "LOW"))
  }
  
  return(binned)
}
