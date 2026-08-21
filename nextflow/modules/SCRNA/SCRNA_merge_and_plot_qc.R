# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(ggrastr)
  library(scales)
})

# ---- 🛠️ 2. Smart Setup (Find & source UTILS.R) ----

get_modules_dir <- function() {
  # 1. Windows dev machine
  if (.Platform$OS.type == "windows") {
    return("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/Scripts/nextflow/modules")
  }
  # 2. Interactive Linux / macOS (HPC interactive session)
  if (interactive()) {
    return(file.path(getwd(), "modules"))
  }
  # 3. Non-interactive (Nextflow / Rscript)
  initial.options <- commandArgs(trailingOnly = FALSE)
  file_arg        <- grep("--file=", initial.options, value = TRUE)
  if (length(file_arg) == 0) stop("Cannot detect modules directory path!")
  modules_dir      <- dirname(sub("--file=", "", file_arg))
  return(modules_dir)
}

find_utils_dir <- function(start_dir) {
  d <- normalizePath(start_dir, mustWork = FALSE)
  repeat {
    if (file.exists(file.path(d, "UTILS.R"))) return(d)
    parent <- dirname(d)
    if (parent == d) stop("❌ UTILS.R not found searching upward from: ", start_dir)
    d <- parent
  }
}

modules_dir <- get_modules_dir()
utils_path  <- file.path(find_utils_dir(modules_dir), "UTILS.R")
source(utils_path)

# ---- 🧬 3. Function Definition (Always Runs) ----

# ═══════════════════════════════════════════════════════════════════════════
# HELPER: plot_flag_rate_dotplot — generic engine, reusable beyond QC
# ═══════════════════════════════════════════════════════════════════════════
# NOT QC-specific: takes any data frame, a grouping column, and a vector of
# logical (or 0/1) "flag" columns. For each group x flag combination, computes
# the % of rows in that group where the flag is TRUE — encoded as dot size —
# and how that group's rate compares to OTHER groups on the SAME flag only —
# encoded as color (z-score across groups, clipped). NA flag values are
# excluded from the rate calculation (na.rm = TRUE); if every value for a
# group x flag combination is NA (e.g. a classifier that didn't run for a
# given sample's input_type), that dot is simply omitted.
#
# Candidate for promotion to UTILS.R if a second caller shows up later (e.g.
# a cluster-level junk/sparse audit dotplot) — kept local here for now since
# QC is currently the only user.

plot_flag_rate_dotplot <- function(data,
                                   group_var,
                                   flag_cols,
                                   scale_clip  = 2.5,
                                   max_dot_size = 15,
                                   title       = "Flag Rate",
                                   size_legend = "% Flagged",
                                   fill_legend = "Relative Rate\n(z-score)") {

  df_long <- data %>%
    dplyr::select(dplyr::all_of(c(group_var, flag_cols))) %>%
    tidyr::pivot_longer(cols      = dplyr::all_of(flag_cols),
                        names_to  = "flag",
                        values_to = "value")

  rates <- df_long %>%
    dplyr::group_by(.data[[group_var]], flag) %>%
    dplyr::summarise(pct = mean(value, na.rm = TRUE) * 100, .groups = "drop")

  # z-score per flag (row-wise group) — compares each group's rate to other
  # groups on the SAME flag only, so "unusually high MitoRatio failures" is
  # judged relative to other samples' MitoRatio failures, not against
  # unrelated flags. sd == 0 (or a single group) → all zscores set to 0
  # rather than NaN — no variation to flag as an outlier.
  rates <- rates %>%
    dplyr::group_by(flag) %>%
    dplyr::mutate(
      zscore = if (dplyr::n() > 1 && isTRUE(stats::sd(pct, na.rm = TRUE) > 0)) {
        as.numeric(scale(pct))
      } else {
        0
      },
      zscore = dplyr::case_when(
        zscore >  scale_clip ~  scale_clip,
        zscore < -scale_clip ~ -scale_clip,
        TRUE                 ~ zscore
      )
    ) %>%
    dplyr::ungroup()

  rates[[group_var]] <- factor(rates[[group_var]],
                               levels = sort(unique(as.character(rates[[group_var]]))))
  rates$flag          <- factor(rates$flag, levels = rev(flag_cols))

  ggplot2::ggplot(data    = rates,
                  mapping = ggplot2::aes(x    = .data[[group_var]],
                                         y    = flag,
                                         size = pct,
                                         fill = zscore)) +
    ggplot2::geom_point(shape = 21, colour = "black", stroke = 0.25, na.rm = TRUE) +
    ggplot2::scale_fill_distiller(type      = "div",
                                  palette   = "RdYlGn",
                                  direction = -1,
                                  limits    = c(-scale_clip, scale_clip),
                                  name      = fill_legend) +
    ggplot2::scale_size_area(max_size = max_dot_size, limits = c(0, 100), name = size_legend) +
    ggplot2::guides(size = ggplot2::guide_legend(
      override.aes = list(shape = 21, colour = "black", fill = "white", stroke = 0.75))) +
    ggplot2::labs(title = title, x = group_var, y = NULL) +
    theme_publication() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1))
}

# ═══════════════════════════════════════════════════════════════════════════
# HELPER: plot_qc — generates the QC plots from barcode metadata
# ═══════════════════════════════════════════════════════════════════════════

plot_qc <- function(metadata,
                    gene_cutoff    = 250,
                    umi_cutoff     = 500,
                    mito_cutoff    = 0.2,
                    novelty_cutoff = 0.8,
                    ribo_cutoff    = 0.05,
                    output_dir     = NULL) {

  # Ensure required columns exist
  required_cols <- c("Sample", "Barcode_QC", "nUMIs", "nGenes", "MitoRatio", "RiboRatio", "Novelty")
  missing_cols  <- setdiff(required_cols, colnames(metadata))
  if (length(missing_cols) > 0) {
    log_error(sample = "", step = "plot_qc",
              msg = glue::glue("Missing required columns: ",
                               "{paste(missing_cols, collapse = ', ')}"))
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Shared Aesthetics
  # ═══════════════════════════════════════════════════════════════════════════
  # QC levels ordered for consistent legend and fill ordering across all plots.
  # Colors chosen for maximum perceptual distance and colorblind friendliness.

  qc_levels  <- c("Empty Droplet", "Low Quality", "Doublet", "Singlet")
  fill_colors <- c("Empty Droplet" = "#FFC61E",
                   "Low Quality"   = "#D62728",
                   "Doublet"       = "#1F77B4",
                   "Singlet"       = "#2CA02C")

  # Explicitly arrange data frame and set factor levels based on sorted values.
  # This guarantees consistent layouts and tracking lines across all subplots.
  # Renamed to QC locally so downstream plotting code stays unchanged.
  metadata <- metadata %>%
    dplyr::rename(QC = Barcode_QC) %>%
    dplyr::arrange(Sample)
  metadata$Sample <- factor(metadata$Sample, levels = stringr::str_sort(unique(as.character(metadata$Sample))))
  metadata$QC     <- factor(metadata$QC,     levels = qc_levels)

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Helper — Cell Count Bar Chart
  # ═══════════════════════════════════════════════════════════════════════════
  # Log10 y-axis handles the large dynamic range between empty droplets
  # (potentially millions) and singlets (typically thousands).
  # tidyr::complete ensures all QC levels appear for every sample even if
  # absent — prevents bars from silently disappearing for clean samples.

  cell_qc <- function(meta) {

    df <- meta %>%
      dplyr::add_count(Sample, QC) %>%
      dplyr::distinct(Sample, QC, n) %>%
      tidyr::complete(Sample, QC, fill = list(n = 1))

    sample_levels    <- levels(meta$Sample)
    vline_positions  <- seq(1.5, length(sample_levels) - 0.5, by = 1)

    ggplot2::ggplot(data    = df,
                    mapping = ggplot2::aes(x = Sample, y = n, fill = QC)) +
      ggplot2::geom_bar(stat     = "identity",
                        position = ggplot2::position_dodge(0.9)) +
      # y = 1.2 positions the numbers uniformly right above the x-axis line
      ggplot2::geom_text(mapping  = ggplot2::aes(label = n, y = 1.2),
                         position = ggplot2::position_dodge(width = 0.9),
                         hjust    = 0,
                         angle    = 90,
                         size     = 2.5) +
      ggplot2::geom_vline(xintercept = vline_positions,
                          linetype   = "dotted",
                          color      = "gray50") +
      ggplot2::scale_y_log10(breaks = c(10, 100, 1000, 10000, 100000, 1000000),
                             labels = c("10", "100", "1K", "10K", "100K", "1M")) +
      ggplot2::scale_fill_manual(values = fill_colors) +
      ggplot2::coord_cartesian(ylim = c(1, 1e7), clip = "off", expand = FALSE) +
      ggplot2::labs(x = "Sample", y = "Cell Counts", title = "Number of Cells") +
      theme_publication() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1))
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Helper — Violin + Boxplot per Metric
  # ═══════════════════════════════════════════════════════════════════════════
  # Violin shows full distribution shape; overlaid narrow boxplot shows
  # median and IQR. cutoff draws a horizontal reference line matching the
  # threshold applied during calc_qc_and_identify_singlets.

  violin_qc <- function(meta, yvar, ylab, title, cutoff, ylog, ylim) {

    # Remove empty droplets — millions of barcodes with 0 UMIs/genes add no
    # distribution information to QC metrics and cause extreme memory overhead.
    # Empty droplet counts are already captured in the cell_qc bar chart.
    df <- meta %>%
      dplyr::filter(QC != "Empty Droplet") %>%
      dplyr::mutate(QC = base::droplevels(QC))

    sample_levels   <- levels(meta$Sample)
    vline_positions <- seq(1.5, length(sample_levels) - 0.5, by = 1)

    p <- ggplot2::ggplot(data    = df,
                         mapping = ggplot2::aes(x    = Sample,
                                                y    = .data[[yvar]],
                                                fill = QC)) +
      # Use scale = "width" else violins will be scaled by counts and look narrow
      ggplot2::geom_violin(position = ggplot2::position_dodge(0.9),
                           scale    = "width",
                           drop     = FALSE) +
      # Add show.legend = FALSE here to strip boxplot glyphs from the legend
      ggplot2::geom_boxplot(position     = ggplot2::position_dodge(0.9),
                            width        = 0.05,
                            outlier.size = 0.5,
                            show.legend  = FALSE) +
      ggplot2::geom_vline(xintercept = vline_positions,
                          linetype   = "dotted",
                          color      = "gray50") +
      ggplot2::scale_fill_manual(values = fill_colors) +
      ggplot2::labs(x = "Sample", y = ylab, title = title) +
      theme_publication() +
      # Overrides the key glyphs to use clean, solid color rectangles
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(linetype = 0, shape = NA))) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) +
      ggplot2::geom_hline(yintercept = cutoff, linetype   = "dashed") +
      ggplot2::coord_cartesian(clip = "off")

    if (ylog) {
      p <- p + ggplot2::scale_y_log10(limits = ylim, oob = scales::oob_keep)
    } else {
      p <- p + ggplot2::scale_y_continuous(limits = ylim, oob = scales::oob_keep)
    }

    return(p)
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Helper — Genes vs UMIs vs MitoRatio Scatter
  # ═══════════════════════════════════════════════════════════════════════════
  # Faceted by Sample × QC — reveals whether low quality cells cluster
  # distinctly from high quality cells and whether any QC category dominates.
  # Bottom-right quadrant (high UMIs, low genes) flags low-complexity cells
  # such as erythrocytes or dying cells with mitochondrial enrichment.

  gene_umi_mito_qc <- function(meta) {

    # Remove empty droplets — scatter of millions of barcodes with 0 UMIs/genes
    # obscures the real cell population and crashes rendering.
    df <- meta %>%
      dplyr::filter(QC != "Empty Droplet") %>%
      dplyr::mutate(QC = base::droplevels(QC))

    ggplot2::ggplot(data    = df,
                    mapping = ggplot2::aes(x     = nUMIs,
                                           y     = nGenes,
                                           color = MitoRatio)) +
      ggrastr::geom_point_rast(alpha = 0.3, size = 0.25, raster.dpi = 300) +
      ggplot2::geom_vline(xintercept = umi_cutoff,  linetype = "dashed") +
      ggplot2::geom_hline(yintercept = gene_cutoff, linetype = "dashed") +
      ggplot2::stat_smooth(method    = lm,
                           color     = "red",
                           se        = FALSE,
                           linewidth = 0.3) +
      ggplot2::facet_wrap(Sample ~ QC, ncol = 3, drop = FALSE) +
      ggplot2::scale_x_log10(breaks = c(100, 1000, 10000, 100000, 1e6),
                             labels = c("100", "1K", "10K", "100K", "1M")) +
      ggplot2::scale_y_log10(breaks = c(100, 1000, 10000, 100000),
                             labels = c("100", "1K", "10K", "100K")) +
      ggplot2::scale_color_viridis_c(option = "D", limits = c(0, 1)) +
      ggplot2::coord_cartesian(xlim = c(10, 1e6), ylim = c(10, 20000), clip = "off") +
      ggplot2::labs(x     = "Number of UMIs",
                    y     = "Number of Genes",
                    title = "UMIs vs Genes (Colored by MitoRatio)") +
      theme_publication() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1))
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5: Helper — QC Failure Dotplot (all criteria, per sample)
  # ═══════════════════════════════════════════════════════════════════════════
  # Sample-level diagnostic only — does NOT feed back into any filtering
  # decision. A cell that fails is excluded either way; this exists purely to
  # answer "does anything look technically off about a specific sample" by
  # comparing each sample's failure rate on each criterion against the OTHER
  # samples on that same criterion (see plot_flag_rate_dotplot above).
  #
  # RiboRatio is deliberately excluded — it is not currently an exclusion
  # gate (passes_qc in calc_qc_and_identify_singlets never uses ribo_cutoff),
  # only a monitored/reference metric. Including it here would imply it
  # drives exclusion when it doesn't.
  #
  # DropletUtils/CellRanger/scDblFinder "ND" (classifier didn't run for this
  # sample's input_type, or <50 high-quality cells) become NA rather than
  # FALSE — an untested barcode is not the same as a barcode that passed.

  qc_failure_dotplot <- function(meta) {

    df <- meta %>%
      dplyr::filter(QC != "Empty Droplet") %>%
      dplyr::mutate(
        nUMIs_fail         = nUMIs     < umi_cutoff,
        nGenes_fail        = nGenes    < gene_cutoff,
        MitoRatio_fail     = MitoRatio > mito_cutoff,
        Novelty_fail       = Novelty   < novelty_cutoff,
        DropletUtils_fail  = dplyr::case_when(DropletUtils == "ND" ~ NA,
                                              TRUE                 ~ DropletUtils == "Empty Droplet"),
        CellRanger_fail    = dplyr::case_when(CellRanger == "ND"   ~ NA,
                                              TRUE                 ~ CellRanger == "Empty Droplet"),
        scDblFinder_fail   = dplyr::case_when(scDblFinder == "ND"  ~ NA,
                                              TRUE                 ~ scDblFinder == "Doublet"))

    flag_cols <- c("nUMIs_fail", "nGenes_fail", "MitoRatio_fail", "Novelty_fail",
                   "DropletUtils_fail", "CellRanger_fail", "scDblFinder_fail")

    plot_flag_rate_dotplot(data        = df,
                           group_var   = "Sample",
                           flag_cols   = flag_cols,
                           title       = "QC Failure Rate by Sample (excl. Empty Droplets)",
                           size_legend = "% Failing",
                           fill_legend = "Relative Rate\n(z-score,\nvs other samples)")
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 6: Build Plot List
  # ═══════════════════════════════════════════════════════════════════════════
  # Cutoff values from function arguments flow into violin reference lines —
  # user sees exactly where the thresholds were applied relative to data.

  plot_list <- list(
    Cell_Counts                      = function(x) cell_qc(x),
    UMI_Distribution                 = function(x) violin_qc(x, "nUMIs",
                                                             "Number of UMIs",
                                                             "UMI Distribution",
                                                             umi_cutoff, TRUE,
                                                             c(1, 1e6)),
    Gene_Distribution                = function(x) violin_qc(x, "nGenes",
                                                             "Number of Genes",
                                                             "Gene Distribution",
                                                             gene_cutoff, TRUE,
                                                             c(1, 30000)),
    MitoRatio_Distribution           = function(x) violin_qc(x, "MitoRatio",
                                                             "MitoRatio",
                                                             "MitoRatio Distribution",
                                                             mito_cutoff, FALSE,
                                                             c(0, 1)),
    RiboRatio_Distribution           = function(x) violin_qc(x, "RiboRatio",
                                                             "RiboRatio",
                                                             "RiboRatio Distribution",
                                                             ribo_cutoff, FALSE,
                                                             c(0, 1)),
    Novelty_Distribution             = function(x) violin_qc(x, "Novelty",
                                                             "Novelty",
                                                             "Novelty Distribution",
                                                             novelty_cutoff, FALSE,
                                                             c(0, 1)),
    Genes_UMIs_MitoRatio             = function(x) gene_umi_mito_qc(x),
    QC_Failure_Dotplot               = function(x) qc_failure_dotplot(x)
  )

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 7: Generate and Save Plots
  # ═══════════════════════════════════════════════════════════════════════════

  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  }

  n_samples <- length(unique(metadata$Sample))
  for (plot_name in names(plot_list)) {

    p <- plot_list[[plot_name]](metadata)

    width_per_bar            <- 0.4 # 0.4 inches
    width_per_violin         <- 0.4 # 0.4 inches
    width_per_facet          <- 2   # 2 inches
    height_per_facet         <- 2   # 2 inches
    qc_levels_bar            <- 4   # Singlet, Doublet, Low Quality, Empty Droplet
    qc_levels_violin         <- 3   # Singlet, Doublet, Low Quality
    qc_levels_genes_umi_mito <- 3   # Singlet, Doublet, Low Quality
    margin                   <- 2   # 2 inches

    if (plot_name == "Cell_Counts") {
      pdf_width  <- min(n_samples * qc_levels_bar * width_per_bar + margin, 30)
      pdf_height <- 8
    } else if (plot_name == "Genes_UMIs_MitoRatio") {
      pdf_width  <- min(qc_levels_genes_umi_mito * width_per_facet + margin, 30)
      pdf_height <- min(n_samples * height_per_facet + margin, 30)
    } else if (plot_name == "QC_Failure_Dotplot") {
      height_per_flag <- 0.6 # inches per criterion row
      n_flags         <- 7   # nUMIs, nGenes, MitoRatio, Novelty, DropletUtils, CellRanger, scDblFinder
      pdf_width  <- min(n_samples * width_per_violin + margin + 2, 30) # +2 for legend
      pdf_height <- min(n_flags * height_per_flag + margin, 15)
    } else {
      pdf_width  <- min(n_samples * qc_levels_violin * width_per_violin + margin, 50)
      pdf_height <- 8
    }

    if (!is.null(output_dir)) {

      output_file <- file.path(output_dir, paste0("QC_", plot_name, ".pdf"))

      grDevices::cairo_pdf(filename = output_file,
                           width    = pdf_width,
                           height   = pdf_height,
                           onefile  = TRUE)
      print(p)
      grDevices::dev.off()

      log_info(sample = "", step = "plot_qc",
               msg = glue::glue("Saved: '{output_file}'"))
    }
  }

  log_info(sample = "", step = "plot_qc",
           msg = glue::glue("QC plots complete — {length(plot_list)} plots generated."))

  return(invisible(NULL))
}

# ═══════════════════════════════════════════════════════════════════════════
# MAIN FUNCTION: merge_and_plot_qc
# ═══════════════════════════════════════════════════════════════════════════

merge_and_plot_qc <- function(sample_rds_list, metadata_xlsx, output_dir = NULL) {

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Load and Validate Per-Sample Seurat Objects
  # ═══════════════════════════════════════════════════════════════════════════

  # Step 1: Check if it's already an R list (Interactive mode / R loop)
  if (is.list(sample_rds_list)) {
    rds_inputs <- sample_rds_list

    # Step 2: Check if it's a character string (Nextflow CLI mode)
  } else if (is.character(sample_rds_list) && length(sample_rds_list) == 1) {

    # Since your Nextflow script explicitly joins on commas, we split on commas.
    # If it doesn't have a comma (single file), strsplit still returns a vector of length 1.
    rds_inputs <- unlist(strsplit(sample_rds_list, ","))
    rds_inputs <- trimws(rds_inputs)

    # Step 3: Fallback for standard character vectors (e.g., c("path1", "path2"))
  } else {
    rds_inputs <- sample_rds_list
  }

  if (length(rds_inputs) == 0) {
    log_error(sample = "", step = "merge_and_plot_qc",
              msg = "No RDS files or Seurat objects provided.")
  }

  # IMPORTANT: Sorting cells inside seurat object is NOT straightforward
  # It is best to sort before object creation
  rds_inputs <- sort(rds_inputs)

  log_info(sample = "", step = "merge_and_plot_qc",
           msg = glue::glue("Processing {length(rds_inputs)} sample input(s)."))

  # Pre-allocate for metadata (using a list is faster than growing a dataframe)
  metadata_list <- list()
  seurat_list   <- list()

  # Loop using indices to safely handle both path strings and Seurat objects
  for (i in seq_along(rds_inputs)) {
    item <- rds_inputs[[i]]

    # load_smart handles either the path string or passes the object right through
    sample_seurat <- load_smart(input_path = item)

    # Harvest metadata BEFORE subsetting and extract sample ID
    sample_id <- as.character(unique(sample_seurat@meta.data$Sample))
    metadata_list[[sample_id]] <- sample_seurat@meta.data

    # Subset immediately to free up memory from the "Empty Droplets"
    sample_seurat <- subset(sample_seurat, Barcode_QC != "Empty Droplet")

    # Store the pruned version
    seurat_list[[sample_id]] <- sample_seurat

    log_info(sample = sample_id, step = "merge_and_plot_qc",
             msg = glue::glue("Loaded: {ncol(sample_seurat)} barcodes (all Barcode_QC categories)."))
  }

  # Combine metadata after the loop finishes
  pre_filtered_metadata <- dplyr::bind_rows(metadata_list)
  rm(metadata_list)
  gc()

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Merge Seurat Objects
  # ═══════════════════════════════════════════════════════════════════════════
  # add.cell.ids prefixes each barcode with its sample name to prevent
  # barcode collisions — Cell Ranger reuses the same barcode sequences across
  # samples so without prefixing, barcodes from different samples would
  # overwrite each other after merging.

  samples    <- names(seurat_list)
  sample_ids <- gsub(pattern = "\\.Spatial.*", replacement = "", x = samples)

  merged_seurat <- tryCatch({
    merge(x            = seurat_list[[1]],
          y            = seurat_list[-1],
          add.cell.ids = sample_ids,
          merge.data   = FALSE)
  }, error = function(e) {
    log_error(sample = "", step = "merge_and_plot_qc",
              msg = glue::glue("Merge failed: {e$message}"))
  })

  log_info(sample = "", step = "merge_and_plot_qc",
           msg = glue::glue("Merged object: {ncol(merged_seurat)} total barcodes across ",
                            "{length(samples)} samples."))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Remove HTO Assay
  # ═══════════════════════════════════════════════════════════════════════════
  # HTO assay is only relevant for demultiplexing — already captured in
  # HTO_Final metadata column. Remove before integration to avoid downstream
  # tools attempting to process it alongside RNA.

  if ("HTO" %in% Seurat::Assays(merged_seurat)) {
    merged_seurat[["HTO"]] <- NULL
    log_info(sample = "", step = "merge_and_plot_qc",
             msg = "HTO assay removed.")
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Load and Join Extra Metadata
  # ═══════════════════════════════════════════════════════════════════════════
  # Sample-level metadata (e.g. treatment, genotype, batch) is stored in a
  # user-provided xlsx/csv and joined here — BEFORE the raw object is saved,
  # so both pre_filtered_seurat.rds and filtered_seurat.rds carry the same annotations.

  extra_metadata <- tryCatch({
    load_smart(metadata_xlsx) %>%
      dplyr::select(-dplyr::any_of(c("Comments")))
  }, error = function(e) {
    log_warn(sample = "", step = "merge_and_plot_qc",
             msg = glue::glue("Failed to load metadata file: {e$message}. Skipping join."))
    data.frame()
  })

  if (nrow(extra_metadata) > 0) {

    meta_data <- merged_seurat@meta.data %>%
      tibble::rownames_to_column(var = "temp_barcode_id")

    # Check if the object contains spatial images instead of checking DefaultAssay
    is_spatial <- length(merged_seurat@images) > 0

    if (!is_spatial) {
      # For multiplexed samples, join on Sample + HTO combination (Sample_ID).
      # For non-multiplexed, join on Sample alone.
      meta_data <- meta_data %>%
        dplyr::mutate(
          Sample_ID = dplyr::case_when(
            "HTO_Final" %in% colnames(.) & HTO_Final != "ND" ~ paste0(Sample, "_", HTO_Final),
            TRUE              ~ Sample
          )
        ) %>%
        dplyr::left_join(extra_metadata, by = "Sample_ID")

    } else {
      # Spatial assays join by Slide — one metadata row per capture area
      meta_data <- meta_data %>%
        dplyr::left_join(extra_metadata, by = c("Sample" = "Slide"))
    }

    # Drop columns that are entirely NA — arise when metadata file has columns
    # not applicable to all samples (e.g. HTO columns for non-multiplexed runs)
    meta_data <- meta_data %>%
      dplyr::select(
        where(~ !all(is.na(.))),
        -dplyr::any_of(c("Initial filename", "Comments"))
      )

    merged_seurat@meta.data <- meta_data %>%
      tibble::column_to_rownames(var = "temp_barcode_id")

    log_info(sample = "", step = "merge_and_plot_qc",
             msg = "Extra metadata joined successfully.")

  } else {

    log_info(sample = "", step = "merge_and_plot_qc",
             msg = "No external metadata loaded — skipping join.")
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5: Filter to Singlets
  # ═══════════════════════════════════════════════════════════════════════════
  # Filtered object — singlet barcodes from every sample, all
  # classification columns intact. This is the object for downstream analysis

  filtered_seurat <- base::subset(x = merged_seurat, Barcode_QC == "Singlet")

  log_info(sample = "", step = "merge_and_plot_qc",
           msg = glue::glue("Filtered object: {ncol(filtered_seurat)} singlets retained ",
                            "out of {ncol(merged_seurat)} total barcodes."))

  # QC cutoffs used during CALC_QC_AND_IDENTIFY_SINGLETS — carried in @misc so
  # plot_qc always draws reference lines matching the actual thresholds applied
  merged_seurat@misc$qc_params       <- seurat_list[[1]]@misc$qc_params
  merged_seurat@misc$full_metadata   <- pre_filtered_metadata
  filtered_seurat@misc$qc_params     <- seurat_list[[1]]@misc$qc_params
  filtered_seurat@misc$full_metadata <- pre_filtered_metadata

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 6: Save Pre-Filtered and Filtered Seurat Object
  # ═══════════════════════════════════════════════════════════════════════════
  # Pre-Filtered object — non-empty barcodes from every sample, all
  # classification columns intact. This is the object to come back to if
  # cutoffs need revisiting, or to inspect empty/low-quality/doublet barcodes
  # directly (e.g. ambient RNA diagnostics) without rerunning anything upstream.

  if (!is.null(output_dir)) {

    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    pre_filtered_rds_file <- file.path(output_dir, "pre_filtered_seurat.rds")
    saveRDS(object = merged_seurat, file = pre_filtered_rds_file)

    log_info(sample = "", step = "merge_and_plot_qc",
             msg = glue::glue("Non-empty Seurat object saved to: '{pre_filtered_rds_file}'"))

    filtered_rds_file <- file.path(output_dir, "filtered_seurat.rds")
    saveRDS(object = filtered_seurat, file = filtered_rds_file)

    log_info(sample = "", step = "merge_and_plot_qc",
         msg = glue::glue("Filtered Seurat object saved to: '{filtered_rds_file}'"))

  } else {

    log_info(sample = "", step = "merge_and_plot_qc",
             msg = "No output_dir provided — skipping raw object save.")
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 7: Generate QC Plots
  # ═══════════════════════════════════════════════════════════════════════════
  # Read straight from merged_seurat@meta.data — full barcode-level data for
  # every sample is right here, no @misc reconstruction needed.

  qc_params <- seurat_list[[1]]@misc$qc_params

  plot_qc(metadata       = pre_filtered_metadata,
          gene_cutoff    = qc_params$gene_cutoff,
          umi_cutoff     = qc_params$umi_cutoff,
          mito_cutoff    = qc_params$mito_cutoff,
          novelty_cutoff = qc_params$novelty_cutoff,
          ribo_cutoff    = qc_params$ribo_cutoff,
          output_dir     = output_dir)

  log_info(sample = "", step = "merge_and_plot_qc",
           msg = "Merge and plot QC completed successfully.")

  return(invisible(list(pre_filtered_seurat = merged_seurat, filtered_seurat = filtered_seurat)))
}

# ---- 🚀 4. Smart Execution (Nextflow Only) ----

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)

  get_arg <- function(idx, default = NULL) {
    if (idx > length(args)) return(default)
    val <- args[idx]
    if (is.na(val) || val == "" || val == "null" || val == "NULL") return(default)
    return(val)
  }

  merge_and_plot_qc(
    sample_rds_list = get_arg(1),
    metadata_xlsx   = get_arg(2),
    output_dir    = get_arg(3, default = ".")
  )
}
