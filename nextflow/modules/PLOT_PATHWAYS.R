#!/usr/bin/env Rscript

# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(openxlsx)
  library(tidyr)
  library(ggplot2)
  library(stringr)   # str_wrap for description labels
  library(cowplot)   # plot_grid for multi-panel PDF pages
  library(grid)      # grid.draw / grid.newpage for heatmaps
  library(grDevices) # cairo_pdf / dev.off
  library(scales)
})

# ---- 🛠️ 2. Smart Setup (Find & source UTILS.R) ----

get_utils_path <- function() {
  # 1. Windows dev machine
  if (.Platform$OS.type == "windows") {
    return("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/Scripts/nextflow/modules/UTILS.R")
  }
  # 2. Interactive Linux / macOS (HPC interactive session)
  if (interactive()) {
    return(file.path(getwd(), "modules", "UTILS.R"))
  }
  # 3. Non-interactive (Nextflow / Rscript)
  initial.options <- commandArgs(trailingOnly = FALSE)
  file_arg        <- grep("--file=", initial.options, value = TRUE)
  if (length(file_arg) == 0) stop("Cannot detect script path for UTILS.R!")
  script_dir      <- dirname(sub("--file=", "", file_arg))
  return(file.path(script_dir, "UTILS.R"))
}

utils_path <- get_utils_path()
if (!file.exists(utils_path)) stop(paste("❌ UTILS.R not found at:", utils_path))
source(utils_path)

# ---- 🧬 3. Function Definition (Always Runs) ----

plot_pathways <- function(contrast, pathway_xlsx, output_dir,
                          vst_nonblind    = NULL,
                          metadata        = NULL,
                          col_annotations = NULL,
                          col_cluster_by  = "all",
                          n_hits           = 10) {

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Load Data
  # ═══════════════════════════════════════════════════════════════════════════
  # Each sheet = one analysis (e.g. one gene set database like KEGG, Hallmark).
  # Methods sheet is metadata only — no plottable data, always excluded.

  expr_mat  <- load_smart(input_path = vst_nonblind)
  metadata  <- load_smart(input_path = metadata)

  all_sheets    <- setdiff(openxlsx::getSheetNames(pathway_xlsx), "Methods")
  analysis_list <- lapply(all_sheets, function(s) load_smart(pathway_xlsx, sheet = s))
  names(analysis_list) <- all_sheets

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Prep expr_mat
  # ═══════════════════════════════════════════════════════════════════════════
  # expr_mat comes in as a data.frame with SYMBOL + gene_id + gene_biotype cols.
  # We need it subsetted to only the samples relevant to this contrast.
  # expr_mat is only used for pathway-level heatmaps — skip if not provided.
  # Note: column_to_rownames() is no longer needed here — plot_heatmap() now
  # accepts a data.frame directly with row_id = "SYMBOL".

  if (!is.null(expr_mat) && !is.null(metadata)) {

    # Parse contrast string ("GroupA - GroupB") into group names
    # all.vars() handles any valid R expression, not just " - " format
    target_groups <- all.vars(base::parse(text = contrast))

    # Find which metadata column contains these group labels
    # (could be "Group", "Condition", "Treatment" etc — user-agnostic)
    group_col <- names(metadata)[
      sapply(metadata, function(x) any(target_groups %in% x))
    ][1]

    if (is.na(group_col)) {
      log_error(sample = "", step = "plot_pathways",
                msg = "Could not find contrast groups in any metadata column.")
    }

    relevant_samples <- metadata %>%
      dplyr::filter(.data[[group_col]] %in% target_groups) %>%
      dplyr::pull(Sample_ID)

    # Subset to contrast-relevant samples — keep SYMBOL + sample columns only
    expr_mat <- expr_mat %>%
      dplyr::select(dplyr::any_of(c("SYMBOL", relevant_samples)))

    metadata <- metadata %>% dplyr::filter(Sample_ID %in% relevant_samples)
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Per-Sheet Loop  —  detect → filter → plot → save
  # ═══════════════════════════════════════════════════════════════════════════

  # Set global ggplot theme once so all plots (bar, dot) inherit it.
  # Heatmap font handled separately via par() around plot_heatmap() call.
  ggplot2::theme_set(theme_publication())

  # Defined outside loop so log_info at the end always has a valid value
  safe_contrast <- gsub("[^[:alnum:]_-]", "_", contrast)
  output_path   <- file.path(output_dir, safe_contrast)

  for (analysis_name in names(analysis_list)) {

    df <- analysis_list[[analysis_name]]
    if (is.null(df) || nrow(df) == 0) next

    # ── 3a. Detect analysis type ─────────────────────────────────────────────
    # Column signature distinguishes tool used:
    #   NES        → GSEA  (pre-ranked, all genes used, positive + negative scores)
    #   GeneRatio  → ORA   (over-representation, only DEGs used, always positive)
    # This drives x-axis label, count column, and gene list column downstream.

    metric <- intersect(colnames(df), c("NES", "GeneRatio"))[1]

    if (is.na(metric)) {
      log_error(sample = "", step = "plot_pathways",
                msg = glue::glue("Sheet '{analysis_name}': no 'NES' or 'GeneRatio' column found."))
      next
    }

    x_label  <- base::switch(metric,
                              "NES"       = "Normalized Enrichment Score (NES)",
                              "GeneRatio" = "Gene Ratio")

    # counts = number of genes driving the enrichment signal
    counts   <- base::switch(metric,
                              "NES"       = "leading_edge_size",
                              "GeneRatio" = "k")

    # gene_col = column holding the actual gene symbols for heatmap plotting
    gene_col <- base::switch(metric,
                              "NES"       = "leadingEdge",
                              "GeneRatio" = "overlapGenes")

    #plot_colors <- c("Upregulated" = "#E69F00", "Downregulated" = "#56B4E9")

    # ── 3b. Filter + get top hits ────────────────────────────────────────────
    # Why pval and not padj for the significance filter?
    # padj is computed by BH correction across ALL pathways in a collection.
    # For small gene set collections (e.g. custom or niche databases), the
    # number of tested pathways can be so small that BH correction sets padj
    # to 1 for everything — even genuinely enriched pathways. Using padj here
    # would silently produce zero rows and skip the entire collection without
    # any warning. pval <= 0.05 is a looser filter that ensures we always show
    # something meaningful while still excluding noise.
    # FDR significance is still communicated visually via the alpha aesthetic
    # (padj <= 0.05 → fully opaque, nominal only → faded).
    #
    # Minimum 3 genes required — fewer is not biologically meaningful.
    # Top N taken per Collection × Direction to keep plots readable.
    #
    # NOTE: str_wrap applied at plot time only (not here) so Description
    # stays unmodified for gene lookups in the heatmap section below.

    top_df <- df %>%
      dplyr::filter(!is.na(.data[[metric]]),
                    .data[[counts]] >= 3,
                    pval <= 0.05) %>%
      dplyr::mutate(Direction = dplyr::case_when(
        .data[[metric]] > 0 ~ "Upregulated",
        .data[[metric]] < 0 ~ "Downregulated",
        TRUE                ~ NA_character_
      )) %>%
      dplyr::filter(!is.na(Direction)) %>%
      dplyr::group_by(Collection, Direction) %>%
      dplyr::arrange(dplyr::desc(abs(.data[[metric]])), .by_group = TRUE) %>%
      dplyr::slice_head(n = n_hits) %>%
      dplyr::mutate(rank = dplyr::row_number()) %>%
      dplyr::ungroup()

    if (nrow(top_df) == 0) {
      
      if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
      
      # Define the output filename clearly
      no_results_file <- file.path(output_path, paste0("No_Significant_Pathways_Found(", analysis_name, ").pdf"))
      
      # Open a PDF device
      grDevices::cairo_pdf(no_results_file, width = 8, height = 5)
      
      # Create a blank page with a clear message
      grid::grid.newpage()
      grid::grid.text(label = paste0("Analysis: ", analysis_name, "\n\n",
          "No pathways met the significance criteria:\n",
          "• p-value <= 0.05\n",
          "• Gene count >= 3"),
        gp = grid::gpar(fontsize = 14, col = "grey30"))
      
      # Close the device
      grDevices::dev.off()
      
      next
    }

    # ── 3c. Build skeleton ───────────────────────────────────────────────────
    
    # ~~ Wrangling ~~
    # Full grid of Collection × Direction × rank ensures every panel has the same
    # number of rows. NA slots are handled by plot_bar() which replaces them with
    # " " so ggplot retains the row as blank space rather than collapsing the panel.
    skeleton <- tidyr::expand_grid(Collection = unique(top_df$Collection),
      Direction  = c("Upregulated", "Downregulated"),
      rank       = seq_len(n_hits))
    
    plot_df <- skeleton %>%
      dplyr::left_join(top_df, by = c("Collection", "Direction", "rank")) %>%
      dplyr::mutate(PlotLabel = ifelse(is.na(Description),
                                       NA,
                                       paste0(stringr::str_wrap(gsub("_", " ", Description), width = 30), " (", .data[[counts]], ")"))) %>%
      dplyr::group_by(Collection) %>%
      # Sort lowest → highest so negative bars go left and NA rows are pushed
      # to the top where they render as blank space.
      dplyr::arrange(.data[[metric]], .by_group = TRUE) %>%
      dplyr::mutate(rank      = seq(from = n_hits * 2, to = 1, by = -1),
                    PlotLabel = factor(PlotLabel, levels = unique(PlotLabel))) %>%
      dplyr::ungroup()

    # ── 3d. Plot per Collection ──────────────────────────────────────────────

    bar_plots     <- list()
    dot_plots     <- list()
    heatmap_plots <- list()

    for (collection in unique(plot_df$Collection)) {

      coll_df <- plot_df %>% 
        dplyr::filter(Collection == collection) %>%
        dplyr::mutate(Significant = as.character(!is.na(padj) & padj <= 0.05))

      # nrow is always n_hits*2 (skeleton) — check if any real pathways exist
      if (all(is.na(coll_df$Description))) next

      # ── 3d-i. Bar Plot ────────────────────────────────────────────────────
      
      bar_plots[[collection]] <- plot_bar(df                  = coll_df,
                                          x_var               = metric,
                                          y_var               = "PlotLabel",
                                          x_label             = x_label,
                                          fill_var            = "Direction",
                                          fill_legend_title   = "Direction",
                                          alpha_var           = "Significant",
                                          alpha_values        = c("TRUE" = 1, "FALSE" = 0.35),
                                          alpha_legend_title  = "Significance",
                                          alpha_legend_labels = c("TRUE" = "FDR <= 0.05", "FALSE" = "Nominal p <= 0.05"),
                                          title               = collection)

      # ── 3d-ii. Dot Plot ───────────────────────────────────────────────────
      # Dot size = gene count (k or leading_edge_size) — conveys both
      # enrichment magnitude (x) and evidence strength (size) simultaneously.
      
      dot_plots[[collection]] <- plot_dot(df                  = coll_df,
                                          x_var               = metric,
                                          y_var               = "PlotLabel",
                                          x_label             = x_label,
                                          color_var           = "Direction",
                                          fill_var            = "Direction",
                                          size_var            = counts,
                                          size_legend_title   = "Gene Count",
                                          alpha_var           = "Significant",
                                          alpha_values        = c("TRUE" = 1, "FALSE" = 0.35),
                                          alpha_legend_title  = "Significance",
                                          alpha_legend_labels = c("TRUE" = "FDR <= 0.05", "FALSE" = "Nominal p <= 0.05"),
                                          title               = collection)

      # ── 3d-iii. Heatmap ───────────────────────────────────────────────────
      # One heatmap per significant pathway showing leading edge / overlap genes.
      # Skipped entirely if expr_mat not provided.

      if (is.null(expr_mat)) next

      metadata_col <- NULL
      if (!is.null(col_annotations) && !is.null(metadata)) {
        metadata_col <- metadata %>%
          dplyr::select(dplyr::any_of(col_annotations), Sample_ID)
      }

      # Use original unmodified Description (not WrappedLabel) for gene lookup
      for (pathway in unique(coll_df$Description[!is.na(coll_df$Description)])) {

        # ~~ Wrangling ~~
        # Parse gene list — stored as "/" separated string (clusterProfiler format)
        plot_genes <- coll_df %>%
          dplyr::filter(Description == pathway) %>%
          dplyr::pull(.data[[gene_col]]) %>%
          base::strsplit("/") %>%
          base::unlist() %>%
          base::trimws() %>%
          stats::na.omit() %>%
          .[. != ""] %>%
          unique() %>%
          intersect(expr_mat$SYMBOL)   # keep only genes present in expr_mat

        if (length(plot_genes) < 2) next  # pheatmap needs at least 2 rows

        # ~~ Plotting ~~
        # Why save and restore par()? plot_heatmap() internally calls par() to
        # set graphical parameters (font family, margins) but does NOT restore
        # them afterward. Without old_par / par(old_par), the font change leaks
        # into every plot drawn after this heatmap in the same PDF — bar and dot
        # plots would silently inherit the wrong font. par(family="Helvetica")
        # sets the font BEFORE plot_heatmap() so the heatmap matches the rest;
        # par(old_par) restores whatever was set before so nothing leaks out.
        old_par <- par(family = "Helvetica")
        ph <- plot_heatmap(
          df               = expr_mat %>% dplyr::filter(SYMBOL %in% plot_genes),
          row_id           = "SYMBOL",
          col_id           = "Sample_ID",
          data_type        = "counts_log",
          label_genes      = NULL,
          metadata_col     = metadata_col,
          metadata_row     = NULL,
          col_annotations  = col_annotations,
          row_annotations  = NULL,
          col_gap_by       = NULL,
          row_gap_by       = NULL,
          col_cluster_by   = col_cluster_by,
          row_cluster_by   = "all",
          plot_title       = stringr::str_wrap(pathway, width = 30),
          heatmap_palette  = "rdbu",
          border_color     = NA,
          show_expr_legend = TRUE
        )
        par(old_par)

        if (!is.null(ph)) {
          heatmap_plots[[paste0(collection, "_", pathway)]] <- ph$ph$gtable
        }
      }
    }

    # ── 3e. Save ─────────────────────────────────────────────────────────────

    if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

    # Bar and Dot: all collections on one page, 3 columns
    for (type in c("Bar", "Dot")) {

      plot_list <- if (type == "Bar") bar_plots else dot_plots
      if (length(plot_list) == 0) next

      n_rows      <- ceiling(length(plot_list) / 3)
      output_file <- file.path(output_path, paste0(type, "_Plot_Pathways(", analysis_name, ").pdf"))

      ggplot2::ggsave(
        filename = output_file,
        plot     = cowplot::plot_grid(plotlist = plot_list, ncol = 3, align = "hv"),
        device   = grDevices::cairo_pdf,
        width    = 3 * 6,
        height   = n_rows * 6,    # 15 points * 20 rows / 72 points ~ 4.16 inches
        bg       = "white"
      )
    }

    # Heatmaps: one per page in a multi-page PDF
    if (length(heatmap_plots) > 0) {
      output_file <- file.path(output_path, paste0("Heatmap_Plot_Pathways(", analysis_name, ").pdf"))
      grDevices::cairo_pdf(output_file, width = 8, height = 11.5, onefile = TRUE)
      for (ht in heatmap_plots) {
        grid::grid.newpage()
        grid::grid.draw(ht)
      }
      grDevices::dev.off()
    }
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Log and Return
  # ═══════════════════════════════════════════════════════════════════════════

  log_info(sample = "", step = "plot_pathways",
           msg    = glue::glue("Pathway plots saved to: '{output_path}'"))

  return(invisible(NULL))
}

# ---- 🚀 4. Smart Execution (Nextflow Only) ----

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  get_arg <- function(idx, default = NULL) {
    if (idx > length(args)) return(default)
    val <- args[idx]
    if (is.na(val) || val == "" || val == "null" || val == "NULL") return(default)
    return(val)
  }

  # Split comma-separated col_annotations string into a character vector
  raw_col_ann  <- get_arg(6)
  col_ann_list <- if (!is.null(raw_col_ann)) trimws(strsplit(raw_col_ann, ",")[[1]]) else NULL

  plot_pathways(
    contrast        = get_arg(1),
    pathway_xlsx    = get_arg(2),
    output_dir      = get_arg(3, "."),
    vst_nonblind    = get_arg(4),
    metadata        = get_arg(5),
    col_annotations = col_ann_list,
    col_cluster_by  = get_arg(7, "all"),
    n_hits           = as.numeric(get_arg(8, 10))
  )
}
