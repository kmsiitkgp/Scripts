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

plot_tfs <- function(contrast, tf_xlsx,
                     metadata        = NULL,
                     col_annotations = NULL,
                     col_cluster_by  = "all",
                     n_hits          = 20,
                     output_dir = NULL) {

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Load Data
  # ═══════════════════════════════════════════════════════════════════════════
  # Each sheet in tf_xlsx is one analysis type. We load all at once and name
  # by sheet so we can dispatch by content later.

  metadata      <- load_smart(input_path = metadata)
  all_sheets    <- openxlsx::getSheetNames(tf_xlsx)
  analysis_list <- lapply(all_sheets, function(s) load_smart(tf_xlsx, sheet = s))
  names(analysis_list) <- all_sheets

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Select Top Hits  —  done ONCE, shared across all sheets
  # ═══════════════════════════════════════════════════════════════════════════
  # Top hits are selected from the activity score sheets only (DEG or Sample).
  # All other sheets (evidence dot plots) are filtered to these same top hits
  # so every plot in the output highlights the same TFs / pathways.
  #
  # Priority: DEG-level (condition == "t") > Sample-level (condition == sample IDs)
  # Reason  : DEG-level gives a single ranked list; sample-level needs variance.

  get_top_hits <- function(df, n_hits) {
    if ("condition" %in% colnames(df) && all(df$condition == "stat", na.rm = TRUE)) {
      # DEG-level: rank by absolute activity score within each direction
      df %>%
        dplyr::mutate(Direction = dplyr::case_when(score > 0 ~ "Upregulated",
                                                   score < 0 ~ "Downregulated")) %>%
        dplyr::filter(!is.na(Direction)) %>%
        dplyr::group_by(statistic, Direction) %>%
        dplyr::arrange(dplyr::desc(abs(score)), .by_group = TRUE) %>%
        dplyr::slice_head(n = n_hits) %>%
        dplyr::mutate(rank = dplyr::row_number()) %>%
        dplyr::ungroup()
    } else {
      # Sample-level: rank by variability across samples
      # (no single score to rank by, so SD across conditions is the best proxy)
      df %>%
        dplyr::group_by(statistic, source) %>%
        dplyr::summarise(std = stats::sd(score, na.rm = TRUE), .groups = "drop") %>%
        dplyr::slice_max(order_by = std, n = n_hits, with_ties = FALSE) %>%
        dplyr::select(statistic, source)
    }
  }

  # %||% is a null-coalescing operator (defined in UTILS.R): returns the left
  # side if non-NULL, otherwise the right side. DEG sheet takes priority —
  # if both DEG and Sample sheets exist, we rank by DEG-level scores because
  # they give a single ranked list per TF. Sample-level sheets are the fallback
  # for experiments where no formal contrast was run (e.g. exploratory cohorts).
  tf_df      <- analysis_list[["DEG_TF_Activity"]]      %||% analysis_list[["Sample_TF_Activity"]]
  pathway_df <- analysis_list[["DEG_Pathway_Activity"]] %||% analysis_list[["Sample_Pathway_Activity"]]

  top_tf_hits      <- if (!is.null(tf_df))      get_top_hits(tf_df,      n_hits) else NULL
  top_pathway_hits <- if (!is.null(pathway_df)) get_top_hits(pathway_df, n_hits) else NULL

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Per-Sheet Loop  —  classify → prep → plot → save
  # ═══════════════════════════════════════════════════════════════════════════

  # Set global theme once so all ggplot-based plots (bar, dot) inherit it.
  # Heatmap font is controlled separately via par() around plot_heatmap() calls.
  ggplot2::theme_set(theme_publication())

  # Defined outside loop so log_info at the end always has a valid value
  safe_contrast <- gsub("[^[:alnum:]_-]", "_", contrast)
  output_dir    <- file.path(output_dir %||% ".", safe_contrast)

  for (analysis_name in names(analysis_list)) {

    df <- analysis_list[[analysis_name]]
    if (is.null(df) || nrow(df) == 0) next

    # ── 3a. Classify sheet ──────────────────────────────────────────────────
    # We sniff columns (not sheet names) because column structure is controlled
    # by decoupleR and is guaranteed, whereas sheet names are user-defined.
    #
    # Column signatures:
    #   mor    column  → TF network weights       → TF evidence dot plot
    #   weight column  → Pathway network weights  → Pathway evidence dot plot
    #   condition == "t"    → single DEG t-stat   → Bar plot
    #   condition != "t"    → per-sample scores   → Heatmap

    has_cols <- function(...) all(c(...) %in% colnames(df))

    if (has_cols("source", "target", "stat", "mor")) {
      type         <- "Dot"
      is_pathway   <- FALSE
      metric       <- "mor"
      title_prefix <- "TF Evidence: "
      plot_colors  <- c("Concordant" = "#2196A6", "Discordant" = "#E07B39", "Neutral" = "#9E9E9E")

    } else if (has_cols("source", "target", "stat", "weight")) {
      type         <- "Dot"
      is_pathway   <- TRUE
      metric       <- "weight"
      title_prefix <- "Pathway Evidence: "
      plot_colors  <- c("Concordant" = "#2196A6", "Discordant" = "#E07B39", "Neutral" = "#9E9E9E")

    } else if (has_cols("source", "statistic", "condition", "score", "n_targets", "padj") &&
               all(df$condition == "stat", na.rm = TRUE)) {
      type         <- "Bar"
      is_pathway   <- grepl("pathway", analysis_name, ignore.case = TRUE)
      title_prefix <- "Top Hits "
      plot_colors  <- c("Upregulated" = "#E69F00", "Downregulated" = "#56B4E9")

    } else if (has_cols("source", "statistic", "condition", "score") &&
               !all(df$condition == "t", na.rm = TRUE)) {
      type         <- "Heatmap"
      is_pathway   <- grepl("pathway", analysis_name, ignore.case = TRUE)
      title_prefix <- "Top Hits "

    } else {
      next  # unrecognised sheet structure — skip silently
    }

    # ── 3b. Resolve top hits for this sheet ─────────────────────────────────
    top_df <- if (is_pathway) top_pathway_hits else top_tf_hits
    if (is.null(top_df) || nrow(top_df) == 0) next

    # ── 3c-i. BAR PLOT ───────────────────────────────────────────────────────
    if (type == "Bar") {

      # ~~ Wrangling ~~
      # Full grid of statistic × Direction × rank ensures every panel has the same
      # number of rows. NA slots are handled by plot_bar() which replaces them with
      # " " so ggplot retains the row as blank space rather than collapsing the panel.
      skeleton <- tidyr::expand_grid(statistic = unique(top_df$statistic),
                                     #statistic = unique(df$statistic),
                                     Direction = c("Upregulated", "Downregulated"),
                                     rank      = seq_len(n_hits))

      plot_df <- skeleton %>%
        dplyr::left_join(top_df, by = c("statistic", "Direction", "rank")) %>%
        #dplyr::left_join(top_bar_df, by = c("statistic", "Direction", "rank")) %>%
        dplyr::mutate(PlotLabel = ifelse(is.na(source),
                                         NA,
                                         paste0(stringr::str_wrap(gsub("_", " ", source), width = 30), " (", n_targets, ")"))) %>%
        dplyr::group_by(statistic) %>%
        # Sort lowest → highest so negative bars go left and NA rows are pushed
        # to the top where they render as blank space.
        dplyr::arrange(score, .by_group = TRUE) %>%
        dplyr::mutate(rank      = seq(from = n_hits * 2, to = 1, by = -1),
                      PlotLabel = factor(PlotLabel, levels = unique(PlotLabel))) %>%
        dplyr::ungroup()

      # ~~ Plotting ~~
      bar_plots <- list()
      for (stat in unique(plot_df$statistic)) {

        stat_df <- plot_df %>%
          dplyr::filter(statistic == stat) %>%
          dplyr::mutate(Significant = as.character(!is.na(padj) & padj <= 0.05))

        if (all(is.na(stat_df$source))) next   # skip if truly no data for this method

        bar_plots[[stat]] <- plot_bar(df                  = stat_df,
                                      x_var               = "score",
                                      y_var               = "PlotLabel",
                                      x_label             = "Activity Score",
                                      fill_var            = "Direction",
                                      fill_legend_title   = "Direction",
                                      alpha_var           = "Significant",
                                      alpha_values        = c("TRUE" = 1, "FALSE" = 0.35),
                                      alpha_legend_title  = "Significance",
                                      alpha_legend_labels = c("TRUE" = "FDR <= 0.05", "FALSE" = "Nominal p <= 0.05"),
                                      title               = paste0(title_prefix, " (", stat, ")"))
      }

      plots_out <- bar_plots

    # ── 3c-ii. HEATMAP ───────────────────────────────────────────────────────
    } else if (type == "Heatmap") {

      # ~~ Wrangling ~~
      # Pivot to wide: rows = sources (TFs/pathways), cols = samples, values = scores
      # Keep the "source" column intact as the row ID for plot_heatmap().
      # column_to_rownames() is no longer called here — plot_heatmap() handles
      # the conversion internally when row_id = "source" is passed.
      wide_df <- df %>%
        dplyr::semi_join(top_df %>% dplyr::select(statistic, source),
                         by = c("statistic", "source")) %>%
        tidyr::pivot_wider(id_cols     = c("statistic", "source"),
                           names_from  = "condition",
                           values_from = "score")

      # Filter metadata to only samples present in this matrix
      metadata_col <- NULL
      if (!is.null(col_annotations) && !is.null(metadata)) {
        metadata_col <- metadata %>%
          dplyr::filter(Sample_ID %in% colnames(wide_df)) %>%
          dplyr::select(dplyr::any_of(col_annotations), Sample_ID)
      }

      # ~~ Plotting ~~
      heatmap_plots <- list()
      for (stat in unique(wide_df$statistic)) {

        stat_df <- wide_df %>% dplyr::filter(statistic == stat) %>% dplyr::select(-statistic)

        if (nrow(stat_df) < 2) next   # pheatmap needs at least 2 rows to cluster

        # Why save and restore par()? plot_heatmap() internally calls par() to
        # set graphical parameters (font family, margins) but does NOT restore
        # them afterward. Without old_par / par(old_par), the font change leaks
        # into every subsequent plot in the same PDF.
        old_par <- par(family = "Helvetica")
        ph <- plot_heatmap(
          df               = stat_df,
          row_id           = "source",
          col_id           = "Sample_ID",
          data_type        = "auto",
          label_genes      = NULL,
          metadata_col     = metadata_col,
          metadata_row     = NULL,
          col_annotations  = col_annotations,
          row_annotations  = NULL,
          col_gap_by       = NULL,
          row_gap_by       = NULL,
          col_cluster_by   = col_cluster_by,
          row_cluster_by   = "all",
          plot_title       = paste0(title_prefix, " (", stat, ")"),
          heatmap_palette  = "rdbu",
          border_color     = NA,
          show_expr_legend = TRUE
        )
        par(old_par)

        heatmap_plots[[stat]] <- ph$ph$gtable
      }

      plots_out <- heatmap_plots

    # ── 3c-iii. DOT PLOT ─────────────────────────────────────────────────────
    } else if (type == "Dot") {

      # ~~ Wrangling ~~
      # Direction: does a target's network weight agree with its expression change?
      #   mor/weight > 0 and stat > 0  → gene activated as expected   → Concordant
      #   mor/weight > 0 and stat < 0  → gene suppressed despite pos weight → Discordant
      x_var <- if ("log2FoldChange" %in% colnames(df)) "log2FoldChange" else "stat"
      y_var <- if ("padj"  %in% colnames(df)) "padj"  else "abs_stat"

      # dot plot df has no statistic col — match on source only
      top_sources <- unique(top_df$source)

      plot_df <- df %>%
        dplyr::filter(source %in% top_sources) %>%
        dplyr::mutate(Direction = dplyr::case_when(
          .data[[metric]] * stat > 0 ~ "Concordant",
          .data[[metric]] * stat < 0 ~ "Discordant",
          TRUE                       ~ "Neutral"
        )) %>%
        dplyr::group_by(source) %>%
        dplyr::slice_max(abs(stat), n = 50, with_ties = FALSE) %>%   # top 50 targets per TF/pathway
        dplyr::ungroup()

      # Transform y-axis: -log10(padj) makes significant hits appear at the top.
      # Why +1e-10 offset? padj can be exactly 0 for very significant hits —
      # log10(0) is -Inf, which breaks axis scaling. The tiny offset avoids this
      # without meaningfully affecting the plot since -log10(1e-10) = 10,
      # which simply becomes the axis ceiling for the most significant genes.
      if (y_var == "padj")     plot_df$padj     <- -log10(plot_df$padj + 1e-10)
      if (y_var == "abs_stat") plot_df$abs_stat <- abs(plot_df$stat)

      # ~~ Plotting ~~
      dot_plots <- list()
      for (hit in unique(plot_df$source)) {

        hit_df <- plot_df %>% dplyr::filter(source == hit)

        dot_plots[[hit]] <- ggplot2::ggplot(
          data    = hit_df,
          mapping = aes(x = .data[[x_var]], y = .data[[y_var]], color = Direction)
        ) +
          ggplot2::geom_point(
            aes(size = abs(.data[[metric]])),
            position = ggplot2::position_jitter(width = 0.1, height = 0.1),
            alpha = 0.6
          ) +
          ggrepel::geom_text_repel(aes(label = target), size = 3, max.overlaps = 15) +
          ggplot2::geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.5) +
          ggplot2::geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.5) +
          ggplot2::scale_color_manual(values = plot_colors, na.translate = FALSE) +
          ggplot2::scale_size_continuous(range = c(0.5, 3)) +
          ggplot2::labs(
            title    = paste0(title_prefix, hit),
            subtitle = "Top 50 targets by |stat|",
            x        = x_var,
            y        = y_var,
            size     = "Weight / MOR"
          )
      }

      plots_out <- dot_plots
    }

    # ── 3d. Save ─────────────────────────────────────────────────────────────

    if (length(plots_out) == 0) next

    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    output_file <- file.path(output_dir, paste0(type, "_Plot_", analysis_name, ".pdf"))

    if (type == "Bar") {
      # All methods (ulm, mlm, viper …) side-by-side on one page
      ggplot2::ggsave(
        filename = output_file,
        plot     = cowplot::plot_grid(plotlist = plots_out,
                                      ncol = length(plots_out), nrow = 1,
                                      align = "hv"),
        device   = grDevices::cairo_pdf,
        width    = length(plots_out) * 6,
        height   = 9,    # ~15 rows × 2 directions × (15pt / 72pt) ≈ 8.3 in
        bg       = "white"
      )
    } else {
      # Heatmaps and dot plots: one panel per page in a multi-page PDF.
      # Heatmaps are stored as gtable objects (pheatmap's internal format) and
      # must be drawn with grid::grid.draw() — they cannot be passed to print().
      # ggplot objects (dot plots) use print(). inherits() check dispatches correctly.
      grDevices::cairo_pdf(output_file, width = 8, height = 11.5, onefile = TRUE)
      for (p in plots_out) {
        if (inherits(p, "gtable")) { grid::grid.newpage(); grid::grid.draw(p) } else print(p)
      }
      grDevices::dev.off()
    }
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Log and Return
  # ═══════════════════════════════════════════════════════════════════════════

  log_info(sample = "", step = "plot_tfs",
           msg    = glue::glue("TF plots saved to: '{output_dir}'"))

  return(invisible(NULL))
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

  # Split comma-separated col_annotations string into a character vector
  raw_col_ann  <- get_arg(4)
  col_ann_list <- if (!is.null(raw_col_ann)) trimws(strsplit(raw_col_ann, ",")[[1]]) else NULL

  plot_tfs(
    contrast        = get_arg(1),
    tf_xlsx         = get_arg(2),
    metadata        = get_arg(3),
    col_annotations = col_ann_list,
    col_cluster_by  = get_arg(5, default = "all"),
    n_hits          = as.numeric(get_arg(6, 20)),
    output_dir      = get_arg(7, default = ".")
  )
}
