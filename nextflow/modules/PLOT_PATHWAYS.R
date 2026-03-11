#!/usr/bin/env Rscript

# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)      # Provides the %>% operator
  library(tibble)
  library(openxlsx)
  library(tidyr)
  library(ggplot2)    # For ggplot, aes, theme_classic, etc.
  library(stringr)    # For str_wrap (used in your Description cleaning)
  library(cowplot)    # For plot_grid (used to stitch the PDF together)
  library(grid)       # For grid.draw and grid.newpage (used for Heatmaps)
  library(grDevices)  # For cairo_pdf and dev.off
})

# ---- 🛠️ 2. Smart Setup (Find & source UTILS.R) ----

get_utils_path <- function() {
  # 1. Windows dev machine
  if (.Platform$OS.type == "windows") {
    return("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/Scripts/nextflow/modules/UTILS.R")
  }
  
  # 2. Interactive Linux / macOS (HPC interactive session)
  if (interactive()) {
    # Assume project root is current working directory
    return(file.path(getwd(), "modules", "UTILS.R"))
  }
  
  # 3. Non-interactive (Nextflow / Rscript)
  initial.options <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", initial.options, value = TRUE)
  if (length(file_arg) == 0) stop("Cannot detect script path for UTILS.R!")
  script_dir <- dirname(sub("--file=", "", file_arg))
  return(file.path(script_dir, "UTILS.R"))
}

utils_path <- get_utils_path()
if (!file.exists(utils_path)) stop(paste("❌ UTILS.R not found at:", utils_path))
source(utils_path)

# ---- 🧬 3. Function Definition (Always Runs) ----

plot_pathways <- function(contrast, pathway_df, output_dir, expr_mat = NULL, 
                          metadata = NULL, col_annotations = NULL, 
                          col_cluster_by = "all", top_n = 10){
  
  # 1. Load Data
  gsea_data <- load_smart(input_path = pathway_df, sheet = "GSEA")
  ora_data  <- load_smart(input_path = pathway_df, sheet = "ORA")
  expr_mat  <- load_smart(input_path = expr_mat)
  metadata  <- load_smart(input_path = metadata)
  
  # Create a named list to track which is which
  analysis_list <- list("GSEA" = gsea_data, "ORA" = ora_data)
  
  for (analysis_name in names(analysis_list)) {
    
    current_pathway_df <- analysis_list[[analysis_name]]
    
    # Skip if the sheet was empty or not found
    if (is.null(current_pathway_df) || nrow(current_pathway_df) == 0) next
    
    # ---- 📊 Define Map Logic ----
    
    # Deduce the enrichment metric based on columns present in current_pathway_df
    metric <- base::intersect(colnames(current_pathway_df), c("NES", "GeneRatio"))[1]
    
    # Safety Check: If neither exists, metric will be NA or NULL
    if (is.na(metric) || is.null(metric)) {
      log_error(sample = "", step = "plot_pathways", 
                msg = "Execution halted : Could not find 'GeneRatio' or 'NES' in the data.")
    }
    
    # Set the x-axis label based on the deduced metric
    x_label <- base::switch(EXPR        = metric,
                            "GeneRatio" = "Gene Ratio",
                            "NES"       = "Normalized Enrichment Score (NES)")
    
    # Set the count column based on the deduced metric
    counts  <- base::switch(EXPR        = metric,
                            "GeneRatio" = "k",
                            "NES"       = "leading_edge_size")
    
    # Set the gene_col based on the deduced metric
    gene_col  <- base::switch(EXPR        = metric,
                              "GeneRatio" = "overlapGenes",
                              "NES"       = "leadingEdge")
    
    sig_col     <- "pval"   # Sometimes, we get 0 pathways with padj <= 0.05
    color_col   <- "Direction"
    plot_colors <- c("Upregulated" = "#E69F00", "Downregulated" = "#56B4E9")
    
    # ---- 🔄 Iterate Through Pathway Collections ----
    
    dot_plots     <- list()
    bar_plots     <- list()
    heatmap_plots <- list()
    
    # 1. Pre-process and Filter Data 
    processed_df <- current_pathway_df %>%
      # Cleaning
      dplyr::mutate(Description = base::gsub("_", " ", Description),
                    Description = stringr::str_wrap(Description, width = 30)) %>%
      # Filtering
      dplyr::filter(!is.na(.data[[metric]]), 
                    .data[[counts]] >= 3, 
                    .data[[sig_col]] <= 0.05) %>%
      # Grouping
      dplyr::group_by(Collection, Direction) %>%
      # SORTING INSIDE THE GROUP (This is the critical step)
      # .by_group = TRUE ensures the sort stays within the Collection/Direction bucket
      dplyr::arrange(dplyr::desc(abs(.data[[metric]])), .by_group = TRUE) %>%
      # SLICING (Take the top N that are now at the top of the group)
      dplyr::slice_head(n = top_n) %>% 
      # RANKING (row_number now strictly counts 1 to 10 within the group)
      dplyr::mutate(rank = dplyr::row_number()) %>% 
      # CLEANUP
      dplyr::ungroup()
    
    # 2. Create the Master Skeleton
    # This creates a row for every possible combination of Collection, Direction, and Rank
    skeleton_df <- tidyr::expand(processed_df, 
                                 Collection, 
                                 Direction = c("Upregulated", "Downregulated"), 
                                 rank = seq_len(top_n))
    
    # 3. Join Actual Data to Skeleton
    master_plot_df <- skeleton_df %>%
      dplyr::left_join(processed_df, by = c("Collection", "Direction", "rank")) %>%
      # Create a 'CleanDescription' that uses unique placeholder for NAs.
      # This prevents ggplot from collapsing all NA rows into a single bar.
      dplyr::mutate(PlotLabel = ifelse(is.na(Description), 
                                       paste0("empty_", Collection, "_", Direction, "_", rank), 
                                       Description)) %>%
      dplyr::arrange(Collection, Direction, dplyr::desc(rank)) %>% # Rank 1 at top
      dplyr::mutate(PlotLabel = factor(PlotLabel, levels = unique(PlotLabel)))
    
    for (collection in unique(master_plot_df$Collection)) {
      
      # Get top pathways for each collection 
      collection_df <- master_plot_df %>% 
        dplyr::filter(Collection == collection)
      
      # 🚨 NEW SAFETY CHECK
      # Since nrow is always 20, we check if all 'Description' values are NA
      if (all(is.na(collection_df$Description))) next
      
      # Dynamic Y-axis text sizing
      max_label_len <- base::max(base::nchar(collection_df$Description), na.rm = TRUE)
      y_text_size <- dplyr::case_when(max_label_len > 50 ~ 6,
                                      max_label_len > 35 ~ 7,
                                      max_label_len > 25 ~ 8,
                                      TRUE ~ 10)
      
      # Calculate limits dynamically per collection
      x_min <- if (metric == "NES") { base::pmin(0, floor(min(collection_df[[metric]], na.rm = TRUE))) } else  { 0 }
      x_limits <- c(x_min, NA)
      
      # ---- 📈 Generate Bar Plot ---- 
      
      bar_p <- ggplot2::ggplot(data = collection_df,
                               mapping = aes(x     = .data[[metric]],
                                             y     = PlotLabel,
                                             fill  = .data[[color_col]],
                                             # Binary alpha for clear communication
                                             alpha = padj <= 0.05)) +
        ggplot2::geom_col(width = 0.75, na.rm = TRUE) +
        ggplot2::labs(x = x_label, 
                      y = "", 
                      title = collection,
                      subtitle = "Solid bars = FDR < 0.05 | Faded bars = nominal p < 0.05",
                      fill = color_col) +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = y_text_size)) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::scale_x_continuous(limits = x_limits, expand = ggplot2::expansion(mult = c(0, 0.05))) +
        ggplot2::scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.35), 
                                    name = "Significance",
                                    labels = c("TRUE" = "FDR <= 0.05", "FALSE" = "Nominal p ≤ 0.05")) +
        # USE YOUR FORMATTED LABELS FOR THE Y-AXIS
        ggplot2::scale_y_discrete(labels = function(x) ifelse(startsWith(x, "empty_"), "", x)) +
        ggplot2::scale_fill_manual(values = plot_colors,
                                   na.translate = FALSE) +  # This hides "NA" from the legend
        guides(fill  = guide_legend(override.aes = list(shape = 22, size = 6)),
               color = guide_legend(override.aes = list(shape = 22, size = 6)),
               alpha = guide_legend(override.aes = list(shape = 22, size = 6))) +
        ggplot2::geom_text(aes(label = .data[[counts]]), x = 0, hjust = -0.5, size = 3, show.legend = FALSE)
      
      bar_plots[[collection]] <- bar_p
      
      # ---- 🟢 Generate Dot Plot ---- 
      
      size_vals <- collection_df[[counts]][!base::is.na(collection_df[[counts]])]
      #breaks <- as.vector(floor(stats::quantile(size_vals, na.rm = TRUE) / 10) * 10)
      
      # Filter: Keep only integers AND stay within the actual data range
      my_breaks <- base::pretty(size_vals, n = 4)
      my_breaks <- my_breaks[my_breaks >= min(size_vals) & 
                             my_breaks <= max(size_vals) & 
                             my_breaks %% 1 == 0]       # Checking x %% 1 == 0 ensures it's a whole number
      
      dot_p <- ggplot2::ggplot(data = collection_df,
                               mapping = aes(x     = .data[[metric]],
                                             y     = PlotLabel,
                                             fill  = .data[[color_col]],
                                             # Binary alpha for clear communication
                                             alpha = padj <= 0.05,
                                             color = .data[[color_col]],
                                             size  = .data[[counts]])) +
        ggplot2::geom_point(na.rm = TRUE) +
        ggplot2::labs(x = x_label, 
                      y = "", 
                      title = collection, 
                      color = color_col, 
                      size = "Counts") +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = y_text_size)) +
        ggplot2::coord_cartesian(clip = "off") + 
        ggplot2::scale_x_continuous(limits = x_limits, expand = ggplot2::expansion(mult = c(0, 0.05))) +
        ggplot2::scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.35), 
                                    name = "Significance",
                                    labels = c("TRUE" = "FDR <= 0.05", "FALSE" = "Nominal p ≤ 0.05")) +
        # USE YOUR FORMATTED LABELS FOR THE Y-AXIS
        ggplot2::scale_y_discrete(labels = function(x) ifelse(startsWith(x, "empty_"), "", x)) +
        ggplot2::scale_color_manual(values = plot_colors,
                                    na.translate = FALSE) +  # This hides "NA" from the legend
        ggplot2::scale_fill_manual(values = plot_colors,     # need for coloring the legend
                                   na.translate = FALSE) +   # This hides "NA" from the legend  
        ggplot2::scale_size(range = c(2, 6), breaks = unique(my_breaks)) 
      
      dot_plots[[collection]] <- dot_p
      
      # ---- 🔥 ️ Generate Heatmap ----
      
      if (is.null(expr_mat)) next
      
      #   gene_pathway_mapping <- collection_df %>%
      #     dplyr::filter(!is.na(.data[[gene_col]])) %>%
      #     dplyr::select(Description, .data[[gene_col]], .data[[metric]]) %>%
      #     tidyr::separate_rows(.data[[gene_col]], sep = "/") %>%
      #     dplyr::rename(GeneID = .data[[gene_col]]) %>%
      #     dplyr::mutate(GeneID = base::trimws(GeneID)) %>%
      #     # Keep only genes present in your expression matrix
      #     dplyr::filter(GeneID %in% base::rownames(expr_mat))
      
      #   metadata_row <- gene_pathway_mapping %>%
      #     dplyr::group_by(GeneID) %>%
      #     dplyr::slice_max(order_by = abs(.data[[metric]]), n = 1, with_ties = FALSE) %>%
      #     dplyr::ungroup() %>%
      #     dplyr::select(GeneID, Description) %>%
      #     dplyr::rename(Pathways = Description, SYMBOL = GeneID)
      
      metadata_col <- NULL
      if (!is.null(col_annotations)) {
        # Only select if col_annotations is not NULL
        metadata_col <- metadata %>% 
          dplyr::select(dplyr::any_of(col_annotations), Sample_ID)
      }
      
      for (pathway in unique(collection_df$Description)) {
        
        # Extract genes for this specific pathway from the long-format dataframe
        plot_genes <- collection_df %>%
          dplyr::filter(Description == pathway) %>%
          dplyr::pull(.data[[gene_col]]) %>%
          # Split the string by the slash into a vector
          base::strsplit(split = "/") %>%
          base::unlist() %>%
          # Clean up whitespace and remove empty/NA values
          base::trimws() %>%
          stats::na.omit() %>%
          .[. != ""] %>%
          # Keep only unique genes that actually exist in your expression matrix
          base::unique() %>%
          base::intersect(rownames(expr_mat))
        
        # Skip plotting if less than 2 genes
        if (base::length(plot_genes) < 2) next
        
        # Plot heatmap
        ph_output <- plot_heatmap(expr_mat            = expr_mat[plot_genes, ,drop = FALSE],
                                  data_type           = "counts_log",
                                  label_genes         = NULL,
                                  metadata_col        = metadata_col,
                                  metadata_row        = NULL,
                                  col_annotations     = col_annotations,
                                  row_annotations     = NULL,
                                  col_gap_by          = NULL,
                                  row_gap_by          = NULL,
                                  col_cluster_by      = col_cluster_by,
                                  row_cluster_by      = "all",
                                  plot_title          = stringr::str_wrap(string = pathway, width = 30),
                                  heatmap_palette     = "rdbu",
                                  border_color        = NA,
                                  show_expr_legend    = TRUE)
        if (!is.null(ph_output)) {
          heatmap_plots[[paste0(collection, "_", pathway)]] <- ph_output$ph$gtable
        }
      }
    } 
    
    # ---- 💾 Save Consolidated Summary Plots ----
    
    summary_plots <- base::list(Bar = bar_plots, Dot = dot_plots, Heatmap = heatmap_plots)
    
    for (type in names(summary_plots)) {
        
        # Skip if there is nothing to plot
        plot_list <- summary_plots[[type]]
        if (length(plot_list) == 0) next
        
        # Setup Paths
        safe_contrast <- gsub("[^[:alnum:]_-]", "_", contrast)
        output_path   <- file.path(output_dir, safe_contrast)
        if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
        output_file   <- file.path(output_path, paste0(type, "_Plot_Pathways.pdf"))
        
        # Save by Type
        if (type %in% c("Bar", "Dot")) {
          
          # Calculate rows (3 plots per row)
          n_rows <- ceiling(length(plot_list) / 3)
          
          ggplot2::ggsave(filename = output_file,
                          plot     = cowplot::plot_grid(plotlist = plot_list, ncol = 3, align = "hv"),
                          device   = grDevices::cairo_pdf,
                          width    = 3 * 6, 
                          height   = n_rows * 6,
                          bg       = "white")
        } else {
          # Heatmaps: Multi-page PDF
          grDevices::cairo_pdf(output_file, width = 8, height = 11.5, onefile = TRUE)
          for (ht in plot_list) {
            grid::grid.newpage()
            grid::grid.draw(ht)
          }
          grDevices::dev.off()
        }
    }
  }
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
  
  # --- PROCESS BATCH VARS SAFELY ---
  raw_col_ann <- get_arg(6)
  
  # split by comma, then trim leading/trailing whitespace from each element
  col_ann_list <- if (!is.null(raw_col_ann)) {
    trimws(strsplit(raw_col_ann, ",")[[1]]) 
  } else {
    NULL
  }
  
  plot_pathways(
    contrast        = get_arg(1),
    pathway_df      = get_arg(2),
    output_dir      = get_arg(3, "."),
    expr_mat        = get_arg(4),              
    metadata        = get_arg(5), 
    col_annotations = col_ann_list,
    col_cluster_by  = get_arg(7, "all"),
    top_n           = as.numeric(get_arg(8, 10))
    )
  
}