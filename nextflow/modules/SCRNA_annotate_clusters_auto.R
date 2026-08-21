# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)
  library(glue)
  library(openxlsx)
  library(Seurat)
  library(SeuratObject)
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

modules_dir <- get_modules_dir()
utils_path          <- file.path(modules_dir, "UTILS.R")
if (!file.exists(utils_path)) stop(paste("❌ UTILS.R not found at:", modules_dir))
source(utils_path)

# ---- 🧬 3. Function Definition (Always Runs) ----

annotate_clusters_auto <- function(rds_file,
                               output_dir          = NULL,
                               AVG_PROB_THRESHOLD  = 0.3,
                               DELTA_THRESHOLD     = 0.1,
                               STABILITY_THRESHOLD = 0.9) {

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Load Seurat Object
  # ═══════════════════════════════════════════════════════════════════════════
  # Expects scored_seurat.rds produced by CALC_MODULE_SCORES.
  # That object contains:
  #   - SCT assay with normalized expression
  #   - Seurat module score columns in metadata (e.g. T.Cell1, B.Cell2...)
  #   - UCell module score columns in metadata (e.g. T.Cell_UCell...)
  #   - Cluster columns in metadata (harmony_res0.2 ... harmony_res1.2)
  #   - @misc$scoring_params with exact column names and signatures
  #   - @misc$clustering_params and @misc$integration_params from upstream

  integrated_seurat <- load_smart(input_path = rds_file)

  log_info(sample = "", step = "annotate_clusters_auto",
           msg = glue("Loaded: {ncol(integrated_seurat)} cells, ",
                      "{nrow(integrated_seurat)} features."))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Resolve Parameters from @misc
  # ═══════════════════════════════════════════════════════════════════════════
  # Read exact column names written by CALC_MODULE_SCORES.
  # Using exact stored names rather than regex pattern matching avoids the
  # fragility of detecting columns by suffix — regex on "[0-9]+$" also matches
  # cluster columns (harmony_res0.8) requiring hacky prefix exclusion.
  #
  # cluster_cols       : e.g. c("harmony_res0.2", "harmony_res0.4"... "harmony_res1.2")
  #                      All clustering resolutions for Harmony —
  #                      annotation runs across ALL resolutions equally.
  #                      No single resolution is trusted over others;
  #                      the multi-resolution consensus (mode vote) determines
  #                      final assignment. This avoids brittleness of any
  #                      single resolution choice.
  #
  # seurat_module_cols : e.g. c("T.Cell1", "B.Cell2", "Epithelial3"...)
  # ucell_module_cols  : e.g. c("T.Cell_UCell", "B.Cell_UCell"...)

  # Read exact column names written by upstream modules.
  # cluster_cols: written by CLUSTER_AND_FINDMARKERS → @misc$clustering_params$cluster_cols
  # score cols:  written by CALC_MODULE_SCORES → @misc$scoring_params
  # SC3 optimal_resolutions no longer exists — resolution selection is handled
  # by the multi-resolution consensus in Section 5/6, not by a single chosen
  # resolution. Each resolution is an equal vote in the consensus.

  seurat_module_cols <- integrated_seurat@misc$scoring_params$seurat_module_cols %||% {
    log_error(sample = "", step = "annotate_clusters_auto",
              msg = "seurat_module_cols missing from @misc. Run CALC_MODULE_SCORES first.")
  }

  ucell_module_cols <- integrated_seurat@misc$scoring_params$ucell_module_cols %||% {
    log_error(sample = "", step = "annotate_clusters_auto",
              msg = "ucell_module_cols missing from @misc. Run CALC_MODULE_SCORES first.")
  }

  # Read cluster cols from CLUSTER_AND_FINDMARKERS @misc first,
  # fall back to CALC_MODULE_SCORES @misc, then derive from scratch.
  cluster_cols <- integrated_seurat@misc$clustering_params$cluster_cols %||%
    integrated_seurat@misc$scoring_params$cluster_cols %||% {
      resolutions <- integrated_seurat@misc$clustering_params$resolutions %||% c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
      as.vector(outer("harmony", resolutions, function(m, r) paste0(m, "_res", r)))
    }
  cluster_cols <- cluster_cols[cluster_cols %in% colnames(integrated_seurat@meta.data)]

  if (length(cluster_cols) == 0) {
    log_error(sample = "", step = "annotate_clusters_auto",
              msg = "No cluster columns found in metadata. Run CLUSTER_AND_FINDMARKERS first.")
  }

  log_info(sample = "", step = "annotate_clusters_auto",
           msg = glue("Found {length(cluster_cols)} cluster columns, ",
                      "{length(seurat_module_cols)} Seurat score columns, ",
                      "{length(ucell_module_cols)} UCell score columns."))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Z-score Transformation of Module Scores
  # ═══════════════════════════════════════════════════════════════════════════
  # WHY NOT STORED FROM CALC_MODULE_SCORES:
  # Z-score computation on a 100k × 15 matrix takes ~2 seconds — negligible
  # cost that doesn't justify storing 4 large matrices in @misc permanently.
  # Raw scores are already in metadata; zscores are recomputed here on the fly.
  #
  # STEP 1 — Clamp negative Seurat scores to zero:
  #   Seurat AddModuleScore can produce negative values when a gene set is
  #   expressed below background. In single cell data, absence of expression
  #   is unreliable due to dropout — a cell may fail to detect a gene not
  #   because it's truly unexpressed but because of capture failure.
  #   Analogous to flow cytometry: we annotate on positive expression,
  #   not on absence. So negative = "no confident signal" = 0.
  #   UCell scores are always >= 0 by construction (clamping is redundant
  #   but applied for consistency).
  #
  # STEP 2 — Row-wise Z-score standardisation:
  #   For each cell, standardise scores ACROSS cell types (not across cells).
  #   This answers: "relative to all other cell type scores for THIS cell,
  #   how does each score compare?"
  #   Achieved by transposing (so cell types become rows), scaling (Z-score),
  #   then transposing back.
  #   Result: each cell has a Z-score per cell type with mean=0, SD=1 across
  #   its own scores.
  #
  # STEP 3 — Safe scaling for zero-variance rows:
  #   If a cell has identical scores across ALL cell types (e.g. all zero
  #   after clamping — likely ambient RNA or very low quality cell that
  #   slipped QC), standard scale() produces NaN because SD=0.
  #   safe_scale() detects this and returns zeros instead, so these cells
  #   get a flat probability distribution after softmax → correctly receive
  #   low confidence scores → filtered by stability threshold.

  safe_scale <- function(x) {
    # If all values identical (including all-zero), SD=0 → return zeros
    # rather than NaN. These cells will get uniform softmax distribution
    # → low confidence → filtered downstream by STABILITY_THRESHOLD.
    if (stats::sd(x, na.rm = TRUE) == 0) {
      return(rep(0, length(x)))
    } else {
      return(as.numeric(scale(x)))
    }
  }

  # Extract raw score matrices from metadata
  seurat_score_matrix <- integrated_seurat@meta.data[, seurat_module_cols, drop = FALSE] %>%
    as.matrix()
  ucell_score_matrix  <- integrated_seurat@meta.data[, ucell_module_cols,  drop = FALSE] %>%
    as.matrix()

  # Clamp negatives to zero — treat as no signal (see rationale above)
  seurat_score_matrix[seurat_score_matrix < 0] <- 0
  ucell_score_matrix[ucell_score_matrix   < 0] <- 0    # redundant for UCell, applied for consistency

  # Row-wise Z-score — standardise across cell types per cell
  seurat_zscore_matrix <- t(apply(seurat_score_matrix, 1, safe_scale))
  ucell_zscore_matrix  <- t(apply(ucell_score_matrix,  1, safe_scale))

  # Restore column names lost by apply()
  colnames(seurat_zscore_matrix) <- seurat_module_cols
  colnames(ucell_zscore_matrix)  <- ucell_module_cols

  log_info(sample = "", step = "annotate_clusters_auto",
           msg = "Z-score transformation completed.")

  # Free raw score matrices immediately — no longer needed
  rm(seurat_score_matrix, ucell_score_matrix); gc()

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Softmax Transformation — Z-scores to Pseudo-Probabilities
  # ═══════════════════════════════════════════════════════════════════════════
  # Softmax converts Z-scores into a distribution summing to 1, making
  # values comparable across cells and interpretable as relative confidence.
  #
  # WHY SOFTMAX IS APPROPRIATE HERE:
  #   - Marker gene sets are broad and orthogonal (T cell vs Epithelial vs
  #     Erythrocyte) — genuine cell types produce one dominant high Z-score
  #     with all others near zero. Softmax correctly amplifies this contrast.
  #   - For fine-grained subtypes (CD4 vs CD8) softmax would be overconfident
  #     on small differences, but we are NOT doing subtype annotation here.
  #   - Cluster-level averaging (Section 5) absorbs per-cell softmax noise —
  #     individual cells with ambiguous distributions are outvoted by the
  #     cluster majority.
  #
  # NUMERICAL STABILITY:
  #   Subtracting max(x) before exp() prevents floating point overflow
  #   when Z-scores are large. Does not change the output distribution.
  #
  # NOTE: Softmax always produces a winner even for low-quality cells.
  # This is handled downstream by three cluster-level filters (Section 6):
  #   AVG_PROB_THRESHOLD, delta, and entropy — not at the per-cell level.

  softmax <- function(x) {
    e_x <- exp(x - max(x, na.rm = TRUE))   # subtract max for numerical stability
    e_x / sum(e_x, na.rm = TRUE)
  }

  # Apply row-wise softmax — each cell gets a probability distribution
  seurat_prob_matrix <- t(apply(seurat_zscore_matrix, 1, softmax))
  ucell_prob_matrix  <- t(apply(ucell_zscore_matrix,  1, softmax))

  # Restore column names
  colnames(seurat_prob_matrix) <- seurat_module_cols
  colnames(ucell_prob_matrix)  <- ucell_module_cols

  # Handle any residual NaN (should be caught by safe_scale but defensive check)
  seurat_prob_matrix[is.nan(seurat_prob_matrix)] <- 0
  ucell_prob_matrix[is.nan(ucell_prob_matrix)]   <- 0

  log_info(sample = "", step = "annotate_clusters_auto",
           msg = "Softmax transformation completed.")

  # Free zscore matrices immediately — no longer needed
  rm(seurat_zscore_matrix, ucell_zscore_matrix); gc()

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5: Weighted Cluster Consensus — Core Annotation Logic
  # ═══════════════════════════════════════════════════════════════════════════
  # For each scoring method (Seurat, UCell) × each clustering resolution,
  # compute cluster-level average probabilities and assign cell types.
  #
  # WHY CLUSTER-LEVEL NOT PER-CELL:
  #   Single cell expression is noisy — dropout, ambient RNA, and capture
  #   variability mean per-cell scores are unreliable for individual cells.
  #   Clusters share transcriptional programs; averaging within a cluster
  #   smooths noise and gives a robust estimate of the dominant cell type.
  #
  # TRIMMED MEAN (trim = 0.1):
  #   Drops bottom 10% and top 10% of cells by probability before averaging.
  #   Rationale: every cluster contains a small fraction of outlier cells —
  #   boundary cells, ambient RNA contaminated cells, or rare doublets that
  #   passed QC. These outliers pull the mean toward zero or toward a wrong
  #   cell type. Trimmed mean is robust against these without discarding data.
  #   trim=0.1 means for a cluster of 100 cells, 10 cells dropped from each
  #   tail — aggressive enough to remove outliers, conservative enough to
  #   retain cluster biology.
  #
  # THREE CONFIDENCE FILTERS applied to cluster average probabilities:
  #
  #   1. AVG_PROB_THRESHOLD (default 0.3):
  #      Is the winning cell type's average probability strong enough?
  #      With 15 cell types the uniform (random) baseline is 1/15 = 0.067.
  #      A threshold of 0.3 requires the winner to be ~4.5x above baseline.
  #      Clusters below this have no meaningful dominant signal.
  #
  #   2. DELTA_THRESHOLD (default 0.1):
  #      Is the winner clearly ahead of the runner-up?
  #      delta = prob(rank1) - prob(rank2)
  #      Catches cases where two cell types are nearly tied — e.g.
  #      T.Cell=0.35, B.Cell=0.34 both pass AVG_PROB_THRESHOLD but
  #      delta=0.01 reveals the assignment is unreliable.
  #      A threshold of 0.1 requires the winner to lead by 10 probability
  #      points — meaningful separation for broad cell types.
  #
  #   3. ENTROPY (stored as metadata, not a hard filter):
  #      Shannon entropy measures how spread the probability distribution is.
  #      Low entropy = confident (probability concentrated on one cell type).
  #      High entropy = uncertain (probability spread across many types).
  #      Normalised by log(n_celltypes) so range is always [0, 1].
  #      Stored per cluster and mapped to cells for downstream QC —
  #      lets researchers apply their own threshold rather than hardcoding one.
  #      High entropy clusters that pass AVG_PROB and delta filters likely
  #      represent cell types not in the marker list or genuinely mixed clusters.
  #
  # TIE HANDLING:
  #   If two cell types share exactly equal maximum probability after
  #   AVG_PROB and delta filters, the cluster is assigned "Ambiguous"
  #   rather than pasting both names (e.g. "Tcell,Bcell").
  #   Pasted compound labels would be treated as unique classes by the
  #   mode calculation in Section 6, destabilising the consensus.
  #
  # LABEL CLEANING:
  #   AddModuleScore appends numeric suffix (T.Cell1 → T.Cell).
  #   UCell appends _UCell (T.Cell_UCell → T.Cell).
  #   Both suffixes are stripped so Seurat and UCell assignments for the
  #   same cell type produce identical labels — essential for the cross-method
  #   mode calculation in Section 6.

  score_cols_list  <- list(Seurat = seurat_module_cols,
                           UCell  = ucell_module_cols)
  prob_matrix_list <- list(Seurat = seurat_prob_matrix,
                           UCell  = ucell_prob_matrix)

  # Initialise list to store per-method per-resolution cluster summaries
  cluster_summary_list <- list()

  for (method in c("Seurat", "UCell")) {

    score_cols  <- score_cols_list[[method]]
    prob_matrix <- prob_matrix_list[[method]]

    # Build working dataframe: probability columns + all cluster columns
    temp_df <- integrated_seurat@meta.data[, cluster_cols, drop = FALSE]
    temp_df[, score_cols] <- prob_matrix

    for (col in cluster_cols) {

      # ---- Compute trimmed mean probabilities per cluster ----
      # trim=0.1 drops top and bottom 10% of cells per cluster per cell type
      # before computing mean — robust against outlier cells in each cluster
      avg_probs <- temp_df %>%
        dplyr::group_by(.data[[col]]) %>%
        dplyr::summarise(
          dplyr::across(
            .cols = dplyr::all_of(score_cols),
            .fns  = ~ mean(.x, trim = 0.1, na.rm = TRUE)   # trimmed mean
          ),
          .groups = "drop"
        )

      # ---- Compute entropy per cluster ----
      # Shannon entropy normalised by log(n_celltypes) → range [0, 1]
      # 0 = perfectly certain (one cell type dominates)
      # 1 = completely uncertain (all cell types equally probable)
      # Stored as metadata — not used as hard filter — lets researcher
      # identify clusters where no cell type in the marker list fits well
      n_types <- length(score_cols)
      entropy_vals <- apply(avg_probs[, score_cols], 1, function(p) {
        p <- p[p > 0]                              # ignore zero probabilities
        H <- -sum(p * log(p))                      # Shannon entropy
        H / log(n_types)                           # normalise to [0,1]
      })
      avg_probs$Annotation_Entropy <- round(entropy_vals, digits = 3)

      # ---- Apply confidence filters ----
      # Filter 1: AVG_PROB_THRESHOLD — absolute strength of winner
      # Filter 2: DELTA_THRESHOLD — relative strength vs runner-up
      # Both must pass for a confident assignment; otherwise → Unknown

      avg_probs$Assigned_Type <- apply(avg_probs[, score_cols], 1, function(row) {

        max_prob <- max(row, na.rm = TRUE)

        # Filter 1: winner must exceed absolute probability threshold
        if (max_prob < AVG_PROB_THRESHOLD) return("Unknown")

        # Filter 2: winner must be clearly ahead of runner-up
        sorted_probs <- sort(row, decreasing = TRUE)
        delta        <- sorted_probs[1] - sorted_probs[2]
        if (delta < DELTA_THRESHOLD) return("Unknown")

        # Tie handling: if two cell types share exactly equal maximum
        # probability, assign Ambiguous rather than pasting both names
        # which would create spurious unique labels in mode calculation
        winning_types <- names(row)[row == max_prob]
        if (length(winning_types) > 1) return("Ambiguous")

        return(winning_types)
      })

      # Store max probability for each cluster (confidence indicator)
      max_vals       <- do.call(pmax, c(avg_probs[, score_cols], na.rm = TRUE))
      avg_probs$Max_P <- round(max_vals, digits = 3)

      # ---- Clean suffixes from assigned labels ----
      # Remove numeric suffix added by AddModuleScore: "T.Cell1" → "T.Cell"
      # Remove _UCell suffix added by UCell: "T.Cell_UCell" → "T.Cell"
      # This ensures Seurat and UCell produce identical label strings for
      # the same cell type — critical for cross-method mode calculation
      cleaned_types <- gsub(pattern = "_UCell$", replacement = "", x = avg_probs$Assigned_Type)
      cleaned_types <- gsub(pattern = "[0-9]+$",  replacement = "", x = cleaned_types)
      avg_probs$Assigned_Type <- cleaned_types

      # ---- Map cluster assignments back to individual cells ----
      # Convert to named vectors: names = cluster IDs, values = assignment
      # Lookup by cell's cluster ID returns that cluster's assignment
      assigned_vec  <- stats::setNames(avg_probs$Assigned_Type,   avg_probs[[col]])
      meanprob_vec  <- stats::setNames(avg_probs$Max_P,            avg_probs[[col]])
      entropy_vec   <- stats::setNames(avg_probs$Annotation_Entropy, avg_probs[[col]])

      cell_clusters       <- integrated_seurat@meta.data[[col]]
      cluster_assignment  <- assigned_vec[as.character(cell_clusters)]
      cluster_meanprob    <- meanprob_vec[as.character(cell_clusters)]
      cluster_entropy     <- entropy_vec[as.character(cell_clusters)]

      # Restore barcodes as names (lost during vector lookup)
      cluster_assignment <- stats::setNames(cluster_assignment, rownames(integrated_seurat@meta.data))
      cluster_meanprob   <- stats::setNames(cluster_meanprob,   rownames(integrated_seurat@meta.data))
      cluster_entropy    <- stats::setNames(cluster_entropy,    rownames(integrated_seurat@meta.data))

      # Store in Seurat metadata
      # Column naming: "{Method}_{resolution}" e.g. "Seurat_rpca0.4", "UCell_harmony0.8"
      integrated_seurat[[paste(method, col, sep = "_")]]              <- cluster_assignment
      integrated_seurat[[paste(method, col, "MeanProb", sep = "_")]]  <- cluster_meanprob
      integrated_seurat[[paste(method, col, "Entropy",  sep = "_")]]  <- cluster_entropy

      # Store cluster-level summary for Excel export
      cluster_summary_list[[paste(method, col, sep = "_")]] <- avg_probs %>%
        dplyr::select(dplyr::any_of(col), Max_P, Annotation_Entropy, Assigned_Type)
    }
  }

  # Save cluster-level summary to Excel for manual review
  # One sheet with all methods and resolutions stacked — includes
  # cluster ID, winning probability, entropy, and assigned type
  if (!is.null(output_dir)) {
    cluster_summary_df <- dplyr::bind_rows(cluster_summary_list, .id = "Method_Resolution")
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Cluster_Summary")
    openxlsx::writeData(wb, "Cluster_Summary", cluster_summary_df)
    openxlsx::saveWorkbook(wb,
                           file      = file.path(output_dir, "Cluster_Annotation_Summary.xlsx"),
                           overwrite = TRUE)
    log_info(sample = "", step = "annotate_clusters_auto",
             msg = "Cluster annotation summary saved to Cluster_Annotation_Summary.xlsx")
  }

  log_info(sample = "", step = "annotate_clusters_auto",
           msg = glue("Cluster consensus completed across ",
                      "{length(cluster_cols)} resolutions × 2 methods ",
                      "= {length(cluster_cols) * 2} annotation columns."))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 6: Stability Score and Final CellType Assignment
  # ═══════════════════════════════════════════════════════════════════════════
  # After Section 5, each cell has multiple annotations — one per
  # method × resolution combination (e.g. 2 methods × 6 resolutions = 12).
  # These are now consolidated into a single final CellType and a
  # stability score reflecting annotation confidence.
  #
  # FINAL CELLTYPE — mode across all method × resolution combinations:
  #   The most frequently assigned label across all 12 combinations is
  #   taken as the final CellType. This ensemble approach is analogous to
  #   random forest majority voting — each resolution is a different "view"
  #   of the data and the mode is robust against any single view being wrong.
  #
  #   tabulate() + match() used for speed — much faster than table() on
  #   100k cells across many combinations. Critical for large datasets.
  #
  # STABILITY SCORE — fraction of combinations agreeing with the mode:
  #   stability = (n combinations matching mode) / (total combinations)
  #   Range: [0, 1]. Higher = more consistent annotation across methods
  #   and resolutions.
  #
  #   Interpretation:
  #   1.0  = all 12 combinations agree — highly confident
  #   0.9  = 11/12 agree — confident (default clean subset threshold)
  #   0.5  = 6/12 agree — uncertain, likely transitional or boundary cell
  #   <0.5 = majority disagree — unclassifiable with current markers
  #
  #   NOTE: STABILITY_THRESHOLD (default 0.9) is used to create a clean
  #   subset for downstream DE analysis where pure populations are needed.
  #   ALL cells are retained in the full annotated object with stability
  #   score as metadata — researchers can apply their own threshold for
  #   specific analyses (e.g. trajectory analysis may accept lower stability).
  #
  # WHY NOT WEIGHT ANY SINGLE RESOLUTION:
  #   All resolutions contribute equally to the consensus vote. Weighting
  #   one resolution higher would reintroduce single-resolution brittleness
  #   that multi-resolution consensus was designed to avoid.

  # Identify all annotation columns produced in Section 5
  # Pattern: starts with "Seurat_" or "UCell_" followed by a cluster column name
  consensus_cols <- colnames(integrated_seurat@meta.data)
  consensus_cols <- consensus_cols[
    grepl(
      pattern = paste0("^(Seurat|UCell)_(", paste(cluster_cols, collapse = "|"), ")$"),
      x       = consensus_cols
    )
  ]

  consensus_mat <- integrated_seurat@meta.data[, consensus_cols, drop = FALSE]

  # Final CellType — mode across all method × resolution combinations
  # tabulate(match()) is fast for large vectors — avoids slow table() calls
  cell_mode_assignment <- apply(consensus_mat, 1, function(x) {
    ux     <- stats::na.omit(x)
    if (length(ux) == 0) return(NA)
    ux_tab <- tabulate(match(ux, unique(ux)))
    unique(ux)[which.max(ux_tab)]
  })

  # Clean suffixes from mode labels
  # (should already be clean from Section 5 but applied defensively)
  cell_mode_assignment <- gsub("_UCell$", "", cell_mode_assignment)
  cell_mode_assignment <- gsub("[0-9]+$",  "", cell_mode_assignment)
  integrated_seurat$CellType <- cell_mode_assignment

  # Stability Score — fraction of combinations matching the mode
  # Mode matrix: replicate the mode assignment across all columns for comparison
  mode_matrix <- matrix(
    data   = cell_mode_assignment,
    nrow   = nrow(consensus_mat),
    ncol   = ncol(consensus_mat),
    byrow  = FALSE
  )

  stability_score <- rowSums(consensus_mat == mode_matrix, na.rm = TRUE) / ncol(consensus_mat)
  integrated_seurat$Stability_Score <- round(stability_score, digits = 3)

  log_info(sample = "", step = "annotate_clusters_auto",
           msg = glue("CellType and Stability_Score computed. ",
                      "Median stability: {round(median(stability_score), 3)}. ",
                      "Cells >= {STABILITY_THRESHOLD} stability: ",
                      "{sum(stability_score >= STABILITY_THRESHOLD)} ",
                      "({round(mean(stability_score >= STABILITY_THRESHOLD)*100, 1)}%)."))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 7: Store Annotation Parameters in @misc
  # ═══════════════════════════════════════════════════════════════════════════
  # Documents all parameters and outputs for reproducibility.
  # Downstream functions (marker finding, compositional analysis) can
  # read celltype_levels from here rather than re-deriving from metadata.

  integrated_seurat@misc$annotation_params <- list(
    AVG_PROB_THRESHOLD  = AVG_PROB_THRESHOLD,
    DELTA_THRESHOLD     = DELTA_THRESHOLD,
    STABILITY_THRESHOLD = STABILITY_THRESHOLD,
    consensus_cols      = consensus_cols,
    celltype_levels     = sort(unique(cell_mode_assignment[!is.na(cell_mode_assignment)])),
    n_unknown           = sum(cell_mode_assignment == "Unknown",   na.rm = TRUE),
    n_ambiguous         = sum(cell_mode_assignment == "Ambiguous", na.rm = TRUE)
  )

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 8: UMAP Plots
  # ═══════════════════════════════════════════════════════════════════════════
  # Three sets of plots saved:
  #   1. CellType and Stability_Score — main annotation output
  #   2. Stability_Score split by CellType — reveals which cell types are
  #      consistently annotated vs ambiguous across methods and resolutions
  #   3. One PDF per resolution — clusters alongside CellType for visual
  #      validation that cluster boundaries align with cell type boundaries.
  #      Separate PDFs are easier to open and compare than one large file.

  if (!is.null(output_dir)) {

    # Main annotation — CellType and Stability_Score
    plot_seurat(integrated_seurat,
                features     = c("CellType", "Stability_Score"),
                feature_type = "metadata",
                reduction    = "umap_harmony",
                filename     = "UMAP_Annotation_CellType",
                output_dir   = output_dir)

    # Stability split by cell type — reveals which cell types are consistently
    # annotated vs ambiguous across methods and resolutions
    plot_seurat(integrated_seurat,
                features     = "Stability_Score",
                feature_type = "metadata",
                reduction    = "umap_harmony",
                filename     = "UMAP_Stability_ByCellType",
                split_col    = "CellType",
                output_dir   = output_dir)

    # One PDF per resolution — clusters alongside CellType for visual validation
    # Useful for checking whether cluster boundaries align with cell type boundaries
    for (col in cluster_cols) {
      plot_seurat(integrated_seurat,
                  features     = c(col, "CellType"),
                  feature_type = "metadata",
                  reduction    = "umap_harmony",
                  filename     = glue("UMAP_Annotation_{col}"),
                  output_dir   = output_dir)
    }

    log_info(sample = "", step = "annotate_clusters_auto",
             msg = glue("Annotation UMAP plots saved — ",
                        "CellType + Stability + {length(cluster_cols)} resolution plots."))
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 9: Save Seurat Objects
  # ═══════════════════════════════════════════════════════════════════════════
  # Two objects saved:
  #
  #   annotated_seurat_auto.rds — full object, all cells, all metadata
  #     Used for: visualisation, trajectory analysis, compositional analysis,
  #     any analysis where you want to retain transitional/uncertain cells
  #     with their stability scores as confidence metadata.
  #
  #   annotated_seurat_auto_clean.rds — high confidence cells only
  #     Filter: Stability_Score >= STABILITY_THRESHOLD (default 0.9)
  #     Used for: differential expression, marker finding, analyses that
  #     require pure, confidently annotated populations.
  #     NOTE: 0.9 means 90% of method × resolution combinations agreed —
  #     strict but appropriate for DE where contaminating cells inflate
  #     false positives. Make STABILITY_THRESHOLD a parameter so researchers
  #     can relax this for tissues with inherently lower stability
  #     (e.g. tumour microenvironment, developmental datasets).

  if (!is.null(output_dir)) {

    # Full annotated object
    rds_full  <- file.path(output_dir, "annotated_seurat_auto.rds")
    saveRDS(object = integrated_seurat, file = rds_full)
    log_info(sample = "", step = "annotate_clusters_auto",
             msg = glue("Full annotated object saved: '{rds_full}' ",
                        "({ncol(integrated_seurat)} cells)."))

    # Clean high-confidence subset
    clean_seurat <- base::subset(integrated_seurat,
                                 subset = Stability_Score >= STABILITY_THRESHOLD)
    rds_clean <- file.path(output_dir, "annotated_seurat_auto_clean.rds")
    saveRDS(object = clean_seurat, file = rds_clean)
    log_info(sample = "", step = "annotate_clusters_auto",
             msg = glue("Clean annotated object saved: '{rds_clean}' ",
                        "({ncol(clean_seurat)} cells, ",
                        "stability >= {STABILITY_THRESHOLD})."))

  } else {
    log_info(sample = "", step = "annotate_clusters_auto",
             msg = "No output_dir provided — skipping save.")
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 10: Log Summary and Return
  # ═══════════════════════════════════════════════════════════════════════════

  celltype_table <- sort(table(integrated_seurat$CellType), decreasing = TRUE)
  log_info(sample = "", step = "annotate_clusters_auto",
           msg = "CellType distribution:")
  for (ct in names(celltype_table)) {
    pct <- round(celltype_table[ct] / ncol(integrated_seurat) * 100, 1)
    log_info(sample = "", step = "annotate_clusters_auto",
             msg = glue("  {ct}: {celltype_table[ct]} cells ({pct}%)"))
  }

  log_info(sample = "", step = "annotate_clusters_auto",
           msg = glue("Annotation completed. ",
                      "New metadata: CellType, Stability_Score, ",
                      "{length(consensus_cols)} consensus columns, ",
                      "Annotation_Entropy per resolution."))

  return(invisible(integrated_seurat))
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

  # AVG_PROB_THRESHOLD and DELTA_THRESHOLD passed as numeric from Nextflow
  annotate_clusters_auto(
    rds_file            = get_arg(1),
    output_dir          = get_arg(2, default = "."),
    AVG_PROB_THRESHOLD  = as.numeric(get_arg(3, default = "0.3")),
    DELTA_THRESHOLD     = as.numeric(get_arg(4, default = "0.1")),
    STABILITY_THRESHOLD = as.numeric(get_arg(5, default = "0.9"))
  )
}
