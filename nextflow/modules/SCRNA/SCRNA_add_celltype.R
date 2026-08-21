# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)
  library(glue)
  library(ggplot2)
  library(Seurat)
  library(SeuratObject)
  library(openxlsx)
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


# FUNCTION: add_celltype
# Purpose: Assembles the FINAL CellType for every cell in the dataset, by
#          combining Whole's coarse Lineage call with each subcluster branch's
#          own fine-grained re-annotation, then produces the final clean
#          object, final marker table, a per-CellType marker review xlsx, and
#          final UMAP for the whole project.
#
# Design (per project discussion — NOT the old majority-vote/Ambiguous
# design):
#   - Whole's Lineage (coarse: Epithelial/Mesenchymal/Hematopoietic/...) comes
#     in already named "Lineage" (SCRNA_ANNOTATE_CLUSTERS writes it directly
#     to that column for the "Whole" branch — no rename happens here any
#     more) — kept as a permanent reference column, but no longer
#     authoritative on its own.
#   - Each subset branch's own CellType (fine-grained, e.g. "Goblet"/"Basal",
#     from that branch's own SCRNA_ANNOTATE_CLUSTERS run — including that
#     branch's own "Exclude" calls) OVERWRITES CellType for its cells.
#   - Cells never covered by any subset branch this run (lineage not yet
#     subclustered, or Lineage was already "Exclude" and so never entered
#     any branch) simply keep Lineage as their final CellType.
#   - NO majority-vote/consensus check against the original Lineage call —
#     a subset branch's marker-informed reclassification is trusted outright,
#     not flagged as "Ambiguous" for disagreeing with the coarser Whole-level
#     clustering (that coarse call has less information, not more — see
#     project discussion re: cells landing in the wrong Whole-level cluster
#     despite having clear real markers for a different, correctly-identified
#     compartment at the subset level).
#   - CellType_Source tracks provenance per cell: "Lineage" (never
#     subclustered / fell back) or the subset branch label that provided
#     the final call — informational only, for debugging/transparency.
#   - Clean object (Exclude removed, by whatever the FINAL CellType is) is
#     computed HERE for the first time in the whole pipeline — CellType
#     isn't genuinely final until this point.
#   - On the clean object: NO reclustering (CellType is final, a fresh
#     Seurat cluster ID would be redundant). Just:
#       - FindAllMarkers(group.by = "CellType") once — the one authoritative
#         final marker table for the whole project.
#       - A per-CellType top-markers review xlsx built straight off that same
#         result (no master marker file lookup, no "new vs known type" check
#         — every final CellType gets a top-10-marker column to review by
#         eye; genes get added to the master marker file by hand afterward).
#         This replaced the old per-branch "candidate_markers_review.xlsx"
#         that used to live in SCRNA_ANNOTATE_CLUSTERS.
#       - RunUMAP-only recompute (no FindNeighbors — nothing downstream
#         consumes a neighbor graph here) on the `reduction` embedding,
#         which Seurat's own subset() already restricts to just the clean
#         cells automatically. Gives a manifold that reflects only the
#         final population, not distorted by cells since excluded.

add_celltype <- function(whole_rds,
                         subset_labels,
                         subset_rds_list,
                         reduction  = "integ_harmony",
                         output_dir = NULL) {
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Load Whole Object and Validate
  # ═══════════════════════════════════════════════════════════════════════════
  # whole_rds is the FULL annotated_seurat.rds from the "Whole" branch. Auto-
  # detected junk/sparse cells are already physically absent by this point
  # (removed upstream in cluster_and_findmarkers) -- any "Exclude" value in
  # Lineage here is a manual template call on a whole cluster, not an auto-flag.
  
  whole_seurat <- load_smart(input_path = whole_rds)
  whole_ledger <- whole_seurat@misc$excluded_barcodes   # NULL if nothing was ever dropped upstream of Whole
  
  log_info(sample = "", step = "add_celltype",
           msg = glue::glue("Loaded Whole: {ncol(whole_seurat)} cells, {nrow(whole_seurat)} features."))
  
  if (!"Lineage" %in% colnames(whole_seurat@meta.data)) {
    log_error(sample = "", step = "add_celltype",
              msg = "Whole object metadata missing required column: Lineage.")
  }
  
  if (!reduction %in% names(whole_seurat@reductions)) {
    log_error(sample = "", step = "add_celltype",
              msg = glue::glue("Whole object missing reduction '{reduction}' — required for the ",
                               "final UMAP recompute. Run INTEGRATE first."))
  }
  
  if (!is.null(output_dir) && !dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Parse Subset Branch Arguments
  # ═══════════════════════════════════════════════════════════════════════════
  # Both arrive as comma-separated strings (Nextflow CLI convention) — split
  # and validate they're the same length before pairing them up. Zero subset
  # branches (e.g. first-ever run, nothing subclustered yet) is valid — every
  # cell just keeps Lineage as its final CellType in that case.
  
  labels    <- if (is.null(subset_labels) || subset_labels == "") {
    character(0)
  } else {
    trimws(strsplit(subset_labels, ",")[[1]])
  }
  rds_files <- if (is.null(subset_rds_list) || subset_rds_list == "") {
    character(0)
  } else {
    trimws(strsplit(subset_rds_list, ",")[[1]])
  }
  
  if (length(labels) != length(rds_files)) {
    log_error(sample = "", step = "add_celltype",
              msg = glue::glue("subset_labels ({length(labels)}) and subset_rds_list ",
                               "({length(rds_files)}) have mismatched lengths."))
  }
  
  log_info(sample = "", step = "add_celltype",
           msg = if (length(labels) == 0) {
             "No subset branches provided this run — every cell keeps Lineage as its final CellType."
           } else {
             glue::glue("Processing {length(labels)} subset branch(es): {paste(labels, collapse=', ')}")
           })
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Initialize CellType / CellType_Source / Exclude_Reason Defaults
  # ═══════════════════════════════════════════════════════════════════════════
  # Lineage: already present on the Whole object exactly as SCRNA_ANNOTATE_CLUSTERS
  #          wrote it — permanent record, never overwritten again.
  # CellType: working column, starts as a copy of Lineage, then gets
  #          overwritten per subset branch below. This is the FINAL call.
  # CellType_Source: provenance — "Lineage" (default/fallback) or the subset
  #          branch label that provided the final call for that cell.
  # Exclude_Reason: NA by default; set to "Manual_Exclude" in the branch loop
  #          below, or to the ledger's own reason for upstream-dropped cells.
  #          Initialized here (not after the loop) so the loop's own writes
  #          survive.
  
  whole_seurat$CellType        <- as.character(whole_seurat$Lineage)
  whole_seurat$CellType_Source <- "Lineage"
  whole_seurat$Exclude_Reason  <- NA_character_
  whole_seurat$Exclude_Reason[whole_seurat$Lineage == "Exclude"] <- "Lineage_Excluded_Manual"

# ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Overlay Each Subset Branch's Own CellType
  # ═══════════════════════════════════════════════════════════════════════════

  branch_summary <- list()
  seen_barcodes  <- character(0)   # tracks cells already assigned, to catch overlap across branches
  all_excluded   <- c()            # accumulates every branch's own drops + each branch's inherited ledger

  for (i in seq_along(labels)) {

    label         <- labels[i]
    subset_seurat <- load_smart(input_path = rds_files[i])

    # Validate required metadata BEFORE anything else touches CellType
    if (!"CellType" %in% colnames(subset_seurat@meta.data)) {
      log_error(sample = "", step = "add_celltype",
                msg = glue::glue("Subset branch '{label}' metadata missing required column: CellType."))
    }

    n_total    <- nrow(subset_seurat@meta.data)
    n_excluded <- nrow(subset_seurat@meta.data %>% dplyr::filter(CellType == "Exclude"))

    cells_keep <- colnames(subset_seurat)[subset_seurat@meta.data[["CellType"]] != "Exclude"]

    # Validate that every subset cell exists in Whole
    missing_barcodes <- setdiff(cells_keep, colnames(whole_seurat))
    if (length(missing_barcodes) > 0) {
      log_error(sample = "", step = "add_celltype",
                msg = glue::glue("Subset branch '{label}': {length(missing_barcodes)} barcode(s) ",
                                 "not found in Whole object. Branches must be derived from Whole via ",
                                 "SCRNA_SUBSET_CELLS — barcodes should always match exactly."))
    }

    # Validate that barcodes do not overlap between subset branches
    overlap_barcodes <- intersect(cells_keep, seen_barcodes)
    if (length(overlap_barcodes) > 0) {
      log_info(sample = "", step = "add_celltype",
               msg = glue::glue("WARNING: {length(overlap_barcodes)} barcode(s) in branch '{label}' ",
                                "already assigned by a previous branch — will be overwritten. Check ",
                                "subcluster_lineages for overlapping CellType filters."))
    }

    # Track all barcodes assigned a real CellType across all branches
    seen_barcodes <- union(seen_barcodes, cells_keep)

    # Assign cell types directly to whole object
    celltype <- setNames(as.character(subset_seurat$CellType), cells_keep)
    whole_seurat$CellType[cells_keep]        <- celltype
    whole_seurat$CellType_Source[cells_keep] <- label

    # This branch's own "Exclude" calls (contaminants identified within this branch)
    dropped_barcodes <- setdiff(colnames(subset_seurat), cells_keep)
    new_excluded <- NULL
    if (length(dropped_barcodes) > 0) {
      new_excluded <- stats::setNames(rep("Contaminant", length(dropped_barcodes)), dropped_barcodes)
    }

    # Carry forward this branch's own inherited ledger (e.g. from SCT_TRANSFORM's
    # <50-cell sample drop, or an earlier SUBSET_CELLS call), then combine with
    # this branch's own drops and fold into the running total across all branches.
    prior_excluded    <- subset_seurat@misc$excluded_barcodes   # NULL, or a named vector
    combined_excluded <- c(prior_excluded, new_excluded)
    all_excluded      <- c(all_excluded, combined_excluded)

    # Store branch summary
    branch_summary[[label]] <- list(n_cells    = n_total,
                                    n_excluded = n_excluded)

    # Log the execution step
    log_info(sample = "", step = "add_celltype",
             msg = glue::glue("Branch '{label}': {n_total} cells — ",
                              "{n_excluded} labeled 'Exclude' by this branch's own annotation, ",
                              "{n_total - n_excluded} assigned a real CellType."))
    rm(subset_seurat); gc()
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4b: Recover Cells Dropped Upstream (Whole/Initial/SCT/Branch Ledgers)
  # ═══════════════════════════════════════════════════════════════════════════
  # Barcodes physically absent from every branch's Final object, but present in
  # Whole's ledger or any branch's own carried-forward + own-Exclude ledger
  # (all_excluded, built in Section 4), were dropped somewhere upstream —
  # AutoQC at Whole, manual Exclude at Initial, SCTransform's <50-cell sample
  # removal, or a branch's own contaminant call — NOT genuinely "never
  # covered." These must be forced to Exclude rather than falling back to raw
  # Lineage. Guarded for the case of the combined ledger being completely
  # empty (e.g. first-ever run before upstream scripts wrote any ledger).

  combined_vector <- c(whole_ledger, all_excluded)

  if (length(combined_vector) > 0) {

    # Defensive: keep the earliest recorded reason if a barcode appears in multiple ledgers
    combined_vector <- combined_vector[!duplicated(names(combined_vector), fromLast = FALSE)]

    # Only barcodes NOT physically present in seen_barcodes (cells that successfully
    # passed all filters and got a real CellType in Section 4) need to be forced to Exclude.
    final_exclude_barcodes <- setdiff(names(combined_vector), seen_barcodes)

    if (length(final_exclude_barcodes) > 0) {
      final_reasons <- combined_vector[final_exclude_barcodes]

      # Stamp final metadata updates for all excluded cells
      whole_seurat@meta.data[final_exclude_barcodes, "CellType"]        <- "Exclude"
      whole_seurat@meta.data[final_exclude_barcodes, "CellType_Source"] <- "Excluded_Upstream"
      whole_seurat@meta.data[final_exclude_barcodes, "Exclude_Reason"]  <- final_reasons

      log_info(sample = "", step = "add_celltype",
               msg = glue::glue("{length(final_exclude_barcodes)} total cell(s) finalized as 'Exclude' ",
                                "via unified exclusion ledgers."))
    }
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5: Save Full Final Object
  # ═══════════════════════════════════════════════════════════════════════════
  
  if (!is.null(output_dir)) {
    full_out <- file.path(output_dir, "annotated_seurat_final.rds")
    saveRDS(object = whole_seurat, file = full_out)
    log_info(sample = "", step = "add_celltype",
             msg = glue::glue("Saved full final object: '{full_out}' ({ncol(whole_seurat)} cells, ",
                              "includes Exclude-labeled cells)."))
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 6: Build Clean Object
  # ═══════════════════════════════════════════════════════════════════════════
  # First point in the whole pipeline where CellType is genuinely final —
  # this is the only place "clean" (Exclude removed) is computed.
  
  clean_seurat <- base::subset(whole_seurat, subset = CellType != "Exclude")
  
  log_info(sample = "", step = "add_celltype",
           msg = glue::glue("Clean object: {ncol(clean_seurat)} / {ncol(whole_seurat)} cells retained ",
                            "({round(ncol(clean_seurat)/ncol(whole_seurat)*100, 1)}%) after removing ",
                            "Exclude-labeled cells."))
  
  if (ncol(clean_seurat) == 0) {
    log_error(sample = "", step = "add_celltype",
              msg = "Clean object contains 0 cells — every cell was labeled Exclude. Check upstream annotation.")
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 7: Final Marker Table — FindAllMarkers(group.by = "CellType")
  # ═══════════════════════════════════════════════════════════════════════════
  # No reclustering — CellType is the final grouping variable directly, no
  # Seurat cluster ID needed or generated. One authoritative marker table for
  # the whole project.
  
  Seurat::DefaultAssay(clean_seurat) <- "RNA"
  Seurat::Idents(clean_seurat)       <- "CellType"
  
  final_markers <- tryCatch({
    Seurat::FindAllMarkers(clean_seurat, only.pos = TRUE, min.pct = 0.1,
                           logfc.threshold = 0.25, verbose = FALSE)
  }, error = function(e) {
    log_warn(sample = "", step = "add_celltype",
             msg = glue::glue("FindAllMarkers(group.by='CellType') failed: {e$message}"))
    NULL
  })
  
  sig_markers <- NULL
  
  if (!is.null(final_markers) && nrow(final_markers) > 0 && !is.null(output_dir)) {
    
    sig_markers <- final_markers %>%
      dplyr::mutate(
        pct.1 = dplyr::if_else(pct.1 == 0, 0.001, pct.1),
        pct.2 = dplyr::if_else(pct.2 == 0, 0.001, pct.2),
        ratio = (pct.1 + 1e-4) / (pct.2 + 1e-4)) %>%
      dplyr::filter(p_val_adj <= 0.05, ratio > 1) %>%
      dplyr::arrange(cluster, dplyr::desc(pct.1))
    
    markers_file <- file.path(output_dir, "Markers.CellType_Final.xlsx")
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Sig_Markers")
    openxlsx::writeData(wb, "Sig_Markers", sig_markers)
    openxlsx::saveWorkbook(wb, file = markers_file, overwrite = TRUE)
    
    log_info(sample = "", step = "add_celltype",
             msg = glue::glue("Saved final marker table: '{markers_file}' — {nrow(sig_markers)} sig ",
                              "markers, {length(unique(sig_markers$cluster))} CellType(s)."))
  } else {
    log_warn(sample = "", step = "add_celltype",
             msg = "No significant final markers found or FindAllMarkers failed — skipping Markers.CellType_Final.xlsx.")
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 7b: Per-CellType Marker Review xlsx
  # ═══════════════════════════════════════════════════════════════════════════
  # Replaces the old per-branch candidate_markers_review.xlsx that used to be
  # generated in SCRNA_ANNOTATE_CLUSTERS. No master-marker-file lookup and no
  # "new vs known type" distinction any more — every final CellType (the
  # clean object's, Exclude already removed) gets a top-10-marker column,
  # straight off the sig_markers computed just above, for manual review.
  # Genes get added to the master celltype_marker_xlsx by hand afterward.
  # Ranking: p_val_adj asc, pct.1 desc, avg_log2FC desc — same cascade used
  # everywhere else in the pipeline for "top markers per group."
  
  if (!is.null(sig_markers) && nrow(sig_markers) > 0 && !is.null(output_dir)) {
    
    sig_markers$cluster <- as.character(sig_markers$cluster)
    
    review_list <- list()
    for (ct in sort(unique(sig_markers$cluster))) {
      review_list[[ct]] <- sig_markers %>%
        dplyr::filter(cluster == ct) %>%
        dplyr::arrange(p_val_adj, dplyr::desc(pct.1), dplyr::desc(avg_log2FC)) %>%
        dplyr::slice_head(n = 10) %>%
        dplyr::pull(gene) %>%
        unique()
    }
    
    max_len      <- max(lengths(review_list))
    padded       <- lapply(review_list, function(x) c(x, rep(NA_character_, max_len - length(x))))
    review_df    <- as.data.frame(padded, check.names = FALSE)
    
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Top_Markers")
    openxlsx::writeData(wb, "Top_Markers", review_df)
    review_file <- file.path(output_dir, "CellType_Markers_Review.xlsx")
    openxlsx::saveWorkbook(wb, file = review_file, overwrite = TRUE)
    
    log_info(sample = "", step = "add_celltype",
             msg = glue::glue("Saved marker review file: '{review_file}' — top 10 markers for ",
                              "{length(review_list)} final CellType(s). Review and add new genes to ",
                              "the master celltype_marker_xlsx by hand."))
  } else {
    log_warn(sample = "", step = "add_celltype",
             msg = "No sig_markers available — skipping CellType_Markers_Review.xlsx.")
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 8: Final UMAP — RunUMAP Only (No FindNeighbors, No Reclustering)
  # ═══════════════════════════════════════════════════════════════════════════
  # `reduction` (default "integ_harmony") was already automatically restricted
  # to just the clean cells by base::subset() above — Seurat subsets existing
  # reductions along with everything else, no manual embedding extraction
  # needed. RunUMAP directly on that gives a fresh manifold reflecting only
  # the final population, without the influence of cells since excluded.
  # No FindNeighbors: nothing downstream consumes an SNN graph at this stage.
  
  n_dims <- ncol(Seurat::Embeddings(clean_seurat, reduction = reduction))
  
  clean_seurat <- Seurat::RunUMAP(clean_seurat,
                                  reduction      = reduction,
                                  dims           = 1:n_dims,
                                  reduction.name = "umap_final",
                                  verbose        = FALSE)
  
  log_info(sample = "", step = "add_celltype",
           msg = glue::glue("Final UMAP computed on '{reduction}' ({n_dims} dims) restricted to ",
                            "{ncol(clean_seurat)} clean cells — reduction.name = 'umap_final'."))
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 9: Final Plots
  # ═══════════════════════════════════════════════════════════════════════════
  
  if (!is.null(output_dir)) {
    
    plot_seurat(clean_seurat,
                features        = "CellType",
                feature_type    = "metadata",
                reduction       = "umap_final",
                plot_unlabelled = TRUE,
                filename        = "UMAP_Final_CellType",
                output_dir      = output_dir)
    
    plot_seurat(clean_seurat,
                features        = "Lineage",
                feature_type    = "metadata",
                reduction       = "umap_final",
                plot_unlabelled = TRUE,
                filename        = "UMAP_Final_Lineage",
                output_dir      = output_dir)
    
    plot_seurat(clean_seurat,
                features        = "Sample",
                feature_type    = "metadata",
                reduction       = "umap_final",
                split_col       = "Sample",
                plot_unlabelled = TRUE,
                filename        = "UMAP_Final_Sample",
                output_dir      = output_dir)
    
    log_info(sample = "", step = "add_celltype",
             msg = "Final plots saved: UMAP_Final_CellType, UMAP_Final_Lineage, UMAP_Final_Sample.")
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 10: Save Clean Object
  # ═══════════════════════════════════════════════════════════════════════════
  
  if (!is.null(output_dir)) {
    clean_out <- file.path(output_dir, "annotated_seurat_final_clean.rds")
    saveRDS(object = clean_seurat, file = clean_out)
    log_info(sample = "", step = "add_celltype",
             msg = glue::glue("Saved clean final object: '{clean_out}' ({ncol(clean_seurat)} cells) — ",
                              "use this for downstream DE / composition / pseudobulk analysis."))
  } else {
    log_info(sample = "", step = "add_celltype", msg = "No output_dir provided — skipping save.")
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 11: Summary
  # ═══════════════════════════════════════════════════════════════════════════
  
  celltype_table <- sort(table(clean_seurat$CellType), decreasing = TRUE)
  log_info(sample = "", step = "add_celltype", msg = "Final CellType distribution (clean object):")
  for (ct in names(celltype_table)) {
    pct <- round(celltype_table[ct] / ncol(clean_seurat) * 100, 1)
    log_info(sample = "", step = "add_celltype",
             msg = glue::glue("  {ct}: {celltype_table[ct]} cells ({pct}%)"))
  }
  
  log_info(sample = "", step = "add_celltype", msg = "Completed successfully.")
  
  return(invisible(list(full = whole_seurat, clean = clean_seurat)))
}

# ---- 🚀 4. Smart Execution (Nextflow Only) ----
# CLI arg order matches SCRNA_add_celltype.nf's Rscript call exactly:
#   whole_rds, subset_labels, subset_rds_list, reduction, output_dir

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)

  get_arg <- function(idx, default = NULL) {
    if (idx > length(args)) return(default)
    val <- args[idx]
    if (is.na(val) || val == "" || val == "null" || val == "NULL") return(default)
    return(val)
  }

  add_celltype(
    whole_rds       = get_arg(1),
    subset_labels   = get_arg(2, default = ""),
    subset_rds_list = get_arg(3, default = ""),
    reduction       = get_arg(4, default = "integ_harmony"),
    output_dir      = get_arg(5, default = ".")
  )
}
