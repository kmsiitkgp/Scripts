# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)
  library(glue)
  library(ggplot2)
  library(openxlsx)
  library(Seurat)
  library(SeuratObject)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(presto)
})

# ---- 🛠️ 2. Smart Setup (Find & source UTILS.R) ----

get_modules_dir <- function() {
  if (.Platform$OS.type == "windows") {
    return("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/Scripts/nextflow/modules")
  }
  if (interactive()) {
    return(file.path(getwd(), "modules"))
  }
  initial.options <- commandArgs(trailingOnly = FALSE)
  file_arg        <- grep("--file=", initial.options, value = TRUE)
  if (length(file_arg) == 0) stop("Cannot detect modules directory path!")
  modules_dir     <- dirname(sub("--file=", "", file_arg))
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

# Fixed reference list of the 6 top-level lineages used across the whole
# project (matches params.scrna.seurat.subcluster_lineages in the pipeline
# yaml). Never varies by dataset -- hardcoded here rather than threaded
# through as a CLI arg since it's a permanent project constant, not a
# per-run parameter.
LINEAGE_LEVELS <- c("Epithelial", "Mesenchymal", "Hematopoietic",
                    "Endothelial", "Neural", "Germline")

# ---- 🧬 3. Function Definition (Always Runs) ----

# ─────────────────────────────────────────────────────────────────────────────
# MANUAL REVIEW TRICK: using call_distribution to resolve Uncertain/Mixed
# ─────────────────────────────────────────────────────────────────────────────
# When final_call (primary resolution, e.g. res 1.2) is "Uncertain" or "Mixed",
# don't just guess -- check that cluster's call_distribution in reasons.xlsx.
# Bucket order is always Retain / Exclude / Uncertain / Mixed (4 numbers,
# sum to n_res = 6), e.g. "3/0/2/1" means 3 of the 6 resolutions called this
# cluster's cells Retain, 0 called Exclude, 2 Uncertain, 1 Mixed.
#
# Rule of thumb (validated by hand against real Neural_Initial output):
#   - If the Retain bucket is the largest nonzero count (even if not an
#     outright majority, e.g. "3/0/2/1" or "4/0/2/0") -> lean RETAIN.
#     The primary resolution alone was underpowered/ambiguous (often because
#     n_expected_hits = 0 at that one resolution), but other resolutions with
#     different cluster boundaries did find real lineage signal.
#   - If the Retain bucket is 0 across ALL 6 resolutions (e.g. "0/1/5/0",
#     "0/2/0/4", "0/4/2/0") -> lean EXCLUDE. Not one resolution, at any
#     granularity, ever supported keeping this cluster.
#   - Cross-check against UMAP position as a final sanity check either way.
# This is a manual-review heuristic only -- not automated in final_call
# itself, to avoid a second silent decision-maker overriding the primary
# resolution's verdict (see design discussion: only the primary resolution's
# call is ever supposed to drive CellType/Exclude downstream).
# ─────────────────────────────────────────────────────────────────────────────


# ---- HELPER 1: run_clustering() ----

# Purely mechanical: FindNeighbors + FindClusters (all resolutions) + RunUMAP.
# No sparse/junk decisions live here — those are QC judgments that belong in
# identify_junk_and_sparse(). Always produces reduction "umap_harmony" and
# cluster columns "harmony_res{X}" — one call per object (Whole or a subset
# branch), no iterative re-clustering.
#
# Returns: seurat_obj with graphs, cluster columns, UMAP; cluster_cols vector.

run_clustering <- function(seurat_obj,
                           resolutions,
                           step_label = "CLUSTER") {
  
  set.seed(1234)
  
  n_cells <- ncol(seurat_obj)
  
  harmony_reduction <- seurat_obj@misc$integration_params$harmony_reduction %||% "integ_harmony"
  max_dims          <- seurat_obj@misc$integration_params$n_dims             %||% 30
  
  dims_available <- min(max_dims,
                        ncol(Seurat::Embeddings(seurat_obj, reduction = harmony_reduction)))
  
  # ── FindNeighbors ───────────────────────────────────────────────────────────
  # annoy (approximate nearest neighbors) is dramatically faster than exact kNN
  # at large cell counts with negligible accuracy loss for graph construction.
  # k.param = 30 balances local structure resolution with computational cost —
  # standard for datasets of this scale.
  # graph.name uses "harmony" prefix so downstream code can identify which
  # integration method produced these graphs.
  
  log_info(sample = "", step = step_label,
           msg = glue::glue("FindNeighbors: {n_cells} cells, {dims_available} dims, annoy."))
  
  seurat_obj <- Seurat::FindNeighbors(
    object       = seurat_obj,
    reduction    = harmony_reduction,
    dims         = 1:dims_available,
    k.param      = 30,
    nn.method    = "annoy",
    annoy.metric = "euclidean",
    graph.name   = c("harmony.nn", "harmony.snn"),
    verbose      = FALSE
  )
  
  # ── FindClusters at all resolutions ─────────────────────────────────────────
  # Leiden (algorithm = 4) over Louvain — Leiden guarantees well-connected
  # communities and is the current scRNAseq standard.
  #
  # n.start and n.iter reduced from defaults (10, 10):
  #   Default = 10 starts × 10 iterations × 6 resolutions = 600 Leiden passes.
  #   At 500k+ cells this takes hours. At our scale the graph structure is
  #   rich enough that 3 starts × 5 iterations converges to equivalent
  #   biological results. This is not a shortcut — it is appropriate scaling
  #   of stochastic initialization effort to dataset size.
  
  cluster_cols <- c()
  
  for (res in resolutions) {
    
    cluster_name <- paste0("harmony_res", res)
    
    seurat_obj <- Seurat::FindClusters(
      object       = seurat_obj,
      graph.name   = "harmony.snn",
      cluster.name = cluster_name,
      resolution   = res,
      algorithm    = 4,
      n.start      = 3,
      n.iter       = 5,
      random.seed  = 1,
      verbose      = FALSE
    )
    
    cluster_cols <- c(cluster_cols, cluster_name)
    
    log_info(sample = "", step = step_label,
             msg = glue::glue("Leiden clustering done — res {res} -> '{cluster_name}': ",
                              "{length(unique(seurat_obj@meta.data[[cluster_name]]))} clusters."))
  }
  
  # ── RunUMAP ──────────────────────────────────────────────────
  # Run on full Harmony embedding. n.neighbors matches k.param in FindNeighbors
  # for consistency — the UMAP neighborhood structure mirrors the clustering
  # neighborhood structure.
  
  log_info(sample = "", step = step_label, msg = "Running UMAP -> 'umap_harmony'.")
  
  seurat_obj <- Seurat::RunUMAP(
    object         = seurat_obj,
    reduction      = harmony_reduction,
    dims           = 1:dims_available,
    n.neighbors    = 30L,
    reduction.name = "umap_harmony",
    return.model   = FALSE,
    verbose        = FALSE
  )
  
  gc()
  
  # SCT as default assay for all downstream analysis — Pearson residuals are
  # appropriate for visualization and module scoring. RNA is used only for
  # FindAllMarkers (log-normalized counts, set explicitly there).
  Seurat::DefaultAssay(seurat_obj)     <- "SCT"
  Seurat::VariableFeatures(seurat_obj) <- Seurat::VariableFeatures(seurat_obj, assay = "SCT")
  
  return(list(seurat_obj = seurat_obj, cluster_cols = cluster_cols))
}

# ---- HELPER 2: find_markers() ----

# Runs FindAllMarkers at every resolution on the FULL cluster set — sparse
# micro-clusters are intentionally left in (removing them first would change
# the pct.2/ratio denominator for every other cluster's specificity stats,
# an unforced change to the marker statistics themselves, not just
# bookkeeping). Saves one xlsx per resolution containing only significant
# markers (padj <= 0.05) — no ranking, no matrix, no unfiltered top-N sheet.
#
# Ranking by logFC/pct.1/ratio was deliberately dropped: these metrics
# penalize real markers of tightly related subclusters (e.g. Krt genes
# across 5 epithelial subclusters have high pct.1/pct.2 but modest ratio/
# logFC because they're shared across the family, not exclusive to one
# subcluster) and equally over-reward genes with high logFC but low
# specificity (e.g. a neuronal marker leaking into an epithelial cluster).
# Neither ranking axis reliably tracks "real marker" — raw significant
# markers are left for the caller (identify_junk_and_sparse or manual
# review) to interpret without a pre-baked, possibly misleading sort order.
#
# Returns: named list of all_markers data frames (padj <= 0.05), one per
# resolution (column name = cluster column name).

compute_seurat_style_fc <- function(mat, groups, pseudocount = 1) {
  groups   <- as.character(groups)
  clusters <- sort(unique(groups))
  
  mat_exp    <- mat
  mat_exp@x  <- expm1(mat_exp@x)
  
  total_sum   <- Matrix::rowSums(mat_exp)
  total_ncell <- ncol(mat_exp)
  
  fc_list <- lapply(clusters, function(cl) {
    idx   <- which(groups == cl)
    n_in  <- length(idx)
    n_out <- total_ncell - n_in
    
    sum_in  <- Matrix::rowSums(mat_exp[, idx, drop = FALSE])
    sum_out <- total_sum - sum_in
    
    data.frame(
      gene       = rownames(mat_exp),
      cluster    = cl,
      avg_log2FC = log2((sum_in + pseudocount) / n_in) -
        log2((sum_out + pseudocount) / n_out),
      stringsAsFactors = FALSE
    )
  })
  
  dplyr::bind_rows(fc_list)
}

find_markers <- function(seurat_obj,
                         label,
                         cluster_cols,
                         output_dir = NULL,
                         step_label = "FINDMARKERS") {
  
  # Switch to RNA for FindAllMarkers.
  # RNA log-normalized counts are appropriate for marker finding across multiple
  # samples — SCT Pearson residuals are not comparable across samples because
  # SCT fits are sample-specific. RNA data layer is uniform across all samples.
  Seurat::DefaultAssay(seurat_obj) <- "RNA"
  
  # ── Confirm presto is available ─────────────────────────────────────────
  presto_available <- requireNamespace("presto", quietly = TRUE)
  log_info(sample = "", step = step_label,
           msg = glue::glue("presto available for accelerated Wilcoxon: {presto_available}"))
  if (!presto_available) {
    log_error(sample = "", step = step_label,
              msg = "presto not found — FindAllMarkers will use slow base-R Wilcoxon test.")
  }
  
  mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
  
  markers_list <- list()
  
  for (col in cluster_cols) {
    
    log_info(sample = "", step = step_label,
             msg = glue::glue("Finding markers for: '{col}'."))
    
    Seurat::Idents(seurat_obj) <- seurat_obj@meta.data[[col]]
    groups <- seurat_obj@meta.data[[col]]
    
    # ── FindAllMarkers ─────────────────────────────────────────────────────────
    
    all_markers <- tryCatch({
      raw <- presto::wilcoxauc(X = mat, y = groups)
      fc  <- compute_seurat_style_fc(mat, groups)
      
      raw %>%
        dplyr::rename(gene = feature, cluster = group,
                      p_val = pval, p_val_adj = padj,
                      pct.1 = pct_in, pct.2 = pct_out) %>%
        dplyr::mutate(pct.1 = pct.1 / 100, pct.2 = pct.2 / 100) %>%
        dplyr::group_by(cluster) %>%
        dplyr::mutate(p_val_adj = p.adjust(p_val, method = "bonferroni", n = nrow(mat))) %>%
        dplyr::ungroup() %>%
        dplyr::inner_join(fc, by = c("gene", "cluster")) %>%
        # pre-filters matching find_markers()'s FindAllMarkers call
        dplyr::filter(avg_log2FC > 0, pmax(pct.1, pct.2) >= 0.1, avg_log2FC >= 0.25)
    }, error = function(e) {
      log_error(sample = "", step = step_label,
                msg = glue::glue("Marker finding (presto) failed for '{col}': {e$message}"))
      return(NULL)
    })
    
    # only.pos = TRUE — we want defining markers, not repressed genes.
    # logfc.threshold = 0.25 and min.pct = 0.1 pre-filter to reduce runtime.
    # These are permissive thresholds — we do not want to miss real markers.
    # Stricter filtering (padj <= 0.05) applied post-hoc, no other filtering.
    # all_markers <- tryCatch({
    #   Seurat::FindAllMarkers(
    #     object          = seurat_obj,
    #     assay           = "RNA",
    #     layer           = "data",
    #     logfc.threshold = 0.25,
    #     test.use        = "wilcox",
    #     min.pct         = 0.1,
    #     min.diff.pct    = -Inf,
    #     only.pos        = TRUE,
    #     verbose         = FALSE)
    # }, error = function(e) {
    #   log_error(sample = "", step = step_label,
    #             msg = glue::glue("FindAllMarkers failed for '{col}': {e$message}"))
    # })
    
    # Compute ratio = pct.1 / pct.2 — fraction of cells in cluster expressing
    # the gene divided by fraction outside. ratio > 1 means the gene is more
    # specific to this cluster than to the rest of the dataset.
    # Small epsilon added to avoid division by zero on rare all-or-nothing genes.
    # NOTE: ratio is computed and retained for reference/manual review, but is
    # NOT used to rank or filter here — see rationale above.
    if (is.null(all_markers) || nrow(all_markers) == 0) {
      log_warn(sample = "", step = step_label,
               msg = glue::glue("No markers found for '{col}'; skipping formatting."))
      # Keep the same column structure as the normal case (just 0 rows) so
      # downstream code (identify_junk_and_sparse's dplyr::filter on cluster/
      # ratio/avg_log2FC/pct.1) doesn't break on a missing column — a bare
      # data.frame() with zero columns would error there instead of just
      # returning zero matching rows.
      sig_markers <- data.frame(
        p_val            = numeric(0),
        avg_log2FC       = numeric(0),
        pct.1            = numeric(0),
        pct.2            = numeric(0),
        p_val_adj        = numeric(0),
        cluster          = character(0),
        gene             = character(0),
        ratio            = numeric(0),
        stringsAsFactors = FALSE
      )
    } else {
      sig_markers <- all_markers %>%
        dplyr::mutate(
          pct.1 = dplyr::if_else(pct.1 == 0, 0.001, pct.1),
          pct.2 = dplyr::if_else(pct.2 == 0, 0.001, pct.2),
          ratio = (pct.1 + 1e-4) / (pct.2 + 1e-4)) %>%
        dplyr::filter(p_val_adj <= 0.05, base::round(ratio, 2) > 1) %>%
        dplyr::select(p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster, gene, ratio) %>%
        dplyr::arrange(cluster, dplyr::desc(pct.1))
      
      if (!is.null(output_dir)) {
        if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
        
        file_name <- file.path(output_dir, glue::glue("Markers_{col}_{label}.xlsx"))
        wb        <- openxlsx::createWorkbook()
        openxlsx::addWorksheet(wb, "Sig_Markers")
        openxlsx::writeData(wb, "Sig_Markers", sig_markers)
        openxlsx::saveWorkbook(wb, file = file_name, overwrite = TRUE)
        
        log_info(sample = "", step = step_label,
                 msg = glue::glue("Saved: '{file_name}' — {nrow(sig_markers)} sig markers, ",
                                  "{length(unique(sig_markers$cluster))} clusters."))
      }
    }
    
    markers_list[[col]] <- sig_markers
  }
  
  return(markers_list)
}

# ---- HELPER 3: identify_junk_and_sparse() ----

# Both QC judgments live here, side by side, since they're the same category
# of decision (what to exclude and why) even though the mechanism differs:
#   - Sparse: <=sparse_threshold cells in a cluster at ANY resolution tested
#     (union across resolutions) — graph-connectivity issue, typically forms
#     at high resolutions.
#   - Junk: fewer than MIN_MARKERS genes simultaneously pass ratio >= RATIO_CUTOFF
#     (1.1), avg_log2FC >= LFC_CUTOFF (0.58), pct.1 >= PCT1_MIN (0.30) at a given
#     resolution, intersected across ALL resolutions (a cell is only junk if
#     consistently homeless everywhere) — lack of transcriptional identity.
# Kept distinguishable in the returned labels (sparse vs junk) since they're
# different failure modes, even though both result in exclusion.
#
# CURRENT LIMITATION (flagged for discussion, not yet resolved): the
# ratio/logFC/pct.1 junk criteria can misfire on real markers shared across
# a family of closely related subclusters (e.g. Krt genes across 5
# epithelial subclusters — high pct.1/pct.2 but modest ratio/logFC since
# the gene isn't exclusive to one subcluster). pct.1 alone may be a more
# reliable signal in that scenario than ratio/logFC. Thresholds intentionally
# left as originally validated (see project notes) pending further review.
#
# Returns: list with
#  $markerless_barcodes : character vector, consensus markerless (marker-based)
#  $sparse_barcodes     : character vector, sparse-cluster cells (any resolution)
#  $excluded_barcodes   : union of the two
#  $cluster_audit_df    : one row per cluster per resolution — is_junk,
#                         is_sparse, n_cells, top gene info, QC metrics

identify_junk_and_sparse <- function(seurat_obj,
                                     cluster_cols,
                                     markers_list,
                                     sparse_threshold,
                                     output_dir = NULL,
                                     step_label = "JUNK_SPARSE",
                                     umi_col    = "nUMIs",
                                     gene_col   = "nGenes") {
  
  RATIO_CUTOFF <- 1.1
  LFC_CUTOFF   <- 0.58
  PCT1_MIN     <- 0.30
  MIN_MARKERS  <- 2
  
  has_qc <- all(c(umi_col, gene_col) %in% colnames(seurat_obj@meta.data))
  
  junk_sets          <- list()
  sparse_cells        <- c()
  cluster_audit_rows <- list()
  
  for (col in cluster_cols) {
    
    sig_markers  <- markers_list[[col]]
    all_clusters <- unique(as.character(seurat_obj@meta.data[[col]]))
    
    # ── Sparse check ────────────────────────────────────────────────────────
    cluster_counts <- table(seurat_obj@meta.data[[col]])
    sparse_ids     <- names(cluster_counts[cluster_counts <= sparse_threshold])
    if (length(sparse_ids) > 0) {
      sparse_here <- rownames(seurat_obj@meta.data)[
        as.character(seurat_obj@meta.data[[col]]) %in% sparse_ids]
      sparse_cells <- union(sparse_cells, sparse_here)
    }
    
    # ── Junk check per cluster ─────────────────────────────────────────────
    junk_clusters_this_res <- c()
    
    for (cl in all_clusters) {
      
      cl_markers <- sig_markers %>% dplyr::filter(as.character(cluster) == cl)
      
      strong_markers <- cl_markers %>%
        dplyr::filter(ratio >= RATIO_CUTOFF, avg_log2FC >= LFC_CUTOFF, pct.1 >= PCT1_MIN)
      
      is_junk <- nrow(strong_markers) < MIN_MARKERS
      if (is_junk) junk_clusters_this_res <- c(junk_clusters_this_res, cl)
      
      # Prefer the best-supported strong marker if one exists; otherwise fall
      # back to the best candidate by ratio among all sig markers, so even
      # "no strong marker" clusters (like cluster 1) show something informative
      # rather than a generic housekeeping gene
      top_gene_row <- if (nrow(strong_markers) > 0) {
        strong_markers %>% dplyr::arrange(dplyr::desc(ratio)) %>% dplyr::slice(1)
      } else {
        cl_markers %>% dplyr::arrange(dplyr::desc(ratio)) %>% dplyr::slice(1)
      }
      
      cl_bc   <- rownames(seurat_obj@meta.data)[as.character(seurat_obj@meta.data[[col]]) == cl]
      cl_meta <- seurat_obj@meta.data[cl_bc, , drop = FALSE]
      
      cluster_audit_rows[[paste(col, cl, sep = "_")]] <- data.frame(
        cluster_col      = col,
        cluster          = cl,
        n_cells          = length(cl_bc),
        is_sparse        = cl %in% sparse_ids,
        is_junk          = is_junk,
        n_strong_markers = nrow(strong_markers),
        top_gene         = if (nrow(top_gene_row) > 0) top_gene_row$gene else NA_character_,
        top_gene_pct1    = if (nrow(top_gene_row) > 0) round(top_gene_row$pct.1, 3) else NA_real_,
        top_gene_ratio   = if (nrow(top_gene_row) > 0) round(top_gene_row$ratio, 3) else NA_real_,
        top_gene_log2FC  = if (nrow(top_gene_row) > 0) round(top_gene_row$avg_log2FC, 3) else NA_real_,
        mean_nUMI        = if (has_qc) round(mean(cl_meta[[umi_col]]),  1) else NA_real_,
        mean_nGenes      = if (has_qc) round(mean(cl_meta[[gene_col]]), 1) else NA_real_,
        stringsAsFactors = FALSE
      )
      
      if (is_junk) {
        top_strong <- if (nrow(strong_markers) > 0) {
          strong_markers %>% dplyr::arrange(dplyr::desc(ratio)) %>% dplyr::slice(1) %>%
            glue::glue_data("{gene}(ratio={round(ratio,2)},pct1={pct.1},lFC={round(avg_log2FC,2)})")
        } else "none"
        
        log_info(sample = "", step = step_label,
                 msg = glue::glue("Junk cluster: '{cl}' at {col} | sig_markers={nrow(cl_markers)} | ",
                                  "strong_markers={nrow(strong_markers)} | best_strong={top_strong}"))
      }
    }
    
    # Always record an entry for this resolution — even an empty one. A resolution
    # with ZERO junk clusters means "nothing is junk at this resolution," which
    # must veto markerless status for every barcode once intersected below.
    # Previously this only ran `if (length(junk_clusters_this_res) > 0)`, which
    # silently dropped clean resolutions from junk_sets entirely instead of
    # contributing an empty set — letting Reduce(intersect, ...) skip them and
    # over-include barcodes as "markerless" that a clean resolution should have
    # vetoed.
    junk_sets[[col]] <- rownames(seurat_obj@meta.data)[
      as.character(seurat_obj@meta.data[[col]]) %in% junk_clusters_this_res]
  }
  
  markerless_barcodes <- if (length(junk_sets) == 0) character(0) else Reduce(intersect, junk_sets)
  excluded_barcodes   <- union(markerless_barcodes, sparse_cells)
  cluster_audit_df    <- dplyr::bind_rows(cluster_audit_rows)
  
  cluster_audit_df$in_consensus_junk <- NA_integer_
  for (i in seq_len(nrow(cluster_audit_df))) {
    col   <- cluster_audit_df$cluster_col[i]
    cl    <- cluster_audit_df$cluster[i]
    cl_bc <- rownames(seurat_obj@meta.data)[as.character(seurat_obj@meta.data[[col]]) == cl]
    cluster_audit_df$in_consensus_junk[i] <- sum(cl_bc %in% markerless_barcodes)
  }
  
  if (length(markerless_barcodes) > 0 && has_qc) {
    clean_bc  <- setdiff(colnames(seurat_obj), markerless_barcodes)
    umi_ratio <- mean(seurat_obj@meta.data[markerless_barcodes, umi_col]) /
      mean(seurat_obj@meta.data[clean_bc, umi_col])
    
    log_info(sample = "", step = step_label,
             msg = glue::glue("Markerless QC audit | n_markerless={length(markerless_barcodes)} | ",
                              "umi_ratio={round(umi_ratio, 2)}"))
    
    if (umi_ratio > 0.8) {
      log_warn(sample = "", step = step_label,
               msg = glue::glue("WARNING: markerless barcodes have similar UMIs to clean cells ",
                                "(ratio={round(umi_ratio, 2)}). Review thresholds."))
    }
  }
  
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Cluster_Junk_Status")
    openxlsx::writeData(wb, "Cluster_Junk_Status",
                        cluster_audit_df %>%
                          dplyr::arrange(cluster_col, dplyr::desc(is_junk), dplyr::desc(is_sparse), dplyr::desc(n_cells)))
    openxlsx::saveWorkbook(wb, file = file.path(output_dir, "Junk_Cluster_Audit.xlsx"), overwrite = TRUE)
  }
  
  log_info(sample = "", step = step_label,
           msg = glue::glue("{length(markerless_barcodes)} markerless + {length(sparse_cells)} sparse ",
                            "= {length(excluded_barcodes)} unique cells excluded."))
  
  return(list(
    markerless_barcodes  = markerless_barcodes,
    sparse_barcodes      = sparse_cells,
    excluded_barcodes    = excluded_barcodes,
    cluster_audit_df     = cluster_audit_df
  ))
}


# ---- HELPER 4: generate_annotation_template() ----

# Builds one manual-annotation template xlsx per cluster resolution, entirely
# from objects already computed in this script — no dependency on module
# scoring or automatic annotation. Sheet 1 ("Annotation") is what gets filled
# in by hand; Sheet 2 ("Reference_Lists") is a copy-paste reference only and
# is ALWAYS written with two independent columns (no row-to-row mapping
# between them, just two side-by-side reference lists, padded with blanks
# to the longer column's length):
#   - Lineage  : the fixed 6-value project constant (LINEAGE_LEVELS).
#   - CellType : column headers of celltype_marker_xlsx (empty if that file
#                isn't supplied).
#
# Sheet 1's editable label column is named by `label`: "Lineage" when
# label == "Whole", "CellType" for every other (subset) label -- matching
# what SCRNA_ANNOTATE_CLUSTERS will validate/write. The output filename also
# encodes label: Annotation_{cluster_col}_{label}.xlsx.
#
# Marker ranking per cluster (3-criteria cascade, since large cell counts
# routinely produce many genes tied at p_val_adj == 0):
#   1. p_val_adj  ascending
#   2. pct.1      descending
#   3. avg_log2FC descending
# Top 10 genes per cluster after that ordering, comma-separated into one
# "top_markers" cell — 1 row per cluster, not 1 row per gene.

generate_annotation_template <- function(seurat_obj,
                                         label,
                                         cluster_cols,
                                         markers_list,
                                         celltype_marker_xlsx = NULL,
                                         output_dir = NULL,
                                         step_label = "ANNOTATION_TEMPLATE") {
  
  if (is.null(output_dir)) {
    log_info(sample = "", step = step_label, msg = "No output_dir — skipping template generation.")
    return(invisible(NULL))
  }
  
  template_dir <- file.path(output_dir, "Annotation_Templates")
  if (!dir.exists(template_dir)) dir.create(template_dir, recursive = TRUE)
  
  # Sheet 1's editable column: Lineage for Whole, CellType for every subset.
  target_col <- if (identical(label, "Whole")) "Lineage" else "CellType"
  
  # Sheet 2 — always both reference lists, independent of label.
  known_celltypes <- character(0)
  if (!is.null(celltype_marker_xlsx) && file.exists(celltype_marker_xlsx)) {
    marker_df       <- openxlsx::read.xlsx(celltype_marker_xlsx)
    known_celltypes <- colnames(marker_df)
  }
  max_len   <- max(length(LINEAGE_LEVELS), length(known_celltypes))
  sheet2 <- tibble::tibble(
    Lineage  = c(LINEAGE_LEVELS,  rep(NA_character_, max_len - length(LINEAGE_LEVELS))),
    CellType = c(known_celltypes, rep(NA_character_, max_len - length(known_celltypes)))
  )
  
  for (col in cluster_cols) {
    
    cluster_vals <- seurat_obj@meta.data[[col]]
    n_cells_tbl  <- table(cluster_vals)
    
    top_markers_df <- markers_list[[col]] %>%
      dplyr::filter(!grepl("^[Mm][Tt]-", gene)) %>%
      dplyr::arrange(dplyr::desc(pct.1), dplyr::desc(avg_log2FC)) %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_head(n = 100) %>%
      dplyr::summarise(top_markers = paste(gene, collapse = ", "), .groups = "drop")
    
    sheet1 <- tibble::tibble(
      cluster_id    = names(n_cells_tbl),
      n_cells       = as.integer(n_cells_tbl),
      pct_of_total  = round(as.integer(n_cells_tbl) / sum(n_cells_tbl) * 100, 1)
    ) %>%
      dplyr::left_join(top_markers_df, by = c("cluster_id" = "cluster")) %>%
      dplyr::mutate(!!target_col := NA_character_) %>%
      dplyr::relocate(dplyr::all_of(target_col), .after = cluster_id) %>%
      dplyr::arrange(as.numeric(as.character(cluster_id)))
    
    file_name <- file.path(template_dir, glue::glue("Annotation_{col}_{label}.xlsx"))
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Annotation")
    openxlsx::writeData(wb, "Annotation", sheet1)
    openxlsx::addWorksheet(wb, "Reference_Lists")
    openxlsx::writeData(wb, "Reference_Lists", sheet2)
    openxlsx::saveWorkbook(wb, file = file_name, overwrite = TRUE)
    
    log_info(sample = "", step = step_label,
             msg = glue::glue("Saved: '{file_name}' — {nrow(sheet1)} clusters, ",
                              "editable column '{target_col}', {length(known_celltypes)} known cell type(s) ",
                              "in Reference_Lists."))
  }
}


# ---- DATA: lineage_panels ----

# Universal mammalian major-lineage marker panels, one entry per value in
# LINEAGE_LEVELS. Used only by score_cluster_lineages() / auto_annotate_lineage()
# below. A cluster does NOT need to express every gene in a panel — a coherent
# signal of a few genes from one lineage, with low competing-lineage scores,
# is treated as evidence for that lineage. Panels are intentionally broad and
# redundant; subtype-specific genes are retained when they contribute useful
# cross-tissue lineage evidence.

lineage_panels <- list(
  
  Epithelial = list(
    genes = c(
      # pan-epithelial structural / junctional
      "Epcam","Cdh1","Krt8","Krt18","Krt19","Krt7","Krt17","Krt5","Krt14",
      # epithelial junction / polarity
      "Cldn3","Cldn4","Cldn7","Cldn10","Cldn18","Ocln","Tjp1","Msl1",
      # epithelial transcriptional / differentiation program
      "Elf3","Grhl1","Grhl2","Klf5","Klf4","Perp",
      # epithelial mucosal / secretory programs
      "Muc1","Muc13","Muc20","Agr2","Spdef","Klf3",
      "Muc2","Clca1","Tff3","Fcgbp","Spink4","Reg4",
      "Gcg","Pyy","Dclk1","Aqp8","Slc26a3","Lgr5",
      # epithelial differentiation / tissue-associated
      "Hnf4g","Krt20",
      # keratinization / stratified epithelial programs
      "Krt1","Krt2","Krt10","Krt16","Krt23"
    ),
    master_regulators = c()
  ),
  
  Mesenchymal = list(
    genes = c(
      # fibroblast / ECM core
      "Col1a1","Col1a2","Col3a1","Col5a1","Col5a2","Col6a1","Col6a2","Col6a3",
      "Dcn","Lum","Fbln1","Fbln2","Fbn1","Lox","Loxl1","Pdgfra",
      # stromal / fibroblast-associated
      "Pdgfrb","Col14a1","Pi16","C7","C1s","C1r","Dpt","Mfap4","Sparc","Tcf21",
      # contractile / smooth-muscle / mural mesenchyme
      "Acta2","Tagln","Myh11","Cnn1","Des","Myl9","Mylk","Rgs5","Notch3","Mcam",
      "Prkg1",
      # general mesenchymal-associated
      "Vim","Postn","Thy1","Cd34","S100a10"
    ),
    master_regulators = c()
  ),
  
  Hematopoietic = list(
    genes = c(
      # broad hematopoietic
      "Ptprc","Cd52","Hcls1","Lsp1",
      # lymphoid
      "Cd3d","Cd3e","Cd3g","Trbc1","Trbc2","Lck","Il7r","Nkg7",
      # B-cell-associated
      "Cd79a","Cd79b","Cd19","Ms4a1","Cd37","Bank1","H2-Aa","Ighm","Pax5",
      # myeloid
      "Lyz2","Lyz1","Csf1r","Ctss","Tyrobp","Aif1","Fcerg1","Fcer1g","Adgre1","Cd68",
      # macrophage / mononuclear phagocyte
      "C1qa","C1qb","C1qc","Mrc1","Cd163","F13a1",
      # mast-cell-associated
      "Ms4a2","Cpa3","Kit",
      # plasma / antibody-producing
      "Jchain","Mzb1","Sdc1","Xbp1","Igkc","Iglc2"
    ),
    master_regulators = c("Ptprc")
  ),
  
  Endothelial = list(
    genes = c(
      "Pecam1","Cdh5","Kdr","Flt1","Tek","Vwf","Egfl7","Emcn",
      "Cldn5","Esam","Eng","Robo4","Klf2","Esm1",
      "Ets1","Klf4","Sox17","Gja5","Sema3g","Bmx",
      "Mmrn1","Cd93","Plvap","Ramp2","Aplnr","Jam3",
      "Flt4","Prox1","Lyve1","Ccl21a",
      "Ackr1","Madcam1","Selp","Fabp4"
    ),
    master_regulators = c()
  ),
  
  Neural = list(
    genes = c(
      "Tubb3","Elavl3","Elavl4","Rbfox3","Map2","Stmn2","Uchl1","Syt1",
      "Snap25","Syp","Dcx","Tbr1","Neurod1","Neurod2","Rbfox1",
      "Nefm","Nefl","Nrg1","Cntn4","Cntnap2","Elavl2",
      "Sox10","Cdh19","Gfra1","Ret","Phox2b",
      "Plp1","S100b","Fabp7","Gfap","Mpz","Mbp",
      "Vip","Scn3a"
    ),
    master_regulators = c()
  ),
  
  Germline = list(
    genes = c(
      "Ddx4","Dazl","Nanos3","Piwil2","Piwil1",
      "Dppa3","Nanos2","Tfcp2l1",
      "Sycp1","Sycp2","Sycp3","Syce1","Spo11","Stra8",
      "Mael","Mov10l1","Tdrd1","Tdrd5","Tex14",
      "Zp3"
    ),
    master_regulators = c()
  )
)

# ---- HELPER 5: classify_gene() ----

# Tags a gene symbol into a junk class (mito/ribosomal/predicted/pseudogene/
# housekeeping) or "informative". Used both for the QC-junk diagnostic and to
# exclude junk classes before specificity filtering in score_cluster_lineages().
classify_gene <- function(gene) {
  dplyr::case_when(
    stringr::str_detect(gene, "^[Mm][Tt]-")                     ~ "mito",
    stringr::str_detect(gene, "^Rp[sl][0-9]")                    ~ "ribosomal",
    stringr::str_detect(gene, "^Rp[sl]-")                        ~ "ribosomal",
    stringr::str_detect(gene, "^Gm[0-9]+$")                      ~ "predicted",
    stringr::str_detect(gene, "Rik$")                            ~ "predicted",
    stringr::str_detect(gene, "-ps[0-9]*$")                      ~ "pseudogene",
    stringr::str_detect(gene, "^(Actb|Actg1|Tpt1|Eef1a1|Ptma|Fth1|Ftl1|Tmsb4x|Hsp90ab1|Calm1|Myl6|Cox[0-9a-z]*|Ndufa[0-9]*|Ndufb[0-9]*|Atp5[a-z0-9]*)$") ~ "housekeeping",
    TRUE                                                         ~ "informative"
  )
}

# ---- HELPER 5b: get_pct_expressed_matrix() ----

# Computes, for EVERY gene in the assay (not just FindAllMarkers-significant
# ones) and EVERY cluster at one resolution, the fraction of cells in that
# cluster with nonzero counts -- the same quantity FindAllMarkers calls
# pct.1, just unfiltered by significance. This is required for tau (below):
# tau computed from FindAllMarkers output alone is unreliable, because
# genes missing from a cluster's marker rows (not significant there) get
# silently treated as absent, which manufactures artificial specificity for
# genes that are actually broadly present but only "significant" in a few
# clusters (validated on real data: this exact bug initially made Muc2/Fcgbp
# look "highly specific" instead of ambient).
#
# Genes with zero detection in EVERY cluster are dropped before returning --
# tau is undefined/degenerate for them (no real expression pattern exists to
# measure "flatness" of), not because they're safe to treat as ambient.
#
# Run once per resolution (cluster_col changes each call).

get_pct_expressed_matrix <- function(seurat_obj, cluster_col, assay = "RNA", layer = "counts") {
  
  mat <- SeuratObject::GetAssayData(seurat_obj, assay = assay, layer = layer)
  clusters <- as.character(seurat_obj@meta.data[[cluster_col]])
  
  cluster_levels <- sort(unique(clusters))
  pct_mat <- vapply(cluster_levels, function(cl) {
    idx <- which(clusters == cl)
    Matrix::rowSums(mat[, idx, drop = FALSE] > 0) / length(idx)
  }, numeric(nrow(mat)))
  
  colnames(pct_mat) <- cluster_levels
  rownames(pct_mat) <- rownames(mat)
  
  keep <- Matrix::rowSums(pct_mat) > 0
  pct_mat <- pct_mat[keep, , drop = FALSE]
  
  as.data.frame(pct_mat) %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(-gene, names_to = "cluster", values_to = "pct.1")
}


# ---- HELPER 5c: compute_tau_specificity() ----

# Tau tissue-specificity index (Yanai et al. 2005), repurposed here with
# clusters standing in for "tissues". tau = sum(1 - x_hat_i) / (n-1), where
# x_hat_i = x_i / max(x) and x = pct.1 across clusters at one resolution.
#   tau near 0 -> gene detected roughly equally everywhere (ambient candidate)
#   tau near 1 -> gene is essentially exclusive to one cluster (real marker)
# Literature-standard cutoffs: tau > 0.8 = "specific", tau < 0.2 = "broadly
# expressed" -- used as-is below, not re-derived per dataset.
#
# VALIDATED (real project data, 6/6 lineage subsets + Whole, all 6
# resolutions): with get_pct_expressed_matrix()'s zero-everywhere genes
# already dropped, ambient_candidate (tau < 0.2) consistently isolates a
# small (~15-20 gene), biologically coherent set: Malat1, mt-* genes,
# Eef1a1/Tpt1 (independently rediscovering the classify_gene() housekeeping
# list from data alone), plus tissue-specific ambient leakage genes
# (Muc2/Fcgbp/Clca1/etc. in gut-adjacent tissue). No extra magnitude
# threshold (e.g. mean_pct1 > X) is needed on top of tau < 0.2, once the
# zero-everywhere genes are excluded upstream -- the remaining low-tau genes
# were consistently high-magnitude (mean_pct1 mostly > 0.75) with no fuzzy
# in-between cases.
#
# CAVEAT: a gene central to the very lineage being scored (e.g. Epcam within
# an all-Epithelial subset, Pecam1/Cdh5 within an all-Endothelial subset)
# can ALSO show low tau -- not because it's ambient, but because it's
# uniformly high across that lineage's OWN subclusters by definition, so it
# doesn't "discriminate" within this restricted comparison set. Confirmed on
# real data in both directions (Epcam flagged within Epithelial; Pecam1/Cdh5
# NOT flagged when scored against Neural/Mesenchymal/Hematopoietic, where
# they're genuinely rare). This is why ambient exclusion in
# score_cluster_lineages() below always exempts any gene appearing in ANY
# lineage panel, regardless of which lineage is being verified/classified.

compute_tau_specificity <- function(markers_df, value_col = "pct.1") {
  
  all_clusters <- sort(unique(markers_df$cluster))
  n_clusters   <- length(all_clusters)
  if (n_clusters < 2) {
    return(tibble::tibble(gene = character(), tau = numeric(),
                          n_clusters_present = integer(), n_clusters_total = integer(),
                          max_pct1 = numeric(), cluster_of_max = character(),
                          mean_pct1 = numeric(), ambient_candidate = logical(),
                          highly_specific = logical()))
  }
  
  wide <- markers_df %>%
    dplyr::select(gene, cluster, value = dplyr::all_of(value_col)) %>%
    dplyr::group_by(gene, cluster) %>%
    dplyr::summarise(value = max(value, na.rm = TRUE), .groups = "drop") %>%
    tidyr::complete(gene, cluster = all_clusters, fill = list(value = 0)) %>%
    tidyr::pivot_wider(names_from = cluster, values_from = value)
  
  cluster_cols <- setdiff(colnames(wide), "gene")
  m <- as.matrix(wide[, cluster_cols])
  m[is.na(m)] <- 0
  
  max_per_gene <- apply(m, 1, max)
  max_per_gene[max_per_gene == 0] <- NA   # guard only; get_pct_expressed_matrix()
  # should already have dropped these
  
  x_hat <- m / max_per_gene
  tau_val <- rowSums(1 - x_hat, na.rm = TRUE) / (n_clusters - 1)
  cluster_of_max <- cluster_cols[apply(m, 1, which.max)]
  n_present <- rowSums(m > 0)
  
  tibble::tibble(
    gene               = wide$gene,
    tau                = round(tau_val, 4),
    n_clusters_present = n_present,
    n_clusters_total   = n_clusters,
    max_pct1           = round(max_per_gene, 4),
    cluster_of_max      = cluster_of_max,
    mean_pct1          = round(rowMeans(m), 4),
    ambient_candidate   = tau_val < 0.2,
    highly_specific     = tau_val > 0.8
  ) %>%
    dplyr::arrange(tau)
}

# ---- HELPER 6: score_cluster_lineages() ----

# Shared scoring core, called once per resolution by auto_annotate_lineage().
# Not specific to "Whole" (classify) vs "*_Initial" (verify) mode — it simply
# scores every cluster against ALL 6 lineage panels simultaneously and returns
# a tidy long table; the caller decides how to turn that into a verdict.
#
# SECTION 1 (QC diagnostic): mirrors the ORIGINAL naive ranking (pct.1 desc,
#   avg_log2FC desc), no gene exclusion — asks what fraction of the top
#   qc_top_n genes by that ranking are junk/housekeeping classes.
# SECTION 2 (specificity filter): excludes junk gene classes AND tau-ambient
#   genes (panel-exempted -- see active_panel_genes above), then filters to
#   genes that are both prevalent (pct.1 > pct1_floor) and specific
#   (ratio > ratio_thresh). No top-N truncation — a cap here would bury real,
#   high-prevalence markers under rare high-ratio noise before panel-matching
#   ever sees them.
# SECTION 3 (lineage scoring): for each of the 6 lineages, count how many of
#   its panel genes survived Section 2 for each cluster (n_hits), their mean
#   prevalence (avg_pct1), total evidence mass (sum of pct.1), whether a
#   master-regulator gene is among the hits, and which genes those were.
#
# ambient_genes: character vector of gene symbols flagged tau-ambient for
#   THIS resolution/subset (from compute_tau_specificity()'s ambient_candidate
#   column). Pass character(0) to disable ambient filtering entirely (e.g.
#   if the caller didn't/couldn't compute tau -- keeps this function usable
#   standalone against just a markers xlsx, same as before).
#
# Returns: list(qc = <one row per cluster>, scores = <one row per cluster
#   per lineage — NOT collapsed to a single winner/verdict>).

score_cluster_lineages <- function(markers_df, ambient_genes = character(0),
                                   active_panel_genes = character(0),
                                   ratio_thresh = 1.1, pct1_floor = 0.1,
                                   qc_top_n = 50, qc_junk_thresh = 0.5) {
  
  # ---- Section 1: QC diagnostic (naive top-N junk fraction) ----------------
  qc <- markers_df %>%
    dplyr::mutate(gene_class = classify_gene(gene)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(dplyr::desc(pct.1), dplyr::desc(avg_log2FC), .by_group = TRUE) %>%
    dplyr::slice_head(n = qc_top_n) %>%
    dplyr::summarise(
      n_top_checked   = dplyr::n(),
      n_junk          = sum(gene_class %in% c("mito","ribosomal","predicted","pseudogene","housekeeping")),
      frac_junk       = n_junk / n_top_checked,
      qc_flag         = frac_junk > qc_junk_thresh,
      top_genes_naive = paste(gene, collapse = ", "),
      .groups = "drop"
    )
  
  # ---- Section 2: specificity-filtered candidate markers --------------------
  # ambient exclusion: drop tau-flagged genes UNLESS they belong to any
  # lineage panel (pan-lineage identity genes are expected to look "flat"
  # within their own lineage's subclusters -- that's not ambient contamination).
  ambient_to_drop <- setdiff(ambient_genes, active_panel_genes)
  
  candidates <- markers_df %>%
    dplyr::mutate(gene_class = classify_gene(gene)) %>%
    dplyr::filter(gene_class == "informative") %>%
    dplyr::filter(!gene %in% ambient_to_drop) %>%
    dplyr::filter(ratio > ratio_thresh, pct.1 > pct1_floor) %>%
    dplyr::ungroup()
  
  # ---- Section 3: score every cluster against every lineage -----------------
  scores <- purrr::map_dfr(names(lineage_panels), function(lin) {
    panel <- lineage_panels[[lin]]
    candidates %>%
      dplyr::filter(gene %in% panel$genes) %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(
        lineage       = lin,
        n_hits        = dplyr::n(),
        avg_pct1      = mean(pct.1),
        evidence_mass = sum(pct.1),
        has_master    = any(gene %in% panel$master_regulators),
        hit_genes     = paste(gene, collapse = ", "),
        .groups = "drop"
      )
  })
  
  # Every cluster gets a row for every lineage, even 0-hit ones, so downstream
  # argmax/verify logic never has to special-case a missing lineage row.
  all_clusters <- unique(markers_df$cluster)
  full_grid <- tidyr::expand_grid(cluster = all_clusters, lineage = names(lineage_panels))
  scores <- full_grid %>%
    dplyr::left_join(scores, by = c("cluster", "lineage")) %>%
    dplyr::mutate(
      n_hits        = dplyr::coalesce(n_hits, 0L),
      avg_pct1      = dplyr::if_else(n_hits == 0, NA_real_, avg_pct1),
      evidence_mass = dplyr::coalesce(evidence_mass, 0),
      has_master    = dplyr::coalesce(has_master, FALSE),
      hit_genes     = dplyr::coalesce(hit_genes, "")
    )
  
  list(qc = qc, scores = scores)
}


# ---- HELPER 7: auto_annotate_lineage() ----

# Auto-annotation orchestrator. Runs ONLY for label == "Whole" or labels
# matching "*_Initial" — never for "*_Final" (no auto-annotation there; Final
# only re-clusters/re-marks/re-templates on the already-Exclude-filtered set).
#
# Consumes objects already in memory from the current run — no re-reading of
# any xlsx from disk, no separate script:
#   - markers_list        : named list, one sig-markers df per resolution
#                            (already produced by find_markers() this pass)
#   - seurat_obj           : used only to pull @meta.data[, cluster_cols],
#                            i.e. a barcode x resolution cluster-ID table —
#                            never passed into any scoring function itself.
#
# MODE (derived from `label`, same convention as generate_annotation_template()):
#   label == "Whole"        -> classify mode: argmax lineage across all 6
#                               panels per cluster -> Lineage column (one of
#                               LINEAGE_LEVELS).
#   label matches "*_Initial" -> verify mode: expected_lineage (parsed as the
#                               prefix of `label` before "_") vs best
#                               contaminant -> CellType column, one of
#                               Retain / Exclude / Uncertain / Mixed.
#
# For BOTH modes:
#   - Score every resolution once (score_cluster_lineages()) — all 6 stay in
#     memory simultaneously (mirrors find_markers()'s markers_list pattern).
#   - The PRIMARY resolution (primary_resolution, default "harmony_res1.2")
#     alone decides the written call — this is the only thing that ever
#     drives CellType/Lineage downstream (SCRNA_ANNOTATE_CLUSTERS etc.).
#     The other 5 resolutions are a confidence signal only, never a second
#     decision-maker.
#   - call_distribution / n_res / disagreeing_resolutions: for every barcode,
#     look up its cluster's call at EACH of the 6 resolutions (using the
#     same verdict rule as the primary resolution, mode-appropriate), then
#     roll up to the cluster (at primary resolution) by averaging over its
#     member barcodes. n_res is always 6. Bucket order is fixed and spelled
#     out in the header comment above each mode's writer below.
#
# Two files written per call, both under output_dir/Annotation_Templates/:
#   Annotation_{primary_resolution}_{label}_auto.xlsx  — same shape as the
#     manual Annotation_{col}_{label}.xlsx template (cluster_id, n_cells,
#     pct_of_total, top_markers, and the target column pre-filled), so it
#     can be diffed against / used directly as the manual template.
#   Annotation_{primary_resolution}_{label}_reasons.xlsx — one row per
#     cluster: winning/expected call + its hit genes, runner-up/contaminant
#     + its hit genes, margin, QC flag, and the cross-resolution columns.

auto_annotate_lineage <- function(seurat_obj,
                                  label,
                                  cluster_cols,
                                  markers_list,
                                  primary_resolution = "harmony_res1.2",
                                  ratio_thresh = 1.1, pct1_floor = 0.1,
                                  hit_bar = 3, hit_margin_min = 2, avg_margin_min = 0.10,
                                  qc_top_n = 50, qc_junk_thresh = 0.5,
                                  output_dir = NULL,
                                  step_label = "AUTO_ANNOTATE") {
  
  is_whole <- identical(label, "Whole")
  is_initial <- grepl("_Initial$", label)
  
  if (!is_whole && !is_initial) {
    log_info(sample = "", step = step_label,
             msg = glue::glue("label '{label}' is neither 'Whole' nor '*_Initial' — skipping auto-annotation."))
    return(invisible(NULL))
  }
  if (!primary_resolution %in% cluster_cols) {
    log_error(sample = "", step = step_label,
              msg = glue::glue("primary_resolution '{primary_resolution}' not in cluster_cols."))
  }
  if (is.null(output_dir)) {
    log_info(sample = "", step = step_label, msg = "No output_dir — skipping auto-annotation.")
    return(invisible(NULL))
  }
  
  expected_lineage <- if (is_initial) sub("_Initial$", "", label) else NA_character_
  target_col <- if (is_whole) "Lineage" else "CellType"

  # ---- active_panel_genes (data) ----
  # Used as the ambient-gene exemption set: a gene never gets dropped by the 
  # tau ambient filter
  if (label == "Whole") {
    active_panel_genes <- unique(unlist(lapply(lineage_panels, function(p) p$genes)))
  } else if (grepl("_Initial", label)) {
    active_panel_genes <- lineage_panels[[expected_lineage]]$genes
  } else {
    active_panel_genes <- NULL
  }


  # ---- Ambient gene detection (tau), once per resolution --------------------
  # Uses the FULL assay expression matrix (get_pct_expressed_matrix()), not
  # markers_list -- tau computed from FindAllMarkers output alone is unreliable
  # (see compute_tau_specificity() header comment). ambient_genes_by_res feeds
  # directly into score_cluster_lineages() below, which exempts any gene in
  # active_panel_genes regardless of tau.
  
  ambient_audit_list <- list()
  ambient_genes_by_res <- list()
  
  for (col in cluster_cols) {
    pct_long <- get_pct_expressed_matrix(seurat_obj, cluster_col = col)
    tau_tbl  <- compute_tau_specificity(pct_long)  # tau_tbl is created here per iteration
    
    audit_df <- tau_tbl %>%
      dplyr::mutate(
        resolution         = col,
        protected_by_panel = if (is.null(active_panel_genes)) {
          FALSE
        } else {
          gene %in% active_panel_genes
        },
        action             = dplyr::case_when(
          ambient_candidate & !protected_by_panel ~ "DROPPED (Ambient)",
          ambient_candidate & protected_by_panel  ~ "RETAINED (Protected by Panel)",
          TRUE                                    ~ "RETAINED (Informative/Specific)"
        )
      )
    
    ambient_audit_list[[col]] <- audit_df
    ambient_genes_by_res[[col]] <- tau_tbl$gene[tau_tbl$ambient_candidate]
  }
  
  # ── Write ONE Master Audit File along with Annotation templates ──
  if (!is.null(output_dir)) {
    audit_dir <- file.path(output_dir, "Annotation_Templates")
    if (!dir.exists(audit_dir)) dir.create(audit_dir, recursive = TRUE)
    
    master_audit_df <- dplyr::bind_rows(ambient_audit_list) %>%
      dplyr::arrange(gene, tau)
    
    # Simple deduplicated list of unique ambient candidates across any resolution
    simple_ambient_list <- master_audit_df %>%
      dplyr::filter(ambient_candidate) %>%
      dplyr::distinct(gene, .keep_all = TRUE) %>%
      dplyr::select(gene, tau, resolution, protected_by_panel, action) %>%
      dplyr::arrange(gene)
    
    master_file <- file.path(audit_dir, glue::glue("Master_Ambient_Audit_{label}.xlsx"))
    wb_master <- openxlsx::createWorkbook()
    
    openxlsx::addWorksheet(wb_master, "Ambient_Genes_List")
    openxlsx::writeData(wb_master, "Ambient_Genes_List", simple_ambient_list)
    
    openxlsx::addWorksheet(wb_master, "All_Resolutions_Detail")
    openxlsx::writeData(wb_master, "All_Resolutions_Detail", master_audit_df)
    
    openxlsx::saveWorkbook(wb_master, file = master_file, overwrite = TRUE)
    
    log_info(sample = "", step = step_label,
             msg = glue::glue("Saved ambient audit: '{master_file}'"))
  }
  
  # ---- Score every resolution once, keep all 6 in memory --------------------
  scored_by_res <- purrr::map(cluster_cols, function(col) {
    score_cluster_lineages(markers_list[[col]],
                           ambient_genes = ambient_genes_by_res[[col]],
                           active_panel_genes = active_panel_genes,
                           ratio_thresh = ratio_thresh, pct1_floor = pct1_floor,
                           qc_top_n = qc_top_n, qc_junk_thresh = qc_junk_thresh)
  })
  names(scored_by_res) <- cluster_cols
  
  # ---- Per-resolution verdict function (mode-aware), shared by primary calc
  # and by the cross-resolution comparison below, so both use identical logic.
  verdict_for_resolution <- function(res_col) {
    
    sc <- scored_by_res[[res_col]]
    
    if (is_whole) {
      # classify mode: argmax lineage per cluster
      out <- sc$scores %>%
        dplyr::group_by(cluster) %>%
        dplyr::arrange(dplyr::desc(n_hits), dplyr::desc(evidence_mass), .by_group = TRUE) %>%
        dplyr::mutate(rank = dplyr::row_number()) %>%
        dplyr::summarise(
          winner_lineage          = lineage[rank == 1],
          winner_hits             = n_hits[rank == 1],
          winner_avg_pct1         = avg_pct1[rank == 1],
          winner_evidence_mass    = evidence_mass[rank == 1],
          winner_genes            = hit_genes[rank == 1],
          winner_has_master       = has_master[rank == 1],
          runnerup_lineage        = if (dplyr::n() >= 2) lineage[rank == 2] else NA_character_,
          runnerup_hits           = if (dplyr::n() >= 2) n_hits[rank == 2] else 0L,
          runnerup_avg_pct1       = if (dplyr::n() >= 2) avg_pct1[rank == 2] else NA_real_,
          runnerup_evidence_mass  = if (dplyr::n() >= 2) evidence_mass[rank == 2] else 0,
          runnerup_genes          = if (dplyr::n() >= 2) hit_genes[rank == 2] else "",
          margin                  = winner_hits - runnerup_hits,
          .groups = "drop"
        ) %>%
        dplyr::left_join(sc$qc, by = "cluster") %>%
        dplyr::mutate(
          call = dplyr::case_when(
            winner_has_master                              ~ winner_lineage,
            margin >= hit_margin_min & winner_hits > 0       ~ winner_lineage,
            winner_hits > 0 & winner_hits >= hit_bar         ~ winner_lineage,
            TRUE                                            ~ "Uncertain"
          ),
          call = dplyr::if_else(dplyr::coalesce(qc_flag, FALSE), "Uncertain", call)
        )
      return(out)
    }
    
    # verify mode: expected_lineage vs best competing lineage
    exp_rows   <- sc$scores %>% dplyr::filter(lineage == expected_lineage)
    other_best <- sc$scores %>%
      dplyr::filter(lineage != expected_lineage) %>%
      dplyr::group_by(cluster) %>%
      dplyr::arrange(dplyr::desc(n_hits), .by_group = TRUE) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::ungroup()
    
    out <- exp_rows %>%
      dplyr::select(cluster, exp_hits = n_hits, exp_avg_pct1 = avg_pct1,
                    exp_mass = evidence_mass, exp_master = has_master, exp_genes = hit_genes) %>%
      dplyr::left_join(
        other_best %>% dplyr::select(cluster, contaminant_lineage = lineage,
                                     contam_hits = n_hits, contam_avg_pct1 = avg_pct1,
                                     contam_mass = evidence_mass, contam_master = has_master,
                                     contam_genes = hit_genes),
        by = "cluster"
      ) %>%
      dplyr::left_join(sc$qc, by = "cluster") %>%
      dplyr::mutate(
        hit_diff      = exp_hits - contam_hits,
        hits_decisive = abs(hit_diff) >= hit_margin_min,
        hits_dir      = dplyr::case_when(hit_diff > 0 ~ "exp", hit_diff < 0 ~ "contam", TRUE ~ "tie"),
        avg_diff      = exp_avg_pct1 - contam_avg_pct1,
        avg_decisive  = !is.na(avg_diff) & abs(avg_diff) >= avg_margin_min,
        avg_dir       = dplyr::case_when(is.na(avg_diff) ~ "tie", avg_diff > 0 ~ "exp",
                                         avg_diff < 0 ~ "contam", TRUE ~ "tie"),
        call = dplyr::case_when(
          exp_master                                                                  ~ "Retain",
          contam_master                                                               ~ "Exclude",
          hits_decisive & avg_decisive & hits_dir == avg_dir & hits_dir == "exp"       ~ "Retain",
          hits_decisive & avg_decisive & hits_dir == avg_dir & hits_dir == "contam"    ~ "Exclude",
          hits_decisive & avg_decisive & hits_dir != avg_dir                          ~ "Mixed",
          hits_decisive & !avg_decisive & hits_dir == "exp"                           ~ "Retain",
          hits_decisive & !avg_decisive & hits_dir == "contam"                        ~ "Exclude",
          !hits_decisive & avg_decisive & avg_dir == "exp" & exp_hits >= 2             ~ "Retain",
          !hits_decisive & avg_decisive & avg_dir == "contam" & contam_hits >= 2       ~ "Exclude",
          !hits_decisive & !avg_decisive & (exp_hits >= hit_bar | contam_hits >= hit_bar) ~ "Mixed",
          TRUE                                                                        ~ "Uncertain"
        ),
        call = dplyr::if_else(dplyr::coalesce(qc_flag, FALSE), "Uncertain", call)
      )
    out
  }
  
  verdicts_by_res <- purrr::map(cluster_cols, verdict_for_resolution)
  names(verdicts_by_res) <- cluster_cols
  primary <- verdicts_by_res[[primary_resolution]]
  
  # ---- Cross-resolution call_distribution / n_res / disagreeing_resolutions -
  # Bucket order fixed here so counts always line up with the header comment
  # in the written xlsx:
  #   classify mode buckets = LINEAGE_LEVELS, in that order (6-wide)
  #   verify mode buckets   = Retain, Exclude, Uncertain, Mixed (4-wide)
  bucket_levels <- if (is_whole) LINEAGE_LEVELS else c("Retain", "Exclude", "Uncertain", "Mixed")
  
  bc_map <- seurat_obj@meta.data[, cluster_cols, drop = FALSE] %>%
    tibble::rownames_to_column("barcode")
  
  # barcode -> call at each resolution, via that resolution's cluster id
  bc_calls <- bc_map
  for (col in cluster_cols) {
    lut <- stats::setNames(verdicts_by_res[[col]]$call, verdicts_by_res[[col]]$cluster)
    bc_calls[[paste0("call_", col)]] <- unname(lut[as.character(bc_map[[col]])])
  }
  
  # For every barcode: does each resolution's call for its cluster match the
  # PRIMARY resolution's call for that same barcode's primary cluster? Rolled
  # up to primary-resolution clusters below by averaging over member barcodes.
  primary_call_lut <- stats::setNames(primary$call, primary$cluster)
  bc_calls$primary_call <- unname(primary_call_lut[as.character(bc_calls[[primary_resolution]])])
  
  n_res_total <- length(cluster_cols)
  
  cross_res <- purrr::map_dfr(unique(bc_calls[[primary_resolution]]), function(cl) {
    rows <- bc_calls[bc_calls[[primary_resolution]] == cl, , drop = FALSE]
    n_bc <- nrow(rows)
    
    # Mean count-per-resolution landing in each bucket, across this cluster's
    # barcodes -- sums to n_res_total by construction (each barcode contributes
    # exactly one call per resolution, every resolution counted).
    bucket_counts <- vapply(bucket_levels, function(b) {
      per_res_frac <- vapply(paste0("call_", cluster_cols), function(cc) {
        mean(rows[[cc]] == b, na.rm = TRUE)
      }, numeric(1))
      sum(per_res_frac)
    }, numeric(1))
    dist_str <- paste(round(bucket_counts), collapse = "/")
    
    # A resolution "disagrees" if, on average, that resolution's call for
    # this cluster's barcodes doesn't match the primary resolution's call.
    disagree_frac <- vapply(cluster_cols, function(cc) {
      mean(rows[[paste0("call_", cc)]] != rows$primary_call, na.rm = TRUE)
    }, numeric(1))
    disagreeing <- paste(cluster_cols[disagree_frac > 0.5], collapse = ", ")
    
    tibble::tibble(cluster = cl, call_distribution = dist_str,
                   n_res = n_res_total, disagreeing_resolutions = disagreeing)
  })
  
  # ---- Build the two output tables ------------------------------------------
  cluster_vals <- seurat_obj@meta.data[[primary_resolution]]
  n_cells_tbl  <- table(cluster_vals)
  
  top_markers_df <- markers_list[[primary_resolution]] %>%
    dplyr::filter(!grepl("^[Mm][Tt]-", gene)) %>%
    dplyr::arrange(dplyr::desc(pct.1), dplyr::desc(avg_log2FC)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_head(n = 100) %>%
    dplyr::summarise(top_markers = paste(gene, collapse = ", "), .groups = "drop")
  
  if (is_whole) {
    call_df <- primary %>% dplyr::transmute(cluster, !!target_col := call)
    reasons <- primary %>%
      dplyr::left_join(cross_res, by = "cluster") %>%
      dplyr::transmute(
        cluster, Lineage = call,
        winner_lineage, winner_hits, winner_avg_pct1, winner_evidence_mass, winner_genes,
        runnerup_lineage, runnerup_hits, runnerup_avg_pct1, runnerup_evidence_mass, runnerup_genes,
        margin, qc_flag, frac_junk, call_distribution, n_res, disagreeing_resolutions
      )
  } else {
    call_df <- primary %>% dplyr::transmute(cluster, !!target_col := call)
    reasons <- primary %>%
      dplyr::left_join(cross_res, by = "cluster") %>%
      dplyr::transmute(
        cluster, CellType = call,
        expected_lineage = expected_lineage,
        n_expected_hits = exp_hits, expected_avg_pct1 = exp_avg_pct1,
        expected_evidence_mass = exp_mass, expected_hit_genes = exp_genes,
        contaminant_lineage, n_contaminant_hits = contam_hits, contaminant_avg_pct1 = contam_avg_pct1,
        contaminant_evidence_mass = contam_mass, contaminant_genes = contam_genes,
        qc_flag, frac_junk, call_distribution, n_res, disagreeing_resolutions
      )
  }
  
  sheet1 <- tibble::tibble(
    cluster_id   = names(n_cells_tbl),
    n_cells      = as.integer(n_cells_tbl),
    pct_of_total = round(as.integer(n_cells_tbl) / sum(n_cells_tbl) * 100, 1)
  ) %>%
    dplyr::left_join(top_markers_df, by = c("cluster_id" = "cluster")) %>%
    dplyr::left_join(call_df, by = c("cluster_id" = "cluster")) %>%
    dplyr::relocate(dplyr::all_of(target_col), .after = cluster_id) %>%
    dplyr::arrange(as.numeric(as.character(cluster_id)))
  
  template_dir <- file.path(output_dir, "Annotation_Templates")
  if (!dir.exists(template_dir)) dir.create(template_dir, recursive = TRUE)
  
  auto_file <- file.path(template_dir, glue::glue("Annotation_{primary_resolution}_{label}_auto.xlsx"))
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Annotation")
  openxlsx::writeData(wb, "Annotation", sheet1)
  openxlsx::saveWorkbook(wb, file = auto_file, overwrite = TRUE)
  
  reasons_file <- file.path(template_dir, glue::glue("Annotation_{primary_resolution}_{label}_reasons.xlsx"))
  wb2 <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb2, "Reasons")
  openxlsx::writeData(wb2, "Reasons", reasons)
  # Bucket order key, so call_distribution counts are decodable without
  # re-reading this function's source.
  openxlsx::writeData(wb2, "Reasons",
                      glue::glue("call_distribution bucket order: {paste(bucket_levels, collapse=' / ')}"),
                      startRow = nrow(reasons) + 3)
  openxlsx::saveWorkbook(wb2, file = reasons_file, overwrite = TRUE)
  
  log_info(sample = "", step = step_label,
           msg = glue::glue("Saved: '{auto_file}' and '{reasons_file}' — {nrow(sheet1)} clusters, ",
                            "mode={if (is_whole) 'classify(Whole)' else glue::glue('verify(Initial, expected={expected_lineage})')}."))
  
  invisible(list(auto = sheet1, reasons = reasons))
}

# ---- MAIN FUNCTION: cluster_and_findmarkers() ----

# Flow:
#   - Non-"Whole" labels: unchanged single pass, exactly as before — cluster
#     all cells once, find markers once, flag junk/sparse via AutoQC_Flag/
#     CellType, save everyone (nothing physically removed).
#   - "Whole" label: iterative passes. Each pass clusters the CURRENT
#     surviving cells, finds markers, and runs identify_junk_and_sparse().
#     If (n_markerless + n_sparse) / ORIGINAL_N > junk_threshold_pct, those
#     cells are physically subset out, their barcodes+reasons are folded
#     into a running ledger, and a new pass starts on the reduced object —
#     nothing from a non-final pass (clustering columns, markers, plots) is
#     kept. Iteration stops when a pass's junk fraction drops at/below
#     threshold, or max_passes is hit, whichever comes first. Only the
#     FINAL pass's clustering/markers/plots/templates are saved; the saved
#     object also carries @misc$excluded_barcodes recording every cell
#     dropped across all prior passes (reason = "AutoQC_Markerless",
#     "AutoQC_Sparse", or "AutoQC_Markerless+Sparse").
#
#   junk_threshold_pct is evaluated against the ORIGINAL cell count (fixed
#   once, before pass 1) — not the shrinking current-pass count — since a
#   fixed absolute bar converges faster and doesn't get harder to satisfy
#   as cells are removed (see project discussion).
#
#   Setting junk_threshold_pct = 100 (or any value >= 100) disables
#   iteration entirely — junk/sparse fraction can never exceed 100%, so the
#   loop always stops after pass 1, reproducing the old single-pass-only
#   behavior (nothing physically removed here; AutoQC_Flag/CellType
#   pre-filled, physical removal deferred to SCRNA_ADD_CELLTYPE) even for
#   label == "Whole".
#
# KNOWN TRADE-OFF (non-"Whole" labels, and "Whole" when iteration is off):
# junk cells are flagged rather than removed-and-reclustered, so cluster
# boundaries for retained cells were computed while junk was still
# physically present in the graph. See prior discussion — deliberate
# trade-off for speed/simplicity on branches that are typically much
# cleaner than Whole to begin with.

cluster_and_findmarkers <- function(seurat_rds,
                                    label,
                                    resolutions           = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2),
                                    sparse_threshold       = 5,
                                    #junk_threshold_pct     = 2,     # % of ORIGINAL cells; >=100 disables iteration
                                    #max_passes             = 3,     # safety cap regardless of convergence
                                    celltype_marker_xlsx   = NULL,
                                    output_dir            = NULL) {
  
  set.seed(1234)
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Load and Validate
  # ═══════════════════════════════════════════════════════════════════════════
  junk_threshold_pct     <- 2
  max_passes             <- 3
  seurat_obj <- load_smart(input_path = seurat_rds)
  original_n <- ncol(seurat_obj)   # fixed denominator for junk_threshold_pct across all passes
  
  log_info(sample = "", step = "CLUSTER_AND_FINDMARKERS",
           msg = glue::glue("Loaded: {original_n} cells across ",
                            "{length(unique(seurat_obj@meta.data$Sample))} samples."))
  
  if (!"integ_harmony" %in% names(seurat_obj@reductions)) {
    log_error(sample = "", step = "CLUSTER_AND_FINDMARKERS",
              msg = "'integ_harmony' reduction missing. Run INTEGRATE first.")
  }
  if (!"RNA" %in% Seurat::Assays(seurat_obj)) {
    log_error(sample = "", step = "CLUSTER_AND_FINDMARKERS",
              msg = "RNA assay missing — required for FindAllMarkers.")
  }
  
  iterate        <-  junk_threshold_pct < 100
  prior_excluded <- seurat_obj@misc$excluded_barcodes   # NULL, or inherited ledger from upstream
  ledger         <- c()   # accumulates THIS function's own drops across passes only
  
  log_info(sample = "", step = "CLUSTER_AND_FINDMARKERS",
           msg = if (iterate) {
             glue::glue("Iterative mode ON for label '{label}': threshold={junk_threshold_pct}% of ",
                        "original {original_n} cells (={round(original_n*junk_threshold_pct/100)} cells), ",
                        "max_passes={max_passes}.")
           } else {
             glue::glue("Single-pass mode for label '{label}' (junk_threshold_pct={junk_threshold_pct}).")
           })
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Cluster / Markers / Junk-Sparse — one pass, possibly looped
  # ═══════════════════════════════════════════════════════════════════════════
  
  clean_obj <- seurat_obj
  pass <- 1L
  
  repeat {
    
    log_info(sample = "", step = "CLUSTER_AND_FINDMARKERS",
             msg = glue::glue("--- Pass {pass}: {ncol(clean_obj)} cells entering clustering ---"))
    
    clust_result <- run_clustering(clean_obj, resolutions = resolutions,
                                   step_label = glue::glue("CLUSTER_P{pass}"))
    clean_obj   <- clust_result$seurat_obj
    cluster_cols <- clust_result$cluster_cols
    
    if (length(SeuratObject::Layers(clean_obj, assay = "RNA")) > 1) {
      log_info(sample = "", step = "CLUSTER_AND_FINDMARKERS",
               msg = glue::glue("Pass {pass}: joining split RNA layers."))
      clean_obj <- SeuratObject::JoinLayers(clean_obj, assay = "RNA")
    }
    
    # Only the FINAL pass's markers get written to disk — intermediate passes
    # are pure discard-after-use, so markers_dir is NULL until we know this
    # pass is the last one. We still need markers_list in-memory every pass
    # to run identify_junk_and_sparse().
    markers_list <- find_markers(clean_obj,
                                 label = label,
                                 cluster_cols = cluster_cols,
                                 output_dir = NULL,
                                 step_label = glue::glue("FINDMARKERS_P{pass}"))
    
    junk_sparse <- identify_junk_and_sparse(clean_obj, cluster_cols = cluster_cols,
                                            markers_list = markers_list,
                                            sparse_threshold = sparse_threshold,
                                            output_dir = NULL,
                                            step_label = glue::glue("JUNKSPARSE_P{pass}"))
    
    n_flagged   <- length(junk_sparse$excluded_barcodes)
    pct_flagged <- 100 * n_flagged / original_n
    
    log_info(sample = "", step = "CLUSTER_AND_FINDMARKERS",
             msg = glue::glue("Pass {pass}: {n_flagged} markerless+sparse cells ",
                              "({round(pct_flagged, 2)}% of original {original_n})."))
    
    # Remove this pass's junk/sparse, ledger it
    dropped_barcodes <- junk_sparse$excluded_barcodes
    reasons <- dplyr::case_when(
      dropped_barcodes %in% intersect(junk_sparse$markerless_barcodes,
                                      junk_sparse$sparse_barcodes) ~ glue::glue("{label}_Markerless+Sparse"),
      dropped_barcodes %in% junk_sparse$markerless_barcodes        ~ glue::glue("{label}_Markerless"),
      TRUE                                                         ~ glue::glue("{label}_Sparse")
    )
    ledger <- c(ledger, stats::setNames(reasons, dropped_barcodes))
    
    keep_cells <- setdiff(colnames(clean_obj), dropped_barcodes)
    clean_obj <- subset(clean_obj, cells = keep_cells)
    
    stop_now <- (!iterate) || (pct_flagged <= junk_threshold_pct) || (pass >= max_passes)
    
    if (stop_now) {
      if (iterate && pct_flagged > junk_threshold_pct && pass >= max_passes) {
        log_warn(sample = "", step = "CLUSTER_AND_FINDMARKERS",
                 msg = glue::glue("Stopping at max_passes={max_passes} — this pass's ",
                                  "{round(pct_flagged,2)}% flagged (threshold={junk_threshold_pct}%) was still ",
                                  "removed before saving; convergence just took longer than max_passes allowed."))
      }
      break
    }
    
    # Strip this pass's clustering artifacts before looping — the whole point
    # of another pass is a clean re-cluster, not carrying stale columns/graphs.
    for (col in cluster_cols) clean_obj@meta.data[[col]] <- NULL
    for (r in c("umap_harmony")) if (r %in% SeuratObject::Reductions(clean_obj)) clean_obj[[r]] <- NULL
    clean_obj@graphs <- list()
    
    log_info(sample = "", step = "CLUSTER_AND_FINDMARKERS",
             msg = glue::glue("Pass {pass}: removed {length(dropped_barcodes)} cells, ",
                              "{ncol(clean_obj)} remain — starting pass {pass + 1}."))
    
    pass <- pass + 1L
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Finalize — save FINAL pass's markers, flag remaining junk/sparse
  # ═══════════════════════════════════════════════════════════════════════════
  # Only the final pass's marker tables get written to disk. Re-running
  # find_markers a second time here (rather than caching pass N's in-memory
  # result) keeps this identical in shape/behavior to the pre-existing
  # single-pass code path — the small re-compute cost buys simplicity and
  # avoids a separate "did we already write these" branch.
  
  markers_dir  <- if (!is.null(output_dir)) file.path(output_dir, "FindAllMarkers") else NULL
  
  markers_list <- find_markers(clean_obj,
                               label = label,
                               cluster_cols = cluster_cols,
                               output_dir = markers_dir,
                               step_label = "FINDMARKERS")
  
  excluded_vec <- c(prior_excluded, ledger)   # named: names = barcodes, values = reason
  seurat_obj$AutoQC_Flag <- ifelse(colnames(seurat_obj) %in% names(excluded_vec),
                                   excluded_vec[colnames(seurat_obj)],
                                   "Retained")
  
  # Fold this function's own multi-pass ledger together with whatever ledger
  # the object already carried in from upstream (SUBSET_CELLS, SCT_TRANSFORM,
  # etc.) — same accumulation pattern used throughout the pipeline.
  clean_obj@misc$excluded_barcodes <- c(prior_excluded, ledger)
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Annotation Template
  # ═══════════════════════════════════════════════════════════════════════════
  
  generate_annotation_template(clean_obj,
                               label                 = label,
                               cluster_cols          = cluster_cols,
                               markers_list          = markers_list,
                               celltype_marker_xlsx  = celltype_marker_xlsx,
                               output_dir            = output_dir,
                               step_label            = "ANNOTATION_TEMPLATE")
  
  # Auto-annotation: runs only for label == "Whole" or "*_Initial" (no-op,
  # early-return for "*_Final" — see auto_annotate_lineage() header comment).
  auto_annotate_lineage(clean_obj,
                        label         = label,
                        cluster_cols  = cluster_cols,
                        markers_list  = markers_list,
                        output_dir    = output_dir,
                        step_label    = "AUTO_ANNOTATE")
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5: Plots
  # ═══════════════════════════════════════════════════════════════════════════
  # Plots reflect only the FINAL pass's surviving cells — earlier passes'
  # physically-removed junk is gone and can't appear here. Everyone still
  # present (including this pass's own flagged-but-not-yet-removed junk, if
  # iteration stopped at threshold/max_passes with some still flagged) stays
  # visible, same rationale as before: template numbering must match, and
  # the AutoQC_Flag audit view needs the flagged cells present to review.
  
  if (!is.null(output_dir)) {
    
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    for (col in cluster_cols) {
      
      plot_seurat(clean_obj,
                  features = col, feature_type = "metadata",
                  reduction = "umap_harmony",
                  filename = glue::glue("UMAP_Clusters_{col}_{label}"),
                  output_dir = output_dir)
      
      plot_seurat(clean_obj,
                  features = col, feature_type = "metadata",
                  reduction = "umap_harmony",
                  filename = glue::glue("UMAP_Clusters_{col}_{label}_split"),
                  split_col = col, output_dir = output_dir)
    }
    
    plot_seurat(clean_obj,
                features = c("nUMIs", "nGenes", "MitoRatio", "S.Score", "G2M.Score", "CC.Score"),
                feature_type = "metadata", reduction = "umap_harmony",
                filename = glue::glue("UMAP_Clustering_QC_Metrics_{label}"),
                output_dir = output_dir)
    
    plot_seurat(clean_obj,
                features = "Sample", feature_type = "metadata",
                reduction = "umap_harmony",
                filename = glue::glue("UMAP_Clustering_Sample_{label}"),
                split_col = "Sample", output_dir = output_dir)
    
    plot_seurat(seurat_obj,
                features = "AutoQC_Flag", feature_type = "metadata",
                reduction = "umap_harmony_qc",
                filename = glue::glue("UMAP_AutoQC_Flag_{label}"),
                output_dir = output_dir)
    
    log_info(sample = "", step = "CLUSTER_AND_FINDMARKERS", msg = "Post-clustering plots saved.")
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 6: Save
  # ═══════════════════════════════════════════════════════════════════════════
  
  if (!is.null(output_dir)) {
    
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    rds_out <- file.path(output_dir, "clustered_seurat.rds")
    saveRDS(object = clean_obj, file = rds_out)
    log_info(sample = "", step = "CLUSTER_AND_FINDMARKERS",
             msg = glue::glue("Saved: '{rds_out}' ({ncol(clean_obj)} cells)."))
    
  } else {
    log_info(sample = "", step = "CLUSTER_AND_FINDMARKERS", msg = "No output_dir — skipping save.")
  }
  
  log_info(sample = "", step = "CLUSTER_AND_FINDMARKERS", msg = "Completed successfully.")
  
  return(invisible(clean_obj))
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
  
  raw_res          <- get_arg(3, default = "0.2,0.4,0.6,0.8,1.0,1.2")
  resolutions_list <- as.numeric(trimws(strsplit(raw_res, ",")[[1]]))
  
  cluster_and_findmarkers(
    seurat_rds            = get_arg(1),
    label                 = get_arg(2),
    resolutions           = resolutions_list,
    sparse_threshold      = as.integer(get_arg(4, default = "5")),
    celltype_marker_xlsx  = get_arg(5, default = NULL),
    output_dir            = get_arg(6, default = ".")
  )
}
