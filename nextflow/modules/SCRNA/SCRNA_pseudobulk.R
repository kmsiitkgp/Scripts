# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)
  library(glue)
  library(Matrix)
  library(Seurat)
  library(SeuratObject)
  library(openxlsx)
  library(readr)
  library(tibble)
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

# ---- 🧬 3. Function Definition (Always Runs) ----

# MIN_TOTAL_SAMPLES — structural floor for RNA_CREATE_DDS's DESeq2::DESeq() call
# to even attempt a model fit on this group. NOT a tunable statistical/scientific
# choice like min_cells_per_sample is — it's a fixed fact about the tool (DESeq()
# fits dispersions on the WHOLE group, once, before any contrast is extracted;
# with <2 total samples that fit can't run at all, independent of which specific
# contrast might later prove viable — that's RNA_extract_deseq2_results.R's own
# per-contrast check, a separate and later concern). Hardcoded on purpose rather
# than exposed as a yaml parameter.
MIN_TOTAL_SAMPLES <- 2

# =============================================================================
# FUNCTION: run_pseudobulk
# =============================================================================
# Purpose: Aggregates raw RNA counts per Sample, for ONE analyses entry (one
#          CellType/Lineage selector), into a pseudobulk matrix ready for
#          RNA_CREATE_DDS.
#
# Design (per project discussion — fan-out-per-entry, not single-task-loops-all):
#   - Nextflow fans `params.scrna.pseudobulk.analyses` out into one task PER
#     ENTRY (scrnaseq.nf: Channel.fromList(pseudobulk_analyses)), mirroring
#     exactly how bulk RNAseq fans contrasts out via Channel.fromList(contrasts).
#     Each task here handles exactly one entry — no looping over multiple
#     entries inside one R run, no JSON/jsonlite needed to hand R a nested
#     structure, since selector_values is just one flat comma-joined string
#     (same convention as batch_vars elsewhere). Nextflow's channels — not a
#     string encoding — are what keep entries separate.
#   - Trade-off, explicit: each task independently reloads seurat_rds. For a
#     large whole-tissue object this means N reloads instead of one shared
#     load — accepted in exchange for incremental downstream pipelining (each
#     group's RNA_CREATE_DDS can start the moment ITS task finishes, not after
#     every entry has finished) and dropping the JSON dependency entirely.
#   - group_name is NOT derived here — it's computed once in scrnaseq.nf (via
#     the same label_for()/sanitize_celltype() used for subcluster labels) and
#     passed straight through as an argument. Avoids having two independent
#     sanitization implementations (Groovy + R) that would need to stay in
#     sync purely to make a downstream directory-name join work.
#   - contrasts are NOT an input to this script at all — this step only ever
#     used them for audit-log bookkeeping, never for aggregation logic itself.
#     That bookkeeping is gone (see below), so there's nothing left in this
#     script that has any use for contrast information.
#   - No audit file. Every reason a sample or an entire entry gets excluded
#     reduces to exactly one criterion — too few cells of this population from
#     that sample (explicitly filtered below min_cells_per_sample, or trivially
#     zero to begin with) — and that's already written to this task's own log
#     via log_info()/log_warn() below. A separate xlsx would only reformat
#     information already in plain text, for a single-entry-per-task script
#     where "no output produced" already IS the skip signal Nextflow needs.
#   - MIN_TOTAL_SAMPLES (hardcoded, see above): if fewer than this many samples
#     survive min_cells_per_sample filtering, this task produces NO files at
#     all and exits 0 — a graceful skip, not an error — because RNA_CREATE_DDS's
#     DESeq() fit can't run on fewer than that regardless of condition/group
#     membership. All 3 outputs are declared optional: true in
#     SCRNA_pseudobulk.nf for exactly this reason.
#   - Sample-level metadata is pulled directly from seurat_obj@meta.data
#     (deduplicated to one row per Sample), NOT re-read from metadata_xlsx —
#     that file was already joined onto the object back in MERGE_AND_PLOT_QC.
#   - Raw count aggregation uses a sparse matrix cross-product (Matrix::fac2sparse
#     + matrix multiply), NOT Seurat::AggregateExpression() — see inline
#     comments below for the full rationale (unchanged from prior design).
# =============================================================================

run_pseudobulk <- function(seurat_rds,
                           selector_col,
                           selector_values,
                           group_name,
                           min_cells_per_sample = 100,
                           output_dir            = ".") {

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Load and Validate
  # ═══════════════════════════════════════════════════════════════════════════

  seurat_obj <- load_smart(input_path = seurat_rds)

  log_info(sample = group_name, step = "pseudobulk",
           msg = glue::glue("Loaded: {ncol(seurat_obj)} cells, {nrow(seurat_obj)} features."))

  for (req_col in c("CellType", "Lineage", "Sample")) {
    if (!req_col %in% colnames(seurat_obj@meta.data)) {
      log_error(sample = group_name, step = "pseudobulk",
                msg = glue::glue("Required metadata column missing: {req_col}"))
    }
  }

  if (!selector_col %in% c("CellType", "Lineage")) {
    log_error(sample = group_name, step = "pseudobulk",
              msg = glue::glue("selector_col must be 'CellType' or 'Lineage', got: {selector_col}"))
  }

  if (!"RNA" %in% SeuratObject::Assays(seurat_obj)) {
    log_error(sample = group_name, step = "pseudobulk",
              msg = "RNA assay missing — required for raw count aggregation.")
  }

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # selector_values arrives as one comma-joined string (same convention as
  # batch_vars elsewhere) — split back into a character vector here.
  selector_values_vec <- trimws(strsplit(selector_values, ",")[[1]])

  log_info(sample = group_name, step = "pseudobulk",
           msg = glue::glue("Selecting {selector_col} %in% {paste(selector_values_vec, collapse = ', ')}"))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Select Cells, Filter Samples, Aggregate
  # ═══════════════════════════════════════════════════════════════════════════

  cells_in_group <- colnames(seurat_obj)[as.character(seurat_obj@meta.data[[selector_col]]) %in% selector_values_vec]

  if (length(cells_in_group) == 0) {
    log_warn(sample = group_name, step = "pseudobulk",
             msg = glue::glue("SKIP '{group_name}': 0 cells match {selector_col} %in% ",
                              "{paste(selector_values_vec, collapse = ', ')}. No output produced."))
    return(invisible(NULL))
  }

  subset_obj <- subset(seurat_obj, cells = cells_in_group)

  # ── Per-sample cell counts (before aggregation) — used for min_cells_per_sample ──
  raw_samples             <- as.character(subset_obj$Sample)
  names(raw_samples)      <- colnames(subset_obj)
  sample_key_map          <- stats::setNames(make.names(unique(raw_samples)), unique(raw_samples))
  cell_counts_per_sample  <- table(sample_key_map[raw_samples])

  surviving_samples <- names(cell_counts_per_sample)[cell_counts_per_sample >= min_cells_per_sample]
  dropped_samples   <- setdiff(names(cell_counts_per_sample), surviving_samples)

  if (length(dropped_samples) > 0) {
    log_info(sample = group_name, step = "pseudobulk",
             msg = glue::glue("'{group_name}': dropping {length(dropped_samples)} sample(s) below ",
                              "min_cells_per_sample ({min_cells_per_sample}): ",
                              "{paste(glue::glue('{dropped_samples}({cell_counts_per_sample[dropped_samples]})'), collapse=', ')}"))
  }

  # ── MIN_TOTAL_SAMPLES structural gate — before any aggregation work ──────
  if (length(surviving_samples) < MIN_TOTAL_SAMPLES) {
    log_warn(sample = group_name, step = "pseudobulk",
             msg = glue::glue("SKIP '{group_name}': only {length(surviving_samples)} sample(s) survived ",
                              "min_cells_per_sample filtering — need >= {MIN_TOTAL_SAMPLES} for RNA_CREATE_DDS's ",
                              "DESeq2::DESeq() fit to run at all. No output produced."))
    return(invisible(NULL))
  }

  # ═══════════════════════════════════════════════════════════════════════
  # Aggregate raw RNA counts per Sample — sparse matrix cross-product
  # ═══════════════════════════════════════════════════════════════════════
  # JoinLayers defensively: Seurat v5 stores counts as separate per-sample
  # layers (e.g. counts.Sample1, counts.Sample2, ...) unless already joined
  # upstream. Called here, on the (already cell-subsetted, so cheaper) group
  # subset, regardless of whether an earlier pipeline step already did this.
  subset_obj <- SeuratObject::JoinLayers(subset_obj, assay = "RNA")

  counts_mat <- SeuratObject::LayerData(subset_obj, assay = "RNA", layer = "counts")

  keep_cells <- colnames(subset_obj)[sample_key_map[raw_samples] %in% surviving_samples]
  counts_mat <- counts_mat[, keep_cells, drop = FALSE]

  # Sparse indicator matrix: one row per surviving sample, one column per
  # kept cell, 1 where that cell belongs to that sample. Column ORDER (not
  # names) must match counts_mat's column order exactly — both built from the
  # same keep_cells vector, in the same order.
  sample_factor       <- factor(sample_key_map[raw_samples[keep_cells]], levels = surviving_samples)
  indicator           <- Matrix::fac2sparse(sample_factor)          # samples x cells, sparse 0/1
  expr_mat_sparse      <- counts_mat %*% Matrix::t(indicator)        # genes x samples, one BLAS matmul
  expr_mat             <- as.data.frame(as.matrix(expr_mat_sparse))
  colnames(expr_mat)   <- surviving_samples

  # ═══════════════════════════════════════════════════════════════════════
  # Build sample-level metadata from seurat_obj@meta.data (NOT metadata_xlsx)
  # ═══════════════════════════════════════════════════════════════════════

  sample_metadata <- subset_obj@meta.data %>%
    dplyr::mutate(Sample_ID = sample_key_map[as.character(Sample)]) %>%
    dplyr::filter(Sample_ID %in% colnames(expr_mat)) %>%
    dplyr::distinct(Sample_ID, .keep_all = TRUE) %>%
    dplyr::select(Sample_ID, dplyr::any_of(setdiff(colnames(subset_obj@meta.data), c(
      "CellType", "Lineage", "CellType_Source", "AutoQC_Flag",
      "Stability_Score", "orig.ident", "Sample"))))

  # ═══════════════════════════════════════════════════════════════════════
  # Save expr_mat.csv, metadata.xlsx, tx2gene_dummy.csv — flat into output_dir.
  # No group_name subdirectory needed here: this task IS already scoped to
  # one group, and SCRNA_pseudobulk.nf's publishDir places these under
  # 05.Pseudobulk/{group_name}/ on the Nextflow side.
  # ═══════════════════════════════════════════════════════════════════════

  gene_ids           <- make.names(rownames(expr_mat))
  rownames(expr_mat) <- gene_ids

  readr::write_csv(expr_mat %>% tibble::rownames_to_column("gene_id"),
                   file.path(output_dir, "expr_mat.csv"))

  openxlsx::write.xlsx(sample_metadata, file.path(output_dir, "metadata.xlsx"), overwrite = TRUE)

  # tx2gene_dummy — identity mapping. Column 1 is a placeholder (never read by
  # add_annotation(), which selects columns 2/3/4 positionally). Column 4
  # (gene_biotype) is never functionally consumed either — filled with
  # "protein_coding" as a harmless default.
  tx2gene_dummy <- data.frame(
    transcript_id = gene_ids,
    gene_id       = gene_ids,
    gene_name     = gene_ids,
    gene_biotype  = "protein_coding"
  )
  readr::write_csv(tx2gene_dummy, file.path(output_dir, "tx2gene_dummy.csv"))

  log_info(sample = group_name, step = "pseudobulk",
           msg = glue::glue("Saved '{group_name}': {ncol(expr_mat)} samples, {nrow(expr_mat)} genes."))

  log_info(sample = group_name, step = "pseudobulk", msg = "Completed successfully.")

  return(invisible(NULL))
}

# ---- 🚀 4. Smart Execution (Nextflow Only) ----
# CLI arg order matches SCRNA_pseudobulk.nf's Rscript call exactly:
#   seurat_rds, selector_col, selector_values, group_name, min_cells_per_sample, output_dir

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)

  get_arg <- function(idx, default = NULL) {
    if (idx > length(args)) return(default)
    val <- args[idx]
    if (is.na(val) || val == "" || val == "null" || val == "NULL") return(default)
    return(val)
  }

  run_pseudobulk(
    seurat_rds             = get_arg(1),
    selector_col           = get_arg(2),
    selector_values        = get_arg(3),
    group_name             = get_arg(4),
    min_cells_per_sample   = as.integer(get_arg(5, default = "100")),
    output_dir             = get_arg(6, default = ".")
  )
}
