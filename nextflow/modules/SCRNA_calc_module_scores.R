# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)
  library(glue)
  library(openxlsx)
  library(Seurat)
  library(SeuratObject)
  library(UCell)
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

calc_module_scores <- function(rds_file, marker_file, output_dir = NULL) {

  # Set seed for reproducible stochastic processes
  # AddModuleScore uses random background gene sampling — seed ensures
  # identical scores on reruns, critical for pipeline reproducibility.
  set.seed(1234)

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Load Seurat Object
  # ═══════════════════════════════════════════════════════════════════════════
  # Expects clustered_seurat.rds produced by CLUSTER_AND_FINDMARKERS.
  # That object contains:
  #   - SCT assay with normalized expression
  #   - Integration reductions (integ_harmony)
  #   - UMAP reductions (umap_harmony)
  #   - Cluster columns (harmony_res0.2 ... harmony_res1.2)
  #   - @misc$integration_params and @misc$clustering_params

  integrated_seurat <- load_smart(input_path = rds_file)

  log_info(sample = "", step = "calc_module_scores",
           msg = glue("Loaded: {ncol(integrated_seurat)} cells across ",
                      "{length(unique(integrated_seurat@meta.data$Sample))} samples."))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Validate Prerequisites
  # ═══════════════════════════════════════════════════════════════════════════
  # SCT assay must exist — scoring always runs on SCT normalized data.
  # Marker file must exist and be readable.

  if (!"SCT" %in% Seurat::Assays(integrated_seurat)) {
    log_error(sample = "", step = "calc_module_scores",
              msg = "SCT assay missing. Run SCT_TRANSFORM before CALC_MODULE_SCORES.")
  }

  if (!file.exists(marker_file)) {
    log_error(sample = "", step = "calc_module_scores",
              msg = glue("Marker file not found: '{marker_file}'."))
  }

# ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Resolve Parameters from @misc
  # ═══════════════════════════════════════════════════════════════════════════
  # Read clustering parameters written to @misc by CLUSTER_AND_FINDMARKERS.
  # Using stored column names rather than re-deriving from resolutions avoids
  # fragility — if resolutions change between runs the stored names remain
  # accurate. Fall back to derivation only for interactive / standalone use.
  #
  # cluster_cols : exact clustering columns present in metadata
  #               e.g. c("harmony_res0.2", "harmony_res0.4"... "harmony_res1.2")
  #               Stored here so ANNOTATE_CLUSTERS reads them directly from
  #               @misc rather than re-deriving by regex.
  cluster_cols <- integrated_seurat@misc$clustering_params$cluster_cols %||% {
    resolutions <- integrated_seurat@misc$clustering_params$resolutions %||% c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
    as.vector(outer("harmony", resolutions, function(m, r) paste0(m, "_res", r)))
  }
  cluster_cols <- cluster_cols[cluster_cols %in% colnames(integrated_seurat@meta.data)]

  log_info(sample = "", step = "calc_module_scores",
           msg = glue("Cluster columns found: {length(cluster_cols)}."))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Set Active Assay
  # ═══════════════════════════════════════════════════════════════════════════
  # Scoring always runs on SCT assay — this step always follows CLUSTER_AND_FINDMARKERS
  # which itself requires SCT. SCT data layer contains log-normalised, variance-
  # stabilised expression values appropriate for module scoring.
  #
  # NOTE: AddModuleScore and UCell both use the normalised 'data' layer, not
  # raw counts. SCT 'data' layer is more appropriate than RNA 'data' for
  # integrated multi-sample datasets because SCT corrects for sequencing depth
  # differences across samples before integration.

  active_assay <- "SCT"
  Seurat::DefaultAssay(integrated_seurat) <- active_assay

  log_info(sample = "", step = "calc_module_scores",
           msg = glue("Active assay set to: '{active_assay}'."))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5: Build Validated Gene Set List from Marker File
  # ═══════════════════════════════════════════════════════════════════════════
  # Marker file is an Excel file where each column = one cell type.
  # Column headers become module names. Each column contains marker genes
  # for that cell type, with NAs padding shorter columns to equal length.
  #
  # Validation steps per gene set:
  #   1. Remove NAs
  #   2. Case-insensitive match against genes present in the SCT assay
  #      (handles human/mouse capitalisation differences e.g. Hba-a1 vs HBA1)
  #   3. Calculate coverage = detected / expected genes
  #   4. Require >= 3 matched genes (below this AddModuleScore is unreliable)
  #
  # Gene sets with < 3 matched genes are skipped with a warning.
  # Coverage is logged so low-coverage sets can be identified and improved.
  #
  # NOTE: We use make.names() on column headers to ensure module names are
  # valid R variable names — spaces and special characters are replaced with
  # dots. This matters because AddModuleScore appends a numeric suffix to
  # create metadata column names, and invalid names cause errors.

  marker_df <- tryCatch({
    openxlsx::read.xlsx(xlsxFile = marker_file)
  }, error = function(e) {
    log_error(sample = "", step = "calc_module_scores",
              msg = glue("Failed to read marker file: '{e$message}'."))
  })

  # Get all genes present in the SCT assay data layer
  present_features <- rownames(
    SeuratObject::GetAssayData(object   = integrated_seurat,
                               assay    = active_assay,
                               layer    = "data")
  )

  # Build validated signatures list
  signatures_list <- list()

  for (i in base::seq_len(ncol(marker_df))) {

    module_name <- make.names(colnames(marker_df)[i])

    # Extract non-NA genes from this column
    xlsx_features <- marker_df[[i]] %>%
      stats::na.omit() %>%
      as.vector()

    # Case-insensitive matching — handles Hba-a1 (mouse) vs HBA1 (human)
    # present_features retains original capitalisation for correct gene lookup
    features <- present_features[tolower(present_features) %in% tolower(xlsx_features)]

    # Calculate and log coverage
    coverage <- round(length(features) / length(xlsx_features), digits = 2)

    log_info(sample = "", step = "calc_module_scores",
             msg = glue("Module '{module_name}': ",
                        "{length(features)}/{length(xlsx_features)} genes detected ",
                        "(coverage = {coverage})."))

    if (coverage < 0.3) {
      log_warn(sample = "", step = "calc_module_scores",
               msg = glue("Low coverage for '{module_name}' ({coverage}). ",
                          "Consider updating marker list for this dataset."))
    }

    # Require minimum 3 matched genes
    # Below this threshold AddModuleScore background comparison is unreliable
    # and UCell rank statistics have insufficient signal
    if (length(features) >= 3) {
      signatures_list[[module_name]] <- sort(features)
    } else {
      log_warn(sample = "", step = "calc_module_scores",
               msg = glue("Skipping '{module_name}' — fewer than 3 matching markers found. ",
                          "({length(features)} detected)."))
    }
  }

  if (length(signatures_list) < 2) {
    log_error(sample = "", step = "calc_module_scores",
              msg = "Fewer than 2 valid gene sets found. Cannot proceed with scoring.")
  }

  log_info(sample = "", step = "calc_module_scores",
           msg = glue("{length(signatures_list)} valid gene sets ready for scoring."))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 6: Calculate Module Scores — Seurat AddModuleScore
  # ═══════════════════════════════════════════════════════════════════════════
  # AddModuleScore computes per-cell enrichment scores by comparing mean
  # expression of the gene set against randomly sampled control genes drawn
  # from the same expression bins. This controls for expression level bias —
  # highly expressed genes don't inflate scores simply by being highly expressed.
  #
  # Output: one metadata column per gene set, named "{module_name}1"
  # (Seurat appends a numeric index starting at 1 as suffix).
  # Scores can be negative — cells where the gene set is expressed BELOW
  # background level get negative scores.
  #
  # Negative scores are clamped to 0 in ANNOTATE_CLUSTERS because:
  #   - Single cell dropout means absence of expression is unreliable
  #   - A negative score reflects "expressed below background" not
  #     "actively not this cell type" — analogous to flow cytometry
  #     where you annotate on positive expression only, not absence
  #   - For broad orthogonal cell types (T cell vs Epithelial) a cell
  #     that is not a T cell will simply have zero T cell markers detected,
  #     not negative expression — the negative score is noise not signal

  integrated_seurat <- tryCatch({
    Seurat::AddModuleScore(object   = integrated_seurat,
                           features = signatures_list,
                           assay    = active_assay,
                           layer    = "data",
                           name     = names(signatures_list))
  }, error = function(e) {
    log_error(sample = "", step = "calc_module_scores",
              msg = glue("AddModuleScore failed: {e$message}"))
  })

  log_info(sample = "", step = "calc_module_scores",
           msg = "Seurat AddModuleScore completed.")

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 7: Calculate Module Scores — UCell
  # ═══════════════════════════════════════════════════════════════════════════
  # UCell computes per-cell enrichment scores using the Mann-Whitney U
  # statistic on per-cell gene rankings. This is fundamentally different
  # from AddModuleScore:
  #
  #   AddModuleScore : mean expression vs random background
  #                   → sensitive to dataset composition, expression depth,
  #                     gene set size, background gene selection
  #
  #   UCell          : rank-based, no background needed
  #                   → robust to batch effects, sequencing depth variation,
  #                     dataset composition, gene set size
  #                   → always produces scores in [0, 1]
  #                   → particularly important for small gene sets
  #                     (LymphaticEndothelial=5, Myofibroblasts=6 genes)
  #                     where AddModuleScore background is less stable
  #
  # UCell is notably slow on large datasets (20-40 min for 100k cells) —
  # this is the primary reason CALC_MODULE_SCORES is kept as a separate
  # function from ANNOTATE_CLUSTERS. The scored RDS is cached so annotation
  # parameters can be tuned without rerunning UCell.
  #
  # Output: one metadata column per gene set, named "{module_name}_UCell"

  integrated_seurat <- tryCatch({
    UCell::AddModuleScore_UCell(obj      = integrated_seurat,
                                features = signatures_list,
                                assay    = active_assay,
                                slot     = "data",
                                name     = "_UCell")
  }, error = function(e) {
    log_error(sample = "", step = "calc_module_scores",
              msg = glue("UCell scoring failed: {e$message}"))
  })

  log_info(sample = "", step = "calc_module_scores",
           msg = "UCell scoring completed.")

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 8: Identify and Validate Module Score Columns
  # ═══════════════════════════════════════════════════════════════════════════
  # After scoring, we identify the exact metadata columns added by each method.
  # We derive expected names from our signatures_list rather than using regex
  # pattern matching on all metadata columns — this is more robust because:
  #   - Regex on "[0-9]+$" also matches cluster columns (harmony0.4)
  #     requiring additional filtering with bad_prefixes
  #   - Explicit derivation is exact — no risk of capturing wrong columns
  #
  # AddModuleScore naming: appends sequential integers starting at 1
  #   e.g. T.Cell → "T.Cell1", B.Cell → "B.Cell2"
  #   NOTE: the integer suffix is the position index in signatures_list,
  #   not always "1" — hence we use seq_along() not hardcoded "1"
  #
  # UCell naming: appends "_UCell" suffix
  #   e.g. T.Cell → "T.Cell_UCell"

  module_names      <- names(signatures_list)
  seurat_module_cols <- paste0(module_names, seq_along(module_names))
  ucell_module_cols  <- paste0(module_names, "_UCell")

  # Verify columns actually exist in metadata — catch any scoring failures
  seurat_module_cols <- seurat_module_cols[seurat_module_cols %in% colnames(integrated_seurat@meta.data)]
  ucell_module_cols  <- ucell_module_cols[ucell_module_cols  %in% colnames(integrated_seurat@meta.data)]

  if (length(seurat_module_cols) < 2) {
    log_error(sample = "", step = "calc_module_scores",
              msg = glue("Fewer than 2 Seurat module score columns found in metadata. ",
                         "AddModuleScore may have failed silently."))
  }

  if (length(ucell_module_cols) < 2) {
    log_error(sample = "", step = "calc_module_scores",
              msg = glue("Fewer than 2 UCell module score columns found in metadata. ",
                         "UCell scoring may have failed silently."))
  }

  log_info(sample = "", step = "calc_module_scores",
           msg = glue("Identified {length(seurat_module_cols)} Seurat score columns ",
                      "and {length(ucell_module_cols)} UCell score columns."))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 9: Store Scoring Parameters in @misc
  # ═══════════════════════════════════════════════════════════════════════════
  # Stored at object level so ANNOTATE_CLUSTERS reads exact column names
  # directly from @misc rather than re-deriving by regex.
  #
  # seurat_module_cols : exact metadata column names for Seurat scores
  #                      e.g. c("T.Cell1", "B.Cell2", "Epithelial3"...)
  #                      Used by ANNOTATE_CLUSTERS to build score matrix
  #                      Used by plot functions for FeaturePlot
  #
  # ucell_module_cols  : exact metadata column names for UCell scores
  #                      e.g. c("T.Cell_UCell", "B.Cell_UCell"...)
  #
  # signatures_list    : validated gene sets actually used for scoring
  #                      (post case-insensitive matching, post coverage filter)
  #                      Documents exactly which genes drove the scores —
  #                      important for reproducibility and manuscript methods
  #
  # cluster_cols       : exact cluster column names derived from clustering
  #                      parameters written by CLUSTER_AND_FINDMARKERS —
  #                      passed forward so ANNOTATE_CLUSTERS doesn't need to
  #                      re-derive them
  #
  # active_assay       : assay used for scoring — always SCT in this pipeline
  #                      but stored explicitly for audit trail

  integrated_seurat@misc$scoring_params <- list(
    active_assay       = active_assay,
    signatures_list    = signatures_list,
    seurat_module_cols = seurat_module_cols,
    ucell_module_cols  = ucell_module_cols,
    cluster_cols       = cluster_cols
  )

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 10: UMAP Plots
  # ═══════════════════════════════════════════════════════════════════════════
  # Two visualisations — one per scoring method:
  #   1. Seurat AddModuleScore — one panel per cell type, global color scale
  #   2. UCell scores — one panel per cell type, global color scale
  #
  # global_scale = TRUE so all cell type panels share the same color scale —
  # allows direct visual comparison of score intensity across cell types.
  # A T cell score of 0.5 maps to the same color as an Epithelial score of 0.5.
  #
  # Plotted before annotation runs — validates that marker genes produce
  # meaningful spatial patterns on the UMAP. If scores look wrong (e.g. T cell
  # markers scattered uniformly) → fix marker file and rerun before annotation.

  if (!is.null(output_dir)) {

    plot_seurat(integrated_seurat,
                features      = seurat_module_cols,
                feature_type  = "module_score",
                reduction     = "umap_harmony",
                filename      = "UMAP_Seurat_Module_Scores",
                output_dir    = output_dir)

    plot_seurat(integrated_seurat,
                features      = ucell_module_cols,
                feature_type  = "module_score",
                reduction     = "umap_harmony",
                filename      = "UMAP_UCell_Module_Scores",
                output_dir    = output_dir)
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 11: Save and Return
  # ═══════════════════════════════════════════════════════════════════════════
  # Saved as scored_seurat.rds — fixed filename, no parameter needed.
  #
  # WHY WE SAVE HERE AND NOT IN ANNOTATE_CLUSTERS:
  # UCell scoring on 100k cells takes 20-40 minutes on HPC.
  # Saving the scored object as a checkpoint means annotation parameters
  # (AVG_PROB_THRESHOLD, stability_threshold, delta threshold) can be
  # tuned by rerunning only ANNOTATE_CLUSTERS without repeating UCell.
  # This is the primary reason these two functions are kept separate.

  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    rds_out <- file.path(output_dir, "scored_seurat.rds")
    saveRDS(object = integrated_seurat, file = rds_out)
    log_info(sample = "", step = "calc_module_scores",
             msg = glue("Scored Seurat object saved to: '{rds_out}'"))
  } else {
    log_info(sample = "", step = "calc_module_scores",
             msg = "No output_dir provided — skipping save.")
  }

  log_info(sample = "", step = "calc_module_scores",
           msg = glue("Module scoring completed successfully. ",
                      "{length(seurat_module_cols)} Seurat + ",
                      "{length(ucell_module_cols)} UCell score columns added to metadata."))

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

  calc_module_scores(
    rds_file    = get_arg(1),
    marker_file = get_arg(2),
    output_dir  = get_arg(3, default = ".")
  )
}
