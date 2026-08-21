# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(Seurat)
  library(SeuratObject)
  library(DropletUtils)
  library(SingleCellExperiment)
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

load_cellranger_and_classify_droplets <- function(sample_id, sample_dir, gene_column = 2, output_dir = NULL) {

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Detect Input Type
  # ═══════════════════════════════════════════════════════════════════════════
  # Cell Ranger writes two matrix directories per sample:
  #   raw_feature_bc_matrix/      : all barcodes including empty droplets
  #   filtered_feature_bc_matrix/ : Cell Ranger's filtered cells only
  #
  # We always prefer raw — it preserves all barcodes for DropletUtils/emptyDrops
  # and lets the CellRanger classifier compare against it. If raw is absent
  # (user only has filtered output), fall back to filtered as the main object.
  # input_type drives which empty droplet classifiers run below.

  raw_dir  <- file.path(sample_dir, "raw_feature_bc_matrix")
  filt_dir <- file.path(sample_dir, "filtered_feature_bc_matrix")
  raw_h5   <- file.path(sample_dir, "raw_feature_bc_matrix.h5")
  filt_h5  <- file.path(sample_dir, "filtered_feature_bc_matrix.h5")

  has_raw_dir  <- dir.exists(raw_dir)
  has_filt_dir <- dir.exists(filt_dir)
  has_raw_h5   <- file.exists(raw_h5)
  has_filt_h5  <- file.exists(filt_h5)

  has_raw  <- has_raw_dir  || has_raw_h5
  has_filt <- has_filt_dir || has_filt_h5

  if (!has_raw && !has_filt) {
    log_error(sample = sample_id, step = "load_cellranger_and_classify_droplets",
              msg = glue::glue("Neither raw nor filtered data found under: '{sample_dir}'"))
  }

  input_type <- dplyr::case_when(
    has_raw  & has_filt  ~ "raw_and_filtered",
    has_raw  & !has_filt ~ "raw_only",
    !has_raw & has_filt  ~ "filtered_only"
  )

  log_info(sample = sample_id, step = "load_cellranger_and_classify_droplets",
           msg = glue::glue("Input type detected: '{input_type}'"))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Helper — Read a 10x Matrix (h5 preferred, MEX fallback)
  # ═══════════════════════════════════════════════════════════════════════════
  # h5 is a single file, faster to read, and less prone to partial write
  # corruption than the MEX (matrix.mtx.gz) directory format. Both formats
  # contain identical data. Shared here since both raw and filtered reads
  # (Sections 3 and 5) need identical logic.

  read_10x_matrix <- function(h5_path, dir_path, label) {

    if (!is.null(h5_path) && file.exists(h5_path)) {

      log_info(sample = sample_id, step = "load_cellranger_and_classify_droplets",
               msg = glue::glue("Reading HDF5 ({label}): '{basename(h5_path)}'"))

      tryCatch({
        result <- Seurat::Read10X_h5(filename        = h5_path,
                                     use.names       = TRUE,
                                     unique.features = TRUE)
        if (is.list(result)) {
          log_error(sample = sample_id, step = "load_cellranger_and_classify_droplets",
                    msg = "HDF5 file returned a list — multimodal data detected. Check file content.")
        }
        result
      }, error = function(e) {
        log_error(sample = sample_id, step = "load_cellranger_and_classify_droplets",
                  msg = glue::glue("Failed to read HDF5 ({label}): {e$message}"))
      })

    } else if (!is.null(dir_path) && dir.exists(dir_path)) {

      log_info(sample = sample_id, step = "load_cellranger_and_classify_droplets",
               msg = glue::glue("Reading MEX ({label}): '{basename(dir_path)}'"))

      tryCatch({
        Seurat::Read10X(data.dir        = dir_path,
                        gene.column     = gene_column,
                        cell.column     = 1,
                        unique.features = TRUE,
                        strip.suffix    = FALSE)
      }, error = function(e) {
        log_error(sample = sample_id, step = "load_cellranger_and_classify_droplets",
                  msg = glue::glue("Failed to read MEX ({label}): {e$message}"))
      })

    } else {
      NULL
    }
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Load Main Seurat Object
  # ═══════════════════════════════════════════════════════════════════════════
  # min.cells = 0, min.features = 0 — no filtering at this stage. All barcodes
  # must be retained for emptyDrops to model the ambient RNA profile correctly.

  if (has_raw) {
    counts <- read_10x_matrix(h5_path  = if (has_raw_h5) raw_h5 else NULL,
                              dir_path = if (has_raw_dir) raw_dir else NULL,
                              label    = "raw")
  } else {
    log_warn(sample = sample_id, step = "load_cellranger_and_classify_droplets",
             msg = "Raw matrix absent — loading filtered matrix as the main object.")
    counts <- read_10x_matrix(h5_path  = if (has_filt_h5) filt_h5 else NULL,
                              dir_path = if (has_filt_dir) filt_dir else NULL,
                              label    = "filtered")
  }

  sample_seurat <- SeuratObject::CreateSeuratObject(
    counts       = counts,
    project      = sample_id,
    assay        = "RNA",
    names.field  = 1,
    names.delim  = "_",
    min.cells    = 0,
    min.features = 0
  )

  sample_seurat@misc$input_type <- input_type
  all_barcodes                  <- colnames(sample_seurat)

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: DropletUtils Classification (emptyDrops)
  # ═══════════════════════════════════════════════════════════════════════════
  # emptyDrops tests each barcode against an ambient RNA profile estimated from
  # barcodes with very low total counts. Requires the full unfiltered matrix —
  # skipped when input is filtered_only, since the ambient pool is absent.

  set.seed(100)

  if (input_type %in% c("raw_and_filtered", "raw_only")) {

    sce    <- Seurat::as.SingleCellExperiment(x = sample_seurat)
    niters <- 10000

    df <- tryCatch({

      e.out <- DropletUtils::emptyDrops(m      = SingleCellExperiment::counts(sce),
                                        niters = niters)
      df    <- as.data.frame(e.out)

      # If Limited == TRUE and FDR > 0.05, those droplets need more iterations
      # to achieve a reliable FDR estimate. Iterate until none remain, capped
      # at MAX_NITERS to avoid looping forever on droplets that never resolve
      # below the FDR threshold (e.g. genuinely ambiguous low-count barcodes).
      MAX_NITERS <- 100000
      n_improve  <- nrow(dplyr::filter(df, Limited == TRUE, FDR > 0.05))

      while (n_improve > 0 && niters < MAX_NITERS) {
        niters    <- niters + 10000
        e.out     <- DropletUtils::emptyDrops(m      = SingleCellExperiment::counts(sce),
                                              niters = niters)
        df        <- as.data.frame(e.out)
        n_improve <- nrow(dplyr::filter(df, Limited == TRUE, FDR > 0.05))

        log_info(sample = sample_id, step = "load_cellranger_and_classify_droplets",
                 msg = glue::glue("{n_improve} droplets need more iterations. niters = {niters}"))
      }

      if (n_improve > 0) {
        log_warn(sample = sample_id, step = "load_cellranger_and_classify_droplets",
                 msg = glue::glue("Reached MAX_NITERS ({MAX_NITERS}) with {n_improve} droplet(s) ",
                                  "still Limited & FDR > 0.05. Proceeding with the last computed ",
                                  "result rather than looping indefinitely."))
      }

      df

    }, error = function(e) {

      # emptyDrops fails when there are no low-count barcodes to estimate the
      # ambient profile. Mark all barcodes as non-empty rather than failing.
      if (grepl("no counts available to estimate the ambient profile", e$message)) {
        log_warn(sample = sample_id, step = "load_cellranger_and_classify_droplets",
                 msg = "emptyDrops failed: no ambient profile. Marking all barcodes as non-empty.")

        data.frame(Total   = colSums(SingleCellExperiment::counts(sce)),
                   LogProb = NA,
                   PValue  = NA,
                   FDR     = 0,
                   Limited = FALSE,
                   row.names = colnames(sce))
      } else {
        log_error(sample = sample_id, step = "load_cellranger_and_classify_droplets",
                  msg = glue::glue("emptyDrops error: {e$message}"))
      }
    })

    # FDR <= 0.05 — standard significance threshold for non-empty classification.
    # NA barcodes were not tested (below emptyDrops minimum count threshold).
    true_barcodes <- rownames(dplyr::filter(df, FDR <= 0.05, !is.na(FDR)))

    dropletutils_calls <- dplyr::case_when(
      all_barcodes %in% true_barcodes ~ "Non-Empty Droplet",
      TRUE                            ~ "Empty Droplet"
    )
    names(dropletutils_calls) <- all_barcodes

    sample_seurat <- SeuratObject::AddMetaData(object   = sample_seurat,
                                               metadata = dropletutils_calls,
                                               col.name = "DropletUtils")

    log_info(sample = sample_id, step = "load_cellranger_and_classify_droplets",
             msg = glue::glue("DropletUtils: {sum(dropletutils_calls == 'Non-Empty Droplet')} ",
                              "non-empty droplets identified."))

  } else {

    log_info(sample = sample_id, step = "load_cellranger_and_classify_droplets",
             msg = "Skipping DropletUtils — filtered_only input has no ambient pool.")
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5: CellRanger Classification
  # ═══════════════════════════════════════════════════════════════════════════
  # CellRanger's filtered_feature_bc_matrix contains only barcodes it called as
  # cells. We check which barcodes in the full matrix appear in the filtered
  # set. Only meaningful when BOTH raw and filtered are present — reads just
  # the filtered matrix here (not a full Seurat object) to avoid the overhead
  # of building and discarding a second Seurat object for a barcode lookup.

  if (input_type == "raw_and_filtered") {

    filt_counts       <- read_10x_matrix(h5_path  = if (has_filt_h5) filt_h5 else NULL,
                                         dir_path = if (has_filt_dir) filt_dir else NULL,
                                         label    = "filtered")
    filtered_barcodes <- colnames(filt_counts)

    cellranger_calls <- dplyr::case_when(
      all_barcodes %in% filtered_barcodes ~ "Non-Empty Droplet",
      TRUE                                ~ "Empty Droplet"
    )
    names(cellranger_calls) <- all_barcodes

    sample_seurat <- SeuratObject::AddMetaData(object   = sample_seurat,
                                               metadata = cellranger_calls,
                                               col.name = "CellRanger")

    log_info(sample = sample_id, step = "load_cellranger_and_classify_droplets",
             msg = glue::glue("CellRanger: {sum(cellranger_calls == 'Non-Empty Droplet')} ",
                              "non-empty droplets identified."))

  } else {

    log_info(sample = sample_id, step = "load_cellranger_and_classify_droplets",
             msg = glue::glue("Skipping CellRanger classifier — ",
                              "filtered_feature_bc_matrix absent for input_type: '{input_type}'"))
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 6: Save and Return
  # ═══════════════════════════════════════════════════════════════════════════

  if (!is.null(output_dir)) {

    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    output_file <- file.path(output_dir, paste0(sample_id, ".LOAD_CLASSIFY.rds"))
    saveRDS(object = sample_seurat, file = output_file)

    log_info(sample = sample_id, step = "load_cellranger_and_classify_droplets",
             msg = glue::glue("Seurat object saved to: '{output_file}'"))

  } else {

    log_info(sample = sample_id, step = "load_cellranger_and_classify_droplets",
             msg = "No output_dir provided — skipping save.")
  }

  return(invisible(sample_seurat))
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

  load_cellranger_and_classify_droplets(
    sample_id   = get_arg(1),
    sample_dir  = get_arg(2),
    gene_column = as.integer(get_arg(3, default = 2)),
    output_dir  = get_arg(4, default = ".")
  )
}
