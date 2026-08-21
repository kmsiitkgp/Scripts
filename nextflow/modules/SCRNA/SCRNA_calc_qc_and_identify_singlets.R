# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(Seurat)
  library(SeuratObject)
  library(scDblFinder)
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

calc_qc_and_identify_singlets <- function(sample_id, sample_rds,
                                          gene_cutoff    = 250,
                                          umi_cutoff     = 500,
                                          mito_cutoff    = 0.2,
                                          novelty_cutoff = 0.8,
                                          ribo_cutoff    = 0.05,
                                          multiplet_rate = 0.008,
                                          output_dir     = NULL) {

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Load Seurat Object
  # ═══════════════════════════════════════════════════════════════════════════

  sample_seurat <- load_smart(input_path = sample_rds)
  all_barcodes  <- colnames(sample_seurat)

  # Check if the object contains spatial images instead of checking DefaultAssay
  is_spatial <- length(sample_seurat@images) > 0
  assay      <- if (is_spatial) "Spatial" else SeuratObject::DefaultAssay(sample_seurat)

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Populate Missing Optional Classifier Columns
  # ═══════════════════════════════════════════════════════════════════════════
  # DropletUtils, CellRanger, HTO_Final may be absent depending on input_type
  # or experimental design. Fill with "ND" for safe downstream case_when logic.

  optional_cols <- c("DropletUtils", "CellRanger", "HTO_Final")
  for (col in optional_cols) {
    if (!col %in% colnames(sample_seurat@meta.data)) {
      sample_seurat <- SeuratObject::AddMetaData(object   = sample_seurat,
                                                 metadata = "ND",
                                                 col.name = col)
    }
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Compute Percentage Feature Sets
  # ═══════════════════════════════════════════════════════════════════════════
  # Patterns are case-insensitive to handle both human (MT-) and mouse (mt-)
  # gene naming conventions without separate logic per species.

  # Mitochondrial genes — high mito ratio indicates damaged/dying cells that
  # have lost cytoplasmic RNA but retained mitochondrial transcripts
  sample_seurat <- Seurat::PercentageFeatureSet(object   = sample_seurat,
                                                pattern  = "^[Mm][Tt]-",
                                                col.name = "MitoPercent",
                                                assay    = assay)

  # Ribosomal genes — very high ribo ratio can indicate low-complexity cells
  sample_seurat <- Seurat::PercentageFeatureSet(object   = sample_seurat,
                                                pattern  = "^[Rr][Pp][SsLl]",
                                                col.name = "RiboPercent",
                                                assay    = assay)

  # Hemoglobin genes — presence indicates red blood cell contamination.
  # No trailing hyphen: human symbols have none (HBA1, HBA2, HBB, HBD), while
  # mouse symbols do (Hba-a1, Hba-a2, Hbb-bs, Hbb-bt). Requiring a hyphen here
  # would only ever match mouse, leaving HemePercent at 0 for human samples.
  sample_seurat <- Seurat::PercentageFeatureSet(object   = sample_seurat,
                                                pattern  = "^[Hh][Bb][AaBb]",
                                                col.name = "HemePercent",
                                                assay    = assay)

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Build QC Metadata Table
  # ═══════════════════════════════════════════════════════════════════════════

  sample_metadata <- sample_seurat@meta.data %>%
    tibble::rownames_to_column(var = "barcode") %>%
    dplyr::mutate(
      Cell      = paste0(orig.ident, "_", barcode),
      Sample    = orig.ident,
      nUMIs     = .data[[paste0("nCount_",   assay)]],
      nGenes    = .data[[paste0("nFeature_", assay)]],
      MitoRatio = MitoPercent / 100,
      RiboRatio = RiboPercent / 100,
      HemeRatio = HemePercent / 100,
      # Novelty: log10(genes) / log10(UMIs) — measures transcriptional complexity.
      # Low novelty (many UMIs but few genes) flags low-complexity cells like
      # red blood cells or droplets capturing a single dominant transcript.
      Novelty   = log10(nGenes) / log10(nUMIs)
    )

  # Rename HTO columns if present — only relevant for multiplexed experiments
  if ("nCount_HTO" %in% colnames(sample_metadata)) {
    sample_metadata <- sample_metadata %>%
      dplyr::rename(nHTO_UMIs = nCount_HTO,
                    nHTOs     = nFeature_HTO)
  }

  # ---- 📥 Add Spatial Coordinates (if available) ----

  if (length(sample_seurat@images) > 0){

    # Collect spatial coordinates into a single data frame
    df_coords <- data.frame()

    for (image_name in names(sample_seurat@images)){
      df <- data.frame(barcode = sample_seurat@images[[image_name]]@boundaries$centroids@cells,
                       X = sample_seurat@images[[image_name]]@boundaries$centroids@coords[,1],
                       Y = sample_seurat@images[[image_name]]@boundaries$centroids@coords[,2])
      df_coords <- dplyr::bind_rows(df_coords, df)
    }

    # Join coordinates based on barcode
    sample_metadata <- sample_metadata %>%
      dplyr::left_join(df_coords, by=c("barcode"="barcode"))
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5: Apply QC Cutoffs and Classify Quality
  # ═══════════════════════════════════════════════════════════════════════════
  # Quality classification hierarchy:
  #   Empty Droplet : BOTH DropletUtils AND CellRanger agree — high specificity
  #                   ensures true empty droplets are removed without losing
  #                   real low-UMI cells. "ND" columns (absent classifiers) are
  #                   never equal to "Empty Droplet" so they don't trigger this.
  #   Low Quality   : fails any QC metric cutoff
  #   High Quality  : passes all cutoffs and not flagged as empty

  log_info(sample = sample_id, step = "calc_qc_and_identify_singlets",
           msg = glue::glue("QC cutoffs — nGenes >= {gene_cutoff} | ",
                            "nUMIs >= {umi_cutoff} | ",
                            "MitoRatio <= {mito_cutoff} | ",
                            "Novelty >= {novelty_cutoff}"))

  sample_metadata <- sample_metadata %>%
    dplyr::mutate(
      # NOTE: ribo_cutoff is intentionally NOT included here. RiboRatio is
      # computed and reported (and ribo_cutoff is recorded in qc_params below)
      # for informational/downstream-review purposes only - high ribosomal
      # content alone is not treated as a hard QC failure by this pipeline.
      passes_qc = nGenes    >= gene_cutoff    &
        nUMIs     >= umi_cutoff     &
        MitoRatio <= mito_cutoff    &
        Novelty   >= novelty_cutoff,
      Quality = dplyr::case_when(
        # Both ran and agree — high specificity
        DropletUtils == "Empty Droplet" & CellRanger == "Empty Droplet" ~ "Empty Droplet",
        # Only DropletUtils ran (raw_only) — trust it alone
        DropletUtils == "Empty Droplet" & CellRanger == "ND"            ~ "Empty Droplet",
        # Only CellRanger ran (filtered_only) — trust it alone
        DropletUtils == "ND"            & CellRanger == "Empty Droplet" ~ "Empty Droplet",
        !passes_qc                                                      ~ "Low Quality",
        TRUE                                                            ~ "High Quality"
      )
    ) %>%
    dplyr::select(-passes_qc) %>%
    tibble::column_to_rownames(var = "barcode")

  # Store QC params — MERGE_AND_PLOT_QC reads cutoffs from here, not CLI args
  sample_seurat@misc$qc_params <- list(
    gene_cutoff    = gene_cutoff,
    umi_cutoff     = umi_cutoff,
    mito_cutoff    = mito_cutoff,
    novelty_cutoff = novelty_cutoff,
    ribo_cutoff    = ribo_cutoff
  )

  sample_seurat <- SeuratObject::AddMetaData(object   = sample_seurat,
                                             metadata = sample_metadata)

  n_high  <- sum(sample_metadata$Quality == "High Quality",  na.rm = TRUE)
  n_low   <- sum(sample_metadata$Quality == "Low Quality",   na.rm = TRUE)
  n_empty <- sum(sample_metadata$Quality == "Empty Droplet", na.rm = TRUE)

  log_info(sample = sample_id, step = "calc_qc_and_identify_singlets",
           msg = glue::glue("Quality summary — High: {n_high} | ",
                            "Low: {n_low} | ",
                            "Empty: {n_empty}"))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 6: Subset High Quality Cells for Doublet Detection
  # ═══════════════════════════════════════════════════════════════════════════
  # scDblFinder must run on high quality cells only — empty droplets and low
  # quality cells would distort the doublet score distribution and produce
  # unreliable classifications. Spatial data is skipped entirely — a spot
  # capturing multiple cells is expected by design, not a doublet artifact.

  if (is_spatial) {

    log_info(sample = sample_id, step = "calc_qc_and_identify_singlets",
             msg = "Spatial data detected — skipping doublet detection.")

    sample_seurat <- SeuratObject::AddMetaData(object   = sample_seurat,
                                               metadata = rep("ND", ncol(sample_seurat)),
                                               col.name = "scDblFinder")

  } else {

    subset_seurat         <- base::subset(x = sample_seurat, subset = Quality == "High Quality")
    high_quality_barcodes <- colnames(subset_seurat)

    # scDblFinder needs a minimum number of cells to build a reliable kNN graph.
    # Below 50 cells, skip doublet detection and treat all as singlets.
    if (ncol(subset_seurat) < 50) {

      log_warn(sample = sample_id, step = "calc_qc_and_identify_singlets",
               msg = glue::glue("Only {ncol(subset_seurat)} high quality cells — ",
                                "skipping doublet detection, treating all as Singlets."))

      scdblfinder_calls                        <- rep("ND", length(all_barcodes))
      names(scdblfinder_calls)                 <- all_barcodes
      scdblfinder_calls[high_quality_barcodes] <- "Singlet"

    } else {

      # ═══════════════════════════════════════════════════════════════════════
      # SECTION 7: scDblFinder Classification
      # ═══════════════════════════════════════════════════════════════════════
      # scDblFinder runs directly on raw counts — no PCA/normalization needed,
      # it handles preprocessing internally.

      set.seed(100)

      sce <- Seurat::as.SingleCellExperiment(x = subset_seurat)
      sce <- scDblFinder::scDblFinder(sce       = sce,
                                      clusters  = NULL,
                                      samples   = NULL,
                                      dbr.per1k = multiplet_rate)

      scdbl_calls            <- base::tolower(as.character(sce$scDblFinder.class))
      doublet_barcodes_scdbl <- colnames(sce)[scdbl_calls == "doublet"]

      scdblfinder_calls                         <- rep("ND", length(all_barcodes))
      names(scdblfinder_calls)                  <- all_barcodes
      scdblfinder_calls[high_quality_barcodes]  <- "Singlet"
      scdblfinder_calls[doublet_barcodes_scdbl] <- "Doublet"

      log_info(sample = sample_id, step = "calc_qc_and_identify_singlets",
               msg = glue::glue("scDblFinder: {sum(scdblfinder_calls == 'Doublet')} doublets identified."))
    }

    sample_seurat <- SeuratObject::AddMetaData(object   = sample_seurat,
                                               metadata = scdblfinder_calls,
                                               col.name = "scDblFinder")
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 8: Final Barcode-Level QC Classification
  # ═══════════════════════════════════════════════════════════════════════════
  # No cells are removed here — Barcode_QC is a classification column only.
  # Actual subsetting to singlets happens later, in MERGE_AND_PLOT_QC, so that
  # the full unfiltered per-sample object (and its QC plots) are never lost.
  #
  #   Empty Droplet : from Quality
  #   Low Quality   : from Quality
  #   Doublet       : scDblFinder calls doublet
  #   Singlet       : passes all above

  metadata <- sample_seurat@meta.data %>%
    dplyr::mutate(
      Barcode_QC = dplyr::case_when(
        Quality      == "Empty Droplet" ~ "Empty Droplet",
        Quality      == "Low Quality"   ~ "Low Quality",
        scDblFinder  == "Doublet"       ~ "Doublet",
        TRUE                            ~ "Singlet"
      )
    )

  # Retain only relevant columns — keep metadata lean for the merge step.
  # No barcodes are dropped by this select() — only unused columns.
  keep_cols <- c("Cell", "Sample", "Barcode_QC", "nUMIs", "nGenes",
                 "MitoRatio", "RiboRatio", "HemeRatio", "Novelty", "Quality",
                 "DropletUtils", "CellRanger", "scDblFinder",
                 "nHTO_UMIs", "nHTOs", "HTO_Final", "X", "Y")

  metadata <- metadata %>%
    dplyr::select(dplyr::any_of(keep_cols))

  sample_seurat@meta.data <- metadata

  n_singlets <- sum(metadata$Barcode_QC == "Singlet", na.rm = TRUE)
  n_doublets <- sum(metadata$Barcode_QC == "Doublet", na.rm = TRUE)
  n_low      <- sum(metadata$Barcode_QC == "Low Quality", na.rm = TRUE)
  n_empty    <- sum(metadata$Barcode_QC == "Empty Droplet", na.rm = TRUE)

  log_info(sample = sample_id, step = "calc_qc_and_identify_singlets",
           msg = glue::glue("Barcode_QC — Singlets: {n_singlets} | ",
                            "Doublets: {n_doublets} | ",
                            "Low Quality: {n_low} | ",
                            "Empty: {n_empty}"))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 9: Save and Return
  # ═══════════════════════════════════════════════════════════════════════════
  # Full, unfiltered object saved — no cells removed at this stage.

  if (!is.null(output_dir)) {

    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    output_file <- file.path(output_dir, paste0(sample_id, ".CALC_QC_IDENTIFY.rds"))
    saveRDS(object = sample_seurat, file = output_file)

    log_info(sample = sample_id, step = "calc_qc_and_identify_singlets",
             msg = glue::glue("Seurat object saved to: '{output_file}'"))

  } else {

    log_info(sample = sample_id, step = "calc_qc_and_identify_singlets",
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

  calc_qc_and_identify_singlets(
    sample_id      = get_arg(1),
    sample_rds     = get_arg(2),
    gene_cutoff    = as.numeric(get_arg(3, default = 250)),
    umi_cutoff     = as.numeric(get_arg(4, default = 500)),
    mito_cutoff    = as.numeric(get_arg(5, default = 0.2)),
    novelty_cutoff = as.numeric(get_arg(6, default = 0.8)),
    ribo_cutoff    = as.numeric(get_arg(7, default = 0.05)),
    multiplet_rate = as.numeric(get_arg(8, default = 0.008)),
    output_dir     = get_arg(9, default = ".")
  )
}
