# =============================================================================
# 01_classification_association_survival.R
#
# Purpose: Using BostonGeneMFP TME classifications and ssGSEA scores,
#          test association between NPEPPS/ERAP1/ERAP2 expression and
#          TME class (Infiltrated/Desert/Fibrotic), run pathway correlations,
#          and perform survival analysis (KM + Cox).
#
# Scope of analysis:
#   - Wilcoxon (gene vs TME class)      : per-dataset only (40 cohorts)
#   - Gene vs pathway correlation       : per-dataset only (40 cohorts)
#   - Survival (KM + Cox)               : all 40 cohorts, single run,
#                                          faceted + per-facet cutoff by
#                                          Project_ID (functionally per-dataset)
#   - Heatmap                           : per-dataset (40) + one PanCancer
#                                          (33 TCGA pooled) heatmap
#
#   Gene expression (NPEPPS/ERAP1/ERAP2) is NOT pooled across cohorts for
#   Wilcoxon/correlation testing. DESeq2 log-norm values are only calibrated
#   within their own cohort's size factors — pooling raw values across 33
#   different cancer types confounds true TME association with cancer-type-
#   driven baseline expression differences. Pathway scores (ssGSEA, rescaled
#   0-1 per plot) don't have this issue, so PanCancer pooling is retained
#   only for the heatmap.
#
# Inputs:
#   combined_classification.tsv   — sample_id, TME (4-class), cohort
#   combined_ssgsea_scores.tsv    — sample_id, cohort, 29 ssGSEA scores
#   {Project_ID}_metadata.xlsx    — per-cohort clinical/survival data
#   {Project_ID}_Log_Norm_Counts.tsv — per-cohort log-norm expression
#
# Outputs (all under RESULTS_DIR):
#   master_sample_info.xlsx         — merged sample info + TME + scores
#   Association/
#     Gene_Infiltration_Association.xlsx
#       Wilcoxon_TME              — pairwise Wilcoxon per gene per dataset
#       Gene_Pathway_Correlation  — gene vs 22 pathway scores per dataset
#   Survival_Curves/
#     All_Dataset_{endpoint}_{variable}.pdf/.xlsx   (faceted by cohort)
#     Survival_Summary.xlsx
#   Heatmaps/
#     TME_Heatmap_PanCancer.pdf     — 33 TCGA pooled
#     TME_Heatmap_{cohort}.pdf      — one per dataset (40 total)
#
# Workflow:
#   S1 : Libraries, paths, parameters
#   S2 : Read classifications, map 4-class -> 3-class
#   S3 : Read ssGSEA scores
#   S4 : Read metadata (per-cohort, stack)
#   S5 : Extract NPEPPS/ERAP1/ERAP2 expression (per-cohort, stack)
#   S6 : Merge into master dataframe
#   S7 : Pairwise Wilcoxon (gene vs TME class) — per-dataset only
#   S8 : Gene vs 22 pathway correlations — per-dataset only
#   S9 : Survival analysis (KM + Cox) — all datasets, faceted
#   S10: Heatmap figures — per-dataset + PanCancer
# =============================================================================

# =============================================================================
# S1: LIBRARIES, PATHS, PARAMETERS
# =============================================================================

options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(openxlsx)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
})

source("/home/kailasamms/scripts/nextflow/modules/UTILS.R")

# --- Key parameters — modify here only ---
GENES_OF_INTEREST <- c("NPEPPS", "ERAP1", "ERAP2")

# TME groups — remove "Fibrotic" from this vector to drop it from all analyses
TME_GROUPS <- c("Infiltrated", "Desert", "Fibrotic")

TME_COLORS <- c(
  "Infiltrated" = "#C10020",
  "Desert"      = "#00538A",
  "Fibrotic"    = "#8B7355"
)

# 22 BostonGene signatures used for heatmap and pathway correlations
HEATMAP_22 <- c(
  "MHCI", "MHCII", "Coactivation_molecules", "Effector_cells",
  "T_cell_traffic", "NK_cells", "T_cells", "B_cells",
  "M1_signatures", "Th1_signature", "Antitumor_cytokines",
  "Checkpoint_inhibition", "Treg", "T_reg_traffic",
  "Neutrophil_signature", "Granulocyte_traffic", "MDSC",
  "MDSC_traffic", "Macrophages", "Macrophage_DC_traffic",
  "Th2_signature", "Protumor_cytokines"
)

# Survival endpoints
endpoints <- list(
  OS  = list(time_col = "OS.months",  status_col = "OS")
  #DSS = list(time_col = "DSS.months", status_col = "DSS"),
  #PFI = list(time_col = "PFI.months", status_col = "PFI")
)

MIN_N_SURVIVAL     <- 10
GENE_CUTOFF_METHOD <- "optimal"
MV_COVARIATES      <- c("Age", "Stage")

# --- Paths ---
BASE_DIR        <- "/home/kailasamms/scripts/nextflow/projects/Custom/Saswat_NPEPPS"
CLASS_DIR       <- file.path(BASE_DIR, "01.Classification")
RESULTS_DIR     <- file.path(BASE_DIR, "02.Results")
METADATA_DIR    <- "/home/kailasamms/scripts/nextflow/resources/Datasets/01.Metadata"
COUNTS_DIR      <- "/home/kailasamms/scripts/nextflow/resources/Datasets/03.Normalized"

ASSOC_DIR       <- file.path(RESULTS_DIR, "Association")
SURV_DIR        <- file.path(RESULTS_DIR, "Survival_Curves")
FIG_DIR         <- file.path(RESULTS_DIR, "Heatmaps")

for (p in c(RESULTS_DIR, ASSOC_DIR, SURV_DIR, FIG_DIR)) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

cat("S1: Parameters and paths set\n")
cat("  TME groups :", paste(TME_GROUPS, collapse = ", "), "\n")
cat("  Genes      :", paste(GENES_OF_INTEREST, collapse = ", "), "\n")

# =============================================================================
# S2: READ CLASSIFICATIONS, MAP 4-CLASS -> 3-CLASS
# =============================================================================

cat("\nS2: Reading TME classifications...\n")

# Standardize to Sample_ID throughout
class_df <- read.table(
  file.path(CLASS_DIR, "combined_classification.tsv"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
) %>%
  dplyr::rename(Sample_ID = sample_id)

cat("  Raw classifications:", nrow(class_df), "samples\n")
cat("  4-class distribution:\n")
print(table(class_df$TME))

# Map 4-class -> 3-class
class_df <- class_df %>%
  dplyr::mutate(
    TME_3class = dplyr::case_when(
      TME %in% c("IE", "IE/F") ~ "Infiltrated",
      TME == "D"               ~ "Desert",
      TME == "F"               ~ "Fibrotic",
      TRUE                     ~ NA_character_
    )
  ) %>%
  dplyr::filter(TME_3class %in% TME_GROUPS)

cat("  3-class distribution:\n")
print(table(class_df$TME_3class))
cat("  Samples after TME_GROUPS filter:", nrow(class_df), "\n")

# =============================================================================
# S3: READ SSGSEA SCORES
# =============================================================================

cat("\nS3: Reading ssGSEA scores...\n")

ssgsea_df <- read.table(
  file.path(CLASS_DIR, "combined_ssgsea_scores.tsv"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
) %>%
  dplyr::rename(Sample_ID = sample_id)

cat("  ssGSEA scores:", nrow(ssgsea_df), "samples x",
    ncol(ssgsea_df) - 2, "signatures\n")

# Keep only 22 heatmap signatures + Sample_ID + cohort
ssgsea_22 <- ssgsea_df %>%
  dplyr::select(Sample_ID, cohort, dplyr::all_of(HEATMAP_22))

cat("  Subset to 22 signatures for analysis\n")

# =============================================================================
# S4: READ METADATA (per-cohort, stack)
# =============================================================================

cat("\nS4: Reading metadata...\n")

meta_cols_keep <- c(
  "Sample_ID", "Project_ID", "Age", "Sex", "Stage",
  "OS", "OS.months"
)

meta_files <- list.files(METADATA_DIR, pattern = "_metadata\\.xlsx$",
                         full.names = TRUE)

# Exclude a pancancer file by name, if one ever exists — PanCancer here is
# always derived by filtering already-loaded per-cohort data (see S10), never
# read from a separate combined file.
meta_files <- meta_files[!grepl("pancancer", meta_files, ignore.case = TRUE)]

cat("  Found", length(meta_files), "metadata files\n")

meta_list <- lapply(meta_files, function(f) {
  df <- openxlsx::read.xlsx(f, sheet = 1, colNames = TRUE)

  if ("Age" %in% colnames(df)) {
    df$Age <- as.numeric(df$Age)
  }

  if ("pT" %in% colnames(df)) {
    df$pT <- as.character(df$pT)
  }

  # Ensure Stage exists and is character
  if ("Stage" %in% colnames(df)) {
    df$Stage <- as.character(df$Stage)

    # Apply standard cleaning inside the loop
    df$Stage <- stringr::str_to_upper(stringr::str_trim(df$Stage))

    df$Stage <- dplyr::case_when(
      df$Stage %in% c("STAGE I", "STAGE IA", "STAGE IB") ~ "I",
      df$Stage %in% c("STAGE II", "STAGE IIA", "STAGE IIB", "STAGE IIC") ~ "II",
      df$Stage %in% c("STAGE III", "STAGE IIIA", "STAGE IIIB", "STAGE IIIC") ~ "III",
      df$Stage %in% c("STAGE IV", "STAGE IVA", "STAGE IVB", "STAGE IVC") ~ "IV",
      df$Stage %in% c("STAGE 0", "IS") ~ "0",
      df$Stage %in% c("[DISCREPANCY]", "STAGE X", "I/II NOS") ~ NA_character_,
      TRUE ~ df$Stage
    )
  }

  # IMvigor210 specific subsetting
  if (grepl(pattern = "IMvigor210", x = f)){
    df <- df %>%
      dplyr::filter(Tissue == "bladder", `Received.platinum` == "Y",
                    `Sample.collected.pre-platinum` == "N")
  }
  # IMvigor010 specific subsetting
  if (grepl(pattern = "IMvigor010", x = f)){
    df <- df %>%
      dplyr::filter(ARM == "Atezolizumab (MPDL3280A) 1200 mg",
                    prior_neoadjuvant_chemotherapy == "YES")
  }

  # keep only columns that exist in this file
  keep <- intersect(meta_cols_keep, colnames(df))
  df <- df[, keep, drop = FALSE]

  cat(f, ":", nrow(df), "samples\n")

  return(df)
})

metadata <- dplyr::bind_rows(meta_list)

cat("  Stacked metadata:", nrow(metadata), "samples\n")

# =============================================================================
# S5: EXTRACT GENE EXPRESSION (per-cohort, stack)
# =============================================================================

cat("\nS5: Extracting", paste(GENES_OF_INTEREST, collapse = "/"),
    "expression...\n")

counts_files <- list.files(COUNTS_DIR,
                           pattern = "_Log_Norm_Counts\\.tsv$",
                           full.names = TRUE)
# Exclude a pancancer file by name, if one ever exists — same rationale as S4
counts_files <- counts_files[!grepl("pancancer", counts_files, ignore.case = TRUE)]

cat("  Found", length(counts_files), "log-norm count files\n")

expr_list <- lapply(counts_files, function(f) {

  cohort <- sub("_Log_Norm_Counts\\.tsv$", "", basename(f))

  # Read only the rows for genes of interest — avoids loading full matrix
  full <- data.table::fread(f, sep = "\t", header = TRUE, data.table = FALSE)

  df <- full %>%
    dplyr::filter(SYMBOL %in% GENES_OF_INTEREST) %>%
    dplyr::select(-dplyr::any_of(c("ID", "gene_type"))) %>%
    tibble::column_to_rownames("SYMBOL")

  if (nrow(df) == 0) return(NULL)

  t(df) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Sample_ID") %>%
    dplyr::mutate(cohort = cohort)
})

expr_df <- dplyr::bind_rows(expr_list) %>%
  dplyr::rename_with(
    .fn   = ~ paste0(.x, "_LogNorm"),
    .cols = dplyr::any_of(GENES_OF_INTEREST)
  )

gene_cols <- paste0(GENES_OF_INTEREST, "_LogNorm")

cat("  Expression matrix:", nrow(expr_df), "samples x",
    length(gene_cols), "genes\n")

# =============================================================================
# S6: MERGE INTO MASTER DATAFRAME
# =============================================================================

cat("\nS6: Merging into master dataframe...\n")

master_df <- class_df %>%
  dplyr::inner_join(metadata,  by = "Sample_ID") %>%
  dplyr::left_join(expr_df %>% dplyr::select(Sample_ID, dplyr::all_of(gene_cols)),
                   by = "Sample_ID") %>%
  dplyr::left_join(ssgsea_22 %>% dplyr::select(-cohort),
                   by = "Sample_ID")

cat("  Master dataframe:", nrow(master_df), "samples x",
    ncol(master_df), "columns\n")
cat("  TME_3class distribution:\n")
print(table(master_df$TME_3class, useNA = "ifany"))

# Save master sample info
master_path <- file.path(RESULTS_DIR, "master_sample_info.xlsx")
openxlsx::write.xlsx(master_df, master_path, overwrite = TRUE)
cat("  Saved:", master_path, "\n")

dataset_ids <- sort(unique(master_df$cohort))
cat("  Datasets:", length(dataset_ids), "\n")

# PanCancer = filtered view of the 33 TCGA rows already in master_df.
# Used ONLY for the heatmap (S10) — NOT for Wilcoxon/correlation, since
# pooling raw log-norm gene expression across 33 cancer types confounds
# TME association with cancer-type-driven baseline expression (see header).
tcga_df <- master_df %>% dplyr::filter(grepl("^TCGA-", cohort))
cat("  PanCancer (TCGA) samples for heatmap use:", nrow(tcga_df), "\n")

# =============================================================================
# S7: PAIRWISE WILCOXON (gene vs TME class) — per-dataset only
# =============================================================================

cat("\nS7: Pairwise Wilcoxon association tests (per-dataset)...\n")

# All pairs from TME_GROUPS
tme_pairs <- combn(TME_GROUPS, 2, simplify = FALSE)
cat("  TME pairs:", length(tme_pairs), "\n")

wilcox_results <- list()

run_wilcoxon <- function(df, dataset_label) {
  rows <- list()
  for (gene in GENES_OF_INTEREST) {
    gcol <- paste0(gene, "_LogNorm")
    if (!gcol %in% colnames(df)) next
    for (pair in tme_pairs) {
      g1 <- pair[1]; g2 <- pair[2]
      v1 <- df[[gcol]][df$TME_3class == g1 & !is.na(df[[gcol]])]
      v2 <- df[[gcol]][df$TME_3class == g2 & !is.na(df[[gcol]])]
      if (length(v1) < 3 || length(v2) < 3) next
      wt <- tryCatch(wilcox.test(v1, v2), error = function(e) NULL)
      if (is.null(wt)) next
      rows[[length(rows) + 1]] <- data.frame(
        Dataset   = dataset_label,
        Gene      = gene,
        Group1    = g1,
        Group2    = g2,
        n1        = length(v1),
        n2        = length(v2),
        median1   = median(v1),
        median2   = median(v2),
        W         = as.numeric(wt$statistic),
        pval      = wt$p.value,
        stringsAsFactors = FALSE
      )
    }
  }
  dplyr::bind_rows(rows)
}

# Per-dataset only (40 cohorts)
for (dataset in dataset_ids) {
  dataset_df <- master_df %>% dplyr::filter(cohort == dataset)
  res   <- run_wilcoxon(dataset_df, dataset)
  if (!is.null(res) && nrow(res) > 0) wilcox_results[[dataset]] <- res
}

# Combine and FDR correct
wilcox_df <- dplyr::bind_rows(wilcox_results) %>%
  dplyr::group_by(Gene, Group1, Group2) %>%
  dplyr::mutate(FDR = p.adjust(pval, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Gene, Group1, Group2, FDR)

cat("  Wilcoxon results:", nrow(wilcox_df), "tests\n")

# Save
wb_assoc <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb_assoc, "Wilcoxon_TME")
openxlsx::writeData(wb_assoc, "Wilcoxon_TME", wilcox_df)

# =============================================================================
# S8: GENE VS 22 PATHWAY CORRELATIONS — per-dataset only
# =============================================================================

cat("\nS8: Gene vs 22 pathway correlations (per-dataset)...\n")

corr_results <- list()

run_correlations <- function(df, dataset_label) {
  rows <- list()
  for (gene in GENES_OF_INTEREST) {
    gcol <- paste0(gene, "_LogNorm")
    if (!gcol %in% colnames(df)) next
    for (pw in HEATMAP_22) {
      if (!pw %in% colnames(df)) next
      vals <- df %>%
        dplyr::filter(!is.na(.data[[gcol]]), !is.na(.data[[pw]])) %>%
        dplyr::select(dplyr::all_of(c(gcol, pw)))
      if (nrow(vals) < 10) next
      ct <- tryCatch(
        cor.test(vals[[gcol]], vals[[pw]], method = "spearman"),
        error = function(e) NULL
      )
      if (is.null(ct)) next
      rows[[length(rows) + 1]] <- data.frame(
        Dataset   = dataset_label,
        Gene      = gene,
        Pathway   = pw,
        n         = nrow(vals),
        rho       = as.numeric(ct$estimate),
        pval      = ct$p.value,
        stringsAsFactors = FALSE
      )
    }
  }
  dplyr::bind_rows(rows)
}

# Per-dataset only (40 cohorts)
for (dataset in dataset_ids) {
  dataset_df <- master_df %>% dplyr::filter(cohort == dataset)
  res   <- run_correlations(dataset_df, dataset)
  if (!is.null(res) && nrow(res) > 0) corr_results[[dataset]] <- res
}

corr_df <- dplyr::bind_rows(corr_results) %>%
  dplyr::group_by(Gene, Pathway) %>%
  dplyr::mutate(FDR = p.adjust(pval, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Gene, Pathway, FDR)

cat("  Correlation results:", nrow(corr_df), "tests\n")

openxlsx::addWorksheet(wb_assoc, "Gene_Pathway_Correlation")
openxlsx::writeData(wb_assoc, "Gene_Pathway_Correlation", corr_df)

# Save association workbook
assoc_path <- file.path(ASSOC_DIR, "Gene_Infiltration_Association.xlsx")
openxlsx::saveWorkbook(wb_assoc, assoc_path, overwrite = TRUE)
cat("  Saved:", assoc_path, "\n")

# =============================================================================
# S9: SURVIVAL ANALYSIS (KM + Cox) — all datasets, faceted
# =============================================================================
# One call per endpoint/variable, faceted by cohort with global_cutoff = FALSE
# — this fits a separate optimal cutoff WITHIN each cohort before pooling for
# display, so it is functionally per-dataset even though it's one call and
# one output file per variable (40 facets in a single plot/workbook).

cat("\nS9: Survival analysis (all datasets, faceted by cohort)...\n")

all_cox_stats    <- list()
all_mv_cox_stats <- list()

# Create a vector of all continuous variables you want to test
surv_features <- c(paste0(GENES_OF_INTEREST, "_LogNorm"), "Effector_cells")

for (endpoint in names(endpoints)) {

  ep      <- endpoints[[endpoint]]
  n_valid <- sum(!is.na(master_df[[ep$time_col]]) &
                   !is.na(master_df[[ep$status_col]]))
  if (n_valid < MIN_N_SURVIVAL) next

  # (a) TME_3class
  tme_df <- master_df %>%
    dplyr::filter(TME_3class %in% TME_GROUPS) %>%
    dplyr::mutate(TME_3class = factor(TME_3class, levels = TME_GROUPS))

  if (nrow(tme_df) >= MIN_N_SURVIVAL) {
    result <- tryCatch(
      plot_survival(
        metadata      = tme_df,
        stratify_var  = "TME_3class",
        time_col      = ep$time_col,
        status_col    = ep$status_col,
        facet_var     = "cohort",
        output_dir    = SURV_DIR,
        filename      = paste0("All_Dataset_", endpoint, "_TME")
      ),
      error = function(e) { cat("  WARNING:", endpoint, "TME failed —", conditionMessage(e), "\n"); NULL }
    )
    if (!is.null(result$cox_stats))
      all_cox_stats[[paste0("All_Dataset_", endpoint, "_TME")]] <- result$cox_stats %>%
        dplyr::mutate(Run = paste0("All_Dataset_", endpoint, "_TME"))
  }

  # (b) Features — univariate + multivariate
  for (feature in surv_features) {
    feature_data <- master_df %>% dplyr::filter(!is.na(.data[[feature]]))
    if (nrow(feature_data) < MIN_N_SURVIVAL) next

    expr_data <- feature_data %>%
      dplyr::select(Sample_ID, .data[[feature]]) %>%
      dplyr::distinct(Sample_ID, .keep_all = TRUE) %>%
      tidyr::pivot_wider(names_from = Sample_ID, values_from = .data[[feature]]) %>%
      dplyr::mutate(SYMBOL = feature) %>%
      dplyr::select(SYMBOL, everything())

    # Drop the feature column from the metadata argument itself — feature_data
    # (= master_df subset) already contains this column, which otherwise
    # collides with expr_data's SYMBOL match and trips plot_survival()'s
    # ambiguity check (Section 1h), silently aborting via tryCatch.
    feature_meta <- feature_data %>% dplyr::select(-dplyr::all_of(feature))

    available_covs <- intersect(MV_COVARIATES, colnames(feature_meta))

    available_covs <- available_covs[
      sapply(available_covs, function(cv) sum(!is.na(feature_meta[[cv]])) >= MIN_N_SURVIVAL)
    ]

    # Univariate
    result <- tryCatch(
      plot_survival(
        metadata      = feature_meta,
        stratify_var  = feature,
        time_col      = ep$time_col,
        status_col    = ep$status_col,
        expr_data     = expr_data,
        facet_var     = "cohort",
        global_cutoff = FALSE,
        cutoff_method = GENE_CUTOFF_METHOD,
        output_dir    = SURV_DIR,
        filename      = paste0("All_Dataset_", endpoint, "_", feature)
      ),
      error = function(e) { cat("  WARNING:", endpoint, feature, "failed —", conditionMessage(e), "\n"); NULL }
    )
    if (!is.null(result$cox_stats))
      all_cox_stats[[paste0("All_Dataset_", endpoint, "_", feature)]] <- result$cox_stats %>%
        dplyr::mutate(Run = paste0("All_Dataset_", endpoint, "_", feature))

    # Multivariate
    if (length(available_covs) > 0) {
      result <- tryCatch(
        plot_survival(
          metadata       = feature_meta,
          stratify_var   = feature,
          time_col       = ep$time_col,
          status_col     = ep$status_col,
          expr_data      = expr_data,
          facet_var      = "cohort",
          global_cutoff  = FALSE,
          cutoff_method  = GENE_CUTOFF_METHOD,
          covariate_cols = available_covs,
          output_dir     = SURV_DIR,
          filename       = paste0("All_Dataset_", endpoint, "_", feature, "_MV")
        ),
        error = function(e) { cat("  WARNING:", endpoint, feature, "MV failed —", conditionMessage(e), "\n"); NULL }
      )
      if (!is.null(result$mv_cox_stats))
        all_mv_cox_stats[[paste0("All_Dataset_", endpoint, "_", feature, "_MV")]] <- result$mv_cox_stats %>%
          dplyr::mutate(Run = paste0("All_Dataset_", endpoint, "_", feature, "_MV"))
    }
  }
}

# ── Combine cox stats ─────────────────────────────────────────────────────────
wb_surv <- openxlsx::createWorkbook()

if (length(all_cox_stats) > 0) {
  cox_summary <- dplyr::bind_rows(all_cox_stats) %>%
    dplyr::arrange(Run)
  openxlsx::addWorksheet(wb_surv, "Cox_Summary")
  openxlsx::writeData(wb_surv, "Cox_Summary", cox_summary)
  cat("  Cox_Summary:", nrow(cox_summary), "rows\n")
}

if (length(all_mv_cox_stats) > 0) {
  mv_regex  <- paste0("^(", paste(c(GENES_OF_INTEREST, "Effector_cells"), collapse = "|"), ")")
  mv_summary <- dplyr::bind_rows(all_mv_cox_stats) %>%
    dplyr::filter(grepl(mv_regex, Term)) %>%
    dplyr::arrange(Run)
  openxlsx::addWorksheet(wb_surv, "MV_Cox_Genes")
  openxlsx::writeData(wb_surv, "MV_Cox_Genes", mv_summary)
  cat("  MV_Cox_Genes:", nrow(mv_summary), "rows\n")
}

if (length(all_cox_stats) > 0 || length(all_mv_cox_stats) > 0) {
  openxlsx::saveWorkbook(wb_surv,
                         file.path(SURV_DIR, "Survival_Summary.xlsx"), overwrite = TRUE)
  cat("  Saved: Survival_Summary.xlsx\n")
}

# =============================================================================
# S10: HEATMAP FIGURES — per-dataset (40) + PanCancer (33 TCGA pooled)
# =============================================================================

cat("\nS10: Building TME heatmap figures...\n")

rescale_01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0.5, length(x)))
  (x - rng[1]) / diff(rng)
}

# One function, reused for every heatmap (per-dataset and PanCancer alike).
# Rescaling is computed WITHIN whatever sample set is passed in (per point 1
# from the discussion — each heatmap's own min/max, not a shared global scale).
make_tme_heatmap <- function(df, label, out_dir) {

  heatmap_samples <- df %>%
    dplyr::filter(!is.na(TME_3class), TME_3class %in% TME_GROUPS) %>%
    dplyr::select(Sample_ID, TME_3class, dplyr::all_of(HEATMAP_22)) %>%
    tidyr::drop_na()

  if (nrow(heatmap_samples) < MIN_N_SURVIVAL) {
    cat("  SKIPPING", label, "— only", nrow(heatmap_samples), "samples with complete data\n")
    return(invisible(NULL))
  }

  heatmap_mat <- t(as.matrix(heatmap_samples[, HEATMAP_22]))
  colnames(heatmap_mat) <- make.names(heatmap_samples$Sample_ID)

  heatmap_mat_scaled <- t(apply(heatmap_mat, 1, rescale_01))

  mean_score <- colMeans(heatmap_mat_scaled, na.rm = TRUE)
  tme_vec    <- heatmap_samples$TME_3class
  names(tme_vec) <- make.names(heatmap_samples$Sample_ID)

  col_order <- unlist(lapply(TME_GROUPS, function(grp) {
    samps <- names(tme_vec)[tme_vec == grp]
    samps[order(mean_score[samps])]
  }))

  heatmap_mat_scaled <- heatmap_mat_scaled[, col_order, drop = FALSE]
  tme_annotation_vec <- tme_vec[col_order]

  col_fun <- circlize::colorRamp2(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    colors = c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B")
  )

  top_ann <- ComplexHeatmap::HeatmapAnnotation(
    TME_Class = tme_annotation_vec,
    col       = list(TME_Class = TME_COLORS[names(TME_COLORS) %in% TME_GROUPS]),
    annotation_height    = unit(4, "mm"),
    annotation_name_side = "left",
    annotation_name_gp   = gpar(fontsize = 9),
    annotation_legend_param = list(
      TME_Class = list(
        title    = "TME Class",
        title_gp = gpar(fontsize = 9, fontface = "bold"),
        labels_gp = gpar(fontsize = 8)
      )
    )
  )

  ht <- ComplexHeatmap::Heatmap(
    matrix            = heatmap_mat_scaled,
    name              = "Rescaled\nPathway Score",
    col               = col_fun,
    row_names_side    = "left",
    row_names_gp      = gpar(fontsize = 8),
    cluster_rows      = FALSE,
    show_row_dend     = FALSE,
    cluster_columns   = FALSE,
    show_column_names = FALSE,
    show_column_dend  = FALSE,
    column_title      = paste0("BostonGene TME Classification — ", label),
    column_title_gp   = gpar(fontsize = 14, fontface = "bold"),
    top_annotation    = top_ann,
    border            = FALSE,
    rect_gp           = gpar(col = NA),
    heatmap_legend_param = list(
      title      = "Rescaled\nPathway Score",
      title_gp   = gpar(fontsize = 9, fontface = "bold"),
      labels_gp  = gpar(fontsize = 8),
      at         = c(0, 0.25, 0.5, 0.75, 1),
      labels     = c("0 (min)", "0.25", "0.5", "0.75", "1 (max)")
    )
  )

  heatmap_pdf <- file.path(out_dir, paste0("TME_Heatmap_", label, ".pdf"))
  plot_width  <- max(8, min(40, ncol(heatmap_mat_scaled) * 0.05))

  grDevices::cairo_pdf(filename = heatmap_pdf, width = plot_width, height = 8)
  ComplexHeatmap::draw(ht, merge_legend = TRUE,
                       heatmap_legend_side    = "right",
                       annotation_legend_side = "right")
  grDevices::dev.off()

  cat("  Saved:", basename(heatmap_pdf), "(", nrow(heatmap_samples), "samples )\n")
}

# PanCancer (33 TCGA pooled)
make_tme_heatmap(tcga_df, "PanCancer", FIG_DIR)

# One per individual dataset (40 total: TCGA + GEO + IMvigor)
for (dataset in dataset_ids) {
  dataset_df <- master_df %>% dplyr::filter(cohort == dataset)
  make_tme_heatmap(dataset_df, dataset, FIG_DIR)
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("COMPLETE\n")
cat(strrep("=", 60), "\n")
cat("Datasets          :", length(dataset_ids), "\n")
cat("Total samples     :", nrow(master_df), "\n")
print(table(master_df$TME_3class))
cat("Results directory :", RESULTS_DIR, "\n")
cat(strrep("=", 60), "\n")
