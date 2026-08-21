#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════════════════════
# 📄 Script   : 08b_biomodal_dmr_analysis.R
# 📌 Purpose  : Full DMR analysis pipeline for Carotuximab cfDNA study.
#               Pre vs Post treatment comparison using 5hmC and 5mC marks.
# “We used CpG-bias–corrected, baseline-normalized gene-level epigenetic effect 
# sizes to perform pathway analysis, reducing gene-architecture–driven 
# enrichment artifacts.CpG correction was applied to prevent enrichment bias 
# driven by gene length and CpG density, which is a known confounder in 
# methylation-based GSEA"
# ===============================================================================
# 🎯 ANALYSIS DESIGN
# Selected configuration: paired_no_od_cov (from 08a_biomodal_dmr_screen.R)
#
# 1. DESIGN   [paired]   : Paired Pre vs Post per patient — maximises power
#                          by using each patient as their own control.
# 2. OD       [no_od]    : Overdispersion correction too conservative for this
#                          dataset — wiped out 5hmC signal entirely.
# 3. COVARIATES [cov]    : Retained to protect against batch effects.
# 4. COHORTS             : Short and Long only — Full cohort excluded.
#
# BIOLOGICAL FILTERS:
#   5hmC → Gene_body regions only
#   5mC  → Promoter regions only
#
# GSEA RANKING METRIC: CpG-corrected per-patient z-score
#   Rationale: Raw mod_lfc is biased toward large genes with more CpGs.
#   Z-scores normalise each gene relative to its own Pre baseline per patient.
#   CpG count bias is removed by regressing out log10(total_CpGs).
#   5mC z-score is negated for GSEA so positive rank = de-repression.
#
# ORA HIT SELECTION:
#   ORA1/2 : padj <= 0.05 AND n_CpGs >= 10 AND |mean_Z| >= 1
#   ORA3   : intersect of ORA1/ORA2 concordant genes (|mean_Z| >= 1)
#   ORA4/5 : |mean_Z| >= 1 AND percent_consistent >= 0.75
#   ORA6   : intersect of ORA4/ORA5 concordant genes
#
# ===============================================================================
# STRUCTURE
#
# ── BLOCK 1: DATA PREPARATION (always run) ────────────────────────────────────
#   (A) QC           — verify Biomodal output integrity
#   (B) Read DMR     — read TSVs, rename columns, biological filter, compute LFC
#   (C) Z-scores     — per-patient z-scores, CpG bias correction
#   (D) Master DMR   — join z-scores + patient consistency to DMR data
#   (E) Save outputs — Excel (4 files × 4 sheets) + key objects summary
#
# ── BLOCK 2: ANALYSIS MODULES (independent, use Block 1 objects) ──────────────
#   (F) DMR funnel   — print filter counts from master
#   (G) Venn         — 6 Venns (2 marks × 3 filter types)
#   (H) GSEA         — 2 cohorts × 2 marks = 4 GSEA runs
#   (I) ORA          — 2 cohorts × 6 ORA types = 12 ORA runs
#
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(openxlsx)
  library(matrixStats)
  library(VennDiagram)
  library(ggplot2)
})

# ── Paths ──────────────────────────────────────────────────────────────────────

HOME_DIR     <- "/home/kailasamms"
BASE_DIR     <- file.path(HOME_DIR, "analysis/Biomodal/Kevin_Carotuximab")
MODALITY_DIR <- file.path(BASE_DIR, "03.modality")
DMR_DIR      <- file.path(BASE_DIR, "04.dmr")
OUTPUT_DIR   <- file.path(BASE_DIR, "05.analysis")
GMT_DIR      <- "/home/kailasamms/resources/signatures/Human"
PROJ_DIR     <- file.path(HOME_DIR, "scripts/nextflow/projects/Biomodal/Kevin_Carotuximab")

source(file.path(HOME_DIR, "scripts/nextflow/modules/UTILS.R"))

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Modality output file paths ─────────────────────────────────────────────────

COUNT_TOTAL_C <- file.path(MODALITY_DIR, "get_count",         "count_total_c.tsv")
FRAC_HMC      <- file.path(MODALITY_DIR, "get_regional_frac", "frac_hmc.tsv")
FRAC_MC       <- file.path(MODALITY_DIR, "get_regional_frac", "frac_mc.tsv")
SUM_HMC       <- file.path(MODALITY_DIR, "get_sum",           "sum_hmc.tsv")
SUM_MC        <- file.path(MODALITY_DIR, "get_sum",           "sum_mc.tsv")
SUM_TOTAL_C   <- file.path(MODALITY_DIR, "get_sum",           "sum_total_c.tsv")

# ── DMR file paths (short and long only) ──────────────────────────────────────

dmr_subdirs <- list(
  short_paired = file.path(DMR_DIR, "dmr_short_paired_no_od_cov"),
  long_paired  = file.path(DMR_DIR, "dmr_long_paired_no_od_cov")
)

COHORTS <- c("short_paired", "long_paired")
MARKS   <- c("hmc", "mc")

# ── Metadata paths (one per cohort) ───────────────────────────────────────────

meta_paths <- list(
  short_paired = file.path(PROJ_DIR, "Kevin_Carotuximab_metadata_short_paired.csv"),
  long_paired  = file.path(PROJ_DIR, "Kevin_Carotuximab_metadata_long_paired.csv")
)

# ── Analysis constants ─────────────────────────────────────────────────────────

META_COLS          <- c("Chromosome", "Start", "End", "Annotation", "Name")
PADJ_CUTOFF        <- 0.05
N_CPGS_CUTOFF      <- 10
LFC_CUTOFF         <- 0.58
ZSCORE_CUTOFF      <- 1.0
CONSISTENCY_CUTOFF <- 0.75    # 75% of patients must agree on direction

# ═══════════════════════════════════════════════════════════════════════════════
# ██████████████████████████████████████████████████████████████████████████████
# BLOCK 1 — DATA PREPARATION
# ██████████████████████████████████████████████████████████████████████████████
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n")
cat("██████████████████████████████████████████████████████\n")
cat(" BLOCK 1 — DATA PREPARATION\n")
cat("██████████████████████████████████████████████████████\n\n")

# ═══════════════════════════════════════════════════════════════════════════════
# (A) QC — Cross-check Biomodal output for internal consistency
# ═══════════════════════════════════════════════════════════════════════════════

cat("══════════════════════════════════════════════════════\n")
cat(" SECTION A — QC\n")
cat("══════════════════════════════════════════════════════\n")

# ── (A1) CpG site counts ───────────────────────────────────────────────────────

cat("\n── (A1): CpG count consistency ──\n")
cat("⚠️  Skipping A1: count_hmc.tsv and count_mc.tsv not yet generated.\n")
cat("   Uncomment once all 3 count files are available.\n")

# ── (A2) Fraction math check ──────────────────────────────────────────────────
# sum(hmc) / sum(total_c) should equal frac_hmc — residual should be 0

cat("\n── (A2): Fraction math check (residual should be 0) ──\n")

frac_hmc_qc <- readr::read_tsv(FRAC_HMC,    show_col_types = FALSE)
frac_mc_qc  <- readr::read_tsv(FRAC_MC,     show_col_types = FALSE)
sum_hmc     <- readr::read_tsv(SUM_HMC,     show_col_types = FALSE)
sum_mc      <- readr::read_tsv(SUM_MC,      show_col_types = FALSE)
sum_total_c <- readr::read_tsv(SUM_TOTAL_C, show_col_types = FALSE)

sample_cols <- setdiff(colnames(sum_hmc), META_COLS)
sample_idx  <- which(colnames(sum_hmc) %in% sample_cols)

cat("Sample columns found:", length(sample_cols), "\n")
cat("Column range used:   ", sample_cols[1], "→", sample_cols[length(sample_cols)], "\n\n")

residual_hmc <- sum(
  round(sum_hmc[sample_idx] / sum_total_c[sample_idx] - frac_hmc_qc[sample_idx], 2),
  na.rm = TRUE
)
residual_mc <- sum(
  round(sum_mc[sample_idx] / sum_total_c[sample_idx] - frac_mc_qc[sample_idx], 2),
  na.rm = TRUE
)

cat("5hmC residual (should be 0):", residual_hmc, "\n")
cat("5mC  residual (should be 0):", residual_mc,  "\n")

if (residual_hmc == 0 && residual_mc == 0) {
  cat("✅ A2 PASSED: Fraction math is internally consistent.\n")
} else {
  cat("❌ A2 FAILED: Residuals non-zero — check pipeline outputs.\n")
}

rm(frac_hmc_qc, frac_mc_qc, sum_hmc, sum_mc, sum_total_c)

# ── (A3) DMR file discovery ────────────────────────────────────────────────────

cat("\n── (A3): DMR file discovery ──\n")

find_dmr_file <- function(subdir, mark) {
  pattern <- paste0("DMR_", mark, "_Pre__Post.*\\.tsv$")
  hits    <- list.files(subdir, pattern = pattern, full.names = TRUE)
  if (length(hits) == 0) return(NA_character_)
  if (length(hits) >  1) {
    warning("Multiple matches for ", mark, " in ", subdir, " — using most recent.")
    hits <- sort(hits, decreasing = TRUE)
  }
  hits[1]
}

dmr_paths <- lapply(names(dmr_subdirs), function(label) {
  subdir <- dmr_subdirs[[label]]
  list(
    label = label,
    hmc   = find_dmr_file(subdir, "hmc"),
    mc    = find_dmr_file(subdir, "mc")
  )
})
names(dmr_paths) <- names(dmr_subdirs)

all_found <- TRUE
for (label in names(dmr_paths)) {
  entry  <- dmr_paths[[label]]
  hmc_ok <- !is.na(entry$hmc) && file.exists(entry$hmc)
  mc_ok  <- !is.na(entry$mc)  && file.exists(entry$mc)
  cat(sprintf("  %-20s  hmc: %s   mc: %s\n",
              label,
              ifelse(hmc_ok, "✅", "❌"),
              ifelse(mc_ok,  "✅", "❌")))
  if (!hmc_ok || !mc_ok) all_found <- FALSE
}

if (all_found) {
  cat("✅ A3 PASSED: All DMR TSV files found.\n")
} else {
  cat("❌ A3 FAILED: One or more DMR files missing.\n")
  stop("Aborting: missing DMR input files.")
}

cat("\n══════════════════════════════════════════════════════\n")
cat(" Section A complete.\n")
cat("══════════════════════════════════════════════════════\n\n")

# ═══════════════════════════════════════════════════════════════════════════════
# (B) Read and prepare DMR data
#
# For each cohort × mark:
#   1. Read DMR TSV
#   2. Rename columns to DESeq2 convention
#   3. Compute log2FoldChange from group means with epsilon guard
#      (Biomodal mod_fold_change has inf/NA when Pre=0 — recompute from means)
#   4. Apply biological filter:
#        5hmC → Gene_body only
#        5mC  → Promoter only
#      AND gene-annotated regions only (!is.na(Name), n_CpGs >= N_CPGS_CUTOFF)
#
# Result: dmr_raw[[cohort]][[mark]] — filtered DMR data, no z-scores yet
# ═══════════════════════════════════════════════════════════════════════════════

cat("══════════════════════════════════════════════════════\n")
cat(" SECTION B — Read and Prepare DMR Data\n")
cat("══════════════════════════════════════════════════════\n\n")

read_dmr <- function(tsv_path, mark) {

  annotation_type <- ifelse(mark == "hmc", "Gene_body", "Promoter")

  readr::read_tsv(tsv_path, show_col_types = FALSE) %>%
    # ── Rename to DESeq2 convention ──────────────────────────────────────────
    dplyr::rename(
      n_CpGs        = num_contexts,
      baseMean_Pre  = mean_mod_group_1,
      baseMean_Post = mean_mod_group_2,
      foldChange    = mod_fold_change,
      pvalue        = dmr_pvalue,
      padj          = dmr_qvalue,
      stat          = test_statistic
    ) %>%
    # ── Compute log2FoldChange with epsilon guard ─────────────────────────────
    # Biomodal mod_fold_change has inf/NA when Pre or Post mean = 0
    # Recompute from group means with epsilon to avoid log2(0) = -Inf
    dplyr::mutate(
      baseMean_Pre   = ifelse(baseMean_Pre  == 0 | is.na(baseMean_Pre),  1e-6, baseMean_Pre),
      baseMean_Post  = ifelse(baseMean_Post == 0 | is.na(baseMean_Post), 1e-6, baseMean_Post),
      foldChange     = baseMean_Post / baseMean_Pre,
      log2FoldChange = log2(foldChange),
      Direction      = ifelse(log2FoldChange > 0, "UP", "DOWN")
    ) %>%
    # ── Biological filter ─────────────────────────────────────────────────────
    # Gene-annotated regions only, correct feature type, minimum CpG coverage
    dplyr::filter(
      !is.na(Name), Name != "",
      Annotation == annotation_type,
      n_CpGs >= N_CPGS_CUTOFF
    )
}

dmr_raw <- lapply(COHORTS, function(cohort) {
  lapply(setNames(MARKS, MARKS), function(mark) {
    cat(sprintf("  Reading: %s / %s\n", cohort, mark))
    df <- read_dmr(dmr_paths[[cohort]][[mark]], mark)
    cat(sprintf("    Regions after biological filter: %d\n", nrow(df)))
    df
  })
})
names(dmr_raw) <- COHORTS

cat("\n✅ Section B complete.\n\n")

# ═══════════════════════════════════════════════════════════════════════════════
# (C) Compute z-scores
#
# Per-patient z-scores normalise each gene relative to its own Pre baseline,
# removing inter-patient variability. CpG bias correction removes the inflation
# of z-scores for large genes with many CpGs.
#
# Steps:
#   C1. Load fraction matrices (frac_hmc.tsv, frac_mc.tsv)
#   C2. Load metadata — identifies Pre (C1) and Post samples per cohort
#   C3. Prepare per-gene methylation matrices
#        5hmC: raw fraction space (sparse, near-zero — logit would distort)
#        5mC:  logit transform (bimodal beta → M-value space)
#        Collapse multiple windows per gene by mean
#   C4. Compute per-patient z-scores
#        z = (Post_fraction - mean(Pre)) / sd(Pre)
#        sd_floor = 5th percentile of non-zero SDs (data-driven, per mark)
#   C5. CpG bias correction
#        Regress log10(total_CpGs) out of z-scores per Post sample
#        Average corrected z across Post samples → mean_Z (GSEA ranking metric)
#
# Result: z_hmc_corrected, z_mc_corrected — Name + mean_Z per gene
# ═══════════════════════════════════════════════════════════════════════════════

cat("══════════════════════════════════════════════════════\n")
cat(" SECTION C — Compute Z-scores\n")
cat("══════════════════════════════════════════════════════\n\n")

# ── (C1) Load fraction matrices ───────────────────────────────────────────────

cat("── (C1): Loading fraction matrices ──\n")

frac_hmc_raw <- readr::read_tsv(FRAC_HMC, show_col_types = FALSE)
frac_mc_raw  <- readr::read_tsv(FRAC_MC,  show_col_types = FALSE)

# Strip assay suffix so column names match Sample_ID in metadata exactly
# e.g. "361C1_num_hmc_region_frac" → "361C1"
colnames(frac_hmc_raw) <- gsub("_num_hmc_region_frac", "", colnames(frac_hmc_raw))
colnames(frac_mc_raw)  <- gsub("_num_mc_region_frac",  "", colnames(frac_mc_raw))

cat(sprintf("  frac_hmc: %d regions, %d columns\n", nrow(frac_hmc_raw), ncol(frac_hmc_raw)))
cat(sprintf("  frac_mc:  %d regions, %d columns\n", nrow(frac_mc_raw),  ncol(frac_mc_raw)))

# ── (C2) Load metadata ────────────────────────────────────────────────────────

cat("\n── (C2): Loading metadata ──\n")

meta_all <- dplyr::bind_rows(
  lapply(names(meta_paths), function(cohort) {
    readr::read_csv(meta_paths[[cohort]], show_col_types = FALSE) %>%
      dplyr::select(Sample_ID, Patient_ID, Condition, Cycle) %>%
      dplyr::mutate(cohort = cohort)
  })
) %>%
  dplyr::distinct(Sample_ID, .keep_all = TRUE)

cat(sprintf("  Total samples: %d (Pre: %d, Post: %d)\n",
            nrow(meta_all),
            sum(meta_all$Condition == "Pre"),
            sum(meta_all$Condition == "Post")))

# ── (C3) Prepare per-gene methylation matrices ────────────────────────────────
# Filter to correct annotation type, replace NA with 0, apply logit for 5mC,
# collapse multiple windows per gene by mean

cat("\n── (C3): Preparing gene methylation matrices ──\n")

prepare_gene_matrix <- function(frac_df, annotation_type, transform_logit = FALSE) {

  sample_cols <- setdiff(colnames(frac_df), META_COLS)

  df <- frac_df %>%
    dplyr::filter(
      !is.na(Name), Name != "",
      Annotation == annotation_type
    ) %>%
    dplyr::mutate(across(all_of(sample_cols), ~ replace(.x, is.na(.x), 0)))

  # Logit transform for 5mC only — promoter methylation is bimodal/beta-distributed
  # epsilon guard prevents log2(0) and log2(Inf)
  if (transform_logit) {
    eps <- 0.001
    df <- df %>%
      dplyr::mutate(across(
        all_of(sample_cols),
        ~ log2((.x + eps) / (1 - .x + eps))
      ))
  }

  # Collapse to one row per gene by averaging across windows
  df %>%
    dplyr::group_by(Name) %>%
    dplyr::summarise(
      across(all_of(sample_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )
}

hmc_gene_mat <- prepare_gene_matrix(frac_hmc_raw, "Gene_body", transform_logit = FALSE)
mc_gene_mat  <- prepare_gene_matrix(frac_mc_raw,  "Promoter",  transform_logit = TRUE)

cat(sprintf("  hmc gene matrix: %d genes × %d samples\n",
            nrow(hmc_gene_mat), ncol(hmc_gene_mat) - 1))
cat(sprintf("  mc  gene matrix: %d genes × %d samples\n",
            nrow(mc_gene_mat),  ncol(mc_gene_mat)  - 1))

# ── (C4) Compute per-patient z-scores ─────────────────────────────────────────
# z = (Post_fraction - mean(Pre)) / sd(Pre)
# sd_floor = 5th percentile of non-zero SDs (data-driven per mark)
# Using all samples from both cohorts as baseline (maximises Pre sample n)

cat("\n── (C4): Computing per-patient z-scores ──\n")

compute_zscores <- function(gene_mat, meta_df) {

  all_samples   <- setdiff(colnames(gene_mat), "Name")
  valid_samples <- intersect(all_samples, meta_df$Sample_ID)

  pre_samples  <- meta_df %>%
    dplyr::filter(Condition == "Pre",  Sample_ID %in% valid_samples) %>%
    dplyr::pull(Sample_ID)
  post_samples <- meta_df %>%
    dplyr::filter(Condition == "Post", Sample_ID %in% valid_samples) %>%
    dplyr::pull(Sample_ID)

  cat(sprintf("    Pre samples: %d | Post samples: %d\n",
              length(pre_samples), length(post_samples)))

  pre_matrix  <- as.matrix(gene_mat %>% dplyr::select(all_of(pre_samples)))
  post_matrix <- as.matrix(gene_mat %>% dplyr::select(all_of(post_samples)))

  baseline_mean <- matrixStats::rowMeans2(pre_matrix, na.rm = TRUE)
  baseline_sd   <- matrixStats::rowSds(pre_matrix,   na.rm = TRUE)

  # Data-driven sd_floor: 5th percentile of observed non-zero SDs
  # Floors only genuinely zero-variance genes without deflating low-but-real variance
  sd_floor <- quantile(baseline_sd[baseline_sd > 0], 0.05, na.rm = TRUE)
  n_sd0    <- sum(baseline_sd == 0 | is.na(baseline_sd))
  cat(sprintf("    SD floor (5th pct of non-zero SDs): %.4f\n", sd_floor))
  cat(sprintf("    SD=0 genes floored: %d\n", n_sd0))

  baseline_sd[baseline_sd == 0 | is.na(baseline_sd)] <- sd_floor

  z_mat <- sweep(post_matrix, 1, baseline_mean, "-") %>%
    sweep(1, baseline_sd, "/")

  dplyr::bind_cols(Name = gene_mat$Name, as.data.frame(z_mat))
}

cat("  5hmC z-scores:\n")
z_hmc <- compute_zscores(hmc_gene_mat, meta_all)
cat("  5mC z-scores:\n")
z_mc  <- compute_zscores(mc_gene_mat,  meta_all)

# ── (C5) CpG bias correction ──────────────────────────────────────────────────
# Large genes accumulate more CpGs → inflated z-scores by chance
# Fix: regress log10(total_CpGs) out of each Post sample's z-score
# total_CpGs = sum of n_CpGs across all DMR regions for that gene
# Use full cohort DMR data (both cohorts combined) for most stable CpG estimates

cat("\n── (C5): CpG bias correction ──\n")

# Build total CpG count per gene from all cohort DMR data combined
build_ncpg <- function(mark) {
  dplyr::bind_rows(lapply(COHORTS, function(cohort) dmr_raw[[cohort]][[mark]])) %>%
    dplyr::group_by(Name) %>%
    dplyr::summarise(total_CpGs = sum(n_CpGs, na.rm = TRUE), .groups = "drop")
}

ncpg_hmc <- build_ncpg("hmc")
ncpg_mc  <- build_ncpg("mc")

correct_cpg_bias <- function(z_df, ncpg_df) {

  post_cols <- setdiff(colnames(z_df), "Name")

  z_corrected <- z_df %>%
    dplyr::left_join(ncpg_df, by = "Name") %>%
    dplyr::mutate(across(
      .cols = all_of(post_cols),
      .fns  = ~ {
        y    <- .x
        cpg  <- total_CpGs
        keep <- !is.na(y) & !is.na(cpg) & cpg > 0
        result <- rep(NA_real_, length(y))
        if (sum(keep) >= 3) {
          lm_fit       <- lm(y[keep] ~ log10(cpg[keep]))
          result[keep] <- residuals(lm_fit)
        }
        result
      }
    )) %>%
    dplyr::select(-total_CpGs)

  # Average corrected z across all Post samples → one score per gene
  post_vals <- z_corrected %>% dplyr::select(all_of(post_cols))
  z_corrected %>%
    dplyr::mutate(
      mean_Z = rowMeans(as.matrix(post_vals), na.rm = TRUE),
      mean_Z = ifelse(is.nan(mean_Z), NA_real_, mean_Z)
    ) %>%
    dplyr::select(Name, all_of(post_cols), mean_Z)
}

z_hmc_corrected <- correct_cpg_bias(z_hmc, ncpg_hmc)
z_mc_corrected  <- correct_cpg_bias(z_mc,  ncpg_mc)

cat(sprintf("  Corrected z_hmc: %d genes, NA mean_Z: %d\n",
            nrow(z_hmc_corrected), sum(is.na(z_hmc_corrected$mean_Z))))
cat(sprintf("  Corrected z_mc:  %d genes, NA mean_Z: %d\n",
            nrow(z_mc_corrected),  sum(is.na(z_mc_corrected$mean_Z))))

cat("\n✅ Section C complete.\n\n")

# ═══════════════════════════════════════════════════════════════════════════════
# (D) Build master DMR dataframes
#
# Join z-scores and patient consistency metrics to DMR data.
# One master per cohort × mark (4 masters total).
#
# Patient consistency:
#   For each gene, count Post patients where Z > 0 vs Z < 0.
#   percent_consistent = max(n_positive, n_negative) / n_post_patients
#   patient_consistent = "n_agree/n_total" e.g. "3/4"
#   This is pure patient-level majority voting — no reference to group mean.
#
# Master columns:
#   Genomic     : Chromosome, Start, End, Annotation, Name
#   DMR stats   : padj, pvalue, stat, n_CpGs, baseMean_Pre, baseMean_Post,
#                 foldChange, log2FoldChange, Direction
#   Z-score     : mean_Z (CpG-corrected average post-treatment z-score)
#   Consistency : percent_consistent, patient_consistent
# ═══════════════════════════════════════════════════════════════════════════════

cat("══════════════════════════════════════════════════════\n")
cat(" SECTION D — Build Master DMR Dataframes\n")
cat("══════════════════════════════════════════════════════\n\n")

build_master <- function(cohort, mark) {

  # Post sample IDs for this cohort
  post_samples <- meta_all %>%
    dplyr::filter(cohort == !!cohort, Condition == "Post") %>%
    dplyr::pull(Sample_ID)

  n_post <- length(post_samples)

  # Select correct z-score table
  z_corrected <- if (mark == "hmc") z_hmc_corrected else z_mc_corrected

  # Post-sample z-score columns available for this cohort
  z_post_cols <- intersect(post_samples, colnames(z_corrected))

  # Compute patient consistency from individual post-sample z-scores
  consistency_df <- z_corrected %>%
    dplyr::select(Name, all_of(z_post_cols)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      n_positive         = sum(dplyr::c_across(all_of(z_post_cols)) > 0, na.rm = TRUE),
      n_negative         = sum(dplyr::c_across(all_of(z_post_cols)) < 0, na.rm = TRUE),
      n_agree            = max(n_positive, n_negative),
      percent_consistent = n_agree / n_post,
      patient_consistent = paste0(n_agree, "/", n_post)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(Name, percent_consistent, patient_consistent)

  # Join z-score (mean_Z only) and consistency to DMR data
  dmr_raw[[cohort]][[mark]] %>%
    dplyr::left_join(
      z_corrected %>% dplyr::select(Name, mean_Z),
      by = "Name"
    ) %>%
    dplyr::left_join(consistency_df, by = "Name")
}

master <- lapply(COHORTS, function(cohort) {
  lapply(setNames(MARKS, MARKS), function(mark) {
    cat(sprintf("  Building master: %s / %s\n", cohort, mark))
    df <- build_master(cohort, mark)
    cat(sprintf("    Regions: %d | with mean_Z: %d\n",
                nrow(df), sum(!is.na(df$mean_Z))))
    df
  })
})
names(master) <- COHORTS

cat("\n✅ Section D complete.\n\n")

# ═══════════════════════════════════════════════════════════════════════════════
# (E) Save outputs
#
# Excel: 4 files (short/long × hmc/mc), each with 4 sheets:
#   All_regions      — full master (biological filter + z-score + consistency)
#   Hits_padj        — padj <= PADJ_CUTOFF AND n_CpGs >= N_CPGS_CUTOFF
#   Hits_padj_zscore — padj <= PADJ_CUTOFF AND n_CpGs >= N_CPGS_CUTOFF
#                      AND |mean_Z| >= ZSCORE_CUTOFF
#   Hits_padj_lfc    — padj <= PADJ_CUTOFF AND n_CpGs >= N_CPGS_CUTOFF
#                      AND |log2FoldChange| >= LFC_CUTOFF
# ═══════════════════════════════════════════════════════════════════════════════

cat("══════════════════════════════════════════════════════\n")
cat(" SECTION E — Save Outputs\n")
cat("══════════════════════════════════════════════════════\n\n")

for (cohort in COHORTS) {
  for (mark in MARKS) {

    df       <- master[[cohort]][[mark]]
    filename <- sprintf("DMR_%s_%s.xlsx", cohort, mark)
    wb       <- openxlsx::createWorkbook()

    # Sheet 1: All biologically filtered regions
    openxlsx::addWorksheet(wb, "All_regions")
    openxlsx::writeData(wb, "All_regions", df)

    # Sheet 2: Statistical significance filter
    openxlsx::addWorksheet(wb, "Hits_padj")
    openxlsx::writeData(wb, "Hits_padj",
      df %>% dplyr::filter(padj <= PADJ_CUTOFF))

    # Sheet 3: Statistical + z-score effect size filter (used for ORA1/2)
    openxlsx::addWorksheet(wb, "Hits_padj_zscore")
    openxlsx::writeData(wb, "Hits_padj_zscore",
      df %>% dplyr::filter(padj <= PADJ_CUTOFF, !is.na(mean_Z),
                           abs(mean_Z) >= ZSCORE_CUTOFF))

    # Sheet 4: Statistical + LFC effect size filter (used for Venn1)
    openxlsx::addWorksheet(wb, "Hits_padj_lfc")
    openxlsx::writeData(wb, "Hits_padj_lfc",
      df %>% dplyr::filter(padj <= PADJ_CUTOFF,
                           abs(log2FoldChange) >= LFC_CUTOFF))

    openxlsx::saveWorkbook(wb, file.path(OUTPUT_DIR, filename), overwrite = TRUE)
    cat(sprintf("  ✅ Saved: %s\n", filename))
  }
}

cat("\n✅ Section E complete.\n\n")

# ═══════════════════════════════════════════════════════════════════════════════
# ██████████████████████████████████████████████████████████████████████████████
# BLOCK 2 — ANALYSIS MODULES
# Each module uses objects from Block 1 already in memory.
# To rerun a single module interactively, run Block 1 first (~1-2 min).
# ██████████████████████████████████████████████████████████████████████████████
# ═══════════════════════════════════════════════════════════════════════════════

cat("██████████████████████████████████████████████████████\n")
cat(" BLOCK 2 — ANALYSIS MODULES\n")
cat("██████████████████████████████████████████████████████\n\n")

# ═══════════════════════════════════════════════════════════════════════════════
# (F) DMR Funnel — print filter counts from master
#
# Prints how many regions pass each successive filter:
#   All regions (biological filter already applied in Block 1)
#   → padj <= PADJ_CUTOFF
#   → padj + |mean_Z| >= ZSCORE_CUTOFF
#   → padj + |log2FoldChange| >= LFC_CUTOFF
#   → padj + |mean_Z| + percent_consistent >= CONSISTENCY_CUTOFF
# ═══════════════════════════════════════════════════════════════════════════════

cat("══════════════════════════════════════════════════════\n")
cat(" SECTION F — DMR Funnel\n")
cat("══════════════════════════════════════════════════════\n\n")

print_dmr_funnel <- function(cohort, mark) {

  df         <- master[[cohort]][[mark]]
  mark_label <- ifelse(mark == "hmc", "5hmC Gene_body", "5mC Promoter")

  cat(sprintf("\n────────────────────────────────────────────────────\n"))
  cat(sprintf(" %s — %s\n", cohort, mark_label))
  cat(sprintf("────────────────────────────────────────────────────\n"))

  step0 <- nrow(df)
  step1 <- df %>% dplyr::filter(padj <= PADJ_CUTOFF) %>% nrow()
  step2 <- df %>% dplyr::filter(padj <= PADJ_CUTOFF,
                                 !is.na(mean_Z), abs(mean_Z) >= ZSCORE_CUTOFF) %>% nrow()
  step3 <- df %>% dplyr::filter(padj <= PADJ_CUTOFF,
                                 abs(log2FoldChange) >= LFC_CUTOFF) %>% nrow()
  step4 <- df %>% dplyr::filter(padj <= PADJ_CUTOFF,
                                 !is.na(mean_Z), abs(mean_Z) >= ZSCORE_CUTOFF,
                                 percent_consistent >= CONSISTENCY_CUTOFF) %>% nrow()

  # Direction breakdown for padj + zscore hits
  hits_z <- df %>% dplyr::filter(padj <= PADJ_CUTOFF,
                                  !is.na(mean_Z), abs(mean_Z) >= ZSCORE_CUTOFF)

  cat(sprintf("All regions (biological filter):          %d\n", step0))
  cat(sprintf("padj <= %.2f:                             %d\n", PADJ_CUTOFF, step1))
  cat(sprintf("padj + |mean_Z| >= %.1f:                  %d\n", ZSCORE_CUTOFF, step2))
  cat(sprintf("  UP:                                     %d\n",
              sum(hits_z$Direction == "UP")))
  cat(sprintf("  DOWN:                                   %d\n",
              sum(hits_z$Direction == "DOWN")))
  cat(sprintf("padj + |LFC| >= %.2f:                     %d\n", LFC_CUTOFF, step3))
  cat(sprintf("padj + |mean_Z| + consistent >= %.0f%%:    %d\n",
              CONSISTENCY_CUTOFF * 100, step4))
}

for (cohort in COHORTS) {
  for (mark in MARKS) {
    print_dmr_funnel(cohort, mark)
  }
}

cat("\n══════════════════════════════════════════════════════\n")
cat(" Section F complete.\n")
cat("══════════════════════════════════════════════════════\n\n")

# ═══════════════════════════════════════════════════════════════════════════════
# (G) Venn Diagrams — Short vs Long cohort overlap
#
# 6 Venns total: 2 marks × 3 filter types
#   Venn type 1 (padj_lfc)      : padj + |LFC| >= LFC_CUTOFF
#   Venn type 2 (padj_zscore)   : padj + |mean_Z| >= ZSCORE_CUTOFF
#   Venn type 3 (patient_zscore): |mean_Z| >= ZSCORE_CUTOFF
#                                  + percent_consistent >= CONSISTENCY_CUTOFF
#
# Each Venn is 4-way: Short_UP | Short_DOWN | Long_UP | Long_DOWN
# ═══════════════════════════════════════════════════════════════════════════════

cat("══════════════════════════════════════════════════════\n")
cat(" SECTION G — Venn Diagrams\n")
cat("══════════════════════════════════════════════════════\n\n")

# Helper: extract gene names by cohort, mark, direction, filter type
get_venn_genes <- function(cohort, mark, direction, filter_type) {

  df <- master[[cohort]][[mark]]

  hits <- switch(filter_type,
    "padj_lfc" = df %>%
      dplyr::filter(padj <= PADJ_CUTOFF,
                    abs(log2FoldChange) >= LFC_CUTOFF),
    "padj_zscore" = df %>%
      dplyr::filter(padj <= PADJ_CUTOFF,
                    !is.na(mean_Z), abs(mean_Z) >= ZSCORE_CUTOFF),
    "patient_zscore" = df %>%
      dplyr::filter(!is.na(mean_Z), abs(mean_Z) >= ZSCORE_CUTOFF,
                    percent_consistent >= CONSISTENCY_CUTOFF)
  )

  hits %>%
    dplyr::filter(Direction == direction) %>%
    dplyr::pull(Name) %>%
    na.omit() %>%
    unique() %>%
    as.character()
}

# Helper: build and save one 4-way Venn
save_venn <- function(mark, filter_type) {

  mark_label   <- ifelse(mark == "hmc", "GeneBody_5hmC", "Promoter_5mC")
  venn_name    <- paste0(mark_label, "_", filter_type, "_Short_vs_Long")

  s_up   <- get_venn_genes("short_paired", mark, "UP",   filter_type)
  s_down <- get_venn_genes("short_paired", mark, "DOWN", filter_type)
  l_up   <- get_venn_genes("long_paired",  mark, "UP",   filter_type)
  l_down <- get_venn_genes("long_paired",  mark, "DOWN", filter_type)

  max_len <- max(length(s_up), length(s_down), length(l_up), length(l_down), 1)

  venn_df <- data.frame(
    Short_UP   = c(s_up,   rep(NA, max_len - length(s_up))),
    Short_DOWN = c(s_down, rep(NA, max_len - length(s_down))),
    Long_UP    = c(l_up,   rep(NA, max_len - length(l_up))),
    Long_DOWN  = c(l_down, rep(NA, max_len - length(l_down)))
  )

  plot_venn(venn_df, venn_name, OUTPUT_DIR)
  cat(sprintf("  ✅ Venn saved: %s\n", venn_name))
}

for (mark in MARKS) {
  for (filter_type in c("padj_lfc", "padj_zscore", "patient_zscore")) {
    save_venn(mark, filter_type)
  }
}

cat("\n══════════════════════════════════════════════════════\n")
cat(" Section G complete.\n")
cat("══════════════════════════════════════════════════════\n\n")

# ═══════════════════════════════════════════════════════════════════════════════
# (H) GSEA — 2 cohorts × 2 marks = 4 GSEA runs
#
# Ranking metric: CpG-corrected mean_Z
#   5hmC: mean_Z as-is       (positive = 5hmC gained post-treatment)
#   5mC:  -mean_Z (negated)  (5mC LOSS becomes positive rank = de-repression)
#
# Input: SYMBOL + score (renamed log2FoldChange for analyze_pathways detection)
#        all gene-annotated regions from master (full background, only_deg=FALSE)
# ═══════════════════════════════════════════════════════════════════════════════

cat("══════════════════════════════════════════════════════\n")
cat(" SECTION H — GSEA\n")
cat("══════════════════════════════════════════════════════\n\n")

biomodal_gsea <- function(cohort, mark) {

  mark_label  <- ifelse(mark == "hmc", "5hmC", "5mC")
  contrast    <- paste0(cohort, "_", mark_label, "_gsea")

  # For 5mC negate mean_Z so 5mC loss = positive rank = epigenetic de-repression
  sign_mult <- ifelse(mark == "hmc", 1, -1)

  gsea_input <- master[[cohort]][[mark]] %>%
    dplyr::filter(!is.na(mean_Z)) %>%
    dplyr::mutate(SYMBOL = Name) %>%
    # One entry per gene: keep highest |mean_Z| if gene appears in multiple regions
    dplyr::group_by(SYMBOL) %>%
    dplyr::slice_max(order_by = abs(mean_Z), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(log2FoldChange = sign_mult * mean_Z) %>%
    dplyr::select(SYMBOL, log2FoldChange) %>%
    dplyr::arrange(dplyr::desc(log2FoldChange))

  cat(sprintf("  %s / %s — background genes: %d\n",
              cohort, mark_label, nrow(gsea_input)))

  analyze_pathways(
    contrast   = contrast,
    input_data = gsea_input,
    only_deg   = FALSE,
    gmt_dir    = GMT_DIR,
    output_dir = OUTPUT_DIR,
    padj_cutoff = 0.05, minsize = 15, maxsize = 500
  )

  plot_pathways(
    contrast     = contrast,
    pathway_xlsx = file.path(OUTPUT_DIR, contrast, "Pathways.xlsx"),
    output_dir   = OUTPUT_DIR,
    vst_nonblind = NULL,
    metadata     = NULL
  )

  cat(sprintf("  ✅ GSEA done: %s\n", contrast))
}

for (cohort in COHORTS) {
  for (mark in MARKS) {
    biomodal_gsea(cohort, mark)
  }
}

cat("\n══════════════════════════════════════════════════════\n")
cat(" Section H complete.\n")
cat("══════════════════════════════════════════════════════\n\n")

# ═══════════════════════════════════════════════════════════════════════════════
# (I) ORA — 2 cohorts × 6 ORA types = 12 ORA runs
#
# ORA1 — 5hmC hits   : padj + n_CpGs + |mean_Z| >= 1  → SYMBOL + mean_Z
# ORA2 — 5mC hits    : padj + n_CpGs + |mean_Z| >= 1  → SYMBOL + (-mean_Z)
# ORA3 — concordant  : intersect(ORA1, ORA2) opposite directions
#                      → SYMBOL + dummy LFC (+1 reactivation, -1 silencing)
# ORA4 — 5hmC patient: |mean_Z| >= 1 + percent_consistent >= 0.75
#                      → SYMBOL + mean_Z
# ORA5 — 5mC patient : |mean_Z| >= 1 + percent_consistent >= 0.75
#                      → SYMBOL + (-mean_Z)
# ORA6 — concordant  : intersect(ORA4, ORA5) opposite directions
#                      → SYMBOL + dummy LFC (+1 reactivation, -1 silencing)
#
# All ORAs pass SYMBOL + log2FoldChange to analyze_pathways() with only_deg=TRUE
# analyze_pathways() auto-splits UP/DOWN by sign of log2FoldChange
# ═══════════════════════════════════════════════════════════════════════════════

cat("══════════════════════════════════════════════════════\n")
cat(" SECTION I — ORA\n")
cat("══════════════════════════════════════════════════════\n\n")

# Helper: run analyze_pathways + plot_pathways for one ORA contrast
biomodal_ora <- function(contrast, ora_input, cohort) {

  if (nrow(ora_input) < 5) {
    cat(sprintf("  ⚠️  %s skipped — only %d genes (need >= 5)\n",
                contrast, nrow(ora_input)))
    return(invisible(NULL))
  }

  cat(sprintf("  Running ORA: %s (%d genes)\n", contrast, nrow(ora_input)))

  analyze_pathways(
    contrast   = contrast,
    input_data = ora_input,
    only_deg   = TRUE,
    gmt_dir    = GMT_DIR,
    output_dir = OUTPUT_DIR
  )

  plot_pathways(
    contrast     = contrast,
    pathway_xlsx = file.path(OUTPUT_DIR, contrast, "Pathways.xlsx"),
    output_dir   = OUTPUT_DIR,
    vst_nonblind = NULL,
    metadata     = NULL
  )

  cat(sprintf("  ✅ ORA done: %s\n", contrast))
}

# Helper: build concordant gene input from two hit sets (ORA3 and ORA6)
# Keep only genes where the two marks show opposite directions (concordant)
# Reactivation: 5hmC UP + 5mC DOWN → dummy LFC = +1
# Silencing:    5hmC DOWN + 5mC UP → dummy LFC = -1
build_concordant_input <- function(hmc_up, hmc_down, mc_up, mc_down) {
  
  reactivation <- intersect(hmc_up,   mc_down)
  silencing    <- intersect(hmc_down, mc_up)
  
  # Guard against empty intersects
  rows <- list()
  if (length(reactivation) > 0) {
    rows[[1]] <- data.frame(SYMBOL = reactivation, log2FoldChange =  1)
  }
  if (length(silencing) > 0) {
    rows[[2]] <- data.frame(SYMBOL = silencing, log2FoldChange = -1)
  }
  
  if (length(rows) == 0) {
    return(data.frame(SYMBOL = character(0), log2FoldChange = numeric(0)))
  }
  
  dplyr::bind_rows(rows) %>%
    dplyr::distinct(SYMBOL, .keep_all = TRUE)
}

# Run all 6 ORA types per cohort
run_all_oras <- function(cohort) {

  cat(sprintf("\n────────────────────────────────────────────────────\n"))
  cat(sprintf(" ORA cohort: %s\n", cohort))
  cat(sprintf("────────────────────────────────────────────────────\n"))

  hmc_df <- master[[cohort]][["hmc"]]
  mc_df  <- master[[cohort]][["mc"]]

  # ── ORA1: 5hmC padj + z-score hits ────────────────────────────────────────
  ora1_input <- hmc_df %>%
    dplyr::filter(padj <= PADJ_CUTOFF, !is.na(mean_Z),
                  abs(mean_Z) >= ZSCORE_CUTOFF) %>%
    dplyr::mutate(SYMBOL = Name, log2FoldChange = mean_Z) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(SYMBOL, log2FoldChange)

  biomodal_ora(paste0(cohort, "_5hmC_ora1"), ora1_input, cohort)

  # ── ORA2: 5mC padj + z-score hits (negated: 5mC loss = positive) ──────────
  ora2_input <- mc_df %>%
    dplyr::filter(padj <= PADJ_CUTOFF, !is.na(mean_Z),
                  abs(mean_Z) >= ZSCORE_CUTOFF) %>%
    dplyr::mutate(SYMBOL = Name, log2FoldChange = -mean_Z) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(SYMBOL, log2FoldChange)

  biomodal_ora(paste0(cohort, "_5mC_ora2"), ora2_input, cohort)

  # ── ORA3: concordant — intersect ORA1 and ORA2 opposite directions ─────────
  ora1_up   <- ora1_input %>% dplyr::filter(log2FoldChange >  0) %>% dplyr::pull(SYMBOL)
  ora1_down <- ora1_input %>% dplyr::filter(log2FoldChange <  0) %>% dplyr::pull(SYMBOL)
  ora2_up   <- ora2_input %>% dplyr::filter(log2FoldChange >  0) %>% dplyr::pull(SYMBOL)
  ora2_down <- ora2_input %>% dplyr::filter(log2FoldChange <  0) %>% dplyr::pull(SYMBOL)
  # Note: ora2 is negated so ora2_up = 5mC lost, ora2_down = 5mC gained

  ora3_input <- build_concordant_input(ora1_up, ora1_down, ora2_up, ora2_down)
  biomodal_ora(paste0(cohort, "_concordant_ora3"), ora3_input, cohort)

  cat(sprintf("    ORA3 concordant: %d reactivation + %d silencing genes\n",
              sum(ora3_input$log2FoldChange > 0),
              sum(ora3_input$log2FoldChange < 0)))

  # ── ORA4: 5hmC patient-consistent hits ────────────────────────────────────
  ora4_input <- hmc_df %>%
    dplyr::filter(!is.na(mean_Z), abs(mean_Z) >= ZSCORE_CUTOFF,
                  percent_consistent >= CONSISTENCY_CUTOFF) %>%
    dplyr::mutate(SYMBOL = Name, log2FoldChange = mean_Z) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(SYMBOL, log2FoldChange)

  biomodal_ora(paste0(cohort, "_5hmC_ora4"), ora4_input, cohort)

  # ── ORA5: 5mC patient-consistent hits (negated) ───────────────────────────
  ora5_input <- mc_df %>%
    dplyr::filter(!is.na(mean_Z), abs(mean_Z) >= ZSCORE_CUTOFF,
                  percent_consistent >= CONSISTENCY_CUTOFF) %>%
    dplyr::mutate(SYMBOL = Name, log2FoldChange = -mean_Z) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(SYMBOL, log2FoldChange)

  biomodal_ora(paste0(cohort, "_5mC_ora5"), ora5_input, cohort)

  # ── ORA6: concordant patient-consistent ───────────────────────────────────
  ora4_up   <- ora4_input %>% dplyr::filter(log2FoldChange >  0) %>% dplyr::pull(SYMBOL)
  ora4_down <- ora4_input %>% dplyr::filter(log2FoldChange <  0) %>% dplyr::pull(SYMBOL)
  ora5_up   <- ora5_input %>% dplyr::filter(log2FoldChange >  0) %>% dplyr::pull(SYMBOL)
  ora5_down <- ora5_input %>% dplyr::filter(log2FoldChange <  0) %>% dplyr::pull(SYMBOL)

  ora6_input <- build_concordant_input(ora4_up, ora4_down, ora5_up, ora5_down)
  biomodal_ora(paste0(cohort, "_concordant_ora6"), ora6_input, cohort)

  cat(sprintf("    ORA6 concordant: %d reactivation + %d silencing genes\n",
              sum(ora6_input$log2FoldChange > 0),
              sum(ora6_input$log2FoldChange < 0)))
}

for (cohort in COHORTS) {
  run_all_oras(cohort)
}

cat("\n══════════════════════════════════════════════════════\n")
cat(" Section I complete.\n")
cat("══════════════════════════════════════════════════════\n\n")

cat("========================================================\n")
cat(sprintf("End: %s\n", Sys.time()))
cat("========================================================\n")
