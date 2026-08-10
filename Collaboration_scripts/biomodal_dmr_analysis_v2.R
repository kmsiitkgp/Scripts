# ═══════════════════════════════════════════════════════════════════════════════
# Biomodal DMR Analysis Pipeline
#
# Goal: Identify tumor-specific differentially methylated regions (DMRs)
#       in cfDNA comparing Pre vs Post treatment, using both 5hmC and 5mC marks,
#       and cross-reference with bulk RNASeq DESeq2 results.
#
# Analysis sections:
#   (A) QC           — verify Biomodal output integrity
#   (B) Normal ref   — build baseline methylation reference from healthy samples
#   (C) Group DMRs   — identify significant DMRs (Pre vs Post, group-level)
#                      output: sig_hmc, sig_mc (with Type, Direction columns)
#   (D) Sample trends— patient-wise directional consistency (individual level)
#                      output: final_hmc_out, final_mc_out
#   (E) Concordance  — genes changing consistently in both 5hmC AND 5mC
#   (F) DESeq2       — cross-reference with bulk RNASeq at group + patient level
#
# Note on treatment context:
#   Biomodal: cfDNA Pre (on Apalutamide) vs Post (combo Apa+Carotuximab)
#   DESeq2:   AC_A (combo vs Apa alone) is closest biological match
#             All 5 comparisons run for completeness
# ═══════════════════════════════════════════════════════════════════════════════

#install.packages("emmeans")
library(readr)
library(dplyr)
library(tibble)
library(openxlsx)
library(GenomicRanges)
library(rtracklayer)
library(purrr)
library(S4Vectors)
library(ggplot2)
library(ggrepel)


# ── Paths ─────────────────────────────────────────────────────────────────────
# Uncomment appropriate path for your environment

# Local (Windows)
# path    <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Biomodal"
# gmt_dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/Signatures/Human"

# HPC
path    <- "/hpc/home/kailasamms/scratch/Biomodal"
gmt_dir <- "/hpc/home/kailasamms/NGSTools/Signatures/Human/"
source("/hpc/home/kailasamms/scripts/nextflow/modules/UTILS.R")

output_dir <- path

# ── DMR timestamp — update when new modality run is produced ──────────────────
# Pattern: DMR_{TIMESTAMP}_DMR_hmc_Pre__Post_{TIMESTAMP}.tsv
DMR_TIMESTAMP <- "20260331_123140"

# ── DESeq2 file paths ─────────────────────────────────────────────────────────
RNASEQ_DIR <- "/hpc/home/kailasamms/scratch/RNASeq_Xenograft_Kevin_22Rv1_Xeno (A1,C1,U3 removed)/Human/08.DESeq2"

deseq2_paths <- list(
  AC_A = file.path(RNASEQ_DIR, "Apalutamide_Carotuximab-Apalutamide/DEGs.xlsx"),
  AC_C = file.path(RNASEQ_DIR, "Apalutamide_Carotuximab-Carotuximab/DEGs.xlsx"),
  AC_V = file.path(RNASEQ_DIR, "Apalutamide_Carotuximab-Vehicle/DEGs.xlsx")
  #A_V  = file.path(RNASEQ_DIR, "Apalutamide-Vehicle/DEGs.xlsx"),
  #C_V  = file.path(RNASEQ_DIR, "Carotuximab-Vehicle/DEGs.xlsx")
)


# ═══════════════════════════════════════════════════════════════════════════════
# (A) QC — Cross-check Biomodal output for internal consistency
# ═══════════════════════════════════════════════════════════════════════════════

# (A1) CpG site counts must be identical across all 3 modification types
#      because they all interrogate the same set of CpG positions
count_hmc     <- readr::read_tsv(file.path(path, "biomodal_output", "count_hmc.tsv"),     show_col_types = FALSE)
count_mc      <- readr::read_tsv(file.path(path, "biomodal_output", "count_mc.tsv"),      show_col_types = FALSE)
count_total_c <- readr::read_tsv(file.path(path, "biomodal_output", "count_total_c.tsv"), show_col_types = FALSE)

cat("\n── QC (A1): CpG count consistency ──\n")
cat("hmc vs mc:      ", all.equal(as.data.frame(count_hmc), as.data.frame(count_mc),      check.attributes = FALSE), "\n")
cat("hmc vs total_c: ", all.equal(as.data.frame(count_hmc), as.data.frame(count_total_c), check.attributes = FALSE), "\n")

# (A2) Verify fractions: sum(hmc) / sum(total_c) should equal frac_hmc
#      A sum of zero would indicate systematic computation errors
frac_hmc_qc <- readr::read_tsv(file.path(path, "biomodal_output", "frac_hmc.tsv"),     show_col_types = FALSE)
frac_mc_qc  <- readr::read_tsv(file.path(path, "biomodal_output", "frac_mc.tsv"),      show_col_types = FALSE)
sum_hmc     <- readr::read_tsv(file.path(path, "biomodal_output", "sum_hmc.tsv"),      show_col_types = FALSE)
sum_mc      <- readr::read_tsv(file.path(path, "biomodal_output", "sum_mc.tsv"),       show_col_types = FALSE)
sum_total_c <- readr::read_tsv(file.path(path, "biomodal_output", "sum_total_c.tsv"),  show_col_types = FALSE)

cat("\n── QC (A2): Fraction math check (should be 0) ──\n")
cat("5hmC residual:", sum(round(sum_hmc[4:24] / sum_total_c[4:24] - frac_hmc_qc[4:24], 2), na.rm = TRUE), "\n")
cat("5mC  residual:", sum(round(sum_mc[4:24]  / sum_total_c[4:24] - frac_mc_qc[4:24],  2), na.rm = TRUE), "\n")

rm(list = intersect(ls(), c("count_hmc", "count_mc", "count_total_c",
                            "frac_mc_qc", "frac_hmc_qc",
                            "sum_hmc", "sum_mc", "sum_total_c")))


# ═══════════════════════════════════════════════════════════════════════════════
# (B) Build Normal Methylation Reference
# ═══════════════════════════════════════════════════════════════════════════════
#
# Goal: identify genes with baseline 5hmC/5mC in normal samples so we can
#       subtract them from cfDNA results to isolate tumor-specific signal.
#
# 5hmC reference: benign prostate tissue (GSE144530, n=5, narrowPeak, hg19)
# 5mC  reference: healthy plasma cfDNA  (n=29, BED, hg19)
#
# Since Biomodal data is hg38 and references are hg19, direct coordinate
# comparison is not feasible. Instead we annotate both against gene models
# and compare at gene name + region level (Promoter / Gene body).
#
# Threshold: gene must appear in >=80% of normal samples (majority rule).
# Strict intersection was tested but rejected — one low-coverage 5hmC sample
# (2604 genes vs 3173-5701 in others) artificially halved the reference.

# ── (B1) Helper functions ─────────────────────────────────────────────────────

get_collapsed_genes <- function(peaks, hits, subject_gr) {
  mcols_data <- S4Vectors::mcols(subject_gr)
  col_to_use <- if ("gene_name" %in% colnames(mcols_data)) "gene_name" else "gene_id"
  hits_by_query <- S4Vectors::splitAsList(
    x = S4Vectors::subjectHits(hits),
    f = factor(S4Vectors::queryHits(hits), levels = seq_along(peaks))
  )
  vapply(hits_by_query, function(idx) {
    if (length(idx) == 0) return(NA_character_)
    raw_names    <- as.character(mcols_data[[col_to_use]][idx])
    unique_names <- unique(raw_names[!is.na(raw_names)])
    if (length(unique_names) == 0) return(NA_character_)
    paste(sort(unique_names), collapse = ";")
  }, character(1))
}

extract_unique_genes <- function(df_list, col_name, logical_col) {
  lapply(df_list, function(df) {
    mask  <- df[[logical_col]] & !is.na(df[[col_name]])
    genes <- df[[col_name]][mask]
    if (length(genes) == 0) return(character(0))
    unique(unlist(strsplit(as.character(genes), ";", fixed = TRUE)))
  })
}

get_majority_genes <- function(gene_lists, threshold = 0.8) {
  if (length(gene_lists) == 0) return(character(0))
  n         <- length(gene_lists)
  min_count <- ceiling(threshold * n)
  gene_freq <- table(unlist(gene_lists))
  if (length(gene_freq) == 0) return(character(0))
  names(gene_freq[gene_freq >= min_count])
}

# ── (B2) Load genome annotation (hg19 — to match reference files) ────────────
gtf_file      <- file.path(path, "biomodal_output", "hg19.refGene.gtf.gz")
gtf           <- rtracklayer::import(con = gtf_file, format = "gtf")
transcripts   <- gtf[S4Vectors::mcols(gtf)$type == "transcript"]
gene_body_gr  <- transcripts
promoter_gr   <- GenomicRanges::promoters(transcripts, upstream = 1000, downstream = 0)
tss_region_gr <- GenomicRanges::promoters(transcripts, upstream = 200,  downstream = 200)

# ── (B3) Annotate each normal reference file against gene regions ─────────────
annotate_peaks <- function(peak_file, format) {
  peaks      <- rtracklayer::import(con = peak_file, format = format)
  peak_names <- S4Vectors::mcols(peaks)$name
  if (is.null(peak_names) || all(is.na(peak_names))) {
    peak_names <- paste0("peak_", seq_along(peaks))
  }
  promoter_hits  <- GenomicRanges::findOverlaps(peaks, promoter_gr)
  tss_hits       <- GenomicRanges::findOverlaps(peaks, tss_region_gr)
  gene_body_hits <- GenomicRanges::findOverlaps(peaks, gene_body_gr)
  data.frame(
    peak_name            = as.character(peak_names),
    peak_chr             = as.character(GenomicRanges::seqnames(peaks)),
    peak_start           = GenomicRanges::start(peaks),
    peak_end             = GenomicRanges::end(peaks),
    in_promoter          = GenomicRanges::countOverlaps(peaks, promoter_gr)   > 0,
    in_tss               = GenomicRanges::countOverlaps(peaks, tss_region_gr) > 0,
    in_gene_body         = GenomicRanges::countOverlaps(peaks, gene_body_gr)  > 0,
    promoter_gene_names  = get_collapsed_genes(peaks, promoter_hits,  promoter_gr),
    tss_gene_names       = get_collapsed_genes(peaks, tss_hits,       tss_region_gr),
    gene_body_gene_names = get_collapsed_genes(peaks, gene_body_hits, gene_body_gr)
  )
}

normal_5hmc_files <- list.files(file.path(path, "Benign_prostate_5hmc"),
                                pattern = "\\.narrowPeak.gz$", full.names = TRUE)
normal_5mc_files  <- list.files(file.path(path, "Healthy_plasma_5mc"),
                                pattern = "\\.bed.gz$", ignore.case = TRUE, full.names = TRUE)

annotated_list <- list(
  normal_5hmc = lapply(normal_5hmc_files, function(f) annotate_peaks(f, "narrowPeak")),
  normal_5mc  = lapply(normal_5mc_files,  function(f) annotate_peaks(f, "BED"))
)

# ── (B4) Apply majority rule and save normal reference ───────────────────────
for (i in names(annotated_list)) {
  promoter_genes_list  <- extract_unique_genes(annotated_list[[i]], "promoter_gene_names",  "in_promoter")
  tss_genes_list       <- extract_unique_genes(annotated_list[[i]], "tss_gene_names",       "in_tss")
  gene_body_genes_list <- extract_unique_genes(annotated_list[[i]], "gene_body_gene_names", "in_gene_body")

  common_promoter_genes  <- get_majority_genes(promoter_genes_list)
  common_tss_genes       <- get_majority_genes(tss_genes_list)
  common_gene_body_genes <- get_majority_genes(gene_body_genes_list)

  cat("\n── Normal reference:", i, "──\n")
  cat("Promoter  | Count:", length(common_promoter_genes),  "\n")
  cat("TSS       | Count:", length(common_tss_genes),       "\n")
  cat("Gene Body | Count:", length(common_gene_body_genes), "\n")

  all_genes_df <- dplyr::bind_rows(
    data.frame(Name = common_promoter_genes,  Annotation = "Promoter"),
    data.frame(Name = common_tss_genes,       Annotation = "TSS"),
    data.frame(Name = common_gene_body_genes, Annotation = "Gene")
  ) %>% dplyr::distinct(Name, .keep_all = TRUE)

  if (nrow(all_genes_df) > 0) {
    output_file <- file.path(path, paste0(i, "_summary.xlsx"))
    openxlsx::write.xlsx(x = all_genes_df, file = output_file,
                         sheetName = "Summary", overwrite = TRUE)
    cat("✅ Saved:", output_file, "\n")
  } else {
    cat("⚠️  No genes met threshold for:", i, "\n")
  }
}

# Load saved references for downstream use
ref_hmc <- openxlsx::read.xlsx(file.path(path, "normal_5hmc_summary.xlsx"))
ref_mc  <- openxlsx::read.xlsx(file.path(path, "normal_5mc_summary.xlsx"))


# ═══════════════════════════════════════════════════════════════════════════════
# (C) Group-Level DMR Funnel — produces sig_hmc and sig_mc
# ═══════════════════════════════════════════════════════════════════════════════
#
# Sequential filters applied to raw Biomodal DMR output:
#   Step 0: Remove CpG island/shore/shelf regions (!is.na(Name))
#   Step 1: padj <= 0.05 (statistical significance, DESeq2-style naming)
#   Step 2: n_CpGs >= 10 (minimum CpG quality per region)
#   Step 3: |log2FoldChange| >= 0.58 (~1.5x fold change)
#           → Direction column added: UP / DOWN
#   Step 4: Subtract normal reference
#           → Type column added: Tumor only / Normal
#
# Column naming follows DESeq2 convention for easy cross-comparison:
#   dmr_qvalue      → padj
#   dmr_pvalue      → pvalue
#   mod_lfc         → log2FoldChange
#   mod_fold_change → foldChange
#   num_contexts    → n_CpGs
#   mean_mod_group_1→ baseMean_Pre
#   mean_mod_group_2→ baseMean_Post
#   test_statistic  → stat
#
# Normal subtraction strategy:
#   5hmC: strict (Name only) — hg19 vs hg38 annotation mismatch means region
#         labels can differ for same gene; only ~5% extra genes removed vs loose
#   5mC:  loose (Name + Annotation) — strict removed 54% of hits (too aggressive)
#         given 5mC is ubiquitous in normal cells
#         in_strict_ref flag added for sensitivity analysis
#
# Note: WBC-derived 5hmC not accounted for (benign prostate reference only).
#       Some hits may originate from leukocytes (~80-90% of cfDNA).
#       Adding blood cell 5hmC reference (future work) will reduce false positives.

# ── Load raw DMR output ───────────────────────────────────────────────────────
dmr_hmc <- readr::read_tsv(
  file.path(path, "biomodal_output",
            paste0("DMR_", DMR_TIMESTAMP, "_DMR_hmc_Pre__Post_", DMR_TIMESTAMP, ".tsv")),
  show_col_types = FALSE)

dmr_mc <- readr::read_tsv(
  file.path(path, "biomodal_output",
            paste0("DMR_", DMR_TIMESTAMP, "_DMR_mc_Pre__Post_", DMR_TIMESTAMP, ".tsv")),
  show_col_types = FALSE)

# ── Group DMR funnel function ─────────────────────────────────────────────────
run_group_funnel <- function(raw_df, ref_df, ref_join_cols, mark_label) {

  cat("\n════════════════════════════════════\n")
  cat(" Group-level DMR Funnel —", mark_label, "\n")
  cat("════════════════════════════════════\n")

  # Step 0: Remove CpG islands/shores/shelves — keep gene-annotated regions only
  step0 <- raw_df %>% dplyr::filter(!is.na(Name))
  cat("Step 0 | Gene regions only:           ", nrow(step0), "\n")

  # Step 1: Statistical significance
  step1 <- step0 %>% dplyr::filter(dmr_qvalue <= 0.05)
  cat("Step 1 | padj <= 0.05:                ", nrow(step1), "\n")

  # Step 2: Minimum CpG quality per region
  step2 <- step1 %>% dplyr::filter(num_contexts >= 10)
  cat("Step 2 | n_CpGs >= 10:                ", nrow(step2), "\n")

  # Step 3: Effect size filter with epsilon guard on zero fractions
  #         to avoid log2(0) = -Inf corrupting fold change calculations
  step2 <- step2 %>%
    dplyr::mutate(
      mean_mod_group_1 = ifelse(mean_mod_group_1 == 0, 1e-6, mean_mod_group_1),
      mean_mod_group_2 = ifelse(mean_mod_group_2 == 0, 1e-6, mean_mod_group_2),
      mod_fold_change  = mean_mod_group_2 / mean_mod_group_1,
      mod_lfc          = log2(mod_fold_change)
    )

  step3 <- step2 %>%
    dplyr::filter(abs(mod_lfc) >= 0.58) %>%
    dplyr::mutate(Direction = ifelse(mod_lfc > 0, "UP", "DOWN"))

  cat("Step 3 | |log2FoldChange| >= 0.58:    ", nrow(step3), "\n")
  cat("        UP:                            ", sum(step3$Direction == "UP"),   "\n")
  cat("        DOWN:                          ", sum(step3$Direction == "DOWN"), "\n")

  # Step 4: Subtract normal reference — label Tumor only vs Normal
  tumor_names <- step3 %>%
    dplyr::anti_join(ref_df, by = ref_join_cols) %>%
    dplyr::mutate(key = paste(Chromosome, Start, End, Name, Annotation, sep = "_")) %>%
    dplyr::pull(key)

  sig_all <- step3 %>%
    dplyr::mutate(
      key  = paste(Chromosome, Start, End, Name, Annotation, sep = "_"),
      Type = ifelse(key %in% tumor_names, "Tumor only", "Normal")
    )

  cat("Step 4 | After normal subtraction:\n")
  cat("        Tumor only:                    ", sum(sig_all$Type == "Tumor only"), "\n")
  cat("        Normal:                        ", sum(sig_all$Type == "Normal"),     "\n")

  # Rename columns to match DESeq2 convention for easy cross-comparison
  sig_all <- sig_all %>%
    dplyr::rename(
      log2FoldChange = mod_lfc,
      foldChange     = mod_fold_change,
      pvalue         = dmr_pvalue,
      padj           = dmr_qvalue,
      stat           = test_statistic,
      baseMean_Pre   = mean_mod_group_1,
      baseMean_Post  = mean_mod_group_2,
      n_CpGs         = num_contexts
    )

  # Add in_strict_ref flag for 5mC sensitivity analysis
  # (marks genes that would be removed if strict Name-only subtraction were used)
  if ("Name" %in% ref_join_cols && "Annotation" %in% ref_join_cols) {
    sig_all <- sig_all %>%
      dplyr::mutate(in_strict_ref = Name %in% ref_df$Name)
  }

  sig_all
}

# ── Run funnel for 5hmC and 5mC ──────────────────────────────────────────────

# 5hmC: strict subtraction by Name only (hg19 vs hg38 annotation mismatch)
sig_hmc <- run_group_funnel(dmr_hmc, ref_hmc,
                            ref_join_cols = "Name",
                            mark_label    = "5hmC")

# 5mC: loose subtraction by Name + Annotation (strict removes 54% — too aggressive)
sig_mc  <- run_group_funnel(dmr_mc,  ref_mc,
                            ref_join_cols = c("Name", "Annotation"),
                            mark_label    = "5mC")

# ── Save group-level outputs ──────────────────────────────────────────────────
openxlsx::write.xlsx(sig_hmc, file = file.path(path, "sig_hmc.xlsx"), overwrite = TRUE)
openxlsx::write.xlsx(sig_mc,  file = file.path(path, "sig_mc.xlsx"),  overwrite = TRUE)
cat("\n✅ Group-level DMR files saved: sig_hmc.xlsx, sig_mc.xlsx\n")


# ═══════════════════════════════════════════════════════════════════════════════
# (D) Sample-Level Trend Analysis — patient-wise directional consistency
# ═══════════════════════════════════════════════════════════════════════════════
#
# For each patient, compare post-treatment timepoints (C2, C3, EOT) to their
# own baseline (C1). Count how many patients show UP vs DOWN methylation per
# region. This adds individual-level consistency on top of group significance.
#
# A region must pass BOTH filters:
#   Group level:  padj <= 0.05, |log2FoldChange| >= 0.58 (from funnel above)
#   Sample level: |n_Diff| >= 2 (at least 2 more patients in same direction)
#
# Sample IDs follow pattern: {PatientID}{Timepoint} e.g. "283C1", "283C2"

# ── (D1) Build patient C1 → post-treatment pairs ─────────────────────────────
metadata <- openxlsx::read.xlsx(file.path(path, "biomodal_output", "Metadata.xlsx"))
samples  <- unique(metadata$Sample_ID)

combns     <- utils::combn(x = samples, m = 2)
p_controls <- c()
p_expts    <- c()

for (i in 1:ncol(combns)) {
  id_a <- gsub("C1|C2|C3|EOT", "", combns[1, i])
  id_b <- gsub("C1|C2|C3|EOT", "", combns[2, i])

  # Only pair samples from the same patient, with C1 as control
  if (id_a == id_b) {
    if (grepl("C1", combns[1, i]) & !grepl("C1", combns[2, i])) {
      p_controls <- c(p_controls, combns[1, i])
      p_expts    <- c(p_expts,    combns[2, i])
    }
  }
}

patient_comparisons <- list(control = p_controls, expt = p_expts)

# ── (D2) Load per-sample methylation fraction matrices ────────────────────────
frac_hmc_raw <- readr::read_tsv(file.path(path, "biomodal_output", "frac_hmc.tsv"), show_col_types = FALSE)
frac_mc_raw  <- readr::read_tsv(file.path(path, "biomodal_output", "frac_mc.tsv"),  show_col_types = FALSE)

# Clean column names, replace NA fractions with 0
# NA = region not covered in that sample = treated as no methylation
# Only applied to sample columns to avoid corrupting coordinate columns
process_fracs <- function(df, suffix, sample_ids) {
  colnames(df) <- gsub(suffix, "", colnames(df))
  sample_cols  <- intersect(sample_ids, colnames(df))
  df %>%
    dplyr::mutate(across(all_of(sample_cols), ~replace(.x, is.na(.x), 0))) %>%
    dplyr::mutate(n_UP = 0, n_DOWN = 0, n_EQUAL = 0)
}

frac_hmc <- process_fracs(frac_hmc_raw, "_num_hmc_region_frac", samples)
frac_mc  <- process_fracs(frac_mc_raw,  "_num_mc_region_frac",  samples)

# Keep only pairs where both samples exist in the fraction matrices
valid_rows <- (patient_comparisons$control %in% colnames(frac_mc)) &
              (patient_comparisons$expt    %in% colnames(frac_mc))

patient_comparisons <- list(
  control = patient_comparisons$control[valid_rows],
  expt    = patient_comparisons$expt[valid_rows]
)

cat("\n── Patient comparisons (C1 → post-treatment) ──\n")
print(data.frame(Control    = patient_comparisons$control,
                 Experiment = patient_comparisons$expt))

# ── (D3) Count directional changes per region across all patient pairs ────────
for (i in seq_along(patient_comparisons$control)) {
  ctrl_id <- patient_comparisons$control[i]
  expt_id <- patient_comparisons$expt[i]

  frac_hmc <- frac_hmc %>%
    dplyr::mutate(
      n_UP    = n_UP    + as.integer(.data[[expt_id]] >  .data[[ctrl_id]]),
      n_DOWN  = n_DOWN  + as.integer(.data[[expt_id]] <  .data[[ctrl_id]]),
      n_EQUAL = n_EQUAL + as.integer(.data[[expt_id]] == .data[[ctrl_id]])
    )

  frac_mc <- frac_mc %>%
    dplyr::mutate(
      n_UP    = n_UP    + as.integer(.data[[expt_id]] >  .data[[ctrl_id]]),
      n_DOWN  = n_DOWN  + as.integer(.data[[expt_id]] <  .data[[ctrl_id]]),
      n_EQUAL = n_EQUAL + as.integer(.data[[expt_id]] == .data[[ctrl_id]])
    )
}

# ── (D4) Sample-level funnel ──────────────────────────────────────────────────
#
# Step 5: |n_Diff| >= 2 — patient-level consistency filter
#         inner_join with sig ensures only group-significant regions kept
#         Type and Direction already encoded in sig from group funnel

run_sample_funnel <- function(frac_df, sig_df, mark_label) {

  cat("\n════════════════════════════════════\n")
  cat(" Sample-level Trend Funnel —", mark_label, "\n")
  cat("════════════════════════════════════\n")

  stopifnot("log2FoldChange" %in% colnames(sig_df))
  stopifnot("padj"           %in% colnames(sig_df))
  stopifnot("Type"           %in% colnames(sig_df))
  stopifnot("Direction"      %in% colnames(sig_df))

  # Join patient trend counts with group-significant DMRs
  # Type and Direction carried over from sig_df
  base <- frac_df %>%
    dplyr::filter(!is.na(Name)) %>%
    dplyr::inner_join(sig_df, by = c("Chromosome", "Start", "End", "Name", "Annotation")) %>%
    dplyr::mutate(n_Diff = n_UP - n_DOWN)

  # Step 5: Patient consistency — at least 2 more patients in same direction
  final <- base %>% dplyr::filter(abs(n_Diff) >= 2)

  cat("Step 5 | Patient consistent |n_Diff| >= 2:\n")
  cat("        (input already filtered: padj<=0.05,\n")
  cat("         |log2FoldChange|>=0.58, n_CpGs>=10)\n")
  cat("        UP   Tumor only:                ", sum(final$n_Diff >  0 & final$Type == "Tumor only"), "\n")
  cat("        DOWN Tumor only:                ", sum(final$n_Diff <  0 & final$Type == "Tumor only"), "\n")
  cat("        UP   Normal:                    ", sum(final$n_Diff >  0 & final$Type == "Normal"),     "\n")
  cat("        DOWN Normal:                    ", sum(final$n_Diff <  0 & final$Type == "Normal"),     "\n")
  cat("        ─────────────────────────────────────────\n")
  cat("        TOTAL:                          ", nrow(final),            "\n")
  cat("        Regions lost to n_Diff filter:  ", nrow(base) - nrow(final), "\n")

  # Format output columns
  final %>%
    dplyr::select(Chromosome, Start, End, Name, Annotation, Type, Direction,
                  n_UP, n_DOWN, n_EQUAL, n_Diff,
                  n_CpGs, log2FoldChange, padj,
                  everything())
}

final_hmc_out <- run_sample_funnel(frac_hmc, sig_hmc, "5hmC")
final_mc_out  <- run_sample_funnel(frac_mc,  sig_mc,  "5mC")

# ── Save sample-level outputs ─────────────────────────────────────────────────
openxlsx::write.xlsx(
  list("5hmC_Trends" = final_hmc_out,
       "5mC_Trends"  = final_mc_out),
  file      = file.path(path, "Biomodal_Patient_Trend_Analysis.xlsx"),
  overwrite = TRUE
)
cat("\n✅ Patient trend analysis saved\n")

# 
# # ═══════════════════════════════════════════════════════════════════════════════
# # (E) Concordance — genes changing consistently in BOTH 5hmC AND 5mC
# # ═══════════════════════════════════════════════════════════════════════════════
# #
# # Uses patient-consistent tumor-specific genes from Step 5 (final outputs)
# # Not raw sig — ensures concordance based on most stringently filtered set
# #
# # Activated = 5hmC UP (gene body + promoter gaining hmc) + 5mC DOWN (promoter losing mc)
# #             → classic gene reactivation signature post-treatment
# # Repressed  = 5hmC DOWN + 5mC UP
# #             → gene silencing signature
# 
# cat("\n════════════════════════════════════\n")
# cat(" Concordance — 5hmC + 5mC\n")
# cat("════════════════════════════════════\n")
# 
# hmc_up_genes   <- final_hmc_out %>% dplyr::filter(Type == "Tumor only", log2FoldChange > 0) %>% dplyr::pull(Name)
# hmc_down_genes <- final_hmc_out %>% dplyr::filter(Type == "Tumor only", log2FoldChange < 0) %>% dplyr::pull(Name)
# mc_up_genes    <- final_mc_out  %>% dplyr::filter(Type == "Tumor only", log2FoldChange > 0, Annotation == "Promoter") %>% dplyr::pull(Name)
# mc_down_genes  <- final_mc_out  %>% dplyr::filter(Type == "Tumor only", log2FoldChange < 0, Annotation == "Promoter") %>% dplyr::pull(Name)
# 
# activated_genes <- intersect(hmc_up_genes,   mc_down_genes)
# repressed_genes <- intersect(hmc_down_genes, mc_up_genes)
# 
# cat("Concordant activated (5hmC↑ + 5mC↓):", length(activated_genes), "\n")
# cat("Concordant repressed (5hmC↓ + 5mC↑):", length(repressed_genes), "\n")
# 
# concordant_df <- bind_rows(
#   tibble(Gene = activated_genes, Direction = "Activated"),
#   tibble(Gene = repressed_genes, Direction = "Repressed")
# )
# 
# openxlsx::write.xlsx(concordant_df,
#                      file      = file.path(output_dir, "Concordant_genes_hmc_mc.xlsx"),
#                      sheetName = "Genes",
#                      overwrite = TRUE)
# cat("✅ Concordant genes saved\n")
# 
# 
# # ═══════════════════════════════════════════════════════════════════════════════
# # (F) DESeq2 Cross-reference
# # ═══════════════════════════════════════════════════════════════════════════════
# #
# # Compare Biomodal DMRs against bulk RNASeq DESeq2 results at two levels:
# #   (F1) Group level  — sig_hmc / sig_mc (all significant DMRs)
# #   (F2) Sample level — final_hmc_out / final_mc_out (patient consistent)
# #
# # Treatment context:
# #   Biomodal Pre  = patients on Apalutamide (baseline)
# #   Biomodal Post = patients on Apa + Carotuximab (combo)
# #   → AC_A (combo vs Apa alone) is closest biological match
# #   → All 5 comparisons run for completeness
# #
# # Overlap logic:
# #   5hmC UP   + RNA UP   → active transcription (gene activation)
# #   5hmC DOWN + RNA DOWN → transcriptional repression
# #   5mC DOWN  + RNA UP   → promoter demethylation → gene reactivation (KEY)
# #   5mC UP    + RNA DOWN → promoter methylation → gene silencing
# #
# # Note: 5mC comparison uses Promoter regions only (directly controls expression)
# #       Gene_body 5mC has less clear functional link to transcription
# 
# # ── (F1) Load DESeq2 results ──────────────────────────────────────────────────
# deseq2_dfs <- lapply(deseq2_paths, openxlsx::read.xlsx)
# 
# # ── (F2) Comparison function ──────────────────────────────────────────────────
# compare_biomodal_deseq2 <- function(biomodal_df,
#                                     deseq2_df,
#                                     comparison_label,
#                                     biomodal_label,
#                                     promoter_only = FALSE) {
#   
#   # Infer mark from biomodal_label — no need to pass separately
#   mark <- ifelse(grepl("5mC", biomodal_label) & !grepl("5hmC", biomodal_label), "5mC", "5hmC")
#   
#   cat("\n── Biomodal", biomodal_label, "vs DESeq2:", comparison_label, "──\n")
#   
#   # DESeq2 gene sets
#   deg      <- deseq2_df %>% dplyr::filter(padj <= 0.05, abs(log2FoldChange) >= 0.58)
#   deg_up   <- deg %>% dplyr::filter(log2FoldChange >  0) %>% dplyr::pull(SYMBOL)
#   deg_down <- deg %>% dplyr::filter(log2FoldChange <  0) %>% dplyr::pull(SYMBOL)
#   
#   # Biomodal gene sets — tumor specific only
#   b_df <- biomodal_df %>% dplyr::filter(Type == "Tumor only")
#   if (promoter_only) b_df <- b_df %>% dplyr::filter(Annotation == "Promoter")
#   
#   b_up   <- b_df %>% dplyr::filter(Direction == "UP")   %>% dplyr::pull(Name) %>% unique()
#   b_down <- b_df %>% dplyr::filter(Direction == "DOWN") %>% dplyr::pull(Name) %>% unique()
#   
#   # Overlaps — same 4 combinations, but labels reflect biological meaning per mark
#   results <- list(
#     label          = comparison_label,
#     biomodal_label = biomodal_label,
#     up_rna_up      = intersect(b_up,   deg_up),
#     down_rna_down  = intersect(b_down, deg_down),
#     down_rna_up    = intersect(b_down, deg_up),
#     up_rna_down    = intersect(b_up,   deg_down)
#   )
#   
#   if (mark == "5hmC") {
#     # 5hmC: UP mark + UP RNA = concordant (active transcription)
#     #       DOWN mark + DOWN RNA = concordant (loss of active mark)
#     cat(sprintf("  %-50s %d\n", "5hmC UP   + RNA UP   (concordant — activation):",    length(results$up_rna_up)))
#     cat(sprintf("  %-50s %d\n", "5hmC DOWN + RNA DOWN (concordant — repression):",    length(results$down_rna_down)))
#     cat(sprintf("  %-50s %d\n", "5hmC UP   + RNA DOWN (discordant):",                 length(results$up_rna_down)))
#     cat(sprintf("  %-50s %d\n", "5hmC DOWN + RNA UP   (discordant):",                 length(results$down_rna_up)))
#   } else {
#     # 5mC: DOWN mark + UP RNA = concordant (promoter demethylation → gene ON)
#     #      UP mark + DOWN RNA = concordant (promoter methylation → gene OFF)
#     cat(sprintf("  %-50s %d\n", "5mC DOWN + RNA UP   (concordant — reactivation):",   length(results$down_rna_up)))
#     cat(sprintf("  %-50s %d\n", "5mC UP   + RNA DOWN (concordant — silencing):",      length(results$up_rna_down)))
#     cat(sprintf("  %-50s %d\n", "5mC DOWN + RNA DOWN (discordant):",                  length(results$down_rna_down)))
#     cat(sprintf("  %-50s %d\n", "5mC UP   + RNA UP   (discordant):",                  length(results$up_rna_up)))
#   }
#   
#   results
# }
# 
# # ── (F3) Run comparisons at GROUP level ──────────────────────────────────────
# cat("\n════════════════════════════════════════════════════════\n")
# cat(" (F1) GROUP LEVEL — sig_hmc / sig_mc vs DESeq2\n")
# cat("════════════════════════════════════════════════════════\n")
# 
# group_overlaps_hmc <- lapply(names(deseq2_dfs), function(name) {
#   compare_biomodal_deseq2(
#     biomodal_df      = sig_hmc,
#     deseq2_df        = deseq2_dfs[[name]],
#     comparison_label = name,
#     biomodal_label   = "5hmC",
#     promoter_only    = FALSE   # 5hmC is gene body mark
#   )
# })
# names(group_overlaps_hmc) <- names(deseq2_dfs)
# 
# group_overlaps_mc <- lapply(names(deseq2_dfs), function(name) {
#   compare_biomodal_deseq2(
#     biomodal_df      = sig_mc,
#     deseq2_df        = deseq2_dfs[[name]],
#     comparison_label = name,
#     biomodal_label   = "5mC",
#     promoter_only    = TRUE    # 5mC Promoter only for functional relevance
#   )
# })
# names(group_overlaps_mc) <- names(deseq2_dfs)
# 
# # ── (F4) Run comparisons at SAMPLE (patient-consistent) level ─────────────────
# cat("\n════════════════════════════════════════════════════════\n")
# cat(" (F2) SAMPLE LEVEL — final_hmc_out / final_mc_out vs DESeq2\n")
# cat("════════════════════════════════════════════════════════\n")
# 
# sample_overlaps_hmc <- lapply(names(deseq2_dfs), function(name) {
#   compare_biomodal_deseq2(
#     biomodal_df      = final_hmc_out,
#     deseq2_df        = deseq2_dfs[[name]],
#     comparison_label = name,
#     biomodal_label   = "5hmC (patient consistent)",
#     promoter_only    = FALSE
#   )
# })
# names(sample_overlaps_hmc) <- names(deseq2_dfs)
# 
# sample_overlaps_mc <- lapply(names(deseq2_dfs), function(name) {
#   compare_biomodal_deseq2(
#     biomodal_df      = final_mc_out,
#     deseq2_df        = deseq2_dfs[[name]],
#     comparison_label = name,
#     biomodal_label   = "5mC (patient consistent)",
#     promoter_only    = TRUE
#   )
# })
# names(sample_overlaps_mc) <- names(deseq2_dfs)
# 
# # ── (F5) Save overlap gene lists ──────────────────────────────────────────────
# # Save key finding: 5mC DOWN at Promoter + RNA UP (reactivation)
# # for each DESeq2 comparison at both group and sample level
# 
# reactivation_genes <- dplyr::bind_rows(lapply(names(deseq2_dfs), function(name) {
#   dplyr::bind_rows(
#     data.frame(
#       Gene       = group_overlaps_mc[[name]]$down_rna_up,
#       Level      = "Group",
#       Comparison = name,
#       Finding    = "5mC_DOWN_RNA_UP"
#     ),
#     data.frame(
#       Gene       = sample_overlaps_mc[[name]]$down_rna_up,
#       Level      = "Sample",
#       Comparison = name,
#       Finding    = "5mC_DOWN_RNA_UP"
#     )
#   )
# }))
# 
# openxlsx::write.xlsx(reactivation_genes,
#                      file      = file.path(output_dir, "Reactivation_genes_5mC_RNA.xlsx"),
#                      sheetName = "Genes",
#                      overwrite = TRUE)
# cat("\n✅ DESeq2 overlap results saved\n")
# 
# 
# # ═══════════════════════════════════════════════════════════════════════════════
# # Final Summary
# # ═══════════════════════════════════════════════════════════════════════════════
# 
# cat("\n\n════════════════════════════════════════════════════════\n")
# cat(" FINAL SUMMARY\n")
# cat("════════════════════════════════════════════════════════\n")
# cat(sprintf("%-45s %8s %8s\n", "Step", "5hmC", "5mC"))
# cat(sprintf("%-45s %8d %8d\n", "Raw DMR regions",
#             nrow(dmr_hmc), nrow(dmr_mc)))
# cat(sprintf("%-45s %8d %8d\n", "Step 0: Gene regions only",
#             nrow(dmr_hmc %>% dplyr::filter(!is.na(Name))),
#             nrow(dmr_mc  %>% dplyr::filter(!is.na(Name)))))
# cat(sprintf("%-45s %8d %8d\n", "Step 1: padj <= 0.05",
#             nrow(sig_hmc),
#             nrow(sig_mc)))
# cat(sprintf("%-45s %8d %8d\n", "Step 3: Tumor only",
#             sum(sig_hmc$Type == "Tumor only"),
#             sum(sig_mc$Type  == "Tumor only")))
# cat(sprintf("%-45s %8d %8d\n", "Step 5: Patient |n_Diff| >= 2",
#             nrow(final_hmc_out),
#             nrow(final_mc_out)))
# cat(sprintf("%-45s %8d %8d\n", "Concordant activated (5hmC↑ + 5mC↓)",
#             length(activated_genes), length(activated_genes)))
# cat(sprintf("%-45s %8d %8d\n", "Concordant repressed (5hmC↓ + 5mC↑)",
#             length(repressed_genes), length(repressed_genes)))
# cat("════════════════════════════════════════════════════════\n")
# cat("\nOutput files:\n")
# cat("  sig_hmc.xlsx\n")
# cat("  sig_mc.xlsx\n")
# cat("  Biomodal_Patient_Trend_Analysis.xlsx\n")
# cat("  Concordant_genes_hmc_mc.xlsx\n")
# cat("  Reactivation_genes_5mC_RNA.xlsx\n")

# ═══════════════════════════════════════════════════════════════════════════════
# cfDNA Epigenomic Analysis: 5mC + 5hmC → Composite Z → Quadrant Plot
# Flow: Diagnostics → Hypergeometric → Composite Z → Quadrant Plot
# ═══════════════════════════════════════════════════════════════════════════════

# ── Define sample groups ──────────────────────────────────────────────────────
c1_samples  <- c("361", "663", "676", "657", "700", "704", "692", "283", "772")
c2_samples  <- c("361", "676", "700", "704", "692", "283", "772", "557")
c3_samples  <- c("692", "557")
eot_samples <- c("663", "657")

# ── Helper: get columns that exist in df ──────────────────────────────────────
get_cols <- function(df, patients, suffix) {
  intersect(paste0(patients, suffix), colnames(df))
}

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1: DIAGNOSTICS
# Purpose: justify min_c1 = 6 coverage threshold and sd_floor choices
# ═══════════════════════════════════════════════════════════════════════════════

frac_hmc_raw <- readr::read_tsv(file.path(path, "biomodal_output", "frac_hmc.tsv"), show_col_types = FALSE)
frac_mc_raw  <- readr::read_tsv(file.path(path, "biomodal_output", "frac_mc.tsv"),  show_col_types = FALSE)

c1_cols_hmc <- grep("C1_", colnames(frac_hmc_raw), value = TRUE)
c1_cols_mc  <- grep("C1_", colnames(frac_mc_raw),  value = TRUE)

# ── 1a: Windows per gene (confirms averaging is needed) ───────────────────────
cat("\n=== 5hmC: Gene_body windows per gene ===\n")
frac_hmc_raw %>%
  dplyr::filter(Annotation == "Gene_body", !is.na(Name), Name != "") %>%
  dplyr::count(Name) %>%
  dplyr::count(n, name = "num_genes") %>%
  print()

cat("\n=== 5mC: Promoter windows per gene ===\n")
frac_mc_raw %>%
  dplyr::filter(Annotation == "Promoter", !is.na(Name), Name != "") %>%
  dplyr::count(Name) %>%
  dplyr::count(n, name = "num_genes") %>%
  print()

# ── 1b: Coverage threshold analysis (justifies min_c1 = 6) ───────────────────
check_coverage <- function(df, c1_cols, annotation, mark_name) {
  cat("\n===", mark_name, annotation, "coverage threshold analysis ===\n")
  
  df_filt <- df %>%
    dplyr::filter(Annotation == annotation, !is.na(Name), Name != "")
  
  total_windows <- nrow(df_filt)
  total_genes   <- n_distinct(df_filt$Name)
  cat("Total windows:", total_windows, "\n")
  cat("Total genes:  ", total_genes,   "\n\n")
  
  df_filt <- df_filt %>%
    dplyr::mutate(
      n_c1_covered = rowSums(
        !is.na(dplyr::select(., all_of(c1_cols))) &
          dplyr::select(., all_of(c1_cols)) != ""
      )
    )
  
  cat("Distribution of C1 samples covered per window:\n")
  print(table(df_filt$n_c1_covered))
  
  cat("\nWindows and genes surviving each threshold:\n")
  for (thresh in 1:length(c1_cols)) {
    surv <- df_filt %>% dplyr::filter(n_c1_covered >= thresh)
    cat(sprintf("  >= %d C1 samples: %d windows (%d%%), %d genes (%d%%)\n",
                thresh,
                nrow(surv),
                round(100 * nrow(surv) / total_windows),
                n_distinct(surv$Name),
                round(100 * n_distinct(surv$Name) / total_genes)))
  }
}

check_coverage(frac_hmc_raw, c1_cols_hmc, "Gene_body", "5hmC")
check_coverage(frac_mc_raw,  c1_cols_mc,  "Promoter",  "5mC")
# → Inspect output and confirm min_c1 = 6 is appropriate before proceeding


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2: HYPERGEOMETRIC TEST
# Purpose: validate that 5hmC (not 5mC) signal reflects tumor transcription
# Result: justifies using only tumor_hmc genes for quadrant plot filtering
# ═══════════════════════════════════════════════════════════════════════════════

# ── Load significant cfDNA gene sets (output from upstream analysis) ──────────
sig_hmc <- read.xlsx(file.path(output_dir, "sig_hmc.xlsx"))
sig_mc  <- read.xlsx(file.path(output_dir, "sig_mc.xlsx"))

# ── Load RNA-seq DESeq2 results ───────────────────────────────────────────────
deseq2_dfs <- lapply(deseq2_paths, openxlsx::read.xlsx)

# ── Define universes: genes measured in BOTH platforms ───────────────────────
universe_mc  <- intersect(unique(sig_mc$Name),  unique(deseq2_dfs[["AC_A"]]$SYMBOL))
universe_hmc <- intersect(unique(sig_hmc$Name), unique(deseq2_dfs[["AC_A"]]$SYMBOL))

cat("\nUniverse (5mC  ∩ RNA-seq):", length(universe_mc),  "\n")
cat("Universe (5hmC ∩ RNA-seq):", length(universe_hmc), "\n")

# ── Define gene sets ──────────────────────────────────────────────────────────

final_mc_out  <- read.xlsx(file.path(output_dir, "Biomodal_Patient_Trend_Analysis.xlsx"), sheet = "5mC_Trends")
final_hmc_out <- read.xlsx(file.path(output_dir, "Biomodal_Patient_Trend_Analysis.xlsx"), sheet = "5hmC_Trends")

# 5mC promoter
mc_down_promoter_all   <- final_mc_out %>% filter(Direction == "DOWN", Annotation == "Promoter")               %>% pull(Name) %>% unique()
mc_down_promoter_tumor <- final_mc_out %>% filter(Direction == "DOWN", Annotation == "Promoter", Type == "Tumor only") %>% pull(Name) %>% unique()
mc_up_promoter_all     <- final_mc_out %>% filter(Direction == "UP",   Annotation == "Promoter")               %>% pull(Name) %>% unique()
mc_up_promoter_tumor   <- final_mc_out %>% filter(Direction == "UP",   Annotation == "Promoter", Type == "Tumor only") %>% pull(Name) %>% unique()

# 5hmC (all annotations)
hmc_up_tumor   <- final_hmc_out %>% filter(Direction == "UP",   Type == "Tumor only") %>% pull(Name) %>% unique()
hmc_down_tumor <- final_hmc_out %>% filter(Direction == "DOWN", Type == "Tumor only") %>% pull(Name) %>% unique()
hmc_up_all     <- final_hmc_out %>% filter(Direction == "UP")                         %>% pull(Name) %>% unique()
hmc_down_all   <- final_hmc_out %>% filter(Direction == "DOWN")                       %>% pull(Name) %>% unique()

# 5hmC (gene body only)
hmc_up_genebody_tumor   <- final_hmc_out %>% filter(Direction == "UP",   Annotation == "Gene_body", Type == "Tumor only") %>% pull(Name) %>% unique()
hmc_down_genebody_tumor <- final_hmc_out %>% filter(Direction == "DOWN", Annotation == "Gene_body", Type == "Tumor only") %>% pull(Name) %>% unique()
hmc_up_genebody_all     <- final_hmc_out %>% filter(Direction == "UP",   Annotation == "Gene_body")                       %>% pull(Name) %>% unique()
hmc_down_genebody_all   <- final_hmc_out %>% filter(Direction == "DOWN", Annotation == "Gene_body")                       %>% pull(Name) %>% unique()

# RNA-seq thresholds
rna_up     <- deseq2_dfs[["AC_A"]] %>% filter(padj <= 0.05, log2FoldChange >=  0.58) %>% pull(SYMBOL) %>% unique()
rna_down   <- deseq2_dfs[["AC_A"]] %>% filter(padj <= 0.05, log2FoldChange <= -0.58) %>% pull(SYMBOL) %>% unique()
rna_up_L   <- deseq2_dfs[["AC_A"]] %>% filter(padj <= 0.05, log2FoldChange >   0)   %>% pull(SYMBOL) %>% unique()
rna_down_L <- deseq2_dfs[["AC_A"]] %>% filter(padj <= 0.05, log2FoldChange <   0)   %>% pull(SYMBOL) %>% unique()

# ── Hypergeometric test function ──────────────────────────────────────────────
run_hyper_test <- function(set_a, set_b, universe, label_a, label_b, category) {
  set_b_clean <- intersect(set_b, universe)
  set_a_clean <- intersect(set_a, universe)
  overlap     <- intersect(set_a_clean, set_b_clean)
  
  N   <- length(universe)
  m   <- length(set_b_clean)   # successes in universe
  n   <- N - m                 # failures in universe
  k   <- length(set_a_clean)   # sample size
  obs <- length(overlap)
  exp <- (k * m) / N
  
  p_val <- phyper(q = obs - 1, m = m, n = n, k = k, lower.tail = FALSE)
  
  tibble(
    Category    = category,
    Comparison  = paste(label_a, "vs", label_b),
    N_Universe  = N,
    n_SetA      = k,
    n_SetB      = m,
    Observed    = obs,
    Expected    = round(exp, 2),
    Fold_Enrich = ifelse(exp > 0, round(obs / exp, 2), NA),
    P_Value     = p_val
  )
}

# ── Run tests ─────────────────────────────────────────────────────────────────
# 5mC: all combinations (promoter, strict/lenient, tumor/all) — included to
# demonstrate null result across all comparisons
# 5hmC: concordant and discordant tests as internal specificity controls
tasks <- list(
  
  # --- 5mC Promoter — Tumor only ---
  list(a = mc_down_promoter_tumor, b = rna_up_L,   u = universe_mc, la = "5mC DOWN (Prom Tumor)", lb = "RNA UP (Lenient)",   cat = "5mC Prom_T_Lenient"),
  list(a = mc_up_promoter_tumor,   b = rna_down_L, u = universe_mc, la = "5mC UP (Prom Tumor)",   lb = "RNA DOWN (Lenient)", cat = "5mC Prom_T_Lenient"),
  list(a = mc_down_promoter_tumor, b = rna_up,     u = universe_mc, la = "5mC DOWN (Prom Tumor)", lb = "RNA UP (Strict)",    cat = "5mC Prom_T_Strict"),
  list(a = mc_up_promoter_tumor,   b = rna_down,   u = universe_mc, la = "5mC UP (Prom Tumor)",   lb = "RNA DOWN (Strict)",  cat = "5mC Prom_T_Strict"),
  
  # # --- 5mC Promoter — Tumor + Normal ---
  # list(a = mc_down_promoter_all,   b = rna_up_L,   u = universe_mc, la = "5mC DOWN (Prom All)",   lb = "RNA UP (Lenient)",   cat = "5mC Prom_TN_Lenient"),
  # list(a = mc_up_promoter_all,     b = rna_down_L, u = universe_mc, la = "5mC UP (Prom All)",     lb = "RNA DOWN (Lenient)", cat = "5mC Prom_TN_Lenient"),
  # list(a = mc_down_promoter_all,   b = rna_up,     u = universe_mc, la = "5mC DOWN (Prom All)",   lb = "RNA UP (Strict)",    cat = "5mC Prom_TN_Strict"),
  # list(a = mc_up_promoter_all,     b = rna_down,   u = universe_mc, la = "5mC UP (Prom All)",     lb = "RNA DOWN (Strict)",  cat = "5mC Prom_TN_Strict"),
  
  # --- 5hmC All Annotations — Tumor only ---
  list(a = hmc_up_tumor,   b = rna_up,     u = universe_hmc, la = "5hmC UP (Any Tumor)",   lb = "RNA UP (Strict)",    cat = "5hmC T_Strict"),
  #list(a = hmc_down_tumor, b = rna_down,   u = universe_hmc, la = "5hmC DOWN (Any Tumor)", lb = "RNA DOWN (Strict)",  cat = "5hmC T_Strict"),
  list(a = hmc_up_tumor,   b = rna_up_L,   u = universe_hmc, la = "5hmC UP (Any Tumor)",   lb = "RNA UP (Lenient)",   cat = "5hmC T_Lenient"),
  #list(a = hmc_down_tumor, b = rna_down_L, u = universe_hmc, la = "5hmC DOWN (Any Tumor)", lb = "RNA DOWN (Lenient)", cat = "5hmC T_Lenient"),
  
  # # --- 5hmC All Annotations — Tumor + Normal ---
  # list(a = hmc_up_all,   b = rna_up,     u = universe_hmc, la = "5hmC UP (Any All)",   lb = "RNA UP (Strict)",    cat = "5hmC TN_Strict"),
  # list(a = hmc_down_all, b = rna_down,   u = universe_hmc, la = "5hmC DOWN (Any All)", lb = "RNA DOWN (Strict)",  cat = "5hmC TN_Strict"),
  # list(a = hmc_up_all,   b = rna_up_L,   u = universe_hmc, la = "5hmC UP (Any All)",   lb = "RNA UP (Lenient)",   cat = "5hmC TN_Lenient"),
  # list(a = hmc_down_all, b = rna_down_L, u = universe_hmc, la = "5hmC DOWN (Any All)", lb = "RNA DOWN (Lenient)", cat = "5hmC TN_Lenient"),
  # 
  # # --- 5hmC Gene Body Only — Tumor only ---
  # list(a = hmc_up_genebody_tumor,   b = rna_up,     u = universe_hmc, la = "5hmC GB UP (Tumor)",   lb = "RNA UP (Strict)",    cat = "5hmC GB_T_Strict"),
  # list(a = hmc_down_genebody_tumor, b = rna_down,   u = universe_hmc, la = "5hmC GB DOWN (Tumor)", lb = "RNA DOWN (Strict)",  cat = "5hmC GB_T_Strict"),
  # list(a = hmc_up_genebody_tumor,   b = rna_up_L,   u = universe_hmc, la = "5hmC GB UP (Tumor)",   lb = "RNA UP (Lenient)",   cat = "5hmC GB_T_Lenient"),
  # list(a = hmc_down_genebody_tumor, b = rna_down_L, u = universe_hmc, la = "5hmC GB DOWN (Tumor)", lb = "RNA DOWN (Lenient)", cat = "5hmC GB_T_Lenient"),
  # 
  # # --- 5hmC Gene Body Only — Tumor + Normal ---
  # list(a = hmc_up_genebody_all,   b = rna_up,     u = universe_hmc, la = "5hmC GB UP (All)",   lb = "RNA UP (Strict)",    cat = "5hmC GB_TN_Strict"),
  # list(a = hmc_down_genebody_all, b = rna_down,   u = universe_hmc, la = "5hmC GB DOWN (All)", lb = "RNA DOWN (Strict)",  cat = "5hmC GB_TN_Strict"),
  # list(a = hmc_up_genebody_all,   b = rna_up_L,   u = universe_hmc, la = "5hmC GB UP (All)",   lb = "RNA UP (Lenient)",   cat = "5hmC GB_TN_Lenient"),
  # list(a = hmc_down_genebody_all, b = rna_down_L, u = universe_hmc, la = "5hmC GB DOWN (All)", lb = "RNA DOWN (Lenient)", cat = "5hmC GB_TN_Lenient"),
  
  # --- Discordant controls (specificity check) ---
  list(a = hmc_up_tumor, b = rna_down,   u = universe_hmc, la = "5hmC UP (Any Tumor)", lb = "RNA DOWN (Strict)",  cat = "5hmC T_Strict(Discordant)"),
  list(a = hmc_up_tumor, b = rna_down_L, u = universe_hmc, la = "5hmC UP (Any Tumor)", lb = "RNA DOWN (Lenient)", cat = "5hmC T_Lenient(Discordant)")
)

results_summary <- purrr::map_df(tasks, ~run_hyper_test(.x$a, .x$b, .x$u, .x$la, .x$lb, .x$cat)) %>%
  mutate(
    FDR = p.adjust(P_Value, method = "BH"),
    Sig = case_when(FDR < 0.001 ~ "***", FDR < 0.01 ~ "**", FDR < 0.05 ~ "*", TRUE ~ "ns")
  )

cat("\n=== Hypergeometric Test Results ===\n")
results_summary %>%
  arrange(P_Value) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3: COMPOSITE Z-SCORE
# Formula: composite_Z = (Z_hmc - Z_mc) / 2
# High score = 5hmC gained AND 5mC lost = epigenetic activation
# ═══════════════════════════════════════════════════════════════════════════════

# ── Prepare marks ─────────────────────────────────────────────────────────────
# min_c1 = 6: justified by coverage diagnostic above (Step 1b)
# transform = FALSE for 5hmC: sparse, near-zero, gamma-like distribution
# transform = TRUE  for 5mC:  bimodal beta → logit (M-value) space
prepare_mark <- function(df, annotation, sample_suffix, min_c1 = 6, transform = FALSE) {
  
  c1_cols   <- get_cols(df, c1_samples, paste0("C1_", sample_suffix))
  meta_cols <- c("Chromosome", "Start", "End", "Annotation", "Name", "n_c1_covered")
  
  df <- df %>%
    dplyr::filter(Annotation == annotation, !is.na(Name), Name != "") %>%
    dplyr::mutate(
      n_c1_covered = rowSums(
        !is.na(dplyr::select(., all_of(c1_cols))) &
          dplyr::select(., all_of(c1_cols)) != ""
      )
    ) %>%
    dplyr::filter(n_c1_covered >= min_c1)
  
  sample_cols <- setdiff(colnames(df), meta_cols)
  
  df <- df %>%
    dplyr::mutate(across(all_of(sample_cols), ~ as.numeric(.)))
  
  if (transform) {
    df <- df %>%
      dplyr::mutate(across(all_of(sample_cols),
                           ~ log2((. + 0.001) / (1 - . + 0.001))))
  }
  
  # Average across multi-window genes
  df %>%
    dplyr::group_by(Name) %>%
    dplyr::summarise(
      across(all_of(sample_cols), ~ {
        val <- mean(., na.rm = TRUE)
        if (is.nan(val)) NA_real_ else val
      }),
      .groups = "drop"
    )
}

hmc_gene <- prepare_mark(frac_hmc_raw, "Gene_body", "num_hmc_region_frac", transform = FALSE)
mc_gene  <- prepare_mark(frac_mc_raw,  "Promoter",  "num_mc_region_frac",  transform = TRUE)

cat("\nNaN in hmc_gene:", sum(sapply(hmc_gene, function(x) sum(is.nan(x)))), "\n")
cat("NaN in mc_gene: ", sum(sapply(mc_gene,  function(x) sum(is.nan(x)))), "\n")

# ── Define column groups ──────────────────────────────────────────────────────
hmc_c1_cols   <- get_cols(hmc_gene, c1_samples,  "C1_num_hmc_region_frac")
hmc_c2_cols   <- get_cols(hmc_gene, c2_samples,  "C2_num_hmc_region_frac")
hmc_c3_cols   <- get_cols(hmc_gene, c3_samples,  "C3_num_hmc_region_frac")
hmc_eot_cols  <- get_cols(hmc_gene, eot_samples, "EOT_num_hmc_region_frac")
hmc_post_cols <- c(hmc_c2_cols, hmc_c3_cols, hmc_eot_cols)

mc_c1_cols    <- get_cols(mc_gene, c1_samples,  "C1_num_mc_region_frac")
mc_c2_cols    <- get_cols(mc_gene, c2_samples,  "C2_num_mc_region_frac")
mc_c3_cols    <- get_cols(mc_gene, c3_samples,  "C3_num_mc_region_frac")
mc_eot_cols   <- get_cols(mc_gene, eot_samples, "EOT_num_mc_region_frac")
mc_post_cols  <- c(mc_c2_cols, mc_c3_cols, mc_eot_cols)

# ── Identify SD=0 genes (justifies sd_floor values) ──────────────────────────
# sd_floor = 0.01 for 5hmC: many genes have zero baseline signal (fraction space)
# sd_floor = 0.5  for 5mC:  logit/M-value space has wider dynamic range
sd0_hmc_genes <- hmc_gene$Name[apply(hmc_gene %>% select(all_of(hmc_c1_cols)), 1, sd, na.rm = TRUE) == 0]
sd0_mc_genes  <- mc_gene$Name[ apply(mc_gene  %>% select(all_of(mc_c1_cols)),  1, sd, na.rm = TRUE) == 0]

cat("\nSD=0 genes in 5hmC:", length(sd0_hmc_genes), "\n")
cat("SD=0 genes in 5mC: ", length(sd0_mc_genes),  "\n")

# ── Compute z-scores ──────────────────────────────────────────────────────────
compute_zscores <- function(gene_df, c1_cols, post_cols, sd_floor = 0.01) {
  c1_mat   <- gene_df %>% dplyr::select(all_of(c1_cols))
  post_mat <- gene_df %>% dplyr::select(all_of(post_cols))
  
  baseline_mean <- rowMeans(c1_mat, na.rm = TRUE)
  baseline_sd   <- apply(c1_mat, 1, sd, na.rm = TRUE)
  
  cat("  SD=0 genes at baseline:", sum(baseline_sd == 0, na.rm = TRUE), "\n")
  cat("  All-NA C1 genes:       ", sum(is.na(baseline_mean)), "\n")
  
  baseline_sd[baseline_sd == 0] <- sd_floor   # floor preserves de novo activation signal
  
  z_mat <- sweep(post_mat, 1, baseline_mean, "-") %>%
    sweep(1, baseline_sd, "/")
  
  bind_cols(gene = gene_df$Name, z_mat)
}

z_hmc_all <- compute_zscores(hmc_gene, hmc_c1_cols, hmc_post_cols,                sd_floor = 0.01)
z_hmc_c2  <- compute_zscores(hmc_gene, hmc_c1_cols, hmc_c2_cols,                  sd_floor = 0.01)
z_hmc_c3  <- compute_zscores(hmc_gene, hmc_c1_cols, c(hmc_c3_cols, hmc_eot_cols), sd_floor = 0.01)

z_mc_all  <- compute_zscores(mc_gene, mc_c1_cols, mc_post_cols,                   sd_floor = 0.5)
z_mc_c2   <- compute_zscores(mc_gene, mc_c1_cols, mc_c2_cols,                     sd_floor = 0.5)
z_mc_c3   <- compute_zscores(mc_gene, mc_c1_cols, c(mc_c3_cols, mc_eot_cols),     sd_floor = 0.5)

# ── Compute composite Z ───────────────────────────────────────────────────────
compute_composite <- function(z_hmc, z_mc, comparison_name,
                              sd0_hmc_genes, sd0_mc_genes,
                              ncpg_genebody, ncpg_promoter) {
  
  # ── Per-sample CpG correction for 5hmC BEFORE averaging ───────────────
  z_hmc_corrected <- z_hmc %>%
    dplyr::left_join(ncpg_genebody, by = c("gene" = "Name")) %>%
    dplyr::mutate(across(
      .cols = -c(gene, total_CpGs_hmc),
      .fns  = ~ {
        y    <- .x
        cpg  <- total_CpGs_hmc
        keep <- !is.na(y) & !is.na(cpg) & cpg > 0
        result <- rep(NA_real_, length(y))
        if (sum(keep) > 2) {
          df_lm <- data.frame(y = y[keep], cpg = cpg[keep])
          result[keep] <- residuals(lm(y ~ log10(cpg), data = df_lm))
        }
        result
      }
    )) %>%
    dplyr::select(-total_CpGs_hmc)
  
  # ── Mean Z per gene (post per-sample correction) ───────────────────────
  hmc_mean <- z_hmc_corrected %>%
    dplyr::mutate(Z_hmc = {
      val <- rowMeans(dplyr::select(., -gene), na.rm = TRUE)
      ifelse(is.nan(val), NA_real_, val)
    }) %>%
    dplyr::select(gene, Z_hmc)
  
  mc_mean <- z_mc %>%
    dplyr::mutate(Z_mc = {
      val <- rowMeans(dplyr::select(., -gene), na.rm = TRUE)
      ifelse(is.nan(val), NA_real_, val)
    }) %>%
    dplyr::select(gene, Z_mc)
  
  # ── Join diagnostics ───────────────────────────────────────────────────
  cat("\n===", comparison_name, "- join diagnostics ===\n")
  cat("Genes in 5hmC but not 5mC:", nrow(dplyr::anti_join(hmc_mean, mc_mean, by = "gene")), "\n")
  cat("Genes in 5mC but not 5hmC:", nrow(dplyr::anti_join(mc_mean, hmc_mean, by = "gene")), "\n")
  
  # ── 5mC CpG correction after averaging (r = -0.05, kept for symmetry) ─
  mc_mean <- mc_mean %>%
    dplyr::left_join(ncpg_promoter, by = c("gene" = "Name")) %>%
    dplyr::mutate(Z_mc = {
      keep <- !is.na(Z_mc) & !is.na(total_CpGs_mc) & total_CpGs_mc > 0
      result <- rep(NA_real_, n())
      if (sum(keep) > 2) {
        result[keep] <- residuals(
          lm(Z_mc ~ log10(total_CpGs_mc),
             data = dplyr::pick(everything())[keep, ])
        )
      }
      result
    }) %>%
    dplyr::select(gene, Z_mc)
  
  # ── Combine with variance normalization ───────────────────────────────
  composite <- hmc_mean %>%
    dplyr::inner_join(mc_mean, by = "gene") %>%
    dplyr::mutate(
      # Variance-normalize before combining: accounts for sd(Z_hmc) ~ 2x sd(Z_mc)
      Z_hmc_scaled = Z_hmc / sd(Z_hmc, na.rm = TRUE),
      Z_mc_scaled  = Z_mc  / sd(Z_mc,  na.rm = TRUE),
      composite_Z  = Z_hmc_scaled - Z_mc_scaled,
      signal_class = dplyr::case_when(
        Z_hmc_scaled >  1 & Z_mc_scaled < -1 ~ "concordant",
        Z_hmc_scaled >  1 & Z_mc_scaled > -1 ~ "hmc_driven",
        Z_hmc_scaled <  1 & Z_mc_scaled < -1 ~ "mc_driven",
        TRUE                                  ~ "weak"),
      sd0_hmc      = gene %in% sd0_hmc_genes,
      sd0_mc       = gene %in% sd0_mc_genes,
      sd0_any      = sd0_hmc | sd0_mc
    ) %>%
    dplyr::rename(SYMBOL = gene) %>%
    dplyr::arrange(desc(composite_Z))

  cat("Genes in composite:  ", nrow(composite), "\n")
  cat("Z-score range:       ", round(min(composite$composite_Z, na.rm = TRUE), 2),
      "to", round(max(composite$composite_Z, na.rm = TRUE), 2), "\n")
  cat("NA in composite_Z:   ", sum(is.na(composite$composite_Z)), "\n")
  cat("SD=0 flagged genes:  ", sum(composite$sd0_any), "\n")
  
  composite
}

# ── Calls ──────────────────────────────────────────────────────────────────────

ncpg_genebody <- dmr_hmc %>%
  filter(!is.na(Name), Name != "", Annotation == "Gene_body") %>%
  group_by(Name) %>%
  summarise(total_CpGs_hmc = sum(num_contexts), .groups = "drop")

ncpg_promoter <- dmr_hmc %>%
  filter(!is.na(Name), Name != "", Annotation == "Promoter") %>%
  group_by(Name) %>%
  summarise(total_CpGs_mc = sum(num_contexts), .groups = "drop")

composite_all <- compute_composite(z_hmc_all, z_mc_all, "All post vs C1",
                                   sd0_hmc_genes, sd0_mc_genes,
                                   ncpg_genebody, ncpg_promoter)
composite_c2  <- compute_composite(z_hmc_c2,  z_mc_c2,  "C2 vs C1 (early)",
                                   sd0_hmc_genes, sd0_mc_genes,
                                   ncpg_genebody, ncpg_promoter)
composite_c3  <- compute_composite(z_hmc_c3,  z_mc_c3,  "C3/EOT vs C1 (late)",
                                   sd0_hmc_genes, sd0_mc_genes,
                                   ncpg_genebody, ncpg_promoter)

openxlsx::write.xlsx(composite_all, file = file.path(path, "Composite_Z_Scores_C3C2_vs_C1.xlsx"), overwrite = TRUE)
openxlsx::write.xlsx(composite_c2,  file = file.path(path, "Composite_Z_Scores_C2_vs_C1.xlsx"),   overwrite = TRUE)
openxlsx::write.xlsx(composite_c3,  file = file.path(path, "Composite_Z_Scores_C3_vs_C1.xlsx"),   overwrite = TRUE)

# ── Full bias audit ───────────────────────────────────────────────────────────
# Tracks CpG context correlation through every transformation step
audit_cpg_bias <- function(label, 
                           hmc_gene, mc_gene,
                           z_hmc, z_mc,
                           composite_df,
                           hmc_c1_cols, mc_c1_cols,
                           ncpg_genebody, ncpg_promoter) {
  
  cat("\n════════════════════════════════════════════════════\n")
  cat(" CpG Bias Audit —", label, "\n")
  cat("════════════════════════════════════════════════════\n")
  
  # ── 1. Raw 5hmC fraction (baseline mean) vs Gene_body CpGs ─────────────
  dmr_hmc <- hmc_gene %>%
    dplyr::mutate(mean_raw = rowMeans(
      dplyr::select(., all_of(hmc_c1_cols)), na.rm = TRUE)) %>%
    dplyr::select(Name, mean_raw) %>%
    dplyr::inner_join(ncpg_genebody, by = "Name")
  
  rho_dmr_hmc <- cor.test(dmr_hmc$total_CpGs_hmc, 
                          dmr_hmc$mean_raw,
                          method = "spearman")$estimate
  
  # ── 2. Raw 5mC fraction (baseline mean) vs Promoter CpGs ───────────────
  dmr_mc <- mc_gene %>%
    dplyr::mutate(mean_raw = rowMeans(
      dplyr::select(., all_of(mc_c1_cols)), na.rm = TRUE)) %>%
    dplyr::select(Name, mean_raw) %>%
    dplyr::inner_join(ncpg_promoter, by = "Name")
  
  rho_dmr_mc <- cor.test(dmr_mc$total_CpGs_mc,
                         dmr_mc$mean_raw,
                         method = "spearman")$estimate
  
  # ── 3. Uncorrected Z_hmc mean vs Gene_body CpGs ─────────────────────────
  z_hmc_mean <- z_hmc %>%
    dplyr::mutate(Z_hmc_uncorrected = rowMeans(
      dplyr::select(., -gene), na.rm = TRUE)) %>%
    dplyr::select(gene, Z_hmc_uncorrected) %>%
    dplyr::inner_join(ncpg_genebody, by = c("gene" = "Name"))
  
  rho_z_hmc_uncorr <- cor.test(z_hmc_mean$total_CpGs_hmc,
                               z_hmc_mean$Z_hmc_uncorrected,
                               method = "spearman")$estimate
  
  # ── 4. Uncorrected Z_mc mean vs Promoter CpGs ───────────────────────────
  z_mc_mean <- z_mc %>%
    dplyr::mutate(Z_mc_uncorrected = rowMeans(
      dplyr::select(., -gene), na.rm = TRUE)) %>%
    dplyr::select(gene, Z_mc_uncorrected) %>%
    dplyr::inner_join(ncpg_promoter, by = c("gene" = "Name"))
  
  rho_z_mc_uncorr <- cor.test(z_mc_mean$total_CpGs_mc,
                              z_mc_mean$Z_mc_uncorrected,
                              method = "spearman")$estimate
  
  # ── 5. Corrected Z_hmc vs Gene_body CpGs ────────────────────────────────
  corr_hmc <- composite_df %>%
    dplyr::inner_join(ncpg_genebody, by = c("SYMBOL"="Name"))
  
  rho_z_hmc_corr <- cor.test(corr_hmc$total_CpGs_hmc,
                             corr_hmc$Z_hmc,
                             method = "spearman")$estimate
  
  # ── 6. Corrected Z_mc vs Promoter CpGs ──────────────────────────────────
  corr_mc <- composite_df %>%
    dplyr::inner_join(ncpg_promoter, by = c("SYMBOL"="Name"))
  
  rho_z_mc_corr <- cor.test(corr_mc$total_CpGs_mc,
                            corr_mc$Z_mc,
                            method = "spearman")$estimate
  
  # ── 7. composite_Z vs Gene_body CpGs ────────────────────────────────────
  rho_composite <- cor.test(corr_hmc$total_CpGs_hmc,
                            corr_hmc$composite_Z,
                            method = "spearman")$estimate
  
  # ── Print summary table ──────────────────────────────────────────────────
  cat(sprintf("%-45s rho = %6.3f\n", "Raw 5hmC fraction vs Gene_body CpGs:",    rho_dmr_hmc))
  cat(sprintf("%-45s rho = %6.3f\n", "Raw 5mC  fraction vs Promoter  CpGs:",    rho_dmr_mc))
  cat(sprintf("%-45s rho = %6.3f\n", "Uncorrected Z_hmc vs Gene_body CpGs:",    rho_z_hmc_uncorr))
  cat(sprintf("%-45s rho = %6.3f\n", "Uncorrected Z_mc  vs Promoter  CpGs:",    rho_z_mc_uncorr))
  cat(sprintf("%-45s rho = %6.3f\n", "Corrected Z_hmc   vs Gene_body CpGs:",    rho_z_hmc_corr))
  cat(sprintf("%-45s rho = %6.3f\n", "Corrected Z_mc    vs Promoter  CpGs:",    rho_z_mc_corr))
  cat(sprintf("%-45s rho = %6.3f\n", "composite_Z       vs Gene_body CpGs:",    rho_composite))
}

audit_cpg_bias("All post vs C1",
               hmc_gene, mc_gene,
               z_hmc_all, z_mc_all,
               composite_all,
               hmc_c1_cols, mc_c1_cols,
               ncpg_genebody, ncpg_promoter)

audit_cpg_bias("C2 vs C1 early",
               hmc_gene, mc_gene,
               z_hmc_c2, z_mc_c2,
               composite_c2,
               hmc_c1_cols, mc_c1_cols,
               ncpg_genebody, ncpg_promoter)

audit_cpg_bias("C3/EOT vs C1 late",
               hmc_gene, mc_gene,
               z_hmc_c3, z_mc_c3,
               composite_c3,
               hmc_c1_cols, mc_c1_cols,
               ncpg_genebody, ncpg_promoter)

# ═══════════════════════════════════════════════════════════════════════════════
# STEPS 4-6 COMMENTED OUT
# Reason: Spearman correlation between composite_Z and RNA-seq log2FC = 0.0005
# (p = 0.98) within tumor-specific 5hmC genes — no gene-level concordance.
# Hypergeometric test also non-significant (FE = 1.12, p = 0.086).
# Species mismatch (human cfDNA vs mouse RNA-seq) and cell-type heterogeneity
# in cfDNA likely explain absence of gene-level correlation.
# Pathway-level analysis (GSEA) pursued instead — see Step 7.
# ═══════════════════════════════════════════════════════════════════════════════

# # ═══════════════════════════════════════════════════════════════════════════════
# # STEP 4: INTEGRATE WITH RNA-SEQ + ASSIGN EPIGENETIC STATES
# # ═══════════════════════════════════════════════════════════════════════════════
# 
# integrate_rnaseq <- function(composite, rnaseq, comparison_name) {
#   
#   integrated <- composite %>%
#     dplyr::inner_join(rnaseq %>% dplyr::select(SYMBOL, log2FoldChange, padj), by = "SYMBOL") %>%
#     dplyr::filter(!is.na(composite_Z), !is.na(log2FoldChange)) %>%
#     mutate(Epigenetic_State = case_when(
#       composite_Z >  2 & log2FoldChange >  1  ~ "Activated",
#       composite_Z >  2 & log2FoldChange <= 1  ~ "Epigenetically Primed",
#       composite_Z < -2 & log2FoldChange < -1  ~ "Repressed",
#       composite_Z >= -2 & composite_Z <= 2 & log2FoldChange > 1 ~ "Expression Only",
#       TRUE ~ "Stable"
#     ))
#   
#   cat("\n===", comparison_name, "integrated ===\n")
#   cat("Genes after integration:", nrow(integrated), "\n")
#   
#   cor_result <- cor.test(integrated$composite_Z, integrated$log2FoldChange, method = "spearman")
#   cat("Spearman r:", round(cor_result$estimate, 3),
#       "(p =", round(cor_result$p.value, 4), ")\n")
#   print(table(integrated$Epigenetic_State))
#   
#   integrated
# }
# 
# integrated_all <- integrate_rnaseq(composite_all, deseq2_dfs[["AC_A"]], "All post vs C1")
# integrated_c2  <- integrate_rnaseq(composite_c2,  deseq2_dfs[["AC_A"]], "C2 vs C1 early")
# integrated_c3  <- integrate_rnaseq(composite_c3,  deseq2_dfs[["AC_A"]], "C3/EOT vs C1 late")
# 
# 
# # ═══════════════════════════════════════════════════════════════════════════════
# # STEP 5: FILTER TO TUMOR-SPECIFIC 5hmC GENES
# # Justification: hypergeometric test showed only 5hmC tumor genes have
# # transcriptional concordance; 5mC promoter genes excluded (all tests ns)
# # ═══════════════════════════════════════════════════════════════════════════════
# 
# tumor_hmc <- read.xlsx(file.path(output_dir, "sig_hmc.xlsx")) %>%
#   dplyr::filter(Type == "Tumor only")
# # tumor_mc excluded: 5mC promoter changes showed no transcriptional concordance
# # across all hypergeometric comparisons (strict, lenient, tumor-only, tumor+normal)
# 
# tumor_all <- integrated_all %>% dplyr::filter(SYMBOL %in% tumor_hmc$Name)
# tumor_c2  <- integrated_c2  %>% dplyr::filter(SYMBOL %in% tumor_hmc$Name)
# tumor_c3  <- integrated_c3  %>% dplyr::filter(SYMBOL %in% tumor_hmc$Name)
# 
# cat("\nTumor-specific 5hmC genes in integrated_all:", nrow(tumor_all), "\n")
# cat("Tumor-specific 5hmC genes in integrated_c2: ", nrow(tumor_c2),  "\n")
# cat("Tumor-specific 5hmC genes in integrated_c3: ", nrow(tumor_c3),  "\n")
# 
# 
# # ═══════════════════════════════════════════════════════════════════════════════
# # STEP 6: QUADRANT PLOT
# # X-axis: transcriptional change (RNA-seq log2FC, combination vs apalutamide)
# # Y-axis: composite epigenetic remodeling score (5hmC gain - 5mC loss)
# # Gene selection: tumor-specific 5hmC genes only (justified by hypergeometric)
# # ═══════════════════════════════════════════════════════════════════════════════
# 
# fc_cut <- 1
# z_cut  <- 2
# 
# make_quadrant_plot <- function(plot_data, label_data, title_suffix) {
#   ggplot(plot_data, aes(x = log2FoldChange, y = composite_Z)) +
#     
#     geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
#     geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
#     geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dotted", color = "firebrick", linewidth = 0.6) +
#     geom_hline(yintercept = c(-z_cut,  z_cut),  linetype = "dotted", color = "blue4",     linewidth = 0.6) +
#     
#     geom_point(aes(color = Epigenetic_State), alpha = 0.4, size = 1.5) +
#     
#     geom_text_repel(
#       data = label_data,
#       aes(label = SYMBOL),
#       size = 3.5, fontface = "italic", max.overlaps = 15, box.padding = 0.5
#     ) +
#     
#     scale_color_manual(values = c(
#       "Activated"             = "#e31a1c",
#       "Epigenetically Primed" = "#ff7f00",
#       "Expression Only"       = "#33a02c",
#       "Repressed"             = "#1f78b4",
#       "Stable"                = "gray90"
#     )) +
#     
#     scale_x_continuous(breaks = seq(-10, 10, by = 1)) +
#     scale_y_continuous(breaks = seq(-10, 30, by = 2)) +
#     
#     labs(
#       title    = paste("Epigenetic Remodeling vs. Transcriptional Output —", title_suffix),
#       subtitle = "Gene selection: tumor-specific 5hmC | Y-axis: composite Z (5hmC gain − 5mC loss)\nThresholds: |Log2FC| > 1, |Composite Z| > 2",
#       x        = "Transcriptional Change (RNA log2FC: Combination vs Apalutamide)",
#       y        = "Epigenetic Remodeling Score (Composite Z)",
#       color    = "Biological State"
#     ) +
#     theme_minimal() +
#     theme(
#       panel.grid.minor = element_blank(),
#       legend.position  = "right",
#       plot.title       = element_text(face = "bold")
#     )
# }
# 
# # All post-treatment timepoints
# p_all <- make_quadrant_plot(
#   plot_data  = tumor_all,
#   label_data = subset(tumor_all, composite_Z > 10 | (composite_Z > 5 & log2FoldChange > 3) | composite_Z < -5),
#   title_suffix = "All Post vs C1"
# )
# 
# # Early (C2 vs C1)
# p_c2 <- make_quadrant_plot(
#   plot_data  = tumor_c2,
#   label_data = subset(tumor_c2, composite_Z > 10 | (composite_Z > 5 & log2FoldChange > 3) | composite_Z < -5),
#   title_suffix = "C2 vs C1 (Early)"
# )
# 
# # Late (C3/EOT vs C1)
# p_c3 <- make_quadrant_plot(
#   plot_data  = tumor_c3,
#   label_data = subset(tumor_c3, composite_Z > 10 | (composite_Z > 5 & log2FoldChange > 3) | composite_Z < -5),
#   title_suffix = "C3/EOT vs C1 (Late)"
# )
# 
# print(p_all)
# print(p_c2)
# print(p_c3)
# 
# # ── Save plots ────────────────────────────────────────────────────────────────
# ggsave(file.path(output_dir, "quadrant_all.pdf"), p_all, width = 10, height = 8)
# ggsave(file.path(output_dir, "quadrant_c2.pdf"),  p_c2,  width = 10, height = 8)
# ggsave(file.path(output_dir, "quadrant_c3.pdf"),  p_c3,  width = 10, height = 8)

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 7: GSEA
# ═══════════════════════════════════════════════════════════════════════════════

# ── GSEA using composite_Z ────────────────────────────────────────────────────
input_all <- read.xlsx(file.path(path, "Composite_Z_Scores_C3C2_vs_C1.xlsx")) %>%
  dplyr::select(SYMBOL, log2FoldChange = composite_Z)

input_c2 <- read.xlsx(file.path(path, "Composite_Z_Scores_C2_vs_C1.xlsx")) %>%
  dplyr::select(SYMBOL, log2FoldChange = composite_Z)

input_c3 <- read.xlsx(file.path(path, "Composite_Z_Scores_C3_vs_C1.xlsx")) %>%
  dplyr::select(SYMBOL, log2FoldChange = composite_Z)

gsea_all <- analyze_pathways(contrast = "composite_all", input_data = input_all, gmt_dir = gmt_dir, output_dir = output_dir)
gsea_c2  <- analyze_pathways(contrast = "composite_c2",  input_data = input_c2,  gmt_dir = gmt_dir, output_dir = output_dir)
gsea_c3  <- analyze_pathways(contrast = "composite_c3",  input_data = input_c3,  gmt_dir = gmt_dir, output_dir = output_dir)

# ── GSEA using Z_hmc only ─────────────────────────────────────────────────────
input_all_hmc <- read.xlsx(file.path(path, "Composite_Z_Scores_C3C2_vs_C1.xlsx")) %>%
  dplyr::select(SYMBOL, log2FoldChange = Z_hmc)

input_c2_hmc <- read.xlsx(file.path(path, "Composite_Z_Scores_C2_vs_C1.xlsx")) %>%
  dplyr::select(SYMBOL, log2FoldChange = Z_hmc)

input_c3_hmc <- read.xlsx(file.path(path, "Composite_Z_Scores_C3_vs_C1.xlsx")) %>%
  dplyr::select(SYMBOL, log2FoldChange = Z_hmc)

gsea_all_hmc <- analyze_pathways(contrast = "composite_all_hmc", input_data = input_all_hmc, gmt_dir = gmt_dir, output_dir = output_dir)
gsea_c2_hmc  <- analyze_pathways(contrast = "composite_c2_hmc",  input_data = input_c2_hmc,  gmt_dir = gmt_dir, output_dir = output_dir)
gsea_c3_hmc  <- analyze_pathways(contrast = "composite_c3_hmc",  input_data = input_c3_hmc,  gmt_dir = gmt_dir, output_dir = output_dir)

# ── GSEA using Z_mc only (flipped: high Z_mc = repression) ───────────────────
input_all_mc <- read.xlsx(file.path(path, "Composite_Z_Scores_C3C2_vs_C1.xlsx")) %>%
  dplyr::select(SYMBOL, log2FoldChange = Z_mc) %>%
  dplyr::mutate(log2FoldChange = -1 * log2FoldChange)

input_c2_mc <- read.xlsx(file.path(path, "Composite_Z_Scores_C2_vs_C1.xlsx")) %>%
  dplyr::select(SYMBOL, log2FoldChange = Z_mc) %>%
  dplyr::mutate(log2FoldChange = -1 * log2FoldChange)

input_c3_mc <- read.xlsx(file.path(path, "Composite_Z_Scores_C3_vs_C1.xlsx")) %>%
  dplyr::select(SYMBOL, log2FoldChange = Z_mc) %>%
  dplyr::mutate(log2FoldChange = -1 * log2FoldChange)

gsea_all_mc <- analyze_pathways(contrast = "composite_all_mc", input_data = input_all_mc, gmt_dir = gmt_dir, output_dir = output_dir)
gsea_c2_mc  <- analyze_pathways(contrast = "composite_c2_mc",  input_data = input_c2_mc,  gmt_dir = gmt_dir, output_dir = output_dir)
gsea_c3_mc  <- analyze_pathways(contrast = "composite_c3_mc",  input_data = input_c3_mc,  gmt_dir = gmt_dir, output_dir = output_dir)

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 8: Bulk vs Biomodal
# ═══════════════════════════════════════════════════════════════════════════════

bulk_gsea <- read.xlsx("/hpc/home/kailasamms/scratch/RNASeq_Xenograft_Kevin_22Rv1_Xeno (A1,C1,U3 removed)/Human/09.Pathways/Apalutamide_Carotuximab-Apalutamide/Pathways.xlsx")

biomodal_all_gsea  <- gsea_all[[1]]
biomodal_c2_gsea   <- gsea_c2[[1]]
biomodal_c3_gsea   <- gsea_c3[[1]]

biomodal_all_gsea_hmc  <- gsea_all_hmc[[1]]
biomodal_c2_gsea_hmc   <- gsea_c2_hmc[[1]]
biomodal_c3_gsea_hmc   <- gsea_c3_hmc[[1]]

biomodal_all_gsea_mc  <- gsea_all_mc[[1]]
biomodal_c2_gsea_mc   <- gsea_c2_mc[[1]]
biomodal_c3_gsea_mc   <- gsea_c3_mc[[1]]

# tumor_hmc_genes <- sig_hmc %>%
#   filter(Type == "Tumor only") %>%
#   pull(Name)
# 
# tumor_mc_genes <- sig_mc %>%
#   filter(Type == "Tumor only") %>%
#   pull(Name)

# Define the comparison function
get_pathway_overlap <- function(bulk_df, biomodal_df, comparison_name){
                                #tumor_hmc_genes, tumor_mc_genes) {
  
  # Ensure we are only looking at columns we need to avoid clutter
  bulk_clean <- bulk_df %>% 
    select(Collection, Description, NES_Bulk = NES, pval_Bulk = pval, padj_Bulk = padj, Edge_Bulk = leadingEdge)
  
  biomodal_clean <- biomodal_df %>% 
    select(Collection, Description, NES_Biomodal = NES, pval_Biomodal = pval, padj_Biomodal = padj, Edge_Biomodal = leadingEdge)
  
  # Merge on Collection and Description
  comparison <- inner_join(bulk_clean, biomodal_clean, by = c("Collection", "Description")) %>%
    mutate(
      # ── parse once safely ─────────────────────────────
      bulk_genes = strsplit(ifelse(is.na(Edge_Bulk), "", as.character(Edge_Bulk)), "/"),
      bio_genes  = strsplit(ifelse(is.na(Edge_Biomodal), "", as.character(Edge_Biomodal)), "/"),
      
      # ── shared genes ──────────────────────────────────
      Common_Genes = mapply(function(x, y) {
        common <- intersect(x, y)
        if(length(common) > 0) paste(common, collapse = "/") else NA_character_
      }, bulk_genes, bio_genes),
      
      Common_Gene_Count = mapply(function(x, y) {
        length(intersect(x, y))
      }, bulk_genes, bio_genes),
      
      # # ── tumor-specific enrichment (NEW) ───────────────
      # tumor_hmc_n = mapply(function(x, y) {
      #   length(intersect(intersect(x, y), tumor_hmc_genes))
      # }, bulk_genes, bio_genes),
      # 
      # tumor_mc_n = mapply(function(x, y) {
      #   length(intersect(intersect(x, y), tumor_mc_genes))
      # }, bulk_genes, bio_genes),
      # 
      # tumor_total_n = tumor_hmc_n + tumor_mc_n,
      
      # ── directionality ────────────────────────────────
      Status = ifelse(sign(NES_Bulk) == sign(NES_Biomodal),
                      "Concordant", "Discordant"),
      
      Direction = case_when(
        NES_Bulk > 0 & NES_Biomodal > 0 ~ "Both Up (Opening)",
        NES_Bulk < 0 & NES_Biomodal < 0 ~ "Both Down (Closing)",
        NES_Bulk > 0 & NES_Biomodal < 0 ~ "Bulk Up / Biomodal Down",
        NES_Bulk < 0 & NES_Biomodal > 0 ~ "Bulk Down / Biomodal Up"
      ),
      
      Comparison_Group = comparison_name
    ) %>%
    # Remove unwanted columns
    select(-Edge_Bulk, -Edge_Biomodal, -bulk_genes, -bio_genes)

  return(comparison)
}

results_all <- get_pathway_overlap(bulk_gsea, biomodal_all_gsea, "All")#, tumor_hmc_genes, tumor_mc_genes)
results_c2  <- get_pathway_overlap(bulk_gsea, biomodal_c2_gsea,  "C2")#, tumor_hmc_genes, tumor_mc_genes)
results_c3  <- get_pathway_overlap(bulk_gsea, biomodal_c3_gsea,  "C3")#, tumor_hmc_genes, tumor_mc_genes)

all_comparisons <- bind_rows(results_all, results_c2, results_c3)
openxlsx::write.xlsx(all_comparisons, file = file.path(path, "Pathway_Comparisons.xlsx"), overwrite = TRUE)

results_all_hmc <- get_pathway_overlap(bulk_gsea, biomodal_all_gsea_hmc, "All")#, tumor_hmc_genes, tumor_mc_genes)
results_c2_hmc  <- get_pathway_overlap(bulk_gsea, biomodal_c2_gsea_hmc,  "C2")#, tumor_hmc_genes, tumor_mc_genes)
results_c3_hmc  <- get_pathway_overlap(bulk_gsea, biomodal_c3_gsea_hmc,  "C3")#, tumor_hmc_genes, tumor_mc_genes)

all_comparisons_hmc <- bind_rows(results_all_hmc, results_c2_hmc, results_c3_hmc)
openxlsx::write.xlsx(all_comparisons_hmc, file = file.path(path, "Pathway_Comparisons_hmc.xlsx"), overwrite = TRUE)

results_all_mc <- get_pathway_overlap(bulk_gsea, biomodal_all_gsea_mc, "All")#, tumor_hmc_genes, tumor_mc_genes)
results_c2_mc  <- get_pathway_overlap(bulk_gsea, biomodal_c2_gsea_mc,  "C2")#, tumor_hmc_genes, tumor_mc_genes)
results_c3_mc  <- get_pathway_overlap(bulk_gsea, biomodal_c3_gsea_mc,  "C3")#, tumor_hmc_genes, tumor_mc_genes)

all_comparisons_mc <- bind_rows(results_all_mc, results_c2_mc, results_c3_mc)
openxlsx::write.xlsx(all_comparisons_mc, file = file.path(path, "Pathway_Comparisons_mc.xlsx"), overwrite = TRUE)

# Check if pathways Neil is interested are relevant
keywords <- c("androgen", "protein.*fold", "fold.*protein", 
              "chaperone", "apoptosis", "caspase", "heat.*shock")

check_pathways <- function(results_df, label) {
  results_df %>%
    dplyr::filter(grepl(paste(keywords, collapse = "|"), 
                        Pathway, ignore.case = TRUE)) %>%
    dplyr::select(Pathway, Timepoint, N_genes, Median_Z, Background_Z, 
                  Shift, Direction, P_Value, FDR) %>%
    dplyr::arrange(P_Value) %>%
    dplyr::mutate(Source = label)
}

keyword_hits <- bind_rows(
  check_pathways(wilcoxon_results %>% filter(Timepoint == "All"), "All"),
  check_pathways(wilcoxon_results %>% filter(Timepoint == "C2"),  "C2"),
  check_pathways(wilcoxon_results %>% filter(Timepoint == "C3"),  "C3")
)

# Print to console
keyword_hits %>% 
  select(Pathway, Timepoint, Shift, Direction, P_Value, FDR) %>%
  print(n = 50)

# Save separately
openxlsx::write.xlsx(keyword_hits,
                     file = file.path(path, "Keyword_Pathway_Hits.xlsx"),
                     overwrite = TRUE)

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 9: Wilcoxon shift test
# ═══════════════════════════════════════════════════════════════════════════════

composite_all <- read.xlsx(file.path(path, "Composite_Z_Scores_C3C2_vs_C1.xlsx"))
composite_c2  <- read.xlsx(file.path(path, "Composite_Z_Scores_C2_vs_C1.xlsx"))
composite_c3  <- read.xlsx(file.path(path, "Composite_Z_Scores_C3_vs_C1.xlsx"))

# ── Wilcoxon shift test for all GMT pathways ──────────────────────────────────
run_wilcoxon_gsea <- function(composite_df, gmt_files, score_col, contrast_name) {
  
  cat("\n=== Wilcoxon GSEA —", contrast_name, "===\n")
  
  # Read all GMT files
  gene_sets <- list()
  for (gmt_file in gmt_files) {
    lines <- readLines(gmt_file)
    for (line in lines) {
      parts <- strsplit(line, "\t")[[1]]
      if (length(parts) < 3) next
      gs_name <- parts[1]
      genes   <- parts[3:length(parts)]
      genes   <- genes[genes != "" & !is.na(genes)]
      gene_sets[[gs_name]] <- genes
    }
  }
  
  cat("Total gene sets loaded:", length(gene_sets), "\n")
  
  # Get score vector
  scores     <- composite_df[[score_col]]
  symbols    <- composite_df$SYMBOL
  background <- scores[!is.na(scores)]
  
  # Test each pathway
  results <- purrr::map_df(names(gene_sets), function(gs_name) {
    
    p_genes    <- gene_sets[[gs_name]]
    in_pathway <- scores[symbols %in% p_genes & !is.na(scores)]
    bg         <- scores[!symbols %in% p_genes & !is.na(scores)]
    
    if (length(in_pathway) < 5) {
      return(tibble(
        Pathway    = gs_name,
        N_genes    = length(in_pathway),
        Median_Z   = NA_real_,
        Background = NA_real_,
        Shift      = NA_real_,
        P_Value    = NA_real_,
        Direction  = NA_character_
      ))
    }
    
    wt <- wilcox.test(in_pathway, bg, alternative = "two.sided")
    
    tibble(
      Pathway    = gs_name,
      N_genes    = length(in_pathway),
      Median_Z   = round(median(in_pathway), 3),
      Background = round(median(bg), 3),
      Shift      = round(median(in_pathway) - median(bg), 3),
      P_Value    = signif(wt$p.value, 3),
      Direction  = ifelse(median(in_pathway) > median(bg), "UP", "DOWN")
    )
  })
  
  # FDR correction
  results <- results %>%
    dplyr::filter(!is.na(P_Value)) %>%
    dplyr::mutate(
      FDR = p.adjust(P_Value, method = "BH"),
      Sig = case_when(
        FDR < 0.001 ~ "***",
        FDR < 0.01  ~ "**",
        FDR < 0.05  ~ "*",
        P_Value < 0.05 ~ "nominal",
        TRUE ~ "ns"
      )
    ) %>%
    dplyr::arrange(P_Value)
  
  cat("Significant (FDR<0.05):", sum(results$FDR < 0.05, na.rm = TRUE), "\n")
  cat("Nominal (p<0.05):      ", sum(results$P_Value < 0.05, na.rm = TRUE), "\n")
  
  results
}

# ── Run for composite_Z ───────────────────────────────────────────────────────
gmt_files <- list.files(gmt_dir, full.names = TRUE, pattern = "\\.gmt$")

wilcox_all <- run_wilcoxon_gsea(composite_all, gmt_files, "composite_Z", "All post vs C1")
wilcox_c2  <- run_wilcoxon_gsea(composite_c2,  gmt_files, "composite_Z", "C2 vs C1 early")
wilcox_c3  <- run_wilcoxon_gsea(composite_c3,  gmt_files, "composite_Z", "C3/EOT vs C1 late")

# ── Run for Z_hmc only ────────────────────────────────────────────────────────
wilcox_all_hmc <- run_wilcoxon_gsea(composite_all, gmt_files, "Z_hmc", "All post vs C1 — 5hmC")
wilcox_c2_hmc  <- run_wilcoxon_gsea(composite_c2,  gmt_files, "Z_hmc", "C2 vs C1 — 5hmC")
wilcox_c3_hmc  <- run_wilcoxon_gsea(composite_c3,  gmt_files, "Z_hmc", "C3/EOT vs C1 — 5hmC")

# ── Run for Z_mc only (flipped) ───────────────────────────────────────────────
composite_all_mc_flip <- composite_all %>% dplyr::mutate(Z_mc_flip = -1 * Z_mc)
composite_c2_mc_flip  <- composite_c2  %>% dplyr::mutate(Z_mc_flip = -1 * Z_mc)
composite_c3_mc_flip  <- composite_c3  %>% dplyr::mutate(Z_mc_flip = -1 * Z_mc)

wilcox_all_mc <- run_wilcoxon_gsea(composite_all_mc_flip, gmt_files, "Z_mc_flip", "All post vs C1 — 5mC")
wilcox_c2_mc  <- run_wilcoxon_gsea(composite_c2_mc_flip,  gmt_files, "Z_mc_flip", "C2 vs C1 — 5mC")
wilcox_c3_mc  <- run_wilcoxon_gsea(composite_c3_mc_flip,  gmt_files, "Z_mc_flip", "C3/EOT vs C1 — 5mC")

# ── Save all results ──────────────────────────────────────────────────────────
openxlsx::write.xlsx(
  list("All_composite"  = wilcox_all,
       "C2_composite"   = wilcox_c2,
       "C3_composite"   = wilcox_c3,
       "All_5hmC"       = wilcox_all_hmc,
       "C2_5hmC"        = wilcox_c2_hmc,
       "C3_5hmC"        = wilcox_c3_hmc,
       "All_5mC"        = wilcox_all_mc,
       "C2_5mC"         = wilcox_c2_mc,
       "C3_5mC"         = wilcox_c3_mc),
  file      = file.path(output_dir, "Wilcoxon_Pathway_Results.xlsx"),
  overwrite = TRUE
)

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 10: Targeted Wilcoxon shift test
# For each bulk RNA-seq pathway: test if leading edge genes show
# composite_Z shift in cfDNA, across all three timepoint comparisons
# ═══════════════════════════════════════════════════════════════════════════════

bulk_gsea <- read.xlsx("/hpc/home/kailasamms/scratch/RNASeq_Xenograft_Kevin_22Rv1_Xeno (A1,C1,U3 removed)/Human/09.Pathways/Apalutamide_Carotuximab-Apalutamide/Pathways.xlsx")

# ── Named list of composite scores by timepoint ───────────────────────────────
composite_list <- list(
  All = composite_all,
  C2  = composite_c2,
  C3  = composite_c3
)

# ── Run Wilcoxon test for every pathway × timepoint combination ───────────────
wilcoxon_results <- bulk_gsea %>%
  dplyr::mutate(
    pathway_id  = paste0(Collection, "_", Description),
    leading_genes = strsplit(ifelse(is.na(leadingEdge), "", as.character(leadingEdge)), "/")
  ) %>%
  dplyr::select(pathway_id, leading_genes) %>%
  tidyr::crossing(timepoint = names(composite_list)) %>%
  dplyr::mutate(
    result = purrr::pmap(list(pathway_id, leading_genes, timepoint), function(pathway, genes, tp) {
      
      composite_df <- composite_list[[tp]]
      in_pathway   <- composite_df$composite_Z[composite_df$SYMBOL %in% genes]
      background   <- composite_df$composite_Z[!composite_df$SYMBOL %in% genes]
      
      in_pathway <- in_pathway[!is.na(in_pathway)]
      background <- background[!is.na(background)]
      
      if (length(in_pathway) < 3) return(NULL)
      
      wt <- wilcox.test(in_pathway, background, alternative = "two.sided")
      
      tibble(
        Pathway        = pathway,
        Timepoint      = tp,
        N_genes        = length(in_pathway),
        Median_Z       = round(median(in_pathway),  3),
        Background_Z   = round(median(background),  3),
        Shift          = round(median(in_pathway) - median(background), 3),
        Direction      = ifelse(median(in_pathway) > median(background),
                                "Epigenetically Opening",
                                "Epigenetically Closing"),
        P_Value        = signif(wt$p.value, 3)
      )
    })
  ) %>%
  dplyr::select(result) %>%
  tidyr::unnest(result) %>%
  dplyr::group_by(Timepoint) %>%
  dplyr::mutate(FDR = p.adjust(P_Value, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Timepoint, P_Value)

openxlsx::write.xlsx(wilcoxon_results, 
                     file = file.path(path, "Targeted_Pathway_Comparisons.xlsx"), 
                     overwrite = TRUE)

# Check results look sensible
wilcoxon_results %>%
  group_by(Timepoint) %>%
  summarise(
    n_pathways    = n(),
    n_sig_p05     = sum(P_Value <= 0.05),
    n_sig_fdr05   = sum(FDR <= 0.05),
    n_opening     = sum(P_Value <= 0.05 & Direction == "Epigenetically Opening"),
    n_closing     = sum(P_Value <= 0.05 & Direction == "Epigenetically Closing")
  )

# Pathways significant in Wilcoxon (p<0.05)
wilcoxon_sig <- wilcoxon_results %>%
  filter(P_Value <= 0.05) %>%
  mutate(Pathway_clean = gsub("^[^_]+_", "", Pathway))

# Check if those same pathways appear concordantly in Step 8
df <- all_comparisons %>%
  filter(Description %in% wilcoxon_sig$Pathway_clean,
         Status == "Concordant") %>%
  select(Description, Comparison_Group, NES_Bulk, NES_Biomodal, Status)
openxlsx::write.xlsx(df, 
                     file = file.path(path, "Wilcoxon_Sig_GSEA_Concordant.xlsx"), 
                     overwrite = TRUE)

df <- all_comparisons %>%
  filter(Description %in% wilcoxon_sig$Pathway_clean,
         Status == "Discordant") %>%
  select(Description, Comparison_Group, NES_Bulk, NES_Biomodal, Status)
openxlsx::write.xlsx(df, 
                     file = file.path(path, "Wilcoxon_Sig_GSEA_Discordant.xlsx"), 
                     overwrite = TRUE)

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 11: Overlap with Felix Fang ATAC Seq + ORA
# ═══════════════════════════════════════════════════════════════════════════════

# bulk_gsea <- read.xlsx("/hpc/home/kailasamms/scratch/RNASeq_Xenograft_Kevin_22Rv1_Xeno (A1,C1,U3 removed)/Human/09.Pathways/Apalutamide_Carotuximab-Apalutamide/Pathways.xlsx")
# gmt_files <- list.files(gmt_dir, full.names = TRUE, pattern = "\\.gmt$")
# 
# protein_folding_genes <- bulk_gsea %>%
#   dplyr::filter(grepl("fold", Description, ignore.case = TRUE),
#                 grepl("protein", Description, ignore.case = TRUE)) %>%
#   dplyr::mutate(
#     pathway_id  = paste0(Collection, "_", Description),
#     leading_genes = strsplit(ifelse(is.na(leadingEdge), "", as.character(leadingEdge)), "/")
#   ) %>%
#   #dplyr::pull(pathway_id) %>%
#   dplyr::pull(leading_genes) %>%
#   unlist() %>%
#   unique()
# 
# protein_folding_genes_hmc <- dmr_hmc %>% 
#   dplyr::filter(Name %in% protein_folding_genes, dmr_qvalue <= 0.05)
# 
# openxlsx::write.xlsx(protein_folding_genes_hmc, 
#                      file = file.path(path, "protein_folding_genes_hmc.xlsx"), 
#                      overwrite = TRUE)


#sig_hmc <- read.xlsx(file.path(path, "sig_hmc.xlsx"))
dmr_hmc <- readr::read_tsv(
  file.path(path, "biomodal_output",
            paste0("DMR_", DMR_TIMESTAMP, "_DMR_hmc_Pre__Post_", DMR_TIMESTAMP, ".tsv")),
  show_col_types = FALSE)

felix_atac <- read.xlsx(file.path(path, "DOI_10.1158_0008-5472.CAN-24-0890.xlsx"), 
                        sheet = "supplementary_table_s3", startRow = 2) %>%
  dplyr::mutate(Annotation = dplyr::case_when(feature == "Promoter" ~ "Promoter",
                                              feature %in% c("Exon", "Intron", "3' UTR", "5' UTR") ~ "Gene_body",
                                              TRUE ~ "Distal_Intergenic" # Captures Intergenic and Downstream
  ))

# Robust intersection helper
get_int <- function(df_atac, df_dmr, feature, direction) {
  atac_genes <- df_atac %>% filter(Annotation == feature) %>% pull(gene)
  
  if(direction == "up") {
    dmr_genes <- df_dmr %>% filter(Annotation == feature, mod_fold_change > 1.5, dmr_pvalue < 0.05) %>% pull(Name)
  } else {
    dmr_genes <- df_dmr %>% filter(Annotation == feature, mod_fold_change < 0.67, dmr_pvalue < 0.05) %>% pull(Name)
  }
  
  return(intersect(atac_genes, dmr_genes))
}

# Then just call it:
up_p    <- get_int(felix_atac, dmr_hmc, "Promoter", "up")
down_p  <- get_int(felix_atac, dmr_hmc, "Promoter", "down")
up_gb   <- get_int(felix_atac, dmr_hmc, "Gene_body", "up")
down_gb <- get_int(felix_atac, dmr_hmc, "Gene_body", "down")

up <- unique(c(up_p, up_gb))
down <- unique(c(down_p, down_gb))
gene_lists <- list(up = up, down = down, all = NULL)
gmt_files <- list.files(gmt_dir, full.names = TRUE, pattern = "\\.gmt$")
universe <- unique(dmr_hmc$Name)
ora_results <- run_ora(gene_lists = gene_lists,
                       universe   = universe,
                       gmt_files  = gmt_files)
ora_results %>% filter(grepl("fold", Description, ignore.case = TRUE),
                       grepl("protein", Description, ignore.case = TRUE)) %>% select(Description, pval, Direction) %>% arrange(pval)
ora_results %>% filter(grepl("androgen", Description, ignore.case = TRUE)) %>% select(Description, pval, Direction) %>% arrange(pval)
ora_results %>% filter(grepl("apoptosis", Description, ignore.case = TRUE)) %>% select(Description, pval, Direction) %>% arrange(pval) %>% head()

up <- c(up_p)
down <- c(down_p)
gene_lists <- list(up = up, down = down, all = NULL)
gmt_files <- list.files(gmt_dir, full.names = TRUE, pattern = "\\.gmt$")
universe <- unique(dmr_hmc$Name)
ora_results <- run_ora(gene_lists = gene_lists,
                       universe   = universe,
                       gmt_files  = gmt_files)
ora_results %>% filter(grepl("fold", Description, ignore.case = TRUE),
                       grepl("protein", Description, ignore.case = TRUE)) %>% select(Description, pval, Direction) %>% arrange(pval)
ora_results %>% filter(grepl("androgen", Description, ignore.case = TRUE)) %>% select(Description, pval, Direction) %>% arrange(pval)
ora_results %>% filter(grepl("apoptosis", Description, ignore.case = TRUE)) %>% select(Description, pval, Direction) %>% arrange(pval) %>% head()

up <- c(up_gb)
down <- c(down_gb)
gene_lists <- list(up = up, down = down, all = NULL)
gmt_files <- list.files(gmt_dir, full.names = TRUE, pattern = "\\.gmt$")
universe <- unique(dmr_hmc$Name)
ora_results <- run_ora(gene_lists = gene_lists,
                       universe   = universe,
                       gmt_files  = gmt_files)
ora_results %>% filter(grepl("fold", Description, ignore.case = TRUE),
                       grepl("protein", Description, ignore.case = TRUE)) %>% select(Description, pval, Direction) %>% arrange(pval)
ora_results %>% filter(grepl("androgen", Description, ignore.case = TRUE)) %>% select(Description, pval, Direction) %>% arrange(pval)
ora_results %>% filter(grepl("apoptosis", Description, ignore.case = TRUE)) %>% select(Description, pval, Direction) %>% arrange(pval) %>% head()