#!/usr/bin/env Rscript
# ===============================================================================
# 📄 Script   : 08a_biomodal_dmr_screen.R
# 📌 Purpose  : Phase 1 QC — screen all 18 DMR analyses (36 TSVs)
#               to identify which runs produce signal before deeper analysis.
#
# 📊 Per TSV summary:
#     total        — DMRs with annotated Name (non-NA)
#     q_0.05       — hits at dmr_qvalue <= 0.05
#     q_0.10       — hits at dmr_qvalue <= 0.10
#     min_qval     — minimum q-value (signal strength indicator)
#     median_qval  — median q-value
#
# 📁 Outputs:
#     dmr_phase1_summary.csv   — full 36-row table
#     dmr_phase1_summary.tsv   — same, tab-delimited
#     dmr_phase1_heatmap.pdf   — heatmap of q_0.05 hits across all runs
#
# ⚙️  Edit the paths below before running.
# ===============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(stringr)
})

# ===============================================================================
# ⚙️ USER SETTINGS
# ===============================================================================

INPUT_DIR  <- "/home/kailasamms/analysis/Biomodal/Kevin_Carotuximab/04.dmr"
OUTPUT_DIR <- "/home/kailasamms/analysis/Biomodal/Kevin_Carotuximab/05.dmr_qc"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
}

# ===============================================================================
# 📋 Build manifest of all expected TSVs
# ===============================================================================

groups   <- c("long", "short", "full")
variants <- c(
  "unpaired_no_od_no_cov",
  "unpaired_od_no_cov",
  "unpaired_no_od_cov",
  "unpaired_od_cov",
  "paired_no_od_cov",
  "paired_od_cov"
)
marks <- c("hmc", "mc")

manifest <- tidyr::expand_grid(group = groups, variant = variants, mark = marks) %>%
  dplyr::mutate(
    label     = paste0(group, "_", variant),
    dmr_dir   = file.path(INPUT_DIR, paste0("dmr_", label)),
    design    = dplyr::if_else(str_starts(variant, "paired"), "paired", "unpaired"),
    od        = dplyr::if_else(str_detect(variant, "_no_od_"), "no", "yes"),
    covariate = dplyr::if_else(str_detect(variant, "_no_cov$"), "no", "yes")
  )

cat(sprintf("📋 Manifest built: %d expected TSV files\n", nrow(manifest)))

# ===============================================================================
# 🔍 Helper: find the correct TSV for a given dir + mark
# ===============================================================================

find_tsv <- function(dmr_dir, label, mark) {

  if (!dir.exists(dmr_dir)) return(NA_character_)

  tsvs <- list.files(dmr_dir, pattern = "\\.tsv$", full.names = TRUE)

  if (length(tsvs) == 0) return(NA_character_)

  if (mark == "hmc") {
    matched <- tsvs[str_detect(basename(tsvs), "hmc")]
  } else {
    matched <- tsvs[str_detect(basename(tsvs), "mc") & !str_detect(basename(tsvs), "hmc")]
  }

  if (length(matched) == 0) return(NA_character_)
  if (length(matched) > 1) {
    warning(sprintf("Multiple TSVs matched for %s / %s — using first:\n  %s",
                    label, mark, paste(matched, collapse = "\n  ")))
  }
  matched[1]
}

# ===============================================================================
# 📊 QC function: summarise one TSV
# ===============================================================================

summarise_dmr <- function(tsv_path, label, mark) {

  if (is.na(tsv_path) || !file.exists(tsv_path)) {
    return(tibble(
      tsv_found   = FALSE,
      total       = NA_integer_,
      q_0.05      = NA_integer_,
      q_0.10      = NA_integer_,
      min_qval    = NA_real_,
      median_qval = NA_real_
    ))
  }

  raw <- tryCatch(
    readr::read_tsv(tsv_path, show_col_types = FALSE),
    error = function(e) {
      warning(sprintf("Failed to read %s: %s", tsv_path, e$message))
      return(NULL)
    }
  )

  if (is.null(raw)) {
    return(tibble(
      tsv_found   = TRUE,
      total       = NA_integer_,
      q_0.05      = NA_integer_,
      q_0.10      = NA_integer_,
      min_qval    = NA_real_,
      median_qval = NA_real_
    ))
  }

  if (!"dmr_qvalue" %in% colnames(raw)) {
    warning(sprintf("Column 'dmr_qvalue' not found in %s\n  Columns: %s",
                    tsv_path, paste(colnames(raw), collapse = ", ")))
    return(tibble(
      tsv_found   = TRUE,
      total       = NA_integer_,
      q_0.05      = NA_integer_,
      q_0.10      = NA_integer_,
      min_qval    = NA_real_,
      median_qval = NA_real_
    ))
  }

  if (!"Name" %in% colnames(raw)) {
    warning(sprintf("Column 'Name' not found in %s\n  Columns: %s",
                    tsv_path, paste(colnames(raw), collapse = ", ")))
  }

  filtered <- if ("Name" %in% colnames(raw)) {
    raw %>% dplyr::filter(!is.na(Name), Name != "")
  } else {
    raw
  }

  filtered %>%
    dplyr::summarise(
      tsv_found   = TRUE,
      total       = n(),
      q_0.05      = sum(dmr_qvalue <= 0.05, na.rm = TRUE),
      q_0.10      = sum(dmr_qvalue <= 0.10, na.rm = TRUE),
      min_qval    = min(dmr_qvalue,    na.rm = TRUE),
      median_qval = median(dmr_qvalue, na.rm = TRUE)
    )
}

# ===============================================================================
# 🔁 Loop over manifest and collect results
# ===============================================================================

cat("\n🔁 Processing TSVs...\n")

results <- manifest %>%
  dplyr::mutate(
    tsv_path = purrr::map2_chr(dmr_dir, mark, ~ find_tsv(.x, label, .y))
  )

qc_list <- vector("list", nrow(results))
for (i in seq_len(nrow(results))) {
  row <- results[i, ]
  cat(sprintf("  [%02d/36] %s / %s ... ", i, row$label, row$mark))
  qc <- summarise_dmr(row$tsv_path, row$label, row$mark)
  qc_list[[i]] <- qc
  if (!qc$tsv_found) {
    cat("⚠️  TSV not found\n")
  } else if (is.na(qc$total)) {
    cat("❌ Read error\n")
  } else {
    cat(sprintf("total=%d  q≤0.05=%d  q≤0.10=%d  min_q=%.4f\n",
                qc$total, qc$q_0.05, qc$q_0.10, qc$min_qval))
  }
}

summary_tbl <- results %>%
  dplyr::select(group, variant, design, od, covariate, mark, label, tsv_path) %>%
  dplyr::bind_cols(dplyr::bind_rows(qc_list)) %>%
  dplyr::mutate(
    od_zero_flag = od == "yes" & !is.na(q_0.05) & q_0.05 == 0
  ) %>%
  dplyr::arrange(group, design, od, covariate, mark)

# ===============================================================================
# 💾 Save outputs
# ===============================================================================

csv_out <- file.path(OUTPUT_DIR, "dmr_phase1_summary.csv")
tsv_out <- file.path(OUTPUT_DIR, "dmr_phase1_summary.tsv")

readr::write_csv(summary_tbl, csv_out)
readr::write_tsv(summary_tbl, tsv_out)

cat(sprintf("\n💾 Summary saved:\n  %s\n  %s\n", csv_out, tsv_out))

# ===============================================================================
# 🖨️ Print clean console table
# ===============================================================================

cat("\n========================================================\n")
cat("📊 PHASE 1 QC SUMMARY — all 36 TSVs\n")
cat("========================================================\n")

print_tbl <- summary_tbl %>%
  dplyr::select(group, design, od, covariate, mark,
                total, q_0.05, q_0.10, min_qval, median_qval,
                tsv_found, od_zero_flag) %>%
  dplyr::mutate(
    min_qval    = round(min_qval,    4),
    median_qval = round(median_qval, 4),
    flag        = dplyr::case_when(
      !tsv_found                 ~ "⚠️ MISSING",
      is.na(total)               ~ "❌ READ ERROR",
      od_zero_flag               ~ "🚨 OD=0 HITS",
      q_0.05 == 0 & q_0.10 == 0 ~ "⬜ NO HITS",
      q_0.05 == 0 & q_0.10  > 0 ~ "🟡 q≤0.10 only",
      q_0.05  > 0               ~ "✅ HAS HITS",
      TRUE                       ~ ""
    )
  ) %>%
  dplyr::select(-od_zero_flag, -tsv_found)

print(as.data.frame(print_tbl), row.names = FALSE)

# ===============================================================================
# 📊 Heatmap: q_0.05 hits across all 36 runs
# ===============================================================================

cat("\n📊 Generating heatmap...\n")

heatmap_tbl <- summary_tbl %>%
  dplyr::mutate(
    q_0.05_plot = dplyr::if_else(is.na(q_0.05), -1L, q_0.05)
  )

p <- ggplot(heatmap_tbl,
            aes(x = interaction(design, od, covariate, sep = "\n"),
                y = interaction(group, mark, sep = " / "),
                fill = q_0.05_plot)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = dplyr::if_else(q_0.05_plot < 0, "NA", as.character(q_0.05_plot))),
            size = 3, color = "black") +
  scale_fill_gradient2(
    low      = "#d73027",
    mid      = "#ffffbf",
    high     = "#1a9850",
    midpoint = 5,
    na.value = "grey80",
    name     = "Hits\n(q≤0.05)"
  ) +
  labs(
    title    = "DMR Phase 1 QC — Hits at q≤0.05",
    subtitle = "Rows = group / mark | Columns = design × OD × covariate",
    x        = "Analysis variant",
    y        = "Group / Mark"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y     = element_text(size = 9),
    plot.title      = element_text(face = "bold"),
    legend.position = "right"
  )

pdf_out <- file.path(OUTPUT_DIR, "dmr_phase1_heatmap.pdf")
ggsave(pdf_out, plot = p, width = 14, height = 6)
cat(sprintf("💾 Heatmap saved: %s\n", pdf_out))

# ===============================================================================
# 📋 Quick summary: which runs have hits?
# ===============================================================================

cat("\n========================================================\n")
cat("✅ RUNS WITH HITS AT q≤0.05\n")
cat("========================================================\n")
hits <- summary_tbl %>%
  dplyr::filter(!is.na(q_0.05), q_0.05 > 0) %>%
  dplyr::select(group, design, od, covariate, mark, total, q_0.05, q_0.10, min_qval)
if (nrow(hits) == 0) {
  cat("  ⚠️  No runs produced hits at q≤0.05\n")
} else {
  print(as.data.frame(hits), row.names = FALSE)
}

cat("\n========================================================\n")
cat("🚨 OD RUNS WITH ZERO HITS (collapsed models)\n")
cat("========================================================\n")
od_zeros <- summary_tbl %>%
  dplyr::filter(od == "yes", !is.na(q_0.05), q_0.05 == 0) %>%
  dplyr::select(group, design, od, covariate, mark, total, q_0.05, min_qval)
if (nrow(od_zeros) == 0) {
  cat("  ✅ No OD runs collapsed to zero hits\n")
} else {
  print(as.data.frame(od_zeros), row.names = FALSE)
}

cat("\n========================================================\n")
cat(sprintf("End: %s\n", Sys.time()))
cat("========================================================\n")
