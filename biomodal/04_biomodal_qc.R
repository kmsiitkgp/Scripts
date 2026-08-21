#!/usr/bin/env Rscript

#===============================================================================
# 📦 Libraries
#===============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
  library(dplyr)
})

#===============================================================================
# ⚙️ Parameters — adjust these as needed
#===============================================================================

BASE_DIR   <- "/home/kailasamms/analysis/Biomodal/Kevin_Carotuximab"
OUTPUT_DIR <- file.path(BASE_DIR, "02.QC")

# QC Thresholds
THRESH_DUP_RATE        <- 0.95   # FAIL if duplication rate ABOVE this
THRESH_MEAN_COV        <- 1.0    # FAIL if mean coverage (excl Ns) BELOW this
THRESH_MODC_SENS       <- 0.95   # FAIL if modC sensitivity (lambda) BELOW this
THRESH_MODC_SPEC       <- 0.95   # FAIL if modC specificity (pUC19) BELOW this
THRESH_READ_RESOLUTION <- 0.90   # FAIL if % input reads that resolve BELOW this

#===============================================================================
# 🔧 Helper: read one Summary.xlsx → two columns (Metric, SampleID)
#===============================================================================

read_summary_xlsx <- function(filepath) {

  # Read entire sheet with no column names
  raw <- read_excel(filepath, col_names = FALSE, .name_repair = "minimal")

  # Find row where first cell == "Sample ID"
  start_row <- which(trimws(as.character(raw[[1]])) == "Sample ID")
  if (length(start_row) == 0) {
    warning(paste("Could not find 'Sample ID' row in:", filepath))
    return(NULL)
  }

  # Slice from Sample ID row onwards, first two columns only
  dat <- raw[start_row:nrow(raw), 1:2]
  colnames(dat) <- c("Metric", "Value")

  # Get sample name from first row (the Sample ID row)
  sample_id <- trimws(as.character(dat$Value[1]))

  # Rename Value column to the sample ID
  colnames(dat)[2] <- sample_id

  return(dat)
}

#===============================================================================
# 🔎 Find all Summary.xlsx files
#===============================================================================

cat("========================================================\n")
cat("🔎 Searching for Summary.xlsx files...\n")

xlsx_files <- list.files(
  path       = BASE_DIR,
  pattern    = "*Summary.xlsx",
  recursive  = TRUE,
  full.names = TRUE
)

cat(sprintf("Found %d Summary.xlsx files\n", length(xlsx_files)))
if (length(xlsx_files) == 0) stop("❌ No Summary.xlsx files found. Check BASE_DIR.")

#===============================================================================
# 📊 Parse all samples and merge by column binding
#===============================================================================

cat("⚙️  Parsing samples...\n")

all_dfs <- lapply(xlsx_files, function(f) {
  cat(sprintf("  → %s\n", basename(f)))
  tryCatch(read_summary_xlsx(f), error = function(e) {
    warning(paste("Failed to parse:", f, "-", e$message))
    return(NULL)
  })
})

# Remove NULLs
all_dfs <- Filter(Negate(is.null), all_dfs)

# Use first sample as base (has Metric column)
merged_df <- all_dfs[[1]]

# Bind remaining samples by joining on Metric column
for (i in seq(2, length(all_dfs))) {
  merged_df <- full_join(merged_df, all_dfs[[i]], by = "Metric")
}

# Sort sample columns alphabetically, keep Metric first
sample_cols <- sort(setdiff(colnames(merged_df), "Metric"))
merged_df   <- merged_df[, c("Metric", sample_cols)]

cat(sprintf("✅ Merged %d samples, %d rows\n", length(sample_cols), nrow(merged_df)))

#===============================================================================
# 🚦 PASS/FAIL per sample
#===============================================================================

cat("🚦 Applying QC thresholds...\n")
cat(sprintf("   Duplication rate     > %.0f%%  → FAIL\n", THRESH_DUP_RATE * 100))
cat(sprintf("   Mean coverage        < %.1fx   → FAIL\n", THRESH_MEAN_COV))
cat(sprintf("   modC sensitivity     < %.0f%%  → FAIL\n", THRESH_MODC_SENS * 100))
cat(sprintf("   modC specificity     < %.0f%%  → FAIL\n", THRESH_MODC_SPEC * 100))
cat(sprintf("   Read resolution      < %.0f%%  → FAIL\n", THRESH_READ_RESOLUTION * 100))

# Helper: extract numeric value for a metric label for one sample column
get_metric <- function(sample_col, metric_label) {
  row <- merged_df[trimws(merged_df$Metric) == trimws(metric_label), sample_col, drop = TRUE]
  if (length(row) == 0 || all(is.na(row))) return(NA)
  suppressWarnings(as.numeric(row[1]))
}

qc_results <- data.frame(
  Sample_ID    = sample_cols,
  QC_Status    = NA_character_,
  Fail_Reasons = NA_character_,
  stringsAsFactors = FALSE
)

for (i in seq_along(sample_cols)) {
  s       <- sample_cols[i]
  reasons <- c()

  dup  <- get_metric(s, "Genome-mapped duplication rate")
  cov  <- get_metric(s, "Mean coverage excluding Ns")
  sens <- get_metric(s, "Methylated lambda modC sensitivity")
  spec <- get_metric(s, "Non-methylated pUC19 control modC specificity")
  res  <- get_metric(s, "Percentage of input reads that resolve")

  if (!is.na(dup)  && dup  > THRESH_DUP_RATE)       reasons <- c(reasons, sprintf("Dup rate %.1f%%",  dup  * 100))
  if (!is.na(cov)  && cov  < THRESH_MEAN_COV)        reasons <- c(reasons, sprintf("Coverage %.2fx",   cov))
  if (!is.na(sens) && sens < THRESH_MODC_SENS)        reasons <- c(reasons, sprintf("modC sens %.1f%%", sens * 100))
  if (!is.na(spec) && spec < THRESH_MODC_SPEC)        reasons <- c(reasons, sprintf("modC spec %.1f%%", spec * 100))
  if (!is.na(res)  && res  < THRESH_READ_RESOLUTION)  reasons <- c(reasons, sprintf("Read res %.1f%%",  res  * 100))

  qc_results$QC_Status[i]    <- ifelse(length(reasons) == 0, "PASS", "FAIL")
  qc_results$Fail_Reasons[i] <- ifelse(length(reasons) == 0, "", paste(reasons, collapse = "; "))
}

#===============================================================================
# 💾 Save outputs
#===============================================================================

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("💾 Writing merged_QC_summary.xlsx...\n")

wb <- createWorkbook()

# --- Sheet 1: All Metrics (Metric | 283C1 | 283C2 | ...) ---
addWorksheet(wb, "All_Metrics")
writeData(wb, "All_Metrics", merged_df)

header_style <- createStyle(
  fontColour = "#FFFFFF", fgFill = "#2C3E50",
  halign = "CENTER", textDecoration = "Bold"
)
addStyle(wb, "All_Metrics", header_style,
         rows = 1, cols = 1:ncol(merged_df), gridExpand = TRUE)
freezePane(wb, "All_Metrics", firstRow = TRUE)
setColWidths(wb, "All_Metrics", cols = 1:ncol(merged_df), widths = "auto")

# --- Sheet 2: QC Summary (one row per sample, PASS/FAIL + reasons) ---
addWorksheet(wb, "QC_Summary")
writeData(wb, "QC_Summary", qc_results)

addStyle(wb, "QC_Summary", header_style,
         rows = 1, cols = 1:ncol(qc_results), gridExpand = TRUE)

pass_rows <- which(qc_results$QC_Status == "PASS") + 1
fail_rows <- which(qc_results$QC_Status == "FAIL") + 1

if (length(pass_rows) > 0) {
  addStyle(wb, "QC_Summary", createStyle(fgFill = "#D5F5E3"),
           rows = pass_rows, cols = 1:ncol(qc_results), gridExpand = TRUE)
}
if (length(fail_rows) > 0) {
  addStyle(wb, "QC_Summary", createStyle(fgFill = "#FADBD8"),
           rows = fail_rows, cols = 1:ncol(qc_results), gridExpand = TRUE)
}

freezePane(wb, "QC_Summary", firstRow = TRUE)
setColWidths(wb, "QC_Summary", cols = 1:ncol(qc_results), widths = "auto")

saveWorkbook(wb, file.path(OUTPUT_DIR, "merged_QC_summary.xlsx"), overwrite = TRUE)

# --- Text files ---
pass_samples <- qc_results$Sample_ID[qc_results$QC_Status == "PASS"]
fail_samples <- qc_results[qc_results$QC_Status == "FAIL", c("Sample_ID", "Fail_Reasons")]

writeLines(pass_samples, file.path(OUTPUT_DIR, "QC_pass_samples.txt"))
write.table(fail_samples, file.path(OUTPUT_DIR, "QC_fail_samples.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

#===============================================================================
# 📋 Print summary
#===============================================================================

cat("\n========================================================\n")
cat(sprintf("✅ PASS: %d samples\n", length(pass_samples)))
cat(sprintf("❌ FAIL: %d samples\n", nrow(fail_samples)))
cat("\nFailed samples:\n")
if (nrow(fail_samples) > 0) {
  for (i in seq_len(nrow(fail_samples))) {
    cat(sprintf("  %-15s → %s\n", fail_samples$Sample_ID[i], fail_samples$Fail_Reasons[i]))
  }
} else {
  cat("  None\n")
}
cat("========================================================\n")
cat(sprintf("📁 Outputs saved to: %s\n", OUTPUT_DIR))
cat("   merged_QC_summary.xlsx  (Sheet 1: All_Metrics | Sheet 2: QC_Summary)\n")
cat("   QC_pass_samples.txt\n")
cat("   QC_fail_samples.txt\n")
cat("========================================================\n")
