# =============================================================================
# 00_create_gene_lengths.R
#
# Purpose: Compute a gene-level length table from a GTF, for use when TPM
#          must be derived from gene-level raw counts that have no per-
#          transcript effective-length information (e.g. IMvigor — counts
#          only, no Salmon/kallisto quantification available).
#
# Method : union-of-exons length per gene (merge overlapping exons across
#          all transcripts of a gene, sum resulting non-overlapping widths).
#          This is a standard fallback length for TPM when only gene-level
#          counts are available. It is NOT equivalent to Salmon's effective
#          length (which also accounts for the fragment-size distribution of
#          the sequencing library) — expect IMvigor TPM computed with this
#          length table to differ somewhat in methodology from TCGA (GDC/
#          GENCODE length) and GEO (depositor's own pipeline).
#
# Run once per genome build/release. Re-run only if the GTF version changes.
# Designed to be reused for any future gene-level-counts-only dataset, not
# just IMvigor.
#
# Output : gene_lengths_{build}.{release}.csv   (columns: gene_id, length)
#          saved alongside tx2gene_{build}.{release}.csv
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(rtracklayer)
  library(GenomicRanges)
})

# =============================================================================
# INPUTS / OUTPUTS  — HPC paths
# =============================================================================

GTF_PATH   <- "/home/kailasamms/resources/genomes/Homo_sapiens.GRCh38.115.gtf"
OUTPUT_DIR <- "/home/kailasamms/scripts/nextflow/resources/Datasets"

# Derive output filename from GTF filename, same convention as
# salmon_create_tx2gene() — e.g. "Homo_sapiens.GRCh38.115.gtf" -> "GRCh38.115"
base_name <- basename(GTF_PATH)
ver_name  <- gsub("^[^.]+\\.|\\.gtf$", "", base_name)
csv_name  <- paste0("gene_lengths_", ver_name, ".csv")
out_path  <- file.path(OUTPUT_DIR, csv_name)

cat("GTF_PATH   :", GTF_PATH, "\n")
cat("Output     :", out_path, "\n")

if (!file.exists(GTF_PATH)) {
  stop("GTF not found: ", GTF_PATH)
}

# =============================================================================
# STEP 1: Import exon records from GTF
# =============================================================================
# Same import approach as salmon_create_tx2gene() — read only exon-level
# rows directly, skip building a full TxDb (not needed just for lengths).

cat("\nImporting exon records from GTF...\n")
t0 <- Sys.time()

exon_data <- rtracklayer::import(GTF_PATH,
                                 feature.type = "exon",
                                 colnames     = c("type", "gene_id"))

cat("  Imported", length(exon_data), "exon records in",
    round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

# =============================================================================
# STEP 2: Union-of-exons length per gene
# =============================================================================
# split() groups exon ranges per gene. reduce() merges any overlapping exons
# within a gene (from alternative transcripts) into a non-redundant set of
# intervals. sum(width()) gives total exonic footprint.

cat("\nComputing union-of-exons length per gene...\n")

exons_by_gene <- GenomicRanges::split(exon_data, exon_data$gene_id)
gene_lengths  <- sum(GenomicRanges::width(GenomicRanges::reduce(exons_by_gene)))

gene_length_df <- data.frame(
  gene_id = names(gene_lengths),
  length  = as.numeric(gene_lengths),
  stringsAsFactors = FALSE
)

cat("  Genes with computed length:", nrow(gene_length_df), "\n")
cat("  Length range (bp): [", min(gene_length_df$length), ",",
    max(gene_length_df$length), "]\n")

# =============================================================================
# STEP 3: Save
# =============================================================================

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

readr::write_csv(gene_length_df, out_path)

cat("\nSaved:", out_path, "\n")
cat("Done.\n")
