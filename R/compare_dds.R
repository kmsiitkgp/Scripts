
# ---- Full DESeqDataSet comparison ----

compare_dds <- function(dds1_path, dds2_path) {
  cat("\n============================================================\n")
  cat("  ULTIMATE DESEQ2 COMPARISON: STATS & STRUCTURE\n")
  cat("============================================================\n")
  
  # 1. Load Files
  dds1 <- readRDS(dds1_path)
  dds2 <- readRDS(dds2_path)
  
  # 1. Metadata & Design
  v1 <- metadata(dds1)$version; v2 <- metadata(dds2)$version
  d1 <- deparse(design(dds1)); d2 <- deparse(design(dds2))
  cat(sprintf("[1] PIPELINE VERSIONS:\n    dds1=%s | dds2=%s\n", v1, v2))
  cat(sprintf("    Designs Match: %s\n", identical(d1, d2)))
  
  # 2. colData Classes
  c1 <- sapply(colData(dds1), class); c2 <- sapply(colData(dds2), class)
  mismatches <- names(c1)[c1 != c2]
  if(length(mismatches) > 0) {
    cat("\n[2] !!! COLDATA CLASS MISMATCHES (Found variables with different types):\n")
    print(data.frame(column=mismatches, dds1=c1[mismatches], dds2=c2[mismatches]))
  } else cat("\n[2] COLDATA: All metadata columns have matching classes.\n")
  
  # 3. Statistical Drift (The "Does it matter?" check)
  res_cols <- grep("WaldStatistic", colnames(mcols(dds1)), value = TRUE)
  if (length(res_cols) > 0) {
    # Use the first one but label it clearly
    stat_diff <- abs(mcols(dds1)[[res_cols[1]]] - mcols(dds2)[[res_cols[1]]])
    cat(sprintf("\n[3] STATISTICAL DRIFT (Difference in Wald Z-scores):\n"))
    cat(sprintf("    (Checking column: %s)\n", res_cols[1]))
    print(round(quantile(stat_diff, probs=c(0, 0.5, 0.9, 0.99, 1), na.rm=TRUE), 5))
    cat("    Note: A difference > 2.0 at the 99% quantile means your top genes are moving.\n")
  }
  
  # 4. Normalization Factors
  nf1 <- normalizationFactors(dds1)
  nf2 <- normalizationFactors(dds2)
  
  if(!is.null(nf1) && !is.null(nf2)) {
    # Find genes present in both runs
    common_genes <- intersect(rownames(nf1), rownames(nf2))
    
    # Check if we lost any genes
    lost_in_run2 <- setdiff(rownames(nf1), rownames(nf2))
    lost_in_run1 <- setdiff(rownames(nf2), rownames(nf1))
    
    if(length(lost_in_run1) > 0 || length(lost_in_run2) > 0) {
      cat("\n[!] WARNING: Run 1 and Run 2 have different gene counts!")
      cat("\n    Genes only in Run 1:", length(lost_in_run2))
      cat("\n    Genes only in Run 2:", length(lost_in_run1))
    }
    
    # Subset to common genes to allow subtraction
    nf1_common <- nf1[common_genes, ]
    nf2_common <- nf2[common_genes, ]
    
    val_to_check <- abs(nf1_common - nf2_common)
    cat("\n[4] NORMALIZATION FACTOR DIFFERENCES (Common Genes Quantiles):\n")
    print(round(quantile(val_to_check, na.rm=TRUE), 6))
    
  } else {
    # Fallback to Size Factors if Norm Factors don't exist
    sf1 <- sizeFactors(dds1)
    sf2 <- sizeFactors(dds2)
    if(!is.null(sf1) && !is.null(sf2)) {
      cat("\n[4] SIZE FACTOR DIFFERENCES (Quantiles):\n")
      print(round(quantile(abs(sf1 - sf2), na.rm=TRUE), 6))
    }
  }
  
  # 5. Outlier Mismatches
  out1 <- mcols(dds1)$dispOutlier; out2 <- mcols(dds2)$dispOutlier
  mismatch_idx <- which((out1 != out2) | (is.na(out1) != is.na(out2)))
  cat(sprintf("\n[5] GENES WITH DIFFERENT OUTLIER STATUS: %d\n", length(mismatch_idx)))
  if(length(mismatch_idx) > 0) {
    cat("    First 5 mismatched genes:\n")
    print(head(rownames(dds1)[mismatch_idx], 5))
  }
  
  cat("\n============================================================\n")
}

compare_txi <- function(txi1_path, txi2_path) {
  cat("\n============================================================\n")
  cat("  TXIMPORT (TXI) INPUT DIVE + QUANTILES\n")
  cat("============================================================\n")
  
  # 1. Load Files
  txi1 <- readRDS(txi1_path)
  txi2 <- readRDS(txi2_path)
  
  # 1. Transcript Length Quantiles (The "Smoking Gun")
  len_diff <- abs(txi1$length - txi2$length)
  cat("[1] TRANSCRIPT LENGTH DIFFERENCES (Base Pairs):\n")
  print(round(quantile(len_diff, probs=c(0, 0.5, 0.9, 0.99, 1)), 4))
  
  # 2. Count Quantiles
  cnt_diff <- abs(txi1$counts - txi2$counts)
  cat("\n[2] RAW COUNT DIFFERENCES (Reads):\n")
  print(round(quantile(cnt_diff, probs=c(0, 0.5, 0.9, 0.99, 1)), 4))
  
  # 3. Root Cause Gene
  max_gene <- rownames(txi1$length)[which.max(rowMeans(len_diff))]
  cat(sprintf("\n[3] MOST AFFECTED GENE (Largest Length Change):\n"))
  cat(sprintf("    Gene: %s (Average change of %.2f bp)\n", 
              max_gene, max(rowMeans(len_diff))))
  
  cat("\n============================================================\n")
}

compare_tx2gene <- function(tx2gene1_path, tx2gene2_path) {
  cat("\n============================================================\n")
  cat("  TX2GENE MAPPING FILE COMPARISON\n")
  cat("============================================================\n")
  
  # 1. Load Files
  df1 <- read.csv(tx2gene1_path, stringsAsFactors = FALSE)
  df2 <- read.csv(tx2gene2_path, stringsAsFactors = FALSE)
  
  cat(sprintf("[1] File Stats:\n"))
  cat(sprintf("    File 1: %s (%d rows)\n", basename(tx2gene1_path), nrow(df1)))
  cat(sprintf("    File 2: %s (%d rows)\n", basename(tx2gene2_path), nrow(df2)))
  
  # 2. Check for ID Overlap
  tx1 <- unique(df1$transcript_id)
  tx2 <- unique(df2$transcript_id)
  
  common_tx <- intersect(tx1, tx2)
  only_in_1 <- setdiff(tx1, tx2)
  only_in_2 <- setdiff(tx2, tx1)
  
  cat(sprintf("\n[2] Transcript ID Overlap:\n"))
  cat(sprintf("    Shared Transcripts: %d\n", length(common_tx)))
  cat(sprintf("    Unique to File 1:   %d\n", length(only_in_1)))
  cat(sprintf("    Unique to File 2:   %d\n", length(only_in_2)))
  
  # 3. Check for Mapping Shifts (The "Smoking Gun")
  # Join them on transcript_id to see if they point to the same gene
  merged <- merge(df1[,c("transcript_id", "gene_id")], 
                  df2[,c("transcript_id", "gene_id")], 
                  by = "transcript_id", suffixes = c("_f1", "_f2"))
  
  mismatched_mappings <- merged[merged$gene_id_f1 != merged$gene_id_f2, ]
  
  cat(sprintf("\n[3] Mapping Consistency:\n"))
  if (nrow(mismatched_mappings) == 0) {
    cat("    SUCCESS: All shared transcripts map to the same Gene IDs.\n")
  } else {
    cat(sprintf("    !!! ALERT: %d transcripts map to DIFFERENT Gene IDs between files.\n", 
                nrow(mismatched_mappings)))
    cat("    Showing first 5 mismatches:\n")
    print(head(mismatched_mappings, 5))
  }
  
  # 4. Check for Weighted Average Impact
  # Count how many transcripts are assigned to each gene
  counts1 <- as.data.frame(table(df1$gene_id))
  counts2 <- as.data.frame(table(df2$gene_id))
  colnames(counts1) <- c("gene_id", "n_tx_f1")
  colnames(counts2) <- c("gene_id", "n_tx_f2")
  
  count_merge <- merge(counts1, counts2, by = "gene_id")
  diff_tx_per_gene <- count_merge[count_merge$n_tx_f1 != count_merge$n_tx_f2, ]
  
  cat(sprintf("\n[4] Gene Architecture (Transcripts per Gene):\n"))
  cat(sprintf("    Genes with different transcript counts: %d\n", nrow(diff_tx_per_gene)))
  if (nrow(diff_tx_per_gene) > 0) {
    cat("    (This is why your average gene lengths shifted!)\n")
    cat("    Showing genes with largest change in transcript number:\n")
    diff_tx_per_gene$delta <- abs(diff_tx_per_gene$n_tx_f1 - diff_tx_per_gene$n_tx_f2)
    print(head(diff_tx_per_gene[order(-diff_tx_per_gene$delta), ], 5))
  }
  
  cat("\n============================================================\n")
}


# ---- Run full comparison ----

dds1_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/DESeq2_dds.rds"
dds2_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/RNASeq_Xenograft_Manish_22Rv1_Xeno/Human_split/08.DESeq2/DESeq2_dds.rds"

txi1_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Salmon_txi.rds"
txi2_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/RNASeq_Xenograft_Manish_22Rv1_Xeno/Human_split/03.Salmon/Salmon_txi_Human_split.rds"

tx2gene1_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/tx2gene_GRCh38.115.csv"
tx2gene2_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/RNASeq_Xenograft_Manish_22Rv1_Xeno/tx2gene_GRCh38.115.csv"

compare_dds(dds1_path, dds2_path)
compare_txi(txi1_path, txi2_path)
compare_tx2gene(tx2gene1_path, tx2gene2_path)