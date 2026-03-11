df1 <- read.xlsx("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/LH_9033517_20260212_RP_Pos_CD_Results.xlsx")
df2 <- read.xlsx("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/LH_9033517_20260212_HILIC_Neg_CD_Results.xlsx")

analyze_metabolites <- function(df, filename) {
  
  # ---- STEP 1: Preliminary Filtering (Junk Removal) ----
  # Identify Columns using dots (Standard R format for Compound Discoverer headers)
  nl_cols     <- grep("Area:.*_NL",    colnames(df), value = TRUE)
  y_cols      <- grep("Area:.*_Y",     colnames(df), value = TRUE)
  qc_cols     <- grep("Area:.*_QC",    colnames(df), value = TRUE)
  blank_cols  <- grep("Area:.*_Blank", colnames(df), value = TRUE)
  all_samples <- c(nl_cols, y_cols, qc_cols) 
  
  # Identify Gap Fill columns (used to detect "missing" peaks)
  nl_gap_cols <- grep("Gap.Fill.Status:.*_NL", colnames(df), value = TRUE)
  y_gap_cols  <- grep("Gap.Fill.Status:.*_Y",  colnames(df), value = TRUE)
  
  # Calculate mean of blanks and max of samples for S/B ratio
  mean_blank <- rowMeans(df[, blank_cols, drop = FALSE], na.rm = TRUE)
  max_sample <- apply(df[, c(nl_cols, y_cols)], 1, max, na.rm = TRUE)
  
  # Frequency Filter: Calculate % of 'real' (non-gap-filled) peaks
  # Logic: CD uses 8 or 32 to indicate a peak was "filled" (missing/noise)
  is_bad_nl     <- matrix(as.matrix(df[, nl_gap_cols]) %in% c(8, 32, "8", "32"), nrow = nrow(df))
  is_bad_y      <- matrix(as.matrix(df[, y_gap_cols])  %in% c(8, 32, "8", "32"), nrow = nrow(df))
  prop_found_nl <- rowMeans(!is_bad_nl)
  prop_found_y  <- rowMeans(!is_bad_y)
  
  # Apply Filter: 3x Blank Signal AND Found in 80% of at least one group
  keep_features <- max_sample > (3 * mean_blank) & (prop_found_nl >= 0.8 | prop_found_y >= 0.8)
  df_pre <- df[keep_features, ]
  
  # ---- STEP 2: Redundancy Curation (Grouping & De-duplicating) ----
  # We group by Name/Formula to compare peaks and remove adducts/duplicates
  curated_df <- df_pre %>%
    # Temporary math to pick the "strongest" version of a peak
    mutate(Max_Area = rowSums(select(., all_of(c(nl_cols, y_cols))), na.rm = TRUE)) %>%
    # Create a grouping key: Name > Formula > RowID
    mutate(TempGroup = case_when(
      !is.na(Name) & Name != "" ~ Name,
      !is.na(Formula) & Formula != "" ~ Formula,
      TRUE ~ as.character(row_number())
    )) %>%
    # Sort: Best mzCloud score first, then highest intensity
    arrange(TempGroup, desc(`mzCloud.Best.Match.Confidence`), desc(Max_Area)) %>%
    group_by(TempGroup) %>%
    # REDUNDANCY FILTER: Keep the best hit, OR keep if it's a separate isomer (RT > 0.2 min)
    filter(row_number() == 1 | abs(`RT.[min]` - first(`RT.[min]`)) > 0.2) %>%
    ungroup()
  
  # ---- STEP 3: Naming & Unique ID Creation (Non-Redundant) ----
  # Now that we have only the rows we want, we create the final unique labels
  curated_df <- curated_df %>%
    # 1. Fill in missing names for Unknowns
    mutate(Name = ifelse(!is.na(Name) & Name != "", 
                         Name, 
                         paste0(Formula, "_", formatC(as.numeric(`m/z`), format = "f", digits = 4), 
                                "@", formatC(as.numeric(`RT.[min]`), format = "f", digits = 2)))) %>%
    # 2. Tag Isomers (Surviving duplicates of the same Name/Formula)
    group_by(Name) %>%
    mutate(Name = if(n() > 1) paste0(Name, " (Isomer ", row_number(), ")") else Name) %>%
    ungroup() %>%
    # 3. Create FeatureID (The unique computer key for the matrix)
    mutate(FeatureID = make.unique(Name)) %>%
    select(-TempGroup, -Max_Area) # Clean up temp columns
  
  # ---- STEP 4: Normalization & CV Filtering ----
  mat <- as.matrix(curated_df[, all_samples])
  rownames(mat) <- curated_df$FeatureID
  
  # Median scaling normalization
  sample_medians <- apply(mat, 2, median, na.rm = TRUE)
  mat_norm <- sweep(mat, 2, sample_medians / mean(sample_medians), "/")
  
  # QC CV Filtering: Fixed for 2 QCs
  calc_cv <- function(x) {
    x_clean <- x[!is.na(x) & x > 0]
    if(length(x_clean) < 2) return(999) 
    mean_val <- mean(x_clean)
    if(mean_val <= 0) return(999)
    (sd(x_clean) / mean_val) * 100
  }
  
  qc_cv_vals <- apply(mat_norm[, qc_cols, drop = FALSE], 1, calc_cv)
  
  # Keep only reproducible features (CV < 30%)
  df_final <- curated_df[qc_cv_vals <= 30, ]
  mat_final_norm <- mat_norm[qc_cv_vals <= 30, ]
  
  # ---- STEP 5: Export for MetaboAnalyst (Features in Rows) ----
  # 1. Create the data matrix for ONLY the samples (NL and Y)
  # Ensure the matrix only contains the sample columns, not QC columns
  export_matrix <- as.data.frame(mat_final_norm[, c(nl_cols, y_cols)])
  
  # 2. Create the group labels (Must be exactly 10 labels for 10 columns)
  groups <- c(rep("NL", length(nl_cols)), rep("Y", length(y_cols)))
  
  # 3. Create a 1-row dataframe for the labels
  # We set the colnames of this 1-row df to match our matrix exactly
  group_row <- as.data.frame(t(groups))
  colnames(group_row) <- colnames(export_matrix)
  rownames(group_row) <- "Label"
  
  # 4. Use rbind to put the Label row on top of the Data rows
  final_export <- rbind(group_row, export_matrix)
  
  # Write to CSV
  write.csv(final_export, paste0("MetaboAnalyst_Input_", filename, ".csv"), row.names = TRUE)

  # ---- STEP 6: PCA Preparation ----
  # Log transformation (standard for metabolomics PCA)
  mat_log <- log2(mat_final_norm[, c(nl_cols, y_cols)] + 1)
  pca_res <- prcomp(t(mat_log), center = TRUE, scale. = TRUE)
  
  cat("Success! Kept", nrow(df_final), "high-quality features.\n")
  return(list(data = df_final, matrix = mat_log, pca = pca_res))
}

rp_clean    <- analyze_metabolites(df1, "RP")
hilic_clean <- analyze_metabolites(df2, "HILIC")
