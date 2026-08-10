# =============================================================================
# PRE-PROCESSING STEP: GET PUBCHEM CIDs FOR COMPOUND NAME MAPPING
# =============================================================================
# Before running this R script, you need to create name_map CSV files for
# both RP_Pos and HILIC_Neg datasets using the following workflow:
#
# STEP 1: EXTRACT UNIQUE COMPOUND NAMES
#   - Open the CD Results xlsx files
#   - Copy the 'Name' column, remove blanks, remove duplicates
#   - Save as a text file (one name per line)
#
# STEP 2: GET PUBCHEM CIDs
#   - Go to: http://pubchem.ncbi.nlm.nih.gov/idexchange
#   - Input format:    Synonyms
#   - Operator type:   Same Connectivity
#   - Paste your compound names (one per line)
#   - Output type:     CID
#   - Output method:   Two column file showing input-output correspondence
#   - Compression:     None
#   - Click Submit Job and download the result
#
# STEP 3: CLEAN THE CID OUTPUT IN EXCEL
#   - Open the downloaded txt file in Excel
#   - Sort by Column A (Name) A->Z, then Column B (CID) smallest to largest
#     NOTE: Ensure CID column is formatted as NUMBER not text before sorting
#   - Remove duplicates based on Column A (Name) only
#     This keeps the lowest CID per name (most curated/canonical entry)
#   - You now have one CID per compound name
#
# STEP 4: GET MAPPED NAMES VIA METABOANALYST
#   - Go to: https://www.metaboanalyst.ca/MetaboAnalyst/upload/ConvertView.xhtml
#   - PASS 1: Input type = PubChem CID, paste all CIDs
#             Note which CIDs fail (not recognized)
#   - PASS 2: For failed CIDs, go back to original names
#             Input type = Common Name, paste failed compound names directly
#             Note which names still fail
#   - PASS 3: For remaining failures, try HMDB ID if known
#   - Download the result table after each pass
#     This gives you: Query | Match | HMDB | PubChem | ChEBI | KEGG | SMILES
#
# STEP 5: COMBINE ALL PASSES INTO ONE name_map CSV
#   - Merge all successful mappings from Pass 1, 2, and 3
#   - Final CSV must have columns: Query, Match, HMDB, PubChem, ChEBI, KEGG, SMILES
#   - Original = original messy name from CD Results
#   - Query = renamed names to get matches from Metaboanalyst
#   - Match = proper name recognized by MetaboAnalyst
#   - Save as: name_map_RP.csv and name_map_HILIC.csv
#
# NOTES:
#   - Compounds not mapped after all 3 passes are excluded from MetaboAnalyst
#     export but KEPT for QC filtering and PCA (they are real detected peaks)
#   - Greek letters must be spelled out in MetaboAnalyst (alpha, beta, gamma)
#   - Lipids with complex names (PE, PS, PG etc.) may need LipidMaps separately
#   - Placeholder names (e.g. VDNGMDHARYXCBF-UHFFFAOYSA-N, FL1AAKGM0001_a)
#     and non-human compounds (pesticides, fungal natural products, industrial
#     chemicals) should be excluded - they will not map in any database
# =============================================================================

library(openxlsx)
library(dplyr)
library(stringr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(tibble)
library(readr)

in_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Past/Boopati/input"
out_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Past/Boopati/output"

# ---- Load Data ----
df1 <- read.xlsx(file.path(in_path, "Natural_LOY_vs_YPos_RP.xlsx"))
df2 <- read.xlsx(file.path(in_path, "Natural_LOY_vs_YPos_HILIC.xlsx"))

# ---- Load Name Maps ----
nm_rp    <- read.csv(file.path(in_path, "Natural_LOY_vs_YPos_RP_name_map.csv"),    stringsAsFactors = FALSE)
nm_hilic <- read.csv(file.path(in_path, "Natural_LOY_vs_YPos_HILIC_name_map.csv"), stringsAsFactors = FALSE)

# ---- Main Processing Function ----
analyze_metabolites <- function(df, filename) {
  
  # ---- STEP 1: Preliminary Filtering ----
  nl_cols     <- grep("Area:.*_NL",    colnames(df), value = TRUE)
  y_cols      <- grep("Area:.*_Y",     colnames(df), value = TRUE)
  qc_cols     <- grep("Area:.*_QC",    colnames(df), value = TRUE)
  blank_cols  <- grep("Area:.*_Blank", colnames(df), value = TRUE)
  all_samples <- c(nl_cols, y_cols, qc_cols)
  
  # Safety checks
  stopifnot(length(nl_cols) > 0, length(y_cols) > 0,
            length(blank_cols) > 0, length(qc_cols) > 0)
  
  nl_gap_cols <- grep("Gap.Fill.Status:.*_NL", colnames(df), value = TRUE)
  y_gap_cols  <- grep("Gap.Fill.Status:.*_Y",  colnames(df), value = TRUE)
  
  mean_blank <- rowMeans(df[, blank_cols, drop = FALSE], na.rm = TRUE)
  max_sample <- apply(df[, c(nl_cols, y_cols)], 1, max, na.rm = TRUE)
  
  is_bad_nl <- matrix(as.matrix(df[, nl_gap_cols]) %in% 
                        #c(8, 32, 8192, "8", "32", "8192"), nrow = nrow(df))
                        c(8, 32, "8", "32"), nrow = nrow(df))
  is_bad_y  <- matrix(as.matrix(df[, y_gap_cols])  %in% 
                        #c(8, 32, 8192, "8", "32", "8192"), nrow = nrow(df))
                        c(8, 32, "8", "32"), nrow = nrow(df))
  prop_found_nl <- rowMeans(!is_bad_nl)
  prop_found_y  <- rowMeans(!is_bad_y)
  
  keep_features <- max_sample > (3 * mean_blank) & 
    (prop_found_nl >= 0.8 | prop_found_y >= 0.8)
  df_pre <- df[keep_features, ]
  
  # ---- STEP 2: Redundancy Curation ----
  curated_df <- df_pre %>%
    mutate(Max_Area = rowSums(select(., all_of(c(nl_cols, y_cols))), na.rm = TRUE)) %>%
    mutate(TempGroup = case_when(
      !is.na(Name) & Name != "" ~ Name,
      !is.na(Formula) & Formula != "" ~ Formula,
      TRUE ~ as.character(row_number())
    )) %>%
    arrange(TempGroup,
            desc(Tags == "Confirmed ID (High Confidence)"),
            desc(`mzCloud.Best.Match.Confidence`),
            desc(Max_Area)) %>%
    group_by(TempGroup) %>%
    filter(row_number() == 1 | abs(`RT.[min]` - first(`RT.[min]`)) > 0.2) %>%
    ungroup()
  
  # ---- STEP 3: Naming & Unique ID Creation ----
  curated_df <- curated_df %>%
    mutate(Name = ifelse(
      !is.na(`mzCloud.Best.Match.Confidence`) &
        !is.na(`mzCloud.Best.Match`) &
        (`mzCloud.Best.Match.Confidence` < 80 & `mzCloud.Best.Match` < 80),
      paste0(Name, " (Low Conf)"), Name)) %>%
    mutate(Name = ifelse(
      !is.na(Name) & Name != "" & !is.na(Formula),
      Name,
      paste0(ifelse(is.na(Formula) | Formula == "", "Unknown", Formula),
             "_", formatC(as.numeric(`m/z`), format = "f", digits = 4),
             "@", formatC(as.numeric(`RT.[min]`), format = "f", digits = 2)))) %>%
    group_by(Name) %>%
    mutate(Name = if(n() > 1) paste0(Name, " (Isomer ", row_number(), ")") 
           else Name) %>%
    ungroup() %>%
    mutate(FeatureID = make.unique(Name)) %>%
    select(-TempGroup, -Max_Area)
  
  # ---- STEP 4: QC CV Filter (on raw areas) ----
  mat <- as.matrix(curated_df[, all_samples])
  rownames(mat) <- curated_df$FeatureID
  
  calc_cv <- function(x) {
    x_clean <- x[!is.na(x) & x > 0]
    if(length(x_clean) < 2) return(999)
    mean_val <- mean(x_clean)
    if(mean_val <= 0) return(999)
    (sd(x_clean) / mean_val) * 100
  }
  
  qc_cv_vals <- apply(mat[, qc_cols, drop = FALSE], 1, calc_cv)
  df_final   <- curated_df[qc_cv_vals <= 30, ]
  
  cat("Success! Kept", nrow(df_final), "high-quality features.\n")
  return(df_final)
}

# ---- Name Replacement Function ----
replace_names <- function(df, name_map) {
  # 1. Create the lookup mapping
  lookup <- setNames(name_map$Match, name_map$Original)
  
  # 2. Preserve original name and apply new mapping
  df$Name_Original <- df$Name 
  df$Name <- ifelse(df$Name %in% names(lookup), lookup[df$Name], df$Name)
  df$Name_Mapped <- df$Name %in% name_map$Match  # TRUE = successfully mapped
  
  # 3. SPLIT: Mapped vs Unmapped
  # We only want to remove duplicates for the things we successfully named
  df_mapped   <- df[df$Name_Mapped == TRUE, ]
  df_unmapped <- df[df$Name_Mapped == FALSE, ]
  
  # 4. COLLAPSE MAPPED ENTRIES
  # We pick the one "Representative" row with the highest total area
  area_cols <- grep("Area:", colnames(df), value = TRUE)
  
  df_mapped_clean <- df_mapped %>%
    mutate(Temp_Total = rowSums(select(., all_of(area_cols)), na.rm = TRUE)) %>%
    group_by(Name) %>%
    arrange(desc(Temp_Total)) %>%
    slice(1) %>%         # Only the best row per name survives
    ungroup() %>%
    select(-Temp_Total)
  
  df_unmapped_clean <- df_unmapped %>%
    mutate(Temp_Total = rowSums(select(., all_of(area_cols)), na.rm = TRUE)) %>%
    group_by(Name) %>%
    arrange(desc(Temp_Total)) %>%
    slice(1) %>%         # Only the best row per name survives
    ungroup() %>%
    select(-Temp_Total)
  
  # 5. RECOMBINE
  # This puts the unique mapped names back together with all your unknowns
  final_df <- rbind(df_mapped_clean, df_unmapped_clean)
  
  cat("Consolidated", nrow(df_mapped), "rows into", nrow(df_mapped_clean), 
      "unique metabolites. Kept", nrow(df_unmapped_clean), "unmapped features.\n")
  
  return(final_df)
}

# ---- Export for MetaboAnalyst (RAW areas, mapped names only) ----
export_for_metaboanalyst <- function(df_final, filename) {
  
  nl_cols    <- grep("Area:.*_NL", colnames(df_final), value = TRUE)
  y_cols     <- grep("Area:.*_Y",  colnames(df_final), value = TRUE)
  area_cols  <- c(nl_cols, y_cols)
  
  # Drop unmapped compounds before export
  df_mapped <- df_final[df_final$Name_Mapped == TRUE, ]
  
  # Use Name (proper metabolite name from metaboanalyst) as row identifier
  export_mat <- as.data.frame(df_mapped[, area_cols])
  #rownames(export_mat) <- df_mapped$FeatureID
  rownames(export_mat) <- df_mapped$Name
  
  # Replace 0/NA with 1 (MetaboAnalyst handles imputation)
  export_mat[is.na(export_mat) | export_mat == 0] <- 1
  
  # Label row
  group_labels <- ifelse(grepl("_NL_", colnames(export_mat)), "NL", "Y")
  label_row    <- as.data.frame(t(group_labels))
  colnames(label_row)  <- colnames(export_mat)
  rownames(label_row)  <- "Label"
  
  final_export <- rbind(label_row, export_mat)
  
  write.csv(final_export,
            file.path(out_path, paste0("MetaboAnalyst_Ready_", filename, ".csv")),
            row.names = TRUE)
  
  cat("Exported", nrow(export_mat), "mapped compounds to MetaboAnalyst_Ready_",
      filename, ".csv\n")
}

#rp_clean    <- analyze_metabolites(df1, "RP")
#hilic_clean <- analyze_metabolites(df2, "HILIC")

rp_clean <- df1
hilic_clean <- df2
rp_metaboanalyst <- replace_names(rp_clean, nm_rp)
hilic_metaboanalyst <- replace_names(hilic_clean, nm_hilic)

export_for_metaboanalyst(rp_metaboanalyst,    "RP_Pos")
export_for_metaboanalyst(hilic_metaboanalyst, "HILIC_Neg")

# ---- PCA (log2 of raw areas for visualization only) ----
plot_pca <- function(df_final, title) {
  
  library(ggplot2)
  
  nl_cols <- grep("Area:.*_NL", colnames(df_final), value = TRUE)
  y_cols  <- grep("Area:.*_Y",  colnames(df_final), value = TRUE)
  mat     <- as.matrix(df_final[, c(nl_cols, y_cols)])
  rownames(mat) <- df_final$FeatureID
  
  # Log2 transform for PCA only
  mat_log <- log2(mat + 1)
  mat_log <- mat_log[apply(mat_log, 1, var) > 0, ]
  pca_res <- prcomp(t(mat_log), center = TRUE, scale. = TRUE)
  
  df_pca        <- as.data.frame(pca_res$x)
  df_pca$Group  <- ifelse(grepl("_NL_", rownames(df_pca)), "NL", "Y")
  vars          <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
  
  ggplot(df_pca, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 5, alpha = 0.8) +
    stat_ellipse(level = 0.95) +
    theme_minimal() +
    scale_color_manual(values = c("NL" = "#2C7BB6", "Y" = "#D7191C")) +
    labs(title = title,
         subtitle = paste(nrow(mat_log), "features"),
         x = paste0("PC1 (", vars[1], "%)"),
         y = paste0("PC2 (", vars[2], "%)"))
}

library(gridExtra)
p1 <- plot_pca(rp_metaboanalyst,    "PCA: RP-Pos")
p2 <- plot_pca(hilic_metaboanalyst, "PCA: HILIC-Neg")
grid.arrange(p1, p2, ncol = 2)


#========
# Merge HILIC and RP normalized values
#=====

# 1. Define the file paths
rp_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/Downloads/rp_data_normalized.csv"
hilic_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/Downloads/hilic_data_normalized.csv"

# ---- Load Name Maps ----
in_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Past/Boopati/input"
nm_rp    <- read.csv(file.path(in_path, "Natural_LOY_vs_YPos_RP_name_map.csv"),    stringsAsFactors = FALSE)
nm_hilic <- read.csv(file.path(in_path, "Natural_LOY_vs_YPos_HILIC_name_map.csv"), stringsAsFactors = FALSE)

# 2. Read the CSV files
rp_norm <- read.csv(rp_path)
hilic_norm <- read.csv(hilic_path)
colnames(rp_norm)[1] <- "Name"
colnames(hilic_norm)[1] <- "Name"

label_row <- rp_norm %>% filter(Name == "Label") %>% rename(HMDB = Name)

# 3. Combine the datasets row-wise (stacking the features)
combined_data <- bind_rows(rp_norm, hilic_norm)

# Dynamically get the name of the first column (e.g., "Metabolite" or "Feature")
metabolite_col <- colnames(combined_data)[1]

# 4. Calculate row-wise SD and filter out duplicates with the lower SD
cleaned_data <- combined_data %>%
  rowwise() %>%
  # Calculate SD across all sample columns (excludes the first name column)
  mutate(row_sd = sd(c_across(-1), na.rm = TRUE)) %>%
  ungroup() %>%
  # Group by the metabolite name column
  group_by(across(all_of(metabolite_col))) %>%
  # Keep only the single row per group that has the maximum SD
  slice_max(row_sd, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  # Drop the temporary SD column
  select(-row_sd)

# 6. Replace names with HMDB or KEGG ID from name map
nm_full <- dplyr::bind_rows(nm_rp, nm_hilic) %>% 
  dplyr::select(Match, HMDB) %>% #, KEGG, PubChem) %>%
  dplyr::distinct()

cleaned_data <- cleaned_data %>% 
  dplyr::left_join(nm_full, by = c("Name" = "Match")) %>%
  dplyr::select(HMDB, everything(), -Name) %>%
  dplyr::arrange(HMDB) %>%
  dplyr::filter(!is.na(HMDB) & HMDB != "")

# 7. Put Label row at 2nd row manually
cleaned_data <- dplyr::bind_rows(label_row, cleaned_data)

# 8. Export the clean, merged dataset back to your Downloads folder
output_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/Downloads/merged_hilic_rp_norm.csv"
write_csv(cleaned_data, output_path)

cat("Success! Merged dataset saved to:", output_path, "\n")

# ==============================================================================
# METABOANALYST STATISTICAL WORKFLOW SUMMARY
# DATASET: HILIC POSITIVE MODE (Cleaned via 80% Rule & 30% QC CV Filter)
# ==============================================================================

## 1. DATA UPLOAD & FILTERING
# - Module: Statistical Analysis [1-factor]
# - Integrity Check: Passed
# - Filters: ALL DISABLED (0%) to preserve pre-cleaned R metabolites.

## 2. NORMALIZATION SETTINGS
# - Sample Normalization: Normalization by Median
# - Data Transformation: Log 2 Transformation
# - Data Scaling: Pareto Scaling
# - Result: Bell-shaped density and aligned sample medians

## 3. DIFFERENTIAL EXPRESSION (Limma)
# - Fold Change: Limma Method, Threshold = 1 (Full list export)
# - T-Test: Limma Moderated Statistics, P-value Threshold = 1 (Full list export)
# - Output: Merged XLSX with Log2FC, Raw P-value, and FDR.

## 4. VISUALIZATION PARAMETERS
# - Volcano Plot: Classic style, 1.5 FC threshold, 0.05 P-value (Raw)
# - PCA: 2D Scores Plot with 95% Confidence Regions enabled.
# - Heatmap: Ward clustering, Euclidean distance, Autoscaled, top 100, 
# - Export Format: PDF, 300 DPI, 12x12 inches for publication quality.

## 5. PATHWAY SETTINGS
# Get normalized values from download.zip.
# Convert names to HMDID and upload to 
# - Module: Pathway Analysis


# ==============================================================================
# QUADRANT PLOT
# ==============================================================================


### COmpare results of human and mouse
rp <- openxlsx::read.xlsx(file.path(out_path, "RP_Yneg_Ypos/rp_DMA.xlsx"))
hilic <- openxlsx::read.xlsx(file.path(out_path, "HILIC_Yneg_Ypos/hilic_DMA.xlsx"))
mb49 <- openxlsx::read.xlsx(file.path(out_path, "MB49_Yneg_Ypos/mb49_DMA.xlsx"))

# 1. Combine them and label them
human <- dplyr::bind_rows(rp %>% mutate(Method = "RP"),
                          hilic %>% mutate(Method = "HILIC")) %>%
  # 2. Group by metabolite name
  group_by(SYMBOL) %>%
  # 3. Sort by p.value (lowest first)
  arrange(p.value, .by_group = TRUE) %>%
  # 4. Keep only the first (most significant) row for each metabolite
  slice(1) %>%
  ungroup()

# Merge Human and Mouse by SYMBOL
# inner_join keeps only metabolites found in both species
# full_join would keep everything, even if missing in one
comparison <- inner_join(human, 
                         mb49, 
                         by = "SYMBOL", 
                         suffix = c("_Human", "_Mouse"))

# Add a column to categorize the "Conserved" status
comparison <- comparison %>%
  mutate(Status = case_when(`log2(FC)_Human` > 0 & `log2(FC)_Mouse` > 0 ~ "Conserved Up",
                            `log2(FC)_Human` < 0 & `log2(FC)_Mouse` < 0 ~ "Conserved Down",
                            TRUE ~ "Divergent"))

# See a quick summary of the counts
table(comparison$Status)

library(ggplot2)
library(ggrepel)

# 1. First, create the label list for the conserved metabolites
# Since you have 12 conserved (8 down, 4 up), we can label all of them
label_metabolites <- comparison %>% 
  filter(Status != "Divergent") %>% 
  pull(SYMBOL)

# 2. Define your specific color palette
status_colors <- c(
  "Conserved Up"   = "#C0392B", 
  "Conserved Down" = "#1A6FA8", 
  "Divergent"      = "grey82"
)

# 3. Create the Plot
p_final <- ggplot(comparison, aes(x = `log2(FC)_Human`, y = `log2(FC)_Mouse`)) +
  # Quadrant lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.35) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.35) +
  
  # Layer 1: Background (Divergent/Ambiguous)
  geom_point(data  = dplyr::filter(comparison, Status == "Divergent"),
             color = "grey82", size = 1.5, alpha = 0.4) +
  
  # Layer 2: Foreground (Conserved)
  geom_point(data  = dplyr::filter(comparison, Status != "Divergent"),
             aes(color = Status), size = 2.5, alpha = 0.9) +
  
  # Layer 3: Labels
  ggrepel::geom_text_repel(
    data         = dplyr::filter(comparison, SYMBOL %in% label_metabolites),
    aes(label    = SYMBOL, color = Status),
    size         = 3, fontface = "bold.italic", max.overlaps = Inf,
    box.padding  = 0.5, point.padding = 0.3,
    segment.size = 0.25, segment.alpha = 0.5, 
    show.legend  = FALSE
  ) +
  
  # Annotations (adjust x/y based on your actual data range)
  annotate("text", x =  2, y =  2.5, label = "Conserved UP",
           color = "#C0392B", size = 3.5, fontface = "bold", hjust = 1) +
  annotate("text", x = -2, y = -2.5, label = "Conserved DOWN",
           color = "#1A6FA8", size = 3.5, fontface = "bold", hjust = 0) +
  
  scale_color_manual(values = status_colors) +
  # Typical metabolomics range is -3 to 3; adjust if your data is wider
  scale_x_continuous(limits = c(-3, 3)) + 
  scale_y_continuous(limits = c(-3, 3)) +
  
  labs(x     = "Log\u2082FC (Human Patient Samples)",
       y     = "Log\u2082FC (Mouse MB49 Model)",
       color = NULL,
       title = "Metabolic Concordance: Human vs. Mouse Model") +
  
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(size = 12, face = "bold"),
    axis.title      = element_text(size = 10),
    legend.position = "bottom",
    panel.grid      = element_blank()
  )

# 4. Save as PDF
ggsave("Human_Mouse_Metabolomics_Concordance.pdf", p_final, width = 6, height = 6.5)

# ==============================================================================
# GLUCODE GENES
# ==============================================================================

# Filter for Glycolysis and keep only the Human/Mouse mapping columns
glucose_metabolites <- c("D-Glucose", 
                         "D-Fructose 1-6-bisphosphate", 
                         "D-Glyceraldehyde 3-phosphate",
                         "1-3-Bisphosphoglycerate",
                         "2/3-Phospho-D-glycerate",
                         "Phosphoenolpyruvate",
                         "Pyruvate",
                         "Lactate",
                         "Maltose",
                         "NAD+",
                         "NADH",
                         "Succinate",
                         "Fumarate",
                         "Malate",
                         "Adenosine")

# Function to remove all non-alphanumeric characters and lowercase
clean_name <- function(x) {
  x <- gsub("[^[:alnum:]]", "", x) # Removes -, /, +, spaces, etc.
  x <- tolower(x)
  return(x)
}

# Create the cleaned reference list
#cleaned_ref <- clean_name(glucose_metabolites)

human_glucose <- human %>%
  dplyr::filter(SYMBOL %in% glucose_metabolites) %>%
  select(SYMBOL, logFC_Human = `log2(FC)`) # Ensure the column name matches your data (e.g., logFC)

# 3. Process Mouse Data (MB49)
mouse_glucose <- mb49 %>%
  dplyr::filter(SYMBOL %in% glucose_metabolites) %>%
  select(SYMBOL, logFC_Mouse = `log2(FC)`)

combined_glucose <- full_join(human_glucose, mouse_glucose, by = "SYMBOL")

# 2. Prepare the matrix for the heatmap
# We convert the SYMBOL column to row names
plot_matrix <- combined_glucose %>%
  # If there are duplicate symbols, keep the first one
  distinct(SYMBOL, .keep_all = TRUE) %>% 
  remove_rownames() %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

# 1. Convert any Infinite values to NA (log of 0 creates -Inf)
plot_matrix[is.infinite(plot_matrix)] <- NA

# 2. Identify and remove rows that are ALL NA
# Clustering cannot handle a row where every single column is missing
rows_to_keep <- rowSums(is.na(plot_matrix)) < ncol(plot_matrix)
filtered_matrix <- plot_matrix[rows_to_keep, , drop = FALSE]

# 3. Check: If you still have rows with SOME NA, pheatmap's default clustering 
# might still complain. You have two choices:

# Choice A: Turn off row clustering (Easiest/Fastest)
pheatmap(filtered_matrix, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         na_col = "grey95",
         display_numbers = TRUE,
         main = "Metabolites (No Clustering)")

# Choice B: Keep clustering but handle the NAs (Best for patterns)
# We replace NAs with 0 just for the clustering calculation
pheatmap(filtered_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",
         na_col = "grey95",
         display_numbers = TRUE,
         main = "Metabolites (Clustered)")

####################

















path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/LH_9033517_20260212 (56 naturals)-selected/input"
df1 <- read.xlsx(file.path(path, "LH_9033517_20260212_RP_Pos_CD_Results.xlsx"))
df2 <- read.xlsx(file.path(path, "LH_9033517_20260212_HILIC_Neg_CD_Results.xlsx"))

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
  # Logic: CD uses 8, 32, or 8192 to indicate a peak was "filled" (missing/noise/low-quality)
  is_bad_nl <- matrix(as.matrix(df[, nl_gap_cols]) %in% c(8, 32, 8192, "8", "32", "8192"), nrow = nrow(df))
  is_bad_y  <- matrix(as.matrix(df[, y_gap_cols])  %in% c(8, 32, 8192, "8", "32", "8192"), nrow = nrow(df))
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
    # Prioritize Confirmed IDs first, then mzCloud confidence, then intensity
    arrange(TempGroup, 
            desc(Tags == "Confirmed ID (High Confidence)"), 
            desc(`mzCloud.Best.Match.Confidence`), 
            desc(Max_Area)) %>% 
    group_by(TempGroup) %>%
    # REDUNDANCY FILTER: Keep the best hit, OR keep if it's a separate isomer (RT > 0.2 min)
    filter(row_number() == 1 | abs(`RT.[min]` - first(`RT.[min]`)) > 0.2) %>%
    ungroup()

  # ---- STEP 3: Naming & Unique ID Creation (Non-Redundant) ----
  # Now that we have only the rows we want, we create the final unique labels
  curated_df <- curated_df %>%
    # 1. Flag Low Confidence hits (based on the 80% rule in your documentation)
    mutate(Name = ifelse(!is.na(`mzCloud.Best.Match.Confidence`) & 
                           !is.na(`mzCloud.Best.Match`) &
                           (`mzCloud.Best.Match.Confidence` < 80 | `mzCloud.Best.Match` < 80), 
                         paste0(Name, " (Low Conf)"), 
                         Name)) %>%
    # 2. Fill in missing names for Unknowns
    mutate(Name = ifelse(!is.na(Name) & Name != "" & !is.na(Formula), 
                         Name, 
                         paste0(ifelse(is.na(Formula) | Formula == "", "Unknown", Formula), 
                                "_", formatC(as.numeric(`m/z`), format = "f", digits = 4), 
                                "@", formatC(as.numeric(`RT.[min]`), format = "f", digits = 2)))) %>%
    # 3. Tag Isomers (Surviving duplicates of the same Name/Formula)
    group_by(Name) %>%
    mutate(Name = if(n() > 1) paste0(Name, " (Isomer ", row_number(), ")") else Name) %>%
    ungroup() %>%
    # 4. Create FeatureID (The unique computer key for the matrix)
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
  write.csv(final_export, paste0("Processed_", filename, ".csv"), row.names = TRUE)

  # ---- STEP 6: PCA Preparation ----
  # Log transformation (standard for metabolomics PCA)
  mat_log <- log2(mat_final_norm[, c(nl_cols, y_cols)] + 1)
  # Remove rows with zero variance to prevent PCA crash
  mat_log <- mat_log[apply(mat_log, 1, var) > 0, ]
  pca_res <- prcomp(t(mat_log), center = TRUE, scale. = TRUE)
  
  cat("Success! Kept", nrow(df_final), "high-quality features.\n")
  return(list(data = df_final, matrix = mat_log, pca = pca_res))
}

rp_clean    <- analyze_metabolites(df1, "RP")
hilic_clean <- analyze_metabolites(df2, "HILIC")


# Function to plot the PCA object already sitting in your list
plot_existing_pca <- function(processed_list, title) {
  
  library(ggplot2)
  library(gridExtra)
  
  # 1. Extract the PCA results and the log matrix
  pca_res <- processed_list$pca
  mat_log <- processed_list$matrix
  
  # 2. Prepare plotting data (PC1 and PC2)
  df_pca <- as.data.frame(pca_res$x)
  
  # 3. Assign Groups based on sample names (NL vs Y)
  # This looks at the row names of your matrix (the sample names)
  df_pca$Group <- ifelse(grepl("_NL_", rownames(df_pca)), "NL", "Y")
  
  # 4. Calculate % variance for the axes
  vars <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
  
  # 5. Create the ggplot
  ggplot(df_pca, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 5, alpha = 0.8) +
    stat_ellipse(level = 0.95) +
    theme_minimal() +
    scale_color_manual(values = c("NL" = "#2C7BB6", "Y" = "#D7191C")) +
    labs(title = title,
         subtitle = paste(nrow(mat_log), "High-Quality Features"),
         x = paste0("PC1 (", vars[1], "%)"),
         y = paste0("PC2 (", vars[2], "%)"))
}

# Use the lists you already created
p1 <- plot_existing_pca(rp_clean, "PCA: RP-Pos")
p2 <- plot_existing_pca(hilic_clean, "PCA: HILIC-Neg")

# Arrange side-by-side
grid.arrange(p1, p2, ncol = 2)


export_for_metaboanalyst <- function(processed_list, filename) {
  # 1. Use the 'data' part (linear scale) instead of 'matrix' (log scale)
  # This ensures we sum areas correctly.
  df_raw <- processed_list$data
  
  # 2. Identify sample columns (Area: ... NL/Y)
  area_cols <- grep("Area:.*(_NL_|_Y_)", colnames(df_raw), value = TRUE)
  
  # 3. Clean Names (Isomers, Tags, Encoding)
  df_raw$CleanName <- df_raw$Name %>%
    # 1. Remove "Low Conf" regardless of capitalization or extra spaces
    str_replace_all(regex(" ?\\(low conf\\)?", ignore_case = TRUE), "") %>%
    
    # 2. Remove "Isomer X" regardless of capitalization
    str_replace_all(regex(" ?\\(isomer \\d+\\)?", ignore_case = TRUE), "") %>%
    
    # 3. Fix the specific "messy" characters and suffixes
    str_replace_all("Î¥", "gamma") %>%
    str_replace_all("Î²", "beta") %>%
    str_replace_all("Î±", "alpha") %>%
    
    # 4. Final cleanup: Remove any trailing dashes or double spaces
    str_replace_all("\\s+", " ") %>%
    str_trim()
  
  # 4. Consolidate Isomers by SUMMING the linear areas
  df_consolidated <- df_raw %>%
    group_by(CleanName) %>%
    summarise(across(all_of(area_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
    as.data.frame()
  
  # Set the cleaned names as row names and remove the helper column
  rownames(df_consolidated) <- df_consolidated$CleanName
  df_consolidated$CleanName <- NULL
  
  # 5. Normalization & Log Transformation (Standard for MetaboAnalyst)
  # Replace 0s with a small value to allow log
  df_consolidated[df_consolidated == 0] <- 1
  # Log2 transform the summed areas
  mat_final <- log2(df_consolidated)
  
  # 6. Create the 'Label' row
  group_labels <- ifelse(grepl("_NL_", colnames(mat_final)), "NL", "Y")
  label_row <- as.data.frame(t(group_labels))
  colnames(label_row) <- colnames(mat_final)
  rownames(label_row) <- "Label"
  
  # 7. Combine Label row with Data
  final_export <- rbind(label_row, mat_final)
  
  # 8. Save to CSV
  write.csv(final_export, paste0("MetaboAnalyst_Ready_", filename, ".csv"), row.names = TRUE)
  
  cat("Success! Exported", nrow(mat_final), "unique compounds to:", 
      paste0("MetaboAnalyst_Ready_", filename, ".csv\n"))
}

export_for_metaboanalyst(rp_clean, "RP_Pos")
export_for_metaboanalyst(hilic_clean, "HILIC_Neg")