#!/usr/bin/env Rscript

# Read and store variables from command line interface (CLI)
cli <- base::commandArgs(trailingOnly = TRUE) 
args <- base::strsplit(x = cli, split = "=", fixed = TRUE)

for (e in args){
  argname <- e[1]
  argval <- e[2]
  assign(argname, argval)
}

# ==== 🧬 BULK RNA-SEQ WORKFLOW ====

# ⚙️️ Project Setup 

if (.Platform$OS.type == "windows") {
  parent_dir  <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data"
  gmt_dir     <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GSEA_genesets"
  scripts_dir <- NULL
  script_file <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/Scripts/R/Custom_Functions.R"
} else {  # Linux/macOS (e.g., HPC)
  parent_dir  <- "/hpc/home/kailasamms/scratch"
  gmt_dir     <- "/hpc/home/kailasamms/projects/GSEA_genesets"
  scripts_dir <- "/hpc/home/kailasamms/projects/scRNASeq/scripts"
  script_file <- "/hpc/home/kailasamms/projects/scRNASeq/Custom_Functions.R"
}

if (file.exists(script_file)) {
  source(script_file)
} else{
  stop(paste("Custom_Functions.R not found at:", script_file))
}

# Define these 4 variables in project specific script
# proj
# species
# contrasts <- c()
# deseq2.override <- list()
# heatmap.override <- list()
# volcano.override <- list()

proj.params <- setup_project(proj             = proj,
                             species          = species,  #"Mus musculus", "Homo sapiens"
                             contrasts        = contrasts,
                             parent_dir       = parent_dir,
                             gmt_dir          = gmt_dir,
                             scripts_dir      = scripts_dir,
                             deseq2.override  = deseq2.override,
                             heatmap.override = heatmap.override,
                             volcano.override = volcano.override)

proj_name  <- proj.params$proj
species    <- proj.params$species
proj_dir   <- proj.params$proj_dir
gmt_dir    <- proj.params$gmt_dir
counts_dir <- proj.params$counts_dir
salmon_dir <- proj.params$salmon_dir

# ---- 🔖 Fetch Annotations ONCE at the start ----

#ann_list <- get_annotations()

# ---- 🧾 Load Metadata ----

metadata_xlsx <- file.path(proj_dir, paste0(proj_name, "_Metadata.xlsx"))
metadata <- openxlsx::read.xlsx(xlsxFile = metadata_xlsx)

# ---- 🔢 Load raw counts ----

raw_counts_xlsx <- file.path(proj.params$proj_dir, paste0(proj, "_Raw_counts.xlsx"))
if (dir.exists(counts_dir) && length(list.files(counts_dir)) > 0){
  
  # Compile raw counts from  txt files
  raw_counts_mat <- merge_counts(counts_dir = counts_dir,
                                 filename   = proj_name, 
                                 output_dir = proj_dir)
  txi <- NULL
  
} else if (file.exists(raw_counts_xlsx)){
  
  # Read raw counts from excel file
  raw_counts_df <- openxlsx::read.xlsx(xlsxFile = raw_counts_xlsx)
  txi <- NULL
  
  if (!"SYMBOL" %in% colnames(raw_counts_df)) {
    log_warn(sample = "",
             step   = "RNASeq_Workflow.R",
             msg    = "1st column of xlsx file should be named 'SYMBOL'")
  } else {
    raw_counts_mat <- raw_counts_df %>% 
      tibble::column_to_rownames("SYMBOL")
  }
  
} else if (dir.exists(salmon_dir) && length(list.files(salmon_dir)) > 0){
  
  # Read counts from salmon quant files
  txi <- prep_txi(salmon_dir = salmon_dir, 
                  species    = species, 
                  db_version = "113", 
                  filename   = NULL, 
                  output_dir = proj_dir)
  raw_counts_mat <- NULL
  
} else {
  log_error(sample = "",
            step   = "RNASeq_Workflow.R",
            msg    = glue::glue("Please provide raw counts in excel file '{raw_counts_xlsx}' or as txt files in '{counts_dir}'"))
}

# ---- 🚧️ Prepare DESeq2 Input ----

# We use annotation mainly for plotting gene names and pathway analysis.
# If trasncript/gene has no gene name, we keep it for DESeq2 but remove before plotting, pathway analysis
ann_df  <- readr::read_csv("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/RNASeq_Human_Manish_22Rv1_ARCaPM/tx2gene_GRCh38.115.csv") %>%
  # Retain unique gene_id-> gene_name combinations
  select(gene_id, gene_name, gene_biotype) %>% 
  distinct() %>%
  # If protein-coding gene has gene_id but no gene_name, use gene_id instead of gene_name
  mutate(gene_name = ifelse((is.na(gene_name) | gene_name == "") & (gene_biotype == "protein_coding"), gene_id, gene_name))

add_annotation <- function(df) {
  
  # Ensure ID column is character for a clean join
  df <- df %>% dplyr::mutate(ID = as.character(ID))
  
  log_info(sample = "", step = "add_annotation", msg = "Joining with local GTF annotation...")
  
  # Join and apply your "Protein-Coding fallback" logic
  annotated_df <- df %>%
    dplyr::left_join(ann_df, by = c("ID" = "gene_id")) %>%
    # Use the pre-processed gene_name from your ann_df (which already has the fallbacks)
    # Then filter out the remaining NAs (unnamed non-coding)
    dplyr::filter(!is.na(gene_name) & gene_name != "") %>%
    dplyr::select(gene_name, dplyr::everything()) %>%
    dplyr::rename(SYMBOL = gene_name)
  
  return(annotated_df)
}

metadata <- read.xlsx("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/RNASeq_Human_Manish_22Rv1_ARCaPM/input/Manish_22Rv1_ARCaPM_metadata.xlsx")
metadata <- read.xlsx("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/RNASeq_Xenograft_Manish_22Rv1_Xeno/input/Manish_22RV1_Xeno_metadata.xlsx")
metadata <- read.xlsx("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/RNASeq_Xenograft_Manish_22Rv1_Xeno2/input/Manish_22RV1_Xeno2_metadata.xlsx")

txi <- readRDS("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/RNASeq_Human_Manish_22Rv1_ARCaPM/Human_full/03.Salmon/Salmon_txi_Human_full.rds")
txi <- readRDS("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/RNASeq_Xenograft_Manish_22Rv1_Xeno/Human_split/03.Salmon/Salmon_txi_Human_split.rds")
txi <- readRDS("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/RNASeq_Xenograft_Manish_22Rv1_Xeno2/Human_split/03.Salmon/Salmon_txi_Human_split.rds")

design      <- proj.params$deseq2$design
raw_counts_mat <- NULL
deseq2_data <- prepare_deseq2_input(expr_mat = raw_counts_mat,
                                    txi      = txi, 
                                    metadata = metadata,
                                    design   = design)

# ---- 🔍 QC : Initial PCA & Batch Inspection ----

plot_pca(expr_mat    = deseq2_data$expr_mat,
         txi         = deseq2_data$txi,
         metadata    = deseq2_data$metadata,
         filename    = proj_name,
         output_dir  = proj_dir,
         top_n_genes = 500,
         perform_vst = TRUE, 
         skip_plot   = FALSE)

# ---- 🧹 Sample Filtering : Outlier Removal ----

remove_samples <- NULL
#remove_samples <- c("MT2") # Manish_22Rv1_Xeno
# remove_samples <- c("SBQuadFc2", "SBQuadFc4")

if (!is.null(remove_samples)){
  
  raw_counts_mat <- raw_counts_mat[, !(colnames(raw_counts_mat) %in% remove_samples), drop = FALSE]
  
  txi$counts    <- txi$counts[, colnames(txi$counts) != remove_samples]
  txi$abundance <- txi$abundance[, colnames(txi$abundance) != remove_samples]
  txi$length    <- txi$length[, colnames(txi$length) != remove_samples]
  
  # ---- 🚧️ Prepare DESeq2 Input ----
  
  design      <- proj.params$deseq2$design
  deseq2_data <- prepare_deseq2_input(expr_mat = raw_counts_mat,
                                      txi      = txi, 
                                      metadata = metadata,
                                      design   = design)
}

# ---- 🧬 Differential Expression (DESeq2 model) ----

dds <- fit_deseq2_model(expr_mat    = deseq2_data$expr_mat, 
                        txi         = deseq2_data$txi, 
                        metadata    = deseq2_data$metadata,
                        design      = design)

# Extract non-blind VST Counts (for all downstream visualization)
vsd <- DESeq2::vst(dds, blind = FALSE)
vst_counts <- SummarizedExperiment::assay(vsd) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ID") %>%
  add_annotation() %>%
  # Keep ONLY copy with highest expresion
  mutate(row_mean = rowMeans(select(., -ID, -SYMBOL, -gene_biotype))) %>%
  group_by(SYMBOL) %>%
  slice_max(order_by = row_mean, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-ID, -gene_biotype, -row_mean) %>%
  tibble::column_to_rownames("SYMBOL") %>%
  as.matrix()
  
  # add_annotation(ann_list = ann_list) %>%
  # dplyr::group_by(SYMBOL) %>%
  # dplyr::summarize(across(.cols = where(is.numeric), .fns = mean, na.rm = TRUE), .groups = "drop") %>%
  # tibble::column_to_rownames("SYMBOL") %>%
  # as.matrix()

# Extract blind VST Counts (for global heatmap ONLY)
# NOTE: vst counts are affected by design ONLY when blind=FALSE
vsd_blind <- DESeq2::vst(object = dds, blind  = TRUE)
vst_counts_blind <- SummarizedExperiment::assay(vsd_blind) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ID") %>%
  add_annotation() %>%
  # Keep ONLY copy with highest expresion
  mutate(row_mean = rowMeans(select(., -ID, -SYMBOL, -gene_biotype))) %>%
  group_by(SYMBOL) %>%
  slice_max(order_by = row_mean, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-ID, -gene_biotype, -row_mean) %>%
  tibble::column_to_rownames("SYMBOL") %>%
  as.matrix()

  # add_annotation(ann_list = ann_list) %>%
  # dplyr::group_by(SYMBOL) %>%
  # dplyr::summarize(across(.cols = where(is.numeric), .fns = mean, na.rm = TRUE), .groups = "drop") %>%
  # tibble::column_to_rownames("SYMBOL") %>%
  # as.matrix()

# Extract normalized Counts (for violin/box plot of gene expression ONLY)
norm_counts <- DESeq2::counts(dds, normalized = TRUE) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ID") %>%
  add_annotation() %>%
  # Keep ONLY copy with highest expresion
  mutate(row_mean = rowMeans(select(., -ID, -SYMBOL, -gene_biotype))) %>%
  group_by(SYMBOL) %>%
  slice_max(order_by = row_mean, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-ID, -gene_biotype, -row_mean) %>%
  tibble::column_to_rownames("SYMBOL") %>%
  as.matrix()

  # add_annotation(ann_list = ann_list) %>%
  # dplyr::group_by(SYMBOL) %>%
  # dplyr::summarize(across(.cols = where(is.numeric), .fns = mean, na.rm = TRUE), .groups = "drop") %>%
  # tibble::column_to_rownames("SYMBOL") %>%
  # as.matrix()

# Create a new workbook
wb <- createWorkbook()
addWorksheet(wb, "VST_counts_non_blind")
writeData(wb, sheet = "VST_counts_non_blind", vst_counts, rowNames = TRUE)
saveWorkbook(wb, file = file.path(proj_dir, "VST_counts_non_blind.xlsx"), overwrite = TRUE)

wb <- createWorkbook()
addWorksheet(wb, "VST_counts_blind")
writeData(wb, sheet = "VST_counts_blind", vst_counts_blind, rowNames = TRUE)
saveWorkbook(wb, file = file.path(proj_dir, "VST_counts_blind.xlsx"), overwrite = TRUE)

wb <- createWorkbook()
addWorksheet(wb, "Normalized_counts")
writeData(wb, sheet = "Normalized_counts", norm_counts, rowNames = TRUE)
saveWorkbook(wb, file = file.path(proj_dir, "Normalized_counts.xlsx"), overwrite = TRUE)


# ==== 🌐 GLOBAL QUALITY & EXPLORATION ====

# ---- 📉 QC : Dispersion Estimates ----

plot_dispersion(dds        = dds,
                filename   = proj_name,
                output_dir = proj_dir)

# ---- 🔥 Visualization : Heatmap (Global) ----

# Perform LRT
# 'reduced = ~1' tests if the 'design(dds)' significantly improves the model over no groups at all
dds_LRT <- DESeq2::DESeq(dds, test = "LRT", reduced = ~1)
res_LRT <- DESeq2::results(dds_LRT) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ID") %>%
  add_annotation(ann_list = ann_list)

# Identify significant DEGs
sig_genes <- res_LRT %>% 
  dplyr::filter(padj <= 0.05) %>% 
  dplyr::pull(SYMBOL) %>%
  base::intersect(rownames(vst_counts_blind))

# Plot heatmap
ph <- plot_heatmap(expr_mat            = vst_counts_blind[sig_genes, , drop = FALSE], 
                   label_genes         = NULL,
                   filename            = proj_name,
                   output_dir          = proj_dir,
                   metadata_col        = metadata, 
                   metadata_row        = NULL,
                   col_annotations     = proj.params$heatmap$col_annotations,
                   row_annotations     = proj.params$heatmap$row_annotations,
                   col_gap_by          = proj.params$heatmap$col_gap_by,
                   row_gap_by          = proj.params$heatmap$row_gap_by,
                   col_cluster_by      = proj.params$heatmap$col_cluster_by,
                   row_cluster_by      = proj.params$heatmap$row_cluster_by,
                   plot_title          = proj.params$heatmap$plot_title,
                   heatmap_palette     = proj.params$heatmap$heatmap_palette,
                   annotation_palette  = proj.params$heatmap$annotation_palette,
                   border_color        = proj.params$heatmap$border_color,
                   force_log           = proj.params$heatmap$force_log,
                   show_expr_legend    = proj.params$heatmap$show_expr_legend,
                   save_plot           = proj.params$heatmap$save_plot,
                   save_matrix         = proj.params$heatmap$save_matrix)



# ==== 🔁 DIFFERENTIAL EXPRESSION PER CONTRAST ====

contrasts <- proj.params$deseq2$contrasts

for (contrast in contrasts) {
  
  log_info(sample = contrast,
           step   = "RNASeq_Workflow.R",
           msg    = glue::glue("Extracting DESeq2 results for contrast: '{contrast}'"))
  
  # Identify samples for current contrast
  contrast_samples <- metadata %>%
    dplyr::filter(Comparisons %in% all.vars(expr = as.formula(paste0("~", contrast)))) %>%
    dplyr::pull(Sample_ID) %>%
    as.character()
  
  # ---- 🧬 Differential Expression (DESeq2 results) ----
  
  valid_contrast <- contrast |> 
    gsub(pattern = "/", replacement = "_over_") |> 
    gsub(pattern = "//+", replacement = "_plus_") |> 
    #gsub(pattern = "-", replacement = "_vs_") |> 
    gsub(pattern = "[[:space:][:punct:]]+", replacement = "_") |> 
    gsub(pattern = "(^_+|_+$)", replacement = "")
  contrast_dir <- file.path(proj_dir, valid_contrast)
  
  # deseq2_results <- get_deseq2_results(dds         = dds,
  #                                      contrast    = contrast,
  #                                      output_dir  = contrast_dir,
  #                                      lfc_cutoff  = 0, 
  #                                      padj_cutoff = 0.1)
  
  extract_deseq2_results(contrast=contrast, dds=dds, tx2gene_csv_path=NULL, output_dir=NULL, lfc_cutoff = 0, padj_cutoff = 0.1)}
  
  # ---- 🏹 Visualization : MA Plot ----
  
  plot_ma(dds        = deseq2_results$dds,  # results() MUST be called on dds earlier
          filename   = paste(proj_name, contrast, sep = "_"), 
          output_dir = contrast_dir)
  
  # ---- 🌋 Visualization : Volcano Plot ----
  
  plot_volcano(res_df      = deseq2_results$degs, 
               filename    = paste(proj_name, contrast, sep = "_"), 
               output_dir  = contrast_dir, 
               contrast    = contrast,
               label_genes = proj.params$volcano$label_genes,
               top_n       = proj.params$volcano$top_n,
               lfc_cutoff  = proj.params$volcano$lfc_cutoff, 
               padj_cutoff = proj.params$volcano$padj_cutoff)
  
  # ---- 🔥 Visualization : Heatmap ----
  
  # Keep only samples present in vst
  common_samples <- base::intersect(contrast_samples, colnames(vst_counts))

  # Identify significant DEGs
  sig_genes <- deseq2_results$degs %>% 
    dplyr::filter(padj <= 0.05) %>% 
    dplyr::pull(SYMBOL) %>%
    base::intersect(rownames(vst_counts))
  
  # Plot heatmap
  ph <- plot_heatmap(expr_mat            = vst_counts[sig_genes, common_samples, drop = FALSE], 
                     label_genes         = NULL,
                     filename            = paste(contrast, proj_name, sep = "_"),
                     output_dir          = contrast_dir,
                     metadata_col        = metadata, 
                     metadata_row        = NULL,
                     col_annotations     = proj.params$heatmap$col_annotations,
                     row_annotations     = proj.params$heatmap$row_annotations,
                     col_gap_by          = proj.params$heatmap$col_gap_by,
                     row_gap_by          = proj.params$heatmap$row_gap_by,
                     col_cluster_by      = proj.params$heatmap$col_cluster_by,
                     row_cluster_by      = proj.params$heatmap$row_cluster_by,
                     plot_title          = proj.params$heatmap$plot_title,
                     heatmap_palette     = proj.params$heatmap$heatmap_palette,
                     annotation_palette  = proj.params$heatmap$annotation_palette,
                     border_color        = proj.params$heatmap$border_color,
                     force_log           = proj.params$heatmap$force_log,
                     show_expr_legend    = proj.params$heatmap$show_expr_legend,
                     save_plot           = proj.params$heatmap$save_plot,
                     save_matrix         = proj.params$heatmap$save_matrix)
  
  # ---- 🛤️ Pathway Analysis (GSEA & ORA) ----
  
  pathway_dir <- file.path(proj_dir, contrast, "Pathway_Analysis")
  pathway_results <- analyze_pathway(res_df     = deseq2_results$degs,
                                     species    = species, 
                                     gmt_dir    = gmt_dir,
                                     output_dir = pathway_dir,
                                     minsize    = 15, 
                                     maxsize    = 500)
  
  # ---- 🌿 Visualization : Pathway Enrichment ----
  
  # Identify top 10 Up & top 10 Down GSEA pathways for each collection
  top_gsea <- pathway_results$consensus %>%
    dplyr::filter(method != "ORA") %>%
    # Deduplicate by method priority
    dplyr::group_by(Collection, Consensus, Description) %>%
    dplyr::slice_min(order_by = match(method, c("FGSEA", "GSEA", "ORA")), n = 1) %>%
    dplyr::ungroup() %>%
    # Rank based on abs(NES) for each direction
    dplyr::group_by(Collection, Consensus) %>%
    dplyr::slice_max(order_by = abs(NES), n = 10, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  # Identify top 10 Up & top 10 Down ORA pathways for each collection
  top_ora <- pathway_results$consensus %>%
    dplyr::filter(method == "ORA") %>%
    # Rank based on padj for each direction
    dplyr::group_by(Collection, Consensus) %>%
    dplyr::slice_min(order_by = padj, n = 10, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  # Keep only samples present in vst
  common_samples <- base::intersect(contrast_samples, colnames(vst_counts))
  
  if (nrow(top_gsea) > 0) {
    plot_pathway(pathway_df = top_gsea, 
                 method     = "GSEA",
                 expr_mat   = vst_counts[, common_samples], 
                 metadata   = metadata,
                 output_dir = pathway_dir)
  }
  if (nrow(top_ora) > 0) {
    plot_pathway(pathway_df = top_ora, 
                 method     = "ORA",
                 expr_mat   = vst_counts[, common_samples], 
                 metadata   = metadata,
                 output_dir = pathway_dir)
  }
  

  # ---- 📡 Regulatory Network Analysis ----
  
  # Keep only samples present in vst
  common_samples <- base::intersect(contrast_samples, colnames(vst_counts))
  
  tf_dir <- file.path(proj_dir, contrast, "TF_Analysis")
  tf_activity_samples <- analyze_tf(expr_mat   = vst_counts[, common_samples], 
                                    res_df     = NULL, 
                                    species    = species, 
                                    output_dir = tf_dir,
                                    stats      = c("ulm", "mlm", "viper"), 
                                    minsize    = 5, 
                                    top_n      = 500) 
  
  tf_activity_degs <- analyze_tf(expr_mat   = NULL, 
                                 res_df     = deseq2_results$degs, 
                                 species    = species, 
                                 output_dir = tf_dir,
                                 stats      = c("ulm", "mlm", "viper"), 
                                 minsize    = 5, 
                                 top_n      = 500)
  
  # ---- 🌿 Visualization : Regulatory Networks ----
   
  plot_tf(tf_df      = tf_activity_samples$tf, 
          contrast   = contrast, 
          metadata   = metadata, 
          output_dir = tf_dir,
          top_n      = 20)
  
  plot_tf(tf_df      = tf_activity_degs$tf, 
          contrast   = contrast, 
          metadata   = metadata, 
          output_dir = tf_dir,
          top_n      = 20)
}

# ---- THE END ----

path1 <- NULL
#path1 <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/RNASeq_Xenograft_Manish_22Rv1_Xeno"
path2 <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/RNASeq_Xenograft_Manish_22Rv1_Xeno2"
path3 <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/RNASeq_Human_Manish_22Rv1_ARCaPM"

# Find all Pathway_results.xlsx files in the 09.Pathway directory
pathway_files <- list.files(
  path       = c(path1, path2, path3),
  pattern    = "^Pathway_results//.xlsx$", # Matches exactly this name
  recursive  = TRUE, 
  full.names = TRUE
)

tf_files <- list.files(
  path       = c(path1, path2, path3),
  pattern    = "^TF_results//.xlsx$", # Matches exactly this name
  recursive  = TRUE, 
  full.names = TRUE
)


# 2. Read pathway files and combine
pathway_data <- purrr::map_df(pathway_files, function(f) {
  
  # Get Names
  comp_name <- basename(dirname(f))
  proj_name <- basename(dirname(dirname(dirname(dirname(f)))))
  full_id   <- paste(proj_name, comp_name, sep = " | ")
  
  # Read Data
  df <- openxlsx::read.xlsx(f, sheet = "consensus") %>%
    dplyr::filter(padj < 0.05, !is.na(NES)) %>%
    dplyr::select(Collection, Description, NES) %>%
    dplyr::distinct(Collection, Description, .keep_all = TRUE)
  
  # PRINT THE NROW TO CONSOLE
  message(paste0("ID: ", full_id, " | Rows: ", nrow(df)))
  
  # Add the ID column
  df %>% dplyr::mutate(Comparison = full_id)
})

# 2. Read TF files and combine
tf_data <- purrr::map_df(tf_files, function(f) {
  
  # Get Names
  comp_name <- basename(dirname(f))
  proj_name <- basename(dirname(dirname(dirname(dirname(f)))))
  full_id   <- paste(proj_name, comp_name, sep = " | ")
  
  # Read Data
  df <- openxlsx::read.xlsx(f, sheet = "tf") %>%
    #dplyr::filter(p_value < 0.05, statistic == "mlm") %>%
    #dplyr::filter(p_value < 0.05, statistic == "ulm") %>%
    #dplyr::filter(p_value < 0.05, statistic == "viper") %>%
    dplyr::filter(p_value < 0.05, statistic == "consensus") %>%
    dplyr::select(source, score)

  # summary_df <- tf_results %>%
  #   # 1. Focus only on the core methods
  #   filter(statistic %in% c("ulm", "mlm", "viper")) %>%
  #   group_by(statistic) %>%
  #   # 2. Standardize scores (Z-score) within each method to equalize their 'vote'
  #   mutate(std_score = as.vector(scale(score))) %>% 
  #   group_by(source) %>%
  #   summarize(
  #     methods_count = n_distinct(statistic[p_value < 0.05]),
  #     # 3. Simple mean of standardized scores across methods
  #     avg_std_score = mean(std_score, na.rm = TRUE),
  #     # 4. Final Weighted Score: Average Score * Method Agreement
  #     weighted_consensus = avg_std_score * (methods_count / 3),
  #     methods_list = paste(unique(statistic[p_value < 0.05]), collapse = ", "),
  #     consistent_direction = all(score > 0) | all(score < 0)
  #   ) %>%
  #   # 5. Filter for TFs found by at least 1 method and sort by the new Consensus
  #   filter(methods_count > 0) %>%
  #   arrange(desc(abs(weighted_consensus)))
  
  # PRINT THE NROW TO CONSOLE
  message(paste0("ID: ", full_id, " | Rows: ", nrow(df)))
  
  # Add the ID column
  df %>% dplyr::mutate(Comparison = full_id)
})

# Create an NES Matrix for Heatmaps
pathway_comparison_nes <- pathway_data %>%
  #dplyr::filter(!grepl("0Gy_ARCaPM|4Gy_ARCaPM|PRDX1", Comparison)) %>%
  dplyr::filter(grepl("PRDX1_4Gy_22Rv1|NDRG1_4Gy_22Rv1|SPB_IRR", Comparison)) %>%
  dplyr::filter(!grepl("0Gy|Vehicle", Comparison)) %>%
   tidyr::pivot_wider(names_from = Comparison, values_from = NES) %>%
  dplyr::mutate(across(where(is.numeric), ~tidyr::replace_na(., 0)))

# Create an TF activity Matrix for Heatmaps
tf_comparison_score <- tf_data %>%
  #dplyr::filter(!grepl("0Gy_ARCaPM|4Gy_ARCaPM|PRDX1", Comparison)) %>%
  dplyr::filter(grepl("PRDX1_4Gy_22Rv1|NDRG1_4Gy_22Rv1|SPB_IRR", Comparison)) %>%
  dplyr::filter(!grepl("0Gy|Vehicle", Comparison)) %>%
  tidyr::pivot_wider(names_from = Comparison, values_from = score) %>%
  dplyr::mutate(across(where(is.numeric), ~tidyr::replace_na(., 0)))

# 1. Prepare the Annotation Data Frame first
# This maps every pathway to its Collection
row_ann <- pathway_comparison_nes %>%
  mutate(RowID = make.names(Description, unique = TRUE)) %>%
  select(RowID, Collection) %>%
  dplyr::rename(SYMBOL = RowID)

# row_ann <- tf_comparison_score %>%
#   mutate(RowID = make.names(source, unique = TRUE)) %>%
#   select(RowID) %>%
#   dplyr::rename(SYMBOL = RowID)

col_ann <- data.frame(Sample_ID = colnames(pathway_comparison_nes)) %>%
  dplyr::filter(!Sample_ID %in% c("Collection", "Description")) %>%
  dplyr::mutate(ID = Sample_ID, 
                Sample_ID = gsub(pattern="^.*//| ", "", Sample_ID),
                Source = case_when( grepl("Xeno", ID) ~ "Tumor",
                                    grepl("ARCaPM", Sample_ID) ~ "ARCaPM",
                                    grepl("22Rv1", Sample_ID) ~ "22Rv1",
                                    TRUE ~ "CellLine"),
                Treatment = case_when(grepl("NDRG1", Sample_ID) ~ "NDRG1",
                                      grepl("PRDX1", Sample_ID) ~ "PRDX1",
                                      grepl("SPB", Sample_ID) ~ "SPB",
                                      grepl("IRR", Sample_ID) ~ "IRR",
                                      TRUE ~ "WT"),
                Gray = case_when(grepl("4Gy|IRR", Sample_ID) ~ "4Gy",
                                 TRUE ~ "0Gy"),
                )

col_ann <- data.frame(Sample_ID = colnames(tf_comparison_score)) %>%
  dplyr::filter(!Sample_ID %in% c("Collection", "Description")) %>%
  dplyr::mutate(ID = Sample_ID, 
                Sample_ID = gsub(pattern="^.*//| ", "", Sample_ID),
                Source = case_when( grepl("Xeno", ID) ~ "Tumor",
                                    grepl("ARCaPM", Sample_ID) ~ "ARCaPM",
                                    grepl("22Rv1", Sample_ID) ~ "22Rv1",
                                    TRUE ~ "CellLine"),
                Treatment = case_when(grepl("NDRG1", Sample_ID) ~ "NDRG1",
                                      grepl("PRDX1", Sample_ID) ~ "PRDX1",
                                      grepl("SPB", Sample_ID) ~ "SPB",
                                      grepl("IRR", Sample_ID) ~ "IRR",
                                      TRUE ~ "WT"),
                Gray = case_when(grepl("4Gy|IRR", Sample_ID) ~ "4Gy",
                                 TRUE ~ "0Gy"),
  )
                
# Rename after creating col_ann                                
colnames(pathway_comparison_nes) <- gsub(pattern="^.*//| ", "", colnames(pathway_comparison_nes))                                
colnames(tf_comparison_score) <- gsub(pattern="^.*//| ", "", colnames(tf_comparison_score)) 

# 2. Prepare the Numeric Matrix
# We use the same RowID logic to ensure they match
pathway_mat <- pathway_comparison_nes %>%
  mutate(RowID = make.names(Description, unique = TRUE)) %>%
  tibble::column_to_rownames("RowID") %>%
  select(-Collection, -Description) %>%
  as.matrix()

tf_mat <- tf_comparison_score %>%
  mutate(RowID = make.names(source, unique = TRUE)) %>%
  tibble::column_to_rownames("RowID") %>%
  select(-source) %>%
  as.matrix()

# Plot heatmap
ph <- plot_heatmap(expr_mat            = pathway_mat, 
                   label_genes         = NULL,
                   #filename            = proj_name,
                   output_dir          = "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/",
                   metadata_col        = col_ann, 
                   metadata_row        = row_ann,
                   col_annotations     = c("Source", "Treatment", "Gray"),
                   row_annotations     = c("Collection"),
                   #col_gap_by          = c("Gray"),
                   row_gap_by          = c("Collection"),
                   col_cluster_by      = "all", #c("Gray"),
                   row_cluster_by      = c("Collection"),
                   plot_title          = "Pathway NES Comparison",
                   #heatmap_palette     = "rdbu",
                   #annotation_palette  = proj.params$heatmap$annotation_palette,
                   #border_color        = proj.params$heatmap$border_color,
                   #force_log           = proj.params$heatmap$force_log,
                   #show_expr_legend    = proj.params$heatmap$show_expr_legend,
                   save_plot           = TRUE,
                   save_matrix         = TRUE)

ph <- plot_heatmap(expr_mat            = tf_mat, 
                   label_genes         = NULL,
                   #filename            = proj_name,
                   output_dir          = "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/",
                   metadata_col        = col_ann, 
                   #metadata_row        = row_ann,
                   col_annotations     = c("Source", "Treatment", "Gray"),
                   #row_annotations     = c("Collection"),
                   #col_gap_by          = c("Gray"),
                   #row_gap_by          = c("Collection"),
                   col_cluster_by      = "all",
                   #row_cluster_by      = c("Collection"),
                   plot_title          = "TF Activity Comparison",
                   #heatmap_palette     = "rdbu",
                   #annotation_palette  = proj.params$heatmap$annotation_palette,
                   #border_color        = proj.params$heatmap$border_color,
                   #force_log           = proj.params$heatmap$force_log,
                   #show_expr_legend    = proj.params$heatmap$show_expr_legend,
                   save_plot           = TRUE,
                   save_matrix         = TRUE)

# Find common up /down from ulm, mlm, consensus, viper
base_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/RNASeq_Xenograft+CellLine_Manish_Comparisons"
tf_filenames <- list.files(path = base_path, pattern = "^TF.*.xlsx$")
full_paths <- file.path(base_path, tf_filenames)

# 4. Read all files into a list, then combine into one data frame
# setNames helps keep track of which data came from which file
df_list <- lapply(setNames(full_paths, tf_filenames), read.xlsx)

# Combine them into one big data frame
# .id = "source_file" adds a column telling you which file the data came from
df <- bind_rows(df_list, .id = "source_file")
colnames(df)[2] <- "TF"

# Identify which columns are numeric (excluding your source/ID columns)
num_cols <- sapply(df, is.numeric)

# Filter: (All values > 0) OR (All values < 0)
df_filtered <- df[rowSums(df[, num_cols] > 0) == sum(num_cols) | 
                    rowSums(df[, num_cols] < 0) == sum(num_cols), ] %>%
  dplyr::add_count(TF)
write.xlsx(df_filtered, file = file.path(base_path, "COmmon_TFs.xlsx"), overwrite = TRUE)



