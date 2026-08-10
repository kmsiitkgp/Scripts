#******************************************************************************#
#                         SASWAT PROJECT - ANALYSIS PIPELINE                   #
#  Order: Data prep → Normalize → Classify → DEGs → Hits → Survival → ORA    #
#******************************************************************************#

# =============================================================================#
# PATHS (edit once here, used everywhere)
# =============================================================================#
file_path   <- "/home/kailasamms/analysis/ICB_TCGA_Datasets/original/"
data_path   <- "/home/kailasamms/analysis/ICB_TCGA_Datasets/processed/"
output_path <- "/home/kailasamms/analysis/ICB_TCGA_Datasets/results/"
supp_path   <- "/home/kailasamms/analysis/ICB_TCGA_Datasets/supplemental/"

source("/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

#******************************************************************************#
#             SECTION 1: REFORMAT ICB CLINICAL/METADATA                        #
#******************************************************************************#

gse <- c("GSE78220", "GSE91061", "GSE96619", "GSE115821", "GSE131521")

for (proj in gse){

  meta_data     <- read.xlsx(paste0(file_path, proj, ".Metadata.xlsx"))
  clinical_data <- read.xlsx(paste0(file_path, "ICB_Clinical_Data.xlsx"))

  clinical_data <- clinical_data %>%
    dplyr::filter(series_id %in% proj) %>%
    dplyr::mutate(Time   = as.numeric(OS) / 30,
                  Status = OS_FLAG,
                  Sample = toupper(Sample))

  meta_data <- meta_data %>%
    dplyr::rename(Sample_ID = GEO_Accession__exp_) %>%
    dplyr::mutate(Description = toupper(Description),
                  Run         = toupper(Run))

  if (proj == "GSE78220"){
    metadata <- meta_data %>% dplyr::inner_join(clinical_data, by = c("Description" = "Sample"))
  } else {
    metadata <- meta_data %>% dplyr::inner_join(clinical_data, by = c("Run" = "Sample"))
  }

  cat("Missing from meta_data   :", nrow(metadata) - nrow(meta_data),    "\n")
  cat("Missing from clinical_data:", nrow(metadata) - nrow(clinical_data), "\n")

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "Metadata")
  openxlsx::writeData(wb, sheet = "Metadata", x = metadata)
  openxlsx::saveWorkbook(wb, file = paste0(data_path, make.names(proj), ".Metadata.xlsx"),
                         overwrite = TRUE)
}

#******************************************************************************#
#             SECTION 2: CALCULATE DESEQ2 NORMALIZED COUNTS                    #
#******************************************************************************#

gse        <- c("GSE78220", "GSE91061", "GSE96619", "GSE115821", "GSE131521",
                "IMvigor010", "IMvigor210")
species    <- "Homo sapiens"
annotations <- get_annotations(species)

for (proj in gse){

  meta_data <- read.xlsx(paste0(data_path, proj, ".Metadata.xlsx"))
  read_data <- openxlsx::read.xlsx(xlsxFile = paste0(data_path, proj, ".raw.counts.xlsx"))

  if (proj %in% c("GSE78220", "GSE91061", "GSE96619", "GSE115821", "GSE131521")){
    read_data <- read_data %>% dplyr::rename(SYMBOL = GeneID)
  }

  Comparisons <- list(Variable  = c(NA),
                      Target    = c(NA),
                      Reference = c(NA))

  meta_data <- prep_metadata(meta_data, read_data)
  read_data <- prep_readdata(read_data, meta_data)
  l         <- check_data(read_data, meta_data)
  meta_data <- l[[2]]
  read_data <- l[[1]]

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = read_data,
                                        colData   = meta_data,
                                        design    = ~ 1)

  norm_counts <- deseq2_norm_counts(dds, meta_data, annotations)

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "Normalized counts")
  openxlsx::writeData(wb, sheet = "Normalized counts", x = norm_counts)
  openxlsx::saveWorkbook(wb, file = paste0(data_path, make.names(proj), ".Normalized.counts.xlsx"),
                         overwrite = TRUE)
}

#******************************************************************************#
#         SECTION 3: CLASSIFY PATIENTS — IMMUNE ENRICHED vs IMMUNE DEPLETED    #
#******************************************************************************#

# TME pathway groups (from https://doi.org/10.1016/j.ccell.2021.04.014)
angio <- c("Angiogenesis", "Endothelium", "Cancer.associated.fibroblasts",
           "Matrix", "Matrix.remodeling")

pro_tumor_immune_infiltrate <- c(
  "Tumor.associated.Macrophages", "Macrophage.and.DC.traffic",
  "Myeloid.cells.traffic", "Immune.Suppression.by.Myeloid.Cells",
  "Th2.signature", "Treg.and.Th2.traffic", "Treg", "M1.signature",
  "Protumor.cytokines", "Neutrophil.signature", "Granulocyte.traffic")

anti_tumor_immune_infiltrate <- c(
  "MHCII", "Antitumor.cytokines", "Co.activation.molecules",
  "B.cells", "NK.cells", "Checkpoint.molecules", "Effector.cells",
  "T.cells", "Th1.signature", "Effector.cell.traffic", "MHCI")

malignant_cell <- c("EMT.signature", "Tumor.proliferation.rate")

# NOTE: HLA-A etc appear as HLA.A in norm counts — fix with gsub below
sig <- read.xlsx(paste0(file_path, "Saswat_Sig_list.xlsx"))
sig <- sig[-1, -1] %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames(colnames(.)[1]) %>%
  t() %>%
  data.frame() %>%
  dplyr::select(all_of(pro_tumor_immune_infiltrate),
                all_of(anti_tumor_immune_infiltrate)) %>%
  dplyr::mutate(across(.cols = everything(),
                       ~ gsub(pattern = "-", replacement = ".", x = .)))

meta_data_pancancer <- read.xlsx(paste0(data_path, "TCGA.PanAtlas.Metadata.xlsx"))
proj_list <- make.names(unique(meta_data_pancancer$Project_ID))
proj_list <- c(proj_list, "IMvigor010", "IMvigor210",
               "GSE78220", "GSE91061", "GSE96619", "GSE115821", "GSE131521",
               "TCGA.PanAtlas")

for (f in proj_list){

  meta_data <- read.xlsx(paste0(data_path, f, ".Metadata.xlsx"))

  if (f == "TCGA.PanAtlas"){
    read_data <- read.table(paste0(data_path, f, ".Normalized.counts.tsv"), header = TRUE)
    read_data <- read_data[, -c(2:3)]
  } else if (f %in% c("IMvigor010", "IMvigor210")){
    read_data <- read.xlsx(paste0(data_path, f, ".Normalized.counts.xlsx"))
    read_data <- read_data[, -c(2:4)]
  } else {
    read_data <- read.xlsx(paste0(data_path, f, ".Normalized.counts.xlsx"))
    read_data <- read_data[, -c(2:3)]
  }
  colnames(read_data) <- make.names(colnames(read_data))

  # Reformat: deduplicate, remove zero rows, keep only samples in metadata
  norm_counts <- read_data %>%
    base::replace(is.na(.), 0) %>%
    dplyr::mutate(n = rowSums(.[, -1])) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::slice_max(n) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n != 0) %>%
    dplyr::select(-n) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("SYMBOL") %>%
    dplyr::mutate(across(.cols = everything(), .fns = as.numeric)) %>%
    dplyr::select(all_of(intersect(make.names(meta_data$Sample_ID),
                                   make.names(colnames(read_data)))))

  colnames(norm_counts) <- base::make.names(colnames(norm_counts))

  log_norm_counts <- log(1 + norm_counts, base = 2)
  t_med <- base::apply(X = log_norm_counts, MARGIN = 1, FUN = median, na.rm = TRUE)
  log_norm_counts <- base::sweep(x = log_norm_counts, MARGIN = 1, FUN = "-", STATS = t_med)

  if (f == "TCGA.PanAtlas"){
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheetName = "median centered log normalized")
    openxlsx::writeData(wb, sheet = "median centered log normalized",
                        x = log_norm_counts, rowNames = TRUE)
    openxlsx::saveWorkbook(wb,
                           file = paste0(output_path, "TCGA.PanAtlas.median_centered_log_norm_counts.xlsx"),
                           overwrite = TRUE)
  }

  # Build gene set list
  gs <- list()
  for (i in 1:ncol(sig)){
    plot_genes_i <- list(sig[, i][!is.na(sig[, i])])
    names(plot_genes_i) <- colnames(sig)[i]
    gs <- c(gs, plot_genes_i)
  }

  # GSVA scores
  gsvaPar    <- GSVA::gsvaParam(exprData = as.matrix(log_norm_counts), geneSets = as.list(gs))
  gsva.scores <- GSVA::gsva(gsvaPar, verbose = TRUE)

  # ssGSEA scores
  ssgseaPar    <- GSVA::ssgseaParam(exprData = as.matrix(log_norm_counts), geneSets = as.list(gs))
  ssgsea.scores <- GSVA::gsva(ssgseaPar, verbose = TRUE)

  # Z-scores (Levine et al.)
  z.scores <- log_norm_counts %>%
    t() %>%
    data.frame() %>%
    dplyr::select(1) %>%          # dummy column to preserve df structure
    tibble::rownames_to_column("Sample")

  for (i in 1:ncol(sig)){
    plot_genes_i <- sig[, i][!is.na(sig[, i])]
    expr_df_i    <- as.data.frame(advanced_Z(plot_genes_i, log_norm_counts))
    colnames(expr_df_i) <- colnames(sig)[i]
    expr_df_i <- expr_df_i %>% tibble::rownames_to_column("Sample")
    z.scores  <- z.scores %>% dplyr::left_join(expr_df_i, by = "Sample")
  }

  z.scores <- z.scores[, -2] %>%
    tibble::column_to_rownames("Sample") %>%
    t()

  rownames(z.scores) <- make.names(rownames(z.scores))
  colnames(z.scores) <- make.names(colnames(z.scores))

  # ---- Heatmaps ----
  row_ann <- data.frame(
    Signature  = c(pro_tumor_immune_infiltrate, anti_tumor_immune_infiltrate),
    Row_Groups = c(rep("Pro-tumor Immune infiltrate",  11),
                   rep("Anti-tumor Immune infiltrate", 11)),
    Dummy_col  = NA) %>%
    tibble::column_to_rownames("Signature")

  ann_colors <- list(
    Column_Groups = c(`Immune Depleted`  = "#CB181D", `Immune Enriched` = "#A6D854"),
    Row_Groups    = c(`Pro-tumor Immune infiltrate` = "white",
                      `Anti-tumor Immune infiltrate` = "white"))

  for (s in c("gsva.scores", "ssgsea.scores", "z.scores")){

    scores <- get(s)
    mat    <- scores

    col_ann <- scores %>%
      t() %>%
      data.frame() %>%
      dplyr::mutate(Mean          = rowMeans(as.matrix(.)),
                    Column_Groups = dplyr::case_when(Mean > 0 ~ "Immune Enriched",
                                                     TRUE     ~ "Immune Depleted")) %>%
      dplyr::select(Column_Groups, Mean)

    col_order <- rownames(col_ann %>% dplyr::arrange(Mean))

    temp_mat  <- mat
    rowclust  <- hclust(dist(temp_mat))
    row_order <- rownames(temp_mat[rowclust$order, ])

    mat <- mat[row_order, col_order]

    if (max(mat) == 0){
      breaks     <- seq(from = min(mat), to = 0, length.out = 100)
    } else if (min(mat) == 0){
      breaks     <- seq(from = 0, to = max(mat), length.out = 100)
    } else if (min(mat) < -3 | max(mat) > 3){
      breaks     <- c(seq(-1.5, 0, length.out = 50), seq(1.5/100, 1.5, length.out = 50))
    } else {
      breaks     <- c(seq(from = min(mat), to = 0, length.out = 50),
                      seq(from = max(mat)/100, to = max(mat), length.out = 50))
    }

    rdbu <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100)

    gaps_col <- (col_ann %>%
                   dplyr::count(Column_Groups) %>%
                   dplyr::mutate(n = cumsum(n)) %>%
                   dplyr::pull(n))
    gaps_col <- gaps_col[gaps_col < ncol(mat)]

    row_mat <- t(scale(t(mat)))

    pheatmap::pheatmap(
      mat              = row_mat,
      scale            = "none",
      breaks           = breaks,
      color            = rdbu,
      border_color     = NA,
      fontsize         = 10,
      annotation_colors = ann_colors,
      cluster_rows     = FALSE,
      cluster_cols     = FALSE,
      gaps_col         = gaps_col,
      annotation_col   = col_ann %>% dplyr::select(Column_Groups),
      annotation_legend = TRUE,
      show_rownames    = FALSE,
      show_colnames    = FALSE,
      annotation_names_row = FALSE,
      annotation_names_col = FALSE,
      width            = 20,
      filename         = paste0(output_path, f, "_Heatmap_", s, ".tiff"))
  }

  # ---- Merge classifications & NPEPPS expression ----
  npepps_zscore <- log_norm_counts %>%
    tibble::rownames_to_column("SYMBOL") %>%
    dplyr::filter(SYMBOL == "NPEPPS") %>%
    tibble::column_to_rownames("SYMBOL") %>%
    t() %>% data.frame() %>%
    tibble::rownames_to_column("Sample") %>%
    dplyr::mutate(Sample = make.names(Sample))

  npepps <- norm_counts %>%
    tibble::rownames_to_column("SYMBOL") %>%
    dplyr::filter(SYMBOL == "NPEPPS") %>%
    tibble::column_to_rownames("SYMBOL") %>%
    t() %>% data.frame() %>%
    tibble::rownames_to_column("Sample") %>%
    dplyr::mutate(Sample = make.names(Sample))

  summary_df <- data.frame(Sample = colnames(z.scores))

  for (s in c("gsva.scores", "ssgsea.scores", "z.scores")){

    scores  <- get(s)
    col_ann <- scores %>%
      t() %>%
      data.frame() %>%
      dplyr::mutate(Mean  = rowMeans(as.matrix(.)),
                    Class = dplyr::case_when(Mean > 0 ~ "Immune Enriched",
                                             TRUE     ~ "Immune Depleted")) %>%
      tibble::rownames_to_column("Sample") %>%
      dplyr::select(Sample, Class, Mean, Effector.cells) %>%
      dplyr::mutate(Sample = make.names(Sample)) %>%
      dplyr::rename(
        !!rlang::sym(gsub(".scores", ".mean",           rlang::as_string(rlang::sym(s)))) := Mean,
        !!rlang::sym(gsub(".scores", ".classification", rlang::as_string(rlang::sym(s)))) := Class,
        !!rlang::sym(gsub(".scores", ".CTL",            rlang::as_string(rlang::sym(s)))) := Effector.cells)

    summary_df <- summary_df %>% dplyr::left_join(col_ann, by = "Sample")
  }

  summary_df <- summary_df %>%
    dplyr::mutate(
      Comments = dplyr::case_when(
        gsva.classification == ssgsea.classification &
          gsva.classification == z.classification ~ "OK",
        TRUE ~ "AMBIGUOUS CLASSIFICATION"),
      gsva.33 = dplyr::case_when(
        gsva.mean <= quantile(gsva.mean, 0.33) |
          gsva.mean >= quantile(gsva.mean, 0.66) ~ gsva.classification,
        TRUE ~ NA_character_),
      ssgsea.33 = dplyr::case_when(
        ssgsea.mean <= quantile(ssgsea.mean, 0.33) |
          ssgsea.mean >= quantile(ssgsea.mean, 0.66) ~ ssgsea.classification,
        TRUE ~ NA_character_),
      # BUG FIX: was comparing z.mean to quantile(z.scores,...) — wrong object
      z.33 = dplyr::case_when(
        z.mean <= quantile(z.mean, 0.33) |
          z.mean >= quantile(z.mean, 0.66) ~ z.classification,
        TRUE ~ NA_character_)) %>%
    dplyr::left_join(npepps_zscore, by = "Sample") %>%
    dplyr::left_join(npepps,        by = "Sample")

  cat(nrow(summary_df %>% dplyr::filter(Comments != "OK")), "out of",
      nrow(summary_df), "samples ambiguous in", f, "\n")

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "NPEPPS");        openxlsx::writeData(wb, "NPEPPS",        summary_df,         rowNames = FALSE)
  openxlsx::addWorksheet(wb, "gsva.scores");   openxlsx::writeData(wb, "gsva.scores",   t(gsva.scores),     rowNames = TRUE)
  openxlsx::addWorksheet(wb, "ssgsea.scores"); openxlsx::writeData(wb, "ssgsea.scores", t(ssgsea.scores),   rowNames = TRUE)
  openxlsx::addWorksheet(wb, "z.scores");      openxlsx::writeData(wb, "z.scores",      t(z.scores),        rowNames = TRUE)
  openxlsx::addWorksheet(wb, "Heatmap_matrix");openxlsx::writeData(wb, "Heatmap_matrix", row_mat,           rowNames = TRUE)
  openxlsx::saveWorkbook(wb,
                         file = paste0(output_path, f, ".Scores.&.Heatmap.details.xlsx"),
                         overwrite = TRUE)
}

#******************************************************************************#
#         SECTION 4: DEGs — IMMUNE DEPLETED vs IMMUNE ENRICHED                 #
#******************************************************************************#

species     <- "Homo sapiens"
annotations <- get_annotations(species)

meta_data_pancancer <- read.xlsx(paste0(data_path, "TCGA.PanAtlas.Metadata.xlsx"))
proj_list <- c(make.names(unique(meta_data_pancancer$Project_ID)), "IMvigor010", "IMvigor210")

for (f in proj_list){

  meta_data          <- read.xlsx(paste0(data_path, f, ".Metadata.xlsx"))
  classification_data <- read.xlsx(paste0(output_path, f, ".Scores.&.Heatmap.details.xlsx"),
                                   sheet = "NPEPPS")
  meta_data <- meta_data %>%
    dplyr::left_join(classification_data, by = c("Sample_ID" = "Sample"))

  if (f == "TCGA.PanAtlas"){
    read_data <- read.table(paste0(data_path, f, ".raw.counts.tsv"), header = TRUE)
    read_data <- read_data[, -c(2:3)]
  } else if (f %in% c("IMvigor010", "IMvigor210")){
    read_data <- read.xlsx(paste0(data_path, f, ".raw.counts.xlsx"))
    read_data <- read_data[, -c(2:4)]
  } else {
    read_data <- read.xlsx(paste0(data_path, f, ".raw.counts.xlsx"))
    read_data <- read_data[, -2]
  }
  colnames(read_data)[1] <- "SYMBOL"

  Comparisons <- list(
    Variable   = c("gsva.classification", "ssgsea.classification", "z.classification"),
    Target     = c("Immune Depleted",     "Immune Depleted",       "Immune Depleted"),
    Reference  = c("Immune Enriched",     "Immune Enriched",       "Immune Enriched"),
    lfc.cutoff  = 0,
    padj.cutoff = 0.1)

  meta_data <- prep_metadata(meta_data, read_data)
  read_data <- prep_readdata(read_data, meta_data)
  l         <- check_data(read_data, meta_data)
  meta_data <- l[[2]]
  read_data <- l[[1]]

  for (n in seq_along(Comparisons$Variable)){

    meta_data_comp <- meta_data %>%
      dplyr::mutate(id = get(Comparisons$Variable[n]))

    approach <- "DESeq2"
    if (length(unique(meta_data_comp$Batch)) > 1){
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = read_data,
                                            colData   = meta_data_comp,
                                            design    = ~ Batch + id)
    } else {
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = read_data,
                                            colData   = meta_data_comp,
                                            design    = ~ id)
    }

    prefix <- make.names(f)
    dds    <- run_deseq2(dds, meta_data_comp, annotations, Comparisons, n, approach, prefix, output_path)
  }
}

#******************************************************************************#
#       SECTION 5: DEG FREQUENCY ACROSS 33 CANCERS                             #
#******************************************************************************#

meta_data_pancancer <- read.xlsx(paste0(data_path, "TCGA.PanAtlas.Metadata.xlsx"))
proj_list <- make.names(unique(meta_data_pancancer$Project_ID))

wb <- openxlsx::createWorkbook()

for (s in c("gsva", "ssgsea", "z")){

  deg_genes <- c()
  for (f in proj_list){
    degs <- read.xlsx(paste0(output_path, f,
                             ".DEGs.", s, ".classification_Immune Depleted_vs_Immune Enriched_DESeq2_.xlsx")) %>%
      dplyr::filter(padj <= 0.05)
    deg_genes <- c(deg_genes, degs$SYMBOL)
  }

  all_degs <- unique(deg_genes)
  df       <- data.frame(SYMBOL = all_degs, DUMMY = 0)

  for (f in proj_list){
    degs <- read.xlsx(paste0(output_path, f,
                             ".DEGs.", s, ".classification_Immune Depleted_vs_Immune Enriched_DESeq2_.xlsx")) %>%
      dplyr::filter(padj <= 0.05) %>%
      dplyr::select(SYMBOL, log2FoldChange) %>%
      dplyr::rename(!!rlang::sym(f) := "log2FoldChange") %>%
      dplyr::distinct_at("SYMBOL", .keep_all = TRUE)

    df <- df %>% dplyr::left_join(degs, by = "SYMBOL")
  }

  df <- df %>%
    dplyr::select(-DUMMY) %>%
    tibble::column_to_rownames("SYMBOL")

  df_total <- df; df_up <- df; df_down <- df

  df_total[!is.na(df)] <- 1
  df_up[df > 0]        <- 1;  df_up[df <= 0]  <- 0
  df_down[df >= 0]     <- 0;  df_down[df < 0] <- 1

  df <- df %>%
    dplyr::mutate(n_TOTAL = rowSums(df_total, na.rm = TRUE),
                  n_UP    = rowSums(df_up,    na.rm = TRUE),
                  n_DOWN  = rowSums(df_down,  na.rm = TRUE)) %>%
    dplyr::select(n_TOTAL, n_UP, n_DOWN, everything()) %>%
    dplyr::arrange(desc(n_UP))

  openxlsx::addWorksheet(wb, sheetName = paste0("Summary_", s))
  openxlsx::writeData(wb, sheet = paste0("Summary_", s), x = df, rowNames = TRUE)
}
openxlsx::saveWorkbook(wb, file = paste0(output_path, "zSummary_of_DEGs.xlsx"), overwrite = TRUE)

# gsva classification: NPEPPS UP in 13 cancers, DOWN in 1 → use gsva going forward

#******************************************************************************#
#       SECTION 6: 2D UMAP OF PATHWAY SCORES ACROSS 33 CANCERS                 #
#******************************************************************************#

meta_data_pancancer <- read.xlsx(paste0(data_path, "TCGA.PanAtlas.Metadata.xlsx"))
proj_list <- make.names(unique(meta_data_pancancer$Project_ID))
wb <- openxlsx::createWorkbook()
s  <- "gsva"

# Seed empty matrix from first project then bind all
mat <- read.xlsx(paste0(output_path, proj_list[1], ".Scores.&.Heatmap.details.xlsx"),
                 sheet = paste0(s, ".scores"))
colnames(mat)[1] <- "Sample_ID"
mat <- mat[1, ]

for (f in proj_list){
  data <- read.xlsx(paste0(output_path, f, ".Scores.&.Heatmap.details.xlsx"),
                    sheet = paste0(s, ".scores"))
  colnames(data)[1] <- "Sample_ID"
  mat <- dplyr::bind_rows(mat, data)
}

mat <- mat %>%
  dplyr::distinct_at("Sample_ID", .keep_all = TRUE) %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("Sample_ID")

mat[is.na(mat)] <- 0
mat <- mat[rowSums(mat) != 0, ]

umap_res <- umap::umap(d = mat,
                       n_components  = 3,
                       config        = umap::umap.defaults,
                       method        = "naive",
                       preserve.seed = TRUE)

umap_val <- umap_res$layout %>%
  data.frame() %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::left_join(meta_data_pancancer %>%
                     dplyr::select(Sample_ID, Project_ID) %>%
                     dplyr::mutate(Sample_ID = make.names(Sample_ID)),
                   by = c("Sample" = "Sample_ID"))
colnames(umap_val) <- c("Sample", "UMAP1", "UMAP2", "UMAP3", "Project_ID")

pan_colors <- c("#000000","#D9D9D9","#003C30","#beb9db","#1D91C0","#A6CEE3",
                "#50E991","#A6D854","#74C476","#C7E9B4","#00bfa0","#E5F5F9",
                "#EDBF33","#E6D800","#FFF7BC","#ffee65","#C7EAE5",
                "#67001F","#CB181D","#FD8D3C","#FC9272","#EF3B2C","#F16913",
                "#9b19f5","#6A51A3","#762A83","#D4B9DA","#0bb4ff",
                "#E60049","#DC0AB4","#AE017E","#DF65B0","#FDCCE5")

ggplot2::ggplot(data = umap_val, aes(x = UMAP1, y = UMAP2, color = Project_ID, fill = Project_ID)) +
  geom_point(size = 1) +
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(shape = 22, size = 5, color = "black"))) +
  scale_color_manual(values = pan_colors) +
  scale_fill_manual(values  = pan_colors)

ggsave(paste0(output_path, "UMAP_", s, ".tiff"))

openxlsx::addWorksheet(wb, sheetName = paste0("UMAP_", s))
openxlsx::writeData(wb, sheet = paste0("UMAP_", s), x = umap_val, rowNames = FALSE)
openxlsx::saveWorkbook(wb, file = paste0(output_path, "UMAP_Coords.xlsx"), overwrite = TRUE)

#******************************************************************************#
#       SECTION 7: 3D UMAP — IMMUNE ENRICHED/DEPLETED ACROSS 33 CANCERS        #
#******************************************************************************#

meta_data_pancancer <- read.xlsx(paste0(data_path, "TCGA.PanAtlas.Metadata.xlsx"))
proj_list <- make.names(unique(meta_data_pancancer$Project_ID))
wb <- openxlsx::createWorkbook()

for (s in c("gsva", "ssgsea", "z")){

  group_data <- data.frame(Sample = c(""))

  for (f in proj_list){
    df <- read.xlsx(paste0(output_path, f, ".Scores.&.Heatmap.details.xlsx"),
                    sheet = "NPEPPS") %>%
      dplyr::select(Sample,
                    paste0(s, ".classification"),
                    paste0(s, ".mean"),
                    NPEPPS.x, NPEPPS.y)
    colnames(df) <- c("Sample", "Class", "Mean", "NPEPPS.x", "NPEPPS.y")
    group_data   <- dplyr::bind_rows(group_data, df)
  }
  group_data <- group_data[-1, ]

  openxlsx::addWorksheet(wb, sheetName = paste0(s, ".classification"))
  openxlsx::writeData(wb, sheet = paste0(s, ".classification"), x = group_data, rowNames = FALSE)
}
openxlsx::saveWorkbook(wb,
                       file = paste0(output_path, "Two_Group_Classification_from_Individual_cancers.xlsx"),
                       overwrite = TRUE)

# 3D scatter
s     <- "gsva"
mat   <- read.xlsx(paste0(output_path, "UMAP_Coords.xlsx"), sheet = paste0("UMAP_", s))
group <- read.xlsx(paste0(output_path, "Two_Group_Classification_from_Individual_cancers.xlsx"),
                   sheet = paste0(s, ".classification"))

mat <- mat %>%
  dplyr::left_join(group, by = "Sample") %>%
  dplyr::mutate(color = dplyr::case_when(Class == "Immune Depleted" ~ "#CB181D",
                                         TRUE ~ "#A6D854"))

jpeg(file = paste0(output_path, "UMAP_", s, "_3D.tiff"))
par(mar = c(1, 1, 1, 1))
plot3D::scatter3D(x = mat$UMAP1, y = mat$UMAP2, z = mat$UMAP3,
                  theta = 40, phi = 15, bty = "b2",
                  colvar = NULL, col = mat$color,
                  pch = 19, cex = 0.25,
                  colkey = list(mat$Class), clab = "Class")
dev.off()

#******************************************************************************#
#     SECTION 8: HORIZONTAL STACKED BAR — IMMUNE ENRICHED/DEPLETED PROPORTIONS #
#******************************************************************************#

meta_data_pancancer <- read.xlsx(paste0(data_path, "TCGA.PanAtlas.Metadata.xlsx")) %>%
  dplyr::select(Sample_ID, Project_ID) %>%
  dplyr::mutate(Sample_ID = make.names(Sample_ID))

s <- "gsva"

group <- read.xlsx(paste0(output_path, "Two_Group_Classification_from_Individual_cancers.xlsx"),
                   sheet = paste0(s, ".classification")) %>%
  dplyr::left_join(meta_data_pancancer, by = c("Sample" = "Sample_ID")) %>%
  dplyr::group_by(Project_ID) %>%
  dplyr::count(Class) %>%
  dplyr::mutate(percent       = round(100 * n / sum(n, na.rm = TRUE), digits = 2),
                label_percent = paste0(percent, "%"))

ggplot2::ggplot(data = group, aes(x = percent, y = Project_ID, fill = Class)) +
  geom_bar(stat = "identity", position = "stack", width = 0.95) +
  theme_classic() +
  scale_fill_manual(values = c("#CB181D", "#A6D854")) +
  geom_text(aes(x = percent, label = label_percent),
            position = position_stack(vjust = 0.5),
            fontface = "bold", colour = "white", size = 4, check_overlap = TRUE) +
  ggplot2::labs(title = "", fill = "Immune Class",
                x = "Percent composition", y = "")

ggsave(paste0(output_path, "zPercentages_", s, ".tiff"))

# BUG FIX: was saving `group_data` (undefined at this point) — now saves `group`
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = paste0(s, ".classification"))
openxlsx::writeData(wb, sheet = paste0(s, ".classification"), x = group, rowNames = FALSE)
openxlsx::saveWorkbook(wb, file = paste0(output_path, "zPercentages.xlsx"), overwrite = TRUE)

#******************************************************************************#
#       SECTION 9: NPEPPS IN scRNASeq (FEATURE PLOT + DOT PLOT)                #
#******************************************************************************#

source("/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")
source("/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

integrated_seurat <- readRDS(paste0(seurat_results, "Simon_integrated_seurat_snn.rds"))
Seurat::Idents(integrated_seurat) <- "subtype"

Seurat::FeaturePlot(
  object    = integrated_seurat,
  features  = "NPEPPS",
  reduction = "umap",
  cols      = c("grey", viridis::viridis(n = 10, option = "C", direction = -1)),
  min.cutoff = "q10",
  pt.size   = 0.1,
  order     = TRUE,
  label     = FALSE,
  raster    = FALSE,
  combine   = TRUE) + my_theme

ggsave(paste0(output_path, "NPEPPS.Umap.tiff"))

Seurat::DotPlot(
  object   = subset(integrated_seurat, celltype == "Epithelial"),
  features = "NPEPPS",
  assay    = "RNA",
  dot.min  = 0, dot.scale = 6,
  scale    = TRUE, scale.by = "size",
  scale.min = 0,  scale.max = 100) +
  ggplot2::geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.25) +
  ggplot2::scale_colour_distiller(type = "div", palette = "RdYlGn", direction = -1) +
  ggplot2::guides(size = guide_legend(
    override.aes = list(shape = 21, colour = "black", fill = "white", stroke = 0.75))) +
  ggplot2::theme(axis.text.x = element_text(
    family = "sans", face = "plain", colour = "black", size = 10,
    hjust = 0.5, vjust = 0.5, angle = 45))

ggsave(paste0(output_path, "NPEPPS.Dotplot.tiff"))

#******************************************************************************#
#     SECTION 10: DEGs — KRT13 vs CDH12  &  NON-IMMUNE vs IMMUNE (scRNASeq)   #
#******************************************************************************#

source("/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")
source("/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

species     <- "Homo sapiens"
annotations <- get_annotations(species)

integrated_seurat <- readRDS(paste0(seurat_results, "Simon_integrated_seurat_snn.rds"))

for (r in c("cold_hot", "non.immune_immune")){

  if (r == "cold_hot"){
    integrated_seurat@meta.data <- integrated_seurat@meta.data %>%
      dplyr::mutate(Condition = dplyr::case_when(
        subtype == "KRT13_Epithelial" ~ "KRT13_Epithelial",
        subtype == "CDH12_Epithelial" ~ "CDH12_Epithelial",
        TRUE ~ "Other_cell_types"))
    subtypes <- c("CDH12_Epithelial", "KRT13_Epithelial")
  } else {
    integrated_seurat@meta.data <- integrated_seurat@meta.data %>%
      dplyr::mutate(Condition = dplyr::case_when(
        celltype %in% c("Lymphocyte", "Myeloid") ~ "Immune",
        TRUE ~ "Non-Immune"))
    subtypes <- c("Immune", "Non-Immune")
  }

  meta_data_full <- data.frame()
  read_data_full <- data.frame(SYMBOL = rownames(integrated_seurat@assays$RNA$counts))

  for (s in subtypes){

    subset_seurat <- subset(x = integrated_seurat,
                            subset = Condition == s & !is.na(integrated_seurat@meta.data$Sample))

    meta_data <- subset_seurat@meta.data %>%
      dplyr::distinct(Sample, .keep_all = TRUE) %>%
      dplyr::filter(!is.na(Sample)) %>%
      dplyr::mutate(Batch = 1, Sample_ID = paste0(Sample, "_", s)) %>%
      dplyr::select(Sample_ID, Batch, Condition)

    samples <- subset_seurat@meta.data %>%
      dplyr::select(Sample) %>%
      dplyr::filter(!is.na(Sample)) %>%
      unlist(use.names = FALSE) %>%
      unique()

    read_data <- data.frame(matrix(NA,
                                   nrow = nrow(subset_seurat@assays$RNA$counts),
                                   ncol = nrow(meta_data)))
    rownames(read_data) <- rownames(subset_seurat@assays$RNA$counts)
    colnames(read_data) <- samples

    for (i in samples){
      cells_subset   <- rownames(subset_seurat@meta.data %>% dplyr::filter(Sample == i))
      subset_counts  <- data.frame(subset_seurat@assays$RNA$counts[, cells_subset])
      read_data[, i] <- rowSums(subset_counts)
    }

    colnames(read_data) <- paste0(colnames(read_data), "_", s)
    read_data <- read_data %>% tibble::rownames_to_column("SYMBOL")

    meta_data_full <- dplyr::bind_rows(meta_data_full, meta_data)
    read_data_full <- dplyr::left_join(read_data_full, read_data, by = "SYMBOL")
  }

  if (r == "cold_hot"){
    Comparisons <- list(Variable    = "Condition",
                        Target      = "KRT13_Epithelial",
                        Reference   = "CDH12_Epithelial",
                        lfc.cutoff  = 0,
                        padj.cutoff = 0.1)
  } else {
    Comparisons <- list(Variable    = "Condition",
                        Target      = "Non-Immune",
                        Reference   = "Immune",
                        lfc.cutoff  = 0,
                        padj.cutoff = 0.1)
  }

  meta_data <- prep_metadata(meta_data_full, read_data_full)
  read_data <- prep_readdata(read_data_full, meta_data)
  l         <- check_data(read_data, meta_data)
  meta_data <- l[[2]]
  read_data <- l[[1]]

  for (n in seq_along(Comparisons$Variable)){

    meta_data_comp <- meta_data %>% dplyr::mutate(id = get(Comparisons$Variable[n]))

    approach <- "DESeq2"
    if (length(unique(meta_data_comp$Batch)) > 1){
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = read_data,
                                            colData   = meta_data_comp,
                                            design    = ~ Batch + id)
    } else {
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = read_data,
                                            colData   = meta_data_comp,
                                            design    = ~ id)
    }

    prefix <- "zSimon"
    dds    <- run_deseq2(dds, meta_data_comp, annotations, Comparisons, n, approach, prefix, output_path)
  }

  # Seurat FindAllMarkers
  if (r == "cold_hot"){
    subset_seurat <- subset(x = integrated_seurat,
                            subset = Condition %in% c("CDH12_Epithelial", "KRT13_Epithelial"))
  } else {
    subset_seurat <- subset(x = integrated_seurat,
                            subset = Condition %in% c("Non-Immune", "Immune"))
  }

  Seurat::DefaultAssay(subset_seurat) <- "RNA"
  Seurat::Idents(subset_seurat)       <- "Condition"

  all_markers <- Seurat::FindAllMarkers(
    object             = subset_seurat,
    assay              = "RNA",
    features           = NULL,
    logfc.threshold    = 0.25,
    test.use           = "wilcox",
    slot               = "data",
    min.pct            = 0.1,
    min.diff.pct       = 0.1,
    node               = NULL,
    verbose            = TRUE,
    only.pos           = TRUE,
    max.cells.per.ident = Inf,
    random.seed        = 1,
    latent.vars        = NULL,
    min.cells.feature  = 3,
    min.cells.group    = 1,
    pseudocount.use    = 1,
    mean.fxn           = NULL,
    fc.name            = NULL,
    base               = 2,
    return.thresh      = 0.01,
    densify            = FALSE)

  t1 <- length(intersect(all_markers$gene, annotations$ENSEMBL_SYMBOL))
  t2 <- length(intersect(all_markers$gene, annotations$ENTREZ_SYMBOL))

  if (t1 >= t2){
    all_markers <- all_markers %>%
      dplyr::mutate(pct.1  = dplyr::if_else(pct.1 == 0, 0.001, pct.1),
                    pct.2  = dplyr::if_else(pct.2 == 0, 0.001, pct.2),
                    ratio  = pct.1 / pct.2) %>%
      dplyr::left_join(unique(annotations[, c("ENSEMBL_SYMBOL", "CHR", "DESCRIPTION")]),
                       by = c("gene" = "ENSEMBL_SYMBOL")) %>%
      dplyr::relocate(cluster, gene, CHR, avg_log2FC, p_val, p_val_adj, pct.1, pct.2, ratio, DESCRIPTION)
  } else {
    all_markers <- all_markers %>%
      dplyr::mutate(pct.1  = dplyr::if_else(pct.1 == 0, 0.001, pct.1),
                    pct.2  = dplyr::if_else(pct.2 == 0, 0.001, pct.2),
                    ratio  = pct.1 / pct.2) %>%
      dplyr::left_join(unique(annotations[, c("ENTREZ_SYMBOL", "CHR", "DESCRIPTION")]),
                       by = c("gene" = "ENTREZ_SYMBOL")) %>%
      dplyr::relocate(cluster, gene, CHR, avg_log2FC, p_val, p_val_adj, pct.1, pct.2, ratio, DESCRIPTION)
  }

  filename <- dplyr::case_when(
    r == "cold_hot"           ~ "zSimon.DEGs.KRT13_Epithelial_vs_CDH12_Epithelial_Seurat.xlsx",
    r == "non.immune_immune"  ~ "zSimon.DEGs.Non-Immune_vs_Immune_Seurat.xlsx")

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "All_Markers")
  openxlsx::writeData(wb, "All_Markers", x = all_markers)
  openxlsx::saveWorkbook(wb, file = paste0(output_path, filename), overwrite = TRUE)
}

#******************************************************************************#
#       SECTION 11: HIT FILTERING AND IDENTIFICATION                            #
#******************************************************************************#

source("/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

# ---- Reactome immune pathways ----
reactome_all_pathways <- read.xlsx(paste0(supp_path, "ReactomePathways.gmt.xlsx"), colNames = FALSE)
colnames(reactome_all_pathways) <- c("Pathway", "ID",
                                     paste0("Gene", seq_len(ncol(reactome_all_pathways) - 2)))

reactome_search_immune <- read.xlsx(paste0(supp_path, "Reactome_Immune.xlsx"))

custom_set1 <- reactome_all_pathways %>%
  dplyr::filter(ID %in% reactome_search_immune$ID)

my_patterns <- c("Immune","Immuno","Interferon","Interleukin","Cytokine",
                 "Toll","Complement","Inflamma","Antigen","FC",
                 "Co-stimulation","Co-inhibition")

custom_set2 <- reactome_all_pathways %>%
  dplyr::filter(grepl(paste0(my_patterns, collapse = "|"), Pathway, ignore.case = TRUE)) %>%
  dplyr::distinct_at("Pathway", .keep_all = TRUE)

reactome_immune_pathways <- dplyr::bind_rows(custom_set1, custom_set2) %>%
  dplyr::distinct_at("ID", .keep_all = TRUE)

immune_genes <- reactome_immune_pathways %>%
  dplyr::select(-Pathway, -ID) %>%
  unlist(use.names = FALSE) %>%
  unique()

# ---- TCGA DEGs ----
tcga_degs_up <- read.xlsx(paste0(output_path, "zSummary_of_DEGs.xlsx")) %>%
  dplyr::rename(SYMBOL = 1) %>%
  dplyr::filter(n_UP >= 10, n_DOWN <= 1)

tcga_degs_down <- read.xlsx(paste0(output_path, "zSummary_of_DEGs.xlsx")) %>%
  dplyr::rename(SYMBOL = 1) %>%
  dplyr::filter(n_UP < 4, n_DOWN > 1)

# ---- DGIdb druggable genes ----
druggable_genes <- read.xlsx(paste0(supp_path, "DGIdb_categories.xlsx")) %>%
  dplyr::rename(Class = 2) %>%
  dplyr::filter(Class == "DRUGGABLE GENOME")

# ---- Protein coding filter ----
annotations    <- get_annotations("Homo sapiens")
protein_coding <- annotations %>%
  dplyr::filter(ENSEMBL_BIOTYPE == "protein_coding" | ENTREZ_BIOTYPE == "protein-coding")

tcga_degs_up   <- tcga_degs_up %>%
  dplyr::filter(SYMBOL %in% c(protein_coding$ENSEMBL_SYMBOL,
                               protein_coding$ENTREZ_SYMBOL,
                               protein_coding$ENSEMBL_ID))

tcga_degs_down <- tcga_degs_down %>%
  dplyr::filter(SYMBOL %in% c(protein_coding$ENSEMBL_SYMBOL,
                               protein_coding$ENTREZ_SYMBOL,
                               protein_coding$ENSEMBL_ID))

hits_up <- tcga_degs_up[, 1:4] %>%
  dplyr::mutate(
    Pass = "TCGA",
    Pass = dplyr::case_when(SYMBOL %in% immune_genes      & Pass == "TCGA"         ~ "TCGA+Reactome",           TRUE ~ Pass),
    Pass = dplyr::case_when(SYMBOL %in% druggable_genes$SYMBOL & Pass == "TCGA+Reactome" ~ "TCGA+Reactome+Druggable", TRUE ~ Pass))

hits_down <- tcga_degs_down[, 1:4] %>%
  dplyr::mutate(
    Pass = "TCGA",
    Pass = dplyr::case_when(SYMBOL %in% immune_genes      & Pass == "TCGA"         ~ "TCGA+Reactome",           TRUE ~ Pass),
    Pass = dplyr::case_when(SYMBOL %in% druggable_genes$SYMBOL & Pass == "TCGA+Reactome" ~ "TCGA+Reactome+Druggable", TRUE ~ Pass))

#******************************************************************************#
#       SECTION 12: SURVIVAL ANALYSIS — ICB DATASETS (INDIVIDUAL GENES)        #
#******************************************************************************#

source("/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

survival_params <- list(
  plot_by            = NA,
  split_by           = NA,
  split_plot         = FALSE,
  multiple_cutoff    = FALSE,
  stratify_criteria  = "q",
  reference          = "LOW",
  conf_interval      = FALSE,
  plot_curve         = TRUE,
  plot_risk_table    = TRUE,
  legend_title       = "Expression",
  legend_label       = c("High", "Low"),
  color_palette      = c("#d73027", "#0c2c84"),
  plot_all_bins      = FALSE,
  plot_all_quartiles = FALSE,
  gene_sig_score     = FALSE)

gse <- c("GSE78220", "GSE91061", "GSE96619", "GSE115821", "GSE131521",
         "IMvigor210", "IMvigor010")

for (proj in gse){

  meta_data <- read.xlsx(paste0(data_path, proj, ".Metadata.xlsx"))
  read_data <- openxlsx::read.xlsx(xlsxFile = paste0(data_path, proj, ".Normalized.counts.xlsx"))

  meta_data <- meta_data %>%
    dplyr::mutate(Sample_ID = make.names(Sample_ID),
                  Time      = as.numeric(Time)) %>%
    dplyr::filter(Time > 0 & !is.na(Time)) %>%
    dplyr::distinct_at("Sample_ID", .keep_all = TRUE)

  if (proj == "IMvigor210"){
    meta_data <- meta_data %>%
      dplyr::filter(Tissue == "bladder", Received.platinum == "Y",
                    `Sample.collected.pre-platinum` == "N")
  }
  if (proj == "IMvigor010"){
    meta_data <- meta_data %>%
      dplyr::filter(ARM == "Atezolizumab (MPDL3280A) 1200 mg",
                    prior_neoadjuvant_chemotherapy == "YES")
  }

  norm_counts <- read_data %>%
    dplyr::mutate(SYMBOL = make.names(SYMBOL, unique = TRUE)) %>%
    dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
    tibble::column_to_rownames("SYMBOL") %>%
    dplyr::mutate(across(everything(), as.numeric))
  colnames(norm_counts) <- base::make.names(colnames(norm_counts))
  norm_counts <- norm_counts[, intersect(make.names(meta_data$Sample_ID), colnames(norm_counts))]
  norm_counts <- norm_counts[!rowSums(norm_counts, na.rm = TRUE) == 0, ]

  log_norm_counts <- log(1 + norm_counts, base = 2)
  t_med           <- base::apply(log_norm_counts, 1, median, na.rm = TRUE)
  log_norm_counts <- base::sweep(log_norm_counts, 1, t_med, "-")

  hits <- read.xlsx(paste0(output_path, "zFinal_Hits.xlsx"), sheet = "seurat_hits") %>%
    dplyr::filter(stringr::str_detect(Pass, "Reactome"))
  plot_genes <- intersect(hits$SYMBOL, rownames(log_norm_counts))

  expr_df <- prep_expr_df(log_norm_counts, meta_data, plot_genes, survival_params)

  stats <- list(gene = c(), group = c(),
                lower_cutoff = c(), middle_cutoff = c(), upper_cutoff = c(),
                HR = c(), CI_lower = c(), CI_upper = c(),
                logrank = c(), reg_logrank.late = c(),
                Gehan_Breslow.early = c(), Tarone_Ware.early = c(),
                Peto_Peto.early = c(), modified_Peto_Peto = c(),
                Fleming_Harrington = c())

  classification_df <- expr_df %>% dplyr::select(Sample_ID) %>% dplyr::mutate(Dummy_col = 0)

  for (gene in plot_genes){

    prefix  <- paste0(proj, "_", gene)
    summary <- plot_survival(expr_df, gene, survival_params, prefix, output_path)

    class_df <- summary[[1]] %>%
      dplyr::select(Sample_ID, all_of(gene), model) %>%
      dplyr::rename(!!paste0(gene, "_model") := model)

    classification_df <- classification_df %>%
      dplyr::left_join(class_df, by = "Sample_ID")

    stats$gene              <- c(stats$gene,              summary[[2]]$gene)
    stats$group             <- c(stats$group,             summary[[2]]$group)
    stats$lower_cutoff      <- c(stats$lower_cutoff,      summary[[2]]$lower)
    stats$middle_cutoff     <- c(stats$middle_cutoff,     summary[[2]]$middle)
    stats$upper_cutoff      <- c(stats$upper_cutoff,      summary[[2]]$upper)
    stats$HR                <- c(stats$HR,                summary[[2]]$HR)
    stats$CI_lower          <- c(stats$CI_lower,          summary[[2]]$CI_lower)
    stats$CI_upper          <- c(stats$CI_upper,          summary[[2]]$CI_upper)
    stats$logrank           <- c(stats$logrank,           summary[[2]]$logrank)
    stats$reg_logrank.late  <- c(stats$reg_logrank.late,  summary[[2]]$reg_logrank.late)
    stats$Gehan_Breslow.early <- c(stats$Gehan_Breslow.early, summary[[2]]$Gehan_Breslow.early)
    stats$Tarone_Ware.early   <- c(stats$Tarone_Ware.early,   summary[[2]]$Tarone_Ware.early)
    stats$Peto_Peto.early     <- c(stats$Peto_Peto.early,     summary[[2]]$Peto_Peto.early)
    stats$modified_Peto_Peto  <- c(stats$modified_Peto_Peto,  summary[[2]]$modified_Peto_Peto)
    stats$Fleming_Harrington  <- c(stats$Fleming_Harrington,  summary[[2]]$Fleming_Harrington)
  }

  stats_df <- data.frame(stats)

  norm_df <- norm_counts[intersect(plot_genes, rownames(norm_counts)), ] %>%
    t() %>% data.frame() %>% tibble::rownames_to_column("SYMBOL")

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Summary");        openxlsx::writeData(wb, "Summary",        stats_df)
  openxlsx::addWorksheet(wb, "Classification"); openxlsx::writeData(wb, "Classification", classification_df)
  openxlsx::addWorksheet(wb, "Norm_counts");    openxlsx::writeData(wb, "Norm_counts",    norm_df)
  openxlsx::saveWorkbook(wb, file = paste0(output_path, proj, "_Individual_stats.xlsx"), overwrite = TRUE)
}

#******************************************************************************#
#     SECTION 13: SURVIVAL — CTL GENE SIGNATURE SCORE (IMvigor)                #
#******************************************************************************#

source("/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

survival_params <- list(
  plot_by            = NA, split_by = NA, split_plot = FALSE,
  multiple_cutoff    = FALSE, stratify_criteria = "q", reference = "LOW",
  conf_interval      = FALSE, plot_curve = TRUE, plot_risk_table = TRUE,
  legend_title       = "Expression", legend_label = c("High", "Low"),
  color_palette      = c("#d73027", "#0c2c84"),
  plot_all_bins      = FALSE, gene_sig_score = TRUE)

ctl_genes <- c("IFNG","GZMA","GZMB","PRF1","GZMK","ZAP70",
               "GNLY","FASLG","TBX21","EOMES","CD8A","CD8B")

for (proj in c("IMvigor210", "IMvigor010")){

  meta_data <- read.xlsx(paste0(data_path, proj, ".Metadata.xlsx"))
  read_data <- openxlsx::read.xlsx(xlsxFile = paste0(data_path, proj, ".Normalized.counts.xlsx"))

  meta_data <- meta_data %>%
    dplyr::mutate(Sample_ID = make.names(Sample_ID), Time = as.numeric(Time)) %>%
    dplyr::filter(Time > 0 & !is.na(Time)) %>%
    dplyr::distinct_at("Sample_ID", .keep_all = TRUE)

  if (proj == "IMvigor210"){
    meta_data <- meta_data %>%
      dplyr::filter(Tissue == "bladder", Received.platinum == "Y",
                    `Sample.collected.pre-platinum` == "N")
  }
  if (proj == "IMvigor010"){
    meta_data <- meta_data %>%
      dplyr::filter(ARM == "Atezolizumab (MPDL3280A) 1200 mg",
                    prior_neoadjuvant_chemotherapy == "YES")
  }

  norm_counts <- read_data %>%
    dplyr::mutate(SYMBOL = make.names(SYMBOL, unique = TRUE)) %>%
    dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
    tibble::column_to_rownames("SYMBOL") %>%
    dplyr::mutate(across(everything(), as.numeric))
  colnames(norm_counts) <- base::make.names(colnames(norm_counts))
  norm_counts <- norm_counts[, intersect(make.names(meta_data$Sample_ID), colnames(norm_counts))]
  norm_counts <- norm_counts[!rowSums(norm_counts, na.rm = TRUE) == 0, ]

  log_norm_counts <- log(1 + norm_counts, base = 2)
  t_med           <- base::apply(log_norm_counts, 1, median, na.rm = TRUE)
  log_norm_counts <- base::sweep(log_norm_counts, 1, t_med, "-")

  plot_genes <- intersect(ctl_genes, rownames(log_norm_counts))
  expr_df    <- prep_expr_df(log_norm_counts, meta_data, plot_genes, survival_params)

  gene    <- "combined.exp"
  prefix  <- paste0(proj, "_", gene)
  summary <- plot_survival(expr_df, gene, survival_params, prefix, output_path)

  classification_df <- summary[[1]] %>% dplyr::select(Sample_ID, combined.exp, model, everything())
  stats_df          <- as.data.frame(summary[[2]])

  norm_df <- norm_counts[intersect(plot_genes, rownames(norm_counts)), ] %>%
    t() %>% data.frame() %>% tibble::rownames_to_column("SYMBOL")

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Summary");        openxlsx::writeData(wb, "Summary",        stats_df)
  openxlsx::addWorksheet(wb, "Classification"); openxlsx::writeData(wb, "Classification", classification_df)
  openxlsx::addWorksheet(wb, "Norm_counts");    openxlsx::writeData(wb, "Norm_counts",    norm_df)
  openxlsx::saveWorkbook(wb, file = paste0(output_path, proj, "CTL.score_q_Stats.xlsx"), overwrite = TRUE)
}

#******************************************************************************#
#     SECTION 14: SURVIVAL — IMMUNE ENRICHED vs IMMUNE DEPLETED                 #
#******************************************************************************#

source("/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

survival_params <- list(
  plot_by            = NA, split_by = NA, split_plot = FALSE,
  multiple_cutoff    = FALSE, stratify_criteria = "none", reference = "Immune Enriched",
  conf_interval      = FALSE, plot_curve = TRUE, plot_risk_table = TRUE,
  legend_title       = "Immune Status",
  legend_label       = c("Immune Depleted", "Immune Enriched"),
  color_palette      = c("#CB181D", "#A6D854"),
  plot_all_bins      = FALSE, plot_all_quartiles = FALSE, gene_sig_score = FALSE)

meta_data_pancancer <- read.xlsx(paste0(data_path, "TCGA.PanAtlas.Metadata.xlsx"))
gse <- c(make.names(unique(meta_data_pancancer$Project_ID)),
         "GSE78220", "GSE91061", "GSE96619", "GSE115821", "GSE131521",
         "IMvigor210", "IMvigor010")

for (proj in gse){

  meta_data  <- read.xlsx(paste0(data_path, proj, ".Metadata.xlsx"))
  group_data <- read.xlsx(paste0(output_path, proj, ".Scores.&.Heatmap.details.xlsx"),
                          sheet = "NPEPPS") %>%
    dplyr::select(Sample, gsva.classification, gsva.mean)

  meta_data <- meta_data %>%
    dplyr::mutate(Sample_ID = make.names(Sample_ID)) %>%
    dplyr::left_join(group_data, by = c("Sample_ID" = "Sample")) %>%
    dplyr::mutate(model = gsva.classification,
                  Time   = as.numeric(Time),
                  Status = as.numeric(Status)) %>%
    dplyr::filter(Time > 0 & !is.na(Time)) %>%
    dplyr::distinct_at("Sample_ID", .keep_all = TRUE)

  if (proj == "IMvigor210"){
    meta_data <- meta_data %>%
      dplyr::filter(Tissue == "bladder", Received.platinum == "Y",
                    `Sample.collected.pre-platinum` == "N")
  }
  if (proj == "IMvigor010"){
    meta_data <- meta_data %>%
      dplyr::filter(ARM == "Atezolizumab (MPDL3280A) 1200 mg",
                    prior_neoadjuvant_chemotherapy == "YES")
  }

  expr_df <- meta_data
  gene    <- ""
  prefix  <- paste0(proj, "_", gene)
  summary <- plot_survival(expr_df, gene, survival_params, prefix, output_path)

  classification_df <- summary[[1]] %>% dplyr::select(Sample_ID, model, everything())
  stats_df          <- as.data.frame(summary[[2]])

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Summary");        openxlsx::writeData(wb, "Summary",        stats_df)
  openxlsx::addWorksheet(wb, "Classification"); openxlsx::writeData(wb, "Classification", classification_df)
  openxlsx::saveWorkbook(wb, file = paste0(output_path, proj, "Immune_Stats.xlsx"), overwrite = TRUE)
}

#******************************************************************************#
#     SECTION 15: STACKED BAR — DRUG RESPONSE IN IMvigor210                    #
#******************************************************************************#

source("/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

meta_data  <- read.xlsx(paste0(data_path, "IMvigor210.Metadata.xlsx"))
group_data <- read.xlsx(paste0(output_path, "IMvigor210.Scores.&.Heatmap.details.xlsx"),
                        sheet = "NPEPPS") %>%
  dplyr::select(Sample, gsva.classification, gsva.mean)

plot_data <- meta_data %>%
  dplyr::mutate(Sample_ID = make.names(Sample_ID)) %>%
  dplyr::filter(Tissue == "bladder", Received.platinum == "Y",
                `Sample.collected.pre-platinum` == "N") %>%
  dplyr::left_join(group_data, by = c("Sample_ID" = "Sample")) %>%
  dplyr::mutate(model = gsva.classification) %>%
  dplyr::group_by(model) %>%
  dplyr::count(Best.Confirmed.Overall.Response) %>%
  dplyr::filter(Best.Confirmed.Overall.Response != "NE") %>%
  dplyr::mutate(percent       = round(100 * n / sum(n, na.rm = TRUE), digits = 2),
                label_percent = paste0(percent, "%")) %>%
  data.frame()

ggplot2::ggplot(data = plot_data,
                aes(x = factor(model), y = percent,
                    fill = Best.Confirmed.Overall.Response)) +
  geom_bar(stat = "identity", width = 0.65) +
  theme_classic() +
  scale_fill_brewer(palette = "Set1", aesthetics = "fill") +
  geom_text(aes(y = percent, label = label_percent),
            position   = position_stack(vjust = 0.5),
            fontface   = "bold", colour = "white", size = 4, check_overlap = TRUE) +
  ggplot2::labs(title = "", fill = "Immune Response",
                y = "Percent composition", x = "")

ggsave(paste0(output_path, "Stacked_bar_Imvigor210.tiff"))

#******************************************************************************#
#     SECTION 16: CORRELATION — GENE EXPRESSION vs IMMUNE SCORE & CTL SCORE    #
#******************************************************************************#

source("/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

# Load median centered log norm counts
log_norm_counts <- openxlsx::read.xlsx(
  xlsxFile = paste0(output_path, "TCGA.PanAtlas.median_centered_log_norm_counts.xlsx"))
colnames(log_norm_counts)[1] <- "SYMBOL"

meta_data_normal <- openxlsx::read.xlsx(
  xlsxFile = paste0(data_path, "TCGA.PanAtlas.Metadata.xlsx"), sheet = "Normal")
meta_data_tumor  <- openxlsx::read.xlsx(
  xlsxFile = paste0(data_path, "TCGA.PanAtlas.Metadata.xlsx"), sheet = "Tumor")

ctl_genes  <- c("IFNG","GZMA","GZMB","PRF1","GZMK","ZAP70",
                "GNLY","FASLG","TBX21","EOMES","CD8A","CD8B")
plot_genes <- intersect(ctl_genes, log_norm_counts$SYMBOL)

ctl_survival_params <- list(
  plot_by = NA, split_by = NA, split_plot = FALSE, multiple_cutoff = FALSE,
  stratify_criteria = "none", reference = "Immune Enriched",
  conf_interval = FALSE, plot_curve = TRUE, plot_risk_table = TRUE,
  legend_title = "Immune Status",
  legend_label = c("Immune Depleted", "Immune Enriched"),
  color_palette = c("#CB181D", "#A6D854"),
  plot_all_bins = FALSE, plot_all_quartiles = FALSE, gene_sig_score = TRUE)

log_mat_rownames <- log_norm_counts %>% tibble::column_to_rownames("SYMBOL")
expr_df_tumor  <- prep_expr_df(log_mat_rownames, meta_data_tumor,  plot_genes, ctl_survival_params)
expr_df_normal <- prep_expr_df(log_mat_rownames, meta_data_normal, plot_genes, ctl_survival_params)
expr_df_tumor  <- expr_df_tumor  %>% dplyr::rename(CTL.score = combined.exp) %>% dplyr::select(Sample_ID, CTL.score)
expr_df_normal <- expr_df_normal %>% dplyr::rename(CTL.score = combined.exp) %>% dplyr::select(Sample_ID, CTL.score)

classification_df <- openxlsx::read.xlsx(
  xlsxFile = paste0(output_path, "Two_Group_Classification_from_Individual_cancers.xlsx"))
lower_cutoff <- quantile(classification_df$Mean, c(0.33, 0.66))[[1]]
upper_cutoff <- quantile(classification_df$Mean, c(0.33, 0.66))[[2]]

DEGs_df <- openxlsx::read.xlsx(xlsxFile = paste0(output_path, "zSummary_of_DEGs.xlsx"))
colnames(DEGs_df)[1] <- "SYMBOL"
# BUG FIX: was calling get_annotations() with no args — requires species argument
annotations    <- get_annotations("Homo sapiens")
protein_coding <- annotations %>%
  dplyr::filter(stringr::str_detect(ENSEMBL_BIOTYPE, "protein") |
                  stringr::str_detect(ENTREZ_BIOTYPE, "protein"))

DEGs_vec <- DEGs_df %>%
  dplyr::filter(n_UP > 10, n_DOWN <= 1) %>%
  dplyr::pull(SYMBOL) %>%
  intersect(union(protein_coding$ENSEMBL_SYMBOL, protein_coding$ENTREZ_SYMBOL))
DEGs_vec <- c(DEGs_vec, "CD8A", "CD8B")

df_normal <- log_norm_counts %>%
  dplyr::filter(SYMBOL %in% DEGs_vec) %>%
  tibble::column_to_rownames("SYMBOL") %>%
  t() %>% data.frame() %>%
  tibble::rownames_to_column("Sample_ID") %>%
  dplyr::inner_join(meta_data_normal %>% dplyr::select(Sample_ID, Project_ID), by = "Sample_ID") %>%
  dplyr::inner_join(expr_df_normal,  by = "Sample_ID") %>%
  dplyr::select(Project_ID, Sample_ID, CTL.score, everything())

df_tumor <- log_norm_counts %>%
  dplyr::filter(SYMBOL %in% DEGs_vec) %>%
  tibble::column_to_rownames("SYMBOL") %>%
  t() %>% data.frame() %>%
  tibble::rownames_to_column("Sample_ID") %>%
  dplyr::inner_join(meta_data_tumor %>% dplyr::select(Sample_ID, Project_ID), by = "Sample_ID") %>%
  dplyr::left_join(classification_df %>% dplyr::select(Sample, Class, Mean),
                   by = c("Sample_ID" = "Sample")) %>%
  dplyr::inner_join(expr_df_tumor, by = "Sample_ID") %>%
  dplyr::mutate(Class.33 = dplyr::case_when(
    Mean <= lower_cutoff ~ Class,
    Mean >= upper_cutoff ~ Class,
    TRUE ~ "")) %>%
  dplyr::select(Project_ID, Sample_ID, Class, Class.33, Mean, CTL.score, everything())

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Normal"); openxlsx::writeData(wb, "Normal", df_normal)
openxlsx::addWorksheet(wb, "Tumor");  openxlsx::writeData(wb, "Tumor",  df_tumor)
openxlsx::saveWorkbook(wb, file = paste0(output_path, "zHit.Expression.Tumor.Normal.xlsx"), overwrite = TRUE)

# Per-cancer correlation
DEGs <- intersect(make.names(DEGs_vec), make.names(colnames(df_tumor)))
df   <- data.frame()

for (f in unique(meta_data_tumor$Project_ID)){

  df_sub <- df_tumor %>% dplyr::filter(Project_ID == f)

  gene <- spearman.p.immune <- spearman.r.immune <- pearson.p.immune <- pearson.r.immune <-
    spearman.p.ctl <- spearman.r.ctl <- pearson.p.ctl <- pearson.r.ctl <-
    spearman.p.immune.33 <- spearman.r.immune.33 <- pearson.p.immune.33 <- pearson.r.immune.33 <-
    spearman.p.ctl.33 <- spearman.r.ctl.33 <- pearson.p.ctl.33 <- pearson.r.ctl.33 <- c()

  for (i in seq_along(DEGs)){

    gene <- c(gene, DEGs[i])
    cat(i, ":", DEGs[i], "\n")

    run_cor <- function(df_in, x_col, y_col, method){
      res <- cor.test(x = df_in[[x_col]], y = df_in[[y_col]], method = method)
      list(p = res$p.value, r = res$estimate)
    }

    r <- run_cor(df_sub, "Mean",      DEGs[i], "spearman"); spearman.p.immune <- c(spearman.p.immune, r$p); spearman.r.immune <- c(spearman.r.immune, r$r)
    r <- run_cor(df_sub, "Mean",      DEGs[i], "pearson");  pearson.p.immune  <- c(pearson.p.immune,  r$p); pearson.r.immune  <- c(pearson.r.immune,  r$r)
    r <- run_cor(df_sub, "CTL.score", DEGs[i], "spearman"); spearman.p.ctl    <- c(spearman.p.ctl,    r$p); spearman.r.ctl    <- c(spearman.r.ctl,    r$r)
    r <- run_cor(df_sub, "CTL.score", DEGs[i], "pearson");  pearson.p.ctl     <- c(pearson.p.ctl,     r$p); pearson.r.ctl     <- c(pearson.r.ctl,     r$r)

    df_33 <- df_sub %>% dplyr::filter(nchar(Class.33) != 0)
    r <- run_cor(df_33, "Mean",      DEGs[i], "spearman"); spearman.p.immune.33 <- c(spearman.p.immune.33, r$p); spearman.r.immune.33 <- c(spearman.r.immune.33, r$r)
    r <- run_cor(df_33, "Mean",      DEGs[i], "pearson");  pearson.p.immune.33  <- c(pearson.p.immune.33,  r$p); pearson.r.immune.33  <- c(pearson.r.immune.33,  r$r)
    r <- run_cor(df_33, "CTL.score", DEGs[i], "spearman"); spearman.p.ctl.33    <- c(spearman.p.ctl.33,    r$p); spearman.r.ctl.33    <- c(spearman.r.ctl.33,    r$r)
    r <- run_cor(df_33, "CTL.score", DEGs[i], "pearson");  pearson.p.ctl.33     <- c(pearson.p.ctl.33,     r$p); pearson.r.ctl.33     <- c(pearson.r.ctl.33,     r$r)
  }

  df <- dplyr::bind_rows(df, data.frame(
    Cancer = f, Gene = gene,
    spearman.p.immune, spearman.r.immune, pearson.p.immune, pearson.r.immune,
    spearman.p.ctl, spearman.r.ctl, pearson.p.ctl, pearson.r.ctl,
    spearman.p.immune.33, spearman.r.immune.33, pearson.p.immune.33, pearson.r.immune.33,
    spearman.p.ctl.33, spearman.r.ctl.33, pearson.p.ctl.33, pearson.r.ctl.33))
}

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Summary"); openxlsx::writeData(wb, "Summary", df)
openxlsx::saveWorkbook(wb, file = paste0(output_path, "Correlation.xlsx"), overwrite = TRUE)

#******************************************************************************#
#       SECTION 17: ORA PATHWAY ANALYSIS                                        #
#******************************************************************************#

source("/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

species  <- "Homo sapiens"
gmt_dir  <- paste0("/home/kailasamms/projects/RNASeq/GSEA_genesets/Human")
gmt_files <- list.files(gmt_dir, full.names = TRUE)

up_df   <- data.frame()
down_df <- data.frame()

metadata <- read.xlsx(paste0(data_path, "TCGA.PanAtlas.Metadata.xlsx"), sheet = "Tumor")

for (proj in make.names(unique(metadata$Project_ID))){

  df <- read.xlsx(paste0(output_path, proj,
                          ".DEGs.gsva.classification_Immune Depleted_vs_Immune Enriched_DESeq2_.xlsx"))

  input_genes_up   <- df %>% dplyr::filter(padj <= 0.05, log2FoldChange >  0.58) %>% dplyr::pull(ENSEMBL_SYMBOL)
  input_genes_down <- df %>% dplyr::filter(padj <= 0.05, log2FoldChange < -0.58) %>% dplyr::pull(ENSEMBL_SYMBOL)
  universe_genes   <- df %>% dplyr::pull(ENSEMBL_SYMBOL)

  ora_result_up   <- ora(input_genes_up,   universe_genes, gmt_files) %>% dplyr::mutate(ProjectID = proj)
  ora_result_down <- ora(input_genes_down, universe_genes, gmt_files) %>% dplyr::mutate(ProjectID = proj)

  up_df   <- dplyr::bind_rows(up_df,   ora_result_up)
  down_df <- dplyr::bind_rows(down_df, ora_result_down)

  cat(proj, "| up:", nrow(up_df), "| down:", nrow(down_df), "\n")
}

up_df_sig   <- up_df   %>% dplyr::filter(p.adjust <= 0.05) %>% dplyr::add_count(Description)
down_df_sig <- down_df %>% dplyr::filter(p.adjust <= 0.05) %>% dplyr::add_count(Description)

# BUG FIX: common was computed AFTER up_df_unique/down_df_unique in original
# — moved common computation before filtering
common <- intersect(up_df_sig$Description, down_df_sig$Description)

up_df_unique <- up_df_sig %>%
  dplyr::filter(!(Description %in% common)) %>%
  dplyr::select(n, ProjectID, Description, everything()) %>%
  dplyr::arrange(desc(n), Description)

down_df_unique <- down_df_sig %>%
  dplyr::filter(!(Description %in% common)) %>%
  dplyr::select(n, ProjectID, Description, everything()) %>%
  dplyr::arrange(desc(n), Description)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "ORA_UP_UNIQUE");   openxlsx::writeData(wb, "ORA_UP_UNIQUE",   up_df_unique,   rowNames = FALSE)
openxlsx::addWorksheet(wb, "ORA_DOWN_UNIQUE"); openxlsx::writeData(wb, "ORA_DOWN_UNIQUE", down_df_unique, rowNames = FALSE)
openxlsx::addWorksheet(wb, "ORA_UP_ALL");      openxlsx::writeData(wb, "ORA_UP_ALL",      up_df,          rowNames = FALSE)
openxlsx::addWorksheet(wb, "ORA_DOWN_ALL");    openxlsx::writeData(wb, "ORA_DOWN_ALL",    down_df,        rowNames = FALSE)
openxlsx::saveWorkbook(wb,
                       file = paste0(output_path, "zPanCancer_Pathway_Analysis_Results.xlsx"),
                       overwrite = TRUE)

##


