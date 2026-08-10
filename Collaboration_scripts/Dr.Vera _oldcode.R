
#******************************************************************************#
#               SINGLE REPLICATE DIFFERENTIAL EXPRESSION ANALYSIS  (OUTDATED)           #
#******************************************************************************#

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")

data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Past/Dr.Vera/"
species <- "Homo sapiens"
annotations <- get_annotations(species)

# Read the normalized counts
norm_count <- read.xlsx(paste0(data_path, "Norm_counts_DESeq2.xlsx"))

# Calculate log2FC
deg_file <- norm_count %>%
  tibble::column_to_rownames("ENSEMBL_ID") %>%
  dplyr::select(everything(), -c(SYMBOL, ENSEMBL_SYMBOL)) %>%
  dplyr::mutate(across(.cols=everything(), ~ log2(1+.))) %>%
  dplyr::mutate(Line1_FC1 = Line1_GH1-Line1_Control1,
                Line1_FC2 = Line1_GH2-Line1_Control2,
                Line1_FC3 = Line1_GH3-Line1_Control3,
                Line2_FC1 = Line2_GH1-Line2_Control1,
                Line2_FC2 = Line2_GH2-Line2_Control2,
                Line2_FC3 = Line2_GH3-Line2_Control3) %>%
  dplyr::mutate(log2FoldChange_L1 = (Line1_FC1+Line1_FC2+Line1_FC3)/3,
                log2FoldChange_L2 = (Line2_FC1+Line2_FC2+Line2_FC3)/3)

# Calculate pvalues
subset <- deg_file[,13:15]
t1 <- c()
for (i in 1:nrow(subset)){
  t <- stats::t.test(x = subset[i,], alternative = "two.sided", mu = 0, paired = FALSE, var.equal = FALSE)
  t1 <- c(t1,t$p.value)
}

subset <- deg_file[,16:18]
t2 <- c()
for (i in 1:nrow(subset)){
  t <- stats::t.test(x = subset[i,], alternative = "two.sided", mu = 0, paired = FALSE, var.equal = FALSE)
  t2 <- c(t2,t$p.value)
}

stats_df <- data.frame(ID = rownames(subset), pval_L1 = t1, pval_L2 = t2)
stats_df$padj_L1 <- stats::p.adjust(p = stats_df$pval_L1, method = "fdr", n = length(stats_df$pval_L1))
stats_df$padj_L2 <- stats::p.adjust(p = stats_df$pval_L2, method = "fdr", n = length(stats_df$pval_L2))

deg_file <- deg_file %>% 
  tibble::rownames_to_column("ID") %>%
  dplyr::left_join(stats_df, by=c("ID"="ID"))

# Get gene symbols
deg_file <- add_annotation(deg_file, annotations)

# Remove unwanted genes
deg_file <- deg_file %>%
  dplyr::filter(!(SYMBOL %in% c("Y_RNA", "Metazoa_SRP")))

# Identify genes that have consistent trend in expression upon GH treatment
deg_file <- deg_file %>%
  dplyr::mutate(Status = dplyr::case_when(Line1_FC1>0 & Line1_FC2>0 & Line1_FC3>0 & Line2_FC1>0 & Line2_FC2>0 & Line2_FC3>0 ~ "Up",
                                          Line1_FC1<0 & Line1_FC2<0 & Line1_FC3<0 & Line2_FC1<0 & Line2_FC2<0 & Line2_FC3<0 ~ "Down",
                                          TRUE ~ "Ambiguous"))

# Save the results
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName="DEGs")
openxlsx::writeData(wb, sheet="DEGs", x=deg_file, rowNames=FALSE)
openxlsx::saveWorkbook(wb, file=paste0(data_path, "Similar.Trend.Genes.xlsx"), 
                       overwrite=TRUE)

#******************************************************************************#
#                                    HEATMAP                                   #
#******************************************************************************#

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")

proj <- "RNASeq_Vera"
data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Dr.Vera/"

proj <- "RNASeq_Bhowmick"
data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Bhowmick/"

# Read normalized counts
norm_counts <- read.xlsx(paste0(data_path, "Norm_counts_DESeq2.xlsx"))

# Read DEGs
df <- read.xlsx(paste0(data_path, "Similar.Trend.Genes.xlsx")) 
df1 <- read.xlsx(paste0(data_path, "RNASeq_Bhowmick.DEGs.RF_vs_FC_DESeq2.xlsx"))
df2 <- read.xlsx(paste0(data_path, "RNASeq_Bhowmick.DEGs.RFL_vs_FC_DESeq2.xlsx"))

# Heatmap parameters
plot_genes <- df$SYMBOL
disp_genes <- c()
file_suffix <- ""
output_path <- data_path
heatmap_params <- list(anno.row       = NULL,        # NULL, c("Group")
                       anno.column    = c("Treatment"),
                       row.split      = NA,     
                       col.split      = c("Treatment"),
                       row.cluster    = c("all"),           # c("alphabetical", "group", "all")
                       col.cluster    = c("alphabetical"),  # c("alphabetical", "group", "all")
                       discrete_panel = FALSE, 
                       log.transform  = TRUE,
                       scale          = TRUE,
                       border_color   = "white",
                       bar_width      = NA,              # NA , 5
                       bar_height     = NA,              # NA , 5
                       width          = 5,              # NA
                       height         = 5,              # NA 
                       matrix_color   = "rdbu",          # c("vrds", "rdbu")
                       expr_legend    = TRUE,  # set FALSE if it overlaps with annotation legends
                       file_format    = "tiff")

# Read metadata for columns
# MUST have column Sample_ID
# MUST have columns listed in heatmap_params$anno.column
metadata_column <- openxlsx::read.xlsx(xlsxFile=paste0(data_path, proj, "_Metadata.xlsx")) %>%
  dplyr::distinct_at("Sample_ID", .keep_all=TRUE)

# Read metadata for rows
# MUST have column SYMBOL
# MUST have columns listed in heatmap_params$anno.row
metadata_row <- data.frame(SYMBOL = "")
#metadata_row <- data.frame(SYMBOL = "", Gene = "")
#metadata_row <- data.frame(SYMBOL = norm_counts$SYMBOL, Gene = rep(x=c("A", "B"), times=c(700, 799)))

# Plot for L1 (norm_counts)
L1 <- norm_counts %>% 
  dplyr::select(SYMBOL, starts_with("Line1")) %>%
  dplyr::filter(SYMBOL %in% (df %>% dplyr::filter(pval_L1 < 0.05) %>% dplyr::select(SYMBOL) %>% unlist(use.names = FALSE))) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE)  %>%
  dplyr::filter(!(SYMBOL %in% c("Y_RNA", "Metazoa_SRP")))

metadata_row <- df %>% 
  dplyr::mutate(Group = dplyr::case_when((Line1_FC1 > 0 & Line1_FC2 > 0 & Line1_FC3 > 0) | (Line1_FC1 < 0 & Line1_FC2 < 0 & Line1_FC3 < 0) ~ "Consistent",
                                         TRUE ~ "Ambiguous")) %>%
  dplyr::select(SYMBOL, Group) %>%
  dplyr::filter(SYMBOL %in% L1$SYMBOL) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE)

plot_heatmap(L1, metadata_column, metadata_row, heatmap_params,
             plot_genes, disp_genes, "L1", output_path)

RF <- norm_counts %>% 
  dplyr::select(SYMBOL, contains(c("FC1", "FC2", "RF1", "RF2"))) %>%
  dplyr::filter(SYMBOL %in% (df1 %>% dplyr::filter(padj < 0.05) %>% dplyr::select(SYMBOL) %>% unlist(use.names = FALSE))) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE)  %>%
  dplyr::filter(!(SYMBOL %in% c("Y_RNA", "Metazoa_SRP")))

RFL <- norm_counts %>% 
  dplyr::select(SYMBOL, contains(c("FC1", "FC2", "RFL1", "RFL2"))) %>%
  dplyr::filter(SYMBOL %in% (df2 %>% dplyr::filter(padj < 0.05) %>% dplyr::select(SYMBOL) %>% unlist(use.names = FALSE))) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE)  %>%
  dplyr::filter(!(SYMBOL %in% c("Y_RNA", "Metazoa_SRP")))

plot_genes <- RF$SYMBOL
plot_heatmap(RF, metadata_column, metadata_row, heatmap_params,
             plot_genes, disp_genes, "RF", output_path)
plot_genes <- RFL$SYMBOL
plot_heatmap(RFL, metadata_column, metadata_row, heatmap_params,
             plot_genes, disp_genes, "RFL", output_path)

# Plot for L2 (norm_counts)
L2 <- norm_counts %>% 
  dplyr::select(SYMBOL, starts_with("Line2")) %>%
  dplyr::filter(SYMBOL %in% (df %>% dplyr::filter(pval_L2 < 0.05) %>% dplyr::select(SYMBOL) %>% unlist(use.names = FALSE))) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  dplyr::filter(!(SYMBOL %in% c("Y_RNA", "Metazoa_SRP")))

metadata_row <- df %>% 
  dplyr::mutate(Group = dplyr::case_when((Line2_FC1 > 0 & Line2_FC2 > 0 & Line2_FC3 > 0) | (Line2_FC1 < 0 & Line2_FC2 < 0 & Line2_FC3 < 0) ~ "Consistent",
                                         TRUE ~ "Ambiguous")) %>%
  dplyr::select(SYMBOL, Group) %>%
  dplyr::filter(SYMBOL %in% L2$SYMBOL) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE)

plot_heatmap(L2, metadata_column, metadata_row, heatmap_params,
             plot_genes, disp_genes, "L2", output_path)

# Plot for similar trend genes (norm_counts)
similar <- norm_counts %>% 
  dplyr::filter(ENSEMBL_ID %in% (df %>% dplyr::filter(Status != "Ambiguous") %>% dplyr::select(ENSEMBL_ID) %>% unlist(use.names = FALSE))) %>%
  dplyr::select(SYMBOL, starts_with(c("Line1", "Line2"))) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  dplyr::filter(!(SYMBOL %in% c("Y_RNA", "Metazoa_SRP")))

metadata_row <- df %>% 
  dplyr::mutate(Group = Status) %>%
  dplyr::select(SYMBOL, Group) %>%
  dplyr::filter(Group != "Ambiguous",  SYMBOL %in% similar$SYMBOL) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE)

heatmap_params <- list(anno.row       = c("Group"),        # NULL, c("Group")
                       anno.column    = c("Treatment", "Cell.Line"),
                       row.split      = c("Group"),        # c(), NULL, NA     
                       col.split      = c("Treatment"),    # c(), NULL, NA
                       row.cluster    = c("group"),        # c("alphabetical", "group", "all")
                       col.cluster    = c("group"),          # c("alphabetical", "group", "all")
                       log.transform  = TRUE,
                       scale          = TRUE,
                       border_color   = "white",
                       bar_width      = NA,              # NA , 5
                       bar_height     = NA,              # NA , 5
                       width          = 5,              # NA
                       height         = 5,              # NA 
                       matrix_color   = "rdbu",          # c("vrds", "rdbu")
                       expr_legend    = TRUE,            # set FALSE if it overlaps with annotation legends
                       file_format    = "tiff")

plot_heatmap(similar, metadata_column, metadata_row, heatmap_params,
             plot_genes, disp_genes, "similar.trend", output_path)

#******************************************************************************#
#                               PATHWAY ANALYSIS                               #
#******************************************************************************#

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")

species <- "Homo sapiens"
proj <- "RNASeq_Vera"
data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Dr.Vera/"

species <- "Mus musculus"
proj <- "RNASeq_Bhowmick"
data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Bhowmick/"

annotations <- get_annotations(species)

# Read DEGs
df <- read.xlsx(paste0(data_path, "Similar.Trend.Genes.xlsx"))
df1 <- read.xlsx(paste0(data_path, "RNASeq_Bhowmick.DEGs.RF_vs_FC_DESeq2.xlsx"))
df2 <- read.xlsx(paste0(data_path, "RNASeq_Bhowmick.DEGs.RFL_vs_FC_DESeq2.xlsx"))


# Read the gene set files
species <- "Human"
species <- "Mouse"
gmt_dir <- paste0("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/GSEA_genesets/", species, "/")
gmt_files <- list.files(gmt_dir, full.names = TRUE)

#*********************************ORA ANALYSIS*********************************#

# Define input genes (ONLY significant genes)
input_genes_L1 <- df %>% 
  dplyr::filter(pval_L1 <= 0.05) %>%
  dplyr::select(ENSEMBL_SYMBOL) %>%
  unlist(use.names=FALSE)

input_genes_L2 <- df %>% 
  dplyr::filter(pval_L2 <= 0.05) %>%
  dplyr::select(ENSEMBL_SYMBOL) %>%
  unlist(use.names=FALSE)

input_genes_trend <- df %>% 
  dplyr::filter(Status != "Ambiguous") %>%
  dplyr::select(ENSEMBL_SYMBOL) %>%
  unlist(use.names=FALSE)

input_genes_RF <- df1 %>% 
  dplyr::filter(padj < 0.05) %>%
  dplyr::select(ENSEMBL_SYMBOL) %>%
  unlist(use.names=FALSE)

input_genes_RFL <- df2 %>% 
  dplyr::filter(padj < 0.05) %>%
  dplyr::select(ENSEMBL_SYMBOL) %>%
  unlist(use.names=FALSE)

# Define universe genes
universe_genes <- df %>% 
  dplyr::select(ENSEMBL_SYMBOL) %>%
  unlist(use.names=FALSE)

universe_genes_RF <- df1 %>% 
  dplyr::select(ENSEMBL_SYMBOL) %>%
  unlist(use.names=FALSE)

universe_genes_RFL <- df2 %>% 
  dplyr::select(ENSEMBL_SYMBOL) %>%
  unlist(use.names=FALSE)

# Perform ORA analysis
ora_result_L1 <- ora(input_genes_L1, universe_genes, gmt_files)
ora_result_L2 <- ora(input_genes_L2, universe_genes, gmt_files)
ora_result_trend <- ora(input_genes_trend, universe_genes, gmt_files)
ora_result_RF <- ora(input_genes_RF, universe_genes_RF, gmt_files)
ora_result_RFL <- ora(input_genes_RFL, universe_genes_RFL, gmt_files)

#********************************FGSEA ANALYSIS********************************#

# Input dataframe MUST have SYMBOL, padj, log2FoldChange columns
DEGs_df_L1 <- df %>% 
  dplyr::select(ENSEMBL_SYMBOL, pval_L1, log2FoldChange_L1) %>%
  dplyr::rename(SYMBOL = ENSEMBL_SYMBOL, padj = pval_L1, log2FoldChange = log2FoldChange_L1) %>%
  dplyr::select(SYMBOL, padj, log2FoldChange)

DEGs_df_L2 <- df %>% 
  dplyr::select(ENSEMBL_SYMBOL, pval_L2, log2FoldChange_L2) %>%
  dplyr::rename(SYMBOL = ENSEMBL_SYMBOL, padj = pval_L2, log2FoldChange = log2FoldChange_L2) %>%
  dplyr::select(SYMBOL, padj, log2FoldChange)

DEGs_df_RF <- df1 %>%
  dplyr::rename(SYMBOL = ENSEMBL_SYMBOL) %>%
  dplyr::select(SYMBOL, padj, log2FoldChange)

DEGs_df_RFL <- df2 %>%
  dplyr::rename(SYMBOL = ENSEMBL_SYMBOL) %>%
  dplyr::select(SYMBOL, padj, log2FoldChange)

# Perform FGSEA analysis
fgsea_result_L1 <- fgsea(DEGs_df_L1, gmt_files, annotations)
fgsea_result_L2 <- fgsea(DEGs_df_L2, gmt_files, annotations)
fgsea_result_RF <- fgsea(DEGs_df_RF, gmt_files, annotations)
fgsea_result_RFL <- fgsea(DEGs_df_RFL, gmt_files, annotations)

#*********************************SAVE RESULTS*********************************#

# Create workbook to store results
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName=paste0("ORA"))
openxlsx::writeData(wb, sheet=paste0("ORA"), x=ora_result_L1, rowNames = FALSE)
openxlsx::addWorksheet(wb, sheetName=paste0("GSEA"))
openxlsx::writeData(wb, sheet=paste0("GSEA"), x=fgsea_result_L1[[1]], rowNames = FALSE)
openxlsx::addWorksheet(wb, sheetName=paste0("GSEA_genes"))
openxlsx::writeData(wb, sheet=paste0("GSEA_genes"), x=fgsea_result_L1[[4]], rowNames = FALSE)
openxlsx::saveWorkbook(wb, file = paste0(data_path, "Pathway_Analysis_Results_L1.xlsx"),
                       overwrite = TRUE)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName=paste0("ORA"))
openxlsx::writeData(wb, sheet=paste0("ORA"), x=ora_result_L2, rowNames = FALSE)
openxlsx::addWorksheet(wb, sheetName=paste0("GSEA"))
openxlsx::writeData(wb, sheet=paste0("GSEA"), x=fgsea_result_L2[[1]], rowNames = FALSE)
openxlsx::addWorksheet(wb, sheetName=paste0("GSEA_genes"))
openxlsx::writeData(wb, sheet=paste0("GSEA_genes"), x=fgsea_result_L2[[4]], rowNames = FALSE)
openxlsx::saveWorkbook(wb, file = paste0(data_path, "Pathway_Analysis_Results_L2.xlsx"),
                       overwrite = TRUE)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName=paste0("ORA"))
openxlsx::writeData(wb, sheet=paste0("ORA"), x=ora_result_trend, rowNames = FALSE)
openxlsx::saveWorkbook(wb, file = paste0(data_path, "Pathway_Analysis_Results_Trend.xlsx"),
                       overwrite = TRUE)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName=paste0("ORA"))
openxlsx::writeData(wb, sheet=paste0("ORA"), x=ora_result_RF, rowNames = FALSE)
openxlsx::addWorksheet(wb, sheetName=paste0("GSEA"))
openxlsx::writeData(wb, sheet=paste0("GSEA"), x=fgsea_result_RF[[1]], rowNames = FALSE)
openxlsx::addWorksheet(wb, sheetName=paste0("GSEA_genes"))
openxlsx::writeData(wb, sheet=paste0("GSEA_genes"), x=fgsea_result_RF[[4]], rowNames = FALSE)
openxlsx::saveWorkbook(wb, file = paste0(data_path, "Pathway_Analysis_Results_RF.xlsx"),
                       overwrite = TRUE)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName=paste0("ORA"))
openxlsx::writeData(wb, sheet=paste0("ORA"), x=ora_result_RFL, rowNames = FALSE)
openxlsx::addWorksheet(wb, sheetName=paste0("GSEA"))
openxlsx::writeData(wb, sheet=paste0("GSEA"), x=fgsea_result_RFL[[1]], rowNames = FALSE)
openxlsx::addWorksheet(wb, sheetName=paste0("GSEA_genes"))
openxlsx::writeData(wb, sheet=paste0("GSEA_genes"), x=fgsea_result_RFL[[4]], rowNames = FALSE)
openxlsx::saveWorkbook(wb, file = paste0(data_path, "Pathway_Analysis_Results_RFL.xlsx"),
                       overwrite = TRUE)

#******************************************************************************#

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")

proj <- "RNASeq_Vera"
data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Dr.Vera/"
species <- "Homo sapiens"
annotations <- get_annotations(species)
file_suffix <- ""

# List pathways to plot
plot_pathways <- c("REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS",
                   "REACTOME_COLLAGEN_DEGRADATION",
                   "REACTOME_COLLAGEN_CHAIN_TRIMERIZATION",
                   "GOBP_CELL_SUBSTRATE_ADHESION",
                   "GOBP_CELL_MATRIX_ADHESION",
                   "GOBP_COLLAGEN_FIBRIL_ORGANIZATION",
                   "GOBP_CELL_SUBSTRATE_JUNCTION_ORGANIZATION",
                   "GOBP_NEGATIVE_REGULATION_OF_CELL_SUBSTRATE_ADHESION",
                   "GOCC_COLLAGEN_TRIMER",
                   "GOCC_COMPLEX_OF_COLLAGEN_TRIMERS",
                   "GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT_CONFERRING_TENSILE_STRENGTH")

plot_pathways <- c("REACTOME_DNA_DOUBLE_STRAND_BREAK_RESPONSE",
                   "GOBP_NEGATIVE_REGULATION_OF_DNA_REPAIR",
                   "GOBP_REGULATION_OF_DNA_REPAIR")

# Read ORA results file
ora_results <- read.xlsx(paste0(data_path, "Pathway_Analysis_Results_Similar.Trend.xlsx")) %>%
  dplyr::filter(Description %in% plot_pathways)

#**********************Plot enrichment plot of ORA pathways********************#

plot_ora(ora_results, file_suffix, data_path)

#**********************Plot heatmap of genes in pathways**********************#

# Get genes in these pathways
ora_df <- convert_ora_genelist_to_df(ora_results)
ora_list <- ora_df %>% unlist(use.names = FALSE) %>% unique()
ora_list <- ora_list[!is.na(ora_list)]

# Read normalized counts
norm_counts1 <- read.xlsx(paste0(data_path, "Norm_counts_DESeq2.xlsx")) %>%
  dplyr::select(SYMBOL, starts_with(c("Line1"))) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  dplyr::filter(!(SYMBOL %in% c("Y_RNA", "Metazoa_SRP")))

norm_counts2 <- read.xlsx(paste0(data_path, "Norm_counts_DESeq2.xlsx")) %>%
  dplyr::select(SYMBOL, starts_with(c("Line2"))) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  dplyr::filter(!(SYMBOL %in% c("Y_RNA", "Metazoa_SRP")))

# Read DEGs
df1 <- read.xlsx(paste0(data_path, "Similar.Trend.Genes.xlsx")) %>%
  dplyr::filter(pval_L1 < 0.05 | pval_L2 < 0.05)

df2 <- read.xlsx(paste0(data_path, "Similar.Trend.Genes.xlsx")) %>%
  dplyr::filter(pval_L1 < 0.05 | pval_L2 < 0.05)

# Heatmap parameters
plot_genes <- ora_list
disp_genes <- ora_list
file_suffix <- ""
output_path <- data_path
heatmap_params <- list(anno.row       = c("Group"),        # NULL, c("Group")
                       anno.column    = c("Treatment", "Cell.Line"),
                       row.split      = c("Group"),        # c(), NULL, NA     
                       col.split      = c("Treatment"),    # c(), NULL, NA
                       row.cluster    = c("group"),        # c("alphabetical", "group", "all")
                       col.cluster    = c("alphabetical"), # c("alphabetical", "group", "all")
                       log.transform  = TRUE,
                       scale          = TRUE,
                       border_color   = NA, #"white",
                       bar_width      = 10,              # NA , 5
                       bar_height     = 10,              # NA , 5
                       width          = 5,              # NA
                       height         = 5,              # NA 
                       matrix_color   = "rdbu",         # c("vrds", "rdbu")
                       expr_legend    = TRUE,           # set FALSE if it overlaps with annotation legends
                       file_format    = "tiff")

metadata_column <- openxlsx::read.xlsx(xlsxFile=paste0(data_path, proj, "_Metadata.xlsx")) %>%
  dplyr::distinct_at("Sample_ID", .keep_all=TRUE)

metadata_row1 <- df1 %>% 
  dplyr::mutate(Group = Status) %>%
  dplyr::select(SYMBOL, Group) %>%
  dplyr::filter(Group != "Ambiguous",  SYMBOL %in% plot_genes) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE)

metadata_row2 <- df2 %>% 
  dplyr::mutate(Group = Status) %>%
  dplyr::select(SYMBOL, Group) %>%
  dplyr::filter(Group != "Ambiguous",  SYMBOL %in% plot_genes) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE)

plot_heatmap(norm_counts1, metadata_column, metadata_row1, heatmap_params,
             plot_genes, disp_genes, "ORA.Pathway.Genes.L1", output_path)

plot_heatmap(norm_counts2, metadata_column, metadata_row2, heatmap_params,
             plot_genes, disp_genes, "ORA.Pathway.Genes.L2", output_path)

#**********************Plot volcano of genes in pathways**********************#

# # Read DEGs
# df <- read.xlsx(paste0(data_path, "Similar.Trend.Genes.xlsx")) %>%
#   dplyr::filter(SYMBOL %in% ora_list)
# 
# # Define volcano plot parameters
# disp_genes <- c()
# file_suffix <- ""
# output_path <- data_path
# volcano_params <- list(Target          = "Treatment",
#                        Reference       = "Control", 
#                        lfc.cutoff      = 0.58,
#                        padj.cutoff     = 0.05)
# 
# # Input dataframe MUST have SYMBOL, padj, log2FoldChange columns
# L1 <- df %>% 
#   dplyr::select(SYMBOL, pval_L1, log2FoldChange_L1) %>%
#   dplyr::rename(SYMBOL = SYMBOL, padj = pval_L1, log2FoldChange = log2FoldChange_L1) %>%
#   dplyr::select(SYMBOL, padj, log2FoldChange)
# 
# plot_volcano(L1, volcano_params, disp_genes, "ORA.L1", data_path)
# 
# L2 <- df %>% 
#   dplyr::select(ENSEMBL_SYMBOL, pval_L2, log2FoldChange_L2) %>%
#   dplyr::rename(SYMBOL = ENSEMBL_SYMBOL, padj = pval_L2, log2FoldChange = log2FoldChange_L2) %>%
#   dplyr::select(SYMBOL, padj, log2FoldChange)
# 
# plot_volcano(L2, volcano_params, disp_genes, "ORA.L2", data_path)

#******************************************************************************#








# # Plot enrichment plots for significant pathways (padj < 0.05)
# sig_pathways <- gsea_results %>% 
#   dplyr::filter(padj < 0.05) %>% 
#   dplyr::slice_max(abs_NES, n = 12) %>%
#   dplyr::select(pathway) %>% 
#   unlist(use.names=FALSE)

# if (length(sig_pathways) > 0){
#   for (p in sig_pathways){
#     
#     # Create Enrichment plots for all significant pathways
#     fgsea::plotEnrichment(pathway = gsea_results %>%
#                             dplyr::filter(pathway == p) %>%
#                             dplyr::select(leadingEdge) %>%
#                             unlist(., use.names = FALSE),
#                           stats = DEGs_list,
#                           gseaParam = 1,
#                           ticksSize = 0.2)
#     
#     ggplot2::ggsave(filename = paste0("GSEA_NES_Plot_", p, ".tiff"),
#                     plot = last_plot(),
#                     device = "jpeg",
#                     path = data_path,
#                     scale = 1,
#                     width = 6,
#                     height = 7,
#                     units = c("in"),
#                     dpi = 600,
#                     limitsize = TRUE,
#                     bg = NULL)
#     
#     # # Find geneSetID corresponding to pathway in gseaobject from clusterprofiler
#     # enrichplot::gseaplot2(x = gsea,
#     #                       geneSetID = 1,
#     #                       title = p,
#     #                       color = "green",
#     #                       base_size = 11,
#     #                       rel_heights = c(1.5, 0.5, 1),
#     #                       subplots = 1:3,
#     #                       pvalue_table = FALSE,
#     #                       ES_geom = "line")
#   }
#   
#   # Plot bar plots for top 12 significant pathways by NES (padj < 0.05)
#   # It is difficult to control the width of bars in ggplot. Since, we plot
#   # top 12 pathways, we insert dummy entries to make the data frame have 12
#   # pathways "if" the data frame has less than 12 pathways
#   if (length(sig_pathways) < 12){
#     gsea_results <- gsea_results %>% 
#       dplyr::filter(padj < 0.05)
#     nrows <- nrow(gsea_results)
#     dummy_pathway   <- paste0("None.", seq(1:(12-nrows)))
#     dummy_NES       <- rep(0, times=12-nrows)
#     dummy_Direction <- rep("Downregulated", times=12-nrows)
#     dummy_df <- data.frame(pathway = dummy_pathway, NES = dummy_NES, Direction = dummy_Direction)
#     gsea_results <- dplyr::bind_rows(gsea_results, dummy_df)  
#   } else {
#     gsea_results <- gsea_results %>% 
#       dplyr::filter(padj < 0.05) %>% 
#       dplyr::slice_max(abs_NES, n = 12)
#   }
#   
#   # Modify pathway names to make the plot pretty
#   gsea_results_pretty <- gsea_results %>% 
#     data.frame() %>%
#     dplyr::mutate(pathway = gsub("HALLMARK_|SA_|SIG_|NABA_|GOBP_|GOMF_", "", pathway),
#                   pathway = gsub("_", " ", pathway),
#                   pathway = gsub("ENDOPLASMIC RETICULUM", "ER", pathway),
#                   #pathway = stringr::str_trunc(pathway, 45, "right"),
#                   #pathway = stringr::str_to_title(pathway),
#                   #length = stringr::str_length(pathway),
#                   pathway = stringr::str_wrap(pathway, width = 22))
#   
#   ggplot2::ggplot(data = gsea_results_pretty,
#                   aes(x = NES, y = reorder(pathway, NES), fill = Direction)) +
#     # fill = direction means direction will be arranged in alphabetical order.
#     # So, if you had labeled direction as "Upregulated" and "Downregulated",
#     # then first color in scale_fill_manual() will be assigned to
#     # "downregulated" and it will be labeled as "Activated in Males" in the
#     # plot. So, be careful.
#     ggplot2::geom_col(width = 0.75) +
#     ggplot2::theme_classic() +
#     ggplot2::labs(x = "Normalized Enrichment Score(NES)",
#                   y = "",
#                   title = "GSEA",
#                   fill = "") +
#     ggplot2::coord_cartesian(xlim = c(floor(-max(abs(gsea_results$NES), na.rm=TRUE)), ceiling(max(abs(gsea_results$NES), na.rm=TRUE)))) +
#     ggplot2::theme(#aspect.ratio = 2,
#       plot.title =   element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
#       plot.caption = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0),
#       axis.title.x = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
#       axis.title.y = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
#       axis.text.x =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
#       axis.text.y =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 1),
#       legend.title = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
#       legend.text =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
#       #legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid", colour ="darkblue"),
#       legend.position = "bottom",
#       legend.justification = "left",
#       legend.direction = "horizontal",
#       legend.key.height= unit(0.5, 'cm'),
#       legend.key.width= unit(1.25, 'cm')) +
#     ggplot2::scale_fill_manual(labels=c("Upregulated", "Downregulated"),
#                                values = c(RColorBrewer::brewer.pal(11, "RdYlBu")[c(1)],
#                                           RColorBrewer::brewer.pal(11, "RdYlBu")[c(11)]))
#   
#   # Save the plot
#   ggplot2::ggsave(filename = paste0("GSEA_", gmt_name, "_", file_suffix, ".tiff"),
#                   plot = last_plot(),
#                   device = "jpeg",
#                   path = data_path,
#                   scale = 1,
#                   width = 6,
#                   height = 7,
#                   units = c("in"),
#                   dpi = 600,
#                   limitsize = TRUE,
#                   bg = NULL)
# }