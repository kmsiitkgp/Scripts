# Define axis font etc to use in all plots
my_theme <- ggplot2::theme(aspect.ratio = 1,
                           plot.title =   element_text(family="sans", face="bold",  colour="black", size=15, hjust = 0.5),
                           legend.title = element_text(family="sans", face="bold",  colour="black", size=12, hjust = 0,   vjust = 1,   angle = 0),
                           legend.text =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                           axis.title.x = element_text(family="sans", face="bold",  colour="black", size=12, hjust = 0.5, vjust = 0,   angle = 0),
                           axis.title.y = element_text(family="sans", face="bold",  colour="black", size=12, hjust = 0.5, vjust = 1,   angle = 90),
                           axis.text.x =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5, vjust = 0.5, angle = 45),
                           axis.text.y =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5, vjust = 0.5, angle = 0),
                           strip.text.x = element_text(family="sans", face="bold",  colour="black", size=10, hjust = 0.5),
                           #legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid", colour ="darkblue"),
                           legend.position = "right",
                           legend.justification = "left",
                           legend.direction = "vertical",
                           legend.key.height = unit(0.5, 'cm'),
                           legend.key.width  = unit(0.5, 'cm'), 
                           legend.text.align = 0)

# Fig 1B, supplementary

# Define variables for functions from survival script
parent_path <- "C:/Users/KailasammS/OneDrive - Cedars-Sinai Health System/Desktop/Dr.Chao/"
results_path <- parent_path
gse <- "TCGA_BLCA"
subset_group <- NA
subset_value <- NA
plot_by <- "Expression"
split_by <- subset_group
combine_plot <- FALSE 
multiple_cutoff <- FALSE
stratify_criteria <- "o" #o"
reference <- "LOW"
confidence_interval <- FALSE
legend_title <- "Expression"
plot_risk_table <- TRUE
plot_curve <- TRUE
all_quartiles <- FALSE
gene_signature <- TRUE
legend_label <- c("High", "Low")              # Reference is Low
color_palette <- c("#d73027","#0c2c84")
variable_x <- "Time"
variable_y <- "Status" 

# Create gene lists for CAM, Scaffold, CAM_scaffold
c3_genes <- c("CADM2", "CDH12", "CNTNAP2", "CSMD1", "CSMD3", "CTNNA2", "CTNNA3",
              "DLG2", "DMD", "DPP10", "EYS", "GRID2", "LINGO1", "NRG1", "PCDH15",
              "PTPRD", "RASGEF1B", "RBFOX1", "SLC26A3", "TMEM132D")

ne_genes <- c("MSI1", "PLEKHG4B", "GNG4", "PEG10", "RND2", "APLP1", "SOX2",
              "TUBB2B", "CHGA", "CHGB", "NCAM1", "SYP")

scaffold_genes <- c("DLG2", "DLG4", "HOMER1", "HOMER2", "HOMER3", "SHANK1", 
                    "SHANK2", "SHANK3")
cam_genes <- c("CDH2", "CDH12", "LRFN1", "DAG1", "NLGN1", "NLGN2", "NLGN3",
               "LRRTM1", "LRRTM2", "ADGRL2", "ADGRL3", "EPHB2", "EPHB3", 
               "ADGRB1", "ADGRB3") #"NRXN1", "NRXN2", "NRXN3",
combo <- c(cam_genes, scaffold_genes)

ampa <- c("GRIA1", "GRIA2", "GRIA3", "GRIA4")
muscarinergic <- c("CHRM1", "CHRM2", "CHRM3", "CHRM4", "CHRM5")
gaba_b <- c("GABBR1", "GABBR2")
glycinergic <- c("GLRA1", "GLRA2", "GLRA3", "GLRA4", "GLRB")
kainate <- c("GRIK1", "GRIK2", "GRIK3", "GRIK4", "GRIK5")
nmda <- c("GRIN1", "GRIN2A", "GRIN2B", "GRIN2C", "GRIN2D", "GRIN3A", "GRIN3B")
adrenergic <- c("ADRA1A", "ADRA1B", "ADRA1D", "ADRA2A", "ADRA2B", "ADRA2C",
                "ADRB1", "ADRB2", "ADRB3")
ionotrophic <- c("HTR1A", "HTR1B", "HTR1D", "HTR1E", "HTR1F", "HTR2A", "HTR2B", 
                 "HTR2C", "HTR3A", "HTR3B", "HTR3C", "HTR3D", "HTR3E", "HTR4")
gaba_a <- c("GABRA1", "GABRA2", "GABRA3", "GABRA4", "GABRA5", "GABRA6", "GABRB1",
            "GABRB2", "GABRB3", "GABRD", "GABRE", "GABRG1", "GABRG2", "GABRG3", 
            "GABRP", "GABRQ", "GABRR1", "GABRR2", "GABRR3")
glutamate <- c("GRM1", "GRM2", "GRM3", "GRM4", "GRM5", "GRM6", "GRM7", "GRM8")
nicotinergic <- c("CHRNA1", "CHRNA10", "CHRNA2", "CHRNA3", "CHRNA4", "CHRNA5",
                  "CHRNA6", "CHRNA7", "CHRNA9", "CHRNB1", "CHRNB2", "CHRNB3",
                  "CHRNB4", "CHRND", "CHRNE", "CHRNG")

family_four <- c(adrenergic, ampa, muscarinergic, ionotrophic)
family_eleven <- c(adrenergic, ampa, muscarinergic, ionotrophic, gaba_a, gaba_b,
                   glycinergic, kainate, glutamate, nicotinergic, nmda)

# Import TCGA BLCA metadata
meta_data <- openxlsx::read.xlsx(paste0(parent_path,"TCGA_BLCA_Metadata.xlsx"))

meta_data <- meta_data %>% 
  dplyr::mutate(GEO_ID = make.names(names = GEO_ID)) %>%
  dplyr::mutate(Time = as.numeric(Time)) %>%
  dplyr::filter(Time > 0 & !is.na(Time)) %>%
  dplyr::distinct_at("GEO_ID", .keep_all = TRUE)

# Import TCGA BLCA normalized counts
normalized_counts <- read.xlsx(paste0(parent_path, "TCGA_BLCA_Normalized_Counts.xlsx"))
normalized_counts <- normalized_counts[, -c(2,3,4)]

normalized_counts <- normalized_counts %>%
  dplyr::mutate(SYMBOL = make.names(names = SYMBOL, unique = TRUE)) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  tibble::column_to_rownames(var = "SYMBOL") %>%
  dplyr::mutate(across(.cols = everything(), .fns = as.numeric))
colnames(normalized_counts) <- base::make.names(names = colnames(normalized_counts))
normalized_counts <- normalized_counts[,intersect(make.names(meta_data$GEO_ID), colnames(normalized_counts))]
normalized_counts <- normalized_counts[!rowSums(normalized_counts, na.rm=TRUE) == 0,]

# Perform log transformation 
normalized_counts <- log(1+normalized_counts, base=2)

# Perform median centering for each gene across samples
t <- base::apply(X=normalized_counts, MARGIN=1, FUN=median, na.rm=TRUE)
normalized_counts <- base::sweep(x=normalized_counts, MARGIN=1, FUN="-", STATS=t)

#******************************************************************************#

# Fig 1B Survival for CAMs, Scaffolds, CAMs+Scaffolds
# Fig 3A TCGA stage wise calculation

wb <- openxlsx::createWorkbook()
for (plot_genes in c("adrenergic", "ampa", "muscarinergic", "ionotrophic", 
                     "gaba_a", "gaba_b", "glycinergic", "kainate", "glutamate",
                     "nicotinergic", "nmda", "cam_genes", "scaffold_genes", 
                     "combo")){

  # Calculate z-score using the function described in above paper
  expr_df <- as.data.frame(advanced_Z(get(plot_genes), normalized_counts))
  
  # Merge expression data with survival data
  expr_df <- expr_df %>%
    data.frame() %>%
    dplyr::rename(combined.exp = identity(1)) %>%
    tibble::rownames_to_column("GEO_ID") %>%
    dplyr::inner_join(meta_data, by=c("GEO_ID"="GEO_ID"))
  
  prefix <- plot_genes 
  gene <- "combined.exp"
  if (nrow(expr_df > 0)){
    summary <- wrangle_data(expr_df, stratify_criteria,prefix)
  }
  
  surv_df <- summary[[1]] %>%
    dplyr::select(GEO_ID, combined.exp,	model, everything())
  stats_df <- as.data.frame(summary[[2]])
  
  # Save the results
  openxlsx::addWorksheet(wb, sheetName = paste0("Summary_",plot_genes))
  openxlsx::writeData(wb, sheet = paste0("Summary_",plot_genes), x = stats_df, rowNames = FALSE)
  openxlsx::addWorksheet(wb, sheetName = paste0(plot_genes))
  openxlsx::writeData(wb, sheet = paste0(plot_genes), x = surv_df, rowNames = FALSE)
}
openxlsx::saveWorkbook(wb, file = paste0(parent_path, "Survival.xlsx"), overwrite = TRUE)

# Compile data
new_metadata <- meta_data %>% 
  dplyr::select(GEO_ID, Stage) %>%
  dplyr::mutate(Stage = dplyr::case_when(grepl(pattern="T0", x=Stage) ~ 0,
                                         grepl(pattern="T1", x=Stage) ~ 1,
                                         grepl(pattern="T2", x=Stage) ~ 2,
                                         grepl(pattern="T3", x=Stage) ~ 3,
                                         grepl(pattern="T4", x=Stage) ~ 4,
                                         TRUE ~ 0)) %>%
  dplyr::filter(Stage %in% c(1,2,3,4))

for (plot_genes in c("adrenergic", "ampa", "muscarinergic", "ionotrophic", 
                     "gaba_a", "gaba_b", "glycinergic", "kainate", "glutamate",
                     "nicotinergic", "nmda")){
  
  data <- read.xlsx(paste0(parent_path, "Survival.xlsx"), 
                    sheet = plot_genes) %>%
    dplyr::select(GEO_ID, combined.exp) %>%
    dplyr::rename(!!plot_genes := "combined.exp")
  
  new_metadata <- new_metadata %>% 
    dplyr::left_join(data, by=c("GEO_ID"="GEO_ID"))
}

# Save the results
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "prism")
openxlsx::writeData(wb, sheet = "prism", x = new_metadata, rowNames = FALSE)
openxlsx::saveWorkbook(wb, file = paste0(parent_path, "Prism.xlsx"), overwrite = TRUE)

#******************************************************************************#

# Supplementary survival of individual genes
wb <- openxlsx::createWorkbook()
plot_genes <- c(combo)

# Since with optimal cutoff GRIA2 is slightly better than GRIA1, we plot differently
stratify_criteria <- "q" 
plot_genes <- c(ampa)

# Merge expression data with survival data
expr_df <- normalized_counts %>%
  t() %>%
  data.frame() %>%
  tibble::rownames_to_column("GEO_ID") %>%
  dplyr::select(GEO_ID, all_of(plot_genes)) %>%
  dplyr::inner_join(meta_data, by=c("GEO_ID"="GEO_ID"))
  

# Create a list to store survminer cutoffs, coxph stats, etc..
stats <- list("gene" = c(), 
              "group" = c(),
              "lower_cutoff" = c(),
              "middle_cutoff" = c(),
              "upper_cutoff" = c(),
              "HR" = c(), 
              "CI_lower" = c(), 
              "CI_upper" = c(), 
              "logrank" = c(), 
              "reg_logrank.late" = c(), 
              "Gehan_Breslow.early" = c(),
              "Tarone_Ware.early" = c(), 
              "Peto_Peto.early" = c(),  
              "modified_Peto_Peto" = c(), 
              "Fleming_Harrington" = c())

# Plot survival curves
prefix <- "TCGA_BLCA"
for (gene in intersect(unique(plot_genes), rownames(normalized_counts))) {
  summary <- wrangle_data(expr_df, stratify_criteria, prefix)
  stats$gene          <- c(stats$gene, summary[[2]]$gene)
  stats$group         <- c(stats$group, summary[[2]]$group)
  stats$lower_cutoff  <- c(stats$lower_cutoff, summary[[2]]$lower)
  stats$middle_cutoff <- c(stats$middle_cutoff, summary[[2]]$middle)
  stats$upper_cutoff  <- c(stats$upper_cutoff, summary[[2]]$upper)
  stats$HR            <- c(stats$HR, summary[[2]]$HR )
  stats$CI_lower      <- c(stats$CI_lower, summary[[2]]$CI_lower)
  stats$CI_upper      <- c(stats$CI_upper, summary[[2]]$CI_upper)
  stats$logrank       <- c(stats$logrank, summary[[2]]$logrank)
  stats$reg_logrank.late    <- c(stats$reg_logrank.late, summary[[2]]$reg_logrank.late)
  stats$Gehan_Breslow.early <- c(stats$Gehan_Breslow.early, summary[[2]]$Gehan_Breslow.early)
  stats$Tarone_Ware.early   <- c(stats$Tarone_Ware.early, summary[[2]]$Tarone_Ware.early)
  stats$Peto_Peto.early     <- c(stats$Peto_Peto.early, summary[[2]]$Peto_Peto.early)
  stats$modified_Peto_Peto   <- c(stats$modified_Peto_Peto, summary[[2]]$modified_Peto_Peto)
  stats$Fleming_Harrington   <-c(stats$Fleming_Harrington, summary[[2]]$Fleming_Harrington)
}

stats_df <- data.frame(stats)

# Save the results
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Summary")
openxlsx::writeData(wb, sheet = "Summary", x = stats_df, rowNames = FALSE)
openxlsx::saveWorkbook(wb, file = paste0(parent_path, "Individual_gene_stats.xlsx"), overwrite = TRUE)


#******************************************************************************#

# Fig 1A Oncoprint for CAMs, Scaffolds, CAMs+Scaffolds

# Import the dataframe containing scores with TCGA metadata
for (y in c("cam_genes", "scaffold_genes", "combo")){
  metadata <- read.xlsx(paste0(parent_path, "Survival.xlsx"),
                        sheet = y) %>%
    dplyr::select(GEO_ID, model, Age, Sex, Race, ajcc_pathologic_stage, ajcc_pathologic_m, Stage) %>%
    dplyr::rename("Score" = "model",
                  "Metastatic.stage"= "ajcc_pathologic_m",
                  "Tumor.stage"= "Stage",
                  "Stage" = "ajcc_pathologic_stage") %>%
    dplyr::mutate(Age = dplyr::case_when(Age < 65 ~ "< 65 yrs",
                                         Age >= 65 ~ ">= 65 yrs")) %>%
    tibble::column_to_rownames("GEO_ID") %>%
    dplyr::arrange(desc(Score))
  
  transformed_data <- metadata %>% 
    dplyr::mutate(Score = dplyr::case_when(Score == "LOW" ~ 0,
                                           Score == "HIGH" ~ 1),
                  Age = dplyr::case_when(Age == "< 65 yrs" ~ 0,
                                         Age == ">= 65 yrs" ~ 1),
                  Sex = dplyr::case_when(Sex == "Female" ~ 0,
                                         Sex == "Male" ~ 1),
                  Race = dplyr::case_when(Race == "Not Reported" ~ 0,
                                          Race == "Black/African American" ~ 1,
                                          Race == "Asian" ~ 2,
                                          Race == "White" ~ 3),
                  Metastatic.stage = dplyr::case_when(is.na(Metastatic.stage) ~ 0,
                                                      Metastatic.stage == "M0" ~ 1,
                                                      Metastatic.stage == "M1" ~ 2,
                                                      Metastatic.stage == "MX" ~ 3),
                  Tumor.stage = dplyr::case_when(grepl(pattern="T0", x=Tumor.stage) ~ 0,
                                                 grepl(pattern="T1", x=Tumor.stage) ~ 1,
                                                 grepl(pattern="T2", x=Tumor.stage) ~ 2,
                                                 grepl(pattern="T3", x=Tumor.stage) ~ 3,
                                                 grepl(pattern="T4", x=Tumor.stage) ~ 4,
                                                 TRUE ~ 5),
                  Stage = dplyr::case_when(is.na(Stage) ~ 0,
                                           Stage == "Stage I" ~ 1,
                                           Stage == "Stage II" ~ 2,
                                           Stage == "Stage III" ~ 3,
                                           Stage == "Stage IV" ~ 4))
  
  transformed_data <- t(transformed_data) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("SYMBOL") %>%
    dplyr::mutate(SYMBOL = gsub(pattern= "\\.", replacement= " ", x=SYMBOL))
  colnames(transformed_data)[1] <- "SYMBOL"
  
  metadata_column <- metadata %>% 
    tibble::rownames_to_column("Sample")
  colnames(metadata_column) <- gsub(pattern= "\\.", replacement= " ", x=colnames(metadata_column))
  
  metadata_row <- NULL
  plot_genes <- transformed_data$SYMBOL
  disp_genes <- plot_genes
  columns <- "Sample"
  
  perform_log_transform <- FALSE
  perform_scaling <- FALSE
  row_clustering <- FALSE    
  col_clustering <- FALSE    
  row_clustering_alphabetical <- FALSE
  col_clustering_alphabetical <- FALSE
  gaps_in_col <- TRUE
  gap_columns <- "Score"   # Irrelevant if gaps_in_col is FALSE
  gaps_in_row <- FALSE
  gap_rows <- "Pathway"    # Irrelevant if gaps_in_row is FALSE
  col_clustering_within_group <- FALSE
  row_clustering_within_group <- FALSE
  
  # List annotations you want on heatmap
  # NOTE: anno_columns MUST match one of the column names in metadata_column while
  # anno_rows MUST match one of the column names in metadata_row
  anno_columns <- base::rev(setdiff(colnames(metadata_column), "Sample"))
  anno_rows <- NA
  color_by_cols <- TRUE
  color_by_rows <- FALSE
  my_palette <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100)
  file_format <- ".tiff"
  file_suffix <- y
  bar_width <- NA
  bar_height <- NA
  expr_legend <- FALSE
  plot_heatmap(transformed_data, metadata_column, metadata_row, plot_genes, 
               disp_genes, file_suffix, file_format, results_path, 
               bar_width, bar_height, expr_legend)
}

#******************************************************************************#

# Calculate stats on Fig 1A
pvals <- c()
vars <- c()
sheets <- c()
for (y in c("cam_genes", "scaffold_genes", "combo")){
  
  data <- read.xlsx(paste0(parent_path, "Survival.xlsx"),
                    sheet = y)
  
  for (var in c("Sex", "Stage", "Race", "ajcc_pathologic_m", "ajcc_pathologic_stage")) {
    
    df <- as.data.frame.matrix(table(data %>% dplyr::select(all_of(var)) %>% unlist(use.names=FALSE),
                                     data %>% dplyr::select(model) %>% unlist(use.names=FALSE)))
    
    
    results <- chisq.test(data %>% dplyr::select(all_of(var)) %>% unlist(use.names=FALSE),
                          data %>% dplyr::select(model) %>% unlist(use.names=FALSE))
    
    pvals <- c(pvals, results[["p.value"]])
    vars <- c(vars, var)
    sheets <- c(sheets, y)
  }   
}

df <- data.frame(sheets, vars, pvals)

# Save the results
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Summary")
openxlsx::writeData(wb, sheet = "Summary", x = df, rowNames = FALSE)
openxlsx::saveWorkbook(wb, file = paste0(parent_path, "Fig3A_stats.xlsx"), overwrite = TRUE)

#******************************************************************************#

# Fig 1C UMAP of single nuclei data

proj <- "scRNASeq_Simon_from_Xingyu"
source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")
celltype <- NULL
my_palette <- c("#FF1F5B", "#F28522", "#009ADE", "#AF58BA", "#542788", "#00B000", 
                "#FFC61E", "#808080", "#BF812D", "#35978F", "#C51B7D", "#7FBC41")
# Load the seurat object
integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn",
                                    dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))

Seurat::DimPlot(object = integrated_seurat,
                reduction = "umap",
                cols = my_palette,
                group.by = "initial_celltype")

ggsave("Fig 1C.tiff",
       units = "in",
       width = 10,
       height = 10)

#******************************************************************************#

# Supplementary UMAP of CAM, scaffold, all 11 receptor family genes individually

proj <- "scRNASeq_Simon_from_Xingyu"
source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")
celltype <- NULL

# Load the seurat object
integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn",
                                    dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))

DefaultAssay(integrated_seurat) <- "RNA"

for (feature_subset in c("cam_genes", "scaffold_genes", "adrenergic", "ampa", 
                         "muscarinergic", "ionotrophic", "gaba_a", "gaba_b", 
                         "glycinergic", "kainate", "glutamate", "nicotinergic",
                         "nmda")){
  Seurat::FeaturePlot(object = integrated_seurat ,
                      slot = "data",
                      features = get(feature_subset),
                      cols =  c("grey", viridis(n = 10, option = "C", direction = -1)),
                      pt.size = 0.4,
                      order = TRUE,
                      min.cutoff = 'q10',
                      reduction = "umap",
                      raster = FALSE) +
    plot_layout(ncol = 4,
                nrow = 5)
  
  ggsave(paste0("Supplementary_UMAP_", feature_subset, ".tiff"),
         units = "in",
         width = 20,
         height = 30)
}
#******************************************************************************#

# Fig 1D, Supplementary UMAP of module scores for CAM, Scaffold, CAM+scaffold, 
# C3, 11 receptor families

proj <- "scRNASeq_Simon_from_Xingyu"
source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")
celltype <- NULL

# Load the seurat object
integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn",
                                    dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))

for (feature_subset in c("cam_genes", "scaffold_genes", "combo", 
                         "c3_genes", "adrenergic", "ampa", "muscarinergic", 
                         "ionotrophic", "gaba_a", "gaba_b", "glycinergic", 
                         "kainate", "glutamate", "nicotinergic", "nmda")){
  
  features <- list(base::intersect(get(feature_subset), rownames(integrated_seurat@assays$RNA$counts)))
  
  # Calculate module scores using seurat
  x <- Seurat::AddModuleScore(object = integrated_seurat,
                              features = features,
                              assay = "RNA",
                              slot = "data",
                              name = feature_subset)
  
  # Calculate module scores using ucell
  names(features) <- feature_subset
  x <- UCell::AddModuleScore_UCell(obj = x,
                                   features = features,
                                   assay = "RNA",
                                   slot = "data",
                                   name = "_UCell")
  
  # Plot module scores
  DefaultAssay(integrated_seurat) <- "RNA"
  p1 <- Seurat::FeaturePlot(object = x,
                            slot = "data",
                            features = paste0(feature_subset,"1"),
                            cols =  c("#440154FF", viridis(n = 10, option = "C", direction = 1)),
                            pt.size = 0.4,
                            order = TRUE,
                            min.cutoff = 'q10',
                            reduction = "umap",
                            raster = FALSE)
  
  p2 <- Seurat::FeaturePlot(object = x,
                            slot = "data",
                            features = paste0(feature_subset,"_UCell"),
                            cols =  c("#440154FF", viridis(n = 10, option = "C", direction = 1)),
                            pt.size = 0.4,
                            order = TRUE,
                            min.cutoff = 'q10',
                            reduction = "umap",
                            raster = FALSE)
  
  # Save plot
  ggplot2::ggsave(filename = paste0("Module_score_plot_", feature_subset, ".tiff"),
                  plot = p1+p2,
                  device = "jpeg",
                  width = 10,
                  height = 20,
                  units = c("in"),
                  bg = "white")
}

#******************************************************************************#

# Fig 3B, Supplementary dot plots of CAM+scaffold genes, family_four genes, 
# each of the 11 receptor families in 5 epithelial cell types

proj <- "scRNASeq_Simon_from_Xingyu"
source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")
celltype <- "Epithelial"

# Load the seurat object
integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn",
                                    dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))

Idents(integrated_seurat) <- "subtype"

# We re-order the active ident alphabetically
Idents(integrated_seurat) <- base::factor(x = integrated_seurat@active.ident, 
                                          levels = sort(levels(integrated_seurat@active.ident)))

for (feature_subset in c("combo", "family_four", "adrenergic", 
                         "ampa", "muscarinergic", "ionotrophic", "gaba_a", 
                         "gaba_b", "glycinergic", "kainate", "glutamate",
                         "nicotinergic", "nmda")){
  
  Seurat::DotPlot(object = integrated_seurat,
                  assay = "RNA",
                  features = get(feature_subset),
                  #cols = c("blue", "red"),
                  #col.min = -2,
                  #col.max = 2,
                  dot.min = 0,
                  dot.scale = 6,
                  scale = TRUE,
                  scale.by = "size",
                  scale.min = 0,
                  scale.max = 100) +
    coord_flip() +
    ggplot2::geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.25) + #stroke is width of circle
    ggplot2::scale_colour_distiller(type = "div", palette = "RdYlGn", direction = -1) +
    ggplot2::guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white", stroke=0.75))) +
    ggplot2::theme(axis.text.x =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5, vjust = 0.5, angle = 45))
  
  ggplot2::ggsave(filename = paste0("Dot_plot_", feature_subset, ".tiff"),
                  plot = last_plot(),
                  device = "jpeg",
                  width = length(unique(Idents(integrated_seurat))),
                  height = ceiling((length(get(feature_subset))+3)/4),
                  units = c("in"),
                  bg = "white")
}

#******************************************************************************#

# Calculate population frequencies

proj <- "scRNASeq_Simon_from_Xingyu"
source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")
celltype <- NULL

# Load the seurat object
integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn",
                                    dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
x <- integrated_seurat

for (feature_subset in c("cam_genes", "scaffold_genes", "c3_genes")){
  
  features <- list(base::intersect(get(feature_subset), rownames(integrated_seurat@assays$RNA$counts)))
  
  # Calculate module scores using seurat
  x <- Seurat::AddModuleScore(object = x,
                              features = features,
                              assay = "RNA",
                              slot = "data",
                              name = feature_subset)
  
  # Calculate module scores using ucell
  names(features) <- feature_subset
  x <- UCell::AddModuleScore_UCell(obj = x,
                                   features = features,
                                   assay = "RNA",
                                   slot = "data",
                                   name = "_UCell")
}

test <- x@meta.data[,c(25,27,29)]
test <- test %>% 
  dplyr::mutate(cam_genes1 = dplyr::if_else(cam_genes1 > 0, 1, 0), 
                scaffold_genes1 = dplyr::if_else(scaffold_genes1 > 0, 1, 0), 
                c3_genes1 = dplyr::if_else(c3_genes1 > 0, 1, 0))

class <- test %>%
  dplyr::mutate(class = dplyr::case_when(cam_genes1 == 1 & scaffold_genes1 == 1 & c3_genes1 == 1 ~ "csc",
                                         cam_genes1 == 1 & scaffold_genes1 == 0 & c3_genes1 == 0 ~ "cam",
                                         cam_genes1 == 0 & scaffold_genes1 == 1 & c3_genes1 == 0 ~ "scaffold",
                                         cam_genes1 == 0 & scaffold_genes1 == 0 & c3_genes1 == 1 ~ "c3",
                                         cam_genes1 == 1 & scaffold_genes1 == 1 & c3_genes1 == 0 ~ "cam_scaffold",
                                         cam_genes1 == 1 & scaffold_genes1 == 0 & c3_genes1 == 1 ~ "cam_c3",
                                         cam_genes1 == 0 & scaffold_genes1 == 1 & c3_genes1 == 1 ~ "scaffold_c3",
                                         cam_genes1 == 0 & scaffold_genes1 == 0 & c3_genes1 == 0 ~ "negative"))
class %>% 
  dplyr::count(class) %>%
  dplyr::mutate(percent_n = 100*n/sum(n))

################# NG expression in GSE124312

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")

anno <- read.table(file = "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/GSE124312/GSE124312_Vagus_AnnotationTable.txt",
                   header=TRUE, fill=TRUE, row.names=NULL)

# Use fill=FALSE. If some rows have missing data, then error will be thrown.
expr <-  read.table(file = "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/GSE124312/GSE124312_Vagus_ExpressionTable.txt",
                    header=TRUE, fill=FALSE, row.names=NULL)

anno <- anno %>% dplyr::filter(!is.na(nGene))
colnames(anno)[1] <- "Cell"
colnames(expr)[1] <- "SYMBOL"
colnames(expr) <- gsub(pattern="X0", replacement="0", x=colnames(expr))

putative_efferent <- c("Chat", "Mnx1", "Neurog2", "Trpv1", "Chrm1", "Chrm2",
                       "Chrm3", "Chrm4", "Chrm5", "Chrna1", "Chrna10", "Chrna2",
                       "Crna3", "Chrna4", "Chrna5", "Chrna6", "Chrna7", "Chrnb1",
                       "Chrnb2", "Chrnb3", "Chrnb4", "Chrne", "Chrng")

expr <- expr %>% 
  dplyr::filter(SYMBOL %in% putative_efferent) %>%
  tibble::column_to_rownames("SYMBOL") %>%
  t() %>%
  data.frame() %>%
  tibble::rownames_to_column("Cell") %>%
  dplyr::left_join(anno, by=c("Cell"="Cell")) %>%
  dplyr::select(Cell, Sample, Cluster, everything(), -c(nGene, nUMI, percent.mito))

# For single cell GSE145216
# plot_object <- subset(integrated_seurat, seurat_class=="Neurons")
# anno <- plot_object@meta.data %>% dplyr::select(Cell, Sample)
# expr <- plot_object@assays$RNA$data
# expr <- expr %>%
#   data.frame() %>%
#   tibble::rownames_to_column("SYMBOL") %>%
#   dplyr::filter(SYMBOL %in% putative_efferent) %>%
#   tibble::column_to_rownames("SYMBOL") %>%
#   t() %>%
#   data.frame() %>%
#   tibble::rownames_to_column("Cell") %>%
#   dplyr::left_join(anno, by=c("Cell"="Cell"))

mat <- expr %>%
  tibble::column_to_rownames("Cell") %>%
  dplyr::select(everything(), -c(Sample, Cluster))

my_palette <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)

if(max(mat) == 0){
  breaks <- c(seq(from = min(mat), to = 0, length.out = 100))
  my_palette <- my_palette[1:50]
} else if (min(mat) == 0){
  breaks <- c(seq(from = 0, to = max(mat), length.out = 100))
  my_palette <- my_palette[50:100]
} else if(min(mat) < -3 | max(mat) > 3){
  breaks <- c(seq(-1.5, 0, length.out = 50), seq(1.5/100, 1.5, length.out = 50))
} else{
  breaks <- c(seq(from = min(mat), to = 0, length.out = 50), seq(from = max(mat)/100, to = max(mat), length.out = 50))
}

mat <- mat[,sort(colnames(mat))]
pheatmap(mat, 
         scale="column", 
         breaks=breaks, show_rownames = FALSE, cluster_cols = FALSE, 
         filename = "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/GSE124312_Heatmap.tiff")

# Save the results
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Expr")
openxlsx::writeData(wb, sheet = "Expr", x = mat, rowNames = TRUE)
openxlsx::saveWorkbook(wb, file = "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/GSE124312_Heatmap.xlsx", overwrite = TRUE)

################# DRG expression in GSE249746 

expr <-  read.csv(file = "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/GSE249746.csv",
                  header=TRUE, fill=FALSE, row.names=NULL)

colnames(expr)[1] <- "SYMBOL"
colnames(expr) <- gsub(pattern="X0", replacement="0", x=colnames(expr))

expr <- expr %>% 
  tibble::column_to_rownames("SYMBOL")

# Normalize for seq depth and log1p
expr <- sweep(x=expr, MARGIN=2, STATS=colSums(expr), FUN = "/")
expr <- expr*10000
expr <- log1p(expr)


mat <- expr %>% 
  tibble::rownames_to_column("SYMBOL") %>%
  dplyr::filter(SYMBOL %in% toupper(putative_efferent)) %>%
  tibble::column_to_rownames("SYMBOL") %>%
  t()

my_palette <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)

if(max(mat) == 0){
  breaks <- c(seq(from = min(mat), to = 0, length.out = 100))
  my_palette <- my_palette[1:50]
} else if (min(mat) == 0){
  breaks <- c(seq(from = 0, to = max(mat), length.out = 100))
  my_palette <- my_palette[50:100]
} else if(min(mat) < -3 | max(mat) > 3){
  breaks <- c(seq(-1.5, 0, length.out = 50), seq(1.5/100, 1.5, length.out = 50))
} else{
  breaks <- c(seq(from = min(mat), to = 0, length.out = 50), seq(from = max(mat)/100, to = max(mat), length.out = 50))
}

mat <- mat[,sort(colnames(mat))]
pheatmap(mat, 
         scale="column", 
         breaks=breaks, show_rownames = FALSE, cluster_cols = FALSE, 
         filename = "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/GSE249746_Heatmap.tiff")

# Save the results
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Expr")
openxlsx::writeData(wb, sheet = "Expr", x = mat, rowNames = TRUE)
openxlsx::saveWorkbook(wb, file = "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/GSE249746_Heatmap.xlsx", overwrite = TRUE)

######Reviewer 4 Heatmap

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")

heatmap_params <- list(anno.row       = c("Class"),        # NULL, c("Group")
                       anno.column    = c("ajcc_pathologic_stage"),
                       row.split      = c("Class"),     
                       col.split      = NA, #c("ajcc_pathologic_stage"),
                       row.cluster    = c("alphabetical"),           # c("alphabetical", "group", "all")
                       col.cluster    = c("all"),           # c("alphabetical", "group", "all")
                       discrete_panel = TRUE,
                       log.transform  = TRUE,
                       scale          = TRUE,
                       border_color   = "white",
                       bar_width      = NA,              # NA , 5
                       bar_height     = NA,              # NA , 5
                       width          = 7,               # NA
                       height         = 7,               # NA 
                       matrix_color   = "rdbu",          # c("vrds", "rdbu")
                       expr_legend    = TRUE,  # set FALSE if it overlaps with annotation legends
                       file_format    = "tiff")

plot_genes <- c(cam_genes, scaffold_genes)
disp_genes <- plot_genes
output_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/"

metadata_row <- data.frame(plot_genes) %>%
  dplyr::mutate(Class = dplyr::case_when(plot_genes %in% cam_genes ~ "SynCAM",
                                         TRUE ~ "Scaffold")) %>%
  dplyr::rename(SYMBOL = plot_genes)

# Import TCGA BLCA metadata
parent_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/FromBox/Dr.Chao project/"
metadata_column <- openxlsx::read.xlsx(paste0(parent_path,"!TCGA_BLCA_Metadata.xlsx"))

metadata_column <- metadata_column %>% 
  dplyr::mutate(Sample_ID = make.names(names = Sample_id)) %>%
  dplyr::mutate(Time = as.numeric(Time)) %>%
  dplyr::filter(Time > 0 & !is.na(Time)) %>%
  dplyr::distinct_at("GEO_ID", .keep_all = TRUE) %>%
  dplyr::filter(!is.na(get(heatmap_params$anno.column)))

# Import TCGA BLCA normalized counts
normalized_counts <- read.xlsx(paste0(parent_path, "!TCGA_BLCA_Normalized_Counts.xlsx"))
normalized_counts <- normalized_counts[, -c(2,3,4)]

plot_heatmap(normalized_counts, metadata_column, metadata_row, heatmap_params,
             plot_genes, disp_genes, "TCGA.Stage", output_path)

######Reviewer 4 stage adjusted survival curves

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
output_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Folder8/"
parent_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/FromBox/Dr.Chao project/"

survival_params <- list(plot_by             = c(NA), 
                        split_by            = c("combined_stage"), #c("ajcc_pathologic_stage"),
                        split_plot          = TRUE,
                        multiple_cutoff     = TRUE,
                        stratify_criteria   = c("o"),
                        reference           = c("LOW"),
                        conf_interval       = FALSE,
                        plot_curve          = TRUE,
                        plot_risk_table     = TRUE,
                        legend_title        = "Expression",
                        legend_label        = c("High", "Low"),
                        color_palette       = c("#d73027","#0c2c84"),
                        plot_all_bins       = FALSE,
                        plot_all_quartiles  = FALSE,
                        gene_sig_score      = TRUE)

meta_data <- openxlsx::read.xlsx(paste0(parent_path,"TCGA_BLCA_Metadata.xlsx")) %>%
  dplyr::mutate(combined_stage = dplyr::case_when(ajcc_pathologic_stage %in% c("Stage I", "Stage II") ~ "StageI+II",
                                                  ajcc_pathologic_stage %in% c("Stage III", "Stage IV") ~ "StageIII+IV",
                                                  TRUE ~ "delete")) %>%
  dplyr::filter(combined_stage != "delete")
  
read_data <- read.xlsx(paste0(parent_path, "TCGA_BLCA_Normalized.xlsx"))

# Reformat metadata 
meta_data <- meta_data %>% 
  dplyr::mutate(Sample_ID = make.names(names = Sample_ID)) %>%
  dplyr::mutate(Time = as.numeric(Time)) %>%
  dplyr::filter(Time > 0 & !is.na(Time)) %>%
  dplyr::distinct_at("Sample_ID", .keep_all = TRUE)

# Reformat read data
norm_counts <- read_data %>%
  dplyr::mutate(SYMBOL = make.names(names = SYMBOL, unique = TRUE)) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  tibble::column_to_rownames(var = "SYMBOL") %>%
  dplyr::mutate(across(.cols = everything(), .fns = as.numeric))
colnames(norm_counts) <- base::make.names(names = colnames(norm_counts))
norm_counts <- norm_counts[,intersect(make.names(meta_data$Sample_ID), colnames(norm_counts))]
norm_counts <- norm_counts[!rowSums(norm_counts, na.rm=TRUE) == 0,]

log_norm_counts <- log(1+norm_counts, base=2)
t <- base::apply(X=log_norm_counts, MARGIN=1, FUN=median, na.rm=TRUE)
log_norm_counts <- base::sweep(x=log_norm_counts, MARGIN=1, FUN="-", STATS=t)

for (p in c("cam_genes", "scaffold_genes", "combo")){
  
  plot_genes <- get(p)
  
  # Generate expr_df
  expr_df <- prep_expr_df(log_norm_counts, meta_data, plot_genes, survival_params)
  
  # Create a list to store survminer cutoffs, coxph stats, etc..
  stats <- list("gene" = c(),
                "group" = c(),
                "lower_cutoff" = c(),
                "middle_cutoff" = c(),
                "upper_cutoff" = c(),
                "HR" = c(),
                "CI_lower" = c(),
                "CI_upper" = c(),
                "logrank" = c(),  
                "reg_logrank.late" = c(),
                "Gehan_Breslow.early" = c(),
                "Tarone_Ware.early"  = c(),
                "Peto_Peto.early" = c(),
                "modified_Peto_Peto"  = c(),
                "Fleming_Harrington" = c())
  
  # Create a dataframe to classification info
  classification_df <- expr_df %>% 
    dplyr::select(Sample_ID) %>%
    dplyr::mutate(Dummy_col = 0)
  
  if (survival_params$gene_sig_score == TRUE){
    
    gene <- "combined.exp"
    prefix <- p
    summary <- plot_survival(expr_df, gene, survival_params, prefix, output_path)
    
    # Merge classification info to parent dataframe
    classification_df <- summary[[1]] %>%
      dplyr::select(Sample_ID, combined.exp,	model, everything())
    
    # Merge all stats into a dataframe
    stats_df <- as.data.frame(summary[[2]])
    
    # # Create a dataframe of normalized counts of genes plotted
    # norm_df <- norm_counts[intersect(plot_genes, rownames(norm_counts)),] %>%
    #   t() %>%
    #   data.frame() %>%
    #   tibble::rownames_to_column("SYMBOL")
    
    # Save the results
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheetName = "Summary")
    openxlsx::writeData(wb, sheet = "Summary", x = stats_df)
    openxlsx::addWorksheet(wb, sheetName = "Classification")
    openxlsx::writeData(wb, sheet = "Classification", x = classification_df)
    # openxlsx::addWorksheet(wb, sheetName = "Norm_counts")
    # openxlsx::writeData(wb, sheet = "Norm_counts", x = norm_df)
    openxlsx::saveWorkbook(wb, file = paste0(output_path, p, "_Stats.xlsx"), overwrite = TRUE)
  }
}

##### 

remotes::install_github("cit-bioinfo/BLCAsubtyping", build_vignettes = TRUE)
library(BLCAsubtyping)
remotes::install_github("cit-bioinfo/consensusMIBC", build_vignettes = TRUE)
library(consensusMIBC)

#vignette("BLCAsubtyping")
# data(cit)
# cl <- classify(expMat = cit,
#                classification.systems = c("Baylor", "UNC", "MDA", "CIT", "Lund", "TCGA"))

#vignette("consensusMIBC")
# cl <- getConsensusClass(x, minCor = .2, gene_id = c("entrezgene", "ensembl_gene_id", "hgnc_symbol")[2])

# --- 1. Load Required Libraries (Run this first) ---
library(GEOquery)
library(limma) # Needed for exprs()
library(dplyr) # Needed for %>% , group_by, summarize, full_join
library(purrr) # Needed for reduce
library(tibble) # Needed for column_to_rownames
library(sva) # Needed for ComBat

# --- 2. Configuration ---
gse_list <- c("GSE5287", "GSE13507", "GSE31684", "GSE32894", "GSE48075", 
              "GSE52219", "GSE48276", "GSE87304", "GSE281689") #"GSE83586" <- no probeid even

# --- 3. Function to Process a Single GSE (Data Retrieval & Gene Summarization) ---
process_and_summarize_gse <- function(gse_id) {
  message(paste("\n--- Processing GSE:", gse_id, "---"))
  
  # a. Data Retrieval
  gset <- getGEO(gse_id, GSEMatrix = TRUE)
  if (is.list(gset)) {
    gset <- gset[[1]]
  }
  exprs_data <- exprs(gset)
  
  # b. Prefix column names with GSE ID (CRITICAL for merging)
  colnames(exprs_data) <- paste0(gse_id, "_", colnames(exprs_data))
  
  # c. Annotation Retrieval
  gpl_id <- annotation(gset)
  gpl <- getGEO(gpl_id)
  gpl_table <- Table(gpl)
  
  # d. Prepare Annotation Table
  if (any(colnames(gpl_table) %in% c("Symbol", "symbol", "Gene Name", "Gene Symbol"))){
  annot <- gpl_table %>%
    select(PROBEID = ID, GeneSymbol = matches("Symbol|symbol|Gene Name|Gene Symbol")) %>%
    filter(!is.na(GeneSymbol) & GeneSymbol != "" & GeneSymbol != "---") %>%
    mutate(GeneSymbol = gsub("///.*", "", GeneSymbol)) %>% # Clean up multiple symbols
    distinct(PROBEID, .keep_all = TRUE)
  } else{
    #"GSE87304" has gene symbols in ID column
    annot <- gpl_table %>%
      select(PROBEID = ID, GeneSymbol = ID) %>%
      filter(!is.na(GeneSymbol) & GeneSymbol != "" & GeneSymbol != "---") %>%
      mutate(GeneSymbol = gsub("///.*", "", GeneSymbol)) %>% # Clean up multiple symbols
      distinct(PROBEID, .keep_all = TRUE)
  }
  
  # e. Merge and Summarize to Gene Level (using mean)
  exprs_df <- data.frame(PROBEID = rownames(exprs_data), exprs_data)
  merged <- merge(exprs_df, annot, by = "PROBEID")
  
  gene_level <- merged %>%
    group_by(GeneSymbol) %>%
    summarize(across(where(is.numeric) & !starts_with("PROBEID"), mean, na.rm = TRUE), .groups = 'drop') %>%
    select(-any_of("PROBEID")) # Remove PROBEID column after merge
  
  return(gene_level)
}

# --- 4. Execution and Merging ---

# Process all GSEs and store results in a list
results_list <- lapply(gse_list, process_and_summarize_gse)

# Merge all gene-level data frames using GeneSymbol as the key
final_merged_df <- results_list %>%
  purrr::reduce(full_join, by = "GeneSymbol") %>%
  dplyr::rename(SYMBOL = GeneSymbol)

# --- 5. Prepare metadata ---
Sjodahl <- openxlsx::read.xlsx(
  xlsxFile = "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/16th. Oct. 2025/1-s2.0-S0302283819306955-mmc9.xlsx",
  sheet = "Sjodahl2012,2017")

cl_1750 <- openxlsx::read.xlsx(
  xlsxFile = "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/16th. Oct. 2025/1-s2.0-S0302283819306955-mmc9.xlsx",
  startRow = 2) %>%
  dplyr::select("ID", "Dataset", "Consensus_class", "Gender", "Age", "Stage") %>%
  dplyr::mutate(GSE = dplyr::case_when(Dataset == "Dyrskjot" ~ "GSE5287", #"GSE3167",
                                       Dataset == "Kim" ~ "GSE13507",
                                       Dataset == "Riester" ~ "GSE31684",
                                       Dataset == "Sjodahl2012" ~ "GSE32894",
                                       Dataset == "Choi MDA discovery" ~ "GSE48075",
                                       Dataset == "Choi MDA MVAC" ~ "GSE52219",
                                       Dataset == "Choi MDA validation" ~ "GSE48276",
                                       Dataset == "Choi Philadelphia NAC" ~ "GSE48276",
                                       Dataset == "Seiler" ~ "GSE87304")) %>%
  # Dataset == "Sjodahl2017" ~ "GSE83586")) %>%
  dplyr::filter(!is.na(GSE)) %>%
  tidyr::separate(
    col = ID, # The column to split
    into = c("Sample_ID", "GSM_temp"), # The names of the new columns
    sep = "\\.", # The delimiter (the comma)
    remove = FALSE # Keep the original 'Full_Name_Code' column
  ) %>%
  dplyr::mutate(GSM_temp = dplyr::case_when(is.na(GSM_temp) ~ Sample_ID,
                                            TRUE ~ GSM_temp),
                Sample_ID = GSM_temp) %>%
  dplyr::left_join(Sjodahl,by=c("Sample_ID"="Sample_ID")) %>%
  dplyr::mutate(GSM = dplyr::case_when(is.na(GSM) ~ GSM_temp,
                                       TRUE ~ GSM),
                Sample_ID = paste0(GSE, "_", GSM)) %>%
  dplyr::select(everything(), -GSM_temp) %>%
  dplyr::filter(Sample_ID %in% colnames(final_merged_df))



# --- 6. Plot heatmap for each dataset ---

heatmap_plot_list <- list()
for (gse_id in gse_list){
  
  c3_genes <- c("CADM2", "CDH12", "CNTNAP2", "CSMD1", "CSMD3", "CTNNA2", "CTNNA3",
                "DLG2", "DMD", "DPP10", "EYS", "GRID2", "LINGO1", "NRG1", "PCDH15",
                "PTPRD", "RASGEF1B", "RBFOX1", "SLC26A3", "TMEM132D")
  scaffold_genes <- c("DLG2", "DLG4", "HOMER1", "HOMER2", "HOMER3", "SHANK1", 
                      "SHANK2", "SHANK3")
  cam_genes <- c("CDH2", "CDH12", "LRFN1", "DAG1", "NLGN1", "NLGN2", "NLGN3",
                 "LRRTM1", "LRRTM2", "ADGRL2", "ADGRL3", "EPHB2", "EPHB3", 
                 "ADGRB1", "ADGRB3") #"NRXN1", "NRXN2", "NRXN3",
  ne_genes <- c("MSI1", "PLEKHG4B", "GNG4", "PEG10", "RND2", "APLP1", "SOX2",
                "TUBB2B", "CHGA", "CHGB", "NCAM1", "SYP")
  
  combo <- c(cam_genes, scaffold_genes)
   
  
  df_c3 <- data.frame(SYMBOL = c3_genes, Category = "C3")
  df_scaffold <- data.frame(SYMBOL = scaffold_genes, Category = "Scaffold")
  df_cam <- data.frame(SYMBOL = cam_genes, Category = "CAM")
  df_ne <- data.frame(SYMBOL = ne_genes, Category = "Neuroendocrine")
  metadata_row <- rbind(df_c3, df_scaffold, df_cam, df_ne) %>%
    dplyr::distinct(SYMBOL, .keep_all = TRUE)
  
  metadata_col <- cl_1750 %>% dplyr::filter(GSE %in% gse_id) %>%
    dplyr::select(Sample_ID, Consensus_class, Gender, Age, Stage)
  
  norm_counts <- final_merged_df %>% 
    dplyr::select(SYMBOL, metadata_col$Sample_ID) %>%
    dplyr::filter(SYMBOL %in% disp_genes) %>%
    tibble::column_to_rownames("SYMBOL")
  
  disp_genes <- intersect(c(c3_genes, ne_genes, cam_genes, scaffold_genes), final_merged_df$SYMBOL)
  
  # Heatmap overrides
  heatmap.override <- list(
    #force.log        = TRUE,                 # Force log transformation
    col.ann          = c("Consensus_class"),  # Column annotation
    row.ann          = c("Category"),         # Row annotation
    col.gaps         = c("Consensus_class"),                  # Column gaps
    row.gaps         = c("Category"),                  # Row gaps
    col.cluster      = c("Consensus_class"),  # Column clustering
    row.cluster      = c("Category"),         # Row clustering
    palette         = "rdbu",                # Heatmap palette
    ann.palette     = "discrete",             # Annotation palette
    border.color    = NA,                     # Cell border color
    show.expr.legend = TRUE,                  # Show expression legend
    title           = "",                     # Heatmap title
    format           = "tiff"                 # Output file format
  )
  
  parent.dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop"
  gmt.dir    <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/GSEA_genesets"
  scripts.dir <- NULL
  
  proj.params <- setup_project(
    proj               = "proj",
    species = "species",
    contrasts = "contrasts",
    parent.dir = parent.dir,
    gmt.dir = gmt.dir,
    scripts.dir      = scripts.dir,
    #deseq2.override  = deseq2.override,
    heatmap.override = heatmap.override,
    #volcano.override = volcano.override
  )
  
  proj.params$heatmap$title <- gse_id
  ph <- plot_heatmap(norm_counts, proj.params, metadata_col, metadata_row, disp_genes)
  heatmap_plot_list[[length(heatmap_plot_list) + 1]] <- ph$ph$gtable
}

heatmap_output_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop"
if (length(heatmap_plot_list) > 0) {
  pdf(file.path(heatmap_output_path, paste0("Heatmaps_Dr.Chao.pdf")),
      width = 11, height = 11)
  for (ht in heatmap_plot_list) {
    grid::grid.newpage()
    grid::grid.draw(ht)
  }
  dev.off()
}

## bulk RNASeq Xinyi

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


proj <- "Xinyi"
species <- "Homo sapiens"
contrasts <- c()

parent_dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Past"
gmt_dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/GSEA_genesets"

# DESeq2 overrides
deseq2.override <- list(
  contrasts     = contrasts,
  #design        = "Comparisons",            # DESeq2 design formula or column name
  #lfc.cutoff    = 0,                        # Log fold change cutoff for significance
  #padj.cutoff   = 0.1,                      # Adjusted p-value cutoff for significance
  batch.correct = FALSE                     # Boolean, whether to apply batch correction
)

# Heatmap overrides
heatmap.override <- list(
  #force.log        = TRUE,                  # Force log transformation
  col.ann          = NULL,                  # Column annotation
  #row.ann          = NULL,                  # Row annotation
  #col.gaps         = NULL,                  # Column gaps
  #row.gaps         = NULL,                  # Row gaps
  col.cluster      = "all",                 # Column clustering
  #row.cluster      = "all",                 # Row clustering
  #palette         = "rdbu",                # Heatmap palette
  #ann.palette     = "discrete",            # Annotation palette
  #border.color    = NA,                    # Cell border color
  #show.expr.legend = TRUE,                  # Show expression legend
  #title           = "",                    # Heatmap title
  format           = "tiff"                 # Output file format
)

# Volcano plot overrides
volcano.override <- list(
  #lfc.cutoff  = 0.58,                         # Minimum log2 fold-change to highlight
  #padj.cutoff = 0.05,                      # Adjusted p-value cutoff
  #color       = "vrds",                    # Color palette
  #label.genes = c()                         # Genes to label on the plot
)

salmon_dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Past/Xinyi/salmon"
proj_dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Past/Xinyi"
txi <- prep_txi(salmon_dir = salmon_dir, 
                species    = "Homo sapiens", 
                db_version = "113", 
                filename   = NULL, 
                output_dir = proj_dir)

metadata <- openxlsx::read.xlsx(file.path(proj_dir, "Xinyi_Metadata.xlsx"))
raw_counts_mat <- NULL
design      <- "Batch"
deseq2_data <- prepare_deseq2_input(expr_mat = raw_counts_mat,
                                    txi      = txi, 
                                    metadata = metadata,
                                    design   = design)

dds <- DESeqDataSetFromTximport(txi     = deseq2_data$txi,
                                colData = deseq2_data$metadata,
                                design  = ~ Batch)

dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

library(limma)
vsd_batch_corrected <- removeBatchEffect(vsd_mat,
                                         batch = deseq2_data$metadata$Batch)
ann_list <- get_annotations()
vsd_batch_corrected <- vsd_batch_corrected %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ID") %>%
  add_annotation(ann_list = ann_list) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarize(across(.cols = where(is.numeric), .fns = mean, na.rm = TRUE), .groups = "drop") %>%
  tibble::column_to_rownames("SYMBOL") %>%
  as.matrix()

cl <- getConsensusClass(x = vsd_batch_corrected, gene_id = c("hgnc_symbol"))
metadata_col <- dplyr::left_join(x = metadata %>% dplyr::mutate(Sample_ID = make.names(Sample_ID)), 
                                 y = cl %>% tibble::rownames_to_column("Sample_ID"),
                                 by = c("Sample_ID" = "Sample_ID"))

df_c3 <- data.frame(SYMBOL = c3_genes, Category = "C3")
df_scaffold <- data.frame(SYMBOL = scaffold_genes, Category = "Scaffold")
df_cam <- data.frame(SYMBOL = cam_genes, Category = "CAM")
df_ne <- data.frame(SYMBOL = ne_genes, Category = "Neuroendocrine")
metadata_row <- rbind(df_c3, df_scaffold, df_cam, df_ne) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

mat <- vsd_batch_corrected[rownames(vsd_batch_corrected) %in% c(c3_genes, scaffold_genes, cam_genes, ne_genes), ]

ph <- plot_heatmap(expr_mat            = mat, 
                   label_genes         =  metadata_row$SYMBOL,  
                   filename            = "Xinyi",
                   output_dir          = proj_dir,
                   metadata_col        = metadata_col, 
                   metadata_row        = metadata_row,
                   col_annotations     = c("consensusClass"),
                   row_annotations     = c("Category"),
                   col_gap_by          = c("consensusClass"),
                   row_gap_by          = c("Category"),
                   col_cluster_by      = c("consensusClass"),
                   row_cluster_by      = c("Category"),
                   plot_title          = NULL,
                   heatmap_palette     = "rdbu",
                   annotation_palette  = "discrete",
                   border_color        = NA,
                   force_log           = FALSE,
                   show_expr_legend    = TRUE,
                   save_plot           = TRUE,
                   save_matrix         = TRUE)

# Save the results
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Classification")
openxlsx::writeData(wb, sheet = "Classification", x = metadata_col)
openxlsx::saveWorkbook(wb, file = file.path(proj_dir,"Classification.xlsx"), overwrite = TRUE)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "VST counts")
openxlsx::writeData(wb, sheet = "VST counts", x = vsd_batch_corrected, rowNames = TRUE)
openxlsx::saveWorkbook(wb, file = file.path(proj_dir,"Batch corrected VST counts.xlsx"), overwrite = TRUE)



## scRNAseq

PROJ <- c("scRNASeq_BBN_C57BL6", "scRNASeq_BBN_Rag", "scRNASeq_GSE145216", 
          "scRNASeq_GSE164557", "scRNASeq_GSE200139", "scRNASeq_GSE217093",
          "scRNASeq_Jinfen", "scRNASeq_Chen", "scRNASeq_GSE135337", 
          "scRNASeq_GSE137829", "scRNASeq_GSE222315", "scRNASeq_HRA003620")

for (proj in PROJ){
  
  if (proj %in% c("scRNASeq_BBN_C57BL6", "scRNASeq_BBN_Rag", "scRNASeq_GSE145216", 
                  "scRNASeq_GSE164557", "scRNASeq_GSE200139", "scRNASeq_GSE217093",
                  "scRNASeq_Jinfen")){
    species <- "Mus musculus"
  } else if (proj %in% c("scRNASeq_Chen", "scRNASeq_GSE135337", "scRNASeq_GSE137829",
                         "scRNASeq_GSE222315", "scRNASeq_HRA003620")){
    species <- "Homo sapiens"
  }
  
  parent_dir  <- "/hpc/home/kailasamms/scratch"
  gmt_dir     <- "/hpc/home/kailasamms/projects/GSEA_genesets"
  scripts_dir <- "/hpc/home/kailasamms/projects/scRNASeq"
  
  contrasts <- c()
  deseq2.override <- list()
  heatmap.override <- list()
  volcano.override <- list()
  
  proj.params <- setup_project(proj             = proj,
                               species          = species,  #"Mus musculus", "Homo sapiens"
                               contrasts        = contrasts,
                               parent_dir       = parent_dir,
                               gmt_dir          = gmt_dir,
                               scripts_dir      = scripts_dir,
                               deseq2.override  = deseq2.override,
                               heatmap.override = heatmap.override,
                               volcano.override = volcano.override)

  integrated_seurat <- readRDS(file.path(proj.params$seurat_dir, "integrated_seurat_ann.rds"))
  
  # Initialize an empty list to store all cell type signatures
  signatures_list <- list()
  
  scaffold_genes <- c("DLG2", "DLG4", "HOMER1", "HOMER2", "HOMER3", "SHANK1", 
                      "SHANK2", "SHANK3")
  cam_genes <- c("CDH2", "CDH12", "LRFN1", "DAG1", "NLGN1", "NLGN2", "NLGN3",
                 "LRRTM1", "LRRTM2", "ADGRL2", "ADGRL3", "EPHB2", "EPHB3", 
                 "ADGRB1", "ADGRB3") #"NRXN1", "NRXN2", "NRXN3",
  combo <- c(cam_genes, scaffold_genes)
  
  for (feature_subset in c("cam_genes", "scaffold_genes", "combo")){
    
    # Define name of modules
    module_name <- make.names(feature_subset )
    
    # Determine features from the marker_df
    xlsx_features <- get(feature_subset)
    
    # Determine features present in data set
    present_features <- rownames(SeuratObject::GetAssayData(object = integrated_seurat, 
                                                            assay = "RNA", 
                                                            layer = "data"))
    
    # Match features (case-insensitive) and filter to only genes present in the data
    features <- present_features[base::tolower(present_features) %in% base::tolower(xlsx_features)]
    
    # Only add the validated and sorted feature vector to the master list if it contains 2 or more genes
    if (length(features) >= 2) {
      signatures_list[[module_name]] <- sort(features)
    }
  }
  
  # Calculate module scores using seurat
  integrated_seurat <- Seurat::AddModuleScore(object = integrated_seurat,
                                              features = signatures_list,
                                              assay = "RNA",
                                              slot = "data",
                                              name = names(signatures_list))  # "_UCell"
  
  saveRDS(integrated_seurat, file = file.path("/hpc/home/kailasamms", paste0(proj, "_Chao.rds")))
  
  # Plot module scores
  DefaultAssay(integrated_seurat) <- "RNA"
  p1 <- Seurat::FeaturePlot(object = integrated_seurat,
                            features = paste0(names(signatures_list), seq(1:length(signatures_list))),
                            cols =  c("#440154FF", viridis(n = 10, option = "C", direction = 1)),
                            pt.size = 0.4,
                            order = TRUE,
                            min.cutoff = 'q10',
                            reduction = "umap_harmony",
                            raster = FALSE)
  
  # Save plot
  ggplot2::ggsave(filename = paste0("Module_score_plot_", proj, ".tiff"),
                  path = "/hpc/home/kailasamms/",
                  plot = p1,
                  device = "jpeg",
                  width = 20,
                  height = 20,
                  units = c("in"),
                  bg = "white")
  
  plot_umap(integrated_seurat = integrated_seurat, 
            reduction = "umap_harmony",  
            color_col = "CellType",         
            filename = paste0("UMAP_", proj),   
            output_dir = "/hpc/home/kailasamms/",
            split_col = NULL)
  
}







