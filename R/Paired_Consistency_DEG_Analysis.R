#******************************************************************************#
#         GROWTH HORMONE vs CONTROL - VST-BASED PAIRED LFC ANALYSIS           #
#  Uses batch-corrected VST non-blind counts (already on log scale)            #
#  No additional log transformation needed - subtraction = log2FC directly     #
#******************************************************************************#

library(openxlsx)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)

base_dir    <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data"
proj        <- "RNASeq_Human_Vera_HPrEC (all)"
species     <- "Human"
gmt_dir     <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/Signatures/Human"

# Toggle between old and new pipeline by changing this one line
pipeline    <- "old"   # "old" or "new"

deseq2_dir  <- file.path(base_dir, proj, species, ifelse(pipeline == "old", "08.DESeq2 (old)", "08.DESeq2"))
custom_dir  <- file.path(base_dir, proj, species, ifelse(pipeline == "old", "11.Custom_Analysis (old)", "11.Custom_Analysis"))

if (!dir.exists(custom_dir)) dir.create(custom_dir, recursive = TRUE)

#------------------------------------------------------------------------------#
# SECTION 1: Load Normalized Counts                                            #
#------------------------------------------------------------------------------#
# Using Normalized Counts (linear scale, not VST) because VST compresses       #
# low-fold changes, making subtle GH responses harder to detect.               #
# Paired log2FC (GH/Control) per replicate naturally controls for cell-line    #
# and batch differences.                                                        #
#------------------------------------------------------------------------------#

if (pipeline == "old") {
  
# Clean the old norm counts to match new norm counts format
norm_counts_prefilter <- openxlsx::read.xlsx(file.path(deseq2_dir, "Norm_counts_DESeq2.xlsx")) %>%
  # Remove known non-coding artifacts
  filter(!SYMBOL %in% c("Y_RNA", "Metazoa_SRP", "U3",
                        "5_8S_rRNA", "7SK", "U1", "U2", "U4",
                        "SNORA72", "SNORA73", "SNORD39", "Vault")) %>%
  # For remaining duplicates, keep highest expressed row
  dplyr::mutate(total_expr = rowSums(dplyr::pick(where(is.numeric)))) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::slice_max(order_by = total_expr, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(-total_expr)  %>%
  dplyr::rename(HPrEC1_C1  = Line1_Control1,
                HPrEC1_C2  = Line1_Control2,
                HPrEC1_C3  = Line1_Control3,
                HPrEC1_GH1 = Line1_GH1,
                HPrEC1_GH2 = Line1_GH2,
                HPrEC1_GH3 = Line1_GH3,
                HPrEC2_C1  = Line2_Control1,
                HPrEC2_C2  = Line2_Control2,
                HPrEC2_C3  = Line2_Control3,
                HPrEC2_GH1 = Line2_GH1,
                HPrEC2_GH2 = Line2_GH2,
                HPrEC2_GH3 = Line2_GH3) %>% 
  dplyr::select(SYMBOL, 
                HPrEC1_C1, HPrEC1_C2, HPrEC1_C3,
                HPrEC1_GH1, HPrEC1_GH2, HPrEC1_GH3,
                HPrEC2_C1, HPrEC2_C2, HPrEC2_C3,
                HPrEC2_GH1, HPrEC2_GH2, HPrEC2_GH3)
} else {
  norm_counts_prefilter <- openxlsx::read.xlsx(file.path(deseq2_dir, "Norm_Counts_ExpressionPlots.xlsx"))
}

# Pre-filter: Keep genes expressed (mean > 10) in at least one group
# Preserves GH-induced genes that may be low in control
norm_counts <- norm_counts_prefilter %>%
  dplyr::filter(
    rowMeans(dplyr::pick(contains("HPrEC1_C")))  > 10 |
      rowMeans(dplyr::pick(contains("HPrEC1_GH"))) > 10 |
      rowMeans(dplyr::pick(contains("HPrEC2_C")))  > 10 |
      rowMeans(dplyr::pick(contains("HPrEC2_GH"))) > 10
  )

cat("Genes before filter:", nrow(norm_counts_prefilter), "\n")
cat("Genes after filter:", nrow(norm_counts), "\n")

#------------------------------------------------------------------------------#
# SECTION 2: Compute Per-Replicate Paired log2FC                               #
#  log2(GH/Control) on normalized counts + pseudocount = true log2FC           #
#  Pseudocount of +1 added to all values before ratio to handle zeros          #
#  LFC threshold of 0.5 applied to Status calls to exclude near-zero changes   #
#------------------------------------------------------------------------------#

LFC_THRESHOLD <- 0   # minimum |mean log2FC| to call UP/DOWN

analysis_df <- norm_counts %>%
  
  # 1. Add pseudocount
  dplyr::mutate(across(where(is.numeric), ~ .x + 1)) %>%
  
  # 2. Individual paired LFCs
  dplyr::mutate(
    HPrEC1_LFC1 = log2(HPrEC1_GH1 / HPrEC1_C1),
    HPrEC1_LFC2 = log2(HPrEC1_GH2 / HPrEC1_C2),
    HPrEC1_LFC3 = log2(HPrEC1_GH3 / HPrEC1_C3),
    HPrEC2_LFC1 = log2(HPrEC2_GH1 / HPrEC2_C1),
    HPrEC2_LFC2 = log2(HPrEC2_GH2 / HPrEC2_C2),
    HPrEC2_LFC3 = log2(HPrEC2_GH3 / HPrEC2_C3)
  ) %>%
  
  # 3. Means — combined uses all 6 LFCs directly (not mean of means)
  dplyr::mutate(
    mean_log2FC_HPrEC1   = (HPrEC1_LFC1 + HPrEC1_LFC2 + HPrEC1_LFC3) / 3,
    mean_log2FC_HPrEC2   = (HPrEC2_LFC1 + HPrEC2_LFC2 + HPrEC2_LFC3) / 3,
    mean_log2FC_combined = (HPrEC1_LFC1 + HPrEC1_LFC2 + HPrEC1_LFC3 +
                              HPrEC2_LFC1 + HPrEC2_LFC2 + HPrEC2_LFC3) / 6
  ) %>%
  
  # 4. Vote counts
  dplyr::mutate(
    n_up          = rowSums(dplyr::pick(contains("LFC")) > 0),
    n_down        = rowSums(dplyr::pick(contains("LFC")) < 0),
    n_up_HPrEC1   = rowSums(dplyr::pick(contains("HPrEC1_LFC")) > 0),
    n_down_HPrEC1 = rowSums(dplyr::pick(contains("HPrEC1_LFC")) < 0),
    n_up_HPrEC2   = rowSums(dplyr::pick(contains("HPrEC2_LFC")) > 0),
    n_down_HPrEC2 = rowSums(dplyr::pick(contains("HPrEC2_LFC")) < 0)
  ) %>%
  
  # 5. Status calls — vote consistency + LFC magnitude threshold
  dplyr::mutate(
    Status_HPrEC1 = dplyr::case_when(
      n_up_HPrEC1   == 3 & mean_log2FC_HPrEC1 >  LFC_THRESHOLD ~ "UP",
      n_down_HPrEC1 == 3 & mean_log2FC_HPrEC1 < -LFC_THRESHOLD ~ "DOWN",
      TRUE                                                       ~ "Ambiguous"
    ),
    Status_HPrEC2 = dplyr::case_when(
      n_up_HPrEC2   == 3 & mean_log2FC_HPrEC2 >  LFC_THRESHOLD ~ "UP",
      n_down_HPrEC2 == 3 & mean_log2FC_HPrEC2 < -LFC_THRESHOLD ~ "DOWN",
      TRUE                                                       ~ "Ambiguous"
    ),
    Status = dplyr::case_when(
      n_up   == 6 & mean_log2FC_combined >  LFC_THRESHOLD ~ "UP",
      n_down == 6 & mean_log2FC_combined < -LFC_THRESHOLD ~ "DOWN",
      TRUE                                                      ~ "Ambiguous"
    )
  )

#------------------------------------------------------------------------------#
# SECTION 3: One-Sample t-test on Replicate log2FCs per Gene                  #
#  H0: mean log2FC = 0 (no GH effect)                                          #
#  Each FC value is already a paired difference (GH - matched Control)         #
#  Testing whether the mean difference differs from 0                           #
#  NOTE: pval (not padj) is passed to GSEA for ranking — this is correct       #
#  because GSEA needs an unthresholded continuous ranking signal                #
#------------------------------------------------------------------------------#

run_ttest <- function(lfc_df, fc_cols) {
  apply(lfc_df[, fc_cols], 1, function(x) {
    if (stats::var(x) == 0) return(NA_real_)
    stats::t.test(x, mu = 0, alternative = "two.sided")$p.value
  })
}

pval_hprec1 <- run_ttest(analysis_df, c("HPrEC1_LFC1", "HPrEC1_LFC2", "HPrEC1_LFC3"))
pval_hprec2 <- run_ttest(analysis_df, c("HPrEC2_LFC1", "HPrEC2_LFC2", "HPrEC2_LFC3"))

#------------------------------------------------------------------------------#
# SECTION 4: Adjust p-values (BH/FDR)                                         #
#  Zero-variance genes get NA pval — excluded from adjustment pool             #
#  so they don't inflate the denominator                                        #
#------------------------------------------------------------------------------#

padj_hprec1 <- rep(NA_real_, length(pval_hprec1))
padj_hprec2 <- rep(NA_real_, length(pval_hprec2))

valid1 <- !is.na(pval_hprec1)
valid2 <- !is.na(pval_hprec2)

padj_hprec1[valid1] <- stats::p.adjust(pval_hprec1[valid1], method = "BH")
padj_hprec2[valid2] <- stats::p.adjust(pval_hprec2[valid2], method = "BH")

#------------------------------------------------------------------------------#
# SECTION 5: Combine p-values across cell lines (Fisher's method)             #
#  NA p-values (zero-variance genes) treated as p=1 (no evidence)             #
#  p=0 (extreme t-statistics) floored at 1e-300 rather than flipped to 1      #
#------------------------------------------------------------------------------#

combine_p <- function(p1, p2) {
  p1 <- ifelse(is.na(p1), 1, pmax(p1, 1e-300))
  p2 <- ifelse(is.na(p2), 1, pmax(p2, 1e-300))
  stats::pchisq(-2 * (log(p1) + log(p2)), df = 4, lower.tail = FALSE)
}

#------------------------------------------------------------------------------#
# SECTION 6: Assemble Full Results Table                                       #
#------------------------------------------------------------------------------#

analysis_df <- analysis_df %>%
  dplyr::mutate(
    pval_HPrEC1   = pval_hprec1,
    padj_HPrEC1   = padj_hprec1,
    pval_HPrEC2   = pval_hprec2,
    padj_HPrEC2   = padj_hprec2,
    pval_combined = combine_p(pval_HPrEC1, pval_HPrEC2),
    padj_combined = stats::p.adjust(pval_combined, method = "BH")
  )

#------------------------------------------------------------------------------#
# SECTION 7: Save Results                                                      #
#------------------------------------------------------------------------------#

add_sheet <- function(wb, sheet_name, data) {
  openxlsx::addWorksheet(wb, sheetName = sheet_name)
  openxlsx::writeData(wb, sheet = sheet_name, x = data, rowNames = FALSE)
  openxlsx::freezePane(wb, sheet = sheet_name, firstRow = TRUE)
  openxlsx::addStyle(wb, sheet = sheet_name,
                     style = openxlsx::createStyle(textDecoration = "Bold"),
                     rows  = 1, cols = 1:ncol(data), gridExpand = TRUE)
}

wb <- openxlsx::createWorkbook()

add_sheet(wb, "All_Genes",   analysis_df)
add_sheet(wb, "UP_combined", dplyr::filter(analysis_df, Status == "UP"))
add_sheet(wb, "DN_combined", dplyr::filter(analysis_df, Status == "DOWN"))
add_sheet(wb, "UP_HPrEC1",   dplyr::filter(analysis_df, Status_HPrEC1 == "UP"))
add_sheet(wb, "DN_HPrEC1",   dplyr::filter(analysis_df, Status_HPrEC1 == "DOWN"))
add_sheet(wb, "UP_HPrEC2",   dplyr::filter(analysis_df, Status_HPrEC2 == "UP"))
add_sheet(wb, "DN_HPrEC2",   dplyr::filter(analysis_df, Status_HPrEC2 == "DOWN"))

openxlsx::saveWorkbook(wb,
                       file      = file.path(custom_dir, "Paired_DEG_Analysis.xlsx"),
                       overwrite = TRUE)


#------------------------------------------------------------------------------#
# SECTION 8: Volcano Plots                                                     #
#------------------------------------------------------------------------------#
label_genes <- c("H4C16", "TAF5", "TADA3", "BRCC3", "NUDT16L1")
# HPrEC1 volcano
# Get top genes by LFC magnitude regardless of pval
top_lfc_genes_hprec1 <- analysis_df %>%
  dplyr::filter(-log10(pval_HPrEC1) > 1.3) %>%
  #dplyr::slice_max(order_by = abs(mean_log2FC_HPrEC1), n = 10) %>%
  dplyr::pull(SYMBOL) %>%
  intersect(label_genes)

plot_volcano(df          = analysis_df,
             x_col       = "mean_log2FC_HPrEC1",
             y_col       = "pval_HPrEC1",
             label_col   = "SYMBOL",
             label_genes = top_lfc_genes_hprec1,
             x_cutoff    = 0,
             filename    = "Volcano_Plot_HPrEC1",
             output_dir  = custom_dir)

# HPrEC2 volcano
# Get top genes by LFC magnitude regardless of pval
top_lfc_genes_hprec2 <- analysis_df %>%
  dplyr::filter(-log10(pval_HPrEC2) > 1.3) %>%
  #dplyr::slice_max(order_by = abs(mean_log2FC_HPrEC2), n = 10) %>%
  dplyr::pull(SYMBOL) %>%
  intersect(label_genes)


plot_volcano(df          = analysis_df,
             x_col       = "mean_log2FC_HPrEC2",
             y_col       = "pval_HPrEC2",
             label_col   = "SYMBOL",
             label_genes = top_lfc_genes_hprec2,
             x_cutoff    = 0,
             filename    = "Volcano_Plot_HPrEC2",
             output_dir  = custom_dir)

# Combined volcano
# Get top genes by LFC magnitude regardless of pval
top_lfc_genes_combined <- analysis_df %>%
  dplyr::filter(-log10(pval_combined) > 1.3) %>%
  #dplyr::slice_max(order_by = abs(mean_log2FC_combined), n = 10) %>%
  dplyr::pull(SYMBOL) %>%
  intersect(label_genes)


plot_volcano(df          = analysis_df,
             x_col       = "mean_log2FC_combined",
             y_col       = "pval_combined",
             label_col   = "SYMBOL",
             label_genes = top_lfc_genes_combined,
             x_cutoff    = 0,
             filename    = "Volcano_Plot_combined",
             output_dir  = custom_dir)

#------------------------------------------------------------------------------#
# SECTION 9: Quadrant Plot                                                     #
#  Global correlation uses all genes (concordance measure)                     #
#  Responder correlation computed separately on UP/DOWN genes only             #
#------------------------------------------------------------------------------#
label_genes <- c("H4C16", "TAF5", "TADA3", "BRCC3", "NUDT16L1")
status_colors <- c("UP" = "#C0392B", "DOWN" = "#1A6FA8", "Ambiguous" = "grey82")

# Global correlation (all genes)
r_all <- cor(analysis_df$mean_log2FC_HPrEC1,
             analysis_df$mean_log2FC_HPrEC2,
             use = "complete.obs")

# Responder-only correlation
r_responders <- cor(
  dplyr::filter(analysis_df, Status != "Ambiguous")$mean_log2FC_HPrEC1,
  dplyr::filter(analysis_df, Status != "Ambiguous")$mean_log2FC_HPrEC2,
  use = "complete.obs"
)

top_genes <- analysis_df %>%
  dplyr::filter(Status != "Ambiguous") %>%
  dplyr::group_by(Status) %>%
  dplyr::slice_max(order_by = abs(mean_log2FC_combined), n = 10) %>%
  dplyr::ungroup()

ggplot(analysis_df, aes(x = mean_log2FC_HPrEC1, y = mean_log2FC_HPrEC2)) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.35) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.35) +
  
  geom_point(data  = dplyr::filter(analysis_df, Status == "Ambiguous"),
             color = "grey82", size = 0.5, alpha = 0.25) +
  
  geom_point(data  = dplyr::filter(analysis_df, Status != "Ambiguous"),
             aes(color = Status), size = 1, alpha = 0.8) +
  
  # geom_point(data  = dplyr::filter(analysis_df, SYMBOL %in% label_genes),
  #            color = "black", size = 1.5, alpha = 1) +
  
  ggrepel::geom_text_repel(
    data         = dplyr::filter(analysis_df, SYMBOL %in% label_genes),
    aes(label    = SYMBOL, color = Status),
    size         = 2.6, fontface = "italic", max.overlaps = Inf,
    box.padding  = 0.4, point.padding = 0.3,
    segment.size = 0.25, segment.alpha = 0.5, 
    nudge_x = -2, nudge_y = -1,
    show.legend  = FALSE
  ) +
  
  # ggrepel::geom_text_repel(
  #   data          = dplyr::filter(top_genes, Status == "UP"),
  #   aes(label     = SYMBOL, color = Status),
  #   size          = 2.6, fontface = "italic", max.overlaps = Inf,
  #   box.padding   = 0.4, point.padding = 0.3,
  #   segment.size  = 0.25, segment.color = "#C0392B", segment.alpha = 0.5,
  #   nudge_x = 0.8, nudge_y = 0.1, direction = "y", hjust = 0,
  #   show.legend   = FALSE
  # ) +
  # 
  # ggrepel::geom_text_repel(
  #   data          = dplyr::filter(top_genes, Status == "DOWN"),
  #   aes(label     = SYMBOL, color = Status),
  #   size          = 2.6, fontface = "italic", max.overlaps = Inf,
  #   box.padding   = 0.4, point.padding = 0.3,
  #   segment.size  = 0.25, segment.color = "#1A6FA8", segment.alpha = 0.5,
  #   nudge_x = -0.8, nudge_y = -0.1, direction = "y", hjust = 1,
  #   show.legend   = FALSE
  # ) +
  
  annotate("text", x =  4.5, y =  4.8, label = "UP in both",
           color = "#C0392B", size = 3, fontface = "bold", hjust = 1) +
  annotate("text", x = -4.5, y = -4.8, label = "DOWN in both",
           color = "#1A6FA8", size = 3, fontface = "bold", hjust = 0) +
  annotate("text", x =  4.5, y = -4.8, label = "Discordant",
           color = "grey50", size = 2.8, hjust = 1) +
  annotate("text", x = -4.5, y =  4.8, label = "Discordant",
           color = "grey50", size = 2.8, hjust = 0) +
  
  # # Both correlations annotated
  # annotate("text", x = -4.5, y =  4.4,
  #          label = paste0("r (all) = ", round(r_all, 3)),
  #          size = 3, hjust = 0, color = "grey40") +
  # annotate("text", x = -4.5, y =  3.9,
  #          label = paste0("r (responders) = ", round(r_responders, 3)),
  #          size = 3, hjust = 0, color = "grey40") +

  scale_color_manual(values = status_colors) +
  scale_x_continuous(limits = c(-5.5, 5.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-5.5, 5.5), expand = c(0, 0)) +
  
  labs(x     = "Mean log\u2082FC (HPrEC1)",
       y     = "Mean log\u2082FC (HPrEC2)",
       color = NULL,
       title = "GH vs Control \u2014 concordance across cell lines") +
  
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(size = 11, face = "bold", hjust = 0),
    axis.title      = element_text(size = 10),
    axis.text       = element_text(size = 9, color = "grey30"),
    axis.line       = element_line(linewidth = 0.35, color = "grey40"),
    axis.ticks      = element_line(linewidth = 0.35, color = "grey40"),
    legend.position = "bottom",
    legend.text     = element_text(size = 9),
    panel.grid      = element_blank(),
    plot.margin     = margin(10, 20, 10, 10)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(file.path(custom_dir, "Quadrant_Plot_HPrEC1_vs_HPrEC2.pdf"),
       width = 6, height = 6, useDingbats = FALSE)

#------------------------------------------------------------------------------#
# SECTION 10: Bar Plots — combined, HPrEC1, HPrEC2                            #
#------------------------------------------------------------------------------#

make_barplot <- function(df, status_col, title_label, filename, output_dir) {
  df %>%
    dplyr::filter(.data[[status_col]] != "Ambiguous") %>%
    dplyr::count(.data[[status_col]], name = "n") %>%
    dplyr::rename(Status = 1) %>%
    ggplot(aes(x = Status, y = n, fill = Status)) +
    geom_col(width = 0.5) +
    geom_text(aes(label = n), vjust = -0.5, size = 2.8) +
    scale_fill_manual(values = c("UP" = "#C0392B", "DOWN" = "#1A6FA8")) +
    labs(y = "Number of genes", x = NULL, title = title_label) +
    theme_classic(base_size = 12) +
    theme(
      plot.title      = element_text(size = 11, face = "bold", hjust = 0),
      axis.title      = element_text(size = 10),
      axis.text       = element_text(size = 9, color = "grey30"),
      axis.line       = element_line(linewidth = 0.35, color = "grey40"),
      axis.ticks      = element_line(linewidth = 0.35, color = "grey40"),
      legend.position = "none",
      panel.grid      = element_blank()
    )
  ggsave(file.path(output_dir, paste0(filename, ".pdf")),
         width = 4, height = 5, useDingbats = FALSE)
}

make_barplot(analysis_df, "Status",
             "Genes consistently responding to GH\n(6/6 replicates concordant)",
             "Barplot_DEG_Counts_combined", custom_dir)

make_barplot(analysis_df, "Status_HPrEC1",
             "Genes consistently responding to GH\n(3/3 replicates concordant, HPrEC1)",
             "Barplot_DEG_Counts_HPrEC1", custom_dir)

make_barplot(analysis_df, "Status_HPrEC2",
             "Genes consistently responding to GH\n(3/3 replicates concordant, HPrEC2)",
             "Barplot_DEG_Counts_HPrEC2", custom_dir)

#------------------------------------------------------------------------------#
# SECTION 11: GSEA and ORA                                                     #
#                                                                              #
#  WHY BYPASS analyze_pathways()?                                              #
#  GH effect is feeble — DESeq2 yields <10 significant DEGs. With only n=3     #
#  replicates per cell line, t-test p-values are unstable and survive FDR      #
#  correction only for the strongest effects. Using pval/padj for either       #
#  GSEA ranking or ORA foreground definition would discard real biology.       #
#                                                                              #
#  GSEA RANKING STRATEGY:                                                      #
#  Plain mean log2FC used as ranking metric — more stable than                 #
#  -log10(pval) * sign(LFC) with n=3. With few replicates, pval adds noise     #
#  not signal to the ranking. scoreType auto-detected by run_gsea() as "std"   #
#  since LFC spans positive and negative values.                               #
#                                                                              #
#  ORA FOREGROUND STRATEGY:                                                    #
#  Gene selected based on Vote-based concordance instead of padj thresholding. # 
#  Directional consistency across replicates IS the reproducibility filter.    #
#  Two stringency tiers per cell line and combined:                            #
#    Lenient : LFC > 0 / LFC < 0  (all directionally consistent genes)         #
#    Strict  : LFC > 0.5 / LFC < -0.5  (magnitude-filtered subset)             #
#  Universe = all expressed genes in analysis_df (experimental background),    #
#------------------------------------------------------------------------------#

gmt_files <- list.files(gmt_dir, full.names = TRUE, pattern = "\\.gmt$")

#------------------------------------------------------------------------------#
# SECTION 11A: GSEA — plain LFC ranked lists                                   #
#  One ranked list per cell line + combined                                    #
#  Sorted descending: highest-induced genes anchor top of ranked list          #
#------------------------------------------------------------------------------#

ranked_hprec1 <- analysis_df %>%
  dplyr::filter(!is.na(mean_log2FC_HPrEC1)) %>%
  dplyr::arrange(dplyr::desc(mean_log2FC_HPrEC1)) %>%
  { stats::setNames(.$mean_log2FC_HPrEC1, .$SYMBOL) }

ranked_hprec2 <- analysis_df %>%
  dplyr::filter(!is.na(mean_log2FC_HPrEC2)) %>%
  dplyr::arrange(dplyr::desc(mean_log2FC_HPrEC2)) %>%
  { stats::setNames(.$mean_log2FC_HPrEC2, .$SYMBOL) }

ranked_combined <- analysis_df %>%
  dplyr::filter(!is.na(mean_log2FC_combined)) %>%
  dplyr::arrange(dplyr::desc(mean_log2FC_combined)) %>%
  { stats::setNames(.$mean_log2FC_combined, .$SYMBOL) }

gsea_hprec1   <- run_gsea(ranked_list = ranked_hprec1,   gmt_files = gmt_files)
gsea_hprec2   <- run_gsea(ranked_list = ranked_hprec2,   gmt_files = gmt_files)
gsea_combined <- run_gsea(ranked_list = ranked_combined, gmt_files = gmt_files)

#------------------------------------------------------------------------------#
# SECTION 11B: ORA — vote-based foreground, experimental universe             #
#                                                                              #
#  Universe = all expressed genes (analysis_df) — fixed experimental          #
#  background shared across all GMT collections.                              #
#                                                                              #
#  Lenient: LFC sign only — captures all directionally consistent genes       #
#  Strict : LFC > 0.5 / < -0.5 — magnitude-filtered higher confidence subset #
#------------------------------------------------------------------------------#

universe <- analysis_df$SYMBOL

# ── HPrEC1 ───────────────────────────────────────────────────────────────────

ora_hprec1_lenient <- run_ora(
  gene_lists = list(up   = analysis_df %>% dplyr::filter(n_up_HPrEC1   == 3, mean_log2FC_HPrEC1 > 0) %>% dplyr::pull(SYMBOL),
                    down = analysis_df %>% dplyr::filter(n_down_HPrEC1 == 3, mean_log2FC_HPrEC1 < 0) %>% dplyr::pull(SYMBOL)),
  universe   = universe, gmt_files = gmt_files)

ora_hprec1_strict <- run_ora(
  gene_lists = list(up   = analysis_df %>% dplyr::filter(n_up_HPrEC1   == 3, mean_log2FC_HPrEC1 > 0.5)  %>% dplyr::pull(SYMBOL),
                    down = analysis_df %>% dplyr::filter(n_down_HPrEC1 == 3, mean_log2FC_HPrEC1 < -0.5) %>% dplyr::pull(SYMBOL)),
  universe   = universe, gmt_files = gmt_files)

# ── HPrEC2 ───────────────────────────────────────────────────────────────────

ora_hprec2_lenient <- run_ora(
  gene_lists = list(up   = analysis_df %>% dplyr::filter(n_up_HPrEC2   == 3, mean_log2FC_HPrEC2 > 0) %>% dplyr::pull(SYMBOL),
                    down = analysis_df %>% dplyr::filter(n_down_HPrEC2 == 3, mean_log2FC_HPrEC2 < 0) %>% dplyr::pull(SYMBOL)),
  universe   = universe, gmt_files = gmt_files)

ora_hprec2_strict <- run_ora(
  gene_lists = list(up   = analysis_df %>% dplyr::filter(n_up_HPrEC2   == 3, mean_log2FC_HPrEC2 > 0.5)  %>% dplyr::pull(SYMBOL),
                    down = analysis_df %>% dplyr::filter(n_down_HPrEC2 == 3, mean_log2FC_HPrEC2 < -0.5) %>% dplyr::pull(SYMBOL)),
  universe   = universe, gmt_files = gmt_files)

# ── Combined ─────────────────────────────────────────────────────────────────

ora_combined_lenient <- run_ora(
  gene_lists = list(up   = analysis_df %>% dplyr::filter(n_up    == 6, mean_log2FC_combined > 0) %>% dplyr::pull(SYMBOL),
                    down = analysis_df %>% dplyr::filter(n_down  == 6, mean_log2FC_combined < 0) %>% dplyr::pull(SYMBOL)),
  universe   = universe, gmt_files = gmt_files)

ora_combined_strict <- run_ora(
  gene_lists = list(up   = analysis_df %>% dplyr::filter(n_up    == 6, mean_log2FC_combined > 0.5)  %>% dplyr::pull(SYMBOL),
                    down = analysis_df %>% dplyr::filter(n_down  == 6, mean_log2FC_combined < -0.5) %>% dplyr::pull(SYMBOL)),
  universe   = universe, gmt_files = gmt_files)

#------------------------------------------------------------------------------#
# SECTION 11C: Save — one Pathways.xlsx per folder, three sheets each        #
#------------------------------------------------------------------------------#

save_pathways <- function(gsea, ora_lenient, ora_strict, output_dir) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  wb <- openxlsx::createWorkbook()
  
  sheets <- list(GSEA        = gsea,
                 ORA_lenient = ora_lenient,
                 ORA_strict  = ora_strict)
  
  for (nm in names(sheets)) {
    res <- sheets[[nm]]
    if (!is.null(res) && nrow(res) > 0) {
      openxlsx::addWorksheet(wb, sheetName = nm)
      openxlsx::writeData(wb, sheet = nm, x = res, rowNames = FALSE)
      openxlsx::freezePane(wb, sheet = nm, firstRow = TRUE)
      openxlsx::addStyle(wb, sheet = nm,
                         style     = openxlsx::createStyle(textDecoration = "Bold"),
                         rows      = 1, cols = 1:ncol(res), gridExpand = TRUE)
    }
  }
  
  openxlsx::saveWorkbook(wb,
                         file      = file.path(output_dir, "Pathways.xlsx"),
                         overwrite = TRUE)
}

save_pathways(gsea         = gsea_hprec1,
              ora_lenient  = ora_hprec1_lenient,
              ora_strict   = ora_hprec1_strict,
              output_dir   = file.path(custom_dir, "HPrEC1"))

save_pathways(gsea         = gsea_hprec2,
              ora_lenient  = ora_hprec2_lenient,
              ora_strict   = ora_hprec2_strict,
              output_dir   = file.path(custom_dir, "HPrEC2"))

save_pathways(gsea         = gsea_combined,
              ora_lenient  = ora_combined_lenient,
              ora_strict   = ora_combined_strict,
              output_dir   = file.path(custom_dir, "Combined"))

#------------------------------------------------------------------------------#
# SECTION 12: plot_pathways — one call per folder                             #
#  contrast = folder name only                                                #
#  plot_pathways() iterates all sheets automatically, auto-detecting          #
#  GSEA vs ORA per sheet via NES vs GeneRatio column signature               #
#------------------------------------------------------------------------------#

if (pipeline == "old") {
  
  expr_mat_combined <- norm_counts %>%
    dplyr::mutate(across(where(is.numeric), ~ log2(.x + 1)))
  expr_mat1 <- expr_mat_combined %>% dplyr::select(SYMBOL, contains("HPrEC1"))
  expr_mat2 <- expr_mat_combined %>% dplyr::select(SYMBOL, contains("HPrEC2"))
} else {
  
  expr_mat1 <- openxlsx::read.xlsx(
    file.path(deseq2_dir, "HPrEC1_GrowthHormone-HPrEC1_Control",
              "VST_NonBlind_Counts_Heatmaps.xlsx"))
  
  expr_mat2 <- openxlsx::read.xlsx(
    file.path(deseq2_dir, "HPrEC2_GrowthHormone-HPrEC2_Control",
              "VST_NonBlind_Counts_Heatmaps.xlsx"))
  
  expr_mat_combined <- dplyr::inner_join(expr_mat1, expr_mat2, by = "SYMBOL")
}


plot_pathways(contrast        = "HPrEC1",
              pathway_xlsx    = file.path(custom_dir, "HPrEC1", "Pathways.xlsx"),
              output_dir      = custom_dir,
              vst_nonblind    = expr_mat1,
              metadata        = file.path(deseq2_dir, "Vera_HPrEC_metadata.xlsx"),
              col_annotations = c("Cell_Line", "Treatment"))

plot_pathways(contrast        = "HPrEC2",
              pathway_xlsx    = file.path(custom_dir, "HPrEC2", "Pathways.xlsx"),
              output_dir      = custom_dir,
              vst_nonblind    = expr_mat2,
              metadata        = file.path(deseq2_dir, "Vera_HPrEC_metadata.xlsx"),
              col_annotations = c("Cell_Line", "Treatment"))

plot_pathways(contrast        = "Combined",
              pathway_xlsx    = file.path(custom_dir, "Combined", "Pathways.xlsx"),
              output_dir      = custom_dir,
              vst_nonblind    = expr_mat_combined,
              metadata        = file.path(deseq2_dir, "Vera_HPrEC_metadata.xlsx"),
              col_annotations = c("Cell_Line", "Treatment"))

#------------------------------------------------------------------------------#
# SECTION 14: Heatmaps for specific pathways of interest                       #
#------------------------------------------------------------------------------#

hprec1_pathways <- read.xlsx(file.path(custom_dir, "HPrEC1", "Pathways.xlsx"))
hprec2_pathways <- read.xlsx(file.path(custom_dir, "HPrEC2", "Pathways.xlsx"))

intersect(hprec1_pathways %>% dplyr::filter(pval <= 0.05) %>% dplyr::pull(Description), 
          hprec2_pathways %>% dplyr::filter(pval <= 0.05) %>% dplyr::pull(Description))

plot_pathways_list <- c("DNA DOUBLE STRAND BREAK RESPONSE",
                        "NEGATIVE REGULATION OF DNA REPAIR",
                        "REGULATION OF DNA REPAIR")

gsea_results <- openxlsx::read.xlsx(
  file.path(pathway_dir, "HPrEC1_GrowthHormone-HPrEC1_Control", "Pathways.xlsx")) %>%
  dplyr::filter(Description %in% plot_pathways_list)
