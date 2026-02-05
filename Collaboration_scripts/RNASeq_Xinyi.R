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

