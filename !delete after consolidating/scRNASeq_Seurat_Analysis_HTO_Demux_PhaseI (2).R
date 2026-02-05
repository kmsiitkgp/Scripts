#!/usr/bin/env Rscript

# Read and store variables from CLI
cli <- commandArgs(trailingOnly = TRUE) 
args <- strsplit(cli, "=", fixed = TRUE)

for (e in args){
  argname <- e[1]
  argval <- e[2]
  assign(argname, argval)
}

# NOTE: All variables and functions are defined within the file below
source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")

# NOTE: Run Phase I first. Then, check diagnostic plots and run Phase II.

library(demuxmix)

#******************************************************************************#
#                                   PHASE I                                    #
#******************************************************************************#

# NOTE: RENAME the HTO directories to match the RNA directories
# You need to use umicounts from CITE-seq-count output
# https://github.com/Hoohm/CITE-seq-Count/issues/105

# Create a list of samples which will be added to each barcode.
# Since folder names correspond to sample name, we just use list.files()
samples <- list.files(path = hto_matrix_path)

demux_samples <- function(i, margin){
  
  #***************************SETUP THE SEURAT OBJECT**************************#
  
  # Read the GEX data files from each sample folder
  umis <- Seurat::Read10X(data.dir = paste0(feature_matrix_path, i),
                          gene.column = 2,
                          cell.column = 1,
                          unique.features = TRUE,
                          strip.suffix = FALSE)
  
  # Read the HTO data files from corresponding HTO folder
  htos <- Seurat::Read10X(data.dir = paste0(hto_matrix_path, i),
                          gene.column = 1,  #features.tsv.gz has HTO tag name in 1st column
                          cell.column = 1,
                          unique.features = TRUE,
                          strip.suffix = FALSE)
  
  # Check if format of barcodes is same in GEX and HTO
  # GEX barcode: "AAACCCAAGAAACACT-1" while HTO barcode: "AAACCCAAGAAACACT"
  print(head(colnames(umis)))
  print(head(colnames(htos)))
  
  # Remove unmapped row from HTOs  
  # https://github.com/Hoohm/CITE-seq-Count/issues/143
  # Reformat the ways tags are named in rownames for easy visualization in plots
  # Eg: Rename "HTO-D-TGGTGTCATTCTTGA" as "HTO-D" and so on
  htos <- data.frame(htos) %>%
    tibble::rownames_to_column("id") %>%
    dplyr::filter(id != "unmapped") %>%
    dplyr::mutate(id = stringr::str_extract(string = id, pattern = "HTO-[A-Z]")) %>%
    tibble::column_to_rownames(var = "id")
  
  # If all barcodes in GEX have "-1" at the end, add "-1" to all barcodes of HTO.
  # If not, we have to add "-2", "-3" etc to barcodes of HTO. So, it is easier
  # if we remove "-1", "-2" etc from the barcodes of GEX
  if (length(colnames(umis)) == length(colnames(umis)[stringr::str_detect(colnames(umis), "-1")])){
    colnames(htos) <- paste0(colnames(htos), pattern = "-1")
  } else{
    colnames(umis) <- stringr::str_replace(string = colnames(umis), pattern = "-\\d", "")
  }
  
  # https://github.com/satijalab/seurat/issues/2549
  # Ideally, we must remove HTOs that have 0 counts in all (100%) cells. 
  # In reality, some cells get labelled with non-specific HTO. 
  # So, we remove HTOs that have 0 counts for more than 99% of cells
  for (row in rownames(htos)){
    cat("Number of cells with zero count for", row, ":", sum(htos[row,] == 0), " Number of cells :", ncol(htos), "\n")
    if (sum(htos[row,] == 0) > 0.99*ncol(htos)){
      htos <-  htos[!(rownames(htos) == row),]
    }
  }
  
  # Remove cells that have 0 counts for all HTOs
  htos <- htos[, colSums(htos) != 0]
  
  # Select cell barcodes detected by both GEX and HTO.
  common_bcs <- intersect(colnames(umis), colnames(htos))
  print(length(common_bcs))
  
  # Subset GEX and HTO counts by common cell barcodes
  umis <- umis[, common_bcs]
  htos <- htos[, common_bcs]
  
  # Setup Seurat object
  sample.seurat <- SeuratObject::CreateSeuratObject(counts = umis,
                                                    project = i,
                                                    assay = "RNA",
                                                    names.field = 1,
                                                    names.delim = "_",
                                                    meta.data = NULL,
                                                    min.cells = 0,
                                                    min.features = 0) # not 100
  
  # Normalize RNA data using log normalization.
  # We can demultiplex without normalizing RNA assay. We are normalizing here
  # to avoid normalization after subsetting. If we dont normalize here, we MUST
  # normalize after subsetting i.e.sample.seurat.subset or singlet Seurat objects
  # normalization.method = "LogNormalize" always performs normalization across 
  # features within a cell irrespective of margin=1 or margin=2.
  sample.seurat <- Seurat::NormalizeData(object = sample.seurat,
                                         assay = "RNA",
                                         normalization.method = "LogNormalize")
  
  #***********************ADD HTO DATA TO SEURAT OBJECT************************#
  
  # Add HTO data as a new assay independent from GEX
  sample.seurat[["HTO"]] <- SeuratObject::CreateAssayObject(counts = htos)
  
  # Normalize HTO data, using centered log-ratio (CLR) transformation. 
  # NOTE: HTO normalization is not affected by GEX data
  
  # During single cell preparation, ambient/background mRNA as well as HTOs from
  # damaged cells are present in the liquid. Despite all cells from a single 
  # sample being tagged with one HTO, when multiple samples are pooled together,
  # these ambient HTOs also get incorporated into single cell GEMs. This is the
  # reason why more than 1 HTO has reads for each cell. So, we need to normalize
  # across all HTOs within each cell similar to how we perform logNormalize
  # across features within a cell so that the correct HTo has highest expression.
  # For this, use margin=2 and CLR normalization.
  # https://github.com/satijalab/seurat/issues/2954
  # https://github.com/satijalab/seurat/issues/2203
  # https://github.com/satijalab/seurat/issues/3605
  # https://divingintogeneticsandgenomics.com/post/details-in-centered-log-ratio-clr-normalization-for-cite-seq-protein-count-data/
  # VERY IMPORTANT: When normalization.method = "CLR", 
  # margin=2 : normalization across features within each cell 
  # margin=1 : normalization across cells of each feature independently
  
  sample.seurat <- Seurat::NormalizeData(object = sample.seurat,
                                         assay = "HTO",
                                         normalization.method = "CLR",
                                         scale.factor = 10000,
                                         margin = margin,
                                         verbose = TRUE)
  
  #*************DEMULTIPLEX CELLS BASED ON HTO ENRICHMENT**********************#
  
  # HTODemux:
  # (i) We perform a k-medoid clustering on the normalized HTO values, which 
  # initially separates cells into K(# of samples)+1 clusters.
  # (ii) We calculate a ‘negative’ distribution for HTO. For each HTO, we use 
  # the cluster with the lowest average value as the negative group.
  # (iii) For each HTO, we fit a negative binomial distribution to the negative 
  # cluster. We use the 0.99 quantile of this distribution as a threshold. 
  # (iv) Based on these thresholds, each cell is classified as positive or 
  # negative for each HTO. Cells that are positive for more than one HTOs are 
  # annotated as doublets.
  
  # MULTIseqDemux:
  # (i) Raw barcode reads were log2-transformed before barcode abundance 
  # normalization via mean subtraction i.e. "CLR normalization" 
  # (ii) Following normalization, the probability density function (PDF) for 
  # each barcode was defined by applying the ‘approxfun’ R function to the 
  # Gaussian kernel density estimation produced using the ‘bkde’ function from 
  # the ‘KernSmooth’ R package. 
  # (iii) Classify cells according to the assumption that groups of cells that 
  # are positive and negative for each barcode should manifest as local PDF 
  # maxima. To this end, we trimmed the top and bottom 0.1% of data from each
  # barcode set and chose the lowest and highest maxima as initial solutions. 
  # To avoid noisy maxima identification, we then adjusted the low maxima to the
  # maxima with the largest number of associated cells. 
  # (iv) With these positive and negative approximations in hand, we next sought
  # to define barcode-specific thresholds. To find the best inter-maxima 
  # quantile for threshold definition (e.g., an inter-maxima quantile of 0.5 
  # corresponds to the mid-point), we iterated across 0.01 quantile increments 
  # and chose the value that maximized the number of singlet classifications. 
  # Optimal inter-maxima distances vary across different MULTI-seq datasets and
  # likely reflect technical noise resulting from variable cell numbers and 
  # labeling efficiency between samples. 
  # (v) Sample classifications were then made using these barcode-specific 
  # thresholds by discerning which thresholds each cell surpasses, with doublets
  # being defined as cells surpassing >1 threshold and negative cells as cells
  # surpassing 0 thresholds.
  
  # # Classify using demuxmix (NOT RECOMMENDED..BAD DEMUXING)
  # dmm <- demuxmix::demuxmix(hto = as.matrix(sample.seurat@assays$HTO$counts),
  #                           rna = sample.seurat@meta.data$nFeature_RNA)
  # classes <- demuxmix::dmmClassify(dmm)
  
  # HTODemux will add 6 columns to metadata "HTO_maxID", "HTO_secondID",
  # "HTO_margin", "HTO_classification", "HTO_classification.global", "hash.ID"
  o <- try(Seurat::HTODemux(object = sample.seurat,
                            assay = "HTO",
                            positive.quantile = 0.5,
                            nstarts = 100,
                            kfunc = "clara",
                            nsamples = 100), silent=TRUE)
  
  if (class(o) == "Seurat"){
    sample.seurat <- Seurat::HTODemux(object = sample.seurat,
                                      assay = "HTO",
                                      positive.quantile = 0.5,
                                      nstarts = 100,
                                      kfunc = "clara",
                                      nsamples = 100)
  }  else if (grepl(pattern="Error",x = o)){
    sample.seurat@meta.data <- sample.seurat@meta.data %>%
      dplyr::mutate(HTO_maxID = rownames(sample.seurat[["HTO"]])[1],
                    HTO_secondID = rownames(sample.seurat[["HTO"]])[1],
                    HTO_margin = 0,
                    HTO_classification = rownames(sample.seurat[["HTO"]])[1],
                    HTO_classification.global = "Negative",
                    hash.ID = "Negative")
  }
  
  # We will perform MULTIseqDemux with autoThresh=TRUE. This will add 2 more 
  # columns to metadata "MULTI_ID" and "MULTI_classification"
  o <- try(Seurat::MULTIseqDemux(object = sample.seurat,
                                 assay = "HTO",
                                 #quantile = 0.7,  #DO NOT SET MANUALLY
                                 autoThresh = TRUE,
                                 maxiter = 10,     #increased from 5
                                 qrange = seq(from = 0.1, to = 0.9, by = 0.05),
                                 verbose = TRUE))
  
  if (class(o) == "Seurat"){
    sample.seurat <- Seurat::MULTIseqDemux(object = sample.seurat,
                                           assay = "HTO",
                                           #quantile = 0.7,  #DO NOT SET MANUALLY
                                           autoThresh = TRUE,
                                           maxiter = 10,     #increased from 5
                                           qrange = seq(from = 0.3, to = 0.9, by = 0.1),
                                           verbose = TRUE)
  }  else if (grepl(pattern="Error",x = o)){
    sample.seurat@meta.data <- sample.seurat@meta.data %>%
      dplyr::mutate(MULTI_classification = rownames(sample.seurat[["HTO"]])[1],
                    MULTI_ID = "Negative")
  }
  
  # Add "MULTI_classification.global" column to metadata for visualization
  sample.seurat@meta.data <- sample.seurat@meta.data %>%
    dplyr::mutate(MULTI_classification.global = dplyr::case_when(grepl(pattern="HTO-[ABCDEF]", x=MULTI_ID) ~ "Singlet",
                                                                 MULTI_ID == "Doublet" ~ "Doublet",
                                                                 MULTI_ID == "Negative" ~ "Negative"))
  
  #**********************VISUALIZE DEMULTIPLEXING RESULTS**********************#
  
  # Global classification results
  print(table(sample.seurat$hash.ID))
  print(table(sample.seurat$MULTI_ID))
  #print(table(classes$HTO))
        
  # Visualize enrichment for HTOs using ridge plots
  # If you see good separation for HTO-A peak from other HTO's peaks while
  # plotting ridgeplot for HTO-A, then it indicates cells are labelled.
  # If you see overlap of HTO-A along with HTO-B (for instance) while plotting
  # ridgeplot for HTO-A, it indicates cells are labelled with both HTO-A & HTO-B
  ridge <- function(y){
    
    Idents(sample.seurat) <- y
    print(y)
    Seurat::RidgePlot(object = sample.seurat,
                      features = rownames(sample.seurat[["HTO"]]),
                      assay = "HTO",
                      layer = "data",
                      same.y.lims = TRUE,
                      combine = TRUE,
                      fill.by = "feature") +
      NoLegend() +
      ggplot2::labs(title = y)
  }
  
  # Visualize enrichment for HTOs using heatmap
  htoheatmap <- function(y){
    Idents(sample.seurat) <- y
    print(y)
    Seurat::HTOHeatmap(object = sample.seurat,
                       assay = "HTO",
                       #classification = paste0(assay, "_classification"),
                       #global.classification = paste0(assay, "_classification.global"),
                       ncells = 5000,
                       singlet.names = NULL,
                       raster = FALSE) +
      ggplot2::labs(title = y)
  }
  
  # Compare number of UMIs for each HTO
  htoviolin <- function(y){
    Idents(sample.seurat) <- y
    print(y)
    Seurat::VlnPlot(object = sample.seurat,
                    features = "nCount_RNA",
                    pt.size = 0.1,
                    log = TRUE) +
      ggplot2::labs(title = y)
  }
  
  # Compare number of UMIs for singlets, doublets and negative cells
  demuxviolin <- function(y){
    
    Idents(sample.seurat) <- y
    print(y)
    Seurat::VlnPlot(object = sample.seurat,
                    features = "nCount_RNA",
                    pt.size = 0.1,
                    log = TRUE) +
      ggplot2::labs(title = y)
  }
  
  funcs <- c("ridge", "htoviolin", "demuxviolin")
  
  for (n in 1:length(funcs)){
    
    if (funcs[n] %in% c("ridge", "htoheatmap", "htoviolin")){
      x <- c("hash.ID", "MULTI_ID")
    } else if (funcs[n] == "demuxviolin"){
      x <- c("HTO_classification.global", "MULTI_classification.global")
    }
    
    purrr::map(.x = x,
               .f = get(funcs[n])) %>% 
      cowplot::plot_grid(plotlist = .,
                         align = "hv",
                         axis = "tblr",
                         nrow = 2,  
                         ncol = 1, 
                         rel_widths = 1,
                         rel_heights = 1,
                         #labels = "AUTO",
                         label_size = 14,
                         label_fontfamily = NULL,
                         label_fontface = "bold",
                         label_colour = NULL,
                         label_x = 0,
                         label_y = 1,
                         hjust = -0.5,
                         vjust = 1.5,
                         #scale = 1,
                         greedy = TRUE,
                         byrow = TRUE,
                         cols = NULL,
                         rows = NULL)
    
    # Save the plot
    ggplot2::ggsave(filename = paste0(funcs[n], "_", i, ".tiff"),
                    plot = last_plot(),
                    device = "jpeg",
                    path = paste0(diagnostics_path, "margin", margin, "/"),
                    #scale = 1,
                    width = 15,
                    height = 11,
                    units = c("in"),	 
                    dpi = 600,
                    limitsize = TRUE,
                    bg = NULL)
  }
  
  # # Visualize pairs of HTO signals to confirm mutual exclusivity in singlets
  # # You have to use <hto_assay_name>_<hashtag_name>
  # Seurat::FeatureScatter(object = sample.seurat,
  #                        feature1 = paste0("hto_", rownames(htos)[1]),
  #                        feature2 = paste0("hto_", rownames(htos)[2]),
  #                        cells = NULL,
  #                        shuffle = FALSE,
  #                        seed = 1,
  #                        group.by = NULL,
  #                        cols = NULL,
  #                        pt.size = 1,
  #                        shape.by = NULL,
  #                        span = NULL,
  #                        smooth = FALSE,
  #                        combine = TRUE,
  #                        slot = "data",
  #                        plot.cor = TRUE,
  #                        raster = FALSE,
  #                        raster.dpi = c(512, 512),
  #                        jitter = TRUE)
  #
  # feature1 = stringr::str_extract(rownames(htos)[1], "HTO-[A-F]")
  # feature2 = stringr::str_extract(rownames(htos)[2], "HTO-[A-F]")
  #
  # ggsave(filename = paste0(i, "_", feature1, "_vs_", feature2, "_Feature_Scatter.tiff"),
  #        plot = last_plot(),
  #        device = "jpeg",
  #        path = results_path,
  #        scale = 1,
  #        #width = 8.5,
  #        #height = 11,
  #        units = c("in"),
  #        dpi = 600,
  #        limitsize = TRUE,
  #        bg = "white")
  
  #*********************DIAGNOSTICS ON SINGLETS & DOUBLETS*********************#
  
  # MOST IMPORTANT PLOTS TO CHECK IF DEMULTIPLEXING WORKS PROPERLY
  for (y in c("HTO_classification.global", "MULTI_classification.global")){
    
    Idents(sample.seurat) <- y
    sample.seurat.subset <- sample.seurat
    
    # We marked all cells as Negative if HTODemux didnt work. Run only if 
    # HTODemux worked earlier.
    if (nrow(sample.seurat@meta.data %>% dplyr::filter(get(y) == "Negative")) != nrow(sample.seurat@meta.data) &
        nrow(sample.seurat@meta.data %>% dplyr::filter(get(y) == "Negative")) > 0){
      
      # Remove negative cells from the object
      sample.seurat.subset <- base::subset(x = sample.seurat.subset,
                                           idents = c("Negative"),
                                           invert = TRUE)
    }
    
    # We dont use FindVariableFeatures() before using ScaleData(). Primary use
    # of FindVariableFeatures() is to reduce computational complexity by 
    # focusing on limited set of variable genes rather than all 30,000+ genes 
    # in a cell. In the HTO assay, we only have 6 features i.e. 6 HTOs.
    
    # Scale variable features of HTO data
    DefaultAssay(object = sample.seurat.subset) <- "HTO"
    sample.seurat.subset <- Seurat::ScaleData(object = sample.seurat.subset,
                                              features = rownames(sample.seurat.subset))
    
    # Run PCA on the HTO data
    sample.seurat.subset <- Seurat::RunPCA(object = sample.seurat.subset,
                                           assay = "HTO",
                                           features = rownames(sample.seurat.subset))
    
    # Run UMAP on the HTO data
    sample.seurat.subset <- Seurat::RunUMAP(object = sample.seurat.subset,
                                            features = rownames(sample.seurat.subset),
                                            reduction = "pca")
    
    # Plot UMAP and check if singlets and doublets separate nicely
    p1 <- Seurat::DimPlot(object = sample.seurat.subset,
                          reduction = "umap",
                          group.by = y,
                          label = FALSE,
                          raster = FALSE,
                          combine = TRUE)
    
    p2 <- Seurat::DimPlot(object = sample.seurat.subset,
                          reduction = "umap",
                          group.by = dplyr::if_else(y == "HTO_classification.global", "hash.ID", "MULTI_ID"),
                          split.by = y,
                          label = FALSE,
                          raster = FALSE,
                          combine = TRUE)
    
    ggsave(filename = paste0(y, "_Singlet_Doublet_UMAP_", i, ".tiff"),
           plot = p1 + p2,
           device = "jpeg",
           path = paste0(diagnostics_path, "margin", margin, "/"),
           scale = 1,
           width = 15,
           height = 7.5,
           units = c("in"),
           dpi = 600,
           limitsize = TRUE,
           bg = "white")
  }
  
  # Save demuxed seurat objects with appropriate sample name
  saveRDS(sample.seurat, file=paste0(demux_results, "margin", margin, "/", i, ".rds"))
} 

# Run first with margin=2. Move seurat objects in "results_seurat" folder and 
# figures in "diagnostics" folder to new folder named margin2. Repeat with
# margin=1 and move output to margin1 folder. Go through all diagnostic plot
# to figure out if (i) MULTIseqdemux or HTodemux is better (ii) margin=1 or 2 is
# better.
for (margin in c(2,1)){
  for (i in samples){
    demux_samples(i, margin)
  }
}

#!/usr/bin/env Rscript

# Read and store variables from CLI
cli <- commandArgs(trailingOnly = TRUE) 
args <- strsplit(cli, "=", fixed = TRUE)

for (e in args){
  argname <- e[1]
  argval <- e[2]
  assign(argname, argval)
}

# NOTE: All variables and functions are defined within the file below
source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")

# NOTE: Run Phase I first. Then, check diagnostic plots and run Phase II.

#******************************************************************************#
#                                   PHASE II                                   #
#******************************************************************************#

# Create a list of samples which will be added to each barcode.
# Since folder names correspond to sample name, we just use list.files()
samples <- list.files(path = hto_matrix_path)

# NOTE: STOP!!! Check UMAP diagnostic plots from previous step before proceeding
# to next step. Based on diagnostic plots, subset singlets from appropriate 
# seurat object based on appropriate method.

for (s in samples){
  
  #******************EXTRACT SINGLETS & SAVE AS SEURAT OBJECT******************#
  
  # IMPORTANT: Define which margin and method to use for each sample
  if (s %in% c("B1","B3","B4")){
    margin <- "1"
    method <- "MULTI_classification.global"
  } else{
    margin <- "2"
    method <- "MULTI_classification.global"
  }
  cat("\n", s, ":", margin, ":", method, "\n")
  
  # Read demuxed seurat object
  sample.seurat <- readRDS(paste0(demux_results, "margin", margin, "/", s, ".rds"))
  Idents(sample.seurat) <- method
  singlet <- base::subset(x = sample.seurat,
                          idents = c("Singlet"))
  
  # IMPORTANT: Define sample specific subsetting to remove BAD HTOs
  if (s %in% c("B1","B4")){
    singlet <- base::subset(x = singlet,
                            MULTI_ID == "HTO-E", 
                            invert = TRUE)
  }
  
  # Assign proper idents based on demux method used
  y <- dplyr::if_else(method == "MULTI_classification.global", "MULTI_ID", "hash.ID")
  
  # Add a new column "HTO_Final" that contains the finalized classification
  singlet@meta.data <- singlet@meta.data %>% dplyr::mutate(HTO_Final = get(y))
  
  cat("\nNumber of singlets:", nrow(singlet@meta.data), "\n")
  
  # Save the singlets with appropriate sample name
  saveRDS(singlet, file=paste0(demux_results, "singlets/", s, ".rds"))
  
  #***************************DIAGNOSTICS ON SINGLETS**************************#
  
  # Cluster and visualize cells using the usual scRNA-seq workflow, and examine
  # for the potential presence of batch effects. Since RNA assay was normalized, 
  # we don't need to normalize again. Even if use normalize, results will be same
  # as logNormalize is done within cells.
  
  # Select the top 1000 most variable features
  singlet <- Seurat::FindVariableFeatures(object = singlet,
                                          assay = "RNA",
                                          selection.method = "vst",
                                          nfeatures = 2000)
  
  # Scaling RNA data, we only scale the variable features here for efficiency
  singlet <- Seurat::ScaleData(object = singlet,
                               features = VariableFeatures(singlet),
                               assay = "RNA")
  
  # Run PCA
  singlet <- Seurat::RunPCA(object = singlet,
                            assay = "RNA",
                            features = VariableFeatures(singlet))
  
  # We select the top 10 PCs for clustering and UMAP
  singlet <- Seurat::FindNeighbors(object = singlet,
                                   reduction = "pca",
                                   dims = 1:40)
  
  singlet <- Seurat::FindClusters(object = singlet,
                                  resolution = 0.6,
                                  verbose = FALSE)
  
  singlet <- Seurat::RunUMAP(object = singlet,
                             reduction = "pca",
                             dims = 1:40)
  
  Idents(sample.seurat) <- y
  
  # Plot UMAP
  p1 <- Seurat::DimPlot(object = singlet,
                        group.by = y)
  
  ggsave(filename = paste0(s, "_", y, "_Singlet_UMAP.tiff"),
         plot = p1,
         device = "jpeg",
         path = diagnostics_path,
         scale = 1,
         width = 15,
         height = 7.5,
         units = c("in"),
         dpi = 600,
         limitsize = TRUE,
         bg = "white")
}

#******************************************************************************#