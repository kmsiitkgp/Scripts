#!/usr/bin/env Rscript

# Read and store variables from command line interface (CLI)
cli <- base::commandArgs(trailingOnly = TRUE) 
args <- base::strsplit(x = cli, split = "=", fixed = TRUE)

for (e in args){
  argname <- e[1]
  argval <- e[2]
  assign(argname, argval)
}

# NOTE: All variables and functions are defined within the file below
source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")

#******************************************************************************#
#                       STEP 1: SETUP THE SEURAT OBJECT                        #
#******************************************************************************#

#************************IMPORTING DATA FROM h5AD FILE*************************#

if (h5ad_file == TRUE){
  
  # Load h5ad (useful if analyzing collaborator data in h5ad format)
  SeuratDisk::Convert(source = paste0(seurat_results, proj, ".h5ad"),
                      dest = "h5seurat",
                      assay="RNA",
                      overwrite = FALSE)
  
  raw_seurat <- SeuratDisk::LoadH5Seurat(file = paste0(seurat_results, proj, ".h5seurat"))
  
  # View the contents of the seurat object
  dplyr::glimpse(raw_seurat)
}

#*********************IMPORTING GEX (GENE EXPRESSION) DATA*********************#

if (demultiplexed_data == TRUE){
  
  # Create a list of samples that have been demultiplexed already
  files <- list.files(path = paste0(demux_results, "singlets/"),
                      full.names = FALSE)
  samples <- gsub(pattern="\\..*", replacement="", x=files)
  
  # Loop through each of the individual object in demux directory & import data
  for (i in 1:length(files)){
    
    # Read the seurat object containing demultiplexed singlets
    sample.seurat <- readRDS(file = paste0(demux_results, "singlets/", files[i]))
    
    # Assign the seurat object to its corresponding variable
    assign(samples[i], sample.seurat)
    
    # Explore the meta.data slot
    cat("\nFirst few rows of ", i, "\n")
    print(head(sample.seurat@meta.data))
    cat("\nLast few rows of ", i, "\n")
    print(tail(sample.seurat@meta.data))
  }
} 

if (demultiplexed_data == FALSE){
  
  # Create a list of sample names which will be added to each barcode.
  # Since folder names correspond to sample name, we just use list.files()
  samples <- list.files(path = raw_matrix_path)
  
  # Loop through each of the individual folders in parent directory & import data
  for(i in samples){
    for (matrix_path in c(raw_matrix_path, filt_matrix_path)){
   
      # Read the data files from each sample folder
      # NOTE: gene.column=1 imports Ensembl ids from features.tsv. DO NOT DO THIS. 
      # Read gene symbols from features.tsv as we will calculate mitoratio, 
      # riboratio etc using gene names
      dgCMatrix <- Seurat::Read10X(data.dir = paste0(matrix_path, i),
                                   gene.column = 2,  
                                   cell.column = 1,
                                   unique.features = TRUE,
                                   strip.suffix = FALSE)
      
      # Create a seurat object for each dgCMatrix object
      # NOTE: Set min.features=100 to remove low quality cells & reduce object size
      seurat <- SeuratObject::CreateSeuratObject(counts = dgCMatrix,
                                                 project = i,
                                                 assay = "RNA",
                                                 names.field = 1,
                                                 names.delim = "_",
                                                 meta.data = NULL,
                                                 min.cells = 0,
                                                 min.features = 100)
      
      # Assign the seurat object to its corresponding variable
      if (matrix_path == raw_matrix_path){
        assign(paste0(i, ".raw"), seurat)
        cat("DATA IMPORTED FOR ", i, ".raw dataset\n")
      } else{
        assign(paste0(i, ".filt"), seurat)
        cat("DATA IMPORTED FOR ", i, ".filt dataset\n")
      }
    }
  }
}

#******************************************************************************#
#                           STEP 2: QUALITY CONTROL                            #
#******************************************************************************#

#**********************STEP 2A: CALCULATE ALL QC METRICS***********************#

# NOTE:  We are going to do the same QC on each sample. So, we can merge the 
# individual seurat objects into a single seurat object and calculate metrics on
# the merged seurat object. However, if the number of cells exceeds 1.5 million,
# then this step will fail as there are simply too many cells. So, it is better
# to calculate QC metrics on individual samples, perform filtering on individual
# samples and then merge the filtered seurat objects which have fewer cells.

# Initialize an empty dataframe where the class of each column resembles those 
# of raw_metadata. We will use this dataframe for making QC plots later.
raw_metadata <- data.frame(Cell = c(""), 
                           Sample = as.factor(1), 
                           nUMIs = c(0), 
                           nGenes = c(0), 
                           MitoRatio = c(0), 
                           RiboRatio = c(0), 
                           Novelty = c(0))

# Calculate QC metrics for each sample individually using raw matrices
for(i in paste0(samples, ".raw")){
  
  sample.seurat <- get(i)
  
  # Compute percent mito percent
  sample.seurat <- Seurat::PercentageFeatureSet(object = sample.seurat,
                                                pattern = "^[Mm][Tt]-",
                                                features = NULL,
                                                col.name = "MitoPercent",
                                                assay = "RNA")
  
  # Compute percent ribo percent
  sample.seurat <- Seurat::PercentageFeatureSet(object = sample.seurat,
                                                pattern = "^[Rr][Pp][SsLl]", 
                                                features = NULL,
                                                col.name = "RiboPercent",
                                                assay = "RNA")
  
  # Extract metadata
  sample_metadata <- sample.seurat@meta.data
  
  # Rename columns to be more intuitive and add the additional QC metrics:
  # (i)     Cell      : unique identifiers corresponding to each cell i.e. barcodes
  # (ii)    Sample    : sample names
  # (iii)   nUMIs     : number of transcripts per cell
  # (iv)    nGenes    : number of genes per cell
  # (v)     nHTO_UMIs : number of HTO reads per cell
  # (vi)    nHTOs     : number of HTO types per cell
  # (vii)   MitoRatio : MitoPercent/100
  # (viii)	RiboRatio : RiboPercent/100  
  # (ix)    Novelty   : log ratio of genes per UMI
  sample_metadata <- sample_metadata %>% 
    dplyr::mutate(Cell = paste0(orig.ident, "_", rownames(sample_metadata)),
                  Sample = orig.ident,
                  nUMIs = nCount_RNA,
                  nGenes = nFeature_RNA,
                  MitoRatio = MitoPercent/100,
                  RiboRatio = RiboPercent/100,
                  Novelty = log10(nGenes)/log10(nUMIs), .keep="unused")
  
  if (hto_info == TRUE){
    sample_metadata <- sample_metadata %>% 
      dplyr::mutate(nHTO_UMIs = nCount_HTO,
                    nHTOs = nFeature_HTO,
                    HTO_Final = HTO_Final) %>%
      dplyr::select(Cell, Sample, nUMIs, nGenes, nHTO_UMIs, nHTOs, HTO_Final, MitoRatio, RiboRatio, Novelty)
  }
  
  if (hto_info == FALSE){
    sample_metadata <- sample_metadata %>% 
      dplyr::mutate(nHTO_UMIs = 0,
                    nHTOs = 0,
                    HTO_Final = NA) %>%
      dplyr::select(Cell, Sample, nUMIs, nGenes, nHTO_UMIs, nHTOs, HTO_Final, MitoRatio, RiboRatio, Novelty)
  }
  
  # Replace metadata in raw Seurat object with updated column names
  sample.seurat@meta.data <- sample_metadata
  
  # Append raw metadata of each seurat object which will be used for QC plots later
  raw_metadata <- dplyr::bind_rows(raw_metadata, sample_metadata)
  
  # Assign the seurat object to its corresponding variable
  assign(i, sample.seurat)
}

#******************************STEP 2B: PERFORM QC*****************************#

# Perform QC for each sample individually
for(i in paste0(samples, ".raw")){
  
  sample.seurat <- get(i)
  
  gene_cutoff <- 250
  umi_cutoff <- 500
  mito_cutoff <- 0.2
  ribo_cutoff <- 0.05
  novelty_cutoff <- 0.8  	    # use 0.8 as starting point. Maximum 0.9
  
  if (snRNASeq == TRUE){
    mito_cutoff <- 0.05
  }
  
  sample.seurat <- base::subset(x = sample.seurat,
                                subset = (nGenes >= gene_cutoff) &
                                  (nUMIs >= umi_cutoff) &
                                  (MitoRatio <= mito_cutoff) &
                                  # (RiboRatio >= ribo_cutoff) &
                                  (Novelty >= novelty_cutoff))
  
  # Assign the seurat object to its corresponding variable
  assign(i, sample.seurat)
}

# In case of Simon data, the tumors samples were hashtagged but normal samples
# were not. So, perform QC on demultiplexed tumor samples and then on the 
# normal samples. Define samples <- list.files(path = feature_matrix_path) and
# proceed with creating filtered seurat object

# Create a merged Seurat object
# NOTE: Samples will have same barcodes. To keep track of cell identities 
# (i.e. barcodes) coming from each sample after merging, we add a prefix 
# (i.e. sample name) to each barcode using "add.cell.ids"

# This has all barcodes BEFORE removal of empty droplets and AFTER QC
filtered_pre_emptydrops_seurat <- base::merge(x = get(paste0(samples[1], ".raw")),
                                              y = lapply(paste0(samples[2:length(samples)], ".raw"), get),
                                              add.cell.ids = samples,
                                              merge.data = FALSE)

# Create a merged seurat object from Cellranger filtered matrices
# This has all barcodes AFTER removal of empty droplets but WIHTOUT ANY QC
filt_seurat <- base::merge(x = get(paste0(samples[1], ".filt")),
                           y = lapply(paste0(samples[2:length(samples)], ".filt"), get),
                           add.cell.ids = samples,
                           merge.data = FALSE)

# Identify common barcodes between our filtering and Cellranger filtered matrices
common_bc <- c()
for (i in unique(filtered_pre_emptydrops_seurat@meta.data$Sample)){
  
  common_bc <- c(common_bc, 
                 intersect(rownames(filtered_pre_emptydrops_seurat@meta.data %>% dplyr::filter(Sample == i)),
                           rownames(filt_seurat@meta.data %>% dplyr::filter(orig.ident == i))))
}

# Remove barcodes classified as emptydroplets by CellRanger 
# NOTE: We could start analysis directly with filtered matrices but then we wont 
# know how many cells we are losing.
# This has all barcodes AFTER removal of empty droplets and AFTER QC
filtered_post_emptydrops_seurat <- subset(filtered_pre_emptydrops_seurat,
                                          Cell %in% common_bc)

# Remove HTO assay from to avoid complications during integration, etc
if (demultiplexed_data == TRUE){
  filtered_post_emptydrops_seurat[["HTO"]] <- NULL
}

# If whitelist is not needed, add additional metadata
if (whitelist == FALSE){
  
  # Import other meta data associated with data set.
  # NOTE: This xslx file should have column named "Unique_ID" whose values matches 
  # with  column "Unique_ID" of seurat object's metadata.
  extra_metadata <-  openxlsx::read.xlsx(xlsxFile = paste0(scripts_path, metafile))
  
  # Merge imported metadata with existing metadata
  # NOTE: left_join etc will remove rownames. So, add rownames before replacing
  # metadata in Seurat object
  if (hto_info == TRUE){
    filtered_post_emptydrops_seurat@meta.data <- filtered_post_emptydrops_seurat@meta.data %>%
      dplyr::mutate(Unique_ID = paste0(Sample, "_", HTO_Final)) %>%
      dplyr::left_join(extra_metadata, by=("Unique_ID"="Unique_ID")) %>%
      dplyr::mutate(index = Cell) %>%
      tibble::column_to_rownames(var = "index")
  } else{
    filtered_post_emptydrops_seurat@meta.data <- filtered_post_emptydrops_seurat@meta.data %>%
      dplyr::mutate(Unique_ID = paste0(Sample)) %>%
      dplyr::left_join(extra_metadata, by=("Unique_ID"="Unique_ID")) %>%
      dplyr::mutate(index = Cell) %>%
      tibble::column_to_rownames(var = "index")
  }
  
  #****************************STEP 2C: SAVE THE DATA****************************#
  
  # Create .rds object for filtered seurat object to load at any time
  saveRDS(filtered_post_emptydrops_seurat, file=paste0(seurat_results, "filtered_seurat.rds"))
}  

# Plot graphs before demultiplexing or on data that doesnt need demultiplexing
if (demultiplexed_data == FALSE){
  
  #*****************STEP 2D: VISUALIZE DATA BEFORE AND AFTER QC******************#
  
  # Remove dummy first row from raw_metadata
  raw_metadata <- raw_metadata[-1,]
  rownames(raw_metadata) <- raw_metadata$Cell
  
  # Pre-removal of Empty droplets
  pre_metadata <- filtered_pre_emptydrops_seurat@meta.data
  
  # Post-removal of Empty droplets
  filtered_metadata <- filtered_post_emptydrops_seurat@meta.data
  
  # Visualize the number of cell counts per sample
  cell_qc <- function(metadata, tag){
    
    ggplot(data = get(metadata), aes(x=Sample, fill=Sample)) + 
      geom_bar() +              
      theme_classic() +         #display with x and y axis lines and no gridlines
      my_theme +
      labs(x = "Sample", y = "Number of Cells", title = stringr::str_wrap(paste0("Number of Cells ", tag), 30)) +
      coord_cartesian(clip = "off") +
      geom_text(stat="count", aes(label=after_stat(count)), y = 0, hjust = 0, angle = 90)
  }
  
  # Visualize the number of UMIs per cell
  umi_qc <- function(metadata, tag){
    
    ggplot(data = get(metadata), aes(x=Sample, y=nUMIs, fill=Sample)) +
      geom_violin() + 
      geom_boxplot(width=0.1) +
      theme_classic() +  
      my_theme +
      labs(x = "Sample", y = "Number of UMIs", title = stringr::str_wrap(paste0("Distribution of UMIs ", tag),30)) +
      coord_cartesian(ylim = c(100, 1000000), clip = "off") +
      scale_y_log10(breaks = c(100, 1000, 10000, 100000, 1000000)) +  		     #display y axis in log scale
      geom_hline(yintercept = gene_cutoff, linetype = 2)
  }
  
  # Visualize the number of genes per cell
  gene_qc <- function(metadata, tag){
    
    ggplot(data = get(metadata), aes(x=Sample, y=nGenes, fill = Sample)) +
      geom_violin() + 
      geom_boxplot(width=0.1) +
      theme_classic() + 
      my_theme + 
      labs(x = "Sample", y = "Number of Genes", title = stringr::str_wrap(paste0("Distribution of Genes ", tag),30)) +
      coord_cartesian(ylim = c(1, 30000), clip = "off") +
      scale_y_log10() +    	
      geom_hline(yintercept = gene_cutoff, linetype = 2)
  }
  
  # Visualize the MitoRatio of each cell
  mito_qc <- function(metadata, tag){
    
    ggplot(data = get(metadata), aes(x=Sample, y=MitoRatio, fill = Sample)) +
      geom_violin() + 
      geom_boxplot(width=0.1) + 
      theme_classic() + 
      my_theme +
      labs(x = "Sample", y = "MitoRatio", title = stringr::str_wrap(paste0("Distribution of MitoRatio ", tag),30)) +
      coord_cartesian(ylim = c(0, 1), clip = "off") +
      geom_hline(yintercept = mito_cutoff, linetype = 2)
  }
  
  # Visualize the RiboRatio of each cell
  ribo_qc <- function(metadata, tag){
    
    ggplot(data = get(metadata), aes(x=Sample, y=RiboRatio, fill = Sample)) +
      geom_violin() + 
      geom_boxplot(width=0.1) + 
      theme_classic() + 
      my_theme +
      labs(x = "Sample", y = "RiboRatio", title = stringr::str_wrap(paste0("Distribution of RiboRatio ", tag),30)) +
      coord_cartesian(ylim = c(0, 1), clip = "off") +
      geom_hline(yintercept = ribo_cutoff, linetype = 2)
  }
  
  # Visualize the novelty or complexity of each cell
  novelty_qc <- function(metadata, tag){
    
    ggplot(data = get(metadata), aes(x=Sample, y=Novelty, fill = Sample)) +
      geom_violin() + 
      geom_boxplot(width=0.1) + 
      theme_classic() + 
      my_theme + 
      labs(x = "Sample", y = "Novelty Score", title = stringr::str_wrap(paste0("Distribution of Novelty Score ", tag),30)) +
      coord_cartesian(ylim = c(0, 1), clip = "off") +
      geom_hline(yintercept = novelty_cutoff, linetype = 2)
  }
  
  # Visualize number of genes/cell, number of UMIs/cell & MitoRatio together.
  # Bottom left quadrant : Poor quality cells with low genes & UMIs per cell 
  # Top right quadrant   : Good quality cells with high genes & UMIs per cell
  # Bottom right quadrant: Cells with low genes but high UMIs per cell. These 
  # could be dying cells or population of low complexity cells (i.e erythrocytes)
  gene_umi_mito_qc <- function(metadata, tag){
    
    ggplot(data = get(metadata), aes(x=nUMIs, y=nGenes, color = MitoRatio)) +
      geom_point() +
      theme_classic() + 
      my_theme + 
      labs(x = "Number of UMIs", y = "Number of Genes",	 title = paste0("Distribution of UMIs, Genes & MitoRatio ", tag)) +
      coord_cartesian(xlim = c(100, 1000000), ylim = c(100, 20000), clip = "off") +
      scale_x_log10(breaks = c(100, 1000, 10000, 100000, 1000000)) + 
      scale_y_log10() + 
      facet_wrap(.~Sample, nrow = 4) +   #split the plot by X-axis label
      stat_smooth(method=lm, color="yellow") +
      geom_vline(xintercept = umi_cutoff) +    	#draw a vertical line at x=500 i.e.UMIs cutoff
      geom_hline(yintercept = gene_cutoff) +    #draw a horizontal line at y =250 i.e. Genes cutoff
      scale_color_viridis(option = "D", limits = c(0, 1)) 		# limits sets max and min values of gradient 
  }
  
  # Plot all QC metrics before and after QC
  funcs <- c("cell_qc", "umi_qc", "gene_qc", "mito_qc", "ribo_qc", "novelty_qc",
             "gene_umi_mito_qc")
  
  filenames <- c("Cell_Counts", "UMI_Distribution", "Gene_Distribution",
                 "MitoRatio_Distribution", "RiboRatio_Distribution", 
                 "Novelty_Score_Distribution", "Genes_UMI_MitoRatio_Distribution")
  
  for (i in 1:length(funcs)){
    
    # Plot QC metrics
    purrr::map2(.x = c("raw_metadata", "pre_metadata", "filtered_metadata"),
                .y = c("Pre QC", "Post QC (Before Empty Drops)", "Post QC (After Empty Drops)"),
                .f = get(funcs[i])) %>% 
      cowplot::plot_grid(plotlist = .,
                         align = "hv",
                         axis = "tblr",
                         nrow = 3,  
                         ncol = 1, 
                         rel_widths = 1,
                         rel_heights = 1,
                         labels = NULL,
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
    ggplot2::ggsave(filename = paste0("QC_", filenames[i], ".pdf"),
                    plot = last_plot(),
                    device = "pdf",
                    path = seurat_results,
                    #scale = 1,
                    width = dplyr::if_else(i==7, 17,5+length(samples)/2),
                    height = dplyr::if_else(i==7, 22, 11),
                    units = c("in"),	 
                    dpi = 600,
                    limitsize = TRUE,
                    bg = NULL)
  }
}

# Generate whitelist using filtered seurat and exit script.
if (whitelist == TRUE){
  
  # Extract barcodes and split by "_"
  bc <- filtered_post_emptydrops_seurat@meta.data$Cell
  
  # Adjust this based on how your samples are named
  barcodes <- data.frame(stringr::str_split_fixed(bc, "_", 2)) %>%
    dplyr::rename(Batch = identity(1), Barcodes = identity(2)) %>%
    dplyr::mutate(Barcodes = stringr::str_replace(Barcodes, "-1", ""),
                  Batch = gsub(pattern="-GEX.*", replacement="", x=Batch))
  
  # Remove duplicate barcodes within each batch
  barcodes <- barcodes %>% 
    dplyr::group_by(Batch) %>% 
    dplyr::distinct_at("Barcodes", .keep_all=TRUE) %>% 
    as.data.frame()
  
  # Check how many barcodes are present in each batch
  barcodes %>% dplyr::group_by(Batch) %>% dplyr::count()
  
  # Save barcodes to individual csv files
  for (i in unique(barcodes$Batch)){
    whitelist <- barcodes %>%
      dplyr::filter(Batch == i) %>%
      dplyr::select(Barcodes)
    
    write.table(x = whitelist,
                file = paste0(scripts_path, proj, "_", i, "_whitelist.csv"),
                row.names = FALSE,
                col.names = FALSE)
  }
  
  #base::quit(save = "no")
}

#******************************************************************************#
#                       STEP 3: RUN THE STANDARD PIPELINE                      #
#******************************************************************************#

# cell_type has classification based on UMAP clusters
# ucell_class has classification based on UCell scores
# seurat_class has classification based on seurat scores

if (whitelist == FALSE){
  
  # Load rds file of filtered_seurat object
  filtered_seurat <- base::readRDS(paste0(seurat_results, "filtered_seurat.rds"))

  filt <- filtered_seurat
  sct <- sctransform_data(filt)
  
  kweight <- sct@meta.data %>% dplyr::count(Sample) %>% dplyr::select(n) %>% min()/2
  kweight <- min(kweight, 100) 
  integ <- integrate_data(sct, kweight)
  celltype <- NULL
  integ <- cluster_data(integ, celltype)
  # integ <- add_module_scores(integ, "All Markers")
  # integ <- annotate_data_score(integ, celltype)
  save_data(integ, celltype)
  reduc <- "Harmony"
  res <- 0.4
  get_markers(integ, res, reduc, celltype)
  
  # for (reduc in c("CCA", "RPCA", "Harmony", "JointPCA")){
  #   res <- 1.4
  #   plot_conserved_modules(res, reduc, celltype, "All Markers")
  # }
  
  
  #plot_pre_integration(sct)
  # reduc <- "Harmony"
  # res <- 1.4
  # idents <- paste0("cluster.", res, ".", base::tolower(reduc))
  # plot_post_integration(res, reduc, idents, celltype)
}

#******************************************************************************#

# Next, look at "Module_plot(All Markers)__CCA.jpg", 
# "Module_plot(All Markers)__Harmony.jpg", "Module_plot(All Markers)__RPCA.jpg",
# "Module_plot(All Markers)__JointPCA.jpg" plots to decide which reduction gives
# (i) clear separation of clusters 
# (ii) clusters cells together i.e. B cells are not split between clusters.
# In all cases so far, RPCA and Harmony seemed to work best.

#******************************************************************************#

# # Finding number of common barcodes between samples in Simon snRNASeq
# for (i in samples){
#   barcodes <- rownames(get(i)@meta.data)
#   t <- paste0(i, ".cell")
#   assign(t, barcodes)
# }
# 
# samples.cell <- paste0(samples, ".cell")
# # Identify all possible combinations of 2 samples
# combinations <- utils::combn(x=samples.cell, m=2)
# combinations_result <- c()
# total_sample1 <- c()
# total_sample2 <- c()
# for (i in 1:ncol(combinations)){
#   
#   # Find number of common barcodes for any two samples
#   combinations_result <- c(combinations_result, 
#                            length(intersect(get(combinations[1,i]), get(combinations[2,i]))))
#   total_sample1 <- c(total_sample1, length(get(combinations[1,i])))
#   total_sample2 <- c(total_sample2, length(get(combinations[2,i])))
# }
# 
# results <- data.frame(t(combinations), combinations_result, total_sample1, total_sample2)
# colnames(results) <- c("Sample1", "Sample2", "Number of common barcodes", 
#                        "Total barcodes Sample1", "Total barcodes Sample2")
# 
# wb <- openxlsx::createWorkbook()
# openxlsx::addWorksheet(wb = wb, sheetName = "Summary")
# openxlsx::writeData(wb = wb, sheet = "Summary", x = results)
# openxlsx::saveWorkbook(wb = wb, file = "Summary_QC.xlsx", overwrite = TRUE)

#******************************************************************************#

# # Calculate number of cells that remain at different QC settings 
# summary_df <- data.frame(Sample = "NA")
# cutoffs <- list(gene_cutoff  = c(250, 500, 1000, 5000),
#                 umi_cutoff = c(500, 1000, 10000, 100000),
#                 mito_cutoff = c(0.05, 0.1, 0.2, 0.3),
#                 novelty_cutoff = c(0.8, 0.85, 0.9, 0.95))
# meta.data_cols <- c("nGenes", "nUMIs", "MitoRatio", "Novelty")
# 
# # Loop through each cutoff and calculate cells passing QC
# for (i in 1:4){
#   
#   header_df <- data.frame(Sample = "Cutoffs", 
#                           gene_cutoff = as.numeric(cutoffs[[1]][i]),
#                           umi_cutoff = as.numeric(cutoffs[[2]][i]),
#                           mito_cutoff = as.numeric(cutoffs[[3]][i]),
#                           novelty_cutoff = as.numeric(cutoffs[[4]][i]))
#   
#   df <- data.frame(Sample = unique(raw_seurat@meta.data$Sample))
#   
#   for (j in 1:4){
#     
#     if(j == 3){  # for mito_cutoff, we need to keep cells below than cutoff
#       df1 <- raw_seurat@meta.data %>% 
#         dplyr::group_by(Sample) %>% 
#         dplyr::filter(!!rlang::sym(meta.data_cols[j]) <= cutoffs[[j]][i]) %>%
#         dplyr::count() %>%
#         dplyr::rename(Sample = identity(1), !!names(cutoffs)[j] :=  identity(2))
#     } else { # for all other cutoffs, we need to keep cells above than cutoff
#       df1 <- raw_seurat@meta.data %>% 
#         dplyr::group_by(Sample) %>% 
#         dplyr::filter(!!rlang::sym(meta.data_cols[j]) >= cutoffs[[j]][i]) %>%
#         dplyr::count() %>%
#         dplyr::rename(Sample = identity(1), !!names(cutoffs)[j] :=  identity(2))
#     }
#     df <- df %>% dplyr::left_join(df1, by=c("Sample"="Sample"))
#   }
#   summary_df <- dplyr::bind_rows(summary_df, header_df, df)
# }
# 
# # Save the summary
# wb <- openxlsx::createWorkbook()
# openxlsx::addWorksheet(wb = wb, sheetName = paste0("Summary"))
# openxlsx::writeData(wb = wb, sheet = paste0("Summary"), x = summary_df)
# openxlsx::saveWorkbook(wb = wb, file = paste0(seurat_results, "Summary_QC.xlsx"), overwrite = TRUE)

#******************************************************************************#

# ggplot2::ggplot(data = filtered_seurat@meta.data, aes(x=nUMIs, y=Sample, fill= after_stat(x))) +
#   geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
#   viridis::scale_fill_viridis(name = "nUMIs", alpha = 1, begin = 0, end = 1, 
#                               direction = 1, discrete = FALSE, option = "D") +
#   labs(title = 'UMI Distribution') +
#   coord_cartesian(xlim = c(100, 10000)) +
#   scale_x_continuous(breaks = c(100, 1000, 2500, 5000, 10000))+
#   #hrbrthemes::theme_ipsum() +
#   theme(legend.position="right",
#         panel.spacing = unit(0.1, "lines"),
#         strip.text.x = element_text(size = 8))
# 
# ggsave("1.jpg", bg="white")  

#!/usr/bin/env Rscript

# Read and store variables from command line interface (CLI)
cli <- base::commandArgs(trailingOnly = TRUE) 
args <- base::strsplit(x = cli, split = "=", fixed = TRUE)

for (e in args){
  argname <- e[1]
  argval <- e[2]
  assign(argname, argval)
}

# NOTE: All variables and functions are defined within the file below
source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables_new.R")

#******************************************************************************#

# Create a list of sample names which will be added to each barcode.
# Since folder names correspond to sample name, we just use list.files()
samples <- list.files(path = raw_matrix_path) 

# Create empty dataframe to store raw metadata
raw_metadata <- data.frame(Cell = c("")) #Sample = as.factor(1), nUMIs = c(0)

for (s in samples){
  s.obj <- read_cellranger(s, raw_matrix_path)
  s.obj <- mark_emptydroplets_dropletutils(s.obj)
  s.obj <- mark_emptydroplets_cellranger(s.obj)
  s.obj <- doublet_finder(s.obj)
  s.obj <- scdbl_finder(s.obj)
  s.obj <- calc_qc_metrics(s.obj)
  s.obj <- mark_low_quality(s.obj)
  raw_metadata <- generate_plotdata(s.obj, raw_metadata)
  s.obj <- filter_singlets(s.obj)
  assign(s, s.obj)
}

plot_qc(raw_metadata, seurat_results)
#xpectr::suppress_mw(generate_whitelist(filt, seurat_results))
filt              <- format_filtered(samples, seurat_results)
sct               <- sctransform_singlecell(filt, seurat_results)
reference.samples <- NULL
kweight <- min(sct@meta.data %>% dplyr::count(Sample) %>% dplyr::select(n) %>% min()/2, 100) 
integ             <- integrate_singlecell(sct, reference.samples, kweight, seurat_results)
clust             <- cluster_singlecell(integ, seurat_results)
final             <- remove_sparse_clusters(clust, seurat_results)
suffix <- "Full"
plot_metrics_post_integration(final, suffix, diagnostics_path)
resolution <- 0.8
reduction <- "Harmony"
species <- "Mus musculus"
identify_markers(final, resolution, reduction, species, suffix, seurat_results)

# Decide which reduction gives (i) clear separation of clusters 
# (ii) clusters cells together i.e. B cells are not split between clusters.
# In all cases so far, RPCA and Harmony seemed to work best.

#!/usr/bin/env Rscript

# Read and store variables from command line interface (CLI)
cli <- base::commandArgs(trailingOnly = TRUE) 
args <- base::strsplit(x = cli, split = "=", fixed = TRUE)

for (e in args){
  argname <- e[1]
  argval <- e[2]
  assign(argname, argval)
}

#***********************STEP 2A: IDENTIFY EMPTY DROPLETS***********************#

# Empty droplets often contain RNA from the ambient solution, resulting in 
# non-zero counts. These empty droplets MUST be removed before identifying 
# doublets. The emptyDrops function is designed to distinguish between empty 
# droplets and cells by testing each barcode’s expression profile for 
# significant deviation from the ambient profile. The null hypothesis is that
# "transcript molecules are included into droplets by multinomial sampling from
# the ambient profile" i.e. "droplets contain ambient RNA". So, if FDR < 0.05, 
# we reject null hypothesis.

# NOTE: The ouput of emptyDrops() is a dataframe with columns 
# Total   : total UMI count for each barcode (droplet)
# LogProb : log-probability of observing the barcode's count vector under the null model.
# PValue  : the Monte Carlo p-value against the null model
# Limited : whether a lower p-value could be obtained by increasing ‘niters’
# FDR     : Non-empty droplets have FDR < 0.05 (i.e. 5%)
# glimpse(e.out) shows metadata associated with the dataframe also
# ?emptyDrops shows explanation of why there are NA values in FDR

#**************************STEP 2B: IDENTIFY DOUBLETS**************************#

# NOTE: Doublet finding algorithms need good quality data as starting point for
# correct estimation of parameters used in identifying doublets. So, we perform 
# this step after removing empty droplets.

# cell_type has classification based on UMAP clusters
# ucell_class has classification based on UCell scores
# seurat_class has classification based on seurat scores

integ <- add_module_scores(integ, "All Markers")
integ <- annotate_data_score(integ, celltype)
save_data(integ, celltype)


# for (reduc in c("CCA", "RPCA", "Harmony", "JointPCA")){
#   res <- 1.4
#   plot_conserved_modules(res, reduc, celltype, "All Markers")
# }
#plot_pre_integration(sct)

#******************************************************************************#

# # Finding number of common barcodes between samples in Simon snRNASeq
# for (i in samples){
#   barcodes <- rownames(get(i)@meta.data)
#   t <- paste0(i, ".cell")
#   assign(t, barcodes)
# }
# 
# samples.cell <- paste0(samples, ".cell")
# # Identify all possible combinations of 2 samples
# combinations <- utils::combn(x=samples.cell, m=2)
# combinations_result <- c()
# total_sample1 <- c()
# total_sample2 <- c()
# for (i in 1:ncol(combinations)){
#   
#   # Find number of common barcodes for any two samples
#   combinations_result <- c(combinations_result, 
#                            length(intersect(get(combinations[1,i]), get(combinations[2,i]))))
#   total_sample1 <- c(total_sample1, length(get(combinations[1,i])))
#   total_sample2 <- c(total_sample2, length(get(combinations[2,i])))
# }
# 
# results <- data.frame(t(combinations), combinations_result, total_sample1, total_sample2)
# colnames(results) <- c("Sample1", "Sample2", "Number of common barcodes", 
#                        "Total barcodes Sample1", "Total barcodes Sample2")


#******************************************************************************#

# # Check QC metrics for singlets and doublets
# VlnPlot(object = sample.seurat,
#         group.by = "Sample",
#         split.by = "DF.Class",
#         features = c("nUMIs", "nGenes", "MitoRatio", "RiboRatio", "HemeRatio", "Novelty"))

# ggplot2::ggplot(data = filtered_seurat@meta.data, aes(x=nUMIs, y=Sample, fill= after_stat(x))) +
#   geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
#   viridis::scale_fill_viridis(name = "nUMIs", alpha = 1, begin = 0, end = 1, 
#                               direction = 1, discrete = FALSE, option = "D") +
#   labs(title = 'UMI Distribution') +
#   coord_cartesian(xlim = c(100, 10000)) +
#   scale_x_continuous(breaks = c(100, 1000, 2500, 5000, 10000))+
#   #hrbrthemes::theme_ipsum() +
#   theme(legend.position="right",
#         panel.spacing = unit(0.1, "lines"),
#         strip.text.x = element_text(size = 8))
# 
# ggsave("1.jpg", bg="white")  