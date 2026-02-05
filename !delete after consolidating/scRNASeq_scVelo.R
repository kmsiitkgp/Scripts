# Refer https://jef.works/blog/2020/08/25/using-scvelo-in-R-using-reticulate/
# Refer https://smorabit.github.io/tutorials/8_velocyto/
# Refer figures in https://scvelo.readthedocs.io/en/stable/vignettes/
# https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html

# The steps are as follows:
# 1. Get input matrices (unspliced counts, spliced counts) with velocyto
# 2. Get your clusters/reduced dimensions etc. from Seurat
# 3. Merge the information for scVelo to calculate RNA velocity

#****************STEP 1: PREPARE ENVIRONMENT FOR scvelo ANALYSIS***************#

# Create new conda env with latest python and install necessary packages
# NOTE: Sometimes scvelo may not work with latest python version.
# NOTE: Use login(submit) node to install all packages, not computing node.
# Run the commands below in terminal to check scvelo is installed properly
# You can close this terminal once scvelo is installed properly

conda create --name scVelo python=3.9
conda activate scVelo
pip install -U scvelo
pip install Cython                  #needed for velocyto installation
pip install velocyto
conda install -c bioconda samtools  #needed for velocyto to work

#*STEP 2: GENERATE COUNT MATRICES FOR SPLICED AND UNSPLICED RNA USING VELOCYTO*#

# Run "qsub 02c_velocyto.sh" in HPC cluster and once it completes, copy the loom 
# files to "/hpc/home/kailasamms/scratch/scRNASeq_BBN_C57B6/results_velocyto/"
# and then run the code below in python using scVelo environment
python
import scvelo as scv
import anndata as ad
import os

path="/hpc/home/kailasamms/scratch/scRNASeq_BBN_C57B6/results_velocyto/"
files = sorted(os.listdir(path))

# Load loom files for spliced/unspliced matrices for each sample:
for file in files[0:1]:
  ldata = scv.read(path+file, cache=True)

  # Change barcodes to match Seurat object format
  barcodes = ldata.obs_names.tolist()
  barcodes = [bc.replace(":", "_") for bc in barcodes]
  barcodes = [bc.replace("x", "-1") for bc in barcodes]
  ldata.obs_names = barcodes
  
  # Make variable names unique in case some gene names are duplicated
  ldata.var_names_make_unique()
  
  # View dataframe and its dimensions
  ldata.to_df().head()
  ldata.to_df().shape

for file in files[1:]:
  ldata1 = scv.read(path+file, cache=True)

  # Change barcodes to match Seurat object format
  barcodes = ldata1.obs_names.tolist()
  barcodes = [bc.replace(":", "_") for bc in barcodes]
  barcodes = [bc.replace("x", "-1") for bc in barcodes]
  ldata1.obs_names = barcodes
  
  # Make variable names unique in case some gene names are duplicated
  ldata1.var_names_make_unique()
  
  # View dataframe and its dimensions
  ldata1.to_df().head()
  ldata1.to_df().shape
  
  # Concatenate the loom objects
  ldata = ad.concat([ldata, ldata1])

# Check if number of rows in ldata = sum of number of rows of each loom file
ldata.to_df().shape

# Save the file so we can continue in Seurat
ldata.write(filename="/hpc/home/kailasamms/scratch/scRNASeq_BBN_C57B6/results_scvelo/velocyto.h5ad")

#*****************STEP 3: EXPORT SEURAT OBJECT TO ADATA OBJECT*****************#

# NOTE: While SeuratDisk can do this, it doesn't work well. First, SCT and
# integrated assays have counts for only 3000 genes. So, when h5ad file is 
# created only counts for these 3000 genes are retained while velocyto.h5ad has
# info for all 32000+ genes. Merging seurat.h5ad with velocyto.h5ad will lead to
# loss of count data for almost 29000+ genes which will impact velocity analysis.
# Second, if we remove SCT and integrated assay from Seurat object before saving
# as h5ad, then UMAP and PCA information are lost in the resulting Seurat.h5ad.
# So, UMAP has to be created again in scanpy which will differ from Seurat UMAP.
# Also, perform filtering and normalization using scvelo for velocity analysis
# and use the raw counts from RNA assay while creating h5ad file

# To avoid all this nonsense, refer "Step -1" in the following link
# https://smorabit.github.io/tutorials/8_velocyto/

# Run the for loop in R
for (celltype in c("Epithelial Cells", "Myeloid Cells", "Fibroblasts", "T Cells")){
  
  # Load the integrated seurat object
  integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn",
                                      dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  # Save raw count matrices as mtx as saving in csv takes lot of time and space.
  raw <- integrated_seurat@assays$RNA@counts
  Matrix::writeMM(obj = raw, file = paste0(scvelo_results, celltype, "_raw.mtx"))
  
  # Save the metadata
  obs <- integrated_seurat@meta.data
  write.csv(x = obs, file=paste0(scvelo_results, celltype, "_obs.csv"), quote=F, row.names=F)
  
  # Save the genes
  var <- data.frame("Genes" = rownames(integrated_seurat@assays$RNA@counts))
  write.csv(x = var, file=paste0(scvelo_results, celltype, "_var.csv"), quote=F, row.names=F)
  
  # Save PCA embeddings
  pca <- integrated_seurat@reductions$pca@cell.embeddings
  write.csv(x = pca, file=paste0(scvelo_results, celltype, "_pca.csv"), quote=F, row.names=F)
  
  # Save UMAP embeddings
  umap <- integrated_seurat@reductions$umap@cell.embeddings  
  write.csv(x = umap, file=paste0(scvelo_results, celltype, "_umap.csv"), quote=F, row.names=F)
  
}

#****************STEP 4: PERFORM VELOCITY ANALYSIS USING SCVELO****************#

# NOTE: Perform velocity analysis on individual cell types "ALWAYS"
# (i) import adata (containing count matrices from RNA assay)
# (ii) import ldata (containing spliced, unspliced count matrices)
# (iii) merge adata with ldata (this will retain ONLY common cells)
# (iii) perform scvelo analysis (scvelo has different way of filtering, 
# identifying highly variable genes as compared to Seurat)

# Run the code below in python
import scvelo as scv
import scanpy as sc
import anndata as ad
import pandas as pd
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import os

path = "/hpc/home/kailasamms/scratch/scRNASeq_BBN_C57B6/results_scvelo/"

my_palette = ["#FF1F5B", "#F28522", "#009ADE", "#AF58BA", "#00B000", 
              "#FFC61E", "#808080", "#BF812D", "#000000", "#762A83",  
              "#7FBC41", "#C51B7D", "#D6604D", "#4393C3", "#FFFFBF", 
              "#9E0142", "#E41A1C", "#4DAF4A", "#FF7F00", "#FFFF33", 
              "#A65628", "#F781BF", "#66C2A5", "#FC8D62", "#35978F"]

# Load count matrices for spliced and unspliced RNA
ldata = sc.read_h5ad(path+"velocyto.h5ad")

for celltype in ["Epithelial Cells", "Myeloid Cells", "Fibroblasts", "T Cells"]:
  
# Load count matrix for gene expression
X = io.mmread(path+celltype+"_raw.mtx")

# Read obs and var
obs = pd.read_csv(path+celltype+"_obs.csv")
var = pd.read_csv(path+celltype+"_var.csv")

# Read pca and umap co-ordinates
pca = pd.read_csv(path+celltype+"_pca.csv")
umap = pd.read_csv(path+celltype+"_umap.csv")

# Create the anndata object to store normalized counts
# Notice that output of "integrated_seurat@assays$RNA@counts[1:10,1:5]" and
# "X.todense()[0:10,0:5]" are exactly same.
# Anndata stores cells/obs in rows and genes/vars in columns while Seurat stores
# genes in rows and cells in columns. SO, we transpose the matrices while 
# creating Anndata object
adata = ad.AnnData(
  X = X.transpose().tocsr(),
  obs = obs,
  var = var)

# Although we supplied obs and var while creating Anndata object, we need to 
# manually give row and column names to Anndata object
adata.obs_names = obs.loc[:,"Cell"]
adata.var_names = var.loc[:,"Genes"]

# # Create an anndata object to store raw counts
# raw_adata = anndata.AnnData( 
#   X = raw.transpose().tocsr(),
#   obs = obs,
#   var = var)
# raw_adata.obs_names = obs.loc[:,"Cell"]
# raw_adata.var_names = var.loc[:,"Genes"]
# 
# # Add raw counts to our Anndata object  
# adata.raw = raw_adata

# Set pca and umap
# Refer https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.html
# obsm MUST be ndarray i.e. numpy array
pca.index = obs.loc[:,"Cell"]
umap.index = obs.loc[:,"Cell"]
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = umap.to_numpy()

# Merge ldata and adata matrices
# NOTE: This will ONLY keep only cells present in both adata and ldata
adata = scv.utils.merge(adata, ldata)

scv.pl.proportions(adata, save=celltype+".pdf")

# Pre-process the data
# Preprocessing requisites consist of 
# (i) gene selection by detection (with a minimum number of counts in both 
# spliced & unspliced) and (ii) and high variability (dispersion), 
# (iii) normalizing every cell by its total size and (iv) logarithmizing X. 
# Filtering and normalization is applied in the same vein to spliced/unspliced 
# counts and X. Logarithmizing is only applied to X. 
# If X is already preprocessed from former analysis, it will not be touched.

# scv.pp.filter_and_normalize() is a combination of 4 functions
# (i)   scv.pp.filter_genes(adata, min_shared_counts=20)
# (ii)  scv.pp.normalize_per_cell(adata)
# (iii) scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
# (iv)  scv.pp.log1p(adata) i.e. ln1p()
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# compute velocity
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)

# Save the UMAP
scv.pl.velocity_embedding_stream(adata, 
                                 basis="umap", 
                                 color=['integrated_snn_res.0.4'],
                                 palette=my_palette,
                                 save=path+celltype+"_embedding_stream.pdf")


# Save genes to csv
# pd.DataFrame(adata.var_names.to_list()).to_csv("1.csv")
  
# Subsetting based on barcodes
# mycells = obs.loc[:,"Cell"].to_list()
# ldata1= ldata[ldata.obs_names.isin(mycells)]

# adata.to_df()["Rb1cc1"].sort_values(by=['Rb1cc1']).to_csv("1.csv")
# ldata.to_df()["Rb1cc1"].sort_values().to_csv("2.csv")

# You can set adata.obsm['X_umap'] to the array you have, given that the row 
# order is conserved. Also creating a specific matrix is ok, say 
# adata.obsm['X_seurat'], then use that as basis for all visualizations and
# calculations (e.g. scv.pl.scatter(adata, basis='seurat'…))

# # View embedding
# sc.get.obs_df(adata, obsm_keys=[("X_umap",0), ("X_umap",1)])
# # get embedding  
# emb = pd.DataFrame(adata.obsm['X_umap'])
# clusters = adata.obs.seurat_clusters

######AVOID USING SEURAT DISK################################  
# # Load the integrated seurat object
# integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn.rds"))
# 
# # Read the spliced and unspliced matrices from velocyto.h5ad
# SeuratDisk::Convert(source = paste0(velocyto_results, "velocyto.h5ad"), 
#                     dest = "h5seurat", overwrite = TRUE)
# 
# # Get an idea of what is present in velocyto.h5seurat
# # Refer https://mojaveazure.github.io/seurat-disk/articles/h5Seurat-load.html
# hfile <- Connect(paste0(velocyto_results, "velocyto.h5seurat"), mode="r")
# hfile$index()
# hfile$close_all()
# 
# # We need ONLY the count matrices for spliced and unspliced RNA
# velocyto <- SeuratDisk::LoadH5Seurat(file = paste0(velocyto_results, "velocyto.h5seurat"),
#                                       assays = list(spliced = c("counts"),
#                                                     unspliced = c("counts")))
# 
# # Remove cells which are absent in Seurat object from velocyto object
# velocyto@meta.data <- velocyto@meta.data %>% dplyr::mutate(Cell = rownames(.))
# velocyto <- subset(x = velocyto,
#                    Cell %in% rownames(integrated_seurat@meta.data))
# 
# # Get counts of spliced and unspliced matrices
# s_data <- GetAssayData(object = velocyto[["spliced"]], slot = "counts")
# u_data <- GetAssayData(object = velocyto[["unspliced"]], slot = "counts")
# 
# # Create a sparse matrix with 0 counts for cells that are present in Seurat 
# # object but absent in velocyto object
# missing_cells <- setdiff(rownames(integrated_seurat@meta.data), rownames(velocyto@meta.data))
# missing_data <- matrix(0, nrow = nrow(s_data), ncol = length(missing_cells))
# colnames(missing_data) <- missing_cells
# rownames(missing_data) <- rownames(s_data)
# missing_data <- as(missing_data, "sparseMatrix")
# 
# # Merge missing data with spliced and unspliced data
# s_data <- cbind(s_data, missing_data)
# u_data <- cbind(u_data, missing_data)
# 
# # Rearrange the column names to match the column names of Seurat object because 
# # the order of cells MUST be same to be able to add data to new assay
# s_data <- s_data[,rownames(integrated_seurat@meta.data)]
# u_data <- u_data[,rownames(integrated_seurat@meta.data)]
# 
# # Create new assays to store the count matrices for spliced and unspliced RNA
# integrated_seurat[["spliced"]] <- CreateAssayObject(counts = s_data)
# integrated_seurat[["unspliced"]] <- CreateAssayObject(counts = u_data)
# 
# # Remove SCT and integrated assays before saving as h5ad. All layers/assays in
# # h5ad files have same number of genes and cells. Since SCT and integrated 
# # assays have fewer number of genes, it will create problems
# integrated_seurat[["SCT"]] <- NULL
# integrated_seurat[["integrated"]] <- NULL
# SeuratDisk::SaveH5Seurat(integrated_seurat, filename = paste0(scvelo_results, "adata.h5Seurat"))
# SeuratDisk::Convert(source = paste0(scvelo_results, "adata.h5Seurat"), dest = "h5ad")

#####AVOID USING RETICULATE, USE PYTHON######################
# # Open a new terminal, type the codes below
# # Install python also in environment containing R package reticulate so that
# # reticulate can access scvelo package
# # Refer https://rstudio.github.io/reticulate/articles/versions.html
# # Refer https://rstudio.github.io/reticulate/
# 
# conda create --name R python=3.9  # skip if you already created R with Python3.9
# conda activate R
# R
# library(reticulate)
# use_condaenv("scVelo", required = TRUE)
# scv <- import("scvelo")
# scv$logging$print_version()
# 
# scv$settings$verbosity = 3
# scv$settings$set_figure_params('scvelo', facecolor='white', dpi=100, frameon=FALSE)


  # plot velocity of a selected gene
  scv.pl.velocity(adata, var_names=['Ptgds'], color='celltype')
  
  scv.tl.latent_time(adata_subset)
  scv.pl.scatter(adata_subset, color='latent_time', color_map='gnuplot', size=80)
  
  top_genes = adata_subset.var['fit_likelihood'].sort_values(ascending=False).index[:300]
  scv.pl.heatmap(adata_subset, var_names=top_genes, sortby='latent_time', col_color='celltype', n_convolve=100)
  
  top_genes = adata_subset.var['fit_likelihood'].sort_values(ascending=False).index
  scv.pl.scatter(adata_subset, color='celltype', basis=top_genes[:15], ncols=5, frameon=False)
  
  scv.tl.rank_velocity_genes(adata, groupby='celltype', min_corr=.3)
  
  df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
  df.head()  

# Plot same fig using R


rownames(emb) <- names(seurat_clusters) <- adata.obs_names.values

## get clusters, convert to colors
col <- rainbow(length(levels(clusters)), s=0.8, v=0.8)
cell.cols <- col[clusters]
names(cell.cols) <- names(clusters)

## simple plot
plot(emb, col=cell.cols, pch=16,
     xlab='UMAP X', ylab='UMAP Y')
legend(x=-13, y=0, 
       legend=levels(clusters),
       col=col, 
       pch=16)

## top dynamic genes
topgenes <- adata$var["fit_likelihood"]
topgenes_vals <- topgenes[,1]
names(topgenes_vals) <- rownames(topgenes)
topgenes_vals <- sort(topgenes_vals, decreasing=TRUE)
head(topgenes_vals)

scv$pl$scatter(adata, basis=names(topgenes_vals)[1:5], ncols=5, frameon=FALSE)