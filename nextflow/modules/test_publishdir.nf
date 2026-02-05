// =========================================================================================
// PROCESS: CELLRANGER_COUNT
// =========================================================================================
// Purpose: Gene expression quantification for single-cell RNA-seq using Cell Ranger
//
// What it does:
//   - Demultiplexes raw BCL files or processes FASTQ files
//   - Aligns reads to reference transcriptome using STAR
//   - Filters cells from background (empty droplets)
//   - Generates UMI counts per gene per cell
//   - Produces feature-barcode matrices and quality metrics
//   - Creates interactive web summary report
//
// Cell Ranger advantages:
//   - Handles 10x Genomics chemistry (3' and 5' gene expression)
//   - Automatic cell calling and doublet detection
//   - Built-in secondary analysis (PCA, t-SNE, clustering)
//   - Compatible with Loupe Browser for visualization
//   - Generates both raw and filtered count matrices
//
// Typical resources: 8-16 cores, 64GB RAM, 30-90 minutes per sample
//
// For detailed explanation, see: docs/cellranger_count.md
// =========================================================================================

process TEST_PUBLISHDIR {

    tag "Testing publishDir"
    label 'process_low'                            // Cell Ranger requires 64GB RAM for human

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    path(sample_dir)
    path(log)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("*/outs/**"), emit: sample_dir, includeInputs: true
    path("*.log"),        emit: error_log,  includeInputs: true

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:


    """
    echo "Hi"
    """
}

// =========================================================================================
// QUICK REFERENCE
// =========================================================================================
//
// Output directory structure (${sample_id}/outs/):
//   filtered_feature_bc_matrix/: Main count matrix (genes × cells)
//     - barcodes.tsv.gz: Cell barcodes
//     - features.tsv.gz: Gene IDs and names
//     - matrix.mtx.gz: Sparse count matrix
//   raw_feature_bc_matrix/: Unfiltered matrix (includes empty droplets)
//   web_summary.html: QC report with plots and metrics
//   metrics_summary.csv: Key metrics in tabular format
//   cloupe.cloupe: File for Loupe Browser visualization
//   possorted_genome_bam.bam: Aligned reads (optional, if --no-bam not used)
//   analysis/: Secondary analysis (PCA, t-SNE, UMAP, clustering)
//
// FASTQ naming requirements:
//   Cell Ranger expects Illumina naming: [Sample]_S[#]_L[Lane]_[R1/R2/I1/I2]_001.fastq.gz
//   Example: Sample1_S1_L001_R1_001.fastq.gz, Sample1_S1_L001_R2_001.fastq.gz
//   --sample parameter must match the [Sample] prefix
//
// Expected metrics (good quality run):
//   Cells detected: 1,000-10,000 (depends on loading)
//   Median genes per cell: >500 (higher is better)
//   Median UMI counts per cell: >1,000 (depends on depth)
//   Sequencing saturation: >50% (>70% ideal)
//   Reads mapped to genome: >70%
//   Reads mapped confidently to transcriptome: >50%
//
// Common chemistry types:
//   SC3Pv3: 10x Chromium Single Cell 3' v3 (most common)
//   SC3Pv2: 10x Chromium Single Cell 3' v2
//   SC5P-PE: 10x Chromium Single Cell 5' paired-end
//   SC5P-R2: 10x Chromium Single Cell 5' R2-only
//   ARC-v1: 10x Chromium Single Cell Multiome ATAC + Gene Expression
//
// Downstream analysis (R/Seurat):
//   ```r
//   library(Seurat)
//   library(Matrix)
//
//   # Load 10x data
//   data <- Read10X(data.dir = "sample/outs/filtered_feature_bc_matrix")
//   seurat_obj <- CreateSeuratObject(counts = data, project = "sample")
//
//   # Standard workflow
//   seurat_obj <- NormalizeData(seurat_obj)
//   seurat_obj <- FindVariableFeatures(seurat_obj)
//   seurat_obj <- ScaleData(seurat_obj)
//   seurat_obj <- RunPCA(seurat_obj)
//   seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
//   seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
//   seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
//   ```
//
// Common issues:
//   - Low cells detected (<500) → Poor cell viability, too few cells loaded
//   - Low genes/cell (<300) → Low sequencing depth, poor RNA quality
//   - Low mapping rate (<60%) → Wrong reference, contamination
//   - High mitochondrial % (>20%) → Cell stress/death, poor sample quality
//   - Sample name mismatch → Check FASTQ prefix matches --sample parameter
//
// For comprehensive guide and troubleshooting, see: docs/cellranger_count.md
// =========================================================================================
