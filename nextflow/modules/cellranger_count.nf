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

process CELLRANGER_COUNT {

    tag "Quantifying ${sample_id} with Cell Ranger"
    label 'process_high'                            // Cell Ranger requires 64GB RAM for human

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    val(sample_id)                                   // Sample identifier (must match FASTQ prefix)
    path(raw_fastq_dir)                              // Directory containing raw FASTQ files
    path(cellranger_index_dir)                       // Pre-built Cell Ranger reference directory
    val(cellranger_args)                             // Pre-joined Cell Ranger arguments from config
    // Never do ${params.CELLRANGER_ARGS().join(' ')} inside process. It changes hash on
    // every run and resume fails. So, .join() in main.nf and pass as an argument.

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("${sample_id}/outs/raw_feature_bc_matrix"),        emit: raw_matrix_dir       // Unfiltered matrix (all barcodes)
    path("${sample_id}/outs/filtered_feature_bc_matrix"),   emit: filt_matrix_dir      // Filtered matrix (called cells only)
    path("${sample_id}/outs/web_summary.html"),             emit: web_summary          // Interactive QC report
    path("${sample_id}/outs/metrics_summary.csv"),          emit: metric_summary       // Key metrics (CSV format)
    path("${sample_id}/outs/*.h5"),                         emit: h5_files             // HDF5 format matrices
    path("${sample_id}/outs/*.cloupe"),                     emit: cloupe               // Loupe Browser file
    path("${sample_id}/outs/*.json"),                       emit: run_info,            optional: true    // Run metadata
    path("${sample_id}/outs/*.bam*"),                       emit: bam_files,           optional: true    // Aligned BAM + index
    path("${sample_id}/outs/analysis"),                     emit: analysis_dir,        optional: true    // PCA, t-SNE, clustering
    path("${sample_id}.CELLRANGER_COUNT.error.log"),        emit: error_log                              // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG = "${sample_id}.CELLRANGER_COUNT.error.log"

    """
    # Run Cell Ranger count
    # --transcriptome: Reference transcriptome directory (from cellranger mkref)
    # --fastqs: Directory containing FASTQ files
    # --sample: Sample name (must match FASTQ prefix, e.g., Sample1_S1_L001_*)
    # --id: Output directory name (same as sample for consistency)
    # --localcores: Number of CPU cores to use
    # Additional parameters from config (cellranger_args):
    #   --chemistry: Auto-detect or specify (SC3Pv3, SC5P-PE, etc.)
    #   --expect-cells: Expected number of recovered cells
    #   --nosecondary: Skip PCA/clustering (faster, for QC only)
    #   --include-introns: Count intronic reads (nascent RNA, nuclei)

    cellranger count \\
        ${cellranger_args} \\
        --transcriptome "${cellranger_index_dir}" \\
        --fastqs "${raw_fastq_dir}" \\
        --sample "${sample_id}" \\
        --id "${sample_id}" \\
        --localcores "${task.cpus}" \\
        1>> "${LOG}" 2>&1 \\
        || { echo "❌ ERROR: Cell Ranger count failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: Cell Ranger count completed for ${sample_id}" >> "${LOG}"
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
