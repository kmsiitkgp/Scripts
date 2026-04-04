// =============================================================================
// PROCESS: CELLRANGER_COUNT
// =============================================================================
// Purpose: Gene expression quantification for single-cell RNA-seq using Cell Ranger
//
// What it does:
//   - Aligns reads to reference transcriptome using STAR (internal to Cell Ranger)
//   - Filters cells from background (empty droplets)
//   - Generates UMI count matrix per gene per cell barcode
//   - Produces filtered/raw feature-barcode matrices and QC web summary
//
// Note: This process is standalone (not part of the bulk RNA-seq workflow)
//
// Typical resources: 64GB RAM, 8-16 cores, 30-90 minutes per sample
// =============================================================================

process CELLRANGER_COUNT {

    tag "Cell Ranger count: ${sample_id}"
    label 'process_high'                          // 64GB RAM required for human reference

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    val(sample_id)                 // Sample identifier (must match FASTQ filename prefix)
    path(raw_fastq_dir)            // Directory containing raw FASTQ files
    path(cellranger_index_dir)     // Pre-built Cell Ranger reference (from cellranger mkref)
    val(cellranger_args)           // Pre-joined Cell Ranger argument string from workflow
    // Never call params.CELLRANGER_ARGS().join(' ') inside a process — it changes the
    // task hash on every run and breaks -resume caching

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("${sample_id}/outs/**"),                        emit: result_files    // All Cell Ranger output files
    path("${sample_id}.CELLRANGER_COUNT.error.log"),     emit: error_log       // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG = "${sample_id}.CELLRANGER_COUNT.error.log"

    """
    # Run Cell Ranger count
    # --transcriptome : Reference directory (from cellranger mkref)
    # --fastqs        : Directory containing FASTQ files
    # --sample        : Sample name prefix (must match FASTQ filenames)
    # --id            : Output directory name
    # --localcores    : CPU cores
    # cellranger_args (from config) typically includes:
    #   --chemistry auto       : Auto-detect 10x chemistry version
    #   --expect-cells N       : Expected number of recovered cells
    #   --include-introns true : Count intronic reads (needed for nuclei/nascent RNA)
    cellranger count \
        ${cellranger_args} \
        --transcriptome "${cellranger_index_dir}" \
        --fastqs "${raw_fastq_dir}" \
        --sample "${sample_id}" \
        --id "${sample_id}" \
        --localcores "${task.cpus}" \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: Cell Ranger count failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: Cell Ranger count completed for ${sample_id}" >> "${LOG}"
    """
}

// =============================================================================
// QUICK REFERENCE
// =============================================================================
//
// Output directory structure (${sample_id}/outs/):
//   filtered_feature_bc_matrix/ : Main count matrix (cells passing filtering)
//     barcodes.tsv.gz            : Cell barcodes
//     features.tsv.gz            : Gene IDs and names
//     matrix.mtx.gz              : Sparse count matrix
//   raw_feature_bc_matrix/      : Unfiltered matrix (all barcodes)
//   web_summary.html            : Interactive QC report
//   metrics_summary.csv         : Key metrics in tabular format
//   cloupe.cloupe               : Loupe Browser visualization file
//
// FASTQ naming requirement (Illumina format):
//   [Sample]_S[#]_L[Lane]_[R1/R2/I1/I2]_001.fastq.gz
//   --sample must match the [Sample] prefix exactly
//
// Expected metrics (good quality run):
//   Cells detected              : 1,000-10,000 (depends on loading)
//   Median genes per cell       : >500
//   Reads mapped to transcriptome: >50%
//   Sequencing saturation       : >50%
//
// Common issues:
//   - Low cells detected   → Poor cell viability or too few cells loaded
//   - Low genes/cell       → Low sequencing depth or poor RNA quality
//   - Sample name mismatch → --sample must match FASTQ filename prefix exactly
//   - High mitochondrial % → Cell stress or death; filter in downstream Seurat analysis
// =============================================================================
