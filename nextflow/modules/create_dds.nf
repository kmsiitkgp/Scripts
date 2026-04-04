// =============================================================================
// PROCESS: CREATE_DDS
// =============================================================================
// Purpose: Creates a DESeq2 dataset object (DESeqDataSet) from txi and metadata
//
// What it does:
//   - Reads txi object (from CREATE_TXI) and sample metadata
//   - Constructs DESeqDataSet using the specified design formula
//   - Runs DESeq2 normalization and dispersion estimation
//   - Saves dds object as RDS for downstream contrast extraction
//   - Also exports normalized count matrix as Excel
// =============================================================================

process CREATE_DDS {

    tag "Creating DESeq2 dataset: ${species}"
    label 'process_medium'                        // R with DESeq2; ~8GB RAM

    publishDir { "${params.proj_dir()}/${species}/08.DESeq2" },    mode: 'copy',    pattern: "*.{rds,xlsx}"
    publishDir { "${params.proj_dir()}/${species}/07.Logs" },      mode: 'copy',    pattern: "CREATE_DDS.error.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), path(txi), path(tx2gene)
    // species  : "Human" or "Mouse"
    // txi      : Salmon_txi.rds from CREATE_TXI
    // tx2gene  : tx2gene CSV passed through for downstream use
    path(metadata)    // Sample metadata Excel (columns must match design formula variables)
    val(design)       // DESeq2 design formula string (e.g., "~condition" or "~batch+condition")

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), path("DESeq2_dds.rds"), path(tx2gene),    emit: dds          // DESeqDataSet + tx2gene for EXTRACT_DESEQ2_RESULTS
    tuple val(species), path("Norm_Counts*.xlsx"),                emit: norm_counts       // Count matrices
    tuple val(species), path("Log_Norm*.xlsx"),                   emit: log_norm_counts     // DEGs file for this contrast
    tuple val(species), path("VST_Blind*.xlsx"),                  emit: vst_blind_counts
    path(metadata),                                               emit: metadata
    path("CREATE_DDS.error.log"),                                 emit: error_log    // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def script_path = "${workflow.projectDir}/modules"
    def LOG         = "CREATE_DDS.error.log"

    """
    # Arg 1: txi      — path to Salmon_txi.rds
    # Arg 2: metadata — path to sample metadata Excel file
    # Arg 3: tx2gene  — path to tx2gene CSV
    # Arg 4: "."      — output directory (current work dir)
    # Arg 5: design   — DESeq2 design formula string
    Rscript "${script_path}/CREATE_DDS.R" \
        "${txi}" \
        "${metadata}" \
        "${tx2gene}" \
        "." \
        "${design}" > "${LOG}" 2>&1 \
        || { echo "❌ ERROR: CREATE_DDS.R failed for ${species}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: DESeq2 dataset created for ${species}" >> "${LOG}"
    """
}
