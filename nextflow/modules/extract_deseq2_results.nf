// =============================================================================
// PROCESS: EXTRACT_DESEQ2_RESULTS
// =============================================================================
// Purpose: Extracts DESeq2 differential expression results for a single contrast
//
// What it does:
//   - Runs DESeq2 results() for the specified contrast
//   - Applies LFC shrinkage (apeglm)
//   - Filters DEGs by padj and LFC cutoffs (from config)
//   - Handles optional batch correction variables
//   - Saves results as Excel and CSV under a contrast-named subdirectory
//
// Called once per contrast (fan-out via .combine() in rnaseq.nf)
// =============================================================================

process EXTRACT_DESEQ2_RESULTS {

    tag "DESeq2 results: ${species} / ${contrast}"
    label 'process_medium'                        // R with DESeq2; ~8GB RAM

    // get_safe sanitizes the contrast label for use as a directory name
    // e.g. "Treated vs Control" → "Treated_vs_Control"
    def get_safe = { c -> c.replaceAll(/[^a-zA-Z0-9-]/, '_') }

    publishDir { "${params.proj_dir()}/${species}/08.DESeq2" },    mode: 'copy',    pattern: "${get_safe(contrast)}/*.{rds,xlsx}"
    publishDir { "${params.proj_dir()}/${species}/07.Logs" },      mode: 'copy',    pattern: "*.error.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), path(dds), path(tx2gene), val(contrast)
    // species   : "Human" or "Mouse"
    // dds       : DESeq2_dds.rds from CREATE_DDS
    // tx2gene   : tx2gene CSV for gene symbol annotation
    // contrast  : Contrast string (e.g., "condition_Treated_vs_Control")
    val(deseq2)    // DESeq2 params map from config (lfc_cutoff, padj_cutoff, batch_vars)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), val(contrast), path("${get_safe(contrast)}/*.rds"),                 emit: res
    tuple val(species), val(contrast), path("${get_safe(contrast)}/DEGs.xlsx"),             emit: degs         // DEGs file for this contrast
    tuple val(species), val(contrast), path("${get_safe(contrast)}/VST_NonBlind*.xlsx"),    emit: vst
    path("*.error.log"),                                                                    emit: error_log   // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def safe_contrast = get_safe(contrast)
    def script_path   = "${workflow.projectDir}/modules"
    def LOG           = "${species}_${safe_contrast}_EXTRACT_DESEQ2_RESULTS.error.log"

    // Handles all 3 batch_vars cases:
    //   null                   → "" (no batch correction)
    //   "batch"                → "batch"
    //   ["batch", "condition"] → "batch,condition"
    def batch_str = (deseq2.batch_vars instanceof List) \
        ? deseq2.batch_vars.join(',') \
        : (deseq2.batch_vars ?: "")

    """
    # Arg 1: contrast    — contrast string for DESeq2 results()
    # Arg 2: dds         — path to DESeq2_dds.rds
    # Arg 3: tx2gene     — path to tx2gene CSV
    # Arg 4: "."         — output directory (current work dir; R creates safe_contrast subdir internally)
    # Arg 5: batch_str   — comma-separated batch variables (or "" if none)
    # Arg 6: lfc_cutoff  — log2 fold change threshold for DEG filtering
    # Arg 7: padj_cutoff — adjusted p-value threshold for DEG filtering
    Rscript "${script_path}/EXTRACT_DESEQ2_RESULTS.R" \
        "${contrast}" \
        "${dds}" \
        "${tx2gene}" \
        "." \
        "${batch_str}" \
        "${deseq2.lfc_cutoff}" \
        "${deseq2.padj_cutoff}" > "${LOG}" 2>&1 \
        || { echo "❌ ERROR: EXTRACT_DESEQ2_RESULTS.R failed for ${species} / ${contrast}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: DESeq2 results extracted for ${species} / ${contrast}" >> "${LOG}"
    """
}
