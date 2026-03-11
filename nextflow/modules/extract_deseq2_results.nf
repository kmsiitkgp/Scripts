process EXTRACT_DESEQ2_RESULTS {

    tag "Extracting results: ${species} ${contrast}"
    label 'process_medium'                      // STAR indexing requires 30-50GB RAM for human

    // We define the logic in a closure so it is evaluated when the block is called
    // Rscript creates the folder using safe_contrast internally
    def get_safe = { c -> c.replaceAll(/[^a-zA-Z0-9-]/, '_') }

    publishDir { "${params.proj_dir()}/${species}/08.DESeq2" },    mode: 'copy',    pattern: "${get_safe(contrast)}"
    publishDir { "${params.proj_dir()}/${species}/07.Logs" },      mode: 'copy',    pattern: "*.error.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), path(dds), path(tx2gene), val(contrast)
    val(deseq2)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), val(contrast), path("${get_safe(contrast)}"),    emit: deg_dir
    path("*.error.log"),                                                 emit: error_log    // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================

    script:

    // Here we can use a local variable for the bash part
    def safe_contrast = get_safe(contrast)

    // This points to the modules folder relative to your project root
    def script_path = "${workflow.projectDir}/modules"
    def LOG = "${species}_${safe_contrast}_EXTRACT_DESEQ2_RESULTS.error.log"

    // --- THIS HANDLES ALL 3 CASES ---
    // 1. If null -> returns empty string ""
    // 2. If "batch" -> returns "batch"
    // 3. If ["batch", "condition"] -> returns "batch,condition"
    def batch_str = (deseq2.batch_vars instanceof List) ? deseq2.batch_vars.join(',') : (deseq2.batch_vars ?: "")

    """
    # We pass '.' as the 4th arg so the Excel file is saved in the current folder
    Rscript ${script_path}/EXTRACT_DESEQ2_RESULTS.R \
        "${contrast}" \
        "${dds}" \
        "${tx2gene}" \
        "." \
        "${batch_str}" \
        "${deseq2.lfc_cutoff}" \
        "${deseq2.padj_cutoff}" > ${LOG} 2>&1
    """
}