process RNA_DEEP_AUDIT {

    tag "Auditing results"
    label 'process_medium'                      // STAR indexing requires 30-50GB RAM for human
    label 'omics_r'

    publishDir = [
        [path: { "${params.proj_dir()}" }, mode: 'copy', pattern: "*.{csv,pdf}"]
        // [path: { "${params.proj_dir()}/${species}/Logs" }, mode: 'copy', pattern: "DEEP_AUDIT.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    // We don't want to stage the whole project dir (too slow).
    // Instead, we just pass it as a value to the bash scripts.
    val(proj_dir)
    val trigger  // This is the "Gatekeeper"

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("*.{csv,pdf}"),             emit: summary_files
    path("DEEP_AUDIT.log"),    emit: error_log    // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    // This points to the modules folder relative to your project root
    def script_path = "${workflow.projectDir}/modules/RNA"
    def LOG = "DEEP_AUDIT.log"

    """
    Rscript ${script_path}/RNA_DEEP_AUDIT.R "${proj_dir}" > ${LOG} 2>&1
    """
}