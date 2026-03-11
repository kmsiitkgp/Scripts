process DEEP_AUDIT {

    tag "Auditing results"
    label 'process_medium'                      // STAR indexing requires 30-50GB RAM for human

    publishDir { "${params.proj_dir()}" },    mode: 'copy',    pattern: "*.{csv,pdf}"
    //publishDir { "${params.proj_dir()}/${species}/07.Logs" },      mode: 'copy',    pattern: "DEEP_AUDIT.error.log"

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
    path("DEEP_AUDIT.error.log"),    emit: error_log    // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    // This points to the modules folder relative to your project root
    def script_path = "${workflow.projectDir}/modules"
    def LOG = "DEEP_AUDIT.error.log"

    """
    Rscript ${script_path}/DEEP_AUDIT.R "${proj_dir}" > ${LOG} 2>&1
    """
}