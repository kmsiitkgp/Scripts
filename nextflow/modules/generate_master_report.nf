process GENERATE_MASTER_REPORT {

    tag "Generating Master Report"
    label 'process_low'                            // STAR requires 30-50GB RAM for human

    publishDir { "${params.proj_dir()}" },         mode: 'copy',    pattern: "*Master_Summary.txt"
    publishDir { "${params.proj_dir()}/Logs" },    mode: 'copy',    pattern: "*.error.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    // We don't want to stage the whole project dir (too slow).
    // Instead, we just pass it as a value to the bash scripts.
    val(proj_dir)
    val(trigger)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("STAR_Master_Summary.txt"),             emit: star_stats                 // Alignment stats
    path("Salmon_Master_Summary.txt"),           emit: salmon_stats               // Alignment stats
    path("GENERATE_MASTER_REPORT.error.log"),    emit: error_log,    optional: true // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    // This points to the modules folder relative to your project root
    def script_path = "${workflow.projectDir}/modules"
    def LOG = "GENERATE_MASTER_REPORT.error.log"

    """

    bash ${script_path}/star_stats_collector.sh "${proj_dir}" > ${LOG} 2>&1
    bash ${script_path}/salmon_stats_collector.sh "${proj_dir}" >> ${LOG} 2>&1
    """
}