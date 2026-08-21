process SALMON_STATS {

    tag "Salmon stats"
    label 'process_low'

    publishDir { "${params.proj_dir()}" }, mode: 'copy', pattern: "Salmon_Master_Summary.txt"

    // =================================================================================
    // INPUT
    // =================================================================================
    // proj_dir  : project root — passed as a value so the bash script can crawl it.
    // trigger   : collected channel of all SALMON_QUANT output dirs — acts as a
    //             gatekeeper ensuring this process only runs after ALL samples finish.
    //             The value itself is not used inside the script; its presence in the
    //             channel forces Nextflow to wait for all upstream items to complete.

    input:
    val(proj_dir)
    val(trigger)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("Salmon_Master_Summary.txt"), emit: salmon_summary

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:
    def script_path = "${workflow.projectDir}/modules"
    def LOG         = "SALMON_STATS.log"

    """
    bash ${script_path}/salmon_stats_collector.sh "${proj_dir}" > ${LOG} 2>&1
    """
}
