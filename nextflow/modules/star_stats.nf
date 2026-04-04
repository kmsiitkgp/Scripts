process STAR_STATS {

    tag "STAR stats"
    label 'process_low'

    publishDir { "${params.proj_dir()}" }, mode: 'copy', pattern: "STAR_Master_Summary.txt"

    // =================================================================================
    // INPUT
    // =================================================================================
    // proj_dir  : project root — passed as a value so the bash script can crawl it.
    // trigger   : collected channel of all STAR_ALIGN outputs — acts as a gatekeeper
    //             ensuring this process only runs after ALL samples finish aligning.

    input:
    val(proj_dir)
    val(trigger)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("STAR_Master_Summary.txt"), emit: star_summary

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:
    def script_path = "${workflow.projectDir}/modules"
    def LOG         = "STAR_STATS.error.log"

    """
    bash ${script_path}/star_stats_collector.sh "${proj_dir}" > ${LOG} 2>&1
    """
}
