// =============================================================================
// PROCESS: GENERATE_MASTER_REPORT
// =============================================================================
// Purpose: Collects alignment and quantification stats across all samples
//          into summary text files
//
// What it does:
//   - Scans the project output directory for STAR Log.final.out files
//     and aggregates key alignment metrics into STAR_Master_Summary.txt
//   - Scans for Salmon meta_info.json files and aggregates mapping rates
//     into Salmon_Master_Summary.txt
//
// Why val(proj_dir) instead of path(proj_dir)?
//   - Staging the entire project directory as a path input would be very slow
//   - The bash scripts scan the real published output paths directly
//   - val() passes the directory path as a string without staging
//
// Why val(trigger)?
//   - Ensures this process runs only after all upstream processes have completed
//   - The trigger value itself is not used in the script
// =============================================================================

process GENERATE_MASTER_REPORT {

    tag "Generating master QC summary"
    label 'process_low'                           // Bash file scanning; minimal resources

    publishDir { "${params.proj_dir()}" },         mode: 'copy',    pattern: "*Master_Summary.txt"
    publishDir { "${params.proj_dir()}/Logs" },    mode: 'copy',    pattern: "GENERATE_MASTER_REPORT.error.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    val(proj_dir)    // Project output directory path (passed as string, not staged)
    val(trigger)     // Completion signal from upstream processes (value unused)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("STAR_Master_Summary.txt"),             emit: star_summary      // STAR alignment stats across all samples
    path("Salmon_Master_Summary.txt"),           emit: salmon_summary    // Salmon mapping rates across all samples
    path("GENERATE_MASTER_REPORT.error.log"),    emit: error_log,    optional: true    // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def script_path = "${workflow.projectDir}/modules"
    def LOG         = "GENERATE_MASTER_REPORT.error.log"

    """
    bash "${script_path}/star_stats_collector.sh" "${proj_dir}" \
        > "${LOG}" 2>&1 \
        || { echo "❌ ERROR: star_stats_collector.sh failed" | tee -a "${LOG}" >&2; exit 1; }

    bash "${script_path}/salmon_stats_collector.sh" "${proj_dir}" \
        >> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: salmon_stats_collector.sh failed" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: Master report generated" >> "${LOG}"
    """
}
