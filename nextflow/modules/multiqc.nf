// =============================================================================
// PROCESS: MULTIQC
// =============================================================================
// Purpose: Aggregates QC reports from all tools into a single interactive HTML report
//
// What it does:
//   - Scans work directory recursively for recognized output files
//   - Parses metrics from whichever tools ran upstream (e.g. FastQC, STAR,
//     Salmon, RSeQC for rnaseq.nf; FastQC, Cell Ranger for scrnaseq.nf)
//   - Generates interactive HTML report and data directory
//
// This module is shared between rnaseq.nf and scrnaseq.nf - the calling
// pipeline is responsible for assembling/flattening the right set of report
// files/dirs into all_reports before invoking this process.
//
// Typical resources: <1GB RAM, 1-5 minutes
// =============================================================================

process MULTIQC {

    tag "MultiQC: ${species}"
    label 'process_low'                           // Lightweight file parsing only

    publishDir = [
        [path: { "${params.proj_dir()}/${species}/12.MultiQC" }, mode: 'copy', pattern: "*.html"],
        [path: { "${params.proj_dir()}/${species}/12.MultiQC" }, mode: 'copy', pattern: "*_data"],
        [path: { "${params.proj_dir()}/${species}/Logs" },      mode: 'copy', pattern: "MULTIQC.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), path(all_reports)
    // all_reports: Flat list of all QC output files/dirs (FastQC + upstream tool
    //              outputs). Assembled and flattened by the calling pipeline
    //              (rnaseq.nf or scrnaseq.nf) via .mix().groupTuple().map{flatten}

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("*.html"),                    emit: multiqc_html     // Interactive HTML report
    path("*_data"),                    emit: multiqc_data     // Parsed data directory
    path("MULTIQC.log"),         emit: error_log        // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG              = "MULTIQC.log"
    def multiqc_title    = "${species} ${params.project} QC Report"
    def multiqc_filename = "${species}_${params.project}_multiqc"

    """
    # Run MultiQC on current directory (Nextflow stages all input files here)
    # Using "." avoids ARG_MAX errors with large sample counts and handles
    # directory inputs (e.g. Salmon quant dirs) automatically
    #
    # --force    : Overwrite existing reports
    # --clean-up : Remove intermediate files after completion
    # --quiet    : Suppress progress messages
    multiqc \
        --force \
        --clean-up \
        --quiet \
        --title "${multiqc_title}" \
        --filename "${multiqc_filename}" \
        . \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: MultiQC failed for ${species}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: MultiQC report generated for ${species}" >> "${LOG}"
    """
}

// =============================================================================
// QUICK REFERENCE
// =============================================================================
//
// Tools automatically detected and parsed (depends on calling pipeline):
//   FastQC      : *_fastqc.zip                          (rnaseq.nf, scrnaseq.nf)
//   STAR        : *.Log.final.out                       (rnaseq.nf)
//   Salmon      : meta_info.json + quant.sf              (rnaseq.nf)
//   RSeQC       : *.read_distribution.txt, *.junction*   (rnaseq.nf)
//   Cell Ranger : web_summary.html, metrics_summary.csv  (scrnaseq.nf)
//
// Output files:
//   *_multiqc.html        : Main interactive report (open in browser)
//   *_multiqc_data/       : Parsed data directory:
//     multiqc_data.json   : All data in JSON
//     multiqc_general_stats.txt : Summary table (TSV)
//     multiqc_sources.txt : List of all parsed files
//
// Common issues:
//   - "No results found"    → Check file patterns; update MultiQC version
//   - Missing samples       → Check for special characters in sample names
//   - Report > 50MB         → Many samples; expected behavior
// =============================================================================
