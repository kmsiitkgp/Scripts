// =========================================================================================
// PROCESS: MULTIQC
// =========================================================================================
// Purpose: Aggregates QC reports from multiple tools into single interactive HTML report
//
// What it does:
//   - Scans for output files from FastQC, STAR, Salmon, RSeQC
//   - Parses metrics from each tool
//   - Aggregates data across all samples
//   - Generates interactive HTML report + data directory
//
// Key features:
//   - Cross-sample comparison (identify outliers)
//   - Interactive plots (hover, zoom, export)
//   - Summary tables (sortable, filterable)
//   - Export options (PNG, SVG, CSV)
//
// Typical resources: <1GB RAM, ~1-5 minutes
//
// For detailed explanation, see: docs/multiqc.md
// =========================================================================================

process MULTIQC {

    tag "Aggregating all QC reports"
    label 'process_low'                                              // Lightweight (file parsing only)

    publishDir { "${params.proj_dir()}/${species}_${type}/06.MULTIQC" },    mode: 'copy',    pattern: "*.html"
    publishDir { "${params.proj_dir()}/${species}_${type}/06.MULTIQC" },    mode: 'copy',    pattern: "*_data"
    publishDir { "${params.proj_dir()}/${species}_${type}/07.Logs" },       mode: 'copy',    pattern: "*.MULTIQC.error.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), val(type), path(all_reports)        // All QC files from upstream processes

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("*.html"),    emit: multiqc_html  // HTML report
    path("*_data"),    emit: multiqc_dir   // Data directory
    path("MULTIQC.error.log"),                            emit: error_log     // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG = "MULTIQC.error.log"
    def multiqc_title    = "${species} ${type} ${params.project} QC Report"
    def multiqc_filename = "${species}_${type}_${params.project}_multiqc"

    """
    # Run MultiQC to consolidate all QC reports
    # Searches current directory recursively for recognized output files
    #
    # Parameters:
    #   --force: Overwrite existing reports
    #   --clean-up: Remove intermediate files after completion
    #   --quiet: Suppress progress messages (cleaner logs)
    #   --title: Report header title
    #   --filename: Output filename (without .html extension)
    #   . (dot): Search current directory recursively
    #
    # Why "." instead of listing files explicitly?
    #   - Avoids ARG_MAX errors (command line too long)
    #   - Handles directories (e.g., Salmon output) automatically
    #   - MultiQC designed for recursive directory search
    #   - Nextflow stages all input files in working directory

    multiqc \
        --force \
        --clean-up \
        --quiet \
        --title "${multiqc_title}" \
        --filename "${multiqc_filename}" \
        . \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: MultiQC failed" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: MultiQC completed" >> "${LOG}"
    """
}

// =========================================================================================
// QUICK REFERENCE
// =========================================================================================
//
// Tools detected automatically:
//   - FastQC: *_fastqc.zip files
//   - STAR: Log.final.out files
//   - Salmon: meta_info.json + quant.sf
//   - RSeQC: Various *.txt, *.log files
//
// Output files:
//   - ${project}_report.html: Main interactive report
//   - ${project}_report_data/: Data directory with:
//     * multiqc_data.json (all data, JSON format)
//     * multiqc_general_stats.txt (summary table, TSV)
//     * multiqc_sources.txt (list of parsed files)
//     * Tool-specific data files
//
// Report sections:
//   1. General Statistics (cross-tool summary)
//   2. FastQC (quality, adapters, duplication)
//   3. STAR (alignment rates, junctions)
//   4. Salmon (mapping rates, library type)
//   5. RSeQC (read distribution, gene body coverage)
//
// Interactive features:
//   - Hover for exact values
//   - Click to highlight sample across all plots
//   - Export plots (PNG/SVG) and data (CSV)
//   - Filter/rename samples
//
// Common issues:
//   - "No results found" → Check file patterns, update MultiQC
//   - Missing samples → Check special characters in names
//   - Large file (>50MB) → Many samples or high-res plots
//
// For comprehensive guide, see: docs/multiqc.md
// ========================================================================================='