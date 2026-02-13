// =========================================================================================
// PROCESS: FASTQC
// =========================================================================================
// Purpose: Quality control analysis on FASTQ files using FastQC
//
// What it does:
//   - Analyzes read quality, GC content, adapter contamination, duplication levels
//   - Generates HTML reports (human-readable) and ZIP files (for MultiQC)
//   - Runs on both raw and trimmed reads to compare before/after
//
// When to use:
//   - Always run before alignment to catch sequencing issues
//   - Run after trimming to verify improvement
//
// For detailed explanation, see: docs/fastqc.md
// =========================================================================================

process FASTQC {

    tag { "FastQC on ${species} sample ${sample_id} (${read_type})" }
    label 'process_low'                                      // FastQC needs ~250MB RAM, 1-2 CPUs

    // =================================================================================
    // INPUT
    // =================================================================================

    publishDir { "${params.proj_dir()}/${species}_${type}/02.FastQC/${read_type}" },    mode: 'copy',    pattern: "*.html"
    publishDir { "${params.proj_dir()}/${species}_${type}/02.FastQC/${read_type}" },    mode: 'copy',    pattern: "*.zip"
    publishDir { "${params.proj_dir()}/${species}_${type}/07.Logs" },                   mode: 'copy',    pattern: "*.FASTQC.error.log"


    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), val(type), val(sample_id), path(fastq_files), val(read_type)
    // species      : Human or Mouse
    // type         : full or split
    // sample_id    : Sample identifier (e.g., "Sample1")
    // fastq_files  : List of FASTQ files [R1.fq.gz] for SE or [R1.fq.gz, R2.fq.gz] for PE
    // read_type    : "raw" or "trimmed" (determines output directory)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), val(type), path("*_fastqc.zip"), path("*_fastqc.html"),    emit: fastqc_results        // Data for MultiQC aggregation
    path("*.FASTQC.error.log"),                                         emit: error_log         // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG = "${sample_id}.${read_type}.FASTQC.error.log"

    """
    # Run FastQC on all FASTQ files
    # --threads: Use multiple cores for parallel processing of multiple files
    # --quiet: Suppress progress output (cleaner logs)
    # ${fastq_files.join(' ')}: Converts list to space-separated string
    #   SE example: "Sample1_R1.fq.gz"
    #   PE example: "Sample1_R1.fq.gz Sample1_R2.fq.gz"

    fastqc \
        --threads "${task.cpus}" \
        --quiet \
        ${fastq_files.join(' ')} \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: FastQC failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: FastQC completed for ${sample_id}" >> "${LOG}"
    """
}

// =========================================================================================
// QUICK REFERENCE
// =========================================================================================
//
// Output files per FASTQ:
//   - *_fastqc.html: Visual report with plots
//   - *_fastqc.zip: Raw data (parsed by MultiQC)
//
// Typical runtime: 1-2 minutes per FASTQ file
// Memory usage: ~250MB per file
//
// Common issues:
//   - Corrupted FASTQ → Check with: gzip -t file.fq.gz
//   - No output → Check LOG file and FASTQ format
//   - High duplication in RNA-seq → Normal! (highly expressed genes)
//
// For interpretation guide and troubleshooting, see: docs/fastqc.md
// ========================================================================================='