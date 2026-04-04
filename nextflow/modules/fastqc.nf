// =============================================================================
// PROCESS: FASTQC
// =============================================================================
// Purpose: Per-sample quality control on FASTQ files
//
// What it does:
//   - Analyzes read quality, GC content, adapter contamination, duplication
//   - Generates an HTML report (human-readable) and ZIP file (parsed by MultiQC)
//   - Runs on raw reads; can also be aliased for trimmed reads (FASTQC_TRIMMED)
//
// Typical resources: ~250MB RAM, 1-2 CPUs, 1-2 minutes per FASTQ
// =============================================================================

process FASTQC {

    tag "FastQC: ${species} / ${sample_id} (${read_type})"
    label 'process_low'                           // ~250MB RAM, 1-2 CPUs per file

    publishDir { "${params.proj_dir()}/${species}/02.FastQC/${read_type}" },    mode: 'copy',    pattern: "*.html"
    publishDir { "${params.proj_dir()}/${species}/02.FastQC/${read_type}" },    mode: 'copy',    pattern: "*.zip"
    publishDir { "${params.proj_dir()}/${species}/07.Logs" },                   mode: 'copy',    pattern: "*.FASTQC.error.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), val(sample_id), path(fastq_files), val(read_type)
    // species      : "Human", "Mouse", or "Xenograft"
    // sample_id    : Sample identifier (e.g., "Sample1")
    // fastq_files  : [R1.fq.gz] for SE or [R1.fq.gz, R2.fq.gz] for PE
    // read_type    : "raw" or "trimmed" — determines output subdirectory

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), path("*_fastqc.zip"), path("*_fastqc.html"),    emit: fastqc_results    // Parsed by MultiQC
    path("*.FASTQC.error.log"),                                         emit: error_log         // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG = "${sample_id}.${read_type}.FASTQC.error.log"

    """
    # Run FastQC on all FASTQ files for this sample
    # --threads: One thread per file (SE=1, PE=2)
    # --quiet  : Suppress progress output for cleaner logs
    # fastq_files.join(' '): Converts Nextflow list → space-separated string
    #   SE: "Sample1_R1.fq.gz"
    #   PE: "Sample1_R1.fq.gz Sample1_R2.fq.gz"
    fastqc \
        --threads "${task.cpus}" \
        --quiet \
        ${fastq_files.join(' ')} \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: FastQC failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: FastQC completed for ${sample_id} (${read_type})" >> "${LOG}"
    """
}

// =============================================================================
// QUICK REFERENCE
// =============================================================================
//
// Output files per FASTQ:
//   *_fastqc.html : Visual QC report (open in browser)
//   *_fastqc.zip  : Raw data parsed by MultiQC
//
// Key metrics to check in the HTML report:
//   Per base quality    : Should be >Q28 across read length
//   GC content          : Should match expected genome GC (~50% human, ~42% mouse)
//   Adapter content     : Flag if >5% — consider trimming
//   Duplication         : High duplication in RNA-seq is normal (highly expressed genes)
//
// Common issues:
//   - Corrupted FASTQ   → Check with: gzip -t file.fq.gz
//   - No output         → Check LOG file; verify FASTQ format
//   - High duplication  → Normal for RNA-seq; not a problem unless extreme (>80%)
// =============================================================================
