// =============================================================================
// PROCESS: GENERATE_MD5
// =============================================================================
// Purpose: Generates MD5 checksum for each FASTQ file
//
// What it does:
//   - Runs md5sum on a single FASTQ file
//   - Writes checksum to a .md5 sidecar file
//   - Per-file .md5 outputs are merged into a single manifest in rnaseq.nf
//     via .collectFile() → manifest_md5.txt
//
// Why MD5?
//   - Verifies file integrity after transfer or download
//   - Manifest provides an audit trail of all input files
// =============================================================================

process GENERATE_MD5 {

    tag "MD5: ${fastq}"
    label 'process_low'                           // md5sum is I/O bound; minimal CPU/RAM

    //publishDir { "${params.proj_dir()}" },         mode: 'copy',    pattern: "*.md5"
    publishDir { "${params.proj_dir()}/reports" }, mode: 'copy',    pattern: "*.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    path(fastq)    // Single FASTQ file (process runs once per file in parallel)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("${fastq}.md5"),               emit: md5_file     // Per-file MD5 checksum
    path("GENERATE_MD5.log"),     emit: error_log    // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG = "GENERATE_MD5.log"

    """
    md5sum "${fastq}" > "${fastq}.md5" \
        2>> "${LOG}" \
        || { echo "❌ ERROR: md5sum failed for ${fastq}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: MD5 generated for ${fastq}" >> "${LOG}"
    """
}
