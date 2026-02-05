// =========================================================================================
// PROCESS: RENAME_FASTQS
// =========================================================================================
// Purpose: Batch rename FASTQ files using a mapping file
//
// What it does:
//   - Reads a two-column mapping file (old_name → new_name)
//   - Renames all FASTQ files according to the mapping
//   - Handles already-renamed files gracefully
//   - Reports missing files without stopping the entire batch
//   - Provides detailed renaming statistics
//
// Why this matters:
//   - SRA downloads (SRR*) need conversion to sample names
//   - Standardizes naming for downstream processes
//   - Prevents manual errors when renaming many files
//   - Validates all files exist before pipeline continues
//
// Mapping file format (rename_map.txt):
//   - Two columns: old_name new_name (whitespace-separated)
//   - No header row
//   - Unix line endings (LF, not CRLF)
//   - Example:
//       SRR10162346_1.fastq.gz    Sample1_S1_L001_R1_001.fastq.gz
//       SRR10162346_2.fastq.gz    Sample1_S1_L001_R2_001.fastq.gz
//
// For detailed explanation, see: docs/rename_fastqs.md
// =========================================================================================

process GENERATE_MD5 {

    tag "Generating MD5 for ${fastq}"
    label 'process_low'                                      // Simple file operations, minimal resources

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    path(fastq)                                             // FASTQ file

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("${fastq}.md5"),                  emit: md5_file        // Renamed FASTQ files
    path("GENERATE_MD5.error.log"),        emit: error_log       // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG = "GENERATE_MD5.error.log"

    """
    md5sum ${fastq} > ${fastq}.md5 \
        2>> "${LOG}" \
        || { echo "❌ ERROR: md5sum failed" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: md5sum completed" >> "${LOG}"
    """
}

// =========================================================================================
// QUICK REFERENCE
// =========================================================================================
//
// Mapping file requirements:
//   - Format: old_name<whitespace>new_name
//   - No column headers
//   - Unix line endings (LF), not Windows (CRLF)
//   - Can use spaces or tabs as separator
//   - Comments allowed (lines starting with #)
//   - Empty lines ignored
//
// Example mapping file (rename_map.txt):
//   # RNA-seq project samples
//   SRR10162346_1.fastq.gz    Sample1_S1_L001_R1_001.fastq.gz
//   SRR10162346_2.fastq.gz    Sample1_S1_L001_R2_001.fastq.gz
//   SRR10162347_1.fastq.gz    Sample2_S2_L001_R1_001.fastq.gz
//   SRR10162347_2.fastq.gz    Sample2_S2_L001_R2_001.fastq.gz
//
// Common use cases:
//   1. SRA downloads → Project names:
//      SRR12345_1.fq.gz → Patient001_Tumor_R1.fq.gz
//
//   2. Sequencing facility → Standard format:
//      ABC123_ATGC_L001_R1_001.fastq.gz → Sample1_S1_L001_R1_001.fastq.gz
//
//   3. Add tissue/condition labels:
//      Sample1_R1.fq.gz → Sample1_Tumor_R1.fq.gz
//
// Process behavior:
//   - Idempotent: Safe to re-run (skips already-renamed files)
//   - Fail-fast: Exits with error if ANY file is missing
//   - Detailed logging: Reports every rename operation
//   - Validates all files before downstream processing
//
// Troubleshooting:
//   - "File not found" errors → Check mapping file paths match staged files
//   - Line endings issue → Convert CRLF to LF: dos2unix rename_map.txt
//   - Wrong separator → Use spaces or tabs, not both mixed
//   - Missing files → Verify all files listed in map exist in input
//
// Usage in main workflow:
//   RENAME_FASTQS(
//       Channel.fromPath("${params.raw_fastq_dir}/*.fq.gz").collect(),
//       file("${params.proj_dir}/rename_map.txt")
//   )
//
// For comprehensive guide, see: docs/rename_fastqs.md
// =========================================================================================
