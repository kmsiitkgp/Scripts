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

process RENAME_FASTQS {

    tag "Renaming FASTQ files"
    label 'process_low'                                      // Simple file operations, minimal resources

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    path(fastqs)                                             // All FASTQ files to rename (staged in work dir)
    path(map_file)                                           // Mapping file: old_name new_name

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("*.{fastq,fq}.gz"),                emit: renamed_fastqs       // Renamed FASTQ files
    path("RENAME_FASTQS.error.log"),        emit: error_log            // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG = "RENAME_FASTQS.error.log"

    """
    # Temporarily disable "exit on error" to allow batch processing
    # This lets us rename as many files as possible and report all failures
    set +e

    # =============================================================================
    # INITIALIZATION
    # =============================================================================

    echo "=== FASTQ Renaming Process ===" > "${LOG}"
    echo "Start time: \$(date)" >> "${LOG}"
    echo "Mapping file: ${map_file}" >> "${LOG}"
    echo "" >> "${LOG}"

    # Initialize counters
    RENAMED_COUNT=0
    FAILURE_COUNT=0

    # Ensure mapping file ends with newline (prevents last line from being skipped)
    [[ \$(tail -c1 "${map_file}") != "" ]] && echo >> "${map_file}"

    # =============================================================================
    # BATCH RENAME LOOP
    # =============================================================================

    echo "Processing renames..." >> "${LOG}"
    echo "" >> "${LOG}"

    # Read mapping file line by line
    # IFS=\$' \t' allows both space and tab as separators
    # read -r prevents backslash interpretation
    while IFS=\$' \t' read -r OLD_NAME NEW_NAME; do

        # Skip comment lines (starting with #)
        if [[ "\$OLD_NAME" =~ ^# ]]; then
            continue
        fi

        # Skip empty lines
        if [[ -z "\$OLD_NAME" ]]; then
            continue
        fi

        # --- Rename Logic ---
        if [[ -f "\$OLD_NAME" ]]; then
            # Case 1: Original file found → Rename it
            mv "\$OLD_NAME" "\$NEW_NAME"
            echo "✓ Renamed: \$OLD_NAME → \$NEW_NAME" >> "${LOG}"
            ((RENAMED_COUNT++))

        elif [[ -f "\$NEW_NAME" ]]; then
            # Case 2: File already has new name → Skip (idempotent)
            echo "ℹ️  Already renamed: \$NEW_NAME" >> "${LOG}"
            ((RENAMED_COUNT++))

        else
            # Case 3: Neither old nor new name exists → Report error
            echo "❌ ERROR: File not found: \$OLD_NAME" >> "${LOG}"
            ((FAILURE_COUNT++))
        fi

    done < "${map_file}"

    # =============================================================================
    # SUMMARY AND VALIDATION
    # =============================================================================

    echo "" >> "${LOG}"
    echo "========================================================" >> "${LOG}"
    echo "RENAMING SUMMARY" >> "${LOG}"
    echo "========================================================" >> "${LOG}"
    echo "Total files renamed/verified: \$RENAMED_COUNT" >> "${LOG}"
    echo "Total failures: \$FAILURE_COUNT" >> "${LOG}"
    echo "End time: \$(date)" >> "${LOG}"
    echo "========================================================" >> "${LOG}"

    # Critical validation: Fail the process if any files are missing
    # This prevents downstream processes from running with incomplete data
    if [[ \$FAILURE_COUNT -gt 0 ]]; then
        echo "" >> "${LOG}"
        echo "⚠️  WARNING: \$FAILURE_COUNT file(s) could not be renamed!" >> "${LOG}"
        echo "Check the mapping file and ensure all files exist." >> "${LOG}"
        echo "" >> "${LOG}"

        # Output to stderr for Nextflow to capture
        echo "❌ ERROR: RENAME_FASTQS failed - \$FAILURE_COUNT file(s) not found" >&2
        echo "See ${LOG} for details" >&2

        exit 1
    fi

    echo "✅ SUCCESS: All files renamed successfully" >> "${LOG}"

    # Re-enable exit on error for subsequent commands
    set -e
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
