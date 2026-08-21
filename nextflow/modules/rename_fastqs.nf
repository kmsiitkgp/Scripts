// =============================================================================
// PROCESS: RENAME_FASTQS
// =============================================================================
// Purpose: Batch renames FASTQ files using a two-column mapping file
//
// What it does:
//   - Reads old_name → new_name pairs from a mapping file
//   - Renames files in the Nextflow work directory
//   - Skips files already renamed (idempotent — safe to re-run)
//   - Reports missing files without stopping the batch
//   - Fails the process if any file could not be found under either name
//
// Typical use case: SRA downloads (SRR*) → project sample names
//
// Mapping file format (whitespace-separated, no header):
//   SRR10162346_1.fastq.gz    Sample1_R1.fq.gz
//   SRR10162346_2.fastq.gz    Sample1_R2.fq.gz
// =============================================================================

process RENAME_FASTQS {

    tag "Renaming FASTQ files"
    label 'process_low'                           // Simple file operations; minimal resources

    //publishDir { "${params.read_dir}/renamed" },    mode: 'copy',    pattern: "*.{fastq,fq}.gz"
    //publishDir { "${params.proj_dir()}" },          mode: 'copy',    pattern: "RENAME_FASTQS.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    path(fastqs)        // All FASTQ files staged into work directory (via .collect())
    path(map_file)      // Mapping file: old_name<whitespace>new_name

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("*.{fastq,fq}.gz"),            emit: renamed_fastqs    // Renamed FASTQ files
    path("RENAME_FASTQS.log"),    emit: error_log         // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG = "RENAME_FASTQS.log"

    """
    # Allow loop to continue on individual failures so we can report all missing files
    set +e

    echo "=== FASTQ Rename Process ===" > "${LOG}"
    echo "Working Directory: \$(pwd)" >> "${LOG}"
    echo "Files in Directory: \$(ls *.fastq.gz 2>/dev/null | wc -l)" >> "${LOG}"
    echo "Start time: \$(date)" >> "${LOG}"
    echo "Mapping file: ${map_file}" >> "${LOG}"
    echo "" >> "${LOG}"

    RENAMED_COUNT=0
    FAILURE_COUNT=0

    # Ensure mapping file ends with newline (prevents last line being skipped by read)
    [[ \$(tail -c1 "${map_file}") != "" ]] && echo >> "${map_file}"

    # Process each mapping entry
    # IFS=\$' \t' allows both space and tab as column separators
    while IFS=\$' \t' read -r OLD_NAME NEW_NAME; do

        # Skip comment lines and empty lines
        [[ "\$OLD_NAME" =~ ^# ]] && continue
        [[ -z "\$OLD_NAME" ]] && continue
        [[ -z "\$NEW_NAME" ]] && continue

        if [[ -f "\$OLD_NAME" ]]; then
            # Original file found → rename it
            mv "\$OLD_NAME" "\$NEW_NAME"
            echo "✓ Renamed: \$OLD_NAME → \$NEW_NAME" >> "${LOG}"
            ((RENAMED_COUNT++))

        elif [[ -f "\$NEW_NAME" ]]; then
            # File already has new name → idempotent skip
            echo "ℹ️  Already renamed: \$NEW_NAME" >> "${LOG}"
            ((RENAMED_COUNT++))

        else
            # Neither name found → log error, continue batch
            echo "❌ ERROR: File not found: \$OLD_NAME" >> "${LOG}"
            ((FAILURE_COUNT++))
        fi

    done < "${map_file}"

    echo "" >> "${LOG}"
    echo "======================================" >> "${LOG}"
    echo "Files renamed/verified : \$RENAMED_COUNT" >> "${LOG}"
    echo "Files not found        : \$FAILURE_COUNT" >> "${LOG}"
    echo "End time: \$(date)" >> "${LOG}"
    echo "======================================" >> "${LOG}"

    # Fail the process if any files were missing
    # Prevents downstream processes from running with incomplete data
    if [[ \$FAILURE_COUNT -gt 0 ]]; then
        echo "❌ ERROR: \$FAILURE_COUNT file(s) not found — check mapping file" | tee -a "${LOG}" >&2
        exit 1
    fi

    echo "✅ SUCCESS: All files renamed successfully" >> "${LOG}"
    set -e
    """
}
