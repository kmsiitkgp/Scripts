// =============================================================================
// PROCESS: MD5_VALIDATE
// =============================================================================
// Purpose: Validates MD5 checksums of all downloaded TCGA count files
//
// What it does:
//   - Computes md5sum for every downloaded *_star_gene_counts.tsv file
//   - Compares computed md5 against expected md5 from GDC manifest
//   - Prints summary message to log (pass count / total)
//   - If ALL pass  : prints ✅ message, no files created, workflow continues
//   - If ANY fail  : prints ❌ message, writes failed_manifest.txt (same format
//                    as original manifest), workflow exits with error
//
// Input:
//   - Directory containing all flattened *_star_gene_counts.tsv files
//   - GDC manifest file (contains expected md5 per file)
//
// Output:
//   - failed_manifest.txt (only created if failures exist) → published to proj_dir
//   - MD5_VALIDATE.log                                     → published to reports/
//
// Re-running after failures:
//   - Update manifest_file in TCGA_config.yaml to point to failed_manifest.txt
//   - Re-run pipeline with: bash ~/scripts/nextflow/run_nextflow.sh
//   - gdc-client will re-download only the failed files
//
// Typical resources: 4GB RAM, 1 core, 30-60 minutes for ~11k files (I/O bound)
// =============================================================================

process MD5_VALIDATE {

    tag "MD5 Validate: ${manifest.getName()}"
    label 'process_low'

    publishDir { "${params.proj_dir()}" },         mode: 'copy',    pattern: "failed_manifest.txt"
    publishDir { "${params.proj_dir()}/Logs" },    mode: 'copy',    pattern: "*.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    path(count_files)       // All *_star_gene_counts.tsv files collected into workdir
                            // .collect() in tcga.nf sends all files to one process instance
    path(manifest)          // GDC manifest file (id, filename, md5, size, state)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("failed_manifest.txt"),    optional: true,    emit: failed_manifest    // Only exists if failures found
    path("MD5_VALIDATE.log"),                          emit: log_file           // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG = "MD5_VALIDATE.log"

    """
    echo "=========================================================" >> "${LOG}"
    echo "Start date : \$(date)"                                     >> "${LOG}"
    echo "Manifest   : ${manifest}"                                  >> "${LOG}"
    echo "=========================================================" >> "${LOG}"

    # -----------------------------------------------------------------------
    # STEP 1: Build a lookup table of filename → expected_md5 from manifest
    # Manifest format (tab-separated, skip header line 1):
    #   id    filename    md5    size    state
    # We key on filename (col 2) and store md5 (col 3)
    # -----------------------------------------------------------------------
    N_TOTAL=\$(tail -n +2 "${manifest}" | wc -l)
    echo "Total files in manifest: \${N_TOTAL}" >> "${LOG}"

    # -----------------------------------------------------------------------
    # STEP 2: Validate each downloaded file against expected md5
    # -----------------------------------------------------------------------
    N_PASS=0
    N_FAIL=0

    # Write header for failed manifest to a TEMP file — never write directly
    # to "failed_manifest.txt" here, because the INPUT manifest for this run
    # may itself be named "failed_manifest.txt" (e.g. on a retry run). Writing
    # to that name mid-script would truncate the input file before the while
    # loop below finishes reading it
    echo -e "id\\tfilename\\tmd5\\tsize\\tstate" > failed_manifest.txt.tmp

    while IFS=\$'\t' read -r id filename md5 size state || [ -n "\$id" ]; do

        if [ ! -f "\${filename}" ]; then
            # File is missing entirely — treat as failure
            echo "❌ MISSING : \${filename}" >> "${LOG}"
            echo -e "\${id}\\t\${filename}\\t\${md5}\\t\${size}\\t\${state}" >> failed_manifest.txt.tmp
            N_FAIL=\$((N_FAIL + 1))
            continue
        fi

        # Compute actual md5 of downloaded file
        ACTUAL_MD5=\$(md5sum "\${filename}" | awk '{print \$1}')

        if [ "\${ACTUAL_MD5}" == "\${md5}" ]; then
            N_PASS=\$((N_PASS + 1))
        else
            echo "❌ MISMATCH: \${filename}" >> "${LOG}"
            echo "   Expected : \${md5}"     >> "${LOG}"
            echo "   Actual   : \${ACTUAL_MD5}" >> "${LOG}"
            echo -e "\${id}\\t\${filename}\\t\${md5}\\t\${size}\\t\${state}" >> failed_manifest.txt.tmp
            N_FAIL=\$((N_FAIL + 1))
        fi

    done < <(tail -n +2 "${manifest}")

    # -----------------------------------------------------------------------
    # STEP 3: Report results and exit accordingly
    # -----------------------------------------------------------------------
    echo "=========================================================" >> "${LOG}"
    echo "PASS : \${N_PASS} / \${N_TOTAL}"                           >> "${LOG}"
    echo "FAIL : \${N_FAIL} / \${N_TOTAL}"                           >> "${LOG}"
    echo "=========================================================" >> "${LOG}"

    if [ "\${N_FAIL}" -eq 0 ]; then
        # All files passed — remove empty failed_manifest.txt.tmp (optional output)
        rm -f failed_manifest.txt.tmp
        echo "✅ SUCCESS: All \${N_PASS}/\${N_TOTAL} files passed MD5 validation" >> "${LOG}"
    else
        # Promote temp file to the real failed_manifest.txt now that the
        # input manifest has been fully read and is no longer needed
        mv failed_manifest.txt.tmp failed_manifest.txt
        echo "❌ ERROR: \${N_FAIL} files failed MD5 validation" >> "${LOG}"
        echo "   failed_manifest.txt written to \$(pwd)" >> "${LOG}"
        echo "   Re-run pipeline with manifest_file set to failed_manifest.txt" >> "${LOG}"
        exit 1
    fi
    """
}

// =============================================================================
// QUICK REFERENCE
// =============================================================================
//
// How MD5 validation works:
//   manifest col 2 (filename) → key to find the downloaded file
//   manifest col 3 (md5)      → expected checksum
//   md5sum <file>             → actual checksum
//   match → PASS | mismatch → FAIL
//
// failed_manifest.txt format (same as original manifest):
//   id          filename                          md5        size     state
//   744a6d3d    94027f46..._star_gene_counts.tsv  4a917...   4231952  released
//
// Re-run workflow with failed files:
//   1. Update TCGA_config.yaml:
//        manifest_file: "/home/.../TCGA_Human_PanCancer/failed_manifest.txt"
//   2. Run: cd ~/scripts/nextflow/projects/TCGA && bash ~/scripts/nextflow/run_nextflow.sh
//   3. Repeat until no failed_manifest.txt is produced
//
// Common failure causes:
//   - Network interruption mid-download  → partial file, wrong md5
//   - Disk full during download          → truncated file
//   - GDC server error                   → corrupted transfer
//
// =============================================================================
