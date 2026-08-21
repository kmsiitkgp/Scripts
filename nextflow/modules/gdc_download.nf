// =============================================================================
// PROCESS: GDC_DOWNLOAD
// =============================================================================
// Purpose: Downloads TCGA RNA-seq count files from GDC portal using gdc-client
//
// What it does:
//   - Runs gdc-client download using the provided manifest file
//   - Downloads all files into subdirectories named by File ID (GDC default)
//   - Flattens downloaded *_star_gene_counts.tsv files into a single output directory
//
// Input:
//   - GDC manifest file (tab-separated: id, filename, md5, size, state)
//
// Output:
//   - All *_star_gene_counts.tsv files published to 01.Downloads/
//   - GDC_DOWNLOAD.log published to reports/
//
// Notes:
//   - gdc-client creates one subdirectory per File ID during download
//   - mv step flattens: gdc_tmp/{file_id}/*.tsv → 01.Downloads/*.tsv
//   - --no-related-files : skip annotation/index sidecar files (counts only)
//   - --no-annotations   : skip annotation files
//   - --retry-amount     : retry each file N times on network failure
//   - -n                 : number of parallel download connections
//
// Typical resources: 64GB RAM, 8 cores, 2-6 hours for ~11k files
// =============================================================================

process GDC_DOWNLOAD {

    tag "GDC Download: ${manifest.getName()}"
    label 'process_high'

    publishDir { "${params.proj_dir()}/01.Downloads" },    mode: 'copy',    pattern: "*.tsv"
    publishDir { "${params.proj_dir()}/reports" },         mode: 'copy',    pattern: "*.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    path(manifest)          // GDC manifest file (gdc_manifest.*.txt)
    val(gdc_args)           // Pre-joined gdc-client argument string from tcga.nf
                            // Never call params closures inside a process —
                            // it changes the task hash on every run and breaks -resume

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("*.tsv"),              emit: count_files    // Flattened *_star_gene_counts.tsv files
    path("GDC_DOWNLOAD.log"),   emit: log_file       // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG       = "GDC_DOWNLOAD.log"
    def DL_TMPDIR = "gdc_tmp"

    """
    echo "=========================================================" >> "${LOG}"
    echo "Start date : \$(date)"                                     >> "${LOG}"
    echo "Manifest   : ${manifest}"                                  >> "${LOG}"
    echo "GDC args   : ${gdc_args}"                                  >> "${LOG}"
    echo "=========================================================" >> "${LOG}"

    # -----------------------------------------------------------------------
    # STEP 1: Download using gdc-client
    # --dir             : temporary destination directory
    # --manifest        : manifest file listing files to download
    # --no-related-files: skip annotation/index sidecar files
    # --no-annotations  : skip annotation files
    # --retry-amount    : retry each file N times on failure
    # -n                : number of parallel download connections
    # -----------------------------------------------------------------------
    mkdir -p "${DL_TMPDIR}"

    set +e
    gdc-client download \
        --dir "${DL_TMPDIR}" \
        --manifest "${manifest}" \
        --no-related-files \
        --no-annotations \
        ${gdc_args} \
        1>> "${LOG}" 2>&1
    GDC_EXIT=\$?
    set -e

    echo "gdc-client exit code: \${GDC_EXIT}" >> "${LOG}"

    # -----------------------------------------------------------------------
    # STEP 2: Count downloaded files before flattening
    # -----------------------------------------------------------------------
    N_MANIFEST=\$(tail -n +2 "${manifest}" | wc -l)
    N_DOWNLOADED=\$(find "${DL_TMPDIR}" -name "*_star_gene_counts.tsv" | wc -l)

    echo "Files in manifest : \${N_MANIFEST}" >> "${LOG}"
    echo "Files downloaded  : \${N_DOWNLOADED}" >> "${LOG}"

    if [ "\${N_DOWNLOADED}" -eq 0 ]; then
        echo "❌ ERROR: No files downloaded. Check manifest and network." | tee -a "${LOG}" >&2
        exit 1
    fi

    # -----------------------------------------------------------------------
    # STEP 3: Flatten — move all .tsv files from subdirectories into workdir
    # gdc-client creates: gdc_tmp/{file_id}/*_star_gene_counts.tsv
    # publishDir will move them to: 01.Downloads/
    # -----------------------------------------------------------------------
    find "${DL_TMPDIR}" -name "*_star_gene_counts.tsv" -exec mv -t . {} + \
        2>> "${LOG}" || { echo "ERROR: Failed to flatten downloaded files" | tee -a "${LOG}" >&2; exit 1; }

    # Clean up temporary directory
    rm -rf "${DL_TMPDIR}"

    echo "✅ SUCCESS: Downloaded and flattened \${N_DOWNLOADED}/\${N_MANIFEST} files" >> "${LOG}"
    """
}

// =============================================================================
// QUICK REFERENCE
// =============================================================================
//
// GDC manifest format (tab-separated, header on line 1):
//   id          filename                          md5        size     state
//   744a6d3d    94027f46..._star_gene_counts.tsv  4a917...   4231952  released
//
// gdc-client directory structure (before flattening):
//   gdc_tmp/
//   ├── 744a6d3d-b666-49aa-8d26-47f34e3d1eb5/
//   │   └── 94027f46-..._star_gene_counts.tsv
//   └── 4ecc1f1a-8ff4-4552-a5e8-7a9652b6d1d5/
//       └── df45fb41-..._star_gene_counts.tsv
//
// After flattening (publishDir moves to 01.Downloads/):
//   01.Downloads/
//   ├── 94027f46-..._star_gene_counts.tsv
//   ├── df45fb41-..._star_gene_counts.tsv
//   └── ...
//
// Common issues:
//   "Connection refused"      → GDC portal may be down; check https://portal.gdc.cancer.gov/
//   "No space left on device" → Ensure base_dir partition has enough space (~50GB)
//   File count mismatch       → MD5_VALIDATE will catch incomplete downloads
//   gdc-client skips files    → Already downloaded files are skipped on re-run (-resume safe)
//
// =============================================================================
