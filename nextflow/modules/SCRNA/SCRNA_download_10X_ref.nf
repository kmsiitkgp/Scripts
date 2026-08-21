// =============================================================================
// PROCESS: SCRNA_DOWNLOAD_10X_REF
// =============================================================================
// Purpose: Downloads genome FASTA and GTF annotation from Ensembl FTP
//
// What it does:
//   - Downloads primary assembly FASTA and GTF for Human or Mouse
//   - Retries up to 5 times per file to handle transient DNS/network failures
//   - Decompresses .gz files after successful download
//
// Why retry logic?
//   - Ensembl FTP can be slow or briefly unavailable
//   - A single transient failure would otherwise abort the entire pipeline run
//
// Uses storeDir: files are downloaded once and reused across runs
// =============================================================================

process SCRNA_DOWNLOAD_10X_REF {

    tag "Downloading ${species} reference (Ensembl release ${release})"
    label 'process_medium'                        // Network-bound; moderate resources

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), val(fa_name), val(gtf_name), val(genome_version), val(assembly), val(release)
    // species        : "Human" or "Mouse"
    // fa_name        : Full FASTA ref_tarball (used in output block for storeDir checking)
    // gtf_name       : Full GTF ref_tarball (used in output block for storeDir checking)
    // genome_version : Version string (e.g., "GRCh38.115")
    // assembly       : Assembly name (e.g., "GRCh38")
    // release        : Ensembl release number (e.g., "115")

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), path("refdata-gex-*"),                              emit: index_dir
    tuple val(species), val(genome_version), val(assembly), val(release),   emit: ref_tuple

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG = "DOWNLOAD_10X_REF.log"

    """
    set -euo pipefail

    # Map pipeline species name to 10X URL path component
    if [[ "${species}" == "Human" ]]; then
        url="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
    elif [[ "${species}" == "Mouse" ]]; then
        url="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz"
    elif [[ "${species}" == "Xenograft" ]]; then
        url="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38_and_GRCm39-2024-A.tar.gz"
    else
        echo "❌ ERROR: Unsupported species: ${species}" | tee -a "${LOG}" >&2
        exit 1
    fi

    # Extract the tarball name from the URL
    ref_tarball=\$(basename "\${url}")

    # --- Download Reference ---
    # Retry loop handles transient DNS/connection failures
    success=false
    for i in {1..5}; do
        if wget --retry-connrefused --waitretry=10 --read-timeout=30 --tries=10 \
            --output-document "\${ref_tarball}" "\${url}"; then
            success=true
            break
        else
            echo "Reference download attempt \$i failed. Retrying in 30s..." >> "${LOG}"
            sleep 30
        fi
    done
    \$success || { echo "❌ ERROR: Reference download failed after 5 attempts" | tee -a "${LOG}" >&2; exit 1; }

    tar -xzf "\${ref_tarball}"  2>> "${LOG}"
    rm "\${ref_tarball}" # Clean up the huge tarball

    echo "✅ SUCCESS: Reference download completed for ${species} (release ${release})" >> "${LOG}"
    """
}
