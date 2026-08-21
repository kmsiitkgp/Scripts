// =============================================================================
// PROCESS: RNA_DOWNLOAD_ENSEMBL_REF
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

process RNA_DOWNLOAD_ENSEMBL_REF {

    tag "Downloading ${species} reference (Ensembl release ${release})"
    label 'process_medium'                        // Network-bound; moderate resources

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), val(fa_name), val(gtf_name), val(genome_version), val(assembly), val(release)
    // species        : "Human" or "Mouse"
    // fa_name        : Full FASTA filename (used in output block for storeDir checking)
    // gtf_name       : Full GTF filename (used in output block for storeDir checking)
    // genome_version : Version string (e.g., "GRCh38.115")
    // assembly       : Assembly name (e.g., "GRCh38")
    // release        : Ensembl release number (e.g., "115")

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), path(fa_name), path(gtf_name), val(genome_version), val(assembly), val(release),
        emit: ref_tuple

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG = "DOWNLOAD_ENSEMBL_REF.log"

    """
    set -euo pipefail

    # Map pipeline species name to Ensembl URL path component
    if [[ "${species}" == "Human" ]]; then
        ensembl_species="homo_sapiens"
    elif [[ "${species}" == "Mouse" ]]; then
        ensembl_species="mus_musculus"
    else
        echo "❌ ERROR: Unsupported species: ${species}" | tee -a "${LOG}" >&2
        exit 1
    fi

    # --- Download FASTA ---
    # Retry loop handles transient DNS/connection failures
    success=false
    for i in {1..5}; do
        if wget --retry-connrefused --waitretry=10 --read-timeout=30 --tries=10 \
            --output-document "${fa_name}.gz" \
            "https://ftp.ensembl.org/pub/release-${release}/fasta/\${ensembl_species}/dna/${fa_name}.gz"; then
            success=true
            break
        else
            echo "FASTA download attempt \$i failed. Retrying in 30s..." >> "${LOG}"
            sleep 30
        fi
    done
    \$success || { echo "❌ ERROR: FASTA download failed after 5 attempts" | tee -a "${LOG}" >&2; exit 1; }

    gunzip -f "${fa_name}.gz" 2>> "${LOG}"
    echo "✅ SUCCESS: FASTA downloaded and decompressed" >> "${LOG}"

    # --- Download GTF ---
    # Retry loop handles transient DNS/connection failures
    success=false
    for i in {1..5}; do
        if wget --retry-connrefused --waitretry=10 --read-timeout=30 --tries=10 \
            --output-document "${gtf_name}.gz" \
            "https://ftp.ensembl.org/pub/release-${release}/gtf/\${ensembl_species}/${gtf_name}.gz"; then
            success=true
            break
        else
            echo "GTF download attempt \$i failed. Retrying in 30s..." >> "${LOG}"
            sleep 30
        fi
    done
    \$success || { echo "❌ ERROR: GTF download failed after 5 attempts" | tee -a "${LOG}" >&2; exit 1; }

    gunzip -f "${gtf_name}.gz" 2>> "${LOG}"
    echo "✅ SUCCESS: GTF downloaded and decompressed" >> "${LOG}"

    echo "✅ SUCCESS: Reference download completed for ${species} (release ${release})" >> "${LOG}"
    """
}
