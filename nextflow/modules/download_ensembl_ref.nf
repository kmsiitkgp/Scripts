process DOWNLOAD_ENSEMBL_REF {

    tag "Downloading ${species} Ensembl Reference"
    label 'process_medium'

    input:
    val(species)          // Human or Mouse
    val(fa_name)          // Full FASTA filename from metadata in order use in output: for storeDir checking
    val(gtf_name)         // Full GTF filename from metadata in order use in output: for storeDir checking:
    val(release)          // Release from metadata

    output:
    tuple val(species), path(fa_name), path(gtf_name),   emit: ref_tuple

    script:

    def LOG = "DOWNLOAD_ENSEMBL_REF.error.log"

    """
    set -euo pipefail
    echo "Starting reference download for: ${species}" >> "${LOG}"

    # Determine Ensembl species code
    if [[ "${species}" == "Human" ]]; then
        ensembl_species="homo_sapiens"
    elif [[ "${species}" == "Mouse" ]]; then
        ensembl_species="mus_musculus"
    else
        echo "Unsupported species: ${species}" >&2
        exit 1
    fi

    # --- Download FASTA ---
    # We wrap this in a loop to handle temporary DNS failures
    success=false
    for i in {1..5}; do
        if wget --retry-connrefused --waitretry=10 --read-timeout=30 --tries=10 \
            --output-document "${fa_name}.gz" \
            "https://ftp.ensembl.org/pub/release-${release}/fasta/\${ensembl_species}/dna/${fa_name}.gz"; then
            success=true
            break
        else
            echo "Download attempt \$i failed. Retrying in 20s..." >> "${LOG}"
            sleep 30
        fi
    done

    if [ "\$success" = false ]; then
        echo "Failed to download Fasta after 5 attempts." >> "${LOG}"
        exit 1
    fi

    gunzip -f "${fa_name}.gz" 2>> "${LOG}"

    # --- Download GTF ---
    # We wrap this in a loop to handle temporary DNS failures
    success=false
    for i in {1..5}; do
        if wget --retry-connrefused --waitretry=10 --read-timeout=30 --tries=10 \
            --output-document "${gtf_name}.gz" \
            "https://ftp.ensembl.org/pub/release-${release}/gtf/\${ensembl_species}/${gtf_name}.gz"; then
            success=true
            break
        else
            echo "Download attempt \$i failed. Retrying in 20s..." >> "${LOG}"
            sleep 30
        fi
    done

    if [ "\$success" = false ]; then
        echo "Failed to download GTF after 5 attempts." >> "${LOG}"
        exit 1
    fi

    gunzip -f "${gtf_name}.gz" 2>> "${LOG}"

    echo "✅ Download completed for ${species}" >> "${LOG}"
    """
}
