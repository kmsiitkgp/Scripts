process PREPARE_XENOGRAFT_REF {

    tag "Creating Xenograft Reference"
    label 'process_medium'

    input:
    path(human_fa)
    path(human_gtf)
    path(mouse_fa)
    path(mouse_gtf)
    val(fasta_name)  // final merged fasta
    val(gtf_name)    // final merged gtf

    output:
    path("${fasta_name}"), emit: ref_fa
    path("${gtf_name}"),   emit: ref_gtf

    script:

    def LOG = "PREPARE_XENOGRAFT.error.log"

    """
    set -euo pipefail
    echo "Starting Xenograft hybridization..." >> "${LOG}"

    # --- Prefix FASTA headers ---
    awk '{if(\$0 ~ /^>/) print ">h_" substr(\$0,2); else print \$0}' ${human_fa} > h_p.fa 2>> "${LOG}"
    awk '{if(\$0 ~ /^>/) print ">m_" substr(\$0,2); else print \$0}' ${mouse_fa} > m_p.fa 2>> "${LOG}"

    # Merge into final Xenograft FASTA
    cat h_p.fa m_p.fa > "${fasta_name}" 2>> "${LOG}"

    # --- Prefix GTF chromosomes ---
    awk '/^[^#]/ { sub(/^/, "h_"); print }' ${human_gtf} > h_p.gtf 2>> "${LOG}"
    awk '/^[^#]/ { sub(/^/, "m_"); print }' ${mouse_gtf} > m_p.gtf 2>> "${LOG}"

    # Merge into final Xenograft GTF
    cat h_p.gtf m_p.gtf > "${gtf_name}" 2>> "${LOG}"

    # Clean up temporary files
    rm h_p.fa m_p.fa h_p.gtf m_p.gtf

    echo "Xenograft hybridization completed successfully." >> "${LOG}"
    """
}