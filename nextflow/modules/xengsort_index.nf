process XENGSORT_INDEX {

    tag "Indexing ${href_fasta.name} and ${mref_fasta.name}"
    label 'process_high'                      // XENGSORT indexing requires 30-50GB RAM for human

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    path(href_fasta)        // From RESOLVE_REFERENCES.out.ref_fa
    path(mref_fasta)        // From RESOLVE_REFERENCES.out.ref_fa
    val(genome_version)     // From RESOLVE_REFERENCES.out.genome_version

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("xengsort_index_dir_${genome_version}"),        emit: xengsort_index_dir        // XENGSORT genome index
    //path("XENGSORT_INDEX.error.log"),                  emit: error_log             // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG = "XENGSORT_INDEX.error.log"
    def dir_name = "xengsort_index_dir_${genome_version}"

    """
    xengsort index \
    --index "${genome_version}" \
    --host "${mref_fasta}" \
    --graft "${href_fasta}" \
    --kmer 25 \
    --nobjects 4500000000 \
        2>> "${LOG}" \
        || { echo "❌ ERROR: XENGSORT index generation failed" | tee -a "${LOG}" >&2; exit 1; }

    mkdir -p ${dir_name}
    mv "${genome_version}".info "${genome_version}".hash ${dir_name}/

    echo "✅ SUCCESS: XENGSORT index generation completed" >> "${LOG}"
    """
}
