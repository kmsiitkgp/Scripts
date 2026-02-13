process XENGSORT_CLASSIFY {

    tag "Separating human and mouse reads for ${sample_id}"
    label 'process_high'                      // XENGSORT classification requires 30-50GB RAM for human

    publishDir { "${params.proj_dir()}/Human_split/01.FastQ/xengsort" },    mode: 'copy',    pattern: "${sample_id}_graft_R*.gz"
    publishDir { "${params.proj_dir()}/Mouse_split/01.FastQ/xengsort" },    mode: 'copy',    pattern: "${sample_id}_host_R*.gz"
    publishDir { "${params.proj_dir()}/Mouse_split/07.Logs" },              mode: 'copy',    pattern: "*.XENGSORT_CLASSIFY.error.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(sample_id), path(fastq_files)
    path(xengsort_index_dir)
    val(genome_version)
    // sample_id    : Sample identifier (e.g., "Sample1")
    // fastq_files  : List of FASTQ files [R1.fq.gz] for SE or [R1.fq.gz, R2.fq.gz] for PE

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(sample_id), path("${sample_id}_host_R*.gz"),     emit: host_fastqs
    tuple val(sample_id), path("${sample_id}_graft_R*.gz"),    emit: graft_fastqs
    //path("${sample_id}_{both,neither,ambiguous}_R*.gz"),     emit: junk_fastqs,    optional: true
    path("${sample_id}.XENGSORT_CLASSIFY.error.log"),          emit: error_log             // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG = "${sample_id}.XENGSORT_CLASSIFY.error.log"
    def is_paired = fastq_files.size() > 1
    def r1 = fastq_files[0]
    def r2 = is_paired ? fastq_files[1] : ""

    """
    xengsort classify \
        --index "${xengsort_index_dir}/${genome_version}" \
        --fastq "${r1}" \
        ${is_paired ? "--pairs ${r2}" : ""} \
        --prefix "${sample_id}" \
        --compression gz \
        --threads "${task.cpus}" \
        --mode quick \
    1>> "${LOG}" 2>&1 \
    || { echo "❌ ERROR: Xengsort classification  failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: Xengsort classification completed for ${sample_id}" >> "${LOG}"

    # Rename Surgery: MT4-host.1.fq.gz -> MT4_host_R1.fq.gz
    for f in ${sample_id}-*.1.f*q.gz; do
        [ -e "\$f" ] || continue
        # Use double backslashes for the dots so Nextflow doesn't crash
        newname=\$(echo \$f | sed 's/-/_/' | sed 's/\\.1\\./_R1./')
        mv "\$f" "\$newname"
    done

    if [ "${is_paired}" = "true" ]; then
        for f in ${sample_id}-*.2.f*q.gz; do
            [ -e "\$f" ] || continue
            # Use double backslashes here too
            newname=\$(echo \$f | sed 's/-/_/' | sed 's/\\.2\\./_R2./')
            mv "\$f" "\$newname"
        done
    fi

    """
}