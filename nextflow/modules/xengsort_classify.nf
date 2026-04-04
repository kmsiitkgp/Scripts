// =============================================================================
// PROCESS: XENGSORT_CLASSIFY
// =============================================================================
// Purpose: Separates PDX reads into human (graft) and mouse (host) FASTQ files
//
// What it does:
//   - Classifies each read as graft, host, both, neither, or ambiguous
//   - Writes graft (human) and host (mouse) reads to separate FASTQ files
//   - Renames XenGSort output files to pipeline naming convention:
//       Sample-graft.1.fq.gz → Sample_graft_R1.fq.gz
//       Sample-host.2.fq.gz  → Sample_host_R2.fq.gz
//
// Downstream in rnaseq.nf:
//   graft_fastqs → labeled "Human"     → STAR/Salmon with Human index
//   host_fastqs  → labeled "Mouse"     → STAR/Salmon with Mouse index
//   original FASTQs (unseparated)      → labeled "Xenograft" → Human reference
//
// Typical resources: ~30GB RAM, 4-8 cores
// =============================================================================

process XENGSORT_CLASSIFY {

    tag "Classifying graft/host reads for ${sample_id}"
    label 'process_high'                          // K-mer lookup across two genomes; ~30GB RAM

    // Uncomment if you want the fastq files to be saved
    //publishDir { "${params.proj_dir()}/Human/01.FastQ/xengsort" },     mode: 'copy',    pattern: "${sample_id}_graft_R*.gz"
    //publishDir { "${params.proj_dir()}/Mouse/01.FastQ/xengsort" },     mode: 'copy',    pattern: "${sample_id}_host_R*.gz"
    publishDir { "${params.proj_dir()}/Xenograft/07.Logs" },           mode: 'copy',    pattern: "*.XENGSORT_CLASSIFY.error.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(sample_id), path(fastq_files)
    // sample_id   : Sample identifier (e.g., "Sample1")
    // fastq_files : [R1.fq.gz] for SE or [R1.fq.gz, R2.fq.gz] for PE
    path(xengsort_index_dir)                      // Index directory from XENGSORT_INDEX
    val(genome_version)                           // Composite version (e.g., "GRCh38.115_GRCm39.115")

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(sample_id), path("${sample_id}_graft_R*.gz"),    emit: graft_fastqs    // Human reads
    tuple val(sample_id), path("${sample_id}_host_R*.gz"),     emit: host_fastqs     // Mouse reads
    //path("${sample_id}_{both,neither,ambiguous}_R*.gz"),     emit: junk_fastqs,    optional: true
    path("${sample_id}.XENGSORT_CLASSIFY.error.log"),          emit: error_log       // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG       = "${sample_id}.XENGSORT_CLASSIFY.error.log"
    def is_paired = fastq_files.size() > 1
    def r1        = fastq_files[0]
    def r2        = is_paired ? fastq_files[1] : ""

    """
    # Classify reads as graft (human) or host (mouse)
    # --mode quick: Fast classification (recommended for production)
    # --compression gz: Compress output FASTQ files
    xengsort classify \
        --index "${xengsort_index_dir}/${genome_version}" \
        --fastq "${r1}" \
        ${is_paired ? "--pairs ${r2}" : ""} \
        --prefix "${sample_id}" \
        --compression gz \
        --threads "${task.cpus}" \
        --mode quick \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: XenGSort classification failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: XenGSort classification completed for ${sample_id}" >> "${LOG}"

    # Rename R1 files to pipeline convention
    # XenGSort default: Sample-graft.1.fq.gz → pipeline: Sample_graft_R1.fq.gz
    for f in ${sample_id}-*.1.f*q.gz; do
        [ -e "\$f" ] || continue
        newname=\$(echo \$f | sed 's/-/_/' | sed 's/\\.1\\./_R1./')
        mv "\$f" "\$newname"
    done

    if [ "${is_paired}" = "true" ]; then
        # Rename R2 files to pipeline convention
        # XenGSort default: Sample-graft.2.fq.gz → pipeline: Sample_graft_R2.fq.gz
        for f in ${sample_id}-*.2.f*q.gz; do
            [ -e "\$f" ] || continue
            newname=\$(echo \$f | sed 's/-/_/' | sed 's/\\.2\\./_R2./')
            mv "\$f" "\$newname"
        done
    fi
    """
}

// =============================================================================
// QUICK REFERENCE
// =============================================================================
//
// Classification output categories:
//   graft     : Confidently human reads → Human RNA-seq analysis
//   host      : Confidently mouse reads → Mouse RNA-seq analysis
//   both      : Maps equally to both genomes → discarded
//   neither   : Maps to neither genome → discarded
//   ambiguous : Low confidence assignment → discarded
//
// Common issues:
//   - Low graft rate (<30%)  → Verify sample is truly PDX; check index genome versions
//   - Index path error       → Ensure xengsort_index_dir contains <genome_version>.hash
//   - Rename failure         → Check XenGSort version outputs expected filename format
// =============================================================================
