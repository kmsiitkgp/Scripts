// =============================================================================
// PROCESS: RNA_SALMON_QUANT
// =============================================================================
// Purpose: Fast, alignment-free transcript quantification using Salmon
//
// What it does:
//   - Quantifies transcript abundances using selective alignment (k-mer based)
//   - Auto-detects library strand orientation
//   - Applies GC, sequence, and positional bias corrections
//   - Generates TPM values and estimated read counts per transcript
//   - Copies quant.sf to a sample-named file for easy downstream collection
//
// Typical resources: 8-12GB RAM, 2-5 minutes per sample on 4 cores
// =============================================================================

process RNA_SALMON_QUANT {

    tag "Salmon quant: ${species} / ${sample_id}"
    label 'process_high'                        // ~30GB RAM, 4 cores

    publishDir = [
        [path: { "${params.proj_dir()}/${species}/03.Salmon" },             mode: 'copy', pattern: "${sample_id}"],
        [path: { "${params.proj_dir()}/${species}/03.Salmon/quant_files" }, mode: 'copy', pattern: "*.quant.sf"],
        [path: { "${params.proj_dir()}/${species}/Logs" },                  mode: 'copy', pattern: "*.SALMON_QUANT.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), val(sample_id), path(fastq_files), path(salmon_index_dir)
    // species          : "Human", "Mouse", or "Xenograft"
    // sample_id        : Sample identifier (e.g., "Sample1")
    // fastq_files      : [R1.fq.gz] for SE or [R1.fq.gz, R2.fq.gz] for PE
    // salmon_index_dir : Pre-built Salmon index directory from SALMON_INDEX
    val(salmon_args)
    // Pre-joined Salmon argument string from rnaseq.nf
    // Never call params.SALMON_ARGS().join(' ') inside a process — it changes the
    // task hash on every run and breaks -resume caching

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), val(sample_id), path(sample_id),       emit: salmon_quant_dir     // Full Salmon output directory
    tuple val(species), path("${sample_id}.quant.sf"),         emit: salmon_quant_file    // Transcript abundances (TPM + counts)
    path("${sample_id}.SALMON_QUANT.log"),               emit: error_log            // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    // Build read input flag(s) based on SE vs PE
    def MATES_ARGS = fastq_files.size() == 2 \
        ? "--mates1 ${fastq_files[0]} --mates2 ${fastq_files[1]}" \
        : "--unmatedReads ${fastq_files[0]}"

    def LOG = "${sample_id}.SALMON_QUANT.log"

    """
    # Quantify transcripts using Salmon selective alignment
    # --validateMappings : Validate k-mer chains (more accurate than quasi-mapping)
    # --index            : Path to Salmon index
    # --output           : Output directory (named after sample for easy collection)
    # --threads          : CPU cores (use 1 thread for reproducibility)
    # salmon_args (from config) typically includes:
    #   --libType A   : Auto-detect strand orientation
    #   --gcBias      : Correct for GC content bias (~5% accuracy improvement)
    #   --seqBias     : Correct for random hexamer priming bias
    #   --posBias     : Correct for positional fragment distribution bias
    salmon quant \
        --validateMappings \
        --index "${salmon_index_dir}" \
        ${salmon_args} \
        ${MATES_ARGS} \
        --threads ${task.cpus} \
        --output "${sample_id}" \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: Salmon quantification failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    # Copy quant.sf to a sample-named file for easy collection across samples
    # Allows gathering all quant.sf files without complex glob patterns
    cp "${sample_id}/quant.sf" "${sample_id}.quant.sf" \
        || { echo "❌ ERROR: Failed to copy quant.sf for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: Salmon quantification completed for ${sample_id}" >> "${LOG}"
    """
}

// =============================================================================
// QUICK REFERENCE
// =============================================================================
//
// Output directory structure (${sample_id}/):
//   quant.sf             : Main results — Name, Length, EffectiveLength, TPM, NumReads
//   meta_info.json       : Run parameters and mapping rate summary
//   lib_format_counts.json : Library type detection results
//   aux_info/            : Auxiliary data for tximeta/tximport
//
// Library type detection:
//   SE: SF (forward), SR (reverse), U (unstranded)
//   PE: ISF (forward), ISR (reverse), IU (unstranded)
//
// Expected mapping rates:
//   >70% : Good
//   50-70%: Investigate
//   <50%  : Poor — check species, contamination, RNA degradation
//
// Common issues:
//   - Low mapping     → Wrong index version, contamination, degraded RNA
//   - Mixed lib types → Inconsistent library prep across samples
//   - OOM error       → Increase RAM to 16GB+
// =============================================================================
