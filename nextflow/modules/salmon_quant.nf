// =========================================================================================
// PROCESS: SALMON_QUANT
// =========================================================================================
// Purpose: Fast, alignment-free transcript quantification using Salmon
//
// What it does:
//   - Quantifies transcript abundances using k-mer pseudo-alignment
//   - Estimates both transcript-level and gene-level expression
//   - Applies bias corrections (GC, sequence, positional)
//   - Auto-detects library type (strand orientation)
//   - Generates TPM and read count estimates
//
// Salmon advantages:
//   - 10-30x faster than traditional alignment + counting
//   - Better multi-mapping read handling (EM algorithm)
//   - No alignment BAM files needed
//   - Built-in bias correction
//
// Typical resources: 8-12GB RAM, 2-5 minutes per sample on 4 cores
//
// For detailed explanation, see: docs/salmon_quant.md
// =========================================================================================

process SALMON_QUANT {

    tag "Quantifying ${sample_id}"
    label 'process_medium'                          // 4 cores, 12GB RAM typical

    publishDir { "${params.proj_dir()}/${species}_${type}/03.Salmon" },                mode: 'copy',    pattern: "${sample_id}"
    publishDir { "${params.proj_dir()}/${species}_${type}/03.Salmon/quant_files" },    mode: 'copy',    pattern: "*.quant.sf"
    publishDir { "${params.proj_dir()}/${species}_${type}/07.Logs" },                  mode: 'copy',    pattern: "*.SALMON_QUANT.error.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), val(type), val(sample_id), path(fastq_files), path(salmon_index_dir)
    // sample_id        : Sample identifier (e.g., "Sample1")
    // fastq_files      : List of FASTQ files [R1.fq.gz] for SE or [R1.fq.gz, R2.fq.gz] for PE
    // salmon_index_dir : Pre-built Salmon index directory
    val(salmon_args)                                     // Pre-joined SALMON ARGS from config
    // Never do ${params.SALMON_ARGS().join(' ')} inside process. It changes hash on
    // every run and resume fails. So, .join() in main.nf and pass as an argument.

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), val(type), val(sample_id), path(sample_id),   emit: salmon_quant_dir      // Full output directory
    tuple val(species), path("${sample_id}.quant.sf"),     emit: salmon_quant_file     // Transcript abundances (TPM, counts)
    path("${sample_id}.SALMON_QUANT.error.log"),           emit: error_log             // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG = "${sample_id}.SALMON_QUANT.error.log"

    // Build read input arguments based on SE vs PE
    // PE: --mates1 R1.fq.gz --mates2 R2.fq.gz
    // SE: --unmatedReads R1.fq.gz
    def MATES_ARGS = fastq_files.size() == 2 ?
        "--mates1 ${fastq_files[0]} --mates2 ${fastq_files[1]}" :
        "--unmatedReads ${fastq_files[0]}"

    """
    # Quantify transcripts using Salmon
    # --validateMappings: Selective alignment mode (validates k-mer chains, more accurate)
    # --index: Path to Salmon index
    # --output: Directory for results (named after sample)
    # --threads: Number of CPU cores
    # Additional parameters from config (salmon_args):
    #   --libType A: Auto-detect library type (strand orientation)
    #   --gcBias: Correct for GC content bias
    #   --seqBias: Correct for hexamer priming bias
    #   --posBias: Correct for positional coverage bias

    salmon quant \
        --validateMappings \
        --index "${salmon_index_dir}" \
        ${salmon_args} \
        ${MATES_ARGS} \
        --threads "${task.cpus}" \
        --output "${sample_id}" \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: Salmon quantification failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    # Copy quant.sf to sample-named file for easy collection
    # This enables gathering all quant.sf files for tximport without complex glob patterns
    cp "${sample_id}/quant.sf" ${sample_id}.quant.sf \
        || { echo "❌ ERROR: Failed to copy quant.sf for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: Salmon quantification completed for ${sample_id}" >> "${LOG}"
    """
}

// =========================================================================================
// QUICK REFERENCE
// =========================================================================================
//
// Output directory structure (${sample_id}/):
//   quant.sf: Main results (5 columns: Name, Length, EffectiveLength, TPM, NumReads)
//   meta_info.json: Run parameters and version info
//   lib_format_counts.json: Library type detection results
//   aux_info/: Auxiliary data for downstream tools
//   cmd_info.json: Command line that was run
//
// Library type detection results:
//   SE: SF (forward), SR (reverse), U (unstranded)
//   PE: ISF (forward), ISR (reverse), IU (unstranded)
//
// Bias corrections (from config):
//   --gcBias: Corrects PCR amplification GC bias (~5% improvement)
//   --seqBias: Corrects random hexamer priming bias (~3% improvement)
//   --posBias: Corrects fragment distribution bias (~2% improvement)
//
// Expected mapping rates:
//   Good: >70% mapped
//   Acceptable: 50-70%
//   Poor: <50% (check species, contamination, degradation)
//
// Downstream analysis (R/tximport):
//   ```r
//   library(tximport)
//   files <- list.files(pattern = "quant.sf", recursive = TRUE, full.names = TRUE)
//   names(files) <- basename(dirname(files))
//   txi <- tximport(files, type = "salmon", tx2gene = tx2gene_df)
//   dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~condition)
//   ```
//
// Common issues:
//   - Low mapping → Wrong index, contamination, degraded RNA
//   - High duplication → Low input, over-amplification
//   - Mixed library types → Check protocol consistency
//   - OOM error → Increase RAM to 16GB+
//
// For comprehensive guide and troubleshooting, see: docs/salmon_quant.md
// =========================================================================================
