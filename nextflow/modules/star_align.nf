// =========================================================================================
// PROCESS: STAR_ALIGN
// =========================================================================================
// Purpose: Splice-aware alignment of RNA-seq reads to reference genome using STAR
//
// What it does:
//   - Maps reads to genome with splice junction awareness
//   - Discovers novel junctions (two-pass mode)
//   - Generates sorted BAM file for downstream analysis
//   - Produces gene-level counts and junction tables
//
// Typical resources: 30-50GB RAM, ~10-30 minutes per sample on 8 cores
// Output: BAM file, gene counts, splice junctions, alignment stats
//
// For detailed explanation, see: docs/star_align.md
// =========================================================================================

process STAR_ALIGN {

    tag "Aligning fastqs of ${sample_id}"
    label 'process_high'                            // STAR requires 30-50GB RAM for human

    publishDir { "${params.proj_dir()}/${species}_${type}/04.STAR/gene_counts" },        mode: 'copy',    pattern: "*.ReadsPerGene.out.tab"
    publishDir { "${params.proj_dir()}/${species}_${type}/04.STAR/splice_junction" },    mode: 'copy',    pattern: "*.SJ.out.tab"
    publishDir { "${params.proj_dir()}/${species}_${type}/04.STAR/alignment_stats" },    mode: 'copy',    pattern: "*.Log.final.out"
    publishDir { "${params.proj_dir()}/${species}_${type}/07.Logs" },                    mode: 'copy',    pattern: "*.STAR_ALIGN.error.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), val(type), val(sample_id), path(fastq_files), path(star_index_dir)
    // sample_id        : Sample identifier (e.g., "Sample1")
    // fastq_files      : List of FASTQ files [R1.fq.gz] for SE or [R1.fq.gz, R2.fq.gz] for PE
    // star_index_dir   : Pre-built STAR index directory
    val(star_args)                                  // Pre-joined STAR ARGS
    // Never do ${params.STAR_ARGS().join(' ')} inside process. It changes hash on
    // every run and resume fails. So, .join() in main.nf and pass as an argument.

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), val(type),
        val(sample_id),
        path("${sample_id}.Aligned.sortedByCoord.out.bam"),    // Sorted BAM file
        path("${sample_id}.ReadsPerGene.out.tab"),             // Gene counts
        path("${sample_id}.SJ.out.tab"),                       // Splice junctions
        path("${sample_id}.Log.final.out"),                    // Alignment stats
        emit: star_results

    path("${sample_id}.STAR_ALIGN.error.log"),      emit: error_log     // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    // Build read input arguments based on SE vs PE
    // PE: --readFilesIn R1.fq.gz R2.fq.gz
    // SE: --readFilesIn R1.fq.gz
    def MATES_ARGS = fastq_files.size() == 2 ?
        "--readFilesIn ${fastq_files[0]} ${fastq_files[1]}" :
        "--readFilesIn ${fastq_files[0]}"

    def LOG = "${sample_id}.STAR_ALIGN.error.log"

    """
    # Align reads using STAR
    # --genomeDir: Path to STAR index (loaded into RAM)
    # --readFilesIn: Input FASTQ file(s) - SE or PE
    # --outFileNamePrefix: Prefix for all output files
    # --runThreadN: Number of CPU cores to use
    # Additional parameters from config (params.STAR_ARGS):
    #   - Two-pass mode for novel junction discovery
    #   - Gene counting for differential expression
    #   - BAM output sorted by coordinate
    #   - Various filtering and alignment parameters

    STAR \
        --genomeDir "${star_index_dir}" \
        ${MATES_ARGS} \
        ${star_args} \
        --outFileNamePrefix "${sample_id}." \
        --runThreadN "${task.cpus}" \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: STAR alignment failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: STAR alignment completed for ${sample_id}" >> "${LOG}"
    """
}

// =========================================================================================
// QUICK REFERENCE
// =========================================================================================
//
// Key STAR parameters (from config):
//   --twopassMode Basic: Discovers novel junctions, improves mapping 2-5%
//   --quantMode GeneCounts: Generates gene count matrix
//   --outSAMtype BAM SortedByCoordinate: Sorted BAM output
//   --outFilterMultimapNmax 10: Allow up to 10 mapping locations
//
// Output files:
//   ${sample_id}.Aligned.sortedByCoord.out.bam: Main BAM file for IGV/RSeQC
//   ${sample_id}.ReadsPerGene.out.tab: Gene counts (4 columns for strandedness)
//   ${sample_id}.SJ.out.tab: Splice junction table (9 columns)
//   ${sample_id}.Log.final.out: Alignment statistics
//
// Expected mapping rates:
//   Good: >70% uniquely mapped
//   Acceptable: 50-70%
//   Poor: <50% (investigate contamination/wrong reference)
//
// Common issues:
//   - Low mapping → Wrong species/index, contamination, degraded RNA
//   - OOM error → Increase RAM or reduce --limitBAMsortRAM
//   - High multi-mapping → rRNA contamination
//
// For comprehensive parameter guide, see: docs/star_align.md
// ========================================================================================='