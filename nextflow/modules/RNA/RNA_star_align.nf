// =============================================================================
// PROCESS: RNA_STAR_ALIGN
// =============================================================================
// Purpose: Splice-aware alignment of RNA-seq reads to reference genome using STAR
//
// What it does:
//   - Maps reads to the genome with full splice junction awareness
//   - Discovers novel junctions (via two-pass mode configured in star_args)
//   - Produces coordinate-sorted BAM, gene counts, splice junction table, and stats
//
// Typical resources: 30-50GB RAM, 10-30 minutes per sample on 8 cores
// =============================================================================

process RNA_STAR_ALIGN {

    tag "STAR align: ${species} / ${sample_id}"
    label 'process_high'                          // 30-50GB RAM required for human/mouse

    publishDir = [
        [path: { "${params.proj_dir()}/${species}/04.STAR/gene_counts" },     mode: 'copy', pattern: "*.ReadsPerGene.out.tab"],
        [path: { "${params.proj_dir()}/${species}/04.STAR/splice_junction" }, mode: 'copy', pattern: "*.SJ.out.tab"],
        [path: { "${params.proj_dir()}/${species}/04.STAR/alignment_stats" }, mode: 'copy', pattern: "*.Log.final.out"],
        [path: { "${params.proj_dir()}/${species}/Logs" },                    mode: 'copy', pattern: "*.STAR_ALIGN.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), val(sample_id), path(fastq_files), path(star_index_dir)
    // species        : "Human", "Mouse", or "Xenograft"
    // sample_id      : Sample identifier (e.g., "Sample1")
    // fastq_files    : [R1.fq.gz] for SE or [R1.fq.gz, R2.fq.gz] for PE
    // star_index_dir : Pre-built STAR index directory from STAR_INDEX
    val(star_args)
    // Pre-joined STAR argument string from rnaseq.nf
    // Never call params.STAR_ARGS().join(' ') inside a process — it changes the
    // task hash on every run and breaks -resume caching

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species),
        val(sample_id),
        path("${sample_id}.Aligned.sortedByCoord.out.bam"),    // Coordinate-sorted BAM
        path("${sample_id}.ReadsPerGene.out.tab"),             // Gene-level read counts
        path("${sample_id}.SJ.out.tab"),                       // Splice junction table
        path("${sample_id}.Log.final.out"),                    // Alignment summary stats
        emit: star_results
    path("${sample_id}.STAR_ALIGN.log"),    emit: error_log    // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    // Build read input flag(s) based on SE vs PE
    def MATES_ARGS = fastq_files.size() == 2 \
        ? "--readFilesIn ${fastq_files[0]} ${fastq_files[1]}" \
        : "--readFilesIn ${fastq_files[0]}"

    def LOG = "${sample_id}.STAR_ALIGN.log"

    """
    # Align reads to genome using STAR
    # --genomeDir        : Path to STAR index (loaded into shared memory)
    # --outFileNamePrefix: Prefix for all output files
    # --runThreadN       : CPU cores
    # star_args (from config) typically includes:
    #   --twopassMode Basic          : Two-pass mode for improved novel junction detection
    #   --quantMode GeneCounts       : Generate ReadsPerGene.out.tab
    #   --outSAMtype BAM SortedByCoordinate
    #   --readFilesCommand zcat      : Decompress gzipped input on the fly
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

// =============================================================================
// QUICK REFERENCE
// =============================================================================
//
// Output files:
//   *.Aligned.sortedByCoord.out.bam : Sorted BAM for IGV / RSeQC / sambamba
//   *.ReadsPerGene.out.tab          : Gene counts (cols: geneID, unstranded, sense, antisense)
//   *.SJ.out.tab                    : Splice junction table (novel + annotated)
//   *.Log.final.out                 : Alignment rate summary
//
// Expected mapping rates:
//   >70% uniquely mapped : Good
//   50-70%               : Investigate (contamination? wrong reference?)
//   <50%                 : Poor — check species, index version, RNA quality
//
// Common issues:
//   - Low mapping    → Wrong species/index, contamination, degraded RNA
//   - OOM error      → Increase RAM or reduce --limitBAMsortRAM in star_args
//   - High multi-map → rRNA contamination (add --outFilterMultimapNmax 1 to filter)
// =============================================================================
