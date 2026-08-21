// =============================================================================
// PROCESS: HMMCOPY_READCOUNTER
// =============================================================================
// Purpose: Count reads in fixed genomic bins from tumor BAM files using
//          HMMcopy readCounter utility
//
// What it does:
//   - Slides across the genome in fixed-size bins (1Mb)
//   - Counts reads falling in each bin per sample
//   - Filters low mapping quality reads (MAPQ < 20)
//   - Outputs a WIG file per sample for input to ICHORCNA_RUN
//
// Typical resources: 2GB RAM, 5-15 minutes per sample
// =============================================================================

process HMMCOPY_READCOUNTER {

    tag "HMMcopy readCounter: ${sample_id}"
    label 'process_low'

    publishDir { "${params.proj_dir()}/02.ReadCounter" },    mode: 'copy',    pattern: "*.wig"
    publishDir { "${params.proj_dir()}/Logs" },              mode: 'copy',    pattern: "*.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(sample_id), path(bam), path(bai)
    // sample_id : Sample identifier (e.g., "Sample1")
    // bam       : Tumor BAM file (deduped, coordinate-sorted)
    // bai       : BAM index file (.bam.bai) — must be staged alongside bam
    val(bin_size)
    // bin_size  : Genomic bin size in base pairs (e.g. 1000000 = 1Mb)
    //             Must match bin_size used in HMMCOPY_REFWIG and ICHORCNA_RUN

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(sample_id), path("${sample_id}.wig"),    emit: wig
    path("${sample_id}.HMMCOPY_READCOUNTER.log"),      emit: log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:
    def LOG = "${sample_id}.HMMCOPY_READCOUNTER.log"

    """
    # Count reads in genomic bins across all autosomes and chrX
    # --window    : Bin size in base pairs
    # --quality   : Minimum mapping quality (MAPQ >= 20 excludes multimappers)
    # --chromosome: Chromosomes to include; must match chr naming in BAM header
    readCounter \
        --window ${bin_size} \
        --quality 20 \
        --chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,\
chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX" \
        ${bam} \
        > ${sample_id}.wig \
        2>> "${LOG}" \
        || { echo "❌ ERROR: readCounter failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: readCounter completed for ${sample_id}" >> "${LOG}"
    """
}

// =============================================================================
// QUICK REFERENCE
// =============================================================================
//
// Output files:
//   *.wig : Read counts per 1Mb bin across autosomes + chrX
//
// MAPQ filter (--quality 20):
//   Removes reads with mapping quality below 20
//   These are likely multimappers that would create false copy number signal
//
// Chromosome naming:
//   Uses chr prefix (chr1, chr2...) matching GRCh38Decoy BAM headers
//   If BAM uses numeric chromosomes (1, 2...) remove the chr prefix
//
// Common issues:
//   - Empty WIG        → BAM index (.bai) missing or chromosome names don't match
//   - All zeros in WIG → Wrong chromosome naming convention (chr vs no chr)
//   - readCounter OOM  → Rare; increase process_low memory if needed
// =============================================================================