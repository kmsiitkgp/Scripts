// =============================================================================
// PROCESS: RSEQC_BED
// =============================================================================
// Purpose: Converts GTF gene annotation to BED12 format required by RSeQC
//
// What it does:
//   - GTF → genePred → BED12 (two-step conversion via UCSC tools)
//   - Creates a housekeeping gene subset (top 5000 longest transcripts)
//     for faster gene body coverage calculation in RSEQC
//
// Why BED12 instead of GTF?
//   - RSeQC tools only accept BED12 format
//   - BED12 is faster to parse (one line per transcript vs many lines in GTF)
//   - Exon structure encoded compactly in blockCount/blockSizes/blockStarts columns
//
// Typical resources: <2GB RAM, ~10-30 seconds
// Uses storeDir: built once per genome version, reused across runs
// =============================================================================

process RSEQC_BED {

    tag "Converting GTF to BED ${assembly}.${release}"
    label 'process_low'                           // Quick format conversion; minimal resources

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), path(ref_fasta), path(ref_gtf), val(genome_version), val(assembly), val(release)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species),
        path("${ref_gtf.baseName}.bed"),                  // Full BED12 annotation
        path("${ref_gtf.baseName}.housekeeping.bed"),     // Top 5000 longest transcripts
        emit: rseqc_bed_tuple
    //path("RSEQC_BED.error.log"),    emit: log    // storeDir: log not captured

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG             = "RSEQC_BED.error.log"
    def ref_bed         = "${ref_gtf.baseName}.bed"
    def housekeeping_bed = "${ref_gtf.baseName}.housekeeping.bed"

    """
    # Step 1: GTF → genePred (UCSC intermediate format preserving exon structure)
    gtfToGenePred "${ref_gtf}" tmp.genePred \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: gtfToGenePred failed" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: GTF → genePred conversion completed" >> "${LOG}"

    # Step 2: genePred → BED12
    # BED12 encodes exon structure in the final 3 columns:
    #   blockCount (exon count), blockSizes, blockStarts (relative to txStart)
    genePredToBed tmp.genePred "${ref_bed}" \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: genePred → BED12 conversion failed" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: genePred → BED12 conversion completed" >> "${LOG}"

    # Step 3: Create housekeeping subset (top 5000 longest transcripts)
    # Used for gene body coverage to get a stable, high-expression baseline
    # awk adds transcript length (col13 = end - start), sort by length desc, take top 5000
    awk 'BEGIN {OFS="\\t"} {print \$0, \$3-\$2}' "${ref_bed}" | \
        sort -t \$'\\t' -k13,13rn | \
        head -n 5000 | \
        cut -f1-12 > "${housekeeping_bed}" \
        2>> "${LOG}" \
        || { echo "❌ ERROR: Housekeeping BED creation failed" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: Housekeeping BED created (top 5000 longest transcripts)" >> "${LOG}"
    """
}

// =============================================================================
// QUICK REFERENCE
// =============================================================================
//
// Conversion chain:
//   GTF (multi-line per gene) → genePred (tabular) → BED12 (one line per transcript)
//
// BED12 format (12 columns):
//   1-3  : chrom, chromStart, chromEnd
//   4-6  : name, score, strand
//   7-8  : thickStart, thickEnd (CDS boundaries)
//   9    : itemRgb (unused)
//   10-12: blockCount, blockSizes, blockStarts (exon structure)
//
// Housekeeping subset:
//   - Top 5000 longest transcripts by genomic span
//   - ~10x faster than full BED for gene body coverage
//   - Sufficient for detecting 5' to 3' bias
//
// Common issues:
//   - Chromosome name mismatch (BED uses "1", BAM uses "chr1") → Check FASTA headers
//   - Empty BED output → Invalid GTF format; check with gtfToGenePred directly
// =============================================================================
