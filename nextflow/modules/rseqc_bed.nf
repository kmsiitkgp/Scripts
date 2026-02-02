// =========================================================================================
// PROCESS: RSEQC_BED
// =========================================================================================
// Purpose: Converts GTF gene annotations to BED12 format for RSeQC tools
//
// What it does:
//   - Converts GTF → genePred → BED12 (two-step conversion)
//   - Creates housekeeping gene subset for fast gene body coverage
//   - Preserves exon structure and CDS/UTR boundaries
//
// Why needed:
//   - RSeQC requires BED12 format (not GTF)
//   - BED is faster to parse than GTF
//   - One line per transcript vs multiple lines in GTF
//
// Typical resources: <2GB RAM, ~10-30 seconds
//
// For detailed explanation, see: docs/gtf_to_bed.md
// =========================================================================================

process RSEQC_BED {

    tag "Converting ${ref_gtf.name} to BED"
    label 'process_low'                          // Minimal resources (quick conversion)

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    path(ref_gtf)    // Gene annotation GTF file

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("${ref_gtf.baseName}.bed"),                        emit: ref_bed              // Full BED12 file
    path("${ref_gtf.baseName}.housekeeping.bed"),           emit: housekeeping_bed     // Subset for fast QC
    path("RSEQC_BED.error.log"),                            emit: error_log            // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def ref_bed = "${ref_gtf.baseName}.bed"
    def housekeeping_bed = "${ref_gtf.baseName}.housekeeping.bed"
    def LOG = "RSEQC_BED.error.log"

    """
    # Step 1: GTF → genePred (intermediate format)
    # Captures gene models with exon structure
    gtfToGenePred "${ref_gtf}" tmp.genePred \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: gtfToGenePred failed" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: gtfToGenePred completed" >> "${LOG}"

    # Step 2: genePred → BED12
    # BED12 columns: chr, start, end, name, score, strand,
    #                thickStart, thickEnd, itemRgb,
    #                blockCount, blockSizes, blockStarts
    genePredToBed tmp.genePred "${ref_bed}" \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: GTF to BED conversion failed" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: GTF to BED conversion completed" >> "${LOG}"

    # Step 3: Create housekeeping subset (top 5000 longest transcripts)
    # Used for faster gene body coverage calculation
    # Calculates transcript length, sorts by length, takes top 5000
    awk 'BEGIN {OFS="\\t"} {print \$0, \$3-\$2}' "${ref_bed}" | \
    sort -t \$'\\t' -k13,13rn | \
    head -n 5000 | \
    cut -f1-12 > "${housekeeping_bed}" \
        2>> "${LOG}" \
        || { echo "❌ ERROR: Housekeeping BED creation failed" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: Housekeeping BED creation completed" >> "${LOG}"
    """
}

// =========================================================================================
// QUICK REFERENCE
// =========================================================================================
//
// Conversion steps:
//   GTF (multi-line per gene) → genePred (tabular) → BED12 (one line per transcript)
//
// BED12 format (12 columns):
//   1-3: chr, start, end (genomic coordinates)
//   4-6: name, score, strand
//   7-8: thickStart, thickEnd (CDS boundaries)
//   9: itemRgb (color, not used)
//   10-12: blockCount, blockSizes, blockStarts (exon structure)
//
// Housekeeping subset:
//   - Top 5000 longest transcripts
//   - Used for fast gene body coverage (~10x speedup)
//   - Sufficient for detecting 3' bias
//
// Output files:
//   ${gtf_basename}.bed: Complete annotation (~100-500MB)
//   ${gtf_basename}.housekeeping.bed: Subset (~5-10MB)
//
// Common issues:
//   - Chromosome mismatch (BED "chr1" vs BAM "1")
//   - Invalid GTF format
//   - Wrong strand information
//
// For comprehensive guide, see: docs/gtf_to_bed.md
// ========================================================================================='