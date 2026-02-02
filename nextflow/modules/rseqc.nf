// =========================================================================================
// PROCESS: RSEQC
// =========================================================================================
// Purpose: Comprehensive quality control analysis on aligned RNA-seq BAM files
//
// What it does:
//   - Read distribution across genomic features (exons, introns, intergenic)
//   - Insert size distribution (PE only)
//   - Splice junction annotation (known vs novel)
//   - Junction saturation analysis
//   - Gene body coverage uniformity (5' to 3' bias)
//   - Sequencing artifact profiling
//
// Why this matters:
//   - Detects library prep issues (3' bias, RNA degradation)
//   - Identifies contamination (genomic DNA, non-coding RNA)
//   - Verifies sequencing depth and strand-specificity
//
// Typical resources: 8-16GB RAM, 15-30 minutes per sample (with subsampling)
//
// For detailed explanation, see: docs/rseqc.md
// =========================================================================================

process RSEQC {

    tag "RSeQC on ${sample_id}"
    label 'process_medium'                          // 4 cores, 12GB RAM typical

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(sample_id),
        path(bam),
        path(bai),
        path(bam_1M),
        path(bai_1M),
        path(read_len_file)                          // Sample BAM files and metadata
    path(ref_bed)                                    // Gene annotation in BED12 format
    path(housekeeping_bed)                           // Housekeeping genes for gene body coverage
    val(mode)                                        // "SINGLE_END" or "PAIRED_END" (as list from .collect())

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("${sample_id}*.{pdf,jpeg,png,tiff}"),      emit: rseqc_plots,     optional: true      // Visualization plots
    path("${sample_id}*.{txt,log,r,xls}"),          emit: rseqc_logs,      optional: true      // Data files and logs
    path("${sample_id}*.bed"),                      emit: rseqc_beds,      optional: true      // Junction BED files
    path("${sample_id}.RSEQC.error.log"),           emit: error_log                            // Process log
    // Note: optional: true prevents errors when PE-only files don't exist for SE data

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    // Translate Nextflow mode to RSeQC parameter format
    // mode is a list from .collect(), so we access first element
    def SEQUENCING_MODE = (mode[0] == "PAIRED_END") ? "PE" : "SE"

    def LOG = "${sample_id}.RSEQC.error.log"

    """
    # =============================================================================
    # QC STEP 1: Read Distribution Across Genomic Features
    # =============================================================================
    # Determines where reads map (exons, introns, intergenic)
    # Expected: 50-70% CDS, <10% introns, <5% intergenic

    read_distribution.py \
        --input-file "${bam}" \
        --refgene "${ref_bed}" \
        > "${sample_id}.read_distribution.txt" \
        2>> "${LOG}" \
        || { echo "❌ ERROR: Read distribution failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: Read distribution completed for ${sample_id}" >> "${LOG}"


    # =============================================================================
    # QC STEP 2: Inner Distance Between Read Pairs (PE ONLY)
    # =============================================================================
    # Measures insert size minus read lengths
    # Expected: 50-300bp for mRNA libraries
    # --mapq 30: Only use high-quality alignments (99.9% confidence)

    if [[ "${SEQUENCING_MODE}" == "PE" ]]; then
        inner_distance.py \
            --input-file "${bam}" \
            --refgene "${ref_bed}" \
            --mapq 30 \
            --out-prefix "${sample_id}" \
            1>> "${LOG}" 2>&1 \
            || { echo "❌ ERROR: Inner distance failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

        echo "✅ SUCCESS: Inner distance completed for ${sample_id}" >> "${LOG}"
    fi


    # =============================================================================
    # QC STEP 3: Splice Junction Annotation
    # =============================================================================
    # Classifies junctions as known vs novel
    # Expected: >80% known junctions for well-studied organisms
    # --min-intron 50: Gaps <50bp are deletions, not introns

    junction_annotation.py \
        --input-file "${bam}" \
        --refgene "${ref_bed}" \
        --mapq 30 \
        --min-intron 50 \
        --out-prefix "${sample_id}" \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: Junction annotation failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: Junction annotation completed for ${sample_id}" >> "${LOG}"


    # =============================================================================
    # QC STEP 4: Junction Saturation Analysis
    # =============================================================================
    # Tests if sequencing depth is sufficient
    # Uses 1M read subsample for speed (full BAM takes hours)

    junction_saturation.py \
        --input-file "${bam_1M}" \
        --refgene "${ref_bed}" \
        --mapq 30 \
        --min-intron 50 \
        --out-prefix "${sample_id}" \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: Junction saturation failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: Junction saturation completed for ${sample_id}" >> "${LOG}"


    # =============================================================================
    # QC STEP 5: Gene Body Coverage (5' to 3' Bias Detection)
    # =============================================================================
    # Analyzes coverage uniformity across genes
    # Expected: Relatively flat profile (uniform coverage)
    # Uses housekeeping genes for more uniform baseline

    geneBody_coverage.py \
        --input "${bam_1M}" \
        --refgene "${housekeeping_bed}" \
        --out-prefix "${sample_id}" \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: Gene body coverage failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: Gene body coverage completed for ${sample_id}" >> "${LOG}"


    # =============================================================================
    # QC STEP 6: Insertion Profile
    # =============================================================================
    # Analyzes insertion artifact patterns
    # --sequencing: PE or SE (different analysis strategies)

    insertion_profile.py \
        --input-file "${bam}" \
        --sequencing "${SEQUENCING_MODE}" \
        --out-prefix "${sample_id}" \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: Insertion profile failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: Insertion profile completed for ${sample_id}" >> "${LOG}"


    # =============================================================================
    # QC STEP 7: Deletion Profile
    # =============================================================================
    # Analyzes deletion patterns along reads
    # --read-align-length: From read length file calculated upstream

    deletion_profile.py \
        --input "${bam}" \
        --out-prefix "${sample_id}" \
        --read-align-length \$(cat "${read_len_file}") \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: Deletion profile failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: Deletion profile completed for ${sample_id}" >> "${LOG}"


    # =============================================================================
    # QC STEP 8: Clipping Profile
    # =============================================================================
    # Analyzes where reads are soft/hard clipped
    # High end-clipping: adapter issues
    # High internal-clipping: alignment problems

    clipping_profile.py \
        --input-file "${bam}" \
        --sequencing "${SEQUENCING_MODE}" \
        --out-prefix "${sample_id}" \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: Clipping profile failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: Clipping profile completed for ${sample_id}" >> "${LOG}"


    # =============================================================================
    # QC STEP 9: Mismatch Profile
    # =============================================================================
    # Analyzes mismatch patterns to detect systematic errors

    mismatch_profile.py \
        --input "${bam}" \
        --out-prefix "${sample_id}" \
        --read-align-length \$(cat "${read_len_file}") \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: Mismatch profile failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: Mismatch profile completed for ${sample_id}" >> "${LOG}"


    # =============================================================================
    # All QC steps completed successfully
    # =============================================================================
    echo "✅ SUCCESS: All RSeQC analyses completed for ${sample_id}" >> "${LOG}"
    """
}

// =========================================================================================
// QUICK REFERENCE
// =========================================================================================
//
// RSeQC analyses performed:
//   1. Read distribution: Where reads map (exons/introns/intergenic)
//   2. Inner distance: Insert size distribution (PE only)
//   3. Junction annotation: Known vs novel splice junctions
//   4. Junction saturation: Sequencing depth adequacy
//   5. Gene body coverage: 5' to 3' bias detection
//   6-9. Artifact profiling: Insertions, deletions, clipping, mismatches
//
// Key parameters:
//   --mapq 30: Only high-quality alignments (99.9% confidence)
//   --min-intron 50: Gaps <50bp are deletions, not introns
//
// Expected metrics (good RNA-seq library):
//   Read distribution: 50-70% CDS, <10% introns, <5% intergenic
//   Junction annotation: >80% known junctions
//   Gene body coverage: Relatively flat (uniform) profile
//   Inner distance: 50-300bp mean (mRNA libraries)
//
// Output files:
//   ${sample_id}.read_distribution.txt: Feature distribution stats
//   ${sample_id}.junction.txt: Junction classification
//   ${sample_id}.geneBodyCoverage.txt: Coverage profile data
//   ${sample_id}*.pdf: Individual sample plots
//   ${sample_id}*.r: R scripts to regenerate plots
//
// Common issues:
//   High introns (>20%): Genomic DNA contamination or degraded RNA
//   High intergenic (>10%): Wrong annotation or contamination
//   High novel junctions (>30%): Wrong annotation version
//   Strong 3' bias: RNA degradation or poly-A selection issues
//   BED format error: Chromosome name mismatch between BAM and BED
//
// For comprehensive guide and troubleshooting, see: docs/rseqc.md
// =========================================================================================
