// =============================================================================
// PROCESS: RNA_RSEQC_QC
// =============================================================================
// Purpose: Comprehensive post-alignment QC on RNA-seq BAM files
//
// What it does (9 analyses):
//   1. Read distribution   : Where reads map (exons, introns, intergenic)
//   2. Inner distance      : Insert size distribution (PE only)
//   3. Junction annotation : Known vs novel splice junctions
//   4. Junction saturation : Is sequencing depth sufficient?
//   5. Gene body coverage  : 5' to 3' bias detection
//   6. Insertion profile   : Insertion artifact patterns
//   7. Deletion profile    : Deletion artifact patterns
//   8. Clipping profile    : Soft/hard clipping patterns
//   9. Mismatch profile    : Systematic error detection
//
// Uses 1M-read subsample (from SAMBAMBA_PREP) for analyses 4 and 5 for speed
// Typical resources: 8-16GB RAM, 15-30 minutes per sample
// =============================================================================

process RNA_RSEQC_QC {

    tag "RSeQC: ${species} / ${sample_id}"
    label 'process_medium'                        // 8-16GB RAM, 4 cores

    publishDir = [
        [path: { "${params.proj_dir()}/${species}/05.RSEQC/01_read_distribution" },   mode: 'copy', pattern: "*.read_distribution*"],
        [path: { "${params.proj_dir()}/${species}/05.RSEQC/02_inner_distance" },      mode: 'copy', pattern: "*.inner_distance*"],
        [path: { "${params.proj_dir()}/${species}/05.RSEQC/03_junction_annotation" }, mode: 'copy', pattern: "*.{splice*pdf,junction_summary.txt,junction_plot.r,junction.bed,junction.xls}"],
        [path: { "${params.proj_dir()}/${species}/05.RSEQC/04_junction_saturation" }, mode: 'copy', pattern: "*junctionSaturation*"],
        [path: { "${params.proj_dir()}/${species}/05.RSEQC/05_deletion_profile" },    mode: 'copy', pattern: "*.deletion_profile*"],
        [path: { "${params.proj_dir()}/${species}/05.RSEQC/06_mismatch_profile" },    mode: 'copy', pattern: "*.mismatch_profile*"],
        [path: { "${params.proj_dir()}/${species}/05.RSEQC/07_insertion_profile" },   mode: 'copy', pattern: "*.insertion_profile*"],
        [path: { "${params.proj_dir()}/${species}/05.RSEQC/08_clipping_profile" },    mode: 'copy', pattern: "*.clipping_profile*"],
        [path: { "${params.proj_dir()}/${species}/05.RSEQC/09_gene_body_coverage" },  mode: 'copy', pattern: "*.geneBodyCoverage*"],
        [path: { "${params.proj_dir()}/${species}/Logs" },                            mode: 'copy', pattern: "*.RSEQC_QC.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species),
        val(sample_id),
        path(bam),              // Full BAM from SAMBAMBA_PREP
        path(bai),              // Full BAM index
        path(bam_1M),           // 1M-read subsample BAM
        path(bai_1M),           // Subsample index
        path(read_len_file),    // Detected read length (used by deletion/mismatch profile)
        path(ref_bed),          // Gene annotation BED12 (from RSEQC_BED)
        path(housekeeping_bed)  // Top 5000 longest transcripts BED (from RSEQC_BED)
    val(mode)                   // "SINGLE_END" or "PAIRED_END" (as List from .collect())

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), val(sample_id), path("${sample_id}*.{txt,log,r,xls}"),
                                                    emit: rseqc_logs,      optional: true    // Data files (parsed by MultiQC)
    path("${sample_id}*.{pdf,jpeg,png,tiff}"),      emit: rseqc_plots,     optional: true    // QC plots
    path("${sample_id}*.bed"),                      emit: rseqc_beds,      optional: true    // Junction BED files
    path("${sample_id}.RSEQC_QC.log"),           emit: error_log                          // Process log
    // optional: true prevents errors when PE-only outputs are absent for SE data

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    // mode arrives as a List from .collect() in rnaseq.nf; access first element
    def SEQUENCING_MODE = (mode[0] == "PAIRED_END") ? "PE" : "SE"
    def LOG = "${sample_id}.RSEQC_QC.log"

    """
    # ==========================================================================
    # QC 1: Read Distribution — where reads map across genomic features
    # Expected: >50% CDS exons, <10% introns, <5% intergenic
    # ==========================================================================
    read_distribution.py \
        --input-file "${bam}" \
        --refgene "${ref_bed}" \
        > "${sample_id}.read_distribution.txt" \
        2>> "${LOG}" \
        || { echo "❌ ERROR: Read distribution failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
    echo "✅ Read distribution completed" >> "${LOG}"


    # ==========================================================================
    # QC 2: Inner Distance — insert size distribution (PE only)
    # Expected: 50-300bp mean for mRNA libraries
    # --mapq 30: High-quality alignments only (>99.9% confidence)
    # ==========================================================================
    if [[ "${SEQUENCING_MODE}" == "PE" ]]; then
        inner_distance.py \
            --input-file "${bam}" \
            --refgene "${ref_bed}" \
            --mapq 30 \
            --out-prefix "${sample_id}" \
            1>> "${LOG}" 2>&1 \
            || { echo "❌ ERROR: Inner distance failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
        echo "✅ Inner distance completed" >> "${LOG}"
    fi


    # ==========================================================================
    # QC 3: Junction Annotation — known vs novel splice junctions
    # Expected: >80% known junctions for well-annotated organisms
    # --min-intron 50: Gaps <50bp treated as deletions, not introns
    # ==========================================================================
    junction_annotation.py \
        --input-file "${bam}" \
        --refgene "${ref_bed}" \
        --mapq 30 \
        --min-intron 50 \
        --out-prefix "${sample_id}" \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: Junction annotation failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
    echo "✅ Junction annotation completed" >> "${LOG}"


    # ==========================================================================
    # QC 4: Junction Saturation — is sequencing depth sufficient?
    # Uses 1M subsample for speed (full BAM can take hours)
    # ==========================================================================
    junction_saturation.py \
        --input-file "${bam_1M}" \
        --refgene "${ref_bed}" \
        --mapq 30 \
        --min-intron 50 \
        --out-prefix "${sample_id}" \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: Junction saturation failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
    echo "✅ Junction saturation completed" >> "${LOG}"


    # ==========================================================================
    # QC 5: Gene Body Coverage — 5' to 3' bias detection
    # Uses 1M subsample + housekeeping genes for a stable baseline
    # Expected: relatively flat coverage profile across gene body
    # ==========================================================================
    geneBody_coverage.py \
        --input "${bam_1M}" \
        --refgene "${housekeeping_bed}" \
        --out-prefix "${sample_id}" \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: Gene body coverage failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
    echo "✅ Gene body coverage completed" >> "${LOG}"


    # ==========================================================================
    # QC 6: Insertion Profile — insertion artifact patterns along reads
    # ==========================================================================
    insertion_profile.py \
        --input-file "${bam}" \
        --sequencing "${SEQUENCING_MODE}" \
        --out-prefix "${sample_id}" \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: Insertion profile failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
    echo "✅ Insertion profile completed" >> "${LOG}"


    # ==========================================================================
    # QC 7: Deletion Profile — deletion patterns along read positions
    # --read-align-length: Read length detected by SAMBAMBA_PREP
    # ==========================================================================
    deletion_profile.py \
        --input "${bam}" \
        --out-prefix "${sample_id}" \
        --read-align-length \$(cat "${read_len_file}") \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: Deletion profile failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
    echo "✅ Deletion profile completed" >> "${LOG}"


    # ==========================================================================
    # QC 8: Clipping Profile — where reads are soft/hard clipped
    # High end-clipping: adapter contamination
    # High internal-clipping: alignment problems
    # ==========================================================================
    clipping_profile.py \
        --input-file "${bam}" \
        --sequencing "${SEQUENCING_MODE}" \
        --out-prefix "${sample_id}" \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: Clipping profile failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
    echo "✅ Clipping profile completed" >> "${LOG}"


    # ==========================================================================
    # QC 9: Mismatch Profile — systematic sequencing error detection
    # ==========================================================================
    mismatch_profile.py \
        --input "${bam}" \
        --out-prefix "${sample_id}" \
        --read-align-length \$(cat "${read_len_file}") \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: Mismatch profile failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
    echo "✅ Mismatch profile completed" >> "${LOG}"


    echo "✅ SUCCESS: All RSeQC analyses completed for ${sample_id}" >> "${LOG}"
    """
}

// =============================================================================
// QUICK REFERENCE
// =============================================================================
//
// Expected metrics (good quality RNA-seq library):
//   Read distribution : >50% CDS exons, <10% introns, <5% intergenic
//   Junction annotation: >80% known junctions
//   Gene body coverage : Relatively flat profile (uniform 5' to 3')
//   Inner distance     : 50-300bp mean
//
// Common issues:
//   High introns (>20%)       : Genomic DNA contamination or RNA degradation
//   High novel junctions (>30%): Wrong annotation version
//   Strong 3' bias            : RNA degradation or poly-A selection issue
//   BED format error          : Chromosome name mismatch between BAM and BED
// =============================================================================
