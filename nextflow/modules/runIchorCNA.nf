// =============================================================================
// PROCESS: RUN_ICHORCNA
// =============================================================================
// Purpose: Estimate ctDNA tumor fraction and copy number profile from
//          low-coverage tumor WIG file using ichorCNA R package
//
// What it does:
//   - Corrects read counts for GC content and mappability biases
//   - Fits a Hidden Markov Model to segment the genome into copy number states
//   - Estimates tumor fraction as 1 - normal fraction
//   - Outputs tumor fraction estimate, copy number segments, and profile plots
//
// Typical resources: 4GB RAM, 5-10 minutes per sample
// =============================================================================

process RUN_ICHORCNA {

    tag "ichorCNA: ${sample_id}"
    label 'process_low'

    publishDir { "${params.proj_dir()}/03.ichorCNA/${sample_id}" },    mode: 'copy',    pattern: "*.{txt,seg}"
    publishDir { "${params.proj_dir()}/03.ichorCNA/${sample_id}" },    mode: 'copy',    pattern: "${sample_id}"
    publishDir { "${params.proj_dir()}/Logs" },                        mode: 'copy',    pattern: "*.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(sample_id), path(reads_wig), path(gc_wig), path(map_wig)
    // sample_id : Sample identifier (e.g., "Sample1")
    // reads_wig : Per-sample read count WIG from HMMCOPY_READCOUNTER
    // gc_wig    : GC content WIG from HMMCOPY_REFWIG (shared across all samples)
    // map_wig   : Mappability WIG from HMMCOPY_REFWIG (shared across all samples)
    val(ichorcna_params)
    // ichorcna_params : Map of ichorCNA settings from project_info.yaml Section 4g
    //                   contains bin_size, ploidy, normal, max_cn, include_homd

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(sample_id),
        path("${sample_id}.params.txt"),    // Tumor fraction + ploidy estimates
        path("${sample_id}.seg.txt"),        // Copy number segments
        emit: results
    path("${sample_id}"),                   emit: plots_dir      // Copy number profile plots
    path("${sample_id}.RUN_ICHORCNA.log"),  emit: log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:
    def LOG          = "${sample_id}.RUN_ICHORCNA.log"
    def script_path  = "${workflow.projectDir}/modules"
    def ploidy       = ichorcna_params.ploidy       ?: "c(2,3)"
    def normal       = ichorcna_params.normal        ?: "c(0.5,0.6,0.7,0.8,0.9,0.95,0.99)"
    def max_cn       = ichorcna_params.max_cn        ?: 5
    def include_homd = ichorcna_params.include_homd  ?: false

    """
    Rscript "${script_path}/runIchorCNA.R" \
        --id           "${sample_id}" \
        --WIG          "${reads_wig}" \
        --gcWig        "${gc_wig}" \
        --mapWig       "${map_wig}" \
        --ploidy       "${ploidy}" \
        --normal       "${normal}" \
        --maxCN        ${max_cn} \
        --includeHOMD  ${include_homd} \
        --chrs         "c(1:22, \\"X\\")" \
        --chrTrain     "c(1:22)" \
        --chrNormalize "c(1:22)" \
        --genomeStyle  "UCSC" \
        --genomeBuild  "hg38" \
        --estimateNormal        True \
        --estimatePloidy        True \
        --estimateScPrevalence  False \
        --outDir       "." \
        > "${LOG}" 2>&1
    EXIT_CODE=\$?
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "❌ ERROR: runIchorCNA.R failed for ${sample_id}" | tee -a "${LOG}" >&2
        exit \$EXIT_CODE
    fi
    echo "✅ SUCCESS: ichorCNA completed for ${sample_id}" >> "${LOG}"
    """
}

// =============================================================================
// QUICK REFERENCE
// =============================================================================
//
// Output files (per sample in 03.ichorCNA/SampleID/):
//   *.params.txt : Tumor fraction, ploidy, and model parameters
//                  KEY OUTPUT — look for "Tumor Fraction" column
//   *.seg.txt    : Copy number segments across the genome
//   *_*.pdf      : Copy number profile plots (one per chromosome)
//
// Interpreting tumor fraction:
//   Value in params.txt = 1 - normal fraction
//   e.g. normal fraction 0.85 → tumor fraction 0.15 → 15% ctDNA
//   Reliable minimum detection threshold: ~3% tumor fraction
//
// Key parameters:
//   --genomeStyle UCSC    : required for GRCh38 BAMs with chr prefix
//   --genomeBuild hg38    : tells ichorCNA which centromere coords to use
//   --estimateScPrevalence False : no subclone modeling (standard for ctDNA)
//   --chrs / --chrTrain   : autosomes + chrX for analysis; autosomes only for training
//
// Common issues:
//   - Tumor fraction = 0   → no detectable CNAs; sample may be very low ctDNA
//                            or tumor is near-diploid with few copy number changes
//   - Noisy profile plot   → very low coverage (<0.1x); result may be unreliable
//   - All segments neutral → same as tumor fraction = 0 above
// =============================================================================
