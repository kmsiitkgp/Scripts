process SCRNA_CALC_QC_AND_IDENTIFY_SINGLETS {

    tag "Sample: ${sample_id}"
    label 'process_high'
    label 'omics_r'

    publishDir = [
        [path: { "${params.proj_dir()}/${params.species}/Logs" }, mode: 'copy', pattern: "*.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(sample_id), path(sample_rds)
    // sample_id  : Sample identifier (e.g., "SampleA")
    // sample_rds : Seurat object from LOAD_CELLRANGER_AND_CLASSIFY_DROPLETS
    //              contains DropletUtils and/or CellRanger classification columns

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    // IMPORTANT: Output filename MUST include the step name (e.g. sample.CALC_QC_IDENTIFY.rds)
    // Nextflow stages inputs as symlinks with the same basename into each process work dir.
    // If two steps both output "sample.rds", the downstream process gets a symlink named
    // "sample.rds" pointing back to the upstream work dir. When R saves to "./sample.rds"
    // it follows the symlink and overwrites the upstream file — corrupting the cache and
    // causing cache misses on every resume. Unique filenames per step prevent this entirely.
    //
    // NOTE: this Seurat object is NOT filtered — every barcode (Empty Droplet, Low Quality,
    // Doublet, Singlet) is retained. Barcode_QC is a classification column only. Actual
    // subsetting to singlets happens later, in MERGE_AND_PLOT_QC.

    tuple val(sample_id), path("${sample_id}.CALC_QC_IDENTIFY.rds"), emit: sample_rds
    path("${sample_id}.SCRNA_CALC_QC_AND_IDENTIFY_SINGLETS.log"),    emit: log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def script_path    = "${workflow.projectDir}/modules"
    def LOG            = "${sample_id}.SCRNA_CALC_QC_AND_IDENTIFY_SINGLETS.log"
    def gene_cutoff    = params.scrna.qc.gene_cutoff                ?: 250
    def umi_cutoff     = params.scrna.qc.umi_cutoff                 ?: 500
    def mito_cutoff    = params.scrna.qc.mito_cutoff                ?: 0.2
    def novelty_cutoff = params.scrna.qc.novelty_cutoff             ?: 0.8
    def ribo_cutoff    = params.scrna.qc.ribo_cutoff                ?: 0.05
    def multiplet_rate = params.scrna.scdblfinder.multiplet_rate    ?: 0.008

    """
    Rscript "${script_path}/SCRNA/SCRNA_calc_qc_and_identify_singlets.R" \
        "${sample_id}" \
        "${sample_rds}" \
        "${gene_cutoff}" \
        "${umi_cutoff}" \
        "${mito_cutoff}" \
        "${novelty_cutoff}" \
        "${ribo_cutoff}" \
        "${multiplet_rate}" \
        "." \
        > "${LOG}" 2>&1

    EXIT_CODE=\$?
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "❌ ERROR: SCRNA_CALC_QC_AND_IDENTIFY_SINGLETS failed for ${sample_id}" | tee -a "${LOG}" >&2
        exit \$EXIT_CODE
    fi

    echo "✅ SUCCESS: SCRNA_CALC_QC_AND_IDENTIFY_SINGLETS completed for ${sample_id}" >> "${LOG}"
    """
}
