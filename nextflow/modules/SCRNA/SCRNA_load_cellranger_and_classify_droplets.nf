process SCRNA_LOAD_CELLRANGER_AND_CLASSIFY_DROPLETS {

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
    tuple val(sample_id), val(sample_dir)
    // sample_id  : Sample identifier (e.g., "SampleA")
    // sample_dir : absolute path to CellRanger output for this sample.
    //              Passed as val() not path() — if passed as path(), Nextflow would
    //              checksum the entire CellRanger directory tree to validate the cache,
    //              which is expensive and fragile (any timestamp change = cache miss).
    //              As val(), Nextflow compares only the string — fast and stable.
    //              The R script uses it as an absolute path string directly, and reads
    //              both raw_feature_bc_matrix/ and filtered_feature_bc_matrix/ from it
    //              (load + emptyDrops + CellRanger classification all happen in one pass).

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    // IMPORTANT: Output filename MUST include the step name
    // (e.g. sample.LOAD_CLASSIFY.rds)
    // Nextflow stages inputs as symlinks with the same basename into each process work dir.
    // If two steps both output "sample.rds", the downstream process gets a symlink named
    // "sample.rds" pointing back to the upstream work dir. When R saves to "./sample.rds"
    // it follows the symlink and overwrites the upstream file — corrupting the cache and
    // causing cache misses on every resume. Unique filenames per step prevent this entirely.
    //
    // sample_dir is NOT re-emitted downstream — CALC_QC_AND_IDENTIFY_SINGLETS and everything
    // after it only need the Seurat object, not the raw CellRanger directory.

    tuple val(sample_id), path("${sample_id}.LOAD_CLASSIFY.rds"),         emit: sample_rds
    path("${sample_id}.SCRNA_LOAD_CELLRANGER_AND_CLASSIFY_DROPLETS.log"),   emit: log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def script_path = "${workflow.projectDir}/modules"
    def LOG         = "${sample_id}.SCRNA_LOAD_CELLRANGER_AND_CLASSIFY_DROPLETS.log"
    def gene_column = params.scrna.read10x.gene_column ?: 2

    """
    Rscript "${script_path}/SCRNA/SCRNA_load_cellranger_and_classify_droplets.R" \
        "${sample_id}" \
        "${sample_dir}" \
        "${gene_column}" \
        "." \
        > "${LOG}" 2>&1

    EXIT_CODE=\$?
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "❌ ERROR: SCRNA_LOAD_CELLRANGER_AND_CLASSIFY_DROPLETS failed for ${sample_id}" | tee -a "${LOG}" >&2
        exit \$EXIT_CODE
    fi

    echo "✅ SUCCESS: SCRNA_LOAD_CELLRANGER_AND_CLASSIFY_DROPLETS completed for ${sample_id}" >> "${LOG}"
    """
}
