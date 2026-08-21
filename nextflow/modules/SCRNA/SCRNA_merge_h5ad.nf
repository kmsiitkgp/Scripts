// =============================================================================
// PROCESS: SCRNA_MERGE_H5AD
// =============================================================================
// Purpose: Merges all per-sample kb-python h5ad files (04.KB/{sample_id}/*.h5ad)
//          into a single AnnData, with barcodes rewritten to match the Seurat
//          object's convention ({sample_id}_{barcode}-1), so downstream lineage
//          subsetting (SCRNA_PREP_VELOCITY_INPUT) can match against Seurat
//          barcodes directly.
//
// Runs ONCE per pipeline run (all samples -> one merged h5ad), unlike
// SCRNA_PREP_VELOCITY_INPUT which runs once per lineage off this shared output.
// =============================================================================

process SCRNA_MERGE_H5AD {

    tag "Merging ${sample_ids.size()} sample h5ad files"
    label 'process_heavy'
    label 'omics_py'

    publishDir = [
        [path: { "${params.proj_dir()}/${params.species}/06.scVelo/Whole" }, mode: 'copy', pattern: "merged.h5ad"],
        [path: { "${params.proj_dir()}/${params.species}/Logs" },              mode: 'copy', pattern: "*.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(sample_ids), path(h5ad_files)
    // sample_ids : List<String>  — e.g. ["GHRKO24m_1", "GHRKO24m_2", ...]
    // h5ad_files : List<Path>    — matching per-sample .h5ad files, same order as
    //              sample_ids (built via SCRNA_KB_COUNT.out.h5ad_files.collect(flat:false)
    //              in scrnaseq.nf, same "keep paired, sort together" pattern already
    //              used for subset_prepared_ch)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("merged.h5ad"),               emit: merged
    path("SCRNA_MERGE_H5AD.log"),      emit: log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def script_path = "${workflow.projectDir}/modules"
    def LOG          = "SCRNA_MERGE_H5AD.log"

    // Build "SAMPLE_ID=path.h5ad" args, pairing sample_ids[i] with h5ad_files[i]
    def sample_args = [sample_ids, h5ad_files].transpose()
        .collect { sid, f -> "${sid}=${f}" }
        .join(' ')

    """
    python3 "${script_path}/SCRNA/SCRNA_merge_h5ad.py" \
        --output merged.h5ad \
        ${sample_args} > "${LOG}" 2>&1

    EXIT_CODE=\$?
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "❌ ERROR: SCRNA_MERGE_H5AD failed" | tee -a "${LOG}" >&2
        exit \$EXIT_CODE
    fi
    echo "✅ SUCCESS: SCRNA_MERGE_H5AD completed" >> "${LOG}"
    """
}
