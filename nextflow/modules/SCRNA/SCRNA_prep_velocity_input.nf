// =============================================================================
// PROCESS: SCRNA_PREP_VELOCITY_INPUT
// =============================================================================
// Purpose: Subsets the merged kb-python AnnData (all samples) down to one
//          lineage's barcodes, using that lineage's own _Final Seurat
//          barcodes.csv (CellType, cluster, UMAP — written by
//          SCRNA_ANNOTATE_CLUSTERS alongside annotated_seurat.rds),
//          producing a lineage-specific velocity_input.h5ad ready for
//          scVelo/veloVI.
//
// Runs once per lineage (parallel), each task reusing the SAME merged.h5ad
// (broadcast via value-channel semantics in scrnaseq.nf) against its OWN
// lineage's barcodes.csv.
//
// Barcode overlap between Seurat and merged h5ad is expected to be <100%
// (kb-python's --filter bustools and Cell Ranger's cell-calling are
// independent algorithms) — handled via strict inner join inside the script,
// with overlap % logged as a diagnostic.
// =============================================================================

process SCRNA_PREP_VELOCITY_INPUT {

    tag "Prepping velocity input: ${label}"
    label 'process_high'
    label 'omics_py'

    publishDir = [
        [path: { "${params.proj_dir()}/${params.species}/06.scVelo/${label}" }, mode: 'copy', pattern: "velocity_input.h5ad"],
        [path: { "${params.proj_dir()}/${params.species}/Logs" },                 mode: 'copy', pattern: "*.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(label), path(barcodes_csv), path(merged_h5ad)
    // label        : base lineage label, e.g. "Epithelial" (stripped of "_Final"
    //                upstream in scrnaseq.nf)
    // barcodes_csv : SCRNA_ANNOTATE_CLUSTERS.out.barcodes for this lineage's
    //                _Final branch — barcode, CellType, cluster, UMAP_1, UMAP_2
    // merged_h5ad  : SCRNA_MERGE_H5AD.out.merged — the SAME file broadcast to
    //                every lineage task (single-item channel, reused across
    //                all parallel calls, standard Nextflow value-channel behavior)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(label), path("velocity_input.h5ad"),        emit: velocity_input
    path("SCRNA_PREP_VELOCITY_INPUT.${label}.log"),       emit: log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def script_path = "${workflow.projectDir}/modules"
    def LOG          = "SCRNA_PREP_VELOCITY_INPUT.${label}.log"

    """
    python3 "${script_path}/SCRNA/SCRNA_prep_velocity_input.py" \
        --merged "${merged_h5ad}" \
        --barcodes "${barcodes_csv}" \
        --output velocity_input.h5ad > "${LOG}" 2>&1

    EXIT_CODE=\$?
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "❌ ERROR: SCRNA_PREP_VELOCITY_INPUT failed for ${label}" | tee -a "${LOG}" >&2
        exit \$EXIT_CODE
    fi
    echo "✅ SUCCESS: SCRNA_PREP_VELOCITY_INPUT completed for ${label}" >> "${LOG}"
    """
}
