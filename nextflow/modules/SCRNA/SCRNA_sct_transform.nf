process SCRNA_SCT_TRANSFORM {

    tag "SCTransform on merged object: ${label}"
    label 'process_heavy'
    label 'omics_r'

    publishDir = [
        [path: { "${params.proj_dir()}/${params.species}/05.Seurat/${label}" }, mode: 'copy', pattern: "*.rds"],
        [path: { "${params.proj_dir()}/${params.species}/Logs" },               mode: 'copy', pattern: "*.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(label), path(seurat_rds)
    // label      : "Whole" or a cell-type label (e.g. "Epithelial") — bundled with
    //              seurat_rds (not a separate channel) so multiple labels can flow
    //              through this ONE process call concurrently without mispairing —
    //              e.g. when several subcluster branches run in the same execution.
    // seurat_rds : filtered_seurat.rds — merged object of high quality singlets.
    //              "Whole" branch: from MERGE_AND_PLOT_QC. Subcluster branches:
    //              from SCRNA_SUBSET_CELLS. Identical entry-point contract either
    //              way (raw RNA counts + basic QC metadata, nothing computed).
    path(cellcycle_marker_xlsx)
    // cellcycle_marker_xlsx : Cell_Cycle_Markers.xlsx — columns Phase, Mouse_Gene, Human_Gene
    //                 Species filtering is automatic via intersect with present features
    //                 Shared across every label — not bundled into the tuple above.

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(label), path("sct_seurat.rds"),   emit: seurat_rds
    // label echoed back unchanged so downstream processes (INTEGRATE, etc.) can
    // keep pairing the right file with the right label without a separate lookup.
    path("SCRNA_SCT_TRANSFORM.${label}.log"),            emit: log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def script_path = "${workflow.projectDir}/modules"
    def LOG         = "SCRNA_SCT_TRANSFORM.${label}.log"

    """
    Rscript "${script_path}/SCRNA/SCRNA_sct_transform.R" \
        "${seurat_rds}" \
        "${cellcycle_marker_xlsx}" \
        "." > "${LOG}" 2>&1

    EXIT_CODE=\$?
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "❌ ERROR: SCRNA_SCT_TRANSFORM failed for ${label}" | tee -a "${LOG}" >&2
        exit \$EXIT_CODE
    fi

    echo "✅ SUCCESS: SCRNA_SCT_TRANSFORM completed for ${label}" >> "${LOG}"
    """
}
