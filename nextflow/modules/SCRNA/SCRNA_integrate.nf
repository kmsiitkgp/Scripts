process SCRNA_INTEGRATE {
    tag "Harmony integration: ${label}"
    label 'process_heavy'
    label 'omics_r'
    publishDir = [
        [path: { "${params.proj_dir()}/${params.species}/05.Seurat/${label}" },                   mode: 'copy', pattern: "*.rds"],
        [path: { "${params.proj_dir()}/${params.species}/05.Seurat/${label}/Plots/Integration" }, mode: 'copy', pattern: "*.pdf"],
        [path: { "${params.proj_dir()}/${params.species}/Logs" },                                 mode: 'copy', pattern: "*.log"]
    ]
    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(label), path(seurat_rds)
    // label      : "Whole" or a cell-type label — bundled with seurat_rds (not a separate
    //              channel) so multiple labels can flow through this ONE process call
    //              concurrently without mispairing.
    // seurat_rds : sct_seurat.rds from SCRNA_SCT_TRANSFORM — SCTransformed object
    //              with sct_pca reduction and @misc$sct_params populated.
    //              Harmony integration method and n_dims are set internally —
    //              no additional parameters needed from Nextflow.
    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(label), path("integrated_seurat.rds"),   emit: seurat_rds
    // label echoed back unchanged so downstream processes can keep pairing the
    // right file with the right label without a separate lookup.
    path("*.pdf"),  optional: true,                    emit: plots
    path("SCRNA_INTEGRATE.${label}.log"),              emit: log
    // =================================================================================
    // EXECUTION
    // NOTE: Harmony operates on PCA embeddings — fast and memory-safe at 500k+ cells.
    // A QC UMAP (umap_harmony_qc) is saved here for integration validation only.
    // The final UMAP used downstream is produced by SCRNA_CLUSTER_AND_FINDMARKERS
    // on the clean dataset after junk cluster removal.
    // =================================================================================
    script:

    def script_path = "${workflow.projectDir}/modules"
    def LOG         = "SCRNA_INTEGRATE.${label}.log"

    """
    Rscript "${script_path}/SCRNA/SCRNA_integrate.R" \
        "${seurat_rds}" \
        "." > "${LOG}" 2>&1

    EXIT_CODE=\$?
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "❌ ERROR: SCRNA_INTEGRATE failed for ${label}" | tee -a "${LOG}" >&2
        exit \$EXIT_CODE
    fi

    echo "✅ SUCCESS: SCRNA_INTEGRATE completed for ${label}" >> "${LOG}"
    """
}
