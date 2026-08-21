process SCRNA_CALC_MODULE_SCORES {

    tag "Calculating module scores: ${subset_label}"
    label 'process_heavy'
    label 'omics_r'

    publishDir = [
        [path: { "${params.proj_dir()}/${params.species}/04.Seurat/${subset_label}" },    mode: 'copy', pattern: "*.rds"],
        [path: { "${params.proj_dir()}/${params.species}/04.Seurat/${subset_label}/QC" }, mode: 'copy', pattern: "*.pdf"],
        [path: { "${params.proj_dir()}/${params.species}/Logs" },                         mode: 'copy', pattern: "*.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    path(rds_file)
    // rds_file    : clustered_seurat.rds from SCRNA_CLUSTER_AND_FINDMARKERS — clean
    //               object with cluster columns (harmony_res0.2 ... harmony_res1.2)
    //               in metadata and @misc$clustering_params populated.
    //               Module scores are cell-level and resolution-independent —
    //               computed once on the full object, available for annotation
    //               at any resolution.
    path(marker_file)
    // marker_file : Excel file where each column = one broad cell type and each
    //               row = a marker gene for that cell type. Column headers become
    //               module names. NAs pad shorter columns to equal length.
    //               e.g. BPlasma, Dendritic, Endothelial, Epithelial, TNK...
    //               Both human and mouse gene symbols are supported — matching
    //               is case-insensitive to handle capitalisation differences.
    val(subset_label)
    // subset_label : "Whole" or a cell-type label (e.g. "Epithelial") — every label
    //               gets an identical output shape under 04.Seurat/{label}/.
    //               Determines output location only.

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("scored_seurat.rds"),              emit: rds
    path("*.pdf"),  optional: true,         emit: plots
    path("SCRNA_CALC_MODULE_SCORES.log"),   emit: log

    // =================================================================================
    // EXECUTION
    // NOTE: UCell scoring is notably slow on large datasets (20-40 min for 100k cells).
    // scored_seurat.rds is cached by Nextflow so annotation parameters can be tuned
    // by rerunning only SCRNA_ANNOTATE_CLUSTERS without repeating UCell scoring.
    // =================================================================================
    script:

    def script_path = "${workflow.projectDir}/modules"
    def LOG         = "SCRNA_CALC_MODULE_SCORES.log"

    """
    Rscript "${script_path}/SCRNA_calc_module_scores.R" \
        "${rds_file}" \
        "${marker_file}" \
        "." > "${LOG}" 2>&1

    EXIT_CODE=\$?
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "❌ ERROR: SCRNA_CALC_MODULE_SCORES failed" | tee -a "${LOG}" >&2
        exit \$EXIT_CODE
    fi

    echo "✅ SUCCESS: SCRNA_CALC_MODULE_SCORES completed" >> "${LOG}"
    """
}
