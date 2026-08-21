process SCRNA_ANNOTATE_CLUSTERS_AUTO {

    tag "Annotating clusters: ${subset_label}"
    label 'process_heavy'
    label 'omics_r'

    publishDir = [
        [path: { "${params.proj_dir()}/${params.species}/04.Seurat/${subset_label}" },      mode: 'copy', pattern: "*.rds"],
        [path: { "${params.proj_dir()}/${params.species}/04.Seurat/${subset_label}/Auto" }, mode: 'copy', pattern: "*.pdf"],
        [path: { "${params.proj_dir()}/${params.species}/04.Seurat/${subset_label}/Auto" }, mode: 'copy', pattern: "*.xlsx"],
        [path: { "${params.proj_dir()}/${params.species}/Logs" },                           mode: 'copy', pattern: "*.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    path(rds_file)
    // rds_file            : scored_seurat.rds from SCRNA_CALC_MODULE_SCORES — object
    //                       with Seurat module score columns (e.g. T.Cell1, B.Cell2...)
    //                       and UCell score columns (e.g. T.Cell_UCell...) in metadata,
    //                       and @misc$scoring_params with exact column names.
    //                       Annotation runs across ALL resolutions equally — no single
    //                       resolution is trusted over others. The final CellType is the
    //                       mode label across all method × resolution combinations.
    val(avg_prob_threshold)
    // avg_prob_threshold  : minimum average probability for the winning cell type
    //                       across cells in a cluster to be a confident assignment.
    //                       Clusters below this → Unknown.
    //                       Default: 0.3 (winner must be ~4.5x above uniform baseline
    //                       of 1/15 = 0.067 for 15 cell types).
    val(delta_threshold)
    // delta_threshold     : minimum difference between top-1 and top-2 cell type
    //                       probabilities at the cluster level. Catches cases where
    //                       two cell types are nearly tied despite passing
    //                       avg_prob_threshold. Default: 0.1.
    val(stability_threshold)
    // stability_threshold : minimum fraction of method × resolution combinations
    //                       that must agree with the mode CellType for a cell to be
    //                       included in the clean output object.
    //                       Default: 0.9 (90% agreement across 12 combinations).
    //                       Full annotated object retains all cells regardless.
    val(subset_label)
    // subset_label        : "Whole" or a cell-type label (e.g. "Epithelial") — every
    //                       label gets an identical output shape under
    //                       04.Seurat/{label}/. Determines output location only.

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("annotated_seurat_auto.rds"),               emit: rds
    path("annotated_seurat_auto_clean.rds"),         emit: rds_clean
    path("*.pdf"),  optional: true,                  emit: plots
    path("Cluster_Annotation_Summary.xlsx"),         emit: tables
    path("SCRNA_ANNOTATE_CLUSTERS_AUTO.log"),        emit: log

    // =================================================================================
    // EXECUTION
    // NOTE: Annotation logic (Z-score, softmax, cluster averaging, mode voting) is
    // fast compared to UCell scoring. Parameters can be tuned and this process rerun
    // without rerunning SCRNA_CALC_MODULE_SCORES thanks to Nextflow caching.
    // Two objects saved:
    //   annotated_seurat_auto.rds       — all cells, full metadata, for exploration
    //   annotated_seurat_auto_clean.rds — high-confidence cells only (Stability_Score
    //                                >= stability_threshold), for DE analysis
    // =================================================================================
    script:

    def script_path = "${workflow.projectDir}/modules"
    def LOG         = "SCRNA_ANNOTATE_CLUSTERS_AUTO.log"

    """
    Rscript "${script_path}/SCRNA_annotate_clusters_auto.R" \
        "${rds_file}" \
        "." \
        "${avg_prob_threshold}" \
        "${delta_threshold}" \
        "${stability_threshold}" > "${LOG}" 2>&1

    EXIT_CODE=\$?
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "❌ ERROR: SCRNA_ANNOTATE_CLUSTERS_AUTO failed" | tee -a "${LOG}" >&2
        exit \$EXIT_CODE
    fi

    echo "✅ SUCCESS: SCRNA_ANNOTATE_CLUSTERS_AUTO completed" >> "${LOG}"
    """
}
