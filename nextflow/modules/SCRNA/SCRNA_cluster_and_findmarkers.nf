process SCRNA_CLUSTER_AND_FINDMARKERS {

    tag "Clustering, marker finding, junk flagging: ${label}"
    label 'process_heavy'
    label 'omics_r'

    publishDir = [
        [path: { "${params.proj_dir()}/${params.species}/05.Seurat/${label}" },                  mode: 'copy', pattern: "*.rds"],
        [path: { "${params.proj_dir()}/${params.species}/05.Seurat/${label}" },                  mode: 'copy', pattern: "FindAllMarkers"],
        [path: { "${params.proj_dir()}/${params.species}/05.Seurat/${label}" },                  mode: 'copy', pattern: "Annotation_Templates"],
        [path: { "${params.proj_dir()}/${params.species}/05.Seurat/${label}/Plots/Clustering" }, mode: 'copy', pattern: "*.pdf"],
        [path: { "${params.proj_dir()}/${params.species}/Logs" },                                mode: 'copy', pattern: "*.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(label), path(seurat_rds)
    // label      : "Whole" or a cell-type label — bundled with seurat_rds (not a separate
    //              channel) so multiple labels can flow through this ONE process call
    //              concurrently without mispairing.
    // seurat_rds : integrated_seurat.rds from SCRNA_INTEGRATE — object with
    //              integ_harmony reduction and @misc$integration_params populated.
    //              FindNeighbors, FindClusters, RunUMAP, FindAllMarkers all run
    //              internally, exactly ONCE (no iterative junk-removal-and-
    //              recluster loop — junk/sparse cells are flagged via
    //              CellType="Exclude", not physically removed here). Identical
    //              logic for Whole and every subset branch.
    val(resolutions)
    // resolutions      : comma-separated string of Leiden resolutions to sweep.
    //                    e.g. "0.2,0.4,0.6,0.8,1.0,1.2" — same for every label,
    //                    not bundled into the tuple above.
    //                    Markers are found at every resolution — Excel files saved
    //                    per resolution for manual review.
    //                    Passed from params.scrna.seurat.resolutions.
    val(sparse_threshold)
    // sparse_threshold : clusters with <= this many cells at any resolution are
    //                    flagged CellType="Exclude" — catches micro-clusters
    //                    that form at high resolutions and have no biological meaning.
    //                    Default: 5. Passed from params.scrna.seurat.sparse_threshold.
    path(celltype_marker_xlsx)
    // celltype_marker_xlsx : same file as params.scrna.seurat.celltype_marker_file,
    //                    used by SCRNA_CALC_MODULE_SCORES downstream. Used here
    //                    only to build the Reference_Lists sheet (Sheet 2,
    //                    CellType column) in each Annotation_{col}_{label}.xlsx —
    //                    copy-paste reference only, no scoring/annotation logic
    //                    runs here.
    //                    Same file for every label — not bundled into the tuple above.

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(label), path("clustered_seurat.rds"),          emit: seurat_rds
    // label echoed back unchanged so downstream processes can keep pairing the
    // right file with the right label without a separate lookup. ALL cells are
    // in this one object — junk/sparse cells are flagged CellType="Exclude",
    // not removed. No separate full/clean/junk variants any more.
    path("FindAllMarkers/"),                 optional: true,   emit: qc_folders
    path("Annotation_Templates/"),           optional: true,   emit: annotation_templates
    path("*.pdf"),                           optional: true,   emit: plots
    path("SCRNA_CLUSTER_AND_FINDMARKERS.${label}.log"),        emit: log

    // =================================================================================
    // EXECUTION
    // NOTE: One FindNeighbors/FindClusters/FindAllMarkers pass per object — no
    // iterative junk-removal-and-recluster loop (that design was retired: junk
    // cells are now flagged via CellType="Exclude" in the single pass's
    // clustering, not physically removed and reclustered). Identical logic for
    // Whole and every subset branch — no more Whole-only special casing.
    // Physical Exclude removal happens once, downstream, in SCRNA_ADD_CELLTYPE,
    // after manual annotation (SCRNA_ANNOTATE_CLUSTERS) has finalized
    // CellType for every cluster.
    // Outputs are all saved under a single flat FindAllMarkers/ folder:
    //   FindAllMarkers/Markers.{resolution}.xlsx (one per resolution)
    //   FindAllMarkers/Junk_Cluster_Audit.xlsx (per-cluster junk/sparse status summary)
    // Annotation_Templates/Annotation_{cluster_col}_{label}.xlsx (one per
    // resolution) is generated once, after this single clustering/marker pass.
    // Sheet 1's editable column is "Lineage" for label=="Whole", "CellType"
    // for every other label.
    // =================================================================================
    script:

    def script_path = "${workflow.projectDir}/modules"
    def LOG         = "SCRNA_CLUSTER_AND_FINDMARKERS.${label}.log"
    def resolutions_str = (resolutions instanceof List) \
        ? resolutions.join(',') \
        : (resolutions ?: "0.2,0.4,0.6,0.8,1.0,1.2")

    """
    Rscript "${script_path}/SCRNA/SCRNA_cluster_and_findmarkers.R" \
        "${seurat_rds}" \
        "${label}" \
        "${resolutions_str}" \
        "${sparse_threshold}" \
        "${celltype_marker_xlsx}" \
        "." > "${LOG}" 2>&1

    EXIT_CODE=\$?
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "❌ ERROR: SCRNA_CLUSTER_AND_FINDMARKERS failed for ${label}" | tee -a "${LOG}" >&2
        exit \$EXIT_CODE
    fi

    echo "✅ SUCCESS: SCRNA_CLUSTER_AND_FINDMARKERS completed for ${label}" >> "${LOG}"
    """
}
