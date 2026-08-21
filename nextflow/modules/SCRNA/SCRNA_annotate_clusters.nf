process SCRNA_ANNOTATE_CLUSTERS {

    tag "Annotating clusters: ${label}"
    label 'process_heavy'
    label 'omics_r'

    publishDir = [
        [path: { "${params.proj_dir()}/${params.species}/05.Seurat/${label}" },                  mode: 'copy', pattern: "*.rds"],
        [path: { "${params.proj_dir()}/${params.species}/05.Seurat/${label}" },                  mode: 'copy', pattern: "barcodes.csv"],
        [path: { "${params.proj_dir()}/${params.species}/05.Seurat/${label}/Plots/Annotation" }, mode: 'copy', pattern: "*.pdf"],
        [path: { "${params.proj_dir()}/${params.species}/Logs" },                                mode: 'copy', pattern: "*.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(label), path(seurat_rds), path(annotation_xlsx), val(cluster_col)
    // Bundled together (not separate channels) because every one of them varies PER
    // LABEL -- bundling is what lets several labels flow through this ONE process call
    // concurrently without mispairing a label with the wrong file/value.
    // label            : "Whole" or a cell-type label -- every label gets an identical
    //                    output shape under 05.Seurat/{label}/. Also determines which
    //                    metadata column the R script writes: "Whole" -> Lineage,
    //                    anything else -> CellType (see R script Section 3).
    // seurat_rds       : clustered_seurat.rds from SCRNA_CLUSTER_AND_FINDMARKERS --
    //                    ALL cells, CellType already pre-filled to "Exclude" for
    //                    barcodes its automated junk/sparse detection flagged.
    // annotation_xlsx  : Annotation_{cluster_col}_{label}.xlsx from
    //                    SCRNA_CLUSTER_AND_FINDMARKERS, with the Lineage (Whole) or
    //                    CellType (subset) column filled in by hand. Every cluster in
    //                    cluster_col must have a non-blank label, or this process
    //                    fails fast with the list of missing ones. Clusters the user
    //                    judges to be junk/mixed should be labeled with exclude_label
    //                    (default "Exclude").
    //                    Barcodes CLUSTER_AND_FINDMARKERS already pre-flagged
    //                    "Exclude" keep that label regardless of what the
    //                    template says for their cluster -- see R script Section 3.
    // cluster_col      : which resolution the annotation_xlsx was built for,
    //                    e.g. "harmony_res0.8". Must match exactly -- this is
    //                    how cluster labels get mapped back onto cells.
    val(exclude_label)
    // exclude_label         : reserved label string marking junk/mixed clusters
    //                         (and individually pre-flagged auto-junk barcodes).
    //                         Default: "Exclude". Physical removal happens once,
    //                         downstream, in SCRNA_ADD_CELLTYPE -- not here.
    //                         Passed from params.scrna.seurat.exclude_label. Same
    //                         value for every label -- not bundled into the tuple above.

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(label), path("annotated_seurat.rds"),        emit: seurat_rds
    // label echoed back so downstream processes (SCRNA_ADD_CELLTYPE, etc.) can keep
    // pairing the right file with the right label. ALL cells are in this one
    // object -- Exclude-labeled cells (manual + auto-flagged) are NOT removed
    // here. No clean variant: CellType isn't truly final until SCRNA_ADD_CELLTYPE
    // assembles every branch's own CellType calls, so clean is computed there,
    // exactly once, instead of at every annotate call site. Candidate-marker
    // review generation also lives there now, not here (see SCRNA_ADD_CELLTYPE).
    tuple val(label), path("barcodes.csv"),                 emit: barcodes
    // Lightweight barcode/CellType(or Lineage)/cluster/UMAP CSV, written on
    // EVERY call regardless of label. Only "_Final" branches' CSVs are
    // actually consumed downstream (SCRNA_PREP_VELOCITY_INPUT via
    // lineage_ch) -- Whole/_Initial CSVs are simply unused, same as *.pdf
    // being emitted unconditionally below.
    path("*.pdf"),                       optional: true,        emit: plots
    path("SCRNA_ANNOTATE_CLUSTERS.${label}.log"),                        emit: log

    // =================================================================================
    // EXECUTION
    // NOTE: This process is NOT chained automatically after
    // SCRNA_CLUSTER_AND_FINDMARKERS -- annotation_xlsx doesn't exist until the user
    // fills it in by hand. The main workflow gates inclusion of this process on
    // annotation_xlsx existing on disk (checked before the DAG is built), so a first
    // run completes through clustering with this step skipped, and a second run
    // (after the xlsx is filled in and its path added to yaml) picks it up while
    // every upstream step hits cache untouched.
    // One object saved: annotated_seurat.rds -- all cells, including every
    // Exclude-labeled one, for visual audit on UMAP. For "Whole", the annotation
    // lands in metadata column Lineage; for every other label, in CellType.
    // =================================================================================
    script:

    def script_path = "${workflow.projectDir}/modules"
    def LOG         = "SCRNA_ANNOTATE_CLUSTERS.${label}.log"

    """
    Rscript "${script_path}/SCRNA/SCRNA_annotate_clusters.R" \
        "${seurat_rds}" \
        "${label}" \
        "${annotation_xlsx}" \
        "${cluster_col}" \
        "${exclude_label}" \
        "." > "${LOG}" 2>&1

    EXIT_CODE=\$?
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "❌ ERROR: SCRNA_ANNOTATE_CLUSTERS failed for ${label}" | tee -a "${LOG}" >&2
        exit \$EXIT_CODE
    fi

    echo "✅ SUCCESS: SCRNA_ANNOTATE_CLUSTERS completed for ${label}" >> "${LOG}"
    """
}
