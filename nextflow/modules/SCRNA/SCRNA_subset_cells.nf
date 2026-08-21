process SCRNA_SUBSET_CELLS {

    tag "Subsetting cells: ${label}"
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
    path(seurat_rds)
    // seurat_rds      : annotated_seurat.rds from whichever label is upstream of this
    //                   call — "Whole" for the "_Initial" pass (must have a populated
    //                   Lineage column), or a "{label}_Initial" branch's own annotated
    //                   object for the "_Final" pass (must have a populated CellType
    //                   column, values Retain/Exclude). Must have filter_col populated.
    //                   For the "_Initial" alias, ONE object is shared across every
    //                   subset request this run (Whole, broadcast); for the "_Final"
    //                   alias, each request has its OWN source object (that branch's own
    //                   "_Initial" annotation) — either way it's not bundled into the
    //                   tuple below, since filter_value/label genuinely differ per
    //                   request while filter_col is constant for a given process call.
    tuple val(filter_value), val(label)
    // filter_value    : comma-separated values to keep, e.g. "Epithelial" or
    //                   "T Cell,B Cell,Myeloid" to merge several into one branch
    //                   (when filter_col="Lineage"), or "Retain" (when
    //                   filter_col="CellType", second-pass mixed-population cleanup).
    // label           : derived folder/branch name for this subset, e.g. "Epithelial_Initial"
    //                   or "Epithelial_Final". Computed in the calling workflow (see
    //                   label_for closure / "_Initial"/"_Final" suffixing) so the folder
    //                   name always reflects exactly what was filtered and which pass.
    // Bundled together (not separate channels) so multiple subset requests can flow
    // through this ONE process call concurrently without mispairing a label with the
    // wrong filter string.
    val(filter_col)
    // filter_col      : which metadata column to filter on — "Lineage" for the
    //                   Whole -> {label}_Initial call, "CellType" for the
    //                   {label}_Initial -> {label}_Final call (second-pass mixed-
    //                   population removal). Same value for every request in a given
    //                   process call, so broadcast as a plain val, not bundled into the
    //                   tuple above — same convention as exclude_label in
    //                   SCRNA_ANNOTATE_CLUSTERS.

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(label), path("filtered_seurat.rds"), optional: true,   emit: seurat_rds
    // optional: true — a filter_value value with no matching cells in the object
    // SKIPS (produces no output for that label) rather than erroring out. Handled in
    // SCRNA_subset_cells.R Section 2: zero matched_types logs a warning and returns
    // early without writing filtered_seurat.rds, so this optional:true correctly
    // reflects "genuinely absent lineage" as a non-fatal skip.
    path("SCRNA_SUBSET_CELLS.${label}.log"),                                  emit: log

    // =================================================================================
    // EXECUTION
    // NOTE: Output is always named "filtered_seurat.rds" — identical entry-point name
    // across every label's folder (05.Seurat/{label}/filtered_seurat.rds), same
    // convention as MERGE_AND_PLOT_QC produces for "Whole". All Whole-pipeline
    // computed state (SCT assay, PCA/Harmony reductions, cluster columns, module
    // scores, CellType, Stability_Score, @misc provenance) is stripped inside the R
    // script — see SCRNA_subset_cells.R Sections 2 and 4 — so this object is safe to
    // feed directly into SCRNA_SCT_TRANSFORM as if it were a fresh merged object.
    // =================================================================================
    script:

    def script_path = "${workflow.projectDir}/modules"
    def LOG         = "SCRNA_SUBSET_CELLS.${label}.log"

    """
    Rscript "${script_path}/SCRNA/SCRNA_subset_cells.R" \
        "${seurat_rds}" \
        "${filter_col}" \
        "${filter_value}" \
        "${label}" \
        "." > "${LOG}" 2>&1

    EXIT_CODE=\$?
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "❌ ERROR: SCRNA_SUBSET_CELLS failed for ${label}" | tee -a "${LOG}" >&2
        exit \$EXIT_CODE
    fi

    echo "✅ SUCCESS: SCRNA_SUBSET_CELLS completed for ${label}" >> "${LOG}"
    """
}
