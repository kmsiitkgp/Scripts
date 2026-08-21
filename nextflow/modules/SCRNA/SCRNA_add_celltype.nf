// =========================================================================================
// PROCESS: SCRNA_ADD_CELLTYPE
// =========================================================================================
//
// Purpose:
//   Assemble the final CellType annotation by combining:
//     1. Whole-object lineage annotation (coarse)
//     2. Fine-grained annotations from completed subset branches
//
// Workflow position:
//
//   Whole annotated_seurat.rds
//              +
//   Epithelial_Final / Immune_Final / etc.
//              ↓
//   SCRNA_ADD_CELLTYPE
//              ↓
//   Final CellType object + markers + plots
//
// Key behavior:
//   - Whole CellType starts from Lineage.
//   - Subset annotations overwrite CellType for covered cells.
//   - Cells not covered by subsets retain their original Lineage.
//   - Excluded cells are removed only when creating the clean object.
//
// This is a fan-in process: one task receives multiple annotated subset objects.
// =========================================================================================


process SCRNA_ADD_CELLTYPE {

    tag "Assembling final CellType: Whole + subset branches"

    label 'process_heavy'
    label 'omics_r'

    publishDir = [
        [path: { "${params.proj_dir()}/${params.species}/05.Seurat/Whole" },        mode: 'copy', pattern: "*.rds"],
        [path: { "${params.proj_dir()}/${params.species}/05.Seurat/Whole" },        mode: 'copy', pattern: "Markers.CellType_Final.xlsx"],
        [path: { "${params.proj_dir()}/${params.species}/05.Seurat/Whole" },        mode: 'copy', pattern: "CellType_Markers_Review.xlsx"],
        [path: { "${params.proj_dir()}/${params.species}/05.Seurat/Whole/Plots" },  mode: 'copy', pattern: "*.pdf"],
        [path: { "${params.proj_dir()}/${params.species}/Logs" },                   mode: 'copy', pattern: "*.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(whole_label), path(whole_rds)
    // whole_label : always "Whole" — kept as a tuple element for tag/log consistency
    //               with every other process in this chain, not because it varies.
    // whole_rds   : annotated_seurat.rds (the FULL object, not a clean
    //               variant — that no longer exists) from the "Whole" branch.
    //               Lineage here is already Whole's own coarse call (Epithelial/
    //               Mesenchymal/Hematopoietic/...), written directly by
    //               SCRNA_ANNOTATE_CLUSTERS — used only as the fallback for
    //               cells no subset branch covers.
    val(subset_labels)
    // subset_labels   : comma-joined list of every subset branch's label that
    //                   successfully finished annotation this run, e.g.
    //                   "Epithelial,Immune". Order matches subset_rds_list.
    //                   Can be empty (first-ever run, nothing subclustered yet) —
    //                   every cell then just keeps Lineage as its final CellType.
    path(subset_rds_list, stageAs: "subset_input_*/*")
    // subset_rds_list : list of annotated_seurat.rds paths (the FULL object
    //                   per branch, not a clean variant), one per subset branch —
    //                   order matches subset_labels. Using the FULL object (not
    //                   clean) is what lets a subset branch's own "Exclude" calls
    //                   correctly propagate back onto Whole, instead of silently
    //                   vanishing had a pre-filtered clean object been used.
    //                   Declared as path() (not buried in a val()) so Nextflow
    //                   stages every file into the task work dir correctly.
    //                   stageAs: "subset_input_*/*" is REQUIRED here, not optional —
    //                   every branch's SCRNA_ANNOTATE_CLUSTERS output is named
    //                   identically ("annotated_seurat.rds"), so staging 2+ of them
    //                   into the task root at once is a guaranteed Nextflow "input
    //                   file name collision" (task never gets submitted). This
    //                   pattern auto-numbers each into its own subdirectory
    //                   (subset_input_1/annotated_seurat.rds, subset_input_2/..., ...)
    //                   so no two staged paths collide.
    //                   NOTE: this input only resolves once every requested subset
    //                   branch that's actually ready has finished (built via
    //                   .collect() upstream in scrnaseq.nf) — SCRNA_ADD_CELLTYPE is
    //                   a genuine fan-in point, not something that runs per-label
    //                   like the rest of the chain.
    val(reduction)
    // reduction : name of the pre-UMAP embedding to reuse for the final UMAP
    //             recompute (default "integ_harmony" — the Harmony-corrected PCA
    //             space from SCRNA_INTEGRATE, still present on the Whole object).
    //             Same value regardless of label — not bundled into the tuple above.
    //             Passed from params.scrna.seurat.final_umap_reduction (or similar).

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("annotated_seurat_final.rds"),                            emit: seurat_rds
    path("annotated_seurat_final_clean.rds"),                      emit: seurat_rds_clean
    path("Markers.CellType_Final.xlsx"),  optional: true,          emit: final_markers
    path("CellType_Markers_Review.xlsx"), optional: true,          emit: markers_review
    path("*.pdf"),                        optional: true,          emit: plots
    path("SCRNA_ADD_CELLTYPE.log"),                                emit: log

    // =================================================================================
    // EXECUTION
    // NOTE: For every cell in whole_rds:
    //   Lineage   : Whole's original coarse CellType call — permanent record, never
    //               overwritten again.
    //   CellType  : starts as a copy of Lineage, then overwritten per subset branch
    //               with THAT branch's own fine-grained CellType call (including its
    //               own "Exclude" calls). Cells never covered by any subset branch
    //               this run keep Lineage as their final CellType. NO majority-vote/
    //               consensus check against Lineage — a subset branch's marker-
    //               informed reclassification is trusted outright (see project
    //               discussion: the coarser Whole-level clustering has LESS
    //               information than a subset branch's own re-annotation, not more —
    //               disagreement isn't "ambiguous," it's often a correction).
    //   CellType_Source : provenance per cell — "Lineage" (fallback) or the subset
    //               branch label that provided the final call. Informational only.
    // Clean object (Exclude removed by the FINAL CellType) is computed HERE for the
    // first time in the whole pipeline. On the clean object: NO reclustering
    // (CellType is the final grouping, a fresh cluster ID would be redundant) — just
    // one FindAllMarkers(group.by="CellType") pass, and a RunUMAP-only recompute (no
    // FindNeighbors) on `reduction`, which base::subset() already restricted to just
    // the clean cells automatically. Three final plots: CellType, Lineage, Sample —
    // all on the new "umap_final" reduction.
    // =================================================================================
    script:

    def script_path    = "${workflow.projectDir}/modules"
    def LOG             = "SCRNA_ADD_CELLTYPE.log"
    def subset_rds_str  = (subset_rds_list instanceof List) \
        ? subset_rds_list.collect { it.toString() }.join(',') \
        : subset_rds_list.toString()

    """
    Rscript "${script_path}/SCRNA/SCRNA_add_celltype.R" \
        "${whole_rds}" \
        "${subset_labels}" \
        "${subset_rds_str}" \
        "${reduction}" \
        "." > "${LOG}" 2>&1

    EXIT_CODE=\$?
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "❌ ERROR: SCRNA_ADD_CELLTYPE failed" | tee -a "${LOG}" >&2
        exit \$EXIT_CODE
    fi

    echo "✅ SUCCESS: SCRNA_ADD_CELLTYPE completed" >> "${LOG}"
    """
}
