process SCRNA_PSEUDOBULK {

    tag "Pseudobulk aggregation: ${species} / ${group_name}"
    label 'process_heavy'
    label 'omics_r'

    publishDir = [
        [path: { "${params.proj_dir()}/${species}/07.Pseudobulk/${group_name}" }, mode: 'copy', pattern: "*.{csv,xlsx}"],
        [path: { "${params.proj_dir()}/${species}/Logs" },                        mode: 'copy', pattern: "*.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    // Fans out ONE task PER analyses entry — mirrors exactly how bulk RNAseq fans
    // contrasts out via Channel.fromList(params.deseq2.contrasts).combine(dds_ch).
    // Each task handles exactly one entry: no looping over multiple entries inside one
    // R run, no JSON needed to hand R a nested structure (selector_values is a single
    // flat comma-joined string, same convention as batch_vars elsewhere). Trade-off,
    // explicit: each task independently reloads seurat_rds — accepted in exchange for
    // incremental downstream pipelining (this group's RNA_CREATE_DDS can start the
    // moment THIS task finishes, not after every other entry's task also finishes) and
    // dropping the JSON/jsonlite dependency entirely.
    input:
    tuple val(species), path(seurat_rds), val(selector_col), val(selector_values), val(group_name)
    // seurat_rds      : annotated_seurat_final_clean.rds from SCRNA_ADD_CELLTYPE.
    // selector_col    : "CellType" or "Lineage" — which metadata column this entry selects on.
    // selector_values : comma-joined string, e.g. "Basal,Club" — the value(s) merged into
    //                   this one population. Split back into a vector in R via strsplit.
    // group_name      : precomputed ONCE in scrnaseq.nf via the same label_for()/
    //                   sanitize_celltype() used for subcluster labels — NOT re-derived
    //                   here, so there's only one sanitization implementation to keep
    //                   consistent with downstream directory naming, not two.
    val(min_cells_per_sample)
    // min_cells_per_sample : an individual Sample is dropped from this entry if fewer
    //                   than this many cells contributed to it. Applies to the COMBINED
    //                   cell count across every value in selector_values.
    //                   Passed from params.scrna.pseudobulk.min_cells_per_sample.

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    // All 3 optional: true — this task produces NO files and exits 0 (graceful skip,
    // not an error) if 0 cells matched the selector, or if fewer than MIN_TOTAL_SAMPLES
    // (hardcoded in SCRNA_pseudobulk.R) samples survived min_cells_per_sample filtering.
    // scrnaseq.nf's downstream RNA_CREATE_DDS simply never receives an item for a group
    // that produced nothing here — same "absence = skip" idiom used throughout this
    // pipeline (e.g. RNA_extract_deseq2_results.nf's per-contrast skip).
    tuple val(species), val(group_name), path("expr_mat.csv"), path("metadata.xlsx"), path("tx2gene_dummy.csv"), emit: pseudobulk_files, optional: true
    //tuple val(species), val(group_name), path("expr_mat.csv"),      emit: expr_mat,      optional: true
    //tuple val(species), val(group_name), path("metadata.xlsx"),     emit: metadata,      optional: true
    //tuple val(species), val(group_name), path("tx2gene_dummy.csv"), emit: tx2gene_dummy, optional: true
    path("SCRNA_PSEUDOBULK_${group_name}.log"),                     emit: log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def script_path = "${workflow.projectDir}/modules/SCRNA"
    def LOG          = "SCRNA_PSEUDOBULK_${group_name}.log"

    """
    Rscript "${script_path}/SCRNA_pseudobulk.R" \
        "${seurat_rds}" \
        "${selector_col}" \
        "${selector_values}" \
        "${group_name}" \
        "${min_cells_per_sample}" \
        "." > "${LOG}" 2>&1

    EXIT_CODE=\$?
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "❌ ERROR: SCRNA_PSEUDOBULK failed for ${species} / ${group_name}" | tee -a "${LOG}" >&2
        exit \$EXIT_CODE
    fi

    echo "✅ SUCCESS: SCRNA_PSEUDOBULK completed (or gracefully skipped) for ${species} / ${group_name}" >> "${LOG}"
    """
}
