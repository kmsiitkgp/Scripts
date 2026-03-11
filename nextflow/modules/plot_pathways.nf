process PLOT_PATHWAYS {

    tag "Plotting Pathway analysis: ${species} ${contrast}"
    label 'process_medium'                      // STAR indexing requires 30-50GB RAM for human

    // We define the logic in a closure so it is evaluated when the block is called
    // Rscript creates the folder using safe_contrast internally
    def get_safe = { c -> c.replaceAll(/[^a-zA-Z0-9-]/, '_') }

    publishDir { "${params.proj_dir()}/${species}/09.Pathways" },    mode: 'copy',    pattern: "${get_safe(contrast)}"
    publishDir { "${params.proj_dir()}/${species}/07.Logs" },        mode: 'copy',    pattern: "*.error.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), val(contrast), path(pathways_xlsx), val(expr_mat_xlsx)
    // using val(expr_mat_xlsx) instead of path(expr_mat_xlsx) so R doesnt error
    // NULL is passed as a value to R instead of a character string 
    path(metadata_xlsx)
    val(heatmap_list)
    val(gsea_list)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("${get_safe(contrast)}"),    emit: pathway_plot_dir
    path("*.error.log"),              emit: error_log    // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================

    script:

    // Here we can use a local variable for the bash part
    def safe_contrast = get_safe(contrast)

    // This points to the modules folder relative to your project root
    def script_path = "${workflow.projectDir}/modules"
    def LOG = "${species}_${safe_contrast}_PLOT_PATHWAYS.error.log"

    // ========================================================================
    // BUG FIX: UNWRAP GSEA PARAMETERS
    // ------------------------------------------------------------------------
    // .collect() wraps the 'gsea' Map in a List: [ [padj_cutoff: 0.05, ...] ]
    // If we use "${gsea_list.padj_cutoff}", Bash receives "[0.05]", which 
    // causes R to throw "NAs introduced by coercion" and fail the analysis.
    // We unwrap the first element here to pass clean numbers to R.
    // ========================================================================
    def gsea    = gsea_list[0]
    def heatmap = heatmap_list[0]

    // --- THIS HANDLES ALL 3 CASES ---
    // 1. If null -> returns empty string ""
    // 2. If "batch" -> returns "batch"
    // 3. If ["batch", "condition"] -> returns "batch,condition"
    def col_annotation_str = (heatmap.col_annotations instanceof List) ? heatmap.col_annotations.join(',') : (heatmap.col_annotations ?: "")

    """
    # We pass '.' as the 4th arg so the Excel file is saved in the current folder
    Rscript ${script_path}/PLOT_PATHWAYS.R \
        "${contrast}" \
        "${pathways_xlsx}" \
        "." \
        "${expr_mat_xlsx}" \
        "${metadata_xlsx}" \
        "${col_annotation_str}" \
        "${heatmap.col_cluster_by}" \
        "${gsea.top_n}" > ${LOG} 2>&1
    """
}