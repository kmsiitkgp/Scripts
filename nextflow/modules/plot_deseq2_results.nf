// =============================================================================
// PROCESS: PLOT_DESEQ2_RESULTS
// =============================================================================
// Purpose: Generates pathway plots and optional expression heatmap for one contrast
//
// What it does:
//   - Reads pathway results Excel from ANALYZE_PATHWAYS
//   - Generates dot plots and bar plots summarizing enriched pathways
//   - Optionally overlays gene expression heatmap using expression matrix
//   - All plots saved under a contrast-named subdirectory
//
// Called once per contrast (fan-out via .combine() in rnaseq.nf)
// =============================================================================

process PLOT_DESEQ2_RESULTS {

    tag "Plotting DESeq2 results: ${species} / ${contrast}"
    label 'process_medium'                        // R with ggplot2; ~8GB RAM

    // get_safe sanitizes the contrast label for use as a directory name
    // e.g. "Treated vs Control" → "Treated_vs_Control"
    def get_safe = { c -> c.replaceAll(/[^a-zA-Z0-9-]/, '_') }

    publishDir { "${params.proj_dir()}/${species}/08.DESeq2" },    mode: 'copy',    pattern: "${get_safe(contrast)}/*.pdf"
    publishDir { "${params.proj_dir()}/${species}/07.Logs" },      mode: 'copy',    pattern: "*.error.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), val(contrast), path(res_rds), path(deg_xlsx), path(vst_nonblind_xlsx)
    // species       : "Human" or "Mouse"
    // contrast      : Contrast label (e.g., "condition_Treated_vs_Control")
    // pathway_xlsx : Pathway results from ANALYZE_PATHWAYS
    // vst_nonblind_xlsx : Expression matrix for heatmap overlay; passed as val (not path)
    //                 so R receives the string "NULL" when no matrix is provided
    path(metadata_xlsx)    // Sample metadata for heatmap column annotations
    val(heatmap_list)      // Heatmap params map from config (wrapped in List by .collect())

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("${get_safe(contrast)}/*.pdf"),    emit: deseq2_plots    // PDF plots for this contrast
    path("*.error.log"),                    emit: error_log       // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def safe_contrast = get_safe(contrast)
    def script_path   = "${workflow.projectDir}/modules"
    def LOG           = "${species}_${safe_contrast}_PLOT_DESEQ2_RESULTS.error.log"

    // .collect() wraps Maps in Lists: [ [padj_cutoff: 0.05, ...] ]
    // If we use "${gsea_list.padj_cutoff}", Bash receives "[0.05]", which 
    // causes R to throw "NAs introduced by coercion" and fail the analysis.
    // Unwrap first element to pass clean scalar values to R
    def heatmap = heatmap_list[0]

    // Handles all 3 col_annotations cases:
    //   null                   → "" (no column annotations)
    //   "batch"                → "batch"
    //   ["batch", "condition"] → "batch,condition"
    def col_annotation_str = (heatmap.col_annotations instanceof List) \
        ? heatmap.col_annotations.join(',') \
        : (heatmap.col_annotations ?: "")

    """
    # Arg 1: contrast          — contrast label string
    # Arg 2: deg_xlsx     — path to pathway results Excel
    # Arg 3: "."               — output directory (current work dir)
    # Arg 4: vst_nonblind_xlsx           — path to expression matrix Excel (or "NULL")
    # Arg 5: metadata_xlsx     — path to sample metadata Excel
    # Arg 6: col_annotation_str— comma-separated metadata columns for heatmap annotations
    # Arg 7: col_cluster_by    — column to use for hierarchical clustering
    # Arg 8: top_n             — number of top pathways to plot
    Rscript "${script_path}/PLOT_DESEQ2_RESULTS.R" \
        "${contrast}" \
        "${res_rds}" \
        "${deg_xlsx}" \
        "." \
        "${vst_nonblind_xlsx}" \
        "${metadata_xlsx}" \
        "${col_annotation_str}" \
        "${heatmap.col_cluster_by}" > "${LOG}" 2>&1 \
        || { echo "❌ ERROR: PLOT_DESEQ2_RESULTS.R failed for ${species} / ${contrast}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: DESeq2 plots generated for ${species} / ${contrast}" >> "${LOG}"
    """
}
