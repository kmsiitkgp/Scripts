// =============================================================================
// PROCESS: RNA_PLOT_PATHWAYS
// =============================================================================
// Purpose: Generates pathway enrichment plots for one contrast.
//
// What it does:
//   - Reads pathway results from RNA_ANALYZE_PATHWAYS.
//   - Generates pathway plots using the configured top_n value.
//   - Optionally generates an expression heatmap when a VST non-blind
//     expression matrix and sample metadata are provided.
//   - Saves plots under a contrast-specific subdirectory.
//
// This process is called once per contrast:
//   - Bulk RNAseq: one task per configured contrast.
//   - scRNAseq pseudobulk: one task per (group_name, contrast) pair.
//
// group_name identifies the pseudobulk population and is empty for bulk RNAseq.
// =============================================================================

process RNA_PLOT_PATHWAYS {

    tag "Plotting pathways: ${species}${group_name ? '|' + group_name : ''} | ${contrast}"

    label 'process_medium'
    label 'omics_r'

    // Sanitize the contrast name for use as a filesystem directory name.
    // Example: "Treated vs Control" → "Treated_vs_Control"
    def get_safe = { c -> c.replaceAll(/[^a-zA-Z0-9-]/, '_') }

    publishDir = [
        [path: { group_name
                    ? "${params.proj_dir()}/${species}/07.Pseudobulk/${group_name}/08.Pathways"
                    : "${params.proj_dir()}/${species}/08.Pathways" },
         mode: 'copy',
         pattern: "*/*.pdf",
         optional: true],

        [path: { "${params.proj_dir()}/${species}/Logs" },
         mode: 'copy',
         pattern: "*.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), val(group_name), val(contrast),
          path(pathway_xlsx), val(vst_nonblind_xlsx), path(metadata_xlsx)

    // species            : "Human" or "Mouse"
    // group_name         : pseudobulk group identifier; empty for bulk RNAseq
    // contrast           : DESeq2 contrast string, e.g. "GHRKO_Late-WT_Late"
    // pathway_xlsx       : pathway results Excel from RNA_ANALYZE_PATHWAYS
    // vst_nonblind_xlsx  : VST non-blind expression matrix; passed as val so
    //                      "NULL" can be supplied when no matrix is available
    // metadata_xlsx : sample metadata used for heatmap sample annotations

    val(heatmap_list)
    // heatmap_list : heatmap configuration map from config

    val(gsea_list)
    // gsea_list : GSEA configuration map from config

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:

    // A contrast with no significant pathways may produce no PDF files.
    // Treat this as a valid outcome rather than a process failure.
    path("${get_safe(contrast)}/*.pdf"),
         emit: pathway_plots,
         optional: true

    path("*.log"), emit: error_log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def safe_contrast = get_safe(contrast)
    def script_path   = "${workflow.projectDir}/modules/RNA"

    // Include group_name in pseudobulk log names so parallel groups do not
    // overwrite each other's logs.
    def LOG = group_name
        ? "${species}_${group_name}_${safe_contrast}_PLOT_PATHWAYS.log"
        : "${species}_${safe_contrast}_PLOT_PATHWAYS.log"

    // .collect() wraps each configuration map in a List.
    // Extract the single configuration map before accessing its parameters.
    def gsea    = gsea_list[0]
    def heatmap = heatmap_list[0]

    // Convert col_annotations to the comma-separated format expected by R.
    // Supported configuration forms:
    //   null                   → ""
    //   "batch"                → "batch"
    //   ["batch", "condition"] → "batch,condition"
    def col_annotation_str = (heatmap.col_annotations instanceof List) \
        ? heatmap.col_annotations.join(',') \
        : (heatmap.col_annotations ?: "")

    """
    # Arg 1: contrast               — DESeq2 contrast string
    # Arg 2: pathway_xlsx            — pathway results Excel
    # Arg 3: vst_nonblind_xlsx       — VST non-blind expression matrix or "NULL"
    # Arg 4: metadata_xlsx           — sample metadata Excel
    # Arg 5: col_annotation_str      — comma-separated metadata columns for heatmap
    # Arg 6: col_cluster_by          — metadata column used for sample clustering
    # Arg 7: top_n                   — number of top pathways to plot
    # Arg 8: "."                     — output directory

    Rscript "${script_path}/RNA_plot_pathways.R" \
        "${contrast}" \
        "${pathway_xlsx}" \
        "${vst_nonblind_xlsx}" \
        "${metadata_xlsx}" \
        "${col_annotation_str}" \
        "${heatmap.col_cluster_by}" \
        "${gsea.top_n}" \
        "." > "${LOG}" 2>&1

    EXIT_CODE=\$?
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "❌ ERROR: RNA_plot_pathways.R failed for ${species}${group_name ? '|' + group_name : ''} | ${contrast}" | tee -a "${LOG}" >&2
        exit \$EXIT_CODE
    fi

    echo "✅ SUCCESS: Pathway plots generated for ${species}${group_name ? '|' + group_name : ''} | ${contrast}" >> "${LOG}"
    """
}
