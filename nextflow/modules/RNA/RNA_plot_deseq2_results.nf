// =============================================================================
// PROCESS: RNA_PLOT_DESEQ2_RESULTS
// =============================================================================
// Purpose: Generates DESeq2 result plots for one contrast.
//
// What it does:
//   - Generates an MA plot from the DESeq2 results object.
//   - Generates a volcano plot from the DEG results.
//   - Optionally generates a DEG heatmap using the VST non-blind expression
//     matrix and sample metadata.
//   - Saves all plots under a contrast-specific subdirectory.
//
// This process is called once per contrast:
//   - Bulk RNAseq: one task per configured contrast.
//   - scRNAseq pseudobulk: one task per (group_name, contrast) pair.
//
// group_name identifies the pseudobulk population and is empty for bulk RNAseq.
// =============================================================================

process RNA_PLOT_DESEQ2_RESULTS {

    tag "Plotting DESeq2 results: ${species}${group_name ? '|' + group_name : ''} | ${contrast}"

    label 'process_medium'
    label 'omics_r'

    // Sanitize the contrast name for use as a filesystem directory name.
    // Example: "Treated vs Control" → "Treated_vs_Control"
    def get_safe = { c -> c.replaceAll(/[^a-zA-Z0-9-]/, '_') }

    publishDir = [
        [path: { group_name
                    ? "${params.proj_dir()}/${species}/07.Pseudobulk/${group_name}/07.DESeq2"
                    : "${params.proj_dir()}/${species}/07.DESeq2" },
         mode: 'copy',
         pattern: "*/*.pdf"],

        [path: { "${params.proj_dir()}/${species}/Logs" },
         mode: 'copy',
         pattern: "*.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), val(group_name), val(contrast),
          path(res_rds), path(degs_xlsx), path(vst_nonblind_xlsx), path(metadata_xlsx)
    // species            : "Human" or "Mouse"
    // group_name         : pseudobulk group identifier; empty for bulk RNAseq
    // contrast           : DESeq2 contrast string, e.g. "GHRKO_Late-WT_Late"
    // res_rds            : DESeq2 results RDS from RNA_EXTRACT_DESEQ2_RESULTS
    // degs_xlsx          : DEG results Excel from RNA_EXTRACT_DESEQ2_RESULTS
    // vst_nonblind_xlsx  : VST non-blind expression matrix Excel; used for DEG heatmap
    // metadata_xlsx      : sample metadata used for heatmap sample annotations

    val(heatmap_list)
    // heatmap_list : heatmap configuration map from config

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("${get_safe(contrast)}/*.pdf"), emit: deseq2_plots
    path("*.log"),                       emit: error_log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def safe_contrast = get_safe(contrast)
    def script_path   = "${workflow.projectDir}/modules/RNA"

    // Include group_name in pseudobulk log names so parallel groups do not
    // overwrite each other's logs.
    def LOG = group_name
        ? "${species}_${group_name}_${safe_contrast}_PLOT_DESEQ2_RESULTS.log"
        : "${species}_${safe_contrast}_PLOT_DESEQ2_RESULTS.log"

    // .collect() wraps the heatmap configuration map in a List.
    // Extract the single configuration map before accessing its parameters.
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
    # Arg 2: res_rds                — DESeq2 results RDS
    # Arg 3: degs_xlsx              — DEG results Excel
    # Arg 4: vst_nonblind_xlsx      — VST non-blind expression matrix Excel
    # Arg 5: metadata_xlsx          — sample metadata Excel
    # Arg 6: col_annotation_str     — comma-separated metadata columns for heatmap
    # Arg 7: col_cluster_by         — metadata column used for sample clustering
    # Arg 8: "."                    — output directory

    Rscript "${script_path}/RNA_plot_deseq2_results.R" \
        "${contrast}" \
        "${res_rds}" \
        "${degs_xlsx}" \
        "${vst_nonblind_xlsx}" \
        "${metadata_xlsx}" \
        "${col_annotation_str}" \
        "${heatmap.col_cluster_by}" \
        "." > "${LOG}" 2>&1

    EXIT_CODE=\$?
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "❌ ERROR: RNA_plot_deseq2_results.R failed for ${species}${group_name ? '|' + group_name : ''} | ${contrast}" | tee -a "${LOG}" >&2
        exit \$EXIT_CODE
    fi

    echo "✅ SUCCESS: DESeq2 plots generated for ${species}${group_name ? '|' + group_name : ''} | ${contrast}" >> "${LOG}"
    """
}
