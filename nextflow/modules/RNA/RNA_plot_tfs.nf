// =============================================================================
// PROCESS: RNA_PLOT_TFS
// =============================================================================
// Purpose: Generates transcription-factor activity plots for one contrast.
//
// What it does:
//   - Reads TF activity/results from RNA_ANALYZE_TFS.
//   - Generates TF activity plots for the specified contrast.
//   - Optionally adds sample metadata annotations to the plots.
//   - Saves plots under a contrast-specific subdirectory.
//
// This process is called once per contrast:
//   - Bulk RNAseq: one task per configured contrast.
//   - scRNAseq pseudobulk: one task per (group_name, contrast) pair.
//
// group_name identifies the pseudobulk population and is empty for bulk RNAseq.
// It is used to distinguish pseudobulk output directories and task names.
// =============================================================================

process RNA_PLOT_TFS {

    tag "Plotting TFs: ${species}${group_name ? '|' + group_name : ''} | ${contrast}"

    label 'process_medium'
    label 'omics_r'

    // Sanitize the contrast name for use as a filesystem directory name.
    // Example: "Treated vs Control" → "Treated_vs_Control"
    def get_safe = { c -> c.replaceAll(/[^a-zA-Z0-9-]/, '_') }

    publishDir = [
        [path: { group_name
                    ? "${params.proj_dir()}/${species}/07.Pseudobulk/${group_name}/09.TFs"
                    : "${params.proj_dir()}/${species}/09.TFs" },
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
    tuple val(species), val(group_name), val(contrast), path(tf_xlsx), path(metadata_xlsx)

    // species   : "Human" or "Mouse"
    // group_name: pseudobulk group identifier; empty for bulk RNAseq
    // contrast  : contrast label, e.g. "condition_Treated_vs_Control"
    // tf_xlsx   : TF activity/results Excel file from RNA_ANALYZE_TFS
    // metadata_xlsx : sample metadata used for plot annotations

    val(heatmap_list)
    // Heatmap configuration map from config, wrapped in a List by .collect()

    val(decoupler_list)
    // Decoupler/TF plotting configuration map from config, wrapped in a List by .collect()

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    // Optional because a contrast may contain no TF results that can be plotted.
    // In that case the R script can complete successfully without producing PDFs.
    path("${get_safe(contrast)}/*.pdf"),
        emit: tf_plots,
        optional: true

    path("*.log"),
        emit: error_log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def safe_contrast = get_safe(contrast)
    def script_path   = "${workflow.projectDir}/modules/RNA"

    // Include group_name and contrast in pseudobulk log names so parallel
    // tasks do not overwrite each other's logs before publishDir.
    def LOG = group_name
        ? "${species}_${group_name}_${safe_contrast}_PLOT_TFS.log"
        : "${species}_${safe_contrast}_PLOT_TFS.log"

    // .collect() wraps configuration Maps in Lists:
    //   [ [top_n: 20, ...] ]
    //
    // Unwrap the first element so individual configuration values are passed
    // to R as scalar values rather than strings representing a List.
    def decoupler = decoupler_list[0]
    def heatmap   = heatmap_list[0]

    // Convert col_annotations to the comma-separated format expected by R.
    // Supported configuration forms:
    //   null                   → ""
    //   "batch"                → "batch"
    //   ["batch", "condition"] → "batch,condition"
    def col_annotation_str = (heatmap.col_annotations instanceof List) \
        ? heatmap.col_annotations.join(',') \
        : (heatmap.col_annotations ?: "")

    """
    # Arg 1: contrast               — contrast label string
    # Arg 2: tf_xlsx                — TF activity/results Excel file
    # Arg 3: metadata_xlsx          — sample metadata Excel file
    # Arg 4: col_annotation_str     — comma-separated metadata columns for annotations
    # Arg 5: col_cluster_by         — column used for hierarchical clustering
    # Arg 6: top_n                  — number of top TFs to plot
    # Arg 7: "."                   — output directory (current work directory)

    Rscript "${script_path}/RNA_plot_tfs.R" \
        "${contrast}" \
        "${tf_xlsx}" \
        "${metadata_xlsx}" \
        "${col_annotation_str}" \
        "${heatmap.col_cluster_by}" \
        "${decoupler.top_n}" \
        "." > "${LOG}" 2>&1

    EXIT_CODE=\$?
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "❌ ERROR: RNA_plot_tfs.R failed for ${species}${group_name ? '|' + group_name : ''} | ${contrast}" | tee -a "${LOG}" >&2
        exit \$EXIT_CODE
    fi

    echo "✅ SUCCESS: TF plots generated for ${species}${group_name ? '|' + group_name : ''} | ${contrast}" >> "${LOG}"
    """
}
