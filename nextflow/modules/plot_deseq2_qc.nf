// =============================================================================
// PROCESS: PLOT_DESEQ2_QC
// =============================================================================
// Purpose: Generates pathway plots and optional expression heatmap for one contrast
//
// What it does:
//   - Reads pathway results Excel from ANALYZE_PATHWAYS
//   - Generates dot plots and bar plots summarizing enriched pathways
//   - Optionally overlays gene expression heatmap using expression matrix
//   - All plots saved under a contrast-named subdirectory
//
// Called once per dataset (fan-out via .combine() in rnaseq.nf)
// =============================================================================

process PLOT_DESEQ2_QC {

    tag "Plotting DESeq2 QC: ${species}"
    label 'process_medium'                        // R with ggplot2; ~8GB RAM

    publishDir { "${params.proj_dir()}/${species}/08.DESeq2" },    mode: 'copy',    pattern: "*.pdf"
    publishDir { "${params.proj_dir()}/${species}/07.Logs" },      mode: 'copy',    pattern: "*.error.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), path(dds_rds), path(vst_blind_xlsx)
    // species       : "Human" or "Mouse"
    // pathway_xlsx : Pathway results from ANALYZE_PATHWAYS
    // vst_nonblind_xlsx : Expression matrix for heatmap overlay; passed as val (not path)
    //                 so R receives the string "NULL" when no matrix is provided
    path(metadata_xlsx)    // Sample metadata for heatmap column annotations

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("*.pdf"),        emit: qc_plots    // PDF plots for QC
    path("*.error.log"),  emit: error_log   // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def script_path   = "${workflow.projectDir}/modules"
    def LOG           = "${species}_PLOT_DESEQ2_QC.error.log"

    """
    # Arg 1: contrast          — contrast label string
    # Arg 2: deg_xlsx     — path to pathway results Excel
    # Arg 3: "."               — output directory (current work dir)
    # Arg 4: vst_nonblind_xlsx           — path to expression matrix Excel (or "NULL")
    # Arg 5: metadata_xlsx     — path to sample metadata Excel
    # Arg 6: col_annotation_str— comma-separated metadata columns for heatmap annotations
    # Arg 7: col_cluster_by    — column to use for hierarchical clustering
    # Arg 8: top_n             — number of top pathways to plot
    Rscript "${script_path}/PLOT_DESEQ2_QC.R" \
        "${dds_rds}" \
        "${vst_blind_xlsx}" \
        "${metadata_xlsx}" \
        "." > "${LOG}" 2>&1 \
        || { echo "❌ ERROR: PLOT_DESEQ2_QC.R failed for ${species}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: DESeq2 QC plots generated for ${species}" >> "${LOG}"
    """
}
