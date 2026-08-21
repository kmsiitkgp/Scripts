// =============================================================================
// PROCESS: RNA_EXTRACT_DESEQ2_RESULTS
// =============================================================================
// Purpose: Extracts DESeq2 differential expression results for one contrast.
//
// What it does:
//   - Runs DESeq2 results() for the specified contrast.
//   - Applies LFC shrinkage.
//   - Filters DEGs using the configured LFC and adjusted-p-value cutoffs.
//   - Handles optional batch correction variables.
//   - Saves results under a contrast-specific subdirectory.
//
// This process is called once per contrast:
//   - Bulk RNAseq: one task per configured contrast.
//   - scRNAseq pseudobulk: one task per (group_name, contrast) pair.
//
// group_name identifies the pseudobulk population and is empty for bulk RNAseq.
// It is used only to distinguish pseudobulk output directories and task names.
//
// The R script can gracefully skip a contrast when either referenced group
// does not contain enough samples according to min_samples_per_group.
// =============================================================================

process RNA_EXTRACT_DESEQ2_RESULTS {

    tag "DESeq2 results: ${species}${group_name ? '|' + group_name : ''} | ${contrast}"

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
         pattern: "*/*.{rds,xlsx}"],

        [path: { "${params.proj_dir()}/${species}/Logs" },
         mode: 'copy',
         pattern: "*.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), val(group_name), path(dds_rds), path(tx2gene_csv), 
          path(log_norm_counts_xlsx), path(metadata_xlsx), val(contrast)

    // species              : "Human" or "Mouse"
    // group_name           : pseudobulk group identifier; empty for bulk RNAseq
    // dds_rds              : DESeq2_dds.rds produced by RNA_CREATE_DDS
    // tx2gene_csv          : gene/transcript annotation mapping used for result annotation
    // log_norm_counts_xlsx : 
    // metadata_xlsx        : sample metadata for this group
    // contrast             : DESeq2 contrast string, e.g. "GHRKO_Late-WT_Late"

    val(deseq2)
    // DESeq2 configuration map containing:
    //   - batch_vars
    //   - lfc_cutoff
    //   - padj_cutoff
    //   - min_samples_per_group

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:

    // These outputs are optional because the R script can intentionally skip a
    // contrast when the referenced groups do not meet min_samples_per_group.
    // In that case the process still exits successfully but produces no result files.

    tuple val(species), val(group_name), val(contrast),
          path("${get_safe(contrast)}/res_parser.rds"),
          path("${get_safe(contrast)}/DEGs.xlsx"),
          path("${get_safe(contrast)}/VST_NonBlind_Counts_Heatmaps.xlsx"),
          path(log_norm_counts_xlsx),
          path(metadata_xlsx),
          emit: results,
          optional: true

    path("${get_safe(contrast)}/res_standard.rds"),
         emit: res_standard,
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
        ? "${species}_${group_name}_${safe_contrast}_EXTRACT_DESEQ2_RESULTS.log"
        : "${species}_${safe_contrast}_EXTRACT_DESEQ2_RESULTS.log"

    // Convert batch_vars to the comma-separated format expected by the R script.
    // Supported configuration forms:
    //   null                   → ""
    //   "batch"                → "batch"
    //   ["batch", "condition"] → "batch,condition"
    def batch_str = (deseq2.batch_vars instanceof List) \
        ? deseq2.batch_vars.join(',') \
        : (deseq2.batch_vars ?: "")

    """
    # Arg 1: contrast               — contrast string for DESeq2 results()
    # Arg 2: dds_rds                — path to DESeq2_dds.rds
    # Arg 3: tx2gene_csv            — path to tx2gene CSV
    # Arg 4: batch_str              — comma-separated batch variables (or "" if none)
    # Arg 5: lfc_cutoff             — log2 fold change threshold for DEG filtering
    # Arg 6: padj_cutoff            — adjusted p-value threshold for DEG filtering
    # Arg 7: min_samples_per_group  — per-contrast sample-count floor
    # Arg 8: "."                    — output directory (current work dir; R creates safe_contrast subdir internally)

    Rscript "${script_path}/RNA_extract_deseq2_results.R" \
        "${contrast}" \
        "${dds_rds}" \
        "${tx2gene_csv}" \
        "${batch_str}" \
        "${deseq2.lfc_cutoff}" \
        "${deseq2.padj_cutoff}" \
        "${deseq2.min_samples_per_group ?: 2}" \
        "." > "${LOG}" 2>&1

    EXIT_CODE=\$?
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "❌ ERROR: RNA_extract_deseq2_results.R failed for ${species}${group_name ? '|' + group_name : ''} | ${contrast}" | tee -a "${LOG}" >&2
        exit \$EXIT_CODE
    fi

    echo "✅ SUCCESS: DESeq2 results extracted (or gracefully skipped) for ${species}${group_name ? '|' + group_name : ''} | ${contrast}" >> "${LOG}"
    """
}
