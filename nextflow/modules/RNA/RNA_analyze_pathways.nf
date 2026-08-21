// =============================================================================
// PROCESS: RNA_ANALYZE_PATHWAYS
// =============================================================================
// Purpose: Performs pathway analysis for one DESeq2 contrast.
//
// What it does:
//   - Runs GSEA using the configured GMT gene sets.
//   - Optionally restricts the GSEA input to significant DEGs.
//   - Applies the configured pathway-analysis parameters.
//   - Saves pathway results under a contrast-specific subdirectory.
//
// This process is called once per contrast:
//   - Bulk RNAseq: one task per configured contrast.
//   - scRNAseq pseudobulk: one task per (group_name, contrast) pair.
//
// group_name identifies the pseudobulk population and is empty for bulk RNAseq.
// =============================================================================

process RNA_ANALYZE_PATHWAYS {

    tag "Pathway analysis: ${species}${group_name ? '|' + group_name : ''} | ${contrast}"

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
         pattern: "*/Pathways.xlsx"],

        [path: { "${params.proj_dir()}/${species}/Logs" },
         mode: 'copy',
         pattern: "*.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), val(group_name), val(contrast),
          path(degs_xlsx), path(vst_nonblind_xlsx), path(metadata_xlsx)

    // species           : "Human" or "Mouse"
    // group_name        : pseudobulk group identifier; empty for bulk RNAseq
    // contrast          : DESeq2 contrast string, e.g. "GHRKO_Late-WT_Late"
    // degs_xlsx         : DEG results Excel from RNA_EXTRACT_DESEQ2_RESULTS
    // vst_nonblind_xlsx : VST non-blind expression matrix for the contrast; used for expression heatmaps
    // metadata_xlsx     : sample metadata for this group; used for sample annotations

    path(gmt_dir)
    // gmt_dir : directory containing species-specific GMT gene-set files

    val(gsea_list)
    // gsea_list : GSEA configuration map from config

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:

    tuple val(species), val(group_name), val(contrast),
          path("${get_safe(contrast)}/Pathways.xlsx"),
          path(vst_nonblind_xlsx),
          path(metadata_xlsx),
          emit: pathways

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
        ? "${species}_${group_name}_${safe_contrast}_ANALYZE_PATHWAYS.log"
        : "${species}_${safe_contrast}_ANALYZE_PATHWAYS.log"

    // .collect() wraps the GSEA configuration map in a List.
    // Extract the single configuration map before accessing its parameters.
    def gsea = gsea_list[0]

    """
    # Arg 1: contrast       — DESeq2 contrast string
    # Arg 2: degs_xlsx      — DEG results Excel
    # Arg 3: gmt_dir        — directory containing GMT gene-set files
    # Arg 4: only_deg       — whether to restrict GSEA to significant DEGs
    # Arg 5: padj_cutoff    — adjusted p-value threshold
    # Arg 6: minsize        — minimum gene-set size
    # Arg 7: maxsize        — maximum gene-set size
    # Arg 8: "."            — output directory

    Rscript "${script_path}/RNA_analyze_pathways.R" \
        "${contrast}" \
        "${degs_xlsx}" \
        "${gmt_dir}" \
        "${gsea.only_deg}" \
        "${gsea.padj_cutoff}" \
        "${gsea.minsize}" \
        "${gsea.maxsize}" \
        "." > "${LOG}" 2>&1

    EXIT_CODE=\$?
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "❌ ERROR: RNA_analyze_pathways.R failed for ${species}${group_name ? '|' + group_name : ''} | ${contrast}" | tee -a "${LOG}" >&2
        exit \$EXIT_CODE
    fi

    echo "✅ SUCCESS: Pathway analysis completed for ${species}${group_name ? '|' + group_name : ''} | ${contrast}" >> "${LOG}"
    """
}
