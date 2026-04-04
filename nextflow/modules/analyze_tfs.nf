// =============================================================================
// PROCESS: ANALYZE_PATHWAYS
// =============================================================================
// Purpose: Runs GSEA and ORA pathway analysis on DESeq2 results for one contrast
//
// What it does:
//   - Runs gene set enrichment analysis (GSEA) using fgsea
//   - Runs over-representation analysis (ORA) using clusterProfiler
//   - Filters to DEGs only if configured (gsea.only_deg)
//   - Saves pathway results as Excel under a contrast-named subdirectory
//
// Called once per contrast (fan-out via .combine() in rnaseq.nf)
// =============================================================================

process ANALYZE_TFS {

    tag "Pathway analysis: ${species} / ${contrast}"
    label 'process_medium'                        // R with fgsea/clusterProfiler; ~8GB RAM

    // get_safe sanitizes the contrast label for use as a directory name
    // e.g. "Treated vs Control" → "Treated_vs_Control"
    def get_safe = { c -> c.replaceAll(/[^a-zA-Z0-9-]/, '_') }

    publishDir { "${params.proj_dir()}/${species}/10.TFs" },    mode: 'copy',    pattern: "${get_safe(contrast)}/TF*.xlsx"
    publishDir { "${params.proj_dir()}/${species}/07.Logs" },   mode: 'copy',    pattern: "*.error.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), val(contrast), path(deg_xlsx), path(log_norm_counts_xlsx)
    // species   : "Human" or "Mouse"
    // contrast  : Contrast label (e.g., "condition_Treated_vs_Control")
    // deg_xlsx  : DEGs.xlsx from EXTRACT_DESEQ2_RESULTS
    // gmt_dir   : Directory of GMT gene set files for this species
    path(metadata_xlsx)
    val(decoupler_list)    // GSEA params map from config (wrapped in List by .collect())

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), val(contrast), path("${get_safe(contrast)}/TF*.xlsx"),    emit: tf_results    // Pathway results directory
    path("*.error.log"),                                                  emit: error_log      // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def safe_contrast = get_safe(contrast)
    def script_path   = "${workflow.projectDir}/modules"
    def LOG           = "${species}_${safe_contrast}_ANALYZE_TFS.error.log"

    // .collect() wraps the gsea Map in a List: [ [padj_cutoff: 0.05, ...] ]
    // If we use "${gsea_list.padj_cutoff}", Bash receives "[0.05]", which 
    // causes R to throw "NAs introduced by coercion" and fail the analysis.
    // Unwrap first element to pass clean scalar values to R
    def decoupler = decoupler_list[0]

    """
    # Arg 1: contrast      — contrast label string
    # Arg 2: deg_xlsx          — path to DEGs.xlsx from EXTRACT_DESEQ2_RESULTS
    # Arg 3: only_deg      — boolean; restrict GSEA ranking to DEGs only
    # Arg 4: gmt_dir       — directory containing GMT gene set files
    # Arg 5: "."           — output directory (current work dir)
    # Arg 6: padj_cutoff   — adjusted p-value threshold for pathway significance
    # Arg 7: minsize       — minimum gene set size for GSEA
    # Arg 8: maxsize       — maximum gene set size for GSEA
    Rscript "${script_path}/ANALYZE_TFS.R" \
        "${contrast}" \
        "${deg_xlsx}" \
        "${log_norm_counts_xlsx}" \
        "${metadata_xlsx}" \
        "${species}" \
        "." \
        "${decoupler.methods}" \
        "${decoupler.minsize}" > "${LOG}" 2>&1 \
        || { echo "❌ ERROR: ANALYZE_TFS.R failed for ${species} / ${contrast}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: Pathway analysis completed for ${species} / ${contrast}" >> "${LOG}"
    """
}
