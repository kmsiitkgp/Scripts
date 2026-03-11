process ANALYZE_PATHWAYS {

    tag "Performing Pathway analysis: ${species} ${contrast}"
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
    tuple val(species), val(contrast), path(deg_xlsx), path(gmt_dir)
    val(gsea_list)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), val(contrast), path("${get_safe(contrast)}"),    emit: pathway_dir
    path("*.error.log"),              emit: error_log    // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================

    script:

    // Here we can use a local variable for the bash part
    def safe_contrast = get_safe(contrast)

    // This points to the modules folder relative to your project root
    def script_path = "${workflow.projectDir}/modules"
    def LOG = "${species}_${safe_contrast}_ANALYZE_PATHWAYS.error.log"

    // ========================================================================
    // BUG FIX: UNWRAP GSEA PARAMETERS
    // ------------------------------------------------------------------------
    // .collect() wraps the 'gsea' Map in a List: [ [padj_cutoff: 0.05, ...] ]
    // If we use "${gsea_list.padj_cutoff}", Bash receives "[0.05]", which 
    // causes R to throw "NAs introduced by coercion" and fail the analysis.
    // We unwrap the first element here to pass clean numbers to R.
    // ========================================================================
    def gsea    = gsea_list[0]

    """
    # We pass '.' as the 4th arg so the Excel file is saved in the current folder
    Rscript ${script_path}/ANALYZE_PATHWAYS.R \
        "${contrast}" \
        "${deg_xlsx}" \
        "${gsea.only_deg}" \
        "${gmt_dir}" \
        "." \
        "${gsea.padj_cutoff}" \
        "${gsea.minsize}" \
        "${gsea.maxsize}" > ${LOG} 2>&1
    """
}