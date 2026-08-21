// =============================================================================
// PROCESS: RNA_ANALYZE_TFS
// =============================================================================
// Purpose: Performs transcription-factor activity analysis for one DESeq2 contrast.
//
// What it does:
//   - Uses DESeq2 differential-expression results and log-normalized expression
//     values as input for TF activity analysis.
//   - Runs the configured decoupleR methods.
//   - Uses the species argument to select the appropriate regulon database.
//   - Saves TF activity results as an Excel file under a contrast-specific
//     subdirectory.
//
// This process is called once per contrast:
//   - Bulk RNAseq: one task per configured contrast.
//   - scRNAseq pseudobulk: one task per (group_name, contrast) pair.
//
// group_name identifies the pseudobulk population and is empty for bulk RNAseq.
// It is used to distinguish pseudobulk output directories and task names.
//
// species is a biological argument and must remain the actual species
// ("Human" or "Mouse") because it determines the appropriate TF regulon database.
// =============================================================================

process RNA_ANALYZE_TFS {

    tag "TF analysis: ${species}${group_name ? '|' + group_name : ''} | ${contrast}"

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
         pattern: "*/TF_and_Pathway_Activity.xlsx"],

        [path: { "${params.proj_dir()}/${species}/Logs" },
         mode: 'copy',
         pattern: "*.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), val(group_name), val(contrast),
          path(deg_xlsx), path(log_norm_counts_xlsx), path(metadata_xlsx)

    // species               : "Human" or "Mouse"; selects the appropriate TF regulon database
    // group_name            : pseudobulk group identifier; empty for bulk RNAseq
    // contrast              : contrast label, e.g. "condition_Treated_vs_Control"
    // deg_xlsx              : DEGs.xlsx from RNA_EXTRACT_DESEQ2_RESULTS
    // log_norm_counts_xlsx  : log-normalized expression matrix from RNA_CREATE_DDS
    // metadata_xlsx         : sample metadata used by the TF analysis

    val(decoupler_list)
    // decoupleR configuration map from config, wrapped in a List by .collect()

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), val(group_name), val(contrast),
          path("${get_safe(contrast)}/TF_and_Pathway_Activity.xlsx"),
          path(metadata_xlsx),
          emit: tfs

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
        ? "${species}_${group_name}_${safe_contrast}_ANALYZE_TFS.log"
        : "${species}_${safe_contrast}_ANALYZE_TFS.log"

    // .collect() wraps configuration Maps in Lists:
    //   [ [methods: [...], minsize: 5, ...] ]
    //
    // Unwrap the first element so individual configuration values are passed
    // to R as scalar values rather than strings representing a List.
    def decoupler = decoupler_list[0]

    // Convert the configured TF activity methods to the comma-separated format
    // expected by the R script.
    // Supported configuration forms:
    //   null                   → default methods
    //   "ulm"                  → "ulm"
    //   ["ulm", "mlm", "viper"] → "ulm,mlm,viper"
    def methods_str = (decoupler.methods instanceof List) \
        ? decoupler.methods.join(',') \
        : (decoupler.methods ?: "ulm,mlm,viper")

    """
    # Arg 1: contrast               — contrast label string
    # Arg 2: deg_xlsx               — DESeq2 differential-expression results Excel
    # Arg 3: log_norm_counts_xlsx   — log-normalized expression matrix Excel
    # Arg 4: metadata_xlsx          — sample metadata Excel
    # Arg 5: species                — biological species; selects the TF regulon database
    # Arg 6: methods_str            — comma-separated decoupleR methods
    # Arg 7: minsize                — minimum regulon size
    # Arg 8: "."                    — output directory (current work directory)

    Rscript "${script_path}/RNA_analyze_tfs.R" \
        "${contrast}" \
        "${deg_xlsx}" \
        "${log_norm_counts_xlsx}" \
        "${metadata_xlsx}" \
        "${species}" \
        "${methods_str}" \
        "${decoupler.minsize}" \
        "." > "${LOG}" 2>&1

    EXIT_CODE=\$?
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "❌ ERROR: RNA_analyze_tfs.R failed for ${species}${group_name ? '|' + group_name : ''} | ${contrast}" | tee -a "${LOG}" >&2
        exit \$EXIT_CODE
    fi

    echo "✅ SUCCESS: TF analysis completed for ${species}${group_name ? '|' + group_name : ''} | ${contrast}" >> "${LOG}"
    """
}
