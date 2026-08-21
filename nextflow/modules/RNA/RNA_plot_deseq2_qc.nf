// =============================================================================
// PROCESS: RNA_PLOT_DESEQ2_QC
// =============================================================================
// Purpose: Generates DESeq2 quality-control plots for one dataset/group.
//
// What it does:
//   - Reads the DESeq2 dataset (DDS) for the dataset/group.
//   - Generates a DESeq2 dispersion plot from the DDS.
//   - Reads the VST-blind expression matrix and sample metadata.
//   - Generates a PCA plot using the VST-blind expression matrix and metadata.
//   - Saves all QC plots as PDF files.
//
// Called once per dataset/group.
// For bulk RNAseq, group_name is empty and outputs are published under
// {species}/07.DESeq2/.
// For scRNAseq pseudobulk, group_name identifies the pseudobulk population and
// outputs are published under {species}/07.Pseudobulk/{group_name}/07.DESeq2/.
// =============================================================================

process RNA_PLOT_DESEQ2_QC {

    tag "Plotting DESeq2 QC: ${species}${group_name ? '|' + group_name : ''}"

    label 'process_medium'
    label 'omics_r'

    publishDir = [
        [path: { group_name
                    ? "${params.proj_dir()}/${species}/07.Pseudobulk/${group_name}/07.DESeq2"
                    : "${params.proj_dir()}/${species}/07.DESeq2" },
         mode: 'copy',
         pattern: "*.pdf"],

        [path: { "${params.proj_dir()}/${species}/Logs" },
         mode: 'copy',
         pattern: "*.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), val(group_name), path(dds_rds), path(vst_blind_xlsx), path(metadata_xlsx)

    // species        : "Human" or "Mouse"
    // group_name     : pseudobulk group name; "" for bulk RNAseq
    // dds_rds        : DESeq2_dds.rds from RNA_CREATE_DDS
    // vst_blind_xlsx : VST_Blind*.xlsx from RNA_CREATE_DDS
    // metadata_xlsx  : sample metadata used for heatmap column annotations

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("*.pdf"), emit: qc_plots
    path("*.log"), emit: error_log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def script_path = "${workflow.projectDir}/modules/RNA"

    // Include group_name in pseudobulk log names so parallel groups do not
    // overwrite each other's logs.
    def LOG = group_name
        ? "${species}_${group_name}_PLOT_DESEQ2_QC.log"
        : "${species}_PLOT_DESEQ2_QC.log"

    """
    # Arg 1: dds_rds                — path to DESeq2_dds.rds
    # Arg 2: vst_nonblind_xlsx      — path to expression matrix Excel (or "NULL")
    # Arg 3: metadata_xlsx          — path to sample metadata Excel
    # Arg 4: "."                    — output directory (current work dir)

    Rscript "${script_path}/RNA_plot_deseq2_qc.R" \
        "${dds_rds}" \
        "${vst_blind_xlsx}" \
        "${metadata_xlsx}" \
        "." > "${LOG}" 2>&1

    EXIT_CODE=\$?
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "❌ ERROR: RNA_plot_deseq2_qc.R failed for ${species}${group_name ? '|' + group_name : ''}" | tee -a "${LOG}" >&2
        exit \$EXIT_CODE
    fi

    echo "✅ SUCCESS: DESeq2 QC plots generated for ${species}${group_name ? '|' + group_name : ''}" >> "${LOG}"
    """
}
