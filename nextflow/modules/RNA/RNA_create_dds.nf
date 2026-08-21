// =============================================================================
// PROCESS: RNA_CREATE_DDS
// =============================================================================
// Purpose: Creates a DESeq2 dataset (DESeqDataSet) from count data and
// sample-level metadata.
//
// What it does:
//   - Reads either:
//       1. a Salmon tximport object for bulk RNAseq, or
//       2. a raw pseudobulk count matrix for scRNAseq pseudobulk analysis.
//   - Reads sample metadata and the tx2gene annotation mapping.
//   - Constructs the DESeqDataSet using the configured design formula.
//   - Performs DESeq2 normalization and dispersion estimation.
//   - Saves the DESeqDataSet as an RDS file.
//   - Exports normalized, log-normalized, and VST-blind count matrices.
//
// This process is shared between bulk RNAseq and scRNAseq pseudobulk DE.
// group_name is empty for bulk RNAseq and identifies the pseudobulk
// population for scRNAseq pseudobulk analysis.
// =============================================================================

process RNA_CREATE_DDS {

    tag "Creating DESeq2 dataset: ${species}${group_name ? '|' + group_name : ''}"

    label 'process_medium'
    label 'omics_r'

    publishDir = [
        [path: { group_name
                    ? "${params.proj_dir()}/${species}/07.Pseudobulk/${group_name}"
                    : "${params.proj_dir()}/${species}/07.DESeq2" },
         mode: 'copy',
         pattern: "*.{rds,xlsx}"],

        [path: { "${params.proj_dir()}/${species}/Logs" },
         mode: 'copy',
         pattern: "*.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), val(group_name), path(count_data), path(metadata_xlsx), path(tx2gene_csv), val(design)

    // species       : "Human" or "Mouse"
    // group_name    : pseudobulk group identifier; empty for bulk RNAseq
    // count_data    : Salmon tximport RDS for bulk RNAseq, or raw pseudobulk
    //                 count matrix for scRNAseq pseudobulk
    // metadata_xlsx : sample-level metadata Excel file; columns must match
    //                 variables referenced by the DESeq2 design formula
    // tx2gene_csv   : transcript/gene annotation mapping used for annotation
    //                 and duplicate-gene handling

    // design        : DESeq2 design formula string, e.g. "~condition"
    //                 or "~batch+condition"

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:

    tuple val(species), val(group_name),
      path("DESeq2_dds.rds"),
      path(tx2gene_csv),
      path("Norm_Counts*.xlsx"),
      path("Log_Norm*.xlsx"),
      path("VST_Blind*.xlsx"),
      path(metadata_xlsx),
      emit: dds

    path("*.log"), emit: error_log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def script_path = "${workflow.projectDir}/modules/RNA"

    // Include group_name in pseudobulk log names so parallel groups do not
    // overwrite each other's logs.
    def LOG = group_name
        ? "${species}_${group_name}_CREATE_DDS.log"
        : "${species}_CREATE_DDS.log"

    """
    # Arg 1: count_data     — Salmon tximport RDS or raw pseudobulk count matrix
    # Arg 2: metadata_xlsx  — sample metadata Excel file
    # Arg 3: tx2gene_csv    — transcript/gene annotation CSV
    # Arg 4: design         — DESeq2 design formula
    # Arg 5: "."            — output directory (current work directory)

    Rscript "${script_path}/RNA_create_dds.R" \
        "${count_data}" \
        "${metadata_xlsx}" \
        "${tx2gene_csv}" \
        "${design}" \
        "." > "${LOG}" 2>&1

    EXIT_CODE=\$?

    if [ \$EXIT_CODE -ne 0 ]; then
        echo "❌ ERROR: RNA_create_dds.R failed for ${species}${group_name ? '|' + group_name : ''}" | tee -a "${LOG}" >&2
        exit \$EXIT_CODE
    fi

    echo "✅ SUCCESS: DESeq2 dataset created for ${species}${group_name ? '|' + group_name : ''}" >> "${LOG}"
    """
}
