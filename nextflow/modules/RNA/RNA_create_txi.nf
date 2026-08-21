// =============================================================================
// PROCESS: RNA_CREATE_TXI
// =============================================================================
// Purpose: Imports Salmon quantification results into R using tximeta/tximport
//
// What it does:
//   - Reads all sample quant.sf files (staged into work directory)
//   - Links transcripts to genes using the tx2gene mapping
//   - Produces a txi object (Salmon_txi.rds) for DESeq2 input
//   - Also exports gene-level count and TPM tables as plain text
// =============================================================================

process RNA_CREATE_TXI {

    tag "Creating txi: ${species}"
    label 'process_medium'                        // R with tximeta; ~8GB RAM
    label 'omics_r'

    publishDir = [
        [path: { "${params.proj_dir()}/${species}/03.Salmon" }, mode: 'copy', pattern: "*.{rds,txt}"],
        [path: { "${params.proj_dir()}/${species}/Logs" },      mode: 'copy', pattern: "RNA_CREATE_TXI.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), path(quant_files), path(tx2gene_csv)
    // species     : "Human" or "Mouse"
    // quant_files : All sample quant.sf files staged into work directory
    // tx2gene_csv : Transcript-to-gene mapping from CREATE_TX2GENE

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), path("Salmon_txi.rds"), path("${tx2gene_csv}"),    emit: txi             // txi object + tx2gene for CREATE_DDS
    tuple val(species), path("Salmon_Gene_Counts.txt"),                    emit: gene_counts     // Gene-level read counts
    tuple val(species), path("Salmon_TPM_Values.txt"),                     emit: tpm             // Gene-level TPM values
    path("RNA_CREATE_TXI.log"),                                          emit: error_log       // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def script_path = "${workflow.projectDir}/modules/RNA"
    def LOG         = "RNA_CREATE_TXI.log"

    """
    # Arg 1: species      — used to label output files
    # Arg 2: "."          — directory containing staged quant.sf files (current work dir)
    # Arg 3: "."          — directory to save output files (current work dir)
    # Arg 4: tx2gene_csv  — transcript-to-gene mapping file
    Rscript "${script_path}/RNA_create_txi.R" \
        "${species}" \
        "." \
        "." \
        "${tx2gene_csv}" > "${LOG}" 2>&1

    EXIT_CODE=\$?
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "❌ ERROR: RNA_CREATE_TXI failed for ${species}" | tee -a "${LOG}" >&2
        exit \$EXIT_CODE
    fi

    echo "✅ SUCCESS: txi created for ${species}" >> "${LOG}"
    """
}
