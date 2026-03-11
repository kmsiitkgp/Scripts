process CREATE_TXI {

    tag "Creating txi: ${species}"
    label 'process_medium'                      // STAR indexing requires 30-50GB RAM for human

    publishDir { "${params.proj_dir()}/${species}/03.Salmon" },    mode: 'copy',    pattern: "*.{rds,txt}"
    publishDir { "${params.proj_dir()}/${species}/07.Logs" },      mode: 'copy',    pattern: "CREATE_TXI.error.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), path(quant_files), path(tx2gene_csv)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), path("Salmon_txi.rds"),  path("${tx2gene_csv}"),            emit: txi
    tuple val(species), path("Salmon_Gene_Counts.txt"),    emit: gene_counts
    tuple val(species), path("Salmon_TPM_Values.txt"),     emit: tpm
    path("CREATE_TXI.error.log"),                                     emit: error_log    // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    // This points to the modules folder relative to your project root
    def script_path = "${workflow.projectDir}/modules"
    def LOG = "CREATE_TXI.error.log"

    """
    # We pass '.' because Nextflow staged all 'gene_counts' into the current folder
    # We pass '.' as the second arg so the Excel file is saved in the current folder
    Rscript ${script_path}/CREATE_TXI.R \
        "${species}" \
        "." \
        "." \
        "${tx2gene_csv}" > ${LOG} 2>&1
    """
}