process CREATE_TXI {

    tag "Creating txi: ${species}_${type}"
    label 'process_medium'                      // STAR indexing requires 30-50GB RAM for human

    publishDir { "${params.proj_dir()}/${species}_${type}/03.Salmon" },    mode: 'copy',    pattern: "*.{rds,txt}"
    publishDir { "${params.proj_dir()}/${species}_${type}/07.Logs" },      mode: 'copy',    pattern: "CREATE_TXI.error.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), val(type), path(quant_files), path(tx2gene_csv)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), val(type), path("Salmon_txi.rds"),            emit: txi
    tuple val(species), val(type), path("Salmon_Gene_Counts.txt"),    emit: gene_counts
    tuple val(species), val(type), path("Salmon_TPM_Values.txt"),     emit: tpm
    path("CREATE_TXI.error.log"),                                     emit: error_log    // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    // This points to the modules folder relative to your project root
    def script_path = "${workflow.projectDir}/modules"

    """
    # We pass '.' because Nextflow staged all 'gene_counts' into the current folder
    # We pass '.' as the second arg so the Excel file is saved in the current folder
    Rscript ${script_path}/CREATE_TXI.R ${species} . . ${tx2gene_csv} > CREATE_TXI.error.log 2>&1
    """
}