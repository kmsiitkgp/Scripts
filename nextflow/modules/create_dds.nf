process CREATE_DDS {

    tag "Creating dds: ${species}"
    label 'process_medium'                      // STAR indexing requires 30-50GB RAM for human

    publishDir { "${params.proj_dir()}/${species}/08.DESeq2" },    mode: 'copy',    pattern: "*.{rds,xlsx}"
    publishDir { "${params.proj_dir()}/${species}/07.Logs" },      mode: 'copy',    pattern: "CREATE_DDS.error.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), path(txi), path(tx2gene)
    path(metadata)
    val(design)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), path("DESeq2_dds.rds"), path(tx2gene),    emit: dds
    path("*.xlsx"),                                               emit: count_xlsx
    path("CREATE_DDS.error.log"),                                 emit: error_log    // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    // This points to the modules folder relative to your project root
    def script_path = "${workflow.projectDir}/modules"
    def LOG = "CREATE_DDS.error.log"

    """
    # We pass '.' because Nextflow staged all 'gene_counts' into the current folder
    # We pass '.' as the second arg so the Excel file is saved in the current folder
    Rscript ${script_path}/CREATE_DDS.R \
        "${txi}" \
        "${metadata}" \
        "${tx2gene}" \
        "." \
        "${design}" > ${LOG} 2>&1
    """
}