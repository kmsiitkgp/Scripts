process CREATE_DDS {

    tag "Creating dds: ${species}_${type}"
    label 'process_medium'                      // STAR indexing requires 30-50GB RAM for human

    publishDir { "${params.proj_dir()}/${species}_${type}/08.DESeq2" },    mode: 'copy',    pattern: "DESeq2_dds.rds"
    publishDir { "${params.proj_dir()}/${species}_${type}/07.Logs" },      mode: 'copy',    pattern: "CREATE_DDS.error.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), val(type), path(txi)
    path(metadata)
    val(design)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("DESeq2_dds.rds"),            emit: dds
    path("CREATE_DDS.error.log"),      emit: error_log    // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    // This points to the modules folder relative to your project root
    def script_path = "${workflow.projectDir}/modules"

    """
    # We pass '.' because Nextflow staged all 'gene_counts' into the current folder
    # We pass '.' as the second arg so the Excel file is saved in the current folder
    Rscript ${script_path}/CREATE_DDS.R ${txi} ${metadata} . ${design} > CREATE_DDS.error.log 2>&1
    """
}