process MERGE_STAR_COUNTS {

    tag "Merging counts: ${species}"
    label 'process_low'                            // STAR requires 30-50GB RAM for human

    publishDir { "${params.proj_dir()}/${species}/04.STAR/" },    mode: 'copy',    pattern: "STAR_Gene_counts.xlsx"
    publishDir { "${params.proj_dir()}/${species}/07.Logs" },     mode: 'copy',    pattern: "MERGE_STAR_COUNTS.error.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), path(gene_counts)
    // sample_id        : Sample identifier (e.g., "Sample1")
    // gene_counts      : list of ReadsPerGene.out.tab files

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), path("STAR_Gene_counts.xlsx"),    emit: star_counts                 // Alignment stats
    path("MERGE_STAR_COUNTS.error.log"),                             emit: error_log,    optional: true // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    // This points to the modules folder relative to your project root
    def script_path = "${workflow.projectDir}/modules"
    def LOG = "MERGE_STAR_COUNTS.error.log"

    """
    # We pass '.' because Nextflow staged all 'gene_counts' into the current folder
    # We pass '.' as the second arg so the Excel file is saved in the current folder
    Rscript ${script_path}/MERGE_STAR_COUNTS.R \
        "." \
        "." > ${LOG} 2>&1
    """
}