process CREATE_TX2GENE {

    tag "Creating tx2gene ${species}"
    label 'process_high'                      // STAR indexing requires 30-50GB RAM for human

    publishDir { "${params.proj_dir()}" },                               mode: 'copy',    pattern: "*.csv"
    publishDir { "${params.proj_dir()}/${species}/07.Logs" },    mode: 'copy',    pattern: "CREATE_TXI.error.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), path(ref_fasta), path(ref_gtf), val(genome_version), val(assembly), val(release)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    // WARNING: Do NOT use a generic wildcard like path("*.csv") here if using storeDir.
    // If multiple species files (Human/Mouse) exist in the storeDir, .combine() in 
    // the workflow will stage BOTH species' CSVs into the same work directory, 
    // leading to species-mismatch errors in downstream R scripts.
    // Use a specific naming pattern to ensure a 1-to-1 species match.
    // genome_version matches with suffix extracted by R from gtf name

    tuple val(species), path("tx2gene_*.csv"),    emit: tx2gene_tuple
    //path("CREATE_TX2GENE.error.log"),           emit: error_log         // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    // This points to the modules folder relative to your project root
    def script_path = "${workflow.projectDir}/modules"
    def LOG = "CREATE_TX2GENE.error.log"

    """
    # Arg 1: "." (Save in current working directory)
    # Arg 2: "${ref_gtf}" (The path to the staged GTF)
    echo ${ref_gtf}
    Rscript ${script_path}/CREATE_TX2GENE.R \
        "." \
        "${ref_gtf}" \
        "${species}" \
        "${assembly}" \
        "${release}" > ${LOG} 2>&1
    """
}