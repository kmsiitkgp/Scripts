// =============================================================================
// PROCESS: CREATE_TX2GENE
// =============================================================================
// Purpose: Creates a transcript-to-gene mapping CSV from GTF annotation
//
// What it does:
//   - Parses the GTF file to extract transcript_id → gene_id mappings
//   - Saves result as tx2gene_<assembly>.<release>.csv
//   - Output passed to CREATE_TXI for tximeta/tximport gene-level summarization
//
// Why a specific filename pattern (tx2gene_*.csv)?
//   - Using storeDir means multiple species CSVs may coexist in the cache
//   - A specific pattern ensures .combine() in rnaseq.nf matches each species
//     to its own tx2gene file (not another species' file)
//
// Uses storeDir: built once per genome version, reused across runs
// =============================================================================

process CREATE_TX2GENE {

    tag "Creating tx2gene: ${species} (${assembly}.${release})"
    label 'process_low'                           // R parsing GTF; ~2-4GB RAM

    publishDir { "${params.proj_dir()}" },                         mode: 'copy',    pattern: "*.csv"
    publishDir { "${params.proj_dir()}/${species}/07.Logs" },      mode: 'copy',    pattern: "CREATE_TX2GENE.error.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), path(ref_fasta), path(ref_gtf), val(genome_version), val(assembly), val(release)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), path("tx2gene_*.csv"),    emit: tx2gene_tuple
    //path("CREATE_TX2GENE.error.log"),           emit: log    // storeDir: log not captured

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def script_path = "${workflow.projectDir}/modules"
    def LOG         = "CREATE_TX2GENE.error.log"

    """
    # Arg 1: "."          — save output CSV in current work directory
    # Arg 2: ref_gtf      — staged GTF file path
    # Arg 3: species      — used to label output file (e.g. tx2gene_Human_GRCh38.115.csv)
    # Arg 4: assembly     — genome assembly name (e.g. GRCh38)
    # Arg 5: release      — Ensembl release number (e.g. 115)
    Rscript "${script_path}/CREATE_TX2GENE.R" \
        "." \
        "${ref_gtf}" \
        "${species}" \
        "${assembly}" \
        "${release}" > "${LOG}" 2>&1 \
        || { echo "❌ ERROR: CREATE_TX2GENE.R failed for ${species}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: tx2gene created for ${species} (${assembly}.${release})" >> "${LOG}"
    """
}
