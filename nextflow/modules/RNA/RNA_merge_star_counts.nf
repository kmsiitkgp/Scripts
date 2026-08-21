// =============================================================================
// PROCESS: RNA_MERGE_STAR_COUNTS
// =============================================================================
// Purpose: Merges per-sample STAR gene count files into a single Excel workbook
//
// What it does:
//   - Receives all ReadsPerGene.out.tab files for a species (staged into work dir)
//   - Calls R script to read, merge, and format them into one Excel file
//   - Output used for differential expression analysis downstream
// =============================================================================

process RNA_MERGE_STAR_COUNTS {

    tag "Merging STAR counts: ${species}"
    label 'process_low'                           // R with file I/O only; minimal compute
    label 'omics_r'

    publishDir = [
        [path: { "${params.proj_dir()}/${species}/04.STAR/" }, mode: 'copy', pattern: "STAR_Gene_counts.xlsx"],
        [path: { "${params.proj_dir()}/${species}/Logs" },     mode: 'copy', pattern: "RNA_MERGE_STAR_COUNTS.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), path(gene_counts)
    // species     : "Human" or "Mouse"
    // gene_counts : List of all ReadsPerGene.out.tab files staged into work directory

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), path("STAR_Gene_counts.xlsx"),    emit: star_counts              // Merged count matrix
    path("RNA_MERGE_STAR_COUNTS.log"),                  emit: error_log,  optional: true // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def script_path = "${workflow.projectDir}/modules/RNA"
    def LOG         = "RNA_MERGE_STAR_COUNTS.log"

    """
    # Both args are ".": Nextflow staged all gene_counts files into the current
    # work directory, and the R script saves its Excel output there too
    Rscript "${script_path}/RNA_merge_star_counts.R" \
        "." \
        "." > "${LOG}" 2>&1

    EXIT_CODE=\$?
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "❌ ERROR: RNA_MERGE_STAR_COUNTS failed for ${species}" | tee -a "${LOG}" >&2
        exit \$EXIT_CODE
    fi

    echo "✅ SUCCESS: STAR counts merged for ${species}" >> "${LOG}"
    """
}
