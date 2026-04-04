// =============================================================================
// PROCESS: MERGE_STAR_COUNTS
// =============================================================================
// Purpose: Merges per-sample STAR gene count files into a single Excel workbook
//
// What it does:
//   - Receives all ReadsPerGene.out.tab files for a species (staged into work dir)
//   - Calls R script to read, merge, and format them into one Excel file
//   - Output used for differential expression analysis downstream
// =============================================================================

process MERGE_STAR_COUNTS {

    tag "Merging STAR counts: ${species}"
    label 'process_low'                           // R with file I/O only; minimal compute

    publishDir { "${params.proj_dir()}/${species}/04.STAR/" },    mode: 'copy',    pattern: "STAR_Gene_counts.xlsx"
    publishDir { "${params.proj_dir()}/${species}/07.Logs" },     mode: 'copy',    pattern: "MERGE_STAR_COUNTS.error.log"

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
    path("MERGE_STAR_COUNTS.error.log"),                  emit: error_log,  optional: true // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def script_path = "${workflow.projectDir}/modules"
    def LOG         = "MERGE_STAR_COUNTS.error.log"

    """
    # Both args are ".": Nextflow staged all gene_counts files into the current
    # work directory, and the R script saves its Excel output there too
    Rscript "${script_path}/MERGE_STAR_COUNTS.R" \
        "." \
        "." > "${LOG}" 2>&1 \
        || { echo "❌ ERROR: MERGE_STAR_COUNTS.R failed for ${species}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: STAR counts merged for ${species}" >> "${LOG}"
    """
}
