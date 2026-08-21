process SCRNA_MERGE_AND_PLOT_QC {

    tag "Merging all samples"
    label 'process_heavy'
    label 'omics_r'

    publishDir = [
        [path: { "${params.proj_dir()}/${params.species}/05.Seurat/Whole" },          mode: 'copy', pattern: "*.rds"],
        [path: { "${params.proj_dir()}/${params.species}/05.Seurat/Whole/Plots/QC" }, mode: 'copy', pattern: "*.pdf"],
        [path: { "${params.proj_dir()}/${params.species}/Logs" },                     mode: 'copy', pattern: "*.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    path(sample_rds_list)
    // sample_rds_list : All per-sample Seurat .rds files collected from
    //                 CALC_QC_AND_IDENTIFY_SINGLETS — one per sample, UNFILTERED
    //                 (every barcode, all Barcode_QC categories retained)
    path(metadata_xlsx)
    // metadata_xlsx : params.metadata_file — sample-level annotations (xlsx)
    //                 joined to merged Seurat object by Sample or Unique_ID

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    // Two checkpoints now instead of one:
    //   pre_filtered_seurat.rds  : non-empty barcodes from every sample, all classification columns
    //                         intact — nothing removed. Kept for QC diagnostics and so
    //                         cutoffs can be revisited without rerunning upstream steps.
    //   filtered_seurat.rds : pre_filtered_seurat.rds subset to Barcode_QC == "Singlet" — this is
    //                         the only file downstream (SCT_TRANSFORM etc.) actually reads.

    path("pre_filtered_seurat.rds"),      emit: pre_filtered_rds
    path("filtered_seurat.rds"),          emit: filtered_rds
    path("*.pdf"),                        emit: plots
    path("SCRNA_MERGE_AND_PLOT_QC.log"),  emit: log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def script_path = "${workflow.projectDir}/modules"
    def LOG         = "SCRNA_MERGE_AND_PLOT_QC.log"

    // Handles both single and multiple rds files
    // Comma-separated to avoid issues with spaces in paths
    def sample_rds_str = sample_rds_list instanceof List ? sample_rds_list.join(',') : sample_rds_list.toString()

    """
    Rscript "${script_path}/SCRNA/SCRNA_merge_and_plot_qc.R" \
    "${sample_rds_str}" \
    "${metadata_xlsx}" \
    "." > "${LOG}" 2>&1

    EXIT_CODE=\$?
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "❌ ERROR: SCRNA_MERGE_AND_PLOT_QC failed" | tee -a "${LOG}" >&2
        exit \$EXIT_CODE
    fi

    echo "✅ SUCCESS: SCRNA_MERGE_AND_PLOT_QC completed" >> "${LOG}"
    """
}
