#!/usr/bin/env nextflow
nextflow.enable.dsl=2        // Enable DSL2 syntax (modern Nextflow with explicit workflow blocks)

// =========================================================================================
// IMPORT PROCESS MODULES
// =========================================================================================
// Modular architecture: Each process in separate file for maintainability

include { VALIDATE_INPUT }               from '../modules/validate_input.nf'
include { FASTQC as FASTQC_RAW }         from '../modules/fastqc.nf'
include { CELLRANGER_COUNT }             from '../modules/cellranger_count.nf'
include { MULTIQC }                      from '../modules/multiqc.nf'
//include { TEST_INDEX }                   from '../modules/test_index.nf'  // Debugging utility

// =========================================================================================
// PIPELINE HEADER
// =========================================================================================
// Displays configuration for user verification before execution

log.info """
    ===========================================
    SCRNA-SEQ PIPELINE
    ===========================================
    Project              : ${params.project}
    Species              : ${params.species}
    Genome               : ${params.genome_version}
    Fasta File           : ${params.ref_fasta()}
    GTF File             : ${params.ref_gtf()}
    BED File             : ${params.ref_bed()}
    Housekeeping         : ${params.housekeeping_bed()}

    PATHS:
    Reference Dir        : ${params.ref_dir()}
    Project Dir          : ${params.proj_dir()}
    Input (FastQ)        : ${params.fastq_dir()}
    Input (Raw FastQ)    : ${params.raw_fastq_dir()}
    Output (FastQC)      : ${params.fastqc_dir()}
    Output (CellRanger)    : ${params.cellranger_dir()}
    Output (MultiQC)     : ${params.multiqc_dir()}
    Logs                 : ${params.log_dir()}
    ===========================================
    """

// =========================================================================================
// MAIN WORKFLOW
// =========================================================================================
// Data flows through processes via channels
// Nextflow automatically parallelizes based on channel cardinality

workflow SCRNASEQ {

    // =====================================================================================
    // STEP 1: VALIDATE INPUT FASTQ FILES
    // =====================================================================================
    // Validates naming conventions and creates organized sample channels

    VALIDATE_INPUT(params.raw_fastq_dir())
    mode_ch         = VALIDATE_INPUT.out.mode.collect()         // [mode] → collected for RSeQC
    sample_fastq_ch = VALIDATE_INPUT.out.grouped_samples_ch     // [sample_id, [R1, R2]]
    sample_ch       = VALIDATE_INPUT.out.sample_names_ch        // [sample_id]

    // =====================================================================================
    // STEP 2: QUALITY CONTROL ON RAW READS
    // =====================================================================================
    // FastQC analyzes read quality before processing

    // Add read_type tag to channel: [sample_id, [fastqs], "raw"]
    fastqc_ch = sample_fastq_ch.map { sample_id, fastqs -> tuple(sample_id, fastqs, "raw") }
    FASTQC_RAW(fastqc_ch)

    // =====================================================================================
    // STEP 3: CELLRANGER COUNT
    // =====================================================================================

    // CRITICAL: Pre-join arguments to prevent cache invalidation
    // If params.CELLRANGER_ARGS().join(' ') called inside process → hash changes → resume fails
    cellranger_args     = params.CELLRANGER_ARGS().join(' ')
    raw_fastq_ch        = Channel.value(file(params.raw_fastq_dir(),        type: 'dir', checkIfExists: true))
    cellranger_index_ch = Channel.value(file(params.cellranger_index_dir(), type: 'dir', checkIfExists: true))
    CELLRANGER_COUNT(sample_ch, raw_fastq_ch, cellranger_index_ch, cellranger_args)

    // =====================================================================================
    // STEP 4: AGGREGATE QC REPORTS (MULTIQC)
    // =====================================================================================
    // Combines all QC outputs into single interactive HTML report

    multiqc_ch = Channel.empty()
        .mix(FASTQC_RAW.out.fastqc_zip)                     // FastQC reports
        .mix(CELLRANGER_COUNT.out.web_summary)
        .mix(CELLRANGER_COUNT.out.metric_summary)
        .collect()                                          // Wait for all samples

    // CRITICAL: Convert closures to string and pass into process to  prevent cache invalidation
    multiqc_title = params.multiqc_titlename()
    multiqc_file  = params.multiqc_filename()
    MULTIQC(multiqc_ch, multiqc_title, multiqc_file)
}

// =========================================================================================
// PIPELINE OVERVIEW
// =========================================================================================
//
// This pipeline performs comprehensive RNA-seq analysis:
//   1. Input validation and quality control (FastQC)
//   2. CellRanger count (STAR, Salmon, BED conversion)
//   3. Report aggregation (MultiQC - single HTML report)
//
// OUTPUT STRUCTURE:
// ${proj_dir}/
// ├── 01.FastQ/raw/                      # Raw input FASTQs
// ├── 02.FastQC/raw/                     # QC on raw reads
// ├── 03.CellRanger/                     # Transcript quantification
// │   ├── Sample1/                       # Full CellRanger output
// │   └── Sample2/                       # Full CellRanger output
// ├── 06.MultiQC/
// │   ├── Project_MultiQC_Report.html
// │   └── Project_MultiQC_Report_data/
// └── 07.Logs/                           # All error logs + reports
//     ├── trace.txt
//     ├── report.html
//     └── timeline.html
//
// DOWNSTREAM ANALYSIS:
//   Differential Expression:
//     - Option 1: Use STAR gene counts with DESeq2/edgeR directly
//     - Option 2: Import Salmon quant.sf files with tximport → DESeq2/edgeR
//
//   Visualization:
//     - Load BAM files into IGV for gene-level inspection
//     - Use RSeQC or deepTools for coverage plots
//
//   Splice Analysis:
//     - Use STAR SJ.out.tab files with rMATS or LeafCutter
//
// QUALITY CONTROL THRESHOLDS (check MultiQC report):
//   - Mapping rate: >70% (good), 60-70% (acceptable), <60% (investigate)
//   - Duplication: <50% (good), 50-70% (acceptable), >70% (low complexity)
//   - Read distribution: >50% CDS exons, <10% introns, <5% intergenic
//   - 3' bias: <3 (good), 3-5 (acceptable), >5 (degraded RNA)
//
// RESUMING FAILED RUNS:
//   Nextflow caches completed processes in work/ directory
//   To resume: bash run_nextflow.sh (script includes -resume flag)
//   Cache is invalidated if: process code changes, input files change, or work/ deleted
//
/// CLEANING UP:
//   After successful run and verification:
//     rm -rf ${work_dir}/*              # Free disk space
//   Warning: Cannot resume after deleting work directory!
//
// COMMON ISSUES:
//   "No such file or directory":
//     → Check paths in project_info.yaml
//     → Verify Singularity bind mounts (run_nextflow.sh prints mappings)
//
//   Out of memory:
//     → Increase memory in process labels (nextflow.config)
//     → Reduce parallel jobs with maxForks in nextflow.config
//
//   Process won't resume:
//     → Check if params.FUNCTION().join(' ') called inside process (should be in workflow)
//     → Verify work/ directory not deleted
//     → Review .nextflow.log for cache invalidation reasons
//
// EXTENDING THE PIPELINE:
//   1. Create new process module in modules/ directory
//   2. Add include statement above (top of file)
//   3. Call process in workflow block (add to appropriate step)
//   4. Configure publishDir in nextflow.config
//   5. Add outputs to multiqc_ch if process generates QC data
//
// For detailed documentation on each process:
//   docs/validate_input.md, docs/star_align.md, docs/salmon_quant.md, etc.
//
// =========================================================================================
