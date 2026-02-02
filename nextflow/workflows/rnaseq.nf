#!/usr/bin/env nextflow
nextflow.enable.dsl=2        // Enable DSL2 syntax (modern Nextflow with explicit workflow blocks)

// =========================================================================================
// IMPORT PROCESS MODULES
// =========================================================================================
// Modular architecture: Each process in separate file for maintainability

include { VALIDATE_INPUT }               from '../modules/validate_input.nf'
include { FASTQC as FASTQC_RAW }         from '../modules/fastqc.nf'
include { FASTQC as FASTQC_TRIMMED }     from '../modules/fastqc.nf'     // Reusable process with alias
include { STAR_INDEX }                   from '../modules/star_index.nf'
include { EXTRACT_GENTROME }             from '../modules/salmon_index.nf'
include { SALMON_INDEX }                 from '../modules/salmon_index.nf'
include { RSEQC_BED }                    from '../modules/rseqc_bed.nf'
include { SALMON_QUANT }                 from '../modules/salmon_quant.nf'
include { STAR_ALIGN }                   from '../modules/star_align.nf'
include { SAMBAMBA_PREP }                from '../modules/sambamba_prep.nf'
include { RSEQC }                        from '../modules/rseqc.nf'
include { MULTIQC }                      from '../modules/multiqc.nf'
//include { TEST_INDEX }                   from '../modules/test_index.nf'  // Debugging utility

// =========================================================================================
// PIPELINE HEADER
// =========================================================================================
// Displays configuration for user verification before execution

log.info """
    ===========================================
    RNA-SEQ PIPELINE
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
    STAR Index           : ${params.star_index_dir()}
    Salmon Index         : ${params.salmon_index_dir()}
    Project Dir          : ${params.proj_dir()}
    Input (FastQ)        : ${params.fastq_dir()}
    Input (Raw FastQ)    : ${params.raw_fastq_dir()}
    Output (FastQC)      : ${params.fastqc_dir()}
    Output (Salmon)      : ${params.salmon_dir()}
    Output (STAR)        : ${params.star_dir()}
    Output (RSeQC)       : ${params.rseqc_dir()}
    Output (MultiQC)     : ${params.multiqc_dir()}
    Logs                 : ${params.log_dir()}
    ===========================================
    """

// =========================================================================================
// MAIN WORKFLOW
// =========================================================================================
// Data flows through processes via channels
// Nextflow automatically parallelizes based on channel cardinality

workflow RNASEQ {

    // =====================================================================================
    // STEP 1: VALIDATE INPUT FASTQ FILES
    // =====================================================================================
    // Validates naming conventions and creates organized sample channels

    VALIDATE_INPUT(params.raw_fastq_dir())
    mode_ch         = VALIDATE_INPUT.out.mode.collect()          // [mode] → collected for RSeQC
    sample_fastq_ch = VALIDATE_INPUT.out.grouped_samples_ch      // [sample_id, [R1, R2]]
    sample_ch        = VALIDATE_INPUT.out.sample_names_ch        // [sample_id]

    // =====================================================================================
    // STEP 2: QUALITY CONTROL ON RAW READS
    // =====================================================================================
    // FastQC analyzes read quality before processing

    // Add read_type tag to channel: [sample_id, [fastqs], "raw"]
    fastqc_ch = sample_fastq_ch.map { sample_id, fastqs -> tuple(sample_id, fastqs, "raw") }
    FASTQC_RAW(fastqc_ch)

    // =====================================================================================
    // STEP 3: BUILD REFERENCE INDEXES
    // =====================================================================================

    // Create value channels for reference files (singleton channels, reusable)
    // checkIfExists: true → pipeline fails immediately if files missing
    ref_fasta_ch = Channel.value(file(params.ref_fasta(), checkIfExists: true))
    ref_gtf_ch   = Channel.value(file(params.ref_gtf(), checkIfExists: true))

    // STAR index for genome alignment
    STAR_INDEX(ref_fasta_ch, ref_gtf_ch)
    star_index_ch = STAR_INDEX.out.star_index_dir.collect()

    // Salmon index for transcript quantification
    EXTRACT_GENTROME(ref_fasta_ch, ref_gtf_ch)
    SALMON_INDEX(EXTRACT_GENTROME.out.gentrome.collect(), EXTRACT_GENTROME.out.decoy.collect())
    salmon_index_ch = SALMON_INDEX.out.salmon_index_dir.collect()

    // BED files for RSeQC analysis
    RSEQC_BED(ref_gtf_ch)
    ref_bed_ch          = RSEQC_BED.out.ref_bed.collect()
    housekeeping_bed_ch = RSEQC_BED.out.housekeeping_bed.collect()

    // =====================================================================================
    // STEP 4: TRANSCRIPT QUANTIFICATION (SALMON)
    // =====================================================================================
    // Fast, alignment-free quantification
    // Runs in parallel with STAR (independent processes)

    // CRITICAL: Pre-join arguments to prevent cache invalidation
    // If params.SALMON_ARGS().join(' ') called inside process → hash changes → resume fails
    salmon_args = params.SALMON_ARGS().join(' ')
    SALMON_QUANT(sample_fastq_ch, salmon_index_ch, salmon_args)

    // =====================================================================================
    // STEP 5: GENOME ALIGNMENT (STAR)
    // =====================================================================================
    // Splice-aware alignment, generates BAM files for visualization and QC

    // CRITICAL: Pre-join arguments to prevent cache invalidation
    star_args = params.STAR_ARGS().join(' ')
    STAR_ALIGN(sample_fastq_ch, star_index_ch, star_args)
    sample_unindexed_bam_ch = STAR_ALIGN.out.bam_unindexed  // [sample_id, bam]

    // Index BAM files and create subsampled versions for faster RSeQC
    SAMBAMBA_PREP(sample_unindexed_bam_ch)
    sample_indexed_bam_ch = SAMBAMBA_PREP.out.bam_indexed  // [sample_id, bam, bai, 1M.bam, 1M.bai, read_len]

    // =====================================================================================
    // STEP 6: ALIGNMENT QUALITY CONTROL (RSEQC)
    // =====================================================================================
    // Comprehensive QC: read distribution, gene body coverage, junction analysis, etc.

    RSEQC(sample_indexed_bam_ch, ref_bed_ch, housekeeping_bed_ch, mode_ch)

    // =====================================================================================
    // STEP 7: AGGREGATE QC REPORTS (MULTIQC)
    // =====================================================================================
    // Combines all QC outputs into single interactive HTML report

    multiqc_ch = Channel.empty()
        .mix(FASTQC_RAW.out.fastqc_zip)                         // FastQC reports
        .mix(SALMON_QUANT.out.salmon_quant_dir.map { it[1] })   // Salmon dirs (extract from tuple)
        .mix(STAR_ALIGN.out.gene_counts)                        // ReadsPerGene.out.tab
        .mix(STAR_ALIGN.out.sj_tab)                             // Splice junctions
        .mix(STAR_ALIGN.out.star_log)                           // Alignment stats
        .mix(RSEQC.out.rseqc_logs)                              // RSeQC outputs
        .collect()                                              // Wait for all samples

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
//   2. Reference genome indexing (STAR, Salmon, BED conversion)
//   3. Transcript quantification (Salmon - fast, alignment-free)
//   4. Read alignment (STAR - splice-aware, generates BAM)
//   5. Alignment QC (RSeQC - detects biases and issues)
//   6. Report aggregation (MultiQC - single HTML report)
//
// OUTPUT STRUCTURE:
// ${proj_dir}/
// ├── 01.FastQ/raw/                      # Raw input FASTQs
// ├── 02.FastQC/raw/                     # QC on raw reads
// ├── 03.Salmon/                         # Transcript quantification
// │   ├── Sample1/                       # Full Salmon output
// │   └── quant_files/                   # Collected quant.sf files
// ├── 04.STAR/                           # Alignment outputs
// │   ├── gene_counts/                   # ReadsPerGene.out.tab files
// │   ├── splice_junction/               # SJ.out.tab files
// │   ├── alignment_stats/               # Log.final.out files
// │   └── bam                              # BAM + BAI files
// ├── 05.RSEQC/                          # Organized by analysis type
// │   ├── 01_read_distribution/
// │   ├── 02_inner_distance/
// │   ├── 03_junction_annotation/
// │   └── ... (09 subdirectories total)
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
// CLEANING UP:
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
