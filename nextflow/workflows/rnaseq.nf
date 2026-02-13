#!/usr/bin/env nextflow
nextflow.enable.dsl=2        // Enable DSL2 syntax (modern Nextflow with explicit workflow blocks)

// =========================================================================================
// IMPORT PROCESS MODULES
// =========================================================================================
// Modular architecture: Each process in separate file for maintainability

include { RENAME_FASTQS }                from '../modules/rename_fastq.nf'
include { GENERATE_MD5 }                 from '../modules/generate_md5.nf'
include { VALIDATE_INPUT }               from '../modules/validate_input.nf'
include { FASTQC as FASTQC_RAW }         from '../modules/fastqc.nf'
include { FASTQC as FASTQC_TRIMMED }     from '../modules/fastqc.nf'           // Reusable process with alias
include { FETCH_ENSEMBL_VERSIONS }       from '../modules/fetch_ensembl_versions.nf'
include { DOWNLOAD_ENSEMBL_REF }         from '../modules/download_ensembl_ref.nf'
//include { PREPARE_XENOGRAFT_REF }        from '../modules/prepare_xenograft_ref.nf'
include { XENGSORT_INDEX }               from '../modules/xengsort_index.nf'
include { XENGSORT_CLASSIFY }            from '../modules/xengsort_classify.nf'
include { STAR_INDEX }                   from '../modules/star_index.nf'
include { EXTRACT_GENTROME }             from '../modules/extract_gentrome.nf'
include { SALMON_INDEX }                 from '../modules/salmon_index.nf'
include { RSEQC_BED }                    from '../modules/rseqc_bed.nf'
include { SALMON_QUANT }                 from '../modules/salmon_quant.nf'
include { STAR_ALIGN }                   from '../modules/star_align.nf'
include { SAMBAMBA_PREP }                from '../modules/sambamba_prep.nf'
include { RSEQC }                        from '../modules/rseqc.nf'
include { MULTIQC }                      from '../modules/multiqc.nf'
//include { TEST_INDEX }                   from '../modules/test_index.nf'       // Debugging utility
//include { TEST_PUBLISHDIR }              from '../modules/test_publishdir.nf'  // Debugging utility

// =========================================================================================
// MAIN WORKFLOW
// =========================================================================================
// Data flows through processes via channels
// Nextflow automatically parallelizes based on channel cardinality

workflow RNASEQ {

    // =========================================================================================
    // PIPELINE HEADER
    // =========================================================================================
    // Displays configuration for user verification before execution

    log.info """
        ===========================================
        RNA-SEQ PIPELINE
        ===========================================
        Experiment               : ${params.expt}
        Species                  : ${params.species}
        Project                  : ${params.project}
        Mapping File             : ${params.map_file}
        Input FastQ Dir          : ${params.read_dir}

        Reference Dir            : ${params.ref_dir}
        ===========================================
        """.stripIndent()

    // =====================================================================================
    // STEP -1: RENAME FASTQ FILES
    // =====================================================================================

    // Check if the map file actually exists on disk
    map_file_path = params.map_file ? file(params.map_file) : null

    if (map_file_path && map_file_path.exists()) {

        map_ch        = Channel.value(map_file_path)
        srr_fastq_ch  = Channel.fromPath("${params.proj_dir()}/downloads/*.{fastq,fq}.gz")

        // CRITICAL: Use .collect() so all files go to ONE rename process
        RENAME_FASTQS(srr_fastq_ch.collect(), map_ch)

        // Use .flatten() here to break the list into 48 individual items
        raw_fastq_ch = RENAME_FASTQS.out.renamed_fastqs.flatten()

    } else {
        log.info "No map file found. Passing raw FASTQs directly to validation."
        raw_fastq_ch  = Channel.fromPath("${params.read_dir}/*.{fastq,fq}.gz")
    }
    if (params.stop_after == 'RENAME_FASTQS') {
        log.info "Stopping pipeline after RENAME_FASTQS as requested."
        return
    }

    // =====================================================================================
    // STEP 0: GENERATE MD5
    // =====================================================================================
    // Validates naming conventions and creates organized sample channels

    GENERATE_MD5(raw_fastq_ch)

    // Merge all MD5 into one file and save it to your project dir
    manifest_ch = GENERATE_MD5.out.md5_file.collectFile(
            name: 'manifest_md5.txt',
            newLine: true,
            storeDir: "${params.proj_dir()}"
        )
    if (params.stop_after == 'GENERATE_MD5') {
        manifest_ch.subscribe {
            log.info "✅ Manifest created in ${params.proj_dir()}"
            log.info "Stopping pipeline after GENERATE_MD5 as requested."
        }
        return
    }

    // =====================================================================================
    // STEP 1: VALIDATE INPUT FASTQ FILES
    // =====================================================================================
    // Validates naming conventions and creates organized sample channels

    // Workflow can accept path strings like VALIDATE_INPUT(params.raw_fastq_dir())
    // Here we pass the fastqs as channel instead
    VALIDATE_INPUT(raw_fastq_ch.collect())
    mode_ch          = VALIDATE_INPUT.out.mode.collect()          // [mode] → collected for RSeQC
    sample_fastq_ch  = VALIDATE_INPUT.out.grouped_samples_ch      // [sample_id, [R1, R2]]
    sample_ch        = VALIDATE_INPUT.out.sample_names_ch         // [sample_id]
    if (params.stop_after == 'VALIDATE_INPUT') {
        log.info "Stopping pipeline after VALIDATE_INPUT as requested."
        return
    }

    // Create value channels for reference files (singleton channels, reusable)
    // Without .collect(), output will be a 'Queue Channel'.
    // If you have 10 samples, only 1st sample will find an index; the other 9 will wait forever.
    // With .collect(), output will be a 'Value Channel'.
    // If you have 1 index and 10 samples, all 10 samples can use it simultaneously.

    // =================================================================================
    // STEP 2: REFERENCE METADATA FETCHING
    // =================================================================================

    // Determine which species metadata to fetch based on user parameters
    if (params.species == 'Xenograft') {
        // Xenograft requires both parent species AND the combined naming convention
        target_species_list = ['Human', 'Mouse', 'Xenograft']
    } else {
        // Standard run only requires the single requested species
        target_species_list = [params.species]
    }

    species_list_ch = Channel.fromList(target_species_list)
    FETCH_ENSEMBL_VERSIONS(species_list_ch)

    // Parse the pipe-delimited metadata file into an accessible list
    // Structure: [Species, FA_Name, GTF_Name, Version, Assembly, Release]
    parsed_metadata_ch = FETCH_ENSEMBL_VERSIONS.out.meta
        .splitText()
        .map { it.trim().split('\\|') }


    // =================================================================================
    // STEP 3: REFERENCE FILE DOWNLOADS
    // =================================================================================

    // Filter out "Xenograft" because it is a virtual entry; physical files only exist
    // for Human and Mouse.
    download_queue_ch = parsed_metadata_ch.filter { it[0] != "Xenograft" }

    DOWNLOAD_ENSEMBL_REF(
        download_queue_ch.map { it[0] }, // species_name
        download_queue_ch.map { it[1] }, // fasta_filename
        download_queue_ch.map { it[2] }, // gtf_filename
        download_queue_ch.map { it[5] }  // ensembl_release
    )

    // =================================================================================
    // STEP 4: MASTER REFERENCE TUPLE CONSTRUCTION
    // =================================================================================

    // Join the metadata attributes with the actual downloaded file paths.
    // Pair every fasta and gtf with its specific species
    // Using 'by: 0' means "match where the first element is the same, in our case `species`
    ref_ch = download_queue_ch
        .map { it -> [ it[0], it[3] ] } // [Species, Version]
        .combine(DOWNLOAD_ENSEMBL_REF.out.ref_tuple, by: 0)
        .map { species, version, fasta, gtf -> tuple(species, fasta, gtf, version) }

    // ref_ch result:
    // [Species, FA_File, GTF_File, Full_Version]

    // =====================================================================================
    // STEP 5: BUILD REFERENCE INDEXES
    // =====================================================================================

    // STAR index for genome alignment
    STAR_INDEX(ref_ch)
    if (params.stop_after == 'STAR_INDEX') {
        log.info "Stopping pipeline after STAR_INDEX as requested."
        return
    }

    // Salmon index for transcript quantification
    EXTRACT_GENTROME(ref_ch)
    decoy_gentrome_ch = EXTRACT_GENTROME.out.decoy_gentrome_tuple

    SALMON_INDEX(decoy_gentrome_ch)
    if (params.stop_after == 'SALMON_INDEX') {
        log.info "Stopping pipeline after SALMON_INDEX as requested."
        return
    }

    // BED files for RSeQC analysis
    RSEQC_BED(ref_ch)
    if (params.stop_after == 'RSEQC_BED') {
        log.info "Stopping pipeline after RSEQC_BED as requested."
        return
    }

    // =================================================================================
    // STEP 6: SPECIES SEPARATION & MAPPING PREPARATION
    // =================================================================================

    if (params.species == 'Xenograft') {

        // --- XENOGRAFT WORKFLOW BRANCH ---

        // Extract physical FASTA paths for the index builder
        human_fasta_ch = ref_ch.filter { it[0] == "Human" }.map { it[1] }
        mouse_fasta_ch = ref_ch.filter { it[0] == "Mouse" }.map { it[1] }

        // Extract the composite version string for the Xenograft (e.g., GRCh38.115_GRCm39.115)
        xeno_version_ch = parsed_metadata_ch.filter { it[0] == "Xenograft" }.map { it[3] }.first()

        // Generate the K-mer index for species separation
        XENGSORT_INDEX(human_fasta_ch, mouse_fasta_ch, xeno_version_ch)

        // Classify raw reads into Graft (Human) and Host (Mouse) bins
        XENGSORT_CLASSIFY(sample_fastq_ch, XENGSORT_INDEX.out.xengsort_index_dir.collect(), xeno_version_ch)

        // Manually add the "Human" and "Mouse" labels to the outputs and
        graft_labeled_ch = XENGSORT_CLASSIFY.out.graft_fastqs
            .map { sample_id, fastqs -> tuple("Human", "split", sample_id, fastqs) }

        host_labeled_ch = XENGSORT_CLASSIFY.out.host_fastqs
            .map { sample_id, fastqs -> tuple("Mouse", "split", sample_id, fastqs) }

        full_fastq_ch = sample_fastq_ch
            .map { sample_id, fastqs -> tuple("Human", "full", sample_id, fastqs) }

        // Merge both streams into a single downstream mapping channel
        all_fastq_ch = graft_labeled_ch
            .concat(host_labeled_ch)
            .concat(full_fastq_ch)

        //  Bridge the separated reads back to their specific species
        sample_fastq_metadata_ch = all_fastq_ch
            .combine(ref_ch, by: 0)

    } else {

        // --- STANDARD WORKFLOW BRANCH ---
        // For non-xenograft runs, we map the input reads directly to the single reference
        sample_fastq_metadata_ch = sample_fastq_ch
            .combine(ref_ch)
            .map { sample_id, fastqs, species, fasta, gtf, version ->
                tuple(species, "full", sample_id, fastqs, fasta, gtf, version) }
    }

    // -------------------------------------------------------------------------------------
    // Resulting sample_fastq_metadata_ch structure for both branches:
    // [species, type, sample_id, [R1, R2], fasta, gtf, version]
    // type indicates if we split fastq or use full fastq
    // -------------------------------------------------------------------------------------

    // =====================================================================================
    // STEP 7: QUALITY CONTROL ON RAW (UNTRIMMED) READS
    // =====================================================================================
    // FastQC analyzes read quality before processing

    // Add read_type tag to channel: [sample_id, [fastqs], species, "raw"]
    fastqc_input_ch = sample_fastq_metadata_ch
        .map { species, type, sample_id, fastqs, fasta, gtf, version ->
            tuple(species, type, sample_id, fastqs, "raw")
        }

    FASTQC_RAW(fastqc_input_ch)
    if (params.stop_after == 'FASTQC_RAW') {
        log.info "Stopping pipeline after FASTQC_RAW as requested."
        return
    }

    // =====================================================================================
    // STEP 8: TRANSCRIPT QUANTIFICATION (SALMON)
    // =====================================================================================
    // Fast, alignment-free quantification
    // Runs in parallel with STAR (independent processes)

    // CRITICAL: Pre-join arguments to prevent cache invalidation
    // If params.SALMON_ARGS().join(' ') called inside process → hash changes → resume fails

    salmon_args           = params.SALMON_ARGS().join(' ')
    salmon_quant_input_ch = sample_fastq_metadata_ch
        .map { species, type, sample_id, fastqs, fasta, gtf, version ->
            tuple(species, type, sample_id, fastqs)
        }
        .combine(SALMON_INDEX.out.salmon_index_tuple, by: 0)

    SALMON_QUANT(salmon_quant_input_ch, salmon_args)

    if (params.stop_after == 'SALMON_QUANT') {
        log.info "Stopping pipeline after SALMON_QUANT as requested."
        return
    }

    // =====================================================================================
    // STEP 9: GENOME ALIGNMENT (STAR)
    // =====================================================================================
    // Splice-aware alignment, generates BAM files for visualization and QC

    // CRITICAL: Pre-join arguments to prevent cache invalidation
    // If params.STAR_ARGS().join(' ') called inside process → hash changes → resume fails

    star_args               = params.STAR_ARGS().join(' ')
    star_alignment_input_ch = sample_fastq_metadata_ch
        .map { species, type, sample_id, fastqs, fasta, gtf, version ->
            tuple(species, type, sample_id, fastqs)
        }
        .combine(STAR_INDEX.out.star_index_tuple, by: 0)

    STAR_ALIGN(star_alignment_input_ch, star_args)

    if (params.stop_after == 'STAR_ALIGN') {
        log.info "Stopping pipeline after STAR_ALIGN as requested."
        return
    }

    // Index BAM files and create subsampled versions for faster RSeQC
    sambamba_input_ch = STAR_ALIGN.out.star_results
        .map { species, type, sample_id, bam, gene_counts, sj_out, log ->
            tuple(species, type, sample_id, bam)
            }
    SAMBAMBA_PREP(sambamba_input_ch)

    if (params.stop_after == 'SAMBAMBA_PREP') {
        log.info "Stopping pipeline after SAMBAMBA_PREP as requested."
        return
    }

    // =====================================================================================
    // STEP 10: ALIGNMENT QUALITY CONTROL (RSEQC)
    // =====================================================================================
    // Comprehensive QC: read distribution, gene body coverage, junction analysis, etc.

    rseqc_input_ch = SAMBAMBA_PREP.out.bam_indexed          // [species, sample_id, bam, bai, 1M.bam, 1M.bai, read_len]
        .combine(RSEQC_BED.out.rseqc_bed_tuple, by: 0)      // [species, ref_bed, housekeeping_bed]
    RSEQC(rseqc_input_ch, mode_ch)

    if (params.stop_after == 'RSEQC') {
        log.info "Stopping pipeline after RSEQC as requested."
        return
    }

    // =====================================================================================
    // STEP 11: AGGREGATE QC REPORTS (MULTIQC)
    // =====================================================================================
    // Combines all QC outputs into single interactive HTML report

    // Group files by species first
    multiqc_input_ch = Channel.empty()
        .mix(FASTQC_RAW.out.fastqc_results.map     { species, type, zip, html                                -> tuple(species, type, [zip]) })
        .mix(SALMON_QUANT.out.salmon_quant_dir.map { species, type, sample_id, sample_dir                    -> tuple(species, type, [sample_dir]) })
        .mix(STAR_ALIGN.out.star_results.map       { species, type, sample_id, bam, gene_counts, sj_out, log -> tuple(species, type, [gene_counts, sj_out, log]) })
        .mix(RSEQC.out.rseqc_logs.map              { species, type, sample_id, logs                          -> tuple(species, type, [logs]) })
        .groupTuple(by: [0,1])                            // [ "Human", type, [zip], [sample_dir], [gene_counts, sj_out, log],[logs] ] and [ "Mouse", type, [zip], [sample_dir], [gene_counts, sj_out, log],[logs] ]
        .map { species, type, report_lists ->
            tuple(species, type, report_lists.flatten())    // [ "Human", type, [zip, sample_dir, gene_counts, sj_out, log, logs] ] and [ "Mouse", type, [zip, sample_dir, gene_counts, sj_out, log, logs] ]
    }

    MULTIQC(multiqc_input_ch)

    if (params.stop_after == 'MULTIQC') {
        log.info "Stopping pipeline after MULTIQC as requested."
        return
    }
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
//${base_dir}
//└──${expt}/
//   └──${project}/
//      ├── 01.FastQ/              //      Raw input FASTQs
//      │   ├── srr/
//      │   ├── raw/
//      │   └── trimmed/
//      ├── 02.FastQC/             //      QC on raw reads
//      │   ├── raw
//      │   └── trimmed
//      ├── 03.Salmon/                 //      Transcript quantification
//      │   ├── Sample1/               //      Full Salmon output
//      │   ├── Sample2/               //      Full Salmon output
//      │   └── quant_files/           //      Collected quant.sf files
//      ├── 04.STAR/                   //      Alignment outputs
//      │   ├── gene_counts/           //      ReadsPerGene.out.tab files
//      │   ├── splice_junction/       //      SJ.out.tab files
//      │   ├── alignment_stats/       //      Log.final.out files
//      │   └── bam/                   //      BAM + BAI files
//      ├── 05.RSEQC/                  //      Organized by analysis type
//      │   ├── 01_read_distribution/
//      │   ├── 02_inner_distance/
//      │   ├── 03_junction_annotation/
//      │   └── ... (09 subdirectories total)
//      ├── 06.MultiQC/
//      │   ├── Project_MultiQC_Report.html
//      │   └── Project_MultiQC_Report_data/
//      └── 07.Logs/                   //      All error logs + reports
//          ├── trace.txt
//          ├── report.html
//          └── timeline.html
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
