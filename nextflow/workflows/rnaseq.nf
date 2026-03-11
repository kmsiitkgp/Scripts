#!/usr/bin/env nextflow
nextflow.enable.dsl=2        // DSL2: enables modular workflows with explicit process blocks

// =========================================================================================
// IMPORT PROCESS MODULES
// =========================================================================================

include { RENAME_FASTQS }                from '../modules/rename_fastq.nf'
include { GENERATE_MD5 }                 from '../modules/generate_md5.nf'
include { VALIDATE_INPUT }               from '../modules/validate_input.nf'

include { FASTQC as FASTQC_RAW }         from '../modules/fastqc.nf'
include { FASTQC as FASTQC_TRIMMED }     from '../modules/fastqc.nf'           // Same process, different alias

include { FETCH_ENSEMBL_VERSIONS }       from '../modules/fetch_ensembl_versions.nf'
include { DOWNLOAD_ENSEMBL_REF }         from '../modules/download_ensembl_ref.nf'

include { XENGSORT_INDEX }               from '../modules/xengsort_index.nf'
include { XENGSORT_CLASSIFY }            from '../modules/xengsort_classify.nf'

include { STAR_INDEX }                   from '../modules/star_index.nf'
include { EXTRACT_GENTROME }             from '../modules/extract_gentrome.nf'
include { SALMON_INDEX }                 from '../modules/salmon_index.nf'
include { RSEQC_BED }                    from '../modules/rseqc_bed.nf'

include { CREATE_TX2GENE }               from '../modules/create_tx2gene.nf'
include { CREATE_TXI }                   from '../modules/create_txi.nf'
include { SALMON_QUANT }                 from '../modules/salmon_quant.nf'

include { STAR_ALIGN }                   from '../modules/star_align.nf'
include { SAMBAMBA_PREP }                from '../modules/sambamba_prep.nf'

include { RSEQC }                        from '../modules/rseqc.nf'
include { MULTIQC }                      from '../modules/multiqc.nf'

include { MERGE_STAR_COUNTS }            from '../modules/merge_star_counts.nf'
include { CREATE_DDS }                   from '../modules/create_dds.nf'

include { EXTRACT_DESEQ2_RESULTS }       from '../modules/extract_deseq2_results.nf'
include { ANALYZE_PATHWAYS }             from '../modules/analyze_pathways.nf'
include { PLOT_PATHWAYS }                from '../modules/plot_pathways.nf'

// =========================================================================================
// MAIN WORKFLOW
// =========================================================================================

workflow RNASEQ {

    // Declare channel variables upfront; populated conditionally below
    def reference_ch          = Channel.empty()
    def create_tx2gene_out_ch = Channel.empty()
    def salmon_quant_out_ch   = Channel.empty()
    def star_align_out_ch     = Channel.empty()

    // =====================================================================================
    // MODE 1: BULK RNA-SEQ RUN
    // =====================================================================================

    if (params.run_mode == 'bulk' || params.run_mode == 'bulk+de') {

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
        // STEP 1: RENAME FASTQ FILES
        // =====================================================================================

        map_file_path = params.map_file ? file(params.map_file) : null

        if (map_file_path && map_file_path.exists()) {

            map_ch       = Channel.value(map_file_path)
            srr_fastq_ch = Channel.fromPath("${params.proj_dir()}/downloads/*.{fastq,fq}.gz")

            // .collect() sends ALL files to ONE process instance (not one per file)
            RENAME_FASTQS(srr_fastq_ch.collect(), map_ch)

            // .flatten() breaks the output list back into individual file items
            raw_fastq_ch = RENAME_FASTQS.out.renamed_fastqs.flatten()

        } else {
            log.info "No map file found. Passing raw FASTQs directly to validation."
            raw_fastq_ch = Channel.fromPath("${params.read_dir}/*.{fastq,fq}.gz")
        }

        if (params.stop_after == 'RENAME_FASTQS') {
            log.info "Stopping pipeline after RENAME_FASTQS as requested."
            return
        }

        // =====================================================================================
        // STEP 2: GENERATE MD5
        // =====================================================================================

        GENERATE_MD5(raw_fastq_ch)

        // Merge all per-file MD5s into one manifest file
        manifest_ch = GENERATE_MD5.out.md5_file.collectFile(
            name:     'manifest_md5.txt',
            newLine:  true,
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
        // STEP 3: VALIDATE INPUT FASTQ FILES
        // =====================================================================================
        // Checks naming conventions; emits grouped samples and sequencing mode

        VALIDATE_INPUT(raw_fastq_ch.collect())
        mode_ch         = VALIDATE_INPUT.out.mode.collect()      // collected → Value Channel so RSeQC can reuse it
        sample_fastq_ch = VALIDATE_INPUT.out.grouped_samples_ch  // [sample_id, [R1, R2]]
        sample_ch       = VALIDATE_INPUT.out.sample_names_ch     // [sample_id]

        // NOTE on Queue vs Value Channel:
        // Queue Channel (no .collect()): 10 samples compete for 1 index → 9 hang forever.
        // Value Channel (.collect()):    all 10 samples read the same index simultaneously.

        if (params.stop_after == 'VALIDATE_INPUT') {
            log.info "Stopping pipeline after VALIDATE_INPUT as requested."
            return
        }

        // =================================================================================
        // STEP 4: FETCH ENSEMBL REFERENCE METADATA
        // =================================================================================

        // Xenograft: need Human + Mouse metadata individually, plus a synthetic combined entry
        if (params.species == 'Xenograft') {
            target_species_list = ['Human', 'Mouse', 'Xenograft']
        } else {
            target_species_list = [params.species]
        }

        species_list_ch = Channel.fromList(target_species_list)
        FETCH_ENSEMBL_VERSIONS(species_list_ch)

        fetch_ensembl_versions_out_ch = FETCH_ENSEMBL_VERSIONS.out.meta
            .splitText()
            .map { it.trim() }
            .filter { it }
            .map { line ->
                def parts = line.trim().split('\\|').collect { it.trim() }
                // splitText() emits a Java List, not a Nextflow Tuple.
                // Index-based access .map { it[3] } works fine, but named deconstruction
                // .map { a, b -> } triggers 'Invalid method invocation' on nested lists.
                // Use tuple() constructor to produce a proper Nextflow Tuple.
                if (parts.size() >= 6) {
                    return tuple(parts[0], parts[1], parts[2], parts[3], parts[4], parts[5])
                }
            }
            .filter { it != null }

        // [Species, FA_Name, GTF_Name, Version, Ensembl_Assembly, Ensembl_Release]
        // ["Human",     "Homo_sapiens.GRCh38.dna.primary_assembly.fa", "Homo_sapiens.GRCh38.115.gtf", "GRCh38.115",            "GRCh38", "115"]
        // ["Mouse",     "Mus_musculus.GRCm39.dna.primary_assembly.fa", "Mus_musculus.GRCm39.115.gtf", "GRCm39.115",            "GRCm39", "115"]
        // ["Xenograft", "Homo_sapiens.GRCh38.dna.primary_assembly.fa", "Homo_sapiens.GRCh38.115.gtf", "GRCh38.115_GRCm39.115", "GRCh38", "115"]  ← composite version; no real files on Ensembl

        // =================================================================================
        // STEP 5: DOWNLOAD REFERENCE FILES
        // =================================================================================

        // BUG FIX: Xenograft has no FA/GTF on Ensembl — filter it out before downloading.
        // The Xenograft reference entry is synthetic (Human FA+GTF paths, composite version tag);
        // it is reconstructed below by combining the downloaded Human paths with Xenograft metadata.
        download_ensembl_ref_in_ch = fetch_ensembl_versions_out_ch
            .filter { it[0] != "Xenograft" }

        // [Species, FA_Name, GTF_Name, Version, Ensembl_Assembly, Ensembl_Release]
        // ["Human", "Homo_sapiens.GRCh38.dna.primary_assembly.fa", "Homo_sapiens.GRCh38.115.gtf", "GRCh38.115", "GRCh38", "115"]
        // ["Mouse", "Mus_musculus.GRCm39.dna.primary_assembly.fa", "Mus_musculus.GRCm39.115.gtf", "GRCm39.115", "GRCm39", "115"]

        DOWNLOAD_ENSEMBL_REF(download_ensembl_ref_in_ch)
        reference_ch = DOWNLOAD_ENSEMBL_REF.out.ref_tuple

        // [Species, FA_Path, GTF_Path, Version, Ensembl_Assembly, Ensembl_Release]
        // ["Human", "/ref/GRCh38.115/Homo_sapiens.GRCh38.dna.primary_assembly.fa", "/ref/GRCh38.115/Homo_sapiens.GRCh38.115.gtf", "GRCh38.115", "GRCh38", "115"]
        // ["Mouse", "/ref/GRCm39.115/Mus_musculus.GRCm39.dna.primary_assembly.fa", "/ref/GRCm39.115/Mus_musculus.GRCm39.115.gtf", "GRCm39.115", "GRCm39", "115"]

        if (params.species == 'Xenograft') {

            // Build synthetic Xenograft row: Human FA+GTF paths, but with the composite version tag.
            // This lets STAR/Salmon index steps produce a Xenograft-labelled index directory
            // that maps unseparated reads against the human reference.

            xeno_meta_ch = fetch_ensembl_versions_out_ch
                .filter { it[0] == "Xenograft" }
                .map { species, fa, gtf, version, assembly, release -> tuple(species, version, assembly, release) }
            // ["Xenograft", "GRCh38.115_GRCm39.115", "GRCh38", "115"]

            human_ref_ch = reference_ch
                .filter { it[0] == "Human" }
                .map { species, fa, gtf, version, assembly, release -> tuple(fa, gtf) }
            // ["/ref/GRCh38.115/...fa", "/ref/GRCh38.115/...gtf"]

            xeno_ref_ch = xeno_meta_ch
                .combine(human_ref_ch)
                .map { species, version, assembly, release, fa, gtf -> tuple(species, fa, gtf, version, assembly, release) }
            // ["Xenograft", "/ref/GRCh38.115/...fa", "/ref/GRCh38.115/...gtf", "GRCh38.115_GRCm39.115", "GRCh38", "115"]

            // Append Xenograft row so all downstream index steps receive all 3 species
            reference_ch = reference_ch.concat(xeno_ref_ch)
            // ["Human",     "/ref/GRCh38.115/...fa", "/ref/GRCh38.115/...gtf", "GRCh38.115",            "GRCh38", "115"]
            // ["Mouse",     "/ref/GRCm39.115/...fa", "/ref/GRCm39.115/...gtf", "GRCm39.115",            "GRCm39", "115"]
            // ["Xenograft", "/ref/GRCh38.115/...fa", "/ref/GRCh38.115/...gtf", "GRCh38.115_GRCm39.115", "GRCh38", "115"]
        }

        // =====================================================================================
        // STEP 6: BUILD REFERENCE INDEXES
        // =====================================================================================
        // Each index process receives all 3 rows (Human, Mouse, Xenograft) and runs in parallel

        STAR_INDEX(reference_ch)
        if (params.stop_after == 'STAR_INDEX') {
            log.info "Stopping pipeline after STAR_INDEX as requested."
            return
        }

        EXTRACT_GENTROME(reference_ch)
        decoy_gentrome_ch = EXTRACT_GENTROME.out.decoy_gentrome_tuple

        SALMON_INDEX(decoy_gentrome_ch)
        if (params.stop_after == 'SALMON_INDEX') {
            log.info "Stopping pipeline after SALMON_INDEX as requested."
            return
        }

        RSEQC_BED(reference_ch)
        rseqc_bed_out_ch = RSEQC_BED.out.rseqc_bed_tuple

        // [Species, ref_bed, housekeeping_bed]
        // ["Human",     "/bed/GRCh38_ref.bed", "/bed/GRCh38_hk.bed"]
        // ["Mouse",     "/bed/GRCm39_ref.bed",  "/bed/GRCm39_hk.bed"]
        // ["Xenograft", "/bed/GRCh38_ref.bed", "/bed/GRCh38_hk.bed"]  ← same beds as Human

        if (params.stop_after == 'RSEQC_BED') {
            log.info "Stopping pipeline after RSEQC_BED as requested."
            return
        }

        // =================================================================================
        // STEP 7: ASSIGN SPECIES LABELS & PREPARE MAPPING CHANNEL
        // =================================================================================

        if (params.species == 'Xenograft') {

            // --- XENOGRAFT BRANCH ---
            // XenGSort separates mixed-species reads into Graft (Human) and Host (Mouse) bins.
            // The full (unseparated) reads are also kept and mapped under the "Xenograft" label.

            human_fasta_ch = reference_ch
                .filter { it[0] == "Human" }
                .map { species, fa, gtf, version, assembly, release -> fa }
            // "/ref/GRCh38.115/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

            mouse_fasta_ch = reference_ch
                .filter { it[0] == "Mouse" }
                .map { species, fa, gtf, version, assembly, release -> fa }
            // "/ref/GRCm39.115/Mus_musculus.GRCm39.dna.primary_assembly.fa"

            // .first() not .collect() — .collect() wraps the string in a list,
            // so XENGSORT_INDEX would receive ["GRCh38.115_GRCm39.115"] instead of "GRCh38.115_GRCm39.115"
            xeno_version_ch = reference_ch
                .filter { it[0] == "Xenograft" }
                .map { species, fa, gtf, version, assembly, release -> version }
                .first()
            // "GRCh38.115_GRCm39.115"

            XENGSORT_INDEX(human_fasta_ch, mouse_fasta_ch, xeno_version_ch)
            xengsort_index_out_ch = XENGSORT_INDEX.out.xengsort_index_dir.collect()

            XENGSORT_CLASSIFY(sample_fastq_ch, xengsort_index_out_ch, xeno_version_ch)

            graft_fastq_ch = XENGSORT_CLASSIFY.out.graft_fastqs
                .map { sample_id, fastqs -> tuple("Human", sample_id, fastqs) }
            // ["Human", "H_S1", [H_S1_graft_R1.fq, H_S1_graft_R2.fq]]
            // ["Human", "H_S2", [H_S2_graft_R1.fq, H_S2_graft_R2.fq]]

            host_fastq_ch = XENGSORT_CLASSIFY.out.host_fastqs
                .map { sample_id, fastqs -> tuple("Mouse", sample_id, fastqs) }
            // ["Mouse", "M_S1", [M_S1_host_R1.fq, M_S1_host_R2.fq]]
            // ["Mouse", "M_S2", [M_S2_host_R1.fq, M_S2_host_R2.fq]]

            full_fastq_ch = sample_fastq_ch
                .map { sample_id, fastqs -> tuple("Xenograft", sample_id, fastqs) }
            // ["Xenograft", "S1", [S1_R1.fq, S1_R2.fq]]
            // ["Xenograft", "S2", [S2_R1.fq, S2_R2.fq]]

            // Merge all three labeled streams into one channel
            all_fastq_ch = graft_fastq_ch
                .concat(host_fastq_ch)
                .concat(full_fastq_ch)
            // ["Human",     "H_S1", [H_S1_graft_R1.fq, H_S1_graft_R2.fq]]
            // ["Human",     "H_S2", [H_S2_graft_R1.fq, H_S2_graft_R2.fq]]
            // ["Mouse",     "M_S1", [M_S1_host_R1.fq,  M_S1_host_R2.fq]]
            // ["Mouse",     "M_S2", [M_S2_host_R1.fq,  M_S2_host_R2.fq]]
            // ["Xenograft", "S1",   [S1_R1.fq, S1_R2.fq]]
            // ["Xenograft", "S2",   [S2_R1.fq, S2_R2.fq]]

            // Join each labeled sample to its matching reference row by Species key
            sample_fastq_metadata_ch = all_fastq_ch
                .combine(reference_ch, by: 0)

        } else {

            // --- STANDARD BRANCH ---
            // Single species: attach the one reference row to every sample
            sample_fastq_metadata_ch = sample_fastq_ch
                .combine(reference_ch)
                .map { sample_id, fastqs, species, fa, gtf, version, assembly, release ->
                    tuple(species, sample_id, fastqs, fa, gtf, version, assembly, release) }
        }

        // [Species, sample_id, [R1,R2], FA_path, GTF_path, Version, Assembly, Release]
        //
        // Standard (Human):
        // ["Human", "S1", [S1_R1.fq, S1_R2.fq], "/ref/GRCh38.115/...fa", "/ref/GRCh38.115/...gtf", "GRCh38.115", "GRCh38", "115"]
        // ["Human", "S2", [S2_R1.fq, S2_R2.fq], "/ref/GRCh38.115/...fa", "/ref/GRCh38.115/...gtf", "GRCh38.115", "GRCh38", "115"]
        //
        // Xenograft:
        // ["Human",     "H_S1", [H_S1_graft_R1.fq, ...], "/ref/GRCh38.115/...fa", "/ref/GRCh38.115/...gtf", "GRCh38.115",            "GRCh38", "115"]
        // ["Human",     "H_S2", [H_S2_graft_R1.fq, ...], "/ref/GRCh38.115/...fa", "/ref/GRCh38.115/...gtf", "GRCh38.115",            "GRCh38", "115"]
        // ["Mouse",     "M_S1", [M_S1_host_R1.fq,  ...], "/ref/GRCm39.115/...fa", "/ref/GRCm39.115/...gtf", "GRCm39.115",            "GRCm39", "115"]
        // ["Mouse",     "M_S2", [M_S2_host_R1.fq,  ...], "/ref/GRCm39.115/...fa", "/ref/GRCm39.115/...gtf", "GRCm39.115",            "GRCm39", "115"]
        // ["Xenograft", "S1",   [S1_R1.fq, ...],         "/ref/GRCh38.115/...fa", "/ref/GRCh38.115/...gtf", "GRCh38.115_GRCm39.115", "GRCh38", "115"]
        // ["Xenograft", "S2",   [S2_R1.fq, ...],         "/ref/GRCh38.115/...fa", "/ref/GRCh38.115/...gtf", "GRCh38.115_GRCm39.115", "GRCh38", "115"]

        // =====================================================================================
        // STEP 8: FASTQC ON RAW READS
        // =====================================================================================

        fastqc_in_ch = sample_fastq_metadata_ch
            .map { species, sample_id, fastqs, fa, gtf, version, assembly, release ->
                tuple(species, sample_id, fastqs, "raw")
            }

        // [Species, sample_id, [R1,R2], fastq_type]
        // ["Human",     "H_S1", [H_S1_graft_R1.fq, H_S1_graft_R2.fq], "raw"]
        // ["Mouse",     "M_S1", [M_S1_host_R1.fq,  M_S1_host_R2.fq],  "raw"]
        // ["Xenograft", "S1",   [S1_R1.fq, S1_R2.fq],                  "raw"]

        FASTQC_RAW(fastqc_in_ch)
        if (params.stop_after == 'FASTQC_RAW') {
            log.info "Stopping pipeline after FASTQC_RAW as requested."
            return
        }

        // =====================================================================================
        // STEP 9: TRANSCRIPT QUANTIFICATION (SALMON)
        // =====================================================================================

        // Pre-join args string here — calling .join() inside the process script
        // generates a new hash every run and breaks -resume cache
        salmon_args = (params.SALMON_ARGS() instanceof List) ? params.SALMON_ARGS().join(' ') : params.SALMON_ARGS()

        salmon_quant_in_ch = sample_fastq_metadata_ch
            .map { species, sample_id, fastqs, fa, gtf, version, assembly, release ->
                tuple(species, sample_id, fastqs)
            }
            .combine(SALMON_INDEX.out.salmon_index_tuple, by: 0)

        // [Species, sample_id, [R1,R2], salmon_index_dir]  — joined on Species
        // ["Human",     "H_S1", [H_S1_graft_R1.fq, H_S1_graft_R2.fq], "/idx/salmon_GRCh38.115/"]
        // ["Mouse",     "M_S1", [M_S1_host_R1.fq,  M_S1_host_R2.fq],  "/idx/salmon_GRCm39.115/"]
        // ["Xenograft", "S1",   [S1_R1.fq, S1_R2.fq],                  "/idx/salmon_GRCh38.115_GRCm39.115/"]

        SALMON_QUANT(salmon_quant_in_ch, salmon_args)
        salmon_quant_out_ch = SALMON_QUANT.out.salmon_quant_dir

        // [Species, sample_id, quant_dir]
        // ["Human",     "H_S1", "/03.Salmon/H_S1/"]
        // ["Human",     "H_S2", "/03.Salmon/H_S2/"]
        // ["Mouse",     "M_S1", "/03.Salmon/M_S1/"]
        // ["Mouse",     "M_S2", "/03.Salmon/M_S2/"]
        // ["Xenograft", "S1",   "/03.Salmon/S1/"]
        // ["Xenograft", "S2",   "/03.Salmon/S2/"]

        if (params.stop_after == 'SALMON_QUANT') {
            log.info "Stopping pipeline after SALMON_QUANT as requested."
            return
        }

        // =====================================================================================
        // STEP 10: GENOME ALIGNMENT (STAR)
        // =====================================================================================

        // Pre-join args string here — same reason as salmon_args above
        star_args = (params.STAR_ARGS() instanceof List) ? params.STAR_ARGS().join(' ') : params.STAR_ARGS()

        star_align_in_ch = sample_fastq_metadata_ch
            .map { species, sample_id, fastqs, fa, gtf, version, assembly, release ->
                tuple(species, sample_id, fastqs)
            }
            .combine(STAR_INDEX.out.star_index_tuple, by: 0)

        // [Species, sample_id, [R1,R2], star_index_dir]  — joined on Species
        // ["Human",     "H_S1", [H_S1_graft_R1.fq, H_S1_graft_R2.fq], "/idx/star_GRCh38.115/"]
        // ["Mouse",     "M_S1", [M_S1_host_R1.fq,  M_S1_host_R2.fq],  "/idx/star_GRCm39.115/"]
        // ["Xenograft", "S1",   [S1_R1.fq, S1_R2.fq],                  "/idx/star_GRCh38.115_GRCm39.115/"]

        STAR_ALIGN(star_align_in_ch, star_args)
        star_align_out_ch = STAR_ALIGN.out.star_results

        // [Species, sample_id, bam, gene_counts, sj_out, log]
        // ["Human",     "H_S1", "H_S1.bam", "H_S1_ReadsPerGene.tab", "H_S1_SJ.out.tab", "H_S1_Log.final.out"]
        // ["Mouse",     "M_S1", "M_S1.bam", "M_S1_ReadsPerGene.tab", "M_S1_SJ.out.tab", "M_S1_Log.final.out"]
        // ["Xenograft", "S1",   "S1.bam",   "S1_ReadsPerGene.tab",   "S1_SJ.out.tab",   "S1_Log.final.out"]

        if (params.stop_after == 'STAR_ALIGN') {
            log.info "Stopping pipeline after STAR_ALIGN as requested."
            return
        }

        // Sort + index BAM; also creates a 1M-read subsample for faster RSeQC
        sambamba_prep_in_ch = star_align_out_ch
            .map { species, sample_id, bam, gene_counts, sj_out, log ->
                tuple(species, sample_id, bam)
            }

        // [Species, sample_id, bam]
        // ["Human",     "H_S1", "H_S1.bam"]
        // ["Mouse",     "M_S1", "M_S1.bam"]
        // ["Xenograft", "S1",   "S1.bam"]

        SAMBAMBA_PREP(sambamba_prep_in_ch)
        sambamba_prep_out_ch = SAMBAMBA_PREP.out.bam_indexed

        // [Species, sample_id, bam, bai, 1M.bam, 1M.bai, read_len]
        // ["Human",     "H_S1", "H_S1.bam", "H_S1.bai", "H_S1_1M.bam", "H_S1_1M.bai", 150]
        // ["Mouse",     "M_S1", "M_S1.bam", "M_S1.bai", "M_S1_1M.bam", "M_S1_1M.bai", 150]
        // ["Xenograft", "S1",   "S1.bam",   "S1.bai",   "S1_1M.bam",   "S1_1M.bai",   150]

        if (params.stop_after == 'SAMBAMBA_PREP') {
            log.info "Stopping pipeline after SAMBAMBA_PREP as requested."
            return
        }

        // =====================================================================================
        // STEP 11: ALIGNMENT QC (RSEQC)
        // =====================================================================================

        rseqc_in_ch = sambamba_prep_out_ch
            // [Species, sample_id, bam, bai, 1M.bam, 1M.bai, read_len]
            .combine(rseqc_bed_out_ch, by: 0)
            // attaches [ref_bed, housekeeping_bed] matched on Species

        // [Species, sample_id, bam, bai, 1M.bam, 1M.bai, read_len, ref_bed, housekeeping_bed]
        // ["Human",     "H_S1", "H_S1.bam", "H_S1.bai", "H_S1_1M.bam", "H_S1_1M.bai", 150, "/bed/GRCh38_ref.bed", "/bed/GRCh38_hk.bed"]
        // ["Mouse",     "M_S1", "M_S1.bam", "M_S1.bai", "M_S1_1M.bam", "M_S1_1M.bai", 150, "/bed/GRCm39_ref.bed",  "/bed/GRCm39_hk.bed"]
        // ["Xenograft", "S1",   "S1.bam",   "S1.bai",   "S1_1M.bam",   "S1_1M.bai",   150, "/bed/GRCh38_ref.bed", "/bed/GRCh38_hk.bed"]

        RSEQC(rseqc_in_ch, mode_ch)
        if (params.stop_after == 'RSEQC') {
            log.info "Stopping pipeline after RSEQC as requested."
            return
        }

        // =====================================================================================
        // STEP 12: AGGREGATE QC REPORT (MULTIQC)
        // =====================================================================================
        // Mix all QC outputs → group by Species → flatten into one file list per species

        multiqc_in_ch = Channel.empty()
            .mix(FASTQC_RAW.out.fastqc_results.map     { species, zip, html                                -> tuple(species, [zip]) })
            .mix(SALMON_QUANT.out.salmon_quant_dir.map { species, sample_id, sample_dir                    -> tuple(species, [sample_dir]) })
            .mix(STAR_ALIGN.out.star_results.map       { species, sample_id, bam, gene_counts, sj_out, log -> tuple(species, [gene_counts, sj_out, log]) })
            .mix(RSEQC.out.rseqc_logs.map              { species, sample_id, logs                          -> tuple(species, [logs]) })
            .groupTuple(by: 0)
            .map { species, report_lists ->
                tuple(species, report_lists.flatten())
            }
            // After groupTuple: one row per Species; each file field is now a list-of-lists
            // ["Human",     [[zip1],[zip2]], [[dir1],[dir2]], [[gc1,sj1,log1],[gc2,sj2,log2]], [[logs1],[logs2]]]
            // ["Mouse",     [[zip1],[zip2]], [[dir1],[dir2]], [[gc1,sj1,log1],[gc2,sj2,log2]], [[logs1],[logs2]]]
            // ["Xenograft", [[zip1]],        [[dir1]],        [[gc1,sj1,log1]],                [[logs1]]]
            //
            // Splat (...report_lists) captures all fields after species;
            // flatten() collapses the nested list-of-lists into a single flat file list

            // After flatten: clean flat list ready for MultiQC
            // ["Human",     [zip1, zip2, dir1, dir2, gc1, sj1, log1, gc2, sj2, log2, logs1, logs2]]
            // ["Mouse",     [zip1, zip2, dir1, dir2, gc1, sj1, log1, gc2, sj2, log2, logs1, logs2]]
            // ["Xenograft", [zip1, dir1, gc1, sj1, log1, logs1]]

        MULTIQC(multiqc_in_ch)
        if (params.stop_after == 'MULTIQC') {
            log.info "Stopping pipeline after MULTIQC as requested."
            return
        }
    }

    // =====================================================================================
    // MODE 2: DIFFERENTIAL EXPRESSION RUN
    // =====================================================================================

    if (params.run_mode == 'de' || params.run_mode == 'bulk+de') {

        // 'de' mode:      user provides pre-computed counts/quant dirs + GTF.
        // 'bulk+de' mode: reference_ch / star_align_out_ch / salmon_quant_out_ch already populated above.
        if (params.run_mode == 'de') {

            // For Xenograft in 'de' mode: run DE on separated Human reads only
            def target_species = (params.species == 'Xenograft') ? 'Human' : params.species

            def gtf_file = params.gtf ? file(params.gtf) : null
            def fa_file  = null
            reference_ch = Channel.of(tuple(target_species, fa_file, gtf_file, params.genome_version, params.assembly, params.release))
            // ["Human", null, "/ref/Homo_sapiens.GRCh38.115.gtf", "GRCh38.115", "GRCh38", "115"]

            // User-supplied STAR count files (glob e.g. "counts/*.tab")
            star_align_out_ch = Channel.fromPath(params.gene_counts)
                .map { file -> tuple(target_species, file.baseName, null, file, null, null) }
            // [Species, sample_id, bam=null, gene_counts, sj_out=null, log=null]
            // ["Human", "H_S1", null, "/counts/H_S1_ReadsPerGene.tab", null, null]
            // ["Human", "H_S2", null, "/counts/H_S2_ReadsPerGene.tab", null, null]

            // User-supplied Salmon output dirs (glob e.g. "03.Salmon/*/")
            salmon_quant_out_ch = Channel.fromPath(params.salmon_dir)
                .map { sample_dir -> tuple(target_species, sample_dir.name, sample_dir) }
            // [Species, sample_id, quant_dir]
            // ["Human", "H_S1", "/03.Salmon/H_S1/"]
            // ["Human", "H_S2", "/03.Salmon/H_S2/"]
        }

        // =====================================================================================
        // STEP DE-1: CREATE TX2GENE MAPPING
        // =====================================================================================

        CREATE_TX2GENE(reference_ch)
        create_tx2gene_out_ch = CREATE_TX2GENE.out.tx2gene_tuple

        // [Species, tx2gene_csv]
        // ["Human", "/path/human_tx2gene.csv"]
        // ["Mouse", "/path/mouse_tx2gene.csv"]  ← only present in bulk+de Xenograft runs

        if (params.stop_after == 'CREATE_TX2GENE') {
            log.info "Stopping pipeline after CREATE_TX2GENE as requested."
            return
        }

        // =====================================================================================
        // STEP DE-2: MERGE STAR COUNTS
        // =====================================================================================

        merge_star_counts_in_ch = star_align_out_ch
            .map { species, sample_id, bam, gene_counts, sj_out, log -> tuple(species, gene_counts) }
            // [Species, gene_counts]
            // ["Human", "/counts/H_S1_ReadsPerGene.tab"]
            // ["Human", "/counts/H_S2_ReadsPerGene.tab"]
            // ["Mouse", "/counts/M_S1_ReadsPerGene.tab"]
            // ["Mouse", "/counts/M_S2_ReadsPerGene.tab"]
            .groupTuple(by: 0)
            // [Species, [gene_counts...]]
            // ["Human", ["/counts/H_S1_ReadsPerGene.tab", "/counts/H_S2_ReadsPerGene.tab"]]
            // ["Mouse", ["/counts/M_S1_ReadsPerGene.tab", "/counts/M_S2_ReadsPerGene.tab"]]

        MERGE_STAR_COUNTS(merge_star_counts_in_ch)

        // =====================================================================================
        // STEP DE-3: CREATE TXI OBJECT (tximeta/tximport)
        // =====================================================================================

        create_txi_in_ch = salmon_quant_out_ch
            .map { species, sample_id, sample_dir -> tuple(species, sample_dir) }
            // ["Human", "/03.Salmon/H_S1/"]
            // ["Human", "/03.Salmon/H_S2/"]
            // ["Mouse", "/03.Salmon/M_S1/"]
            // ["Mouse", "/03.Salmon/M_S2/"]
            .groupTuple(by: 0)
            // ["Human", ["/03.Salmon/H_S1/", "/03.Salmon/H_S2/"]]
            // ["Mouse", ["/03.Salmon/M_S1/", "/03.Salmon/M_S2/"]]
            .combine(create_tx2gene_out_ch, by: 0)
            // [Species, [quant_dirs...], tx2gene_csv]
            // ["Human", ["/03.Salmon/H_S1/", "/03.Salmon/H_S2/"], "/path/human_tx2gene.csv"]
            // ["Mouse", ["/03.Salmon/M_S1/", "/03.Salmon/M_S2/"], "/path/mouse_tx2gene.csv"]

        CREATE_TXI(create_txi_in_ch)
        create_txi_out_ch = CREATE_TXI.out.txi

        // [Species, txi_rds, tx2gene_csv]
        // ["Human", "/path/human_txi.rds", "/path/human_tx2gene.csv"]
        // ["Mouse", "/path/mouse_txi.rds", "/path/mouse_tx2gene.csv"]

        // =====================================================================================
        // STEP DE-4: CREATE DESEQ2 DATASET
        // =====================================================================================

        metadata_ch = Channel.fromPath(params.metadata_file, checkIfExists: true).collect()

        CREATE_DDS(create_txi_out_ch, metadata_ch, params.deseq2.design)
        create_dds_out_ch = CREATE_DDS.out.dds

        // [Species, dds_rds, tx2gene_csv]
        // ["Human", "/path/human_dds.rds", "/path/human_tx2gene.csv"]
        // ["Mouse", "/path/mouse_dds.rds", "/path/mouse_tx2gene.csv"]

        // =====================================================================================
        // STEP DE-5: EXTRACT DESEQ2 RESULTS
        // =====================================================================================

        extract_deseq2_results_in_ch = create_dds_out_ch
            .combine(Channel.fromList(params.deseq2.contrasts))

        // [Species, dds_rds, tx2gene_csv, contrast]
        // ["Human", "/path/human_dds.rds", "/path/human_tx2gene.csv", ["condition","Treated","Control"]]
        // ["Mouse", "/path/mouse_dds.rds", "/path/mouse_tx2gene.csv", ["condition","Treated","Control"]]

        deseq2_ch = Channel.value(params.deseq2)
        EXTRACT_DESEQ2_RESULTS(extract_deseq2_results_in_ch, deseq2_ch)

        gsea_ch = Channel.fromList([params.gsea]).collect()
        gmt_dir_ch = Channel.from(["Human", "Mouse"])
            .map { species -> tuple(species, file("${params.gmt_dir}/${species}", checkIfExists: true)) }

        analyze_pathways_in_ch = EXTRACT_DESEQ2_RESULTS.out.deg_dir
           .map { species, contrast, dir -> tuple(species, contrast, file("${dir}/DEGs.xlsx", checkIfExists: true)) }
           .combine(gmt_dir_ch, by: 0)
        ANALYZE_PATHWAYS(analyze_pathways_in_ch, gsea_ch)

        heatmap_ch = Channel.fromList([params.heatmap]).collect()
        plot_pathways_in_ch = ANALYZE_PATHWAYS.out.pathway_dir
            .map { species, contrast, dir -> tuple(species, contrast, file("${dir}/Pathways.xlsx", checkIfExists: true), "NULL") }
        PLOT_PATHWAYS(plot_pathways_in_ch, metadata_ch, heatmap_ch, gsea_ch )

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
// ${base_dir}/
// └── ${expt}/
//     └── ${project}/
//         ├── 01.FastQ/
//         │   ├── srr/
//         │   ├── raw/
//         │   └── trimmed/
//         ├── 02.FastQC/
//         │   ├── raw/
//         │   └── trimmed/
//         ├── 03.Salmon/
//         │   ├── H_S1/              ← full Salmon output per sample
//         │   ├── H_S2/
//         │   └── quant_files/       ← collected quant.sf files
//         ├── 04.STAR/
//         │   ├── gene_counts/       ← ReadsPerGene.out.tab files
//         │   ├── splice_junction/   ← SJ.out.tab files
//         │   ├── alignment_stats/   ← Log.final.out files
//         │   └── bam/               ← BAM + BAI files
//         ├── 05.RSEQC/
//         │   ├── 01_read_distribution/
//         │   ├── 02_inner_distance/
//         │   ├── 03_junction_annotation/
//         │   └── ...  (09 subdirectories total)
//         ├── 06.MultiQC/
//         │   ├── Project_MultiQC_Report.html
//         │   └── Project_MultiQC_Report_data/
//         └── 07.Logs/
//             ├── trace.txt
//             ├── report.html
//             └── timeline.html
//
// QUALITY CONTROL THRESHOLDS (check MultiQC report):
//   Mapping rate:      >70% good  |  60–70% acceptable  |  <60% investigate
//   Duplication:       <50% good  |  50–70% acceptable  |  >70% low complexity
//   Read distribution: >50% CDS exons, <10% introns, <5% intergenic
//   3' bias:           <3 good    |  3–5 acceptable      |  >5 degraded RNA
//
// RESUMING FAILED RUNS:
//   Nextflow caches completed processes in work/
//   To resume: bash run_nextflow.sh  (script includes -resume flag)
//   Cache invalidated if: process code changes, input files change, or work/ deleted
//
// CLEANING UP:
//   rm -rf ${work_dir}/*   ← frees disk space; WARNING: cannot resume after this
//
// COMMON ISSUES:
//   "No such file or directory"
//     → Check paths in project_info.yaml
//     → Verify Singularity bind mounts (run_nextflow.sh prints mappings)
//
//   Out of memory
//     → Increase memory in process labels (nextflow.config)
//     → Reduce parallel jobs with maxForks in nextflow.config
//
//   Process won't resume
//     → Ensure params.FUNCTION().join(' ') is called in workflow block, not inside process
//     → Verify work/ directory not deleted
//     → Review .nextflow.log for cache invalidation reasons
//
// DOWNSTREAM ANALYSIS:
//   DE:            STAR gene counts or Salmon quant.sf → tximport → DESeq2/edgeR
//   Visualization: BAM files → IGV; coverage plots → RSeQC or deepTools
//   Splice:        STAR SJ.out.tab → rMATS or LeafCutter
//
// =========================================================================================
