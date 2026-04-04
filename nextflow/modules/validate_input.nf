// =============================================================================
// WORKFLOW: VALIDATE_INPUT
// =============================================================================
// Purpose: Validates FASTQ naming conventions, detects SE/PE mode, and
//          creates organized sample channels for all downstream processes
//
// What it does:
//   1. Validates FASTQ filenames against expected naming pattern
//   2. Detects single-end vs paired-end sequencing
//   3. Enforces consistent file extension (.fq.gz or .fastq.gz, not mixed)
//   4. Pairs R1/R2 files correctly and verifies each R2 exists
//   5. Emits grouped sample tuples and metadata channels
//
// Why fail-fast validation here?
//   - Catches naming errors before wasting hours of compute time
//   - Provides clear, actionable error messages instead of cryptic downstream failures
// =============================================================================

workflow VALIDATE_INPUT {

    take:
    fastq_files_ch    // Channel of Path objects (individual FASTQ files)

    main:
    log.info "==> Validating file naming conventions"

    // .collect() on a channel acts as a synchronization barrier:
    // waits for all upstream items, then emits them as a single List.
    // .map{} on the collected List runs synchronously in Groovy memory — no channel semantics.
    results_ch = fastq_files_ch.collect().map { fastq_list ->

        def fastq_files = fastq_list.sort { it.name }

        // =================================================================================
        // 1. DEFINE VALID NAMING PATTERN
        // =================================================================================

        def VALID_PATTERN
        if (params.expt == "RNASeq") {
            VALID_PATTERN = ~/.*(_Tumor|_Normal)?.*(_[Rr][12]).*\.f(q|astq)\.gz/
        } else if (params.expt == "scRNASeq") {
            VALID_PATTERN = ~/^([A-Za-z0-9-]+)_S\d+_L\d{3}_(R[12]|I[12])_001\.fastq\.gz$/
        } else {
            error "❌ Unknown experiment type: ${params.expt}. Options: 'RNASeq', 'scRNASeq'"
        }

        def valid_files   = fastq_files.findAll { it.name ==~ VALID_PATTERN }
        def invalid_files = fastq_files.findAll { !(it.name ==~ VALID_PATTERN) }

        // =================================================================================
        // 2. DETECT FILE EXTENSION
        // =================================================================================

        def fq_gz_files    = valid_files.findAll { it.name.endsWith(".fq.gz") }
        def fastq_gz_files = valid_files.findAll { it.name.endsWith(".fastq.gz") }

        def FILE_FORMAT = ""
        if (fq_gz_files.size() == valid_files.size()) {
            FILE_FORMAT = "fq.gz"
        } else if (fastq_gz_files.size() == valid_files.size()) {
            FILE_FORMAT = "fastq.gz"
        }

        // =================================================================================
        // 3. DETECT SEQUENCING MODE (SE vs PE)
        // =================================================================================

        def r1_files = valid_files.findAll { it.name.contains("_r1") }
        def r2_files = valid_files.findAll { it.name.contains("_r2") }
        def R1_files = valid_files.findAll { it.name.contains("_R1") }
        def R2_files = valid_files.findAll { it.name.contains("_R2") }

        def MODE = ""
        if (r1_files.size() == valid_files.size() || R1_files.size() == valid_files.size()) {
            MODE = "SINGLE_END"
        } else if (r1_files.size() > 0 && r1_files.size() == r2_files.size() && r1_files.size() * 2 == valid_files.size()) {
            MODE = "PAIRED_END"
        } else if (R1_files.size() > 0 && R1_files.size() == R2_files.size() && R1_files.size() * 2 == valid_files.size()) {
            MODE = "PAIRED_END"
        }

        // =================================================================================
        // 4. DETECT READ TAGS (case: _r1/_r2 or _R1/_R2)
        // =================================================================================

        def READ1_TAG = ""
        def READ2_TAG = ""

        if      (MODE == "SINGLE_END" && r1_files.size() == valid_files.size()) { READ1_TAG = "_r1" }
        else if (MODE == "SINGLE_END" && R1_files.size() == valid_files.size()) { READ1_TAG = "_R1" }
        else if (MODE == "PAIRED_END" && r1_files.size() * 2 == valid_files.size()) { READ1_TAG = "_r1" }
        else if (MODE == "PAIRED_END" && R1_files.size() * 2 == valid_files.size()) { READ1_TAG = "_R1" }

        if (MODE == "PAIRED_END") {
            if      (r2_files.size() * 2 == valid_files.size()) { READ2_TAG = "_r2" }
            else if (R2_files.size() * 2 == valid_files.size()) { READ2_TAG = "_R2" }
        }

        // =================================================================================
        // 5. VALIDATION CHECKS (fail fast with clear messages)
        // =================================================================================

        if (valid_files.size() == 0) {
            log.error """
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     ERROR: NO VALID FASTQ FILES FOUND
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             Search Path : ${params.raw_fastq_dir()}
             Files found : ${fastq_files.size()}
             Possible reasons:
               1. Wrong directory path in config
               2. Files not gzipped (.gz extension missing)
               3. Files don't match naming convention
               4. Directory is empty
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            """.stripIndent()
            error "Aborting: Zero files matched the pattern in the directory."
        }

        if (invalid_files.size() > 0) {
            log.error """
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ERROR: ${invalid_files.size()} INVALID FILE(S) DETECTED
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             These files don't match the expected naming pattern:
             Expected: *_R1.fq.gz / *_R2.fq.gz (or _r1/_r2)

               ✗ ${invalid_files.collect{ it.name }.join('\n   ✗ ')}

             Common fixes:
               1. Rename _1.fq.gz → _R1.fq.gz (add 'R')
               2. Add .gz extension if missing
               3. Ensure consistent case (_R1 not _r1 mixed)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            """.stripIndent()
            error "Aborting: Please rename the files listed above."
        }

        if (FILE_FORMAT == "") {
            log.error """
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              FORMAT ERROR: Mixed FASTQ extensions detected
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             .fq.gz files    : ${fq_gz_files.size()}
             .fastq.gz files : ${fastq_gz_files.size()}
             All files must use the SAME extension.
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            """.stripIndent()
            error "Aborting: Standardize all files to either .fq.gz or .fastq.gz"
        }

        if (MODE == "") {
            log.error """
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   ERROR: INCONSISTENT READ PAIR NAMING
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             File counts by tag:
               Lowercase → r1: ${r1_files.size()}, r2: ${r2_files.size()}
               Uppercase → R1: ${R1_files.size()}, R2: ${R2_files.size()}
               Total files   : ${valid_files.size()}
             Solution:
               - Use ONLY _R1/_R2 (or ONLY _r1/_r2)
               - Ensure every R1 has a matching R2 with the same sample name
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            """.stripIndent()
            error "Mixed case or mismatched pairs detected. Please standardize naming."
        }

        // =================================================================================
        // 6. COUNT SAMPLES
        // =================================================================================

        def tumor_files  = valid_files.findAll { it.name.contains("_Tumor") }
        def normal_files = valid_files.findAll { it.name.contains("_Normal") }

        def N_SAMPLES        = (MODE == "PAIRED_END") ? valid_files.size().intdiv(2)  : valid_files.size()
        def N_TUMOR_SAMPLES  = (MODE == "PAIRED_END") ? tumor_files.size().intdiv(2)  : tumor_files.size()
        def N_NORMAL_SAMPLES = (MODE == "PAIRED_END") ? normal_files.size().intdiv(2) : normal_files.size()

        // =================================================================================
        // 7. BUILD GROUPED SAMPLE TUPLES: [sample_id, [R1] or [R1, R2]]
        // =================================================================================

        // .collect{} on a List is a synchronous Groovy loop (not a Nextflow channel operation)
        def all_r1_files = (r1_files + R1_files).sort()
        def grouped_samples = all_r1_files.collect { r1 ->

            def sample_id = ""
            if (params.expt == "scRNASeq") {
                def matcher = (r1.name =~ VALID_PATTERN)
                sample_id = matcher ? matcher[0][1] : r1.simpleName
            } else {
                def idx = r1.name.lastIndexOf(READ1_TAG)
                sample_id = (idx != -1) ? r1.name.take(idx) : r1.simpleName
            }

            if (MODE == "PAIRED_END") {
                // Find R2 mate: reverse string, replace first (rightmost) READ1_TAG, reverse back
                def r1_name = r1.name
                def r2_name = r1_name.reverse().replaceFirst(READ1_TAG.reverse(), READ2_TAG.reverse()).reverse()
                def r2 = r1.parent.resolve(r2_name)

                if (!r2.exists()) {
                    log.error """
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                  ERROR: MISSING R2 PAIR
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     R1 file  : ${r1.name}
                     Expected : ${r2_name}
                     Location : ${r1.parent}
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    """.stripIndent()
                    error "Aborting: Missing paired-end mate for ${r1.name}"
                }
                return [ sample_id, [r1, r2] ]
            } else {
                return [ sample_id, [r1] ]
            }
        }

        // =================================================================================
        // 8. COLLECT UNIQUE SAMPLE NAMES
        // =================================================================================

        // .unique() handles multi-lane samples (same sample_id across lanes)
        def sample_names = grouped_samples.collect { id, fastqs -> id }.unique()

        // =================================================================================
        // 9. LOG SUMMARY
        // =================================================================================

        log.info """
             ============================================================
             FASTQ INPUT VALIDATION SUMMARY
             ============================================================
             FILE FORMAT         : $FILE_FORMAT
              - .fq.gz           : ${fq_gz_files.size()}
              - .fastq.gz        : ${fastq_gz_files.size()}
            ------------------------------------------------------------
             SEQUENCING MODE     : $MODE
             READ TAGS USED      :
              - R1 tag           : $READ1_TAG
              - R2 tag           : ${MODE == "PAIRED_END" ? READ2_TAG : 'N/A'}
             TAG DISTRIBUTION    :
              - Uppercase (R1/R2): ${R1_files.size()} / ${R2_files.size()}
              - Lowercase (r1/r2): ${r1_files.size()} / ${r2_files.size()}
            ------------------------------------------------------------
             SAMPLE SUMMARY      :
              - Total samples    : $N_SAMPLES
              - Tumor samples    : $N_TUMOR_SAMPLES
              - Normal samples   : $N_NORMAL_SAMPLES
              - Other samples    : ${N_SAMPLES - N_TUMOR_SAMPLES - N_NORMAL_SAMPLES}
             TOTAL FASTQ FILES   : ${valid_files.size()}
            ============================================================
        """.stripIndent()

        return [
            mode             : MODE,
            read1_tag        : READ1_TAG,
            read2_tag        : READ2_TAG,
            file_format      : FILE_FORMAT,
            total_samples    : N_SAMPLES,
            n_tumor_samples  : N_TUMOR_SAMPLES,
            n_normal_samples : N_NORMAL_SAMPLES,
            valid_files      : valid_files,
            tumor_files      : tumor_files,
            normal_files     : normal_files,
            sample_names     : sample_names,
            grouped_samples  : grouped_samples
        ]
    }

    emit:
    mode             = results_ch.map { it.mode }
    read1_tag        = results_ch.map { it.read1_tag }
    read2_tag        = results_ch.map { it.read2_tag }
    file_format      = results_ch.map { it.file_format }
    total_samples    = results_ch.map { it.total_samples }
    n_tumor_samples  = results_ch.map { it.n_tumor_samples }
    n_normal_samples = results_ch.map { it.n_normal_samples }
    all_fastq_ch         = results_ch.flatMap { it.valid_files }
    tumor_fastq_ch       = results_ch.flatMap { it.tumor_files }
    normal_fastq_ch      = results_ch.flatMap { it.normal_files }
    grouped_samples_ch   = results_ch.flatMap { it.grouped_samples }
    sample_names_ch      = results_ch.flatMap { it.sample_names }
}
