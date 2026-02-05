// =========================================================================================
// WORKFLOW: VALIDATE_INPUT
// =========================================================================================
// Purpose: Validates FASTQ naming, detects SE/PE mode, creates sample channels
//
// What it does:
//   1. Validates FASTQ file naming conventions
//   2. Detects single-end vs paired-end sequencing
//   3. Ensures consistent file format (.fq.gz vs .fastq.gz)
//   4. Pairs R1/R2 files correctly
//   5. Creates organized sample channels for downstream processes
//
// Why this matters:
//   - Catches naming errors BEFORE wasting compute time
//   - Prevents cryptic pairing errors later in pipeline
//   - Provides clear, actionable error messages
//
// For detailed explanation, see: docs/validate_input.md
// =========================================================================================

workflow VALIDATE_INPUT {

    take:
    fastq_files_ch                         // This is now a list of file objects

    main:
    log.info "==> Validating file naming conventions"

    // =================================================================================
    // 1. VALIDATE FASTQ FILES
    // =================================================================================

    // Using .collect() on a channel (e.g., fastq_files_ch) acts as a synchronization
    // barrier. It waits for the complete signal from the upstream channel, then gathers
    // all emitted items (e.g., Path objects) into a single new event: a List object.
    // The resulting channel emits this List exactly once.
    // Using .map{} on a channel (e.g., fastq_files_ch) is a 1-to-1 transformation of the
    // live data stream. As each item is emitted by the upstream channel, the code within
    // {} is applied to it immediately. The result is put onto a new channel, maintaining
    // the same number of items (events) as the original.

    results_ch = fastq_files_ch.collect().map { fastq_list ->

        // Sort fastq files by name
        def fastq_files = fastq_list.sort { it.name }

        // Define validation regex pattern
        // Expected: *_R1.fq.gz, *_R2.fq.gz (or _r1/_r2)
        // Optional: _Tumor or _Normal designation
        def VALID_PATTERN
        if (params.expt == "RNASeq") {
            //VALID_PATTERN = ~/.*((_Tumor|_Normal))?.*((_R|_r)[12]).*\.f(q|astq)\.gz/
            VALID_PATTERN = ~/.*(_Tumor|_Normal)?.*(_[Rr][12]).*\.f(q|astq)\.gz/
        }
        else if (params.expt == "scRNASeq") {
            VALID_PATTERN = ~/^([A-Za-z0-9-]+)_S\d+_L\d{3}_(R[12]|I[12])_001\.fastq\.gz$/
        } else {
            error "❌ Unknown experiment type: ${params.expt}. Options: 'RNASeq', 'scRNASeq'"
        }

        // Separate valid from invalid files
        def valid_files   = fastq_files.findAll { it.name ==~ VALID_PATTERN }
        def invalid_files = fastq_files.findAll { !(it.name ==~ VALID_PATTERN) }

        // =================================================================================
        // 2. DETECT FILE FORMAT
        // =================================================================================

        def fq_gz_files     = valid_files.findAll { it.name.endsWith(".fq.gz") }
        def fastq_gz_files  = valid_files.findAll { it.name.endsWith(".fastq.gz") }

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

        // SE: All files are R1, no R2
        def MODE = ""
        if (r1_files.size() == valid_files.size() ||  R1_files.size() == valid_files.size() ) {
            MODE = "SINGLE_END"
        }
        // PE (lowercase): Equal R1/R2 files
        else if (r1_files.size() > 0 && r1_files.size() == r2_files.size() && r1_files.size() * 2 == valid_files.size()) {
            MODE = "PAIRED_END"
        }
        // PE (uppercase): Equal R1/R2 files
        else if (R1_files.size() > 0 && R1_files.size() == R2_files.size() && R1_files.size() * 2 == valid_files.size()) {
            MODE = "PAIRED_END"
        }

        // =================================================================================
        // 4. DETECT READ TAGS
        // =================================================================================

        def Read1_TAG = ""
        def Read2_TAG = ""
        if (MODE == "SINGLE_END" && r1_files.size() == valid_files.size()) {
            Read1_TAG = "_r1"
        } else if (MODE == "SINGLE_END" && R1_files.size() == valid_files.size()) {
            Read1_TAG = "_R1"
        } else if (MODE == "PAIRED_END" && r1_files.size() * 2 == valid_files.size()){
            Read1_TAG = "_r1"
        } else if (MODE == "PAIRED_END" && R1_files.size() * 2 == valid_files.size()){
            Read1_TAG = "_R1"
        }

        if (MODE == "PAIRED_END") {
            if (r2_files.size() * 2 == valid_files.size()) {
                Read2_TAG = "_r2"
            } else if (R2_files.size() * 2 == valid_files.size()) {
                Read2_TAG = "_R2"
            }
        }

        // =================================================================================
        // 5. VALIDATION CHECKS (Fail Fast)
        // =================================================================================

        // CHECK A: At least one valid file found
        if (valid_files.size() == 0) {
            log.error """
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     ERROR: NO VALID FASTQ FILES FOUND
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             Search Path : ${params.raw_fastq_dir()}
             Pattern     : *.f*q.gz
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

        // CHECK B: Report invalid file names
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

        // CHECK C: Ensure consistent file format
        if (FILE_FORMAT == "") {
            log.error """
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              FORMAT ERROR: Mixed FASTQ extensions detected
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             .fq.gz files    : ${fq_gz_files.size()}
             .fastq.gz files : ${fastq_gz_files.size()}

             All files must use the SAME extension.
             Choose one and rename all files:
               Option 1: Rename all to .fq.gz
               Option 2: Rename all to .fastq.gz
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            """.stripIndent()
            error "Aborting: Please standardize all files to either .fq.gz or .fastq.gz"
        }

        // CHECK D: Verify consistent read pair tags
        if (MODE == ""){
            log.error """
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   ERROR: INCONSISTENT READ PAIR NAMING
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             File counts by tag:
               Lowercase → r1: ${r1_files.size()}, r2: ${r2_files.size()}
               Uppercase → R1: ${R1_files.size()}, R2: ${R2_files.size()}
               Total files   : ${valid_files.size()}

             Possible issues:
               1. Mixed case (_r1 and _R1 in same dataset)
               2. Mismatched pairs (different number of R1 vs R2)
               3. Missing R2 files for some samples

             Solution:
               - Use ONLY _R1/_R2 (or ONLY _r1/_r2)
               - Ensure every R1 has matching R2 with same sample name
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            """.stripIndent()
            error "Mixed case or mismatched pairs detected. Please standardize naming."
        }

        // =================================================================================
        // 6. COUNT SAMPLES (Tumor/Normal classification)
        // =================================================================================

        def tumor_files  = valid_files.findAll { it.name.contains("_Tumor") }
        def normal_files = valid_files.findAll { it.name.contains("_Normal") }

        def N_SAMPLES
        def N_TUMOR_SAMPLES
        def N_NORMAL_SAMPLES

        N_SAMPLES        = (MODE == "PAIRED_END") ? valid_files.size().intdiv(2)  : valid_files.size()
        N_TUMOR_SAMPLES  = (MODE == "PAIRED_END") ? tumor_files.size().intdiv(2)  : tumor_files.size()
        N_NORMAL_SAMPLES = (MODE == "PAIRED_END") ? normal_files.size().intdiv(2) : normal_files.size()

        // =================================================================================
        // 7. CREATE SAMPLE-FASTQ CHANNELS
        // =================================================================================

        // Using .collect{} on a List object (e.g., all_r1_files) is a synchronous loop within
        // the script's memory. It iterates through every element in the existing List, applies
        // the code within {} to each, and returns a new List object containing the transformed
        // items. This happens entirely within the current process/block and has no "channel"
        // properties.

        def all_r1_files = (r1_files + R1_files).sort()
        def grouped_samples = all_r1_files.collect { r1 ->

            def sample_id = ""
            if (params.expt == "scRNASeq") {
                // Use the regex capture group to get exactly what's in the first parentheses
                def matcher = (r1.name =~ VALID_PATTERN)
                sample_id = matcher ? matcher[0][1] : r1.simpleName
            } else {
                // Extract sample ID by removing Read1_TAG
                def idx = r1.name.lastIndexOf(Read1_TAG)
                sample_id = (idx != -1) ? r1.name.take(idx) : r1.simpleName
            }

            if (MODE == "PAIRED_END") {
                // Find R2 mate using reverse-replace trick
                // Reverses string, replaces first (=last) occurrence, reverses back
                def r1_name = r1.name
                def r2_name = r1_name.reverse().replaceFirst(Read1_TAG.reverse(), Read2_TAG.reverse()).reverse()
                def r2 = r1.parent.resolve(r2_name)

                // Verify R2 exists
                if (!r2.exists()) {
                    log.error """
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                  ERROR: MISSING R2 PAIR
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     R1 file   : ${r1.name}
                     Expected  : ${r2_name}
                     Location  : ${r1.parent}

                     This R1 file has no matching R2 file.
                     Check for:
                       1. Typo in R2 filename
                       2. Missing file (incomplete download/transfer)
                       3. Case mismatch (R1 vs r1)
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    """.stripIndent()
                    error "Aborting: Missing paired-end mate for ${r1.name}"
                }
                return [ sample_id, [r1, r2] ]
            } else {
                // Single-end: return R1 only (as list for consistency)
                return [ sample_id, [r1] ]
            }
        }

        // =================================================================================
        // 8. CREATE SAMPLE CHANNELS
        // =================================================================================

        // .unique() ensures if a sample has multiple lanes, you only get the name once
        def sample_names = grouped_samples.collect { id, fastqs -> id }.unique()

        // =================================================================================
        // 9. PRINT SUMMARY
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
              - R1 tag           : $Read1_TAG
              - R2 tag           : ${MODE == "PAIRED_END" ? Read2_TAG : 'N/A'}
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

        // Return everything to the results channel
        return [
            // Metadata outputs
            mode: MODE,
            read1_tag        : Read1_TAG,
            read2_tag        : Read2_TAG,
            file_format      : FILE_FORMAT,
            total_samples    : N_SAMPLES,
            n_tumor_samples  : N_TUMOR_SAMPLES,
            n_normal_samples : N_NORMAL_SAMPLES,
            // Raw file lists
            valid_files      : valid_files,
            tumor_files      : tumor_files,
            normal_files     : normal_files,
            // Sample name list [sample_id]
            sample_names     : sample_names,
            // Grouped sample tuples [sample_id, [R1] or [R1, R2]]
            grouped_samples  : grouped_samples
        ]
    }

    emit:
    // Metadata outputs (As Single Values)
    mode             = results_ch.map { it.mode }
    read1_tag        = results_ch.map { it.read1_tag }
    read2_tag        = results_ch.map { it.read2_tag }
    file_format      = results_ch.map { it.file_format }
    total_samples    = results_ch.map { it.total_samples }
    n_tumor_samples  = results_ch.map { it.n_tumor_samples }
    n_normal_samples = results_ch.map { it.n_normal_samples }

    // File channels (Using flatMap to turn lists back into streams)
    all_fastq_ch     = results_ch.flatMap { it.valid_files }
    tumor_fastq_ch   = results_ch.flatMap { it.tumor_files }
    normal_fastq_ch  = results_ch.flatMap { it.normal_files }

    // Grouped sample tuples: [sample_id, [R1, R2]]
    grouped_samples_ch = results_ch.flatMap { it.grouped_samples }

    // Unique sample names
    sample_names_ch    = results_ch.flatMap { it.sample_names }
}

// =========================================================================================
// QUICK REFERENCE
// =========================================================================================
//
// Valid naming patterns:
//   Sample1_R1.fastq.gz, Sample1_R2.fastq.gz
//   Sample1_Tumor_r1.fq.gz, Sample1_Tumor_r2.fq.gz
//   Patient001_Normal_R1.fq.gz, Patient001_Normal_R2.fq.gz
//
// Invalid patterns:
//   Sample1_1.fq.gz (missing "R")
//   Sample1.fastq (not gzipped)
//   Sample1_R3.fq.gz (should be R1 or R2)
//
// Output channels:
//   samples: [sample_id, [fastq_files]] - main output for downstream
//   mode: "SINGLE_END" or "PAIRED_END"
//   fastq_ch: All valid FASTQ files
//
// For comprehensive guide, see: docs/validate_input.md
// ========================================================================================='