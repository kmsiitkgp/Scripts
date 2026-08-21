// =============================================================================
// PROCESS: RNA_RSEQC_PREP_BAM
// =============================================================================
// Purpose: Prepares BAM files for downstream QC and analysis
//
// What it does:
//   1. Renames BAM to standardized sample_id.bam format
//   2. Indexes full BAM (creates .bam.bai for random access)
//   3. Detects read length from first 10K reads
//   4. Creates a 1M-read subsample for fast RSeQC gene body coverage
//   5. Indexes the subsampled BAM
//
// Why subsample?
//   - RSeQC gene body coverage on full BAM can take 30+ minutes
//   - 1M reads is sufficient to detect 3' bias (~2-3 min runtime)
//   - ~10-15x speedup with negligible accuracy loss
//
// Typical resources: 4-8GB RAM, 5-10 minutes per sample
// =============================================================================

process RNA_RSEQC_PREP_BAM {

    tag "BAM prep: ${species} / ${sample_id}"
    label 'process_medium'                        // 4-8GB RAM, 4 cores

    publishDir = [
        //[path: { "${params.star_dir()}/bam" },           mode: 'copy', pattern: "*.{bam,bam.bai}",    saveAs: { filename -> filename.contains('.1M') ? null : filename }],
        [path: { "${params.proj_dir()}/${species}/Logs" }, mode: 'copy', pattern: "*.RSEQC_PREP.log"  ]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), val(sample_id), path(raw_bam)
    // raw_bam: Coordinate-sorted BAM from STAR_ALIGN (long filename e.g. *.Aligned.sortedByCoord.out.bam)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species),
        val(sample_id),
        path("${sample_id}.bam"),                  // Renamed full BAM
        path("${sample_id}.bam.bai"),              // Full BAM index
        path("${sample_id}.1M.bam"),               // 1M-read subsample
        path("${sample_id}.1M.bam.bai"),           // Subsample index
        path("${sample_id}.read_length.txt"),      // Detected read length (used by RSeQC)
        emit: bam_indexed
    path("${sample_id}.RSEQC_PREP.log"),    emit: error_log    // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG = "${sample_id}.RSEQC_PREP.log"

    """
    # Step 1: Standardize BAM filename
    # STAR outputs long names (e.g. Sample.Aligned.sortedByCoord.out.bam)
    mv "${raw_bam}" "${sample_id}.bam"

    # Step 2: Index full BAM
    # sambamba is faster than samtools for indexing (multi-threaded)
    sambamba index -t "${task.cpus}" "${sample_id}.bam" \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: BAM indexing failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: BAM indexed for ${sample_id}" >> "${LOG}"

    # Step 3: Detect read length from first 10K reads
    # Extract CIGAR field (col 10) length; take most frequent value
    # Used by RSeQC deletion_profile and mismatch_profile
    sambamba view "${sample_id}.bam" | \
        head -n 10000 | \
        awk '{print length(\$10)}' | \
        sort | uniq -c | sort -rn | head -n 1 | \
        awk '{print \$2}' > "${sample_id}.read_length.txt"

    # Step 4: Calculate subsampling fraction to target ~1M reads
    # If total reads <= 1M, use fraction 1.0 (keep all reads)
    TOTAL=\$(sambamba view -c -t "${task.cpus}" "${sample_id}.bam")
    FRACTION=\$(awk -v total=\$TOTAL 'BEGIN {
        if (total <= 1000000) print 1.0;
        else print 1000000/total
    }')

    # Step 5: Subsample BAM to ~1M reads
    # --subsampling-seed 42: Reproducible random sampling across runs
    sambamba view \
        -t "${task.cpus}" \
        -f bam \
        -s \$FRACTION \
        --subsampling-seed 42 \
        "${sample_id}.bam" \
        -o "${sample_id}.1M.bam" \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: BAM subsampling failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: BAM subsampled to ~1M reads for ${sample_id}" >> "${LOG}"

    # Step 6: Index subsampled BAM
    sambamba index -t "${task.cpus}" "${sample_id}.1M.bam" \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: Subsample BAM indexing failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: Subsample BAM indexed for ${sample_id}" >> "${LOG}"
    """
}

// =============================================================================
// QUICK REFERENCE
// =============================================================================
//
// Outputs passed to RSEQC:
//   full BAM + .bai         : For read_distribution, junction_annotation, artifact profiles
//   1M subsample + .bai     : For gene body coverage and junction saturation (speed)
//   read_length.txt         : For deletion_profile and mismatch_profile --read-align-length
//
// Why sambamba over samtools?
//   - Multi-threaded indexing and view (faster at scale)
//   - Compatible with samtools BAM format
//
// Common issues:
//   - "truncated file" error → Upstream STAR BAM may be incomplete; check STAR log
//   - read_length.txt empty  → BAM may have too few reads; check alignment rate
// =============================================================================
