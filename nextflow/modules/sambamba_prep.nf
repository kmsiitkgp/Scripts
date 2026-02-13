// =========================================================================================
// PROCESS: SAMBAMBA_PREP
// =========================================================================================
// Purpose: Prepares BAM files for downstream QC analysis
//
// What it does:
//   1. Renames BAM to standardized format
//   2. Indexes full BAM file (creates .bai)
//   3. Calculates read length from BAM
//   4. Creates 1M read subsample for fast QC
//   5. Indexes subsampled BAM
//
// Why subsample?
//   - RSeQC gene body coverage is slow on full BAM (~30 min)
//   - 1M reads sufficient for detecting 3' bias (~2-3 min)
//   - 10-15x speedup with minimal accuracy loss
//
// Typical resources: 4-8GB RAM, ~5-10 minutes per sample
//
// For detailed explanation, see: docs/bam_preparation.md
// =========================================================================================

process SAMBAMBA_PREP {

    tag "Indexing and subsampling ${raw_bam.name}"
    label 'process_medium'                  // Moderate resources: 4-8GB RAM, 4 cores

    //publishDir { "${params.proj_dir()}/${species}_${type}/04.STAR/bam" },    mode: 'copy',    pattern: "*.{bam,bam.bai}",    saveAs: { filename -> filename.contains('.1M') ? null : filename }
    publishDir { "${params.proj_dir()}/${species}_${type}/07.Logs" },        mode: 'copy',    pattern: "*.SAMBAMBA_PREP.error.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), val(type), val(sample_id), path(raw_bam)  // [species, sample_id, sorted_bam_from_STAR]

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), val(type),
        val(sample_id),
        path("${sample_id}.bam"),                                        // Renamed full BAM
        path("${sample_id}.bam.bai"),                                    // Full BAM index
        path("${sample_id}.1M.bam"),                                     // Subsampled BAM (1M reads)
        path("${sample_id}.1M.bam.bai"),                                 // Subsample index
        path("${sample_id}.read_length.txt"),                            // Detected read length
        emit: bam_indexed
    path("${sample_id}.SAMBAMBA_PREP.error.log"),     emit: error_log    // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG = "${sample_id}.SAMBAMBA_PREP.error.log"

    """
    # Step 1: Standardize BAM filename
    # STAR outputs long names like "Sample.Aligned.sortedByCoord.out.bam"
    # Rename to simple "Sample.bam" for consistency
    mv "${raw_bam}" "${sample_id}.bam"

    # Step 2: Index full BAM
    # Creates ${sample_id}.bam.bai for random access
    # sambamba is faster than samtools (multi-threaded)
    sambamba index -t "${task.cpus}" "${sample_id}.bam" \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: BAM indexing failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: BAM indexed for ${sample_id}" >> "${LOG}"

    # Step 3: Calculate read length
    # Extract from first 10K reads, find most common length
    # Used by RSeQC scripts that need read length parameter
    sambamba view "${sample_id}.bam" | \
        head -n 10000 | \
        awk '{print length(\$10)}' | \
        sort | uniq -c | sort -rn | head -n 1 | \
        awk '{print \$2}' > "${sample_id}.read_length.txt"

    # Step 4: Calculate subsampling fraction
    # Target: 1M reads (sufficient for gene body coverage QC)
    # If BAM has <1M reads, use all (fraction = 1.0)
    TOTAL=\$(sambamba view -c -t ${task.cpus} "${sample_id}.bam")
    FRACTION=\$(awk -v total=\$TOTAL 'BEGIN {
        if (total <= 1000000) print 1.0;
        else print 1000000/total
    }')

    # Step 5: Subsample BAM
    # --subsampling-seed=42: Reproducible random sampling
    sambamba view \
        -t "${task.cpus}" \
        -f bam \
        -s \$FRACTION \
        --subsampling-seed=42 \
        "${sample_id}.bam" \
        -o "${sample_id}.1M.bam" \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: Subsampling failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: Subsampled to ~1M reads for ${sample_id}" >> "${LOG}"

    # Step 6: Index subsampled BAM
    sambamba index -t "${task.cpus}" "${sample_id}.1M.bam" \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: Subsample index failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: Subsample indexed for ${sample_id}" >> "${LOG}"
    """
}

// =========================================================================================
// QUICK REFERENCE
// =========================================================================================
//
// Input: Raw BAM from aligner (STAR, HISAT2, etc.)
// Outputs:
//   - Full BAM + index (for visualization, variant calling)
//   - 1M subsample + index (for fast QC)
//   - Read length file (for RSeQC parameters)
//
// Why sambamba?
//   - Multi-threaded (faster than samtools)
//   - Compatible with samtools format
//   - Efficient subsampling
//
// Subsampling strategy:
//   - Target: 1M reads
//   - Seed: 42 (reproducible)
//   - Random sampling (not first N reads)
//
// Typical speedup:
//   - Gene body coverage: 10-15x faster
//   - 30 minutes → 2-3 minutes
//
// For comprehensive guide, see: docs/bam_preparation.md
// ========================================================================================='