// =============================================================================
// PROCESS: CREATE_GENTROME
// =============================================================================
// Purpose: Builds gentrome (transcriptome + genome) for Salmon decoy-aware indexing
//
// What it does:
//   1. Extracts chromosome names from FASTA headers → decoy list
//   2. Extracts transcript sequences from genome using GTF (gffread)
//   3. Concatenates transcripts + genome into gentrome (order critical)
//
// Why gentrome?
//   - Salmon decoy-aware mode reduces false mappings to genomic sequence
//   - Transcripts must come first; genome sequences serve as decoys only
//   - Without decoys, reads from unannotated regions can map to wrong transcripts
//
// Typical resources: ~8GB RAM, 4 cores
// Uses storeDir: rebuilt only when genome version changes
// =============================================================================

process CREATE_GENTROME {

    tag "Creating gentrome for ${assembly}.${release}"
    label 'process_medium'                        // gffread + cat; ~8GB RAM, 4 cores

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), path(ref_fasta), path(ref_gtf), val(genome_version), val(assembly), val(release)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species),
        path("decoy_${assembly}.${release}.txt"),      // Chromosome names for Salmon decoy marking
        path("gentrome_${assembly}.${release}.fa"),    // Combined transcriptome + genome FASTA
        val(assembly), val(release),
        emit: gentrome_tuple
    //path("CREATE_GENTROME.error.log"),    emit: log    // storeDir: log not captured

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG     = "CREATE_GENTROME.error.log"
    def DECOY   = "decoy_${assembly}.${release}.txt"
    def GENTROME = "gentrome_${assembly}.${release}.fa"

    """
    # Step 1: Extract chromosome names from FASTA headers → decoy list
    # Example: ">1 dna:chromosome chromosome..." → "1"
    grep "^>" "${ref_fasta}" | cut -d " " -f1 | sed 's/>//' > ${DECOY} \
        2>> "${LOG}" \
        || { echo "❌ ERROR: Decoy list extraction failed" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: Decoy list created" >> "${LOG}"

    # Step 2: Extract transcript sequences using GTF coordinates
    # gffread -g: genome FASTA input, -w: write transcript sequences to file
    gffread "${ref_gtf}" \
        -g "${ref_fasta}" \
        -w transcriptome.fa \
        2>> "${LOG}" \
        || { echo "❌ ERROR: gffread transcript extraction failed" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: Transcript sequences extracted" >> "${LOG}"

    # Step 3: Concatenate transcriptome + genome → gentrome
    # Order is critical: Salmon requires transcripts first, genome decoys second
    cat transcriptome.fa "${ref_fasta}" > ${GENTROME} \
        2>> "${LOG}" \
        || { echo "❌ ERROR: Gentrome concatenation failed" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: Gentrome created (${assembly}.${release})" >> "${LOG}"
    """
}
