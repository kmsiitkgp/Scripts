// =============================================================================
// PROCESS: EXTRACT_GENTROME
// =============================================================================
// Purpose: Creates gentrome (transcriptome + genome) for Salmon indexing
//
// What it does:
//   1. Extracts chromosome names as decoy list
//   2. Extracts transcript sequences from genome using GTF
//   3. Concatenates transcripts + genome (order matters!)
//
// For detailed explanation of gentrome strategy, see: docs/salmon_index.md
// =============================================================================

process EXTRACT_GENTROME {

    tag "Extracting gentrome from ${ref_fasta.name}"
    label 'process_high'                  // Moderate resources: ~8GB RAM, 4+ cores

    input:
    tuple val(species), path(ref_fasta), path(ref_gtf), val(genome_version)

    output:
    tuple val(species),
    path("decoy_${genome_version}.txt"),      // Chromosome names for decoy marking
    path("gentrome_${genome_version}.fa"),    // Combined transcriptome + genome
    val(genome_version),
    emit: decoy_gentrome_tuple
    //path("EXTRACT_GENTROME.error.log"),        emit: error_log  // Process log

    script:

    def LOG = "EXTRACT_GENTROME.error.log"
    def DECOY = "decoy_${genome_version}.txt"
    def GENTROME = "gentrome_${genome_version}.fa"

    """
    # Step 1: Extract chromosome names from FASTA headers
    # Example: ">1 dna:chromosome..." → "1"
    grep "^>" "${ref_fasta}" | cut -d " " -f1 | sed 's/>//' > ${DECOY}

    # Step 2: Extract transcript sequences using GTF coordinates
    # gffread -g: genome FASTA, -w: write transcript sequences
    gffread "${ref_gtf}" \
        -g "${ref_fasta}" \
        -w transcriptome.fa

    # Step 3: Create gentrome (transcripts FIRST, then genome)
    # Order is critical: Salmon uses position to distinguish targets vs decoys
    cat transcriptome.fa "${ref_fasta}" > ${GENTROME} \
        2>> "${LOG}" \
        || { echo "❌ ERROR: Gentrome creation failed" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: Gentrome creation completed" >> "${LOG}"
    """
}