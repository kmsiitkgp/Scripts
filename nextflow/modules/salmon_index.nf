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
    path(ref_fasta)  // Genome FASTA
    path(ref_gtf)    // Gene annotation GTF

    output:
    path("decoy.txt"),                          emit: decoy      // Chromosome names for decoy marking
    path("gentrome.fa"),                        emit: gentrome   // Combined transcriptome + genome
    path("EXTRACT_GENTROME.error.log"),         emit: error_log  // Process log

    script:

    def LOG = "EXTRACT_GENTROME.error.log"

    """
    # Step 1: Extract chromosome names from FASTA headers
    # Example: ">1 dna:chromosome..." → "1"
    grep "^>" "${ref_fasta}" | cut -d " " -f1 | sed 's/>//' > decoy.txt

    # Step 2: Extract transcript sequences using GTF coordinates
    # gffread -g: genome FASTA, -w: write transcript sequences
    gffread "${ref_gtf}" \
        -g "${ref_fasta}" \
        -w transcriptome.fa

    # Step 3: Create gentrome (transcripts FIRST, then genome)
    # Order is critical: Salmon uses position to distinguish targets vs decoys
    cat transcriptome.fa "${ref_fasta}" > gentrome.fa \
        2>> "${LOG}" \
        || { echo "❌ ERROR: Gentrome creation failed" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: Gentrome creation completed" >> "${LOG}"
    """
}


// =============================================================================
// PROCESS: SALMON_INDEX
// =============================================================================
// Purpose: Builds Salmon k-mer index for fast transcript quantification
//
// What it does:
//   - Creates k-mer hash index from gentrome
//   - Enables decoy-aware selective alignment
//   - Built once per genome/annotation version
//
// Typical resources: 16-32GB RAM, ~30-40 minutes on 8 cores
// Index size: ~5-8GB for human transcriptome
//
// For detailed explanation, see: docs/salmon_index.md
// =============================================================================

process SALMON_INDEX {

    tag "Indexing ${gentrome.name}"
    label 'process_high'  // Substantial RAM required: 16-32GB for human

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    path(gentrome)  // Gentrome FASTA (transcriptome + genome as decoys)
    path(decoy)     // Text file with chromosome names (one per line)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("salmon_index_dir"),               emit: salmon_index_dir  // K-mer index directory
    path("SALMON_INDEX.error.log"),            emit: log               // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG = "SALMON_INDEX.error.log"

    """
    # Build Salmon index with decoy-aware mode
    # --transcripts: Gentrome FASTA (targets + decoys)
    # --decoys: List of decoy sequence names (chromosome names)
    # --kmerLen 31: K-mer size (default for 75-200bp reads)
    # --threads: Parallel indexing

    salmon index \
        --transcripts ${gentrome} \
        --decoys ${decoy} \
        --threads "${task.cpus}" \
        --kmerLen 31 \
        --index salmon_index_dir \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: SALMON index generation failed" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: SALMON index generation completed" >> "${LOG}"
    """
}

// =========================================================================================
// QUICK REFERENCE
// =========================================================================================
//
// When to rebuild index:
//   - New genome assembly version (e.g., GRCh37 → GRCh38)
//   - Updated gene annotation (new GTF version)
//   - Salmon major version upgrade
//   - Different k-mer length needed
//
// Index contents (salmon_index_dir/):
//   - hash.bin: K-mer hash table
//   - txpInfo.bin: Transcript metadata
//   - decoys.txt: Copy of decoy list
//   - duplicate_clusters.tsv: Transcript groups
//
// Common issues:
//   - OOM error → Increase to 32GB+ RAM
//   - Low mapping rates later → Wrong index or corrupted build
//   - Decoy format error → Check chromosome names match FASTA exactly
//
// For comprehensive guide, see: docs/salmon_index.md
// ========================================================================================='