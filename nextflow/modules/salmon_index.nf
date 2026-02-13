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
    tuple val(species), path(decoy), path(gentrome), val(genome_version)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), path("salmon_index_dir_${genome_version}"),    emit: salmon_index_tuple  // K-mer index directory
    //path("SALMON_INDEX.error.log"),                                  emit: log               // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG = "SALMON_INDEX.error.log"
    def index_dir = "salmon_index_dir_${genome_version}"

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
        --index "${index_dir}" \
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