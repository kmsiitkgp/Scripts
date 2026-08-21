// =============================================================================
// PROCESS: RNA_SALMON_INDEX
// =============================================================================
// Purpose: Builds Salmon k-mer index for fast decoy-aware transcript quantification
//
// What it does:
//   - Creates k-mer hash index from gentrome (transcriptome + genome decoys)
//   - Enables selective alignment mode for more accurate quantification
//   - Built once per genome/annotation version
//
// Typical resources: 16-32GB RAM, ~30-40 minutes on 8 cores
// Index size: ~5-8GB for human transcriptome
// Uses storeDir: built once per genome version, reused across runs
// =============================================================================

process RNA_SALMON_INDEX {

    tag "Building Salmon index ${assembly}.${release}"
    label 'process_high'                          // 16-32GB RAM required for human/mouse

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), path(decoy), path(gentrome), val(assembly), val(release)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), path("salmon_index_dir_${assembly}.${release}"),
        emit: salmon_index_tuple
    //path("SALMON_INDEX.log"),    emit: log    // storeDir: log not captured

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG       = "SALMON_INDEX.log"
    def index_dir = "salmon_index_dir_${assembly}.${release}"

    """
    # Build Salmon index with decoy-aware selective alignment
    # --transcripts : Gentrome FASTA (transcripts + genome decoy sequences)
    # --decoys      : List of decoy sequence names (chromosome names from FASTA headers)
    # --kmerLen 31  : K-mer length; default 31 works well for 75-200bp reads
    # --index       : Output index directory
    # --threads     : Parallel indexing cores
    salmon index \
        --transcripts "${gentrome}" \
        --decoys "${decoy}" \
        --kmerLen 31 \
        --threads "${task.cpus}" \
        --index "${index_dir}" \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: Salmon index build failed" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: Salmon index built (${assembly}.${release})" >> "${LOG}"
    """
}

// =============================================================================
// QUICK REFERENCE
// =============================================================================
//
// When to rebuild:
//   - New genome assembly (GRCh37 → GRCh38)
//   - Updated gene annotation (new GTF version)
//   - Salmon major version upgrade
//
// Output files in salmon_index_dir_<assembly>.<release>/:
//   hash.bin                 : K-mer hash table
//   txpInfo.bin              : Transcript metadata
//   decoys.txt               : Copy of decoy list
//   duplicate_clusters.tsv   : Transcripts with identical sequences
//
// Common issues:
//   - OOM error          → Increase RAM to 32GB+
//   - Low mapping later  → Wrong index version or corrupt build
//   - Decoy format error → Chromosome names in decoy file must match FASTA headers exactly
// =============================================================================
