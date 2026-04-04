// =============================================================================
// PROCESS: XENGSORT_INDEX
// =============================================================================
// Purpose: Builds XenGSort k-mer index for separating human (graft) and mouse (host) reads
//
// What it does:
//   - Takes Human and Mouse reference FASTAs
//   - Builds a k-mer hash that can distinguish graft vs host reads at classification time
//   - Moves index files into a named directory for clean Nextflow staging
//
// Why XenGSort?
//   - PDX samples contain a mix of human tumor and mouse stroma reads
//   - Reads must be separated before species-specific alignment and quantification
//   - XenGSort is fast and memory-efficient compared to alignment-based separation
//
// Typical resources: 30-50GB RAM, 1-2 hours
// Uses storeDir: built once per genome version combination, reused across runs
// =============================================================================

process XENGSORT_INDEX {

    tag "Building XenGSort index (${genome_version})"
    label 'process_high'                          // Large k-mer hash; 30-50GB RAM required

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    path(href_fasta)        // Human reference FASTA (graft genome)
    path(mref_fasta)        // Mouse reference FASTA (host genome)
    val(genome_version)     // Composite version string (e.g., "GRCh38.115_GRCm39.115")

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("xengsort_index_dir_${genome_version}"),
        emit: xengsort_index_dir
    //path("XENGSORT_INDEX.error.log"),    emit: log    // storeDir: log not captured

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG       = "XENGSORT_INDEX.error.log"
    def index_dir = "xengsort_index_dir_${genome_version}"

    """
    # Build XenGSort k-mer hash index
    # --host   : Mouse genome (reads assigned here are host)
    # --graft  : Human genome (reads assigned here are graft)
    # --kmer 25: K-mer length for read classification
    # --nobjects 4500000000: Hash table size; sufficient for human + mouse combined
    xengsort index \
        --index "${genome_version}" \
        --host "${mref_fasta}" \
        --graft "${href_fasta}" \
        --kmer 25 \
        --weakthreads "${task.cpus}" \
        --nobjects 4500000000 \
        2>> "${LOG}" \
        || { echo "❌ ERROR: XenGSort index build failed" | tee -a "${LOG}" >&2; exit 1; }

    # Collect index files into a named directory for clean Nextflow output staging
    mkdir -p "${index_dir}"
    mv "${genome_version}".info "${genome_version}".hash "${index_dir}/" \
        2>> "${LOG}" \
        || { echo "❌ ERROR: Failed to move index files to ${index_dir}" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: XenGSort index built (${genome_version})" >> "${LOG}"
    """
}

// =============================================================================
// QUICK REFERENCE
// =============================================================================
//
// Called once per pipeline run; index reused for all samples via .collect()
//
// Output files in xengsort_index_dir_<version>/:
//   <version>.hash : K-mer hash table
//   <version>.info : Index metadata (k-mer size, genome sizes, etc.)
//
// When to rebuild:
//   - New Human or Mouse genome assembly version
//   - XenGSort major version change
//
// Common issues:
//   - OOM error          → Increase RAM; try reducing --nobjects
//   - Low classification → Check genome versions match the sequenced species
// =============================================================================
