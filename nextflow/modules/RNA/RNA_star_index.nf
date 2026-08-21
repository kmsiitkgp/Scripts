// =============================================================================
// PROCESS: RNA_STAR_INDEX
// =============================================================================
// Purpose: Builds STAR genome index for splice-aware RNA-seq alignment
//
// What it does:
//   - Builds suffix array from genome FASTA sequence
//   - Extracts and indexes splice junctions from GTF annotation
//   - Creates all searchable index structures in index_dir/
//
// Typical resources: 30-50GB RAM, 1-2 hours on 8 cores (human genome)
// Index size: ~25-30GB for human
// Uses storeDir: built once per genome version, reused across runs
// =============================================================================

process RNA_STAR_INDEX {

    tag "Building STAR index ${assembly}.${release}"
    label 'process_high'                          // 30-50GB RAM required for human/mouse

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), path(ref_fasta), path(ref_gtf), val(genome_version), val(assembly), val(release)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), path("star_index_dir_${assembly}.${release}"),
        emit: star_index_tuple
    //path("STAR_INDEX.log"),    emit: log    // storeDir: log not captured

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG       = "STAR_INDEX.log"
    def index_dir = "star_index_dir_${assembly}.${release}"

    """
    mkdir -p "${index_dir}"

    # Build STAR genome index
    # --runMode genomeGenerate : Index creation mode (not alignment)
    # --genomeDir              : Output directory for index files
    # --genomeFastaFiles       : Input genome FASTA
    # --sjdbGTFfile            : Gene annotation for splice junction database
    # --sjdbOverhang 100       : ReadLength - 1; 100 works well for 75-150bp reads
    # --genomeSAindexNbases 14 : Suffix array sparsity; optimal for human/mouse genome size
    # --runThreadN             : Parallel indexing cores
    STAR \
        --runMode genomeGenerate \
        --runThreadN "${task.cpus}" \
        --genomeDir "${index_dir}" \
        --genomeFastaFiles "${ref_fasta}" \
        --sjdbGTFfile "${ref_gtf}" \
        --sjdbOverhang 100 \
        --genomeSAindexNbases 14 \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: STAR index build failed" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: STAR index built (${assembly}.${release})" >> "${LOG}"
    """
}

// =============================================================================
// QUICK REFERENCE
// =============================================================================
//
// When to rebuild:
//   - New genome assembly (GRCh37 → GRCh38)
//   - Major GTF update (new gene models)
//   - STAR major version change
//
// Key parameters:
//   sjdbOverhang  : ReadLength - 1 (100 works for 75-150bp reads)
//   genomeSAindexNbases : 14 for human/mouse; reduce for small genomes
//
// Output files in star_index_dir_<assembly>.<release>/:
//   SA             : Suffix array (~20-25GB for human)
//   Genome         : Packed genome sequence
//   sjdbList.out.tab : Splice junctions from GTF
//   chrName.txt    : Chromosome names
//
// Common issues:
//   - OOM error          → Increase RAM allocation in nextflow.config
//   - "SA size error"    → Reduce genomeSAindexNbases
//   - GTF parsing error  → Verify chromosome names match between FASTA and GTF
// =============================================================================
