// =========================================================================================
// PROCESS: STAR_INDEX
// =========================================================================================
// Purpose: Builds STAR genome index for splice-aware RNA-seq alignment
//
// What it does:
//   - Builds suffix array from genome sequence
//   - Extracts splice junctions from GTF annotation
//   - Creates searchable index structures
//   - Stores everything in star_index_dir/
//
// Typical resources: 30-50GB RAM, 1-2 hours on 8 cores (human genome)
// Index size: ~25-30GB for human
//
// For detailed explanation, see: docs/star_index.md
// =========================================================================================

process STAR_INDEX {

    tag "Indexing ${ref_fasta.name}"
    label 'process_high'                      // STAR indexing requires 30-50GB RAM for human

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), path(ref_fasta), path(ref_gtf), val(genome_version)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), path("star_index_dir_${genome_version}"),    emit: star_index_tuple    // STAR genome index
    //path("STAR_INDEX.error.log"),              emit: error_log         // Process log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG = "STAR_INDEX.error.log"
    def index_dir = "star_index_dir_${genome_version}"

    """
    # Create output directory
    mkdir -p "${index_dir}"

    # Build STAR index
    # --runMode genomeGenerate: Index creation mode (not alignment)
    # --runThreadN: Use multiple cores for parallel indexing
    # --genomeDir: Output directory for index files
    # --genomeFastaFiles: Input genome sequence
    # --sjdbGTFfile: Gene annotations for splice junction database
    # --sjdbOverhang 100: Optimal for 75-150bp reads (ReadLength - 1)
    # --genomeSAindexNbases 14: Suffix array sparsity (optimal for human/mouse)

    STAR --runMode genomeGenerate \
        --runThreadN "${task.cpus}" \
        --genomeDir "${index_dir}" \
        --genomeFastaFiles "${ref_fasta}" \
        --sjdbGTFfile "${ref_gtf}" \
        --sjdbOverhang 100 \
        --genomeSAindexNbases 14 \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: STAR index generation failed" | tee -a "${LOG}" >&2; exit 1; }

    echo "✅ SUCCESS: STAR index generation completed" >> "${LOG}"
    """
}

// =========================================================================================
// QUICK REFERENCE
// =========================================================================================
//
// When to rebuild index:
//   - New genome assembly (GRCh37 → GRCh38)
//   - Major GTF update (new gene models)
//   - STAR major version change
//   - Significantly different read lengths
//
// Key parameters:
//   sjdbOverhang: ReadLength - 1 (100 works for 75-150bp reads)
//   genomeSAindexNbases: 14 for human/mouse, smaller for tiny genomes
//
// Output files in star_index_dir/:
//   SA: Suffix array (~20-25GB for human)
//   Genome: Packed genome sequence
//   sjdbList.out.tab: Splice junctions from GTF
//   chrName.txt, chrLength.txt: Chromosome metadata
//
// Common issues:
//   - OOM error → Increase RAM allocation
//   - "SA size error" → Reduce genomeSAindexNbases
//   - GTF parsing error → Check chromosome name match with FASTA
//
// For detailed guide, see: docs/star_index.md
// ========================================================================================='