// =============================================================================
// PROCESS: SCRNA_KB_REF
// =============================================================================
// Purpose: Builds a kallisto|bustools RNA-velocity index (spliced + unspliced)
//          from Ensembl FASTA/GTF, for use by kb-python's "nac" workflow.
//
// What it does:
//   - Runs `kb ref --workflow nac` to build the nascent/mature index used by
//     `kb count --workflow nac` for spliced/unspliced quantification
//   - Emits the index + the t2g files kb count needs downstream
//
// Why separate from SCRNA_DOWNLOAD_10X_REF?
//   - Cell Ranger uses 10X's own pre-built, pre-filtered refdata-gex-* tarball
//     (not a plain FASTA/GTF) — wrong shape for kb-python, which needs to
//     build its own kallisto index from raw Ensembl FASTA/GTF
//   - Input here is RNA_DOWNLOAD_ENSEMBL_REF.out.ref_tuple (same tuple shape
//     already used by the bulk RNAseq workflow) — reused as-is, not refetched
//
// Uses storeDir: index is built once per genome_version and reused across runs,
// same caching pattern as SCRNA_DOWNLOAD_10X_REF.
// =============================================================================

process SCRNA_KB_REF {

    tag "Building kb-python nac index for ${species} (${genome_version})"
    label 'process_high'                          // Index build is CPU/RAM heavy
    label 'omics_py'

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), path(fa), path(gtf), val(genome_version), val(assembly), val(release)
    // species        : "Human", "Mouse", "Xenograft"
    // fa             : Ensembl primary assembly FASTA (RNA_DOWNLOAD_ENSEMBL_REF output)
    // gtf            : Ensembl GTF annotation (RNA_DOWNLOAD_ENSEMBL_REF output)
    // genome_version : Version string (e.g., "GRCm39.115") — storeDir cache key

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(species), path("kb_index_dir_${assembly}.${release}"),
        emit: kb_index_tuple

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG = "KB_REF.${species}.${genome_version}.log"
    def index_dir = "kb_index_dir_${assembly}.${release}"

    """
    mkdir -p "${index_dir}"

    # --workflow nac : builds the nascent (spliced/unspliced) index — the only
    #   workflow this pipeline uses kb-python for, since it exists purely for
    #   velocity input. Never "standard" (gene counts only, no velocity signal).
    # -i / -g        : output index / t2g map
    # -f1/-f2        : cDNA (spliced) / nascent (unspliced) FASTA kb ref extracts
    # -c1/-c2        : cDNA / nascent transcript-to-capture lists (t2c), needed
    #                   by `kb count --workflow nac` to split spliced/unspliced

    kb ref \
        --workflow nac \
        -i index.idx \
        -g t2g.txt \
        -f1 cdna.fa \
        -f2 nascent.fa \
        -c1 cdna_t2c.txt \
        -c2 nascent_t2c.txt \
        -t "${task.cpus}" \
        "${fa}" "${gtf}" \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: kb ref failed for ${species} (${genome_version})" | tee -a "${LOG}" >&2; exit 1; }

    # Move all the generated index files into the index directory
    mv index.idx t2g.txt cdna.fa nascent.fa cdna_t2c.txt nascent_t2c.txt "${index_dir}/"

    echo "✅ SUCCESS: kb ref index built for ${species} (${genome_version})" >> "${LOG}"
    """
}
