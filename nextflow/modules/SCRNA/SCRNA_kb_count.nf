// =============================================================================
// PROCESS: SCRNA_KB_COUNT
// =============================================================================
// Purpose: Spliced/unspliced UMI quantification per sample using kb-python,
//          for RNA velocity input (scVelo / veloVI). Runs directly on FASTQs —
//          no BAM required, independent of SCRNA_CELLRANGER_COUNT.
//
// What it does:
//   - Runs `kb count --workflow nac` against the kb_ref index for this species
//   - Emits an .h5ad file per sample with spliced/unspliced/ambiguous layers
//
// Note: This process does NOT depend on SCRNA_CELLRANGER_COUNT — it consumes
// the same raw FASTQ channel (sample_fastq_metadata_ch) directly, so it can
// run in parallel with the Cell Ranger branch.
//
// Typical resources: scales with kb ref's index build; count itself is lighter
// than Cell Ranger (no STAR alignment) — process_high still safe default.
// =============================================================================

process SCRNA_KB_COUNT {

    tag "kb count (nac): ${sample_id}"
    label 'process_high'
    label 'omics_py'

    publishDir = [
        [path: { "${params.proj_dir()}/${params.species}/04.KB/${sample_id}" }, mode: 'copy', pattern: "*.h5ad"],
        [path: { "${params.proj_dir()}/${params.species}/Logs" },               mode: 'copy', pattern: "*.KB_COUNT.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(species), val(sample_id), path(fastq_files), path(kb_index_dir)
    // species       : "Human", "Mouse", "Xenograft"
    // sample_id     : Sample identifier — MUST match the sample_id string used
    //                 upstream (CreateSeuratObject project=) so barcode prefixing
    //                 lines up during the later merge-with-Seurat step
    // fastq_files   : [R1, R2 across all lanes] — same FASTQs SCRNA_CELLRANGER_COUNT uses
    // kb_index_dir  : SCRNA_KB_REF's bundled index dir for this species
    //                 (kb_index_dir_${assembly}.${release}/ containing index.idx,
    //                 t2g.txt, cdna_t2c.txt, nascent_t2c.txt) — combined by species
    //                 in scrnaseq.nf before this process is called
    val(kb_technology)
    // e.g. "10XV4" — pre-resolved in scrnaseq.nf from yaml params.kb.technology
    // (kb-python has no "auto" chemistry detection, unlike Cell Ranger)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(sample_id), path("${sample_id}.h5ad"),  emit: h5ad_files
    path("${sample_id}.KB_COUNT.log"),                emit: error_log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    def LOG = "${sample_id}.KB_COUNT.log"

    """
    set -euo pipefail

    # kb count has NO Cell Ranger-style filename auto-matching — it takes the
    # FASTQ list literally as ordered R1 R2 R1 R2... pairs (R1 = barcode/UMI,
    # R2 = cDNA read). Nextflow's channel/glob ordering is not guaranteed to
    # already be lane-sorted R1-before-R2, so sort explicitly here rather than
    # trust upstream ordering.
    R1_FILES=(\$(ls *_R1_*.fastq.gz | sort))
    R2_FILES=(\$(ls *_R2_*.fastq.gz | sort))

    if [[ \${#R1_FILES[@]} -ne \${#R2_FILES[@]} ]]; then
        echo "❌ ERROR: R1/R2 count mismatch for ${sample_id} (\${#R1_FILES[@]} vs \${#R2_FILES[@]})" | tee -a "${LOG}" >&2
        exit 1
    fi

    FASTQ_ARGS=()
    for i in "\${!R1_FILES[@]}"; do
        FASTQ_ARGS+=("\${R1_FILES[\$i]}" "\${R2_FILES[\$i]}")
    done

    # --workflow nac    : matches SCRNA_KB_REF's index — spliced/unspliced/ambiguous
    #                     quantification for velocity, never "standard"
    # -x                : chemistry/technology string (no auto-detect in kb-python;
    #                     must match the sample's actual chemistry, e.g. 10XV4)
    # -m                : cap memory explicitly (kb-python defaults to only 4G for
    #                     non-standard workflows, too low at this scale — mirrors
    #                     Cell Ranger's --localmem)
    # -c1/-c2           : cDNA / nascent transcript-to-capture lists from kb ref,
    #                     required by the nac workflow to split spliced/unspliced
    # --filter bustools : knee-point empty-droplet filtering on the final matrix.
    #                     Reduces barcode count from millions to a realistic cell count.
    # --delete-bus      : intermediate BUS files (~30-40GB/sample) aren't needed 
    #                     downstream once the count matrix is built
    # --gene-names      : group counts by gene symbol instead of Ensembl ID, to
    #                     match Seurat's gene naming for the later merge step.
    # --h5ad            : used instead of --loom. loompy's write_loom() densifies
    #                     the sparse count matrix before writing (calls .toarray()),
    #                     which OOM'd on raw (millions-of-empty-droplets) FASTQ input
    #                     — tried to allocate 1.58 TiB for a 78348 x 2765824 dense array. 
    #                     AnnData's write_h5ad() keeps the matrix sparse (CSR) throughout,
    #                     avoiding that blowup entirely.

    kb count \
        --workflow nac \
        -i "${kb_index_dir}/index.idx" \
        -g "${kb_index_dir}/t2g.txt" \
        -x "${kb_technology}" \
        -c1 "${kb_index_dir}/cdna_t2c.txt" \
        -c2 "${kb_index_dir}/nascent_t2c.txt" \
        -o "${sample_id}" \
        -t "${task.cpus}" \
        --h5ad \
        --filter bustools \
        --delete-bus \
        --gene-names \
        "\${FASTQ_ARGS[@]}" \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: kb count failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

    # -m "${task.memory.toGiga()}G" \

    # kb-python writes the h5ad inside counts_filtered/ (or counts_unfiltered/ if
    # --filter fails to find a knee point) — surface it at the top level with a
    # sample-prefixed name for a clean, unambiguous merge input
    find "${sample_id}" -name "*.h5ad" -exec mv {} "${sample_id}.h5ad" \\;

    echo "✅ SUCCESS: kb count completed for ${sample_id}" >> "${LOG}"
    """
}
