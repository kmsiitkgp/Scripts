// =============================================================================
// PROCESS: HMMCOPY_REFWIG
// =============================================================================
// Purpose: Generate reusable GC content and mappability WIG files from
//          reference fasta using HMMcopy utilities gcCounter and mapCounter
//
// What it does:
//   - Computes GC content in fixed-size bins across the reference genome
//   - Computes mappability scores in fixed-size bins across the reference genome
//   - Outputs are genome-build specific and bin-size specific
//   - Run once per genome build; stored in ref_dir via storeDir and reused
//     across all ctDNA projects
//   - Handles both .fa and .fa.gz input — decompresses on the fly if needed
//
// Typical resources: 4GB RAM, 30-60 minutes on 4 cores
// =============================================================================

process HMMCOPY_REFWIG {

    tag "HMMcopy refwig: ${fasta.name}"
    label 'process_medium'

    storeDir { "${params.ref_dir}/ichorcna" }
    publishDir { "${params.proj_dir()}/Logs" }, mode: 'copy', pattern: "*.log"

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    path(fasta)
    // fasta    : Reference genome FASTA file — accepts .fa or .fa.gz
    //            Must be the same reference used for BAM alignment to ensure
    //            chromosome names match exactly (critical for readCounter)
    val(bin_size)
    // bin_size : Genomic bin size in base pairs (e.g. 1000000 = 1Mb)
    //            Must match bin_size used in HMMCOPY_READCOUNTER and ICHORCNA_RUN

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    path("ichorcna_${fasta.name.replaceAll(/\.fa(\.gz)?$/, '').replaceAll(/\.fasta(\.gz)?$/, '')}_${bin_size}"),    emit: refwig_dir
    path("HMMCOPY_REFWIG.log"),                                                                                      emit: log

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:
    def genome_name = fasta.name.replaceAll(/\.fa(\.gz)?$/, "").replaceAll(/\.fasta(\.gz)?$/, "")
    def fa_name     = fasta.name.replaceAll(/\.gz$/, "")   // strips .gz suffix if present
    def out_dir     = "ichorcna_${genome_name}_${bin_size}"
    def GC_WIG      = "${out_dir}/gc.wig"
    def MAP_WIG     = "${out_dir}/map.wig"
    def LOG         = "HMMCOPY_REFWIG.log"

    """
    mkdir -p ${out_dir}

    # Decompress fasta if gzipped — keeps original .fa.gz untouched
    # gunzip -c writes to stdout so the source file is never modified
    if [[ ${fasta} == *.gz ]]; then
        echo "Decompressing ${fasta}..." >> ${LOG}
        gunzip -c ${fasta} > ${fa_name}
        FASTA=${fa_name}
    else
        FASTA=${fasta}
    fi

    # Generate GC content WIG file
    # gcCounter slides across the genome in fixed bins and computes
    # the fraction of G+C bases in each bin
    gcCounter \
        -w ${bin_size} \
        \$FASTA \
        > ${GC_WIG} \
        2>> ${LOG} \
        || { echo "❌ ERROR: gcCounter failed for \$FASTA" | tee -a ${LOG} >&2; exit 1; }

    echo "✅ SUCCESS: gcCounter completed — ${GC_WIG}" >> ${LOG}

    # Generate mappability WIG file
    # mapCounter slides across the genome in fixed bins and computes
    # the mappability score in each bin from the reference fasta
    mapCounter \
        -w ${bin_size} \
        \$FASTA \
        > ${MAP_WIG} \
        2>> ${LOG} \
        || { echo "❌ ERROR: mapCounter failed for \$FASTA" | tee -a ${LOG} >&2; exit 1; }

    echo "✅ SUCCESS: mapCounter completed — ${MAP_WIG}" >> ${LOG}
    """
}

// =============================================================================
// QUICK REFERENCE
// =============================================================================
//
// Output directory (stored in ref_dir/ichorcna/):
//   ichorcna_{genome_name}_{bin_size}/
//     gc.wig  : GC content per bin across genome
//     map.wig : Mappability score per bin across genome
//
// Example with GRCh38Decoy-ss-ctrls-v24.fa.gz at bin_size 1000000:
//   ichorcna_GRCh38Decoy-ss-ctrls-v24_1000000/
//     gc.wig
//     map.wig
//
// Reusing across projects:
//   After first run, provide paths directly in project_info.yaml:
//     gc_wig:  /home/kailasamms/resources/biomodal/ichorcna/ichorcna_GRCh38Decoy-ss-ctrls-v24_1000000/gc.wig
//     map_wig: /home/kailasamms/resources/biomodal/ichorcna/ichorcna_GRCh38Decoy-ss-ctrls-v24_1000000/map.wig
//   This skips HMMCOPY_REFWIG entirely on subsequent projects.
//
// Input fasta:
//   Must be the SAME fasta used for BAM alignment — chromosome names must match
//   BAM headers exactly or readCounter will produce empty/wrong WIG files.
//   For Biomodal DUET data: use GRCh38Decoy-ss-ctrls-v24.fa.gz
//
// Common issues:
//   - gcCounter/mapCounter not found → check hmmcopy container is loaded
//   - Empty WIG file                 → chromosome names don't match BAM headers
//   - Process skipped on rerun       → storeDir found existing output (expected behavior)
// =============================================================================
