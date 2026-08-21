#!/usr/bin/env nextflow
nextflow.enable.dsl=2        // Enable DSL2 syntax (modern Nextflow with explicit workflow blocks)

// ===============================
// INCLUDE WORKFLOWS
// ===============================
include { RNASEQ }       from './workflows/rnaseq.nf'
include { SCRNASEQ }     from './workflows/scrnaseq.nf'
include { CTDNA }        from './workflows/ctdna.nf'
include { TCGA }         from './workflows/tcga.nf'
include { RENAME }       from './workflows/rename.nf'

// ===============================
// DISPATCH
// ===============================
workflow {

    // 1. If the user explicitly requested JUST a rename operation, run it and stop
    if (params.run_mode == "rename") {
        RENAME()
    }
    // 2. Otherwise, route to the full project pipelines (which call renaming internally)
    else if (params.expt == "RNASeq") {
        RNASEQ()
    }
    else if (params.expt == "scRNASeq") {
        SCRNASEQ()
    }
    else if (params.expt == "TCGA") {
        TCGA()
    }
    else if (params.expt == "ctDNA") {
        CTDNA()
    }
    else {
        error "Unknown experiment type: ${params.expt}"
    }
}
