#!/usr/bin/env nextflow
nextflow.enable.dsl=2        // Enable DSL2 syntax (modern Nextflow with explicit workflow blocks)

// ===============================
// INCLUDE WORKFLOWS
// ===============================
include { RNASEQ }       from './workflows/rnaseq'
include { SCRNASEQ }     from './workflows/scrnaseq'

// ===============================
// DISPATCH
// ===============================
workflow {

    if (params.expt == 'RNASeq') {
        RNASEQ()
    }
    else if (params.expt == 'scRNASeq') {
        SCRNASEQ()
    }
    else {
        error "Unknown experiment type: ${params.expt}"
    }
}
