// =============================================================================
// PROCESS: SCRNA_RUN_VELOCITY
// =============================================================================
// Purpose: Runs RNA velocity on one lineage's velocity_input.h5ad using BOTH
//          scVelo (dynamical model) and veloVI (scvi-tools), sharing one
//          preprocessing pass. Each method is independently wrapped in
//          run_velocity.py — a veloVI failure (heavier/more failure-prone:
//          GPU/memory, convergence) never blocks scVelo's already-written
//          output, and either method can be commented out in the .py file
//          with NO change needed here (both outputs are `optional: true`).
//
// Runs once per lineage (parallel), independent per lineage.
// =============================================================================

process SCRNA_RUN_VELOCITY {

    tag "RNA velocity (scVelo + veloVI): ${label}"
    label 'process_heavy'
    label 'omics_py'

    publishDir = [
        [path: { "${params.proj_dir()}/${params.species}/06.scVelo/${label}" },        mode: 'copy', pattern: "*_scvelo.h5ad"],
        [path: { "${params.proj_dir()}/${params.species}/06.scVelo/${label}" },        mode: 'copy', pattern: "*_velovi.h5ad"],
        [path: { "${params.proj_dir()}/${params.species}/06.scVelo/${label}/Plots" },  mode: 'copy', pattern: "*.png"],
        [path: { "${params.proj_dir()}/${params.species}/Logs" },                        mode: 'copy', pattern: "*.log"]
    ]

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(label), path(velocity_input_h5ad)
    // label               : base lineage label, e.g. "Epithelial"
    // velocity_input_h5ad : SCRNA_PREP_VELOCITY_INPUT.out.velocity_input for
    //                       this lineage — barcode-subset, CellType/cluster/
    //                       UMAP already attached

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    tuple val(label), path("${label}_scvelo.h5ad"),  optional: true,  emit: scvelo_h5ad
    tuple val(label), path("${label}_velovi.h5ad"),  optional: true,  emit: velovi_h5ad
    // Both optional: true — if either method fails, or is commented out in
    // run_velocity.py, the process still succeeds and publishes whatever DID
    // get written. No .nf change needed to disable one method.
    path("*.png"),                    optional: true,                 emit: plots
    path("SCRNA_RUN_VELOCITY.${label}.log"),                               emit: log

    // =================================================================================
    // EXECUTION
    // NOTE: run_velocity.py runs shared preprocessing (filter/normalize +
    // moments) once, then scVelo and veloVI each independently in their own
    // try/except block — see script docstring for the exact rationale.
    // =================================================================================
    script:

    def script_path = "${workflow.projectDir}/modules"
    def LOG          = "SCRNA_RUN_VELOCITY.${label}.log"

    """
    python3 "${script_path}/SCRNA/SCRNA_run_velocity.py" \
        --input "${velocity_input_h5ad}" \
        --label "${label}" \
        --output-dir "." > "${LOG}" 2>&1

    EXIT_CODE=\$?
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "❌ ERROR: SCRNA_RUN_VELOCITY failed for ${label}" | tee -a "${LOG}" >&2
        exit \$EXIT_CODE
    fi
    echo "✅ SUCCESS: SCRNA_RUN_VELOCITY completed for ${label}" >> "${LOG}"
    """
}
