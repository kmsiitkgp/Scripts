// =============================================================================
// PROCESS: SCRNA_FETCH_10X_VERSIONS
// =============================================================================
// Purpose: Provides the assembly/release metadata used to build 10X Genomics
//          reference filenames and version tags for each supported species.
//
// What it does:
//   - Returns the hardcoded genome assembly name (e.g., GRCh38, GRCm39) and
//     10X reference release (e.g., 2024-A) for the requested species
//   - Constructs expected FASTA and GTF filenames from those values
//   - For Xenograft: 10X Genomics publishes its OWN dedicated combined
//     human+mouse reference package (it does NOT reuse the plain Human
//     reference) - the composite version tag here just names that
//     dedicated reference distinctly from the standalone Human/Mouse ones
//   - Writes pipe-delimited metadata to metadata.txt for Nextflow to capture
//
// Why hardcoded instead of a live API/website lookup?
//   - 10X Genomics reference releases change infrequently and are not
//     exposed via a stable REST API - these values are maintained by hand
//   - When 10X publishes a new reference version, update H_ASM/M_ASM/X_ASM
//     and H_REL/M_REL/X_REL below (and the corresponding tarball URL/name
//     in SCRNA_DOWNLOAD_10X_REF.nf, which is updated independently)
//
// Species allow-list:
//   - Only "Human", "Mouse", "Xenograft" are supported. This is also
//     validated up front in scrnaseq.nf before any process is launched;
//     the else-branch below is a defensive second check in case this
//     process is ever invoked directly with an unexpected value.
// =============================================================================

process SCRNA_FETCH_10X_VERSIONS {

    tag "Fetching 10X reference metadata for ${species}"
    label 'process_low'

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    val(species)
    // "Human", "Mouse", or "Xenograft"

    // =================================================================================
    // OUTPUT
    // =================================================================================
    // Values are computed in Bash and written to a file so Nextflow can
    // capture and emit them as channel items.
    output:
    path("metadata.txt"),    emit: meta

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    """
    set -euo pipefail

    # Hardcoded assembly and release based on latest 10X Genomics reference packages.
    # Update by hand when 10X publishes a new release.
    H_ASM="GRCh38"
    M_ASM="GRCm39"
    X_ASM="GRCh38_GRCm39"
    H_REL="2024-A"
    M_REL="2024-A"
    X_REL="2024-A"

    # Clear the file at the start to ensure a fresh output for this task
    > metadata.txt

    if [[ "${species}" == "Human" ]]; then
        FASTA="Homo_sapiens.\${H_ASM}.dna.primary_assembly.fa"
        GTF="Homo_sapiens.\${H_ASM}.\${H_REL}.gtf"
        ASSEMBLY=\${H_ASM}
        RELEASE=\${H_REL}
        VERSION=\${H_ASM}.\${H_REL}
    elif [[ "${species}" == "Mouse" ]]; then
        FASTA="Mus_musculus.\${M_ASM}.dna.primary_assembly.fa"
        GTF="Mus_musculus.\${M_ASM}.\${M_REL}.gtf"
        ASSEMBLY=\${M_ASM}
        RELEASE=\${M_REL}
        VERSION=\${M_ASM}.\${M_REL}
    elif [[ "${species}" == "Xenograft" ]]; then
        # Dedicated 10X combined human+mouse reference - NOT a reuse of the
        # plain Human FASTA/GTF above.
        FASTA="GRCh38_and_GRCm39.\${X_ASM}.dna.primary_assembly.fa"
        GTF="GRCh38_and_GRCm39.\${X_ASM}.\${X_REL}.gtf"
        ASSEMBLY=\${X_ASM}
        RELEASE=\${X_REL}
        VERSION=\${X_ASM}.\${X_REL}
    else
        echo "❌ ERROR: Unsupported species '${species}'. Must be one of: Human, Mouse, Xenograft." >&2
        exit 1
    fi

    # Write pipe-delimited metadata; parsed in scrnaseq.nf via .splitText()
    echo "${species}|\$FASTA|\$GTF|\$VERSION|\$ASSEMBLY|\$RELEASE" >> metadata.txt

    """
}
