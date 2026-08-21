// =============================================================================
// PROCESS: RNA_FETCH_ENSEMBL_VERSIONS
// =============================================================================
// Purpose: Queries Ensembl REST API to get the latest assembly and release metadata
//
// What it does:
//   - Fetches current genome assembly name (e.g., GRCh38, GRCm39) via REST API
//   - Fetches current Ensembl release number (e.g., 115)
//   - Constructs expected FASTA and GTF filenames from those values
//   - For Xenograft: produces a synthetic entry reusing Human files but with a
//     composite version tag (GRCh38.115_GRCm39.115) to name the index distinctly
//   - Writes pipe-delimited metadata to metadata.txt for Nextflow to capture
//
// Why runtime fetch instead of hardcoding?
//   - Ensures the pipeline always uses the latest Ensembl release automatically
//   - No stale version strings to update manually between pipeline runs

// =============================================================================

process RNA_FETCH_ENSEMBL_VERSIONS {

    tag "Fetching latest Ensembl metadata for ${species}"
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
    // Runtime values computed in Bash (via curl) must be written to a file
    // so Nextflow can capture and emit them as channel items
    output:
    path("metadata.txt"),    emit: meta

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    """
    # Fetch assembly name and release number for both species
    # Both are always fetched; used selectively based on species below
    H_ASM=\$(curl -s --connect-timeout 10 --max-time 20 \
        https://rest.ensembl.org/info/assembly/homo_sapiens?content-type=application/json \
        | grep -o '"default_coord_system_version":"[^"]*"' | sed 's/.*:"//;s/"//')

    M_ASM=\$(curl -s --connect-timeout 10 --max-time 20 \
        https://rest.ensembl.org/info/assembly/mus_musculus?content-type=application/json \
        | grep -o '"default_coord_system_version":"[^"]*"' | sed 's/.*:"//;s/"//')

    H_REL=\$(curl -s https://rest.ensembl.org/info/software?content-type=application/json \
        | grep -o '"release":[0-9]*' | head -1 | grep -o '[0-9]\\+')

    M_REL=\$(curl -s https://rest.ensembl.org/info/software?content-type=application/json \
        | grep -o '"release":[0-9]*' | head -1 | grep -o '[0-9]\\+')

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
        # Xenograft uses Human FASTA + GTF for alignment, but gets a composite
        # version tag so its index directory is named distinctly from Human's
        FASTA="Homo_sapiens.\${H_ASM}.dna.primary_assembly.fa"
        GTF="Homo_sapiens.\${H_ASM}.\${H_REL}.gtf"
        ASSEMBLY=\${H_ASM}
        RELEASE=\${H_REL}
        VERSION=\${H_ASM}.\${H_REL}_\${M_ASM}.\${M_REL}
    fi

    # Write pipe-delimited metadata; parsed in rnaseq.nf via .splitText()
    echo "${species}|\$FASTA|\$GTF|\$VERSION|\$ASSEMBLY|\$RELEASE" > metadata.txt
    """
}
