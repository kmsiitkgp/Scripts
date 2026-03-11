process FETCH_ENSEMBL_VERSIONS {

    tag "Fetching latest ENSEMBL version"
    executor 'local'

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    val(species)

    // =================================================================================
    // OUTPUT
    // =================================================================================
    output:
    // Runtime values are computed in Bash (via curl).
    // They must be written to a file so Nextflow can capture and emit them.
    path("metadata.txt"),                         emit: meta

    // =================================================================================
    // EXECUTION
    // =================================================================================

    script:

    """
    # Fetch Assembly and Release for human and mouse
    H_ASM=\$(curl -s --connect-timeout 10 --max-time 20 https://rest.ensembl.org/info/assembly/homo_sapiens?content-type=application/json | grep -o '"default_coord_system_version":"[^"]*"' | sed 's/.*:"//;s/"//')
    M_ASM=\$(curl -s --connect-timeout 10 --max-time 20 https://rest.ensembl.org/info/assembly/mus_musculus?content-type=application/json | grep -o '"default_coord_system_version":"[^"]*"' | sed 's/.*:"//;s/"//')
    H_REL=\$(curl -s https://rest.ensembl.org/info/software?content-type=application/json | grep -o '"release":[0-9]*' | head -1 | grep -o '[0-9]\\+')
    M_REL=\$(curl -s https://rest.ensembl.org/info/software?content-type=application/json | grep -o '"release":[0-9]*' | head -1 | grep -o '[0-9]\\+')

    # We use a trick here: output the values to txt file or stdout and capture them

    if [[ "${species}" == "Human" ]]; then
        GTF="Homo_sapiens.\${H_ASM}.\${H_REL}.gtf"
        FASTA="Homo_sapiens.\${H_ASM}.dna.primary_assembly.fa"
        ASSEMBLY=\${H_ASM}
        RELEASE=\${H_REL}
        VERSION=\${H_ASM}.\${H_REL}
    elif [[ "${species}" == "Mouse" ]]; then
        GTF="Mus_musculus.\${M_ASM}.\${M_REL}.gtf"
        FASTA="Mus_musculus.\${M_ASM}.dna.primary_assembly.fa"
        ASSEMBLY=\${M_ASM}
        RELEASE=\${M_REL}
        VERSION=\${M_ASM}.\${M_REL}
    elif [[ "${species}" == "Xenograft" ]]; then
        GTF="Homo_sapiens.\${H_ASM}.\${H_REL}.gtf"
        FASTA="Homo_sapiens.\${H_ASM}.dna.primary_assembly.fa"
        ASSEMBLY=\${H_ASM}
        RELEASE=\${H_REL}
        # Create a composite version string but use Human GTF,FASTA,ASSEMBLY,RELEASE
        VERSION=\${H_ASM}.\${H_REL}_\${M_ASM}.\${M_REL}
    fi

    # Write to a single metadata file
    echo "${species}|\$FASTA|\$GTF|\$VERSION|\$ASSEMBLY|\$RELEASE" > metadata.txt

    """
}