



def genome_indexor(gtf_file, genome_file ):

    if gtf_file.endswith("gtf"):     # if file ends with gtf, use open()
        file_1 = open(gtf_file, "r")

    line = file_1.readline().decode("ascii").strip()
    while line.startswith("#"):
        line = file_1.readline().decode("ascii").strip()



    if genome_file.endswith("fa"):     # if file ends with fa, use open()
        file_2 = open(genome_file, "r")

    dna_id = []
    dna_seq = []

    line = file_2.readline().decode("ascii").strip()
    while line.startswith(">"):
        dna_id = dna_id.append(line)
        seq = ""
        line = file_2.readline().decode("ascii").strip()
        while not line.startswith(">"):
            seq = seq + line
            line = file_2.readline().decode("ascii").strip()
        dna_seq = dna_seq.append(seq)

    file_1.close()
    file_2.close()

    return (dna_id, dna_seq)


if f_read_file.endswith("fastq.gz") or f_read_file.endswith("fq.gz"):  # if file ends with .gz, use gzip.open()
    file_3 = gzip.open(f_read_file, "r")

    f_id = file_1.readline().decode("ascii").strip()  # read 1st line from file_1, store it as f_id

    while f_id.startswith("@"):  # loop will end once "end of file" is reached
        f_seq = file_1.readline().decode("ascii").strip()  # read 2nd line from file_1, store it as f_seq
        file_1.readline().decode("ascii").strip()
        f_qual = file_1.readline().decode("ascii").strip()  # read 4th line from file_1, store it as f_qual

        count = 0
        for i in range(0, len(dna_id)):
            if f_seq in dna_seq[i]:
                count = count+1

