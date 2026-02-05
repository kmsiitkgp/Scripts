# https://genestack.com/blog/2016/07/12/choosing-a-reference-genome/
# https://dnatech.genomecenter.ucdavis.edu/faqs/when-should-i-trim-my-illumina-reads-and-how-should-i-do-it/
#https://sourceforge.net/p/mageck/wiki/Home/#:~:text=code%20on%20Latch!-,Model%2Dbased%20Analysis%20of%20Genome%2Dwide%20CRISPR%2DCas9%20Knockout,screens%20(or%20GeCKO)%20technology.


# Check Linux distro from os-release file or using lsb-release command
cat /etc/os-release
lsb_release -a

# View content between two line numbers
zcat  ~/common/scRNASeq_Simon/B6A-GEX13_S3_L002_R2_001.fastq.gz | awk 'NR >= 508461070 && NR <= 508461080'  

# To find read length i.e. number of bases in each read. If read length is 150, it prints 151
zcat ~/scratch/Neeraj/raw_reads/NA13Ca1R_1.fq.gz | head -2 | awk 'END {print}'| wc -c

# Print first 5 reads alone. Set head -5*4
zcat ~/scratch/WES_Lena/cutadapt_results/SCR1_CKDN240000876-1A_H2YYGDSXC_L3_1.fq.gz | head -20 | awk 'NR % 4 == 2'

# Count number of reads in fastq file. Unlike sam and bam files which have headers, fastq file contain 1 read every 4 lines
zcat ~/scratch/WES_Lena/cutadapt_results/SCR1_CKDN240000876-1A_H2YYGDSXC_L3_1.fq.gz | wc -l
# Divide the output by 4 to get number of reads

# Count number of reads in sam file. sambamba can ONLY count in bam files.
samtools view --count ~/scratch/WES_Lena/alignment_results/SCR1_CKDN240000876-1A_H2YYGDSXC_L3_1.sam

# Count number of reads in bam file.
sambamba view --count ~/scratch/WES_Lena/alignment_results/SCR1_CKDN240000876-1A_H2YYGDSXC_L3_1.sorted.bam

# Extract exons from GTF file and save it to annotation.txt
cat $HOME/scratch/Hany/ref_genome/Mus_musculus.GRCm38.102.gtf | awk '{if($3 == "exon") print $1,$3,$4,$5,$7,$20,$24}' | sed 's/"//g' | sed 's/;//g' > annotation.txt

# To  print reads that have adapters:
zcat ~/scratch/Neeraj/raw_reads/NA13R181_1.fq.gz | head -1000000 | awk /AGATCGGAAGAGCACACGTCTGAACTCCAGTCA/ | head -10
zcat ~/scratch/Neeraj/raw_reads/NA13R181_1.fq.gz | awk /AGATCGGAAGAGCACACGTCTGAACTCCAGTCA/ | wc -l
# zcat $HOME/common/$proj/*.fastq.gz | grep -n TTTGTTGTCTCGCTCA | head -15
zcat /common/theodorescudlab/Sequencing_Data/CRISPR_Prince/DTKPA1/*R1*.gz | awk /GTATCCCTTGGAGAACCACCTTGTTG/ | wc -l

zcat /common/theodorescudlab/CRISPR_Jinfen/Day0_2D_Rep1_1.fq.gz | awk /CTCCAAGGGATACTTATAGTCTCAAAACACACAATTACTTTACAGTTAGGGTGAG/ | wc -l

# To check if trimming was successful and adapters removed in reads
zcat ~/scratch/Neeraj/trimmed_reads/P_NA13R181_1.fq.gz | awk /AGATCGGAAGAGCACACGTCTGAACTCCAGTCA/ | wc -l

# To read the sequence of particular readID if cpp program gives segmentation fault etc
zcat ~/scratch/Neeraj/raw_reads/NA13Ca1R_1.fq.gz | awk '/A00261:107:HHKJ2DMXX:2:1101:27163:30405/{print;getline;print;}'

# To find first 2 line numbers where string is present
cat ~/scratch/Neeraj/ref_genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa | awk '/Na/{print NR}' | head -2

# To get sequence between 2 lines and blast them to verify if seq matches with gene in gtf file
cat -n ~/scratch/Neeraj/ref_genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa | awk '{if((NR>53558)&&(NR<53563)) print}'

# I am pretty sure that your guess is correct and HTseq counts the total number of multimapping alignments rather the reads.
# You can check it by counting the number of multimapping lines yourself, e.g. with
# awk 'substr($1,1,1)!="@" && substr($12,6)>1 {n++} END {print n}' ~/scratch/RNASeq/Mukta_NA13_MB49/STAR_alignment_results/MB49_OE1Aligned.sortedByCoord.out.bam
# This should be equal to the HTseq number for "alignment_not_unique".

# Count number of chr instances before replacement
cat ~/projects/WES/S07604514_hg38/*Covered_list | grep -o chr | wc -l

# Replace chr in a file
cat ~/projects/WES/S07604514_hg38/*Covered_list | awk '{gsub(/chr/,""); print}' > Covered_list

# Count number of chr instances after replacement
cat Covered_list | grep -o chr | wc -l
cat ~/projects/WES/S07604514_hg38/*Covered_list | grep -o chr | wc -l

# It is HIGHLY RECOMMENDED to run one qsub at a time although it is ok to qsub 04,05 and 06 simultaneously by sh (running) this script
# Below all jobs will be submitted one after another WITHOUT waiting for previous job to finish.
# This is fine for FastQC, FAT_Trimming and STAR_Genome_Indexing as they dont depend on each other.
# STAR_Alignment however, depends on results of STAR_Genome_Indexing and FAT_Trimming. So, best to run it separately on console.
# Similarly, HTSeq depends on results of STAR_Alignment. So, best to run it separately on console. 

#qsub ~/projects/RNASeq/05_RNASeq_FAT_Trimming.sh  #TRIMMING READS NOT REQUIRED ANYMORE
qsub ~/projects/RNASeq/02_RNASeq_FastQC.sh
qsub ~/projects/RNASeq/03_RNASeq_STAR_Genome_Indexing.sh
qsub ~/projects/RNASeq/04_RNASeq_STAR_Alignment.sh
qsub ~/projects/RNASeq/05_RNASeq_HTSeq_Read_Counting.sh
# Run 09_RNA_Seq_DESeq2.R in RStudio, not UNIX.

# Check number of reads in bam file and other statistics
~/NGSTools/sambamba-0.8.2-linux-amd64-static  flagstat ~/scratch/ChIPSeq/STAR_alignment_results/ChIP_AR_SRR11467745.fastqAligned.sortedByCoord.out.bam

# The output of the above command is as follows:
# 29642046 + 0 in total (QC-passed reads + QC-failed reads)  <----TOTAL NUMBER OF ALIGNMENTS (PRIMARY +SECONDARY) 
# 4519256 + 0 secondary                                      <----NUMBER OF SECONDARY ALIGNMENTS
# 0 + 0 supplementary
# 0 + 0 duplicates
# 29642046 + 0 mapped (100.00%:N/A)
# 0 + 0 paired in sequencing
# 0 + 0 read1
# 0 + 0 read2
# 0 + 0 properly paired (N/A:N/A)
# 0 + 0 with itself and mate mapped
# 0 + 0 singletons (N/A:N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

# PRIMARY ALIGNMENT for a read is the alignment with best MAPQ score. So, there will be only 1 PRIMARY ALIGNMENT for every mapped read irrespective of unique or multimapped read.
# SECONDARY ALIGNMENT for a read consists of all possible alignments. Unique mapped reads will have NO SECONDARY ALIGNMENT. Multimapped reads will have 1 or more SECONDARY ALIGNMENT.
# NUMBER OF PRIMARY ALIGNMENTS = 1 best alignment for each unique mapped reads + 1 best aligment for each multimapped read.
# In our case, NUMBER OF PRIMARY ALIGNMENTS = 29642046 - 4519256 = 25122790 = 23485281 unique mapped reads (see STAR output below) + 1637509 best aligned read for multimapped reads (see STAR output below)
# There are multiple ways to extract these UNIQUE READS. MAPQ=255 or NH=1 can be used to extract these unique reads from bam file.
# By default, STAR doesnt output unampped reads

# The STAR output from "SRR11467745.fastqLog.final.out" file is as follows:
                                 # Started job on |	Feb 10 17:16:54
                             # Started mapping on |	Feb 10 17:21:08
                                    # Finished on |	Feb 10 17:52:02
       # Mapping speed, Million of reads per hour |	53.33

                          # Number of input reads |	27463731       <------TOTAL READS
                      # Average input read length |	50
                                    # UNIQUE READS:
                   # Uniquely mapped reads number |	23485281  	   <------UNIQUELY MAPPED READS   
                        # Uniquely mapped reads % |	85.51%
                          # Average mapped length |	49.91
                       # Number of splices: Total |	49
            # Number of splices: Annotated (sjdb) |	49
                       # Number of splices: GT/AG |	43
                       # Number of splices: GC/AG |	4
                       # Number of splices: AT/AC |	0
               # Number of splices: Non-canonical |	2
                      # Mismatch rate per base, % |	0.15%
                         # Deletion rate per base |	0.00%
                        # Deletion average length |	1.00
                        # Insertion rate per base |	0.00%
                       # Insertion average length |	1.42
                             # MULTI-MAPPING READS:
        # Number of reads mapped to multiple loci |	1637509		   <------MULTI MAPPED READS
             # % of reads mapped to multiple loci |	5.96%
        # Number of reads mapped to too many loci |	470540         <------MULTI MAPPED READS
             # % of reads mapped to too many loci |	1.71%
                                  # UNMAPPED READS:
  # Number of reads unmapped: too many mismatches |	0
       # % of reads unmapped: too many mismatches |	0.00%
            # Number of reads unmapped: too short |	932513         <------UNMAPPED READS (These reads are not output in SAM or BAM files by default)
                 # % of reads unmapped: too short |	3.40%
                # Number of reads unmapped: other |	937888         <------UNMAPPED READS (These reads are not output in SAM or BAM files by default)
                     # % of reads unmapped: other |	3.42%
                                  # CHIMERIC READS:
                       # Number of chimeric reads |	0
                            # % of chimeric reads |	0.00%

# This doesnt correctly identify unique mapped reads. It identifies all primary alignments which is NOT useful for ChIPSeq analysis.
~/NGSTools/sambamba-0.8.2-linux-amd64-static view --filter "not (unmapped or duplicate or chimeric or supplementary or secondary_alignment)" --format bam --with-header --valid --show-progress --output-filename ~/test.bam ~/scratch/ChIPSeq/STAR_alignment_results/ChIP_AR_SRR11467745.fastqAligned.sortedByCoord.out.bam

~/NGSTools/sambamba-0.8.2-linux-amd64-static  flagstat ~/test.bam

#This correctly identifies unique mapped reads
~/NGSTools/sambamba-0.8.2-linux-amd64-static view --filter "[NH] == 1" --format bam --with-header --valid --show-progress --output-filename ~/test1.bam ~/scratch/ChIPSeq/STAR_alignment_results/ChIP_AR_SRR11467745.fastqAligned.sortedByCoord.out.bam

~/NGSTools/sambamba-0.8.2-linux-amd64-static  flagstat ~/test1.bam

#This correctly identifies unique mapped reads
~/NGSTools/sambamba-0.8.2-linux-amd64-static view --filter "mapping_quality == 255" --format bam --with-header --valid --show-progress --output-filename ~/test2.bam ~/scratch/ChIPSeq/STAR_alignment_results/ChIP_AR_SRR11467745.fastqAligned.sortedByCoord.out.bam

~/NGSTools/sambamba-0.8.2-linux-amd64-static  flagstat ~/test2.bam

samtools flagstat ~/scratch/ChIPSeq/STAR_alignment_results/ChIP_AR_SRR11467745.fastqAligned.sortedByCoord.out.bam
