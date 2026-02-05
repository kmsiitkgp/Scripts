#!/bin/bash -l

#$ -N CRISPR_align   # Set Job Name
#$ -l mem_free=96G            # Request memory
#$ -j y                       # Merge standard output and standard error
#$ -cwd                       # Set current working directory
#$ -t 1-110                   # Submit an array job with n tasks; n=number of samples = number of fastq files/2 if paired end

proj="CRISPR_Lena"
libraries=(CDH12KO CDH12ACT ImmuneKO ImmuneACT)
INPUT_DIR=$HOME/scratch/$proj/cutadapt_results
REF_DIR_BWAMEM2=$HOME/scratch/$proj/index_results/bwamem2
REF_DIR_BOWTIE2=$HOME/scratch/$proj/index_results/bowtie2
REF_DIR_BOWTIE=$HOME/scratch/$proj/index_results/bowtie

for library in ${libraries[*]}
do 

OUTPUT_DIR_BWAMEM2=$HOME/scratch/$proj/alignment_results/bwamem2/$library 
OUTPUT_DIR_BOWTIE2=$HOME/scratch/$proj/alignment_results/bowtie2/$library 
OUTPUT_DIR_BOWTIE=$HOME/scratch/$proj/alignment_results/bowtie/$library 

# Create an array of input files
input=($INPUT_DIR/*.gz)

# Use the SGE_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value
index=$(($SGE_TASK_ID-1))
taskinput=${input[$index]}
taskoutput_bwamem2=${taskinput//$INPUT_DIR/$OUTPUT_DIR_BWAMEM2}
taskoutput_bwamem2=${taskoutput_bwamem2//.*/}
taskoutput_bwamem2_sam=$taskoutput_bwamem2.bwamem2.sam
taskoutput_bwamem2_bam=$taskoutput_bwamem2.bwamem2.bam
taskoutput_bowtie2=${taskinput//$INPUT_DIR/$OUTPUT_DIR_BOWTIE2}
taskoutput_bowtie2=${taskoutput_bowtie2//.*/}
taskoutput_bowtie2_sam=$taskoutput_bowtie2.bowtie2.sam
taskoutput_bowtie2_bam=$taskoutput_bowtie2.bowtie2.bam
taskoutput_bowtie=${taskinput//$INPUT_DIR/$OUTPUT_DIR_BOWTIE}
taskoutput_bowtie=${taskoutput_bowtie//.*/}
taskoutput_bowtie_sam=$taskoutput_bowtie.bowtie.sam
taskoutput_bowtie_bam=$taskoutput_bowtie.bowtie.bam

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date    : $(date)"
echo "Job name      : $JOB_NAME"
echo "Job ID        : $JOB_ID"  
echo "SGE TASK ID   : $SGE_TASK_ID"
echo "Index         : $index"
echo "Task input 1  : $taskinput"
echo "Task input 2  : $taskinput2"
echo "Task output   : $taskoutput_bwamem2_sam"
echo "Task output   : $taskoutput_bowtie2_sam"
echo "Task output   : $taskoutput_bowtie_sam"
echo "Project       : $proj"
echo "Output Folder : $OUTPUT_DIR_BWAMEM2"
echo "Output Folder : $OUTPUT_DIR_BOWTIE2"
echo "Output Folder : $OUTPUT_DIR_BOWTIE"
echo "=========================================================="

# You need to add path for Conda
PATH=$HOME/miniconda3/bin/:$PATH
source activate MAGECK
PATH=$HOME/NGSTools/sambamba-1.0.1-linux-amd64-static/:$PATH

# Although paired end, in CRISPR screens, the read1 or read2 usually 
# contains the guide sequence. So, we do not need use both read1 and read2. 
# In Jinfen CRISPR, the guides are present in read1 only.

printf "\n\n***********Aligning reads using BWAMEM2******************\n\n"
cd $REF_DIR_BWAMEM2

SEED_LENGTH=5
INDEX_BASE_NAME=$proj.$library

# bwa-mem2 mem \
# -o $taskoutput_bwamem2_sam \
# -k $SEED_LENGTH \
# $INDEX_BASE_NAME \
# $taskinput

# bamtools stats -in $taskoutput_bwamem2_bam 
# sambamba flagstat $taskoutput_bwamem2_bam

printf "\n\n***********Aligning reads using BOWTIE2******************\n\n"
cd $REF_DIR_BOWTIE2

TRIM_5_PRIME=0
TRIM_3_PRIME=20
MAX_MISMATCH_IN_SEED=0
SEED_LENGTH=5
INDEX_BASE_NAME=$proj.$library
TEST_RUN_READS=10000
#-q indicates input are fastq files

# bowtie2 \
# --qupto $TEST_RUN_READS \
# -q \
# -5 $TRIM_5_PRIME \
# -3 $TRIM_3_PRIME \
# --phred33 \
# -N $MAX_MISMATCH_IN_SEED \
# -L $SEED_LENGTH \
# --no-1mm-upfront \
# --norc \
# --local \
# --very-sensitive \
# --seed 42 \
# -x $INDEX_BASE_NAME \
# $taskinput \
# -S $taskoutput_bowtie2_sam

# DO NOT USE --nosq option as MAGECK wont be able to count

bowtie2 \
-q \
-5 $TRIM_5_PRIME \
-3 $TRIM_3_PRIME \
--phred33 \
-N $MAX_MISMATCH_IN_SEED \
-L $SEED_LENGTH \
--no-1mm-upfront \
--norc \
--end-to-end \
--very-sensitive \
--seed 42 \
-x $INDEX_BASE_NAME \
$taskinput \
-S $taskoutput_bowtie2_sam
#--qupto $TEST_RUN_READS \

printf "\nConverting from SAM to BAM\n"
sambamba view \
--sam-input \
--format=bam \
--output-filename=$taskoutput_bowtie2_bam \
$taskoutput_bowtie2_sam

printf "\n\n***********Aligning reads using BOWTIE******************\n\n"
cd $REF_DIR_BOWTIE

TRIM_5_PRIME=0
TRIM_3_PRIME=20
MAX_MISMATCH_IN_SEED=0 #3
MAX_MISMATCH_IN_ALIGNMENT=0 #3
SEED_LENGTH=5
QUALITY_CUTOFF=70
INDEX_BASE_NAME=$proj.$library
#TEST_RUN_READS=100000
#-q indicates input are fastq files

# printf "\nAligning using bowtie n mode\n"
# bowtie \
# -q \
# -5 $TRIM_5_PRIME \
# -3 $TRIM_3_PRIME \
# --phred33-quals \
# -n $MAX_MISMATCH_IN_SEED \
# -l $SEED_LENGTH \
# -e $QUALITY_CUTOFF \
# --all \
# --best \
# --strata \
# -m 1 \
# --norc \
# --sam \
# -x $INDEX_BASE_NAME \
# $taskinput \
# $taskoutput_bowtie_sam 
# --qupto $TEST_RUN_READS \
# 2> bowtie-out.stderr

# DO NOT USE --nosq option as MAGECK wont be able to count

printf "\nAligning using bowtie v mode\n"
bowtie \
-q \
-5 $TRIM_5_PRIME \
-3 $TRIM_3_PRIME \
--phred33-quals \
-v MAX_MISMATCH_IN_ALIGNMENT \
--all \
--best \
--strata \
-m 1 \
--norc \
--sam \
-x $INDEX_BASE_NAME \
$taskinput \
$taskoutput_bowtie_sam 
#--qupto $TEST_RUN_READS \
#2> bowtie-out.stderr

printf "\nConverting from SAM to BAM\n"
sambamba view \
--sam-input \
--format=bam \
--output-filename=$taskoutput_bowtie_bam \
$taskoutput_bowtie_sam

done

# printf "\nMarking duplicates\n"
# samtools collate -O -u $taskoutput_bowtie_bam | samtools fixmate  -m -u - - | samtools sort  -u - | samtools markdup  - $taskoutput_bowtie_bam

# printf "\nStats Before Filtering\n"
# sambamba flagstat $taskoutput_bowtie_bam

# printf "\nStats After Filtering\n"
# sambamba flagstat $taskoutput_bowtie_bam

# printf "\nGenerating Index File\n"
# sambamba index $taskoutput_bowtie_bam

# NOTE: 
# Bowtie2 supports (i) gapped, ungapped (ii) local, global and (iii) single-end,paired-end alignment modes.
# Bowtie1 supports only ungapped alignment i.e. global alignment.
# Bowtie2 works best for reads that are greater than 50 bp.
# Bowtie1 works best for reads that are lesser than 50bp.
# bwamem2 is slower than bowtie2 but more accurate. It can work with 70bp to 1Mbp sequences.

# Read about bwamem2, bowtie2, bowtie
# https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
# https://bowtie-bio.sourceforge.net/manual.shtml
# https://www.youtube.com/watch?app=desktop&v=q0_9UAXkVr4
# https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#how-is-bowtie-2-different-from-bowtie-1
# https://hbctraining.github.io/Intro-to-ChIPseq-flipped/lessons/04_alignment_using_bowtie2.html
# By default, Bowtie2 will perform a global end-to-end read alignment, which aligns from the first to the last base of the read. This alignment is best for reads that have already been trimmed for quality and adapters (e.g. reads where nucleotide bases of poor quality or matching adapter sequences have been removed from the ends of the reads prior to alignment). However, Bowtie2 also has a local alignment mode, which, in contrast to end-to-end alignment, ignores portions at the ends of the reads that do match well to the reference genome. This is referred to as soft-clipping and allows for a more accurate alignment. The procedure can carry a small penalty for each soft-clipped base, but amounts to a significantly smaller penalty than mismatching bases. In contrast to trimming, which removes the unwanted sequence (hard-clipping), soft-clipping retains the soft-clipped base in the sequence and simply marks it. If you did not trim reads, use --local when runing Bowtie2. This will perform “soft-clipping” to ignore parts of the reads that are of low quality.

# To perform the Bowtie2 alignment, a genome index is required. The index is analagous to the index in a book. By indexing the genome, we have organized it in a manner that now allows for efficient search and retrieval of matches of the query (sequence read) to the genome. Bowtie2 indexes the genome with an FM Index based on the Burrows-Wheeler Transform method to keep memory requirements low for the alignment process.

# ALIGNMENT MODE:
# -n : number of mismatches in alignment between seed and reference. In this alignment mode, the sum of quality scores at mismatch positions in entire alignment must be < 70 (default). This is more lenient than -v alignment mode.
# -v : number of mismatches in alignment between read and reference. There is no seed region used and quality scores are ignored (-l and -e). Only use this mode if all reads are high quality.

# N ALIGNMENT MODE SPECIFIC PARAMETERS:
# -l : seed length i.e. number of bases at 5' end of read that is used as seed. Used only in -n aligment mode.
# -e : Maximum permitted total of quality values at all mismatched read positions throughout the entire alignment, not just in the “seed” (default 70). Used only in -n aligment mode.

# REPORTING ALIGNMENTS:
# -k : number of valid alignments to report for each read (default 1). Each read can have many valid alignments, some with 1 mismatch, some with 2 mismatches etc. If you do not specify --best, the 1st valid alignment will be reported which will likely not be the best alignment.
# -a : all valid alignments for each read will be reported
# --best : reported alignment(s) are “best” in terms of the number of mismatches and are reported in best-to-worst order
# --strata : number of mismatches in the “seed” region (-n alignment mode) or entire alignment (-v alignment mode). If --strata is specified, --best must also be specified. A read can have 5 valid alignments: 2 alignments with 1 mismatch and 3 alignments with 2 mismatches. --strata will report only the 2 alignments with 1 mismatch as they have minimum mismatch among valid alignments
# -m : if read has more than m "reportable" alignments, none of the alignments are reported

# Scenarios:
# -l=5 : The first 5 bases of the read is the seed region. 
# -n=0 : If the seed region doesnt have 100% match anywhere in the reference, the read will be unmapped.
# -n=2 : All alignments where (i) number of mismatch within seed region <=2 and (ii) sum of quality score at mismatch positions in entire alignment < 70 will be valid alignments. There can be more than 2 mismatch in entire alignment provided (ii) is valid
# -v=2 : There can be no more than 2 mismatch in the entire alignment.


# First, in our CRISPR library, the guides are of different length (17 to 25 bases). 
# Second, the adapter sequence is not identical in all reads. So, we cannot remove adapter sequence from 100% of reads.
# Cutadapt has default error rate of 10% (e=0.1) i.e. if length of adapter=10, then it can find adapters with upto 1 mismatch.
# Also, reads will be of different length after adapter trimming.
# Third, even if we trim all reads to same size (say 30bp) using -cut option, the starting position of the guides in each
# read are different. So, the seed region in bowtie/bowtie2 needs to be minimum (lowest -l in bowtie is 5) and allow for
# maximum mismatch (highest n in bowtie is 3) so that we are able to map most reads. So, if we set -n=3, -l=5, we will be able to map almost all reads and by setting --best and --strata, we will be able to keep the top alignments with the least mismatch in the seed. We can also specify -m=1 to keep only reads having 1 top alignment and ignore reads that have 2 or more valid alignments with least mismatch. Theoretically, we cannot use -v option as it allows only 3 mismatches maximum in whole read and if most of our reads are 25 bases while many guides are less than 20 bases, the number of mismatches for whole read will be grater than 3. HOWEVER, I noticed that -v option gives slightly better alignment percentage as compared to -n option !!! Also, bowtie2 end-to-end works way better than bowtie2 local mode, however, bowtie2 end-to-end mode isnt as good as bowtie -v alignment mode.