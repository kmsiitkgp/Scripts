#!/bin/bash -l

#$ -N CRISPR_count              # Set Job Name
#$ -l mem_free=64G              # Request memory
#$ -j y                         # Merge standard output and standard error
#$ -cwd                         # Set current working directory
#$ -t 1-1                       # Submit an array job with 1 task.

# IMPORTANT:  If you align using bowtie, bowtie2 or bwamem2, and use MAGeCK to
# calculate counts from bam file, you will notice that the total alignments for 
# each sample from these aligners DO NOT MATCH with total counts for each sample
# calculated by MAGeCK. I have used R to calculate counts from sam files directly 
# and the results are identical to counts from MAGeCK. So, the error is with these
# aligners, not with MAGeCK.

# IMPORTANT: If you compare the count.txt output of bowtie, bowtie2, MAGeCK and KMS,
# you will notice bowtie and bowtie2 fail to find reads for some guides while MAGeCK and 
# KMS are able to find map such reads. So, RECOMMENDED using MAGeCK or KMS for aligning 
# guides instead of bowtie, bowtie2 or bwamem2.

# IMPORTANT: You dont need to perform any normalization at count step. The normalized counts
# generated at this step are not used in the next step. The next step uses the raw counts.
# The only use of generating normalized counts in this step is to decide which normalization
# is better for the next step.

# You need to add path for Python3 for MAGeCK to run, if MAGeCK was installed without Conda
# PATH=$HOME/Python3/bin/:$PATH

# You need to add path for Conda
# use "source activate <env_name>" instead of "conda activate <env_name>"
PATH=$HOME/miniconda3/bin/:$PATH
source activate R

# NOTE: Although paired end, in CRISPR screens, either read1 or read2 usually 
# contain the guide sequence. So, we do not need use both read1 and read2. 
# In Jinfen CRISPR, the reverse complement of the guides are present in read1 only.

printf "\n\n***********Counting directly from fastq files******************\n\n"

proj="CRISPR_Jinfen"
libraries=(DTKP)

#proj="CRISPR_Lena"
#libraries=(CDH12KO CDH12ACT ImmuneKO ImmuneACT)

#proj="CRISPR_Prince"
#libraries=(DTKPA1 DTKPA2 DTKPB1)

for library in ${libraries[*]}
do 

CORRECTED_GUIDE=$HOME/projects/CRISPR/$proj.$library.corrected.csv
CONTROL_GUIDE=$HOME/projects/CRISPR/$proj.$library.control_sgRNAs.txt		# RECOMMENDED to use ONLY safe sgRNAs as control
FASTQ_DIR=$HOME/scratch/$proj/cutadapt_results/$library
OUTPUT_DIR=$HOME/scratch/$proj/count_results/$library
mkdir $OUTPUT_DIR
cd $OUTPUT_DIR

methods=(mageck)
for method in ${methods[*]}
do

norms=(median) #(total median control)

for norm in ${norms[*]}
do	

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date    : $(date)"
echo "Job name      : $JOB_NAME"
echo "Job ID        : $JOB_ID"  
echo "SGE TASK ID   : $SGE_TASK_ID"
echo "Project       : $proj"
echo "Library		    : $library"
echo "Method        : $method"
echo "Normalization : $norm"
echo "Fastq Folder  : $FASTQ_DIR"
echo "Output Folder : $OUTPUT_DIR"
echo "=========================================================="

if [ $library == CDH12KO ]
then
mageck count \
--list-seq $CORRECTED_GUIDE \
--fastq $FASTQ_DIR/PCR-01_S1_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-01DO_S52_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-02_S2_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-03_S3_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-04_S4_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-05_S5_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-05DO_S56_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-06_S6_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-07_S7_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-08_S8_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-09_S9_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-09DO_S61_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-10_S10_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-11_S11_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-12_S12_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-13_S13_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-14_S14_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-15_S15_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-16_S16_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-17_S17_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-18_S18_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-19_S19_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-20_S20_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-21_S21_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-22_S22_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-23_S23_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-24_S24_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-25_S25_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-26_S26_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-27_S27_R1_trimmed.fastq.gz \
--sample-label PCR-01,PCR-01DO,PCR-02,PCR-03,PCR-04,PCR-05,PCR-05DO,PCR-06,PCR-07,PCR-08,\
PCR-09,PCR-09DO,PCR-10,PCR-11,PCR-12,PCR-13,PCR-14,PCR-15,PCR-16,PCR-17,PCR-18,PCR-19,\
PCR-20,PCR-21,PCR-22,PCR-23,PCR-24,PCR-25,PCR-26,PCR-27 \
--norm-method=$norm \
--output-prefix=$method.$norm \
--trim-5=0 \
--pdf-report

elif [ $library == CDH12ACT ]
then
mageck count \
--list-seq $CORRECTED_GUIDE \
--fastq $FASTQ_DIR/PCR-02DO_S53_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-06-2DO_S58_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-06DO_S57_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-10DO_S62_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-28_S1_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-29_S2_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-30_S3_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-31_S4_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-32_S5_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-33_S6_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-34_S7_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-35_S8_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-36_S9_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-37_S10_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-38_S11_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-39_S12_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-40_S13_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-41_S14_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-42_S15_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-43_S16_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-44_S17_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-45_S18_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-46_S19_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-47_S20_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-48_S21_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-49_S22_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-50_S23_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-51_S24_R1_trimmed.fastq.gz \
--sample-label PCR-02DO,PCR-06-2DO,PCR-06DO,PCR-10DO,PCR-28,PCR-29,PCR-30,PCR-31,PCR-32,\
PCR-33,PCR-34,PCR-35,PCR-36,PCR-37,PCR-38,PCR-39,PCR-40,PCR-41,PCR-42,PCR-43,PCR-44,\
PCR-45,PCR-46,PCR-47,PCR-48,PCR-49,PCR-50,PCR-51 \
--norm-method=$norm \
--output-prefix=$method.$norm \
--trim-5=0 \
--pdf-report

elif [ $library == ImmuneKO ]
then
mageck count \
--list-seq $CORRECTED_GUIDE \
--fastq $FASTQ_DIR/PCR-03DO_S54_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-07DO_S59_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-11DO_S63_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-52_S28_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-53_S29_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-54_S30_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-55_S31_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-56_S32_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-57_S33_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-58_S34_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-59_S35_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-60_S36_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-61_S37_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-62_S38_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-63_S39_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-64_S40_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-65_S41_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-66_S42_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-67_S43_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-68_S44_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-69_S45_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-70_S46_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-71_S47_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-72_S48_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-73_S49_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-74_S50_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-75_S51_R1_trimmed.fastq.gz \
--sample-label PCR-03DO,PCR-07DO,PCR-11DO,PCR-52,PCR-53,PCR-54,PCR-55,PCR-56,PCR-57,\
PCR-58,PCR-59,PCR-60,PCR-61,PCR-62,PCR-63,PCR-64,PCR-65,PCR-66,PCR-67,PCR-68,PCR-69,\
PCR-70,PCR-71,PCR-72,PCR-73,PCR-74,PCR-75 \
--norm-method=$norm \
--output-prefix=$method.$norm \
--trim-5=0 \
--pdf-report

elif [ $library == ImmuneACT ]
then
mageck count \
--list-seq $CORRECTED_GUIDE \
--fastq $FASTQ_DIR/PCR-04DO_S55_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-08DO_S60_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-12DO_S64_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-76_S25_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-77_S26_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-78_S27_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-79_S28_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-80_S29_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-81_S30_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-82_S31_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-83_S32_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-84_S33_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-85_S34_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-86_S35_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-87_S36_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-88_S37_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-89_S38_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-90_S39_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-91_S40_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-92_S41_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-93_S42_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-94_S43_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-95_S44_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-96_S45_R1_trimmed.fastq.gz \
$FASTQ_DIR/PCR-97_S46_R1_trimmed.fastq.gz \
--sample-label PCR-04DO,PCR-08DO,PCR-12DO,PCR-76,PCR-77,PCR-78,PCR-79,PCR-80,PCR-81,\
PCR-82,PCR-83,PCR-84,PCR-85,PCR-86,PCR-87,PCR-88,PCR-89,PCR-90,PCR-91,PCR-92,PCR-93,\
PCR-94,PCR-95,PCR-96,PCR-97 \
--norm-method=$norm \
--output-prefix=$method.$norm \
--trim-5=0 \
--pdf-report

elif [ $library == DTKPA1 ] || [ $library == DTKPA2 ]
then
mageck count \
--list-seq $CORRECTED_GUIDE \
--control-sgrna $CONTROL_GUIDE \
--fastq $FASTQ_DIR/AP10_R1.fastq.gz $FASTQ_DIR/AP12_R1.fastq.gz \
$FASTQ_DIR/AP2_R1.fastq.gz $FASTQ_DIR/AP4_R1.fastq.gz \
$FASTQ_DIR/AP6_R1.fastq.gz $FASTQ_DIR/AP8_R1.fastq.gz \
$FASTQ_DIR/AP11_R1.fastq.gz $FASTQ_DIR/AP1_R1.fastq.gz \
$FASTQ_DIR/AP3_R1.fastq.gz $FASTQ_DIR/AP5_R1.fastq.gz \
$FASTQ_DIR/AP7_R1.fastq.gz $FASTQ_DIR/AP9_R1.fastq.gz \
--sample-label AP10,AP12,AP2,AP4,AP6,AP8,AP11,AP1,\
AP3,AP5,AP7,AP9 \
--norm-method=$norm \
--output-prefix=$method.$norm \
--trim-5=0 \
--pdf-report

elif [ $library == DTKPB1 ]
then
mageck count \
--list-seq $CORRECTED_GUIDE \
--control-sgrna $CONTROL_GUIDE \
--fastq $FASTQ_DIR/AP10_R2.fastq.gz $FASTQ_DIR/AP12_R2.fastq.gz \
$FASTQ_DIR/AP2_R2.fastq.gz $FASTQ_DIR/AP4_R2.fastq.gz \
$FASTQ_DIR/AP6_R2.fastq.gz $FASTQ_DIR/AP8_R2.fastq.gz \
$FASTQ_DIR/AP11_R2.fastq.gz $FASTQ_DIR/AP1_R2.fastq.gz \
$FASTQ_DIR/AP3_R2.fastq.gz $FASTQ_DIR/AP5_R2.fastq.gz \
$FASTQ_DIR/AP7_R2.fastq.gz $FASTQ_DIR/AP9_R2.fastq.gz \
--sample-label AP10,AP12,AP2,AP4,AP6,AP8,AP11,AP1,\
AP3,AP5,AP7,AP9 \
--norm-method=$norm \
--output-prefix=$method.$norm \
--trim-5=0 \
--reverse-complement \
--pdf-report

elif [ $library == DTKP ]
then
mageck count \
--list-seq $CORRECTED_GUIDE \
--control-sgrna $CONTROL_GUIDE \
--fastq $FASTQ_DIR/Day0_2D_Rep1_1.fq.gz $FASTQ_DIR/Day0_2D_Rep2_1.fq.gz \
$FASTQ_DIR/Day0_3D_Rep1_1.fq.gz $FASTQ_DIR/Day0_3D_Rep2_1.fq.gz \
$FASTQ_DIR/Day14_2D_Rep1_1.fq.gz $FASTQ_DIR/Day14_2D_Rep2_1.fq.gz \
$FASTQ_DIR/Day14_3D_Rep1_1.fq.gz $FASTQ_DIR/Day14_3D_Rep2_1.fq.gz \
$FASTQ_DIR/Day7_2D_Rep1_1.fq.gz $FASTQ_DIR/Day7_2D_Rep2_1.fq.gz \
$FASTQ_DIR/Day7_3D_Rep1_1.fq.gz $FASTQ_DIR/Day7_3D_Rep2_1.fq.gz \
--sample-label Day0_2D_Rep1,Day0_2D_Rep2,Day0_3D_Rep1,Day0_3D_Rep2,\
Day14_2D_Rep1,Day14_2D_Rep2,Day14_3D_Rep1,Day14_3D_Rep2,\
Day7_2D_Rep1,Day7_2D_Rep2,Day7_3D_Rep1,Day7_3D_Rep2 \
--norm-method=$norm \
--output-prefix=$method.$norm \
--trim-5=0 \
--pdf-report \

else
echo $library
fi

done
done
done

#
#--day0-label 
#--reverse-complement 
#--test-run  
# Use --test-run to see if reverse complement of guide is present in reads. If so, 
# use the --reverse-complement and repeat --test-run to confirm and then finally run
# the script  

# mageck count \
# --list-seq $CORRECTED_GUIDE \
# --fastq $FASTQ_DIR/PCR-01_S1_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-01DO_S52_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-02_S2_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-02DO_S53_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-03_S3_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-03DO_S54_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-04_S4_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-04DO_S55_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-05_S5_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-05DO_S56_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-06_S6_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-06-2DO_S58_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-06DO_S57_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-07_S7_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-07DO_S59_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-08_S8_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-08DO_S60_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-09_S9_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-09DO_S61_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-10_S10_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-10DO_S62_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-11_S11_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-11DO_S63_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-12_S12_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-12DO_S64_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-13_S13_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-14_S14_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-15_S15_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-16_S16_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-17_S17_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-18_S18_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-19_S19_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-20_S20_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-21_S21_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-22_S22_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-23_S23_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-24_S24_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-25_S25_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-26_S26_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-27_S27_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-28_S1_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-29_S2_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-30_S3_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-31_S4_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-32_S5_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-33_S6_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-34_S7_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-35_S8_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-36_S9_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-37_S10_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-38_S11_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-39_S12_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-40_S13_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-41_S14_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-42_S15_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-43_S16_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-44_S17_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-45_S18_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-46_S19_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-47_S20_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-48_S21_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-49_S22_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-50_S23_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-51_S24_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-52_S28_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-53_S29_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-54_S30_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-55_S31_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-56_S32_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-57_S33_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-58_S34_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-59_S35_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-60_S36_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-61_S37_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-62_S38_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-63_S39_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-64_S40_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-65_S41_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-66_S42_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-67_S43_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-68_S44_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-69_S45_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-70_S46_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-71_S47_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-72_S48_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-73_S49_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-74_S50_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-75_S51_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-76_S25_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-77_S26_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-78_S27_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-79_S28_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-80_S29_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-81_S30_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-82_S31_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-83_S32_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-84_S33_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-85_S34_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-86_S35_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-87_S36_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-88_S37_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-89_S38_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-90_S39_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-91_S40_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-92_S41_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-93_S42_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-94_S43_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-95_S44_R1_trimmed.fastq.gz \
# $FASTQ_DIR/PCR-96_S45_R1_trimmed.fastq.gz $FASTQ_DIR/PCR-97_S46_R1_trimmed.fastq.gz \
# --sample-label PCR-01,PCR-01DO,PCR-02,PCR-02DO,PCR-03,PCR-03DO,PCR-04,PCR-04DO,PCR-05,\
# PCR-05DO,PCR-06,PCR-06-2DO,PCR-06DO,PCR-07,PCR-07DO,PCR-08,PCR-08DO,PCR-09,PCR-09DO,\
# PCR-10,PCR-10DO,PCR-11,PCR-11DO,PCR-12,PCR-12DO,PCR-13,PCR-14,PCR-15,PCR-16,PCR-17,\
# PCR-18,PCR-19,PCR-20,PCR-21,PCR-22,PCR-23,PCR-24,PCR-25,PCR-26,PCR-27,PCR-28,PCR-29,\
# PCR-30,PCR-31,PCR-32,PCR-33,PCR-34,PCR-35,PCR-36,PCR-37,PCR-38,PCR-39,PCR-40,PCR-41,\
# PCR-42,PCR-43,PCR-44,PCR-45,PCR-46,PCR-47,PCR-48,PCR-49,PCR-50,PCR-51,PCR-52,PCR-53,\
# PCR-54,PCR-55,PCR-56,PCR-57,PCR-58,PCR-59,PCR-60,PCR-61,PCR-62,PCR-63,PCR-64,PCR-65,\
# PCR-66,PCR-67,PCR-68,PCR-69,PCR-70,PCR-71,PCR-72,PCR-73,PCR-74,PCR-75,PCR-76,PCR-77,\
# PCR-78,PCR-79,PCR-80,PCR-81,PCR-82,PCR-83,PCR-84,PCR-85,PCR-86,PCR-87,PCR-88,PCR-89,\
# PCR-90,PCR-91,PCR-92,PCR-93,PCR-94,PCR-95,PCR-96,PCR-97 \
# --norm-method=$norm \
# --output-prefix=$method.$norm \
# --trim-5=0 \
# --pdf-report




# mageck count \
# --list-seq $CORRECTED_GUIDE \
# --fastq $FASTQ_DIR/Day0_2D_Rep1_1.fq.gz $FASTQ_DIR/Day0_2D_Rep2_1.fq.gz \
# $FASTQ_DIR/Day7_2D_Rep1_1.fq.gz $FASTQ_DIR/Day7_2D_Rep2_1.fq.gz \
# $FASTQ_DIR/Day14_2D_Rep1_1.fq.gz $FASTQ_DIR/Day14_2D_Rep2_1.fq.gz \
# $FASTQ_DIR/Day0_3D_Rep1_1.fq.gz $FASTQ_DIR/Day0_3D_Rep2_1.fq.gz \
# $FASTQ_DIR/Day7_3D_Rep1_1.fq.gz $FASTQ_DIR/Day7_3D_Rep2_1.fq.gz \
# $FASTQ_DIR/Day14_3D_Rep1_1.fq.gz $FASTQ_DIR/Day14_3D_Rep2_1.fq.gz \
# --norm-method=control \
# --control-sgrna $CONTROL_GUIDE \
# --sample-label Day0_2D_Rep1,Day0_2D_Rep2,Day7_2D_Rep1,Day7_2D_Rep2,\
# Day14_2D_Rep1,Day14_2D_Rep2,Day0_3D_Rep1,Day0_3D_Rep2,Day7_3D_Rep1,\
# Day7_3D_Rep2,Day14_3D_Rep1,Day14_3D_Rep2 \
# --output-prefix=$PREFIX \
# --trim-5=AUTO \
# --day0-label Day0_2D_Rep1,Day0_2D_Rep2,Day0_3D_Rep1,Day0_3D_Rep2 \
# --pdf-report


# methods=(bowtie2 bowtie)
# for method in ${methods[*]}
# do

# printf "%s " ***********Counting from $method aligned bam files******************
# printf "\n"

# FASTQ_DIR=$HOME/scratch/$proj/alignment_results/$method
# PREFIX=$method
# cd $OUTPUT_DIR/$method

# mageck count \
# --list-seq $CORRECTED_GUIDE \
# --fastq $FASTQ_DIR/Day0_2D_Rep1_1.$method.bam $FASTQ_DIR/Day0_2D_Rep2_1.$method.bam \
# $FASTQ_DIR/Day7_2D_Rep1_1.$method.bam $FASTQ_DIR/Day7_2D_Rep2_1.$method.bam \
# $FASTQ_DIR/Day14_2D_Rep1_1.$method.bam $FASTQ_DIR/Day14_2D_Rep2_1.$method.bam \
# $FASTQ_DIR/Day0_3D_Rep1_1.$method.bam $FASTQ_DIR/Day0_3D_Rep2_1.$method.bam \
# $FASTQ_DIR/Day7_3D_Rep1_1.$method.bam $FASTQ_DIR/Day7_3D_Rep2_1.$method.bam \
# $FASTQ_DIR/Day14_3D_Rep1_1.$method.bam $FASTQ_DIR/Day14_3D_Rep2_1.$method.bam \
# --norm-method=control \
# --control-sgrna $CONTROL_GUIDE \
# --sample-label Day0_2D_Rep1,Day0_2D_Rep2,Day7_2D_Rep1,Day7_2D_Rep2,\
# Day14_2D_Rep1,Day14_2D_Rep2,Day0_3D_Rep1,Day0_3D_Rep2,Day7_3D_Rep1,\
# Day7_3D_Rep2,Day14_3D_Rep1,Day14_3D_Rep2 \
# --output-prefix=$PREFIX \
# --trim-5=AUTO \
# --day0-label Day0_2D_Rep1,Day0_2D_Rep2,Day0_3D_Rep1,Day0_3D_Rep2 \
# --pdf-report

# done

# Read about MAGeCK here 
# https://sourceforge.net/p/mageck/wiki/Home/
# https://sourceforge.net/p/mageck/wiki/output/
# mageck count --help
# mageck test --help

# countsummary.txt has following columns. Recommended values are shown in brackets.
# After mageck count is complete, open countsummary.txt in excel.
# Column		Content
# File			The fastq (or the count table) file used.
# Label			The label of that fastq file assigned.
# Reads			Total number reads in the fastq file. 
# Mapped		Total number of reads mapped to library (Recommended: 300 times the number of sgRNAs)
# Percentage	Mapped percentage, calculated as Mapped/Reads (Recommended: at least 60%)
# TotalsgRNAs	Total number of sgRNAs in the library
# Zerocounts	Total number of sgRNAs that have 0 counts (Recommended: no more than 1%)
# GiniIndex		The Gini Index of the read count distribution. A smaller value indicates more eveness of the count distribution. 
# 				(Recommended: around 0.1 for plasmid/initial state samples, around 0.2-0.3 for negative selection samples)

# --fastq : use space to separate biological replicates and , to separate technical replicates

# --norm-method=control is RECOMMENDED to create accurate null distribution.
# Using --control-sgrna or --control-gene is RECOMMENDED as it improves FDR. 
# By providing the corresponding sgRNA IDs in --control-sgrna or sgRNA gene in --control-gene, 
# MAGeCK will have a better estimation of p values.
# Use safe sgRNAs instead of non-targeting sgRNAs as control (https://doi.org/10.1038/ncomms15178)

# --day0-label doesnt affect count results. It ONLY helps in calculating the
# following columns in the countsummary.txt file.
# NegSelQC						The Enrichment Score (ES) of GSEA
# NegSelQCPval					The p value of the GSEA analysis (Recommended: smaller than 1e-10)
# NegSelQCPvalPermutation		The permutation p value
# NegSelQCPvalPermutationFDR	The FDR of the permutation p value
# NegSelQCGene					The number of essential genes found in the library that are evaluated for GSEA analysis. 