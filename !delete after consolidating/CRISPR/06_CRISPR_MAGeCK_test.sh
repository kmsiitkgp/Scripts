#!/bin/bash -l

#$ -N CRISPR_DEGs               # Set Job Name
#$ -l mem_free=32G              # Request memory
#$ -j y                         # Merge standard output and standard error
#$ -cwd                         # Set current working directory
#$ -t 1-1                       # Submit an array job with 1 task.

# https://doi.org/10.1186/s13059-020-01972-x
# This article recommends MAGeCK RRA over MAGeCK MLE

# IMPORTANT: Sometimes the count table has extra " in sgRNA id and Gene columns. This is because your guides.csv had this problem.
# Read 03_CRISPR_bwamem2_bowtie2_bowtie_index.sh to addresss this issue properly. 
# So, open the count table and replace all " with nothing. Save and re-open in excel to make sure count table is proper.
# sgRNA	Gene	Ctrl_A_1	Ctrl_A_2	Day07_2D_1	Day07_2D_2	Day14_2D_1	Day14_2D_2	Day07_3D_1	Day07_3D_2	Ctrl_B_1	Ctrl_B_2	Day14_3D_1	Day14_3D_2
#"2663	CACNB3"	2142	2959	2542	3192	3759	3277	1339	1152	1340	2252	850	1236
#"4595	CPT2"	210	339	292	300	368	284	85	44	214	281	131	123

# You need to add path for Python3 for MAGeCK to run, if MAGeCK was installed without Conda
# PATH=$HOME/Python3/bin/:$PATH

# You need to add path for Conda
# use "source activate <env_name>" instead of "conda activate <env_name>"
PATH=$HOME/miniconda3/bin/:$PATH
source activate R

proj="CRISPR_Jinfen"
libraries=(DTKP)

#proj="CRISPR_Lena"
#libraries=(CDH12KO CDH12ACT ImmuneKO ImmuneACT)

#proj="CRISPR_Prince"
#libraries=(DTKPA1 DTKPA2 DTKPB1)

COUNT_DIR=$HOME/scratch/$proj/count_results
OUTPUT_DIR=$HOME/scratch/$proj/DEG_results

for library in ${libraries[*]}
do 

mkdir $OUTPUT_DIR/$library
CONTROL_GUIDE=$HOME/projects/CRISPR/$proj.$library.control_sgRNAs.txt		# RECOMMENDED to use ONLY safe sgRNAs as control

methods=(mageck)   #(bowtie2 bowtie mageck kms)
norms=(median)     #(total median control)
lfcs=(alphamedian) #(median alphamedian secondbest)

for method in ${methods[*]}
do

for norm in ${norms[*]}
do	

for lfc in ${lfcs[*]}
do

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date    : $(date)"
echo "Job name      : $JOB_NAME"
echo "Job ID        : $JOB_ID"  
echo "SGE TASK ID   : $SGE_TASK_ID"
echo "Project       : $proj"
echo "Library		: $library"
echo "Method        : $method"
echo "Normalization : $norm"
echo "LFC Method    : $lfc"
echo "Output Folder : $OUTPUT_DIR"
echo "=========================================================="

printf " %s %s" ***********Identifying Differential Guides using "$method" and "$norm" ******************
printf "\n"

if [ $library == CDH12KO ]
then

cd $OUTPUT_DIR/$library

mageck test \
--count-table $COUNT_DIR/$library/$method.$norm.count.txt \
--treatment-id PCR-01,PCR-02,PCR-03,PCR-04,\
PCR-05,PCR-06,PCR-07,PCR-08,PCR-09,PCR-10 \
--control-id PCR-01DO \
--norm-method=$norm \
--gene-test-fdr-threshold=0.25 \
--adjust-method=fdr \
--sort-criteria=neg \
--remove-zero=both \
--remove-zero-threshold=0 \
--gene-lfc-method=$lfc \
--output-prefix=MB49c5.$library.$method.$norm.$lfc \
--pdf-report
#--day0-label \
#--cnv-norm <matrix>\
#--cell-line=NA13 \
#--cnv-est <bed file>
#--control-sgrna $CONTROL_GUIDE \

mageck test \
--count-table $COUNT_DIR/$library/$method.count.txt \
--treatment-id PCR-11,PCR-12,PCR-13,PCR-14,\
PCR-15,PCR-16,PCR-17 \
--control-id PCR-05DO \
--norm-method=$norm \
--gene-test-fdr-threshold=0.25 \
--adjust-method=fdr \
--sort-criteria=neg \
--remove-zero=both \
--remove-zero-threshold=0 \
--gene-lfc-method=$lfc \
--output-prefix=MB83F.$library.$method.$norm.$lfc \
--pdf-report

mageck test \
--count-table $COUNT_DIR/$library/$method.count.txt \
--treatment-id PCR-18,PCR-19,PCR-20,PCR-21,PCR-22,\
PCR-23,PCR-24,PCR-25,PCR-26,PCR-27 \
--control-id PCR-09DO \
--norm-method=$norm \
--gene-test-fdr-threshold=0.25 \
--adjust-method=fdr \
--sort-criteria=neg \
--remove-zero=both \
--remove-zero-threshold=0 \
--gene-lfc-method=$lfc \
--output-prefix=MB81F.$library.$method.$norm.$lfc \
--pdf-report


elif [ $library == CDH12ACT ]
then

cd $OUTPUT_DIR/$library

mageck test \
--count-table $COUNT_DIR/$library/$method.count.txt \
--treatment-id PCR-28,PCR-29,PCR-30,PCR-31,\
PCR-32,PCR-33,PCR-34 \
--control-id PCR-10DO \
--norm-method=$norm \
--gene-test-fdr-threshold=0.25 \
--adjust-method=fdr \
--sort-criteria=neg \
--remove-zero=both \
--remove-zero-threshold=0 \
--gene-lfc-method=$lfc \
--output-prefix=MB81F.$library.$method.$norm.$lfc \
--pdf-report

mageck test \
--count-table $COUNT_DIR/$library/$method.count.txt \
--treatment-id PCR-35,PCR-36,PCR-37,PCR-38,\
PCR-39,PCR-40,PCR-41,PCR-42,PCR-43,PCR-44 \
--control-id PCR-02DO \
--norm-method=$norm \
--gene-test-fdr-threshold=0.25 \
--adjust-method=fdr \
--sort-criteria=neg \
--remove-zero=both \
--remove-zero-threshold=0 \
--gene-lfc-method=$lfc \
--output-prefix=MB49c5.$library.$method.$norm.$lfc \
--pdf-report

mageck test \
--count-table $COUNT_DIR/$library/$method.count.txt \
--treatment-id PCR-46,PCR-47 \
--control-id PCR-06DO \
--norm-method=$norm \
--gene-test-fdr-threshold=0.25 \
--adjust-method=fdr \
--sort-criteria=neg \
--remove-zero=both \
--remove-zero-threshold=0 \
--gene-lfc-method=$lfc \
--output-prefix=MB83F.$library._1.$method.$norm.$lfc \
--pdf-report

mageck test \
--count-table $COUNT_DIR/$library/$method.count.txt \
--treatment-id PCR-48,PCR-50,PCR-51 \
--control-id PCR-06-2DO \
--norm-method=$norm \
--gene-test-fdr-threshold=0.25 \
--adjust-method=fdr \
--sort-criteria=neg \
--remove-zero=both \
--remove-zero-threshold=0 \
--gene-lfc-method=$lfc \
--output-prefix=MB83F.$library._2.$method.$norm.$lfc \
--pdf-report

elif [ $library == ImmuneKO ]
then

cd $OUTPUT_DIR/$library

mageck test \
--count-table $COUNT_DIR/$library/$method.count.txt \
--treatment-id PCR-52,PCR-53,PCR-54 \
--control-id PCR-11DO \
--norm-method=$norm \
--gene-test-fdr-threshold=0.25 \
--adjust-method=fdr \
--sort-criteria=neg \
--remove-zero=both \
--remove-zero-threshold=0 \
--gene-lfc-method=$lfc \
--output-prefix=MB81F.$library.$method.$norm.$lfc \
--pdf-report

mageck test \
--count-table $COUNT_DIR/$library/$method.count.txt \
--treatment-id PCR-58,PCR-59,PCR-61,PCR-62,\
PCR-63,PCR-64,PCR-65,PCR-66 \
--control-id PCR-03DO \
--norm-method=$norm \
--gene-test-fdr-threshold=0.25 \
--adjust-method=fdr \
--sort-criteria=neg \
--remove-zero=both \
--remove-zero-threshold=0 \
--gene-lfc-method=$lfc \
--output-prefix=MB49c5.$library.$method.$norm.$lfc \
--pdf-report

mageck test \
--count-table $COUNT_DIR/$library/$method.count.txt \
--treatment-id PCR-67,PCR-69,PCR-70,\
PCR-71,PCR-72,PCR-73,PCR-74,PCR-75 \
--control-id PCR-07DO \
--norm-method=$norm \
--gene-test-fdr-threshold=0.25 \
--adjust-method=fdr \
--sort-criteria=neg \
--remove-zero=both \
--remove-zero-threshold=0 \
--gene-lfc-method=$lfc \
--output-prefix=MB83F.$library.$method.$norm.$lfc \
--pdf-report

elif [ $library == ImmuneACT ]
then

cd $OUTPUT_DIR/$library

mageck test \
--count-table $COUNT_DIR/$library/$method.count.txt \
--treatment-id PCR-76,PCR-77,PCR-78,PCR-79,PCR-80 \
--control-id PCR-12DO \
--norm-method=$norm \
--gene-test-fdr-threshold=0.25 \
--adjust-method=fdr \
--sort-criteria=neg \
--remove-zero=both \
--remove-zero-threshold=0 \
--gene-lfc-method=$lfc \
--output-prefix=MB81F.$library.$method.$norm.$lfc \
--pdf-report

mageck test \
--count-table $COUNT_DIR/$library/$method.count.txt \
--treatment-id PCR-81,PCR-82,PCR-83,PCR-84,\
PCR-85,PCR-86,PCR-87,PCR-88,PCR-89,PCR-90 \
--control-id PCR-04DO \
--norm-method=$norm \
--gene-test-fdr-threshold=0.25 \
--adjust-method=fdr \
--sort-criteria=neg \
--remove-zero=both \
--remove-zero-threshold=0 \
--gene-lfc-method=$lfc \
--output-prefix=MB49c5.$library.$method.$norm.$lfc \
--pdf-report

mageck test \
--count-table $COUNT_DIR/$library/$method.count.txt \
--treatment-id PCR-91,PCR-93,PCR-94,PCR-95,PCR-97 \
--control-id PCR-08DO \
--norm-method=$norm \
--gene-test-fdr-threshold=0.25 \
--adjust-method=fdr \
--sort-criteria=neg \
--remove-zero=both \
--remove-zero-threshold=0 \
--gene-lfc-method=$lfc \
--output-prefix=MB83F.$library.$method.$norm.$lfc \
--pdf-report

elif [ $library == DTKPA1 ]  || [ $library == DTKPA2 ] || [ $library == DTKPB1 ]
then

cd $OUTPUT_DIR/$library

mageck test \
--count-table $COUNT_DIR/$library/$method.$norm.count.txt \
--treatment-id AP7,AP8,AP9 \
--control-id AP1,AP2,AP3 \
--norm-method=$norm \
--gene-test-fdr-threshold=0.25 \
--adjust-method=fdr \
--sort-criteria=neg \
--remove-zero=both \
--remove-zero-threshold=0 \
--gene-lfc-method=$lfc \
--control-sgrna $CONTROL_GUIDE \
--output-prefix=T14vsT0_YKO.$library.$method.$norm.$lfc \
--pdf-report

mageck test \
--count-table $COUNT_DIR/$library/$method.$norm.count.txt \
--treatment-id AP10,AP11,AP12 \
--control-id AP4,AP5,AP6 \
--norm-method=$norm \
--gene-test-fdr-threshold=0.25 \
--adjust-method=fdr \
--sort-criteria=neg \
--remove-zero=both \
--remove-zero-threshold=0 \
--gene-lfc-method=$lfc \
--control-sgrna $CONTROL_GUIDE \
--output-prefix=T14vsT0_SCR.$library.$method.$norm.$lfc \
--pdf-report

elif [ $library == DTKP ]
then

cd $OUTPUT_DIR/$library

mageck test \
--count-table $COUNT_DIR/$library/$method.$norm.count.txt \
--treatment-id Day7_2D_Rep1,Day7_2D_Rep2 \
--control-id Day0_2D_Rep1,Day0_2D_Rep2 \
--norm-method=$norm \
--gene-test-fdr-threshold=0.25 \
--adjust-method=fdr \
--sort-criteria=neg \
--remove-zero=both \
--remove-zero-threshold=0 \
--gene-lfc-method=$lfc \
--control-sgrna $CONTROL_GUIDE \
--output-prefix=D7vsD0_2D.$library.$method.$norm.$lfc \
--pdf-report

mageck test \
--count-table $COUNT_DIR/$library/$method.$norm.count.txt \
--treatment-id Day14_2D_Rep1,Day14_2D_Rep2 \
--control-id Day0_2D_Rep1,Day0_2D_Rep2 \
--norm-method=$norm \
--gene-test-fdr-threshold=0.25 \
--adjust-method=fdr \
--sort-criteria=neg \
--remove-zero=both \
--remove-zero-threshold=0 \
--gene-lfc-method=$lfc \
--control-sgrna $CONTROL_GUIDE \
--output-prefix=D14vsD0_2D.$library.$method.$norm.$lfc \
--pdf-report

mageck test \
--count-table $COUNT_DIR/$library/$method.$norm.count.txt \
--treatment-id Day7_3D_Rep1,Day7_3D_Rep2 \
--control-id Day0_3D_Rep1,Day0_3D_Rep2 \
--norm-method=$norm \
--gene-test-fdr-threshold=0.25 \
--adjust-method=fdr \
--sort-criteria=neg \
--remove-zero=both \
--remove-zero-threshold=0 \
--gene-lfc-method=$lfc \
--control-sgrna $CONTROL_GUIDE \
--output-prefix=D7vsD0_3D.$library.$method.$norm.$lfc \
--pdf-report

mageck test \
--count-table $COUNT_DIR/$library/$method.$norm.count.txt \
--treatment-id Day14_3D_Rep1,Day14_3D_Rep2 \
--control-id Day0_3D_Rep1,Day0_3D_Rep2 \
--norm-method=$norm \
--gene-test-fdr-threshold=0.25 \
--adjust-method=fdr \
--sort-criteria=neg \
--remove-zero=both \
--remove-zero-threshold=0 \
--gene-lfc-method=$lfc \
--control-sgrna $CONTROL_GUIDE \
--output-prefix=D14vsD0_3D.$library.$method.$norm.$lfc \
--pdf-report

else
echo $library
fi

done
done
done
done

## CRISPR Jinfen comparisons

# mageck test \
# --count-table $COUNT_DIR/$method/$method.count.txt \
# --treatment-id Day7_2D_Rep1,Day7_2D_Rep2 \
# --control-id Day0_2D_Rep1,Day0_2D_Rep2 \
# --norm-method=control \
# --control-sgrna $CONTROL_GUIDE \
# --gene-test-fdr-threshold=0.25 \
# --adjust-method=fdr \
# --sort-criteria=neg \
# --remove-zero=both \
# --remove-zero-threshold=0 \
# --gene-lfc-method=median \
# --output-prefix=Day7vs0_2D \
# --pdf-report
# #--day0-label \
# #--cnv-norm <matrix>\
# #--cell-line=NA13 \
# #--cnv-est <bed file>

# mageck test \
# --count-table $COUNT_DIR/$method/$method.count.txt \
# --treatment-id Day14_2D_Rep1,Day14_2D_Rep2 \
# --control-id Day0_2D_Rep1,Day0_2D_Rep2 \
# --norm-method=control \
# --control-sgrna $CONTROL_GUIDE \
# --gene-test-fdr-threshold=0.25 \
# --adjust-method=fdr \
# --sort-criteria=neg \
# --remove-zero=both \
# --remove-zero-threshold=0 \
# --gene-lfc-method=median \
# --output-prefix=Day14vs0_2D \
# --pdf-report

# # From PCA plot, Day0_3D_Rep2 and Day14_3D_Rep2 doesnt look right. So, we remove them from analysis.
# mageck test \
# --count-table $COUNT_DIR/$method/$method.count.txt \
# --treatment-id Day7_3D_Rep1,Day7_3D_Rep2 \
# --control-id Day0_3D_Rep1 \
# --norm-method=control \
# --control-sgrna $CONTROL_GUIDE \
# --gene-test-fdr-threshold=0.25 \
# --adjust-method=fdr \
# --sort-criteria=neg \
# --remove-zero=both \
# --remove-zero-threshold=0 \
# --gene-lfc-method=median \
# --output-prefix=Day7vs0_3D \
# --pdf-report

# # From PCA plot, Day0_3D_Rep2 and Day14_3D_Rep2 doesnt look right. So, we remove them from analysis.
# mageck test \
# --count-table $COUNT_DIR/$method/$method.count.txt \
# --treatment-id Day14_3D_Rep1 \
# --control-id Day0_3D_Rep1 \
# --norm-method=control \
# --control-sgrna $CONTROL_GUIDE \
# --gene-test-fdr-threshold=0.25 \
# --adjust-method=fdr \
# --sort-criteria=neg \
# --remove-zero=both \
# --remove-zero-threshold=0 \
# --gene-lfc-method=median \
# --output-prefix=Day14vs0_3D \
# --pdf-report

# # From PCA plot, Day0_3D_Rep2 and Day14_3D_Rep2 doesnt look right. So, we remove them from analysis.
# mageck test \
# --count-table $COUNT_DIR/$method/$method.count.txt \
# --treatment-id Day0_3D_Rep1 \
# --control-id Day0_2D_Rep1,Day0_2D_Rep2 \
# --norm-method=control \
# --control-sgrna $CONTROL_GUIDE \
# --gene-test-fdr-threshold=0.25 \
# --adjust-method=fdr \
# --sort-criteria=neg \
# --remove-zero=both \
# --remove-zero-threshold=0 \
# --gene-lfc-method=median \
# --output-prefix=Day0_3Dvs2D \
# --pdf-report

# mageck test \
# --count-table $COUNT_DIR/$method/$method.count.txt \
# --treatment-id Day7_3D_Rep1,Day7_3D_Rep2 \
# --control-id Day7_2D_Rep1,Day7_2D_Rep2 \
# --norm-method=control \
# --control-sgrna $CONTROL_GUIDE \
# --gene-test-fdr-threshold=0.25 \
# --adjust-method=fdr \
# --sort-criteria=neg \
# --remove-zero=both \
# --remove-zero-threshold=0 \
# --gene-lfc-method=median \
# --output-prefix=Day7_3Dvs2D \
# --pdf-report

# # From PCA plot, Day0_3D_Rep2 and Day14_3D_Rep2 doesnt look right. So, we remove them from analysis.
# mageck test \
# --count-table $COUNT_DIR/$method/$method.count.txt \
# --treatment-id Day14_3D_Rep1 \
# --control-id Day14_2D_Rep1,Day14_2D_Rep2 \
# --norm-method=control \
# --control-sgrna $CONTROL_GUIDE \
# --gene-test-fdr-threshold=0.25 \
# --adjust-method=fdr \
# --sort-criteria=neg \
# --remove-zero=both \
# --remove-zero-threshold=0 \
# --gene-lfc-method=median \
# --output-prefix=Day14_3Dvs2D \
# --pdf-report

##*************************************

# Read about MAGeCK here 
# https://sourceforge.net/p/mageck/wiki/Home/
# https://sourceforge.net/p/mageck/wiki/output/
# mageck count --help
# mageck test --help

# --count-table: 	output of mageck count
# --treatment-id:	Labels we assigned to fastq file while running mageck count
# --control-id:		Labels we assigned to fastq file while running mageck count
# --norm-method: 	Use "total" or "control". If using "control", use SAFE sgRNAs as "control", not Non-targeting sgRNAs
# --output-prefix:	The prefix to be added o output of mageck test
# --day0-label:		Use this ONLY if you dont specify --control-id. DO NOT use both parameters simultaneously.
# 					If you want to compare more than 2 sets of treated samples (say Day 7 and Day 14) to the same set of control samples (Day 0), then 
# 					use this parameter. If you want to compare only 2 sets of samples: "Day 7 vs Day 0"/"Day 14 vs Day 0"/"Day 7 vs Day 14", use --control id. 
# 					Turns on negative selection QC. Usually day 0 or plasmid sample is given this label. For every other sample label,
#					the negative selection QC will compare it with day0 sample, and estimate the degree of negative selections in essential genes
# --paired: 		Only use if each treatment sample has ONLY one corresponding control sample 
# --control-sgrna: 	Only use if --norm-method is control. Use SAFE sgRNAs
# --control-gene: 	Only use if --norm-method is control. Use genes corresponding to SAFE sgRNAs
# --cnv-norm: 		A matrix of copy number variation data across cell lines to normalize CNV-biased sgRNA scores prior to gene ranking
# --cell-line:		Cell line to be used for cnv normalization. Must match one of the column names in the file provided by --cnv-norm