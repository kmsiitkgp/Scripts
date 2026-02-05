#!/bin/bash -l
# Run this ONLY ONCE each for human and mouse

#$ -N Germline_Ref              # Set Job Name
#$ -l mem_free=96G              # Request memory
#$ -j y                         # Merge standard output and standard error
#$ -cwd                         # Set current working directory
#$ -t 1-1                       # Submit an array job with 1 task

# In order to accurately eliminate germline mutations, we need
# (i) a germline reference (OR)
# (ii) a matched normal sample for every tumor sample

# Human germline reference are available from gnomad (?? figure out how to download)
# Mouse germline reference are avaiable from "The Mouse Genome Project" https://www.mousegenomes.org/snps-indels/

# Mouse:
# Download the latest version of SNPs and indels with rsIDs
# rsID stands for "reference SNP cluster ID". They are provided by NCBI and 
# allows for easy identification and comparison of genetic variants across studies.
# SNP file: mgp_REL2021_snps.rsID.vcf.gz
# INDEL file: mgp_REL2021_indels.rsID.vcf.gz

# We need to combine the SNPs and INDELs.
# Make sure the SNP and INDEL vcf files have the same column names in same ordering using zcat and head.
# Both SNP and INDEL vcf files need index files. Create them before concatenating.
SNP_vcf=$HOME/NGSTools/WES/mgp_REL2021_snps.rsID.vcf.gz
INDEL_vcf=$HOME/NGSTools/WES/mgp_REL2021_indels.rsID.vcf.gz
SNP_INDEL_vcf=$HOME/NGSTools/WES/mgp_REL2021_snps.indels.rsID.vcf.gz
C57BL6_SNP_INDEL_vcf=$HOME/NGSTools/WES/mgp_REL2021_C57BL6_snps.indels.rsID.vcf.gz
C57BL6_GERMLINE_vcf=$HOME/NGSTools/WES/mgp_REL2021_C57BL6_snps.indels.rsID.AF.vcf.gz
C57BL6_GERMLINE_REF_AF_ONLY_vcf=$HOME/NGSTools/WES/mgp_REL2021_C57BL6_snps.indels.rsID.AFonly.vcf.gz

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date    		: $(date)"
echo "Job name      		: $JOB_NAME"
echo "Job ID        		: $JOB_ID"  
echo "SGE TASK ID   		: $SGE_TASK_ID"
echo "SNP VCF			  	: $SNP_vcf"
echo "INDEL VCF			  	: $INDEL_vcf"
echo "Germline Reference  	: $C57BL6_GERMLINE_vcf"
echo "=========================================================="

# You need to load bcftools module from the cluster
module load bcftools/1.20

cd $HOME/NGSTools/WES/
bcftools index $SNP_vcf
bcftools index $INDEL_vcf
bcftools concat --allow-overlaps --output $SNP_INDEL_vcf --output-type z $SNP_vcf $INDEL_vcf  

# The vcf file has variant info for several mouse strains.
# We filter out variant info for unwanted strains and retain only C57BL6 which we want
bcftools view --samples C57BL_6NJ --output $C57BL6_SNP_INDEL_vcf --output-type z $SNP_INDEL_vcf 

# The germline reference for GATK needs population allele frequencies (AF) in the INFO field
# for every variant. The header of the vcf files do not mention AF tag. So, we create them.
# Also, we need to create an index file for the germline reference. 
# The file extension of index created using GATK (.tbi) and bcftools (.csi) are different.
# GATK needs inedx files to be in  .tbi format. So, bcftools index wont work with GATK
bcftools plugin fill-tags --output $C57BL6_GERMLINE_vcf --output-type z $C57BL6_SNP_INDEL_vcf -- --tags AF

bcftools annotate \
--remove FORMAT,C57BL_6NJ,^INFO/AF \
--output $C57BL6_GERMLINE_REF_AF_ONLY_vcf \
--output-type z \
--write-index=tbi \
$C57BL6_GERMLINE_vcf

# bcftools index $C57BL6_GERMLINE_vcf
~/NGSTools/gatk-4.6.1.0/gatk IndexFeatureFile \
--input $C57BL6_GERMLINE_vcf