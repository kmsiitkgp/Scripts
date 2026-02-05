#!/bin/bash -l

#$ -N mutect2                   # Set Job Name
#$ -l mem_free=96G              # Request memory
#$ -j y                         # Merge standard output and standard error
#$ -cwd                         # Set current working directory
#$ -t 1-2                       # Submit an array job with n tasks; n=number of tumor samples

species="Mouse"
proj="WES_Hany"
INPUT_DIR=$HOME/scratch/$proj/alignment_results  
OUTPUT_DIR=$HOME/scratch/$proj/variant_results
REF_FILE=$HOME/NGSTools/Reference_Genomes/$species/*.fa
PoN_VCF=$OUTPUT_DIR/PoN.vcf.gz
INTERVAL_BED=$HOME/NGSTools/WES/S0276129_SureSelect_Mouse_AllExon_V1_mm39.bed
GERMLINE_REF=$HOME/NGSTools/WES/mgp_REL2021_C57BL6_snps.indels.rsID.AF.vcf.gz
GERMLINE_REF_AF_ONLY=$HOME/NGSTools/WES/mgp_REL2021_C57BL6_snps.indels.rsID.AFonly.vcf.gz

# Create an array of input files
NORMAL=(MB49.coord.sorted.bam)
TUMOR=(MBN.coord.sorted.bam MBP.coord.sorted.bam)
NORMAL=(${NORMAL[*]/#/$INPUT_DIR/})
TUMOR=(${TUMOR[*]/#/$INPUT_DIR/})

# Use the SGE_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value
index=$(($SGE_TASK_ID-1))
tumor_bam=${TUMOR[$index]}   								# tumor sample with path
normal_bam=${NORMAL[0]}										# normal sample with path
normal_sample=${normal_bam//.coord.sorted.bam/}     		# This MUST match the SM field in RG you indicated in bwamem during alignment
normal_sample=$(basename $normal_sample)					

# Create output filenames
unfiltered_vcf=${tumor_bam//.coord.sorted.bam/_unfiltered.vcf.gz}
unfiltered_vcf=${unfiltered_vcf//$INPUT_DIR/$OUTPUT_DIR}				# unfiltered vcf output file with path
filtered_vcf=${tumor_bam//.coord.sorted.bam/_filtered.vcf.gz}
filtered_vcf=${filtered_vcf//$INPUT_DIR/$OUTPUT_DIR}					# filtered vcf output file with path

f1r2_file=${tumor_bam//.coord.sorted.bam/_f1r2.tar.gz}
f1r2_file=${f1r2_file//$INPUT_DIR/$OUTPUT_DIR}							# f1r2 output file with path

read_orientation_file=${tumor_bam//.coord.sorted.bam/_read_orientation.tar.gz}
read_orientation_file=${read_orientation_file//$INPUT_DIR/$OUTPUT_DIR}	# read orientation output file with path

tumor_pile_up_summary=${tumor_bam//.coord.sorted.bam/_pileup.table}
tumor_pile_up_summary=${tumor_pile_up_summary//$INPUT_DIR/$OUTPUT_DIR}	# pileup table output file with path

contamination_table=${tumor_bam//.coord.sorted.bam/_contamination.table}
contamination_table=${contamination_table//$INPUT_DIR/$OUTPUT_DIR}		# contamination table output file with path

segmentation_table=${tumor_bam//.coord.sorted.bam/_segmentation.table}
segmentation_table=${segmentation_table//$INPUT_DIR/$OUTPUT_DIR}		# segmentation table output file with path

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date    		: $(date)"
echo "Job name      		: $JOB_NAME"
echo "Job ID       			: $JOB_ID"  
echo "SGE TASK ID   		: $SGE_TASK_ID"
echo "Project       		: $proj"
echo "Index         		: $index"
echo "Normal BAM			: $normal_bam"
echo "Tumor BAM				: $tumor_bam"
echo "Normal Sample			: $normal_sample"
echo "Reference	Fasta		: $REF_FILE"
echo "Germline VCF			: $GERMLINE_REF"
echo "Panel of Normals VCF	: $PoN_VCF"
echo "Raw VCF				: $unfiltered_vcf"
echo "Filtered VCF 			: $filtered_vcf"
echo "F1R2 					: $f1r2_file"
echo "Read orientation		: $read_orientation_file"
echo "Tumor Pile Up 		: $tumor_pile_up_summary"
echo "Contamination Table 	: $contamination_table"
echo "Segmentation Table 	: $segmentation_table"
echo "=========================================================="

# You need to add path for Java
PATH=$HOME/NGSTools/jdk-23.0.1/bin/:$PATH

# # Find somatic mutations using Mutect2
# ~/NGSTools/gatk-4.6.1.0/gatk Mutect2 \
# --reference $REF_FILE \
# --input $tumor_bam \
# --germline-resource $GERMLINE_REF \
# --panel-of-normals $PoN_VCF \
# --f1r2-tar-gz $f1r2_file \
# --annotation ClippingRankSumTest \
# --annotation DepthPerSampleHC \
# --annotation MappingQualityRankSumTest \
# --annotation MappingQualityZero \
# --annotation QualByDepth \
# --annotation ReadPosRankSumTest \
# --annotation RMSMappingQuality \
# --annotation FisherStrand \
# --annotation MappingQuality \
# --annotation DepthPerAlleleBySample \
# --annotation Coverage \
# --output $unfiltered_vcf

#--input $normal_bam \
#-normal normal_sample_name \
#--intervals $INTERVAL_BED

# # Pass this raw data to LearnReadOrientationModel
# ~/NGSTools/gatk-4.6.1.0/gatk LearnReadOrientationModel \
# --input $f1r2_file \
# --output $read_orientation_file

# Run GetPileupSummaries to summarize read support for a set number of known variant sites.
# --intervals argument is usuallty the same as --variant argument
# However, the mouse germline reference was obtained from WGS and we are doing WES analysis
# So, I think using the padded.bed is sufficient and will save time (double check results..to be done)
~/NGSTools/gatk-4.6.1.0/gatk GetPileupSummaries \
--input $tumor_bam \
--variant $GERMLINE_REF_AF_ONLY \
--intervals $GERMLINE_REF_AF_ONLY \
--output $tumor_pile_up_summary
#--intervals $INTERVAL_BED
#-matched normal-pileups.table

# Calculate the fraction of reads coming from cross-sample contamination
~/NGSTools/gatk-4.6.1.0/gatk CalculateContamination \
--input $tumor_pile_up_summary \
--tumor-segmentation $segmentation_table
--output $contamination_table
#--matched-normal $normal_pileup_summary \

# Pass the learned read orientation model to FilterMutectCalls
# with the -ob-priors argument
~/NGSTools/gatk-4.6.1.0/gatk FilterMutectCalls \
--variant $unfiltered_vcf \
--reference $REF_FILE \
--contamination-table $contamination_table \
--tumor-segmentation $segmentation_table \
--orientation-bias-artifact-priors $read_orientation_file \
--f-score-beta 1 \
--output $filtered_vcf

#--interval-padding,
#--intervals