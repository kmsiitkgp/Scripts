#!/bin/bash -l

#$ -N GATK_PoN                  # Set Job Name
#$ -l mem_free=96G              # Request memory
#$ -j y                         # Merge standard output and standard error
#$ -cwd                         # Set current working directory
#$ -t 1-1                       # Submit an array job with 1 tasks

species="Mouse"
proj="WES_Hany"
INPUT_DIR=$HOME/scratch/$proj/alignment_results  
OUTPUT_DIR=$HOME/scratch/$proj/variant_results

REF_FILE=$HOME/NGSTools/Reference_Genomes/$species/*.fa
INTERVAL_BED=$HOME/NGSTools/WES/S0276129_SureSelect_Mouse_AllExon_V1_mm39.bed

VCF_DB=$OUTPUT_DIR/VCF_DB	# Workspace for GenomicsDB. GATK will create it. You MUST NOT create it.
PoN_VCF=$OUTPUT_DIR/PoN.vcf.gz

# Create an array of normal samples that will be used to generate normal vcfs
# and Panel of Normals (PoN)
NORMAL_SAMPLES=(MB49.coord.sorted.bam)
NORMAL_SAMPLES=(${NORMAL_SAMPLES[*]/#/$INPUT_DIR/})
NORMAL_VCFS=(${NORMAL_SAMPLES[*]//.coord.sorted.bam/.vcf.gz})
NORMAL_VCFS=(${NORMAL_VCFS[*]//$INPUT_DIR/$OUTPUT_DIR})

# PoN is used to eliminate sequencing and other technical artifacts from the data
# matched normal or a germline resource is used to eliminate normal germline calls from your data 

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date    : $(date)"
echo "Job name      : $JOB_NAME"
echo "Job ID        : $JOB_ID"  
echo "SGE TASK ID   : $SGE_TASK_ID"
echo "Project       : $proj"
echo "Reference		: $REF_FILE"
echo "=========================================================="

# You need to add path for Java
PATH=$HOME/NGSTools/jdk-23.0.1/bin/:$PATH

### STEP 1: Run Mutect2 in tumor only mode for each normal sample
for index in ${!NORMAL_SAMPLES[*]}
do
normal_sample=${NORMAL_SAMPLES[$index]}
normal_vcf=${NORMAL_VCFS[$index]}

~/NGSTools/gatk-4.6.1.0/gatk Mutect2 \
--reference $REF_FILE \
--max-mnp-distance 0 \
--input $normal_sample \
--output $normal_vcf

echo "=========================================================="
echo "Normal Sample	: $normal_sample"
echo "Normal VCF	: $normal_vcf"
echo "=========================================================="
done

### STEP 2: Import normal vcfs into GenomicsDB workspace
# --genomicsdb-workspace-path MUST point to a "non-existent or empty directory"
# -L argument is MANDATORY 
# If you have more than 1 normal, add all of them like
# --variant ${NORMAL_VCFS[1]} \
# --variant ${NORMAL_VCFS[2]} \
# --variant ${NORMAL_VCFS[3]}

# For WES, use interval list from Agilent website. 
# https://earray.chem.agilent.com/suredesign/search/entity.htm
# Login --> Find Designs --> SureSelectDNA --> Agilent Catalog --> 
# Filter species --> Download appropriate kit that was used for WES
# https://earray.chem.agilent.com/suredesign/help/Design_analysis_using_tracks.htm
# The intervals are provided in padded.bed file from Agilent but it uses mm9 genome annotation. 
# However, we performed alignment using mm39 genome annotation.
# So, convert the annotations in bed file using https://genome.ucsc.edu/cgi-bin/hgLiftOver 
# We perform sequential lift mm9 -> mm10 -> mm39 since mm9 to mm39 is not possible in one step.
# Open the padded.bed in excel and copy the 1st 3 columns to UCSC tool and convert to mm10.
# Download the converted intervals and repeat for mm10 --> mm39
# You can load into IGV, the mm9 padded.bed, mm39 padded.bed and use mm9 reference.
# You will notice the exons alignment perfectly for mm9 padded.bed but not mm39 padded.bed
# Similarly, you can load into IGV, the mm9 padded.bed, mm39 padded.bed and use mm39 ref.
# You will notice the exons alignment perfectly for mm39 padded.bed but not mm9 padded.bed

# Formating the bed file
# In the bed file from UCSC, chromosomes are listed as chr1, chr2, chrX, chrY, chrM.
# However, the reference files we used from ENSEMBL lists them as 1,2,X,Y,MT.
# So, remove the chr from the bed file and rename chrM as MT before using in GATK.

~/NGSTools/gatk-4.6.1.0/gatk GenomicsDBImport \
--genomicsdb-workspace-path $VCF_DB \
--reference $REF_FILE \
--intervals $INTERVAL_BED \
--interval-padding 100 \
--merge-input-intervals true \
--variant ${NORMAL_VCFS[0]}
# --variant ${NORMAL_VCFS[1]}
# --variant ${NORMAL_VCFS[2]}

### STEP 3: Combine the normal calls to create a panel of normals
# germline resoure is not available for mouse. So, skip --germline-resource

~/NGSTools/gatk-4.6.1.0/gatk CreateSomaticPanelOfNormals \
--reference $REF_FILE \
--variant gendb://$VCF_DB \
--output $PoN_VCF