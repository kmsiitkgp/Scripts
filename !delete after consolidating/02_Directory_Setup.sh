#!/bin/bash -l

#proj="scRNASeq_BBN_C57B6"      
#proj="scRNASeq_BBN_Rag"
#proj="scRNASeq_Chen"
#proj="scRNASeq_GSE164557"
#proj="scRNASeq_GSE217093"
#proj="scRNASeq_GSE222315"
#proj="scRNASeq_HRA003620"
#proj="scRNASeq_Jinfen"
#proj="scRNASeq_NA13_CDH12_C57B6"
#proj="scRNASeq_Jyoti"
proj="scRNASeq_Krizia"


#proj="scRNASeq_NA13_CDH12_Nude"
#proj="scRNASeq_GSE129845"
#proj="scRNASeq_Simon"

#proj="RNASeq_GSE75192"
#proj="EGAD00001007575"
#proj="EGAD00001003977"
#proj="RNASeq_Hany_Antigen"
#proj="RNASeq_Hany_CRISPR_LOY"
#proj="RNASeq_Vera"
#proj="RNASeq_Sandrine.Supriya"
#proj="RNASeq_Manish_22RV1_ARCaPM"
#proj="RNASeq_Manish_22RV1_Xenograft"
proj="RNASeq_Manish_"

#proj="CRISPR_Jinfen"
#proj="CRISPR_Lena"
#proj="CRISPR_Prince"

#proj="WES_Hany"

proj="Spatial_VisiumHD"
proj="Spatial_Bhowmick"

OUTPUT_DIR=$HOME/scratch/$proj

#**********SCRNA SEQ DIRECTORIES**********#
mkdir $OUTPUT_DIR
mkdir $OUTPUT_DIR/logs
mkdir $OUTPUT_DIR/diagnostics
mkdir $OUTPUT_DIR/cellranger
mkdir $OUTPUT_DIR/filt_feature_bc_matrix
mkdir $OUTPUT_DIR/raw_feature_bc_matrix
mkdir $OUTPUT_DIR/seurat
#mkdir $OUTPUT_DIR/raw_reads
#mkdir $OUTPUT_DIR/fastqc_results
#mkdir $OUTPUT_DIR/CITESeq_results
#mkdir $OUTPUT_DIR/raw_hto_bc_matrix
#mkdir $OUTPUT_DIR/diagnostics/margin1
#mkdir $OUTPUT_DIR/diagnostics/margin2
#mkdir $OUTPUT_DIR/results_demux
#mkdir $OUTPUT_DIR/results_demux/margin1
#mkdir $OUTPUT_DIR/results_demux/margin2
#mkdir $OUTPUT_DIR/results_demux/singlets
#mkdir $OUTPUT_DIR/results_cellphonedb
#mkdir $OUTPUT_DIR/results_pyscenic
#mkdir $OUTPUT_DIR/results_scvelo
#mkdir $OUTPUT_DIR/results_velocyto

#**********SPATIAL SEQ DIRECTORIES**********#
# NOTE: raw matrix has reads from whole slide while 
# filtered matrix has reads from areas under tissue ONLY.
# So, we ONLY focus on filtered matrix
mkdir $OUTPUT_DIR
mkdir $OUTPUT_DIR/logs
mkdir $OUTPUT_DIR/diagnostics
mkdir $OUTPUT_DIR/spaceranger_results
mkdir $OUTPUT_DIR/filt_feature_bc_matrix
mkdir $OUTPUT_DIR/results_seurat

# To run this script, use "sh $HOME/projects/scRNASeq/00_scRNASeq_Directory_Setup.sh"
# Alternatively, you can make the script executable using "chmod u+x ./projects/scRNASeq/00_scRNASeq_Directory_Setup.sh" and
# then use "./projects/scRNASeq/00_scRNASeq_Directory_Setup.sh" to run your script

#**********RNA SEQ DIRECTORIES**********#
# DO NOT create a indexed_ref_genome directory. STAR will create it.
# If you create it before STAR, STAR will give error.
mkdir $OUTPUT_DIR
mkdir $OUTPUT_DIR/logs
#mkdir $OUTPUT_DIR/raw_reads
mkdir $OUTPUT_DIR/fastqc_results
#mkdir $OUTPUT_DIR/cutadapt_results
mkdir $OUTPUT_DIR/alignment_results
mkdir $OUTPUT_DIR/count_results

#**********CHIP SEQ DIRECTORIES**********#
# DO NOT create a indexed_ref_genome directory. STAR will create it.
# If you create it before STAR, STAR will give error.
mkdir $OUTPUT_DIR
mkdir $OUTPUT_DIR/logs
#mkdir $OUTPUT_DIR/raw_reads
mkdir $OUTPUT_DIR/fastqc_results
mkdir $OUTPUT_DIR/alignment_results
mkdir $OUTPUT_DIR/count_results
mkdir $OUTPUT_DIR/Sambamba_results
mkdir $OUTPUT_DIR/MACS2_results
mkdir $OUTPUT_DIR/ChIPQC_results

#**********CRISPR DIRECTORIES**********#
mkdir $OUTPUT_DIR
mkdir $OUTPUT_DIR/logs
#mkdir $OUTPUT_DIR/raw_reads
mkdir $OUTPUT_DIR/fastqc_results
mkdir $OUTPUT_DIR/cutadapt_results
mkdir $OUTPUT_DIR/alignment_results
mkdir $OUTPUT_DIR/count_results
mkdir $OUTPUT_DIR/DEG_results

#**********WES DIRECTORIES**********#
mkdir $OUTPUT_DIR
mkdir $OUTPUT_DIR/logs
#mkdir $OUTPUT_DIR/raw_reads
mkdir $OUTPUT_DIR/fastqc_results
mkdir $OUTPUT_DIR/cutadapt_results
mkdir $OUTPUT_DIR/index_results
mkdir $OUTPUT_DIR/index_results/bowtie2
mkdir $OUTPUT_DIR/index_results/bwamem2
mkdir $OUTPUT_DIR/alignment_results
mkdir $OUTPUT_DIR/alignment_results/bowtie2
mkdir $OUTPUT_DIR/alignment_results/bwamem2