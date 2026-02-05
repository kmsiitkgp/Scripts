#!/bin/bash -l

#$ -N multiqc                   # Set Job Name
#$ -l mem_free=40G              # Request memory
#$ -j y                         # Merge standard output and standard error
#$ -cwd                         # Set current working directory
#$ -t 1-1                       # Submit an array job with n tasks; n=number of samples = number of fastq files/2 if paired end

proj="WES_Hany"
species="Mouse" 
INPUT_DIR=$HOME/scratch/$proj/
OUTPUT_DIR=$HOME/scratch/$proj/logs/

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date    : $(date)"
echo "Job name      : $JOB_NAME"
echo "Job ID        : $JOB_ID"  
echo "SGE TASK ID   : $SGE_TASK_ID"
echo "Index         : $index"
echo "Project       : $proj"
echo "Input Folder  : $INPUT_DIR"
echo "Output Folder : $OUTPUT_DIR"
echo "=========================================================="

# You need to add path for Conda
PATH=$HOME/miniconda3/bin/:$PATH
source activate MAGECK

# Run multiqc once after fastqc is complete on all samples
multiqc --force \
--filename MULTIQC_Report \
--outdir=$OUTPUT_DIR \
$INPUT_DIR

#--fullnames \































# # Find which Agilent kit was used for exon capture. In Hany/Lena case, SureSelect Human All Exon V6 was used.
# # Go to https://earray.chem.agilent.com/suredesign/ ; create account and log in. 
# # Find Designs -> SureSelect DNA --> Download the files corresponding to SureSelect Human All Exon V6 r2 (Cat # S07604514).
# # This has 4 bed files and Targets.text file.

# # https://earray.chem.agilent.com/suredesign/help/Target_enrichment_design_files_available_for_download.htm
# # AllTracks.bed : contains the following tracks:
# # The Target Regions track is identical to the track in the Regions BED file.
# # The Covered probes track is identical to the track in the Covered BED file.
# # The Padded Covered track contains the covered regions (from either the Covered or Covered_partial BED file) with 100 bp of padding added on each side of all regions.
# # The No Probes track contains any regions from the Target Regions track that are not included in the Covered probes track.

# # Regions.bed   : contains a single track of the target regions of interest that SureDesign used to select the probes.
# # Covered.bed 	: contains a single track of the genomic regions that are covered by one or more probes in the design. 
# # The fourth column of the file contains annotation information. You can use this file for assessing coverage metrics.
# # Targets.txt 	: contains a list of the target identifiers that you entered when creating the design.

# #********************RUN ONLY ONCE********************

# # # Dict file will be created in same location as .fa file but with .dict extension.
# # java -jar ~/NGSTools/picard.jar CreateSequenceDictionary \
# # --REFERENCE $REFERENCE_DIR/$REFERENCE

# # # Ensembl reference doesnt have chr prefix. So, the bam files also dont have chr prefix.
# # # The bed files from Agilent have chr prefix to each chromosome so we remove chr prefix.
# # cat $BED_DIR/$BAIT_BED | grep -o chr | wc -l#		# Number of lines having chr
# # cat $BED_DIR/$BAIT_BED | wc							# Number of lines, words, characters in before correction
# # cat $BED_DIR/$BAIT_BED | awk '{gsub(/chr/,""); print}' > $BED_DIR/$BAIT_BED.corrected
# # cat $BED_DIR/$BAIT_BED.corrected | wc				# Number of lines, words, characters in after correction
# # cat $BED_DIR/$TARGET_BED | grep -o chr | wc -l		# Number of lines having chr
# # cat $BED_DIR/$TARGET_BED | wc						# Number of lines, words, characters in before correction
# # cat $BED_DIR/$TARGET_BED | awk '{gsub(/chr/,""); print}' > $BED_DIR/$TARGET_BED.corrected
# # cat $BED_DIR/$TARGET_BED.corrected | wc				# Number of lines, words, characters in after correction

# # java -jar ~/NGSTools/picard.jar BedToIntervalList \
# # --INPUT $BED_DIR/$BAIT_BED.corrected \
# # --OUTPUT $BED_DIR/$BAIT_LIST \
# # --SEQUENCE_DICTIONARY $REFERENCE_DIR/$DICT

# # java -jar ~/NGSTools/picard.jar BedToIntervalList \
# # --INPUT $BED_DIR/$TARGET_BED.corrected \
# # --OUTPUT $BED_DIR/$TARGET_LIST \
# # --SEQUENCE_DICTIONARY $REFERENCE_DIR/$DICT

# # # index file will be created in same location as .fa file but with .fai extension 
# # samtools faidx $REFERENCE_DIR/$REFERENCE

# #********************RUN ONLY ONCE UNTIL PREVIOUS LINE********************

# java -jar $HOME/NGSTools/picard.jar CollectHsMetrics \
# --INPUT $taskinput \
# --OUTPUT $taskoutput \
# --INCLUDE_INDELS false \
# --NEAR_DISTANCE 250 \
# --PER_TARGET_COVERAGE $taskoutput_coverage \
# --REFERENCE_SEQUENCE $REFERENCE_DIR/$REFERENCE \
# --BAIT_INTERVALS $BED_DIR/$BAIT_LIST \
# --TARGET_INTERVALS $BED_DIR/$TARGET_LIST

# # Get exon length for each chromosome from gtf file
# # Column number	Column name		Details
# # 1				seqname			name of the chromosome or scaffold; chromosome names can be given with or without the ‘chr’ prefix.
# # 2				source			name of the program that generated this feature, or the data source (database or project name)
# # 3				feature			feature type name, e.g. Gene, Exon
# # 4				start			Start position of the feature, with sequence numbering starting at 1.
# # 5				end	End 		position of the feature, with sequence numbering starting at 1.
# # 6				score			A floating point value.
# # 7				strand			defined as + (forward) or - (reverse).
# # 8				frame			One of ‘0’, ‘1’ or ‘2’. ‘0’ indicates that the first base of the feature is the first base of a codon, ‘1’ that the second base is the first base of a codon, and so on..
# # 9				attribute		A semicolon-separated list of tag-value pairs, providing additional information about each feature.
# # cat ~/NGSTools/Reference_Genomes/Human/*.gtf | awk '$3=="exon" && $1!~/[A-W]/' | cut -f 1,4,5 | awk '{$4=$3-$2}1' | awk '{sum+=$4} END {print sum+0}'

# # # 1st column: Chr; 2nd column: genomic co-ordinates; 3rd column: coverage
# # awk '{sum+=$3} END { print "Average = ",sum/NR}' $taskoutput_depth
# # awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' $taskoutput_depth
# # awk '$3==0 {sum++} END {print sum+0}' $taskoutput_depth
# # awk '$3>=4 {sum++} END {print "Fraction_of_target_covered_with_at_least_4x =",sum+0}' $taskoutput_depth
# # awk '$3>=10 {sum++} END {print "Fraction_of_target_covered_with_at_least_10x =",sum+0}' $taskoutput_depth
# # awk '$3>=20 {sum++} END {print "Fraction_of_target_covered_with_at_least_20x =",sum+0}' $taskoutput_depth
# # awk '$3>=50 {sum++} END {print "Fraction_of_target_covered_with_at_least_50x =",sum+0}' $taskoutput_depth
# # awk '$3>=100 {sum++} END {print "Fraction_of_target_covered_with_at_least_100x =",sum+0}' $taskoutput_depth

