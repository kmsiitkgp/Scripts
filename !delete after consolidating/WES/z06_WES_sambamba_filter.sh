#!/bin/bash -l

#$ -N sambamba                   # Set Job Name
#$ -l mem_free=40G              # Request memory
#$ -j y                         # Merge standard output and standard error
#$ -cwd                         # Set current working directory
#$ -t 1-2                       # Submit an array job with n tasks; n=number of samples = number of fastq files/2 if paired end

proj="WES_Lena"
INPUT_DIR=$HOME/scratch/$proj/alignment_results  
OUTPUT_DIR=$INPUT_DIR

# Create an array of input files
input=($INPUT_DIR/*.sam)

# Use the SGE_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value
index=$(($SGE_TASK_ID-1))
taskinput=${input[$index]}
taskoutput=${taskinput//$INPUT_DIR/$OUTPUT_DIR}
taskoutput=${taskoutput//.sam/}
taskoutput_sam=$taskoutput.sam
taskoutput_bam=$taskoutput.bam
taskoutput_collate_bam=$taskoutput.collate.bam
taskoutput_fixmate_bam=$taskoutput.fixmate.bam
taskoutput_pos_sorted_bam=$taskoutput.pos.sorted.bam
taskoutput_pos_sorted_no_dup_bam=$taskoutput.pos.sorted.no.dup.bam
taskoutput_pos_sorted_no_dup_filtered_bam=$taskoutput.pos.sorted.no.dup.filtered.bam
taskoutput_coverage=$taskoutput.coverage
taskoutput_depth=$taskoutput.depth
taskoutput_bed=$taskoutput.bed
taskoutput_pos_sorted_no_dup_filtered_Y_bam=$taskoutput.pos.sorted.no.dup.filtered.Y.bam
taskoutput_bed_Y=$taskoutput.bed.Y

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date    : $(date)"
echo "Job name      : $JOB_NAME"
echo "Job ID        : $JOB_ID"  
echo "SGE TASK ID   : $SGE_TASK_ID"
echo "Index         : $index"
echo "Task input 1  : $taskinput"
echo "Task output   : $taskoutput"
echo "Project       : $proj"
echo "Output Folder : $OUTPUT_DIR"
echo "=========================================================="

# You need to add path for Conda
PATH=$HOME/miniconda3/bin/:$PATH
source activate MAGECK
PATH=$HOME/NGSTools/sambamba-1.0.1-linux-amd64-static/:$PATH

# Convert to bam
sambamba view \
--sam-input \
--format=bam \
--output-filename=$taskoutput_bam \
$taskoutput_sam

# Mark duplicate reads. sambamba markdup gets stuck. so, use samtools.
# https://www.htslib.org/algorithms/duplicate.html
# The - symbol in samtools it to tell the samtools program to take the input from pipe
# Step (1) groups the reads by read name in the bam file. This puts the read pairs close together so that fixmate can work. 
# Step (2) repairs mate information for read pairs. For markdup the most important parts are the MC and ms tags. 
# MC adds the mate cigar string and is used to calculate the mate position for duplicate matching. 
# The ms tag is the mate quality score and is used to calculate originality. 
# Step (3) sorts the reads into position order to make it reads for markdup. 
# Step (4) is the duplicate marking.
samtools collate -o $taskoutput_collate_bam $taskoutput_bam
samtools fixmate -m $taskoutput_collate_bam $taskoutput_fixmate_bam
samtools sort -o $taskoutput_pos_sorted_bam $taskoutput_fixmate_bam
samtools markdup $taskoutput_pos_sorted_bam $taskoutput_pos_sorted_no_dup_bam

# Get stats pre-filtering
echo "Stats Before Filtering"
sambamba flagstat $taskoutput_pos_sorted_no_dup_bam

# Filter unmapped reads, multi-mapped reads
# https://github.com/biod/sambamba/wiki/%5Bsambamba-view%5D-Filter-expression-syntax
# https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
# https://github.com/hbctraining/Intro-to-ChIPseq-flipped/blob/main/lessons/05_filtering_BAM_files.md
# https://avrilomics.blogspot.com/2013/05/bam-and-sam-flags.html
# If you read bowtie2 link above XS is assigned only if read has multiple alignments   --filter="[XS] == null" \

# Filter and retain ONLY concordant reads (unique + multiple)
# If aligned using bowtie2, you can use --filter="[YT] == 'CP'" 
# If another aligner was used, you can use --filter="proper_pair"
# METHOD 1:
# sambamba view  \
# --with-header \
# --format=bam \
# --filter="[YT] == 'CP'" \
# --output-filename=$taskoutput_sorted_filtered_bam \
# $taskoutput_sorted_bam

# METHOD 2: If you used an aligner other than bowtie2
# samtools view --bam -f 83 $taskoutput_sorted_bam > test1.bam
# samtools view --bam -f 99 $taskoutput_sorted_bam > test2.bam
# samtools view --bam -f 163 $taskoutput_sorted_bam > test3.bam
# samtools view --bam -f 147 $taskoutput_sorted_bam > test4.bam
# sambamba merge $taskoutput_sorted_filtered_bam test1.bam test2.bam test3.bam test4.bam

# Filter and retain concordant reads (unique + multiple) and discordant reads (unique)
# If aligned using bowtie2, you can use --filter="[YT] == 'CP' or [YT] == 'DP'"  
# sambamba view \
# --with-header \
# --format=bam \
# --filter="[YT] == 'CP' or [YT] == 'DP'" \
# --output-filename=$taskoutput_sorted_filtered_bam \
# $taskoutput_sorted_bam

# Filter and retain concordant reads (unique + multiple) and discordant reads (unique + multiple). 
# All unmapped and singly mapped reads and duplicates are ONLY removed in this case.
sambamba view \
--with-header \
--format=bam \
--filter="([YT] == 'CP' or [YT] == 'DP' or [YT] == 'UP') and not (unmapped or mate_is_unmapped or duplicate)" \
--output-filename=$taskoutput_pos_sorted_no_dup_filtered_bam \
$taskoutput_pos_sorted_no_dup_bam

# Get stats post filtering
echo "Stats After Filtering"
sambamba flagstat $taskoutput_pos_sorted_no_dup_filtered_bam

# Generate index files
sambamba index $taskoutput_pos_sorted_no_dup_filtered_bam

# Get reads mapped to each chr
# NOTE: bam file MUST have its index file before running this step.
# NOTE: output is TAB-delimited with each line consisting of reference sequence name, sequence length, # mapped 
# read-segments and # unmapped read-segments (singletons). This may count reads multiple times if they are mapped # more than once or in multiple fragments.
samtools idxstats $taskoutput_pos_sorted_no_dup_filtered_bam

# Calculate coverage at each chromosome
samtools coverage \
--output=$taskoutput_coverage \
$taskoutput_pos_sorted_no_dup_filtered_bam

# Calculate depth at each chromosome
samtools depth \
-o $taskoutput_depth \
$taskoutput_pos_sorted_no_dup_filtered_bam

# Convert bam to bed
bedtools bamtobed \
-i $taskoutput_pos_sorted_no_dup_filtered_bam > $taskoutput_bed

# Y specific files.
# If you check samtools idxstats, chr 1 is corresponds to line 0 i.e ref_id =0.
# Y chr is on 25th line for human and 22nd line for mouse. 
# So, ref_id=24 in human and ref_id=21 in mouse for Y chr. 
# After filtering, convert the bam to sam and view it to confirm.
# samtools view -h -o out.sam in.bam
# cat out.sam | head -300
sambamba view \
--with-header \
--format=bam \
--filter="(ref_id == 24)" \
--output-filename=$taskoutput_pos_sorted_no_dup_filtered_Y_bam \
$taskoutput_pos_sorted_no_dup_filtered_bam

bedtools bamtobed \
-i $taskoutput_pos_sorted_no_dup_filtered_Y_bam > $taskoutput_bed_Y



