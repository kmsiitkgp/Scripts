#!/bin/bash
#SBATCH --job-name=Biomodal_analyze
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --output=%x_%A_%a.log
#SBATCH --array=1-8  # Run in batches of 8

#===============================================================================
# 🚀 Environment
#===============================================================================

source "/home/kailasamms/miniconda3/etc/profile.d/conda.sh"
conda activate omics
module load biomodal/2.0.0
module load apptainer/1.0.1

#===============================================================================
# 🔎 Detect input files
#===============================================================================

# FASTQ MUST be stored within nf-input folder within $IN_DIR
IN_DIR="/common/bhowmicknlab/Biomodal/Kevin_Carotuximab/input/"
BAM_DIR="/home/kailasamms/analysis/Biomodal/Kevin_Carotuximab/02.bams/"

# --- Detect raw R1 FASTQ files ---
echo "⚙️ Looking for R1 FASTQ files in ${IN_DIR}"
shopt -s nullglob
R1_FASTQ_ARRAY=("${IN_DIR%/}/"*_R1_*.f*q.gz)
shopt -u nullglob
echo "Found ${#R1_FASTQ_ARRAY[@]} R1 FASTQ files"

# --- Extract sample names and prefixes ---
INITIAL_SAMPLE_NAMES=()
INITIAL_PREFIXES=()
for f in "${R1_FASTQ_ARRAY[@]}"; do
    # Extract the base filename (remove directory path)
    base=$(basename "$f")
    # Optionally, strip common FASTQ extensions
    base=${base%.fastq.gz}
    base=${base%.fq.gz}
    # PREFIX = remove last two fields: _R1_001
    prefix="${base%_R1*}"
    # Extract SAMPLE_NAME = everything before first "_S"
    sample="${prefix%%_S*}"
    INITIAL_PREFIXES+=("$prefix")
    INITIAL_SAMPLE_NAMES+=("$sample")
done

# Remove duplicates while preserving order
SAMPLE_NAMES=($(printf "%s\n" "${INITIAL_SAMPLE_NAMES[@]}" | awk '!seen[$0]++'))
PREFIXES=($(printf "%s\n" "${INITIAL_PREFIXES[@]}" | awk '!seen[$0]++'))

# --- Select sample for this array task ---
ARRAY_INDEX=$((SLURM_ARRAY_TASK_ID - 1))      # SLURM_TASK_ID is 1-based, bash array index is 0-based
SAMPLE_NAME="${SAMPLE_NAMES[$ARRAY_INDEX]}"
PREFIX="${PREFIXES[$ARRAY_INDEX]}"

if [ -z "${SAMPLE_NAME}" ]; then
    echo "❌ No sample found for task ${SLURM_ARRAY_TASK_ID}. Aborting."
    exit 0
fi

#===============================================================================
# ⚙️ Arguments
#===============================================================================

BIOMODAL_INSTANCE_DIRECTORY="/home/kailasamms/resources/biomodal/"
BIOMODAL_WORK_DIRECTORY="/home/kailasamms/scratch/biomodal/"
INPUT_DIRECTORY="/common/bhowmicknlab/Biomodal/Kevin_Carotuximab/input/"
OUTPUT_DIRECTORY="/home/kailasamms/analysis/Biomodal/Kevin_Carotuximab/01.duet"
MODE="6bp"
DEPTH="super_seq"

#===============================================================================
# 📋 Job Info
#===============================================================================

echo "=========================================================="
printf "Start date      : %s\n" "$(date)"
printf "Job name        : %s\n" "${SLURM_JOB_NAME}"
printf "Job ID          : %s\n" "${SLURM_JOB_ID}"
printf "Array Task ID   : %s\n" "${SLURM_ARRAY_TASK_ID}"
printf "Cores           : %s\n" "${SLURM_CPUS_PER_TASK}"
echo "=========================================================="
printf "Input dir       : %s\n" "${INPUT_DIRECTORY}"
printf "Output dir      : %s\n" "${OUTPUT_DIRECTORY}"
printf "Sample          : %s\n" "${SAMPLE_NAME}"
printf "Prefix          : %s\n" "${PREFIX}"
printf "Mode            : %s\n" "${MODE}"
printf "Depth           : %s\n" "${DEPTH}"
echo "=========================================================="

#===============================================================================
# 🚀 Run Analysis
#===============================================================================

biomodal run duet \
  --instance-directory "${BIOMODAL_INSTANCE_DIRECTORY%/}" \
  --work-dir "${BIOMODAL_WORK_DIRECTORY%/}" \
  --input-path "${INPUT_DIRECTORY%/}" \
  --output-path "${OUTPUT_DIRECTORY%/}" \
  --mode "${MODE}" \
  --tag "${SAMPLE_NAME}" \
  --run-name "${SAMPLE_NAME}" \
  --chg-chh-contexts \
  --additional-params "run_by_default.fastqc.raw=true" \
  --additional-params "publish_by_default.prelude=false" \
  --additional-params "epiquant.primary_chromosomes_only=true" \
  --additional-params "call_somatic_variants=true" \
  --additional-params "lib_prefix=${PREFIX}"

#  --additional-profile "${DEPTH}" \

#===============================================================================
# Symlink bams to new folder for easy access
#===============================================================================
mkdir -p "${BAM_DIR%/}"

find "${OUTPUT_DIRECTORY%/}" -name "*genome.GRCh38Decoy_primary_assembly.dedup.bam" | while read bam; do
    ln -sf "$bam" "${BAM_DIR%/}/$(basename "$bam")"
done

find "${OUTPUT_DIRECTORY%/}" -name "*genome.GRCh38Decoy_primary_assembly.dedup.bam.bai" | while read bai; do
    ln -sf "$bai" "${BAM_DIR%/}/$(basename "$bai")"
done