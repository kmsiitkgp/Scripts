#!/bin/bash

# =========================================================================================
# NEXTFLOW PIPELINE LAUNCHER
# =========================================================================================
# This script:
#   1. Reads paths from project_info.yaml
#   2. Auto-detects physical paths (resolving symlinks)
#   3. Builds Singularity bind mounts dynamically
#   4. Launches the Nextflow pipeline with proper configuration
#
# Usage:
#   bash run_nextflow.sh
#
# Prerequisites:
#   - project_info.yaml configured with correct paths
#   - Nextflow and Singularity modules available
#   - Access to HPC cluster (if using -profile sge)
#
# Resume a failed run:
#   Just run the script again - Nextflow automatically resumes from last checkpoint
# =========================================================================================

# -----------------------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------------------

#project="BBN_C57BL6"
#project="Xinyi"                  # Project name (used in output directory naming)
project="Manish_22Rv1_ARCaPM"
#project="Manish_22Rv1_Xeno"
#project="Manish_22Rv1_Xeno2"

# Nextflow directory (where main.nf, nextflow.config, modules/, workflows/, projects/ are located)
NF_DIR="$HOME/scripts/nextflow"

# YAML file with project specific settings in projects folder
YAML="${NF_DIR}/projects/${project}_project_info.yaml"

# -----------------------------------------------------------------------------------------
# YAML PREPROCESSING
# -----------------------------------------------------------------------------------------
# YAML parsers don't handle tabs well - convert to spaces
# This prevents "mapping values are not allowed here" errors

# Replace tabs with 4 spaces
sed -i 's/\t/    /g' "${YAML}"

# Verify no tabs remain (should print nothing if clean)
grep $'\t' "${YAML}"

# -----------------------------------------------------------------------------------------
# EXTRACT PATHS FROM YAML
# -----------------------------------------------------------------------------------------
# Parse YAML to get directory paths
# These may be logical paths (symlinks) that need to be resolved for Singularity

REF_DIR=$(grep "^[[:space:]]*ref_dir[[:space:]]*:" $YAML | cut -d':' -f2- | xargs | tr -d '"')
READ_DIR=$(grep "^[[:space:]]*read_dir[[:space:]]*:" $YAML | cut -d':' -f2- | xargs | tr -d '"')
BASE_DIR=$(grep "^[[:space:]]*base_dir[[:space:]]*:" $YAML | cut -d':' -f2- | xargs | tr -d '"')
WORK_DIR=$(grep "^[[:space:]]*work_dir[[:space:]]*:" $YAML | cut -d':' -f2- | xargs | tr -d '"')
CACHE_DIR=$(grep "^[[:space:]]*cache_dir[[:space:]]*:" $YAML | cut -d':' -f2- | xargs | tr -d '"')

# -----------------------------------------------------------------------------------------
# RESOLVE PHYSICAL PATHS
# -----------------------------------------------------------------------------------------
# Singularity needs physical paths (not symlinks) for bind mounts
# Example: /scratch -> /gpfs/scratch01 (actual physical location)

P_REF=$(readlink -f "$REF_DIR")
P_READ=$(readlink -f "$READ_DIR")
P_BASE=$(readlink -f "$BASE_DIR")
P_WORK=$(readlink -f "$WORK_DIR")
P_CACHE=$(readlink -f "$CACHE_DIR")

# -----------------------------------------------------------------------------------------
# BUILD SINGULARITY BIND MOUNTS
# -----------------------------------------------------------------------------------------
# Strategy: Bind root directories (e.g., /hpc, /scratch) instead of individual paths
# This covers all subdirectories and avoids complex bind mount conflicts
#
# Why this matters:
#   - If /scratch is a symlink to /gpfs/scratch01, both need to be accessible
#   - Binding parent directories ensures all child paths work
#   - Prevents "No such file or directory" errors inside containers

# Collect all paths (logical and physical)
PATHS_TO_BIND="$REF_DIR $READ_DIR $BASE_DIR $WORK_DIR $CACHE_DIR $P_REF $P_READ $P_BASE $P_WORK $P_CACHE"

# Extract unique root directories (first 2 levels: /hpc, /scratch, /gpfs, etc.)
UNIQUE_ROOTS=$(for p in $PATHS_TO_BIND; do echo "$p" | cut -d'/' -f1-2; done | sort -u | grep '^/')

# Use root directories as bind paths (simplest and most robust approach)
ALL_BIND_PATHS="$UNIQUE_ROOTS"

# Alternative (if you need more specific binds):
# ALL_BIND_PATHS="$UNIQUE_ROOTS $PATHS_TO_BIND"

# Build Singularity bind flags: --bind /hpc --bind /scratch --bind /gpfs
BIND_FLAGS=""
for path in $ALL_BIND_PATHS; do
    if [ -d "$path" ]; then
        BIND_FLAGS+="--bind $path "
    fi
done

echo "Singularity bind flags: $BIND_FLAGS"

# -----------------------------------------------------------------------------------------
# DISPLAY PATH MAPPINGS
# -----------------------------------------------------------------------------------------
# Show how paths map between host (physical) and container (logical)
# Useful for debugging "file not found" issues

echo "---------------------------------------------------------------------------------------------------------------"
printf "%-16s : %-60s -> %s\n" "Mapping:" "PHYSICAL (Host)" "LOGICAL (Container)"
echo "---------------------------------------------------------------------------------------------------------------"
printf "%-16s : %-60s -> %s\n" "Ref Genomes"    "$P_REF"  "$REF_DIR"
printf "%-16s : %-60s -> %s\n" "Raw Reads"      "$P_READ" "$READ_DIR"
printf "%-16s : %-60s -> %s\n" "Project Base"   "$P_BASE" "$BASE_DIR"
printf "%-16s : %-60s -> %s\n" "Temp Work Dir"  "$P_WORK" "$WORK_DIR"
printf "%-16s : %-60s -> %s\n" "Image Cache"    "$P_CACHE" "$CACHE_DIR"
echo "---------------------------------------------------------------------------------------------------------------"

# -----------------------------------------------------------------------------------------
# CLEANUP NEXTFLOW FILES
# -----------------------------------------------------------------------------------------
# Remove problematic characters that cause parsing errors

# Enable nullglob so empty directories don't cause errors
shopt -s nullglob

# Define file list once to keep it DRY
NF_FILES=(
  "${NF_DIR}"/*.nf
  "${NF_DIR}"/*.config
  "${NF_DIR}"/*.txt
  "${NF_DIR}"/*.sh
  "${NF_DIR}"/modules/*.nf
  "${NF_DIR}"/workflows/*.nf
  "${NF_DIR}"/projects/*.yaml
)

# 1. Remove UTF-8 BOM
# 2. Strip trailing whitespace
# 3. Convert tabs to 4 standard spaces
# Using a loop handles the file list more reliably across different OS versions
for file in "${NF_FILES[@]}"; do

    # # 1. Convert to UTF-8 (and remove ANSI/ISO encoding)
    # # iconv -c ignores invalid characters to prevent crashing
    # if command -v iconv >/dev/null; then
        # iconv -f WINDOWS-1252 -t UTF-8 "$file" > "${file}.tmp" 2>/dev/null && mv "${file}.tmp" "$file"
    # fi

    # 2. Remove UTF-8 BOM (Byte Order Mark) if present
    # BOM appears as: 0xEF 0xBB 0xBF at file start (from some text editors)
    # Causes error: "Invalid character at start of file"
    sed -i '1s/^\xEF\xBB\xBF//' "$file"

    # # 3. Convert Windows Line Endings (CRLF) to Linux (LF)
    # # Using 'dos2unix' if available, otherwise 'sed'
    # if command -v dos2unix >/dev/null; then
        # dos2unix -q "$file"
    # else
        # sed -i 's/\r$//' "$file"
    # fi

    # 4. Clean up whitespace: Remove trailing spaces & convert tabs to 4 spaces
    # - Strip trailing whitespace (prevents unnecessary git diffs)
    # - Convert tabs to 4 spaces (consistent indentation)
    sed -i 's/[[:space:]]*$//' "$file"
    sed -i 's/\t/    /g' "$file"
done

echo "All files are now Unix-friendly (UTF-8, LF, No BOM)."

# Disable nullglob
shopt -u nullglob

# -----------------------------------------------------------------------------------------
# EXECUTION SETUP
# -----------------------------------------------------------------------------------------

# Load Nextflow and Singularity from HPC environment modules
# Versions may vary by cluster - adjust as needed
module load nextflow/24.10.5
module load singularity-apptainer/1.1.8

# Define a project-specific workspace within the scratch directory.
# This isolation is CRITICAL for:
#   1. Running multiple projects simultaneously without "Database Lock" errors.
#   2. Keeping logs (.nextflow.log) and metadata (.nextflow/) organized by project.
#   3. Preventing different runs from interfering with each other's cache/resume.
mkdir -p "${WORK_DIR}/${project}"

# Move into the project workspace.
# By changing directory here, Nextflow will automatically create its internal
# './work' folder and session database inside this unique project path.
cd "${WORK_DIR}/${project}"

# Remove any stale lock files from runs to prevent "Run name already used" error
[ -f ".nextflow/lock" ] && rm -f ".nextflow/lock"

# Create a unique run name using a timestamp.
# This avoids the "Run name already used" error if a previous session crashed
# or didn't release its lock. It allows for seamless -resume attempts while
# keeping a unique record for every try in the 'nextflow log'.
RUN_NAME="${project}_$(date +%Y-%m-%d_%H-%M-%S)"

echo "--------------------------------------------------------------------------------"
echo "Lauch Dir : $(pwd)"
echo "Work Dir  : ${WORK_DIR}/${project}/work"
echo "Project   : ${project}"
echo "YAML      : ${YAML}"
echo "Run Name  : ${RUN_NAME}"
echo "--------------------------------------------------------------------------------"

# -----------------------------------------------------------------------------------------
# PIPELINE EXECUTION
# -----------------------------------------------------------------------------------------
# -params-file    : Injects YAML settings (FASTQ paths, species, etc.)
# -name           : Labels the run in 'nextflow log' and SGE (highly recommended)
# -work-dir       : Explicitly sets the heavy data directory. We define this for SAFETY
#                   to ensure task data stays in scratch even if the config changes.
# -profile        : Sets executor to 'sge' (Sun Grid Engine)
# -resume         : Uses cached results to skip successfully completed steps
# --dynamic_binds : Passes calculated physical paths to Singularity containers
# --stop.after    : RENAME_FASTQS, GENERATE_MD5, VALIDATE_INPUT, FASTQC_RAW, STAR_INDEX,
#                 : EXTRACT_GENTROME, SALMON_INDEX, RSEQC_BED, SALMON_QUANT,
#                 : STAR_ALIGN, SAMBAMBA_PREP, RSEQC, MULTIQC, CELLRANGER_COUNT

#screen -S "${project}" -d -m \
nextflow \
    -log "${WORK_DIR}/${project}.nextflow.log" \
    run "${NF_DIR}/main.nf" \
    -params-file "${YAML}" \
    -name "${RUN_NAME}" \
    -work-dir "${WORK_DIR}/${project}/work" \
    -profile sge \
    -resume \
    --dynamic_binds "$BIND_FLAGS" \
    --stop_after "END"
#    --stop_after "GENERATE_MD5"
#-preview

# See logs live
# tail -f "${WORK_DIR}/${project}.nextflow.log"

# Reattach screen session if needed
# screen -r "${project}"

# Detach screen session if needed
# Press Ctrl+A fist, release and press D

# -----------------------------------------------------------------------------------------
# DELETE WORD DIRS FROM FAILED RUNS (RUN AFTER SUCCESSFUL COMPLETION OF PROJECT)
# -----------------------------------------------------------------------------------------

# Dry run first — see what would be deleted
#nextflow log -q -f status | grep -v 'OK' | cut -f1 | xargs -n 1 -I {} echo nextflow clean {} -f

# If looks good, actually delete
#nextflow log -q -f status | grep -v 'OK' | cut -f1 | xargs -n 1 -I {} nextflow clean {} -f


# -----------------------------------------------------------------------------------------
# NOTES ON SINGULARITY BIND MOUNT CHALLENGES
# -----------------------------------------------------------------------------------------
# Problem: Nextflow needs to pass bind mounts to Singularity, but has limited options
#
# Approaches that DIDN'T work:
#   ✗ export NXF_SINGULARITY_RUNOPS="$BIND_FLAGS"
#   ✗ --singularity.runOptions "${BIND_FLAGS}"
#   ✗ --singularity.runOptions "'${BIND_FLAGS}'"
#   ✗ Hardcoding in nextflow.config: runOptions = "--bind /scratch ..."
#
# Solution that WORKS:
#   ✓ Pass bind flags as custom parameter: --dynamic_binds "$BIND_FLAGS"
#   ✓ In nextflow.config: runOptions = "${params.dynamic_binds}"
#
# Why this works:
#   - Custom parameters are evaluated at runtime
#   - Allows dynamic path detection instead of hardcoding
#   - Script calculates binds, passes to Nextflow, Nextflow passes to Singularity

# -----------------------------------------------------------------------------------------
# DEBUGGING TIPS
# -----------------------------------------------------------------------------------------

# Check for UTF-8 BOM in files:
# head -c 3 "${NF_DIR}/main.nf" | od -c
# If output shows: 0000000 357 273 277 # ! /
#   → BOM present, needs removal
# If output shows: 0000000 # ! /
#   → No BOM, file is clean

# View which work directories correspond to which processes:
# nextflow log <run_name> -f name,workdir,status,duration

# List all runs:
# nextflow log

# Clean up old run metadata:
# nextflow clean -f

# Remove temporary Nextflow config directories (safe after pipeline completes):
# rm -rf ~/_nf_config_*

# -----------------------------------------------------------------------------------------
# COMMON ISSUES AND SOLUTIONS
# -----------------------------------------------------------------------------------------

# Issue: "No such file or directory" inside container
# Solution: Check path mappings printed above, ensure physical paths are bound

# Issue: "mapping values are not allowed here" in YAML
# Solution: Run this script - it fixes tabs in YAML automatically

# Issue: Pipeline won't resume after failure
# Solution: Don't run 'nextflow clean' - it deletes cache needed for -resume

# Issue: Out of disk space
# Solution: Clear old work directories: rm -rf ${WORK_DIR}/*
#           Warning: This prevents -resume for those runs!

# Issue: Singularity can't download images
# Solution: Check internet connection, verify cache_dir is writable