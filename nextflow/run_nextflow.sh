#!/bin/bash

#SBATCH --job-name=NextFlow_Launcher
#SBATCH --chdir=./
#SBATCH --export=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G                        # Gives the "Master" plenty of room
#SBATCH --time=96:00:00
#SBATCH --output=nf_launcher_%j.log

# # Use in SGE
# #$ -N NextFlow_Launcher        # Job Name
# #$ -cwd                        # Run in current directory
# #$ -V                          # Export all environment variables
# #$ -j y                        # Merge standard output and error logs
# #$ -l mem_free=6G              # Need atleast 6GB to start job
# #$ -l h_vmem=6G                # Kill job if it exceeds 6G
# #$ -l h_rt=96:00:00            # Wall clock limit (equivalent to SLURM's --time)
# #$ -pe smp 1                   # Request 1 CPU slot (equivalent to SLURM's --cpus-per-task)
# #$ -o nf_launcher_$JOB_ID.log  # Name of the output log

# 1. Cleaner logging for HPC files (Important for sbatch logs!)
export NXF_ANSI_LOG=false

# 2. Safety limits for Java (Keep -Xmx slightly below the --mem above)
export NXF_OPTS="-Xms1g -Xmx5g"
# export NXF_OPTS="-Xms1g -Xmx4g -XX:MaxMetaspaceSize=512m -XX:ActiveProcessorCount=2 -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2"

# =========================================================================================
# NEXTFLOW PIPELINE LAUNCHER
# =========================================================================================
# This script:
#   1. Reads paths from config.yaml
#   2. Auto-detects physical paths (resolving symlinks)
#   3. Builds Singularity bind mounts dynamically
#   4. Launches the Nextflow pipeline with proper configuration
#
# Usage: cd to directory containing yaml file
#    bash run_nextflow.sh
#
# Prerequisites:
#    - config.yaml configured with correct paths
#    - Nextflow and Singularity modules available
#    - Access to HPC cluster (SLURM by default; SGE lines commented below)
#
# Resume a failed run:
#    Just run the script again - Nextflow automatically resumes from last checkpoint
# =========================================================================================

# =========================================================================================
# HELPER FUNCTION: get_yaml_val()
# =========================================================================================
# This function extracts a specific value from the YAML file.
#
# Breakdown of the command chain:
# 1. grep "^[[:space:]]*[^#]" : Filters out any lines that start with a '#' (comments).
#                               Matches lines even if they have spaces before the text.
# 2. grep "${1}[[:space:]]*:" : Searches for the specific key (e.g., read_dir) followed by a colon.
# 3. head -n 1                : Takes only the first match (prevents errors if key is duplicated).
# 4. cut -d':' -f2-           : Splits the line at the first colon and takes everything after it.
# 5. sed 's/#.*//'            : Removes any inline comments (everything from # to end of line)
# 6. xargs                    : Automatically trims leading/trailing whitespace.
# 7. tr -d '"' | tr -d "'"    : Removes both double and single quotes from paths.
# =========================================================================================
get_yaml_val() {
    grep "^[[:space:]]*[^#]" "${YAML}" | grep "${1}[[:space:]]*:" | head -n 1 | cut -d':' -f2- | sed 's/#.*//' | xargs | tr -d '"' | tr -d "'"
}

# =========================================================================================
# HELPER FUNCTION: check_path()
# =========================================================================================
# This function checks if a path exists and return a status label
check_path() {
    if [ -z "$1" ]; then
        echo "Optional (Not Defined)"
    elif [ -e "$1" ]; then
        echo "✔ FOUND"
    else
        echo "❌ MISSING"
    fi
}

# -----------------------------------------------------------------------------------------
# 1. AUTO-DETECT CONFIGURATION
# -----------------------------------------------------------------------------------------
# Looks for the first .yaml file in the current directory
YAML=$(ls *.yaml 2>/dev/null | head -n 1)

if [ -z "${YAML}" ]; then
    echo "❌ ERROR: No YAML configuration file (.yaml) found in this directory."
    exit 1
fi

# RESOLVE FULL PATH IMMEDIATELY
YAML_FILE=$(readlink -f "${YAML}")

# -----------------------------------------------------------------------------------------
# 2. CLEAN YAML CONFIG
# -----------------------------------------------------------------------------------------
# Sanitizes the YAML for Linux/HPC environments (Windows editors are the usual culprit)

# Define your paths
MODULE_DIR="/home/kailasamms/scripts/nextflow/modules"

# # 1. CLEAN R MODULES & NEXTFLOW FILES [ ONLY ONCE is enough ]
# # Removes Windows Line Endings (^M), BOM, and sets Executable permissions
# for FILE in "${MODULE_DIR}"/*.{R,nf}; do
    # if [ -f "$FILE" ]; then
        # # Remove UTF-8 BOM
        # sed -i '1s/^\xEF\xBB\xBF//' "$FILE"
        # # Convert CRLF (Windows) to LF (Linux) - Fixes Exit 127
        # sed -i 's/\r//g' "$FILE"
        # # Ensure scripts are executable
        # chmod +x "$FILE"
        # echo "✅ Sanitized: $(basename "$FILE")"
    # fi
# done

# 2. CLEAN YAML CONFIG (Specific YAML rules)
if [ -f "$YAML" ]; then
    # Remove UTF-8 BOM (Byte Order Mark) from Windows editors
    sed -i '1s/^\xEF\xBB\xBF//' "$YAML"
    # Convert Tabs to Spaces (YAML's biggest enemy)
    sed -i 's/\t/    /g' "$YAML"
    # Strip trailing whitespace
    sed -i 's/[[:space:]]*$//' "$YAML"
    # Final CRLF check for YAML
    sed -i 's/\r//g' "$YAML"
    echo "✅ Configuration file $(basename "$YAML") cleaned."
fi

# -----------------------------------------------------------------------------------------
# 3. EXTRACT PARAMETERS
# -----------------------------------------------------------------------------------------
EXPERIMENT=$(get_yaml_val "expt")
PROJECT=$(get_yaml_val "project")
READ_DIR=$(get_yaml_val "read_dir")


BASE_DIR=$(get_yaml_val "base_dir")
WORK_DIR=$(get_yaml_val "work_dir")
CACHE_DIR=$(get_yaml_val "cache_dir")
REF_DIR=$(get_yaml_val "ref_dir")

# -----------------------------------------------------------------------------------------
# 4. VERIFICATION BLOCK (Launcher-only)
# -----------------------------------------------------------------------------------------
echo "--------------------------------------------------------------------------------"
echo "        LAUNCHER VERIFICATION"
echo "--------------------------------------------------------------------------------"
printf "%-22s : %-15s : %s\n" "PARAMETER" "STATUS" "VALUE"
echo "--------------------------------------------------------------------------------"
printf "%-22s : %-15s : %s\n" "YAML File"       "$(check_path "$YAML_FILE")" "$YAML_FILE"
printf "%-22s : %-15s : %s\n" "Temporary Files" "$(check_path "$WORK_DIR")" "$WORK_DIR"
echo "--------------------------------------------------------------------------------"

# -----------------------------------------------------------------------------------------
# 5. RESOLVE PHYSICAL PATHS (For Singularity Binds)
# -----------------------------------------------------------------------------------------
# Singularity needs physical paths (not symlinks) for bind mounts
P_READ=$(readlink -f "$READ_DIR")
P_BASE=$(readlink -f "$BASE_DIR")
P_WORK=$(readlink -f "$WORK_DIR")
P_CACHE=$(readlink -f "$CACHE_DIR")
P_REF=$(readlink -f "$REF_DIR")

# -----------------------------------------------------------------------------------------
# 6. DISPLAY PATH MAPPINGS
# -----------------------------------------------------------------------------------------
# Show how paths map between host (physical) and container (logical)
# Useful for debugging "file not found" issues

echo "---------------------------------------------------------------------------------------------------------------"
printf "%-16s : %-60s -> %s\n" "Mapping:" "PHYSICAL (Host)" "LOGICAL (Container)"
echo "---------------------------------------------------------------------------------------------------------------"
printf "%-16s : %-60s -> %s\n" "Raw Reads"     "$P_READ" "$READ_DIR"
printf "%-16s : %-60s -> %s\n" "Project Base"   "$P_BASE" "$BASE_DIR"
printf "%-16s : %-60s -> %s\n" "Temp Work Dir"  "$P_WORK" "$WORK_DIR"
printf "%-16s : %-60s -> %s\n" "Image Cache"    "$P_CACHE" "$CACHE_DIR"
printf "%-16s : %-60s -> %s\n" "Ref Genomes"    "$P_REF"  "$REF_DIR"
echo "---------------------------------------------------------------------------------------------------------------"

# -----------------------------------------------------------------------------------------
# 7. BUILD SINGULARITY BIND MOUNTS
# -----------------------------------------------------------------------------------------
# Collect all physical paths and identify unique root directories (e.g., /hpc, /scratch)
if [ "$EXPERIMENT" == "TCGA" ]; then
    ALL_PATHS="$P_BASE $P_WORK $P_CACHE /common/bhowmicknlab/TCGA"
else
    ALL_PATHS="$P_READ $P_BASE $P_WORK $P_CACHE $P_REF"
fi
UNIQUE_ROOTS=$(for p in $ALL_PATHS; do echo "$p" | cut -d'/' -f1-2; done | sort -u | grep '^/')

BIND_FLAGS=""
for path in $UNIQUE_ROOTS; do
    if [ -d "$path" ]; then
        BIND_FLAGS+="--bind $path "
    fi
done

echo "Singularity bind flags: $BIND_FLAGS"

# -----------------------------------------------------------------------------------------
# 8. PIPELINE EXECUTION
# -----------------------------------------------------------------------------------------
# -params-file    : Injects YAML settings (FASTQ paths, species, etc.)
# -name           : Labels the run in 'nextflow log' and SLURM/SGE (highly recommended)
# -work-dir       : Explicitly sets the heavy data directory. We define this for SAFETY
#                   to ensure task data stays in scratch even if the config changes.
# -profile        : Sets executor (slurm below; sge commented alongside)
# -resume         : Uses cached results to skip successfully completed steps
# --dynamic_binds : Passes calculated physical paths to Singularity containers

# Define a project-specific workspace within the scratch directory.
mkdir -p "${WORK_DIR}/${PROJECT}"

# Move into the project workspace.
cd "${WORK_DIR}/${PROJECT}"

# Create a unique run name using a timestamp.
RUN_NAME="${PROJECT}_$(date +%Y-%m-%d_%H-%M-%S)"

# We use absolute path to main.nf so the script works from anywhere
NF_MAIN="/home/kailasamms/scripts/nextflow/main.nf"

# Load Nextflow and Singularity from HPC environment modules
module load nextflow
module load apptainer

# Launch Nextflow
nextflow \
    -log "${WORK_DIR}/${PROJECT}.nextflow.log" \
    run "$NF_MAIN" \
    -params-file "${YAML_FILE}" \
    -name "${RUN_NAME}" \
    -work-dir "${WORK_DIR}/${PROJECT}" \
    -profile slurm \
    -resume \
    --dynamic_binds "$BIND_FLAGS" \
    --stop_after "END"
    # -profile sge \    # Use in SGE environments instead of -profile slurm

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
# COMMON ISSUES AND SOLUTIONS
# -----------------------------------------------------------------------------------------
# "No such file or directory" inside container
#   → Check path mappings printed above, ensure physical paths are bound
#
# "mapping values are not allowed here" in YAML
#   → Run this script - it fixes tabs in YAML automatically
#
# Pipeline won't resume after failure
#   → Don't run 'nextflow clean' - it deletes cache needed for -resume
#
# Out of disk space
#   → Clear old work directories: rm -rf ${WORK_DIR}/*
#   → Warning: This prevents -resume for those runs!
#
# Singularity can't download images
#   → Check internet connection, verify cache_dir is writable

# -----------------------------------------------------------------------------------------
# SESSION MANAGEMENT (tmux/screen) — optional, for interactive monitoring
# -----------------------------------------------------------------------------------------
# tail -f "${WORK_DIR}/${PROJECT}.nextflow.log"   # See logs live

# Start session:
# tmux new -s "${PROJECT}"
# screen -S "${PROJECT}"

# Detach session:
# Press Ctrl+B first, release and press D (for tmux)
# Press Ctrl+A first, release and press D (for screen)

# Reattach session
# tmux attach -t "${PROJECT}"
# screen -r "${PROJECT}"

# -----------------------------------------------------------------------------------------
# DEBUGGING TIPS
# -----------------------------------------------------------------------------------------

# Check for UTF-8 BOM in files:
#   head -c 3 "${NF_DIR}/main.nf" | od -c
#   Output "0000000 357 273 277 # ! /" → BOM present, needs removal
#   Output "0000000 # ! /"             → No BOM, file is clean

# View which work directories correspond to which processes:
#   nextflow log <run_name> -f name,workdir,status,duration

# List all runs:            nextflow log
# Clean up run metadata:     nextflow clean -f
# Remove temp NF configs:    rm -rf ~/_nf_config_*

# -----------------------------------------------------------------------------------------
# DELETE WORK DIRS FROM FAILED RUNS (RUN AFTER SUCCESSFUL COMPLETION OF PROJECT)
# -----------------------------------------------------------------------------------------

# Dry run first — see what would be deleted
#nextflow log -q -f status | grep -v 'OK' | cut -f1 | xargs -n 1 -I {} echo nextflow clean {} -f

# If looks good, actually delete
#nextflow log -q -f status | grep -v 'OK' | cut -f1 | xargs -n 1 -I {} nextflow clean {} -f

# -----------------------------------------------------------------------------------------
# CLEAN THE PIPELINE CODE (Only if needed)
# -----------------------------------------------------------------------------------------
# Tip: run manually after editing files on Windows

# #!/bin/bash

# # Enable nullglob so empty directories don't cause errors
# shopt -s nullglob

# NF_DIR="/home/kailasamms/scripts/nextflow"

# # BOM appears as: 0xEF 0xBB 0xBF at file start (causes "Invalid character at start of file")
# # Strip trailing whitespace (prevents unnecessary git diffs)
# # Convert tabs to 4 spaces (consistent indentation)

# find "${NF_DIR}" \
    # -type d \( \
        # -path "${NF_DIR}/resources" -o \
        # -path "${NF_DIR}/md" -o \
        # -path "${NF_DIR}/modules/junk" -o \
        # -name ".nextflow" \
    # \) -prune -o \
    # -type f \( -name "*.nf" -o -name "*.config" -o -name "*.sh" \
            # -o -name "*.txt" -o -name "*.R" -o -name "*.yaml" \) -print |
# while IFS= read -r file; do
    # sed -i '1s/^\xEF\xBB\xBF//' "$file"      # Remove UTF-8 BOM
    # sed -i 's/\r$//' "$file"                 # CRLF → LF
    # sed -i 's/[[:space:]]*$//' "$file"       # Strip trailing whitespace
    # sed -i 's/\t/    /g' "$file"             # Tabs → 4 spaces
    # echo "✅ Sanitized: $(basename "$file")"
# done

# # Executable bit only for scripts meant to be run directly
# chmod +x "${NF_DIR}"/*.sh 2>/dev/null

# echo "All files are now Unix-friendly (UTF-8, LF, No BOM)."

# # Disable nullglob
# shopt -u nullglob

# # 1. Convert to UTF-8 (and remove ANSI/ISO encoding)
# # iconv -c ignores invalid characters to prevent crashing
# if command -v iconv >/dev/null; then
        # iconv -f WINDOWS-1252 -t UTF-8 "$file" > "${file}.tmp" 2>/dev/null && mv "${file}.tmp" "$file"
# fi