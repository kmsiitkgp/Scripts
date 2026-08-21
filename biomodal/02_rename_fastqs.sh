#!/bin/bash -l

#===============================================================================
# 🔎 Load mapping file
#===============================================================================

# NOTE: mapping file is .txt without column names like below:
# SRR10162346_1.fastq.gz	CRPC1_S01_L001_R1_001.fastq.gz
# SRR10162346_2.fastq.gz	CRPC1_S01_L001_R2_001.fastq.gz
# SRR10162347_1.fastq.gz	CRPC2_S02_L001_R1_001.fastq.gz
# SRR10162347_2.fastq.gz	CRPC2_S02_L001_R2_001.fastq.gz
# NOTE: It MUST have Unix(LF) line endings

# Define the location of your mapping file
MAP_FILE="/home/kailasamms/scripts/nextflow/projects/Biomodal/Kevin_Carotuximab/Kevin_Carotuximab_fastq_map.txt"

# --- Check if map file exists ---
if [[ ! -f "$MAP_FILE" ]]; then
    echo "❌ Error: Mapping file not found at $MAP_FILE"
    exit 1
fi

# --- Check if map file has windows line endings ---
# The -q flag suppresses output. $'\r' is the standard Bash way to represent the Carriage Return character.
if grep -q $'\r' "$MAP_FILE"; then
    echo "❌ ERROR: Mapping file contains Windows line endings (^M)."
    echo "The presence of these invisible characters will cause file renaming to fail."
    echo "Please run: dos2unix $MAP_FILE" 
    exit 1
fi

#===============================================================================
# 🔄 Rename FASTQ files
#===============================================================================

FASTQ_DIR="/common/bhowmicknlab/Biomodal/Kevin_Carotuximab/input"
echo "Starting file renaming process in: $FASTQ_DIR"
echo "--------------------------------------------------------"
RENAMED_COUNT=0
FAILURE_COUNT=0

# Read the mapping file line by line
# IFS= reads lines into two variables, assuming names are separated by whitespace (space or tab)
while IFS=$' \t' read -r OLD_NAME NEW_NAME || [[ -n "${OLD_NAME:-}" ]]; do
    # Skip comments and empty lines
    if [[ "$OLD_NAME" =~ ^# || -z "${OLD_NAME:-}" ]]; then
        continue
    fi

    SOURCE_PATH="${FASTQ_DIR}/${OLD_NAME}"
    DEST_PATH="${FASTQ_DIR}/${NEW_NAME}"

    # --- Safety Check 1: Ensure the source file exists ---
    if [[ ! -f "$SOURCE_PATH" ]]; then
        echo "⚠️ WARNING: Source file not found: $SOURCE_PATH"
        FAILURE_COUNT=$((FAILURE_COUNT + 1))
        continue
    fi
    
    # --- Safety Check 2: Ensure the destination name doesn't already exist ---
    if [[ -f "$DEST_PATH" ]]; then
        echo "❌ ERROR: Destination file already exists: $DEST_PATH. Skipping rename for $OLD_NAME."
        FAILURE_COUNT=$((FAILURE_COUNT + 1))
        continue
    fi

    # --- Execute the rename (mv) command ---
    echo "Renaming: $OLD_NAME  -->  $NEW_NAME"
    mv "$SOURCE_PATH" "$DEST_PATH"
    
    if [ $? -eq 0 ]; then
        RENAMED_COUNT=$((RENAMED_COUNT + 1))
    else
        echo "❌ RENAME FAILED for $OLD_NAME"
        FAILURE_COUNT=$((FAILURE_COUNT + 1))
    fi

done < "$MAP_FILE"

echo "--------------------------------------------------------"
echo "✅ Renaming complete. Total files renamed: $RENAMED_COUNT"
echo "❌ Total failures/skips: $FAILURE_COUNT"
