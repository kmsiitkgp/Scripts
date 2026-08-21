#!/bin/bash
#SBATCH --job-name=Biomodal_dmr
#SBATCH --mem=256G
#SBATCH --cpus-per-task=8
#SBATCH --time=4-00:00:00   # Requests exactly 4 days (96 hours)
#SBATCH --output=%x_%j.log

#===============================================================================
# 📄 Script   : 07_biomodal_dmr.sh
# 📌 Purpose  : Call differentially methylated regions (DMRs) across all
#               valid statistical combinations to determine the best approach
#               for this dataset. Results summarized for PI review.
#
# 🔬 Group designs (3):
#     long  — long duration patients only
#     short — short duration patients only
#     full  — all patients combined
#
# 📊 Statistical options per group:
#     Unpaired (4 variants):
#       unpaired_no_od_no_cov  — no overdispersion, no covariate
#       unpaired_od_no_cov     — overdispersion, no covariate
#       unpaired_no_od_cov     — no overdispersion, Patient_ID covariate
#       unpaired_od_cov        — overdispersion, Patient_ID covariate
#
#     Paired (2 variants — covariate always required for paired):
#       paired_no_od_cov       — no overdispersion, Patient_ID covariate
#       paired_od_cov          — overdispersion, Patient_ID covariate
#
# 🗂️  Total: 3 groups × 6 variants = 18 DMR runs
#            Each run produces 2 TSVs (hmc + mc) = 36 TSV files total
#
# 📁 Output structure under 04.dmr/:
#     merged_long.zarrz
#     merged_short.zarrz
#     merged_full.zarrz
#     dmr_long_unpaired_no_od_no_cov/
#     dmr_long_unpaired_od_no_cov/
#     dmr_long_unpaired_no_od_cov/
#     dmr_long_unpaired_od_cov/
#     dmr_long_paired_no_od_cov/
#     dmr_long_paired_od_cov/
#     dmr_short_unpaired_no_od_no_cov/
#     ... (same pattern for short and full)
#
# 🔁 Depends  : 05_biomodal_annotation.sh must have been run first
#
# 📝 Notes    :
#   - Unpaired sheets include C3 timepoints where available
#   - Paired sheets use C1 and C2 only (balanced 1:1 per patient)
#   - Full paired uses all patients regardless of duration
#   - One merged ZARR per group — sample sheet filters at DMR call time
#   - Coverage filters:
#       --min-coverage 2 : CpG site needs ≥2 reads pooled across group
#       --num-contexts 5 : region needs ≥5 CpGs passing coverage filter
#
# 📊 Directionality:
#   --condition-order Pre Post → Delta = Post - Pre
#   Positive delta_beta = Hypermethylation in Post
#   Negative delta_beta = Hypomethylation in Post
#===============================================================================

#===============================================================================
# 🚀 Environment
#===============================================================================

source "/home/kailasamms/miniconda3/etc/profile.d/conda.sh"
conda activate omics
echo "✅ Environment activated: omics"

#===============================================================================
# ⚙️ Paths
#===============================================================================

HOME_DIR="/home/kailasamms"
DUET_DIR="${HOME_DIR}/analysis/Biomodal/Kevin_Carotuximab/01.duet"
OUTPUT_DIR="${HOME_DIR}/analysis/Biomodal/Kevin_Carotuximab/04.dmr"
PROJ_DIR="${HOME_DIR}/scripts/nextflow/projects/Biomodal/Kevin_Carotuximab"
REF_DIR="${HOME_DIR}/resources/genomes"

# --- Annotation BED (gene-level) ---
BED_GENE="${REF_DIR}/GRCh38_v115_biomodal_annotation_gene.bed.gz"

# --- Sample sheets ---
# Unpaired: includes C3 timepoints where available
# Paired  : C1 and C2 only, balanced 1:1 per patient
SHEET_LONG_UNPAIRED="${PROJ_DIR}/Kevin_Carotuximab_metadata_long.csv"
SHEET_LONG_PAIRED="${PROJ_DIR}/Kevin_Carotuximab_metadata_long_paired.csv"
SHEET_SHORT_UNPAIRED="${PROJ_DIR}/Kevin_Carotuximab_metadata_short.csv"
SHEET_SHORT_PAIRED="${PROJ_DIR}/Kevin_Carotuximab_metadata_short_paired.csv"
SHEET_FULL_UNPAIRED="${PROJ_DIR}/Kevin_Carotuximab_metadata_full.csv"
SHEET_FULL_PAIRED="${PROJ_DIR}/Kevin_Carotuximab_metadata_full_paired.csv"

# --- Merged ZARR paths (one per group, shared across all variants) ---
MERGED_LONG="${OUTPUT_DIR}/merged_long.zarrz"
MERGED_SHORT="${OUTPUT_DIR}/merged_short.zarrz"
MERGED_FULL="${OUTPUT_DIR}/merged_full.zarrz"

# --- Column names ---
CONDITION_COLUMN="Condition"
COVARIATE_COLUMN="Patient_ID"
CONDITION_ORDER=("Pre" "Post")

# --- Coverage filters ---
MIN_COVERAGE=2
NUM_CONTEXTS=5

#===============================================================================
# 📋 Job Info
#===============================================================================

echo "=========================================================="
printf "Start date              : %s\n" "$(date)"
printf "Job name                : %s\n" "${SLURM_JOB_NAME}"
printf "Job ID                  : %s\n" "${SLURM_JOB_ID}"
printf "Cores                   : %s\n" "${SLURM_CPUS_PER_TASK}"
echo "=========================================================="
printf "DUET dir                : %s\n" "${DUET_DIR}"
printf "Output dir              : %s\n" "${OUTPUT_DIR}"
printf "BED gene                : %s\n" "${BED_GENE}"
printf "Sheet Long  (unpaired)  : %s\n" "${SHEET_LONG_UNPAIRED}"
printf "Sheet Long  (paired)    : %s\n" "${SHEET_LONG_PAIRED}"
printf "Sheet Short (unpaired)  : %s\n" "${SHEET_SHORT_UNPAIRED}"
printf "Sheet Short (paired)    : %s\n" "${SHEET_SHORT_PAIRED}"
printf "Sheet Full  (unpaired)  : %s\n" "${SHEET_FULL_UNPAIRED}"
printf "Sheet Full  (paired)    : %s\n" "${SHEET_FULL_PAIRED}"
printf "Condition column        : %s\n" "${CONDITION_COLUMN}"
printf "Covariate column        : %s\n" "${COVARIATE_COLUMN}"
printf "Condition order         : %s\n" "${CONDITION_ORDER[*]}"
printf "Min coverage            : %s\n" "${MIN_COVERAGE}"
printf "Num contexts            : %s\n" "${NUM_CONTEXTS}"
echo "=========================================================="

#===============================================================================
# 🔍 Pre-flight checks
#===============================================================================

echo "Checking required files..."
[[ -f "${BED_GENE}" ]]               || { echo "❌ Gene BED not found: ${BED_GENE}"; exit 1; }
[[ -f "${SHEET_LONG_UNPAIRED}" ]]    || { echo "❌ Not found: ${SHEET_LONG_UNPAIRED}";   exit 1; }
[[ -f "${SHEET_LONG_PAIRED}" ]]      || { echo "❌ Not found: ${SHEET_LONG_PAIRED}";     exit 1; }
[[ -f "${SHEET_SHORT_UNPAIRED}" ]]   || { echo "❌ Not found: ${SHEET_SHORT_UNPAIRED}";  exit 1; }
[[ -f "${SHEET_SHORT_PAIRED}" ]]     || { echo "❌ Not found: ${SHEET_SHORT_PAIRED}";    exit 1; }
[[ -f "${SHEET_FULL_UNPAIRED}" ]]    || { echo "❌ Not found: ${SHEET_FULL_UNPAIRED}";   exit 1; }
[[ -f "${SHEET_FULL_PAIRED}" ]]      || { echo "❌ Not found: ${SHEET_FULL_PAIRED}";     exit 1; }
echo "✅ All required files found."

# --- Create all 18 output directories ---
for GROUP in long short full; do
    for VARIANT in \
        unpaired_no_od_no_cov \
        unpaired_od_no_cov \
        unpaired_no_od_cov \
        unpaired_od_cov \
        paired_no_od_cov \
        paired_od_cov; do
        mkdir -p "${OUTPUT_DIR}/dmr_${GROUP}_${VARIANT}"
    done
done
echo "✅ All 18 output directories created."

#===============================================================================
# 🔎 Helper: build ZARR array filtered to samples in a given sheet
#===============================================================================

build_zarr_array() {
    local sheet="$1"
    local -n result_array="$2"

    local sample_ids
    sample_ids=$(cut -d ',' -f 1 "${sheet}" | tail -n +2)

    result_array=()
    shopt -s nullglob
    for zarr in "${DUET_DIR}"/duet-1.5.0_*/sample_outputs/zarr_store/CG/*.CG.zarrz; do
        local fname
        fname=$(basename "${zarr}")
        local sample_id="${fname%%.*}"
        if echo "${sample_ids}" | grep -qx "${sample_id}"; then
            result_array+=("${zarr}")
        fi
    done
    shopt -u nullglob
}

#===============================================================================
# 🔎 Helper: run a single DMR call
#
# Usage: run_dmr <zarr_path> <sheet> <output_dir> <label> <use_od> <use_cov>
#   use_od  : "yes" or "no"
#   use_cov : "yes" or "no"
#===============================================================================

run_dmr() {
    local zarr_path="$1"
    local sheet="$2"
    local out_dir="$3"
    local label="$4"
    local use_od="$5"
    local use_cov="$6"

    # --- RERUN SAFETY CHECK ---
    # Check if any .bed files already exist in a DMR subfolder.
    # If yes, we skip this entire run block.
    shopt -s nullglob
    local existing_beds=("${out_dir}"/DMR_*/*.bed)
    shopt -u nullglob

    if [ ${#existing_beds[@]} -gt 0 ]; then
        echo "⏩ Skipping DMR call: ${label} (BED files already exist in ${out_dir}!)"
        return 0
    fi

    echo "--- Calling DMRs: ${label}..."

    # Build command as array for clean flag handling
    local cmd=(
        modality dmr call
        --zarr-path  "${zarr_path}"
        --sample-sheet "${sheet}"
        --condition-order "${CONDITION_ORDER[@]}"
        --condition-array-name "${CONDITION_COLUMN}"
        --methylation-contexts mc hmc
        --bedfile "${BED_GENE}"
        --min-coverage "${MIN_COVERAGE}"
        --num-contexts "${NUM_CONTEXTS}"
        --output-dir "${out_dir}"
    )

    [[ "${use_od}"  == "yes" ]] && cmd+=(--overdispersion)
    [[ "${use_cov}" == "yes" ]] && cmd+=(--covariates "${COVARIATE_COLUMN}")

    "${cmd[@]}"

    if [ $? -eq 0 ]; then
        echo "✅ DMR call completed: ${label}"
    else
        echo "❌ DMR call failed: ${label}"
        exit 1
    fi
}

#===============================================================================
# 🔬 STEP 1 — Merge ZARRs (one per group)
#
# Uses the unpaired sheet (largest sample set) to collect all ZARRs per group.
# The sample sheet at DMR call time filters which samples are actually used.
#===============================================================================

echo "=========================================================="
echo "🔬 STEP 1 — Merging ZARRs"
echo "=========================================================="

for GROUP in long short full; do

    echo "--- Merging ${GROUP} ZARRs..."

    case "${GROUP}" in
        long)  SHEET_FOR_MERGE="${SHEET_LONG_UNPAIRED}";  MERGED="${MERGED_LONG}"  ;;
        short) SHEET_FOR_MERGE="${SHEET_SHORT_UNPAIRED}"; MERGED="${MERGED_SHORT}" ;;
        full)  SHEET_FOR_MERGE="${SHEET_FULL_UNPAIRED}";  MERGED="${MERGED_FULL}"  ;;
    esac

    # 1. Check early (OUTSIDE the function)
    if [ -d "${MERGED}" ]; then
        echo "⏩ Skipping ZARR merge for ${GROUP} — ${MERGED} already exists!"
        continue # Instantly jumps to the next group
    fi

    # 2. Only run the helper function if the file is truly missing
    echo "--- Merging ${GROUP} ZARRs..."
    build_zarr_array "${SHEET_FOR_MERGE}" ZARR_ARRAY
    echo "Found ${#ZARR_ARRAY[@]} ZARR files for ${GROUP}:"
    for z in "${ZARR_ARRAY[@]}"; do echo "  → $(basename "${z}")"; done

    if [ ${#ZARR_ARRAY[@]} -eq 0 ]; then
        echo "❌ No ZARR files found for ${GROUP}. Aborting."
        exit 1
    fi

    modality zarr-utils join \
        --zarr-path "${ZARR_ARRAY[@]}" \
        --output-path "${MERGED}"

    if [ $? -eq 0 ]; then
        echo "✅ ${GROUP} ZARR merge completed: ${MERGED}"
    else
        echo "❌ ${GROUP} ZARR merge failed."
        exit 1
    fi

done

#===============================================================================
# 🔬 STEP 2 — DMR calls: LONG (6 variants)
#===============================================================================

echo "=========================================================="
echo "🔬 STEP 2 — DMR calls: LONG"
echo "=========================================================="

# Long unpaired — no OD, no covariate
run_dmr "${MERGED_LONG}" "${SHEET_LONG_UNPAIRED}" \
    "${OUTPUT_DIR}/dmr_long_unpaired_no_od_no_cov" \
    "long_unpaired_no_od_no_cov" "no" "no"

# Long unpaired — OD, no covariate
run_dmr "${MERGED_LONG}" "${SHEET_LONG_UNPAIRED}" \
    "${OUTPUT_DIR}/dmr_long_unpaired_od_no_cov" \
    "long_unpaired_od_no_cov" "yes" "no"

# Long unpaired — no OD, covariate
run_dmr "${MERGED_LONG}" "${SHEET_LONG_UNPAIRED}" \
    "${OUTPUT_DIR}/dmr_long_unpaired_no_od_cov" \
    "long_unpaired_no_od_cov" "no" "yes"

# Long unpaired — OD, covariate
run_dmr "${MERGED_LONG}" "${SHEET_LONG_UNPAIRED}" \
    "${OUTPUT_DIR}/dmr_long_unpaired_od_cov" \
    "long_unpaired_od_cov" "yes" "yes"

# Long paired — no OD, covariate
run_dmr "${MERGED_LONG}" "${SHEET_LONG_PAIRED}" \
    "${OUTPUT_DIR}/dmr_long_paired_no_od_cov" \
    "long_paired_no_od_cov" "no" "yes"

# Long paired — OD, covariate
run_dmr "${MERGED_LONG}" "${SHEET_LONG_PAIRED}" \
    "${OUTPUT_DIR}/dmr_long_paired_od_cov" \
    "long_paired_od_cov" "yes" "yes"

#===============================================================================
# 🔬 STEP 3 — DMR calls: SHORT (6 variants)
#===============================================================================

echo "=========================================================="
echo "🔬 STEP 3 — DMR calls: SHORT"
echo "=========================================================="

# Short unpaired — no OD, no covariate
run_dmr "${MERGED_SHORT}" "${SHEET_SHORT_UNPAIRED}" \
    "${OUTPUT_DIR}/dmr_short_unpaired_no_od_no_cov" \
    "short_unpaired_no_od_no_cov" "no" "no"

# Short unpaired — OD, no covariate
run_dmr "${MERGED_SHORT}" "${SHEET_SHORT_UNPAIRED}" \
    "${OUTPUT_DIR}/dmr_short_unpaired_od_no_cov" \
    "short_unpaired_od_no_cov" "yes" "no"

# Short unpaired — no OD, covariate
run_dmr "${MERGED_SHORT}" "${SHEET_SHORT_UNPAIRED}" \
    "${OUTPUT_DIR}/dmr_short_unpaired_no_od_cov" \
    "short_unpaired_no_od_cov" "no" "yes"

# Short unpaired — OD, covariate
run_dmr "${MERGED_SHORT}" "${SHEET_SHORT_UNPAIRED}" \
    "${OUTPUT_DIR}/dmr_short_unpaired_od_cov" \
    "short_unpaired_od_cov" "yes" "yes"

# Short paired — no OD, covariate
run_dmr "${MERGED_SHORT}" "${SHEET_SHORT_PAIRED}" \
    "${OUTPUT_DIR}/dmr_short_paired_no_od_cov" \
    "short_paired_no_od_cov" "no" "yes"

# Short paired — OD, covariate
run_dmr "${MERGED_SHORT}" "${SHEET_SHORT_PAIRED}" \
    "${OUTPUT_DIR}/dmr_short_paired_od_cov" \
    "short_paired_od_cov" "yes" "yes"

#===============================================================================
# 🔬 STEP 4 — DMR calls: FULL (6 variants)
#===============================================================================

echo "=========================================================="
echo "🔬 STEP 4 — DMR calls: FULL"
echo "=========================================================="

# Full unpaired — no OD, no covariate
run_dmr "${MERGED_FULL}" "${SHEET_FULL_UNPAIRED}" \
    "${OUTPUT_DIR}/dmr_full_unpaired_no_od_no_cov" \
    "full_unpaired_no_od_no_cov" "no" "no"

# Full unpaired — OD, no covariate
run_dmr "${MERGED_FULL}" "${SHEET_FULL_UNPAIRED}" \
    "${OUTPUT_DIR}/dmr_full_unpaired_od_no_cov" \
    "full_unpaired_od_no_cov" "yes" "no"

# Full unpaired — no OD, covariate
run_dmr "${MERGED_FULL}" "${SHEET_FULL_UNPAIRED}" \
    "${OUTPUT_DIR}/dmr_full_unpaired_no_od_cov" \
    "full_unpaired_no_od_cov" "no" "yes"

# Full unpaired — OD, covariate
run_dmr "${MERGED_FULL}" "${SHEET_FULL_UNPAIRED}" \
    "${OUTPUT_DIR}/dmr_full_unpaired_od_cov" \
    "full_unpaired_od_cov" "yes" "yes"

# Full paired — no OD, covariate
run_dmr "${MERGED_FULL}" "${SHEET_FULL_PAIRED}" \
    "${OUTPUT_DIR}/dmr_full_paired_no_od_cov" \
    "full_paired_no_od_cov" "no" "yes"

# Full paired — OD, covariate
run_dmr "${MERGED_FULL}" "${SHEET_FULL_PAIRED}" \
    "${OUTPUT_DIR}/dmr_full_paired_od_cov" \
    "full_paired_od_cov" "yes" "yes"

#===============================================================================
# 📦 Convert DMR BED outputs to TSV
#
# Produces flat TSV files named: <label>_<DMR_dir>_<original_basename>.tsv
# so they are timestamp-agnostic and easy to load in R/Python.
#===============================================================================

echo "=========================================================="
echo "📦 Converting DMR BED outputs to TSV..."
echo "=========================================================="

for GROUP in long short full; do
    for VARIANT in \
        unpaired_no_od_no_cov \
        unpaired_od_no_cov \
        unpaired_no_od_cov \
        unpaired_od_cov \
        paired_no_od_cov \
        paired_od_cov; do

        LABEL="${GROUP}_${VARIANT}"
        DIR="${OUTPUT_DIR}/dmr_${LABEL}"

        shopt -s nullglob
        BED_ARRAY=("${DIR}"/DMR_*/*.bed)
        shopt -u nullglob

        if [ ${#BED_ARRAY[@]} -eq 0 ]; then
            echo "⚠️  No BED files in ${DIR} — skipping."
            continue
        fi

        for BED in "${BED_ARRAY[@]}"; do
            base=$(basename "${BED%.bed}")
            parent=$(basename "$(dirname "${BED}")")
            TSV="${DIR}/${LABEL}_${parent}_${base}.tsv"
            cp "${BED}" "${TSV}"
            echo "✅ ${TSV}"
        done

    done
done

#===============================================================================
# 📋 Summary
#===============================================================================

echo ""
echo "=========================================================="
printf "End date : %s\n" "$(date)"
echo "=========================================================="
echo "✅ All 18 DMR analyses completed."
echo ""
echo "📁 Outputs saved to: ${OUTPUT_DIR}"
echo ""
echo "  18 directories, each containing:"
echo "    *.tsv  — flat DMR results for R/Python"
echo "    *.bed  — raw modality output"
echo "    *.ini  — plot config"
echo "    *.json — run metadata"
echo ""
echo "  Naming convention:"
echo "    dmr_{group}_{paired/unpaired}_{od/no_od}_{cov/no_cov}/"
echo ""
echo "  Groups   : long | short | full"
echo "  Design   : unpaired (4 variants) | paired (2 variants)"
echo "  Total    : 18 runs × 2 marks (hmc + mc) = 36 TSV files"
echo "=========================================================="
