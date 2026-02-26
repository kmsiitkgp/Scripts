#!/bin/bash

FINAL_OUT="Salmon_Master_Summary.txt"

# 1. Header (Expanded with QC metrics)
HEADER="Directory_ID\tSample\tPurity (%)"
HEADER="$HEADER\tTotal (%)\tUnmapped (%)\tUnique (%)"
HEADER="$HEADER\t\tTotal (#)\tUnmapped (#)\tUnique (#)"
HEADER="$HEADER\t\tFrag_Mean\tLibType\tGC_Bias\tSeq_Bias"
HEADER="$HEADER\tDecoy (%)\tDovetail (%)\tTPM Sum"

echo -e "$HEADER" > "$FINAL_OUT"

# Helper Functions
add_commas() { echo "$1" | sed ':a;s/\B[0-9]\{3\}\(>\|$\)/,&/;ta'; }
perc() { awk -v part=$1 -v total=$2 'BEGIN{if(total>0) printf "%.2f", (part/total)*100; else print 0}'; }
parse_json() { grep "$1" "$2" | tr -d '", ' | awk -F':' '{print $2}'; }
clean_bool() { [[ "$1" == "true" ]] && echo "YES" || echo "no"; }

CWD=$(pwd)

# --- PRE-SCAN FOR TOTAL READS ---
declare -A FULL_TOTAL_READS
echo "Scanning for baselines..."

ALL_JSONS=$(find "$CWD" -type f -path "*full*/03.Salmon/*/aux_info/meta_info.json")
for JSON in $ALL_JSONS; do
    SNAME=$(basename $(dirname $(dirname "$JSON")))
    T_VAL=$(parse_json "num_processed" "$JSON")
    FULL_TOTAL_READS["$SNAME"]=$T_VAL
done

# --- PROCESS ALL DIRECTORIES ---
SALMON_DIRS=$(find "$CWD" -type d -name "03.Salmon" | sort -V)

echo "Generating Master Salmon Report (Expanded QC)..."

for S_DIR in $SALMON_DIRS; do
    DIR_ID=$(echo "$S_DIR" | sed "s|$CWD/||g")
    
    # Natural sort MT1, MT2, MT10
    SAMPLES=$(find "$S_DIR" -maxdepth 1 -mindepth 1 -type d | sort -V)

    for SAMP_DIR in $SAMPLES; do
        SAMPLE=$(basename "$SAMP_DIR")
        if [[ "$SAMPLE" == "quant_files" || "$SAMPLE" == "aux_info" || "$SAMPLE" == "logs" ]]; then continue; fi
        
        QFILE="$SAMP_DIR/quant.sf"
        META="$SAMP_DIR/aux_info/meta_info.json"

        if [[ -f "$QFILE" && -f "$META" ]]; then
            # Basic Meta
            LIB_TYPE=$(grep -A 1 "library_types" "$META" | grep -v "library_types" | tr -d '[]", ')
            GC_B=$(clean_bool $(parse_json "gc_bias_correct" "$META"))
            SEQ_B=$(clean_bool $(parse_json "seq_bias_correct" "$META"))
            
            TOTAL=$(parse_json "num_processed" "$META")
            MAPPED=$(parse_json "num_mapped" "$META")
            
            # QC Metrics
            DECOY=$(parse_json "num_decoy_fragments" "$META")
            DOVE=$(parse_json "num_dovetail_fragments" "$META")
            DECOY_P=$(perc $DECOY $TOTAL)
            DOVE_P=$(perc $DOVE $TOTAL)
            
            # Mapping Stats
            MAP_RATE=$(perc $MAPPED $TOTAL)
            UNMAPPED_NUM=$((TOTAL - MAPPED))
            UNMAPPED_PERC=$(perc $UNMAPPED_NUM $TOTAL)
            FRAG_LEN=$(parse_json "frag_length_mean" "$META" | awk '{printf "%.1f", $1}')
            
            # Purity
            FULL_BASE=${FULL_TOTAL_READS["$SAMPLE"]:-$TOTAL}
            PURITY_SCORE=$(perc $TOTAL $FULL_BASE)

            # Sum of Estimated Counts and TPM
            COUNTS_TPM=$(awk 'NR>1{c+=$5; t+=$4} END{printf "%.0f\t%.0f", c, t}' "$QFILE")
            EST_COUNTS=$(echo "$COUNTS_TPM" | cut -f1)
            TPM_SUM=$(echo "$COUNTS_TPM" | cut -f2)

            # --- CONSTRUCT ROW ---
            row="$DIR_ID\t$SAMPLE\t$PURITY_SCORE"
            row="$row\t100.00\t$UNMAPPED_PERC\t$MAP_RATE"
            row="$row\t\t$(add_commas $TOTAL)\t$(add_commas $UNMAPPED_NUM)\t$(add_commas $EST_COUNTS)"
            row="$row\t\t$FRAG_LEN\t$LIB_TYPE\t$GC_B\t$SEQ_B"
            row="$row\t$DECOY_P%\t$DOVE_P%\t$(add_commas $TPM_SUM)"

            echo -e "$row" >> "$FINAL_OUT"
        fi
    done

    # Visual spacer
    spacer=""
    for i in {1..17}; do spacer="$spacer\t"; done
    echo -e "$spacer" >> "$FINAL_OUT"
done

echo "Master Report generated: $FINAL_OUT"