#!/usr/bin/env bash

# --- 1. ARGUMENT HANDLING ---
SEARCH_DIR="${1:-.}"
FINAL_OUT="Salmon_Master_Summary.txt"

# --- 2. HEADER ---
HEADER="Directory_ID\tSample\tPurity (%)"
HEADER="$HEADER\tTotal (%)\tUnmapped (%)\tUnique (%)"
HEADER="$HEADER\t\tTotal (#)\tUnmapped (#)\tUnique (#)"
HEADER="$HEADER\t\tFrag_Mean\tLibType\tGC_Bias\tSeq_Bias"
HEADER="$HEADER\tDecoy (%)\tDovetail (%)\tTPM Sum"

printf "%b\n" "$HEADER" > "$FINAL_OUT"

# --- 3. HELPER FUNCTIONS ---

add_commas() {
    awk -v n="$1" 'BEGIN{
        if(n=="" || n==0){print 0; exit}
        printf "%'\''d\n", n
    }'
}

perc() {
    awk -v p="$1" -v t="$2" 'BEGIN{
        if(t>0) printf "%.2f", (p/t)*100;
        else printf "0.00"
    }'
}

clean_bool() {
    [[ "$1" == "true" ]] && echo "YES" || echo "no"
}

# --- 4. XENOGRAFT BASELINE PRE-SCAN ---
declare -A XENO_BASELINES

echo "Scanning for Xenograft baselines in: $SEARCH_DIR"

while IFS= read -r META; do
    if [[ "$META" =~ /Xenograft/03\.Salmon/([^/]+)/aux_info/meta_info\.json$ ]]; then
        SAMPLE="${BASH_REMATCH[1]}"
        TOTAL=$(awk -F: '/"num_processed"/ {gsub(/[", ]/,"",$2); print $2; exit}' "$META")
        [[ -n "$TOTAL" ]] && XENO_BASELINES["$SAMPLE"]=$TOTAL
        echo "Locked Baseline: $SAMPLE -> $TOTAL"
    fi
done < <(find -L "$SEARCH_DIR" -type f -path "*03.Salmon*/aux_info/meta_info.json")

# --- 5. MAIN REPORT GENERATION ---

while IFS= read -r S_DIR; do

    DIR_ID="${S_DIR#$SEARCH_DIR/}"
    DIR_ID="${DIR_ID#/}"

    while IFS= read -r SAMP_DIR; do

        SAMPLE=$(basename "$SAMP_DIR")
        QFILE="$SAMP_DIR/quant.sf"
        META="$SAMP_DIR/aux_info/meta_info.json"

        [[ ! -f "$QFILE" || ! -f "$META" ]] && continue

        # --- Parse JSON once (Forcing Space Splitting) ---
        read -r TOTAL MAPPED FRAG_MEAN GC_RAW SEQ_RAW DECOY DOVE LIB_TYPE <<< $(
            awk -F: '
                /"num_processed"/ {gsub(/[", ]/,"",$2); np=$2}
                /"num_mapped"/ {gsub(/[", ]/,"",$2); nm=$2}
                /"frag_length_mean"/ {gsub(/[", ]/,"",$2); fl=$2}
                /"gc_bias_correct"/ {gsub(/[", ]/,"",$2); gc=$2}
                /"seq_bias_correct"/ {gsub(/[", ]/,"",$2); sb=$2}
                /"num_decoy_fragments"/ {gsub(/[", ]/,"",$2); nd=$2}
                /"num_dovetail_fragments"/ {gsub(/[", ]/,"",$2); dv=$2}
                /"library_types"/ {getline; gsub(/[\[\]", ]/,""); lt=$0}
                END {print np, nm, fl, gc, sb, nd, dv, lt}
            ' "$META"
        )

        GC_B=$(clean_bool "$GC_RAW")
        SEQ_B=$(clean_bool "$SEQ_RAW")

        UNMAPPED_NUM=$((TOTAL - MAPPED))

        EST_COUNTS=$(awk 'NR>1{sum+=$5} END{printf "%.0f", sum}' "$QFILE")
        TPM_SUM=$(awk 'NR>1{sum+=$4} END{printf "%.0f", sum}' "$QFILE")

        # --- Purity Logic ---
        FINAL_BASE="${XENO_BASELINES[$SAMPLE]:-$TOTAL}"
        PURITY_SCORE=$(perc "$TOTAL" "$FINAL_BASE")

        TOTAL_PCT=$(perc "$TOTAL" "$TOTAL")
        UNMAP_PCT=$(perc "$UNMAPPED_NUM" "$TOTAL")
        MAP_PCT=$(perc "$MAPPED" "$TOTAL")
        DECOY_PCT=$(perc "$DECOY" "$TOTAL")
        DOVE_PCT=$(perc "$DOVE" "$TOTAL")

        # --- Output Row ---
        printf "%s\t%s\t%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t\t%.1f\t%s\t%s\t%s\t%s%%\t%s%%\t%s\n" \
            "$DIR_ID" \
            "$SAMPLE" \
            "$PURITY_SCORE" \
            "$TOTAL_PCT" \
            "$UNMAP_PCT" \
            "$MAP_PCT" \
            "$(add_commas "$TOTAL")" \
            "$(add_commas "$UNMAPPED_NUM")" \
            "$(add_commas "$EST_COUNTS")" \
            "$FRAG_MEAN" \
            "$LIB_TYPE" \
            "$GC_B" \
            "$SEQ_B" \
            "$DECOY_PCT" \
            "$DOVE_PCT" \
            "$(add_commas "$TPM_SUM")" \
            >> "$FINAL_OUT"

    done < <(find -L "$S_DIR" -maxdepth 1 -mindepth 1 -type d \
             ! -name "logs" ! -name "aux_info" ! -name "quant_files" | sort -V)

    printf "%b\n" "$(printf '\t%.0s' {1..17})" >> "$FINAL_OUT"

done < <(find -L "$SEARCH_DIR" -type d -name "03.Salmon" | sort -V)

echo "Done. Salmon Report saved to: $FINAL_OUT"