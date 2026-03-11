#!/usr/bin/env bash

# --- 1. ARGUMENT HANDLING ---
SEARCH_DIR="${1:-.}"
FINAL_OUT="STAR_Master_Summary.txt"

# --- 2. HEADER ---
HEADER="Directory_ID\tSample\tPurity (%)"
HEADER="$HEADER\tTotal (%)\tMulti (%)\tUnmapped (%)\tUnique (%)"
HEADER="$HEADER\t\tFeature U (%)\tAmbiguous U (%)\tNoFeature U (%)"
HEADER="$HEADER\t\tFeature F (%)\tAmbiguous F (%)\tNoFeature F (%)"
HEADER="$HEADER\t\tFeature R (%)\tAmbiguous R (%)\tNoFeature R (%)"
HEADER="$HEADER\t\tTotal (#)\tMulti (#)\tUnmapped (#)\tUnique (#)"
HEADER="$HEADER\t\tFeature U (#)\tAmbiguous U (#)\tNoFeature U (#)"
HEADER="$HEADER\t\tFeature F (#)\tAmbiguous F (#)\tNoFeature F (#)"
HEADER="$HEADER\t\tFeature R (#)\tAmbiguous R (#)\tNoFeature R (#)"

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

# --- 4. XENOGRAFT BASELINE PRE-SCAN (STRUCTURAL) ---

declare -A XENO_BASELINES

echo "Scanning for Xenograft baselines in: $SEARCH_DIR"

while IFS= read -r LOG; do
    if [[ "$LOG" =~ /Xenograft/.*alignment_stats/([^/]+)\.Log\.final\.out$ ]]; then
        SAMPLE="${BASH_REMATCH[1]}"

        TOTAL=$(awk -F'|' '
            /Number of input reads/ {
                gsub(/ /,"",$2);
                print $2+0;
                exit
            }
        ' "$LOG")

        if [[ -n "$TOTAL" ]]; then
            XENO_BASELINES["$SAMPLE"]=$TOTAL
            echo "Locked baseline: $SAMPLE -> $TOTAL"
        fi
    fi
done < <(find -L "$SEARCH_DIR" -type f -name "*.Log.final.out")

# --- 5. MAIN REPORT GENERATION ---

while IFS= read -r BDIR; do

    DIR_ID="${BDIR#$SEARCH_DIR/}"
    DIR_ID="${DIR_ID#/}"

    ALIGN_DIR="$BDIR/alignment_stats"
    GENE_DIR="$BDIR/gene_counts"

    [[ ! -d "$ALIGN_DIR" || ! -d "$GENE_DIR" ]] && continue

    while IFS= read -r GFILE; do

        SAMPLE=$(basename "$GFILE" .ReadsPerGene.out.tab)
        LOG="$ALIGN_DIR/${SAMPLE}.Log.final.out"

        [[ ! -f "$LOG" ]] && continue

        # --- Parse STAR Log Once ---
        read TOTAL UNI M1 M2 U1 U2 < <(
            awk -F'|' '
                /Number of input reads/ {gsub(/ /,"",$2); t=$2}
                /Uniquely mapped reads number/ {gsub(/ /,"",$2); u=$2}
                /Number of reads mapped to multiple loci/ {gsub(/ /,"",$2); m1=$2}
                /Number of reads mapped to too many loci/ {gsub(/ /,"",$2); m2=$2}
                /Number of reads unmapped: too many mismatches/ {gsub(/ /,"",$2); u1=$2}
                /Number of reads unmapped: too short/ {gsub(/ /,"",$2); u2=$2}
                END {print t+0, u+0, m1+0, m2+0, u1+0, u2+0}
            ' "$LOG"
        )

        MULTI=$((M1 + M2))
        UNMAP=$((U1 + U2))

        # --- Parse Gene Counts Once ---
        read GUN GFW GRV AMB_UN AMB_FW AMB_RV NOF_UN NOF_FW NOF_RV < <(
            awk '
                NR==3 {n1=$2; n2=$3; n3=$4}
                NR==4 {a1=$2; a2=$3; a3=$4}
                NR>4 {u1+=$2; u2+=$3; u3+=$4}
                END {
                    print u1+0, u2+0, u3+0,
                          a1+0, a2+0, a3+0,
                          n1+0, n2+0, n3+0
                }
            ' "$GFILE"
        )

        # --- Purity Calculation ---
        FINAL_BASE="${XENO_BASELINES[$SAMPLE]:-$TOTAL}"
        PURITY_SCORE=$(perc "$TOTAL" "$FINAL_BASE")

        # --- Percentage Values ---
        TOTAL_P=$(perc "$TOTAL" "$TOTAL")
        MULTI_P=$(perc "$MULTI" "$TOTAL")
        UNMAP_P=$(perc "$UNMAP" "$TOTAL")
        UNI_P=$(perc "$UNI" "$TOTAL")

        GUN_P=$(perc "$GUN" "$TOTAL")
        AMB_UN_P=$(perc "$AMB_UN" "$TOTAL")
        NOF_UN_P=$(perc "$NOF_UN" "$TOTAL")

        GFW_P=$(perc "$GFW" "$TOTAL")
        AMB_FW_P=$(perc "$AMB_FW" "$TOTAL")
        NOF_FW_P=$(perc "$NOF_FW" "$TOTAL")

        GRV_P=$(perc "$GRV" "$TOTAL")
        AMB_RV_P=$(perc "$AMB_RV" "$TOTAL")
        NOF_RV_P=$(perc "$NOF_RV" "$TOTAL")

        # --- Output Row ---
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t\t%s\t%s\t%s\t\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t\t%s\t%s\t%s\t\t%s\t%s\t%s\n" \
            "$DIR_ID" "$SAMPLE" "$PURITY_SCORE" \
            "$TOTAL_P" "$MULTI_P" "$UNMAP_P" "$UNI_P" \
            "$GUN_P" "$AMB_UN_P" "$NOF_UN_P" \
            "$GFW_P" "$AMB_FW_P" "$NOF_FW_P" \
            "$GRV_P" "$AMB_RV_P" "$NOF_RV_P" \
            "$(add_commas "$TOTAL")" \
            "$(add_commas "$MULTI")" \
            "$(add_commas "$UNMAP")" \
            "$(add_commas "$UNI")" \
            "$(add_commas "$GUN")" \
            "$(add_commas "$AMB_UN")" \
            "$(add_commas "$NOF_UN")" \
            "$(add_commas "$GFW")" \
            "$(add_commas "$AMB_FW")" \
            "$(add_commas "$NOF_FW")" \
            "$(add_commas "$GRV")" \
            "$(add_commas "$AMB_RV")" \
            "$(add_commas "$NOF_RV")" \
            >> "$FINAL_OUT"

    done < <(find "$GENE_DIR" -type f -name "*.ReadsPerGene.out.tab" | sort -V)

    printf "%b\n" "$(printf '\t%.0s' {1..33})" >> "$FINAL_OUT"

done < <(find -L "$SEARCH_DIR" -type d -name "gene_counts" | sed 's/\/gene_counts$//' | sort -u)

echo "Done. STAR Report saved to: $FINAL_OUT"