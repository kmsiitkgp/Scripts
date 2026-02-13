#!/bin/bash

FINAL_OUT="STAR_Summary.txt"

# 1. Header
HEADER="Directory_ID\tSample\tPurity (%)"
HEADER="$HEADER\tTotal (%)\tMulti (%)\tUnmapped (%)\tUnique (%)"
HEADER="$HEADER\t\tFeature U (%)\tAmbiguous U (%)\tNoFeature U (%)"
HEADER="$HEADER\t\tFeature F (%)\tAmbiguous F (%)\tNoFeature F (%)"
HEADER="$HEADER\t\tFeature R (%)\tAmbiguous R (%)\tNoFeature R (%)"
HEADER="$HEADER\t\tTotal (#)\tMulti (#)\tUnmapped (#)\tUnique (#)"
HEADER="$HEADER\t\tFeature U (#)\tAmbiguous U (#)\tNoFeature U (#)"
HEADER="$HEADER\t\tFeature F (#)\tAmbiguous F (#)\tNoFeature F (#)"
HEADER="$HEADER\t\tFeature R (#)\tAmbiguous R (#)\tNoFeature R (#)"

echo -e "$HEADER" > "$FINAL_OUT"

# Helper Functions
add_commas() { echo "$1" | sed ':a;s/\B[0-9]\{3\}\(>\|$\)/,&/;ta'; }
perc() { awk -v part=$1 -v total=$2 'BEGIN{if(total>0) printf "%.2f", (part/total)*100; else print 0}'; }

CWD=$(pwd)
PROJECT_DIRS=$(find "$CWD" -type d -name "gene_counts" | sed 's/\/gene_counts//g' | sort -u)

# --- PRE-SCAN FOR TOTAL READS (FULL) ---
declare -A FULL_TOTAL_READS

for BDIR in $PROJECT_DIRS; do
    LOWER=$(echo "$BDIR" | tr '[:upper:]' '[:lower:]')
    # We only grab the 'Total' baseline from the 'full' directory
    if [[ "$LOWER" == *"full"* ]]; then
        ALIGN_DIR="$BDIR/alignment_stats"
        for LOG in "$ALIGN_DIR"/*.Log.final.out; do
            [ -f "$LOG" ] || continue
            SNAME=$(basename "$LOG" .Log.final.out)
            # Grab the 'Number of input reads'
            T_VAL=$(grep "Number of input reads" "$LOG" | awk -F'|' '{print $2+0}')
            FULL_TOTAL_READS["$SNAME"]=$T_VAL
        done
    fi
done

echo "Generating report with Purity (Total Split / Total Full)..."

for BDIR in $PROJECT_DIRS; do
    DIR_ID=$(echo "$BDIR" | sed "s|$CWD/||g")
    ALIGN_DIR="$BDIR/alignment_stats"
    GENE_DIR="$BDIR/gene_counts"

    if [[ ! -d "$ALIGN_DIR" || ! -d "$GENE_DIR" ]]; then continue; fi

    SAMPLES=($(ls "$GENE_DIR"/*.ReadsPerGene.out.tab | xargs -n 1 basename | sed 's/.ReadsPerGene.out.tab//g' | sort -V))

    for SAMPLE in "${SAMPLES[@]}"; do
        LOG="$ALIGN_DIR/${SAMPLE}.Log.final.out"
        GFILE="$GENE_DIR/${SAMPLE}.ReadsPerGene.out.tab"
        if [[ ! -f "$LOG" ]]; then continue; fi

        # --- EXTRACT DATA ---
        TOTAL=$(grep "Number of input reads" "$LOG" | awk -F'|' '{print $2+0}')
        UNI=$(grep "Uniquely mapped reads number" "$LOG" | awk -F'|' '{print $2+0}')

        M1=$(grep "Number of reads mapped to multiple loci" "$LOG" | awk -F'|' '{print $2+0}')
        M2=$(grep "Number of reads mapped to too many loci" "$LOG" | awk -F'|' '{print $2+0}')
        MULTI=$((M1+M2))
        U1=$(grep "Number of reads unmapped: too many mismatches" "$LOG" | awk -F'|' '{print $2+0}')
        U2=$(grep "Number of reads unmapped: too short" "$LOG" | awk -F'|' '{gsub(/ /,"",$2); print $2+0}')
        UNMAP=$((U1+U2))

        # Gene logic
        GUN=$(awk 'NR>4{sum+=$2} END{print sum+0}' "$GFILE")
        GFW=$(awk 'NR>4{sum+=$3} END{print sum+0}' "$GFILE")
        GRV=$(awk 'NR>4{sum+=$4} END{print sum+0}' "$GFILE")
        AMB_UN=$(awk 'NR==4{print $2+0}' "$GFILE")
        AMB_FW=$(awk 'NR==4{print $3+0}' "$GFILE")
        AMB_RV=$(awk 'NR==4{print $4+0}' "$GFILE")
        NOF_UN=$(awk 'NR==3{print $2+0}' "$GFILE")
        NOF_FW=$(awk 'NR==3{print $3+0}' "$GFILE")
        NOF_RV=$(awk 'NR==3{print $4+0}' "$GFILE")

        # --- CALCULATE PURITY (YOUR WAY) ---
        # Total Reads this row / Total Reads from Full Library
        FULL_BASE=${FULL_TOTAL_READS["$SAMPLE"]:-$TOTAL}
        PURITY_SCORE=$(perc $TOTAL $FULL_BASE)

        # --- CONSTRUCT ROW ---
        row="$DIR_ID\t$SAMPLE\t$PURITY_SCORE"

        # Percentages
        row="$row\t$(perc $TOTAL $TOTAL)\t$(perc $MULTI $TOTAL)\t$(perc $UNMAP $TOTAL)\t$(perc $UNI $TOTAL)"
        row="$row\t\t$(perc $GUN $TOTAL)\t$(perc $AMB_UN $TOTAL)\t$(perc $NOF_UN $TOTAL)"
        row="$row\t\t$(perc $GFW $TOTAL)\t$(perc $AMB_FW $TOTAL)\t$(perc $NOF_FW $TOTAL)"
        row="$row\t\t$(perc $GRV $TOTAL)\t$(perc $AMB_RV $TOTAL)\t$(perc $NOF_RV $TOTAL)"

        # Raw
        row="$row\t\t$(add_commas $TOTAL)\t$(add_commas $MULTI)\t$(add_commas $UNMAP)\t$(add_commas $UNI)"
        row="$row\t\t$(add_commas $GUN)\t$(add_commas $AMB_UN)\t$(add_commas $NOF_UN)"
        row="$row\t\t$(add_commas $GFW)\t$(add_commas $AMB_FW)\t$(add_commas $NOF_FW)"
        row="$row\t\t$(add_commas $GRV)\t$(add_commas $AMB_RV)\t$(add_commas $NOF_RV)"

        echo -e "$row" >> "$FINAL_OUT"
    done

    # Visual spacer
    spacer=""
    for i in {1..33}; do spacer="$spacer\t"; done
    echo -e "$spacer" >> "$FINAL_OUT"
done

echo "Report generated! Purity column = (Total Reads / Full Library Total Reads) * 100"