#!/bin/bash

OUT_FILE="Salmon_Final_Report.txt"

# Standardized Header in your specific order
HEADER="Sample\tLibrary\tStrand Bias Score\tTotal\tDiscarded (Score)\tDiscarded (Decoys)\tDiscarded (Dovetail)\tMapped\tMapping Rate\tGCBias\tSeqBias\tPosBias\tSalmon version"
echo -e "$HEADER" > $OUT_FILE

echo "Generating detailed report from Salmon outputs..."

for dir in */; do
    SAMPLE=${dir%/}

    # Paths to the key files
    CMD="${dir}/cmd_info.json"
    LIB="${dir}/lib_format_counts.json"
    META="${dir}/aux_info/meta_info.json"
    LOG="${dir}/logs/salmon_quant.log"

    if [[ -f "$LIB" && -f "$META" ]]; then

        # 1. Geometry (from lib_format_counts.json)
        LIBRARY=$(grep "expected_format" "$LIB" | cut -d':' -f2 | tr -d '", ')
        STRAND_BIAS=$(grep "strand_mapping_bias" "$LIB" | cut -d':' -f2 | tr -d '", ')

        # 2. Results & Bias (from meta_info.json)
        TOTAL=$(grep "num_processed" "$META" | cut -d':' -f2 | tr -d '", ')
        MAPPED=$(grep "num_mapped" "$META" | cut -d':' -f2 | tr -d '", ')
        MAPPING_RATE=$(grep "percent_mapped" "$META" | cut -d':' -f2 | tr -d '", ')

        GCB=$(grep "gc_bias_correct" "$META" | cut -d':' -f2 | tr -d '", ' | tr '[:lower:]' '[:upper:]')
        SEQB=$(grep "seq_bias_correct" "$META" | cut -d':' -f2 | tr -d '", ' | tr '[:lower:]' '[:upper:]')

        # 3. Deep Log Analysis (using exact strings provided)
        if [[ -f "$LOG" ]]; then
            # Entirely discarded because of alignment score
            DISCARD_SCORE=$(grep "Number of fragments entirely discarded because of alignment score" "$LOG" | awk -F': ' '{print $2}' | tr -d ',' | head -n 1)

            # Best-mapped to decoys
            DISCARD_DECOY=$(grep "Number of fragments discarded because they are best-mapped to decoys" "$LOG" | awk -F': ' '{print $2}' | tr -d ',' | head -n 1)

            # Dovetail (discordant)
            DISCARD_DOVE=$(grep "Number of fragments discarded because they have only dovetail (discordant) mappings to valid targets" "$LOG" | awk -F': ' '{print $2}' | tr -d ',' | head -n 1)

            # PosBias Logic (Checking the LOG as JSON often omits it)
            if grep -q "Computed positional bias model" "$LOG"; then
                POSB="TRUE"
            else
                POSB="FALSE"
            fi
        else
            DISCARD_SCORE="0"
            DISCARD_DECOY="0"
            DISCARD_DOVE="0"
            POSB="FALSE"
        fi

        # 4. Version (from cmd_info.json)
        VERSION=$(grep "salmon_version" "$CMD" | cut -d':' -f2 | tr -d '", ')

        # Write data row in requested order
        echo -e "$SAMPLE\t$LIBRARY\t$STRAND_BIAS\t$TOTAL\t$DISCARD_SCORE\t$DISCARD_DECOY\t$DISCARD_DOVE\t$MAPPED\t$MAPPING_RATE\t$GCB\t$SEQB\t$POSB\t$VERSION" >> $OUT_FILE
        echo "Processed: $SAMPLE"
    else
        echo "Skipping: $SAMPLE (Required files missing)"
    fi
done

echo "-----------------------------------------------"
echo "Report saved to: $OUT_FILE"