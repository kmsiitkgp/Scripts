#!/bin/bash

# =============================================================================
# Copy ichorCNA reference files (GC + mappability) to shared resource directory
# =============================================================================
#
# WHY WE USE PREBUILT FILES INSTEAD OF HMMCOPY_REFWIG:
#
# - HMMCOPY_REFWIG was initially designed to generate BOTH gc.wig and map.wig
#   from a FASTA using gcCounter + mapCounter.
#
# - HOWEVER:
#     • gcCounter = valid and deterministic (can be safely recomputed)
#     • mapCounter = NOT FASTA-based; requires a precomputed mappability BigWig
#       (which we do NOT generate in this pipeline)
#
# - Therefore, HMMCOPY_REFWIG is incomplete for mappability generation and
#   should NOT be used as a source of map.wig.
#
# - Instead, we use prebuilt, validated ichorCNA reference tracks:
#     gc_hg38_1000kb.wig
#     map_hg38_1000kb.wig
#
# =============================================================================
# HOW THESE PREBUILT FILES WERE OBTAINED:
#
# - These files are distributed with the ichorCNA R package (extdata)
#   installed via Bioconductor/conda environment:
#
#   /home/kailasamms/miniconda3/envs/omics/lib/R/library/ichorCNA/extdata
#
# - They are:
#     • hg38-based reference tracks
#     • precomputed at 1Mb bin resolution
#     • generated using standardized genome + mappability assumptions
#
# - These are the SAME reference files used in published ichorCNA workflows
#   and are intended to be reused across projects (not regenerated per sample).
#
# =============================================================================
# REASON FOR CENTRALIZING IN /resources:
#
# - Ensures reproducibility across pipelines and projects
# - Avoids recomputation or accidental mismatch of GC/MAP references
# - Guarantees consistent CNV normalization across all samples
# =============================================================================

mkdir -p /home/kailasamms/resources/genomes/ichorCNA/

cp /home/kailasamms/miniconda3/envs/omics/lib/R/library/ichorCNA/extdata/gc_hg38_1000kb.wig \
   /home/kailasamms/resources/genomes/ichorCNA/

cp /home/kailasamms/miniconda3/envs/omics/lib/R/library/ichorCNA/extdata/map_hg38_1000kb.wig \
   /home/kailasamms/resources/genomes/ichorCNA/

# Verify copied reference files
ls -lh /home/kailasamms/resources/genomes/ichorCNA/