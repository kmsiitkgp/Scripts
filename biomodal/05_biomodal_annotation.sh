#!/bin/bash
#SBATCH --job-name=Biomodal_annotation
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4
#SBATCH --output=%x_%j.log

#===============================================================================
# 📄 Script   : 05_biomodal_annotation.sh
# 📌 Purpose  : One-time preparation of genomic annotation BED files for use
#               with modality get and modality dmr call commands.
#
#               Generates TWO annotation BED files from:
#                 - Ensembl GTF (GRCh38 v115) for gene features
#                 - UCSC CpG islands (hg38) for CpG features
#
#               File 1 — Subgene: Promoter, 5UTR, Exon, Intron, 3UTR,
#                                 CpG_island, CpG_shore, CpG_shelf
#               File 2 — Gene:    Promoter, Gene_body,
#                                 CpG_island, CpG_shore, CpG_shelf
#
# 🔁 Usage    : Run ONCE before 06_biomodal_get.sh and 07_biomodal_dmr.sh
#               bash 05_biomodal_annotation.sh
#
# ⚠️  Requires : bedtools (available in omics conda env)
#===============================================================================

#===============================================================================
# 🚀 Environment
#===============================================================================

source "/home/kailasamms/miniconda3/etc/profile.d/conda.sh"
conda activate omics
echo "✅ Environment activated: omics"

set -euo pipefail    # Exit on error, undefined variable, or pipe failure
export LC_ALL=C      # Consistent sort order across systems

#===============================================================================
# ⚙️ Paths
#===============================================================================

REF_DIR="/home/kailasamms/resources/genomes"
BIO_DIR="${REF_DIR}/biomodal"

# --- Input files ---
# Ensembl GTF v115 — source of all gene/transcript/exon/UTR features
GTF="${REF_DIR}/Homo_sapiens.GRCh38.115.gtf"

# UCSC CpG islands (hg38) — source for CpG_island, CpG_shore, CpG_shelf
# Download from : https://genome.ucsc.edu/cgi-bin/hgTables
#   Genome       : Human
#   Assembly     : hg38
#   Group        : Regulation
#   Track        : CpG Islands
#   Table        : cpgIslandExt
#   Region       : genome
#   Output format: BED
#   Save as      : hg38_CpG_islands.bed.gz
CPG_ISLANDS="${BIO_DIR}/hg38_CpG_islands.bed.gz"

# --- Output files ---
BED_SUBGENE="${REF_DIR}/GRCh38_v115_biomodal_annotation_subgene.bed"
BED_GENE="${REF_DIR}/GRCh38_v115_biomodal_annotation_gene.bed"

# --- Chromosome filter ---
# Keep only canonical chromosomes: chr1-chr22, chrX, chrY, chrM
# This strips alt contigs, random scaffolds, unplaced contigs etc.
# Applied consistently in EVERY awk block after chr name is constructed
CHR_REGEX='^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$'

#===============================================================================
# 📋 Job Info
#===============================================================================

echo "=========================================================="
printf "Start date      : %s\n" "$(date)"
echo "=========================================================="
printf "REF_DIR         : %s\n" "${REF_DIR}"
printf "BIO_DIR         : %s\n" "${BIO_DIR}"
printf "GTF             : %s\n" "${GTF}"
printf "CpG islands     : %s\n" "${CPG_ISLANDS}"
printf "Output subgene  : %s\n" "${BED_SUBGENE}.gz"
printf "Output gene     : %s\n" "${BED_GENE}.gz"
echo "=========================================================="

#===============================================================================
# 🔍 Check dependencies and inputs
#===============================================================================

echo "Checking dependencies..."
command -v bedtools >/dev/null || { echo "❌ bedtools not found"; exit 1; }
echo "✅ bedtools $(bedtools --version)"

echo "Checking input files..."
[[ -f "${GTF}" ]]         || { echo "❌ GTF not found: ${GTF}"; exit 1; }
[[ -f "${CPG_ISLANDS}" ]] || { echo "❌ CpG islands not found: ${CPG_ISLANDS}"; exit 1; }
echo "✅ GTF: ${GTF}"
echo "✅ CpG islands: ${CPG_ISLANDS}"

# --- Skip if output files already exist ---
if [[ -f "${BED_SUBGENE}.gz" && -f "${BED_GENE}.gz" ]]; then
    echo "✅ Annotation BED files already exist — skipping regeneration."
    echo "   ${BED_SUBGENE}.gz"
    echo "   ${BED_GENE}.gz"
    echo "   Delete these files and rerun to regenerate."
    exit 0
fi

#===============================================================================
# 🧬 CpG Features
#
# Built once and shared by BOTH output BED files.
#
# Mutual exclusivity is guaranteed by correct subtraction order:
#   Step A — islands : raw (highest priority, nothing subtracted)
#   Step B — shores  : raw_shores - islands
#   Step C — shelves : raw_shelves - shores_CLEAN - islands
#
# Why subtract shores_CLEAN (not raw shores) for shelves?
# If raw shores were used, any island bases still lurking inside raw shores
# would re-introduce overlap into shelves. Using cleaned shores guarantees
# full transitivity: shelf ∩ shore = ∅  and  shelf ∩ island = ∅
#===============================================================================

CPG_TMP="${BIO_DIR}/cpg_tmp.bed"

echo "--- CpG: Decompressing and filtering CpG islands..."
zcat "${CPG_ISLANDS}" | awk -v regex="${CHR_REGEX}" '
    BEGIN { OFS="\t" }
    {
        chr = $1
        if (chr !~ /^chr/) chr = "chr" chr    # add chr prefix if missing
        if (chr == "chrMT") chr = "chrM"       # Ensembl MT → UCSC M
        if (chr !~ regex) next                 # drop alt/random contigs
        print chr, $2, $3
    }' | sort -k1,1 -k2,2n > "${CPG_TMP}"
echo "✅ CpG islands (filtered + sorted): $(wc -l < "${CPG_TMP}") regions"

# Step A — CpG Islands (raw, highest priority — nothing subtracted)
echo "--- CpG: Step A — Extracting islands..."
awk 'BEGIN{OFS="\t"} { print $1, $2, $3, "CpG_island", "NA" }' \
    "${CPG_TMP}" > "${BIO_DIR}/cpg_islands.bed"
echo "✅ CpG islands: $(wc -l < "${BIO_DIR}/cpg_islands.bed") regions"

# Step B — CpG Shores: 2kb flanks on each side of island, minus islands
echo "--- CpG: Step B — Deriving shores (islands subtracted)..."
awk 'BEGIN{OFS="\t"} {
    s = $2 - 2000; if (s < 0) s = 0
    if (s < $2) print $1, s,  $2,      "CpG_shore", "NA"   # upstream flank
    print         $1, $3, $3+2000,     "CpG_shore", "NA"   # downstream flank
}' "${CPG_TMP}" | \
    sort -k1,1 -k2,2n | \
    bedtools subtract -a stdin -b "${BIO_DIR}/cpg_islands.bed" \
    > "${BIO_DIR}/cpg_shores.bed"
echo "✅ CpG shores (islands removed): $(wc -l < "${BIO_DIR}/cpg_shores.bed") regions"

# Step C — CpG Shelves: 2-4kb flanks, minus shores_CLEAN then minus islands
echo "--- CpG: Step C — Deriving shelves (shores + islands subtracted)..."
awk 'BEGIN{OFS="\t"} {
    s = $2 - 4000; if (s < 0) s = 0
    e = $2 - 2000; if (e < 0) e = 0
    if (s < e) print $1, s,      e,      "CpG_shelf", "NA"   # upstream flank
    print         $1, $3+2000,  $3+4000, "CpG_shelf", "NA"   # downstream flank
}' "${CPG_TMP}" | \
    sort -k1,1 -k2,2n | \
    bedtools subtract -a stdin -b "${BIO_DIR}/cpg_shores.bed"  | \
    bedtools subtract -a stdin -b "${BIO_DIR}/cpg_islands.bed" \
    > "${BIO_DIR}/cpg_shelves.bed"
echo "✅ CpG shelves (shores + islands removed): $(wc -l < "${BIO_DIR}/cpg_shelves.bed") regions"

rm -f "${CPG_TMP}"

#===============================================================================
# 🧬 Subgene Features (File 1)
#
# Transcript-level features extracted from protein-coding genes only.
# One entry per isoform — captures all known splice variants.
#
# Features: Promoter, 5UTR, Exon (CDS), Intron, 3UTR
#
# Promoter definition: 1kb upstream of TSS (strand-aware)
#   + strand: TSS - 1000 to TSS
#   - strand: TES to TES + 1000
# NOTE: end coordinate is TSS itself ($4 for + strand) so the TSS CpG is
# included in the promoter window — this is the biologically correct behavior
#
# Introns are DERIVED (not extracted directly from GTF) as:
#   Gene_span - (CDS + 5UTR + 3UTR)
# This is standard practice for methylation — avoids isoform complexity
#===============================================================================

echo "--- SUBGENE: Extracting transcript-level promoters..."
awk -v regex="${CHR_REGEX}" '
    BEGIN{OFS="\t"}
    /^#/ {next}
    $3=="transcript" &&
    $0~/gene_biotype "protein_coding"/ {
        match($0, /gene_name "([^"]+)"/, arr)
        chr = ($1=="MT") ? "chrM" : "chr"$1
        if (chr !~ regex) next
        if ($7=="+") { start=$4-1-1000; end=$4-1 }
        else          { start=$5;     end=$5+1000 }
        if (start<0) start=0
        print chr, start, end, "Promoter", arr[1]
    }' "${GTF}" | sort -k1,1 -k2,2n > "${BIO_DIR}/sg_promoters.bed"
echo "✅ Transcript-level promoters: $(wc -l < "${BIO_DIR}/sg_promoters.bed") regions"

echo "--- SUBGENE: Extracting 5' UTRs..."
awk -v regex="${CHR_REGEX}" '
    BEGIN{OFS="\t"}
    /^#/ {next}
    $3=="five_prime_utr" &&
    $0~/gene_biotype "protein_coding"/ {
        match($0, /gene_name "([^"]+)"/, arr)
        chr = ($1=="MT") ? "chrM" : "chr"$1
        if (chr !~ regex) next
        start=$4-1; if(start<0) start=0
        print chr, start, $5, "5UTR", arr[1]
    }' "${GTF}" | sort -k1,1 -k2,2n > "${BIO_DIR}/sg_5utr.bed"
echo "✅ 5' UTRs: $(wc -l < "${BIO_DIR}/sg_5utr.bed") regions"

echo "--- SUBGENE: Extracting CDS exons..."
awk -v regex="${CHR_REGEX}" '
    BEGIN{OFS="\t"}
    /^#/ {next}
    $3=="CDS" &&
    $0~/gene_biotype "protein_coding"/ {
        match($0, /gene_name "([^"]+)"/, arr)
        chr = ($1=="MT") ? "chrM" : "chr"$1
        if (chr !~ regex) next
        start=$4-1; if(start<0) start=0
        print chr, start, $5, "Exon", arr[1]
    }' "${GTF}" | sort -k1,1 -k2,2n > "${BIO_DIR}/sg_exons.bed"
echo "✅ CDS exons: $(wc -l < "${BIO_DIR}/sg_exons.bed") regions"

echo "--- SUBGENE: Extracting 3' UTRs..."
awk -v regex="${CHR_REGEX}" '
    BEGIN{OFS="\t"}
    /^#/ {next}
    $3=="three_prime_utr" &&
    $0~/gene_biotype "protein_coding"/ {
        match($0, /gene_name "([^"]+)"/, arr)
        chr = ($1=="MT") ? "chrM" : "chr"$1
        if (chr !~ regex) next
        start=$4-1; if(start<0) start=0
        print chr, start, $5, "3UTR", arr[1]
    }' "${GTF}" | sort -k1,1 -k2,2n > "${BIO_DIR}/sg_3utr.bed"
echo "✅ 3' UTRs: $(wc -l < "${BIO_DIR}/sg_3utr.bed") regions"

echo "--- SUBGENE: Extracting gene spans (for intron derivation)..."
awk -v regex="${CHR_REGEX}" '
    BEGIN{OFS="\t"}
    /^#/ {next}
    $3=="gene" &&
    $0~/gene_biotype "protein_coding"/ {
        match($0, /gene_name "([^"]+)"/, arr)
        chr = ($1=="MT") ? "chrM" : "chr"$1
        if (chr !~ regex) next
        start=$4-1; if(start<0) start=0
        print chr, start, $5, "Gene", arr[1]
    }' "${GTF}" | sort -k1,1 -k2,2n > "${BIO_DIR}/sg_genes_temp.bed"

echo "--- SUBGENE: Deriving introns (gene span minus CDS + 5UTR + 3UTR)..."
cat "${BIO_DIR}/sg_exons.bed" \
    "${BIO_DIR}/sg_5utr.bed" \
    "${BIO_DIR}/sg_3utr.bed" | \
    sort -k1,1 -k2,2n > "${BIO_DIR}/sg_coding_temp.bed"

bedtools subtract \
    -a "${BIO_DIR}/sg_genes_temp.bed" \
    -b "${BIO_DIR}/sg_coding_temp.bed" | \
    awk 'BEGIN{OFS="\t"} { print $1, $2, $3, "Intron", $5 }' \
    > "${BIO_DIR}/sg_introns.bed"
echo "✅ Introns: $(wc -l < "${BIO_DIR}/sg_introns.bed") regions"

#===============================================================================
# 🧬 Gene-level Features (File 2)
#
# Gene-level (not transcript-level) — one entry per gene.
# Simpler than subgene — useful for broad regional methylation summaries.
#
# Promoter definition: 1kb upstream of gene TSS (strand-aware)
# Gene_body: full gene span (start to end including introns)
#===============================================================================

echo "--- GENE: Extracting gene bodies and promoters..."
awk -v regex="${CHR_REGEX}" \
    -v body_file="${BIO_DIR}/gl_gene_body_unsorted.bed" \
    -v prom_file="${BIO_DIR}/gl_promoters_unsorted.bed" '
    BEGIN{OFS="\t"}
    /^#/ {next}
    $3=="gene" &&
    $0~/gene_biotype "protein_coding"/ {
        match($0, /gene_name "([^"]+)"/, arr)
        gene   = arr[1]
        chr    = ($1=="MT") ? "chrM" : "chr"$1
        if (chr !~ regex) next
        start  = $4-1; if (start<0) start=0
        end    = $5
        strand = $7
        print chr, start, end, "Gene_body", gene > body_file
        if (strand=="+") { ps=start-1000; if(ps<0) ps=0; pe=start }
        else              { ps=end; pe=end+1000 }
        print chr, ps, pe, "Promoter", gene > prom_file
    }' "${GTF}"

sort -k1,1 -k2,2n "${BIO_DIR}/gl_gene_body_unsorted.bed" \
    > "${BIO_DIR}/gl_gene_body.bed"
sort -k1,1 -k2,2n "${BIO_DIR}/gl_promoters_unsorted.bed" \
    > "${BIO_DIR}/gl_promoters.bed"
rm -f "${BIO_DIR}/gl_gene_body_unsorted.bed" \
      "${BIO_DIR}/gl_promoters_unsorted.bed"

echo "✅ Gene bodies:          $(wc -l < "${BIO_DIR}/gl_gene_body.bed") regions"
echo "✅ Gene-level promoters: $(wc -l < "${BIO_DIR}/gl_promoters.bed") regions"

#===============================================================================
# 📦 Merge into two final BED files
#===============================================================================

# --- File 1: Subgene-level ---
echo "Merging subgene-level BED..."
echo -e "Chromosome\tStart\tEnd\tAnnotation\tName" > "${BED_SUBGENE}"
cat \
    "${BIO_DIR}/sg_promoters.bed" \
    "${BIO_DIR}/sg_5utr.bed" \
    "${BIO_DIR}/sg_exons.bed" \
    "${BIO_DIR}/sg_introns.bed" \
    "${BIO_DIR}/sg_3utr.bed" \
    "${BIO_DIR}/cpg_islands.bed" \
    "${BIO_DIR}/cpg_shores.bed" \
    "${BIO_DIR}/cpg_shelves.bed" | \
    sort -k1,1 -k2,2n >> "${BED_SUBGENE}"
echo "✅ Subgene BED written: $(tail -n +2 "${BED_SUBGENE}" | wc -l) regions"

# --- File 2: Gene-level ---
echo "Merging gene-level BED..."
echo -e "Chromosome\tStart\tEnd\tAnnotation\tName" > "${BED_GENE}"
cat \
    "${BIO_DIR}/gl_promoters.bed" \
    "${BIO_DIR}/gl_gene_body.bed" \
    "${BIO_DIR}/cpg_islands.bed" \
    "${BIO_DIR}/cpg_shores.bed" \
    "${BIO_DIR}/cpg_shelves.bed" | \
    sort -k1,1 -k2,2n >> "${BED_GENE}"
echo "✅ Gene-level BED written: $(tail -n +2 "${BED_GENE}" | wc -l) regions"

#===============================================================================
# 🗜️ Compress outputs
#===============================================================================

echo "Compressing..."
gzip "${BED_SUBGENE}"
gzip "${BED_GENE}"

#===============================================================================
# 🧹 Cleanup intermediate files
#===============================================================================

echo "Cleaning up temporary files..."
rm -f \
    "${BIO_DIR}/sg_genes_temp.bed" \
    "${BIO_DIR}/sg_coding_temp.bed"

# Uncomment below to also remove intermediate per-feature BEDs after confirming
# outputs look correct — kept by default for troubleshooting:
rm -f "${BIO_DIR}/sg_promoters.bed"  "${BIO_DIR}/sg_5utr.bed"   \
       "${BIO_DIR}/sg_exons.bed"      "${BIO_DIR}/sg_introns.bed" \
       "${BIO_DIR}/sg_3utr.bed"       "${BIO_DIR}/gl_gene_body.bed" \
       "${BIO_DIR}/gl_promoters.bed"  "${BIO_DIR}/cpg_islands.bed" \
       "${BIO_DIR}/cpg_shores.bed"    "${BIO_DIR}/cpg_shelves.bed"

#===============================================================================
# 📋 Summary
#===============================================================================

echo ""
echo "=========================================================="
printf "End date        : %s\n" "$(date)"
echo "=========================================================="
echo "✅ Done! Annotation BED files written to: ${BIO_DIR}"
echo ""
echo "📄 File 1 — Subgene: ${BED_SUBGENE}.gz"
echo "   Total regions: $(zcat "${BED_SUBGENE}.gz" | tail -n +2 | wc -l)"
echo "   Breakdown:"
zcat "${BED_SUBGENE}.gz" | awk 'NR>1 {print $4}' | sort | uniq -c | sort -rn
echo ""
echo "📄 File 2 — Gene: ${BED_GENE}.gz"
echo "   Total regions: $(zcat "${BED_GENE}.gz" | tail -n +2 | wc -l)"
echo "   Breakdown:"
zcat "${BED_GENE}.gz" | awk 'NR>1 {print $4}' | sort | uniq -c | sort -rn
echo "=========================================================="
