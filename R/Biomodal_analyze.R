library(readr)
library(dplyr)
library(openxlsx)

path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Biomodal"

# ---- Cross check biomodal outputs for accuracy ----

# (i) Number of CpG sites per region ( all 3 files MUST be identical)
count_hmc <- readr::read_tsv(file.path(path, "biomodal_output", "count_hmc.tsv"), show_col_types = FALSE)
count_mc <- readr::read_tsv(file.path(path, "biomodal_output", "count_mc.tsv"), show_col_types = FALSE)
count_total_c <- readr::read_tsv(file.path(path, "biomodal_output", "count_total_c.tsv"), show_col_types = FALSE)

all.equal(as.data.frame(count_hmc),  as.data.frame(count_mc), check.attributes = FALSE)
all.equal(as.data.frame(count_hmc),  as.data.frame(count_total_c),  check.attributes = FALSE)

# (ii) Fraction of 5mc and 5hmc
frac_hmc <- readr::read_tsv(file.path(path, "biomodal_output", "frac_hmc.tsv"), show_col_types = FALSE)
frac_mc <- readr::read_tsv(file.path(path, "biomodal_output", "frac_mc.tsv"), show_col_types = FALSE)

# Number of CpG sites * reads per region
sum_hmc <- readr::read_tsv(file.path(path, "biomodal_output", "sum_hmc.tsv"), show_col_types = FALSE)
sum_mc <- readr::read_tsv(file.path(path, "biomodal_output", "sum_mc.tsv"), show_col_types = FALSE)
sum_total_c <- readr::read_tsv(file.path(path, "biomodal_output", "sum_total_c.tsv"), show_col_types = FALSE)

sum(round(sum_hmc[4:27]/sum_total_c[4:27] - frac_hmc[4:27],2), na.rm=TRUE)
sum(round(sum_mc[4:27]/sum_total_c[4:27] - frac_mc[4:27],2), na.rm=TRUE)

rm(list = intersect(ls(), c("count_hmc", "count_mc", "count_total_c", "frac_mc", "frac_hmc", "sum_hmc", "sum_mc", "sum_total_c")))

# ---- Get genes that already have 5mc and 5hmc in normal samples ----

# http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/

#  0️⃣ Setup
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(openxlsx)

# Input files
normal_5hmc_files <- list.files(file.path(path, "Benign_prostate_5hmc"), full.names = TRUE)
normal_5mc_files <- list.files(file.path(path, "Healthy_plasma_5mc"), full.names = TRUE)

# 1️⃣ Import GTF
gtf_file <- file.path(path, "biomodal_output", "hg19.refGene.gtf.gz")         # RefSeq GTF
gtf <- rtracklayer::import(con = gtf_file, format = "gtf")

# 2️⃣ Define regions to match with Biomodal 
transcripts <- gtf[GenomicRanges::mcols(gtf)$type == "transcript"]

# GRange object for gene body
gene_body_gr <- transcripts

# GRange object for promoters
promoter_gr <- GenomicRanges::promoters(x = transcripts, 
                                        upstream = 1000, 
                                        downstream = 0)

# GRange object around TSS
tss_region_gr <- GenomicRanges::promoters(x = transcripts, 
                                          upstream = 200, 
                                          downstream = 200)

# # GRange object at TSS
# tss_pos <- dplyr::if_else(condition = as.logical(GenomicRanges::strand(transcripts) == "+"),
#                           true = GenomicRanges::start(transcripts),
#                           false = GenomicRanges::end(transcripts))
# tss_gr <- GenomicRanges::GRanges(seqnames = GenomicRanges::seqnames(transcripts),
#                                  ranges = IRanges::IRanges(start = tss_pos, end = tss_pos),
#                                  strand = GenomicRanges::strand(transcripts),
#                                  gene_id = transcripts$gene_id)

# 3️⃣ Annotate peaks for each file 

# Helper function to get collapsed gene names for a single Hits object
get_collapsed_genes <- function(peaks, hits, subject_gr) {
  
  # 1. Create a factor for ALL possible peak indices (1 to N)
  # This ensures that even peaks with 0 hits get an entry in the list
  query_factor <- factor(S4Vectors::queryHits(hits), levels = seq_along(peaks))
  
  # 2. Split subject hits based on that factor
  # This avoids the 'drop' error and keeps the list length = length(peaks)
  hits_by_query <- S4Vectors::splitAsList(S4Vectors::subjectHits(hits), query_factor)
  
  # 3. Extract gene names for all peaks
  sapply(hits_by_query, function(idx) {
    # If a peak had no hits, idx will be an empty integer vector
    if (length(idx) == 0) return("NA")
    
    unique_names <- unique(S4Vectors::mcols(subject_gr)$gene_name[idx])
    
    if (length(unique_names) == 0) {
      return("NA")
    } else {
      return(paste(unique_names, collapse = ";"))
    }
  })
}

annotate_peaks <- function(peak_file, format) {
  
  # Read the peak file
  peaks <- rtracklayer::import(con = peak_file, format = format)
  
  # Overlaps
  promoter_hits   <- GenomicRanges::findOverlaps(query = peaks, subject = promoter_gr)
  tss_hits        <- GenomicRanges::findOverlaps(query = peaks, subject = tss_region_gr)
  gene_body_hits <- GenomicRanges::findOverlaps(query = peaks, subject = gene_body_gr)
  
  # Build annotation table
  peak_df <- data.frame(peak_name     = GenomicRanges::mcols(peaks)$name,
                        peak_chr      = GenomicRanges::seqnames(peaks),
                        peak_start    = GenomicRanges::start(peaks),
                        peak_end      = GenomicRanges::end(peaks),
                        in_promoter   = GenomicRanges::countOverlaps(query = peaks, subject = promoter_gr) > 0,
                        in_tss        = GenomicRanges::countOverlaps(query = peaks, subject = tss_region_gr) > 0,
                        in_gene_body  = GenomicRanges::countOverlaps(query = peaks, subject = gene_body_gr) > 0, 
                        promoter_gene_names  = get_collapsed_genes(peaks, promoter_hits, promoter_gr),
                        tss_gene_names       = get_collapsed_genes(peaks, tss_hits, tss_region_gr),
                        gene_body_gene_names = get_collapsed_genes(peaks, gene_body_hits, gene_body_gr)
  )
  return(peak_df)
}

normal_5hmc <- lapply(normal_5hmc_files, function(f) {annotate_peaks(f, "narrowPeak")})
normal_5mc <- lapply(normal_5mc_files, function(f) {annotate_peaks(f, "BED")})

annotated_list <- list()
annotated_list[["normal_5hmc"]] <- normal_5hmc 
annotated_list[["normal_5mc"]] <- normal_5mc

all_promoter_genes <- list()
all_tss_genes      <- list()
all_gene_body_genes <- list()

for (i in names(annotated_list)){
  
  # 4️⃣ Merge genes across samples
  
  # Get all gene IDs in promoter/TSS/gene body for each sample
  promoter_genes <- lapply(annotated_list[[i]], function(df) {
    df %>%
      dplyr::filter(in_promoter & promoter_gene_names != "NA") %>%
      tidyr::separate_longer_delim(cols = promoter_gene_names, delim = ";") %>%
      dplyr::pull(promoter_gene_names) %>%
      unique()
  })
  
  tss_genes <- lapply(annotated_list[[i]], function(df) {
    df %>%
      dplyr::filter(in_tss & tss_gene_names != "NA") %>%
      tidyr::separate_longer_delim(cols = tss_gene_names, delim = ";") %>%
      dplyr::pull(tss_gene_names) %>%
      unique()
  })
  
  gene_body_genes <- lapply(annotated_list[[i]], function(df) {
    df %>%
      dplyr::filter(in_gene_body & gene_body_gene_names != "NA") %>%
      tidyr::separate_longer_delim(cols = gene_body_gene_names, delim = ";") %>%
      dplyr::pull(gene_body_gene_names) %>%
      unique()
  })
  
  # Save for diagnostics BEFORE intersecting
  all_promoter_genes[[i]]  <- promoter_genes
  all_tss_genes[[i]]       <- tss_genes
  all_gene_body_genes[[i]] <- gene_body_genes
  
  # Genes present in ALL normal samples
  common_promoter_genes <- purrr::reduce(.x = promoter_genes, .f = intersect)
  common_tss_genes      <- purrr::reduce(.x = tss_genes, .f = intersect)
  common_gene_body_genes <- purrr::reduce(.x = gene_body_genes, .f = intersect)
  
  # 5️⃣ Save results
  
  # Create data frames for each region with a column indicating the region
  # Using Name, Annotation, Promoter, TSS, Gene to match biomodal DMR xlsx
  promoter_df   <- data.frame(Name = common_promoter_genes, Annotation = "Promoter")
  tss_df        <- data.frame(Name = common_tss_genes, Annotation = "TSS")
  gene_body_df  <- data.frame(Name = common_gene_body_genes, Annotation = "Gene")
  
  # Merge all into one data frame
  all_genes_df <- bind_rows(promoter_df, tss_df, gene_body_df) %>%
    distinct()  # Remove duplicates if a gene is in multiple regions
  
  # Write to Excel
  wb <- createWorkbook()
  addWorksheet(wb, i)
  writeData(wb, sheet = i, all_genes_df)
  saveWorkbook(wb, file.path(path, paste0(i, ".xlsx")), overwrite = TRUE)
}


# ---- Identify DMR in tumor ----

# 1. Reference Data (Normal DNA)
ref_hmc <- read.xlsx(file.path(path, "normal_5hmc.xlsx"))
ref_mc  <- read.xlsx(file.path(path, "normal_5mc.xlsx"))

# 2. Raw cfDNA Data ( = ctDNA + Normal DNA)
# raw_hmc <- readr::read_tsv(file.path(path, "biomodal_output", "DMR_20251124_130348_DMR_hmc_Pre__Post_20251124_130348.tsv"), show_col_types = FALSE)
# raw_mc  <- readr::read_tsv(file.path(path, "biomodal_output", "DMR_20251124_130348_DMR_mc_Pre__Post_20251124_130348.tsv"), show_col_types = FALSE)

raw_hmc <- readr::read_tsv(file.path(path, "biomodal_output", "DMR_20251208_144543_DMR_hmc_Pre__Post_20251208_144543.tsv"), show_col_types = FALSE)
raw_mc  <- readr::read_tsv(file.path(path, "biomodal_output", "DMR_20251208_144543_DMR_mc_Pre__Post_20251208_144543.tsv"), show_col_types = FALSE)

# 3. Filtered Significant cfDNA DMRs
sig_hmc <- raw_hmc %>% 
  dplyr::mutate(mean_mod_group_1 = ifelse(mean_mod_group_1 == 0, 1e-6, mean_mod_group_1),
                mean_mod_group_2 = ifelse(mean_mod_group_2 == 0, 1e-6, mean_mod_group_2),
                mod_fold_change  = mean_mod_group_2 / mean_mod_group_1,
                mod_lfc          = log2(mod_fold_change)) %>%
  dplyr::filter(dmr_qvalue <= 0.05, num_contexts >= 10, abs(mod_lfc) >= 0.58, !is.na(Name))

sig_mc <- raw_mc %>% 
  dplyr::mutate(mean_mod_group_1 = ifelse(mean_mod_group_1 == 0, 1e-6, mean_mod_group_1),
                mean_mod_group_2 = ifelse(mean_mod_group_2 == 0, 1e-6, mean_mod_group_2),
                mod_fold_change  = mean_mod_group_2 / mean_mod_group_1,
                mod_lfc          = log2(mod_fold_change)) %>% 
  dplyr::filter(dmr_qvalue <= 0.05, num_contexts >= 10, abs(mod_lfc) >= 0.58, !is.na(Name))

# 4. Final Tumor-Specific DMRs (Normal-subtracted)
# Remove any DMR whose gene name appears in the normal reference (ref_hmc/ref_mc),
# regardless of annotation label (Promoter/TSS/Gene). This gene-level subtraction
# is intentionally conservative: since the normal reference (hg19) and Biomodal 
# output (hg38) use different genome builds and platforms, annotation labels for 
# the same gene can differ due to coordinate shifts at region boundaries. Matching 
# on gene name alone avoids retaining false positives caused by this label 
# discordance. The cost is ~5% fewer hits, which is acceptable given the 
# specificity gain.
tumor_hmc <- sig_hmc %>%
  anti_join(ref_hmc, by = c("Name")) %>%
  arrange(Annotation, Name) %>%
  dplyr::mutate(key = paste(Chromosome, Start, End, Name, Annotation, sep = "_"))

tumor_mc <- sig_mc %>%
  anti_join(ref_mc, by = c("Name", "Annotation")) %>%
  arrange(Annotation, Name) %>%
  dplyr::mutate(key = paste(Chromosome, Start, End, Name, Annotation, sep = "_"),
                # Flag genes that would be removed under strict filtering
                in_strict_ref = Name %in% ref_mc$Name)

write.xlsx(sig_hmc, 
           file = file.path(path, "Tumor+Normal_5hmC_DMRs.xlsx"), 
           overwrite = TRUE)

write.xlsx(sig_mc, 
           file = file.path(path, "Tumor+Normal_5mC_DMRs.xlsx"), 
           overwrite = TRUE)

write.xlsx(tumor_hmc, 
           file = file.path(path, "Tumor_Specific_5hmC_DMRs.xlsx"), 
           overwrite = TRUE)

write.xlsx(tumor_mc, 
           file = file.path(path, "Tumor_Specific_5mC_DMRs.xlsx"), 
           overwrite = TRUE)


# ---- Patient wise analysis ----

# 1. Load Metadata and Define Patient Pairs
metadata <- read.xlsx(file.path(path, "Metadata.xlsx"))
samples  <- unique(metadata$Sample_ID)

# 2. Generate all pairs and find C1 vs. Post-baseline matches for the same patient
combns      <- utils::combn(x = samples, m = 2)
p_controls  <- c()
p_expts     <- c()

for (i in 1:ncol(combns)) {
  # Strip timepoint labels to see if IDs match (e.g., Patient_A_C1 -> Patient_A)
  id_a <- gsub("C1|C2|C3|EOT", "", combns[1, i])
  id_b <- gsub("C1|C2|C3|EOT", "", combns[2, i])
  
  if (id_a == id_b) {
    # Ensure we always treat C1 as the baseline/control
    if (grepl("C1", combns[1, i]) & !grepl("C1", combns[2, i])) {
      p_controls <- c(p_controls, combns[1, i])
      p_expts    <- c(p_expts, combns[2, i])
    } 
  }
}

patient_comparisons <- list(control = p_controls, expt = p_expts)

# 3. Process Methylation Fractions
# Load raw fractions
frac_hmc_raw <- readr::read_tsv(file.path(path, "biomodal_output", "frac_hmc.tsv"), show_col_types = FALSE)
frac_mc_raw  <- readr::read_tsv(file.path(path, "biomodal_output", "frac_mc.tsv"), show_col_types = FALSE)

process_fracs <- function(df, suffix) {
  # 1. Clean column names using the specific suffix (hmc or mc)
  colnames(df) <- gsub(suffix, "", colnames(df))
  
  # 2. Fix NAs and add counters
  df %>%
    dplyr::mutate(across(where(is.numeric), ~replace(.x, is.na(.x), 0))) %>%
    dplyr::mutate(n_UP = 0, n_DOWN = 0)
}

frac_hmc <- process_fracs(frac_hmc_raw, "_num_hmc_region_frac")
frac_mc  <- process_fracs(frac_mc_raw,  "_num_mc_region_frac")

# 4. Calculate Patient-Wise Trends (UP/DOWN counts)
for (i in seq_along(patient_comparisons$control)) {
  ctrl_id <- patient_comparisons$control[i]
  expt_id <- patient_comparisons$expt[i]
  
  # Increment counters if Experiment > Control (UP) or Experiment < Control (DOWN)
  frac_hmc <- frac_hmc %>%
    dplyr::mutate(n_UP   = n_UP   + as.integer(.data[[expt_id]] > .data[[ctrl_id]]),
                  n_DOWN = n_DOWN + as.integer(.data[[expt_id]] < .data[[ctrl_id]]))
  
  frac_mc <- frac_mc %>%
    dplyr::mutate(n_UP   = n_UP   + as.integer(.data[[expt_id]] > .data[[ctrl_id]]),
                  n_DOWN = n_DOWN + as.integer(.data[[expt_id]] < .data[[ctrl_id]]))
}

# 5. Reformat and filter results
finalize_dmr_results <- function(frac_df, raw_df, tumor_df) {
  
  frac_df %>%
    # 1. Filter for regions with a consistent trend and valid names
    dplyr::filter(n_UP != n_DOWN, !is.na(Name)) %>%
    
    # 2. Join with significant DMR statistics
    dplyr::inner_join(raw_df, by = c("Chromosome", "Start", "End", "Name", "Annotation")) %>%
    
    # 3. Add the 'Type' (Tumor only vs Normal) using the key
    dplyr::mutate(
      key = paste(Chromosome, Start, End, Name, Annotation, sep = "_"),
      Type = ifelse(key %in% tumor_df$key, "Tumor only", "Normal"),
      n_Diff = n_UP - n_DOWN
    ) %>%
    
    # 4. Select and rename columns for the final report
    dplyr::select(
      Chromosome, Start, End, Name, Annotation, Type, 
      n_UP, n_DOWN, n_Diff, 
      n_CpGs = num_contexts, 
      Log2FC = mod_lfc,    # Standardized mapping
      FDR    = dmr_qvalue, 
      everything()
    ) %>%
    
    # Remove the temporary key before returning
    dplyr::select(-key)
}

sig_hmc   <- read.xlsx(file.path(path, "Tumor+Normal_5hmC_DMRs.xlsx"))
sig_mc    <- read.xlsx(file.path(path, "Tumor+Normal_5mC_DMRs.xlsx"))
tumor_hmc <- read.xlsx(file.path(path, "Tumor_Specific_5hmC_DMRs.xlsx"))
tumor_mc  <- read.xlsx(file.path(path, "Tumor_Specific_5mC_DMRs.xlsx"))

final_hmc_out <- finalize_dmr_results(frac_hmc, sig_hmc, tumor_hmc)
final_mc_out  <- finalize_dmr_results(frac_mc, sig_mc, tumor_mc)

# 6. Write to Excel
# Create a named list where names = Sheet Names and values = DataFrames
output_list <- list("5hmC_Trends" = final_hmc_out, 
                    "5mC_Trends"  = final_mc_out)

# Write everything in one go
write.xlsx(output_list, file = file.path(path, "Biomodal_Patient_Trend_Analysis.xlsx"), overwrite = TRUE)


# ---- Overlap with Felix Fang data ----
path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Biomodal"
ff_data <- read.xlsx(file.path(path, "can-24-0890_supplementary_tables_suppst.xlsx"),
                     sheet = "supplementary_table_s3",
                     startRow = 2) %>%
  dplyr::mutate(bin = case_when(feature %in% c("Promoter") ~ "Promoter",
                                feature %in% c("5' UTR", "Exon", "Intron", "3' UTR") ~ "GeneBody",
                                TRUE ~ "Other"))

hmc_up <- sig_dmr_hmc %>%
  dplyr::mutate(bin = case_when(Annotation %in% c("Promoter", "TSS") ~ "Promoter",
                                Annotation %in% c("Gene") ~ "GeneBody",
                                TRUE ~ "Other")) %>%
  dplyr::filter(mod_lfc > 0)

mc_up <- sig_dmr_mc %>%
  dplyr::mutate(bin = case_when(Annotation %in% c("Promoter", "TSS") ~ "Promoter",
                                Annotation %in% c("Gene") ~ "GeneBody",
                                TRUE ~ "Other")) %>%
  dplyr::filter(mod_lfc > 0)

hmc_down <- sig_dmr_hmc %>%
  dplyr::mutate(bin = case_when(Annotation %in% c("Promoter", "TSS") ~ "Promoter",
                                Annotation %in% c("Gene") ~ "GeneBody",
                                TRUE ~ "Other")) %>%
  dplyr::filter(mod_lfc < 0)

mc_down <- sig_dmr_mc %>%
  dplyr::mutate(bin = case_when(Annotation %in% c("Promoter", "TSS") ~ "Promoter",
                                Annotation %in% c("Gene") ~ "GeneBody",
                                TRUE ~ "Other")) %>%
  dplyr::filter(mod_lfc < 0)

genes_up_5hmc <- intersect(ff_data$gene, hmc_up$Name)
genes_down_5hmc <- intersect(ff_data$gene, hmc_down$Name)

genes_up_5mc <- intersect(ff_data$gene, mc_up$Name)
genes_down_5mc <- intersect(ff_data$gene, mc_down$Name)

genes_up <- union(genes_up_5hmc, genes_down_5mc)
genes_down <- union(genes_down_5hmc, genes_up_5mc)

common <- intersect(genes_up, genes_down)
uniq_genes_up <- setdiff(genes_up, common)
uniq_genes_down <- setdiff(genes_down, common)

gmt_dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GSEA_genesets"
species <- "Homo sapiens"
gmt_files <- list.files(file.path(gmt_dir, species), full.names = TRUE)
# Initialize result dataframes 
fgsea_df         <- data.frame()
gsea_df          <- data.frame()
ora_df_up        <- data.frame()
ora_df_down      <- data.frame()
concise_fgsea_df <- data.frame()
for (gmt_file in gmt_files) {
  
  # Format gene sets for fgsea and keep only genes present in ranked_list
  gmt <- fgsea::gmtPathways(gmt_file)
 
  # Format gene sets for clusterProfiler and keep only genes present in ranked_list
  pathway_gene_df <- utils::stack(x = gmt)
  colnames(pathway_gene_df) <- c("genes", "pathways")
  pathway_gene_df <- pathway_gene_df[, c("pathways", "genes")] # Reorder
  
  ora_res_up <- clusterProfiler::enricher(gene          = uniq_genes_up,
                                          universe      = NULL,
                                          TERM2GENE     = pathway_gene_df,
                                          minGSSize     = 1, #reduced from 15
                                          maxGSSize     = 500,
                                          pvalueCutoff  = 0.05,
                                          pAdjustMethod = "BH",
                                          qvalueCutoff  = 0.2)
  
  ora_res_down <- clusterProfiler::enricher(gene          = uniq_genes_down,
                                          universe      = NULL,
                                          TERM2GENE     = pathway_gene_df,
                                          minGSSize     = 1, #reduced from 15
                                          maxGSSize     = 500,
                                          pvalueCutoff  = 0.05,
                                          pAdjustMethod = "BH",
                                          qvalueCutoff  = 0.2)
  
  if (!is.null(ora_res_up))        { ora_df_up        <- dplyr::bind_rows(ora_df_up,        ora_res_up@result) }
  if (!is.null(ora_res_down))      { ora_df_down      <- dplyr::bind_rows(ora_df_down,      ora_res_down@result) }
}

# Rename columns consistently across different methods
lookup <- c(pathway = "ID", 
            geneID = "leadingEdge", geneID = "core_enrichment", 
            K = "size", K = "setSize", 
            padj = "p.adjust", 
            pval = "pvalue")

# Put your data frames in a named list
dfs <- list(fgsea_df      = fgsea_df,
            gsea_df       = gsea_df,
            ora_df_up     = ora_df_up,
            ora_df_down   = ora_df_down)
df_names <- names(dfs)

dfs <- lapply(X = df_names, FUN = function(df_name) {
  
  # Extract specific df
  df <- dfs[[df_name]]
  if (nrow(df) == 0) return(df)  # skip empty data frames
  
  # Add Direction column
  if (df_name == "ora_df_up"){
    df$Direction <- "Upregulated"
  } else if (df_name == "ora_df_down"){
    df$Direction <- "Downregulated"
  } else {
    df <- df %>%
      dplyr::mutate(Direction = dplyr::case_when(NES > 0 ~ "Upregulated",
                                                 NES < 0 ~ "Downregulated",
                                                 TRUE    ~ "No change"))
  }
  
  # ORA df specific formatting
  # k <- # overlap between pathway and input (sig_genes_up/sig_genes_down/gmt/ranked_list)
  # n <- # overlap between collection and input (sig_genes_up/sig_genes_down/gmt/ranked_list)
  # K <- # overlap between pathway and universe
  # N <- # overlap between collection and universe
  if (df_name %in% c("ora_df_up", "ora_df_down")){
    df <- df %>%
      tidyr::separate(col = GeneRatio, into = c("k", "n")) %>%
      tidyr::separate(col = BgRatio,   into = c("K", "N")) %>%
      dplyr::mutate(across(.cols = c(k, n, K, N), .fns = as.numeric)) %>%
      dplyr::mutate(GeneRatio       = k / n,
                    BackgroundRatio = K / N,
                    EnrichmentRatio = GeneRatio / BackgroundRatio,
                    combined_score  = GeneRatio * -log10(p.adjust),
                    NES             = NA_integer_)
  }
  
  # Standardize column names based on look up table
  df <- df %>%
    dplyr::rename(any_of(lookup))
  
  # Format results
  # IMPORTANT: gsea_df and ora_df store geneID as string "CXCL11/CCL2/..."
  # fgsea_df stores geneID as list of vectors which we convert to string "CXCL11/CCL2/..."
  df <- df %>%
    #dplyr::filter(padj <= 0.05) %>%
    tibble::remove_rownames() %>%
    tidyr::separate(col = pathway, into = c("Collection", "Description"), sep = "_", extra = "merge") %>%
    dplyr::mutate(Description = base::gsub(pattern = "_", replacement = " ", x = Description),
                  geneID = base::sapply(X = geneID, FUN = paste, collapse = "/")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(leading_edge_size = length(unlist(stringr::str_split(geneID, "/")))) %>%
    dplyr::ungroup() %>%
    as.data.frame() %>%
    dplyr::select(Collection, Description, leading_edge_size, K, padj, NES, Direction, everything(), -geneID, geneID)
  
  # Generate separate column for each gene in enriched pathways
  max_len <- max(df$leading_edge_size, na.rm = TRUE)
  if (is.finite(max_len) & max_len > 0) {
    df <- df %>%
      tidyr::separate(col = geneID, into = paste0("gene", 1:max_len), sep = "/", remove = TRUE, fill = "right")
  }
  
  return(df)
})

# Restore names
dfs <- stats::setNames(dfs, df_names)
ora_df   <- dplyr::bind_rows(dfs$ora_df_up, dfs$ora_df_down)

consensus_df <- dplyr::bind_rows(fgsea_df %>% dplyr::mutate(method = "FGSEA"), 
                                 gsea_df %>% dplyr::mutate(method = "GSEA"), 
                                 ora_df %>% dplyr::mutate(method = "ORA")) %>%
  dplyr::add_count(Collection, Description, Direction, name = "n_methods") %>%
  #dplyr::filter(n_methods > 1) %>%
  dplyr::mutate(Consensus = Direction) %>%
  dplyr::arrange(Collection, Description, desc(NES)) %>%
  dplyr::select(n_methods, method, Consensus, Collection, Description, 
                leading_edge_size, K, padj, NES, Direction, everything(), 
                -starts_with("gene",ignore.case = FALSE),
                starts_with("gene", ignore.case = FALSE)) 

# Identify top 10 Up & top 10 Down ORA pathways for each collection
top_ora <- consensus_df %>%
  dplyr::filter(method == "ORA", pval <= 0.05) %>%
  # Rank based on padj for each direction
  dplyr::group_by(Collection, Consensus) %>%
  dplyr::slice_min(order_by = padj, n = 10, with_ties = FALSE) %>%
  dplyr::ungroup()

pathway_results_list <- list(ora       = top_ora)

# Save as excel
file_name <- file.path(path, "Pathway_results.xlsx")
wb <- openxlsx::createWorkbook()
for (i in seq_along(pathway_results_list)) {
  openxlsx::addWorksheet(wb, sheetName = names(pathway_results_list)[i])
  openxlsx::writeData(wb, sheet = names(pathway_results_list)[i], x = pathway_results_list[[i]], rowNames = FALSE)
}
openxlsx::saveWorkbook(wb, file_name, overwrite = TRUE)


plot_pathway <- function(pathway_df, expr_mat, metadata, method, output_dir){
  
  
  # pathway_df$Collection, pathway_df$leading_edge_size,
  # pathway_df$Direction, pathway_df$padj, pathway_df$Description
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(expr_mat = expr_mat, metadata = metadata, output_dir = output_dir)  
  
  if (is.null(method) || !method %in% c("GSEA", "ORA")) {
    log_error(sample = "",
              step   = "plot_pathway",
              msg    = glue::glue("`method` '{method}' is invalid. Must be 'GSEA' or 'ORA'."))
  }
  
  # ---- 📊 Define Plotting Mapping Logic ----
  
  # Map method -> plotting aesthetics
  x_axis_map    <- c(ORA = "GeneRatio", GSEA = "NES")
  x_label_map   <- c(ORA = "Gene Ratio", GSEA = "Normalized Enrichment Score (NES)")
  
  x_axis  <- x_axis_map[[method]]
  x_label <- x_label_map[[method]]
  
  size_col      <- "leading_edge_size"
  color_col     <- "Direction"
  alpha_col     <- "pval"
  
  plot_colors <- c("Upregulated" = "#E69F00", "Downregulated" = "#56B4E9")
  
  # Pre-process descriptions for better wrapping in plots
  pathway_df <- pathway_df %>%
    dplyr::mutate(Description = base::gsub(pattern = "_", replacement = " ", x = Description),
                  Description = stringr::str_wrap(string = Description, width = 30))
  
  collections   <- unique(pathway_df$Collection)
  dot_plots     <- list()
  bar_plots     <- list()
  
  # ---- 🔄 Iterate Through Pathway Collections ----
  
  for (collection in collections) {
    
    # Subset and rank pathways based on score_var
    plot_df <- pathway_df %>% 
      dplyr::filter(Collection == collection) %>% 
      dplyr::filter(!is.na(.data[[x_axis]])) %>%
      dplyr::arrange(dplyr::desc(.data[[x_axis]]))
    
    if (nrow(plot_df) == 0) next
    
    # Pad with empty rows if fewer than 20 pathways for consistent plot scaling
    n_missing <- 20 - nrow(plot_df)
    if (n_missing > 0) {
      plot_df <- dplyr::bind_rows(plot_df, 
                                  data.frame(Description = as.character(base::seq_len(n_missing))))
    }
    
    # Dynamic Y-axis text sizing
    max_label_len <- base::max(base::nchar(plot_df$Description), na.rm = TRUE)
    y_text_size <- dplyr::case_when(max_label_len > 50 ~ 6,
                                    max_label_len > 35 ~ 7,
                                    max_label_len > 25 ~ 8,
                                    TRUE ~ 10)
    
    # Calculate limits dynamically per collection
    x_min <- if (method == "GSEA") { base::pmin(0, floor(min(plot_df[[x_axis]], na.rm = TRUE))) } else  { 0 }
    x_limits <- c(x_min, NA)
    
    # ---- 📈 Generate Bar Plot ---- 
    
    bar_p <- ggplot2::ggplot(data = plot_df,
                             mapping = aes(x     = .data[[x_axis]],
                                           y     = stats::reorder(Description, .data[[x_axis]]),
                                           fill  = .data[[color_col]],
                                           alpha = -log10(.data[[alpha_col]]))) +
      ggplot2::geom_col(width = 0.75, na.rm = TRUE) +
      ggplot2::labs(x = x_label, 
                    y = "", 
                    title = base::paste("Top", collection, "Pathways"), 
                    fill = "Direction") +
      ggplot2::theme_classic() +
      custom_theme +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = y_text_size)) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::scale_x_continuous(limits = x_limits, expand = ggplot2::expansion(mult = c(0, 0.05))) +
      ggplot2::scale_alpha_continuous(range = c(0.5, 1)) +
      ggplot2::scale_fill_manual(values = plot_colors) +
      guides(fill  = guide_legend(override.aes = list(shape = 22, size = 6)),
             color = guide_legend(override.aes = list(shape = 22, size = 6)),
             alpha = guide_legend(override.aes = list(shape = 22, size = 6))) +
      ggplot2::geom_text(aes(label = .data[[size_col]]), x = 0, hjust = -0.5, size = 3, show.legend = FALSE)
    
    bar_plots[[collection]] <- bar_p
    
    # ---- 🟢 Generate Dot Plot ---- 
    
    #size_vals <- c(min(plot_df[[size_col]], na.rm = TRUE), max(plot_df[[size_col]], na.rm = TRUE))
    size_vals <- plot_df$leading_edge_size[!base::is.na(plot_df$leading_edge_size)]
    breaks <- as.vector(floor(stats::quantile(size_vals, na.rm = TRUE) / 10) * 10)
    
    dot_p <- ggplot2::ggplot(data = plot_df,
                             mapping = aes(x     = .data[[x_axis]],
                                           y     = stats::reorder(Description, .data[[x_axis]]),
                                           fill  = .data[[color_col]],
                                           alpha = -log10(.data[[alpha_col]]),
                                           color = .data[[color_col]],
                                           size  = .data[[size_col]])) +
      ggplot2::geom_point(na.rm = TRUE) +
      ggplot2::labs(x = x_label, 
                    y = "", 
                    title = paste("Top", collection, "Pathways"), 
                    color = "Direction", 
                    size = "Counts") +
      ggplot2::theme_classic() +
      custom_theme +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = y_text_size)) +
      ggplot2::coord_cartesian(clip = "off") + 
      ggplot2::scale_x_continuous(limits = x_limits, expand = ggplot2::expansion(mult = c(0, 0.05))) +
      ggplot2::scale_alpha_continuous(range = c(0.5, 1)) +
      ggplot2::scale_color_manual(values = plot_colors) +
      ggplot2::scale_fill_manual(values = plot_colors) +  # need for coloring the legend
      ggplot2::scale_size(breaks = unique(breaks)) 
    
    dot_plots[[collection]] <- dot_p
    
    # ---- 🔥 ️ Generate Heatmap Multi-page PDF ----
    
    if (!is.null(expr_mat)) {
      
      heatmap_plots <- list()
      pathways <- plot_df %>% 
        dplyr::filter(!is.na(Direction)) %>% 
        dplyr::pull(Description) %>% 
        base::unique()
      
      for (pathway in pathways) {
        
        # Extract genes for this specific pathway from the long-format dataframe
        plot_genes <- pathway_df %>%
          dplyr::filter(Collection == collection, Description == pathway) %>%
          tidyr::pivot_longer(cols = dplyr::matches("^(gene|Gene)[0-9]+$"), 
                              values_to = "gene") %>%
          dplyr::pull(gene) %>%
          base::trimws() %>%
          stats::na.omit() %>%
          .[. != ""] %>%          # remove empty strings
          base::unique() %>%
          base::intersect(rownames(expr_mat))
        
        # Skip plotting if less than 2 genes
        if (base::length(plot_genes) < 2) next
        
        plot_title <- stringr::str_wrap(string = pathway, width = 30)
        
        # Plot heatmap
        ph <- plot_heatmap(expr_mat            = expr_mat[plot_genes, ,drop = FALSE], 
                           label_genes         = NULL,
                           filename            = NULL,
                           output_dir          = NULL,
                           metadata_col        = metadata, 
                           metadata_row        = NULL,
                           col_annotations     = proj.params$heatmap$col_annotations,
                           row_annotations     = proj.params$heatmap$row_annotations,
                           col_gap_by          = proj.params$heatmap$col_gap_by,
                           row_gap_by          = proj.params$heatmap$row_gap_by,
                           col_cluster_by      = proj.params$heatmap$col_cluster_by,
                           row_cluster_by      = proj.params$heatmap$row_cluster_by,
                           plot_title          = plot_title,
                           heatmap_palette     = proj.params$heatmap$heatmap_palette,
                           annotation_palette  = proj.params$heatmap$annotation_palette,
                           border_color        = proj.params$heatmap$border_color,
                           force_log           = proj.params$heatmap$force_log,
                           show_expr_legend    = proj.params$heatmap$show_expr_legend,
                           save_plot           = FALSE,
                           save_matrix         = FALSE)
        
        heatmap_plots[[pathway]] <- ph$ph$gtable
      }
      
      # Save stored heatmaps as pdf
      if (length(heatmap_plots) > 0) {
        
        file_extension <- ".pdf"
        file_name <- file.path(output_dir, paste0("Heatmap_", collection, "_", method, file_extension))
        
        # Open multi-page PDF
        grDevices::cairo_pdf(filename = file_name, width = 8, height = 11.5, onefile = TRUE)  
        
        for (ht in heatmap_plots) {
          grid::grid.newpage()
          grid::grid.draw(ht)
        }
        grDevices::dev.off() 
      }
    }
  } 
  
  # ---- 💾 Save Consolidated Summary Plots ----
  
  summary_plots <- base::list(Bar = bar_plots, Dot = dot_plots)
  
  for (type in names(summary_plots)) {
    
    file_extension <- ".pdf"
    file_name <- file.path(output_dir, paste0(type, "_plot_pathways_", method, file_extension))
    
    ggplot2::ggsave(filename = file_name,
                    plot     = cowplot::plot_grid(plotlist = summary_plots[[type]], ncol = 3, align = "hv"),
                    device   = grDevices::cairo_pdf,
                    width    = 3 * 6, 
                    height   = ceiling(length(summary_plots[[type]]) / 3) * 6, 
                    units    = "in",
                    dpi      = 300,
                    bg       = "white")
  }
  
  # ---- 🪵 Log Output and Return ----
  
  log_info(sample = "", 
           step   = "plot_pathway", 
           msg    = glue::glue("Successfully generated {method} visualizations in {output_dir}"))
  
  return(invisible(NULL))
}

plot_pathway(pathway_df = top_ora, 
             method     = "ORA",
             expr_mat   = NULL, 
             metadata   = NULL,
             output_dir = "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop")








# dmr_summary <- df %>%
#   dplyr::mutate(Direction = dplyr::case_when(mod_difference > 0 ~ "Up in Carotuximab",
#                                              mod_difference < 0 ~ "Down in Carotuximab",
#                                              TRUE ~ "No Change")) %>% # Should be rare after filtering
#   # Count unique genes by Annotation and Direction
#   dplyr::group_by(Annotation, Direction) %>%
#   dplyr::summarise(Unique_Genes = n_distinct(Name),
#                    .groups = "drop") %>%
#   # Pivot wider for a clean summary table format
#   tidyr::pivot_wider(names_from = Direction,
#                      values_from = Unique_Genes,
#                      values_fill = 0)
# 
# # Display the summary table
# print(dmr_summary)
# 
# #  Comparison directionof 5mc and 5hmc
# comparison_df <- dplyr::inner_join(x = tumor_DMR_hmc %>% dplyr::rename(hmc_diff = mod_difference, hmc_logFC = mod_fold_change) %>% dplyr::select(Name, Annotation, hmc_diff, hmc_logFC ),
#                                    y = tumor_DMR_mc %>% dplyr::rename(mc_diff = mod_difference, mc_logFC = mod_fold_change) %>% dplyr::select(Name, Annotation, mc_diff, mc_logFC ),
#                                    by = c("Name"= "Name", "Annotation"="Annotation")) %>%
#   dplyr::mutate(Trend = dplyr::case_when((mc_diff > 0 & hmc_diff < 0) | (mc_diff < 0 & hmc_diff > 0) ~ "Inverse",
#                                          (mc_diff > 0 & hmc_diff > 0) | (mc_diff < 0 & hmc_diff < 0) ~ "Concordant",
#                                          TRUE ~ "Unknown")) %>%
#   dplyr::arrange(Trend, Annotation, Name) %>%
#   dplyr::select(Name, Annotation, hmc_diff, mc_diff, hmc_logFC, mc_logFC, Trend)
# 
# comparison_df %>% dplyr::count(Trend)


# ---- UMAP ----

metadata <- read.xlsx(file.path(path, "Metadata.xlsx"))

# Some genes like CD99, MKSS have more than 1 TSS. So, multiple entries exist.
# Collapse into single entry using max signal.
# logit transform the fraction data
hmc_mat <- frac_hmc[,4:ncol(frac_hmc)] %>% 
  dplyr::group_by(Name, Annotation) %>%
  dplyr::summarize(across(.cols = everything(), .fns = max), .groups = "drop") %>%
  dplyr::inner_join(sig_dmr_hmc %>% dplyr::select(Name, Annotation),
                    by = c("Name", "Annotation")) %>%
  dplyr::mutate(ID = paste0(Name, ".", Annotation)) %>%
  dplyr::select(-Name, -Annotation) %>%
  tibble::column_to_rownames("ID") %>%
  dplyr::rename_with(.fn = function(x){ gsub("_num_hmc_region_frac", "", x)})

# Replace NAs and 0/1 extremes and logit transform
epsilon <- 1e-6
hmc_mat[is.na(hmc_mat)] <- epsilon
hmc_mat[hmc_mat == 0] <- epsilon
hmc_mat[hmc_mat == 1] <- 1 - epsilon
logit_hmc <- log(hmc_mat / (1 - hmc_mat)) %>% as.matrix()

mc_mat <- frac_mc[,4:ncol(frac_mc)] %>% 
  dplyr::group_by(Name, Annotation) %>%
  dplyr::summarize(across(.cols = everything(), .fns = max), .groups = "drop") %>%
  dplyr::inner_join(sig_dmr_mc %>% dplyr::select(Name, Annotation),
                    by = c("Name", "Annotation")) %>%
  dplyr::mutate(ID = paste0(Name, ".", Annotation)) %>%
  dplyr::select(-Name, -Annotation) %>%
  tibble::column_to_rownames("ID") %>%
  dplyr::rename_with(.fn = function(x){ gsub("_num_mc_region_frac", "", x)})

# Replace NAs and 0/1 extremes and logit transform
epsilon <- 1e-6
mc_mat[is.na(mc_mat)] <- epsilon
mc_mat[mc_mat == 0] <- epsilon
mc_mat[mc_mat == 1] <- 1 - epsilon
logit_mc <- log(mc_mat / (1 - mc_mat)) %>% as.matrix()


plot_pca(expr_mat = logit_hmc, txi = NULL, metadata = metadata, top_n_genes = 5000, skip_plot = FALSE, filename = "PCA_hmc", output_dir = path)
plot_pca(expr_mat = logit_mc, txi = NULL, metadata = metadata, top_n_genes = 5000, skip_plot = FALSE, filename = "PCA_mc", output_dir = path)

filename <- "UMAP_hmc"
plot_umap(expr_mat = logit_hmc, metadata, n_pcs = 50, n_neighbors = NULL, filename, output_dir = path)
filename <- "UMAP_mc"
plot_umap(expr_mat = logit_mc, metadata, n_pcs = 50, n_neighbors = NULL, filename, output_dir = path)

# ---- Pie Chart ----

sig_dmr_hmc_all <- dmr_hmc_21 %>% 
  dplyr::mutate(mean_mod_group_1 = ifelse(mean_mod_group_1 == 0, 1e-6, mean_mod_group_1),
                mean_mod_group_2 = ifelse(mean_mod_group_2 == 0, 1e-6, mean_mod_group_2),
                mod_fold_change = mean_mod_group_2 / mean_mod_group_1,
                mod_fold_change = log2(mod_fold_change)) %>%
  dplyr::filter(dmr_qvalue <= 0.05, num_contexts >= 10, abs(mod_fold_change) >= 0.58)

sig_dmr_mc_all <- dmr_mc_21 %>% 
  dplyr::mutate(mean_mod_group_1 = ifelse(mean_mod_group_1 == 0, 1e-6, mean_mod_group_1),
                mean_mod_group_2 = ifelse(mean_mod_group_2 == 0, 1e-6, mean_mod_group_2),
                mod_fold_change = mean_mod_group_2 / mean_mod_group_1,
                mod_fold_change = log2(mod_fold_change)) %>% 
  dplyr::filter(dmr_qvalue <= 0.05, num_contexts >= 10, abs(mod_fold_change) >= 0.58)


plot_piechart(metadata = sig_dmr_hmc_all, segment_col = "Annotation", filename = "5hmc_all", output_dir = path, split_col = NULL)
plot_piechart(metadata = sig_dmr_mc_all, segment_col = "Annotation", filename = "5mc_all", output_dir = path, split_col = NULL)


