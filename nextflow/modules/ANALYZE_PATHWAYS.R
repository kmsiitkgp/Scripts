#!/usr/bin/env Rscript

# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)      # Provides the %>% operator
  library(readr)
  library(tibble)
  library(openxlsx)
  library(fgsea)
  library(tidyr)
})

# ---- 🛠️ 2. Smart Setup (Find & source UTILS.R) ----

get_utils_path <- function() {
  # 1. Windows dev machine
  if (.Platform$OS.type == "windows") {
    return("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/Scripts/nextflow/modules/UTILS.R")
  }
  
  # 2. Interactive Linux / macOS (HPC interactive session)
  if (interactive()) {
    # Assume project root is current working directory
    return(file.path(getwd(), "modules", "UTILS.R"))
  }
  
  # 3. Non-interactive (Nextflow / Rscript)
  initial.options <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", initial.options, value = TRUE)
  if (length(file_arg) == 0) stop("Cannot detect script path for UTILS.R!")
  script_dir <- dirname(sub("--file=", "", file_arg))
  return(file.path(script_dir, "UTILS.R"))
}

utils_path <- get_utils_path()
if (!file.exists(utils_path)) stop(paste("❌ UTILS.R not found at:", utils_path))
source(utils_path)

# ---- 🧬 3. Function Definition (Always Runs) ----

# ---- run_gsea() ----
#
# Pure R function - no file I/O
# Loops over GMT files internally, applies per-collection padj correction
#
# Args:
#   ranked_list : Named numeric vector. Names = gene symbols, values = LFC.
#                 Must be sorted descending (high LFC to low LFC).
#   gmt_files   : Character vector of full paths to GMT files
#   minsize     : Minimum pathway size to test (default 15)
#   maxsize     : Maximum pathway size to test (default 500)
#
# Returns: Single formatted dataframe with GSEA results across all collections

run_gsea <- function(ranked_list, gmt_files, minsize = 15, maxsize = 500) {
  
  set.seed(1234)
  
  results <- list()
  
  # Determine score type based on range of ranked list
  # "pos" if all positive, "neg" if all negative, "std" if mixed
  score_type <- dplyr::case_when(min(ranked_list) >= 0 ~ "pos",
                                 max(ranked_list) <= 0 ~ "neg",
                                 TRUE                  ~ "std")
  
  for (gmt_file in gmt_files) {
    
    # Load GMT file into named list of pathway -> gene vectors
    gmt <- fgsea::gmtPathways(gmt_file)
    
    # Filter GMT to only genes present in ranked list
    # IMPORTANT: genes not in ranked list were never measured.
    # Keeping unmeasured genes inflates pathway size and distorts statistics.
    gmt <- lapply(X = gmt, FUN = intersect, y = names(ranked_list))
    
    # Run fgseaMultilevel - adaptive multilevel splitting for accurate p-values
    # nPermSimple = number of permutations for simple estimation (fallback)
    fgsea_res <- fgsea::fgseaMultilevel(
      pathways    = gmt,
      stats       = ranked_list,
      scoreType   = score_type,
      minSize     = minsize,
      maxSize     = maxsize,
      nPermSimple = 10000
    )
    
    # NOTE: collapsePathways() to remove redundant pathways is commented out.
    # Uncomment if needed to reduce output to non-redundant pathways only.
    # collapsed <- fgsea::collapsePathways(fgsea_res, gmt, ranked_list)
    # fgsea_res <- fgsea_res %>% dplyr::filter(pathway %in% collapsed$mainPathways)
    
    results[[gmt_file]] <- fgsea_res
  }
  
  # Combine results across all GMT collections into single dataframe
  combined <- dplyr::bind_rows(results)
  if (is.null(combined) || nrow(combined) == 0) return(NULL)
  
  # ---- Format Results ----
  
  combined <- combined %>%
    
    # Extract Collection and Description from pathway name
    # MSigDB pathway names follow format: COLLECTION_PATHWAY_NAME
    # e.g. HALLMARK_TNFA_SIGNALING_VIA_NFKB -> Collection = HALLMARK
    #                                        -> Description = TNFA SIGNALING VIA NFKB
    tidyr::separate(col   = pathway,
                    into  = c("Collection", "Description"),
                    sep   = "_",
                    extra = "merge") %>%
    dplyr::mutate(
      Description = gsub("_", " ", Description),                      # Replace underscores with spaces for readability
      Direction   = ifelse(NES > 0, "Upregulated", "Downregulated"),  # Direction based on NES sign
      K           = size,                                             # K = pathway size after filtering to measured genes
      N           = length(ranked_list),                              # N = total genes in ranked list (background for GSEA)
      leadingEdge = sapply(leadingEdge, paste, collapse = "/"),       # Collapse leading edge genes to / separated string
      leading_edge_size = lengths(strsplit(leadingEdge, "/"))) %>%
    dplyr::select(Collection, Description, Direction,pval, padj, NES, ES,  # Select and order final columns
                  K, N, log2err,leading_edge_size, leadingEdge)
  
  return(combined)
}


# ---- run_ora() ----
#
# Pure R function - no file I/O
# Loops over GMT files internally, applies per-collection padj correction
# Handles up/down/all gene lists
#
# Args:
#   gene_lists  : Named list with slots: up, down, all (any can be NULL)
#                 up   = upregulated DEGs
#                 down = downregulated DEGs
#                 all  = all DEGs (used when no LFC direction available)
#   universe    : Character vector of background genes.
#                 NULL = derive universe from each GMT collection (per-collection).
#                 NOTE: fora() does NOT accept NULL — handled internally here.
#   gmt_files   : Character vector of full paths to GMT files
#   minsize     : Minimum pathway size to test (default 15)
#   maxsize     : Maximum pathway size to test (default 500)
#
# Returns: Single formatted dataframe with ORA results across all collections

run_ora <- function(gene_lists, universe = NULL, gmt_files, minsize = 15, maxsize = 500) {
  
  # Remove NULL or empty slots from gene_lists before running
  gene_lists <- gene_lists[lengths(gene_lists) > 0]
  
  if (length(gene_lists) == 0) {
    warning("run_ora: all gene lists are NULL or empty. Returning NULL.")
    return(NULL)
  }
  
  results <- list()
  
  for (gmt_file in gmt_files) {
    
    # Load GMT file into named list of pathway -> gene vectors
    gmt <- fgsea::gmtPathways(gmt_file)
    
    # Determine universe for this collection
    # If universe is NULL (no experimental background available),
    # derive it from the GMT collection itself (per-collection universe).
    # NOTE: fora() does NOT accept NULL for universe unlike clusterProfiler::enricher()
    if (is.null(universe)) {
      
      # Per-collection universe = all unique genes in this GMT file
      # N will therefore vary per collection - this is scientifically correct.
      # See Methods sheet N column definition for full explanation.
      universe_to_use <- unique(unlist(gmt))
      universe_source <- "GMT-derived"
      
    } else {
      
      # Use experimental background (all measured genes)
      # Filter GMT to only genes present in universe
      # IMPORTANT: genes not measured in the experiment should not count
      # toward pathway size - keeping them inflates K and distorts BackgroundRatio
      gmt             <- lapply(X = gmt, FUN = intersect, y = universe)
      universe_to_use <- universe
      universe_source <- "Experimental"
    }
    
    # Run ORA for each gene list direction (up, down, or all)
    for (type in names(gene_lists)) {
      
      genes <- gene_lists[[type]]
      
      res <- fgsea::fora(pathways = gmt,
                         genes    = genes,
                         universe = universe_to_use,
                         minSize  = minsize,
                         maxSize  = maxsize)
      
      if (is.null(res) || nrow(res) == 0) next
      
      # Add metadata columns needed for downstream formatting
      res$type            <- type
      res$n               <- length(genes)            # n = number of input genes tested for this direction
      res$N               <- length(universe_to_use)  # N = total universe size for this collection
      res$universe_source <- universe_source
      
      results[[paste(gmt_file, type, sep = "_")]] <- res
    }
  }
  
  # Combine results across all GMT collections and directions
  combined <- dplyr::bind_rows(results)
  if (is.null(combined) || nrow(combined) == 0) return(NULL)
  
  # ---- Format Results ----
  
  combined <- combined %>%
    
    # Extract Collection and Description from pathway name
    # MSigDB pathway names follow format: COLLECTION_PATHWAY_NAME
    tidyr::separate(col   = pathway,
                    into  = c("Collection", "Description"),
                    sep   = "_",
                    extra = "merge") %>%
    dplyr::mutate(
      Description = gsub("_", " ", Description),                      # Replace underscores with spaces for readability
      # Direction label based on which gene list was used
      Direction = dplyr::case_when(type == "up"   ~ "Enriched in UP genes",
                                   type == "down" ~ "Enriched in DOWN genes",
                                   type == "all"  ~ "Enriched in ALL input genes"),
      k = as.integer(overlap),    # k = number of input genes overlapping with the pathway
      K = size,                   # K = pathway size within universe after filtering if experimental background
      # Enrichment ratio calculations
      GeneRatio       = k / n,
      BackgroundRatio = K / N,
      EnrichmentRatio = GeneRatio / BackgroundRatio,
      combined_score  = GeneRatio * -log10(padj),
      overlapGenes = sapply(overlapGenes, paste, collapse = "/")) %>%   # Collapse overlap genes to / separated string
    dplyr::select(Collection, Description, Direction, pval, padj, EnrichmentRatio,   # Select and order final columns
                  k, n, K, N, GeneRatio, BackgroundRatio, combined_score, overlapGenes)
  
  return(combined)
}


# ---- analyze_pathways() ----
#
# Orchestrator function - handles file I/O, input detection, calls run_gsea/
# run_ora, writes Excel output with GSEA, ORA, and Methods sheets
#
# Args:
#   contrast       : Character. Name of the contrast (e.g. "TreatedVsControl").
#                    Used for output subfolder naming.
#   input          : One of:
#                    - File path (csv, xlsx, xls, rds, or plain text gene list)
#                    - data.frame with SYMBOL + padj + log2FoldChange (full DEG results)
#                    - data.frame with SYMBOL + log2FoldChange only (2-col ranked)
#                    - Character vector of gene symbols (gene list)
#   only_deg       : Logical. Only relevant for gene list or 2-col ranked input.
#                    TRUE  = input contains only DEGs (ORA only, GMT universe)
#                    FALSE = input contains all tested genes (GSEA valid)
#                    Ignored when input is a full DEG dataframe (auto-detected).
#   gmt_dir        : Path to directory containing GMT files (.gmt extension)
#   output_dir     : Path to output directory. Results saved in output_dir/contrast/
#   padj : Threshold for defining DEGs from full DEG df (default 0.05)
#   minsize        : Minimum pathway size to test (default 15)
#   maxsize        : Maximum pathway size to test (default 500)
#
# Returns: Invisible list with gsea, ora, and methods dataframes


analyze_pathways <- function(contrast, input_data, only_deg = FALSE, gmt_dir,
                             output_dir, padj_cutoff = 0.05, minsize = 15, maxsize = 500) {
  
  input_data <- load_smart(input_data)

  # ---- STEP 2: Input Type Detection ----
  # Detect which of the 3 supported input types was provided based on column structure
  
  has_full_deg_cols <- is.data.frame(input_data) &&
    all(c("SYMBOL", "padj", "log2FoldChange") %in% colnames(input_data))
  
  has_ranked_cols <- is.data.frame(input_data) &&
    all(c("SYMBOL", "log2FoldChange") %in% colnames(input_data)) &&
    !"padj" %in% colnames(input_data)
  
  is_gene_list <- is.character(input_data) ||
    (is.data.frame(input_data) && ncol(input_data) == 1)
  
  # Set input_type label used in Methods sheet
  if (has_full_deg_cols) {
    input_type <- "full_deg_df"
  } else if (has_ranked_cols) {
    input_type <- "ranked_2col"
  } else if (is_gene_list) {
    input_type <- "gene_list"
    # Extract gene symbols from single-column df or character vector
    if (is.data.frame(input_data)) input_data <- as.character(input_data[[1]])
  } else {
    log_error(sample = "", 
              step = "analyze_pathways",
              msg = glue::glue(
                "Could not detect input type. Provide either: ",
                "(1) df with SYMBOL, padj, log2FoldChange, ",
                "(2) df with SYMBOL, log2FoldChange, or ",
                "(3) a character vector of gene symbols."))
  }
  
  log_info(sample = "", step = "analyze_pathways",
           msg    = paste0("Input type detected: ", input_type))
  
  # ---- STEP 3: Validate GMT Directory ----
  
  if (!dir.exists(gmt_dir)) {
    log_error(sample = "", step = "analyze_pathways",
              msg    = paste0("gmt_dir '", gmt_dir, "' does not exist."))
  }
  
  gmt_files <- list.files(gmt_dir, full.names = TRUE, pattern = "\\.gmt$")
  
  if (length(gmt_files) == 0) {
    log_error(sample = "", step = "analyze_pathways",
              msg    = paste0("No .gmt files found in '", gmt_dir, "'."))
  }
  
  # ---- STEP 4: Build Analysis Inputs ----
  # Decide can_gsea, can_ora, ranked_list, gene_lists, universe
  # based on input type and only_deg flag.
  # Full decision logic is documented in the Methods sheet of the Excel output.
  
  can_gsea         <- FALSE
  can_ora          <- FALSE
  ranked_list      <- NULL
  gene_lists       <- list(up = NULL, down = NULL, all = NULL)
  universe         <- NULL
  universe_label   <- NULL
  gsea_skip_reason <- NULL
  ora_skip_reason  <- NULL
  
  # ---- Case 1: Full DEG dataframe (SYMBOL + padj + log2FoldChange) ----
  if (input_type == "full_deg_df") {
    
    # Clean, deduplicate and sort by LFC descending for GSEA
    ranked_df <- input_data %>%
      dplyr::filter(!is.na(SYMBOL), !is.na(log2FoldChange)) %>%
      dplyr::distinct(SYMBOL, .keep_all = TRUE) %>%
      dplyr::arrange(dplyr::desc(log2FoldChange))
    
    # Check if non-significant genes exist - required for valid GSEA background
    has_non_sig <- any(ranked_df$padj > padj_cutoff, na.rm = TRUE)
    
    # Build ranked list for GSEA (all genes sorted by LFC high to low)
    ranked_list <- stats::setNames(object = as.numeric(ranked_df$log2FoldChange),
                                   nm     = ranked_df$SYMBOL)
    
    # Build ORA gene lists - significant genes split by LFC direction
    sig_up   <- ranked_df %>%
      dplyr::filter(padj <= padj_cutoff, log2FoldChange > 0) %>%
      dplyr::pull(SYMBOL)
    sig_down <- ranked_df %>%
      dplyr::filter(padj <= padj_cutoff, log2FoldChange < 0) %>%
      dplyr::pull(SYMBOL)
    
    gene_lists <- list(up = sig_up, down = sig_down, all = NULL)

    if (has_non_sig) {
      
      # Full background available - both GSEA and ORA valid
      can_gsea       <- TRUE
      can_ora        <- TRUE
      universe       <- ranked_df$SYMBOL
      universe_label <- paste0(
        "Experimental background — all genes in input df (N=", length(universe), ")"
      )
      
    } else {
      
      # No non-significant genes found - sig-only input, GSEA not valid
      can_gsea         <- FALSE
      can_ora          <- TRUE
      gsea_skip_reason <- paste0(
        "All genes in input df have padj <= ", padj_cutoff, ". ",
        "No non-significant genes found. ",
        "GSEA requires a full ranked list including non-DEGs. ",
        "GMT-derived universe used for ORA.")
      # No experimental background available - fall back to GMT-derived universe
      universe       <- NULL
      universe_label <- "GMT-derived (per collection) — no experimental background available ⚠️"
    }
  
  # ---- Case 2: 2-column ranked dataframe (SYMBOL + log2FoldChange, no padj) ----
  } else if (input_type == "ranked_2col") {
    
    ranked_df <- input_data %>%
      dplyr::filter(!is.na(SYMBOL), !is.na(log2FoldChange)) %>%
      dplyr::distinct(SYMBOL, .keep_all = TRUE) %>%
      dplyr::arrange(dplyr::desc(log2FoldChange))
    
    ranked_list <- stats::setNames(object = as.numeric(ranked_df$log2FoldChange),
                                   nm     = ranked_df$SYMBOL)
    
    if (!only_deg) {
      
      # All tested genes are ranked -> GSEA valid
      # ORA not possible - no padj to define which genes are significant
      can_gsea        <- TRUE
      can_ora         <- FALSE
      universe        <- ranked_df$SYMBOL
      universe_label  <- paste0("All genes in ranked list (N=", length(universe), ")")
      ora_skip_reason <- "2-column ranked input with only_deg=FALSE. 
      No padj available to define DEGs for ORA."
      
    } else {
      
      # Ranked but flagged as DEGs only -> ORA only, split by LFC sign
      # GSEA not valid - ranked list contains only DEGs, not all tested genes
      can_gsea         <- FALSE
      can_ora          <- TRUE
      gsea_skip_reason <- "2-column ranked input flagged as DEG-only (only_deg=TRUE).
      GSEA requires all tested genes in the ranked list, not just DEGs."
      universe       <- NULL
      universe_label <- "GMT-derived (per collection) — no experimental background available ⚠️"
      
      sig_up   <- ranked_df %>% dplyr::filter(log2FoldChange > 0) %>% dplyr::pull(SYMBOL)
      sig_down <- ranked_df %>% dplyr::filter(log2FoldChange < 0) %>% dplyr::pull(SYMBOL)
      gene_lists <- list(up = sig_up, down = sig_down, all = NULL)
    }

  # ---- Case 3: Gene list (character vector or single-column df) ----
  } else if (input_type == "gene_list") {
    
    genes <- unique(na.omit(input_data))
    
    if (!only_deg) {
      
      # Not flagged as DEGs - cannot run any meaningful analysis
      can_gsea         <- FALSE
      can_ora          <- FALSE
      gsea_skip_reason <- "Gene list input with only_deg=FALSE. Cannot run GSEA without ranked values."
      ora_skip_reason  <- "Gene list input with only_deg=FALSE. Set only_deg=TRUE if these are DEGs."
      
      log_warn(sample = "", step = "analyze_pathways",
               msg    = "No analysis performed. Gene list provided but only_deg=FALSE.
        Set only_deg=TRUE if your input is a list of DEGs.")
      
    } else {
      
      # DEG gene list - ORA only, no direction (no LFC available)
      can_gsea         <- FALSE
      can_ora          <- TRUE
      gsea_skip_reason <- "Gene list input. GSEA requires a ranked list of all tested genes."
      universe         <- NULL
      universe_label   <- "GMT-derived (per collection) — no experimental background available ⚠️"
      gene_lists       <- list(up = NULL, down = NULL, all = genes)
    }
  } 
  
  # ---- STEP 5: Run Analyses ----
  
  gsea_results <- NULL
  ora_results  <- NULL

  if (can_gsea) {
    log_info(sample = "", step = "analyze_pathways", msg = "Running GSEA...")
    gsea_results <- run_gsea(
      ranked_list = ranked_list,
      gmt_files   = gmt_files,
      minsize     = minsize,
      maxsize     = maxsize
    )
  }
  
  if (can_ora) {
    log_info(sample = "", step = "analyze_pathways", msg = "Running ORA...")
    ora_results <- run_ora(
      gene_lists = gene_lists,
      universe   = universe,
      gmt_files  = gmt_files,
      minsize    = minsize,
      maxsize    = maxsize
    )
  }
  
  # ---- STEP 6: Build Methods Sheet Content ----
  # Three separate dataframes written to the Methods sheet with gap rows between them:
  #   6a. Logic table      - all possible input scenarios with Analyzed flag
  #   6b. Analysis summary - what actually ran in this specific analysis run
  #   6c. Column defs      - explanation of every column in GSEA and ORA sheets
  
  # ---- 6a. Logic Table ----
  # Full decision table showing all supported input scenarios.
  # The Analyzed column flags which row applies to the current analysis.
  
  methods_table <- data.frame(
    Input      = c(
      "Full DEG df — has non-sig genes",
      "Full DEG df — sig-only (no non-sig rows)",
      "2-col ranked (SYMBOL + LFC)",
      "2-col ranked (SYMBOL + LFC)",
      "Gene list",
      "Gene list"
    ),
    only_deg   = c("ignored", "ignored", "FALSE", "TRUE", "TRUE", "FALSE"),
    GSEA       = c("YES", "NO",  "YES", "NO",  "NO",  "NO"),
    ORA        = c("YES (up/down)", "YES (up/down)", "NO", "YES (up/down)", "YES (all)", "NO"),
    Universe   = c(
      "All genes in df",
      "Per-collection GMT genes ⚠️",
      "All genes in ranked list",
      "Per-collection GMT genes ⚠️",
      "Per-collection GMT genes ⚠️",
      "N/A"
    ),
    Filter_GMT = c("YES", "NO", "YES", "NO", "NO", "NO"),
    Reason     = c(
      "Full experimental background available",
      "No background — only DEGs in input",
      "Full measured genes available — no padj for ORA",
      "DEGs only — no experimental background",
      "DEGs only — no background, no LFC for direction",
      "Nothing runs — set only_deg=TRUE if input is DEGs ⚠️"
    ),
    Analyzed   = c(
      ifelse(input_type == "full_deg_df" &&  can_gsea &&  can_ora, "YES <- Your analysis", ""),
      ifelse(input_type == "full_deg_df" && !can_gsea &&  can_ora, "YES <- Your analysis", ""),
      ifelse(input_type == "ranked_2col" && !only_deg,             "YES <- Your analysis", ""),
      ifelse(input_type == "ranked_2col" &&  only_deg,             "YES <- Your analysis", ""),
      ifelse(input_type == "gene_list"   &&  only_deg,             "YES <- Your analysis", ""),
      ifelse(input_type == "gene_list"   && !only_deg,             "YES <- Your analysis", "")
    ),
    stringsAsFactors = FALSE
  )
  
  # ---- 6b. Analysis Summary ----
  # Records exactly what was done in this specific run for full reproducibility
  
  summary_rows <- data.frame(
    Parameter = c(
      "=== ANALYSIS SUMMARY ===",
      "Contrast",
      "Input type detected",
      "only_deg parameter",
      "padj threshold",
      "GSEA performed",
      "GSEA skip reason",
      "ORA performed",
      "ORA skip reason",
      "Universe source",
      "Universe size",
      "Number of UP DEGs tested (ORA)",
      "Number of DOWN DEGs tested (ORA)",
      "Number of ALL DEGs tested (ORA)",
      "GMT directory",
      "GMT files used"
    ),
    Value = c(
      "",
      contrast,
      input_type,
      as.character(only_deg),
      as.character(padj_cutoff),
      ifelse(can_gsea, "YES", "NO"),
      ifelse(is.null(gsea_skip_reason), "N/A", gsea_skip_reason),
      ifelse(can_ora,  "YES", "NO"),
      ifelse(is.null(ora_skip_reason),  "N/A", ora_skip_reason),
      ifelse(is.null(universe_label),   "N/A", universe_label),
      ifelse(is.null(universe), "GMT-derived — varies per collection", as.character(length(universe))),
      ifelse(!is.null(gene_lists$up),   as.character(length(gene_lists$up)),   "N/A"),
      ifelse(!is.null(gene_lists$down), as.character(length(gene_lists$down)), "N/A"),
      ifelse(!is.null(gene_lists$all),  as.character(length(gene_lists$all)),  "N/A"),
      gmt_dir,
      paste(basename(gmt_files), collapse = ", ")
    ),
    stringsAsFactors = FALSE
  )
  
  # ---- 6c. Column Definitions ----
  # Explains every column in the GSEA and ORA output sheets
  
  col_defs <- data.frame(
    Column = c(
      "=== COLUMN DEFINITIONS ===",
      "--- Common Columns (GSEA and ORA) ---",
      "Collection",
      "Description",
      "Direction",
      "pval",
      "padj",
      "K",
      "N",
      "",
      "--- GSEA-Specific Columns ---",
      "NES",
      "ES",
      "log2err",
      "leading_edge_size",
      "leadingEdge",
      "",
      "--- ORA-Specific Columns ---",
      "k",
      "n",
      "GeneRatio",
      "BackgroundRatio",
      "EnrichmentRatio",
      "combined_score",
      "overlapGenes"
    ),
    Definition = c(
      "",
      "",
      "Gene set collection name extracted from pathway name (e.g. HALLMARK, KEGG, REACTOME)",
      "Pathway name with underscores replaced by spaces",
      "For GSEA: Upregulated (NES > 0) or Downregulated (NES < 0). For ORA: Enriched in UP, DOWN, or ALL input genes",
      "Nominal p-value from statistical test",
      "BH-adjusted p-value (FDR) corrected independently per GMT collection",
      "Number of genes in the pathway",
      paste0(
        "Number of genes in the background/universe. ",
        "Total genes in ranked list (GSEA) or total experimental background/GMT-derived genes (ORA). ",
        "When GMT-derived universe is used, N varies per collection as it reflects ",
        "the total unique genes in that specific GMT collection — see Collection column."
      ),
      "",
      "",
      "Normalized Enrichment Score. Positive = pathway upregulated, Negative = pathway downregulated",
      "Raw Enrichment Score before normalization across permutations",
      "Expected error of the p-value estimate. Lower value = more accurate p-value estimate",
      "Number of genes in the leading edge subset — the core genes driving the enrichment signal",
      "/-separated genes in the leading edge subset — the core genes driving enrichment",
      "",
      "",
      "Number of input genes overlapping with the pathway",
      "Number of input genes tested (UP, DOWN, or ALL depending on direction)",
      "k/n — proportion of input genes that belong to this pathway",
      "K/N — proportion of pathway genes in the background. Represents the expected baseline frequency of this pathway in the background",
      "GeneRatio/BackgroundRatio — how much more often this pathway appears in your input genes compared to what would be expected by chance",
      "GeneRatio x -log10(padj) — ranks pathways by both enrichment strength and statistical significance combined",
      "/-separated genes from input list that overlap with the pathway"
    ),
    stringsAsFactors = FALSE
  )
  
  # ---- STEP 7: Save Results to Excel ----
  # Methods sheet uses Option C: 3 separate dataframes written with startRow
  # to avoid column mismatch from binding dataframes of different widths.
  # Two empty rows separate each section for readability.
  
  # Sanitize contrast name for safe file/folder naming
  safe_contrast <- gsub("[^[:alnum:]_-]", "_", contrast)
  output_path   <- file.path(output_dir, safe_contrast)
  
  if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
  
  output_file <- file.path(output_path, "Pathways.xlsx")
  
  # Build workbook manually to control sheet order and Methods layout
  wb <- openxlsx::createWorkbook()
  
  # Add GSEA sheet if results exist
  if (!is.null(gsea_results) && nrow(gsea_results) > 0) {
    openxlsx::addWorksheet(wb, sheetName = "GSEA")
    openxlsx::writeData(wb, sheet = "GSEA", x = gsea_results, rowNames = FALSE)
  }
  
  # Add ORA sheet if results exist
  if (!is.null(ora_results) && nrow(ora_results) > 0) {
    openxlsx::addWorksheet(wb, sheetName = "ORA")
    openxlsx::writeData(wb, sheet = "ORA", x = ora_results, rowNames = FALSE)
  }
  
  # Add Methods sheet - always written regardless of whether analyses ran
  openxlsx::addWorksheet(wb, sheetName = "Methods")
  
  # Section 1: Logic table — starts at row 1
  # +1 for header row, +2 for gap rows between sections
  logic_start   <- 1
  openxlsx::writeData(wb, sheet = "Methods", x = methods_table,
                      startRow = logic_start, rowNames = FALSE)
  
  # Section 2: Analysis summary — starts after logic table + 2 empty gap rows
  summary_start <- logic_start + nrow(methods_table) + 1 + 2
  openxlsx::writeData(wb, sheet = "Methods", x = summary_rows,
                      startRow = summary_start, rowNames = FALSE)
  
  # Section 3: Column definitions — starts after summary + 2 empty gap rows
  coldefs_start <- summary_start + nrow(summary_rows) + 1 + 2
  openxlsx::writeData(wb, sheet = "Methods", x = col_defs,
                      startRow = coldefs_start, rowNames = FALSE)
  
  openxlsx::saveWorkbook(wb, file = output_file, overwrite = TRUE)
  
  # ---- STEP 8: Log and Return ----
  
  log_info(sample = "", step = "analyze_pathways",
           msg    = paste0("Pathway analysis complete. Results saved to '", output_file, "'."))
  
  return(invisible(list(gsea    = gsea_results,
                        ora     = ora_results,
                        methods = list(logic_table      = methods_table,
                                       analysis_summary = summary_rows,
                                       column_defs      = col_defs))))
}

# ---- 🚀 4. Smart Execution (Nextflow Only) ----

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  get_arg <- function(idx, default = NULL) {
    if (idx > length(args)) return(default)
    val <- args[idx]
    if (is.na(val) || val == "" || val == "null" || val == "NULL") return(default)
    return(val)
  }
  
  analyze_pathways(
    contrast    = get_arg(1),
    input_data  = get_arg(2),
    only_deg    = as.logical(toupper(get_arg(3, "FALSE"))),
    gmt_dir     = get_arg(4, "."),
    output_dir  = get_arg(5, "."),
    padj_cutoff = as.numeric(get_arg(6, 0.05)),
    minsize     = as.numeric(get_arg(7, 15)),
    maxsize     = as.numeric(get_arg(8, 500))
  ) 
}
