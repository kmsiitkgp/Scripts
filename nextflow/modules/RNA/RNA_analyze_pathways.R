#!/usr/bin/env Rscript

# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(openxlsx)
  library(fgsea)
  library(tidyr)
})

# ---- 🛠️ 2. Smart Setup (Find & source UTILS.R) ----

get_modules_dir <- function() {
  # 1. Windows dev machine
  if (.Platform$OS.type == "windows") {
    return("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/Scripts/nextflow/modules")
  }
  # 2. Interactive Linux / macOS (HPC interactive session)
  if (interactive()) {
    return(file.path(getwd(), "modules"))
  }
  # 3. Non-interactive (Nextflow / Rscript)
  initial.options <- commandArgs(trailingOnly = FALSE)
  file_arg        <- grep("--file=", initial.options, value = TRUE)
  if (length(file_arg) == 0) stop("Cannot detect modules directory path!")
  modules_dir      <- dirname(sub("--file=", "", file_arg))
  return(modules_dir)
}

find_utils_dir <- function(start_dir) {
  d <- normalizePath(start_dir, mustWork = FALSE)
  repeat {
    if (file.exists(file.path(d, "UTILS.R"))) return(d)
    parent <- dirname(d)
    if (parent == d) stop("❌ UTILS.R not found searching upward from: ", start_dir)
    d <- parent
  }
}

modules_dir <- get_modules_dir()
utils_path  <- file.path(find_utils_dir(modules_dir), "UTILS.R")
source(utils_path)

# ---- 🧬 3. Function Definition (Always Runs) ----

# ═══════════════════════════════════════════════════════════════════════════════
# run_gsea()
#
# Pure function — no file I/O. Runs fgsea across all GMT collections and
# returns a single formatted dataframe. padj is corrected independently per
# GMT collection (not pooled across all collections) so that dense collections
# like REACTOME don't dominate the correction and dilute smaller ones.
#
# Args:
#   ranked_list : Named numeric vector. Names = SYMBOL, values = ranking stat.
#                 Must be sorted descending (highest to lowest) before passing.
#   gmt_files   : Character vector of full paths to .gmt files.
#   minsize     : Min pathway size to test after GMT/ranked_list intersection (default 15).
#   maxsize     : Max pathway size to test (default 500).
#
# Returns: Single formatted dataframe with GSEA results across all collections,
#          or NULL if no results were produced.
# ═══════════════════════════════════════════════════════════════════════════════

run_gsea <- function(ranked_list, gmt_files, minsize = 15, maxsize = 500) {

  # fgseaMultilevel uses permutation-based null distributions — set seed for
  # reproducibility. nPermSimple permutations are used as a fallback estimator
  # when the multilevel adaptive splitting cannot achieve the target precision.
  set.seed(1234)

  results <- list()

  # ── Detect score type ────────────────────────────────────────────────────
  # fgsea's scoreType controls how the enrichment score is normalised:
  #   "std" : standard, bidirectional — ranked list has both positive and negative
  #           values (typical for LFC or t-statistic ranked lists)
  #   "pos" : all values non-negative — e.g. ranking by -log10(pvalue) only
  #   "neg" : all values non-positive
  # Using the wrong scoreType produces incorrect NES normalisation and distorted
  # p-values. Auto-detection here prevents the user from having to specify it.
  score_type <- dplyr::case_when(
    min(ranked_list) >= 0 ~ "pos",
    max(ranked_list) <= 0 ~ "neg",
    TRUE                  ~ "std"
  )

  for (gmt_file in gmt_files) {

    gmt <- fgsea::gmtPathways(gmt_file)

    # Determine the expected prefix based on filename (Human or Mouse)
    file_name <- tolower(basename(gmt_file))

    target_prefix <- dplyr::case_when(
      # Hallmark: h.all (Hs) or mh.all (Mm)
      grepl("^h\\.all|^mh\\.all", file_name)         ~ "HALLMARK",

      # KEGG Pathways
      grepl("kegg", file_name)                       ~ "KEGG",

      # CGP: c2.cgp (Hs) or m2.cgp (Mm)
      grepl("cgp", file_name)                        ~ "CGP",

      # Canonical Pathways
      grepl("biocarta", file_name)                   ~ "BIOCARTA",
      grepl("reactome", file_name)                   ~ "REACTOME",
      grepl("wikipathways", file_name)               ~ "WP",

      # Gene Ontology: c5.go.bp (Hs) or m5.go.bp (Mm)
      grepl("go\\.bp", file_name)                    ~ "GOBP",
      grepl("go\\.cc", file_name)                    ~ "GOCC",
      grepl("go\\.mf", file_name)                    ~ "GOMF",

      # Oncogenic: c6.all (Hs) or m6.all (Mm)
      grepl("c6\\.all|m6\\.all", file_name)          ~ "C6",

      # Catch-all for any other MSigDB collections (C1, C3, C7, etc.)
      grepl("^[cm][0-9]\\.", file_name)              ~ toupper(sub("^([cm][0-9])\\..*", "\\1", file_name)),
      TRUE ~ "OTHER"
    )

    # Only add the prefix if it's NOT already at the start of the names
    if (target_prefix != "") {
      prefix_pattern <- paste0("^", target_prefix, "_")

      # Identify names that DON'T already start with the prefix
      needs_prefix <- !grepl(prefix_pattern, names(gmt))

      if (any(needs_prefix)) {
        names(gmt)[needs_prefix] <- paste0(target_prefix, "_", names(gmt)[needs_prefix])
      }
    }

    # Filter each pathway to only genes present in the ranked list.
    # Genes in the GMT but absent from the ranked list were not measured in
    # this experiment — keeping them inflates effective pathway size (K) and
    # makes enrichment appear weaker than it actually is. Intersecting first
    # ensures minSize/maxSize thresholds reflect only measurable genes.
    gmt <- lapply(X = gmt, FUN = intersect, y = names(ranked_list))

    # fgseaMultilevel uses adaptive multilevel splitting to estimate p-values
    # more accurately than simple permutation, especially for very small p-values
    # (< 1e-4) where simple permutation would need millions of permutations.
    # nPermSimple is the fallback permutation count used when multilevel
    # splitting fails to converge for a pathway — 10,000 is the recommended default.
    fgsea_res <- fgsea::fgseaMultilevel(
      pathways    = gmt,
      stats       = ranked_list,
      scoreType   = score_type,
      minSize     = minsize,
      maxSize     = maxsize,
      #BPPARAM   = BiocParallel::SerialParam(), #uncomment on parallel processing error
      nPermSimple = 10000
    )

    # NOTE: collapsePathways() is available to remove redundant overlapping
    # pathways but is not applied by default — it reduces output volume but
    # may hide biologically distinct pathways that happen to share genes.
    # Uncomment below if output is too large to interpret:
    # collapsed <- fgsea::collapsePathways(fgsea_res, gmt, ranked_list)
    # fgsea_res <- fgsea_res %>% dplyr::filter(pathway %in% collapsed$mainPathways)

    results[[gmt_file]] <- fgsea_res
  }

  combined <- dplyr::bind_rows(results)
  if (is.null(combined) || nrow(combined) == 0) return(NULL)

  # ── Format results ───────────────────────────────────────────────────────
  combined <- combined %>%

    # MSigDB pathway names follow the format: COLLECTION_PATHWAY_NAME
    # e.g. "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
    #       → Collection = "HALLMARK", Description = "TNFA SIGNALING VIA NFKB"
    # extra = "merge" is critical: pathway names contain multiple underscores,
    # so a simple split on "_" would produce many columns. "merge" collapses
    # everything after the first underscore into the Description column.
    tidyr::separate(col   = pathway,
                    into  = c("Collection", "Description"),
                    sep   = "_",
                    extra = "merge") %>%
    dplyr::mutate(
      Description       = gsub("_", " ", Description),
      Direction         = ifelse(NES > 0, "Upregulated", "Downregulated"),
      K                 = size,                                        # pathway size after intersection with ranked_list
      N                 = length(ranked_list),                         # total genes in ranked list (background)
      leadingEdge       = sapply(leadingEdge, paste, collapse = "/"),  # collapse list → /-separated string for Excel
      # Recompute leading_edge_size from the collapsed string rather than the
      # original list column — ensures consistency after paste() above
      leading_edge_size = lengths(strsplit(leadingEdge, "/"))
    ) %>%
    dplyr::select(Collection, Description, Direction, pval, padj, NES, ES,
                  K, N, log2err, leading_edge_size, leadingEdge) %>%
    dplyr::filter(leading_edge_size > 1)

  return(combined)
}


# ═══════════════════════════════════════════════════════════════════════════════
# run_ora()
#
# Pure function — no file I/O. Runs Fisher's exact test (via fgsea::fora) for
# over-representation of up/down/all gene lists across all GMT collections.
# padj is corrected independently per GMT collection — same rationale as run_gsea.
#
# Args:
#   gene_lists  : Named list with slots: up, down, all (any can be NULL).
#                 up   = upregulated DEGs
#                 down = downregulated DEGs
#                 all  = undirected DEG list (used when LFC is unavailable)
#   universe    : Character vector of background genes (all measured genes).
#                 NULL = derive universe from each GMT collection separately.
#                 Statistical implication: a GMT-derived universe is larger and
#                 less specific than the experimental background, which inflates
#                 N, deflates BackgroundRatio, and can produce overly optimistic
#                 enrichment ratios. Use experimental universe whenever possible.
#   gmt_files   : Character vector of full paths to .gmt files.
#   minsize     : Min pathway size to test (default 15).
#   maxsize     : Max pathway size to test (default 500).
#
# Returns: Single formatted dataframe with ORA results across all collections,
#          or NULL if no results were produced.
# ═══════════════════════════════════════════════════════════════════════════════

run_ora <- function(gene_lists, universe = NULL, gmt_files, minsize = 15, maxsize = 500) {

  # Drop NULL or empty gene lists before processing — fora() errors on empty input
  gene_lists <- gene_lists[lengths(gene_lists) > 0]

  if (length(gene_lists) == 0) {
    log_warn(sample = "", step = "run_ora",
             msg    = "All gene lists are NULL or empty. Returning NULL.")
    return(NULL)
  }

  results <- list()

  for (gmt_file in gmt_files) {

    gmt <- fgsea::gmtPathways(gmt_file)

    # Determine the expected prefix based on filename (Human or Mouse)
    file_name <- tolower(basename(gmt_file))

    target_prefix <- dplyr::case_when(
      # Hallmark: h.all (Hs) or mh.all (Mm)
      grepl("^h\\.all|^mh\\.all", file_name)         ~ "HALLMARK",

      # CGP: c2.cgp (Hs) or m2.cgp (Mm)
      grepl("cgp", file_name)                        ~ "CGP",

      # Canonical Pathways
      grepl("biocarta", file_name)                   ~ "BIOCARTA",
      grepl("reactome", file_name)                   ~ "REACTOME",
      grepl("wikipathways", file_name)               ~ "WP",

      # Gene Ontology: c5.go.bp (Hs) or m5.go.bp (Mm)
      grepl("go\\.bp", file_name)                    ~ "GOBP",
      grepl("go\\.cc", file_name)                    ~ "GOCC",
      grepl("go\\.mf", file_name)                    ~ "GOMF",

      # Oncogenic: c6.all (Hs) or m6.all (Mm)
      grepl("c6\\.all|m6\\.all", file_name)          ~ "C6",

      # Catch-all for any other MSigDB collections (C1, C3, C7, etc.)
      grepl("^[cm][0-9]\\.", file_name)              ~ toupper(sub("^([cm][0-9])\\..*", "\\1", file_name)),
      TRUE ~ "OTHER"
    )

    # Only add the prefix if it's NOT already at the start of the names
    if (target_prefix != "") {
      prefix_pattern <- paste0("^", target_prefix, "_")

      # Identify names that DON'T already start with the prefix
      needs_prefix <- !grepl(prefix_pattern, names(gmt))

      if (any(needs_prefix)) {
        names(gmt)[needs_prefix] <- paste0(target_prefix, "_", names(gmt)[needs_prefix])
      }
    }

    # ── Determine universe for this collection ───────────────────────────
    # Two universe strategies with different statistical implications:
    #
    #   Experimental (universe provided):
    #     Filter GMT to genes present in the experiment — genes never measured
    #     cannot be "enriched" and should not count toward pathway size K.
    #     N is fixed across all collections = total measured genes.
    #
    #   GMT-derived (universe = NULL):
    #     Use all genes in the GMT collection as the universe — necessary when
    #     only DEGs are available and no experimental background exists.
    #     N varies per collection (e.g. HALLMARK ~4k genes, REACTOME ~10k).
    #     The ⚠️ warning in the Methods sheet flags this for the user because
    #     GMT-derived universes make cross-collection comparisons less valid.
    if (is.null(universe)) {
      universe_to_use <- unique(unlist(gmt))
      universe_source <- "GMT-derived"

    } else {
      # Filter GMT to only genes present in the experimental background.
      # Genes not measured cannot be enriched — keeping them inflates K
      # and distorts BackgroundRatio.
      gmt             <- lapply(X = gmt, FUN = intersect, y = universe)
      universe_to_use <- universe[!is.na(universe)]
      universe_source <- "Experimental"
    }

    for (type in names(gene_lists)) {

      genes <- gene_lists[[type]]

      res <- fgsea::fora(pathways = gmt,
                         genes    = genes,
                         universe = universe_to_use,
                         minSize  = minsize,
                         maxSize  = maxsize)

      if (is.null(res) || nrow(res) == 0) next

      res$type            <- type
      res$n               <- length(genes)             # input genes tested for this direction
      res$N               <- length(universe_to_use)   # universe size for this collection
      res$universe_source <- universe_source

      results[[paste(gmt_file, type, sep = "_")]] <- res
    }
  }

  combined <- dplyr::bind_rows(results)
  if (is.null(combined) || nrow(combined) == 0) return(NULL)

  # ── Format results ───────────────────────────────────────────────────────
  combined <- combined %>%

    # Same COLLECTION_DESCRIPTION split as run_gsea — see note there for why
    # extra = "merge" is required
    tidyr::separate(col   = pathway,
                    into  = c("Collection", "Description"),
                    sep   = "_",
                    extra = "merge") %>%
    dplyr::mutate(
      Description     = gsub("_", " ", Description),
      Direction       = dplyr::case_when(
        type == "up"   ~ "Enriched in UP genes",
        type == "down" ~ "Enriched in DOWN genes",
        type == "all"  ~ "Enriched in ALL input genes"
      ),
      k               = as.integer(overlap),            # input genes overlapping this pathway
      K               = size,                           # pathway size within universe
      GeneRatio       = k / n,                          # proportion of input genes in pathway
      BackgroundRatio = K / N,                          # baseline frequency of pathway in background
      EnrichmentRatio = GeneRatio / BackgroundRatio,    # fold enrichment over background
      # combined_score combines enrichment magnitude and statistical significance.
      # Equivalent to the scoring used by Enrichr: GeneRatio × -log10(padj).
      # Why padj and not pval? padj is the corrected value — using pval would
      # reward pathways that benefit from lenient multiple testing correction.
      combined_score  = GeneRatio * -log10(padj),
      overlapGenes    = sapply(overlapGenes, paste, collapse = "/")
    ) %>%
    dplyr::select(Collection, Description, Direction, pval, padj, EnrichmentRatio,
                  k, n, K, N, GeneRatio, BackgroundRatio, combined_score, overlapGenes) %>%
    dplyr::filter(k > 1)

  return(combined)
}


# ═══════════════════════════════════════════════════════════════════════════════
# analyze_pathways()
#
# Orchestrator — handles file I/O, input type detection, calls run_gsea/run_ora,
# and writes an Excel output with GSEA, ORA, and Methods sheets.
#
# Args:
#   contrast       : Character. Contrast name — used for output subfolder naming.
#   input_data     : One of:
#                    - File path (csv, xlsx, rds, or plain text gene list)
#                    - data.frame with SYMBOL + padj + log2FoldChange (full DEG results)
#                    - data.frame with SYMBOL + log2FoldChange only (2-col ranked)
#                    - Character vector of gene symbols (gene list)
#   only_deg       : Logical. Relevant only for gene list or 2-col ranked input.
#                    TRUE  = input contains only DEGs (ORA only; GMT universe)
#                    FALSE = input contains all tested genes (GSEA valid)
#                    Ignored when input is a full DEG dataframe (auto-detected).
#   gmt_dir        : Path to directory containing .gmt files.
#   output_dir     : Path to output directory. Results saved in output_dir/contrast/
#   padj_cutoff    : padj threshold for defining DEGs from a full DEG df (default 0.05).
#   minsize        : Minimum pathway size to test (default 15).
#   maxsize        : Maximum pathway size to test (default 500).
#
# Returns: Invisible list with gsea, ora, and methods dataframes.
# ═══════════════════════════════════════════════════════════════════════════════

analyze_pathways <- function(contrast, input_data, gmt_dir, only_deg = FALSE,
                             padj_cutoff = 0.05, minsize = 15, maxsize = 500,
                             output_dir = NULL) {

  input_data <- load_smart(input_data)

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Input Type Detection
  # ═══════════════════════════════════════════════════════════════════════════
  # Three supported input types drive the entire analysis logic:
  #
  #   full_deg_df  : Has SYMBOL + padj + log2FoldChange — the richest input.
  #                  padj lets us separate DEGs from background, enabling both
  #                  GSEA (full ranked list) and ORA (significant genes only).
  #
  #   ranked_2col  : Has SYMBOL + log2FoldChange but NO padj — a pre-ranked list.
  #                  If only_deg = FALSE → all genes present → GSEA valid, ORA not
  #                  (can't define DEGs without padj).
  #                  If only_deg = TRUE  → only DEGs present → ORA valid, GSEA not
  #                  (need full background for GSEA, not just significant genes).
  #
  #   gene_list    : Just gene symbols, no ranking values.
  #                  only_deg = TRUE → ORA only, GMT universe.
  #                  only_deg = FALSE → nothing runs (no ranked values, no DEG flag).

  has_full_deg_cols <- is.data.frame(input_data) &&
    all(c("SYMBOL", "padj", "log2FoldChange") %in% colnames(input_data))

  has_ranked_cols <- is.data.frame(input_data) &&
    all(c("SYMBOL", "log2FoldChange") %in% colnames(input_data)) &&
    !"padj" %in% colnames(input_data)

  is_gene_list <- is.character(input_data) ||
    (is.data.frame(input_data) && ncol(input_data) == 1)

  if (has_full_deg_cols) {
    input_type <- "full_deg_df"

  } else if (has_ranked_cols) {
    input_type <- "ranked_2col"

  } else if (is_gene_list) {
    input_type <- "gene_list"
    # Coerce single-column df to character vector for uniform handling downstream
    if (is.data.frame(input_data)) input_data <- as.character(input_data[[1]])

  } else {
    log_error(sample = "", step = "analyze_pathways",
              msg = glue::glue(
                "Could not detect input type. Provide one of: ",
                "(1) df with SYMBOL + padj + log2FoldChange, ",
                "(2) df with SYMBOL + log2FoldChange, or ",
                "(3) a character vector of gene symbols."))
  }

  log_info(sample = "", step = "analyze_pathways",
           msg    = paste0("Input type detected: ", input_type))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Validate GMT Directory
  # ═══════════════════════════════════════════════════════════════════════════

  if (!dir.exists(gmt_dir)) {
    log_error(sample = "", step = "analyze_pathways",
              msg    = glue::glue("gmt_dir '{gmt_dir}' does not exist."))
  }

  gmt_files <- list.files(gmt_dir, full.names = TRUE, pattern = "\\.gmt$")

  if (length(gmt_files) == 0) {
    log_error(sample = "", step = "analyze_pathways",
              msg    = glue::glue("No .gmt files found in '{gmt_dir}'."))
  }

  log_info(sample = "", step = "analyze_pathways",
           msg = glue::glue("Found {length(gmt_files)} GMT file(s): ",
                            "{paste(basename(gmt_files), collapse = ', ')}"))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Build Analysis Inputs
  # ═══════════════════════════════════════════════════════════════════════════
  # Initialise all flags and inputs to safe defaults. Each case below sets only
  # what it needs — anything not set stays NULL/FALSE, which is the correct
  # "not applicable" state for the Methods sheet and skip-reason logging.
  # The decision logic here is also written to the Methods sheet of the Excel
  # output so the run is fully self-documenting.

  can_gsea         <- FALSE
  can_ora          <- FALSE
  ranked_list      <- NULL
  gene_lists       <- list(up = NULL, down = NULL, all = NULL)
  universe         <- NULL
  universe_label   <- NULL
  gsea_skip_reason <- NULL
  ora_skip_reason  <- NULL

  # ── Case 1: Full DEG dataframe (SYMBOL + padj + log2FoldChange) ─────────
  if (input_type == "full_deg_df") {

    # distinct(SYMBOL) is required because DEG tables can have duplicate gene
    # symbols when multiple Ensembl IDs map to the same symbol. Duplicates
    # would produce duplicate ranked_list names, which fgsea rejects.
    pval_col <- if ("pval" %in% colnames(input_data)) {
      "pval"
    } else if ("padj" %in% colnames(input_data)) {
      "padj"
    } else {
      stop("No p-value column found (expected 'pval' or 'padj')")
    }

    ranked_df <- input_data %>%
      dplyr::filter(!is.na(SYMBOL), !is.na(log2FoldChange), !is.na(.data[[pval_col]])) %>%
      dplyr::mutate(pval_used = ifelse(.data[[pval_col]] == 0, 1e-300, .data[[pval_col]])) %>%
      dplyr::mutate(rank_metric = -log10(pval_used) * sign(log2FoldChange)) %>%
      dplyr::group_by(SYMBOL) %>%
      dplyr::slice_max(order_by = abs(rank_metric), n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(dplyr::desc(rank_metric))

    ranked_list <- stats::setNames(object = as.numeric(ranked_df$rank_metric),
                                   nm     = ranked_df$SYMBOL)

    sig_up   <- ranked_df %>% dplyr::filter(padj <= padj_cutoff, log2FoldChange > 0) %>% dplyr::pull(SYMBOL)
    sig_down <- ranked_df %>% dplyr::filter(padj <= padj_cutoff, log2FoldChange < 0) %>% dplyr::pull(SYMBOL)
    gene_lists <- list(up = sig_up, down = sig_down, all = NULL)

    # GSEA requires a meaningful background of non-significant genes to form a
    # ranked list that spans the full spectrum from activated to repressed.
    # Why not just run GSEA on significant genes only?
    # GSEA tests whether a gene set is non-randomly distributed across the FULL
    # ranked list — if the list contains only DEGs, the "background" is already
    # enriched for signal and the null distribution becomes invalid, producing
    # inflated NES and anti-conservative p-values.
    has_non_sig <- any(ranked_df$padj > padj_cutoff, na.rm = TRUE)

    if (has_non_sig) {
      can_gsea       <- TRUE
      can_ora        <- TRUE
      universe       <- ranked_df$SYMBOL
      universe_label <- paste0("Experimental background — all genes in input df (N=", length(universe), ")")

    } else {
      # All genes are significant — input is effectively a DEG-only list.
      # GSEA is invalid; ORA proceeds with GMT-derived universe.
      can_gsea         <- FALSE
      can_ora          <- TRUE
      gsea_skip_reason <- paste0(
        "All genes in input df have padj <= ", padj_cutoff, ". ",
        "GSEA requires a full ranked list including non-significant genes.")
      universe       <- NULL
      universe_label <- "GMT-derived (per collection) — no experimental background available \u26a0\ufe0f"
    }

  # ── Case 2: 2-column ranked dataframe (SYMBOL + log2FoldChange, no padj) ─
  } else if (input_type == "ranked_2col") {

    ranked_df <- input_data %>%
      dplyr::filter(!is.na(SYMBOL), !is.na(log2FoldChange)) %>%
      dplyr::distinct(SYMBOL, .keep_all = TRUE) %>%
      dplyr::arrange(dplyr::desc(log2FoldChange))

    ranked_list <- stats::setNames(object = as.numeric(ranked_df$log2FoldChange),
                                   nm     = ranked_df$SYMBOL)

    if (!only_deg) {
      # All tested genes present → full ranked background → GSEA valid.
      # ORA skipped: without padj we cannot define which genes are "significant".
      can_gsea        <- TRUE
      can_ora         <- FALSE
      universe        <- ranked_df$SYMBOL
      universe_label  <- paste0("All genes in ranked list (N=", length(universe), ")")
      ora_skip_reason <- paste0("2-column ranked input with only_deg=FALSE. ",
                                "No padj column available to define DEGs for ORA.")

    } else {
      # only_deg = TRUE signals input contains only DEGs → no background → ORA only.
      # LFC sign used to split up/down since padj is unavailable.
      can_gsea         <- FALSE
      can_ora          <- TRUE
      gsea_skip_reason <- paste0("2-column ranked input flagged as DEG-only (only_deg=TRUE). ",
                                 "GSEA requires all tested genes, not just DEGs.")
      universe         <- NULL
      universe_label   <- "GMT-derived (per collection) — no experimental background available \u26a0\ufe0f"
      sig_up   <- ranked_df %>% dplyr::filter(log2FoldChange > 0) %>% dplyr::pull(SYMBOL)
      sig_down <- ranked_df %>% dplyr::filter(log2FoldChange < 0) %>% dplyr::pull(SYMBOL)
      gene_lists <- list(up = sig_up, down = sig_down, all = NULL)
    }

  # ── Case 3: Gene list (character vector or single-column df) ────────────
  } else if (input_type == "gene_list") {

    genes <- unique(na.omit(input_data))

    if (!only_deg) {
      # A gene list with only_deg = FALSE is uninterpretable — we don't know
      # if these are DEGs or background genes, so neither GSEA nor ORA is valid.
      can_gsea         <- FALSE
      can_ora          <- FALSE
      gsea_skip_reason <- "Gene list input with only_deg=FALSE. Cannot run GSEA without ranked values."
      ora_skip_reason  <- "Gene list input with only_deg=FALSE. Set only_deg=TRUE if these are DEGs."
      log_warn(sample = "", step = "analyze_pathways",
               msg    = "No analysis performed. Set only_deg=TRUE if your input is a list of DEGs.")

    } else {
      can_gsea         <- FALSE
      can_ora          <- TRUE
      gsea_skip_reason <- "Gene list input. GSEA requires a ranked list of all tested genes."
      universe         <- NULL
      universe_label   <- "GMT-derived (per collection) — no experimental background available \u26a0\ufe0f"
      gene_lists       <- list(up = NULL, down = NULL, all = genes)
    }
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Run Analyses
  # ═══════════════════════════════════════════════════════════════════════════

  gsea_results <- NULL
  ora_results  <- NULL

  if (can_gsea) {
    log_info(sample = "", step = "analyze_pathways", msg = "Running GSEA...")
    gsea_results <- run_gsea(ranked_list = ranked_list,
                             gmt_files   = gmt_files,
                             minsize     = minsize,
                             maxsize     = maxsize)
  }

  if (can_ora) {
    log_info(sample = "", step = "analyze_pathways", msg = "Running ORA...")
    ora_results <- run_ora(gene_lists = gene_lists,
                           universe   = universe,
                           gmt_files  = gmt_files,
                           minsize    = minsize,
                           maxsize    = maxsize)
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5: Build Methods Sheet
  # ═══════════════════════════════════════════════════════════════════════════
  # Three separate dataframes written to the Methods sheet with gap rows between them:
  #   5a. Logic table      — all possible input scenarios with Analyzed flag
  #   5b. Analysis summary — what actually ran in this specific run
  #   5c. Column defs      — explanation of every output column

  # --- 5a. Logic Table ---

  methods_table <- data.frame(
    Input      = c("Full DEG df — has non-sig genes",
                   "Full DEG df — sig-only (no non-sig rows)",
                   "2-col ranked (SYMBOL + LFC)",
                   "2-col ranked (SYMBOL + LFC)",
                   "Gene list",
                   "Gene list"),
    only_deg   = c("ignored", "ignored", "FALSE", "TRUE", "TRUE", "FALSE"),
    GSEA       = c("YES", "NO",  "YES", "NO",  "NO",  "NO"),
    ORA        = c("YES (up/down)", "YES (up/down)", "NO", "YES (up/down)", "YES (all)", "NO"),
    Universe   = c("All genes in df",
                   "Per-collection GMT genes \u26a0\ufe0f",
                   "All genes in ranked list",
                   "Per-collection GMT genes \u26a0\ufe0f",
                   "Per-collection GMT genes \u26a0\ufe0f",
                   "N/A"),
    Filter_GMT = c("YES", "NO", "YES", "NO", "NO", "NO"),
    Reason     = c("Full experimental background available",
                   "No background — only DEGs in input",
                   "Full measured genes available — no padj for ORA",
                   "DEGs only — no experimental background",
                   "DEGs only — no background, no LFC for direction",
                   "Nothing runs — set only_deg=TRUE if input is DEGs \u26a0\ufe0f"),
    Analyzed   = c(
      ifelse(input_type == "full_deg_df" &&  can_gsea &&  can_ora, "YES <- Your analysis", ""),
      ifelse(input_type == "full_deg_df" && !can_gsea &&  can_ora, "YES <- Your analysis", ""),
      ifelse(input_type == "ranked_2col" && !only_deg,             "YES <- Your analysis", ""),
      ifelse(input_type == "ranked_2col" &&  only_deg,             "YES <- Your analysis", ""),
      ifelse(input_type == "gene_list"   &&  only_deg,             "YES <- Your analysis", ""),
      ifelse(input_type == "gene_list"   && !only_deg,             "YES <- Your analysis", "")),
    stringsAsFactors = FALSE
  )

  # --- 5b. Analysis Summary ---

  summary_rows <- data.frame(
    Parameter = c("=== ANALYSIS SUMMARY ===",
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
                  "GMT files used"),
    Value = c("",
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
              paste(basename(gmt_files), collapse = ", ")),
    stringsAsFactors = FALSE
  )

  # --- 5c. Column Definitions ---

  col_defs <- data.frame(
    Column = c("=== COLUMN DEFINITIONS ===",
               "--- Common Columns (GSEA and ORA) ---",
               "Collection", "Description", "Direction", "pval", "padj", "K", "N", "",
               "--- GSEA-Specific Columns ---",
               "NES", "ES", "log2err", "leading_edge_size", "leadingEdge", "",
               "--- ORA-Specific Columns ---",
               "k", "n", "GeneRatio", "BackgroundRatio", "EnrichmentRatio",
               "combined_score", "overlapGenes"),
    Definition = c(
      "", "",
      "Gene set collection name extracted from pathway name (e.g. HALLMARK, KEGG, REACTOME)",
      "Pathway name with underscores replaced by spaces",
      "For GSEA: Upregulated (NES > 0) or Downregulated (NES < 0). For ORA: Enriched in UP, DOWN, or ALL input genes",
      "Nominal p-value from statistical test",
      "BH-adjusted p-value (FDR) corrected independently per GMT collection",
      "Number of genes in the pathway",
      paste0("Number of background genes. Total genes in ranked list (GSEA) or experimental/GMT-derived background (ORA). ",
             "When GMT-derived, N varies per collection."),
      "",  "",
      "Normalized Enrichment Score. Positive = upregulated, Negative = downregulated",
      "Raw Enrichment Score before normalization",
      "Expected error of the p-value estimate. Lower = more accurate p-value",
      "Number of genes in the leading edge — core genes driving the enrichment signal",
      "/-separated leading edge genes",
      "", "",
      "Number of input genes overlapping with the pathway",
      "Number of input genes tested (UP, DOWN, or ALL depending on direction)",
      "k/n — proportion of input genes that belong to this pathway",
      "K/N — baseline frequency of this pathway in the background",
      "GeneRatio / BackgroundRatio — fold enrichment over background",
      "GeneRatio x -log10(padj) — combined strength and significance score",
      "/-separated input genes overlapping with the pathway"),
    stringsAsFactors = FALSE
  )

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 6: Save Results and Return
  # ═══════════════════════════════════════════════════════════════════════════
  # Methods sheet uses three separate writeData() calls with explicit startRow
  # offsets rather than rbind() — the three sub-tables have different numbers
  # of columns, so binding them would require padding with NAs and would make
  # the sheet unreadable. startRow positions each table independently with
  # two gap rows between sections for visual clarity in Excel.

  safe_contrast <- gsub("[^[:alnum:]_-]", "_", contrast)
  output_dir    <- file.path(output_dir %||% ".", safe_contrast)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  output_file <- file.path(output_dir, "Pathways.xlsx")
  wb          <- openxlsx::createWorkbook()

  if (!is.null(gsea_results) && nrow(gsea_results) > 0) {
    openxlsx::addWorksheet(wb, sheetName = "GSEA")
    openxlsx::writeData(wb, sheet = "GSEA", x = gsea_results, rowNames = FALSE)
  }

  if (!is.null(ora_results) && nrow(ora_results) > 0) {
    openxlsx::addWorksheet(wb, sheetName = "ORA")
    openxlsx::writeData(wb, sheet = "ORA", x = ora_results, rowNames = FALSE)
  }

  openxlsx::addWorksheet(wb, sheetName = "Methods")

  logic_start   <- 1
  summary_start <- logic_start   + nrow(methods_table) + 1 + 2   # +1 header row, +2 gap rows
  coldefs_start <- summary_start + nrow(summary_rows)  + 1 + 2

  openxlsx::writeData(wb, sheet = "Methods", x = methods_table, startRow = logic_start,   rowNames = FALSE)
  openxlsx::writeData(wb, sheet = "Methods", x = summary_rows,  startRow = summary_start, rowNames = FALSE)
  openxlsx::writeData(wb, sheet = "Methods", x = col_defs,      startRow = coldefs_start, rowNames = FALSE)

  openxlsx::saveWorkbook(wb, file = output_file, overwrite = TRUE)

  log_info(sample = "", step = "analyze_pathways",
           msg    = glue::glue("Pathway analysis complete. Results saved to: '{output_file}'"))

  return(invisible(list(gsea    = gsea_results,
                        ora     = ora_results,
                        methods = list(logic_table      = methods_table,
                                       analysis_summary = summary_rows,
                                       column_defs      = col_defs))))
}

# ---- 🚀 4. Smart Execution (Nextflow Only) ----

if (sys.nframe() == 0) {
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
    gmt_dir     = get_arg(3, default = "."),
    only_deg    = as.logical(toupper(get_arg(4, "FALSE"))),
    padj_cutoff = as.numeric(get_arg(5, 0.05)),
    minsize     = as.numeric(get_arg(6, 15)),
    maxsize     = as.numeric(get_arg(7, 500)),
    output_dir  = get_arg(8, default = ".")
  )
}
