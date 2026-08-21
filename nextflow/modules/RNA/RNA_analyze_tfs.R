#!/usr/bin/env Rscript

# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(openxlsx)
  library(decoupleR)
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

analyze_tfs <- function(contrast, species, degs = NULL,
                        log_norm_counts = NULL, metadata = NULL,
                        methods = c("ulm", "mlm", "viper"), minsize = 5,
                        output_dir = NULL) {

  # Why omnipath.offline = TRUE?
  # OmniPath (a dependency of decoupleR) periodically contacts Ensembl and its
  # own servers to check for annotation updates. In HPC / Nextflow environments
  # this causes long hangs or hard failures when outbound internet is blocked.
  # Setting offline = TRUE disables all network checks for this session — the
  # locally cached databases are used instead, which is always what we want in
  # a reproducible pipeline run.
  # We create a temporary log folder inside your existing output_dir
  # to ensure we have write permissions there.
  options(omnipath.offline = TRUE)

  degs_df         <- load_smart(degs)
  log_norm_counts <- load_smart(log_norm_counts)
  metadata        <- load_smart(metadata)

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Resolve Species and Validate Inputs
  # ═══════════════════════════════════════════════════════════════════════════
  # Species string is normalised before any network queries so common
  # abbreviations ("human", "hsap", "Homo_sapiens") all resolve to the formal
  # binomial name that decoupleR's get_collectri() / get_progeny() expect.
  # tolower() + gsub("_", " ") is applied first so matching is case-insensitive
  # and underscore-agnostic.

  clean_species <- gsub(pattern = "_", replacement = " ", x = tolower(species))

  if (grepl(pattern = "human|homo|hsap", x = clean_species)) {
    formal_species <- "Homo sapiens"

  } else if (grepl(pattern = "mouse|mus|mmus", x = clean_species)) {
    formal_species <- "Mus musculus"

  } else if (grepl(pattern = "xeno", x = clean_species)) {
    # Xenograft experiments contain both human (graft) and mouse (host) reads.
    # We annotate against human because the graft signal is what we quantify;
    # mouse reads are filtered out at the alignment stage.
    formal_species <- "Homo sapiens"
    log_warn(sample = species, step = "analyze_tfs",
             msg = "Xeno detected: defaulting to Homo sapiens graft annotation.")

  } else {
    formal_species <- clean_species    # Unknown species — pass through as-is
  }

  log_info(sample = species, step = "analyze_tfs",
           msg = paste("Species resolved to:", formal_species))

  # Hard stop for unsupported species — decoupleR networks only cover human and mouse
  species <- formal_species
  if (!species %in% c("Homo sapiens", "Mus musculus")) {
    log_error(sample = "", step = "analyze_tfs",
              msg    = glue::glue("Unsupported species '{species}'. ",
                                  "Must be 'Homo sapiens' or 'Mus musculus'."))
  }

  # At least one input must be provided:
  #   degs_df         → DEG-level analysis (single ranked list per contrast)
  #   log_norm_counts → Sample-level analysis (per-sample activity scores)
  # Both can be provided — results are stored in separate Excel sheets.
  if (is.null(degs_df) && is.null(log_norm_counts)) {
    log_error(sample = "", step = "analyze_tfs",
              msg    = "Either `degs_df` or `log_norm_counts` must be provided. Both are NULL.")
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Subset log_norm_counts to Contrast-Specific Samples
  # ═══════════════════════════════════════════════════════════════════════════
  # log_norm_counts contains ALL samples from the full experiment. For sample-
  # level TF activity we only want the two groups in this contrast — including
  # other groups would inflate variance estimates and dilute group-specific
  # activity signals.
  # Skipped entirely if log_norm_counts was not provided.

  if (!is.null(log_norm_counts) && !is.null(metadata)) {

    # ── 2a. Identify contrast groups ────────────────────────────────────────
    # all.vars() on the parsed contrast expression extracts group name symbols.
    # e.g. "GroupA - GroupB" → c("GroupA", "GroupB")
    target_groups <- all.vars(base::parse(text = contrast))

    # Find which metadata column contains BOTH group labels.
    # Why all() not any()? any() would match a column that contains only one of
    # the two groups — e.g. a "Batch" column that happens to have a value named
    # "GroupA". all() ensures the column contains every group in the contrast,
    # making the column identification unambiguous.
    group_col <- names(metadata)[
      sapply(metadata, function(x) all(target_groups %in% x))
    ][1]

    if (is.na(group_col)) {
      log_error(sample = "", step = "analyze_tfs",
                msg = glue::glue("Could not find a metadata column containing all ",
                                 "contrast groups: {paste(target_groups, collapse = ', ')}"))
    }

    relevant_samples <- metadata %>%
      dplyr::filter(.data[[group_col]] %in% target_groups) %>%
      dplyr::pull(Sample_ID)

    # ── 2b. Build input matrix ───────────────────────────────────────────────
    # Ensure SYMBOL is an explicit column (may be rownames if loaded from RDS).
    # Why as.data.frame() before checking colnames?
    # Tibbles don't support rownames — coercing first makes the SYMBOL check
    # and rownames_to_column() safe regardless of input class.
    input_df <- log_norm_counts %>% as.data.frame()

    if (!"SYMBOL" %in% colnames(input_df)) {
      input_df <- input_df %>% tibble::rownames_to_column("SYMBOL")
    }

    # Subset to SYMBOL + contrast-relevant samples only.
    # any_of() is used (not all_of()) so missing samples are silently skipped
    # rather than throwing an error — handles cases where a sample was QC-failed
    # and removed from the matrix but is still present in metadata.
    input_df <- input_df %>%
      dplyr::select(SYMBOL, dplyr::any_of(relevant_samples))

    # Move SYMBOL to rownames so the matrix body is purely numeric —
    # decoupleR requires a numeric matrix with gene symbols as rownames.
    # Filter out blank/NA SYMBOLs first to avoid rowname conflicts.
    input_mat <- input_df %>%
      dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
      tibble::column_to_rownames("SYMBOL") %>%
      as.matrix()

    metadata <- metadata %>% dplyr::filter(Sample_ID %in% relevant_samples)
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Load Regulatory Networks
  # ═══════════════════════════════════════════════════════════════════════════

  # ── 3a. Map species to decoupleR organism string ─────────────────────────
  # decoupleR uses "human" / "mouse" rather than binomial names
  organism <- dplyr::case_when(species == "Homo sapiens" ~ "human",
                               species == "Mus musculus" ~ "mouse",
                               TRUE                      ~ "rat")

  # ── 3b. CollecTRI — TF regulatory network ────────────────────────────────
  # Why CollecTRI over DoRothEA?
  # CollecTRI is larger, more recently curated, and has better coverage of
  # experimentally validated TF-target relationships. DoRothEA's confidence
  # levels (A-E) introduce a manual filtering step that can drop many valid TFs.
  #
  # Why split_complexes = FALSE?
  # Some TFs in CollecTRI are annotated as protein complexes (e.g. "RELA_RELB").
  # split_complexes = TRUE would split these into individual subunits, which
  # inflates the number of TF entries and can produce spurious results for
  # subunits that don't independently regulate the target genes. Keeping
  # complexes intact is more conservative and biologically accurate.
  #
  # Note: CollecTRI mor column is binary (-1 or +1) — no continuous edge weights.
  # DoRothEA has continuous weights but they reflect confidence, not biology,
  # so the practical benefit over binary is marginal.
  collectri_db <- decoupleR::get_collectri(organism = organism, split_complexes = FALSE)

  # ── 3c. PROGENy — pathway regulatory network ─────────────────────────────
  # Why top = 500?
  # PROGENy ranks genes by their association with each pathway using a large
  # perturbation dataset. top = 500 takes the 500 most informative genes per
  # pathway. Using all genes adds noise from weak associations; too few (e.g.
  # top = 100) loses sensitivity for pathways with broad gene programs.
  # 500 is the published recommendation for bulk RNA-seq data.

  #progeny_db <- decoupleR::get_progeny(organism = organism, top = 500)
  progeny_db <- NULL

  # Pre-compute network sizes for annotation in output tables
  tf_size      <- collectri_db %>% dplyr::count(source, name = "n_targets")
  #pathway_size <- progeny_db   %>% dplyr::count(source, name = "n_targets")
  pathway_size <- NULL

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Format Input for decoupleR
  # ═══════════════════════════════════════════════════════════════════════════
  # decoupleR expects a numeric matrix: rows = genes (SYMBOL), cols = samples
  # or conditions. The values represent the gene-level statistic used to infer
  # regulatory activity.

  # ── 4a. DEG-level input ──────────────────────────────────────────────────
  # Two sub-cases depending on whether the stat column is available:
  #
  #   stat present   : use directly — it is the Wald test statistic from DESeq2,
  #                    combining effect size and precision in one number. This is
  #                    the ideal input for decoupleR and is always preferred.
  #
  #   stat absent    : stat is dropped by lfcShrink() in some workflows. We
  #                    reconstruct a signed ranking score from pvalue and
  #                    log2FoldChange: t = -log10(pvalue) * sign(LFC).
  #                    Why pvalue and not padj?
  #                    (1) padj has many tied values (BH correction creates
  #                        plateaus) — poor discrimination for ranking.
  #                    (2) padj is NA for genes filtered by DESeq2's independent
  #                        filtering — these genes would be silently dropped.
  #                    (3) pvalue provides a continuous gradient that preserves
  #                        signal from sub-threshold genes that still contribute
  #                        collectively to TF/pathway enrichment.
  #                    Why pmax(pvalue, 1e-300)?
  #                    Some genes have pvalue = 0 due to floating point underflow
  #                    (the true p is below machine precision). log10(0) = -Inf,
  #                    which breaks the matrix. 1e-300 is a conservative floor
  #                    that preserves extreme significance without producing Inf.

  if (!is.null(degs_df)) {

    if ("stat" %in% colnames(degs_df)) {
      input_table <- degs_df %>%
        as.data.frame() %>%
        dplyr::filter(!is.na(stat), !is.na(SYMBOL)) %>%
        dplyr::select(SYMBOL, stat) %>%
        tibble::column_to_rownames("SYMBOL") %>%
        as.matrix()

    } else {
      input_table <- degs_df %>%
        as.data.frame() %>%
        dplyr::mutate(t = -log10(pmax(pvalue, 1e-300)) * sign(log2FoldChange)) %>%
        dplyr::select(SYMBOL, t) %>%
        dplyr::filter(!is.na(t), !is.na(SYMBOL)) %>%
        tibble::column_to_rownames("SYMBOL") %>%
        as.matrix()
    }
  }

  # ── 4b. Sample-level input ───────────────────────────────────────────────
  # Z-score each gene across samples before passing to decoupleR.
  # Why z-score? log_norm_counts are on an absolute expression scale — a
  # highly expressed housekeeping gene dominates activity scores even if it
  # doesn't vary across groups. Z-scoring centres and scales each gene so
  # activity inference reflects RELATIVE expression changes across samples,
  # not absolute expression levels.
  #
  # Why t(scale(t(mat)))?
  # scale() operates column-wise by default (scales each column = each gene
  # across samples). We want to scale row-wise (each gene across samples).
  # t() transposes rows↔cols, scale() operates on the now-column genes,
  # t() transposes back. The double transpose is the standard R idiom for
  # row-wise scaling.
  #
  # Why replace NaN with 0 after scaling?
  # Genes with identical expression across all samples have zero variance —
  # scale() produces NaN (0/0) for these. Replacing with 0 is correct: a
  # gene with no variance contributes no information to activity inference,
  # and a z-score of 0 means "average" which is the appropriate neutral value.

  if (!is.null(log_norm_counts)) {
    input_mat              <- t(scale(t(input_mat)))
    input_mat[is.na(input_mat)] <- 0
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5: Run decoupleR
  # ═══════════════════════════════════════════════════════════════════════════

  # ── 5a. safe_decouple() helper ───────────────────────────────────────────
  # Some decoupleR methods (notably mlm) fail with a collinearity error when
  # the network matrix is rank-deficient for a given input — this happens with
  # small contrasts or very sparse networks. Rather than letting one method
  # crash the entire run, safe_decouple() tests each method individually and
  # skips any that error, running the final decouple() call only on stable methods.
  #
  # Why consensus_score = (length(stable_methods) > 1)?
  # consensus_score computes a z-score across all methods for each TF/pathway,
  # giving a single robust summary score. It requires at least 2 methods to
  # compute a meaningful consensus — with only 1 method there is nothing to
  # aggregate, so we disable it to avoid a spurious single-method "consensus".

  safe_decouple <- function(mat, network, statistics, minsize) {

    stable_methods <- c()

    for (m in statistics) {
      is_stable <- tryCatch({
        decoupleR::decouple(mat = mat, network = network, statistics = m, minsize = minsize)
        TRUE
      }, error = function(e) {
        log_warn(sample = "", step = "analyze_tfs",
                 msg = glue::glue("Method '{m}' failed (likely collinearity). Skipping."))
        FALSE
      })
      if (is_stable) stable_methods <- c(stable_methods, m)
    }

    if (length(stable_methods) == 0) {
      log_warn(sample = "", step = "analyze_tfs",
               msg = "All methods failed for this network/input combination. Returning NULL.")
      return(NULL)
    }

    decoupleR::decouple(mat             = mat,
                        network         = network,
                        statistics      = stable_methods,
                        consensus_score = (length(stable_methods) > 1),
                        minsize         = minsize)
  }

  # ── 5b. run_decoupler() helper ───────────────────────────────────────────
  # Runs safe_decouple() for both TF and pathway networks, then post-processes
  # results: adds BH-adjusted p-values, Direction labels, and network size.
  #
  # Why set.seed(1234)?
  # The viper method uses a permutation-based null distribution — results vary
  # slightly between runs without a fixed seed. Setting the seed here ensures
  # reproducible activity scores and p-values across re-runs.
  #
  # Why group_by(statistic) before p.adjust()?
  # BH correction should be applied within each method independently — pooling
  # p-values across methods (ulm, mlm, viper) would mix p-value distributions
  # from different statistical frameworks, making the correction invalid.
  # Each method produces its own set of TF p-values that should be corrected
  # as a separate family of tests.

  run_decoupler <- function(mat, net_tf, net_pathway) {

    set.seed(1234)

    tf   <- safe_decouple(mat = mat, network = net_tf,       statistics = methods, minsize = minsize)
    #path <- safe_decouple(mat = mat, network = net_pathway,  statistics = methods, minsize = minsize)
    path <- NULL

    tf <- if (!is.null(tf) && nrow(tf) > 0) {
       tf %>%
        dplyr::group_by(statistic) %>%
        dplyr::mutate(
          padj      = stats::p.adjust(p_value, method = "BH"),
          Direction = dplyr::case_when(score > 0 ~ "Upregulated",
                                       score < 0 ~ "Downregulated",
                                       TRUE      ~ "No change")
        ) %>%
        dplyr::ungroup() %>%
        dplyr::left_join(tf_size, by = "source") %>%
        dplyr::filter(Direction != "No change")
    } else NULL

    path <- if (!is.null(path) && nrow(path) > 0) {
      path %>%
        dplyr::group_by(statistic) %>%
        dplyr::mutate(
          padj      = stats::p.adjust(p_value, method = "BH"),
          Direction = dplyr::case_when(score > 0 ~ "Upregulated",
                                       score < 0 ~ "Downregulated",
                                       TRUE      ~ "No change")
        ) %>%
        dplyr::ungroup() %>%
        dplyr::left_join(pathway_size, by = "source") %>%
        dplyr::filter(Direction != "No change")
    } else NULL

    return(list(tf = tf, pathway = path))
  }

  # ── 5c. Run for DEG-level and/or sample-level inputs ────────────────────

  decoupler_results <- list()

  if (!is.null(degs_df)) {
    deg_results <- run_decoupler(input_table, collectri_db, progeny_db)
    decoupler_results[["DEG_TF_Activity"]]      <- deg_results$tf
    decoupler_results[["DEG_Pathway_Activity"]] <- deg_results$pathway
  }

  if (!is.null(log_norm_counts)) {
    mat_results <- run_decoupler(input_mat, collectri_db, progeny_db)
    decoupler_results[["Sample_TF_Activity"]]      <- mat_results$tf
    decoupler_results[["Sample_Pathway_Activity"]] <- mat_results$pathway
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 6: Build Evidence Tables for Plotting
  # ═══════════════════════════════════════════════════════════════════════════
  # The activity results (Sections 5) give a score per TF/pathway. The evidence
  # tables (this section) join the network back to the DEG statistics so that
  # plot_tfs() can show, for each TF, which target genes drove the activity
  # score and whether their expression direction was concordant with the network
  # weight (mor). Only built for DEG-level results — sample-level doesn't have
  # a single DEG stat per gene to join against.

  if (!is.null(degs_df)) {

    # Join network → DEG stats, then filter to only TFs/pathways that survived
    # the activity analysis (i.e. had enough targets, passed minsize filter).
    # TFs absent from deg_results were filtered out by decoupleR — including them
    # in the evidence table would add rows with no corresponding activity score,
    # confusing plot_tfs() which expects a 1:1 match between activity and evidence.

    tf_stats      <- NULL
    pathway_stats <- NULL

    if (!is.null(collectri_db) && !is.null(deg_results$tf)) {
      tf_stats <- collectri_db %>%
        dplyr::add_count(source, name = "n_targets") %>%
        dplyr::left_join(degs_df, by = c("target" = "SYMBOL")) %>%
        dplyr::filter(source %in% unique(deg_results$tf$source))
    }

    if (!is.null(progeny_db) && !is.null(deg_results$pathway)) {
      pathway_stats <- progeny_db %>%
        dplyr::add_count(source, name = "n_targets") %>%
        dplyr::left_join(degs_df, by = c("target" = "SYMBOL")) %>%
        dplyr::filter(source %in% unique(deg_results$pathway$source))
    }

    decoupler_results[["DEG_TF_Evidence"]]      <- tf_stats
    decoupler_results[["DEG_Pathway_Evidence"]] <- pathway_stats
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 7: Save and Return
  # ═══════════════════════════════════════════════════════════════════════════
  # All result sheets are written to a single Excel file — one sheet per
  # analysis type — so plot_tfs() can load them all at once via getSheetNames().

  safe_contrast <- gsub("[^[:alnum:]_-]", "_", contrast)
  output_dir    <- file.path(output_dir %||% ".", safe_contrast)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  output_file <- file.path(output_dir, "TF_and_Pathway_Activity.xlsx")
  openxlsx::write.xlsx(x         = decoupler_results,
                       file      = output_file,
                       overwrite = TRUE)

  log_info(sample = "", step = "analyze_tfs",
           msg    = glue::glue("TF and pathway activity results saved to: '{output_file}'"))

  return(invisible(decoupler_results))
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

  # Split comma-separated method vars into a character vector
  raw_method  <- get_arg(6)
  method_list <- if (!is.null(raw_method)) trimws(strsplit(raw_method, ",")[[1]]) else c("ulm", "mlm", "viper")

  analyze_tfs(
    contrast        = get_arg(1),
    degs            = get_arg(2),
    log_norm_counts = get_arg(3),
    metadata        = get_arg(4),
    species         = get_arg(5),
    methods         = method_list,
    minsize         = as.numeric(get_arg(7, 5)),
    output_dir      = get_arg(8, default = ".")
  )
}
