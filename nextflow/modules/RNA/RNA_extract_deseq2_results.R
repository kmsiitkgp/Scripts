#!/usr/bin/env Rscript

# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(DESeq2)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(limma)
  library(openxlsx)
  library(glue)
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

extract_deseq2_results <- function(contrast, dds, tx2gene, batch_vars = NULL, 
                                   lfc_cutoff = 0, padj_cutoff = 0.1,
                                   min_samples_per_group = 2, output_dir = NULL) {

  dds     <- load_smart(dds)
  tx2gene <- load_smart(tx2gene)

  # ashr shrinkage (Section 4) uses an EM algorithm with random initialisation.
  # Setting a seed here ensures shrinkage estimates are numerically reproducible
  # across re-runs — without this, tiny differences in LFC estimates can appear
  # between otherwise identical runs, making QC comparisons harder.
  set.seed(1234)

  # Sanitise contrast string into a filesystem-safe directory name early —
  # output_dir is used from Section 5 onward so we build it here once.
  safe_contrast <- gsub("[^[:alnum:]_-]", "_", contrast)
  output_dir    <- file.path(output_dir %||% ".", safe_contrast)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Dynamic Contrast Parsing (Parser Method)
  # ═══════════════════════════════════════════════════════════════════════════
  # The parser method converts a human-readable contrast string (e.g. "A - B",
  # "GroupA - GroupB + GroupC") into a numeric contrast vector over the model
  # matrix columns. This is more general than DESeq2's built-in c() syntax,
  # which can only express simple two-group pairwise comparisons.
  #
  # Why not just use c(factor, numerator, denominator) for everything?
  # Because c() is hardcoded to one factor column + two level names. Any
  # contrast with 3+ terms — additive, interaction, or otherwise — cannot be
  # expressed with c() and requires a numeric coefficient vector instead.
  #
  # Strategy:
  #   1. Build a model matrix from the dds design
  #   2. Map each group name to its mean coefficient vector in that matrix
  #   3. Recursively substitute group names in the contrast expression with
  #      their coefficient vectors, then evaluate the resulting arithmetic
  #   4. The output is a numeric vector of length ncol(model_matrix) that
  #      encodes the desired contrast as a linear combination of coefficients

  # all.vars() robustly extracts variable names from any formula structure,
  # including interactions (e.g. ~Batch + Condition * Treatment)
  design_vars <- all.vars(DESeq2::design(dds))

  # Model matrix used for coefficient-vector contrast construction and rank checks.
  # Why model.matrix() here instead of pulling from dds directly?
  # Because we need the actual numeric matrix to compute per-group colMeans —
  # the dds object does not expose this directly.
  mod_mat <- stats::model.matrix(DESeq2::design(dds),
                                 data = SummarizedExperiment::colData(dds))

  # ── 1a. Build Groups column ──────────────────────────────────────────────
  # We create a "Groups" column to map samples to their biological groups so
  # colMeans() can average the correct model matrix rows per group when
  # building group_coef_list in 1b..
  #
  # 1. Why unite all biological variables?
  #    To ensure distinct groups (e.g., 'HPrEC1.Control' vs 'HPrEC2.Control')
  #    aren't accidentally merged into a single 'Control' group.
  #
  # 2. Why exclude batch_vars?
  #    Batch is a technical covariate. Including it would create group names
  #    like 'BatchA.HPrEC1.Control', which won't match the contrast strings
  #    provided by the user (e.g., 'HPrEC1.Control').
  #
  # Batch effects are still handled in the model matrix; we just omit them
  # from the naming convention used for contrast parsing. Section 2 zeroes them
  # from the contrast vector.

  bio_vars <- if (!is.null(batch_vars)) {
    setdiff(design_vars, batch_vars)
  } else {
    design_vars
  }

  if (length(bio_vars) > 0) {
    df_groups <- SummarizedExperiment::colData(dds) %>%
      as.data.frame() %>%
      tidyr::unite(col = "Groups", dplyr::all_of(bio_vars), sep = "_", remove = FALSE)

    SummarizedExperiment::colData(dds)$Groups <- as.factor(df_groups$Groups)

  } else {
    # No design variables — intercept-only model. All samples belong to one group.
    # Must build df_groups as a proper dataframe (not scalar) so that
    # df_groups$Groups is accessible downstream. Chained assignment
    # (df_groups <- colData()$Groups <- ...) would make df_groups a scalar,
    # causing df_groups$Groups to return NULL and breaking group_coef_list.
    df_groups <- SummarizedExperiment::colData(dds) %>%
      as.data.frame() %>%
      dplyr::mutate(Groups = as.factor("All"))

    SummarizedExperiment::colData(dds)$Groups <- df_groups$Groups
  }

  groups <- unique(df_groups$Groups)

  # ── 1b. Map groups to coefficient vectors ────────────────────────────────
  # For each group, average the model matrix rows belonging to that group.
  # Why colMeans() rather than taking a single row?
  # colMeans() is robust to unbalanced designs — a group with many samples
  # gives the same coefficient vector as one with a single sample, because
  # the model matrix encodes group membership, not raw counts.

  group_coef_list <- lapply(groups, function(i) {
    colMeans(as.matrix(mod_mat[df_groups$Groups == i, , drop = FALSE]))
  })
  names(group_coef_list) <- groups

  # ── 1c. Validate groups before parsing ──────────────────────────────────
  # Check for missing groups BEFORE calling replace_symbols() / eval().
  # If a group name in the contrast string is not in group_coef_list,
  # replace_symbols() will leave it as an unevaluated R symbol and eval()
  # will throw a cryptic "object not found" error. Checking here first gives
  # a clear, informative error message instead.

  parsed_expr    <- base::parse(text = contrast)[[1]]
  missing_groups <- setdiff(all.vars(parsed_expr), names(group_coef_list))

  if (length(missing_groups) > 0) {
    log_error(sample = contrast, step = "extract_deseq2_results",
              msg = glue::glue("These groups in contrast do not exist in the design: ",
                               "{paste(missing_groups, collapse = ', ')}"))
  }

  # ── 1c-bis. Validate sample count PER referenced group ───────────────────
  # Distinct from the missing_groups check above: a group can legitimately
  # EXIST in the design (dds has samples belonging to it) while still having
  # too few of them for THIS specific contrast to be meaningful (e.g. n=1 on
  # one side). This is the per-contrast counterpart to SCRNA_pseudobulk.R's
  # MIN_TOTAL_SAMPLES gate: that gate only asks whether the GROUP AS A WHOLE
  # has enough samples for RNA_CREATE_DDS's DESeq2::DESeq() fit to even run;
  # this asks whether THIS contrast's specific referenced groups are each
  # individually viable, once the fit already exists. This script owns
  # contrast-string parsing as the single source of truth, so this is the
  # correct — and only — place this check lives.
  #
  # Skips gracefully (log_warn + early return), NOT log_error()/stop(): one
  # low-n contrast should not take down the whole pipeline, or sibling
  # contrasts sharing the same dds. This task exits 0 with no output files —
  # every declared output in RNA_extract_deseq2_results.nf is optional: true,
  # so downstream .join()s simply have no entry for this (group, contrast)
  # pair, the same "absence = skip" idiom used elsewhere in this pipeline.

  referenced_groups <- intersect(all.vars(parsed_expr), names(group_coef_list))
  samples_per_group <- table(df_groups$Groups)[referenced_groups]
  low_n_groups      <- referenced_groups[samples_per_group < min_samples_per_group]

  if (length(low_n_groups) > 0) {
    log_warn(sample = contrast, step = "extract_deseq2_results",
             msg = glue::glue("SKIP contrast '{contrast}': group(s) ",
                              "{paste(glue::glue('{low_n_groups}(n={samples_per_group[low_n_groups]})'), collapse = ', ')} ",
                              "have fewer than min_samples_per_group ({min_samples_per_group}) sample(s)."))
    return(invisible(NULL))
  }

  # ── 1d. Recursively evaluate contrast expression ─────────────────────────
  # replace_symbols() walks the parsed R expression tree and substitutes each
  # group name symbol with its coefficient vector from group_coef_list.
  # When it encounters an operator like "-" or "+", it returns the node
  # unchanged — these are R language objects (calls), not group names, and
  # R's eval() knows how to apply arithmetic to numeric vectors.
  # The final eval() call then executes the substituted expression, producing
  # a single numeric contrast vector via ordinary vector arithmetic.

  replace_symbols <- function(node) {
    if (is.symbol(node)) {
      nm <- as.character(node)
      if (nm %in% names(group_coef_list)) return(group_coef_list[[nm]])
      return(node)  # Arithmetic operator ("+", "-" etc.) — return as-is for eval()
    } else if (is.call(node)) {
      return(as.call(lapply(node, replace_symbols)))
    } else {
      return(node)
    }
  }

  expr_sub        <- replace_symbols(parsed_expr)
  contrast_parser <- base::eval(expr_sub)

  # ── 1e. Sanity checks on contrast vector ────────────────────────────────

  if (length(contrast_parser) != ncol(mod_mat)) {
    log_error(sample = contrast, step = "extract_deseq2_results",
              msg = glue::glue("Length of contrast_parser ({length(contrast_parser)}) does not ",
                               "match number of model matrix columns ({ncol(mod_mat)})."))
  }

  if (all(contrast_parser == 0)) {
    log_error(sample = contrast, step = "extract_deseq2_results",
              msg = "contrast_parser is all zeros — invalid contrast. Check the contrast string.")
  }

  if (any(is.na(contrast_parser))) {
    log_error(sample = contrast, step = "extract_deseq2_results",
              msg = "Contrast vector contains NA — possible missing or empty group.")
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Zero Out Batch Covariates
  # ═══════════════════════════════════════════════════════════════════════════
  # When batch variables are in the design (e.g. ~Batch + Condition), the
  # parser contrast vector will contain non-zero coefficients for the batch
  # columns. We zero these out so the contrast only captures the biological
  # signal of interest, not the batch effect.
  #
  # Why zero out rather than exclude batch from the design?
  # Keeping batch in the design lets DESeq2 account for it during dispersion
  # estimation and normalisation — we just don't want it contributing to the
  # final contrast direction.

  if (!is.null(batch_vars)) {

    batch_cols <- grep(paste0("^(", paste(batch_vars, collapse = "|"), ")"),
                       names(contrast_parser), value = TRUE)

    if (length(batch_cols) == 0) {
      log_error(sample = contrast, step = "extract_deseq2_results",
                msg = glue::glue("Batch vars {paste(batch_vars, collapse = ', ')} not found in ",
                                 "model matrix columns. Are they in the design formula?"))
    }

    # Warn if a batch variable is fully confounded with Groups — this means
    # every group has only one unique batch level, so zeroing batch coefficients
    # removes information that cannot be separated from the biological signal.
    for (nv in batch_vars) {
      if (nv %in% colnames(SummarizedExperiment::colData(dds))) {
        tab <- table(df_groups$Groups, SummarizedExperiment::colData(dds)[[nv]])
        if (any(rowSums(tab == 0) == ncol(tab) - 1)) {
          log_error(sample = contrast, step = "extract_deseq2_results",
                    msg = glue::glue("Batch variable '{nv}' appears fully confounded with Groups — ",
                                     "batch and biological signal cannot be separated. Results unreliable!"))
        }
      }
    }

    log_info(sample = contrast, step = "extract_deseq2_results",
             msg = glue::glue("Zeroing batch covariate(s) from contrast vector: ",
                              "'{paste(batch_cols, collapse = ', ')}'"))
    contrast_parser[batch_cols] <- 0
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Standard Contrast (for Cross-Validation)
  # ═══════════════════════════════════════════════════════════════════════════
  # The standard contrast uses DESeq2's built-in c(factor, numerator, denominator)
  # syntax as an independent cross-check against the parser contrast vector.
  #
  # Why cross-validate at all?
  # The parser method is more general but also more complex — a bug in
  # replace_symbols() or an unexpected model matrix structure could silently
  # produce a wrong contrast vector. Running the standard method in parallel
  # and comparing results catches any such discrepancies before they propagate
  # to downstream analyses.
  #
  # Why can't we always use c() and skip the parser?
  # DESeq2's c(factor, numerator, denominator) syntax is ONLY valid for simple
  # pairwise comparisons — exactly ONE factor column and exactly TWO group levels.
  # Any contrast with 3+ terms (additive like "A - B + C", interaction, or
  # multi-group) cannot be expressed with c() at all. The parser handles all
  # cases; c() is used here ONLY as a cross-check for the simple case.
  #
  # We detect "simple pairwise" by checking length(terms) == 2:
  # all.vars() on the contrast string returns the group names referenced.
  # Exactly 2 group names → one subtraction → safe to use c() syntax.
  # 3+ group names → complex contrast → fall back to parser vector for both methods.

  # NOTE: design_vars is NOT re-extracted here — it was already computed in
  # Section 1. condition_col uses tail() because DESeq2 convention places the
  # primary condition of interest as the LAST variable in the design formula
  # (e.g. ~Batch + Condition → condition_col = "Condition"). This assumption
  # holds for standard designs; only used in the simple pairwise branch below.
  condition_col <- tail(design_vars, 1)
  terms         <- all.vars(as.formula(paste0("~", contrast)))

  if (length(terms) == 2) {
    # Simple pairwise (e.g. "A - B") — standard c() syntax is valid
    contrast_standard <- c(condition_col, terms[1], terms[2])
  } else {
    # Complex contrast (3+ terms, additive or interaction) — c() cannot express
    # this, so both methods use the same parser vector. Cross-validation still
    # runs to catch any numerical differences from DESeq2's internal handling.
    contrast_standard <- contrast_parser
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: DEG Extraction and LFC Shrinkage
  # ═══════════════════════════════════════════════════════════════════════════
  # Run DESeq2::results() followed by lfcShrink() for both contrast methods.
  #
  # Why ashr shrinkage rather than apeglm?
  # apeglm only supports single model coefficients (coef= argument) — it cannot
  # accept a numeric contrast vector. ashr works with any contrast type and is
  # robust to varied effect size distributions, making it the only viable
  # choice for complex parser contrasts.
  #
  # Why save stat before shrinkage?
  # lfcShrink() drops the stat column from the results object (it modifies LFC
  # so the original test statistic no longer applies). We preserve stat here
  # because it is needed for TF activity analysis downstream (decoupleR uses
  # the t-statistic as the gene-level input).

  # ── 4a. Parser method ────────────────────────────────────────────────────

  res_parser_raw <- DESeq2::results(object               = dds,
                                    contrast             = contrast_parser,
                                    lfcThreshold         = lfc_cutoff,
                                    altHypothesis        = "greaterAbs",
                                    cooksCutoff          = TRUE,
                                    independentFiltering = TRUE,
                                    alpha                = padj_cutoff,
                                    pAdjustMethod        = "BH")

  # Collect stat BEFORE shrinkage — lfcShrink() will drop this column
  stat_parser <- process_degs(res_parser_raw, tx2gene) %>%
    dplyr::select(gene_id, stat)

  res_parser  <- DESeq2::lfcShrink(dds = dds, res = res_parser_raw, type = "ashr", quiet = TRUE)
  degs_parser <- process_degs(res_parser, tx2gene) %>%
    dplyr::left_join(stat_parser, by = "gene_id")

  # ── 4b. Standard method ──────────────────────────────────────────────────

  res_standard_raw <- DESeq2::results(object               = dds,
                                      contrast             = contrast_standard,
                                      lfcThreshold         = lfc_cutoff,
                                      altHypothesis        = "greaterAbs",
                                      cooksCutoff          = TRUE,
                                      independentFiltering = TRUE,
                                      alpha                = padj_cutoff,
                                      pAdjustMethod        = "BH")

  stat_standard <- process_degs(res_standard_raw, tx2gene) %>%
    dplyr::select(gene_id, stat)

  res_standard  <- DESeq2::lfcShrink(dds = dds, res = res_standard_raw, type = "ashr", quiet = TRUE)
  degs_standard <- process_degs(res_standard, tx2gene) %>%
    dplyr::left_join(stat_standard, by = "gene_id")

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5: VST (Contrast-Specific, Design-Aware)
  # ═══════════════════════════════════════════════════════════════════════════
  # Compute VST on only the samples belonging to the two groups in this contrast.
  #
  # Why subset before VST rather than taking columns from a full-dataset VST?
  # VST with blind = FALSE re-estimates dispersions from only the subset samples.
  # A full-dataset VST includes all samples (including unrelated contrasts),
  # inflating dispersion estimates and producing less accurate variance
  # stabilisation for this specific comparison.
  #
  # Why blind = FALSE?
  # blind = TRUE ignores the design formula when estimating dispersions — this
  # is intended for QC/exploration only. For contrast-specific heatmaps used
  # in biological interpretation, blind = FALSE uses the design to give more
  # accurate stabilisation.
  #
  # Why droplevels() before vst()?
  # After subsetting to two groups, the colData factor columns still carry all
  # original levels from the full dds object. Empty levels cause DESeq2 to
  # build a rank-deficient model matrix, breaking dispersion estimation.
  # droplevels() removes unused levels so the design matrix reflects only the
  # samples actually present in the subset.

  group_names <- all.vars(parsed_expr)
  dds_subset  <- dds[, dds$Groups %in% group_names]

  # CRITICAL: drop unused factor levels before vst() — see note above
  SummarizedExperiment::colData(dds_subset) <- droplevels(
    SummarizedExperiment::colData(dds_subset)
  )

  vst_sub_mat <- DESeq2::vst(dds_subset, blind = FALSE) %>%
    SummarizedExperiment::assay()

  meta_sub <- as.data.frame(SummarizedExperiment::colData(dds_subset))

  # ── 5a. Optional limma batch correction on VST matrix ───────────────────
  # removeBatchEffect() removes technical variation (e.g. batch, sex) from the
  # VST matrix while protecting biological signal via a design matrix (~Groups).
  # Applied AFTER VST (not before) because VST expects raw count-scale data —
  # batch correction on log-scale VST values is the correct order of operations.
  # Only applied when batch_vars is provided.

  if (!is.null(batch_vars) && any(batch_vars != "")) {

    design_sub <- stats::model.matrix(~Groups, data = meta_sub)

    b1 <- meta_sub[[batch_vars[1]]]
    b2 <- if (length(batch_vars) > 1) meta_sub[[batch_vars[2]]] else NULL

    log_info(sample = contrast, step = "extract_deseq2_results",
             msg = "Applying limma::removeBatchEffect to VST matrix.")

    vst_sub_mat <- limma::removeBatchEffect(vst_sub_mat,
                                            batch  = b1,
                                            batch2 = b2,
                                            design = design_sub)
  }

  vst_nonblind_counts <- process_counts(vst_sub_mat, tx2gene)

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 6: Cross-Validate Parser vs Standard Results
  # ═══════════════════════════════════════════════════════════════════════════
  # compare_deseq2_results() performs a full column-by-column comparison of
  # the two result objects and returns a list with an any_differences flag.
  # We use this flag to decide what to save — keeping output lean when results
  # are identical (the expected case for simple contrasts).

  comparison <- compare_deseq2_results(
    res_a      = res_parser,
    res_b      = res_standard,
    output_dir = output_dir,
    label_a    = "Parser",
    label_b    = "Standard"
  )

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 7: Save Results
  # ═══════════════════════════════════════════════════════════════════════════
  # Saving strategy mirrors the cross-validation outcome:
  #   Always saved   : res_parser.rds (shrunken), DEGs.xlsx, VST counts
  #   Only if differ : res_standard.rds (shrunken), DEGs_Standard sheet,
  #                    Comparison_Results.xlsx (written by compare_deseq2_results)
  #
  # We save shrunken RDS objects (res_parser / res_standard) rather than raw
  # (res_parser_raw / res_standard_raw) because shrunken LFCs are what all
  # downstream tools (volcano plots, heatmaps, pathway analysis) consume.
  # The raw objects are not needed after this point.

  # Parser RDS and VST counts — always saved
  saveRDS(res_parser, file = file.path(output_dir, "res_parser.rds"))

  openxlsx::write.xlsx(
    x         = list("VST_NonBlind" = vst_nonblind_counts),
    file      = file.path(output_dir, "VST_NonBlind_Counts_Heatmaps.xlsx"),
    overwrite = TRUE
  )

  if (comparison$any_differences) {

    # Results differ — save standard RDS and both DEG sheets for investigation
    saveRDS(res_standard, file = file.path(output_dir, "res_standard.rds"))

    openxlsx::write.xlsx(
      x         = list("DEGs_Parser"   = degs_parser,
                       "DEGs_Standard" = degs_standard),
      file      = file.path(output_dir, "DEGs.xlsx"),
      overwrite = TRUE
    )

  } else {

    # Results identical — parser sheet only, no standard RDS needed
    openxlsx::write.xlsx(
      x         = list("DEGs_Parser" = degs_parser),
      file      = file.path(output_dir, "DEGs.xlsx"),
      overwrite = TRUE
    )
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 8: Log and Return
  # ═══════════════════════════════════════════════════════════════════════════

  log_info(sample = contrast, step = "extract_deseq2_results",
           msg    = glue::glue("DEG results extracted successfully. Saved to: {output_dir}"))

  return(invisible(list(
    degs_parser   = degs_parser,
    degs_standard = degs_standard,
    res_parser    = res_parser,       # shrunken
    res_standard  = res_standard      # shrunken
  )))
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

  # Split comma-separated batch_vars string into a character vector
  raw_batch  <- get_arg(4)
  batch_list <- if (!is.null(raw_batch)) trimws(strsplit(raw_batch, ",")[[1]]) else NULL

  extract_deseq2_results(
    contrast               = get_arg(1),
    dds                    = get_arg(2),
    tx2gene                = get_arg(3),
    batch_vars             = batch_list,
    lfc_cutoff             = as.numeric(get_arg(5, 0)),
    padj_cutoff            = as.numeric(get_arg(6, 0.1)),
    min_samples_per_group  = as.integer(get_arg(7, 2)),
    output_dir             = get_arg(8, default = ".")
  )
}
