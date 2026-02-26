# ---- 📦 Load Required Libraries ----
suppressPackageStartupMessages({
  library(crayon)  # Provides green, yellow, red, blue, magenta
  library(glue)    # Provides glue
  library(rlang)   # Provides the %||% operator
})

# To make colors work even in Nextflow logs, force crayon to be active
options(crayon.enabled = TRUE)

# ---- 📝 LOGGING RELATED FUNCTIONS ----

# Suppress warnings and messages
quiet_msg <- function(expr) {
  # create a temp file and open a connection
  tmp <- tempfile()
  con <- file(tmp, open = "wt")
  
  # sink both output and message streams to the same connection
  sink(con, type = "output")
  sink(con, type = "message")
  
  result <- NULL
  tryCatch({
    result <- expr
  }, finally = {
    # restore normal output in the correct order
    sink(type = "message")
    sink(type = "output")
    close(con)
  })
  
  return(result)
}

# Log info messages (green)
log_info <- function(sample, step, msg) {
  sample <- sample %||% ""  # fallback to empty string if NULL
  prefix <- green(formatC("[INFO]", width = 7, flag = " "))
  message(glue::glue("{prefix} [{sample} | {toupper(step)}] {msg}"))
}

# Log warning messages (yellow)
log_warn <- function(sample, step, msg) {
  sample <- sample %||% ""  # fallback to empty string if NULL
  prefix <- yellow(formatC("[WARN]", width = 7, flag = " "))
  message(glue::glue("{prefix} [{sample} | {toupper(step)}] {msg}"))
}

# Log error messages (red)
log_error <- function(sample, step, msg) {
  sample <- sample %||% ""  # fallback to empty string if NULL
  prefix <- red(formatC("[ERROR]", width = 7, flag = " "))
  message(glue::glue("{prefix} [{sample} | {toupper(step)}] {msg}"))
  stop("Workflow stopped.", call. = FALSE)
}

# Optional: header for sample processing
log_sample_header <- function(sample) {
  cat(blue$bold(glue::glue("\n--- Processing Sample: {sample} ---\n\n")))
}

# Optional: section header
log_section <- function(section_name) {
  cat(magenta$bold(glue::glue("\n[{toupper(section_name)}]\n")))
}


# ---- 📝 ANNOTATION FUNCTION ----

add_annotation <- function(df, ann_df) {
  
  # Ensure ID column is character for a clean join
  df <- df %>% dplyr::mutate(ID = as.character(ID))
  
  # If protein-coding gene has gene_id but no gene_name, use gene_id
  # For all other biotypes, if they have no gene_name, keep them as NA
  ann_df  <- ann_df %>%
    # Retain unique gene_id-> gene_name combinations
    select(gene_id, gene_name, gene_biotype) %>% 
    distinct() %>%
    mutate(gene_name = ifelse((is.na(gene_name) | gene_name == "") & (gene_biotype == "protein_coding"), gene_id, gene_name))
  
  log_info(sample = "", step = "add_annotation", msg = "Joining with local GTF annotation...")
  
  # Join and apply your "Protein-Coding fallback" logic
  annotated_df <- df %>%
    dplyr::left_join(ann_df, by = c("ID" = "gene_id")) %>%
    # Use the pre-processed gene_name from your ann_df (which already has the fallbacks)
    # Then filter out the remaining NAs (unnamed non-coding)
    dplyr::filter(!is.na(gene_name) & gene_name != "") %>%
    dplyr::select(gene_name, dplyr::everything()) %>%
    dplyr::rename(SYMBOL = gene_name)
  
  return(annotated_df)
}