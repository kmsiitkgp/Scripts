# ******************************************************************************
# SURVIVAL & SIGNATURE VALIDATION (sigCheck & Biobase)
# ******************************************************************************

# Create an AnnotatedDataFrame object for meta data
p_data <- meta_data
rownames(p_data) <- make.names(p_data$Sample_ID)
p_data <- Biobase::AnnotatedDataFrame(data = p_data,
                                      varMetadata = data.frame("labelDescription" = colnames(p_data))) 

# Create an AnnotatedDataFrame object for features
f_data <- data.frame("SYMBOL" = rownames(normalized_counts))
rownames(f_data) <- make.names(f_data$SYMBOL)
f_data <- Biobase::AnnotatedDataFrame(data = f_data,
                                      varMetadata = data.frame("labelDescription" = colnames(f_data))) 

# Create an ExpressionSet object for read data
e_data <- as.matrix(normalized_counts)
eset <- Biobase::ExpressionSet(assayData = e_data,
                               phenoData = p_data,
                               featureData = f_data,
                               annotation = "custom")

# NOTE: classes parameter in sigCheck() = Status column of varLabels(eset)
# NOTE: survival parameter in sigCheck() = Time column of varLabels(eset)
varLabels(eset)
eset$Status
eset$Time
# NOTE: annotation parameter in sigCheck() = SYMBOL column of fvarLabels(eset)
# NOTE: plot_genes and genes in SYMBOL column must have some overlap
fvarLabels(eset)   

# Perform anlysis for male and female
# sigCheck is not good at classifying samples into optimal groups. So, we 
# manually classify the samples using survminer and import the classification
# into sigCheck() using scoreMethod and threshold parameters

p <- c()

for (i in 1:100){
  
  plot_genes <- base::sample(x=rownames(normalized_counts), size=64, replace = FALSE)
  
  # Calculate z-score using the function described in above paper
  expr_df <- as.data.frame(advanced_Z(plot_genes, normalized_counts))
  
  # Merge expression data with survival data
  expr_df <- expr_df %>%
    data.frame() %>%
    dplyr::rename(combined.exp = identity(1)) %>%
    tibble::rownames_to_column("Sample_ID") %>%
    dplyr::inner_join(meta_data, by=c("Sample_ID"="Sample_ID"))
  
  gene <- "combined.exp"
  if (nrow(expr_df > 0)){
    summary <- wrangle_data(expr_df, stratify_criteria, prefix)
  }
  
  surv_df <- summary[[1]] %>%
    dplyr::filter(Status == 0 | Status == 1) %>%
    dplyr::mutate(Expression = stringi::stri_replace_all_regex(str = Expression,
                                                               pattern = c("HIGH", "LOW"),
                                                               replacement = c("1", "0"),
                                                               vectorise_all = FALSE)) %>%
    dplyr::filter(Sex == "Male")
  
  sigCheck_score <- function(eset){
    e <- surv_df  %>% 
      dplyr::select(Expression) %>%
      unlist(., use.names = FALSE) %>%
      as.numeric()
    
    return(e)
  }
  
  # Format the object to remove sample not present in surv_df
  # Also, remove samples which have status other than 0 or 1.
  eset_subset <- eset[, eset$Sample_ID %in% surv_df$Sample_ID]
  
  # Create a SigCheck object
  check <- sigCheck(expressionSet = eset_subset, 
                    classes = "Status", 
                    survival = "Time",
                    signature = plot_genes,
                    annotation = "SYMBOL",
                    scoreMethod = sigCheck_score(eset),
                    threshold = sum(sigCheck_score(eset))/length(sigCheck_score(eset)))
  
  p <- c(p, check@survivalPval)
  # sigCheckRandom(check = check,
  #                iterations=100)
}

# ******************************************************************************
# R packages for FlowJo tools
# ******************************************************************************

BiocManager::install(c('flowCore', 'ComplexHeatmap', 'PeacoQC', "FlowSOM", "flowAI"))
install.packages(c('ggplot2', 'FNN', 'igraph', 'Matrix', 'cowsay', 'umap', 
                   'uwot', 'utils', 'devtools', 'data.table', 'dplyr',
                   "pheatmap","png"))

file.path(R.home("bin"), "R")
system("type R")
R.home()

# ******************************************************************************
# Some bulk RNA seq functions (may be useless nowadays)
# ******************************************************************************

### combatseq Batch correction:
# NOTE: This batch correction of known factors is done on raw counts
# The batch corrected raw reads are used in DESeq2
batch_correct_combat <- function(meta_data, read_data, formula_string){ 
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                        colData=meta_data, 
                                        design=~1)
  dds <- DESeq2::estimateSizeFactors(dds) 
  dat  <- DESeq2::counts(dds, normalized = TRUE)
  
  # Full model matrix with the variable of interest
  mod  <- stats::model.matrix(as.formula(formula_string), colData(dds))
  
  if (length(unique(as.numeric(meta_data$Batch))) > 1){
    # Instead of using group & full_mod=TRUE, use covar_mod
    read_data_combat <- sva::ComBat_seq(counts=as.matrix(read_data), 
                                        batch=as.numeric(meta_data$Batch), 
                                        #group=as.numeric(as.factor(meta_data$Condition)),
                                        #full_mod = TRUE,
                                        group = NULL,
                                        covar_mod = mod)
  } else{
    read_data_combat <- read_data
  }
  return(read_data_combat) 
}

### svaseq Batch correction:
# NOTE: This batch correction of unknown factors is done on normalized counts
# NOTE: svaseq() can find n number of surrogate variables. If we model for all 
# of them there could be over correction. Hence, we limit batch correction to
# only the top 3 surrogate variables.
# Here, we just create a new object sva_dds with sva design variables
batch_correct_sva <- function(meta_data, read_data, formula_string){
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                        colData=meta_data, 
                                        design=~1)
  dds <- DESeq2::estimateSizeFactors(dds) 
  dat  <- DESeq2::counts(dds, normalized = TRUE)
  
  # Remove lowly expressed genes before finding surrogate variables
  idx  <- rowMeans(dat) > 1
  dat  <- dat[idx, ]
  
  # Full model matrix with the variable of interest #Reduced/null model matrix with only an intercept term
  mod  <- stats::model.matrix(as.formula(formula_string), colData(dds))
  # Reduced/null model matrix with only an intercept term
  mod0 <- stats::model.matrix(~ 1, colData(dds))
  # Estimate all surrogate variables
  svseq <- sva::svaseq(dat=dat, mod=mod, mod0=mod0)
  
  # If there are more than 3 SV, estimate upto 3 SVs
  # svseq$sv values will change based on n.sv
  if (svseq$n.sv > 3){
    svseq <- sva::svaseq(dat=dat, mod=mod, mod0=mod0, n.sv=3)
  }
  
  # Add these SVs as columns to the DESeqDataSet and then add them to the design
  # ddssva$SV1 <- svseq$sv[,1]
  # ddssva$SV2 <- svseq$sv[,2]
  # design(ddssva) <- ~ SV1 + SV2 + id
  
  ddssva <- dds
  for (i in 1:ncol(svseq$sv)){
    var <- paste0("SV",i)
    ddssva[[var]] <- svseq$sv[,i]
  }
  
  design(ddssva) <- as.formula(paste0("~", 
                                      paste0("SV", seq(1:ncol(svseq$sv)), collapse = "+"), 
                                      "+", 
                                      gsub(pattern="~", replacement="",x=formula_string)))
  
  return(ddssva)
}

# Normalized counts are influenced by sizeFactors.
# sizeFactors are affected by number of samples (all samples vs subset of samples)
# sizeFactors are NOT affected by design formula.
# sizeFactors MUST be estimated first before normalization.
# Normalized counts from dds object are NOT batch corrected. We do this below.
# https://www.biostars.org/p/490181/
norm_counts_DESeq2 <- function(meta_data, read_data, output_path){  
  
  # design doesnt affect size factors. Hence, normalized counts are not affected by design
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                        colData=meta_data, 
                                        design=~ 1)
  # Estimate sizefactors
  dds <- DESeq2::estimateSizeFactors(object = dds)
  
  # EXtract normalized counts from dds object
  normalized_counts <- DESeq2::counts(dds, normalized=TRUE)
  
  # Batch correct
  if (sum(colnames(meta_data) %in% c("Batch")) == 1){
    if (length(unique(meta_data$Batch)) > 1){
      normalized_counts_batch <- limma::removeBatchEffect(x=log2(normalized_counts+1), 
                                                          dds$Batch)
    } else {
      normalized_counts_batch <- normalized_counts
    }
  } else {
    normalized_counts_batch <- normalized_counts
  }
  
  normalized_counts <- normalized_counts %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID")
  
  normalized_counts_batch <- normalized_counts_batch %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID")
  
  # Add gene names
  normalized_counts <- add_annotation(normalized_counts)
  normalized_counts_batch <- add_annotation(normalized_counts_batch)
  
  # Save the results
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="Norm_counts")
  openxlsx::writeData(wb, sheet="Norm_counts", x=normalized_counts, rowNames=FALSE)
  openxlsx::addWorksheet(wb, sheetName="Norm_counts_batch_corrected")
  openxlsx::writeData(wb, sheet="Norm_counts_batch_corrected", x=normalized_counts_batch, rowNames=FALSE)
  openxlsx::saveWorkbook(wb, file=paste0(output_path, "Normalized.counts.DESeq2.xlsx"), 
                         overwrite=TRUE)
  
  return(normalized_counts)
}

norm_counts_combat <- function(meta_data, read_data_combat, output_path){
  
  # design doesnt affect size factors. Hence, normalized counts are not affected by design
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data_combat,
                                        colData=meta_data, 
                                        design=~ 1)
  # Estimate sizefactors
  dds <- DESeq2::estimateSizeFactors(object = dds)
  
  # EXtract normalized counts from dds object
  normalized_counts <- DESeq2::counts(dds, normalized=TRUE)
  
  normalized_counts <- normalized_counts %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID")
  
  # Add gene names
  normalized_counts <- add_annotation(normalized_counts, annotations)
  
  # Save the results
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="Norm_counts")
  openxlsx::writeData(wb, sheet="Norm_counts", x=normalized_counts, rowNames=FALSE)
  openxlsx::saveWorkbook(wb, file=paste0(output_path, "Normalized.counts.combat.xlsx"), 
                         overwrite=TRUE)
  
  return(normalized_counts)
}

# svaseq corrected normalized counts:
# NOTE: ddssva object from svaseq_batch has the top 2 surrogate variables that 
# will be used in DESeq2() but the normalized counts from ddssva object are NOT 
# batch corrected. We do this below.  https://www.biostars.org/p/121489/
# NOTE: Lowly expressed genes are removed before finding surrogate variables.
# So, number of genes is lower than number of DESeq2 normalized counts excel.
norm_counts_sva <- function(meta_data, read_data, formula_string, output_path){
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                        colData=meta_data, 
                                        design=~1)
  dds <- DESeq2::estimateSizeFactors(dds) 
  dat  <- DESeq2::counts(dds, normalized = TRUE)
  
  # Remove lowly expressed genes before finding surrogate variables
  idx  <- rowMeans(dat) > 1
  dat  <- dat[idx, ]
  
  # Full model matrix with the variable of interest #Reduced/null model matrix with only an intercept term
  mod  <- stats::model.matrix(as.formula(formula_string), colData(dds))
  # Reduced/null model matrix with only an intercept term
  mod0 <- stats::model.matrix(~ 1, colData(dds))
  # Estimate all surrogate variables
  svseq <- sva::svaseq(dat=dat, mod=mod, mod0=mod0)
  
  # If there are more than 3 SV, estimate upto 3 SVs
  # svseq$sv values will change based on n.sv
  if (svseq$n.sv > 3){
    svseq <- sva::svaseq(dat=dat, mod=mod, mod0=mod0, n.sv=3)
  }
  
  # %*% indicates Matrix multiplication
  X <- base::cbind(mod, svseq$sv)
  Hat <- base::solve(t(X) %*% X) %*% t(X)
  beta <- (Hat %*% t(dat))
  P <- ncol(mod)
  corrected_data <- dat - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),])
  
  normalized_counts <- corrected_data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID")
  
  normalized_counts <- add_annotation(normalized_counts, annotations) 
  
  # Save batch corrected normalized counts for entire dataset
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="Norm_counts_batch_corrected")
  openxlsx::writeData(wb, sheet="Norm_counts_batch_corrected", x=normalized_counts_batch, rowNames=FALSE)
  openxlsx::saveWorkbook(wb, file=paste0(output_path, "Normalized.counts.SVA.xlsx"), 
                         overwrite=TRUE)
}

#******************************************************************************#
#                              CALCULATE padj AND log2FC
#******************************************************************************#

# Function to calculate pval and log2FoldChange
# norm_counts is matrix with log2 transformed values or non-log transformed values
# DO NOT use log10 transformed values
calc_stats <- function(norm_counts, metadata, Target, Reference, log2_transformed_already){
  
  # Perform t.test
  SYMBOL <- c()
  expt <- c()
  control <- c()
  pval <- c()
  for (j in 1:nrow(norm_counts)){
    
    data <- norm_counts[j,] %>% 
      t() %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "Sample") %>%
      dplyr::left_join(metadata, by = c("Sample" = "Sample"))
    colnames(data) <- c("Sample", "values", "Condition")
    
    # Use F.test() to determine if variances are equal or not.
    # NOTE: p < 0.05 => unequal variance
    
    # NOTE: If Iso <- c(NA, NA, 18) and IP <- c(17, NA, 18),
    # var.test and t.test with throw error since Iso has only 1 value.
    
    # NOTE: If Iso <- c(1,1,1) and IP <- c(1,1,1),
    # t.test will throw error "data are essentially constant".
    
    # NOTE: If Iso <- c(0,0,0) and IP <- c(1,1,1),
    # t.test will throw error "data are essentially constant".
    
    if (sum(!is.na(data[data$Condition == Reference, ]$values)) > 1 &
        sum(!is.na(data[data$Condition == Target, ]$values)) > 1 & 
        (length(unique(data[data$Condition == Reference, ]$values)) + 
         length(unique(data[data$Condition == Target, ]$values)) > 2)){
      #   f_test <- var.test(formula = values ~ Condition, 
      #                      data = data,
      #                      alternative = "two.sided")
      #   
      #   if(!is.na(f_test$p.value)){
      # if (f_test$p.value < 0.05){
      # t_test <- t.test(formula = values ~ Condition,
      #                  data = data,
      #                  alternative = "two.sided",
      #                  var.equal = FALSE)
      # }
      
      # # Remove outliers
      # ref_data <- data[data$Condition == Reference, ]$values
      # low_ref <- quantile(ref_data, na.rm=TRUE)[2] - 1.5*IQR(ref_data, na.rm=TRUE)
      # high_ref <- quantile(ref_data, na.rm=TRUE)[3] + 1.5*IQR(ref_data, na.rm=TRUE)
      # 
      # target_data <- data[data$Condition == Target, ]$values
      # low_tar <- quantile(target_data, na.rm=TRUE)[2] - 1.5*IQR(target_data, na.rm=TRUE)
      # high_tar <- quantile(target_data, na.rm=TRUE)[3] + 1.5*IQR(target_data, na.rm=TRUE)
      # 
      # data <- data %>%
      #   dplyr::filter(!(Condition == Reference & (values > high_ref | values < low_ref))) %>%
      #   dplyr::filter(!(Condition == Target & (values > high_tar | values < low_tar)))
      
      
      # Calculate p values, mean expression
      t_test <- stats::t.test(formula = values ~ Condition, 
                              data = data,
                              alternative = "two.sided",
                              var.equal = FALSE)
      
      if (grepl(Reference, names(t_test$estimate[1]))){
        SYMBOL <- c(SYMBOL, rownames(norm_counts)[j])
        pval <- c(pval, t_test$p.value)
        control <- c(control, t_test$estimate[[1]])
        expt <- c(expt, t_test$estimate[[2]])
      } else if (grepl(Reference, names(t_test$estimate[2]))) {
        SYMBOL <- c(SYMBOL, rownames(norm_counts)[j])
        pval <- c(pval, t_test$p.value)
        control <- c(control, t_test$estimate[[2]])
        expt <- c(expt, t_test$estimate[[1]])
      }
    } else{
      SYMBOL <- c(SYMBOL, rownames(norm_counts)[j])
      pval <- c(pval, 1) # Note: DO NOT SET to NA. It will increase padj.
      control <- c(control, mean(data[data$Condition == Reference, ]$values, na.rm=TRUE))
      expt <- c(expt, mean(data[data$Condition == Target, ]$values, na.rm=TRUE))
    }
  }
  
  stats_df <- data.frame(SYMBOL, expt, control, pval)
  stats_df$padj <- stats::p.adjust(p = stats_df$pval, method = "fdr", n = length(stats_df$pval))
  if(log2_transformed_already){
    stats_df$log2FoldChange <- stats_df$expt - stats_df$control  # if data is already log transformed
  }else{
    stats_df$log2FoldChange <- log(stats_df$expt/stats_df$control, base=2)
  }
  
  result <- norm_counts %>% 
    tibble::rownames_to_column(var = "SYMBOL") %>% 
    dplyr::left_join(stats_df, by = c("SYMBOL" = "SYMBOL"))
  
  return(result)
}

flattenCorrMatrix_pmatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  df <- data.frame(row = rownames(cormat)[row(cormat)[ut]],
                   column = rownames(cormat)[col(cormat)[ut]],
                   cor  =(cormat)[ut],
                   p = pmat[ut])
  
  return(df)
}

flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  df <- data.frame(row = rownames(cormat)[row(cormat)[ut]],
                   column = rownames(cormat)[col(cormat)[ut]],
                   cor  =(cormat)[ut])
  return(df)
}

# Compare two dataframe and output similarity
compare_deg_results <- function(df1, df2, file_suffix, output_path){
  
  df1 <- df1 %>% 
    dplyr::select(SYMBOL, ENSEMBL_ID, ENSEMBL_SYMBOL, log2FoldChange, padj) %>%
    dplyr::filter(padj <= 0.05) %>% 
    dplyr::mutate(padj = round(padj, 2), log2FoldChange = round(log2FoldChange, 2)) 
  
  df2 <- df2 %>% 
    dplyr::select(SYMBOL, ENSEMBL_ID, ENSEMBL_SYMBOL, log2FoldChange, padj) %>%
    dplyr::filter(padj <= 0.05) %>% 
    dplyr::mutate(padj = round(padj, 2), log2FoldChange = round(log2FoldChange, 2)) 
  
  merged_df <- dplyr::full_join(df1, df2,by=c("ENSEMBL_ID"="ENSEMBL_ID")) %>%
    dplyr::mutate(SYMBOL = dplyr::case_when(!is.na(SYMBOL.x) ~ SYMBOL.x,
                                            !is.na(SYMBOL.y) ~ SYMBOL.y,
                                            TRUE ~ ENSEMBL_ID)) %>%
    dplyr::select(SYMBOL, ENSEMBL_ID, log2FoldChange.x,  log2FoldChange.y, padj.x, padj.y) %>%
    dplyr::mutate(Group = dplyr::case_when(padj.x <= 0.05 & padj.y <= 0.05 & log2FoldChange.x >= 0.58 & log2FoldChange.y >= 0.58 ~ "Up in both",
                                           padj.x <= 0.05 & padj.y <= 0.05 & log2FoldChange.x <= -0.58 & log2FoldChange.y <= -0.58 ~ "Down in both",
                                           padj.x <= 0.05 & log2FoldChange.x >= 0.58 & (padj.y > 0.05 | is.na(padj.y)) ~ "Up x",
                                           padj.x <= 0.05 & log2FoldChange.x <= -0.58 & (padj.y > 0.05 | is.na(padj.y)) ~ "Down x",
                                           (padj.x > 0.05 | is.na(padj.x)) & padj.y <= 0.05 & log2FoldChange.y >= 0.58 ~ "Up y",
                                           (padj.x > 0.05 | is.na(padj.x)) & padj.y <= 0.05 & log2FoldChange.y <= -0.58 ~ "Down y",
                                           TRUE ~ "Not Significant in both"))
  
  up.x.y <- nrow(merged_df %>% dplyr::filter(log2FoldChange.x >= 0.58, log2FoldChange.y >= 0.58, padj.x <=0.05, padj.y <= 0.05))
  down.x.y <- nrow(merged_df %>% dplyr::filter(log2FoldChange.x <= -0.58, log2FoldChange.y <= -0.58, padj.x <=0.05, padj.y <= 0.05))
  up.x <- nrow(merged_df %>% dplyr::filter(padj.x <= 0.05 & log2FoldChange.x >= 0.58 & (padj.y > 0.05 | is.na(padj.y))))
  down.x <- nrow(merged_df %>% dplyr::filter(padj.x <= 0.05 & log2FoldChange.x <= -0.58 & (padj.y > 0.05 | is.na(padj.y))))
  up.y <- nrow(merged_df %>% dplyr::filter((padj.x > 0.05 | is.na(padj.x)) & padj.y <= 0.05 & log2FoldChange.y >= 0.58))
  down.y <- nrow(merged_df %>% dplyr::filter((padj.x > 0.05 | is.na(padj.x)) & padj.y <= 0.05 & log2FoldChange.y <= -0.58))               
  
  
  ggplot2::ggplot(data = merged_df, 
                  mapping = aes(x=log2FoldChange.x, y = log2FoldChange.y, color = Group)) +
    geom_point(size=0.75) +
    theme_classic() +
    coord_cartesian(clip = "off") +
    geom_hline(yintercept = c(-0.58, 0.58), color="black", linetype="dashed") +
    geom_vline(xintercept = c(-0.58, 0.58), color="black", linetype="dashed") +
    annotate(geom="text", label=up.x.y, x=10, y=10, col="black", size=5) +
    annotate(geom="text", label=down.x.y , x=-10, y=-10, col="black", size=5) +
    annotate(geom="text", label=up.x, x=10, y=0, col="black", size=5) +
    annotate(geom="text", label=down.x, x=-10, y=0, col="black", size=5) +
    annotate(geom="text", label=up.y, x=0, y=10, col="black", size=5) +
    annotate(geom="text", label=down.y, x=0, y=-10, col="black", size=5)
  
  ggplot2::ggsave(filename = paste0(file_suffix, ".jpg"),
                  plot = last_plot(),
                  device = "tiff",
                  path = output_path,
                  width = 7,
                  height = 7,
                  units = c("in"),
                  dpi = 300,
                  limitsize = TRUE,
                  bg = NULL)
}

#******************************************************************************#
#                         MICROARRAY RELATED FUNCTIONS                         #                       
#******************************************************************************#

# Values should be raw i.e. untransformed. DO NOT use log transformed values etc.
# DO NOT replace NA with 0 etc. Leave NA as they are.
# NA values will be imputed/replaced with average of non-NA values.
# If all values are NA, they will be set to 0.
# NO duplicated genes MUST be present.
impute_with_mean <- function(raw_counts){
  
  # Replace NA with average
  for (j in 1:nrow(raw_counts)){
    
    data <- raw_counts[j,] %>% 
      t() %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "Sample") %>%
      dplyr::left_join(metadata, by = c("Sample" = "Sample"))
    colnames(data) <- c("Sample", "values", "Condition")
    
    # If all values for Reference == NA or Target == NA, set them to 0. 
    # Replace NA with average wherever possible.
    data <- data %>% 
      dplyr::mutate(values = as.numeric(values)) %>%
      dplyr::group_by(Condition) %>% 
      dplyr::mutate(average = mean(values, na.rm=TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(average = dplyr::if_else(is.na(average), 0, average),
                    values = dplyr::if_else(is.na(values), average, values))
    
    data <- data %>% 
      dplyr::select(Sample, values) %>% 
      tibble::column_to_rownames("Sample") %>% 
      t()
    
    if(all(colnames(data) == colnames(raw_counts))){
      raw_counts[j,] <- data[,colnames(raw_counts[j,])] %>% unlist(use.names=FALSE)
    }
  }
  
  imputed_counts <- raw_counts %>% 
    dplyr::mutate(across(.cols=everything(), .fns = as.numeric))
  
  return(imputed_counts)
}

# Perform normalization before imputation so counts can be compared across samples
# NOTE: meta_data MUST have column "Sample" that matches with column names of
# raw_counts 
# NOTE: meta_data MUST have column "Condition" that defines the groups
quantile_norm <- function(raw_counts, meta_data, quant_norm){
  
  # Perform quantile normalization
  if (quant_norm == TRUE){
    
    # https://doi.org/10.1038/s41598-020-72664-6
    # The above paper recommends to quantile normalize within each group rather
    # than whole dataset
    quant_norm_counts <- data.frame(matrix(data=NA, nrow=nrow(raw_counts), ncol=0))
    
    # Subset raw_counts for each group
    for (c in unique(meta_data$Condition)){
      
      samples <- meta_data %>% 
        dplyr::filter(Condition %in% c) %>% 
        dplyr::select(Sample) %>% 
        unlist(use.names=FALSE)
      
      counts <- raw_counts[,samples]
      counts <- as.data.frame(preprocessCore::normalize.quantiles(as.matrix(counts)))
      rownames(counts) <- rownames(raw_counts)
      colnames(counts) <- samples
      
      quant_norm_counts <- dplyr::bind_cols(quant_norm_counts, counts)
    }
  } else {
    quant_norm_counts <- raw_counts
  }
  
  return(quant_norm_counts)
}

#******************************************************************************#
#                             MARKER IDENTIFICATION                            #
#******************************************************************************#

find_markers_across_datasets <- function(datasets){
  
  # Rather than trying to identify clusters that are similar across datasets 
  # and then trying to find marker genes, we try to find which genes are 
  # co-expressed frequently across datasets and then assign cell types based on 
  # such co-expressed gene sets. Since our approach is focused on identifying 
  # markers that are frequently co-expressed across all clusters from all 
  # datasets, it is robust to:
  # (i) sequencing depth (low UMIs dataset vs high UMIs dataset)
  # (ii) cell composition (immune enriched vs whole tumor dataset)
  # (iii) experiment type (single cell vs single nuclei dataset)
  # (iv) purity (high ambient RNA vs low ambient RNA dataset)
  # (v) clustering resolution
  
  datasets <- c("scRNASeq_BBN_C57BL6", "scRNASeq_BBN_Rag", "scRNASeq_GSE217093", 
                "scRNASeq_Jinfen", "scRNASeq_Jyoti", "scRNASeq_GSE164557",
                "scRNASeq_Chen", "scRNASeq_GSE222315", "scRNASeq_HRA003620")
  
  # Get human to mouse ortholog mapping
  ortho <- get_orthologs()
  
  # Create empty dataframe to store markers from all datasets
  markers <- data.frame(cluster = c(""))
  
  # Read markers identified at resolution Harmony 0.8 from each dataset
  for (proj in datasets){
    
    seurat_results <- paste0("/hpc/home/kailasamms/scratch/", proj, "/results_seurat/")
    
    # Genes present repeatedly due to DESCRIPTION column. Remove them.
    df <- read.xlsx(paste0(seurat_results, proj, ".Markers.All.cluster.0.8.harmony.Full.xlsx")) %>%
      dplyr::mutate(Proj = gsub(pattern="scRNASeq_", replacement="", x=proj),
                    Proj.cluster = paste0(Proj, ".", cluster),
                    Avg_Expr = avg_log2FC*pct.1) %>%
      dplyr::distinct_at(c("gene", "cluster", "Proj", "p_val_adj", "avg_log2FC", "pct.1", "pct.2"), .keep_all = TRUE) %>%
      dplyr::select(Proj, cluster, Proj.cluster, gene, p_val_adj, avg_log2FC, pct.1, pct.2, ratio, Avg_Expr)
    
    # Merge markers from all datasets
    markers <- dplyr::bind_rows(markers, df) %>%
      dplyr::filter(!is.na(Proj))
    cat(nrow(markers), "\n")
  }
  
  # Add Human orthologs after removing poor markers
  markers <- markers %>%
    dplyr::filter(p_val_adj <= 0.05) %>% #, pct.1 >= 0.4, ratio >= 2) %>%
    dplyr::left_join(ortho, by=c("gene"="Mouse")) %>%
    dplyr::mutate(avg_log2FC = round(avg_log2FC, 2),
                  ratio = round(ratio, 2),
                  Avg_Expr = round( Avg_Expr, 2),
                  p_val_adj = round(p_val_adj, 2),
                  SYMBOL = dplyr::case_when(is.na(Human) ~ gene,
                                            TRUE ~ Human))
  
  # Get top 100 markers based on avg_log2FC, ratio, Avg_Expr, pct.1 for each cluster
  markers_log2FC <- markers %>%
    dplyr::group_by(Proj.cluster) %>%
    dplyr::slice_max(n=100, order_by = avg_log2FC)
  
  markers_ratio <- markers %>%
    dplyr::group_by(Proj.cluster) %>%
    dplyr::slice_max(n=100, order_by = ratio) %>%
    dplyr::filter(ratio > 1)
  
  markers_pct1 <- markers %>%
    dplyr::group_by(Proj.cluster) %>%
    dplyr::slice_max(n=100, order_by = pct.1) %>%
    # if you exclude this filter, you will get specific but sparsely expressed genes
    dplyr::filter(pct.1 >= 0.4)  
  
  markers_expr <- markers %>%
    dplyr::group_by(Proj.cluster) %>%
    dplyr::slice_max(n=100, order_by = Avg_Expr)
  
  # Merge top 100 markers for each cluster and remove duplicates
  markers_top <- dplyr::bind_rows(markers_log2FC, markers_ratio, 
                                  markers_pct1, markers_expr) %>%
    dplyr::distinct_all(.keep_all = TRUE)
  
  # Get all possible combinations of markers
  marker.pair.df <- data.frame()
  count <- 0
  for (i in unique(markers_top$Proj.cluster)){
    
    # Get markers from each proj and cluster
    markers.subset <- markers_top %>%
      dplyr::filter(Proj.cluster == i)
    
    count <- count+1
    cat(count, ":", nrow(markers.subset), "\t")
    
    if (nrow(markers.subset) >= 2){
      df <- utils::combn(x=markers.subset$SYMBOL, m=2) %>%
        #df <- mixtools::perm(n=length(markers.subset$SYMBOL), r=2, v=markers.subset$SYMBOL)
        t() %>%
        data.frame()
      
      marker.pair.df <- dplyr::bind_rows(marker.pair.df, df)
    }
  }  
  
  # Count all possible combinations of markers from all datasets
  marker.pair.df <- marker.pair.df %>%
    dplyr::rename(PairA = identity(1), PairB = identity(2)) %>%
    dplyr::add_count(PairA, PairB) %>%
    dplyr::distinct_at(c("PairA", "PairB"), .keep_all = TRUE) %>%
    dplyr::rename(n_clusters = n) %>%
    dplyr::filter(n_clusters >=3)   # remove combinations not observed in even 3 clusters
  
  # Calculate overlapping number of genes between PairA and PairB
  marker.pair.df$n_common <- 0  # number of genes commonly coexpressed between A & B
  marker.pair.df$n_PairA <- 0   # number of genes coexpressed with A
  marker.pair.df$n_PairB <- 0   # number of genes coexpressed with B
  marker.pair.df$n_ratio <- 0   # n_common/(n_PairA+n_PairB-n_Common)
  for (i in 1:nrow(marker.pair.df)){
    
    # Find all genes coexpressed with PairA
    coexp_A1 <- marker.pair.df %>% 
      dplyr::filter(PairA == marker.pair.df$PairA[i]) %>% 
      dplyr::select(PairB) %>% unlist(use.names=FALSE)
    coexp_A2 <- marker.pair.df %>% 
      dplyr::filter(PairB ==  marker.pair.df$PairA[i]) %>% 
      dplyr::select(PairA) %>% unlist(use.names=FALSE)
    coexp_A <- unique(coexp_A1, coexp_A2)
    
    # Find all genes coexpressed with PairB
    coexp_B1 <- marker.pair.df %>% 
      dplyr::filter(PairA ==  marker.pair.df$PairB[i]) %>% 
      dplyr::select(PairB) %>% unlist(use.names=FALSE)
    coexp_B2 <- marker.pair.df %>% 
      dplyr::filter(PairB == marker.pair.df$PairB[i]) %>% 
      dplyr::select(PairA) %>% unlist(use.names=FALSE)
    coexp_B <- unique(coexp_B1, coexp_B2)
    
    # Calculate stats
    marker.pair.df$n_common[i] <- length(intersect(coexp_A, coexp_B))
    marker.pair.df$n_PairA[i] <- length(coexp_A)
    marker.pair.df$n_PairB[i] <- length(coexp_B)
    marker.pair.df$n_ratio[i] <- marker.pair.df$n_common[i]/(marker.pair.df$n_PairA[i]+marker.pair.df$n_PairB[i]-marker.pair.df$n_common[i])
    
    cat(i, "\t")
  }
  
  # Remove all pairs that have poor overlap (n_ratio < 0.5) 
  top.marker.pair.df <- marker.pair.df %>% 
    dplyr::filter(n_ratio >= 0.5)
  
  final.markers <- list(T.NK.cell        = c("CD3D"),
                        B.Plasma.cell    = c("CD79A"),
                        Erythrocyte      = c("HBB"),
                        Mast.cell        = c("KIT"),   # tissue resident granule producing cell
                        #Granulocyte      = c(),        # blood resident granule producing cell (Basophil, Eosinophil, Neutrophil)                     
                        Monocyte         = c("GOS2"),  # blood resident phagocyte
                        Macrophage       = c("C1QA"),  # tissue resident phagocyte
                        #Dendritic.cell   = c(),
                        Endothelial.cell = c("VWF"),
                        Myocyte          = c("MYH11"),
                        Neurons          = c("KCNA1"))
  
  # Fill the marker list
  for (i in 1:length(final.markers)){
    
    final.markers[[i]] <- c(final.markers[[i]], marker.pair.df %>% 
                              dplyr::filter(PairA == final.markers[[i]]) %>% 
                              dplyr::select(PairB) %>%
                              unlist(use.names=FALSE))
  }
  
  # Convert to dataframe
  max_l <- max(lengths(final.markers)) 
  final.markers.df <- lapply(X=final.markers, FUN=function(x){c(x, base::rep(x="", times=max_l-length(x)))})
  
  # Save the clustered similarity matrix
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "Markers")
  openxlsx::writeData(wb, sheet = "Markers", x = final.markers.df, rowNames = FALSE)
  openxlsx::addWorksheet(wb, sheetName = "Datasets")
  openxlsx::writeData(wb, sheet = "Datasets", x = data.frame(datasets), rowNames = FALSE)
  openxlsx::addWorksheet(wb, sheetName = "All")
  openxlsx::writeData(wb, sheet = "All", x = marker.pair.df, rowNames = FALSE)
  openxlsx::addWorksheet(wb, sheetName = "Top")
  openxlsx::writeData(wb, sheet = "Top", x = top.marker.pair.df, rowNames = FALSE)
  openxlsx::saveWorkbook(wb = wb, file = "Markers.Compiled.xlsx", overwrite = TRUE)
}

plot_markers_across_datasets <- function(datasets, markers){
  
  # Read integrated seurat object for each dataset
  for (d in datasets){
    seurat_results <- paste0("/hpc/home/kailasamms/scratch/", d, "/results_seurat/")
    integrated.seurat <- readRDS(paste0(seurat_results, "integrated.seurat.rds"))
    Idents(integrated.seurat) <- "cluster.0.8.harmony"
    assign(d, integrated.seurat)
  }
  
  #*****************Plot UMAP at Harmony 0.8 for each dataset******************#
  
  purrr::map(.x = datasets,
             .f = function(x){
               obj <- get(x)
               Seurat::DimPlot(object=obj,
                               reduction="umap.harmony",
                               cols=custom_palette,
                               pt.size=0.2,
                               label.size=1,
                               order = TRUE,  # plot doublets on above rest of cells
                               label = TRUE,
                               raster = FALSE,
                               combine = TRUE) +
                 NoLegend() +
                 ggplot2::labs(title = x,  x="UMAP_1", y="UMAP_2") +
                 custom_theme}) %>%
    cowplot::plot_grid(plotlist=.,
                       align="hv",
                       axis="tblr",
                       #nrow=,
                       ncol=3,
                       rel_widths=1,
                       rel_heights=1,
                       greedy=TRUE,
                       byrow=TRUE)
  
  # Save the plot
  ggplot2::ggsave(filename = "Marker_UMAPs.tiff",
                  plot = last_plot(),
                  device = "jpeg",
                  #path = ,
                  scale = 1,
                  width = 4*3,
                  height = 4*ceiling(length(datasets)/3),
                  units = c("in"),
                  dpi = 300,
                  limitsize = TRUE,
                  bg = "white")
  
  #*****************Plot UMAP of identified markers for each dataset******************#
  
  # Plot these putative markers on UMAP
  for (f in markers){
    
    purrr::map(.x = datasets,
               .f = function(x){
                 obj <- get(x)
                 Seurat::FeaturePlot(object = obj,
                                     features = intersect(rownames(obj@assays$RNA$counts), c(f, stringr::str_to_title(f))),
                                     reduction = "umap.harmony",
                                     cols = c("grey", viridis(n = 10, option = "C", direction = -1)),
                                     pt.size = 0.2,
                                     label.size = 1,
                                     min.cutoff='q10',
                                     order = TRUE,  # plot doublets on above rest of cells
                                     label = TRUE,
                                     raster = FALSE,
                                     combine = TRUE) +
                   NoLegend() +
                   # scale_colour_gradientn(colours=rev(brewer.pal(n=11, name="RdBu"))[5:11])
                   ggplot2::labs(title = x, x="UMAP_1", y="UMAP_2") +
                   custom_theme}) %>% cowplot::plot_grid(plotlist=.,
                                                         align="hv",
                                                         axis="tblr",
                                                         #nrow=,
                                                         ncol=3,
                                                         rel_widths=1,
                                                         rel_heights=1,
                                                         greedy=TRUE,
                                                         byrow=TRUE)
    
    ggplot2::ggsave(filename = paste0(f, ".tiff"),
                    plot = last_plot(),
                    device = "jpeg",
                    #path = ,
                    scale = 1,
                    width = 4*3,
                    height = 4*ceiling(length(datasets)/3),
                    units = c("in"),
                    dpi = 300,
                    limitsize = TRUE,
                    bg = "white")
    
  }
}

### Get human-mouse orthologs
# Output is a dataframe with columns DB.Class.Key, Human, Mouse
# NOTE: Mouse genes (H2-Q10,H2-Q8,H2-Q7,H2-Q6,H2-Q4,H2-Q2,H2-Q1,H2-T23,H2-K1,
# H2-D1) are othologs of the same human gene HLA-A. Similarly, human genes 
# ZNG1A, ZNG1B, ZNG1C, ZNG1E, ZNG1F) are orthologs of the same mouse gene (Zng1).
# So, we make them unique as well syntactically valid (hyphens repalced with .)
# to avoid errors in data analysis.
get_orthologs <- function(){
  
  # This website has a list of human mouse orthologs
  df <- read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
  
  # Get human genes
  df_h <- df %>% 
    dplyr::select(DB.Class.Key, Common.Organism.Name, Symbol) %>%
    dplyr::filter(Common.Organism.Name == "human") %>%
    dplyr::select(DB.Class.Key, Symbol) %>%
    dplyr::distinct_at(c("DB.Class.Key","Symbol"),.keep_all = TRUE) %>%
    dplyr::rename(Human = Symbol)
  
  # Get mouse genes 
  df_m <- df %>% 
    dplyr::select(DB.Class.Key, Common.Organism.Name, Symbol) %>%
    dplyr::filter(Common.Organism.Name == "mouse, laboratory") %>%
    dplyr::select(DB.Class.Key, Symbol) %>%
    dplyr::distinct_at(c("DB.Class.Key","Symbol"),.keep_all = TRUE) %>%
    dplyr::rename(Mouse = Symbol)
  
  # Get human-mouse orthologs & remove genes that dont have orthologs
  df_h_m <- dplyr::full_join(df_h, df_m, by=c("DB.Class.Key"="DB.Class.Key")) %>%
    base::replace(is.na(.), "None") %>%
    dplyr::filter(Human != "None", Mouse != "None")
  #df_h_m[is.na(df_h_m)] <- "None"
  
  # Similar orthologs (mouse and human gene names are identical)
  conf_h_m <- df_h_m %>% dplyr::filter(Human == base::toupper(Mouse))
  
  # Dissimilar orthologs (mouse and human gene names are NOT identical)
  fix_h_m <- df_h_m %>% 
    dplyr::filter(!(Human %in% conf_h_m$Human)) %>%
    dplyr::filter(!(Mouse %in% conf_h_m$Mouse))
  
  # Merge
  df <- dplyr::bind_rows(conf_h_m, fix_h_m) %>%
    dplyr::mutate(Human = make.names(Human, unique=TRUE),
                  Mouse = make.names(Mouse, unique=TRUE))
  
  # Save the excel file  
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb=wb, sheetName="Homologs")
  openxlsx::writeData(wb=wb, sheet="Homologs", x=df)
  openxlsx::saveWorkbook(wb=wb, file="Human.Mouse.Homologs.xlsx", overwrite=TRUE)
  
  return (df)
  
  # # Get mouse and human annotations
  # df <- get_annotations()
  # h <- df[[1]]
  # m <- df[[2]]
  # 
  # # Identify identical genes
  # common <- base::intersect(h$SYMBOL, base::toupper(m$SYMBOL))
  # common_h <- h %>% dplyr::filter(h$SYMBOL %in% common)
  # common_m <- m %>% dplyr::filter(base::toupper(m$SYMBOL) %in% common)
  # 
  # # Identify non-common genes
  # uniq_h <- h %>% dplyr::filter(!(h$SYMBOL %in% common))
  # uniq_m <- m %>% dplyr::filter(!(base::toupper(m$SYMBOL) %in% common))
  # 
  # # Identify similar genes (TIME CONSUMING. SO, save excel and use it in future)
  # # Eg: CD3 and CD3D etc
  # similar_h <- c()
  # similar_m <- c()
  # for (p in unique(uniq_h$SYMBOL)){
  #   
  #   if (sum(grepl(pattern=paste0(p, "."), x=toupper(unique(uniq_m$SYMBOL)))) > 0){
  #     similar_h <- c(similar_h, p)
  #     similar_m <- c(similar_m, unique(uniq_m$SYMBOL)[grepl(pattern=paste0(p, "."), x=toupper(unique(uniq_m$SYMBOL)))][1])
  #                  
  #     #cat(p, ":", unique(uniq_m$SYMBOL)[grepl(pattern=paste0(p, "."), x=toupper(unique(uniq_m$SYMBOL)))][1], "\n")
  #   }
  # }
  # 
  # # Generate final dataframe of human-mouse homolog
  # common_df <- data.frame(HUMAN = unique(sort(common_h$SYMBOL)), 
  #                  MOUSE = unique(sort(common_m$SYMBOL)))
  # similar_df <- data.frame(HUMAN = similar_h,
  #                          MOUSE = similar_m)
  # df1 <- dplyr::bind_rows(common_df, similar_df)
}

### combatseq Batch correction:
# NOTE: This batch correction of known factors is done on raw counts
# The batch corrected raw reads are used in DESeq2
batch_correct_combat <- function(meta_data, read_data, formula_string){ 
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                        colData=meta_data, 
                                        design=~1)
  dds <- DESeq2::estimateSizeFactors(dds) 
  dat  <- DESeq2::counts(dds, normalized = TRUE)
  
  # Full model matrix with the variable of interest
  mod  <- stats::model.matrix(as.formula(formula_string), colData(dds))
  
  if (length(unique(as.numeric(meta_data$Batch))) > 1){
    # Instead of using group & full_mod=TRUE, use covar_mod
    read_data_combat <- sva::ComBat_seq(counts=as.matrix(read_data), 
                                        batch=as.numeric(meta_data$Batch), 
                                        #group=as.numeric(as.factor(meta_data$Condition)),
                                        #full_mod = TRUE,
                                        group = NULL,
                                        covar_mod = mod)
  } else{
    read_data_combat <- read_data
  }
  return(read_data_combat) 
}

### svaseq Batch correction:
# NOTE: This batch correction of unknown factors is done on normalized counts
# NOTE: svaseq() can find n number of surrogate variables. If we model for all 
# of them there could be over correction. Hence, we limit batch correction to
# only the top 3 surrogate variables.
# Here, we just create a new object sva_dds with sva design variables
batch_correct_sva <- function(meta_data, read_data, formula_string){
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                        colData=meta_data, 
                                        design=~1)
  dds <- DESeq2::estimateSizeFactors(dds) 
  dat  <- DESeq2::counts(dds, normalized = TRUE)
  
  # Remove lowly expressed genes before finding surrogate variables
  idx  <- rowMeans(dat) > 1
  dat  <- dat[idx, ]
  
  # Full model matrix with the variable of interest #Reduced/null model matrix with only an intercept term
  mod  <- stats::model.matrix(as.formula(formula_string), colData(dds))
  # Reduced/null model matrix with only an intercept term
  mod0 <- stats::model.matrix(~ 1, colData(dds))
  # Estimate all surrogate variables
  svseq <- sva::svaseq(dat=dat, mod=mod, mod0=mod0)
  
  # If there are more than 3 SV, estimate upto 3 SVs
  # svseq$sv values will change based on n.sv
  if (svseq$n.sv > 3){
    svseq <- sva::svaseq(dat=dat, mod=mod, mod0=mod0, n.sv=3)
  }
  
  # Add these SVs as columns to the DESeqDataSet and then add them to the design
  # ddssva$SV1 <- svseq$sv[,1]
  # ddssva$SV2 <- svseq$sv[,2]
  # design(ddssva) <- ~ SV1 + SV2 + id
  
  ddssva <- dds
  for (i in 1:ncol(svseq$sv)){
    var <- paste0("SV",i)
    ddssva[[var]] <- svseq$sv[,i]
  }
  
  design(ddssva) <- as.formula(paste0("~", 
                                      paste0("SV", seq(1:ncol(svseq$sv)), collapse = "+"), 
                                      "+", 
                                      gsub(pattern="~", replacement="",x=formula_string)))
  
  return(ddssva)
}

norm_counts_combat <- function(meta_data, read_data_combat, output_path){
  
  # design doesnt affect size factors. Hence, normalized counts are not affected by design
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data_combat,
                                        colData=meta_data, 
                                        design=~ 1)
  # Estimate sizefactors
  dds <- DESeq2::estimateSizeFactors(object = dds)
  
  # EXtract normalized counts from dds object
  normalized_counts <- DESeq2::counts(dds, normalized=TRUE)
  
  normalized_counts <- normalized_counts %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID")
  
  # Add gene names
  normalized_counts <- add_annotation(normalized_counts, annotations)
  
  # Save the results
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="Norm_counts")
  openxlsx::writeData(wb, sheet="Norm_counts", x=normalized_counts, rowNames=FALSE)
  openxlsx::saveWorkbook(wb, file=paste0(output_path, "Normalized.counts.combat.xlsx"), 
                         overwrite=TRUE)
  
  return(normalized_counts)
}

# svaseq corrected normalized counts:
# NOTE: ddssva object from svaseq_batch has the top 2 surrogate variables that 
# will be used in DESeq2() but the normalized counts from ddssva object are NOT 
# batch corrected. We do this below.  https://www.biostars.org/p/121489/
# NOTE: Lowly expressed genes are removed before finding surrogate variables.
# So, number of genes is lower than number of DESeq2 normalized counts excel.
norm_counts_sva <- function(meta_data, read_data, formula_string, output_path){
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                        colData=meta_data, 
                                        design=~1)
  dds <- DESeq2::estimateSizeFactors(dds) 
  dat  <- DESeq2::counts(dds, normalized = TRUE)
  
  # Remove lowly expressed genes before finding surrogate variables
  idx  <- rowMeans(dat) > 1
  dat  <- dat[idx, ]
  
  # Full model matrix with the variable of interest #Reduced/null model matrix with only an intercept term
  mod  <- stats::model.matrix(as.formula(formula_string), colData(dds))
  # Reduced/null model matrix with only an intercept term
  mod0 <- stats::model.matrix(~ 1, colData(dds))
  # Estimate all surrogate variables
  svseq <- sva::svaseq(dat=dat, mod=mod, mod0=mod0)
  
  # If there are more than 3 SV, estimate upto 3 SVs
  # svseq$sv values will change based on n.sv
  if (svseq$n.sv > 3){
    svseq <- sva::svaseq(dat=dat, mod=mod, mod0=mod0, n.sv=3)
  }
  
  # %*% indicates Matrix multiplication
  X <- base::cbind(mod, svseq$sv)
  Hat <- base::solve(t(X) %*% X) %*% t(X)
  beta <- (Hat %*% t(dat))
  P <- ncol(mod)
  corrected_data <- dat - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),])
  
  normalized_counts <- corrected_data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID")
  
  normalized_counts <- add_annotation(normalized_counts, annotations) 
  
  # Save batch corrected normalized counts for entire dataset
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="Norm_counts_batch_corrected")
  openxlsx::writeData(wb, sheet="Norm_counts_batch_corrected", x=normalized_counts_batch, rowNames=FALSE)
  openxlsx::saveWorkbook(wb, file=paste0(output_path, "Normalized.counts.SVA.xlsx"), 
                         overwrite=TRUE)
}

#******************************************************************************#
#                              CALCULATE padj AND log2FC
#******************************************************************************#

# Function to calculate pval and log2FoldChange
# norm_counts is matrix with log2 transformed values or non-log transformed values
# DO NOT use log10 transformed values
calc_stats <- function(norm_counts, metadata, Target, Reference, log2_transformed_already){
  
  # Perform t.test
  SYMBOL <- c()
  expt <- c()
  control <- c()
  pval <- c()
  for (j in 1:nrow(norm_counts)){
    
    data <- norm_counts[j,] %>% 
      t() %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "Sample") %>%
      dplyr::left_join(metadata, by = c("Sample" = "Sample"))
    colnames(data) <- c("Sample", "values", "Condition")
    
    # Use F.test() to determine if variances are equal or not.
    # NOTE: p < 0.05 => unequal variance
    
    # NOTE: If Iso <- c(NA, NA, 18) and IP <- c(17, NA, 18),
    # var.test and t.test with throw error since Iso has only 1 value.
    
    # NOTE: If Iso <- c(1,1,1) and IP <- c(1,1,1),
    # t.test will throw error "data are essentially constant".
    
    # NOTE: If Iso <- c(0,0,0) and IP <- c(1,1,1),
    # t.test will throw error "data are essentially constant".
    
    if (sum(!is.na(data[data$Condition == Reference, ]$values)) > 1 &
        sum(!is.na(data[data$Condition == Target, ]$values)) > 1 & 
        (length(unique(data[data$Condition == Reference, ]$values)) + 
         length(unique(data[data$Condition == Target, ]$values)) > 2)){
      #   f_test <- var.test(formula = values ~ Condition, 
      #                      data = data,
      #                      alternative = "two.sided")
      #   
      #   if(!is.na(f_test$p.value)){
      # if (f_test$p.value < 0.05){
      # t_test <- t.test(formula = values ~ Condition,
      #                  data = data,
      #                  alternative = "two.sided",
      #                  var.equal = FALSE)
      # }
      
      # # Remove outliers
      # ref_data <- data[data$Condition == Reference, ]$values
      # low_ref <- quantile(ref_data, na.rm=TRUE)[2] - 1.5*IQR(ref_data, na.rm=TRUE)
      # high_ref <- quantile(ref_data, na.rm=TRUE)[3] + 1.5*IQR(ref_data, na.rm=TRUE)
      # 
      # target_data <- data[data$Condition == Target, ]$values
      # low_tar <- quantile(target_data, na.rm=TRUE)[2] - 1.5*IQR(target_data, na.rm=TRUE)
      # high_tar <- quantile(target_data, na.rm=TRUE)[3] + 1.5*IQR(target_data, na.rm=TRUE)
      # 
      # data <- data %>%
      #   dplyr::filter(!(Condition == Reference & (values > high_ref | values < low_ref))) %>%
      #   dplyr::filter(!(Condition == Target & (values > high_tar | values < low_tar)))
      
      
      # Calculate p values, mean expression
      t_test <- stats::t.test(formula = values ~ Condition, 
                              data = data,
                              alternative = "two.sided",
                              var.equal = FALSE)
      
      if (grepl(Reference, names(t_test$estimate[1]))){
        SYMBOL <- c(SYMBOL, rownames(norm_counts)[j])
        pval <- c(pval, t_test$p.value)
        control <- c(control, t_test$estimate[[1]])
        expt <- c(expt, t_test$estimate[[2]])
      } else if (grepl(Reference, names(t_test$estimate[2]))) {
        SYMBOL <- c(SYMBOL, rownames(norm_counts)[j])
        pval <- c(pval, t_test$p.value)
        control <- c(control, t_test$estimate[[2]])
        expt <- c(expt, t_test$estimate[[1]])
      }
    } else{
      SYMBOL <- c(SYMBOL, rownames(norm_counts)[j])
      pval <- c(pval, 1) # Note: DO NOT SET to NA. It will increase padj.
      control <- c(control, mean(data[data$Condition == Reference, ]$values, na.rm=TRUE))
      expt <- c(expt, mean(data[data$Condition == Target, ]$values, na.rm=TRUE))
    }
  }
  
  stats_df <- data.frame(SYMBOL, expt, control, pval)
  stats_df$padj <- stats::p.adjust(p = stats_df$pval, method = "fdr", n = length(stats_df$pval))
  if(log2_transformed_already){
    stats_df$log2FoldChange <- stats_df$expt - stats_df$control  # if data is already log transformed
  }else{
    stats_df$log2FoldChange <- log(stats_df$expt/stats_df$control, base=2)
  }
  
  result <- norm_counts %>% 
    tibble::rownames_to_column(var = "SYMBOL") %>% 
    dplyr::left_join(stats_df, by = c("SYMBOL" = "SYMBOL"))
  
  return(result)
}

flattenCorrMatrix_pmatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  df <- data.frame(row = rownames(cormat)[row(cormat)[ut]],
                   column = rownames(cormat)[col(cormat)[ut]],
                   cor  =(cormat)[ut],
                   p = pmat[ut])
  
  return(df)
}

flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  df <- data.frame(row = rownames(cormat)[row(cormat)[ut]],
                   column = rownames(cormat)[col(cormat)[ut]],
                   cor  =(cormat)[ut])
  return(df)
}

# Compare two dataframe and output similarity
compare_deg_results <- function(df1, df2, file_suffix, output_path){
  
  df1 <- df1 %>% 
    dplyr::select(SYMBOL, ENSEMBL_ID, ENSEMBL_SYMBOL, log2FoldChange, padj) %>%
    dplyr::filter(padj <= 0.05) %>% 
    dplyr::mutate(padj = round(padj, 2), log2FoldChange = round(log2FoldChange, 2)) 
  
  df2 <- df2 %>% 
    dplyr::select(SYMBOL, ENSEMBL_ID, ENSEMBL_SYMBOL, log2FoldChange, padj) %>%
    dplyr::filter(padj <= 0.05) %>% 
    dplyr::mutate(padj = round(padj, 2), log2FoldChange = round(log2FoldChange, 2)) 
  
  merged_df <- dplyr::full_join(df1, df2,by=c("ENSEMBL_ID"="ENSEMBL_ID")) %>%
    dplyr::mutate(SYMBOL = dplyr::case_when(!is.na(SYMBOL.x) ~ SYMBOL.x,
                                            !is.na(SYMBOL.y) ~ SYMBOL.y,
                                            TRUE ~ ENSEMBL_ID)) %>%
    dplyr::select(SYMBOL, ENSEMBL_ID, log2FoldChange.x,  log2FoldChange.y, padj.x, padj.y) %>%
    dplyr::mutate(Group = dplyr::case_when(padj.x <= 0.05 & padj.y <= 0.05 & log2FoldChange.x >= 0.58 & log2FoldChange.y >= 0.58 ~ "Up in both",
                                           padj.x <= 0.05 & padj.y <= 0.05 & log2FoldChange.x <= -0.58 & log2FoldChange.y <= -0.58 ~ "Down in both",
                                           padj.x <= 0.05 & log2FoldChange.x >= 0.58 & (padj.y > 0.05 | is.na(padj.y)) ~ "Up x",
                                           padj.x <= 0.05 & log2FoldChange.x <= -0.58 & (padj.y > 0.05 | is.na(padj.y)) ~ "Down x",
                                           (padj.x > 0.05 | is.na(padj.x)) & padj.y <= 0.05 & log2FoldChange.y >= 0.58 ~ "Up y",
                                           (padj.x > 0.05 | is.na(padj.x)) & padj.y <= 0.05 & log2FoldChange.y <= -0.58 ~ "Down y",
                                           TRUE ~ "Not Significant in both"))
  
  up.x.y <- nrow(merged_df %>% dplyr::filter(log2FoldChange.x >= 0.58, log2FoldChange.y >= 0.58, padj.x <=0.05, padj.y <= 0.05))
  down.x.y <- nrow(merged_df %>% dplyr::filter(log2FoldChange.x <= -0.58, log2FoldChange.y <= -0.58, padj.x <=0.05, padj.y <= 0.05))
  up.x <- nrow(merged_df %>% dplyr::filter(padj.x <= 0.05 & log2FoldChange.x >= 0.58 & (padj.y > 0.05 | is.na(padj.y))))
  down.x <- nrow(merged_df %>% dplyr::filter(padj.x <= 0.05 & log2FoldChange.x <= -0.58 & (padj.y > 0.05 | is.na(padj.y))))
  up.y <- nrow(merged_df %>% dplyr::filter((padj.x > 0.05 | is.na(padj.x)) & padj.y <= 0.05 & log2FoldChange.y >= 0.58))
  down.y <- nrow(merged_df %>% dplyr::filter((padj.x > 0.05 | is.na(padj.x)) & padj.y <= 0.05 & log2FoldChange.y <= -0.58))               
  
  
  ggplot2::ggplot(data = merged_df, 
                  mapping = aes(x=log2FoldChange.x, y = log2FoldChange.y, color = Group)) +
    geom_point(size=0.75) +
    theme_classic() +
    coord_cartesian(clip = "off") +
    geom_hline(yintercept = c(-0.58, 0.58), color="black", linetype="dashed") +
    geom_vline(xintercept = c(-0.58, 0.58), color="black", linetype="dashed") +
    annotate(geom="text", label=up.x.y, x=10, y=10, col="black", size=5) +
    annotate(geom="text", label=down.x.y , x=-10, y=-10, col="black", size=5) +
    annotate(geom="text", label=up.x, x=10, y=0, col="black", size=5) +
    annotate(geom="text", label=down.x, x=-10, y=0, col="black", size=5) +
    annotate(geom="text", label=up.y, x=0, y=10, col="black", size=5) +
    annotate(geom="text", label=down.y, x=0, y=-10, col="black", size=5)
  
  ggplot2::ggsave(filename = paste0(file_suffix, ".jpg"),
                  plot = last_plot(),
                  device = "tiff",
                  path = output_path,
                  width = 7,
                  height = 7,
                  units = c("in"),
                  dpi = 300,
                  limitsize = TRUE,
                  bg = NULL)
}

#******************************************************************************#
#                         MICROARRAY RELATED FUNCTIONS                         #                       
#******************************************************************************#

# Values should be raw i.e. untransformed. DO NOT use log transformed values etc.
# DO NOT replace NA with 0 etc. Leave NA as they are.
# NA values will be imputed/replaced with average of non-NA values.
# If all values are NA, they will be set to 0.
# NO duplicated genes MUST be present.
impute_with_mean <- function(raw_counts){
  
  # Replace NA with average
  for (j in 1:nrow(raw_counts)){
    
    data <- raw_counts[j,] %>% 
      t() %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "Sample") %>%
      dplyr::left_join(metadata, by = c("Sample" = "Sample"))
    colnames(data) <- c("Sample", "values", "Condition")
    
    # If all values for Reference == NA or Target == NA, set them to 0. 
    # Replace NA with average wherever possible.
    data <- data %>% 
      dplyr::mutate(values = as.numeric(values)) %>%
      dplyr::group_by(Condition) %>% 
      dplyr::mutate(average = mean(values, na.rm=TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(average = dplyr::if_else(is.na(average), 0, average),
                    values = dplyr::if_else(is.na(values), average, values))
    
    data <- data %>% 
      dplyr::select(Sample, values) %>% 
      tibble::column_to_rownames("Sample") %>% 
      t()
    
    if(all(colnames(data) == colnames(raw_counts))){
      raw_counts[j,] <- data[,colnames(raw_counts[j,])] %>% unlist(use.names=FALSE)
    }
  }
  
  imputed_counts <- raw_counts %>% 
    dplyr::mutate(across(.cols=everything(), .fns = as.numeric))
  
  return(imputed_counts)
}

# Perform normalization before imputation so counts can be compared across samples
# NOTE: meta_data MUST have column "Sample" that matches with column names of
# raw_counts 
# NOTE: meta_data MUST have column "Condition" that defines the groups
quantile_norm <- function(raw_counts, meta_data, quant_norm){
  
  # Perform quantile normalization
  if (quant_norm == TRUE){
    
    # https://doi.org/10.1038/s41598-020-72664-6
    # The above paper recommends to quantile normalize within each group rather
    # than whole dataset
    quant_norm_counts <- data.frame(matrix(data=NA, nrow=nrow(raw_counts), ncol=0))
    
    # Subset raw_counts for each group
    for (c in unique(meta_data$Condition)){
      
      samples <- meta_data %>% 
        dplyr::filter(Condition %in% c) %>% 
        dplyr::select(Sample) %>% 
        unlist(use.names=FALSE)
      
      counts <- raw_counts[,samples]
      counts <- as.data.frame(preprocessCore::normalize.quantiles(as.matrix(counts)))
      rownames(counts) <- rownames(raw_counts)
      colnames(counts) <- samples
      
      quant_norm_counts <- dplyr::bind_cols(quant_norm_counts, counts)
    }
  } else {
    quant_norm_counts <- raw_counts
  }
  
  return(quant_norm_counts)
}

#******************************************************************************#
#                         UHMAN - MOUSE HOMOLOGS                              #                       
#******************************************************************************#
# Mouse annotations from ensembl
mouse_full <- read.xlsx("C:/Users/KailasammS/Desktop/Mouse.xlsx", sheet="MOUSE")
# Human annotations from ensembl
human_full <- read.xlsx("C:/Users/KailasammS/Desktop/Mouse.xlsx", sheet="HUMAN")

#Get human_mouse homologs from https://www.informatics.jax.org/downloads/reports/index.html
# Human and Mouse Homology Classes with Sequence information (tab-delimited)
# Separate them into human and mouse on different sheets
mouse_homology <- read.xlsx("C:/Users/KailasammS/Desktop/homology.xlsx", sheet="MOUSE")
human_homology <- read.xlsx("C:/Users/KailasammS/Desktop/homology.xlsx", sheet="HUMAN")

# Filter genes without gene symbols, remove duplicates
mouse_full <- mouse_full %>% 
  dplyr::filter(nchar(MOUSE_SYMBOL) > 0) %>%
  dplyr::distinct(pick(MOUSE_SYMBOL), .keep_all = TRUE)

# Filter genes without gene symbols, remove duplicates
human_full <- human_full %>% 
  dplyr::filter(nchar(HUMAN_SYMBOL) > 0) %>%
  dplyr::distinct(pick(HUMAN_SYMBOL), .keep_all = TRUE)

# Filter duplicated homology pairs
mouse_homology <- mouse_homology %>%
  dplyr::distinct(pick(M_Symbol, DB.Class.Key), .keep_all = TRUE)

# Filter duplicated homology pairs
human_homology <- human_homology %>%
  dplyr::distinct(pick(H_Symbol, DB.Class.Key), .keep_all = TRUE)

# # Find duplicates
# nrow(mouse_full %>% dplyr::count(MOUSE_SYMBOL) %>% dplyr::filter(n>1))
# nrow(mouse_homology %>% dplyr::count(M_Symbol) %>% dplyr::filter(n>1))
# length(unique(intersect(mouse_full$MOUSE_SYMBOL, mouse_homology$M_Symbol)))

mouse_final <- dplyr::inner_join(mouse_full, mouse_homology, by=c("MOUSE_SYMBOL"="M_Symbol")) %>%
  dplyr::select(DB.Class.Key, MOUSE_ID, MOUSE_SYMBOL, M_Genetic.Location, everything()) %>%
  dplyr::rename(MOUSE_GENETIC_LCATION =M_Genetic.Location)

# # Find duplicates
# nrow(human_full %>% dplyr::count(HUMAN_SYMBOL) %>% dplyr::filter(n>1))
# nrow(human_homology %>% dplyr::count(H_Symbol) %>% dplyr::filter(n>1))
# length(unique(intersect(human_full$HUMAN_SYMBOL, human_homology$H_Symbol)))

human_final <- dplyr::inner_join(human_full, human_homology, by=c("HUMAN_SYMBOL"="H_Symbol")) %>%
  dplyr::select(DB.Class.Key, HUMAN_ID, HUMAN_SYMBOL, H_Genetic.Location, everything()) %>%
  dplyr::rename(HUMAN_GENETIC_LCATION =H_Genetic.Location)


homologs1 <- dplyr::inner_join(mouse_homology, human_homology, by=c("DB.Class.Key"="DB.Class.Key")) %>%
  dplyr::select(DB.Class.Key, MOUSE_ID, MOUSE_SYMBOL, HUMAN_ID, HUMAN_SYMBOL, everything(), -c(MOUSE_DESCRIPTION, HUMAN_DESCRIPTION))


#******************************************************************************#
#                                   BAR PLOT                                   #
#******************************************************************************#

# Function to plot bar chart
plot_bar <- function(data, plot_param){
  
  ggplot2::ggplot(data = data,
                  aes(x = !!rlang::sym(plot_param$data_x), 
                      y = reorder(stringr::str_to_upper(string = !!rlang::sym(plot_param$data_y)), desc(!!rlang::sym(plot_param$data_fill))),
                      fill = !!rlang::sym(plot_param$data_fill))) +
    
    # Plot bar plot
    ggplot2::geom_col(width = 0.75) +
    
    # Define the theme of plot
    ggplot2::theme_classic() +
    
    # Define the axis, plot headings
    ggplot2::labs(x = plot_param$title_x,
                  y = plot_param$title_y,
                  fill = plot_param$title_legend_fill,
                  color = plot_param$title_legend_color,
                  size = plot_param$title_legend_size,
                  title = plot_param$title_plot) +
    
    # Define x-axis start and end
    ggplot2::coord_cartesian(xlim = c(0, ceiling(max(abs(data[,plot_param$data_x]))))) +
    
    # Define the color of the bars
    viridis::scale_fill_viridis(discrete = FALSE, option = "D") +
    
    # Adjust font size, style
    my_theme +
    ggplot2::theme(axis.text.y = element_text(family = "sans", face = "plain", colour = "black", size = 10, hjust = 1))
  
  # Save the plot
  ggplot2::ggsave(filename = paste0("Bar_plot_", file_suffix, ".pdf"),
                  plot = last_plot(),
                  device = "pdf",
                  path = results_path,
                  scale = 1,
                  width = 2+nrow(data),
                  height = (2+nrow(data))*aspect_ratio,
                  units = c("in"),
                  dpi = 600,
                  limitsize = TRUE,
                  bg = NULL)
  
  # Create a workbook to store the results
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb = wb, sheetName = "Bar")
  openxlsx::writeData(wb = wb, sheet = "Bar", x = data)
  openxlsx::saveWorkbook(wb, file = paste0(results_path, "Bar_plot_", file_suffix, "_Results.xlsx"),
                         overwrite = TRUE)
  
}

#******************************************************************************#
#                                   DOT PLOT                                   #
#******************************************************************************#

# Function to plot dot plot
plot_dot <- function(data, plot_param){
  
  ggplot2::ggplot(data = data,
                  aes(x = !!rlang::sym(plot_param$data_x), 
                      y = reorder(stringr::str_to_upper(string = !!rlang::sym(plot_param$data_y)), desc(!!rlang::sym(plot_param$data_fill))),
                      size = !!rlang::sym(plot_param$data_size),
                      color = !!rlang::sym(plot_param$data_color))) +
    
    # Plot dot plot
    ggplot2::geom_point() +
    
    # Define the theme of plot
    ggplot2::theme_classic() +
    
    # Define the axis, plot headings
    ggplot2::labs(x = plot_param$title_x,
                  y = plot_param$title_y,
                  color = plot_param$title_legend_color,
                  size = plot_param$title_legend_size,
                  title = plot_param$title_plot) +
    
    # Define x-axis start and end
    ggplot2::coord_cartesian(xlim = c(0, ceiling(max(abs(data[,plot_param$data_x]))))) +
    
    # Define the color of the dots
    viridis::scale_color_viridis(discrete = FALSE, option = "D") +
    
    # Adjust font size, style
    my_theme +
    ggplot2::theme(axis.text.y = element_text(family = "sans", face = "plain", colour = "black", size = 10, hjust = 1))
  
  #ggplot2::scale_size(breaks = sapply(as.vector(quantile(c(min(data$Count), max(data$Count)))), floor))
  
  # Save the plot
  ggplot2::ggsave(filename = paste0("Dot_plot_", file_suffix, ".pdf"),
                  plot = last_plot(),
                  device = "pdf",
                  path = results_path,
                  scale = 1,
                  width = 2+nrow(data),
                  height = (2+nrow(data))*aspect_ratio,
                  units = c("in"),
                  dpi = 600,
                  limitsize = TRUE,
                  bg = NULL)
  
  # Create a workbook to store the results
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb = wb, sheetName = "Dot")
  openxlsx::writeData(wb = wb, sheet = "Dot", x = data)
  openxlsx::saveWorkbook(wb, file = paste0(results_path, "Dot_plot_", file_suffix, "_Results.xlsx"),
                         overwrite = TRUE)
}

#******************************************************************************#
#                               STACKED BAR PLOT                               #
#******************************************************************************#

# Function to plot stacked bar chart
plot_stackedbar <- function(data, plot_param, label_percent){
  
  if (already_percent == FALSE){
    # Calculate percent of each sub type for each sample
    data <- data %>%
      data.frame() %>%
      base::replace(is.na(.), 0) %>%
      dplyr::mutate(across(.cols = c(everything(), -plot_param$data_x), .fns = as.numeric)) %>%
      dplyr::mutate(across(.cols = c(everything(), -plot_param$data_x), .fns = function(x) x*100/rowSums(x = select_if(., is.numeric)))) %>%
      tidyr::pivot_longer(cols = !rlang::sym(plot_param$data_x), names_to = plot_param$data_fill, values_to = "Percent") %>%
      dplyr::mutate(n_gene = gsub(pattern="X", replacement="", x=n_gene))
  }
  
  # Plot stacked bar chart
  p <- ggplot2::ggplot(data = data, 
                       aes(x = !!rlang::sym(plot_param$data_x), 
                           y = Percent, 
                           fill = !!rlang::sym(plot_param$data_fill))) +
    
    # Plot stacked bar plot
    ggplot2::geom_col(position = "stack", width = 0.95) +
    
    # Define the theme of plot
    ggplot2::theme_classic() + 
    
    # Define the axis, plot headings
    ggplot2::labs(x = plot_param$title_x,
                  y = plot_param$title_y,
                  fill = plot_param$title_legend_fill,
                  title = plot_param$title_plot) +
    
    # Define the color of the bars
    ggplot2::scale_fill_manual(values = my_palette) +
    
    # Adjust font size, style
    my_theme
  
  if(label_percent == "TRUE"){
    p <- p +
      ggplot2::geom_text(aes(x = !!rlang::sym(plot_param$data_x), label = round(Percent,1)), 
                         position = position_stack(vjust = 0.5), 
                         fontface = "plain", colour = "white", size = 3, 
                         check_overlap = TRUE)
  }
  
  # Save the plot
  ggplot2::ggsave(filename = paste0("Stacked_Bar_", file_suffix, ".pdf"),
                  #plot = last_plot(),
                  plot = p,
                  device = "pdf",
                  path = results_path,
                  scale = 1,
                  width = dplyr::if_else((2+nrow(data)) < 10, (2+nrow(data)), 10),
                  height = dplyr::if_else((2+nrow(data))*aspect_ratio < 10, (2+nrow(data))*aspect_ratio, 10),
                  units = c("in"),
                  dpi = 600,
                  limitsize = TRUE,
                  bg = NULL)
  
  # Create a workbook to store the results
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb = wb, sheetName = "Stacked_Bar")
  openxlsx::writeData(wb = wb, sheet = "Stacked_Bar", x = data)
  openxlsx::saveWorkbook(wb, file = paste0(results_path, "Stacked_Bar_", file_suffix, "_Results.xlsx"),
                         overwrite = TRUE)
  
}


#******************************************************************************#
#                                 VIOLIN PLOT                                  #
#******************************************************************************#

# Function to plot violin plot
plot_violin <- function(data, plot_param){
  
  p <- ggplot(data = data, 
              aes(x = !!rlang::sym(plot_param$data_x), 
                  y = !!rlang::sym(plot_param$data_y), 
                  fill = !!rlang::sym(plot_param$data_fill))) +
    
    # Plot violin plot with a small box plot within it
    geom_violin(trim=FALSE) + 
    geom_boxplot(width = 0.1) +
    
    # Define the theme of plot
    theme_classic() + 
    
    # Define the axis, plot headings
    labs(x = plot_param$title_x, 
         y = plot_param$title_y, 
         title = plot_param$title_plot)  +
    
    # Define y-axis start and end
    #coord_cartesian(ylim = c(0, ceiling(max(data[,plot_param$data_y])/10)*10)) +
    
    # Define the color of the violins
    ggplot2::scale_fill_manual(values = my_palette) +
    
    # Adjust font size, style
    my_theme 
  
  if (log_scale_y == "TRUE"){
    p <- p +
      scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000)) 	#display y axis in log scale 
  }        
  
  if (show_intercept_y == "TRUE"){
    p <- p +
      geom_hline(yintercept = yintercept_cutoff, linetype = 2)
  }
  
  
  # Save the plot
  ggplot2::ggsave(filename = paste0("Violin_plot_", file_suffix, ".pdf"),
                  #plot = last_plot(),
                  plot = p,
                  device = "pdf",
                  path = results_path,
                  scale = 1,
                  width = 2+length(unique(data[,plot_param$data_x])),
                  height = (2+length(unique(data[,plot_param$data_x])))*aspect_ratio,
                  units = c("in"),
                  dpi = 600,
                  limitsize = TRUE,
                  bg = NULL)
  
}

#******************************************************************************#
#                                   BOX PLOT                                   #
#******************************************************************************#

# Function to plot box plot
plot_box <- function(data, plot_param){
  
  p <- ggplot(data = data, 
              aes(x = !!rlang::sym(plot_param$data_x), 
                  y = !!rlang::sym(plot_param$data_y), 
                  fill = !!rlang::sym(plot_param$data_fill))) +
    
    # Plot box plot
    geom_boxplot() +
    
    # Define the theme of plot
    theme_classic() + 
    
    # Define the axis, plot headings
    labs(x = plot_param$title_x, 
         y = plot_param$title_y, 
         title = plot_param$title_plot)  +
    
    # Define y-axis start and end
    #coord_cartesian(ylim = c(0, ceiling(max(data[,plot_param$data_y])/10)*10)) +
    
    # Define the color of the violins
    ggplot2::scale_fill_manual(values = my_palette) +
    
    # Adjust font size, style
    my_theme 
  
  if (log_scale_y == "TRUE"){
    p <- p +
      scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000)) 	#display y axis in log scale 
  }        
  
  if (show_intercept_y == "TRUE"){
    p <- p +
      geom_hline(yintercept = yintercept_cutoff, linetype = 2)
  }
  
  # Save the plot
  ggplot2::ggsave(filename = paste0("Box_plot_", file_suffix, ".pdf"),
                  #plot = last_plot(),
                  plot = p,
                  device = "pdf",
                  path = results_path,
                  scale = 1,
                  width = 2+length(unique(data[,plot_param$data_x])),
                  height = (2+length(unique(data[,plot_param$data_x])))*aspect_ratio,
                  units = c("in"),
                  dpi = 600,
                  limitsize = TRUE,
                  bg = NULL)
  
}

#******************************************************************************#
#                                HISTOGRAM PLOT                                #
#******************************************************************************#

# Function to plot histogram
plot_histogram <- function(data, plot_param){
  
  ggplot(data = data, 
         aes(x = !!rlang::sym(plot_param$data_x), 
             color = !!rlang::sym(plot_param$data_color),
             fill = !!rlang::sym(plot_param$data_fill))) +
    
    # Plot box plot
    geom_histogram(binwidth = bin_width) +
    
    # Define the theme of plot
    theme_classic() + 
    
    # Define the axis, plot headings
    labs(x = plot_param$title_x, 
         y = plot_param$title_y, 
         title = plot_param$title_plot)  +
    
    # Define y-axis start and end
    #coord_cartesian(ylim = c(0, ceiling(max(data[,plot_param$data_y])/10)*10)) +
    
    # Adjust font size, style
    my_theme 
  
  # Save the plot
  ggplot2::ggsave(filename = paste0("Histogram_plot_", file_suffix, ".pdf"),
                  #plot = last_plot(),
                  plot = p,
                  device = "pdf",
                  path = results_path,
                  scale = 1,
                  width = 2+length(unique(data[,plot_param$data_x])),
                  height = (2+length(unique(data[,plot_param$data_x])))*aspect_ratio,
                  units = c("in"),
                  dpi = 600,
                  limitsize = TRUE,
                  bg = NULL)
}

#******************************************************************************#
#                                 VOLCANO PLOT                                 #
#******************************************************************************#

# Function to calculate pval and log2FC
calc_stats <- function(data_mat, metadata, file_suffix){
  
  genes <- c()
  expt <- c()
  control <- c()
  pval <- c()
  
  for (j in 1:nrow(data_mat)){
    
    data <- data_mat[j,] %>% 
      t() %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "Sample") %>%
      dplyr::left_join(metadata, by = c("Sample" = "Sample"))
    colnames(data) <- c("Sample", "values", "Condition")
    
    # If all values for Reference == NA or Target == NA, set them to 0. 
    # If at least 1 value is not NA, do not make any changes. 
    if (sum(data[data$Condition == Reference, ]$values, na.rm=TRUE) == 0){
      data[data$Condition == Reference, ]$values <- rep(x=0, times=sum(metadata$Condition == Reference))
    }
    if (sum(data[data$Condition == Target, ]$values, na.rm=TRUE) == 0){
      data[data$Condition == Target, ]$values <- rep(x=0, times=sum(metadata$Condition == Target))
    }
    
    # Remove rows which have NA in values.
    # NOTE: This step not needed as f.test/t.test automatically removes entries
    # which have NA
    #data <- data[!is.na(data$values),]
    
    # Use F.test() to determine if variances are equal or not.
    # NOTE: p < 0.05 => unequal variance
    # NOTE: If 3 values for Iso are "NA","NA","18" and 3 values for IP are 
    # "17", "NA", "18", var.test and t.test with throw error since Iso has only
    # 1 value.
    # NOTE: If 3 values for Iso are 1,1,1 and 3 values for IP are 1,1,1, t.test 
    # will throw error "data are essentially constant".
    # NOTE: If 3 values for Iso are 0,0,0 and 3 values for IP are 1,1,1, t.test 
    # will throw error "data are essentially constant".
    
    if (sum(!is.na(data[data$Condition == Reference, ]$values)) > 1 &
        sum(!is.na(data[data$Condition == Target, ]$values)) > 1 & 
        length(unique(data$values)) > 1 &
        (length(unique(data[data$Condition == Reference, ]$values)) + 
         length(unique(data[data$Condition == Target, ]$values)) > 2)){
      #   f_test <- var.test(formula = values ~ Condition, 
      #                      data = data,
      #                      alternative = "two.sided")
      #   
      #   if(!is.na(f_test$p.value)){
      # if (f_test$p.value < 0.05){
      # t_test <- t.test(formula = values ~ Condition,
      #                  data = data,
      #                  alternative = "two.sided",
      #                  var.equal = FALSE)
      # }
      
      # Calculate p values, mean expression
      t_test <- t.test(formula = values ~ Condition, 
                       data = data,
                       alternative = "two.sided",
                       var.equal = TRUE)
      
      if (grepl(Reference, names(t_test$estimate[1]))){
        genes <- c(genes, rownames(data_mat)[j])
        pval <- c(pval, t_test$p.value)
        control <- c(control, t_test$estimate[[1]])
        expt <- c(expt, t_test$estimate[[2]])
      } else if (grepl(Reference, names(t_test$estimate[2]))) {
        genes <- c(genes, rownames(data_mat)[j])
        pval <- c(pval, t_test$p.value)
        control <- c(control, t_test$estimate[[2]])
        expt <- c(expt, t_test$estimate[[1]])
      }
    } else{
      genes <- c(genes, rownames(data_mat)[j])
      pval <- c(pval, 1) # Note: DO NOT SET to NA. It will increase padj.
      control <- c(control, mean(data[data$Condition == Reference, ]$values, na.rm=TRUE))
      expt <- c(expt, mean(data[data$Condition == Target, ]$values, na.rm=TRUE))
    }
  }
  
  stats_df <- data.frame(genes, expt, control, pval)
  stats_df$padj <- stats::p.adjust(p = stats_df$pval, method = "fdr", n = length(stats_df$pval))
  stats_df$log2FC <- stats_df$expt - stats_df$control
  
  result <- data_mat %>% 
    tibble::rownames_to_column(var = "Gene") %>% 
    dplyr::left_join(stats_df, by = c("Gene" = "genes"))
  
  # Save the results
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "1")
  openxlsx::writeData(wb, sheet = "1", x = result, rowNames = FALSE)
  openxlsx::saveWorkbook(wb, 
                         file = paste0(results_path, "Volcano_Results_", file_suffix, ".xlsx"), 
                         overwrite = TRUE)
  
  return(result)
}

# Function to plot volcano plots
plot_volcano <- function(volcano_df, disp_genes, file_suffix){
  
  # Categorize the data points
  volcano_df <- volcano_df %>%
    dplyr::mutate(Direction = dplyr::case_when(padj < padj_cutoff & log2FC > log2_cutoff ~ paste0("Up in ", Source, "_", Target),
                                               padj < padj_cutoff & log2FC < -log2_cutoff ~ paste0("Up in ", Source, "_", Reference),
                                               TRUE ~ "Not Significant"),
                  padj = dplyr::case_when(is.na(padj) ~ 0,
                                          TRUE ~ padj),
                  Significance = dplyr::case_when(abs(log2FC) >= log2_cutoff & padj <= 0.05 & padj > 0.01 ~ "FDR < 0.05",
                                                  abs(log2FC) >= log2_cutoff & padj <= 0.01 & padj > 0.001 ~ "FDR < 0.01",
                                                  abs(log2FC) >= log2_cutoff & padj <= 0.001  ~ "FDR < 0.001",
                                                  TRUE ~ "Not Significant"))
  
  # Define the limits of the x-axis and y-axis
  x_limits <- as.vector(quantile(volcano_df$log2FC, probs = c(0,1)))
  if (is.infinite(max(x_limits))){
    x_limits <- as.vector(quantile(volcano_df$log2FC, probs = c(0,0.999)))
  }
  y_limits <- as.vector(quantile(-log10(volcano_df$padj), na.rm = TRUE, probs = c(0,1)))
  if (is.infinite(max(y_limits))){
    y_limits <- as.vector(quantile(-log10(volcano_df$padj), na.rm = TRUE, probs = c(0,0.999)))
  }
  
  # NOTE: If using labels, sort labels in alphabetical order and then assign 
  # color because R by default will arrange the labels in alphabetical order 
  # first and then match them to colors indicated in values vector and then 
  # color the plot. The coloring in the legend is however dependent on the 
  # order of labels vector and values vector. To understand, create plot first 
  # using the sort and then without the sort(). 
  
  # NOTE: DO NOT USE labels for defining colors due to reasons above. 
  # RECOMMEND using a named vector 
  #volcano_palette <- c("#808080", RColorBrewer::brewer.pal(11, "RdYlBu")[c(11)], RColorBrewer::brewer.pal(11, "RdYlBu")[c(1)])
  if (color_by == "Significance"){
    volcano_palette <- c(viridis(4)[4], viridis(4)[3], viridis(4)[2], viridis(4)[1])
    names(volcano_palette) <- c("Not Significant", "FDR < 0.05", "FDR < 0.01", "FDR < 0.001")
  } else if (color_by == "Direction" & same_color == "TRUE"){
    x <- unique(volcano_df$Direction)
    volcano_palette <- dplyr::case_when(grepl(Target, x) ~ "orange", 
                                        grepl(Reference, x) ~ "purple", 
                                        TRUE ~ "grey")
    names(volcano_palette) <- unique(volcano_df$Direction)
  } else if (color_by == "Direction" & same_color == "FALSE"){
    x <- unique(volcano_df$Direction)
    volcano_palette <- c("grey", my_palette[1:(length(x)-1)])
    names(volcano_palette) <- unique(volcano_df$Direction)
  }
  
  ggplot2::ggplot(data = volcano_df, 
                  aes(x = log2FC, 
                      y = -log10(padj), 
                      label = Gene,
                      color = get(color_by), 
                      shape = Direction,
                      size = 2)) +
    
    # Plot dot plot
    ggplot2::geom_point() +
    
    # Define the theme of plot
    ggplot2::theme_classic() +
    
    # Define the axis, plot headings
    ggplot2::labs(x = expression("log"[2]*"FC"),
                  y = expression("-log"[10]*"padj"),
                  color = plot_param$title_legend_color,
                  size = plot_param$title_legend_size,
                  shape = "Direction",
                  color = color_by,
                  title = paste0(Target, " vs ", Reference)) +   
    
    # Draw line to mark the cutoffs
    geom_vline(xintercept = c(-log2_cutoff,log2_cutoff), color = "black", linetype = "dotted", linewidth = 0.5) +
    geom_hline(yintercept = -log10(padj_cutoff), color = "black", linetype = "dotted", linewidth = 0.5) +
    
    # Define the axis tick marks
    scale_x_continuous(breaks = seq(-ceiling(max(abs(x_limits))), ceiling(max(abs(x_limits))), 1)) +
    #scale_y_continuous(breaks = seq(0, ceiling(y_limits[2]/10)*10, 20)) +
    
    # Define x-axis start and end
    coord_cartesian(xlim = c(-ceiling(max(abs(x_limits))), ceiling(max(abs(x_limits))))) +
    
    # Adjust size of symbols in legend. Since, legend is based on color, we use color = 
    guides(colour = guide_legend(override.aes = list(size = 3)),
           shape = guide_legend(override.aes = list(size = 3))) +
    
    # Define the color of the dots
    #scale_color_viridis_d()+
    scale_color_manual(values = volcano_palette) +
    
    # Adjust font size, style
    my_theme +
    
    # Add gene labels
    geom_text_repel(data = volcano_df %>% dplyr::filter(Gene %in% disp_genes, padj < padj_cutoff),
                    mapping = aes(label = Gene),
                    #size = 2,
                    force = 0.5,
                    point.size = 1,
                    angle = 0,
                    #vjust = 0,
                    #hjust = 0,
                    #direction = "y",
                    box.padding = 1,  # increases line length somehow
                    point.padding = 0.1,
                    max.overlaps = Inf,
                    xlim = c(NA, NA),
                    ylim = c(-Inf,NA),
                    min.segment.length = dplyr::if_else(draw_line == "TRUE", 0.2, Inf),
                    #arrow = arrow(length = unit(0.015, "npc")),
                    position = position_quasirandom())
  
  # Save the plot
  ggplot2::ggsave(filename =  paste0("Volcano_Plot",  file_suffix, ".pdf"),
                  plot = last_plot(),
                  device = "pdf",
                  path = results_path,
                  scale = 1,
                  #width = 8.5,
                  #height = 9,
                  units = c("in"),	 
                  dpi = 300,
                  limitsize = TRUE,
                  bg = "white")
}

#******************************************************************************#
#                                 HEATMAP PLOT                                 #
#******************************************************************************#

# Function to plot heatmap
# NOTE: First column of normalized_counts MUST be "SYMBOL"
# columns parameter defines which column of metadata will be plotted as columns
plot_heatmap <- function(normalized_counts, metadata_column, metadata_row, plot_genes, disp_genes, file_suffix){
  
  #****************************************************************************#
  # Format the matrix for heatmap
  #****************************************************************************#
  
  mat <- normalized_counts %>%
    # Keep only genes that need to be plotted
    dplyr::filter(str_to_upper(SYMBOL) %in% str_to_upper(plot_genes)) %>%
    # If there are duplicated genes, keep only data for highest expressing copy
    dplyr::mutate(n = rowSums(.[,-1])) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::slice_max(n) %>%
    dplyr::ungroup() %>%
    # Duplicated genes with 0 expression in all samples still remain, remove them
    dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
    dplyr::select(everything(), -n) %>%
    # Move gene names to rownames
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("SYMBOL") %>%
    # Make sure all values are numerics
    dplyr::mutate(across(.cols = everything(), .fns = as.numeric))
  # mat[, unlist(lapply(mat, is.numeric))]    #alternative to mutate(across())
  
  # Perform log transform if needed. Count data is usually skewed right or left.
  # So, it will mostly be red or blue. log transform to make it less skewed.
  if (already_log == FALSE){
    mat <- log(1+mat, base = 2)
  }
  
  # Perform scaling for each gene across samples/cells if needed. Since scale() 
  # performs only column scaling, we transpose the dataframe first, so that 
  # genes are on columns & cells are on rows & then perform scaling.
  if (already_scaled == FALSE){
    mat <- mat %>% t() %>% scale() %>% t()
  }
  
  # Replace NA values with 0
  mat[is.na(mat)] <- 0
  
  # Remove rows which have 0 in all samples
  mat <- mat[rowSums(mat) != 0,]
  
  # Keep only genes that need to be plotted
  mat <- mat[intersect(plot_genes, rownames(mat)), ]
  
  # Keep ONLY samples common in metadata_column and mat
  metadata_column <- metadata_column %>% 
    dplyr::filter(make.names(get(columns)) %in% make.names(colnames(mat)))
  
  # Arrange samples in mat in the same order as in metadata_column. 
  # NOTE: This is important because in the next step we assign rownames to 
  # col_annotation assuming identical sample order between metadata_column & mat
  mat <- mat[, metadata_column[,columns]]
  
  #****************************************************************************#
  # Define column and row annotations
  #****************************************************************************#
  
  if (gtools::invalid(anno_columns)){
    col_annotation <- NA
  } else {
    col_annotation <- metadata_column %>% dplyr::select(all_of(anno_columns))
    rownames(col_annotation) <- colnames(mat)
  }
  
  # Define row annotation for genes
  if (gtools::invalid(anno_rows)){
    row_annotation <- NA
  } else {
    row_annotation <- dplyr::left_join(mat %>% as.data.frame() %>% tibble::rownames_to_column("SYMBOL"),
                                       metadata_row,
                                       by=c("SYMBOL"="SYMBOL")) %>%
      dplyr::select(everything(), -colnames(mat)) %>%
      tibble::column_to_rownames("SYMBOL")
  }
  
  #****************************************************************************#
  # Define colors column annotation
  #****************************************************************************#
  
  groups <- unique(unlist(col_annotation, use.names=FALSE))
  
  colors <- c("#BF812D", "#35978F", "#C51B7D", "#7FBC41", "#762A83",
              "#E08214", "#542788", "#D6604D", "#4393C3", "#878787",
              "#1A1A1A", "#FFFFBF", "#377EB8", "#4DAF4A", "#984EA3",
              "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999",
              "#66C2A5", "#FC8D62", "#000000", "#9E0142", "#E41A1C")
  #colors <- c("#74006F", "#FFC606")    # Female: Purple, Male:Gold
  #colors <- c("#F6D2E0", "#C8E7F5")    # Female: Pink,   Male:Blue (BBN paper)
  colors <- colors[1:length(groups)]
  names(colors) <- groups
  
  ann_colors <- list(colors[1:lengths(lapply(col_annotation, unique))[[1]]])
  
  x <- 2
  while (x <= ncol(col_annotation)){ 
    ann_colors <- c(ann_colors, 
                    list(colors[(sum(lengths(ann_colors))+1):(sum(lengths(ann_colors))+lengths(lapply(col_annotation, unique))[[x]])]))
    x <- x+1
  }
  names(ann_colors) <- colnames(col_annotation)
  
  # $Sample
  # FB1       FB2       FB3       FB4       FB5        FC       MB1     
  # "#BF812D" "#35978F" "#C51B7D" "#7FBC41" "#762A83" "#E08214" "#542788"
  # $Sex
  # Female      Male
  # "#9E0142" "#E41A1C"
  # $Condition
  # Tumor    Normal
  # "#377EB8" "#4DAF4A"
  
  #****************************************************************************#
  # Define how samples will be arranged in the heatmap
  #****************************************************************************#
  
  if (row_clustering_alphabetical == TRUE){
    mat <- mat[sort(rownames(mat)),]
  } 
  if (col_clustering_alphabetical == TRUE){
    mat <- mat[,sort(colnames(mat))]
  } 
  #****************************************************************************#
  # Determine breaks for heatmap color scale.
  # breaks correspond to numerical ranges for color palette's bins 
  # i.e. 0 to length(my_palette)
  #****************************************************************************#
  
  if(max(mat) == 0){
    breaks <- c(seq(from = min(mat), to = 0, length.out = 100))
    my_palette <- my_palette[1:100]
  } else if (min(mat) == 0){
    breaks <- c(seq(from = 0, to = max(mat), length.out = 100))
    my_palette <- my_palette[100:200]
  } else if(min(mat) < -3 | max(mat) > 3){
    breaks <- c(seq(-3, 0, length.out = 50), seq(3/100, 3, length.out = 50))
  } else{
    breaks <- c(seq(from = min(mat), to = 0, length.out = 50), seq(from = max(mat)/100, to = max(mat), length.out = 50))
  }
  
  #****************************************************************************#
  # Define vectors indicating positions where you want to have gaps in heatmap
  #****************************************************************************#
  
  if (gaps_in_row == TRUE){
    
    # Count() automatically arranges by alphabetical order unfortunately
    # So, we do this manually.
    element_names <- c()
    element_counts <- c()
    c <- 0
    for (i in 1:nrow(row_annotation)){
      if (!(row_annotation[,anno_rows][i] %in% element_names)){
        element_names <- c(element_names, row_annotation[,anno_rows][i])
        element_counts <- c(element_counts, c)
        c <- 1
      } else{
        c <- c+1
      }
    }
    element_counts <- c(element_counts,c)
    
    gaps_row <- data.frame("Description" = element_names, n = element_counts[-1]) %>%
      dplyr::mutate(n = cumsum(n)) %>%
      dplyr::select(n) %>%
      unlist(use.names = FALSE)
  } else{
    gaps_row <- NULL
  }
  
  if (gaps_in_col == TRUE){
    
    # Count() automatically arranges by alphabetical order unfortunately
    # So, we do this manually.
    element_names <- c()
    element_counts <- c()
    c <- 0
    for (i in 1:nrow(col_annotation)){
      if (!(col_annotation[,anno_columns][i] %in% element_names)){
        element_names <- c(element_names, col_annotation[,anno_columns][i])
        element_counts <- c(element_counts, c)
        c <- 1
      } else{
        c <- c+1
      }
    }
    element_counts <- c(element_counts,c)
    
    gaps_col <- data.frame("Description" = element_names, n = element_counts[-1]) %>%
      dplyr::mutate(n = cumsum(n)) %>%
      dplyr::select(n) %>%
      unlist(use.names = FALSE)
  } else {
    gaps_col <- NULL
  }
  
  #****************************************************************************#
  # Perform row and column clustering
  #****************************************************************************#
  
  if(row_clustering == TRUE){
    # cluster and re-order rows
    rowclust <- hclust(dist(mat))
    reordered <- mat[rowclust$order,]
  } else{
    reordered <- mat
  }
  if(col_clustering == TRUE){
    # cluster and re-order columns
    colclust <- hclust(dist(t(mat)))
    reordered <- reordered[, colclust$order]
  } else{
    reordered <- reordered
  }
  
  #****************************************************************************#
  # List genes and samples you want to display in the plot
  #****************************************************************************#
  
  display_col <- colnames(reordered)
  display_row <- data.frame("gene" = rownames(reordered)) %>%
    dplyr::mutate(gene = dplyr::case_when(gene %in% disp_genes ~ gene, 
                                          TRUE ~ "")) %>% 
    unlist(use.names=FALSE)
  
  # Save the clustered scores in xlsx
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "Heatmap_matrix")
  openxlsx::writeData(wb, sheet = "Heatmap_matrix", x = reordered, rowNames = TRUE)
  openxlsx::saveWorkbook(wb, file = paste0(results_path, "Heatmap_matrix", file_suffix, ".xlsx"), 
                         overwrite = TRUE)
  
  #****************************************************************************#
  # Plot heatmap
  #****************************************************************************#
  pheatmap::pheatmap(mat = as.matrix(reordered),
                     color = my_palette,
                     breaks = breaks, 
                     border_color = "grey90", #"white"
                     cellwidth = dplyr::case_when(ncol(reordered) > 25 ~ 1,
                                                  ncol(reordered) <= 25 ~ 6), 
                     cellheight = dplyr::case_when(nrow(reordered) > 25 ~ 1,
                                                   nrow(reordered) <= 25 ~ 6), 
                     scale = "none",   
                     cluster_rows = row_clustering,   #cluster the rows
                     cluster_cols = col_clustering,   #cluster the columns
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",
                     clustering_method = "complete",
                     legend = TRUE, 
                     legend_breaks = NA,
                     legend_labels = NA, 
                     annotation_row = row_annotation,
                     annotation_col = col_annotation,
                     annotation_colors = ann_colors,
                     annotation_legend = TRUE,
                     annotation_names_row = FALSE,
                     annotation_names_col = TRUE,
                     show_rownames = dplyr::if_else(length(disp_genes) < 80, TRUE, FALSE, missing = NULL), 
                     show_colnames = dplyr::if_else(length(unique(display_col)) < 50, TRUE, FALSE, missing = NULL),
                     fontsize = 5, 
                     fontsize_row = 5, 
                     fontsize_col = 5,
                     gaps_row = gaps_row,
                     gaps_col = gaps_col,
                     angle_col = "45",
                     fontsize_number = 0.8*fontsize, 
                     labels_row = display_row, 
                     labels_col = display_col,
                     width = 8.5,
                     height = 11,
                     filename = paste0(results_path, "Heatmap", file_suffix, ".jpg"))
}

# Note: When you use autoplot, each dot represents the rows of read_data.
# So, if plotting read_data from DESeq2, use tranposed read_data like t(read_data)
# Also, remove genes which have zero expression in all samples to avoid error messages
# Also, use normalized counts not raw counts to replicate PCA plot from DESeq2

iris.pca <- prcomp(data, 
                   center = TRUE, 
                   scale. = TRUE) 


library(ggfortify) 
iris.pca.plot <- autoplot(iris.pca, label=TRUE)


# This R file contains all the user defined functions which can be imported
# to analyze bulk RNA Seq, single cell RNA Seq, plot graphs etc

# Assign colors for UMAP. The current palette supports up to 33 cell types
my_palette <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", 
                "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
                "#FFC61E", "#762A83", "#333333", "#FF1F5B", "#B8E80C",
                "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5",
                "#C49C94", "#F7B6D2", "#C7C7C7", "#DBDB8D", "#9EDAE5",
                "#FFFFBF", "#C51B7D", "#DE77AE", "#F1B6DA", "#FDE0EF", 
                "#F7F7F7", "#E6F5D0", "#B8E186")

#******************************************************************************#
#                                   BAR PLOT                                   #
#******************************************************************************#

# Function to plot bar chart
plot_bar <- function(data, plot_param){
  
  ggplot2::ggplot(data = data,
                  aes(x = !!rlang::sym(plot_param$data_x), 
                      y = reorder(stringr::str_to_upper(string = !!rlang::sym(plot_param$data_y)), desc(!!rlang::sym(plot_param$data_fill))),
                      fill = !!rlang::sym(plot_param$data_fill))) +
    
    # Plot bar plot
    ggplot2::geom_col(width = 0.75) +
    
    # Define the theme of plot
    ggplot2::theme_classic() +
    
    # Define the axis, plot headings
    ggplot2::labs(x = plot_param$title_x,
                  y = plot_param$title_y,
                  fill = plot_param$title_legend_fill,
                  color = plot_param$title_legend_color,
                  size = plot_param$title_legend_size,
                  title = plot_param$title_plot) +
    
    # Define x-axis start and end
    ggplot2::coord_cartesian(xlim = c(0, ceiling(max(abs(data[,plot_param$data_x]))))) +
    
    # Define the color of the bars
    viridis::scale_fill_viridis(discrete = FALSE, option = "D") +
    
    # Adjust font size, style
    my_theme +
    ggplot2::theme(axis.text.y = element_text(family = "sans", face = "plain", colour = "black", size = 10, hjust = 1))
  
  # Save the plot
  ggplot2::ggsave(filename = paste0("Bar_plot_", file_suffix, ".pdf"),
                  plot = last_plot(),
                  device = "pdf",
                  path = results_path,
                  scale = 1,
                  width = 2+nrow(data),
                  height = (2+nrow(data))*aspect_ratio,
                  units = c("in"),
                  dpi = 600,
                  limitsize = TRUE,
                  bg = NULL)
  
  # Create a workbook to store the results
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb = wb, sheetName = "Bar")
  openxlsx::writeData(wb = wb, sheet = "Bar", x = data)
  openxlsx::saveWorkbook(wb, file = paste0(results_path, "Bar_plot_", file_suffix, "_Results.xlsx"),
                         overwrite = TRUE)
  
}

#******************************************************************************#
#                                   DOT PLOT                                   #
#******************************************************************************#

# Function to plot dot plot
plot_dot <- function(data, plot_param){
  
  ggplot2::ggplot(data = data,
                  aes(x = !!rlang::sym(plot_param$data_x), 
                      y = reorder(stringr::str_to_upper(string = !!rlang::sym(plot_param$data_y)), desc(!!rlang::sym(plot_param$data_fill))),
                      size = !!rlang::sym(plot_param$data_size),
                      color = !!rlang::sym(plot_param$data_color))) +
    
    # Plot dot plot
    ggplot2::geom_point() +
    
    # Define the theme of plot
    ggplot2::theme_classic() +
    
    # Define the axis, plot headings
    ggplot2::labs(x = plot_param$title_x,
                  y = plot_param$title_y,
                  color = plot_param$title_legend_color,
                  size = plot_param$title_legend_size,
                  title = plot_param$title_plot) +
    
    # Define x-axis start and end
    ggplot2::coord_cartesian(xlim = c(0, ceiling(max(abs(data[,plot_param$data_x]))))) +
    
    # Define the color of the dots
    viridis::scale_color_viridis(discrete = FALSE, option = "D") +
    
    # Adjust font size, style
    my_theme +
    ggplot2::theme(axis.text.y = element_text(family = "sans", face = "plain", colour = "black", size = 10, hjust = 1))
  
  #ggplot2::scale_size(breaks = sapply(as.vector(quantile(c(min(data$Count), max(data$Count)))), floor))
  
  # Save the plot
  ggplot2::ggsave(filename = paste0("Dot_plot_", file_suffix, ".pdf"),
                  plot = last_plot(),
                  device = "pdf",
                  path = results_path,
                  scale = 1,
                  width = 2+nrow(data),
                  height = (2+nrow(data))*aspect_ratio,
                  units = c("in"),
                  dpi = 600,
                  limitsize = TRUE,
                  bg = NULL)
  
  # Create a workbook to store the results
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb = wb, sheetName = "Dot")
  openxlsx::writeData(wb = wb, sheet = "Dot", x = data)
  openxlsx::saveWorkbook(wb, file = paste0(results_path, "Dot_plot_", file_suffix, "_Results.xlsx"),
                         overwrite = TRUE)
}

#******************************************************************************#
#                               STACKED BAR PLOT                               #
#******************************************************************************#

# Function to plot stacked bar chart
plot_stackedbar <- function(data, plot_param, label_percent){
  
  if (already_percent == FALSE){
    # Calculate percent of each sub type for each sample
    data <- data %>%
      data.frame() %>%
      base::replace(is.na(.), 0) %>%
      dplyr::mutate(across(.cols = c(everything(), -plot_param$data_x), .fns = as.numeric)) %>%
      dplyr::mutate(across(.cols = c(everything(), -plot_param$data_x), .fns = function(x) x*100/rowSums(x = select_if(., is.numeric)))) %>%
      tidyr::pivot_longer(cols = !rlang::sym(plot_param$data_x), names_to = plot_param$data_fill, values_to = "Percent") %>%
      dplyr::mutate(n_gene = gsub(pattern="X", replacement="", x=n_gene))
  }
  
  # Plot stacked bar chart
  p <- ggplot2::ggplot(data = data, 
                       aes(x = !!rlang::sym(plot_param$data_x), 
                           y = Percent, 
                           fill = !!rlang::sym(plot_param$data_fill))) +
    
    # Plot stacked bar plot
    ggplot2::geom_col(position = "stack", width = 0.95) +
    
    # Define the theme of plot
    ggplot2::theme_classic() + 
    
    # Define the axis, plot headings
    ggplot2::labs(x = plot_param$title_x,
                  y = plot_param$title_y,
                  fill = plot_param$title_legend_fill,
                  title = plot_param$title_plot) +
    
    # Define the color of the bars
    ggplot2::scale_fill_manual(values = rep(my_palette, times = ceiling(length(get(plot_param$title_legend_fill))/length(my_palette)))) +
    
    # Adjust font size, style
    my_theme
  
  if(label_percent == "TRUE"){
    p <- p +
      ggplot2::geom_text(aes(x = !!rlang::sym(plot_param$data_x), label = round(Percent,1)), 
                         position = position_stack(vjust = 0.5), 
                         fontface = "plain", colour = "white", size = 3, 
                         check_overlap = TRUE)
  }
  
  # Save the plot
  ggplot2::ggsave(filename = paste0("Stacked_Bar_", file_suffix, ".pdf"),
                  #plot = last_plot(),
                  plot = p,
                  device = "pdf",
                  path = results_path,
                  scale = 1,
                  width = dplyr::if_else((2+nrow(data)) < 10, (2+nrow(data)), 10),
                  height = dplyr::if_else((2+nrow(data))*aspect_ratio < 10, (2+nrow(data))*aspect_ratio, 10),
                  units = c("in"),
                  dpi = 600,
                  limitsize = TRUE,
                  bg = NULL)
  
  # Create a workbook to store the results
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb = wb, sheetName = "Stacked_Bar")
  openxlsx::writeData(wb = wb, sheet = "Stacked_Bar", x = data)
  openxlsx::saveWorkbook(wb, file = paste0(results_path, "Stacked_Bar_", file_suffix, "_Results.xlsx"),
                         overwrite = TRUE)
  
}

#******************************************************************************#
#                                 VIOLIN PLOT                                  #
#******************************************************************************#

# Function to plot violin plot
plot_violin <- function(data, plot_param){
  
  p <- ggplot(data = data, 
              aes(x = !!rlang::sym(plot_param$data_x), 
                  y = !!rlang::sym(plot_param$data_y), 
                  fill = !!rlang::sym(plot_param$data_fill))) +
    
    # Plot violin plot with a small box plot within it
    geom_violin(trim=FALSE) + 
    geom_boxplot(width = 0.1) +
    
    # Define the theme of plot
    theme_classic() + 
    
    # Define the axis, plot headings
    labs(x = plot_param$title_x, 
         y = plot_param$title_y, 
         title = plot_param$title_plot)  +
    
    # Define y-axis start and end
    #coord_cartesian(ylim = c(0, ceiling(max(data[,plot_param$data_y])/10)*10)) +
    
    # Define the color of the violins
    ggplot2::scale_fill_manual(values = my_palette) +
    
    # Adjust font size, style
    my_theme 
  
  if (log_scale_y == "TRUE"){
    p <- p +
      scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000)) 	#display y axis in log scale 
  }        
  
  if (show_intercept_y == "TRUE"){
    p <- p +
      geom_hline(yintercept = yintercept_cutoff, linetype = 2)
  }
  
  
  # Save the plot
  ggplot2::ggsave(filename = paste0("Violin_plot_", file_suffix, ".pdf"),
                  #plot = last_plot(),
                  plot = p,
                  device = "pdf",
                  path = results_path,
                  scale = 1,
                  width = 2+length(unique(data[,plot_param$data_x])),
                  height = (2+length(unique(data[,plot_param$data_x])))*aspect_ratio,
                  units = c("in"),
                  dpi = 600,
                  limitsize = TRUE,
                  bg = NULL)
  
}

#******************************************************************************#
#                                   BOX PLOT                                   #
#******************************************************************************#

# Function to plot box plot
plot_box <- function(data, plot_param){
  
  p <- ggplot(data = data, 
              aes(x = !!rlang::sym(plot_param$data_x), 
                  y = !!rlang::sym(plot_param$data_y), 
                  fill = !!rlang::sym(plot_param$data_fill))) +
    
    # Plot box plot
    geom_boxplot() +
    
    # Define the theme of plot
    theme_classic() + 
    
    # Define the axis, plot headings
    labs(x = plot_param$title_x, 
         y = plot_param$title_y, 
         title = plot_param$title_plot)  +
    
    # Define y-axis start and end
    #coord_cartesian(ylim = c(0, ceiling(max(data[,plot_param$data_y])/10)*10)) +
    
    # Define the color of the violins
    ggplot2::scale_fill_manual(values = my_palette) +
    
    # Adjust font size, style
    my_theme 
  
  if (log_scale_y == "TRUE"){
    p <- p +
      scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000)) 	#display y axis in log scale 
  }        
  
  if (show_intercept_y == "TRUE"){
    p <- p +
      geom_hline(yintercept = yintercept_cutoff, linetype = 2)
  }
  
  # Save the plot
  ggplot2::ggsave(filename = paste0("Box_plot_", file_suffix, ".pdf"),
                  #plot = last_plot(),
                  plot = p,
                  device = "pdf",
                  path = results_path,
                  scale = 1,
                  width = 2+length(unique(data[,plot_param$data_x])),
                  height = (2+length(unique(data[,plot_param$data_x])))*aspect_ratio,
                  units = c("in"),
                  dpi = 600,
                  limitsize = TRUE,
                  bg = NULL)
  
}

#******************************************************************************#
#                                HISTOGRAM PLOT                                #
#******************************************************************************#

# Function to plot histogram
plot_histogram <- function(data, plot_param){
  
  ggplot(data = data, 
         aes(x = !!rlang::sym(plot_param$data_x), 
             color = !!rlang::sym(plot_param$data_color),
             fill = !!rlang::sym(plot_param$data_fill))) +
    
    # Plot box plot
    geom_histogram(binwidth = bin_width) +
    
    # Define the theme of plot
    theme_classic() + 
    
    # Define the axis, plot headings
    labs(x = plot_param$title_x, 
         y = plot_param$title_y, 
         title = plot_param$title_plot)  +
    
    # Define y-axis start and end
    #coord_cartesian(ylim = c(0, ceiling(max(data[,plot_param$data_y])/10)*10)) +
    
    # Adjust font size, style
    my_theme 
  
  # Save the plot
  ggplot2::ggsave(filename = paste0("Histogram_plot_", file_suffix, ".pdf"),
                  #plot = last_plot(),
                  plot = p,
                  device = "pdf",
                  path = results_path,
                  scale = 1,
                  width = 2+length(unique(data[,plot_param$data_x])),
                  height = (2+length(unique(data[,plot_param$data_x])))*aspect_ratio,
                  units = c("in"),
                  dpi = 600,
                  limitsize = TRUE,
                  bg = NULL)
}


#******************************************************************************#
#                           DEFINE GLOBAL PARAMETERS                           #
#******************************************************************************#

parent_path <- "C:/Users/KailasammS/Desktop/"
results_path <- "C:/Users/KailasammS/Desktop/"

# Define any filename you want added to final file
file_suffix <- NULL

# Define cutoff for padj and log2FC
padj_cutoff <- 0.05
log2FC_cutoff <- 0.58

# Define aspect ratio
# NOTE: width of plot has been defined in plot_bar() as height*aspect_ratio.
# So, large the aspect ratio, wider the plot
aspect_ratio <- 1

# Assign colors for UMAP. The current palette supports up to 33 cell types
my_palette <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", 
                "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
                "#FFC61E", "#762A83", "#333333", "#FF1F5B", "#B8E80C",
                "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5",
                "#C49C94", "#F7B6D2", "#C7C7C7", "#DBDB8D", "#9EDAE5",
                "#FFFFBF", "#C51B7D", "#DE77AE", "#F1B6DA", "#FDE0EF", 
                "#F7F7F7", "#E6F5D0", "#B8E186")

# Define axis font etc to use in all plots
my_theme <- ggplot2::theme(aspect.ratio = aspect_ratio,
                           plot.title =   element_text(family = "sans", face = "bold",  colour = "black", size = 15, hjust = 0.5),
                           axis.title.x = element_text(family = "sans", face = "bold",  colour = "black", size = 12, hjust = 0.5, vjust = 0,   angle = 0),
                           axis.title.y = element_text(family = "sans", face = "bold",  colour = "black", size = 12, hjust = 0.5, vjust = 1,   angle = 90),
                           legend.title = element_text(family = "sans", face = "bold",  colour = "black", size = 12, hjust = 0,   vjust = 1,   angle = 0),
                           axis.text.x =  element_text(family = "sans", face = "plain", colour = "black", size = 10, hjust = 0.5, vjust = 0.5, angle = 45),
                           axis.text.y =  element_text(family = "sans", face = "plain", colour = "black", size = 10, hjust = 0.5, vjust = 0.5, angle = 0),
                           legend.text =  element_text(family = "sans", face = "plain", colour = "black", size = 10, hjust = 0.5),
                           strip.text.x = element_text(family = "sans", face = "bold",  colour = "black", size = 10, hjust = 0.5),
                           #legend.background = element_rect(fill = "lightblue", size = 0.5, linetype = "solid", colour = "darkblue"),
                           legend.position = "right",
                           legend.justification = "left",
                           legend.direction = "vertical",
                           legend.key.height = unit(0.5, 'cm'),
                           legend.key.width  = unit(0.5, 'cm'), 
                           legend.text.align = 0)

#******************************************************************************#
#                         DEFINE PARAMETERS FOR GGPLOT2                        #             
#******************************************************************************#

plot_param <- list(
  
  # column to be plotted on x axis
  "data_x" = c("Sample"),                 #"k.K", "Count, "Sample"
  
  # Define title of x axis
  "title_x" = c(),  #"Encrichment Ratio", "Number of Genes"
  
  # column to be plotted on y axis
  "data_y" = c("nUMIs"), #"Description", "nUMIs"
  
  # Define title of y axis
  "title_y" = c("nUMIs"),   #"Percent Composition", "nUMIs"
  
  # column to be used for filling bar colors
  "data_fill" = c("Sample"), #"FDR", "Subgroups", "Sample"
  
  # column to be used for coloring the dots
  "data_color" = c(), #"FDR", "Subgroups"
  
  # column to be used for determining the size of dots
  "data_size" = c(), #"Count"
  
  # Define title of legend
  "title_legend_fill" = c("Sample"),   #"log10(padj)", "Subgroups", "Sample"
  "title_legend_color" = c(),             #"log10(padj)", "Subgroups"
  "title_legend_size" = c(),              #"Number of Genes"
  
  # Define title of plot
  "title_plot" = c("Distribution of UMIs") #"Distribution of UMIs"
)

#******************************************************************************#
#                      IMPORT DATA FOR BAR CHART/DOT PLOT                      #             
#******************************************************************************#

# NOTE: Bar plots gives information only on enrichment ratio & pvalue while dot
# plots give info on number of genes in each pathway additionally.

# NOTE: data MUST have ALL of the columns you define above
# (i) data_y ~ "Description" column containing pathway names
# (ii) data_x ~ "k.K" column containing ratio or total number of genes
# (iii) data_fill/data_color ~ "FDR" column containing log10(padj) values
# (iv) data_size ~ "Count" column containing total number of genes
data <- openxlsx::read.xlsx("C:/Users/KailasammS/Desktop/Data_metascape_DDR2KO.xlsx")

# data <- data %>%
#   dplyr::filter(!!rlang::sym(plot_param$data_fdr) < log10(padj_cutoff)) %>%
#   dplyr::mutate(!!rlang::sym(plot_param$data_pathway) := stringr::str_wrap(get(plot_param$data_pathway), width = 22)) %>%
#   dplyr::slice_min(order_by = get(plot_param$data_fdr), n = 10)

data <- data %>%
  dplyr::filter(FDR < log10(padj_cutoff)) %>%
  dplyr::mutate(Description = stringr::str_wrap(Description, width = 22)) %>%
  dplyr::slice_min(order_by = FDR, n = 10)

# Plot bar chart
plot_bar(data, plot_param)

# Plot dot chart
plot_dot(data, plot_param)

#******************************************************************************#
#                      IMPORT DATA FOR STACKED BAR CHART                       #             
#******************************************************************************#

# Define if you want to label percentages within stacked bar chart
label_percent <- "TRUE"

# Define if input data is already converted to % and can be be plotted as is
already_percent <- "FALSE"

# NOTE: data MUST have ALL of the columns you define above
# (i) data_x ~ "Sample" column containing sample names. This will be on X-axis
# (ii) Rest of columns MUST be the subtype names i.e. sub-categories within each
# sample which will be stacked in each bar which will automatically be named 
# "Subgroups

# Sample  Epithelial-Basal  Epithelial-Luminal
# FB1     350               900
# FB2     700               100

data <- openxlsx::read.xlsx("C:/Users/KailasammS/Desktop/y.xlsx")

# Plot stacked bar chart
plot_stackedbar(data, plot_param, label_percent)

# celltype <- "Epithelial"

# # Load the integrated seurat object
# integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn", 
#                                     dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))

# # Remove unwanted cells and samples
# integrated_seurat <- subset(x = integrated_seurat,
#                             subset = (cell_class %in% c("Mixed", "Unclassified") | (Condition %in% c("Normal"))),
#                             invert = TRUE)

# # Calculate cells per cluster
# data <- integrated_seurat@meta.data %>% 
#   dplyr::group_by(Sample, sub_type) %>%
#   dplyr::count() %>%
#   tidyr::pivot_wider(id_cols = Sample, names_from = sub_type, values_from = n) %>%
#   dplyr::rename_with(make.names)

#******************************************************************************#
#                       IMPORT DATA FOR VIOLIN/BOX PLOT                        #             
#******************************************************************************#

# Define if y axis should be in log scale
log_scale_y <- TRUE

# Define if y axis should have an intercept
show_intercept_y <- TRUE

# Define the value for y axis intercept
yintercept_cutoff <- 5

# NOTE: data MUST have ALL of the columns you define above
# (i) data_x ~ "Sample" column containing sample names. This will be on X-axis
# (ii) data_y ~ "nUMIs" column containing y axis values.

# Sample  nUMIs   
# VC      4.2    
# VC      11.5    
# VC      7.3     
# CJ      5.8     
# CJ      6.4     

data <- openxlsx::read.xlsx("C:/Users/KailasammS/Desktop/x.xlsx")

# Plot violin plot
plot_violin(data, plot_param)

# Plot box plot
plot_box(data, plot_param)

#******************************************************************************#
#                        IMPORT DATA FOR HISTOGRAM PLOT                        #             
#******************************************************************************#

# Define bin_width for histogram
bin_width <- 1

# NOTE: data MUST have ALL of the columns you define above
# (i) data_x ~ "Sample" column containing sample names. This will be on X-axis
# (ii) data_y ~ "nUMIs" column containing y axis values.
data <- openxlsx::read.xlsx("C:/Users/KailasammS/Desktop/x.xlsx")

# Plot histogram
plot_histogram(data, plot_param)

#******************************************************************************#
#               IMPORT DATA & DEFINE PARAMETERS FOR VOLCANO PLOT               #             
#******************************************************************************#

# Import expression data without log2FC and pval
# NOTE: data_mat MUST have 
# (i) "SYMBOL" as rownames, not in 1st column
# (ii) rest of column names MUST match with metadata$Sample
# Values should be raw i.e. untransformed. DO NOT use log transformed values etc.
# DO NOT replace NA with 0 etc. Leave NA as they are.
# NA values will be imputed/replaced with average of non-NA values.
# If all values are NA, they will be set to 0.
# It also performs quantile normalization if specified for mass spec data

data_mat <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "Data_Mass_Spec.xlsx")) %>%
  tibble::column_to_rownames(colnames(.)[1])
data_mat <- log(1+data_mat, base = 2)

# Define/import metadata
# NOTE: metadata MUST have 
# (i) "Sample" as 1st column containing sample names
# (ii) "Condition" as 2nd column which MUST match with "Target" and "Reference"
# variables defined above
# (iii) "Source" as 3rd column which MUST contain the celltype

metadata <- data.frame("Sample" = colnames(data_mat), 
                       "Condition" = c(rep(x = "Control", times = 3), rep(x = "KO", times = 3)),
                       "Source" = c(rep(x = "", times = 6)))

# Define the target and reference groups
# NOTE: These values MUST be present in "Condition" column of metadata
Target <- "KO"
Reference <- "Control"

# Calculate pval and log2FC if not already calculated 
volcano_df <- calc_stats(data_mat, metadata)

# (II) Alternatively, import expression data with log2FC and pval
# NOTE: volcano_df MUST have "Gene", "padj", "log2FC", "Source" columns
# NOTE: Output of calc_stats() can be used as input for plot_volcano()
volcano_df <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "Results_metabolomics_Control_vs_KO.xlsx"))

# (III) Define any genes you want to mark in volcano plot
disp_genes <- c("Guanidinoacetate", "Phosphoserine")
disp_genes <- volcano_df %>%
  dplyr::filter(padj < 0.05 & log2FC > 0) %>%
  dplyr::select(Gene) %>%
  unlist(use.names = FALSE) %>%
  unique()

padj_cutoff <- 0.05
log2_cutoff <- 0
color_by <- "Direction"
same_color <- "FALSE"
draw_line <- "TRUE"
plot_param <- list(
  
  # column to be plotted on x axis
  "data_x" = c(),                 #"k.K", "Count, "Sample"
  
  # Define title of x axis
  "title_x" = c(),  #"Encrichment Ratio", "Number of Genes"
  
  # column to be plotted on y axis
  "data_y" = c(), #"Description", "nUMIs"
  
  # Define title of y axis
  "title_y" = c(),   #"Percent Composition", "nUMIs"
  
  # column to be used for filling bar colors
  "data_fill" = c(), #"FDR", "Subgroups", "Sample"
  
  # column to be used for coloring the dots
  "data_color" = c(), #"FDR", "Subgroups"
  
  # column to be used for determining the size of dots
  "data_size" = c(), #"Count"
  
  # Define title of legend
  "title_legend_fill" = c(),   #"log10(padj)", "Subgroups", "Sample"
  "title_legend_color" = c("log10(padj)"),             #"log10(padj)", "Subgroups"
  "title_legend_size" = c(1),              #"Number of Genes"
  
  # Define title of plot
  "title_plot" = c() #"Distribution of UMIs"
)

# Make volcano plots
plot_volcano(volcano_df, disp_genes,  file_suffix)

#******************************************************************************#
#                 IMPORT DATA & DEFINE PARAMETERS FOR HEATMAP                  #             
#******************************************************************************#

# Define if data is log transformed already
already_log <- FALSE

# Define if data is scaled already
already_scaled <- FALSE

# Define if you want to perform unsupervised row and column clustering
# NOTE: If set to FALSE, samples (columns) and genes(rows) will be arranged 
# in alphabetical order (default) in heatmap. If you want to arranged in 
# specific order, define below.
row_clustering <- TRUE    # Usually TRUE
col_clustering <- TRUE    # Usually FALSE

# Define if you want genes or samples to be arranged in alphabetical order
# NOTE: If set to FALSE, write the plot_genes in order you want in heatmap
# NOTE: If row_clustering == TRUE, then row_clustering_alphabetical is irrelevant
row_clustering_alphabetical <- FALSE
col_clustering_alphabetical <- FALSE

# Define if heatmap columns must have gaps and based on which column of metadata
gaps_in_col <- FALSE
gap_columns <- "Sample"   # Irrelevant if gaps_in_col is FALSE

# Define if heatmap rows must have gaps and based on which column of metadata_row
gaps_in_row <- FALSE
gap_rows <- "Pathway"    # Irrelevant if gaps_in_row is FALSE

# List annotations you want on heatmap
# NOTE: anno_columns MUST match one of the column names in metadata_column while
# anno_rows MUST match one of the column names in metadata_row
anno_columns <- c("Condition")
anno_rows <- c("Pathway")

# Define if you want to make color annotations of rows or columns or both
color_by_cols <- TRUE
color_by_rows <- FALSE

# Define colors for heatmap
my_palette <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100)
# my_palette <- viridis(100)
# my_palette <- c(colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100)[1:49], "#FFFFFF",
#                 colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100)[50:99])

# NOTE: normalized_counts is a dataframe with gene names in 1st column "SYMBOL"
# NOTE: normalized_counts will have genes in rows and samples in columns
# NOTE: metadata is a dataframe with sample names in 1st column "Sample"
# NOTE: metadata_row is a dataframe with gene names in 1st column "SYMBOL"

# NOTE: If there are duplicated gene names, plot_heatmap() will keep only 
# expression data for 1 copy of duplicated gene that has highest expression.
# NOTE: "Error in check.length("fill"):'gpar' element 'fill' must not be length 0"
# This error means sample_annotation dataframe doesnt match with columns of mat
# Try using mat = t(mat) to see if it fixes the error.
# NOTE: If you set scale = none, then you are plotting exact progeny scores
# NOTE: If you set scale = row, then colors do not represent progeny scores 
# (i.e. actual activity) but relative activity

# To plot heat map we need an input matrix of the following format:
#         sample1   sample2
# Gene1    x          z  
# Gene2    y          w
# The rownames of the input matrix will be gene names and colnames will be
# the sample names

# Also, you need to create an annotation matrix to pass to "annotation_col"
# parameter of pheatmap() in the following format:
#            Sex
# Sample1    Male
# Sample2    Female
# The rownames of the annotation matrix will be the sample names

# (I) Read expr data
normalized_counts <- openxlsx::read.xlsx(xlsxFile = 
                                           "C:/Users/KailasammS/Box/Saravana@cedars/10. Ongoing Projects/Boopati project/Yneg_Mass_Spec_Results.xlsx")

normalized_counts <- normalized_counts[, -1]
colnames(normalized_counts)[1] <- "SYMBOL"

# (II) Define metadata for columns if you want to annotate columns
metadata_column <- openxlsx::read.xlsx(xlsxFile = 
                                         "C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/Hany_Y/Metadata.xlsx")

metadata_column <- data.frame("Sample" = colnames(normalized_counts[,-1]), 
                              "Condition" = c("Yneg_IP", "Yneg_IP", "Yneg_iso", "Yneg_iso", "Ypos_IP", "Ypos_IP", "Ypos_iso", "Ypos_iso"))

# (III) Define metadata for rows if you want to annotate rows
metadata_row <- NULL

# (IV) Define genes to plot
plot_genes <- c("Vegfd", "Camk2b", "Mknk2", "Camk1d", "Ncoa1", "Prkd3", 
                "Pik3r5", "Crebbp", "Braf", "Ep300", "Egln3", "Nos2", "Hspa1b",
                "Vegfa", "Adm", "Ldhb", "Rras2", "Mmp2", "Mmp3")
plot_genes <- normalized_counts$SYMBOL

# (V) Genes to display in heatmap
disp_genes <- c("Aldoa", "Pgk1", "Pgam1")
disp_genes <- plot_genes

# (VI) Define which column of metadata to plot as columns in heatmap
columns <- "Sample"

file_suffix <- ""
results_path <- "C:/Users/KailasammS/Desktop/"
# Run function
plot_heatmap(normalized_counts, metadata_column, metadata_row, plot_genes, disp_genes, file_suffix, results_path)

#******************************************************************************
#                 PIE/STACKED BAR CHAR
#******************************************************************************

# Read proteome file if you already prepared it
proteome <- openxlsx::read.xlsx(xlsxFile = paste0(results_path, "Proteome/HPA_Proteome.xlsx"))

# Create list of cells corresponding to import file names
cells <- c("Epithelial_Cells", "Non-Epithelial_Cells", "Fibroblasts", "Myeloid_Cells","B_Cells", "T_Cells" )

# Plot pie chart
for (cell in cells[1:2]){
  
  # Read DEGs file
  DEGs <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "new figs BBN/DEGs_id_Male_vs_Female_", cell, ".xlsx"))
  # DEGs <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "DEGs_id_BBN_vs_Normal_Epithelial Cells.xlsx"))
  # DEGs <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "DEGs_id_BBN_vs_Normal_Fibroblasts.xlsx"))
  DEGs <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "/Fibroblast/", "DEGs_Fibroblast_proj.xlsx"))
  
  # If any of the column is missing column names, you will see
  # "Error in initialize(...):attempt to use zero-length variable name"
  DEGs <- DEGs %>% 
    dplyr::filter(padj < 0.05 & base::abs(log2FoldChange) >= 0.58) %>%
    dplyr::left_join(proteome, by = c("gene_name" = "mouse_gene_name")) %>%
    dplyr::mutate(Class = dplyr::if_else(is.na(Class), "Unknown", Class)) %>%
    dplyr::mutate(regulation = dplyr::if_else(log2FoldChange > 0, "up", "down", "NA"))
  
  # Check which of the secreted proteins are plasma proteins etc
  secreted_DEGs <- DEGs %>% 
    dplyr::filter(Class == "Secreted" | Class == "Secreted & Membrane" | Class == "Membrane") 
  
  secreted_DEGs <- dplyr::left_join(secreted_DEGs, 
                                    secreted_class %>% select("Ensembl", "X.Protein.class."),
                                    by=c("Ensembl"="Ensembl"))
  
  # # Plot pie chart for male and female
  # for (direction in c("up", "down")){
  # 
  #   # Calculate % of each cell type for sample being analyzed
  #   data <- DEGs %>%
  #     dplyr::filter(regulation == direction) %>%
  #     dplyr::count(class) %>%
  #     dplyr::mutate(percent = round(100*n/sum(n, na.rm=TRUE), digits = 0), label_percent = paste0(percent,"%")) %>%
  #     dplyr::arrange(class)
  #   
  #   # If you run ggplot without themevoid(), you will see 0 and 100% dont overlap.
  #   # We need to add the percentage value of 1st cell group  to get accurate label positions
  #   # data <- data %>%
  #   #   mutate(ypos = 100 + data$percent[1]-(cumsum(percent)-percent/2))
  #   
  #   ggplot2::ggplot(data = data, aes(x = "", y = percent, fill = class)) +
  #     geom_bar(stat = "identity", width = 3, color = "white") +
  #     coord_polar(theta = "y", start = 0, direction = -1) +
  #     #geom_label(aes(x = 1.6, y = ypos, label = label_percent), size = 5, label.size = NA, fill = NA) +
  #     #geom_label(aes(x = 2, label = label_percent), position = position_stack(vjust = 0.5), color = "black", size = 5, check_overlap = TRUE) +
  #     geom_text(aes(x=3, label=label_percent), position=position_stack(vjust=0.5), fontface="bold", colour="black", size=5, check_overlap=TRUE) +
  #     theme_void() +        #remove background, grid, numeric labels
  #     scale_fill_brewer(palette = "Set1",
  #                       aesthetics = "fill") +
  #     ggplot2::labs(title = dplyr::if_else(direction == "up", "Male", "Female"),
  #                   fill = "Protein Class",
  #                   x = "",
  #                   y = "") +
  #     ggplot2::theme(plot.title =   element_text(family="sans", face="bold",  colour="black", size=15, hjust = 0.5),
  #                    plot.caption = element_text(family="sans", face="bold",  colour="black", size=10, hjust = 0),
  #                    axis.title.x = element_text(family="sans", face="bold",  colour="black", size=15, hjust = 0.5),
  #                    axis.title.y = element_text(family="sans", face="bold",  colour="black", size=15, hjust = 0.5),
  #                    axis.text.x =  element_blank(),
  #                    axis.text.y =  element_blank(),
  #                    legend.title = element_text(family="sans", face="bold",  colour="black", size=15, hjust = 0),
  #                    legend.text =  element_text(family="sans", face="bold",  colour="black", size=12, hjust = 0.5),
  #                    strip.text.x = element_text(family="sans", face="bold",  colour="black", size=15, hjust = 0.5),
  #                    legend.position = "right",
  #                    legend.direction = "vertical",
  #                    legend.text.align = 0,
  #                    axis.line=element_blank(),
  #                    axis.ticks=element_blank())
  #   
  #   # Save the plot
  #   ggplot2::ggsave(filename = paste0("Pie_chart_Proteome_Class_", cell, "_", direction, "regulated.pdf"),
  #                   plot = last_plot(),
  #                   device = "pdf",
  #                   path = paste0(parent_path, "Proteome/"),
  #                   scale = 1,
  #                   #width = 11,
  #                   #height = 8.5,
  #                   units = c("in"),
  #                   dpi = 600,
  #                   limitsize = TRUE,
  #                   bg = NULL)
  # }
}

# Plot stacked bar chart
data <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "BBN original figs/DEGs_id_Male_vs_Female_", cells[1], ".xlsx"))
data <- data %>% dplyr::mutate(cell_type = cells[1])

for (cell in cells[2:2]){
  
  # Read DEGs file
  DEGs <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "BBN original figs/DEGs_id_Male_vs_Female_", cell, ".xlsx"))
  DEGs <- DEGs %>% dplyr::mutate(cell_type = cell)
  
  
  # Merge all data to be plotted
  data <- rbind(data, DEGs)
  
  # Calculate percentages
  data <- data %>%
    dplyr::filter(padj < 0.05 & base::abs(log2FoldChange) >= 0.58) %>%
    dplyr::left_join(proteome, by = c("gene_name" = "mouse_gene_name")) %>%
    dplyr::mutate(class = dplyr::if_else(is.na(Class), "Unknown", Class)) %>%
    dplyr::mutate(regulation = dplyr::if_else(log2FoldChange > 0, "up", "down", "NA")) %>%
    dplyr::group_by(cell_type, regulation) %>%
    dplyr::count(class) %>%
    dplyr::mutate(percent = round(100*n/sum(n, na.rm=TRUE), digits = 2), label_percent = paste0(percent,"%")) %>%
    dplyr::arrange(class) %>%
    dplyr::mutate(plot.x=paste0(stringr::str_trunc(cell_type, 2, "right", ellipsis = "") ,"_", stringr::str_trunc(regulation, 2, "right", ellipsis = "")))
  
  
  # Plot stacked bar chart
  ggplot2::ggplot(data = data, aes(x = plot.x, y = percent, fill = class)) +
    geom_bar(stat = "identity", position = "stack", width = 0.95) +
    theme_classic() + 
    scale_fill_brewer(palette = "Set1",
                      aesthetics = "fill") +
    geom_text(aes(x=plot.x, label=n), position=position_stack(vjust=0.5), fontface="bold", colour="white", size=6, check_overlap=TRUE) +
    ggplot2::labs(title = "",
                  fill = "Protein Class",
                  x = "",
                  y = "Percent composition")
  
}

#******************************************************************************
#                EXPORT METADATA FROM ADATA
#******************************************************************************
# We export metadata from adata
# We export cell names i.e. barcodes from adata
# We import expr data using SeuratDisk so that we dont have to save expr data in
# csv files

qrsh
conda activate scVelo
conda install -c anaconda openpyxl
python

import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad("/hpc/home/kailasamms/scratch/scRNASeq_GSE132042/scRNASeq_GSE132042.h5ad")

# Info on expr 
adata.X
adata.to_df()

# Info on cells
adata.obs

# Info on genes
adata.var 

# Subsetting adata
# bdata = adata[adata.obs.cell_ontology_class == "bladder urothelial cell"]

# Convert expr data to df and save it. to_excel() is too slow. So, use to_csv()
# Recommend using SeuratDisk rather to import expr data
# df = pd.DataFrame(adata.X.toarray()).transpose()
# df.index = adata.var_names
# df.columns = adata.obs_names
# df.to_csv("/hpc/home/kailasamms/Bladder_droplet_counts.csv", index=True)

# Save metadata to df
adata.obs.to_excel("/hpc/home/kailasamms/scratch/scRNASeq_GSE132042/scRNASeq_GSE132042_Meta.xlsx", index=True)

#******************************************************************************#
#                           EFFECT SIZE AND T-SCORE                            #
#******************************************************************************#

# Dataframe as input
# Column "Gene" MUST be present
# Column Gene MUST have control sgRNAs labelled as "none" and/or "safe"
calc_t_score <- function(data){
  
  # Create a dataframe of control sgRNAs
  data_control <- data %>%
    dplyr::filter(Gene %in% c("none", "safe"))
  
  median_ctrl <- median(data_control$LFC, na.rm=TRUE)
  sd_ctrl <- sd(data_control$LFC, na.rm=TRUE)
  
  # Normalize to control sgRNAs
  data <- data %>%
    dplyr::mutate(pZ = (LFC-median_ctrl)/sd_ctrl)
  
  data_control <- data_control %>%
    dplyr::mutate(pZ = (LFC-median_ctrl)/sd_ctrl)
  
  U_ctrl <- median(data_control$pZ)
  Var_ctrl <- var(data_control$pZ)
  N_ctrl <- mean((data %>% dplyr::count(Gene))$n)
  # Nctrl is the average number of sgRNAs per gene in a given screen
  
  data <- data %>%
    dplyr::group_by(Gene) %>%
    dplyr::mutate(U_gene = median(pZ),
                  Var_gene = var(pZ),
                  N_gene = n(),
                  U_ctrl = U_ctrl,
                  Var_ctrl = Var_ctrl,
                  N_ctrl = N_ctrl,
                  S_gene = (Var_gene*(N_gene-1)) + (Var_ctrl*(N_ctrl-1)),
                  t_score = (U_gene - U_ctrl)/sqrt(S_gene/N_gene + S_gene/N_ctrl),
                  Abs_t_score = abs(t_score)) %>%
    dplyr::select(Gene, U_gene, Var_gene, N_gene, U_ctrl, Var_ctrl, N_ctrl, S_gene, t_score, Abs_t_score) %>%
    dplyr::distinct_at("Gene", .keep_all = TRUE)
  
  return(data)
}

# data is a dataframe output of calc_t_score
# save_path is location to save file
plot_t_score <- function(data, save_path, suffix){
  
  y_cutoff <- sort(data$Abs_t_score, decreasing = TRUE)[100]
  xmin <- floor(min(data$U_gene))
  xmax <- ceiling(max(data$U_gene))
  ymin <- 0
  ymax <- max(data$Abs_t_score)
  
  color_breaks <- c(-20,0,20)
  p <- ggplot2::ggplot(data = data,
                       aes(x = U_gene, 
                           y = Abs_t_score,
                           size = Abs_t_score,
                           #color = pz,
                           fill = U_gene)) +
    # Plot dot plot
    ggplot2::geom_point(col="black", 
                        shape=21,
                        stroke=0.5,
                        position=position_jitter(h=0.01,w=0.01)) +
    # Define the theme of plot
    ggplot2::theme_classic() +
    ggplot2::labs(fill = "U_gene") +
    coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax), clip = "off") +
    #scale_x_continuous(breaks = seq(-5, 5, by = 1)) +
    #scale_y_continuous(breaks = seq(0, 5, by = 1)) +
    ggplot2::guides(size = "none",
                    fill = guide_colourbar(theme = theme(legend.key.width  = unit(0.75, "lines"),
                                                         legend.key.height = unit(10, "lines"),
                                                         legend.ticks = element_blank(),
                                                         legend.frame = element_rect(colour = "Black",
                                                                                     linewidth = 0.5)))) +
    # Define the color of the dots
    ggplot2::scale_fill_viridis_c(option="turbo", limits =c(-5,3))
  #geom_hline(yintercept= y_cutoff, linetype ="dotted")
  
  # scale_fill_gradientn(colors=c("#007ba7", "Black","#FFFF00"), 
  #                      limits=c(-20, 20), 
  #                      values=c(0, scales::rescale(color_breaks, from = range(color_breaks)), 1))
  #scale_fill_gradient2(low="#007ba7", mid="Black", high="Yellow", midpoint = 0, limits=c(-5, 2))
  #scale_fill_continuous_diverging(palette = "Tofino")
  
  ggsave(paste0(save_path, suffix, ".jpg"))
}

#******************************************************************************#
#                          FLATTEN COREELATION MATRIX                          #
#******************************************************************************#

flattenCorrMatrix_pmatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  df <- data.frame(row = rownames(cormat)[row(cormat)[ut]],
                   column = rownames(cormat)[col(cormat)[ut]],
                   cor  =(cormat)[ut],
                   p = pmat[ut])
  
  return(df)
}

flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  df <- data.frame(row = rownames(cormat)[row(cormat)[ut]],
                   column = rownames(cormat)[col(cormat)[ut]],
                   cor  =(cormat)[ut])
  return(df)
}

#******************************************************************************#
#
#******************************************************************************#

# Values should be raw i.e. untransformed. DO NOT use log transformed values etc.
# DO NOT replace NA with 0 etc. Leave NA as they are.
# NA values will be imputed/replaced with average of non-NA values.
# If all values are NA, they will be set to 0.
# NO duplicated genes MUST be present.
impute_with_mean <- function(raw_counts){
  
  # Replace NA with average
  for (j in 1:nrow(raw_counts)){
    
    data <- raw_counts[j,] %>% 
      t() %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "Sample") %>%
      dplyr::left_join(metadata, by = c("Sample" = "Sample"))
    colnames(data) <- c("Sample", "values", "Condition")
    
    # If all values for Reference == NA or Target == NA, set them to 0. 
    # Replace NA with average wherever possible.
    data <- data %>% 
      dplyr::mutate(values = as.numeric(values)) %>%
      dplyr::group_by(Condition) %>% 
      dplyr::mutate(average = mean(values, na.rm=TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(average = dplyr::if_else(is.na(average), 0, average),
                    values = dplyr::if_else(is.na(values), average, values))
    
    data <- data %>% 
      dplyr::select(Sample, values) %>% 
      tibble::column_to_rownames("Sample") %>% 
      t()
    
    if(all(colnames(data) == colnames(raw_counts))){
      raw_counts[j,] <- data[,colnames(raw_counts[j,])] %>% unlist(use.names=FALSE)
    }
  }
  
  imputed_counts <- raw_counts %>% 
    dplyr::mutate(across(.cols=everything(), .fns = as.numeric))
  
  return(imputed_counts)
}

#******************************************************************************#
#
#******************************************************************************#

# Perform normalization before imputation so counts can be compared across samples
# NOTE: meta_data MUST have column "Sample" that matches with column names of
# raw_counts 
# NOTE: meta_data MUST have column "Condition" that defines the groups
quantile_norm <- function(raw_counts, meta_data, quant_norm){
  
  # Perform quantile normalization
  if (quant_norm == TRUE){
    
    # https://doi.org/10.1038/s41598-020-72664-6
    # The above paper recommends to quantile normalize within each group rather
    # than whole dataset
    quant_norm_counts <- data.frame(matrix(data=NA, nrow=nrow(raw_counts), ncol=0))
    
    # Subset raw_counts for each group
    for (c in unique(meta_data$Condition)){
      
      samples <- meta_data %>% 
        dplyr::filter(Condition %in% c) %>% 
        dplyr::select(Sample) %>% 
        unlist(use.names=FALSE)
      
      counts <- raw_counts[,samples]
      counts <- as.data.frame(preprocessCore::normalize.quantiles(as.matrix(counts)))
      rownames(counts) <- rownames(raw_counts)
      colnames(counts) <- samples
      
      quant_norm_counts <- dplyr::bind_cols(quant_norm_counts, counts)
    }
  } else {
    quant_norm_counts <- raw_counts
  }
  
  return(quant_norm_counts)
}

#******************************************************************************#
#                              CALCULATE padj AND log2FC
#******************************************************************************#

# Function to calculate pval and log2FoldChange
# norm_counts is matrix with log2 transformed values or non-log transformed values
# DO NOT use log10 transformed values
calc_stats <- function(norm_counts, metadata, Target, Reference, log2_transformed_already){
  
  # Perform t.test
  SYMBOL <- c()
  expt <- c()
  control <- c()
  pval <- c()
  for (j in 1:nrow(norm_counts)){
    
    data <- norm_counts[j,] %>% 
      t() %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "Sample") %>%
      dplyr::left_join(metadata, by = c("Sample" = "Sample"))
    colnames(data) <- c("Sample", "values", "Condition")
    
    # Use F.test() to determine if variances are equal or not.
    # NOTE: p < 0.05 => unequal variance
    
    # NOTE: If Iso <- c(NA, NA, 18) and IP <- c(17, NA, 18),
    # var.test and t.test with throw error since Iso has only 1 value.
    
    # NOTE: If Iso <- c(1,1,1) and IP <- c(1,1,1),
    # t.test will throw error "data are essentially constant".
    
    # NOTE: If Iso <- c(0,0,0) and IP <- c(1,1,1),
    # t.test will throw error "data are essentially constant".
    
    if (sum(!is.na(data[data$Condition == Reference, ]$values)) > 1 &
        sum(!is.na(data[data$Condition == Target, ]$values)) > 1 & 
        (length(unique(data[data$Condition == Reference, ]$values)) + 
         length(unique(data[data$Condition == Target, ]$values)) > 2)){
      #   f_test <- var.test(formula = values ~ Condition, 
      #                      data = data,
      #                      alternative = "two.sided")
      #   
      #   if(!is.na(f_test$p.value)){
      # if (f_test$p.value < 0.05){
      # t_test <- t.test(formula = values ~ Condition,
      #                  data = data,
      #                  alternative = "two.sided",
      #                  var.equal = FALSE)
      # }
      
      # # Remove outliers
      # ref_data <- data[data$Condition == Reference, ]$values
      # low_ref <- quantile(ref_data, na.rm=TRUE)[2] - 1.5*IQR(ref_data, na.rm=TRUE)
      # high_ref <- quantile(ref_data, na.rm=TRUE)[3] + 1.5*IQR(ref_data, na.rm=TRUE)
      # 
      # target_data <- data[data$Condition == Target, ]$values
      # low_tar <- quantile(target_data, na.rm=TRUE)[2] - 1.5*IQR(target_data, na.rm=TRUE)
      # high_tar <- quantile(target_data, na.rm=TRUE)[3] + 1.5*IQR(target_data, na.rm=TRUE)
      # 
      # data <- data %>%
      #   dplyr::filter(!(Condition == Reference & (values > high_ref | values < low_ref))) %>%
      #   dplyr::filter(!(Condition == Target & (values > high_tar | values < low_tar)))
      
      
      # Calculate p values, mean expression
      t_test <- stats::t.test(formula = values ~ Condition, 
                              data = data,
                              alternative = "two.sided",
                              var.equal = FALSE)
      
      if (grepl(Reference, names(t_test$estimate[1]))){
        SYMBOL <- c(SYMBOL, rownames(norm_counts)[j])
        pval <- c(pval, t_test$p.value)
        control <- c(control, t_test$estimate[[1]])
        expt <- c(expt, t_test$estimate[[2]])
      } else if (grepl(Reference, names(t_test$estimate[2]))) {
        SYMBOL <- c(SYMBOL, rownames(norm_counts)[j])
        pval <- c(pval, t_test$p.value)
        control <- c(control, t_test$estimate[[2]])
        expt <- c(expt, t_test$estimate[[1]])
      }
    } else{
      SYMBOL <- c(SYMBOL, rownames(norm_counts)[j])
      pval <- c(pval, 1) # Note: DO NOT SET to NA. It will increase padj.
      control <- c(control, mean(data[data$Condition == Reference, ]$values, na.rm=TRUE))
      expt <- c(expt, mean(data[data$Condition == Target, ]$values, na.rm=TRUE))
    }
  }
  
  stats_df <- data.frame(SYMBOL, expt, control, pval)
  stats_df$padj <- stats::p.adjust(p = stats_df$pval, method = "fdr", n = length(stats_df$pval))
  if(log2_transformed_already){
    stats_df$log2FoldChange <- stats_df$expt - stats_df$control  # if data is already log transformed
  }else{
    stats_df$log2FoldChange <- log(stats_df$expt/stats_df$control, base=2)
  }
  
  result <- norm_counts %>% 
    tibble::rownames_to_column(var = "SYMBOL") %>% 
    dplyr::left_join(stats_df, by = c("SYMBOL" = "SYMBOL"))
  
  return(result)
}
