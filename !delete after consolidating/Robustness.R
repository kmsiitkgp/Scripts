proj <- "Hany_Y"
#proj <- "Hany_YKO"

# parent directory : directory where input files, results, etc are stored
parent_path <- paste0("C:/Users/KailasammS/Box/Saravana@cedars/10. Ongoing Projects/LOY project/Y Robustness/")
results_path <- paste0("C:/Users/KailasammS/Desktop/")

meta_data1 <- read.xlsx()

#******************************************************************************#
# Perform DESEQ2 on all possible combinations
x <- meta_data1 %>% dplyr::filter(Condition == "scr_KO_centromere") %>% dplyr::select(Sample) %>% unlist(use.names=FALSE)
y <- meta_data1 %>% dplyr::filter(Condition == "YKO_centromere") %>% dplyr::select(Sample) %>% unlist(use.names=FALSE)
count <- 0

for (m in 0:1){
  for (n in 0:1){
    x1 <- utils::combn(x=x, m=m)
    y1 <- utils::combn(x=y, m=n)
    
    for (i in 1:ncol(x1)){
      for (j in 1:ncol(y1)){
        count <- count+1
        new_x <- setdiff(x, x1[,i])
        new_y <- setdiff(y, y1[,j])
        sampless <- c(new_x, new_y)
        #cat("\n[", count, "]", sampless)
        #cat("\n[", count, "]", "\t", length(new_x), "\t", length(new_y))
        
        meta_data <- meta_data1 %>% dplyr::filter(Sample %in% sampless)
        #print(meta_data)
        label <- count
        analyze_DESeq2()
        
      }
    }
  }
}

#******************************************************************************#

# Identify frequently upregulated genes
#common <- annotations$SYMBOL
common <- c()
#files <-list.files("C:/Users/KailasammS/Box/Saravana@cedars/10. Ongoing Projects/LOY project/YKO Robustness/", full.names = TRUE)
files <-list.files("C:/Users/KailasammS/Box/Saravana@cedars/10. Ongoing Projects/LOY project/Y Robustness/", full.names = TRUE)
files <- files[grep(pattern="DEGs", x=files)]

for (f in files){
  
  test <- openxlsx::read.xlsx(f, rowNames = FALSE)
  colnames(test)[1] <- "ID"
  test <- test %>%
    dplyr::filter(padj < 0.05) %>%
    dplyr::select(SYMBOL) %>%
    unlist(use.names=FALSE)
  
  #common <- intersect(common, test)
  common <- c(common, test)
  print(length(common))
}

# Convert to dataframe and find freqeuncy
df <- data.frame("Genes"=common) %>% 
  dplyr::count(Genes) %>% 
  dplyr::arrange(desc(n))

# Calculate number of times each gene has FC>0 or FC<0 in all combinatorial comparisons
df1 <- as.data.frame(df[,1])
colnames(df1) <- "Genes"

for (f in files){
  
  test <- openxlsx::read.xlsx(f, rowNames = FALSE)
  colnames(test)[1] <- "ID"
  test <- test %>%
    dplyr::filter(SYMBOL %in% df$Genes) %>%
    dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
    dplyr::mutate(log2FoldChange = dplyr::case_when(log2FoldChange > 0 & padj < 0.05 ~ 1,
                                                    log2FoldChange < 0 & padj < 0.05 ~ -1,
                                                    TRUE ~ 0)) %>%
    dplyr::select(SYMBOL, log2FoldChange)
  
  df1 <- left_join(df1, test, by=c("Genes"="SYMBOL"))
}

df1 <- df1 %>% 
  dplyr::mutate(totals = rowSums(df1[,-1])) %>%
  dplyr::select(Genes, totals)

#yko_df <- dplyr::left_join(df, df1, by=c("Genes"="Genes"))
y_df <- dplyr::left_join(df, df1, by=c("Genes"="Genes"))

#Add padj ad log2FC from full DESEQ2 results
DEG1 <- openxlsx::read.xlsx("C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/Hany_YKO/Results__id_YKO_centromere_vs_scr_KO_centromere_DEGs.xlsx")
DEG2 <- openxlsx::read.xlsx("C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/Hany_Y/Results__id_Y_neg_vs_Y_pos_DEGs.xlsx")

yko_df <- yko_df %>% 
  dplyr::left_join(DEG1 %>% dplyr::select(SYMBOL, padj, log2FoldChange), by =c("Genes"="SYMBOL"))

y_df <- y_df %>%
  dplyr::left_join(DEG2 %>% dplyr::select(SYMBOL, padj, log2FoldChange), by =c("Genes"="SYMBOL"))

# Save data
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Frequency_YKO")
openxlsx::writeData(wb, sheet = "Frequency_YKO", x = yko_df, rowNames = TRUE)
openxlsx::addWorksheet(wb, sheetName = "Frequency_Y")
openxlsx::writeData(wb, sheet = "Frequency_Y", x = y_df, rowNames = TRUE)

openxlsx::saveWorkbook(wb,
                       file = paste0(results_path, "Frequency_new2.xlsx"),
                       overwrite = TRUE)

