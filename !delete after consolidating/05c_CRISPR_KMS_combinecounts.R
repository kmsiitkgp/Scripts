#!/usr/bin/env Rscript

proj <- "CRISPR_Lena"
libraries <- c("CDH12KO", "CDH12ACT", "ImmuneKO", "ImmuneACT")

# NOTE: All variables and functions are defined within the file below
source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")

for (library in libraries){
  
  # Check parent_path is correct
  parent_path
  count_path <- paste0(parent_path, "count_results/", library, "/")
  
  # Create arrays to store sample names, total reads, mapped reads
  total_reads <- c()
  mapped_reads <- c()
  
  samples <- list.files(count_path, full.names = FALSE)
  samples <- sample_names[grepl(pattern="csv", x=sample_names)]
  samples <- gsub(pattern="\\.csv", replacement = "", x=sample_names)
  
  # Create a  dataframe to store counts from all samples
  count_table <- read.table(file = paste0(count_path, samples[1], ".csv"), header = TRUE, sep = ",")
  count_table <- count_table[,1:5]
  
  # Populate the count table 
  for (i in 1:length(samples)){
    
    # Get total read number from sam file
    read_sam <- read.table(file = paste0(count_path, samples[i], ".sam"), header = TRUE, sep = ",")
    total_reads <- c(total_reads, (nrow(read_sam)-1))
    
    # Read the csv file
    temp_file <- read.table(file = paste0(count_path, samples[i], ".csv"), header = TRUE, sep = ",")
    
    # Get mapped read number from csv file
    mapped_reads <- c(mapped_reads, sum(temp_file$COUNT))
    
    print(total_reads)
    print(mapped_reads)
    
    # Rename the COUNT columnto sample name
    temp_file <- temp_file %>% 
      dplyr::select(sgRNA, COUNT) %>%
      dplyr::rename(!!samples[i] := "COUNT")
    
    # Merge the counts
    count_table <- count_table %>% 
      dplyr::left_join(temp_file, by=c("sgRNA"="sgRNA"))
    
  }
  
  # If first cell has "Byte-Order-Mark" (BOM), correct it
  # count_table[1,1] <- stringr::str_replace(count_table[1,1], "ï»¿","")
  
  # Generate summaries
  read_summary <- data.frame(samples, total_reads, mapped_reads)
  colnames(read_summary) <- c("Label", "Reads", "Mapped")
  read_summary <- read_summary %>%
    dplyr::mutate(Percentage = round(100*Mapped/Reads, 2))
  
  # Save the read summary as txt file
  write.table(x = read_summary, file=paste0(count_path, "kms.countsummary.txt"), sep = "\t",
              row.names = FALSE)
  
  # Save the read table as txt file
  write.table(x = count_table, file=paste0(count_path, "kms.count.txt"), sep = "\t",
              row.names = FALSE)
}