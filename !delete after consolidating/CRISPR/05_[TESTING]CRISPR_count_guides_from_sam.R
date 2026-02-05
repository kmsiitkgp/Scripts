scripts_path <- "/hpc/home/kailasamms/projects/CRISPR/"
results_path <- "/hpc/home/kailasamms/scratch/CRISPR_Jinfen/"
proj <- "CRISPR_Jinfen"

# Save counts in xlsx
wb <- openxlsx::createWorkbook()

# Loop through each method and compile counts
for (method in c("bowtie", "bowtie2")){
  
  # Store sam files with path
  files <- list.files(paste0("/hpc/home/kailasamms/scratch/CRISPR_Jinfen/alignment_results/", method, "/"), full.names = TRUE)
  files <- files[grepl(pattern="\\.sam", x=files)]
  print(files)
  
  # Create a dataframe to store the counts from sam files
  counts <- read.csv(paste0(scripts_path, "CRISPR_Jinfen.corrected.csv"), header=FALSE)
  colnames(counts) <- c("GUIDE", "SEQ", "GENE")
  
  # Loop through each sam file
  for (f in files){
    
    # Read entire sam file
    t <- read.table(f, header=FALSE, sep="\t", fill=TRUE)
    
    # Find last line containing header by searching for "@" in column 1 i.e. V1
    t <- t %>% dplyr::filter(grepl("@", V1))
    
    # Read sam file again but skip the headers
    t <- read.table(f, header=FALSE, skip = nrow(t), sep="\t", fill=TRUE)
    
    # Calculate counts for each guide
    t <- t %>% dplyr::count(V3) %>% dplyr::filter(V3 != "*")
    sample_name <- gsub("^.*/|.sam", replacement="",x=f)
    colnames(t) <- c("GUIDE", sample_name)
    
    # Merge counts from current sam file to previous sam files
    counts <- counts %>% dplyr::left_join(t, by=c("GUIDE"="GUIDE"))
  }
  
  # Remove guide sequence
  counts <- counts %>% 
    dplyr::select(everything(), -c("SEQ"))
  
  # Save to worksheet
  openxlsx::addWorksheet(wb, sheetName = paste0(method, "_counts"))
  openxlsx::writeData(wb, sheet = paste0(method, "_counts"), x = counts, rowNames = TRUE)
}

# save excel file
openxlsx::saveWorkbook(wb, file = paste0("/hpc/home/kailasamms/scratch/CRISPR_Jinfen/Summary.xlsx"), overwrite = TRUE)

