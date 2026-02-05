#!/usr/bin/env Rscript

#******************************************************************************#
#                           LOAD NECESSARY PACKAGES                            #
#******************************************************************************#

# Data wrangling packages
library("openxlsx")
library("dplyr")
library("tibble")
library("stringr")
library("purrr")
library("Matrix")

#******************************************************************************#
#                           DECLARE GLOBAL VARIABLES                           #
#******************************************************************************#

# REQUIREMENTS:
# (i) Excel file with headers Sample Name, Target Name and CT. The headers must be named EXACTLY as specified.
# (ii) The data must be in a sheet named Results within the excel file.
# NOTE: The exported excel file from Applied Biosystems QPCR instruments usually satisfy requirements (i) and  (ii)
# (iii) ONLY the excel files to be analyzed must be stored in a folder named qpcr.
# NOTE: Additional excel files if present in the folder will hinder analysis.
# (iv) Any reference gene MUST be renamed with a prefix "Ref1_" like "Ref1_Actin".
# If multiple reference genes are present, indicate as "Ref1_Actin", "Ref2_Gapdh"
# (v) Any control sample MUST be renamed with a prefix "Control1_"
# If multiple control samples are present, indicate as "Control1_NA13", "Control2_UMUC3"

setwd("C:/Users/KailasammS/Desktop/qpcr")

# Define the pattern present in all reference samples
ref_pattern <- "-R1881"

# Define the pattern present in all experimental samples
non_ref_pattern <-  "+R1881"

# Define the reference gene
ref_gene <- "ACTIN"

#******************************************************************************#
#                         IMPORT DATA FROM EXCEL FILES                         #
#******************************************************************************#  

# Create a list of all xlsx files within qpcr folder that will be analyzed
qpcr_files <- base::list.files(path = "C:/Users/KailasammS/Desktop/qpcr")

# Create an empty dataframe
file <- data.frame()

# Read data from "Results" sheet of each xlsx files stored in working directory
for (i in 1:length(qpcr_files)) {
  temp_file <- openxlsx::read.xlsx(xlsxFile = qpcr_files[1], 
                                   sheet = "Results", 
                                   colNames = FALSE)
  
  #Find the row position where the Sample Name, Target Name and CT values start
  a <- which(temp_file == "Sample Name", arr.ind = TRUE)
  b <- which(temp_file == "Target Name", arr.ind = TRUE)
  c <- which(temp_file == "CT", arr.ind = TRUE)
  
  if (a[[1, 1]] == b[[1, 1]] && a[[1, 1]] == c[[1, 1]]) {
    # Keep ONLY Sample Name, Target Name and CT columns
    # Use dplyr::slice() to remove unwanted info at top of excel file
    # Remove rows containing NA in sample name
    temp_file <- temp_file %>%
      dplyr::select(a[[1, 2]], b[[1, 2]], c[[1, 2]]) %>%
      dplyr::slice((a[[1, 1]]+1):nrow(temp_file)) %>%
      dplyr::filter(!is.na(.[1]))
  } else{
    print("Please check if results are in a proper tabular format")
  }
  
  # Append data from next file if present
  file <- dplyr::bind_rows(file, temp_file)
}

# Rename the 1st, 2nd and 3rd columns
colnames(file) <- c("Sample_Name", "Target_Name", "Ct")

#******************************************************************************#
#                               ANALYZE THE DATA                               #
#******************************************************************************#

# Identify all sample names and target names
sample_names <- file %>% 
  dplyr::distinct(Sample_Name)

target_names <- file %>% 
  dplyr::distinct(Target_Name) %>% 
  dplyr::arrange(Target_Name)

# # Get control sample, treatment sample and reference gene from user
# number_of_controls <- c()
# reference <- c()
# control <- c()
# treatment <- c()
# 
# {
#   print(target_names)
#   r <-
#     readline(prompt = "Enter a number corresponding to the reference gene from the list above (Eg:5)")
#   r <- as.integer(r)
#   #reference <- c(reference,r)
#   reference <- c(r)
#   
#   print(sample_names)
#   c <-
#     readline(prompt = "Enter a number corresponding to the control sample from the list above (Eg:5)")
#   c <- as.integer(c)
#   #control <- c(control,c)
#   control <- c(c)
#   
#   a <-
#     readline(prompt = "Do the rest of samples from the list above belong to treatment group? Type Y for Yes or N for No")
#   a <- as.character(a)
#   if (a == "Y" | a == "y") {
#     treatment <- c(1:nrow(sample_names))[-c]
#   } else if (a == "N" | a == "n") {
#     o <-
#       readline(prompt = "Please go through the list above and indicate the total number of samples in treatment group")
#     o <- as.integer(o)
#     t <-
#       readline(prompt = "Enter a number corresponding to first treatment sample from the list above")
#     t <- as.integer(t)
#     treatment <- c(treatment, t)
#     for (i in 2:o) {
#       t <-
#         readline(prompt = "Enter a number corresponding to next treatment sample from the list above")
#       t <- as.integer(t)
#       treatment <- c(treatment, t)
#     }
#   } else{
#     readline(prompt = "Invalid response. Please re-run the script :D")
#   }
# }

#Replace "Undetermined" Ct values to NA
file$Ct <- dplyr::na_if(x = file$Ct, y = "Undetermined")

# Observe that Ct values are stored as character
sapply(file, class)

# Convert all Ct values are stored as character to numeric
file$Ct <- as.numeric(file$Ct)

# Observe that Ct values are stored as character
sapply(file, class)

# Calculate QC metrics, Average Ct
qc_df <- file %>%
  dplyr::arrange(Sample_Name, Target_Name, Ct) %>%
  dplyr::group_by(Sample_Name, Target_Name) %>%
  dplyr::mutate(Group_size = n(),
                NAs = sum(is.na(Ct)),
                Ct = dplyr::if_else(NAs == Group_size, 0, Ct),
                Max_Diff = max(Ct, na.rm=TRUE) - min(Ct, na.rm=TRUE),
                SD = sd(Ct, na.rm = TRUE),
                Avg_Ct = mean(Ct, na.rm = TRUE),
                Distance = abs(Ct - Avg_Ct),
                Outlier = dplyr::if_else((Max_Diff > 0.5 & Distance == max(Distance, na.rm=TRUE)), TRUE, FALSE), 
                Adj_Ct = dplyr::if_else(Outlier == TRUE, NA, Ct),
                Adj_NAs = sum(is.na(Adj_Ct)),
                Adj_Ct = dplyr::if_else(Adj_NAs == Group_size, 0, Adj_Ct),
                Adj_Max_Diff = max(Adj_Ct, na.rm=TRUE) - min(Adj_Ct, na.rm=TRUE),
                Adj_SD = sd(Adj_Ct, na.rm = TRUE),
                Adj_Avg_Ct = mean(Adj_Ct, na.rm = TRUE),
                Diff_Avg = Avg_Ct - Adj_Avg_Ct,
                Diff_SD = SD - Adj_SD) %>%
  #dplyr::filter(Diff_SD > 0 & Diff_Avg > 0) %>%
  #dplyr::select(Ct, Adj_Ct, Outlier, SD, Adj_SD)
  dplyr::select(Sample_Name, Target_Name, Group_size,
                Ct, Max_Diff, SD, Avg_Ct, Distance, Outlier,   
                Adj_Ct, Adj_Max_Diff, Adj_SD, Adj_Avg_Ct) 

# Create a dataframe for reference gene
ref_gene_df <- qc_df %>%
  dplyr::filter(Target_Name %in% ref_gene) %>%
  dplyr::ungroup() %>%
  dplyr::select(Sample_Name, Adj_Avg_Ct) %>%
  dplyr::distinct_all() %>%
  dplyr::rename(Adj_Avg_Ref_Ct = Adj_Avg_Ct)
  
# Create a dataframe for non-reference genes
non_ref_gene_df <- qc_df %>%
  dplyr::filter(!(Target_Name %in% ref_gene)) %>%
  dplyr::ungroup() %>%
  dplyr::select(Sample_Name, Target_Name, Adj_Avg_Ct) %>%
  dplyr::distinct_all() %>%
  dplyr::left_join(ref_gene_df, by=c("Sample_Name"="Sample_Name")) %>%
  dplyr::mutate(Delta_Ct = Adj_Avg_Ct - Adj_Avg_Ref_Ct)

# Create a dataframe for reference samples
ref_sample_df <- non_ref_gene_df %>%
  dplyr::filter(grepl(pattern = ref_pattern, x = Sample_Name)) %>%
  dplyr::select(Sample_Name, Target_Name, Delta_Ct) %>%
  dplyr::rename(Delta_Ref_Ct = Delta_Ct) %>%
  dplyr::mutate(Sample_Name = gsub(pattern = ref_pattern, replacement = "", x = Sample_Name))

# Create a dataframe for non-reference samples
non_ref_sample_df <- non_ref_gene_df %>%
  dplyr::filter(!(grepl(pattern = ref_pattern, x = Sample_Name))) %>%
  dplyr::mutate(Sample_Name = gsub(pattern = non_ref_pattern, replacement = "", x = Sample_Name)) %>%
  dplyr::left_join(ref_sample_df, by=c("Target_Name"="Target_Name", "Sample_Name"="Sample_Name"))
  




for (i in 1:nrow(temp)) {
  #define index of last replicate for each gene in each sample
  #print(i)
  y <- i + temp$Group_Size[i] - 1
  #condition to make sure Sample_Name, Target_Name are same for first & last element of group being compared
  if (temp$Sample_Name[i] == temp$Sample_Name[y] &&
      temp$Target_Name[i] == temp$Target_Name[y] &&
      y < nrow(temp) &&
      !is.na(temp$Max_Diff[i]) && temp$Max_Diff[i] > 0.6) {
    #determine number of replicates that can be ignored to reduce Max_Diif < 0.6
    for (j in 1:ceiling(temp$Group_Size[i] * 0.5 - 1)) {
      #identify index of most distant replicate
      local_max <-
        match(max(temp$distance[i:y]), temp$distance[i:y]) + i - 1
      #print(local_max)
      #remove the Ct values at index of most distant replicate if Max_Diff <0.6 after removing it
      #set distance at index of most distant replicate to 0
      if (max(as.numeric(temp$Adj_Ct[i:y][!temp$Ct[i:y] %in% temp$Ct[local_max]])) -
          min(as.numeric(temp$Adj_Ct[i:y][!temp$Ct[i:y] %in% temp$Ct[local_max]])) <
          0.6) {
        temp$Adj_Ct[local_max] <- NA
        temp$distance[local_max] <- NA
        break
      } else{
        temp$Adj_Ct[local_max] <- NA
        temp$distance[local_max] <- NA
      }
    }
  }
}

temp1 <- temp %>%
  #sort by Sample_Name, then Target_Name, then Ct
  arrange(Sample_Name, Target_Name, Adj_Ct) %>%
  #group the data based on Sample_Name and Target_Name
  group_by(Sample_Name, Target_Name) %>%
  #compute Maximum Difference in Ct values, Standard Deviation, Group Size for each group
  dplyr::mutate(
    Max_Diff = max(as.numeric(Adj_Ct), na.rm = TRUE) - min(as.numeric(Adj_Ct), na.rm =
                                                                  TRUE),
    SD = sd(Adj_Ct, na.rm = TRUE),
    Group_Size = n()
  ) %>%
  #calculate average
  dplyr::mutate(
    Avg_Ct = mean(as.numeric(Adj_Ct), na.rm = TRUE),
    distance = abs(
      as.numeric(Adj_Ct, na.rm = TRUE) - as.numeric(Avg_Ct, na.rm = TRUE)
    )
  ) %>%
  dplyr::select(
    Sample_Name,
    Target_Name,
    Ct,
    Adj_Ct,
    Avg_Ct,
    Max_Diff,
    SD,
    Group_Size,
    distance
  )

#Adj_Ct <<- temp1
#}

#qpcr_analysis_2 <- function(temp1){
#identify all sample names and target names
#sample_names <- file %>% distinct(Sample_Name)
#target_names <- file %>% distinct(Target_Name) %>% arrange(Target_Name)

#compress the technical triplicates into single row with Avg Ct
temp2 <- temp1 %>%
  summarise(
    Sample_Name = first(Sample_Name),
    Target_Name = first(Target_Name),
    Avg_Ct = first(Avg_Ct)
  ) %>%
  group_by(Sample_Name)

#create a dataframe containing Average Ct of reference gene
temp3 <- temp2 %>%
  dplyr::filter(Target_Name == target_names[r, 1]) %>%
  rename(Avg_Reference_Ct = Avg_Ct)

#add a new column containing Average Reference Ct values by merging the above 2 dataframes
temp4 <-
  left_join(
    temp2,
    temp3,
    by = c("Sample_Name"),
    suffix = c("_original", "_new")
  )

#compute Delta Ct
temp4 <-
  temp4 %>% dplyr::mutate(DeltaCt = Avg_Ct - Avg_Reference_Ct)

#create a dataframe containing delta Ct of control sample
temp5 <- temp4 %>%
  dplyr::filter(Sample_Name == sample_names[c, 1]) %>%
  rename(Control_DeltaCt = DeltaCt)

#add a new column containing Control Delta Ct values by merging control dataframe with our main dataframe
temp6 <-
  right_join(
    temp4,
    temp5,
    by = c("Target_Name_original"),
    suffix = c("_original", "_new")
  ) %>% dplyr::select(1:3, 5, 6, 11) %>%
  dplyr::mutate(
    DeltaDeltaCt = .[[5]] - .[[6]],
    FC = 2 ^ -DeltaDeltaCt,
    Rel_Exp = 2 ^ -.[[5]]
  )

#rename the columns with appropriate column names
detailed_output <- temp6 %>%
  rename(
    Sample_Name = identity(1),
    Target_Name = identity(2),
    Avg_Ct = identity(3),
    Avg_ReferenceCt = identity(4),
    DeltaCt = identity(5),
    DeltaCt_Control = identity(6),
    DeltaDeltaCt = identity(7),
    FoldChange = identity(8),
    Rel_Exp = identity(9)
  )

detailed_results <<- detailed_output
}

concise_rel_exp <- function(detailed_output) {
  #identify all sample names and target names
  sample_names <- detailed_output %>% distinct(Sample_Name)
  target_names <-
    detailed_output %>% distinct(Target_Name) %>% arrange(Target_Name)
  concise_re_output <-
    data.frame(matrix(
      0 * (nrow(sample_names) * nrow(target_names)),
      nrow = nrow(sample_names),
      ncol = nrow(target_names)
    ))
  
  #create a new dataframe by reading Rel_Exp values from detailed output
  z = nrow(detailed_output)
  concise_re_output <- bind_cols(sample_names, concise_re_output)
  oldnames <- colnames(concise_re_output)[-1]
  concise_re_output <-
    concise_re_output %>% rename_at(vars(oldnames),  ~ target_names[[1]])
  for (i in 1:z) {
    concise_re_output[which(concise_re_output == detailed_output[[i, "Sample_Name"]]), detailed_output[[i, "Target_Name"]]] = concise_re_output[which(concise_re_output == detailed_output[[i, "Sample_Name"]]), detailed_output[[i, "Target_Name"]]] + detailed_output[[i, "Rel_Exp"]]
  }
  concise_re_output[concise_re_output == 0] <- NA
  
  concise_rel_exp_results <<- concise_re_output
}

concise_foldchange <- function(detailed_output) {
  #identify all sample names and target names
  sample_names <- detailed_output %>% distinct(Sample_Name)
  target_names <-
    detailed_output %>% distinct(Target_Name) %>% arrange(Target_Name)
  concise_fc_output <-
    data.frame(matrix(
      0 * (nrow(sample_names) * nrow(target_names)),
      nrow = nrow(sample_names),
      ncol = nrow(target_names)
    ))
  
  #create a new dataframe by reading FolfChange values from detailed output
  z = nrow(detailed_output)
  concise_fc_output <- bind_cols(sample_names, concise_fc_output)
  oldnames <- colnames(concise_fc_output)[-1]
  concise_fc_output <-
    concise_fc_output %>% rename_at(vars(oldnames),  ~ target_names[[1]])
  for (i in 1:z) {
    concise_fc_output[which(concise_fc_output == detailed_output[[i, "Sample_Name"]]), detailed_output[[i, "Target_Name"]]] = concise_fc_output[which(concise_fc_output == detailed_output[[i, "Sample_Name"]]), detailed_output[[i, "Target_Name"]]] + detailed_output[[i, "FoldChange"]]
  }
  concise_fc_output[concise_fc_output == 0] <- NA
  
  concise_fold_change_results <<- concise_fc_output
}

#enclose all function within brackets. Otherwise, the userinput prompts from qpcr_analysis() will be skipped
{
  read_files("C:/Users/KailasammS/Desktop/qpcr")
  #read_files("/cloud/project/qpcr")
  
  qpcr_analysis(input_data)
  #qpcr_analysis_2(Adj_Ct)
  concise_rel_exp(detailed_results)
  concise_foldchange(detailed_results)
  
  #export the output files (change KailasammS to your username)
  #write.csv(Adj_Ct, "C:/Users/KailasammS/Desktop/Adj_Ct.csv")
  #write.csv(Adj_Ct, "/cloud/project/Adj_Ct.csv")
  write.csv(detailed_results,
            "C:/Users/KailasammS/Desktop/detailed_results.csv")
  #write.csv(detailed_results, "/cloud/project/detailed_results.csv")
  write.csv(concise_rel_exp_results,
            "C:/Users/KailasammS/Desktop/rel_exp_results.csv")
  #write.csv(concise_rel_exp_results, "/cloud/project/rel_exp_results.csv")
  write.csv(concise_fold_change_results,
            "C:/Users/KailasammS/Desktop/FC_results.csv")
  #write.csv(concise_fold_change_results, "/cloud/project/FC_results.csv")
}