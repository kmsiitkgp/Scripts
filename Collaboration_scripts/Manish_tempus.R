# Configuration
file       <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/EBR_updatedDESeq2_results.xlsx"
gmt_dir    <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/Signatures/Human"
output_dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/"

dir.create(
  file.path(output_dir, "omnipathr-log"),
  recursive = TRUE,
  showWarnings = FALSE
)

# Define the sheet names
contrasts <- c("death", "progression", "recurrence", "metastasis", 
               "metastases_vs_primary", "primary_vs_liver")

# 1. Pathway Processing Loop
for (contrast_name in contrasts) {
  # Load data
  data <- load_smart(input_path = file, sheet = contrast_name)
  colnames(data)[1] <- "SYMBOL"
  
  # Analyze
  analyze_pathways(contrast = contrast_name, 
                   input_data = data, 
                   only_deg = FALSE, 
                   gmt_dir, 
                   output_dir)
}

# 2. Pathway Plotting Loop
for (contrast_name in contrasts) {
  # Construct path dynamically
  pathway_xlsx <- file.path(output_dir, contrast_name, "Pathways.xlsx")
  
  # Plot
  plot_pathways(contrast = contrast_name, 
                pathway_xlsx = pathway_xlsx, 
                output_dir = output_dir)
}

# 3. TF Processing Loop
for (contrast_name in contrasts) {
  # Load data
  data <- load_smart(input_path = file, sheet = contrast_name)
  colnames(data)[1] <- "SYMBOL"
  
  # Analyze
  analyze_tfs(contrast = contrast_name,
              species = "Human",
              output_dir, 
              degs = data)
}

# 4. TF Plotting Loop
for (contrast_name in contrasts) {
  # Construct path dynamically
  tf_xlsx <- file.path(output_dir, contrast_name, "TF_and_Pathway_Activity.xlsx")
  
  # Plot
  plot_tfs(contrast = contrast_name,
         tf_xlsx = tf_xlsx, 
         output_dir)
}
                     
