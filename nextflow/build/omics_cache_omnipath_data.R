# ==============================================================================
# PRE-BUILD NETWORK CACHING SCRIPT
#
# WHY: We are pre-downloading these network files (CollecTRI and PROGENy)
# because the container build process occurs in a network-restricted
# environment (or we want to ensure reproducibility).
#
# INSTRUCTIONS:
# 1. Run this script in the same directory as your 'omics.def' file.
# 2. This will generate four .rds files in this folder.
# 3. After this script completes successfully, run the Apptainer build command:
#    sudo apptainer build --force omics_r4.4.3.sif omics_r4.4.3.def 2>&1 | tee omics.log
#
# The 'omics.def' file is configured to copy these files into the container
# using the %files directive, bypassing the need for internet access during build.
# ==============================================================================

# Load necessary libraries
if (!requireNamespace("OmnipathR", quietly = TRUE)) install.packages("OmnipathR")
if (!requireNamespace("decoupleR", quietly = TRUE)) install.packages("decoupleR")

library(decoupleR)
library(OmnipathR)

message("Downloading networks for offline container build...")

# Save with exact names expected by the Apptainer %files section
saveRDS(get_collectri(organism='human', split_complexes=FALSE), "collectri_human.rds")
saveRDS(get_collectri(organism='mouse', split_complexes=FALSE), "collectri_mouse.rds")
saveRDS(get_progeny(organism='human', top=500), "progeny_human.rds")
saveRDS(get_progeny(organism='mouse', top=500), "progeny_mouse.rds")

message("----------------------------------------------------------------------")
message("Download complete.")
message("Files generated: collectri_human.rds, collectri_mouse.rds, progeny_human.rds, progeny_mouse.rds")
message("You may now proceed to run: apptainer build --force omics.sif omics.def")
message("----------------------------------------------------------------------")