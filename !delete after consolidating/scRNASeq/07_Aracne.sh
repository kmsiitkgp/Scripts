#!/usr/bin/env bash

#!/bin/bash -l

# Set Project Name
# Set Job Name
#$ -N Aracne-AP
# Request memory. You can find maxvmem from qstat -j <jobid> & adjust mem_free accordingly in final run.
#$ -l mem_free=8G
# Merge standard output and standard error
#$ -j y
# Set current working directory
#$ -cwd
# Submit an array job with n tasks. Set n to number of SRR ids in the SRR_Acc_List.txt
#$ -t 1

# Download ARACHNE-AP from https://sourceforge.net/projects/aracne-ap/
# wget https://sourceforge.net/projects/aracne-ap/files/latest/download --no-check-certificate -O $HOME/NGSTools/Aracne.jar

# You need to add path for java as ARACNE-AP needs Java to run
PATH=$HOME/NGSTools/jre1.8.0_271/bin:$PATH

# Download list of mouse transcription factors from  Uniprot (https://www.uniprot.org/uniprot/?query=transcription+factors&sort=score)
# Use Gprofiler to convert Uniprot ID to mouse gene names (https://biit.cs.ut.ee/gprofiler/convert)
# Remove duplicates and paste transcription factors to txt file. This is the "Mouse_TFs.txt" file.

# Go to the sheet containing normalized counts from DESeq2 in excel file
# Format the cells to 2 decimal places. Make sure column names are formatted as text and not numbers
# Save the sheet as txt file and use it as expression matrix.

# increase heap size by 5GB using -Xmx5G
#MATRIX=$HOME/projects/Hany/scripts/Expression_Matrix_B_Cells.txt
MATRIX=$HOME/projects/Hany/scripts/Expression_Matrix_Epithelial_Cells.txt
TRANSCRIPTION_FACTOR=$HOME/projects/Hany/scripts/Mouse_TFs.txt
OUTPUT=$HOME/scratch/Hany/ARACNE_results
P_VAL=1E-8
SEED=1

# Calculate the MI threshold for the dataset
java -Xmx5G -jar $HOME/NGSTools/Aracne.jar -e $MATRIX  -o $OUTPUT --tfs $TRANSCRIPTION_FACTOR --pvalue $P_VAL --seed $SEED --calculateThreshold

# Run ARACNe on bootstraps of the input matrix
for i in {1..5}
do
java -Xmx5G -jar $HOME/NGSTools/Aracne.jar -e $MATRIX  -o $OUTPUT --tfs $TRANSCRIPTION_FACTOR --pvalue $P_VAL  --seed $i
done

# Consolidate, i.e. combine the bootstraps into a final network file
java -Xmx5G -jar $HOME/NGSTools/Aracne.jar -o $OUTPUT --consolidate