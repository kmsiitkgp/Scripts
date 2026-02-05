#!/bin/bash -l

#$ -N pyega3          # Set Job Name
#$ -l mem_free=96G    # Request memory
#$ -j y               # Merge standard output and standard error
#$ -cwd               # Set current working directory
#$ -t 1-1             # Submit an array job with 1 task

# You need to add path for Conda and use source activate instead of conda activate
PATH=$HOME/miniconda3/bin/:$PATH
source activate MAGECK
pyega3 --config-file ~/projects/EGA/EGA_Credentials.json datasets

# # Copy folder names, file names, checksums displayed from files command to excel file
# pyega3 -cf ~/projects/EGA/EGA_Credentials.json files EGAD00001007575
# pyega3 -cf ~/projects/EGA/EGA_Credentials.json files EGAD00001003977
# pyega3 -cf ~/projects/EGA/EGA_Credentials.json files EGAD00001004218

# cd ~/common/EGA
# pyega3 -cf ~/projects/EGA/EGA_Credentials.json fetch EGAD00001011108
# pyega3 -cf ~/projects/EGA/EGA_Credentials.json fetch EGAD00001007574
# pyega3 -cf ~/projects/EGA/EGA_Credentials.json fetch EGAD00001007576
# pyega3 -cf ~/projects/EGA/EGA_Credentials.json fetch EGAD00001007653
# pyega3 -cf ~/projects/EGA/EGA_Credentials.json fetch EGAD00001007654

cd ~/common/EGAD00001004218
# pyega3 --connections 25 -cf ~/projects/EGA/EGA_Credentials.json fetch EGAD00001004218
# cd ~/common/EGA/EGAD00001003977
# pyega3 --connections 25 -cf ~/projects/EGA/EGA_Credentials.json fetch EGAD00001003977
# cd ~/common/EGA/EGAD00001007575
# pyega3 --connections 25 -cf ~/projects/EGA/EGA_Credentials.json fetch EGAD00001007575

# EGAD00001004218 completed download but many folders missing. So, redownload missing folders
files=(EGAF00002053530 EGAF00002053568 EGAF00002053689 EGAF00002053699 EGAF00002053861 EGAF00002054029 \
EGAF00002054030 EGAF00002054033 EGAF00002054034 EGAF00002054035 EGAF00002054036 EGAF00002054037 \
EGAF00002054041 EGAF00002054042 EGAF00002054093 EGAF00002054149 EGAF00002054386 EGAF00002053545 \
EGAF00002053546 EGAF00002053551 EGAF00002053554 EGAF00002053562 EGAF00002053563 EGAF00002053564 \
EGAF00002053565 EGAF00002053566 EGAF00002053599 EGAF00002053621 EGAF00002053659 EGAF00002053712 \
EGAF00002053723 EGAF00002053729 EGAF00002053767 EGAF00002053829 EGAF00002054002 EGAF00002054003 \
EGAF00002054004 EGAF00002054006 EGAF00002054009 EGAF00002054013 EGAF00002054049 EGAF00002054051 \
EGAF00002054175)

for i in ${files[*]}
do
echo $i
pyega3 --connections 25 -cf ~/projects/EGA/EGA_Credentials.json fetch $i
done