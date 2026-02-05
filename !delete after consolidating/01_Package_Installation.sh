#!/usr/bin/env bash

### Use the codes below to install latest versions of softwares from login i.e. submit node, not from computing node. 
### Some commands like unzip dont seem to work in computing node (? check if true)
### After installing softwares, add their path to .bashrc file as described below.
### When you run scripts on computing node, you need to indicate PATH within the script as computing nodes cant access your .bashrc file.

OUTPUT_DIR=/hpc/home/kailasamms/NGSTools
GENOME_DIR=$HOME/NGSTools/Reference_Genomes

# Option 1 (RECOMMENDED): Add PATH to bashrc to automatically invoke every session
# vim ~/.bashrc --> i --> PATH=$HOME/NGSTools/cellranger-8.0.0/bin/:$PATH --> Escape key  --> :wq! --> Enter key --> source .bashrc (to reload bashrc)
# Option 2 (NOT RECOMMENDED): Add directory to your $PATH. This will allow you to invoke the command in this session ONLY.
#export PATH=$OUTPUT_DIR/cellranger-9.0.0/bin/:$PATH

#*****************************************************************#
#***********************DOWNLOAD REFERENCES***********************#
#*****************************************************************#

# Get latest version of "UNMASKED (not sm or rm versions), PRIMARY ASSEMBLY (not top level version)" Human Reference Genome (Fasta file)
wget -P $GENOME_DIR/Human/ https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz --no-check-certificate
gunzip $GENOME_DIR/Human/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
                           
# Download and unpack the latest version of Human Reference Genome (GTF file)
wget -P $GENOME_DIR/Human/ https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz --no-check-certificate
gunzip $GENOME_DIR/Human/Homo_sapiens.GRCh38.113.gtf.gz
                           
# Get latest version of "UNMASKED (not sm or rm versions), PRIMARY ASSEMBLY (not top level version)" Mouse Reference Genome (Fasta file)
wget -P $GENOME_DIR/Mouse/ https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz --no-check-certificate
gunzip $GENOME_DIR/Mouse/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz                         

# Download and unpack the latest version of Mouse Reference Genome (GTF file)
wget -P $GENOME_DIR/Mouse/ https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm39.113.gtf.gz --no-check-certificate
gunzip $GENOME_DIR/Mouse/Mus_musculus.GRCm39.113.gtf.gz
                           
#*****************************************************************#
#***************************CELL RANGER***************************#
#*****************************************************************#

# Download and unpack latest version of CellRanger to NGSTools folder
wget -O $OUTPUT_DIR/cellranger-9.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz?Expires=1744272840&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=OVOzG8KCCrcoeQDgOPvbY3YXsPKjLCFK464hiMvkAF4Eh7obLTLNiFH-Py-3eA7l~VptmRnkoAFXdPJ2Svv1JNIxxw2QYtQxVL0~L1W9mWNLJ-n6OwfcPbnYS48yT1yjYZ4b3DzDYIid0yLIjuEga87nOA~yzWGy4RKt7UZYjALq2X~aNWBXoMVXrD-TfrGsxkjHjX5ARgbL0YlWjQwTZoyZkbLa~Uo42EgYn~GRsepu1eJxi8gFEnGnup5hBss06fqEsFGKpljUrortk3IyzmXs3gV0g8wFzmBKEE2LxMocbjv8~sFY8ekupfinJ~P1vB~1ZEsXDipw0r-Dgz7-tw__"
tar -xzvf $OUTPUT_DIR/cellranger-9.0.1.tar.gz -C $OUTPUT_DIR

# Download and unpack human and mouse reference genome to NGSTools folder
wget -P $OUTPUT_DIR "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
tar -xzvf $OUTPUT_DIR/refdata-gex-GRCh38-2024-A.tar.gz -C $OUTPUT_DIR

wget -P $OUTPUT_DIR "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz"
tar -xzvf $OUTPUT_DIR/refdata-gex-GRCm39-2024-A.tar.gz -C $OUTPUT_DIR

# Add to path
PATH=$HOME/NGSTools/cellranger-9.0.1/bin/:$PATH

# Verify if CellRanger is installed properly
cellranger --help
cellranger --version
#cellranger testrun --id=tiny

#*****************************************************************#
#***************************SPACE RANGER**************************#
#*****************************************************************#

# Download and unpack latest version of SpaceRanger to NGSTools folder
wget -O $OUTPUT_DIR/spaceranger-3.1.3.tar.gz "https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-3.1.3.tar.gz?Expires=1744282819&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=HxrdcyU6iz~MEteUq6LPmQUdINIVHoeZfTj1bBtGI5UkOayCH68bAtnq-6qR39pfX9Zqn~jw8E-oGb9B-v73o7Yq6p1ECeJ1s3836Be3-TclYyY076Jy0F3ALpakfGt5GZEsdQtqhKuT2mdZ35aiXL9MJqgeUEKg3rBrtyfZmleSyjhB89CbdPYimEK1IpFq-deXBMB-9bcMEKmA1tyPOYvJZSII~II35AIaPI4N~yFN3uLqjyA~oWOF4G8n9UHQS4oFRgOIwmDciZO0FT2rhZaoGFiou7byNw7qxoYTVZVQf8m-wLnsZVPapM2JDhF80wx8Gif9X-ZYN3H5AwmrZg__"
tar -xzvf $OUTPUT_DIR/spaceranger-3.1.3.tar.gz -C $OUTPUT_DIR

# Download and unpack human and mouse reference genome to NGSTools folder
wget -P $OUTPUT_DIR "https://cf.10xgenomics.com/supp/spatial-exp/refdata-gex-GRCh38-2020-A.tar.gz"
tar -xzvf $OUTPUT_DIR/refdata-gex-GRCh38-2020-A.tar.gz -C $OUTPUT_DIR

wget -P $OUTPUT_DIR "https://cf.10xgenomics.com/supp/spatial-exp/refdata-gex-mm10-2020-A.tar.gz"
tar -xzvf $OUTPUT_DIR/refdata-gex-mm10-2020-A.tar.gz -C $OUTPUT_DIR

# Add to path
PATH=$HOME/NGSTools/spaceranger-3.1.3/bin/:$PATH

# Verify if Spaceranger is installed properly
spaceranger --help
spaceranger --version

#*****************************************************************#
#******************************GATK*******************************#
#*****************************************************************#

# Download, unpack and compile latest version of GATK to NGSTools folder
# https://github.com/broadinstitute/gatk

wget -P $OUTPUT_DIR https://github.com/broadinstitute/gatk/releases/download/4.6.1.0/gatk-4.6.1.0.zip
unzip $OUTPUT_DIR/gatk-4.6.1.0.zip -d $OUTPUT_DIR

# Add to path
PATH=$HOME/NGSTools/gatk-4.6.1.0/:$PATH

# Verify if GATK is installed properly
gatk --help
gatk --version

#*****************************************************************#
#******************************FASTQC*****************************#
#*****************************************************************#

# NOTE: FastQC needs a JRE (Java Runtime Environment). So, install Java first and then FastQC

# Download and unpack Java to NGSTools folder
wget -P $OUTPUT_DIR https://download.oracle.com/java/24/latest/jdk-24_linux-x64_bin.tar.gz
tar -xzvf $OUTPUT_DIR/jdk-24_linux-x64_bin.tar.gz -C $OUTPUT_DIR

# Add to path
PATH=$HOME/NGSTools/jdk-24/bin/:$PATH

# Verify if Java is installed properly
java --help
java --version

wget -P $OUTPUT_DIR https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip --no-check-certificate
unzip $OUTPUT_DIR/fastqc_v0.12.1.zip -d $OUTPUT_DIR
chmod u+x $OUTPUT_DIR/FastQC/fastqc

# Add to path
PATH=$HOME/NGSTools/FastQC/:$PATH

# Verify if fastqc is installed properly
fastqc --help
fastqc --version

# #*****************************************************************#
# #***************************GDC CLIENT****************************#
# #*****************************************************************#

# # Install from submit i.e.login node. Installation fails from computing node
# module load gcc/11.1.0

# # Download and unpack latest version of GDC Data Transfer Tool Client to NGSTools folder.
# # DO NOT download GDC Data Transfer Tool UI. We will download Ubuntu x64.
# wget -P $OUTPUT_DIR https://gdc.cancer.gov/system/files/public/file/gdc-client_2.3_Ubuntu_x64-py3.8-ubuntu-20.04.zip
# unzip $OUTPUT_DIR/gdc-client_2.3_Ubuntu_x64-py3.8-ubuntu-20.04.zip -d $OUTPUT_DIR
# unzip $OUTPUT_DIR/gdc-client_2.3_Ubuntu_x64.zip -d $OUTPUT_DIR/gdc-client_v2.3_Ubuntu_x64

# # Add to path
# PATH=$HOME/NGSTools/gdc-client_v2.3_Ubuntu_x64/:$PATH

# # Verify if gdc-client is installed properly
# gdc-client --help
# gdc-client download --help

# #*****************************************************************#
# #***************************SRA TOOLKIT***************************#
# #*****************************************************************#

conda activate R
conda install --channel bioconda sra-tools

# Direct installation requires GLIBC2.2.7 for it to work. So, install using conda


#*****************************************************************#
#*************************BOOST LIBRARIES*************************#
#*****************************************************************#

wget -P /hpc/home/kailasamms/NGSTools/ https://archives.boost.io/release/1.88.0/source/boost_1_88_0_rc1.tar.gz --no-check-certificate
tar -xzvf $OUTPUT_DIR/boost_1_88_0_rc1.tar.gz -C $OUTPUT_DIR

#*****************************************************************#
#****************************SAMBAMBA*****************************#
#*****************************************************************#

# Download and unpack latest version of Sambamba to NGSTools folder
wget -P $OUTPUT_DIR https://github.com/biod/sambamba/releases/download/v1.0.1/sambamba-1.0.1-linux-amd64-static.gz
mkdir $OUTPUT_DIR/sambamba-1.0.1-linux-amd64-static
gunzip -c $OUTPUT_DIR/sambamba-1.0.1-linux-amd64-static.gz > $OUTPUT_DIR/sambamba-1.0.1-linux-amd64-static/sambamba
chmod u+x $OUTPUT_DIR/sambamba-1.0.1-linux-amd64-static/sambamba

# Add to path
PATH=$HOME/NGSTools/sambamba-1.0.1-linux-amd64-static/:$PATH

# Verify if sambamba is installed properly
sambamba --help
sambamba --version

#*****************************************************************#
#****************************SAMTOOLS*****************************#
#*****************************************************************#

# Install from submit i.e.login node. Installation fails from computing node
module load gcc/11.1.0

# Download and unpack latest version of Samtools to NGSTools folder
wget -P $OUTPUT_DIR https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2
tar -xvf $OUTPUT_DIR/samtools-1.21.tar.bz2 -C $OUTPUT_DIR
cd $OUTPUT_DIR/samtools-1.21
./configure --prefix=$OUTPUT_DIR/samtools-1.21
make
make install

# Add to path
PATH=$HOME/NGSTools/samtools-1.21/:$PATH

# Verify if samtools is installed properly
samtools --help
samtools --version

#*****************************************************************#
#******************************STAR*******************************#
#*****************************************************************#

# Install from submit i.e.login node. Installation fails from computing node
module load gcc/11.1.0

# Download, unpack and compile latest version of STAR to NGSTools folder
wget -P $OUTPUT_DIR https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz
tar -xvzf $OUTPUT_DIR/2.7.11b.tar.gz -C $OUTPUT_DIR
cd $OUTPUT_DIR/STAR-2.7.11b/source
make STAR

# Add to path
PATH=$HOME/NGSTools/STAR-2.7.11b/bin/Linux_x86_64_static/:$PATH

# Verify if STAR is installed properly
STAR --help
STAR --version

#*****************************************************************#
#******************************CONDA******************************#
#*****************************************************************#

# https://conda.io/projects/conda/en/stable/user-guide/getting-started.html#managing-conda

mkdir $HOME/miniconda3
wget -P $HOME/miniconda3 https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash $HOME/miniconda3/Miniconda3-latest-Linux-x86_64.sh -b -u -p $HOME/miniconda3
# -b to be able to run unattended, which means that all of the agreements are automatically accepted without user prompt
# -u updates any existing Miniconda version in the installation directory if there is one
# -p is the directory into which Miniconda will be installed
# PATH is automatically added for miniconda and all python packages.
# Hence, you do not need to edit .bashrc file 

# Update conda
conda update -n base conda

# Add to path
PATH=$HOME/miniconda3/bin/:$PATH

# Disable automatically activating conda's base environment each time you login
# RECOMMENDED: You should create separate environment and not use base environment
conda config --set auto_activate_base false
source .bashrc

# Verify if conda is installed properly
conda --help
conda --version

#*****************************************************************#
#****************************BIOMODAL*****************************#
#*****************************************************************#

conda create --name biomodal --channel conda-forge python=3.13
conda install --channel conda-forge google-cloud-sdk
conda install --channel conda-forge jq
conda install --channel conda-forge singularity

bash <(curl -s https://app.biomodal.com/cli/download)
# login using username:saravanakumar.kailasammani@cshs.org, pass:usual!usual@

OUTPUT_DIR=/hpc/home/kailasamms/NGSTools
mv biomodal.1.1.3.zip $OUTPUT_DIR
unzip $OUTPUT_DIR/biomodal.1.1.3.zip -d $OUTPUT_DIR

sh $OUTPUT_DIR/biomodal.1.1.3/biomodal-hpc-utils-conda
## Choose relevant geographical location of the biomodal docker registry (1-3)?
# 3) us 
## Current queueSize is: 200. This is the maximum number of concurrent jobs
## that the pipeline will attempt to launch. Would you like to change it? (y/n)
# n 
## Which executor would you like to use? (choose 1-4)
# 3) sge

## Add these to bashrc file after creating a directory "tmp" in $HOME/scratch/
## Make sure to grant permissions, else you will get errors from haplotype caller
chmod og+rwx ~/scratch/tmp
## SINGULARITY_CACHEDIR MUST be same as "libraryDir" defined in $HOME/biomodal/nextflow.config
export TMPDIR=$HOME/scratch/tmp
export SINGULARITY_TMPDIR=$TMPDIR
export SINGULARITY_CACHEDIR=$HOME/biomodal/data_bucket/singularity-images
export NXF_SINGULARITY_CACHEDIR=$SINGULARITY_CACHEDIR
PATH=$HOME/biomodal/:$PATH

## Next, copy these lines to $HOME/biomodal/nextflow.config
singularity {  
	enabled    = true  
	autoMounts = true  
	//libraryDir = "$HOME/biomodal/data_bucket/singularity-images"
	envWhitelist = "TMPDIR,SINGULARITY_TMPDIR,SINGULARITY_CACHEDIR,NXF_SINGULARITY_CACHEDIR"
    runOptions = '--bind "$TMPDIR:/tmp"'
}
params {
  registry = "us-docker.pkg.dev/cegx-releases/us-prod"
}
executor {
  queueSize = "200"
}
report {
    overwrite = true
}
process {  
	process.executor 	= "sge"  
	process.penv 		= "smp"  
	process.queue 		= "all.q"
	clusterOptions 		= "-S /bin/bash" 
//	You can include the below options if needed
	/*process.cpus 		= 2
	clusterOptions 		= "-l h_rt=48:0:0"  
	clusterOptions 		= "-l h_vmem=120G"*/	
}

## Download biomodal pipelines
conda activate biomodal
module load singularity/3.6.0    # although conda has singularity, you need to use admin installed singularity
biomodal init

## Would you like to share events and the pipeline metrics report at the end of all successful analysis runs with biomodal? (1-2)
# 1) Yes
## The duet software location 'init_folder' property is not currently defined. Would you like to set it to /home/kailasamms? (1-2)
# 2) No
## Please enter an alternative existing full path directory:
# /hpc/home/kailasamms/biomodal/
## Select the pipeline error strategy you would like to use for all duet pipeline runs:
- Normal. Retries failed jobs up to 10 times depending on the exit status.
- FailFast. Workflow exits immediately upon an error and does not automatically retry.
# 1) FailFast

## If you get error, it will ask to run line below
$ /home/kailasamms/miniconda3/envs/biomodal/bin/python -m pip install google-crc32c --upgrade --target /hpc/home/kailasamms/miniconda3/envs/biomodal/share/google-cloud-sdk-529.0.0-0/lib/third_party

## Once file downloads are complete, check if installation is working properly using 01_biomodal_test.sh
## Also define the variables again in qsub scripts you run
TMPDIR=$HOME/scratch/tmp
SINGULARITY_TMPDIR=$TMPDIR
SINGULARITY_CACHEDIR=$HOME/biomodal/data_bucket/singularity-images
NXF_SINGULARITY_CACHEDIR=$SINGULARITY_CACHEDIR

biomodal test

conda activate biomodal

#*****************************************************************#
#****************************R (base)*****************************#
#*****************************************************************#

# Installing R base directly is easy but several R packages will need CMake, fftw3, hdf5, etc
# and other dependencies. Fixing each of these is a headache and doesnt work most of the
# time. Several packages like Seurat couldn't be installed if R base is installed directly.

# RECOMMENDED option is to install R through miniconda.
# IMPORTANT: Install python in the environment during creation 
# DO NOT install r-essentials as it has ~100 non-R & ~175 R packages as compared to r-base
# Use conda forge channel rather than default Anaconda channel
# https://stackoverflow.com/questions/39857289/should-conda-or-conda-forge-be-used-for-python-environments

conda create --name R --channel conda-forge python=3.13
conda activate R
module load hdf5/1.8.18                            # this is important. Else, Rcpp wont be installed properly
conda install --channel conda-forge r-base=4.4.3   # install from channel conda forge
# Then type R to enter R and install packages as you normally do in R

conda install bioconda::bioconductor-rhtslib  	   # needed for ensembldb
conda install conda-forge::liblzma-devel           # needed for ensembldb
conda install --channel conda-forge r-textshaping  # needed for CellChat
conda install --channel conda-forge imagemagick    # needed for GSVA
conda install --channel anaconda gmp               # needed for leidenAlg, Banksy
conda install --channel conda-forge xorg-xorgproto # needed for scCustomize
conda install --channel conda-forge scvi-tools     # needed for integration using scVI
conda install --channel conda-forge scapy          # needed for integration using scVI
conda install --channel bioconda gxf2bed           # needed for RSeQC (unable to install in MAGECK env)

# mageck and multiqc needs python 2.7
conda create --name MAGECK
conda activate MAGECK
conda install --channel bioconda mageck
conda install --channel bioconda multiqc
conda install --channel bioconda rseqc

### EXTRA NOTES ON CONDA
# To manually, activate conda's base environment after logging in 
conda activate

# Create and activate new conda environment. Use deactivate to exit
conda create --name R
conda activate R
conda deactivate

# To see all environements you created
conda info --envs

# To see all installed packages in base or user defined environment
conda list

# Check if an environment "R" has package "beautifulsoup4".
# If not install it and see all installed packages in "hany_proj"
conda activate R
conda search beautifulsoup4
conda install beautifulsoup4
conda install pandas=1.4.4  # to install a specific version

# Remove package (in conda list if channel = pypi for a given package, use pip)
conda remove pandas
pip uninstall pandas

#**********************CITE-seq-Count (through conda)**********************#
# Install and upgrade to latest version of CITESeq-count through miniconda
conda activate NGS   #NGS environment was created for all NGS work
conda install pip  #This will install pip to NGS environment. 
# If you do not do this, pip from base environment will be used and pacakge
# will be installed in base environment's bin folder
pip install CITE-seq-Count==1.4.5
pip install CITE-seq-Count --upgrade

# Verify if CITE-seq-Count is installed properly
CITE-seq-Count --help

# PATH is automatically added for miniconda and all python packages.
# Hence, you do not need to edit .bashrc file 

#**********************MACS3 (through conda)**********************#
# Install and upgrade to latest version of MACS3 through miniconda
conda activate NGS   #NGS environment was created for all NGS work
conda install pip  #This will install pip to NGS environment. 
# If you do not do this, pip from base environment will be used and package
# will be installed in base environment's bin folder
pip install macs3

# Verify if macs3 is installed properly 
macs3 --help

# PATH is automatically added for miniconda and all python packages.
# Hence, you do not need to edit .bashrc file 

#**********************leidenalg (through conda)**********************#
# Install leidenalg needed for FindClusters() of Seurat through miniconda
conda activate NGS   #NGS environment was created for all NGS work
conda install pip    #This will install pip to NGS environment. 
# If you do not do this, pip from base environment will be used and package
# will be installed in base environment's bin folder
pip install leidenalg

# PATH is automatically added for miniconda and all python packages.
# Hence, you do not need to edit .bashrc file 


#*****************************************************************#
#*******************************BWA*******************************#
#*****************************************************************#
# Read about BWA from https://github.com/lh3/bwa
# It can be installed via Conda too

conda activate MAGECK
conda install --channel conda-forge bwa
bwa mem

#*****************************************************************#
#****************************bcftools*****************************#
#*****************************************************************#

# Unable to install bcftools. Used module load bcftools

#*****************************************************************#
#*******************************rsem******************************#
#*****************************************************************#

conda activate MAGECK
conda install --channel bioconda rsem

# Install MAGECK, VISPR, bowtie2, bwa-mem2, samtools through conda.
# To generate pdf reports in MAGeCK, you need pdflatex
conda search --channel bioconda --full-name sra-tools
conda install --channel conda-forge --channel bioconda sra-tools=3.2.1
conda install --channel conda-forge pdflatex
conda install --channel bioconda vispr
conda install --channel conda-forge jinja2=3.0.3
conda install --channel bioconda bowtie2
conda install --channel bioconda samtools
conda install --channel bioconda bamtools
conda install --channel bioconda bedtools
conda install --channel bioconda cutadapt
conda install --channel bioconda htseq
conda install --channel bioconda bwa-mem2
conda install --channel conda-forge r-cairo
conda install --channel conda-forge r-xml
conda install --channel conda-forge r-igraph
conda install --channel conda-forge r-leiden
conda install --channel conda-forge r-survminer
conda install --channel conda-forge r-stringi
conda search --channel conda-forge --full-name python
conda search --channel conda-forge --full-name r-igraph
conda search --channel conda-forge --full-name r-leiden
conda install --channel conda-forge r-ragg
conda update mageck
conda update --channel bioconda --all
# NOTE: While installing MAGeCK-VISPR, also install old version of jinja2
# conda install --channel conda-forge jinja2=3.0.3
# The latest version of jinja2 will give following error
# ImportError: cannot import name 'Markup' from 'jinja2'

# Check MAGECK is working properly
mageck --version
vispr --version
bowtie2 --version
samtools --version
multiqc --version
cutadapt --version
bedtools --version
htseq-count --version
java --version
mageck --help
bowtie2 --help
samtools --help
rsem-prepare-reference--help

# Download picard.jar from link provided in below site  
# https://broadinstitute.github.io/picard/
java -jar ~/NGSTools/picard.jar -h
conda install pyega3

########################NOT NECESSARY AS WE USE MINICONDA
# mkdir $HOME/Python3
# wget -O $HOME/Python3/Python-3.11.1.tgz https://www.python.org/ftp/python/3.11.1/Python-3.11.1.tgz
# tar -xzvf $HOME/Python3/Python-3.11.1.tgz -C $HOME/Python3/
# cd $HOME/Python3/Python-3.11.1/
# $HOME/Python3/Python-3.11.1/configure --prefix=$HOME/Python3/    
# make
# make install
# ## Copy and paste PATH to .bashrc file using vim. (vim ~/.bashrc --> i  --> Escape key --> :wq --> Enter key --> source .bashrc (to reload bashrc)
# #PATH=$HOME/Python3/Python-3.9.1/:$PATH
# ## Copy and paste PYTHONPATH to .bashrc file using vim.
# #PYTHONPATH=$HOME/Python3/Python-3.9.1/
# # Verify if Python installed properly using which and version
# which python
# python --version

# ##*******Install latest version of pip*******##
# wget https://bootstrap.pypa.io/get-pip.py -P $HOME/Python3/
# python $HOME/Python3/get-pip.py --user
# #Upgrade pip, setuptools and wheel
# python -m pip install --upgrade pip setuptools wheel
# ## Copy and paste PATH to .bashrc file using vim.	PATH=$HOME/.local/bin/:$PATH
# # Verify if Python installed properly using which and version
# which pip
# pip --version

# ##*******SRAToolKit **************
# # Install from submit i.e.login node. Installation fails from computing node
# module load gcc/11.1.0

# # Download and unpack latest version of SRA Toolkit to NGSTools folder 
# # https://www.ncbi.nlm.nih.gov/sra 
# wget -P $OUTPUT_DIR https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.2.1/sratoolkit.3.2.1-ubuntu64.tar.gz
# tar -xzvf $OUTPUT_DIR/sratoolkit.3.2.1-ubuntu64.tar.gz  -C $OUTPUT_DIR

# # Add to path
# PATH=$HOME/NGSTools/sratoolkit.3.2.1-ubuntu64/bin/:$PATH

# # Verify if fastq-dump is installed properly
# fasterq-dump --help
# which fasterq-dump

# # If you get error asking you to run "vdb-config --interactive", do it.
# # Simply exit by pressing "X" key and now try if fastq-dump is working
# vdb-config --interactive
# fastq-dump --help