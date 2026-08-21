#!/bin/bash

conda activate omics
module load biomodal/2.0.0
module load apptainer/1.0.1

export BIOMODAL_INSTANCE_DIRECTORY=/home/kailasamms/resources/biomodal

# Step 1: Initialize biomodal (only needed once ever)
biomodal init

# Step 2: Test setup
biomodal run duet --test

# Step 3: Install modality XPLR
pip install --extra-index-url https://europe-python.pkg.dev/prj-biomodal-modality/modality-pypi/simple modality