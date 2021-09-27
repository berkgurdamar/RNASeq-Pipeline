##################################################
## Project: RNASeq-Pipeline
## Purpose: Setup file
## Date: September 2021
## Author: Berk Gürdamar
##################################################

#!/bin/bash

echo "Creating conda environment"
conda env create -q -f scripts/depends.yaml

eval "$(conda shell.bash hook)"
conda activate rnaseq-pipeline

echo "Downloading pathfindR"
R CMD BATCH ./scripts/pathfindr_install.R
