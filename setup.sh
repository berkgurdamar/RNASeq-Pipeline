#!/bin/bash

##################################################
## Project: RNASeq-Pipeline
## Purpose: Setup file
## Date: March 2022
## Author: Berk Gürdamar
##################################################

echo "Creating conda environment"
conda env create -q -f env_yaml/rnaseq_pipeline.yml

eval "$(conda shell.bash hook)"
conda activate rnaseq-pipeline

echo "Downloading multiqc"
pip install multiqc

echo "Downloading pathfindR"
Rscript -e "devtools::install_github('egeulgen/pathfindR')"
