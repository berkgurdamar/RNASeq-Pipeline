#!/bin/bash

##################################################
## Project: RNASeq-Pipeline
## Purpose: Setup file
## Date: April 2022
## Author: Berk Gürdamar
##################################################

echo "Creating conda environment"
conda env create -q -f env_yaml/rnaseq_pipeline.yml

eval "$(conda shell.bash hook)"
conda activate rnaseq-pipeline

echo "Downloading pathfindR"
Rscript -e "devtools::install_github('egeulgen/pathfindR', force = T)"
