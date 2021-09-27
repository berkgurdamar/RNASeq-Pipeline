##################################################
## Project: RNASeq-Pipeline
## Purpose: Setup file
## Date: September 2021
## Author: Berk Gürdamar
##################################################

#!/bin/bash

# eval "$(conda shell.bash hook)"
# conda activate rnaseq-pipeline

snakemake --configfile config.yaml --cores 16
