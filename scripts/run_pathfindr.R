##################################################
## Project: RNASeq
## Purpose: automated pathfindR script
## Date: July 2021
## Author: Berk Gürdamar
##################################################

deseq_out <- read.csv2(snakemake@input[["pathfindr_input"]])
gene_Name_converter <- read.csv2(snakemake@params[["gene_name_converter"]], sep = ",")


colnames(deseq_out)[1] <- colnames(gene_Name_converter)[1]


deseq_out <- deseq_out[,c(1,3,7)]
deseq_out[,1] <- sapply(deseq_out[,1], function(x) strsplit(x, "\\.")[[1]][1])

pathfindr_in <- dplyr::left_join(deseq_out, gene_Name_converter, by = "Gene.stable.ID")[, c(4,2,3)]



neoadj_pathfindr_result <- pathfindR::run_pathfindR(pathfindr_in,
                                                    n_processes = snakemake@params[["threads"]],
                                                    p_val_threshold = snakemake@params[["p_val_threshold"]],
                                                    iterations = snakemake@params[["iterations"]],
                                                    output_dir = snakemake@output[["pathfindr_output"]])
