##################################################
## Project: RNASeq
## Purpose: automated pathfindR script
## Date: March 2022
## Author: Berk Gürdamar
##################################################

deseq_out <- read.csv2(snakemake@input[["pathfindr_input"]])

gene_name_converter <- read.csv2(snakemake@params[["gene_name_converter"]], row.names = 1)[,2:3]
gene_name_converter <- unique(gene_name_converter)

colnames(deseq_out)[1] <- colnames(gene_name_converter)[1]

deseq_out <- deseq_out[,c(1,3,7)]

pathfindr_input <- dplyr::left_join(deseq_out, gene_name_converter, by = "ensembl_gene_id")[, c(4,2,3)]

pathfindr_result <- pathfindR::run_pathfindR(pathfindr_input,
                                             plot_enrichment_chart = F,
                                             visualize_enriched_terms = F,
                                             n_processes = snakemake@params[["threads"]],
                                             p_val_threshold = snakemake@params[["p_val_threshold"]],
                                             iterations = snakemake@params[["iterations"]],
                                             output_dir = snakemake@output[["pathfindr_output"]])
