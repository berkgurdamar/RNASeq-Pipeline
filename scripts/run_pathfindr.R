##################################################
## Project: RNASeq
## Purpose: automated pathfindR script
## Date: April 2022
## Author: Berk Gürdamar
##################################################

deseq_out <- read.csv2(snakemake@input[["pathfindr_input"]])

gene_name_converter <- read.csv2(snakemake@params[["gene_name_converter"]], row.names = 1)[,2:3]
gene_name_converter <- unique(gene_name_converter)

colnames(deseq_out)[1] <- colnames(gene_name_converter)[1]

deseq_out <- deseq_out[,c(1,3,7)]

pathfindr_input <- dplyr::left_join(deseq_out, gene_name_converter, by = "ensembl_gene_id")[, c(4,2,3)]



# pathfindR ---------------------------------------------------------------

n_iter <- snakemake@params[["iterations"]]

RA_processed <- pathfindR::input_processing(input = pathfindr_input,
                                            p_val_threshold = snakemake@params[["p_val_threshold"]],
                                            pin_name_path  = "Biogrid",
                                            convert2alias = TRUE)


biocarta_list <- pathfindR::fetch_gene_set(gene_sets = "KEGG",
                                           min_gset_size = 10,
                                           max_gset_size = 300)
biocarta_gsets <- biocarta_list[[1]]
biocarta_descriptions <- biocarta_list[[2]]

dirs <- c()
for (i in base::seq_len(n_iter)) {
  dir_i <- file.path(snakemake@output[["pathfindr_output"]], "active_snw_searches", paste0("Iteration_", i))
  dir.create(dir_i, recursive = TRUE)
  dirs <- c(dirs, dir_i)
}

geneInitProbs <- seq(from = 0.01, to = 0.2, length.out = n_iter)

`%dopar%` <- foreach::`%dopar%`
combined_res <- NULL

cl <- parallel::makeCluster(snakemake@params[["threads"]], setup_strategy = "sequential")
doParallel::registerDoParallel(cl)

combined_res <- foreach::foreach(i = 1:n_iter, .combine = rbind) %dopar% {

  active_snws <- pathfindR::active_snw_search(input_for_search = RA_processed,
                                              pin_name_path = "Biogrid",
                                              score_quan_thr = 0.8,
                                              sig_gene_thr = 0.02,
                                              search_method = "GR",
                                              geneInitProbs = geneInitProbs[i],
                                              snws_file = paste0("active_snws_", i),
                                              dir_for_parallel_run = dirs[i])


  current_res <- pathfindR::enrichment_analyses(snws = active_snws,
                                                sig_genes_vec = RA_processed$GENE,
                                                pin_name_path = "Biogrid",
                                                genes_by_term = biocarta_gsets,
                                                term_descriptions = biocarta_descriptions,
                                                adj_method = "bonferroni",
                                                enrichment_threshold = 0.05,
                                                list_active_snw_genes = TRUE)


  current_res
}

parallel::stopCluster(cl)

summarized_df <- pathfindR::summarize_enrichment_results(combined_res,
                                              list_active_snw_genes = TRUE)


final_res <- pathfindR::annotate_term_genes(result_df = summarized_df,
                                            input_processed = RA_processed,
                                            genes_by_term = biocarta_gsets)


write.csv2(final_res, file.path(snakemake@output[["pathfindr_output"]], "enrichment_result.csv"), row.names = F)
