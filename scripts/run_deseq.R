##################################################
## Project: RNASeq
## Purpose: automated deseq2 script
## Date: March 2022
## Author: Berk Gürdamar
##################################################

mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

transcript_ids <- read.csv2(snakemake@input[["deseq_input"]][1], sep = "\t")[,1]

res <- biomaRt::getBM(attributes = c('ensembl_transcript_id',
                                     'ensembl_gene_id',
                                     'external_transcript_name',
                                     'external_gene_name'),
                      filters = 'ensembl_transcript_id',
                      values = transcript_ids,
                      mart = mart)

final_res <- res[, c(1,2)]
colnames(final_res) <- c("tx_name", "gene_id")


txi <- tximport::tximport(snakemake@input[["deseq_input"]], type = "salmon", tx2gene = final_res)

count_table <- as.data.frame(txi$counts)
colnames(count_table) <- snakemake@input[["deseq_input"]]

write.csv2(res[,c(1,2,4)], snakemake@params[["gene_name_converter"]])

# replicates --------------------------------------------------------------

sample <- snakemake@params[["sample_ids"]]

if(snakemake@params[["replicates"]] == "YES"){

  matched_samples <- snakemake@params[["matched_samples"]]

  sample_list <- list()
  for(i in matched_samples){
    idx <- strsplit(i, "_")[[1]]
    sample_list[[length(sample_list) + 1]] <- idx
  }

  num_vec <- rep(0, length(colnames(count_table)))
  num_vec[sapply(sample, function(x) grep(x, colnames(count_table)))] <- 1
  sample_names <- colnames(count_table)
  count_info <- as.data.frame(cbind(sample_names, num_vec))

  count_info$num_vec <- as.factor(count_info$num_vec)
  count_info$rep_info <- 1:nrow(count_info)

  group_id <- max(count_info$rep_info) + 1
  for(i in sample_list){
    for(j in i){
      count_info$rep_info[grep(j, count_info$sample_names)] <- group_id
    }
    group_id <- group_id + 1
  }

  deseq_obj <- DESeq2::DESeqDataSetFromMatrix(countData = round(count_table), colData = count_info, design = ~ num_vec)
  deseq_obj <- DESeq2::collapseReplicates(deseq_obj, count_info$rep_info, count_info$sample_names)

  write.csv2(deseq_obj@assays@data$counts, snakemake@params[["count_table"]])

}else{
  num_vec <- rep(0, length(colnames(count_table)))
  num_vec[sapply(sample, function(x) grep(x, colnames(count_table)))] <- 1
  sample_names <- colnames(count_table)
  count_info <- as.data.frame(cbind(sample_names, num_vec))

  count_info$num_vec <- as.factor(count_info$num_vec)

  deseq_obj <- DESeq2::DESeqDataSetFromMatrix(countData = round(count_table), colData = count_info, design = ~ num_vec)

  write.csv2(deseq_obj@assays@data$counts, snakemake@params[["count_table"]])
}

run_deseq <- DESeq2::DESeq(deseq_obj)
deseq_result <- as.data.frame(DESeq2::results(run_deseq,
                                              contrast = c("num_vec", "1", "0")))

deseq_result <- deseq_result[order(deseq_result$padj),]
deseq_result <-na.omit(deseq_result)


write.csv2(deseq_result, snakemake@output[["deseq_output"]])
