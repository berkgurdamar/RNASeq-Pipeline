##################################################
## Project: RNASeq
## Purpose: automated deseq2 script
## Date: July 2021
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


write.csv2(count_table, snakemake@params[["count_table"]])



sample <- snakemake@params[["sample_ids"]]


num_vec <- rep(0, length(colnames(count_table)))
num_vec[sapply(sample, function(x) grep(x, colnames(count_table)))] <- 1
sample_names <- colnames(count_table)
count_info <- as.data.frame(cbind(sample_names, num_vec))

count_info$num_vec <- as.factor(count_info$num_vec)

write.csv2(count_info, snakemake@params[["info_table"]])

deseq_obj <- DESeq2::DESeqDataSetFromMatrix(countData = round(count_table), colData = count_info, design = ~ num_vec)
run_deseq <- DESeq2::DESeq(deseq_obj)
deseq_result <- as.data.frame(DESeq2::results(run_deseq,
                                              contrast = c("num_vec", "1", "0")))

deseq_result <- deseq_result[order(deseq_result$padj),]
deseq_result <-na.omit(deseq_result)


write.csv2(deseq_result, snakemake@output[["deseq_output"]])
