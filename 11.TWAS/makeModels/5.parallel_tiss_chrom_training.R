library(parallel)
suppressMessages(library(dplyr))
suppressMessages(library(glmnet))
suppressMessages((library(reshape2)))
suppressMessages(library(methods))

"%&%" <- function(a,b) paste(a,b, sep='')

argv <- commandArgs(trailingOnly = TRUE)
chrom <- argv[1]
tissue <- argv[2]
#tissue <- argv[3]
#celltype <- argv[2]
type <- "eQTL"
prefix <- "Model_training"

setwd("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/TWAS/bulk_eQTL/")
if (!dir.exists(tissue)) {
  dir.create(tissue)
}

gene_annot_file <- "/faststorage/project/cattle_gtexs/CattleGTEx/eQTL_annot.txt"
expression_file <- paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/eQTL/", tissue, "/transformed_expression.txt")
covariates_file <- paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/eQTL/", tissue, "/", tissue, "_LMM.auto_covar")
snp_annot_file <- paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/TWAS/ct_seQTL/", tissue, "/snp_chr/", tissue, "_chr", chrom, ".txt")
genotype_file <- paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/TWAS/ct_seQTL/", tissue, "/geno_chr/", tissue, ".genotype.chr", chrom, ".txt")
prefix <- "Model_training"

get_gene_annotation <- function(gene_annot_file_name, chr, gene_types=c('protein_coding', 'pseudogene', 'lncRNA', 'splicing', 'exon')){
  gene_df <- read.table(gene_annot_file_name, header = TRUE, stringsAsFactors = FALSE) %>%
    filter((chr == chrom) & gene_type %in% gene_types)
  gene_df
}
get_gene_expression <- function(gene_expression_file_name, gene_annot) {
  expr_df <- as.data.frame(t(read.table(gene_expression_file_name, header = T, stringsAsFactors = F, row.names = 1)))
  #expr_df <- expr_df %>% t() %>% as.data.frame()
  expr_df <- expr_df %>% select(one_of(intersect(gene_annot$gene_id, colnames(expr_df))))
  expr_df
}
####
gene_annot <- get_gene_annotation(gene_annot_file, chrom)
expr_df <- get_gene_expression(expression_file, gene_annot)
n_gene <- length(expr_df)


parallel_main = function(i) {
	source("/faststorage/project/cattle_gtexs/script/TWAS/makeModels/6.parallel_nested_cv_elnet.R")
	main(snp_annot_file, gene_annot, genotype_file, expr_df, covariates_file, as.numeric(chrom), prefix, tissue, type, as.numeric(i), null_testing=FALSE)
}


cl = makeCluster(11, type = "FORK")
clusterEvalQ(cl, c("dplyr", "glmnet", "reshape2", "methods"))
clusterExport(cl, c("snp_annot_file", "gene_annot", "genotype_file", "expr_df", "covariates_file", "prefix", "tissue", "type", "chrom"))
parLapply(cl, 1:n_gene, parallel_main)
stopCluster(cl)

