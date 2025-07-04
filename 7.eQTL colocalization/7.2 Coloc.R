options(stringsAsFactors = FALSE)
library(edgeR)
library(preprocessCore)
library(RNOmni)
library(data.table)
library(R.utils)
library(SNPRelate)
library(dplyr)
library(coloc)

path <- "/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/eQTL/"
tissue <- list.dirs(path, full.names = TRUE, recursive = FALSE)
tissue <- sapply(tissue, function(x) unlist(strsplit(x, "\\/"))[9])
tissue <- data.frame(tissue)
tissue <- tissue$tissue

path<-NULL
for(i in tissue){
  path[[i]]<-paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/eQTL/", i, "/ieQTL/")
}

for (i in 3:length(tissue)) {
    setwd(path[[i]])
    ct <- list.dirs(full.names = TRUE, recursive = FALSE)
    ct <- sapply(ct, function(x) unlist(strsplit(x, "\\/"))[2])
    ct <- data.frame(ct)
    ct <- ct$ct
    qtl1 <- fread(paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/Coloc/coloc/Bulk/", tissue[[i]], "_coloc.bed"))
    for (j in 1:length(ct)) {
      if (file.exists(paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/Coloc/coloc/all_cell_components/", tissue[[i]], "_", ct[[j]], "_coloc.bed"))) {
        qtl2 <- fread(paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/Coloc/coloc/all_cell_components/", tissue[[i]], "_", ct[[j]], "_coloc.bed"))
        gene1 <- unique(qtl1$gene_id)
        gene2 <- unique(qtl2$gene_id)
        gene_list <- intersect(gene1,gene2)
        all_summary <- data.frame()
        if (length(gene_list) > 0){
          for (k in 1:length(gene_list)) {
            gene_qtl1 <- qtl1[qtl1$gene_id == gene_list[[k]], ]
            gene_qtl2 <- qtl2[qtl2$gene_id == gene_list[[k]], ]
            input <- merge(gene_qtl1, gene_qtl2, by="rs_id", all=FALSE, suffixes=c("_bulk_eqtl","_celltype_iqtl"))
            if (nrow(input) > 0){
              input <- input[!duplicated(input$rs_id),]
              input$varbeta1 <- input$slope_se_bulk_eqtl ^ 2
              input$varbeta2 <- input$slope_se_celltype_iqtl ^ 2
              #input$beta <- abs(input$beta)
              #input$slope <- abs(input$slope)
              res <- coloc.abf(dataset1=list(snp = input$rs_id, pvalues=input$pval_nominal_bulk_eqtl, beta=input$slope_bulk_eqtl, varbeta=input$varbeta1, N=input$N_bulk_eqtl, type="quant", MAF=input$maf_bulk_eqtl),
                          dataset2=list(snp = input$rs_id, pvalues=input$pval_nominal_celltype_iqtl, beta=input$slope_celltype_iqtl, varbeta=input$varbeta2, N=input$N_celltype_iqtl, type="quant", MAF=input$maf_celltype_iqtl))
              summary <- data.frame(res$summary)
              summary$gene_id <- gene_list[[k]]
              all_summary <- rbind(all_summary, summary)
            }
          }
        }
        setwd("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/Coloc/coloc/Results/Bulk_components/")
        write.csv(all_summary, paste0(tissue[[i]], "_", ct[[j]], "_coloc.csv"))
      }
    }   
}