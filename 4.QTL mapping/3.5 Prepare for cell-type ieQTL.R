options(stringsAsFactors = FALSE)
library(edgeR)
library(preprocessCore)
library(RNOmni)
library(data.table)
library(R.utils)
library(SNPRelate)
library(dplyr)

#extract cell components and corresponding expression matrix
path <- "/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Bulk/Results_DWLS/"
tissue <- list.dirs(path, full.names = TRUE, recursive = FALSE)
tissue <- sapply(tissue, function(x) unlist(strsplit(x, "\\/"))[10])
tissue <- data.frame(tissue)
tissue_sc <- tissue$tissue
#bulk
path <- "/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/eQTL/"
tissue <- list.dirs(path, full.names = TRUE, recursive = FALSE)
tissue <- sapply(tissue, function(x) unlist(strsplit(x, "\\/"))[9])
tissue <- data.frame(tissue)
tissue_bulk <- tissue$tissue
tissue <- intersect(tissue_sc, tissue_bulk)

list1<-NULL
for(i in tissue){
  list1[[i]]<-paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/eQTL/", i, "/expr_tmm_inv.bed.gz")
}
list2<-NULL
for(i in tissue){
  list2[[i]]<-paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/eQTL/", i, "/ieQTL/")
}
list3<-NULL
for(i in tissue){
  list3[[i]]<-paste0("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Bulk/Results_DWLS/", i, "/Results/Predict_", i, "_DWLS.csv")
}

#remove low quality cell type and samples
for (i in 1:length(tissue)) {
  frac <- read.csv(list3[[i]], row.names = 1)
  frac <- na.omit(frac)
  bed <- fread(list1[[i]])
  sample <- colnames(bed)[-c(1:4)]
  frac <- frac[rownames(frac) %in% sample, ]
  zero_ratio <- colSums(frac == 0) / nrow(frac)
  frac1 <- frac[, zero_ratio <= 0.8]
  mean_pro <- colMeans(frac1)
  frac2 <- frac1[, mean_pro >= 0.05]
  ct <- colnames(frac2)
  for (j in 1:length(ct)) {
    setwd(list2[[i]])
    bed <- fread(list1[[i]])
    if (!dir.exists(ct[[j]])) {
        dir.create(ct[[j]])
    }
    setwd(paste0(list2[[i]], ct[[j]]))
    select_ct <- data.frame(frac2[,j])
    rownames(select_ct) <- rownames(frac2)
    names(select_ct)[1] <- "fraction"
    sd <- sd(select_ct$fraction)
    mean <- mean(select_ct$fraction)
    select_ct$select_sample = select_ct$fraction < mean + 3*sd &
          select_ct$fraction > mean - 3*sd
    select_ct <- select_ct[select_ct$select_sample == TRUE,]
    select_frac <- data.frame(select_ct[,-2])
    rownames(select_frac) <- rownames(select_ct)
    names(select_frac)[1] <- "fraction"
    bed1 <- data.frame(bed[,-c(1:4)])
    bed1 <- bed1[,colnames(bed1) %in% rownames(select_frac)]
    bed <- data.frame(bed[,1:4],bed1)
    names(bed)[1] <- "#Chr"
    fwrite(bed, paste0(ct[[j]], "_expr_tmm_inv.bed"))
    write.table(select_frac, paste0(ct[[j]],"_frac.txt"), sep = "\t", col.names = F, row.names = T, quote = F)
  }
}