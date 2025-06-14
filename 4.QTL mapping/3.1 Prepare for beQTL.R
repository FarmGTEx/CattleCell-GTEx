options(stringsAsFactors = FALSE)
library(edgeR)
library(preprocessCore)
library(RNOmni)
library(data.table)
library(R.utils)
library(SNPRelate)
library(dplyr)

##prepare for phenotype
setwd("/faststorage/project/cattle_gtexs/CattleGTEx/Cattle_bulk_exp")
file <- list.files(pattern = ".txt")
list0 <- tools::file_path_sans_ext(file)
list0 <- sapply(list0, function(x) unlist(strsplit(x, "\\_"))[2])
list0 <- data.frame(list0)
tissue <- list0$list0

list1<-NULL
for(i in tissue){
  list1[[i]]<-paste0("/faststorage/project/cattle_gtexs/CattleGTEx/Cattle_bulk_exp/Bulk_", i, ".txt")
}
list2<-NULL
for(i in tissue){
  list2[[i]]<-paste0("/faststorage/project/cattle_gtexs/CattleGTEx/Cattle_bulk_TPM/Bulk_", i, ".txt")
}

for (i in 1:length(tissue)) {
  bulk <- read.csv(list1[[i]], sep = "\t")
  bulk <- na.omit(bulk)
  TPM <- read.csv(list2[[i]], sep = "\t")
  TPM <- na.omit(TPM)
  if (ncol(bulk) >= 40) {
    Counts = counts = bulk
    samids = colnames(Counts) # sample id
    expr_counts = Counts
    expr = DGEList(counts=expr_counts) # counts
    nsamples = length(samids) # sample number
    ngenes = nrow(expr_counts) 
    y = calcNormFactors(expr, method="TMM")
    TMM = cpm(y,normalized.lib.sizes=T)

    count_threshold = 6
    tpm_threshold = 0.1
    sample_frac_threshold = 0.2
    sample_count_threshold = 10
    expr_tpm = TPM[rownames(expr_counts),samids]
    tpm_th = rowSums(expr_tpm >= tpm_threshold)
    count_th = rowSums(expr_counts >= count_threshold)
    ctrl1 = tpm_th >= (sample_frac_threshold * nsamples)
    ctrl2 = count_th >= (sample_frac_threshold * nsamples)
    mask = ctrl1 & ctrl2
    TMM_pass = TMM[mask,]
    rank_qnorm <- function(x) {
      result <- qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
      return(result)
    }
    TMM_inv = t(apply(TMM_pass, MARGIN = 1, FUN = rank_qnorm))

    region_annot <- gtf
    geneid = region_annot$gene_id
    expr_matrix = TMM_inv[rownames(TMM_inv) %in% geneid,]


    bed_annot = region_annot[region_annot$gene_id %in% rownames(expr_matrix),]
    bed = data.frame(bed_annot,expr_matrix[bed_annot$gene_id,])
    bed = bed[bed[,1] %in% as.character(1:30),]
    bed[,1] = as.numeric(bed[,1])
    bed = bed[order(bed[,1],bed[,2]),]
    colnames(bed)[1] = "#Chr"

    bed <- bed[bed$type == "gene",]
    bed <- bed[bed$gene_biotype == "protein_coding" | bed$gene_biotype == "lncRNA",]
    start = bed$start[bed$strand == "-"]
    end = bed$end[bed$strand == "-"]
    bed$start[bed$strand == "-"] = end
    bed$end[bed$strand == "-"] = start
    bed_unique <- bed[!duplicated(bed$gene_id), ]
    rownames(bed_unique) = bed_unique$gene_id
    bed_unique1 <- bed_unique[,-c(1:25)]
    bed_unique <- bed_unique[,-c(4:9,11:25)]
    colnames(bed_unique1) <- colnames(bulk)
    bed <- data.frame(bed_unique[,c(1:4)],bed_unique1)
    names(bed)[1:4] <- c("#Chr","start","end","gene_id")
    
    setwd(paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/eQTL/"))
    if (!dir.exists(tissue[[i]])) {
      dir.create(tissue[[i]])
    }
    setwd(paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/eQTL/", tissue[[i]]))
    fwrite(bed, file = "expr_tmm_inv.bed", sep = "\t")
    system("bgzip expr_tmm_inv.bed")
  }
}

