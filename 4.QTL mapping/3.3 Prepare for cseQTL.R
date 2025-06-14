###cell type expression
path <- "/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/cell_specific/"
tissue <- list.dirs(path, full.names = TRUE, recursive = FALSE)
tissue <- sapply(tissue, function(x) unlist(strsplit(x, "\\/"))[9])
tissue <- data.frame(tissue)
tissue <- tissue$tissue
path<-NULL
for(i in tissue){
  path[[i]]<-paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/cell_specific/", i, "/")
}

for (i in 1:length(tissue)) {
  load(paste0("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Bulk/Result_bMIND/", tissue[[i]], "_bMIND.RData"))
  dat <- data_list$A
  CT <- colnames(dat)
  for (j in 1:dim(dat1)[2]) {  
    select_data <- data.frame(dat[, j, ])
    #select_data2 <- data.frame(dat2[, j, ])
    #select_data1 <- select_data1[rownames(select_data1) %in% rownames(select_data2),]
    #select_data <- cbind(select_data1, select_data2)
    file_name <- paste0("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Bulk/Result_bMIND/CT_exp/", tissue[[i]], "_split_", CT[[j]], ".csv")
    write.csv(select_data, file_name)
  }
}

##construct expression bed file
gtf = rtracklayer::import("/faststorage/project/cattle_gtexs/reference/Bos_taurus.ARS-UCD1.2.110.gtf")
class(gtf)
gtf = as.data.frame(gtf);dim(gtf)
table(gtf$type)
exon = gtf[gtf$type=="exon",
           c("start","end","gene_id")]
gle = lapply(split(exon,exon$gene_id),function(x){
  tmp=apply(x,1,function(y){
    y[1]:y[2]
  })
  length(unique(unlist(tmp)))
})
gle=data.frame(gene_id=names(gle),
               length=as.numeric(gle))


for (i in 1:length(tissue)) {
  setwd("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Bulk/Result_bMIND/CT_exp")
  ct <- list.files(pattern = tissue[[i]])
  ct <- sapply(ct, function(x) unlist(strsplit(x, "\\_"))[3])
  ct <- sapply(ct, function(x) unlist(strsplit(x, "\\."))[1])
  ct <- data.frame(ct)
  ct <- ct$ct
  frac <- read.csv(paste0("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Bulk/Results_DWLS/", tissue[[i]], "/Results/Predict_", tissue[[i]], "_DWLS.csv"), row.names = 1)[-c(601:700),]
  sample_id <- rownames(frac)
  for (j in 1:length(ct)) {
    setwd("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Bulk/Result_bMIND/CT_exp")
    bulk <- read.csv(paste0(tissue[[i]], "_split_", ct[[j]], ".csv"), row.names = 1)
    colnames(bulk) <- rownames(frac)
    Counts = counts = bulk
    le = gle[match(rownames(counts),gle$gene_id),"length"]
    counts$Length <- le
    counts <- na.omit(counts)
    kb <- counts$Length / 1000
    x = ncol(counts)
    countdata <- counts[,1:x-1]
    rpk <- countdata / kb
    rpk = rpk[complete.cases(rpk),]
    tpm <- t(t(rpk)/colSums(rpk) * 1000000)
    TPM <- tpm

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
    
    setwd("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/cell_specific")
    if (!dir.exists(tissue[[i]])) {
      dir.create(tissue[[i]])
    }
    setwd(paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/cell_specific/", tissue[[i]]))
    if (!dir.exists(ct[[j]])) {
      dir.create(ct[[j]])
    }
    setwd(paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/cell_specific/", tissue[[i]], "/", ct[[j]]))
    fwrite(bed, file = "expr_tmm_inv.bed", sep = "\t")
  }
}