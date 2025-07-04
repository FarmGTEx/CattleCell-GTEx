suppressMessages(library(MuSiC))
suppressMessages(library(energy))
suppressMessages(library(dplyr))
suppressMessages(library(SeuratObject))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(remotes))
suppressMessages(library(devtools))
suppressMessages(library(SCDC))
suppressMessages(library(Biobase))
suppressMessages(library(CDSeq))
suppressMessages(library(DeconRNASeq))
suppressMessages(library(Seurat))
suppressMessages(library(BisqueRNA))
suppressMessages(library(FARDEEP))
suppressMessages(library(EpiDISH))
source('/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/MammaryGland/method_benchmark/cibersort.R')

##prepare single cell reference data
#count to tpm
sc <- readRDS("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Cerebral cortex/Cerebral cortex_ref.rds")
Idents(sc) <- "CellType"
##Mammary gland, Ileum
sc <- subset(sc, idents = "Unknown cells", invert = TRUE)
#Heart
#sc <- subset(sc, idents = "Proliferative cells", invert = TRUE)
counts <- data.frame(sc@assays$RNA@counts)


##TPM
gtf = rtracklayer::import("/faststorage/project/cattle_gtexs/reference/Bos_taurus.ARS-UCD1.2.110.gtf")
class(gtf)
gtf = as.data.frame(gtf);dim(gtf)
table(gtf$type)
exon = gtf[gtf$type=="exon",
           c("start","end","gene_name")]
gle = lapply(split(exon,exon$gene_name),function(x){
  tmp=apply(x,1,function(y){
    y[1]:y[2]
  })
  length(unique(unlist(tmp)))
})
gle=data.frame(gene_name=names(gle),
               length=as.numeric(gle))
le = gle[match(rownames(counts),gle$gene_name),"length"]
counts$Length <- le
counts <- na.omit(counts)
kb <- counts$Length / 1000
x = ncol(counts)
countdata <- counts[,1:x-1]
rpk <- countdata / kb
rpk = rpk[complete.cases(rpk),]
tpm <- t(t(rpk)/colSums(rpk) * 1000000)

##gene expression of celltype
tpm <- t(tpm)
tpm <- data.frame(tpm)
meta <- data.frame(sc@meta.data)
meta <- droplevels(meta)
tpm <- tpm %>% mutate(cluster = meta$CellType)
#sc_tpm$cluster <- paste0("cluster", sc_tpm$cluster)
group <- list()
group <- split(tpm,tpm$cluster)
name<- rownames(table(tpm$cluster))
fenleimean <- colMeans(group[[1]][,-ncol(group[[1]])])
for (j in 2:length(group)) {
  flmean<- colMeans(group[[j]][,-ncol(group[[j]])])
  fenleimean <-rbind(fenleimean,flmean)
}
rownames(fenleimean) <- name
fenleimean = t(fenleimean)
tpm <- t(tpm[,-ncol(tpm)])


##CPM
suppressMessages(library(dplyr))
suppressMessages(library(edgeR))  
cpm <- edgeR::cpm(counts)

##gene expression of celltype
cpm <- t(cpm)
cpm <- data.frame(cpm)
meta <- data.frame(sc@meta.data)
meta <- droplevels(meta)
cpm <- cpm %>% mutate(cluster = meta$CellType)
#sc_tpm$cluster <- paste0("cluster", sc_tpm$cluster)
group <- list()
group <- split(cpm,cpm$cluster)
name<- rownames(table(cpm$cluster))
fenleimean <- colMeans(group[[1]][,-ncol(group[[1]])])
for (j in 2:length(group)) {
  flmean<- colMeans(group[[j]][,-ncol(group[[j]])])
  fenleimean <-rbind(fenleimean,flmean)
}
rownames(fenleimean) <- name
fenleimean = t(fenleimean)
cpm <- t(cpm[,-ncol(cpm)])


##counts of celltype
counts <- t(counts)
counts <- data.frame(counts)
meta <- data.frame(sc@meta.data)
meta <- droplevels(meta)
counts <- counts %>% mutate(cluster = meta$CellType)
#sc_tpm$cluster <- paste0("cluster", sc_tpm$cluster)
group <- list()
group <- split(counts,counts$cluster)
name<- rownames(table(counts$cluster))
fenleimean <- colMeans(group[[1]][,-ncol(group[[1]])])
for (j in 2:length(group)) {
  flmean<- colMeans(group[[j]][,-ncol(group[[j]])])
  fenleimean <-rbind(fenleimean,flmean)
}
rownames(fenleimean) <- name
fenleimean = t(fenleimean)
counts <- t(counts[,-ncol(counts)])


#meta data
meta <- data.frame(sc@meta.data)
phenoData <- data.frame(cbind(meta$CellType, meta$sample))
phenoData <- data.frame(phenoData)
rownames(phenoData) <- rownames(meta)
names(phenoData) <- c("CellType","sample")
celltype <- meta$CellType
phenoData$CellType <- celltype
## need to seperate the sample to 3 group (if there are not enough groups)
phenoData$sample[1:ceiling(nrow(phenoData)/3)] <- "sample1"
phenoData$sample[ceiling(nrow(phenoData)/3):ceiling(2*nrow(phenoData)/3)] <- "sample2"
phenoData$sample[ceiling(2*nrow(phenoData)/3):nrow(phenoData)] <- "sample3"
rownames(phenoData) <- gsub("-",".",rownames(phenoData))
rownames(phenoData) <- gsub("/",".",rownames(phenoData))
phenoData <- droplevels(phenoData)

##deconvolution pipeline
setwd("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Heart/method_benchmark/Self-pipeline/Counts/Binorm_Fraction/random_pseudo")
getwd()
file <- list.files(pattern = ".txt")
file <- file[-10]
list0 <- tools::file_path_sans_ext(file)
list1<-NULL
for(i in list0){
  list1[[i]]<-paste0("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Ileum/method_benchmark/Self-pipeline/Counts/Binorm_Fraction/random_pseudo/MuSic/Predict_", i, ".csv")
}
list2<-NULL
for(i in list0){
  list2[[i]]<-paste0("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Ileum/method_benchmark/Self-pipeline/Counts/Binorm_Fraction/random_pseudo/CDSeq/Predict_", i, ".csv")
}
list3<-NULL
for(i in list0){
  list3[[i]]<-paste0("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Ileum/method_benchmark/Self-pipeline/Counts/Binorm_Fraction/random_pseudo/Cibersort/Predict_", i, ".csv")
}
list4<-NULL
for(i in list0){
  list4[[i]]<-paste0("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Ileum/method_benchmark/Self-pipeline/Counts/Binorm_Fraction/random_pseudo/DeconRNASeq/Predict_", i, ".csv")
}
list5<-NULL
for(i in list0){
  list5[[i]]<-paste0("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Ileum/method_benchmark/Self-pipeline/Counts/Binorm_Fraction/random_pseudo/SCDC/Predict_", i, ".csv")
}
list6<-NULL
for(i in list0){
  list6[[i]]<-paste0("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Ileum/method_benchmark/Self-pipeline/Counts/Binorm_Fraction/random_pseudo/MuSic_marker/Predict_", i, ".csv")
}
list7<-NULL
for(i in list0){
  list7[[i]]<-paste0("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Heart/method_benchmark/Self-pipeline/Counts/Binorm_Fraction/random_pseudo/Bisque/Predict_", i, ".csv")
}
list8<-NULL
for(i in list0){
  list8[[i]]<-paste0("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Heart/method_benchmark/Self-pipeline/Counts/Binorm_Fraction/random_pseudo/EpiDISH/Predict_", i, ".csv")
}
list9<-NULL
for(i in list0){
  list9[[i]]<-paste0("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Heart/method_benchmark/Self-pipeline/Counts/Binorm_Fraction/random_pseudo/FARDEEP/Predict_", i, ".csv")
}


sc <- read.table("sc_ref.txt", sep = "\t", row.names = 1, header = T)  ##single cell marker file
marker <- unique(rownames(sc))
for (i in 1:length(list0)) {
  pseudo <- read.table(file[[i]], sep = "\t")
  pseudo1 <- ExpressionSet(as.matrix(pseudo))
  pseudo2 <- exprs(pseudo1)
  
  ##MuSic_ref
  sce <- SingleCellExperiment(list(counts=as.matrix(counts)),colData=DataFrame(cellType=phenoData$CellType,sampleID=phenoData$sample), 
                              metadata=list(phenoData))
  res1 = music_prop(bulk.mtx = pseudo2, # bulk exp
                            sc.sce = sce, # scRNAseq obj
                            clusters = 'cellType',  # cluster column
                            samples = 'sampleID', markers = NULL, normalize = FALSE, verbose = TRUE)$Est.prop.weighted
  write.csv(res1, list1[[i]])

  ##CDSeq
  res2 <- CDSeq(bulk_data = pseudo, reference_gep = fenleimean, cell_type_number = ncol(fenleimean), mcmc_iterations = 1000, block_number = 6, gene_subset_size=15)
  res2 <- t(res2$estProp)
  write.csv(res2, list2[[i]])

  ##Cibersort
  res3 <- CIBERSORT('sc_ref.txt', file[[i]], perm = 1000, QN = T)
  write.csv(res3, list3[[i]])

  ##DeconRNASeq
  res4 <- DeconRNASeq(datasets = pseudo, signatures = sc, proportions = NULL, checksig = FALSE, known.prop = FALSE, use.scale = TRUE, fig = FALSE)
  res4 <- res4$out.all
  rownames(res4) <- paste0("sample", c(1:1000))
  write.csv(res4, list4[[i]])

  ##SCDC
  sce <- ExpressionSet(tpm, phenoData = as(phenoData, "AnnotatedDataFrame"))
  res5 <- t(SCDC::SCDC_prop(bulk.eset = pseudo1, sc.eset = sce, ct.varname = "CellType", sample = "sample", ct.sub = levels(phenoData$CellType), iter.max = 200)$prop.est.mvw)
  res5 <- t(res5)
  write.csv(res5, list5[[i]])
  
  ##MuSic_marker
  sce <- SingleCellExperiment(list(counts=as.matrix(counts)),colData=DataFrame(cellType=phenoData$CellType,sampleID=phenoData$sample), 
                              metadata=list(phenoData))
  res6 = music_prop(bulk.mtx = pseudo2, # bulk exp
                            sc.sce = sce, # scRNAseq obj
                            clusters = 'cellType',  # cluster column
                            samples = 'sampleID', markers = marker, normalize = FALSE, verbose = TRUE)$Est.prop.weighted
  write.csv(res6, list1[[i]])
  
  ##BisqueRNA
  rownames(phenoData) <- gsub("-",".",rownames(phenoData))
  rownames(phenoData) <- gsub("/",".",rownames(phenoData))
  names(phenoData)[1] <- "cellType"
  names(phenoData)[2] <- "SubjectName"
  sce <- ExpressionSet(as.matrix(counts), phenoData = as(phenoData, "AnnotatedDataFrame"))
  res7 <- BisqueRNA::ReferenceBasedDecomposition(pseudo1, sce)$bulk.props
  write.csv(res7, list7[[i]])
  
  ##EpiDISH
  res8 <- t(EpiDISH::epidish(beta.m = pseudo, ref.m = fenleimean, method = "RPC")$estF)
  write.csv(res8, list8[[i]])
  
  ##FARDEEP
  res9 <- fardeep(fenleimean, pseudo, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)
  res9 <- t(res9$abs.beta)
  write.csv(res9, list9[[i]])
}


##normal distribution
setwd("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Ileum/method_benchmark/Self-pipeline/Counts/Normal_Fraction/random_pseudo")
getwd()
file <- list.files(pattern = ".txt")
file <- file[-4]
list0 <- tools::file_path_sans_ext(file)
list1<-NULL
for(i in list0){
  list1[[i]]<-paste0("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Ileum/method_benchmark/Self-pipeline/Counts/Normal_Fraction/random_pseudo/MuSic/Predict_", i, ".csv")
}
list2<-NULL
for(i in list0){
  list2[[i]]<-paste0("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Ileum/method_benchmark/Self-pipeline/Counts/Normal_Fraction/random_pseudo/CDSeq/Predict_", i, ".csv")
}
list3<-NULL
for(i in list0){
  list3[[i]]<-paste0("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Ileum/method_benchmark/Self-pipeline/Counts/Normal_Fraction/random_pseudo/Cibersort/Predict_", i, ".csv")
}
list4<-NULL
for(i in list0){
  list4[[i]]<-paste0("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Ileum/method_benchmark/Self-pipeline/Counts/Normal_Fraction/random_pseudo/DeconRNASeq/Predict_", i, ".csv")
}
list5<-NULL
for(i in list0){
  list5[[i]]<-paste0("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Ileum/method_benchmark/Self-pipeline/Counts/Normal_Fraction/random_pseudo/SCDC/Predict_", i, ".csv")
}
list6<-NULL
for(i in list0){
  list6[[i]]<-paste0("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Ileum/method_benchmark/Self-pipeline/Counts/Normal_Fraction/random_pseudo/MuSic_marker/Predict_", i, ".csv")
}
list7<-NULL
for(i in list0){
  list7[[i]]<-paste0("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Ileum/method_benchmark/Self-pipeline/Counts/Normal_Fraction/random_pseudo/Bisque/Predict_", i, ".csv")
}
list8<-NULL
for(i in list0){
  list8[[i]]<-paste0("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Heart/method_benchmark/Self-pipeline/Counts/Binorm_Fraction/random_pseudo/Bisque/Predict_", i, ".csv")
}

sc <- read.table("sc_ref.txt", sep = "\t", row.names = 1, header = T)

for (i in 1:length(list0)) {
  pseudo <- read.table(file[[i]], sep = "\t")
  pseudo1 <- ExpressionSet(as.matrix(pseudo))
  pseudo2 <- exprs(pseudo1)
  
  ##MuSic_ref
  sce <- SingleCellExperiment(list(counts=as.matrix(counts)),colData=DataFrame(cellType=phenoData$CellType,sampleID=phenoData$sample), 
                              metadata=list(phenoData))
  res1 = music_prop(bulk.mtx = pseudo2, # bulk exp
                            sc.sce = sce, # scRNAseq obj
                            clusters = 'cellType',  # cluster column
                            samples = 'sampleID', markers = NULL, normalize = FALSE, verbose = TRUE)$Est.prop.weighted
  write.csv(res1, list1[[i]])

  ##CDSeq
  res2 <- CDSeq(bulk_data = pseudo, reference_gep = fenleimean, cell_type_number = ncol(fenleimean), mcmc_iterations = 1000, block_number = 6, gene_subset_size=15)
  res2 <- t(res2$estProp)
  write.csv(res2, list2[[i]])

  ##Cibersort
  res3 <- CIBERSORT('sc_ref.txt', file[[i]], perm = 1000, QN = T)
  write.csv(res3, list3[[i]])

  ##DeconRNASeq
  res4 <- DeconRNASeq(datasets = pseudo, signatures = sc, proportions = NULL, checksig = FALSE, known.prop = FALSE, use.scale = TRUE, fig = FALSE)
  res4 <- res4$out.all
  rownames(res4) <- paste0("sample", c(1:1000))
  write.csv(res4, list4[[i]])

  ##SCDC
  sce <- ExpressionSet(tpm, phenoData = as(phenoData, "AnnotatedDataFrame"))
  res5 <- t(SCDC::SCDC_prop(bulk.eset = pseudo1, sc.eset = sce, ct.varname = "CellType", sample = "sample", ct.sub = levels(phenoData$CellType), iter.max = 200)$prop.est.mvw)
  res5 <- t(res5)
  write.csv(res5, list5[[i]])
  
  ##MuSic_marker
  sce <- SingleCellExperiment(list(counts=as.matrix(counts)),colData=DataFrame(cellType=phenoData$CellType,sampleID=phenoData$sample), 
                              metadata=list(phenoData))
  res6 = music_prop(bulk.mtx = pseudo2, # bulk exp
                            sc.sce = sce, # scRNAseq obj
                            clusters = 'cellType',  # cluster column
                            samples = 'sampleID', markers = marker, normalize = FALSE, verbose = TRUE)$Est.prop.weighted
  write.csv(res6, list1[[i]])
  
  ##BisqueRNA
  rownames(phenoData) <- gsub("-",".",rownames(phenoData))
  rownames(phenoData) <- gsub("/",".",rownames(phenoData))
  names(phenoData)[1] <- "cellType"
  names(phenoData)[2] <- "SubjectName"
  sce <- ExpressionSet(as.matrix(counts), phenoData = as(phenoData, "AnnotatedDataFrame"))
  res7 <- BisqueRNA::ReferenceBasedDecomposition(pseudo1, sce)$bulk.props
  write.csv(res7, list7[[i]])
  
  ##EpiDISH
  res8 <- t(EpiDISH::epidish(beta.m = pseudo, ref.m = fenleimean, method = "RPC")$estF)
  write.csv(res8, list8[[i]])
  
  ##FARDEEP
  res9 <- fardeep(fenleimean, pseudo, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)
  res9 <- t(res9$abs.beta)
  write.csv(res9, list9[[i]])
}


##uniform distribution
setwd("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Ileum/method_benchmark/Simbu")
sc <- read.table("sc_ref.txt", sep = "\t", row.names = 1, header = T)
marker <- unique(rownames(sc))
str <- "/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Ileum/method_benchmark/Simbu/"
pseudo <- read.csv("pseudo.txt", sep = "\t", row.names = 1)
#common_genes <- intersect(rownames(counts), rownames(pseudo))
#pseudo <- pseudo[rownames(pseudo) %in% common_genes, ]
#counts <- counts[rownames(counts) %in% common_genes, ]

##MuSic
pseudo1 <- ExpressionSet(as.matrix(pseudo))
pseudo2 <- exprs(pseudo1)
sce <- SingleCellExperiment(list(counts=as.matrix(counts)),colData=DataFrame(cellType=phenoData$CellType,sampleID=phenoData$sample), 
                              metadata=list(phenoData))
res1 = music_prop(bulk.mtx = pseudo2, # bulk exp
                            sc.sce = sce, # scRNAseq obj
                            clusters = 'cellType',  # cluster column
                            samples = 'sampleID', markers = NULL, normalize = FALSE, verbose = TRUE)$Est.prop.weighted
write.csv(res1, paste0(str, "MuSic_ref/Predict_Fraction_pseudo.csv"))

##CDSeq
res2 <- CDSeq(bulk_data = pseudo, reference_gep = fenleimean, cell_type_number = ncol(fenleimean), mcmc_iterations = 1000, block_number = 6, gene_subset_size=15)
res2 <- t(res2$estProp)
write.csv(res2, paste0(str, "CDSeq/Predict_CDSeq.csv"))

##Cibersort
res3 <- CIBERSORT('sc_ref_cpm.txt', "pseudo_simbu_cpm.txt", perm = 1000, QN = T)
write.csv(res3, paste0(str, "Cibersort/Predict_Fraction_pseudo.csv"))

##DeconRNASeq
res4 <- DeconRNASeq(datasets = pseudo, signatures = sc, proportions = NULL, checksig = FALSE, known.prop = FALSE, use.scale = TRUE, fig = FALSE)
res4 <- res4$out.all
rownames(res4) <- paste0("sample", c(1:1000))
write.csv(res4, paste0(str, "DeconRNASeq/Predict_Fraction_pseudo.csv"))

##SCDC
counts <- as.matrix(counts)
sce <- ExpressionSet(counts, phenoData = as(phenoData, "AnnotatedDataFrame"))
res5 <- t(SCDC::SCDC_prop(bulk.eset = pseudo1, sc.eset = sce, ct.varname = "CellType", sample = "sample", ct.sub = levels(phenoData$CellType), iter.max = 200)$prop.est.mvw)
res5 <- t(res5)
write.csv(res5, paste0(str, "SCDC/Predict_SCDC.csv"))

##MuSic_marker
res6 = music_prop(bulk.mtx = pseudo2, # bulk exp
                            sc.sce = sce, # scRNAseq obj
                            clusters = 'cellType',  # cluster column
                            samples = 'sampleID', markers = marker, normalize = FALSE, verbose = TRUE)$Est.prop.weighted
write.csv(res6, paste0(str, "MuSic_marker/Predict_Fraction_pseudo.csv"))

##Bisque
counts <- as.matrix(counts)
rownames(phenoData) <- gsub("-",".",rownames(phenoData))
rownames(phenoData) <- gsub("/",".",rownames(phenoData))
names(phenoData)[1] <- "cellType"
names(phenoData)[2] <- "SubjectName"
colnames(pseudo) <- paste0("sample", c(1:1000))
pseudo1 <- ExpressionSet(as.matrix(pseudo))
sce <- ExpressionSet(counts, phenoData = as(phenoData, "AnnotatedDataFrame"))
res7 <- BisqueRNA::ReferenceBasedDecomposition(pseudo1, sce)$bulk.props
write.csv(res7, paste0(str, "Bisque/Predict_Fraction_pseudo.csv"))

##EpiDISH
res8 <- t(EpiDISH::epidish(beta.m = pseudo, ref.m = fenleimean, method = "RPC")$estF)
write.csv(res8, paste0(str, "EpiDISH/Predict_Fraction_pseudo.csv"))

##FARDEEP
res9 <- fardeep(fenleimean, pseudo, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)
res9 <- t(res9$abs.beta)
write.csv(res9, paste0(str, "FARDEEP/Predict_Fraction_pseudo.csv"))