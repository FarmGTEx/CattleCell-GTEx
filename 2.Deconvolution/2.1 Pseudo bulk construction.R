suppressMessages(library(Seurat))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(dplyr))
suppressMessages(library(FamilyRank))

sc <- readRDS("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Cerebral cortex/Cerebral cortex_pseudo.rds")
Idents(sc) <- "CellType"
#sc <- subset(sc, idents = "Proliferative cells", invert = TRUE)
#sc <- subset(sc, idents = "Unknown cells", invert = TRUE)
counts <- data.frame(sc@assays$RNA@counts)

### TPM  
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


### CPM
library(edgeR)  
cpm <- edgeR::cpm(counts)
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


### Count
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



## Binormal distribution
#cell <- data.frame(table(Idents(sc)))
#real <- cell$Freq / sum(cell$Freq)
#real <- sort(real)
combFunc <- function(...) {
    mapply('bind_cols', ..., SIMPLIFY=FALSE)
}
best_params <- NULL
best_similarity <- -Inf
sd1 <- c(0.1,0.5,0.9)
sd2 <- c(0.1,0.5,0.9)

for (i in 1:length(sd1)) {
  for (j in 1:length(sd2)) {
    fraction <- array(data = NA,dim = c(20,1000))
    for (k in 1:1000) {
      simulated_data <- abs(rbinorm(20, 0, 1, sd1[[i]], sd2[[j]], 0.5))
      simulated_data <- simulated_data / sum(simulated_data)
      fraction[,k] <- simulated_data
    } 
    rownames(fraction) <- paste0("cell", c(1:20))
    colnames(fraction) <- paste0("sample", c(1:1000))
    write.csv(fraction, paste0("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Ileum/method_benchmark/CPM/Binorm_Fraction/Fraction_", sd1[[i]], "_", sd2[[j]], ".csv")) 
  }
}

#produce pseudo bulk expression
setwd("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/MammaryGland/method_benchmark/Self-pipeline/Counts/Binorm_Fraction")
file <- list.files(pattern = ".csv")
list0 <- tools::file_path_sans_ext(file)
list1<-NULL
for(i in list0){
  list1[[i]]<-paste0("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/MammaryGland/method_benchmark/Self-pipeline/Counts/Binorm_Fraction/random_pseudo/", i, "_pseudo.txt")
}

for (i in 1:length(file)) {
  pseudo <- array(data = NA,dim = c(nrow(fenleimean),1000))
  fraction <- read.csv(file[[i]], row.names = 1)
  for (j in 1:1000) {
    for (k in 1:nrow(fenleimean)) {
      gene <- sum(fenleimean[k,] * fraction[,j])
      pseudo[k,j] <- gene
    }
  }
  rownames(pseudo) <- rownames(fenleimean)
  colnames(pseudo) <- paste0("sample", c(1:1000))
  write.table(pseudo, list1[[i]], sep = "\t")
}
      

## normal distribution
combFunc <- function(...) {
    mapply('bind_cols', ..., SIMPLIFY=FALSE)
}
best_params <- NULL
best_similarity <- -Inf
sd <- c(0.1,0.5,0.9)

for (i in 1:length(sd)) {
  fraction <- array(data = NA,dim = c(20,1000))
  for (k in 1:1000) {
    simulated_data <- abs(rnorm(20, 0, sd[[i]]))
    simulated_data <- simulated_data / sum(simulated_data)
    fraction[,k] <- simulated_data
  }
  rownames(fraction) <- paste0("cell", c(1:20))
  colnames(fraction) <- paste0("sample", c(1:1000))
  write.csv(fraction, paste0("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Ileum/method_benchmark/Self-pipeline/Normal_Fraction/Fraction_sd_", sd[[i]], ".csv")) 
}

#produce pseudo bulk expression
setwd("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/MammaryGland/method_benchmark/Self-pipeline/Counts/Normal_Fraction")
file <- list.files(pattern = ".csv")
list0 <- tools::file_path_sans_ext(file)
list1<-NULL
for(i in list0){
  list1[[i]]<-paste0("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/MammaryGland/method_benchmark/Self-pipeline/Counts/Normal_Fraction/random_pseudo/", i, "_pseudo.txt")
}

for (i in 1:length(file)) {
  pseudo <- array(data = NA,dim = c(nrow(fenleimean),1000))
  fraction <- read.csv(file[[i]], row.names = 1)
  for (j in 1:1000) {
    for (k in 1:nrow(fenleimean)) {
      gene <- sum(fenleimean[k,] * fraction[,j])
      pseudo[k,j] <- gene
    }
  }
  rownames(pseudo) <- rownames(fenleimean)
  colnames(pseudo) <- paste0("sample", c(1:1000))
  write.table(pseudo, list1[[i]], sep = "\t")
}



## uniform distribution
combFunc <- function(...) {
    mapply('bind_cols', ..., SIMPLIFY=FALSE)
}
best_params <- NULL
best_similarity <- -Inf
fraction <- array(data = NA,dim = c(20,1000))
for (k in 1:1000) {
  simulated_data <- runif(20, 1, 99)
  simulated_data <- simulated_data / sum(simulated_data)
  fraction[,k] <- simulated_data
}
rownames(fraction) <- paste0("cell", c(1:20))
colnames(fraction) <- paste0("sample", c(1:1000))
write.csv(fraction, "/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Heart/method_benchmark/Self-pipeline/Uniform_Fraction/Farction.csv")

#produce pseudo bulk expression
setwd("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Cerebral cortex/method_benchmark/Self-pipeline/CPM/Uniform_Fraction/")
pseudo <- array(data = NA,dim = c(nrow(fenleimean),1000))
fraction <- read.csv("Fraction.csv", row.names = 1)
for (j in 1:1000) {
  for (k in 1:nrow(fenleimean)) {
    gene <- sum(fenleimean[k,] * fraction[,j])
    pseudo[k,j] <- gene
  }
}
rownames(pseudo) <- rownames(fenleimean)
colnames(pseudo) <- paste0("sample", c(1:1000))
write.table(pseudo, "/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Cerebral cortex/method_benchmark/Self-pipeline/CPM/Uniform_Fraction/random_pseudo/pseudo.txt", sep = "\t")