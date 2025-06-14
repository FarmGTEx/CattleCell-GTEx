suppressMessages(library(energy))
suppressMessages(library(dplyr))
suppressMessages(library(SeuratObject))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(remotes))
suppressMessages(library(devtools))
suppressMessages(library(Biobase))
suppressMessages(library(Seurat))
suppressMessages(library(future.apply))
suppressMessages(library(DWLS))

##prepare single cell reference
sc <- readRDS(list1[[i]])
Idents(sc) <- "CellType"
celltype <- unique(sc$CellType)
if ("Unknown cells" %in% celltype == TRUE) {
  sc <- subset(sc, idents = "Unknown cells", invert = TRUE)   
}
counts <- data.frame(sc@assays$RNA@counts)
meta <- data.frame(sc@meta.data)
phenoData <- data.frame(cbind(meta$CellType, meta$sample))
rownames(phenoData) <- rownames(meta)
names(phenoData) <- c("CellType","sample")
Signature <- buildSignatureMatrixMAST(scdata = counts, id = phenoData[,"CellType"], path = path[[i]], diff.cutoff = 0.5, pval.cutoff = 0.01)

##input bulk data (tpm)
tpm <- fread(list2[[i]])
tpm <- na.omit(tpm)

##run DWLS
res <- data.frame()
for (j in 1:ncol(tpm)) {
  b = setNames(tpm[,j], rownames(tpm))
  tr <- trimData(Signature, b)
  RES <- data.frame(t(solveDampenedWLS(tr$sig, tr$tpm)))
  res <- rbind(res,RES)
}
rownames(res) <- colnames(tpm)
res[res < 10^-5] <- 0
write.csv(res, list3[[i]])