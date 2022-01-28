setwd("F:/Qiaolab/LXX/azoospermia/data")
Nor.TPM <- read.table("NOA_class1_tpm.txt",header = T,sep = "\t")
class <- sapply( strsplit(as.character(colnames(Nor.TPM)), "_"), "[[", 1 )
nosm.tpm <- Nor.TPM[ ,-which(class %in% c("ST","LD","T"))]
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
nbt <- CreateSeuratObject(nosm.tpm, project = "NOA", min.cells = 3, min.features = 200)
nbt
nbt[["percent.mt"]] <- PercentageFeatureSet(nbt, pattern = "^MT-")
nbt <- NormalizeData(nbt, normalization.method = "LogNormalize", scale.factor = 10000)
nbt <- NormalizeData(nbt)
nbt <- FindVariableFeatures(nbt, selection.method = "vst", nfeatures = 2000)

# plot variable features with and without labels
#VariableFeaturePlot(nbt)
length(nbt@assays$RNA@var.features)
var.features <- read.table("F:/Qiaolab/ACE2/human_sperm/azoospermia/Nor_vargene.txt",header = F)
nbt@assays$RNA@var.features <- as.character(var.features$V1)
all.genes <- rownames(nbt)
nbt <- ScaleData(nbt, features = all.genes)
nbt <- RunPCA(nbt, features = VariableFeatures(object = nbt))
DimPlot(nbt, reduction = "pca")
#DimHeatmap(nbt, dims = 1:15, cells = 500, balanced = TRUE)
nbt <- JackStraw(nbt, num.replicate = 100)
nbt <- ScoreJackStraw(nbt, dims = 1:20)
JackStrawPlot(nbt, dims = 1:15)
ElbowPlot(nbt)
nbt <- FindNeighbors(nbt, dims = 1:9)
nbt <- FindClusters(nbt, resolution = 0.3)
nbt <- RunUMAP(nbt, dims = 1:9)
DimPlot(nbt, reduction = "umap",label = T)
#FeaturePlot(nbt, features = c("GFRA1", "KIT", "STRA8", "SPO11", "OVOL2", "NME8", "TXNDC2", "TNP1", "PRM1","WT1","MYH11","CD68"))
umap_out <- data.frame(nbt@reductions$umap@cell.embeddings)
meta_out <- data.frame(nbt@meta.data)
umap_out$Stage=sapply( strsplit(as.character(rownames(umap_out)), "_"), "[[", 1 )
umap_out$cluster <- meta_out$seurat_clusters
umap_out$cell_type <- umap_out$Stage
meta_out <- data.frame(nbt@meta.data)

ggplot(umap_out,aes(x=UMAP_1,y=UMAP_2,color= cell_type))+
  geom_point(size=1.5)+theme_bw()+
  theme(panel.grid=element_blank())+scale_color_manual(values=c("#9a2ff1","#fbc1c8","#3af97e","#fd0d91","#c89420","#fb2f27","#0cfff9","#4377fe","#97fa97","#ccbd6e","#fd69b6","#6bb7fd","#0039ff","#a42831","#9ad030","#FFD92F"),limits=c("SPG1","SPG2","SPG3","L","Z","P","D","SEC","S"),labels=c("SSC","Diff.ing SPG","Diff.ed SPG","L","Z","P","D","SPC7","S"))
#Diffusion PseudoTime  Analysis
library(destiny)
library(ggplot2)# ╪сть destiny...
library(Biobase)
library(Seurat)
ct <-GetAssayData(object = nbt)
ct<-ct[VariableFeatures(nbt),]
ct <- as.ExpressionSet(as.data.frame(t(ct)))
dm <- DiffusionMap(ct,k = 50)
diffusionMap_DATA <- data.frame(DC1=dm$DC1,DC2=dm$DC2)
diffusionMap_DATA$Stage=sapply( strsplit(as.character(rownames(diffusionMap_DATA)), "_"), "[[", 1 )
pdf("../plot/Class1_diffusionMap.pdf",width = 5,height = 3.5)
ggplot(diffusionMap_DATA, aes(DC1, DC2, color = Stage))+geom_point()+theme_bw()+theme(panel.grid=element_blank())+scale_color_manual(values = c("#9E0142","#D8434D","#F67D4A","#FDC070","#F2EA91","#C1E69F","#77C8A4","#388FBA","#5E4FA2"),limits=c("SPG1","SPG2","SPG3","L","Z","P","D","SEC","S"),labels=c("SSC","Diff.ing SPG","Diff.ed SPG","L","Z","P","D","SPC7","S"))
dev.off()
diffusionMap_DATA.use <- diffusionMap_DATA[ ,c(1:2)]
########slingshot
library(slingshot, quietly = FALSE)
nor.TPM.use1 <- as.matrix(nosm.tpm)
sim <- SingleCellExperiment(assays = List(counts = nor.TPM.use1))
geneFilter <- apply(assays(sim)$counts,1,function(x){
  sum(x >= 3) >= 10})
sim <- sim[geneFilter, ]

# full quantile normalization
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(sim)$norm <- FQnorm(assays(sim)$counts)
diffusionMap_DATA.use1 <- as.matrix(diffusionMap_DATA.use)
reducedDims(sim) <- SimpleList(DC =diffusionMap_DATA.use1)

meta_out$cluster  <- as.numeric("1")
meta_out[which(meta_out$orig.ident %in% c("SPG2")), ]$cluster <- as.numeric("2")
meta_out[which(meta_out$orig.ident %in% c("SPG3")), ]$cluster <- as.numeric("3")
meta_out[which(meta_out$orig.ident %in% c("L")), ]$cluster <- as.numeric("4")
meta_out[which(meta_out$orig.ident %in% c("Z")), ]$cluster <- as.numeric("5")
meta_out[which(meta_out$orig.ident %in% c("P")), ]$cluster <- as.numeric("6")
meta_out[which(meta_out$orig.ident %in% c("D")), ]$cluster <- as.numeric("7")
meta_out[which(meta_out$orig.ident %in% c("SEC")), ]$cluster <- as.numeric("8")
meta_out[which(meta_out$orig.ident %in% c("S")), ]$cluster <- as.numeric("9")

colData(sim)$GMM <-meta_out$cluster
sim <- slingshot(sim, clusterLabels = 'GMM', reducedDim = 'DC')
summary(sim$slingPseudotime_1)

library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(9)

plot(reducedDims(sim)$DC, col = colors, pch=16, asp = 1)
lines(SlingshotDataSet(sim), lwd=2, col='black')

plot(reducedDims(sim)$DC, col = rainbow(9)[sim$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(sim), lwd=2, type = 'lineages', col = 'black')
lin1 <- getLineages(nbt@reductions$umap@cell.embeddings, meta_out$cluster, start.clus= '1', end.clus = '9')

lin1 <- getLineages(nbt@reductions$umap@cell.embeddings, meta_out$cluster, start.clus = '1')
## Using full covariance matrix

lin1

crv1 <- getCurves(lin1)
crv1

plot(reducedDims(sim)$UMAP, col = rainbow(9)[sim$GMM], pch=16, asp = 1)
lines(crv1, lwd = 3, col = 'black')

time <- data.frame(sample=sim$slingshot@NAMES,time=sim$slingPseudotime_1,Celltype=sapply( strsplit(as.character(sim$slingshot@NAMES), "_"), "[[", 1 ))
time$UMAP_1=umap_out$UMAP_1
time$UMAP_2=umap_out$UMAP_2
rownames(time) <- time$sample
time.ord <- time[order(time$time), ]
time.ord$num <- c(1:2543)
ggplot(time,aes(x=UMAP_1,y=UMAP_2,color= time))+
  geom_point(size=1.5)+theme_bw()+
  theme(panel.grid=element_blank())#+scale_color_manual(values=c("#9a2ff1","#fbc1c8","#3af97e","#fd0d91","#c89420","#fb2f27","#0cfff9","#4377fe","#97fa97","#ccbd6e","#fd69b6","#6bb7fd","#0039ff","#a42831","#9ad030","#FFD92F"),limits=c("SPG1","SPG2","SPG3","L","Z","P","D","SEC","S1","S2","S3","S4"),labels=c("SSC","Diff.ing SPG","Diff.ed SPG","L","Z","P","D","SPC7","S1","S2","S3","S4"))


library(tradeSeq)
# fit negative binomial GAM
sim <- fitGAM(sim)
# test for dynamic expression
ATres <- associationTest(sim)


topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]
pst.ord <- order(sim$slingPseudotime_1, na.last = NA)
heatdata <- assays(sim)$counts[topgenes, pst.ord]
heatclus <- sim$GMM[pst.ord]

heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])



