setwd("E:/Qiaolab/LXX/azoospermia/data")
file1 <- read.table("NOA_class1_tpm.txt",header = T,sep = "\t")
cell_type=sapply( strsplit(as.character(colnames(file1)), "_"), "[[", 1 )
LD <- file1[ ,which(cell_type %in% c("LD"))]
LD.log <- log2(LD/10+1)
df_pca <- prcomp( t(LD.log) )
plot(df_pca$x[,1], df_pca$x[,2])
df_out <- as.data.frame(df_pca$x)
df_out$sample <- sapply( strsplit(as.character(colnames(LD.log)), "_"), "[[", 1 )
df_out$MYH11 <- t(LD.log["MYH11", ])
df_out$ACTA2 <- t(LD.log["ACTA2", ])
df_out$INSL3 <- t(LD.log["INSL3", ])
df_out$DLK1 <- t(LD.log["DLK1", ])
library(ggplot2)
d <- dist(t(LD.log)) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit # view results
fit_out <- as.data.frame(fit$points)
colnames(fit_out) <- c("Dim1","Dim2")
ggplot(fit_out,aes(x=Dim1,y=Dim2))+geom_point()+theme()
library(ggplot2)
library(grid)
library(gridExtra)

mycolors<-c("red", "blue", "green", "orange", "dark green", "black", "chocolate", "purple", "brown", "dark blue")    

theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

p<-ggplot(df_out,aes(x=PC1,y=PC2,color=MYH11))
p<-p+geom_point()+theme
p


p<-ggplot(df_out,aes(x=PC1,y=PC2,color=ACTA2))
p<-p+geom_point()+theme
p

p<-ggplot(df_out,aes(x=PC1,y=PC2,color=INSL3))
p<-p+geom_point()+theme
p

p<-ggplot(df_out,aes(x=PC1,y=PC2,color=DLK1))
p<-p+geom_point()+theme
p



ensambl_id <- read.table("official_EnSG_genelist.txt",header = F,sep = "\t")
ensambl_id.de <- ensambl_id[!duplicated(ensambl_id$V1), ]
ensambl_id.de1 <- ensambl_id.de[!duplicated(ensambl_id.de$V2), ]
file1.sel <- file1[which(rownames(file1) %in% ensambl_id.de1$V1), ]
dim(file1.sel)
ensambl_id.use <- ensambl_id.de1[which(ensambl_id.de1$V1 %in% rownames(file1.sel)), ]
dim(ensambl_id.use)
seq <- as.vector(ensambl_id.use$V1)
file1.use <- file1.sel[seq, ]
rownames(file1.use) <- ensambl_id.use$V2
file1.use.log <- log2(file1.use/10+1)
write.table(file1.use.log,"NOA_Class1.tpm.new.sel.Ensembl.txt",quote = F,col.names = T,row.names = T,sep="\t")

file1_chart <- data.frame(Cell=colnames(file1.use.log),cell_type=sapply( strsplit(as.character(colnames(file1.use.log)), "_"), "[[", 1 ))
write.table(file1_chart,"NOA_Class1_chart.txt",quote = F,col.names = T,row.names = F,sep="\t")

###########class2
file2 <- read.table("NOA_class2_tpm.txt",header = T,sep = "\t")
colnames(file2)
file2.sel <- file2[which(rownames(file2) %in% ensambl_id.de1$V1), ]
dim(file2.sel)
ensambl_id.use <- ensambl_id.de1[which(ensambl_id.de1$V1 %in% rownames(file2.sel)), ]
dim(ensambl_id.use)
seq <- as.vector(ensambl_id.use$V1)
file2.use <- file2.sel[seq, ]
rownames(file2.use) <- ensambl_id.use$V2
file2.use.log <- log2(file2.use/10+1)
write.table(file2.use.log,"NOA_Class2.tpm.new.sel.Ensembl.txt",quote = F,col.names = T,row.names = T,sep="\t")

file2_chart <- data.frame(Cell=colnames(file2.use.log),cell_type=sapply( strsplit(as.character(colnames(file2.use.log)), "_"), "[[", 1 ))
write.table(file2_chart,"NOA_Class2_chart.txt",quote = F,col.names = T,row.names = F,sep="\t")

###########class3
file3 <- read.table("NOA_class3_tpm.txt",header = T,sep = "\t")
colnames(file3)
file3.sel <- file3[which(rownames(file3) %in% ensambl_id.de1$V1), ]
dim(file3.sel)
ensambl_id.use <- ensambl_id.de1[which(ensambl_id.de1$V1 %in% rownames(file3.sel)), ]
dim(ensambl_id.use)
seq <- as.vector(ensambl_id.use$V1)
file3.use <- file3.sel[seq, ]
rownames(file3.use) <- ensambl_id.use$V2
file3.use.log <- log2(file3.use/10+1)
write.table(file3.use.log,"NOA_class3.tpm.new.sel.Ensembl.txt",quote = F,col.names = T,row.names = T,sep="\t")

file3_chart <- data.frame(Cell=colnames(file3.use.log),cell_type=sapply( strsplit(as.character(colnames(file3.use.log)), "_"), "[[", 1 ))
write.table(file3_chart,"NOA_class3_chart.txt",quote = F,col.names = T,row.names = F,sep="\t")

###########class4
file4 <- read.table("NOA_class4_tpm.txt",header = T,sep = "\t")
colnames(file4)
file4.sel <- file4[which(rownames(file4) %in% ensambl_id.de1$V1), ]
dim(file4.sel)
ensambl_id.use <- ensambl_id.de1[which(ensambl_id.de1$V1 %in% rownames(file4.sel)), ]
dim(ensambl_id.use)
seq <- as.vector(ensambl_id.use$V1)
file4.use <- file4.sel[seq, ]
rownames(file4.use) <- ensambl_id.use$V2
file4.use.log <- log2(file4.use/10+1)
write.table(file4.use.log,"NOA_class4.tpm.new.sel.Ensembl.txt",quote = F,col.names = T,row.names = T,sep="\t")

file4_chart <- data.frame(Cell=colnames(file4.use.log),cell_type=sapply( strsplit(as.character(colnames(file4.use.log)), "_"), "[[", 1 ))
write.table(file4_chart,"NOA_class4_chart.txt",quote = F,col.names = T,row.names = F,sep="\t")
#############Normal
file5 <- read.table("D:/cyd/project/new_tpm.order_sm.tpm.txt",header = T,sep = "\t")
colnames(file5)
file5_chart <- data.frame(cell=colnames(file5),CellType=sapply( strsplit(as.character(colnames(file5)), "_"), "[[", 1 ))
file5_chart$CellType <- gsub("L1","L",file5_chart$CellType)
file5_chart$CellType <- gsub("L2","L",file5_chart$CellType)
DEGsfile5_chart$CellType <- gsub("L3","L",file5_chart$CellType)
file5_chart$CellType <- gsub("S1","S",file5_chart$CellType)
file5_chart$CellType <- gsub("S2","S",file5_chart$CellType)
file5_chart$CellType <- gsub("S3","S",file5_chart$CellType)
file5_chart$CellType <- gsub("S4","S",file5_chart$CellType)
write.table(file5_chart,"D:/cyd/project/CellType.txt",quote = F,col.names = T,row.names = F,sep = "\t")
DEGs <- read.table("E:/Qiaolab/LXX/azoospermia/data/Normal_cellType_DEGs.txt",header = F)
unique(DEGs$V1)
DEGs_tpm <- file5[which(rownames(file5) %in% DEGs$V1), ]
write.table(DEGs_tpm,"E:/Qiaolab/LXX/azoospermia/data/Normal_cellType_DEGs_tpm.txt",quote = F,col.names = T,row.names = T,sep = "\t")
ensambl_id <- read.table("official_EnSG_genelist.txt",header = F,sep = "\t")
ensambl_id.de <- ensambl_id[!duplicated(ensambl_id$V1), ]
ensambl_id.de1 <- ensambl_id.de[!duplicated(ensambl_id.de$V2), ]
file5.sel <- file5[which(rownames(file5) %in% ensambl_id.de1$V1), ]
dim(file5.sel)
ensambl_id.use <- ensambl_id.de1[which(ensambl_id.de1$V1 %in% rownames(file5.sel)), ]
dim(ensambl_id.use)
seq <- as.vector(ensambl_id.use$V1)
file5.use <- file5.sel[seq, ]
rownames(file5.use) <- ensambl_id.use$V2
file5.use.log <- log2(file5.use/10+1)
write.table(file5.use.log,"Normal.tpm.new.sel.Ensembl.txt",quote = F,col.names = T,row.names = T,sep="\t")

file5_chart <- data.frame(Cell=colnames(file5.use.log),cell_type=sapply( strsplit(as.character(colnames(file5.use.log)), "_"), "[[", 1 ))
write.table(file5_chart,"Normal_chart.txt",quote = F,col.names = T,row.names = F,sep="\t")
