setwd("E:/Qiaolab/LXX/azoospermia/data")
######NOA_Class1
file1 <- read.table("NOA_class1_tpm.txt",header = T,sep = "\t")
colnames(file1)
file1.class <- sapply( strsplit(as.character(colnames(file1)), "_"), "[[", 1 )
unique(file1.class)
table(file1.class)
sample <- matrix(1,24153,12)
for (i in 1:12) {
  sample[ ,i] <- apply(file1[ ,which(file1.class %in% unique(file1.class)[i])], 1, mean)
}
file1.use <- as.data.frame(sample)
rownames(file1.use) <- rownames(file1)
colnames(file1.use) <- unique(file1.class)
######NOA_class2
file2 <- read.table("NOA_class2_tpm.txt",header = T,sep = "\t")
colnames(file2)
file2.class <- sapply( strsplit(as.character(colnames(file2)), "_"), "[[", 1 )
unique(file2.class)
table(file2.class)
sample <- matrix(1,24153,12)
for (i in 1:12) {
  sample[ ,i] <- apply(file2[ ,which(file2.class %in% unique(file2.class)[i])], 1, mean)
}
file2.use <- as.data.frame(sample)
rownames(file2.use) <- rownames(file2)
colnames(file2.use) <- unique(file2.class)
######NOA_class3
file3 <- read.table("NOA_class3_tpm.txt",header = T,sep = "\t")
colnames(file3)
file3.class <- sapply( strsplit(as.character(colnames(file3)), "_"), "[[", 1 )
unique(file3.class)
sample <- matrix(1,24153,7)
for (i in 1:7) {
  sample[ ,i] <- apply(file3[ ,which(file3.class %in% unique(file3.class)[i])], 1, mean)
}
file3.use <- as.data.frame(sample)
rownames(file3.use) <- rownames(file3)
colnames(file3.use) <- unique(file3.class)
######NOA_class4
file4 <- read.table("NOA_class4_tpm.txt",header = T,sep = "\t")
colnames(file4)
file4.class <- sapply( strsplit(as.character(colnames(file4)), "_"), "[[", 1 )
unique(file4.class)
sample <- matrix(1,24153,9)
for (i in 1:9) {
  sample[ ,i] <- apply(file4[ ,which(file4.class %in% unique(file4.class)[i])], 1, mean)
}
file4.use <- as.data.frame(sample)
rownames(file4.use) <- rownames(file4)
colnames(file4.use) <- unique(file4.class)
######Normal
file5 <- read.table("D:/cyd/project/new_tpm.order_sm.tpm.txt",header = T,sep = "\t")
colnames(file5) <- gsub("L1","L",colnames(file5))
colnames(file5) <- gsub("L2","L",colnames(file5))
colnames(file5) <- gsub("L3","L",colnames(file5))
colnames(file5) <- gsub("S1","S",colnames(file5))
colnames(file5) <- gsub("S2","S",colnames(file5))
colnames(file5) <- gsub("S3","S",colnames(file5))
colnames(file5) <- gsub("S4","S",colnames(file5))
file5.class <- sapply( strsplit(as.character(colnames(file5)), "_"), "[[", 1 )

unique(file5.class)
sample <- matrix(1,24153,12)
for (i in 1:12) {
  sample[ ,i] <- apply(file5[ ,which(file5.class %in% unique(file5.class)[i])], 1, mean)
}
file5.use <- as.data.frame(sample)
rownames(file5.use) <- rownames(file5)
colnames(file5.use) <- unique(file5.class)

##########OMIM genes
genes <- as.vector(c("DMRT1","TAF4B","DAZL","TEX11","BUB1B","NANOS1","DDX3Y","RPS4Y2","PRDM9","MDC1","SYCP3","SYCE1","MEIOB","FKBP6","HSPA2","BOLL","PIWIL1","STRC","DPY19L2","RSPH1","SLC26A8","CCDC39","CLGN","ESR2","SPATA16","AURKC","NSUN7","CATSPER2","SEPT12","TUBB8","ZMYND15","DNAJB13","TEKT2","KLHL10","CATSPER1","ACE","SMCP","PGR","CAMK4","SPEM1","NT5C1B","NPAS2","SOX9","NR5A1","DCXR","NR0B1","CLCN2","USP9Y","KDM5D","AR","PLP1","PIWIL1","PIWIL2","PIWIL3","PIWIL4"))
file5.use.select <- file5.use[genes, ]
mean <- rowMeans(file5.use.select)
sd_value <- apply(file5.use.select,1,sd)

z_score <- (file5.use.select-mean)/sd_value
z_score[z_score>=1] <- 1
z_score[z_score<= -1] <- -1
z_score<- na.omit(z_score)
#SPG.markers.use.tpm.log=log2(SPG.markers.use.tpm/10+1)
ubj1 <- sapply( strsplit(as.character(colnames(z_score)), "/t"), "[[", 1 )

aaka22 = data.frame(sapply( strsplit(as.character(colnames(z_score)), "_"), "[[", 1 ))

row.names(aaka22)<-ubj1
colnames(aaka22)<-c("class")
aka6 <- list(class=c("SPG1"="#9a2ff1","SPG2"="#fbc1c8","SPG3"="#3af97e","L"="#fd0d91","Z"="#c89420","P"="#fb2f27","D"="#0cfff9","SEC"="#4377fe","S"="#97fa97","ST"="#0039ff","LD"="#a42831","T"="#9ad030"))
library("pheatmap")#
pheatmap( z_score,annotation_col =aaka22,annotation_legend = T,annotation_colors =aka6,gaps_row=c(8,30,42,51),
          color = colorRampPalette(c("#4CC7F2","#7AA1D2","#EBD9ED","#E8A9ED", "#EF3FE3"))(100), cluster_row =F,
          cluster_col =F,show_rownames =T, 
          show_colnames =F, scale = "none",legend =T,border_color =c("white"),
          cellwidth =10,cellheight = 5.5,
          fontsize_row=6,fontsize=8,fontsize_col = 5)
######NOA Pheatmap
#####Class1
file1.use.select <- file1.use[genes, ]
mean <- rowMeans(file1.use.select)
sd_value <- apply(file1.use.select,1,sd)

z_score <- (file1.use.select-mean)/sd_value
z_score[z_score>=1] <- 1
z_score[z_score<= -1] <- -1
z_score<- na.omit(z_score)
#SPG.markers.use.tpm.log=log2(SPG.markers.use.tpm/10+1)
ubj1 <- sapply( strsplit(as.character(colnames(z_score)), "/t"), "[[", 1 )

aaka22 = data.frame(sapply( strsplit(as.character(colnames(z_score)), "_"), "[[", 1 ))

row.names(aaka22)<-ubj1
colnames(aaka22)<-c("class")
aka6 <- list(class=c("SPG1"="#9a2ff1","SPG2"="#fbc1c8","SPG3"="#3af97e","L"="#fd0d91","Z"="#c89420","P"="#fb2f27","D"="#0cfff9","SEC"="#4377fe","S"="#97fa97","ST"="#0039ff","LD"="#a42831","T"="#9ad030"))
library("pheatmap")#
pheatmap( z_score,annotation_col =aaka22,annotation_legend = T,annotation_colors =aka6,gaps_row=c(8,30,42,51),
          color = colorRampPalette(c("#4CC7F2","#7AA1D2","#EBD9ED","#E8A9ED", "#EF3FE3"))(100), cluster_row =F,
          cluster_col =F,show_rownames =T, 
          show_colnames =F, scale = "none",legend =T,border_color =c("white"),
          cellwidth =10,cellheight = 5.5,
          fontsize_row=6,fontsize=8,fontsize_col = 5)

#####Class2
file2.use.select <- file2.use[genes, ]
file2.use.select1 <- file2.use.select[ c(1:6)]
mean <- rowMeans(file2.use.select1)
sd_value <- apply(file2.use.select1,1,sd)

z_score <- (file2.use.select1-mean)/sd_value
z_score[z_score>=1] <- 1
z_score[z_score<= -1] <- -1
z_score<- na.omit(z_score)
#SPG.markers.use.tpm.log=log2(SPG.markers.use.tpm/10+1)
ubj1 <- sapply( strsplit(as.character(colnames(z_score)), "/t"), "[[", 1 )

aaka22 = data.frame(sapply( strsplit(as.character(colnames(z_score)), "_"), "[[", 1 ))

row.names(aaka22)<-ubj1
colnames(aaka22)<-c("class")
aka6 <- list(class=c("SPG1"="#9a2ff1","SPG2"="#fbc1c8","SPG3"="#3af97e","L"="#fd0d91","Z"="#c89420","P"="#fb2f27","D"="#0cfff9","SEC"="#4377fe","S"="#97fa97","ST"="#0039ff","LD"="#a42831","T"="#9ad030"))
library("pheatmap")#
pheatmap( z_score,annotation_col =aaka22,annotation_legend = T,annotation_colors =aka6,gaps_row=c(8,30,42,51),
          color = colorRampPalette(c("#4CC7F2","#7AA1D2","#EBD9ED","#E8A9ED", "#EF3FE3"))(100), cluster_row =F,
          cluster_col =F,show_rownames =T, 
          show_colnames =F, scale = "none",legend =T,border_color =c("white"),
          cellwidth =10,cellheight = 5.5,
          fontsize_row=6,fontsize=8,fontsize_col = 5)

#####Class3
file3.use.select <- file3.use[genes, ]
file3.use.select1 <- file3.use.select[ ,c(1:5)]
mean <- rowMeans(file3.use.select1)
sd_value <- apply(file3.use.select1,1,sd)

z_score <- (file3.use.select1-mean)/sd_value
z_score[z_score>=1] <- 1
z_score[z_score<= -1] <- -1
z_score<- na.omit(z_score)
#SPG.markers.use.tpm.log=log2(SPG.markers.use.tpm/10+1)
ubj1 <- sapply( strsplit(as.character(colnames(z_score)), "/t"), "[[", 1 )

aaka22 = data.frame(sapply( strsplit(as.character(colnames(z_score)), "_"), "[[", 1 ))

row.names(aaka22)<-ubj1
colnames(aaka22)<-c("class")
aka6 <- list(class=c("SPG1"="#9a2ff1","SPG2"="#fbc1c8","SPG3"="#3af97e","L"="#fd0d91","Z"="#c89420","P"="#fb2f27","D"="#0cfff9","SEC"="#4377fe","S"="#97fa97","ST"="#0039ff","LD"="#a42831","T"="#9ad030"))
library("pheatmap")#
pheatmap( z_score,annotation_col =aaka22,annotation_legend = T,annotation_colors =aka6,gaps_row=c(8,30,42,51),
          color = colorRampPalette(c("#4CC7F2","#7AA1D2","#EBD9ED","#E8A9ED", "#EF3FE3"))(100), cluster_row =F,
          cluster_col =F,show_rownames =T, 
          show_colnames =F, scale = "none",legend =T,border_color =c("white"),
          cellwidth =10,cellheight = 5.5,
          fontsize_row=6,fontsize=8,fontsize_col = 5)

#####Class4
file4.use.select <- file4.use[genes, ]
mean <- rowMeans(file4.use.select)
sd_value <- apply(file4.use.select,1,sd)

z_score <- (file4.use.select-mean)/sd_value
z_score[z_score>=1] <- 1
z_score[z_score<= -1] <- -1
z_score<- na.omit(z_score)
#SPG.markers.use.tpm.log=log2(SPG.markers.use.tpm/10+1)
ubj1 <- sapply( strsplit(as.character(colnames(z_score)), "/t"), "[[", 1 )

aaka22 = data.frame(sapply( strsplit(as.character(colnames(z_score)), "_"), "[[", 1 ))

row.names(aaka22)<-ubj1
colnames(aaka22)<-c("class")
aka6 <- list(class=c("SPG1"="#9a2ff1","SPG2"="#fbc1c8","SPG3"="#3af97e","L"="#fd0d91","Z"="#c89420","P"="#fb2f27","D"="#0cfff9","SEC"="#4377fe","S"="#97fa97","ST"="#0039ff","LD"="#a42831","T"="#9ad030"))
library("pheatmap")#
pheatmap( z_score,annotation_col =aaka22,annotation_legend = T,annotation_colors =aka6,gaps_row=c(8,30,42,51),
          color = colorRampPalette(c("#4CC7F2","#7AA1D2","#EBD9ED","#E8A9ED", "#EF3FE3"))(100), cluster_row =F,
          cluster_col =F,show_rownames =T, 
          show_colnames =F, scale = "none",legend =T,border_color =c("white"),
          cellwidth =10,cellheight = 5.5,
          fontsize_row=6,fontsize=8,fontsize_col = 5)
