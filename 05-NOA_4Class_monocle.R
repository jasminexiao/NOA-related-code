setwd("E:/Qiaolab/LXX/azoospermia/data")###R3.3.2
file1 <- read.table("NOA_class1_tpm.txt",header = T,sep = "\t")
file1.class <- sapply( strsplit(as.character(colnames(file1)), "_"), "[[", 1 )
file1.use <- file1[ ,-which(file1.class %in% c("ST","LD","T"))]
colnames(file1.use) <- paste("Class1",colnames(file1.use),sep = "_")

file2 <- read.table("NOA_class2_tpm.txt",header = T,sep = "\t")
file2.class <- sapply( strsplit(as.character(colnames(file2)), "_"), "[[", 1 )
file2.use <- file2[ ,-which(file2.class %in% c("ST","LD","T"))]
colnames(file2.use) <- paste("Class2",colnames(file2.use),sep = "_")

file3 <- read.table("NOA_class3_tpm.txt",header = T,sep = "\t")
file3.class <- sapply( strsplit(as.character(colnames(file3)), "_"), "[[", 1 )
file3.use <- file3[ ,-which(file3.class %in% c("ST","LD","T"))]
colnames(file3.use) <- paste("Class3",colnames(file3.use),sep = "_")

file4 <- read.table("NOA_class4_tpm.txt",header = T,sep = "\t")
file4.class <- sapply( strsplit(as.character(colnames(file4)), "_"), "[[", 1 )
file4.use <- file4[ ,-which(file4.class %in% c("ST","LD","T"))]
colnames(file4.use) <- paste("Class4",colnames(file4.use),sep = "_")

file5 <- read.table("D:/cyd/project/new_tpm.order_sm.tpm.txt",header = T,sep = "\t")
file5.class <- sapply( strsplit(as.character(colnames(file5)), "_"), "[[", 1 )
file5.use <- file5[ ,-which(file5.class %in% c("ST","LD","T"))]
colnames(file5.use) <- paste("Nor",colnames(file5.use),sep = "_")

ALL <- cbind(file5.use,file1.use)
ALL.log <- log2(ALL/10+1)
nor.log <- log2(file5.use/10+1)
##############MONOCLE
library(monocle)
library(ggplot2)
library(devtools)
#install_github("satijalab/seurat@da6cd08")
library(Seurat) 
colne=data.frame(sample=colnames(nor.log),Class=sapply( strsplit(as.character(colnames(nor.log)), "_"), "[[",1 ),CellType=sapply( strsplit(as.character(colnames(nor.log)), "_"), "[[", 2 ))


sample_sheet <- colne
rownames(sample_sheet) <- sample_sheet[,1]
expr=nor.log
gene_annotation <- data.frame(gene=rownames(expr))
rownames(gene_annotation) <- rownames(expr)
expr_matrix <- expr[,as.vector(as.character(sample_sheet[ ,1]))]


pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
# First create a CellDataSet from the relative expression levels
#HSMM <- newCellDataSet(as.matrix(expr_matrix),
#                       phenoData = pd,
#                       featureData = fd)

# Next, use it to estimate RNA counts
#rpc_matrix <- relative2abs(HSMM)

# Now, make a new CellDataSet using the RNA counts
HSMM <- newCellDataSet(as.matrix(expr_matrix),
                       phenoData = pd, 
                       featureData = fd,
                       lowerDetectionLimit=1,
                       expressionFamily=negbinomial())
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 1)
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 100))##old,800
length(expressed_genes)
#holds the identifiers for genes expressed in at least 20 cells of the data set.
print(head(pData(HSMM)))

nbt.monocle <- CreateSeuratObject(raw.data =nor.log, min.cells = 3, min.genes = 2000)

#nbt <- FilterCells(object = nbt, subset.names = c("nGene"), low.thresholds = c(2000))
nbt.monocle <- NormalizeData(object = nbt.monocle, normalization.method = "LogNormalize",scale.factor = 10000)
nbt.monocle <- FindVariableGenes(object = nbt.monocle, mean.function = ExpMean, dispersion.function = LogVMR,
                                 x.low.cutoff = 0.4, y.cutoff = 0.4)######0.1&0.3
length(nbt.monocle@var.genes)

ordering_genes <- intersect(nbt.monocle@var.genes, expressed_genes)
length(ordering_genes)
ordering_genes <- nbt.monocle@var.genes
HSMM <- setOrderingFilter(HSMM, ordering_genes)
#HSMM <- reduceDimension(HSMM, use_irlba=FALSE)
HSMM <- reduceDimension(HSMM)
HSMM <- orderCells(HSMM, num_paths=1, reverse=F)
#HSMM <- orderCells(HSMM, num_paths=13)
#num_paths is determined by the num of cell types
HSMM$CellType <- factor(HSMM$CellType,levels = c("SPG1","SPG2","SPG3","L1","L2","L3","Z","P","D","SEC","S1","S2","S3","S4"))

pdf("../plot/Vivo_monocle_type.pdf",width=5,height = 3.5)
plot_cell_trajectory(HSMM,color_by="CellType",cell_size = 1,show_branch_points = F)+scale_color_manual(values=c("#9a2ff1","#fbc1c8","#3af97e","#fd0d91","#236687","#c1bfc4","#c89420","#fb2f27","#0cfff9","#4377fe","#97fa97","#ccbd6e","#fd69b6","#6bb7fd"),limits=c("SPG1","SPG2","SPG3","L1","L2","L3","Z","P","D","SEC","S1","S2","S3","S4"),labels=c("SSC","Diff.ing SPG","Diff.ed SPG","L1","L2","L3","Z","P","D","SPC7","S1","S2","S3","S4"))
dev.off()

plot_cell_trajectory(HSMM,color_by="Class",cell_size = 1,show_branch_points = F)






pdf("../plot/Vivo_monocle_celltype.pdf",width=8,height = 3.5)
plot_cell_trajectory(HSMM,color_by="Cell_type",cell_size = 1,show_branch_points = F)+scale_color_manual(values=c("S1.TE"="#DE9DD6","S2.MTE"="#3D00E6","S2.PTE"="#FACA82","S3.MTE1"="#FF5EEC","S3.MTE2"="pink","S3.PTE"="#61C2BD"))
dev.off()

#pdf("../plot/monocle_stage.pdf",width=5,height = 3.5)#+scale_color_manual(values = c("#FF5500","#FFD000","#003CFF"),limits=c("S1.TE","S2.TE","S3.TE") )
plot_cell_trajectory(HSMM,color_by="Stage",cell_size = 1,show_branch_points = F)+scale_color_manual(values=c("S1"="#FF5500","S2"="#FFD000","S3"="#003CFF"))


pdf("../plot/Vivo_monocle_location.pdf",width=3.5,height = 3.5)
plot_cell_trajectory(HSMM,color_by="Location",cell_size = 1,show_branch_points = F)#+scale_color_manual(values = c("tomato","cyan"),limits=c("M","P") )
dev.off()
plot_cell_trajectory(HSMM,color_by="Embryo",cell_size = 1,show_branch_points = F)+scale_color_manual(values = mycolors)
plot_cell_trajectory(HSMM,color_by="Pseudotime",cell_size = 1,show_branch_points = F)

mycolors <- c("pink","gray","aquamarine2","green","red","cyan","forestgreen","navy","darkorange2","cornflowerblue","deeppink","skyblue","darkturquoise","deepskyblue4",
              "dimgray","greenyellow","khaki1","darkolivegreen1",
              "black","darkgreen","darkred","darksalmon","darkslateblue","darkturquoise","deepskyblue",
              "purple","black","gray")
plot_cell_trajectory(HSMM,color_by="Embryo",cell_size = 3,show_branch_points = F)+scale_color_manual(values = mycolors )
#pdf("../plot/monocle_type.pdf",width=3.5,height = 3.5)
#plot_cell_trajectory(HSMM,color_by="Type",cell_size = 1,show_branch_points = F)#+scale_color_manual(values = c("#96C4E6","#E560CD","#FAAE4E"),limits=c("TE","MTE","PTE") )
pdf("../plot/monocle_MTE_PTE.pdf",width=3.5,height = 3.5)
plot_cell_trajectory(HSMM,color_by="Type",cell_size = 1,show_branch_points = F)+scale_color_manual(values = c("#6fbce4","#E560CD","#FAAE4E"),limits=c("mixTE","MTE","PTE") )
dev.off()
####################CELL CYCLE
library(ggplot2)
library(ggridges)
library(viridis)
MTE_PTE.CLASS <- sapply( strsplit(as.character(colnames(Stage_TE_TPM)), "_"), "[[", 1 )
MTE_PTE.TPM <- Stage_TE_TPM[ ,which(MTE_PTE.CLASS %in% c("S2.MTE","S3.MTE","S3.PTE"))]
colnames(MTE_PTE.TPM)
MTE_PTE_TPM.log <- log2(MTE_PTE.TPM/10+1)
G1_S_gene <- read.table("G:/tang lab/Preparation/G1_S_gene.txt", header = T, sep = "\t")
G1_S_TPM.log <- MTE_PTE_TPM.log[G1_S_gene$Gene, ]
G1_S_TPM.log_mean <- data.frame(G1_S=colMeans(G1_S_TPM.log))
G1_S_TPM.log_mean$Sample <- rownames(G1_S_TPM.log_mean)
G1_S_TPM.log_mean$Type <- sapply( strsplit(as.character(G1_S_TPM.log_mean$Sample), "_"), "[[", 1 )
ggplot(G1_S_TPM.log_mean, aes(x = G1_S, y = Type, fill = Type, color = Type)) + 
  geom_density_ridges(scale = 4, size = 1) + 
  scale_fill_cyclical(
    values = c("#E560CD", "#FAAE4E"), guide = "legend",
    name = "Type"
  ) +
  scale_color_cyclical(
    values = c("#bd1cf9", "#f45e0b"), guide = "legend",
    labels = c("Fair" = "blue w/ black outline", "Good" = "green w/ yellow outline"),
    name = "Type"
  )

G2_M_gene <- read.table("G:/tang lab/Preparation/G2_M_gene.txt", header = T, sep = "\t")
G2_M_TPM.log <- MTE_PTE_TPM.log[G2_M_gene$Gene,]
G2_M_TPM.log_mean <- data.frame(G2_M=colMeans(G2_M_TPM.log))
G2_M_TPM.log_mean$Sample <- rownames(G2_M_TPM.log_mean)
G2_M_TPM.log_mean$Type <- sapply( strsplit(as.character(G2_M_TPM.log_mean$Sample), "_"), "[[", 1 )
ggplot(G2_M_TPM.log_mean, aes(x = G2_M, y = Type, fill = Type, color = Type)) + 
  geom_density_ridges(scale = 4, size = 1) + 
  scale_fill_cyclical(
    values = c("#E560CD", "#FAAE4E"), guide = "legend",
    name = "Type"
  ) +
  scale_color_cyclical(
    values = c("#bd1cf9", "#f45e0b"), guide = "legend",
    labels = c("Fair" = "blue w/ black outline", "Good" = "green w/ yellow outline"),
    name = "Type"
  )

G1_S_TPM.log_mean$Cell_cycle <- c("G1_S")
colnames(G1_S_TPM.log_mean) <- c("Exp","Sample","Type","Cell_cycle")

G2_M_TPM.log_mean$Cell_cycle <- c("G2_M")
colnames(G2_M_TPM.log_mean) <- c("Exp","Sample","Type","Cell_cycle")

cell_cycle_TPM.log_mean <- rbind(G1_S_TPM.log_mean,G2_M_TPM.log_mean)
cell_cycle_TPM.log_mean$cell_type <- paste(cell_cycle_TPM.log_mean$Cell_cycle,cell_cycle_TPM.log_mean$Type)
cell_cycle_TPM.log_mean$cell_type1 <- NA
cell_cycle_TPM.log_mean[which(cell_cycle_TPM.log_mean$cell_type %in% c("G1_S S2.MTE","G1_S S3.MTE")), ]$cell_type1 <- c("G1_S_MTE")
cell_cycle_TPM.log_mean[which(cell_cycle_TPM.log_mean$cell_type %in% c("G1_S S3.PTE")), ]$cell_type1 <- c("G1_S_PTE")
cell_cycle_TPM.log_mean[which(cell_cycle_TPM.log_mean$cell_type %in% c("G2_M S2.MTE","G2_M S3.MTE")), ]$cell_type1 <- c("G2_M_MTE")
cell_cycle_TPM.log_mean[which(cell_cycle_TPM.log_mean$cell_type %in% c("G2_M S3.PTE")), ]$cell_type1 <- c("G2_M_PTE")

pdf("../plot/MTE_PTE_cell_cycle.pdf",width = 5.3,height = 3.3)
ggplot(cell_cycle_TPM.log_mean, aes(x = `Exp`, y = `cell_type1`, fill = ..x..,colors=cell_type1)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis(name = "Relative Expression", option = "C") +
  labs(title = 'Cell Cycle')+theme_bw()+theme(panel.grid=element_blank())
dev.off()
