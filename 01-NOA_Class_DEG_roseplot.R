###########Class1
file <- read.table("F:/Qiaolab/LXX/azoospermia/data/Nor_Class1_DEG/rose_input.txt",header = T,sep = "\t")
library(reshape2)
rose_long <- melt(file,id.vars="Class",variable.name="CellType",value.name="count")
library(ggplot2)
pdf("E:/Qiaolab/LXX/azoospermia/plot/Class1_DEGs_roseplot.pdf",width = 3.8,height = 3.8)
p <- ggplot(rose_long,aes(x=CellType, y=count, fill=Class))+ geom_bar(stat="identity",color="black")+coord_polar()+scale_fill_manual(values = c("red","blue"),limits=c("Up","Down"))
p
dev.off()
###########Class2
file <- read.table("E:/Qiaolab/LXX/azoospermia/data/Nor_Class2_DEG/rose_input.txt",header = T,sep = "\t")
library(reshape2)
rose_long <- melt(file,id.vars="Class",variable.name="CellType",value.name="count")
library(ggplot2)
pdf("E:/Qiaolab/LXX/azoospermia/plot/Class2_DEGs_roseplot.pdf",width = 3.8,height = 3.8)
p <- ggplot(rose_long,aes(x=CellType, y=count, fill=Class))+ geom_bar(stat="identity",color="black")+coord_polar()+scale_fill_manual(values = c("red","blue"),limits=c("Up","Down"))
p
dev.off()

###########Class3
file <- read.table("E:/Qiaolab/LXX/azoospermia/data/Nor_Class3_DEG/rose_input.txt",header = T,sep = "\t")
library(reshape2)
rose_long <- melt(file,id.vars="Class",variable.name="CellType",value.name="count")
library(ggplot2)
pdf("E:/Qiaolab/LXX/azoospermia/plot/Class3_DEGs_roseplot.pdf",width = 3.8,height = 3.8)
p <- ggplot(rose_long,aes(x=CellType, y=count, fill=Class))+ geom_bar(stat="identity",color="black")+coord_polar()+scale_fill_manual(values = c("red","blue"),limits=c("Up","Down"))
p
dev.off()
###########Class4
file <- read.table("E:/Qiaolab/LXX/azoospermia/data/Nor_Class4_DEG/rose_input.txt",header = T,sep = "\t")
library(reshape2)
rose_long <- melt(file,id.vars="Class",variable.name="CellType",value.name="count")
library(ggplot2)
pdf("E:/Qiaolab/LXX/azoospermia/plot/Class4_DEGs_roseplot.pdf",width = 3.8,height = 3.8)
p <- ggplot(rose_long,aes(x=CellType, y=count, fill=Class))+ geom_bar(stat="identity",color="black")+coord_polar()+scale_fill_manual(values = c("red","blue"),limits=c("Up","Down"))
p
dev.off()
