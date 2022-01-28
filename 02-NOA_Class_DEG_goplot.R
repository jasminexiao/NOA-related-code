###########CLass1 UP
setwd("E:/Qiaolab/LXX/azoospermia/data/Nor_Class1_DEG")
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(reshape2)
data.path="E:/Qiaolab/LXX/azoospermia/data/Nor_Class1_DEG/Class1_upDEGs_GO"

#tissue.list=c("BAT","WAT","Liver", "Kidney", "Aorta", "Skin", "BM")
#class.list="Fig7"
go.df=read.csv(paste(data.path,"\\Enrichment_GO\\GO_AllLists.csv",sep = "") )
go.df$GeneList <- paste("X_LogP",go.df$GeneList,sep = "_")
go.df$GeneList <- gsub(" ",".",go.df$GeneList)
heat.df=read.csv(paste(data.path,"\\Enrichment_heatmap\\HeatmapSelectedGOTop100.csv",sep = "") )
go.list=c("R-HSA-2262752","R-HSA-3371497","GO:1903311","R-HSA-1640170","GO:0033044")
heat.df=heat.df[heat.df$GO %in% go.list,]
heat.df=heat.df[,-c(1)]##,10
heat.df$count=rowSums(as.matrix(heat.df[,2:12]) != 0)
heat.df$sumP=abs(rowSums(as.matrix(heat.df[,2:12]))) 
order.list=as.character(heat.df[order(heat.df$count,heat.df$sumP),'Description']) 
heat.df=heat.df[,-c(13,14)]
heat.df= melt(heat.df,id.vars = c('Description'))
# i=1
 for (i in 1:nrow(heat.df)) {
   if(heat.df[i,]$value != 0){
     
     
     des=as.vector(heat.df[i,]$Description) 
     gro=as.character(heat.df[i,]$variable)
     #gro=substr(gro,8,nchar(gro)) 
     
     count=go.df[go.df$Description==des & go.df$GeneList== gro, ]$X.GeneInGOAndHitList
     
   }else{
     count=0
   }
   if(i==1){
     all.count=count
   }else{
     all.count=c(all.count,count)
   }
   
 }

library(SpectralTAD)
heat.df$count=all.count
heat.df$name=heat.df$variable
heat.df$Description=as.character(heat.df$Description)
#order.list[5]="IGF"
#heat.df[grep("IGF",heat.df$Description),]$Description="IGF"
#order.list[9]="cellular response to growth factor stimulusX"
#heat.df[grep("cellular response to growth factor stimulus",heat.df$Description),]$Description="cellular response to growth factor stimulusX"

heat.df$logP=abs(heat.df$value) 
heat.df$Description=factor(heat.df$Description,levels = order.list)
heat.df$name=sapply( strsplit(as.character(heat.df$name), "_"), "[[", 3 )

plot.df=heat.df
write.table(plot.df,"Class1_Select_go.xls",quote = F,col.names = T,row.names = F,sep = "\t")
# pdf("E:/Qiaolab/LXX/azoospermia/plot/NOA_class1_upgo.pdf",width = 12,height = 4)
# ggplot(plot.df, aes(x=name,y=Description,size=count,colour=logP))+
#   geom_point()+
#   scale_size_continuous(breaks = c(10,50,80),range = c(-1,9),name='Size')+
#   scale_color_gradientn(colors = brewer.pal(9,'Reds') ,breaks=c(0,16,32,48),name='LogP')+  theme_classic()+ylab('Go terms')+xlab('')+##YlOrRd
#   theme(legend.text = element_text(size = 15, colour = "black"),
#         axis.title.x = element_text(size = 15, colour = "black"),
#         axis.title.y = element_text(size = 15, colour = "black"),
#         axis.text.y  = element_text(size = 15,colour = 'black'),
#         axis.text.x = element_text(size = 15,colour = 'black',angle = 45, hjust = 0.5, vjust = 0.5),
#         legend.title = element_text(size = 15),
#         #axis.ticks = element_blank(),
#         legend.position ="right",legend.direction = "vertical")
# dev.off()

###########Class2 UP
setwd("E:/Qiaolab/LXX/azoospermia/data/Nor_Class2_DEG")
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(reshape2)
data.path="E:/Qiaolab/LXX/azoospermia/data/Nor_Class2_DEG/Class2_upDEGs_GO"

#tissue.list=c("BAT","WAT","Liver", "Kidney", "Aorta", "Skin", "BM")
#class.list="Fig7"
go.df=read.csv(paste(data.path,"\\Enrichment_GO\\GO_AllLists.csv",sep = "") )
go.df$GeneList <- paste("X_LogP",go.df$GeneList,sep = "_")
go.df$GeneList <- gsub(" ",".",go.df$GeneList)
heat.df=read.csv(paste(data.path,"\\Enrichment_heatmap\\HeatmapSelectedGOTop100.csv",sep = "") )
go.list=c("GO:0002263","GO:0022618","R-HSA-1640170","GO:0043068","GO:0019058")
heat.df=heat.df[heat.df$GO %in% go.list,]
heat.df=heat.df[,-c(1)]##,10
heat.df$count=rowSums(as.matrix(heat.df[,2:9]) != 0)
heat.df$sumP=abs(rowSums(as.matrix(heat.df[,2:9]))) 
order.list=as.character(heat.df[order(heat.df$count,heat.df$sumP),'Description']) 
heat.df=heat.df[,-c(10,11)]
heat.df= melt(heat.df,id.vars = c('Description'))
# i=1
for (i in 1:nrow(heat.df)) {
  if(heat.df[i,]$value != 0){
    
    
    des=as.vector(heat.df[i,]$Description) 
    gro=as.character(heat.df[i,]$variable)
    #gro=substr(gro,8,nchar(gro)) 
    
    count=go.df[go.df$Description==des & go.df$GeneList== gro, ]$X.GeneInGOAndHitList
    
  }else{
    count=0
  }
  if(i==1){
    all.count=count
  }else{
    all.count=c(all.count,count)
  }
  
}

library(SpectralTAD)
heat.df$count=all.count
heat.df$name=heat.df$variable
heat.df$Description=as.character(heat.df$Description)
#order.list[5]="IGF"
#heat.df[grep("IGF",heat.df$Description),]$Description="IGF"
#order.list[9]="cellular response to growth factor stimulusX"
#heat.df[grep("cellular response to growth factor stimulus",heat.df$Description),]$Description="cellular response to growth factor stimulusX"

heat.df$logP=abs(heat.df$value) 
heat.df$Description=factor(heat.df$Description,levels = order.list)
heat.df$name=sapply( strsplit(as.character(heat.df$name), "_"), "[[", 3 )

plot.df=heat.df
write.table(plot.df,"Class2_Select_go.xls",quote = F,col.names = T,row.names = F,sep = "\t")
summary(plot.df$value)
summary(plot.df$count)
# pdf("E:/Qiaolab/LXX/azoospermia/plot/NOA_Class2_upgo.pdf",width = 12,height = 4)
# ggplot(plot.df, aes(x=name,y=Description,size=count,colour=logP))+
#   geom_point()+
#   scale_size_continuous(breaks = c(10,30,50,70,90),range = c(-1,9),name='Size')+
#   scale_color_gradientn(colors = brewer.pal(9,'Reds') ,breaks=c(0,8,16,24,32,48),name='LogP')+  theme_classic()+ylab('Go terms')+xlab('')+##YlOrRd
#   theme(legend.text = element_text(size = 15, colour = "black"),
#         axis.title.x = element_text(size = 15, colour = "black"),
#         axis.title.y = element_text(size = 15, colour = "black"),
#         axis.text.y  = element_text(size = 15,colour = 'black'),
#         axis.text.x = element_text(size = 15,colour = 'black',angle = 45, hjust = 0.5, vjust = 0.5),
#         legend.title = element_text(size = 15),
#         #axis.ticks = element_blank(),
#         legend.position ="right",legend.direction = "vertical")
# dev.off()

###########Class3 UP
setwd("E:/Qiaolab/LXX/azoospermia/data/Nor_Class3_DEG")
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(reshape2)
data.path="E:/Qiaolab/LXX/azoospermia/data/Nor_Class3_DEG/Class3_upDEGs_GO"

#tissue.list=c("BAT","WAT","Liver", "Kidney", "Aorta", "Skin", "BM")
#class.list="Fig7"
go.df=read.csv(paste(data.path,"\\Enrichment_GO\\GO_AllLists.csv",sep = "") )
go.df$GeneList <- paste("X_LogP",go.df$GeneList,sep = "_")
go.df$GeneList <- gsub(" ",".",go.df$GeneList)
heat.df=read.csv(paste(data.path,"\\Enrichment_heatmap\\HeatmapSelectedGOTop100.csv",sep = "") )
go.list=c("GO:0007017","GO:0016569","GO:0010564","GO:0033044","GO:0006397")
heat.df=heat.df[heat.df$GO %in% go.list,]
heat.df=heat.df[,-c(1)]##,10
heat.df$count=rowSums(as.matrix(heat.df[,2:7]) != 0)
heat.df$sumP=abs(rowSums(as.matrix(heat.df[,2:7]))) 
order.list=as.character(heat.df[order(heat.df$count,heat.df$sumP),'Description']) 
heat.df=heat.df[,-c(8,9)]
heat.df= melt(heat.df,id.vars = c('Description'))
# i=1
for (i in 1:nrow(heat.df)) {
  if(heat.df[i,]$value != 0){
    
    
    des=as.vector(heat.df[i,]$Description) 
    gro=as.character(heat.df[i,]$variable)
    #gro=substr(gro,8,nchar(gro)) 
    
    count=go.df[go.df$Description==des & go.df$GeneList== gro, ]$X.GeneInGOAndHitList
    
  }else{
    count=0
  }
  if(i==1){
    all.count=count
  }else{
    all.count=c(all.count,count)
  }
  
}

library(SpectralTAD)
heat.df$count=all.count
heat.df$name=heat.df$variable
heat.df$Description=as.character(heat.df$Description)
#order.list[5]="IGF"
#heat.df[grep("IGF",heat.df$Description),]$Description="IGF"
#order.list[9]="cellular response to growth factor stimulusX"
#heat.df[grep("cellular response to growth factor stimulus",heat.df$Description),]$Description="cellular response to growth factor stimulusX"

heat.df$logP=abs(heat.df$value) 
heat.df$Description=factor(heat.df$Description,levels = order.list)
heat.df$name=sapply( strsplit(as.character(heat.df$name), "_"), "[[", 3 )

plot.df=heat.df
write.table(plot.df,"Class3_Select_go.xls",quote = F,col.names = T,row.names = F,sep = "\t")

###########Class4 UP
setwd("E:/Qiaolab/LXX/azoospermia/data/Nor_Class4_DEG")
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(reshape2)
data.path="E:/Qiaolab/LXX/azoospermia/data/Nor_Class4_DEG/Class4_upDEGs_GO"

#tissue.list=c("BAT","WAT","Liver", "Kidney", "Aorta", "Skin", "BM")
#class.list="Fig7"
go.df=read.csv(paste(data.path,"\\Enrichment_GO\\GO_AllLists.csv",sep = "") )
go.df$GeneList <- paste("X_LogP",go.df$GeneList,sep = "_")
go.df$GeneList <- gsub(" ",".",go.df$GeneList)
heat.df=read.csv(paste(data.path,"\\Enrichment_heatmap\\HeatmapSelectedGOTop100.csv",sep = "") )
go.list=c("GO:0016569","GO:0007017","GO:0006397","R-HSA-1640170","GO:0006281")
heat.df=heat.df[heat.df$GO %in% go.list,]
heat.df=heat.df[,-c(1)]##,10
heat.df$count=rowSums(as.matrix(heat.df[,2:8]) != 0)
heat.df$sumP=abs(rowSums(as.matrix(heat.df[,2:8]))) 
order.list=as.character(heat.df[order(heat.df$count,heat.df$sumP),'Description']) 
heat.df=heat.df[,-c(9,10)]
heat.df= melt(heat.df,id.vars = c('Description'))
# i=1
for (i in 1:nrow(heat.df)) {
  if(heat.df[i,]$value != 0){
    
    
    des=as.vector(heat.df[i,]$Description) 
    gro=as.character(heat.df[i,]$variable)
    #gro=substr(gro,8,nchar(gro)) 
    
    count=go.df[go.df$Description==des & go.df$GeneList== gro, ]$X.GeneInGOAndHitList
    
  }else{
    count=0
  }
  if(i==1){
    all.count=count
  }else{
    all.count=c(all.count,count)
  }
  
}

library(SpectralTAD)
heat.df$count=all.count
heat.df$name=heat.df$variable
heat.df$Description=as.character(heat.df$Description)
#order.list[5]="IGF"
#heat.df[grep("IGF",heat.df$Description),]$Description="IGF"
#order.list[9]="cellular response to growth factor stimulusX"
#heat.df[grep("cellular response to growth factor stimulus",heat.df$Description),]$Description="cellular response to growth factor stimulusX"

heat.df$logP=abs(heat.df$value) 
heat.df$Description=factor(heat.df$Description,levels = order.list)
heat.df$name=sapply( strsplit(as.character(heat.df$name), "_"), "[[", 3 )

plot.df=heat.df
write.table(plot.df,"Class4_Select_go.xls",quote = F,col.names = T,row.names = F,sep = "\t")

all <- read.table("../NOA_4class_select_UpGO1.txt",header = T,sep = "\t")
summary(all$count)
summary(all$value)
order <- read.table("../NOA_4class_select_UpGO_orderTerm.txt",header = F,sep = "\t")
#all$name <- factor(all$name,levels=rev(c("01.SSC","02.Diff.ing.SPG","03.Diff.ed.SPG","04.L","05.Z","06.P","07.D","08.SPC7","09.S","10.ST","11.MIX")))
all$Description <- factor(all$Description,levels=rev(order$V1))
pdf("E:/Qiaolab/LXX/azoospermia/plot/NOA_4Class_upgo.pdf",width = 12,height = 8)
ggplot(all, aes(x=name,y=Description,size=count,colour=logP))+
  geom_point()+
  scale_size_continuous(breaks = c(10,30,50,70,90,110,130),range = c(-1,9),name='Size')+
  scale_colour_gradientn(colors = brewer.pal(9,'Reds') ,breaks=c(0,5,10,15,20,25,30,35,40,45),name='LogP')+  theme_classic()+ylab('Go terms')+xlab('')+##YlOrRd
  theme(legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        axis.text.y  = element_text(size = 15,colour = 'black'),
        axis.text.x = element_text(size = 15,colour = 'black',angle = 45, hjust = 0.5, vjust = 0.5),
        legend.title = element_text(size = 15),
        #axis.ticks = element_blank(),
        legend.position ="right",legend.direction = "vertical")#+ coord_flip()
dev.off()
##########################################################
########--------------down---------------------###########
#########################################################
###########CLass1 down
setwd("E:/Qiaolab/LXX/azoospermia/data/Nor_Class1_DEG")
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(reshape2)
data.path="E:/Qiaolab/LXX/azoospermia/data/Nor_Class1_DEG/Class1_downDEGs_GO"

#tissue.list=c("BAT","WAT","Liver", "Kidney", "Aorta", "Skin", "BM")
#class.list="Fig7"
go.df=read.csv(paste(data.path,"\\Enrichment_GO\\GO_AllLists.csv",sep = "") )
go.df$GeneList <- paste("X_LogP",go.df$GeneList,sep = "_")
go.df$GeneList <- gsub(" ",".",go.df$GeneList)
heat.df=read.csv(paste(data.path,"\\Enrichment_heatmap\\HeatmapSelectedGOTop100.csv",sep = "") )
go.list=c("GO:0007276","GO:0006119","R-HSA-72766","GO:0003341")
heat.df=heat.df[heat.df$GO %in% go.list,]
heat.df=heat.df[,-c(1)]##,10
heat.df$count=rowSums(as.matrix(heat.df[,2:12]) != 0)
heat.df$sumP=abs(rowSums(as.matrix(heat.df[,2:12]))) 
order.list=as.character(heat.df[order(heat.df$count,heat.df$sumP),'Description']) 
heat.df=heat.df[,-c(13,14)]
heat.df= melt(heat.df,id.vars = c('Description'))
# i=1
for (i in 1:nrow(heat.df)) {
  if(heat.df[i,]$value != 0){
    
    
    des=as.vector(heat.df[i,]$Description) 
    gro=as.character(heat.df[i,]$variable)
    #gro=substr(gro,8,nchar(gro)) 
    
    count=go.df[go.df$Description==des & go.df$GeneList== gro, ]$X.GeneInGOAndHitList
    
  }else{
    count=0
  }
  if(i==1){
    all.count=count
  }else{
    all.count=c(all.count,count)
  }
  
}

library(SpectralTAD)
heat.df$count=all.count
heat.df$name=heat.df$variable
heat.df$Description=as.character(heat.df$Description)
#order.list[5]="IGF"
#heat.df[grep("IGF",heat.df$Description),]$Description="IGF"
#order.list[9]="cellular response to growth factor stimulusX"
#heat.df[grep("cellular response to growth factor stimulus",heat.df$Description),]$Description="cellular response to growth factor stimulusX"

heat.df$logP=abs(heat.df$value) 
heat.df$Description=factor(heat.df$Description,levels = order.list)
heat.df$name=sapply( strsplit(as.character(heat.df$name), "_"), "[[", 3 )

plot.df=heat.df
write.table(plot.df,"Class1_Select_Downgo.xls",quote = F,col.names = T,row.names = F,sep = "\t")
# pdf("E:/Qiaolab/LXX/azoospermia/plot/NOA_class1_downgo.pdf",width = 12,height = 4)
# ggplot(plot.df, aes(x=name,y=Description,size=count,colour=logP))+
#   geom_point()+
#   scale_size_continuous(breaks = c(10,50,80),range = c(-1,9),name='Size')+
#   scale_color_gradientn(colors = brewer.pal(9,'Reds') ,breaks=c(0,16,32,48),name='LogP')+  theme_classic()+ylab('Go terms')+xlab('')+##YlOrRd
#   theme(legend.text = element_text(size = 15, colour = "black"),
#         axis.title.x = element_text(size = 15, colour = "black"),
#         axis.title.y = element_text(size = 15, colour = "black"),
#         axis.text.y  = element_text(size = 15,colour = 'black'),
#         axis.text.x = element_text(size = 15,colour = 'black',angle = 45, hjust = 0.5, vjust = 0.5),
#         legend.title = element_text(size = 15),
#         #axis.ticks = element_blank(),
#         legend.position ="right",legend.direction = "vertical")
# dev.off()

###########Class2 down
setwd("E:/Qiaolab/LXX/azoospermia/data/Nor_Class2_DEG")
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(reshape2)
data.path="E:/Qiaolab/LXX/azoospermia/data/Nor_Class2_DEG/Class2_downDEGs_GO"

#tissue.list=c("BAT","WAT","Liver", "Kidney", "Aorta", "Skin", "BM")
#class.list="Fig7"
go.df=read.csv(paste(data.path,"\\Enrichment_GO\\GO_AllLists.csv",sep = "") )
go.df$GeneList <- paste("X_LogP",go.df$GeneList,sep = "_")
go.df$GeneList <- gsub(" ",".",go.df$GeneList)
heat.df=read.csv(paste(data.path,"\\Enrichment_heatmap\\HeatmapSelectedGOTop100.csv",sep = "") )
go.list=c("GO:0007276","GO:0003006","GO:0009566","GO:0030317","GO:0044782")
heat.df=heat.df[heat.df$GO %in% go.list,]
heat.df=heat.df[,-c(1)]##,10
heat.df$count=rowSums(as.matrix(heat.df[,2:9]) != 0)
heat.df$sumP=abs(rowSums(as.matrix(heat.df[,2:9]))) 
order.list=as.character(heat.df[order(heat.df$count,heat.df$sumP),'Description']) 
heat.df=heat.df[,-c(10,11)]
heat.df= melt(heat.df,id.vars = c('Description'))
# i=1
for (i in 1:nrow(heat.df)) {
  if(heat.df[i,]$value != 0){
    
    
    des=as.vector(heat.df[i,]$Description) 
    gro=as.character(heat.df[i,]$variable)
    #gro=substr(gro,8,nchar(gro)) 
    
    count=go.df[go.df$Description==des & go.df$GeneList== gro, ]$X.GeneInGOAndHitList
    
  }else{
    count=0
  }
  if(i==1){
    all.count=count
  }else{
    all.count=c(all.count,count)
  }
  
}

library(SpectralTAD)
heat.df$count=all.count
heat.df$name=heat.df$variable
heat.df$Description=as.character(heat.df$Description)
#order.list[5]="IGF"
#heat.df[grep("IGF",heat.df$Description),]$Description="IGF"
#order.list[9]="cellular response to growth factor stimulusX"
#heat.df[grep("cellular response to growth factor stimulus",heat.df$Description),]$Description="cellular response to growth factor stimulusX"

heat.df$logP=abs(heat.df$value) 
heat.df$Description=factor(heat.df$Description,levels = order.list)
heat.df$name=sapply( strsplit(as.character(heat.df$name), "_"), "[[", 3 )

plot.df=heat.df
write.table(plot.df,"Class2_Select_Downgo.xls",quote = F,col.names = T,row.names = F,sep = "\t")
summary(plot.df$value)
summary(plot.df$count)
# pdf("E:/Qiaolab/LXX/azoospermia/plot/NOA_Class2_downgo.pdf",width = 12,height = 4)
# ggplot(plot.df, aes(x=name,y=Description,size=count,colour=logP))+
#   geom_point()+
#   scale_size_continuous(breaks = c(10,30,50,70,90),range = c(-1,9),name='Size')+
#   scale_color_gradientn(colors = brewer.pal(9,'Reds') ,breaks=c(0,8,16,24,32,48),name='LogP')+  theme_classic()+ylab('Go terms')+xlab('')+##YlOrRd
#   theme(legend.text = element_text(size = 15, colour = "black"),
#         axis.title.x = element_text(size = 15, colour = "black"),
#         axis.title.y = element_text(size = 15, colour = "black"),
#         axis.text.y  = element_text(size = 15,colour = 'black'),
#         axis.text.x = element_text(size = 15,colour = 'black',angle = 45, hjust = 0.5, vjust = 0.5),
#         legend.title = element_text(size = 15),
#         #axis.ticks = element_blank(),
#         legend.position ="right",legend.direction = "vertical")
# dev.off()

###########Class3 down
setwd("E:/Qiaolab/LXX/azoospermia/data/Nor_Class3_DEG")
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(reshape2)
data.path="E:/Qiaolab/LXX/azoospermia/data/Nor_Class3_DEG/Class3_downDEGs_GO"

#tissue.list=c("BAT","WAT","Liver", "Kidney", "Aorta", "Skin", "BM")
#class.list="Fig7"
go.df=read.csv(paste(data.path,"\\Enrichment_GO\\GO_AllLists.csv",sep = "") )
go.df$GeneList <- paste("X_LogP",go.df$GeneList,sep = "_")
go.df$GeneList <- gsub(" ",".",go.df$GeneList)
heat.df=read.csv(paste(data.path,"\\Enrichment_heatmap\\HeatmapSelectedGOTop100.csv",sep = "") )
go.list=c("R-HSA-72766","hsa00190","GO:0022613","GO:0009060","GO:0048232")
heat.df=heat.df[heat.df$GO %in% go.list,]
heat.df=heat.df[,-c(1)]##,10
heat.df$count=rowSums(as.matrix(heat.df[,2:7]) != 0)
heat.df$sumP=abs(rowSums(as.matrix(heat.df[,2:7]))) 
order.list=as.character(heat.df[order(heat.df$count,heat.df$sumP),'Description']) 
heat.df=heat.df[,-c(8,9)]
heat.df= melt(heat.df,id.vars = c('Description'))
# i=1
for (i in 1:nrow(heat.df)) {
  if(heat.df[i,]$value != 0){
    
    
    des=as.vector(heat.df[i,]$Description) 
    gro=as.character(heat.df[i,]$variable)
    #gro=substr(gro,8,nchar(gro)) 
    
    count=go.df[go.df$Description==des & go.df$GeneList== gro, ]$X.GeneInGOAndHitList
    
  }else{
    count=0
  }
  if(i==1){
    all.count=count
  }else{
    all.count=c(all.count,count)
  }
  
}

library(SpectralTAD)
heat.df$count=all.count
heat.df$name=heat.df$variable
heat.df$Description=as.character(heat.df$Description)
#order.list[5]="IGF"
#heat.df[grep("IGF",heat.df$Description),]$Description="IGF"
#order.list[9]="cellular response to growth factor stimulusX"
#heat.df[grep("cellular response to growth factor stimulus",heat.df$Description),]$Description="cellular response to growth factor stimulusX"

heat.df$logP=abs(heat.df$value) 
heat.df$Description=factor(heat.df$Description,levels = order.list)
heat.df$name=sapply( strsplit(as.character(heat.df$name), "_"), "[[", 3 )

plot.df=heat.df
write.table(plot.df,"Class3_Select_Downgo.xls",quote = F,col.names = T,row.names = F,sep = "\t")

###########Class4 down
setwd("E:/Qiaolab/LXX/azoospermia/data/Nor_Class4_DEG")
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(reshape2)
data.path="E:/Qiaolab/LXX/azoospermia/data/Nor_Class4_DEG/Class4_downDEGs_GO"

#tissue.list=c("BAT","WAT","Liver", "Kidney", "Aorta", "Skin", "BM")
#class.list="Fig7"
go.df=read.csv(paste(data.path,"\\Enrichment_GO\\GO_AllLists.csv",sep = "") )
go.df$GeneList <- paste("X_LogP",go.df$GeneList,sep = "_")
go.df$GeneList <- gsub(" ",".",go.df$GeneList)
heat.df=read.csv(paste(data.path,"\\Enrichment_heatmap\\HeatmapSelectedGOTop100.csv",sep = "") )
go.list=c("R-HSA-72766","GO:0007005","GO:0022613","GO:0048232","GO:0009060")
heat.df=heat.df[heat.df$GO %in% go.list,]
heat.df=heat.df[,-c(1)]##,10
heat.df$count=rowSums(as.matrix(heat.df[,2:8]) != 0)
heat.df$sumP=abs(rowSums(as.matrix(heat.df[,2:8]))) 
order.list=as.character(heat.df[order(heat.df$count,heat.df$sumP),'Description']) 
heat.df=heat.df[,-c(9,10)]
heat.df= melt(heat.df,id.vars = c('Description'))
# i=1
for (i in 1:nrow(heat.df)) {
  if(heat.df[i,]$value != 0){
    
    
    des=as.vector(heat.df[i,]$Description) 
    gro=as.character(heat.df[i,]$variable)
    #gro=substr(gro,8,nchar(gro)) 
    
    count=go.df[go.df$Description==des & go.df$GeneList== gro, ]$X.GeneInGOAndHitList
    
  }else{
    count=0
  }
  if(i==1){
    all.count=count
  }else{
    all.count=c(all.count,count)
  }
  
}

library(SpectralTAD)
heat.df$count=all.count
heat.df$name=heat.df$variable
heat.df$Description=as.character(heat.df$Description)
#order.list[5]="IGF"
#heat.df[grep("IGF",heat.df$Description),]$Description="IGF"
#order.list[9]="cellular response to growth factor stimulusX"
#heat.df[grep("cellular response to growth factor stimulus",heat.df$Description),]$Description="cellular response to growth factor stimulusX"

heat.df$logP=abs(heat.df$value) 
heat.df$Description=factor(heat.df$Description,levels = order.list)
heat.df$name=sapply( strsplit(as.character(heat.df$name), "_"), "[[", 3 )

plot.df=heat.df
write.table(plot.df,"Class4_Select_Downgo.xls",quote = F,col.names = T,row.names = F,sep = "\t")

all <- read.table("../NOA_4class_select_downGO.txt",header = T,sep = "\t")
summary(all$count)
summary(all$value)
order <- read.table("../NOA_4class_select_downGO_orderTerm.txt",header = F,sep = "\t")
#all$name <- factor(all$name,levels=rev(c("01.SSC","02.Diff.ing.SPG","03.Diff.ed.SPG","04.L","05.Z","06.P","07.D","08.SPC7","09.S","10.ST","11.MIX")))
all$Description <- factor(all$Description,levels=rev(order$V1))
pdf("E:/Qiaolab/LXX/azoospermia/plot/NOA_4Class_downgo.pdf",width = 10.8,height = 8)
ggplot(all, aes(x=name,y=Description,size=count,colour=logP))+
  geom_point()+
  scale_size_continuous(breaks = c(10,30,50,70,90,110,130),range = c(-1,9),name='Size')+
  scale_colour_gradientn(colors = brewer.pal(9,'GnBu') ,breaks=c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85),name='LogP')+  theme_classic()+ylab('Go terms')+xlab('')+##YlOrRd
  theme(legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        axis.text.y  = element_text(size = 15,colour = 'black'),
        axis.text.x = element_text(size = 15,colour = 'black',angle = 45, hjust = 0.5, vjust = 0.5),
        legend.title = element_text(size = 15),
        #axis.ticks = element_blank(),
        legend.position ="right",legend.direction = "vertical")#+ coord_flip()
dev.off()
