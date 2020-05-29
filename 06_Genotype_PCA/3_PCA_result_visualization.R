inputpath="/input/path/for/PCA/result"
inputvec="PCA_result_prefix.eigenvec"
inputval="PCA_result_prefix.eigenval"
label="unpruned"
outputpath=inputpath

setwd(inputpath)
pca.x=read.table(inputvec,sep=" ")
pca.val=read.table(inputval,sep="\t")
rownames(pca.x)=pca.x[,1]
pca.x=pca.x[,-c(1,2)]
colnames(pca.x)=paste("PC",seq(1,dim(pca.val)[1]),sep="")

#get race information#
#population information
setwd("/input/path/to/GTEx/phenotype/result")
phenotypename="phs000424.v7.pht002742.v7.p2.c1.GTEx_Subject_Phenotypes.GRU.txt"
phenotype=read.table(phenotypename,sep="\t",skip=10,header=T,fill=TRUE,quote="")
setwd("/input/path/to/GTEx/brain/tissue/sample/annotation")
annotation=read.csv("gtex_v7_brain.csv",header=T)
race=rep(NA,length(rownames(pca.x)))
for (i in 1:length(rownames(pca.x))){
  sampleID=rownames(pca.x)[i]
  race[i]=phenotype[which(phenotype[,"SUBJID"]==as.character(sampleID)),"RACE"]
}
race[which(race==1)]="Asian"
race[which(race==2)]="African American"
race[which(race==3)]="White"
race[which(race==4)]="American Indian"
race[which(race==98)]="Not reported"
race[which(race==99)]="Unknown"


#plot PCA#
library(ggsci)
library(scales)
library(ggplot2)
theme_ydz_scree <- theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
                     plot.title = element_text(size = rel(1.3), vjust = 2, hjust = 0.5, lineheight = 0.8),
                     
                     # axis
                     #axis.title = element_text(size = rel(1.2)),
                     #axis.title = element_blank(),
                     axis.text = element_text(size = rel(1.1)),
                     axis.text.x = element_text(hjust = 0.5,angle = 90, vjust = 0),
                     axis.text.y = element_text(vjust = 0.5, hjust = 0),
                     # axis.line = element_line(colour = "black"),
                     
                     # ticks
                     axis.ticks = element_line(colour = 'darkgrey'),
                     
                     # legend
                     legend.title=element_blank(),
                     legend.text = element_text(size = rel(1.1)),
                     legend.position = "right",
                     legend.direction = 'vertical',
                     legend.background = element_rect(colour = NA, size = 0.5, fill = 'white'),
                     
                     # strip
                     # strip.background=element_rect(fill=NA, color=NA),
                     strip.text=element_text(size = rel(1.3)),
                     
                     # panel
                     panel.background = element_rect(fill = NA, colour = "black", size = 0.7),
                     panel.border = element_rect(colour = "grey", fill=NA, size=1),
                     # panel.grid.major = element_line(colour = "grey", size = 0.05),
                     
                     aspect.ratio=1,
                     complete = T)

theme_ydz_PCA <- theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
                           plot.title = element_text(size = rel(1.3), vjust = 2, hjust = 0.5, lineheight = 0.8),
                           
                           # axis
                           #axis.title = element_text(size = rel(1.2)),
                           #axis.title = element_blank(),
                           axis.text = element_text(size = rel(1.1)),
                           axis.text.x = element_text(hjust = 0.5, vjust = 0),
                           axis.text.y = element_text(vjust = 0.5, hjust = 0),
                           # axis.line = element_line(colour = "black"),
                           
                           # ticks
                           axis.ticks = element_line(colour = 'darkgrey'),
                           
                           # legend
                           legend.title=element_blank(),
                           legend.text = element_text(size = rel(1.1)),
                           legend.position = "right",
                           legend.direction = 'vertical',
                           legend.background = element_rect(colour = NA, size = 0.5, fill = 'white'),
                           
                           # strip
                           # strip.background=element_rect(fill=NA, color=NA),
                           strip.text=element_text(size = rel(1.3)),
                           
                           # panel
                           panel.background = element_rect(fill = NA, colour = "black", size = 0.7),
                           panel.border = element_rect(colour = "grey", fill=NA, size=1),
                           # panel.grid.major = element_line(colour = "grey", size = 0.05),
                           
                           aspect.ratio=1,
                           complete = T)

setwd(outputpath)
pdf(paste("genotype_PCA_plot_",label,".pdf",sep=""))
#scree plot#
eigenvalue=as.numeric(as.matrix(pca.val))
POEV=(abs(eigenvalue)/sum(abs(eigenvalue)))*100   #percentage of explained variance (see https://stats.stackexchange.com/questions/31908/what-is-percentage-of-variance-in-pca)
cutoff=20   #we only print the first 20
df <- data.frame(Percentage.of.explained.variances=POEV[1:cutoff],      
                 PC=paste("PC",seq(1,cutoff),sep=""))
#order the PCs from 1 to cutoff instead of alphabetical order
df$PC <- as.character(df$PC)
#Then turn it back into an ordered factor
df$PC <- factor(df$PC, levels=unique(df$PC))
p<-ggplot(data=df, aes(x=PC, y=Percentage.of.explained.variances)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_ywang_scree
print(p)

#PC plot#
pc1="PC1"
pc2="PC2"
COLOR="race"    #this will control the legend name
temp=as.data.frame(pca.x)
race=as.factor(race)
ggplot(temp,aes_string(x=pc2, y=pc1, color=COLOR)) + geom_point(size=2) + scale_color_npg() + theme_ydz_PCA + 
  labs(x=paste(pc2," (",round(df[which(df[,"PC"]==pc2),"Percentage.of.explained.variances"],2),"%)",sep=""),
       y=paste(pc1," (",round(df[which(df[,"PC"]==pc1),"Percentage.of.explained.variances"],2),"%)",sep=""))

pc1="PC1"
pc2="PC3"
COLOR="race"    #this will control the legend name
temp=as.data.frame(pca.x)
race=as.factor(race)
ggplot(temp,aes_string(x=pc2, y=pc1, color=COLOR)) + geom_point(size=2) + scale_color_npg() + theme_ydz_PCA + 
  labs(x=paste(pc2," (",round(df[which(df[,"PC"]==pc2),"Percentage.of.explained.variances"],2),"%)",sep=""),
       y=paste(pc1," (",round(df[which(df[,"PC"]==pc1),"Percentage.of.explained.variances"],2),"%)",sep=""))

dev.off()




