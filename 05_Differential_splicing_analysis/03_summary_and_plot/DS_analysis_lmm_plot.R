#The purpose of this code:
#summarize the number of age/brainregion/gender dependent exons for SE/A5SS/A3SS and make bar plot of it

library(ggplot2)
library(ggsci)
library(reshape)

splicetypelist=c("SE","A5SS","A3SS")
counttype="JC"
rootinput="./05_Differential_splicing_analysis/02_linear_mixed_model"
rootoutput="./05_Differential_splicing_analysis/03_summary_and_plot/example_output"
fdrcutoff=0.05
library(NMF)
library(RColorBrewer)

signum=matrix(NA,length(splicetypelist),3)
rownames(signum)=splicetypelist
colnames(signum)=c("Brain.region","Age","Gender")

for (splicetype in splicetypelist){
  #read in the DS analysis result
  inputpath=paste(rootinput,splicetype,counttype,sep="/")
  label=strsplit(rootinput,split="/")[[1]][12]
  pvfilename=paste("pvmatrix_",label,".txt",sep="")
  setwd(inputpath)
  pvmatrix=read.table(pvfilename,sep="\t",header=T)
  
  #calculate FDR
  temp=as.numeric(as.matrix(pvmatrix))
  temp[which(is.na(temp))]=1           #assign 1 to missing p values (not converged)
  fdrmatrix=matrix(temp,dim(pvmatrix)[1],dim(pvmatrix)[2])
  rownames(fdrmatrix)=rownames(pvmatrix)
  colnames(fdrmatrix)=colnames(pvmatrix)
  for (i in 1:dim(fdrmatrix)[2]){
    fdrmatrix[,i]=p.adjust(fdrmatrix[,i],"BH")
  }
  
  #get events to plot
  fdrDSbr=which(fdrmatrix[,"pBR"]<fdrcutoff)    #events pass FDR cutoff
  fdrDSage=which(fdrmatrix[,"pAGE"]<fdrcutoff) 
  fdrDSgender=which(fdrmatrix[,"pGENDER"]<fdrcutoff) 
  
  signum[splicetype,"Brain.region"]=length(fdrDSbr)
  signum[splicetype,"Age"]=length(fdrDSage)
  signum[splicetype,"Gender"]=length(fdrDSgender)
}

#generate bar plot
data2plot=melt(signum)
colnames(data2plot)=c("Splice.type","Phenotype","Value")
data2plot$Splice.type=factor(data2plot$Splice.type, levels = c("SE","A5SS","A3SS")) 
data2plot$Phenotype=factor(data2plot$Phenotype, levels = c("Brain.region","Age","Gender")) 
setwd(rootoutput)
p=ggplot(data=data2plot, aes(x=Splice.type, y=Value, fill=Phenotype)) +
  geom_bar(stat="identity",position=position_dodge()) + scale_fill_brewer(palette="Set2") + 
  theme(
    # axis
    axis.title.x = element_text(face="bold", size=16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(face="bold", size=12, margin = margin(t = 0, r = 10, b = 0, l = 0), angle=90),
    axis.text = element_text(size = rel(1.1)),
    axis.text.x = element_text(hjust = 0.5, vjust = 0, size=12, face="bold"),
    axis.text.y = element_text(vjust = 0.5, hjust = 0, size=12, face="bold"),
    axis.line = element_line(colour = "black"),
    #background
    panel.background = element_blank(),
    #legend
    legend.title = element_text(face = "bold"),
    legend.key = element_rect(fill = "white", colour = "black"),
    legend.text = element_text(size = 14, face = "bold"),
    # strip
    strip.text=element_text(size = rel(1.3)),
    aspect.ratio=1,
    complete = T)
pdf("Number_of_significant_event_all_splicetype_all_phenotype.pdf",height=4,width=6)
print(p)
dev.off()



