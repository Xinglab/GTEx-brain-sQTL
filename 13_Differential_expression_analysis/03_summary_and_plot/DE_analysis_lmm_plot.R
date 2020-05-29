expressiontype="allgene"

label="lmm_output_no_interaction_no_warning_one_full_model_original_age_oneforall_original_region_without_SVA"
rootinput="/output/path/of/lmm/result"
rootoutput="/output/path"
fdrcutoff=0.05
library(NMF)
library(RColorBrewer)
library(ggplot2)
library(reshape)

signum=matrix(NA,1,3)
rownames(signum)=expressiontype
colnames(signum)=c("Brain.region","Age","Gender")

#for (expressiontype in expressiontypelist){
setwd(rootinput)
pvfilename=paste("pvmatrix_",label,".txt",sep="")
pvmatrix=read.table(pvfilename,sep="\t",header=T)
rownames(pvmatrix)=rownames(edata)

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
fdrDEbr=which(fdrmatrix[,"pBR"]<fdrcutoff)    #events pass FDR cutoff
fdrDSage=which(fdrmatrix[,"pAGE"]<fdrcutoff) 
fdrDSgender=which(fdrmatrix[,"pGENDER"]<fdrcutoff) 

signum[expressiontype,"Brain.region"]=length(fdrDEbr)
signum[expressiontype,"Age"]=length(fdrDSage)
signum[expressiontype,"Gender"]=length(fdrDSgender)

#generate bar plot
data2plot=melt(signum)
colnames(data2plot)=c("Expression.type","Phenotype","Value")
data2plot$Phenotype=factor(data2plot$Phenotype, levels = c("Brain.region","Age","Gender")) 
setwd(rootoutput)
p=ggplot(data=data2plot, aes(x=Expression.type, y=Value, fill=Phenotype)) +
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
pdf("Number_of_significant_event_all_gene_all_phenotype.pdf",height=4,width=6)
print(p)
dev.off()


