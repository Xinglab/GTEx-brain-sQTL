library(gplots)
library(ggplot2)
library(reshape)
library(RColorBrewer)
require(scales)

outputpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.3_motif_analysis/individual_exon_binding_peaks_scan/DeepBind/gene_expression_per_region"
source("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.3_motif_analysis/individual_exon_binding_peaks_scan/DeepBind/heatmap.3.R")

vectorsplit=function(x,del,pos){
  return(strsplit(x,split=del)[[1]][pos])
}

formalbrainregionlist=c("Amygdala",
                        "Anterior cingulate cortex BA24",
                        "Caudate basal ganglia",
                        "Cerebellar Hemisphere",
                        "Cerebellum",
                        "Cortex",
                        "Frontal Cortex BA9",
                        "Hippocampus",
                        "Hypothalamus",
                        "Nucleus accumbens basal ganglia",
                        "Putamen basal ganglia",
                        "Spinal cord cervical c-1",
                        "Substantia nigra")

#read in the mean TPM in each brain region for all genes (this result is the same with mean_normalized_RBP_expression_per_region.txt, just with more information)
setwd("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/8.5_DE_analysis_oneforall_result_summary_and_plot/result/lmm_output_no_interaction_no_warning_one_full_model_original_age_oneforall_original_region_without_SVA/allgene/fdr=0.01~foldchange=10~deltaTPM=1")
allgenemeanTPM=as.matrix(read.table("DE_BR_allgene_lmm_output_no_interaction_no_warning_one_full_model_original_age_oneforall_original_region_without_SVA.txt",sep="\t",header=T,check.names=F))
rowname=sapply(rownames(allgenemeanTPM),vectorsplit,del="_",pos=2)
rownames(allgenemeanTPM)=rowname
allgenemeanTPM=allgenemeanTPM[,1:13]
colnames(allgenemeanTPM)=formalbrainregionlist


setwd("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/data/RBP_list/RBP_from_NRG")
RBPlist=read.table("summarized_RBP_list.txt",sep="\t",header=T)
RBPlist=as.character(RBPlist[,"gene.name"])
setwd("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/data/splicing_factor_list/SF_from_GR")
SFlist=read.table("summarized_SF_list.txt",sep="\t",header=T)
SFlist=subset(SFlist,SFlist[,"type"]=="RBP/SF")
SFlist=as.character(SFlist[,"HGNC.symbol"])

allRBP=unique(c(RBPlist,SFlist))

#get the RBP list from DeepBind
RBPdb="/u/nobackup/yxing/PROJECT/yidazhan/research/software/deepbind/db/db.tsv"
RBPtable=read.table(RBPdb,sep="\t",header=T)
subRBPtable=subset(RBPtable,RBPtable[,"Species"]=="Homo sapiens")
subRBPtable=subset(subRBPtable,subRBPtable[,"Type"]=="RBP")

#change RBP name (the name of some RBPs in the DeepBind table is different with their name in the gene expression table. They are the same RBP but use different names)
subRBPtable=as.matrix(subRBPtable)
subRBPtable[which(subRBPtable[,"Protein"]=="BRUNOL4"),"Protein"]="CELF4"
subRBPtable[which(subRBPtable[,"Protein"]=="HuR"),"Protein"]="ELAVL1"
subRBPtable[which(subRBPtable[,"Protein"]=="STAR-PAP"),"Protein"]="TUT1"
subRBPtable[which(subRBPtable[,"Protein"]=="BRUNOL5"),"Protein"]="CELF5"
subRBPtable[which(subRBPtable[,"Protein"]=="BRUNOL6"),"Protein"]="CELF6"
subRBPtable[which(subRBPtable[,"Protein"]=="HNRPLL"),"Protein"]="HNRNPLL"

#get the overlap
deepbindRBP=intersect(as.character(subRBPtable[,"Protein"]),union(RBPlist,SFlist))

#add RBFOX2 and RBFOX3 into the plot. Although we don't have Deepbind model for RBFOX2/3, they share very similar motif with RBFOX1. So we still want to include them
#add NOVA1 and NOVA2 into the plot. Although we don't have Deepbind model for NOVA1/2, they have well known motif
deepbindRBP=c(deepbindRBP,"RBFOX2","RBFOX3","NOVA1","NOVA2")
################################################################################

#make the plot 
zscore=function(x){
  return((x-mean(x))/sd(x))     #z-score
}

label1=c("all_RBP","DB_RBP")
label2=c("original_norm_TPM","z_score_norm_TPM")

for (labela in label1){
  for (labelb in label2){
    if (labela=="all_RBP"){
      meanTPM2plot=allgenemeanTPM[intersect(allRBP,rownames(allgenemeanTPM)),]        #"SRSF8"   "AKAP17A" "CPEB3" appears twice in the allgenemeanTPM matrix; VARSL"  "DND1"   "LUC7L2" "NIFK" are not included in the allgenemeanTPM matrix
      plotrowname=NA    #we don't plot row name, i.e., gene symbol since there are too many
    }
    if (labela=="DB_RBP"){
      meanTPM2plot=allgenemeanTPM[deepbindRBP,]
      plotrowname=NULL
    }
    if (labelb=="original_norm_TPM"){
      meanTPM2plot=meanTPM2plot
    }
    if (labelb=="z_score_norm_TPM"){
      meanTPM2plot=t(apply(meanTPM2plot,1,zscore))
    }
    meanTPM2plot[is.na(meanTPM2plot)]=0
    setwd(outputpath)
    pdf(paste("Region_dependent_",labela,"_expression_heatmap_",labelb,".pdf",sep=""),height=10,width=6)
    #ownbreak = c(seq(range(as.numeric(betamatrix))[1],-0.8,length=250),seq(-0.79,0.79,length=501),seq(0.8,range(as.numeric(betamatrix))[2],length=250))
    heatmap.3(meanTPM2plot, 
              col = colorRampPalette(rev(brewer.pal(11,"RdBu")))(1000), 
              key = TRUE,
              Colv = FALSE,
              #breaks = ownbreak,
              labRow = plotrowname)
    dev.off()
  }
}



