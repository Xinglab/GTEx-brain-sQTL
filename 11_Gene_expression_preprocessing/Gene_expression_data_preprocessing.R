inputpath="/input/to/GTEx/processed/gene/expression/data"
filename="GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct"          #gene TPM

##################
#read in raw data#
##################
library(data.table)
setwd(inputpath)
rawEXPR=fread(filename,header=T,data.table=FALSE)

#########################################
#get the subset that is related to brain#
#########################################
#get SRR ID of all samples
setwd("/path/to/PSI/value/from/rMATS")
data=read.table("A3SS_JC_IC_filter.txt",sep="\t",header=T)
SRRlist=colnames(data)

#get ID conversion information
annotationpath="/input/path/to/GTEx/brain/tissue/sample/annotation"
IDconversionpath=annotationpath
IDconversionname="gtex_v7_brain.csv" 
setwd(annotationpath)
IDconversion=read.csv(IDconversionname,header=T)

sampleID=rep(NA,length(SRRlist))
for (i in 1:length(SRRlist)){
  temp=as.character(IDconversion[which(IDconversion[,"Run_s"]==SRRlist[i]),"Sample_Name_s"])[1]
  sampleID[i]=temp
}

commonsamples=intersect(colnames(rawEXPR),sampleID)

subrawEXPR=rawEXPR[,which(colnames(rawEXPR) %in% c("Name","Description",commonsamples))]
setwd(inputpath)
write.table(subrawEXPR,"GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_brain.txt",sep="\t")

IDtable=cbind(SRRlist,sampleID)
write.table(IDtable,"SRRID_2_GTEXID.txt",sep="\t")


###############
#normalization#
###############
#normalization#
medianquartile=function(x){
  percentage=0.5   #this is the median
  #remove 0 from x
  temp=x[x!=0]
  #sort the values in increasing order
  temp_sort=sort(temp)
  #find the 50th percentile value (median)
  normalize_factor=quantile(temp_sort,percentage)      #quantile function will sort the list so we don't actually need to sort it
  #Divide all the expression values by the upper quartile value
  return(x/normalize_factor)
}

normTPM=apply(subrawEXPR[,-c(1,2)],2,medianquartile)
normTPM=cbind(subrawEXPR[,c(1,2)],normTPM)
setwd(inputpath)
write.table(normTPM,"GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_brain_normalized.txt",sep="\t")

##############################################
#get normalized expression of all RBPs/SF/TFs#
##############################################
#read in RBP list
setwd("/path/to/list/of/RBPs")
RBPlist=read.table("summarized_RBP_list.txt",sep="\t",header=T)
expr_RBP_norm=normTPM[which(as.character(normTPM[,"Description"]) %in% as.character(RBPlist[,"gene.name"])),]   #there can be duplications (same gene appears in multiple rows)
#there are 4 RBPs: "VARSL"  "DND1"   "LUC7L2" "NIFK" that are not in the expression data (normTPM)
setwd(inputpath)
write.table(expr_RBP_norm,"all_RBP_normalized_expression.txt",sep="\t")

#read in SF list
setwd("/path/to/list/of/SFs")
SFlist=read.table("summarized_SF_list.txt",sep="\t",header=T)
SFlist=subset(SFlist,SFlist[,"type"]=="RBP/SF")
expr_SF_norm=normTPM[which(as.character(normTPM[,"Description"]) %in% as.character(SFlist[,"HGNC.symbol"])),]    #there can be duplications (same gene appears in multiple rows)
setwd(inputpath)
write.table(expr_SF_norm,"SF_normalized_expression.txt",sep="\t")

