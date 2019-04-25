#The purpose of this code:
#we group the samples into two groups: samples from regions where the sQTL is significant and samples from regions where the sQTL is insignificant
#then we get the expression of those samples in the two groups and perform a differential expression analysis
#also we combine the calculation and plot step into one code

job <- as.numeric(Sys.getenv("SGE_TASK_ID"))
args <- commandArgs(TRUE)
splicetype=args[1]
type=args[2]
chunksize=as.numeric(args[3])

PSItype="original"       #change this to "original" when just making the detailed plot. For other plot, keep it "logit"
counttype="JC"

library(robust)
library(gplots)
library(ggplot2)
library(reshape)
library(RColorBrewer)
require(scales)
library(corrplot)
library(ggsci)

library("data.table")
library(gaston)
library(boot)

'%!in%' <- function(x,y)!('%in%'(x,y))
source("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.3_motif_analysis/individual_exon_binding_peaks_scan/DeepBind/heatmap.3.R")

vectorsplit=function(x,del,pos){
  return(strsplit(x,split=del)[[1]][pos])
}

outputpath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.3_motif_analysis/individual_exon_binding_peaks_scan/DeepBind/deepbind_region_dependent_sQTL_RBP_plot_differential_expression_version",
                 splicetype,type,"duplicated_RBP_version",sep="/")

sqtlrootinput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/sQTL_run/logit/JC",
                    splicetype,sep="/")

brainregionlist=c("Brain-Amygdala",
                  "Brain-AnteriorcingulatecortexBA24",
                  "Brain-Caudatebasalganglia",
                  "Brain-CerebellarHemisphere",
                  "Brain-Cerebellum",
                  "Brain-Cortex",
                  "Brain-FrontalCortexBA9",
                  "Brain-Hippocampus",
                  "Brain-Hypothalamus",
                  "Brain-Nucleusaccumbensbasalganglia",
                  "Brain-Putamenbasalganglia",
                  "Brain-Spinalcordcervicalc-1",
                  "Brain-Substantianigra")

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

#######################
#splicing and genotype#
#######################
summaryinput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/summary/","logit","/",counttype,"/",splicetype,sep="")
rootsqtl="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/sQTL_run"
totalcountinput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/1_get_matrix_from_rMATS/result/",counttype,sep="")
GWASdbpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project/analysis/1_GLMM/sQTL/GWAS_databases"
GWASdbname="gwas_catalog_v1.0.1-associations_e89_r2017-07-31.tsv"
genoplinkprefix="Genotype_V7_plink_maf0.05"
LDpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project/analysis/1_GLMM/sQTL/LD"
sqtlrootinput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/sQTL_run/logit/JC",
                    splicetype,sep="/")
rmatspostpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/data/brain_rMATS_post"

#read in the PSI value
if (PSItype=="original"){   #original PSI value without any correction
  psiinputpath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/1_get_matrix_from_rMATS/result/",counttype,sep="")  
  inputPSI=paste(splicetype,"_",counttype,"_PSI_filter.txt",sep="")
}
if (PSItype=="logit"){     #logit PSI with correction (because we used corrected PSI for all the sQTL calculation, when we plot the result, we also use this but only transform it back to 0-1 scale)
  psiinputpath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/5_correction_for_technical_confounders/result/",counttype,sep="")
  inputPSI=paste(splicetype,"_",counttype,"_logit_corrected_PSI.txt",sep="")
}
inputtotalRC=paste(splicetype,"_",counttype,"_totalRC_filter.txt",sep="")
inputIC=paste(splicetype,"_",counttype,"_IC_filter.txt",sep="")

setwd(psiinputpath)
PSI=read.table(inputPSI,sep="\t",header=T)              #PSI as the normalized inclusion count
setwd(totalcountinput)
totalRC=read.table(inputtotalRC,sep="\t",header=T)
IC=read.table(inputIC,sep="\t",header=T)

#transform logit PSI back to PSI
if (PSItype=="logit"){
  temp=inv.logit(as.numeric(as.matrix(PSI)))
  originalPSI=matrix(temp,dim(PSI)[1],dim(PSI)[2])
  rownames(originalPSI)=rownames(PSI)
  colnames(originalPSI)=colnames(PSI)
  PSI=as.data.frame(originalPSI)
}

#read in exon effective length table
elinput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/1_get_matrix_from_rMATS/result/",counttype,sep="")
setwd(elinput)
effectlengthsubset=try(suppressMessages(read.table(paste("subset_",counttype,".raw.input.",splicetype,".txt",sep=""),sep="\t",header=T)),silent=TRUE) 
if (!(inherits(effectlengthsubset,"try-error"))){       #if we have already have the sutset file (the original file is too big)
  effectlength=effectlengthsubset
}else{      #if we don't, we need to generate that
  #setwd("/u/nobackup/yxing/NOBACKUP/harryyan/gtex_sqtl/gtex_rMATS/post_tissue/brain_rMATS_post")
  setwd(rmatspostpath)
  effectlengthfull=fread(paste(counttype,".raw.input.",splicetype,".txt",sep=""),sep="\t",header=T)
  effectlengthfull=as.data.frame(effectlengthfull)
  
  rownames(effectlengthfull)=effectlengthfull[,"ID"]
  exonID=sapply(rownames(PSI),function(x) return(strsplit(x,split="\\|")[[1]][1]))
  effectlengthsubset=effectlengthfull[exonID,]
  
  setwd(elinput)
  write.table(effectlengthsubset,paste("subset_",counttype,".raw.input.",splicetype,".txt",sep=""),sep="\t")
  effectlength=effectlengthsubset
}
rownames(effectlength)=rownames(PSI)
#read in exon annotation table
setwd(rmatspostpath)
fromGTF=read.table(paste("fromGTF.",splicetype,".txt",sep=""),sep="\t",header=T)

#################
#gene expression#
#################
#read in the joblist
joblistinput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.3_motif_analysis/individual_exon_binding_peaks_scan/DeepBind/deepbind_input/",splicetype,"/",type,sep="")
setwd(joblistinput)
uniquejoblist=read.table(paste(splicetype,"_",type,"_uniquejoblist.txt",sep=""),sep="\t")
#read in the expression table
setwd("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/data/V7_expression_TPM")
normTPM=read.table("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_brain_normalized.txt",sep="\t",header=T,check.names=F)
normTPM=as.matrix(normTPM)

#change RBP name (the name of some RBPs in the DeepBind table is different with their name in the gene expression table. They are the same RBP but use different names)
normTPM[which(normTPM[,2] %in% "CELF4"),2]="BRUNOL4"
normTPM[which(normTPM[,2] %in% "ELAVL1"),2]="HuR"
normTPM[which(normTPM[,2] %in% "TUT1"),2]="STAR-PAP"
normTPM[which(normTPM[,2] %in% "CELF5"),2]="BRUNOL5"
normTPM[which(normTPM[,2] %in% "CELF6"),2]="BRUNOL6"
normTPM[which(normTPM[,2] %in% "HNRNPLL"),2]="HNRPLL"

deepbindRBPlist=as.character(unique(uniquejoblist[,"RBP"]))
#gene ID conversion table
RBPIDconversion=matrix(NA,length(deepbindRBPlist),4)
colnames(RBPIDconversion)=c("DeepBindID","GeneSymbol","EnsemblID","edataID")
RBPIDconversion[,"DeepBindID"]=deepbindRBPlist
RBPIDconversion[,"GeneSymbol"]=sapply(deepbindRBPlist,vectorsplit,del="_",2)
Ensembl2symbol=cbind(as.character(normTPM[,1]),as.character(normTPM[,2]),paste(as.character(normTPM[,1]),as.character(normTPM[,2]),sep="_"))
rownames(Ensembl2symbol)=Ensembl2symbol[,2]
RBPIDconversion[,"EnsemblID"]=Ensembl2symbol[RBPIDconversion[,"GeneSymbol"],1]
RBPIDconversion[,"edataID"]=Ensembl2symbol[RBPIDconversion[,"GeneSymbol"],3]
#######################################################################################################################################################################
#add RBFOX2 and RBFOX3 into the table. Although we don't have Deepbind model for RBFOX2/3, they share very similar motif with RBFOX1. So we still want to include them#
#######################################################################################################################################################################
RBFOX2=c("D00210.001_RBFOX1","RBFOX1","ENSG00000100320.18","ENSG00000100320.18_RBFOX2")
RBFOX3=c("D00210.001_RBFOX1","RBFOX1","ENSG00000167281.14","ENSG00000167281.14_RBFOX3")
RBPIDconversion=rbind(rbind(RBPIDconversion,RBFOX2),RBFOX3)
################################################################################
edata=normTPM
rownames(edata)=paste(as.character(normTPM[,1]),as.character(normTPM[,2]),sep="_")
edata=edata[,-c(1,2)]
temp=matrix(as.numeric(edata),dim(edata)[1],dim(edata)[2])
rownames(temp)=rownames(edata)
colnames(temp)=colnames(edata)
edata=temp

##########################
#brain region information#
##########################
phenotypepath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/2_get_phenotype_information/result"
setwd(phenotypepath)
oribrainregion=read.table("brain_region_table_brain.txt",sep="\t",header=T)
oribrainregion[which(oribrainregion=="Brain - Amygdala")]="Brain-Amygdala"
oribrainregion[which(oribrainregion=="Brain - Anterior cingulate cortex (BA24)")]="Brain-AnteriorcingulatecortexBA24"
oribrainregion[which(oribrainregion=="Brain - Caudate (basal ganglia)")]="Brain-Caudatebasalganglia"
oribrainregion[which(oribrainregion=="Brain - Cerebellar Hemisphere")]="Brain-CerebellarHemisphere"
oribrainregion[which(oribrainregion=="Brain - Cerebellum")]="Brain-Cerebellum"
oribrainregion[which(oribrainregion=="Brain - Cortex")]="Brain-Cortex"
oribrainregion[which(oribrainregion=="Brain - Frontal Cortex (BA9)")]="Brain-FrontalCortexBA9"
oribrainregion[which(oribrainregion=="Brain - Hippocampus")]="Brain-Hippocampus"
oribrainregion[which(oribrainregion=="Brain - Hypothalamus")]="Brain-Hypothalamus"
oribrainregion[which(oribrainregion=="Brain - Nucleus accumbens (basal ganglia)")]="Brain-Nucleusaccumbensbasalganglia"
oribrainregion[which(oribrainregion=="Brain - Putamen (basal ganglia)")]="Brain-Putamenbasalganglia"
oribrainregion[which(oribrainregion=="Brain - Spinal cord (cervical c-1)")]="Brain-Spinalcordcervicalc-1"
oribrainregion[which(oribrainregion=="Brain - Substantia nigra")]="Brain-Substantianigra"
setwd("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/data/V7_expression_TPM")
IDconversion=read.table("SRRID_2_GTEXID.txt",sep="\t",header=T)
rownames(IDconversion)=IDconversion[,"sampleID"]

colnames(edata)=as.character(IDconversion[colnames(edata),"SRRlist"])
edataRBP=as.matrix(edata[RBPIDconversion[,"edataID"],])
rownames(edataRBP)=RBPIDconversion[,"edataID"]
brainregion=oribrainregion[,as.character(colnames(edata))]


#get the full exon information
#get all sQTL exons
rootinput="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/summary/logit/JC"
sQTLexon=c()
inputpath=paste(rootinput,splicetype,sep="/")
setwd(inputpath)
for (i in 1:length(brainregionlist)){
  temp=read.table(paste("all_info_exon_sQTL_GWAS_disease_in_",brainregionlist[i],"_logit_JC_",splicetype,"_",type,".txt",sep=""),sep="\t",header=T)
  sQTLexon=c(sQTLexon,as.character(temp[,"exon_full_ID"]))
}
sQTLexon=unique(sQTLexon)
shortIDlist=rep(NA,length(sQTLexon))
exonsymbol=rep(NA,length(sQTLexon))
for (e in 1:length(sQTLexon)){
  shortIDlist[e]=paste("SE",strsplit(sQTLexon[e],split="\\|")[[1]][1],sep="_")
  exonsymbol[e]=strsplit(sQTLexon[e],split="\\|")[[1]][3]
}
exonIDconversion=cbind(sQTLexon,shortIDlist,exonsymbol)
rownames(exonIDconversion)=shortIDlist
colnames(exonIDconversion)=c("fullID","shortID","gene.symbol")

exonsnppair=unique(uniquejoblist[,c("Exon","SNP")])
#filter out SNPs which are not single base 
allele0=sapply(as.character(exonsnppair[,"SNP"]),vectorsplit,"_",3)
allele1=sapply(as.character(exonsnppair[,"SNP"]),vectorsplit,"_",4)
allele0length=sapply(allele0,nchar)
allele1length=sapply(allele1,nchar)
singlebase=intersect(which(allele0length==1),which(allele1length==1))
singlebaseexonsnppair=exonsnppair[singlebase,]

###################################################################################################################################################
###########
#section 3#
###########
###generate plot of those significant results###
setwd(outputpath)
sigdiffpair=read.table(paste(splicetype,"_",type,"_significant_sQTL_RBP_pairs.txt",sep=""),sep="\t",header=T)

#we make the plot for the following part of the result
start=1+chunksize*(job-1)
end=chunksize*job
if (end>=dim(sigdiffpair)[1]){
  end=dim(sigdiffpair)[1]
}
print(paste(start,end,sep="   -   "))

#outputpath for box plot
boxplotpath=paste(outputpath,"significant_pair_boxplot",sep="/")
command=paste("mkdir -p",boxplotpath)
system(command)

#outputpath for the 13 region plot
regionplotpath=paste(outputpath,"significant_pair_detailed_plot",sep="/")
command=paste("mkdir -p",regionplotpath)
system(command)

#outputpath for the sashimi plot
sashimiplotpath=paste(outputpath,"significant_pair_sashimi_plot",sep="/")
command=paste("mkdir -p",sashimiplotpath)
system(command)


#make box plot for significant pairs
Boxplot=function(SIGREGION,EXPR,sQTLID,RBPID,PV,FC,FDR,path,name){
  sigTPM=as.numeric(as.matrix(EXPR[RBPID,colnames(brainregion)[which(brainregion %in% strsplit(SIGREGION,split=", ")[[1]])]]))
  insigTPM=as.numeric(as.matrix(EXPR[RBPID,setdiff(colnames(brainregion),colnames(brainregion)[which(brainregion %in% strsplit(SIGREGION,split=", ")[[1]])])]))
  exprBR=c(sigTPM,insigTPM)
  label=c(rep("Significant",length(sigTPM)),rep("Insignificant",length(insigTPM)))
  ylim.range=range(exprBR)
  data2plot=data.frame("TPM"=exprBR,"Group"=label)
  
  setwd(path)
  pdf(name,width=4,height=4)
  title=paste(sQTLID,
              RBPID,
              SIGREGION,
              paste("P value:",paste(formatC(PV, format = "e", digits = 3),collapse=",   ")),
              paste("FDR: ",paste(formatC(FDR, format = "e", digits = 3),collapse=",   ")),
              paste("Fold change: ",round(FC, digits = 3),sep=""),
              sep="\n")
  p=ggplot(data2plot, aes(x=Group, y=TPM, color=Group)) + scale_color_brewer(palette="Dark2") +
    geom_boxplot(notch=TRUE) +
    geom_jitter(shape=16, position=position_jitter(0.2)) +
    ggtitle(title) + 
    theme(
      # axis
      axis.title.x = element_text(face="bold", size=12, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(face="bold", size=12, margin = margin(t = 0, r = 10, b = 0, l = 0), angle=90),
      axis.text = element_text(size = rel(1.1)),
      axis.text.x = element_text(hjust = 0.5, vjust = 0, size=14),
      axis.text.y = element_text(vjust = 0.5, hjust = 0, size=14),
      axis.line = element_line(colour = "black"),
      plot.title = element_text(size=8),
      
      #background
      panel.background = element_blank(),
      
      #legend
      legend.position="none",
      
      # strip
      strip.text=element_text(size = rel(1.3)),
      aspect.ratio=1,
      complete = T)
  print(p)
  dev.off()
}

#get events and correponding information to plot
getexon=function(name){
  exoninfo=matrix(NA,length(name),11)
  rownames(exoninfo)=name
  colnames(exoninfo)=c("ID","GeneID","geneSymbol","chr","strand","exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE")
  for (n in 1:length(name)){
    temp=strsplit(name[n],split="\\|")[[1]]
    
    exoninfo[n,"ID"]=temp[1]
    exoninfo[n,"GeneID"]=temp[2]
    exoninfo[n,"geneSymbol"]=temp[3]
    exoninfo[n,"chr"]=strsplit(temp[4],split="chr")[[1]][2]
    exoninfo[n,"strand"]=temp[5]
    exoninfo[n,"exonStart_0base"]=temp[6]
    exoninfo[n,"exonEnd"]=temp[7]
    exoninfo[n,"upstreamES"]=temp[8]
    exoninfo[n,"upstreamEE"]=temp[9]
    exoninfo[n,"downstreamES"]=temp[10]
    exoninfo[n,"downstreamEE"]=temp[11]
  }
  return(exoninfo)
}

#brain region    
cbPalette <- c("#E53537",     
               "#EEA760",     
               "#C9E5C3",     
               "#9ACB3C",     
               "#028B45",     
               "#FFD923",     
               "#F0EC68",     
               "#6D67AF",     
               "#7F388D",     
               "#2999D5",     
               "#4EC3C7",    
               "#CD8ABC",     
               "#F287B7")     

#i=4102 #MAPT-RBFOX2 example

for (i in start:end){     #for each significant pair
  exonshortid=as.character(sigdiffpair[i,"exon.shortID"])
  rbpid=as.character(sigdiffpair[i,"RBP.fullID"])
  snpid=as.character(sigdiffpair[i,"snp"])
  pvalue=as.numeric(sigdiffpair[i,"pvalue"])
  fdr=as.numeric(sigdiffpair[i,"FDR"])
  direction=as.numeric(sigdiffpair[i,"direction"])
  foldchange=as.numeric(sigdiffpair[i,"fold.change"])
  sigregion=as.character(sigdiffpair[i,"sig.region"])
  num.sigregion=as.numeric(sigdiffpair[i,"num.sig.region"])
  exonsymbol=as.character(sigdiffpair[i,"exon.symbol"])
  rbpsymbol=as.character(sigdiffpair[i,"RBP.symbol"])
  exonfullinfo=as.character(sigdiffpair[i,"exon.fullID"])
  dbrbpid=as.character(sigdiffpair[i,"DB.RBPID"])
  
  ##########################################################
  #1. make a box plot of the current wilcoxon rank sum test#
  ##########################################################
  boxplotname=paste(exonshortid,"~",snpid,"~",rbpid,"~",dbrbpid,"~",exonfullinfo,".pdf",sep="")
  Boxplot(sigregion,edataRBP,paste(exonshortid,snpid,sep="~"),rbpid,pvalue,foldchange,fdr,boxplotpath,boxplotname)

  ############################
  #2. make the 13 region plot#
  ############################
  setwd(regionplotpath)
  outfile=paste("Detailed_corplot_all_region_",PSItype,"_",counttype,"_",splicetype,"_",type,"_",
                exonshortid,"_",exonsymbol,"_",snpid,"_",
                rbpid,"_",
                i,".pdf",sep="")
  pdf(outfile,width=16,height=4)
  par(mfcol = c(2, length(brainregionlist))) 
  #layout(matrix(c(1:(2*length(brainregionlist))), 2, length(brainregionlist)), byrow = FALSE)
  par(mar = c(3, 3.5, 1, 0.5))   #margin of each individual plot
  par(mgp = c(1.5, 0.5, 0))
  par(oma = c(0, 0, 7, 1))   #outside margin (keep some room for the main title)
  
  for (br in 1:length(brainregionlist)){
    currentBR=brainregionlist[br]
    formalbrainregion=formalbrainregionlist[br]
    ###get the p value of the current sQTL in the current brain region###
    setwd(paste(sqtlrootinput,currentBR,paste("Glimmps_each_exon_cis_",currentBR,sep=""),sep="/"))
    exonresult=read.table(paste(exonshortid,".asso",sep=""),sep="\t",header=T)
    row=which(exonresult[,"SNPID"] %in% snpid)[1]
    beta=exonresult[row,"Beta"]
    p.value=exonresult[row,"pvals.lm"]
    #get the significant status of the current sQTL in the current brain region (if significant, use the brain region color, otherwise, just black)
    if (type=="pvalue"){
      cutoff=10^-5
    }
    if (type=="permutation"){
      setwd(paste(rootsqtl,"logit",counttype,splicetype,currentBR,sep="/"))
      cutoff=as.numeric(as.matrix(read.table("permutation_p_value_cutoff.txt",sep="\t")))
    }
    if (as.numeric(p.value)<=cutoff){
      #color.title=cbPalette[br]
      color.title="tan3"
    }else{
      color.title="black"
    }
    #get the sample ID of samples in the current brain region
    setwd(paste(strsplit(rootsqtl,split="sQTL_run")[[1]],"/input_splicing/",PSItype,"/",counttype,"/",splicetype,"/",currentBR,sep=""))
    IDs.pheno=as.character(as.matrix(read.table("SRR_ID.txt",sep="\t")))
    
    ###glimmpse plot###
    #get exon information as title
    temp=strsplit(exonfullinfo,split="\\|")[[1]]
    chrom=temp[4]      #the genotype information is by chromosome
    genotypename = paste(rootsqtl,"/","logit","/",counttype,"/",splicetype,"/",currentBR,sep="")
    genoplinkfile= genoplinkprefix
    genofile = genoplinkfile
    coordinate=paste(temp[4],":",temp[6],"-",temp[7],sep="")
    exonsymbol=temp[3]
    #.map file for the current chromosome
    c <-fread(paste(genotypename,"/",chrom,".",genofile, ".map", sep=""), header=F, sep="\t")    # Note: check sep (My code)
    map <- data.frame(c)
    names(map) <- c("Chr","SNPID","cM","Pos")
    #.raw file for the current chromosome
    d <- fread(paste(genotypename,"/",chrom, ".", genofile, ".map_tpose.raw", sep=""), header=T,sep="\t",  na.strings=c("NA"))     # (My code)
    genoplink <- data.frame(d)
    genoplink=genoplink[,-1] 
    IDs.geno <- names(genoplink)           #sample ID of samples with genotype information in the current brain region
    IDs.common <- intersect(IDs.geno,IDs.pheno)
    nsnps <- dim(genoplink)[1] -5                              #remove the first a few lines
    geno <- as.matrix(genoplink[seq(1, nsnps)+5, IDs.common])   #this is the genotype table for all SNPs
    #get genotype information
    exongenoname=snpid    #name of the SNP that correlates with the given exon
    exongeno=geno[which(map[,"SNPID"]==exongenoname),]
    #get exon inclusion level and total read count#
    exonpsi=PSI[exonfullinfo,]
    totalcount=totalRC[exonfullinfo,]
    # only take those individuals with genotypes from the count table, and sort the phenotype to be the same order as in the genotype matrix
    psi=exonpsi[,IDs.common]
    count=totalcount[,IDs.common]
    #5. make the plot
    allele0=strsplit(exongenoname,split="_")[[1]][3]     #reference allele, the value is 0
    allele1=strsplit(exongenoname,split="_")[[1]][4]     #alternative allele, the value is 1
    #0/0 -> 1, 0/1 and 1/0 -> 1, 1/1 -> 2
    Alleles<-c(paste(allele0,allele0,sep="/"),paste(allele0,allele1,sep="/"),paste(allele1,allele1,sep="/"))
    n<- as.numeric(as.matrix(count))
    SNP<- exongeno
    Psi <- as.numeric(as.matrix(psi))
    N <- length(n)  #number of samples
    title=paste(formalbrainregion,
                paste("P value: ",paste(formatC(p.value, format = "e", digits = 2),collapse=",   "),sep=""),
                paste("Beta: ",paste(round(beta, digits = 3),collapse=",   "),sep=""),
                sep="\n")
    ylim.range <- c(0,1) #range(psi,na.rm=T)
    plot(jitter(SNP,factor=0.5), Psi,xlim=c(-0.25,2.25), ylab="",xlab="",xaxt="n",type="n" ,ylim=ylim.range, cex.main=1)
    points(jitter(SNP,factor=0.5)  ,Psi  , pch= 19, cex= log10(n+1)/1 ,col=color.title)
    mtext(text=title, side=3, line=0.2, col=color.title, cex=0.6)
    mtext(text=Alleles, side=1, at= c(0,1,2),cex=0.7,line= 0.3)
    par("new"=T) # add boxplot on top
    boxplot(Psi~SNP, ylab="", xlab="", xaxt="n", yaxt="n",boxwex=0.35, xlim=c(-0.25,2.25), ylim=ylim.range,at=sort(unique(SNP[!is.na(Psi)][!is.na(SNP[!is.na(Psi)])])) ,border=gray(0.45),col=rgb(1,1,1,alpha=0.6),outline=FALSE)  
    
    ###RBP expression plot###
    #get the expression level of samples in the current brain region
    ylim.range=range(edataRBP[rbpid,])
    exprBR=as.numeric(edataRBP[rbpid,which(brainregion %in% currentBR)])
    label=as.integer(rep(1,length(exprBR)))
    plot(jitter(label,factor=0.5), exprBR,xlim=c(0.9,1.1), ylab="",xlab="",xaxt="n",type="n" ,ylim=ylim.range, cex.main=1)
    points(jitter(label,factor=0.5)  ,exprBR  , pch= 19, col=color.title)
    title=paste(formalbrainregion,
                paste("Mean TPM: ",paste(round(mean(exprBR), digits = 3),collapse=",   "),sep=""),
                sep="\n")
    mtext(text=title, side=3, line=0.2, col=color.title, cex=0.6)
    par("new"=T) # add boxplot on top
    boxplot(exprBR~label,
            #ylab="", xlab="", xaxt="n", yaxt="n", xlim=c(0.9,1.1),
            ylim=ylim.range,at=sort(unique(label[!is.na(exprBR)][!is.na(label[!is.na(exprBR)])])),
            border=gray(0.45),col=rgb(1,1,1,alpha=0.6),outline=FALSE)
  }
  main.title1=paste(exonsymbol,snpid,rbpsymbol,sep="      ")
  main.title2=paste(paste("Fold change:",paste(round(foldchange, digits = 3),collapse=",   ")),
                    paste("pval:",paste(formatC(pvalue, format = "e", digits = 2),collapse=",   ")),
                    paste("FDR:",paste(formatC(fdr, format = "e", digits = 2),collapse=",   ")),sep="      ")
  main.title=paste(main.title1,
                   main.title2,
                   sep="\n")
  mtext(main.title, outer = TRUE,side = 3, cex = 0.7, line = 2)
  dev.off()
  
  ######################
  #3. make sashimi plot#
  ######################
  for (sigbr in strsplit(sigregion,split=", ")[[1]]){
    #2. get the sample ID of samples in the current brain region
    setwd(paste(strsplit(rootsqtl,split="sQTL_run")[[1]],"/input_splicing/",PSItype,"/",counttype,"/",splicetype,"/",sigbr,sep=""))
    IDs.pheno=as.character(as.matrix(read.table("SRR_ID.txt",sep="\t")))
    temp=strsplit(exonfullinfo,split="\\|")[[1]]
    chrom=temp[4]      #the genotype information is by chromosome
    genotypename = paste(rootsqtl,"/",PSItype,"/",counttype,"/",splicetype,"/",sigbr,sep="")
    genoplinkfile= genoplinkprefix
    genofile = genoplinkfile
    coordinate=paste(temp[4],":",temp[6],"-",temp[7],sep="")
    genesymbol=temp[3]
    #.map file for the current chromosome
    c <-fread(paste(genotypename,"/",chrom,".",genofile, ".map", sep=""), header=F, sep="\t")    # Note: check sep (My code)
    map <- data.frame(c)
    names(map) <- c("Chr","SNPID","cM","Pos")
    #.raw file for the current chromosome
    d <- fread(paste(genotypename,"/",chrom, ".", genofile, ".map_tpose.raw", sep=""), header=T,sep="\t",  na.strings=c("NA"))     # (My code)
    genoplink <- data.frame(d)
    genoplink=genoplink[,-1] 
    IDs.geno <- names(genoplink)           #sample ID of samples with genotype information in the current brain region
    IDs.common <- intersect(IDs.geno,IDs.pheno)
    nsnps <- dim(genoplink)[1] -5                              #remove the first a few lines
    geno <- as.matrix(genoplink[seq(1, nsnps)+5, IDs.common])   #this is the genotype table for all SNPs
    exongenoname=snpid    #name of the SNP that correlates with the given exon
    exongeno=geno[which(map[,"SNPID"]==exongenoname),]
    #4. get exon inclusion level and total read count#
    exonpsi=PSI[exonfullinfo,]
    totalcount=totalRC[exonfullinfo,]
    # only take those individuals with genotypes from the count table, and sort the phenotype to be the same order as in the genotype matrix
    psi=exonpsi[,IDs.common]
    count=totalcount[,IDs.common]
    allele0=strsplit(exongenoname,split="_")[[1]][3]     #reference allele, the value is 0
    allele1=strsplit(exongenoname,split="_")[[1]][4]     #alternative allele, the value is 1
    #0/0 -> 1, 0/1 and 1/0 -> 1, 1/1 -> 2
    Alleles<-c(paste(allele0,allele0,sep="/"),paste(allele0,allele1,sep="/"),paste(allele1,allele1,sep="/"))
    sampleID0=names(which(exongeno==0))
    sampleID1=names(which(exongeno==1))
    sampleID2=names(which(exongeno==2))
    samplesize=c(length(sampleID0),length(sampleID1),length(sampleID2))
    ###make sashimi plot###
    rootbam="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/data/V7_all_bam_soft_link/"   #the folder containing all BAM files
    
    exoninfo=getexon(exonfullinfo)
    #generate folders for the plot
    setwd(sashimiplotpath)
    outputfolder=paste(PSItype,"_",counttype,"_",splicetype,"_",type,"_",gsub("\\s", "", formalbrainregionlist[which(brainregionlist %in% sigbr)]),"_",paste(strsplit(exonfullinfo,split="\\|")[[1]][-c(2,5)],collapse="_"),"_",exongenoname,sep="")
    if (file.exists(outputfolder)==FALSE){   #if the folder doesn't exist, we create the folder. Otherwise, we don't do anything
      command=paste("mkdir -p",outputfolder)
      system(command)
      
      setwd(paste(sashimiplotpath,outputfolder,sep="/"))
      
      ###prepare input files for sashimi plot###
      sortexongeno=sort(exongeno)
      sample_target=sample_control=c()  #we don't have target and control here but in order to use the previous code, we just refer to the first genotype as target and the rest as control
      sig_group=Alleles[which(samplesize!=0)][1]
      control_group=Alleles[which(samplesize!=0)][-1]
      
      grouping=matrix(NA,length(samplesize[samplesize!=0]),1)      #the file contains the grouping information of samples
      for (g in 1:length(control_group)){   
        group=control_group[g]
        #1. b1 and b2 parameter
        sample_target=names(sortexongeno[sortexongeno==which(Alleles==sig_group)-1])
        temp_control=names(sortexongeno[sortexongeno==(which(Alleles==group)-1)])
        sample_control=c(sample_control,temp_control)
        
        #2. grouping file
        grouping[1,]=paste(sig_group,": ",1,"-",length(sample_target),sep="")
        temp=as.numeric(tail(strsplit(grouping[g,],split="-")[[1]],n=1))
        grouping[g+1,]=paste(group,": ",temp+1,"-",temp+length(temp_control),sep="")
      }
      write.table(grouping,"grouping.gf",sep="\t",row.names = F,col.names = F,quote=F)
      
      b1=paste(paste(rootbam,sample_target,".bam",sep=""),collapse=",")
      b2=paste(paste(rootbam,sample_control,".bam",sep=""),collapse=",")
      
      #3. rMATS format input
      rmats=matrix(NA,dim(exoninfo)[1],23)    #the input file for sashimi plot
      colnames(rmats)=c("ID","GeneID","geneSymbol","chr","strand",
                        colnames(fromGTF)[6:11],
                        #"exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE",
                        "ID",	"IC_SAMPLE_1",	"SC_SAMPLE_1",	"IC_SAMPLE_2",	"SC_SAMPLE_2",	"IncFormLen",	"SkipFormLen",
                        "PValue",	"FDR",	"IncLevel1",	"IncLevel2",	"IncLevelDifference")
      rmats[,1:11]=exoninfo
      rmats[,"PValue"]=rmats[,"FDR"]=0
      rmats[,12]=rmats[,1]
      rownames(rmats)=rownames(exoninfo)  #we need to remove row name when outputting this into file
      
      for(e in 1:dim(rmats)[1]){
        EXON=rownames(exoninfo)[e]
        rmats[EXON,"IC_SAMPLE_1"]=paste(IC[EXON,sample_target],collapse=",")
        rmats[EXON,"SC_SAMPLE_1"]=paste((totalRC[EXON,sample_target]-IC[EXON,sample_target]),collapse=",")
        rmats[EXON,"IC_SAMPLE_2"]=paste(IC[EXON,sample_control],collapse=",")
        rmats[EXON,"SC_SAMPLE_2"]=paste((totalRC[EXON,sample_control]-IC[EXON,sample_control]),collapse=",")
        rmats[EXON,"IncFormLen"]=effectlength[EXON,"IncFormLen"]
        rmats[EXON,"SkipFormLen"]=effectlength[EXON,"SkipFormLen"]
        rmats[EXON,"IncLevel1"]=paste(PSI[EXON,sample_target],collapse=",")
        rmats[EXON,"IncLevel2"]=paste(PSI[EXON,sample_control],collapse=",")
        rmats[EXON,"IncLevelDifference"]=mean(as.matrix(PSI)[EXON,sample_target])-mean(as.matrix(PSI)[EXON,sample_control])
      }
      write.table(rmats,"rmats_events.txt",sep="\t",row.names=F,quote=F)
      
      #run sashimi plot command
      command=paste(#"source /u/local/Modules/default/init/modules.sh",
        #"\n",
        #"module load samtools",       
        #"\n",
        "/u/nobackup/yxing/PROJECT/yidazhan/research/software/anaconda2/bin/rmats2sashimiplot",
        "--b1", b1,
        "--b2", b2,
        "-t", splicetype,
        "-e", "rmats_events.txt",
        "--l1", sig_group,
        "--l2", "other",
        "--exon_s", 1,
        "--intron_s", 1,
        "-o", "events_output",
        "--group-info", "grouping.gf",
        "--min-counts", 0,   
        "--font-size", 8,
        sep=" ")
      system(command)
    }
  }
}




