splicetype="SE"
PSItype="original"    #the calculation is based on the corrected PSI but when we make the plot, we use the original value
counttype="JC"
type="pvalue"
exonshortidlist=c("SE_108190",  #SH3YL1 
                  "SE_108171",  #SH3YL1 
                  "SE_229974",  #CFAP44 
                  "SE_222152",  #Flot1
                  "SE_344917")  #SLC39A13
snpidlist=c("2_242793_G_C_b37",
            "2_224919_A_G_b37",
            "3_113012797_G_A_b37",
            "6_30708695_C_T_b37",
            "11_47434986_G_A_b37")
exonsymbollist=c("SH3YL1",
                 "SH3YL1",
                 "CFAP44",
                 "FLOT1",
                 "SLC39A13")
exonfullinfolist=c("108190|ENSG00000035115.21_3|SH3YL1|chr2|-|242797|242871|234159|234272|247537|247602",
                   "108171|ENSG00000035115.21_3|SH3YL1|chr2|-|224863|224920|218148|219001|229965|230044",
                   "229974|ENSG00000206530.10_3|CFAP44|chr3|-|113012798|113012885|113009699|113010595|113013533|113013668",
                   "222152|ENSG00000137312.14_2|FLOT1|chr6|-|30708966|30709110|30708455|30708575|30709390|30709481",
                   "344917|ENSG00000165915.13_2|SLC39A13|chr11|+|47434950|47435058|47433852|47434018|47435147|47435237")



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

outputpath="/u/nobackup/yxing-BIGDATA/yidazhan/GTEx_V7_analysis/10.6_Glimmpse_plot_for_all_brain_regions_given_RBP_exon_SNP/result"

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
brainregion=oribrainregion


for (i in 1:length(exonshortidlist)){     #for each significant pair
  exonshortid=exonshortidlist[i]
  snpid=snpidlist[i]
  exonsymbol=exonsymbollist[i]
  exonfullinfo=exonfullinfolist[i]

  #########################
  #make the 13 region plot#
  #########################
  setwd(outputpath)
  outfile=paste("Detailed_corplot_all_region_",PSItype,"_",counttype,"_",splicetype,"_",type,"_",
                exonshortid,"_",exonsymbol,"_",snpid,"_",".pdf",sep="")
  pdf(outfile,width=16,height=4)
  par(mfcol = c(2, length(brainregionlist))) 
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
    boxplot(Psi~SNP, ylab="", xlab="", xaxt="n", yaxt="n",boxwex=0.35, xlim=c(-0.25,2.25), ylim=ylim.range,at=sort(unique(SNP[!is.na(Psi)][!is.na(SNP[!is.na(Psi)])])) ,border=gray(0.45),col=rgb(1,1,1,alpha=0.6),outline=FALSE) 
  }
  main.title=paste(exonsymbol,snpid,sep="      ")
  mtext(main.title, outer = TRUE,side = 3, cex = 0.7, line = 2)
  dev.off()
}




