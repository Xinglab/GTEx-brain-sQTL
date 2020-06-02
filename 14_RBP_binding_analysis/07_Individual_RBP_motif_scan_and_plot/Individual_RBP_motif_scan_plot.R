splicetype="SE"      #type of alternative splicing, e.g., SE, A3SS, A5SS, MXE, IR
type="pvalue"       #pvalue or permutation

library(robust)
library(gplots)
library(ggplot2)
library(reshape)
library(RColorBrewer)
require(scales)
library(corrplot)
library(ggsci)
library(ggseqlogo)
library("data.table")
library(gaston)
library(boot)

inputpath=paste("/path/to/DeepBind/deepbind_individual_RBP_motif_scan",
                          splicetype,type,sep="/")

outputpath=paste("/path/to/DeepBind/deepbind_individual_RBP_motif_scan",
                 splicetype,type,sep="/")
outputseqlogo=paste(outputpath,"significant_pair_seqlogo",sep="/")
outputbox=paste(outputpath,"significant_pair_boxplot",sep="/")
outputdetailed=paste(outputpath,"significant_pair_detailed_plot",sep="/")
outputsashimi=paste(outputpath,"significant_pair_sashimi_plot",sep="/")
command=paste("mkdir -p",outputseqlogo)
system(command)
command=paste("mkdir -p",outputbox)
system(command)
command=paste("mkdir -p",outputdetailed)
system(command)
command=paste("mkdir -p",outputsashimi)
system(command)


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

RBP="NOVA"
RBPexpr=c("ENSG00000139910.15_NOVA1","ENSG00000104967.6_NOVA2")

pvaluecutoff=0.05
foldchangecutoff=1.5

###############################
#read in the RBP related sQTLs#
###############################
setwd(inputpath)
rbpsqtllist=read.table(paste(splicetype,"_",type,"_individual_RBP_motif_scan_",RBP,".txt",sep=""),sep="\t")
#there are cases that the no matter with or without the mutation, there is always a motif, we don't want to focus on those cases
motifdiff=rbpsqtllist[,paste(RBP,".motif.num.wt",sep="")]-rbpsqtllist[,paste(RBP,".motif.num.mut",sep="")]
subrbpsqtllist=rbpsqtllist[which(motifdiff!=0),]

#######################
#splicing and genotype#
#######################
PSItype="logit"
counttype="JC"
inputpath="/path/to/corrected/PSI/value"
summaryinput=paste("/path/to/summary/",PSItype,"/",counttype,"/",splicetype,sep="")
rootsqtl="/path/to/sQTL_run"
totalcountinput="./01_Get_PSI_from_rMATS_output/example_output"
GWASdbpath="/path/to/GWAS_catalog/files"
GWASdbname="gwas_catalog_v1.0.1-associations_e89_r2017-07-31.tsv"
genoplinkprefix="genotype_file_name_prefix"
LDpath="/path/to/LD/result"
sqtlrootinput=paste("/path/to/sQTL_run/logit/JC",
                    splicetype,sep="/")
rmatspostpath="/path/to/rMATS/post/result"

#read in the PSI value
inputPSI=paste(splicetype,"_",counttype,"_logit_corrected_PSI.txt",sep="")
inputtotalRC=paste(splicetype,"_",counttype,"_totalRC_filter.txt",sep="")
inputIC=paste(splicetype,"_",counttype,"_IC_filter.txt",sep="")
setwd(inputpath)
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
elinput="./01_Get_PSI_from_rMATS_output/example_output"
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

#################
#read in the TPM#
#################
setwd("/input/to/GTEx/processed/gene/expression/data")
normTPM=read.table("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_brain_normalized.txt",sep="\t",header=T,check.names=F)
edata=normTPM
rownames(edata)=paste(as.character(normTPM[,1]),as.character(normTPM[,2]),sep="_")
edata=edata[,-c(1,2)]

##########################
#brain region information#
##########################
phenotypepath="./03_Get_sample_annotation/example_output"
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
setwd("/input/to/GTEx/processed/gene/expression/data")
IDconversion=read.table("SRRID_2_GTEXID.txt",sep="\t",header=T)
rownames(IDconversion)=IDconversion[,"sampleID"]

colnames(edata)=as.character(IDconversion[colnames(edata),"SRRlist"])
edataRBP=as.matrix(edata[RBPexpr,])
rownames(edataRBP)=RBPexpr
brainregion=oribrainregion[,as.character(colnames(edata))]

#make box plot for significant pairs
Boxplot=function(SIGREGION,EXPR,sQTLID,RBPID,PV,FC,path,name){
  sigTPM=as.numeric(as.matrix(EXPR[RBPID,colnames(brainregion)[which(brainregion %in% SIGREGION)]]))
  insigTPM=as.numeric(as.matrix(EXPR[RBPID,setdiff(colnames(brainregion),colnames(brainregion)[which(brainregion %in% SIGREGION)])]))
  exprBR=c(sigTPM,insigTPM)
  label=c(rep("Significant",length(sigTPM)),rep("Insignificant",length(insigTPM)))
  ylim.range=range(exprBR)
  data2plot=data.frame("TPM"=exprBR,"Group"=label)
  
  setwd(path)
  pdf(name,width=4,height=4)
  title=paste(sQTLID,
              RBPID,
              paste(SIGREGION,collapse=","),
              paste("P value:",paste(formatC(PV, format = "e", digits = 3),collapse=",   ")),
              #paste("FDR: ",paste(formatC(FDR, format = "e", digits = 3),collapse=",   ")),
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
  #the header here doesn't matter. It won't affect any downstream analysis
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
cbPalette <- c("#E53537",     #GTEx Artery color/Amygdala
               "#EEA760",     #GTEx adipose tissue color/ACC
               "#C9E5C3",     #GTEx pituitary color/Caudate
               "#9ACB3C",     #GTEx lung color/CH
               "#028B45",     #GTEx thyroid color/Cerebellum
               "#FFD923",     #GTEx nerve color/cortex
               "#F0EC68",     #GTEx brain color/frontal cortex
               "#6D67AF",     #GTEx muscle color/hippocampus
               "#7F388D",     #GTEx Heart color/Hypothalamus
               "#2999D5",     #GTEx Skin color/Nucleus accumbens
               "#4EC3C7",     #GTEx breast color/Putamen
               "#CD8ABC",     #GTEx LCL color/spinal cord
               "#F287B7")     #GTEx blood color/Substantia nigra

############################################################################
#for each trio, do all the calculation, and make plots for significant ones#
############################################################################
summarystats=matrix(NA,dim(subrbpsqtllist)[1],3*length(RBPexpr))
colnames(summarystats)=c(sort(as.vector(outer(RBPexpr, c("p.value","fold.change","direction"), paste, sep=".")),decreasing=T))

for (i in 1:dim(subrbpsqtllist)[1]){
  print(i)
  exonshortID=as.character(subrbpsqtllist[i,"Exon"])
  exonfullID=as.character(subrbpsqtllist[i,"fullexoninfo"])
  exonsymbol=strsplit(exonfullID,split="\\|")[[1]][3]
  snp=as.character(subrbpsqtllist[i,"SNP"])
  sigregionlist=strsplit(as.character(subrbpsqtllist[i,"Brain.region"]),split=", ")[[1]]
  
  ###########################################
  #generate plot to show the sequence change#
  ###########################################
  wildseq=as.character(subrbpsqtllist[i,paste(RBP,"seq.wild",sep=".")])
  mutseq=as.character(subrbpsqtllist[i,paste(RBP,"seq.mut",sep=".")])
  seq2plot=c(wildseq,mutseq)
  setwd(outputseqlogo)
  outfile=paste("Seqlogo_",PSItype,"_",counttype,"_",splicetype,"_",type,"_",
                exonshortID,"_",exonsymbol,"_",snp,"_",
                RBP,"_",
                i,".pdf",sep="")
  pdf(outfile,width=6,height=4)
  p=ggseqlogo(seq2plot,method = 'prob')
  print(p)
  dev.off()
  
  for (rbp in RBPexpr){
    ################################################
    #1. differential expression analysis + box plot#
    ################################################
    #1. for sQTLs with significant and insignificant region
    if (length(sigregionlist)<length(brainregionlist)){
      #calculate the differential expression between significant and insignificant regions
      sigTPM=as.numeric(as.matrix(edataRBP[rbp,colnames(brainregion)[which(brainregion %in% sigregionlist)]]))
      insigTPM=as.numeric(as.matrix(edataRBP[rbp,setdiff(colnames(brainregion),colnames(brainregion)[which(brainregion %in% sigregionlist)])]))
      difftest=wilcox.test(sigTPM,insigTPM)
      
      #add the p value, fold change and direction to summarystats
      summarystats[i,paste(rbp,"p.value",sep=".")]=difftest$p.value
      if (mean(sigTPM,na.rm=T)>mean(insigTPM,na.rm=T)){
        summarystats[i,paste(rbp,"fold.change",sep=".")]=mean(sigTPM,na.rm=T)/mean(insigTPM,na.rm=T)
        summarystats[i,paste(rbp,"direction",sep=".")]=1
      }else{
        summarystats[i,paste(rbp,"fold.change",sep=".")]=mean(insigTPM,na.rm=T)/mean(sigTPM,na.rm=T)
        summarystats[i,paste(rbp,"direction",sep=".")]=-1
      }
      
      if (summarystats[i,paste(rbp,"p.value",sep=".")]<=pvaluecutoff & summarystats[i,paste(rbp,"fold.change",sep=".")]>=foldchangecutoff) {    #if the p value and fold change is significant
        #generate the box plot
        boxplotname=paste(exonshortID,"~",snp,"~",rbp,"~",exonfullID,".pdf",sep="")
        Boxplot(sigregionlist,edataRBP,paste(exonshortID,snp,sep="~"),rbp,
                summarystats[i,paste(rbp,"p.value",sep=".")],
                summarystats[i,paste(rbp,"fold.change",sep=".")],
                outputbox,boxplotname)
      }
    }
    
    
    ############################
    #2. make the 13 region plot#
    ############################
    setwd(outputdetailed)
    outfile=paste("Detailed_corplot_all_region_",PSItype,"_",counttype,"_",splicetype,"_",type,"_",
                  exonshortID,"_",exonsymbol,"_",snp,"_",
                  rbp,"_",
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
      exonresult=read.table(paste(exonshortID,".asso",sep=""),sep="\t",header=T)
      row=which(exonresult[,"SNPID"] %in% snp)[1]
      beta=exonresult[row,"Beta"]
      p.value=exonresult[row,"pvals.lm"]
      #get the significant status of the current sQTL in the current brain region (if significant, use the brain region color, otherwise, just black)
      if (type=="pvalue"){
        cutoff=10^-5
      }
      if (type=="permutation"){
        setwd(paste(rootsqtl,PSItype,counttype,splicetype,currentBR,sep="/"))
        cutoff=as.numeric(as.matrix(read.table("permutation_p_value_cutoff.txt",sep="\t")))
      }
      if (as.numeric(p.value)<=cutoff){
        color.title=cbPalette[br]
      }else{
        color.title="black"
      }
      #get the sample ID of samples in the current brain region
      setwd(paste(strsplit(rootsqtl,split="sQTL_run")[[1]],"/input_splicing/",PSItype,"/",counttype,"/",splicetype,"/",currentBR,sep=""))
      IDs.pheno=as.character(as.matrix(read.table("SRR_ID.txt",sep="\t")))
      
      ###RBP expression plot###
      #get the expression level of samples in the current brain region
      ylim.range=range(edataRBP[rbp,])
      exprBR=as.numeric(edataRBP[rbp,which(brainregion %in% currentBR)])
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
      
      ###glimmpse plot###
      #get exon information as title
      temp=strsplit(exonfullID,split="\\|")[[1]]
      chrom=temp[4]      #the genotype information is by chromosome
      genotypename = paste(rootsqtl,"/",PSItype,"/",counttype,"/",splicetype,"/",currentBR,sep="")
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
      exongenoname=snp    #name of the SNP that correlates with the given exon
      exongeno=geno[which(map[,"SNPID"]==exongenoname),]
      #get exon inclusion level and total read count#
      exonpsi=PSI[exonfullID,]
      totalcount=totalRC[exonfullID,]
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
      plot(jitter(SNP,factor=0.5), Psi,xlim=c(-0.25,2.75), ylab="",xlab="",xaxt="n",type="n" ,ylim=ylim.range, cex.main=1)
      points(jitter(SNP,factor=0.5)  ,Psi  , pch= 19, cex= log10(n+1)/1 ,col=color.title)
      mtext(text=title, side=3, line=0.2, col=color.title, cex=0.6)
      mtext(text=Alleles, side=1, at= c(0,1,2),cex=0.7,line= 0.3)
      par("new"=T) # add boxplot on top
      boxplot(Psi~SNP, ylab="", xlab="", xaxt="n", yaxt="n",boxwex=0.35, xlim=c(-0.25,2.75), ylim=ylim.range,at=sort(unique(SNP[!is.na(Psi)][!is.na(SNP[!is.na(Psi)])])) ,border=gray(0.45),col=rgb(1,1,1,alpha=0.6),outline=FALSE)  
    }
    main.title1=paste(exonsymbol,snp,strsplit(rbp,split="_")[[1]][2],sep="      ")
    if (length(sigregionlist)<length(brainregionlist)){
      main.title2=paste(paste("Fold change:",paste(round(summarystats[i,paste(rbp,"fold.change",sep=".")], digits = 3),collapse=",   ")),
                        paste("pval:",paste(formatC(summarystats[i,paste(rbp,"p.value",sep=".")], format = "e", digits = 2),collapse=",   ")))
      main.title=paste(main.title1,
                       main.title2,
                       sep="\n")
    }else{
      main.title=main.title1
    }
    mtext(main.title, outer = TRUE,side = 3, cex = 0.7, line = 2)
    dev.off()
    
    ######################
    #3. make sashimi plot#
    ######################
    for (sigbr in sigregionlist){
      #2. get the sample ID of samples in the current brain region
      setwd(paste(strsplit(rootsqtl,split="sQTL_run")[[1]],"/input_splicing/",PSItype,"/",counttype,"/",splicetype,"/",sigbr,sep=""))
      IDs.pheno=as.character(as.matrix(read.table("SRR_ID.txt",sep="\t")))
      temp=strsplit(exonfullID,split="\\|")[[1]]
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
      exongenoname=snp    #name of the SNP that correlates with the given exon
      exongeno=geno[which(map[,"SNPID"]==exongenoname),]
      #4. get exon inclusion level and total read count#
      exonpsi=PSI[exonfullID,]
      totalcount=totalRC[exonfullID,]
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
      rootbam="/path/to/BAM/files/"   #the folder containing all BAM files
      
      exoninfo=getexon(exonfullID)
      #generate folders for the plot
      setwd(outputsashimi)
      outputfolder=paste(PSItype,"_",counttype,"_",splicetype,"_",type,"_",gsub("\\s", "", formalbrainregionlist[which(brainregionlist %in% sigbr)]),"_",paste(strsplit(exonfullID,split="\\|")[[1]][-c(2,5)],collapse="_"),"_",exongenoname,sep="")
      if (file.exists(outputfolder)==FALSE){   #if the folder doesn't exist, we create the folder. Otherwise, we don't do anything
        command=paste("mkdir -p",outputfolder)
        system(command)
        
        setwd(paste(outputsashimi,outputfolder,sep="/"))
        
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
          "/path/to/rmats2sashimiplot",
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
}

summarystats=cbind(subrbpsqtllist,summarystats)
setwd(outputpath)
write.table(summarystats,paste(splicetype,"_",type,"_individual_RBP_motif_scan_",RBP,"_with_stats.txt",sep=""),sep="\t")
