#1.1 get the correlation coefficient and p value of all sQTLs in all 13 brain regions
#1.2 calculate coefficient of variation for each pair
#1.3 classify events into three groups: high CV, medium CV, low CV
#1.4 plot this region dependent sQTL using the beta value
#1.5 for each group, pick one example and show the glimmpse plot for given brain regions

splicetype="SE"
type="pvalue"
PSItype="original"     #the calculation is based on the corrected PSI but when we make the plot, we use the original value
counttype="JC"

library(robust)
library(gplots)
library(ggplot2)
library(reshape)
library(RColorBrewer)
require(scales)
library(corrplot)

'%!in%' <- function(x,y)!('%in%'(x,y))
source("/path/to/heatmap.3.R")

vectorsplit=function(x,del,pos){
  return(strsplit(x,split=del)[[1]][pos])
}

cv=function(numlist){
  return(sd(numlist)/mean(numlist))
}

zscore=function(x){
  return((x-mean(x))/sd(x))     #z-score
}

outputpath=paste("/output/path/DeepBind/deepbind_region_dependent_sQTL_plot_based_on_beta",
                 splicetype,type,sep="/")
command=paste("mkdir -p",outputpath)
system(command)

sqtlrootinput=paste("/path/to/sQTL_run/logit/JC",
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

#read in the joblist
joblistinput=paste("/output/path/DeepBind/deepbind_input/",splicetype,"/",type,sep="")
setwd(joblistinput)
uniquejoblist=read.table(paste(splicetype,"_",type,"_uniquejoblist.txt",sep=""),sep="\t")

#get the full exon information
#get all sQTL exons
rootinput="/path/to/summary/logit/JC"
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

###########
#section 1#
###########
#for each exon, we get all significant SNPs within 300bp
#for each exon-SNP pair, we get its correlation coefficient & p value in all 13 brain regions. 
#we combine this information into a matrix (row: exon-SNP pair. column: region)
#then we did a clustering to get region specific and region ubiquitous sQTLs
exonsnppair=unique(uniquejoblist[,c("Exon","SNP")])
#filter out SNPs which are not single base 
allele0=sapply(as.character(exonsnppair[,"SNP"]),vectorsplit,"_",3)
allele1=sapply(as.character(exonsnppair[,"SNP"]),vectorsplit,"_",4)
allele0length=sapply(allele0,nchar)
allele1length=sapply(allele1,nchar)
singlebase=intersect(which(allele0length==1),which(allele1length==1))
singlebaseexonsnppair=exonsnppair[singlebase,]
betamatrix=pvmatrix=matrix(NA,dim(singlebaseexonsnppair)[1],length(brainregionlist))
rownames(betamatrix)=rownames(pvmatrix)=paste(singlebaseexonsnppair[,"Exon"],singlebaseexonsnppair[,"SNP"],sep="~")
colnames(betamatrix)=colnames(pvmatrix)=brainregionlist

#fill in the coefficient and p value matrix
for (i in 1:dim(singlebaseexonsnppair)[1]){
  exonshortID=as.character(singlebaseexonsnppair[i,"Exon"])
  snpID=as.character(singlebaseexonsnppair[i,"SNP"])
  print(i)
  for (br in 1:length(brainregionlist)){
    currentBR=brainregionlist[br]
    #read in the exon result in the current brain region
    inputpath=paste(sqtlrootinput,currentBR,paste("Glimmps_each_exon_cis_",currentBR,sep=""),sep="/")
    setwd(inputpath)
    exonresult=read.table(paste(exonshortID,".asso",sep=""),sep="\t",header=T)
    row=which(exonresult[,"SNPID"] %in% snpID)[1]
    betamatrix[i,br]=exonresult[row,"Beta"]
    pvmatrix[i,br]=exonresult[row,"pvals.lm"]
  }
}

setwd(outputpath)
write.table(betamatrix,paste("betamatrix_",splicetype,"_",type,".txt",sep=""),sep="\t")
write.table(pvmatrix,paste("pvmatrix_",splicetype,"_",type,".txt",sep=""),sep="\t")


###########
#section 2#
###########
#1. add disease information
diseaselabel=matrix(NA,dim(pvmatrix)[1],2)
rownames(diseaselabel)=rownames(pvmatrix)
exonlistshortID=as.character(singlebaseexonsnppair[,"Exon"])
exonlistfullID=exonIDconversion[exonlistshortID,"fullID"]

#get the disease exon list
setwd("/path/to/disease/exon/list")
disexoninfo=read.table(paste("disease_exon_summary_",splicetype,"_",type,".txt",sep=""),sep="\t",header=T,check.names=F)

for (i in 1:dim(diseaselabel)[1]){
  if (exonlistfullID[i] %in% rownames(disexoninfo)){
    diseaselabel[i,1]=1
    diseaselabel[i,2]=paste(colnames(disexoninfo)[which(disexoninfo[exonlistfullID[i],]==1)],collapse=",")
  }else{
    diseaselabel[i,1]=0
  }
}

#2. calculate CV and add that to the matrix
CVlist=apply(betamatrix,1,cv)
absCVlist=abs(CVlist)

setwd(outputpath)
fullinfomatrix=cbind(pvmatrix,diseaselabel,CVlist,absCVlist,exonlistfullID)
write.table(fullinfomatrix,"fullinfo.txt",sep="\t")

#3. split the matrix into 2 groups based on CV cutoff
#rank the events by absolute CV value
orderedfullinfomatrix=fullinfomatrix[names(sort(absCVlist,decreasing=T)),]
#region specific group: top 150
spergbeta=betamatrix[rownames(orderedfullinfomatrix)[1:150],]
#region ubiquitous group: bottom 150
ubirgbeta=betamatrix[rownames(orderedfullinfomatrix)[(dim(orderedfullinfomatrix)[1]-150):dim(orderedfullinfomatrix)[1]],]

#4. make a heatmap (z score transformed beta value)
for (group in 1:2){
  for (transformation in 1:2){
    if (group==1){
      #region specific group
      data2plot=spergbeta
      label1="Region_specific"
    }
    if (group==2){
      #region ubiquitous group
      data2plot=ubirgbeta
      label1="Region_ubiquitous"
    }
    
    if (transformation==1){
      #original beta value
      zbeta=data2plot
      label2="original_beta"
    }
    if (transformation==2){
      #z score transformation of the beta value
      zbeta=t(apply(data2plot,1,zscore))
      label2="z_score_transformed_beta"
    }
    
    setwd(outputpath)
    pdf(paste(label1,"_sQTL_beta_clustering_heatmap_",splicetype,"_",type,"_",label2,".pdf",sep=""),height=10,width=10)
    ownbreak = c(seq(range(as.numeric(zbeta))[1],-0.8,length=250),seq(-0.79,0.79,length=501),seq(0.8,range(as.numeric(zbeta))[2],length=250))
    heatmap.3(zbeta, 
              col = colorRampPalette(rev(brewer.pal(11,"RdBu")))(1000), 
              key = TRUE,
              breaks = ownbreak,
              Colv=FALSE)
    dev.off()
  }
}


###########
#section 3#
###########
#generate glimmpse plot for the ubiquitous and specific events
singleoutputpath=paste(outputpath,"single_example_plot",sep="/")
command=paste("mkdir -p",singleoutputpath)
system(command)

library("data.table")
library(gaston)
library(boot)

totalcountinput="./01_Get_PSI_from_rMATS_output/example_output"
summaryinput=paste("/path/to/summary",PSItype,counttype,splicetype,sep="/")
rootsqtl="/path/to/sQTL_run"
genoplinkprefix="genotype_file_name_prefix"

singleoutputpath=paste(outputpath,"single_example_plot",sep="/")

if (PSItype=="original"){   #original PSI value without any correction
  psiinputpath="./01_Get_PSI_from_rMATS_output/example_output"
  inputPSI=paste(splicetype,"_",counttype,"_PSI_filter.txt",sep="")
}
if (PSItype=="logit"){     #logit PSI with correction (because we used corrected PSI for all the sQTL calculation, when we plot the result, we also use this but only transform it back to 0-1 scale)
  psiinputpath="/path/to/corrected/PSI/value"
  inputPSI=paste(splicetype,"_",counttype,"_logit_corrected_PSI.txt",sep="")
}
inputtotalRC=paste(splicetype,"_",counttype,"_totalRC_filter.txt",sep="")

setwd(psiinputpath)
PSI=read.table(inputPSI,sep="\t",header=T)              #PSI as the normalized inclusion count
setwd(totalcountinput)
totalRC=read.table(inputtotalRC,sep="\t",header=T)

#transform logit PSI back to PSI
if (PSItype=="logit"){
  temp=inv.logit(as.numeric(as.matrix(PSI)))
  originalPSI=matrix(temp,dim(PSI)[1],dim(PSI)[2])
  rownames(originalPSI)=rownames(PSI)
  colnames(originalPSI)=colnames(PSI)
  PSI=as.data.frame(originalPSI)
}


exonsnp=c("SE_9912~5_176798306_G_A_b37","SE_359526~20_62320968_T_C_b37")                   
brlist=c(6,4)       #the brain region to show for each example         
    

for (es in 1:length(exonsnp)){
  currentpair=exonsnp[es]
  exonshortID=strsplit(currentpair,split="~")[[1]][1]
  exon=exonIDconversion[exonshortID,"fullID"]
  snpid=strsplit(currentpair,split="~")[[1]][2]
  br=brlist[es]                                                                
  brainregion=brainregionlist[br]
  formalbrainregion=formalbrainregionlist[br]

  #2. get the sample ID of samples in the current brain region
  setwd(paste(strsplit(rootsqtl,split="sQTL_run")[[1]],"/input_splicing/",PSItype,"/",counttype,"/",splicetype,"/",brainregion,sep=""))
  IDs.pheno=as.character(as.matrix(read.table("SRR_ID.txt",sep="\t")))
  
  temp=strsplit(exon,split="\\|")[[1]]
  chrom=temp[4]      #the genotype information is by chromosome
  genotypename = paste(rootsqtl,"/","logit","/",counttype,"/",splicetype,"/",brainregion,sep="")
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
  #IDs.common <- IDs.geno[is.element(IDs.geno, IDs.pheno)]     #get the samples with genotype information
  IDs.common <- intersect(IDs.geno,IDs.pheno)
  nsnps <- dim(genoplink)[1] -5                              #remove the first a few lines
  #sub.geno <- match(IDs.common , IDs.geno)          #the position of samples in IDs.geno that are also in IDs.common
  #geno <- as.matrix(genoplink[seq(1, nsnps)+5, sub.geno])    #  I changed here (we basically removed the header and first column from d)
  geno <- as.matrix(genoplink[seq(1, nsnps)+5, IDs.common])   #this is the genotype table for all SNPs
  exongenoname=snpid    #name of the SNP that correlates with the given exon
  exongeno=geno[which(map[,"SNPID"]==exongenoname),]
  
  #4. get exon inclusion level and total read count#
  exonpsi=PSI[exon,]
  totalcount=totalRC[exon,]
  # only take those individuals with genotypes from the count table, and sort the phenotype to be the same order as in the genotype matrix
  psi=exonpsi[,IDs.common]
  count=totalcount[,IDs.common]
  
  #5. make the plot
  setwd(singleoutputpath)
  # This is the same type of plot used in the Glimmps paper. 
  # It is a combination of a boxplot with scatterplot. The 
  # x-axis represents the genotypes and y-axis represents 
  # PSI values.
  allele0=strsplit(exongenoname,split="_")[[1]][3]     #reference allele, the value is 0
  allele1=strsplit(exongenoname,split="_")[[1]][4]     #alternative allele, the value is 1
  #0/0 -> 1, 0/1 and 1/0 -> 1, 1/1 -> 2
  Alleles<-c(paste(allele0,allele0,sep="/"),paste(allele0,allele1,sep="/"),paste(allele1,allele1,sep="/"))
  
  n<- as.numeric(as.matrix(count))
  SNP<- exongeno
  Psi <- as.numeric(as.matrix(psi))
  N <- length(n)  #number of samples
  
  outfile = paste("Glimmpse_plot_",currentpair,"_",brainregion,"_",exon,"_",PSItype,".pdf",sep="")
  if (nchar(outfile)<255){     #if the file name is longer than 255 character, the file cannot be generated and also no need to generate since simple SNP won't have a name longer than 255 character
    pdf(outfile,height=6,width=6)
    
    title=paste(formalbrainregion,
                paste("Gene: ",genesymbol,sep=""),
                paste("Exon coordinate: ",coordinate,sep=""),
                paste("sQTL SNP: ",exongenoname,sep=""),
                sep="\n")
    ylim.range <- c(0,1) #range(psi,na.rm=T)
    plot(jitter(SNP,factor=0.5), Psi,xlim=c(-0.25,2.25), ylab="",xlab="",xaxt="n",type="n" ,ylim=ylim.range, cex.main=1, main=title)
    points(jitter(SNP,factor=0.5)  ,Psi  , pch= 19, cex= log10(n+1) ,col=1)
    mtext(text=Alleles, side=1, at= c(0,1,2),cex=1.5,line= 1)

    par("new"=T) # add boxplot on top
    boxplot(Psi~SNP, ylab="", xlab="", xaxt="n", yaxt="n",boxwex=0.35, xlim=c(-0.25,2.25), ylim=ylim.range,at=sort(unique(SNP[!is.na(SNP)])) ,border=gray(0.45),col=rgb(1,1,1,alpha=0.6),outline=FALSE)
    dev.off()
  }else{
    print(outfile)
  }
}

  







