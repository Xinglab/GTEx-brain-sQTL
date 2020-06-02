#!/usr/bin/env Rscript
# # Note: make sure ID from all input files are vectors of character

#Usage: Rscript sQTLregress.oneexon.cis.R CEU_exon_info.txt CEU_DP.txt CEU_IC.txt CEU_allsnp chr22 1 10 > out.txt
#Usage: Rscript sQTLregress.oneexon.cis.R exon_info.fromGTF.SE.txt GTEx_76_JC.raw.input.SE.txt.DP.txt GTEx_76_JC.raw.input.SE.txt.IC.txt chr21.Genotype_450_maf0.01 chr21 1 10 > out.txt

#############
# sQTL analysis for SNPs within 200kb window of one single exon


## Arguments needed for program listed at the command line
# exoninforfile 
# alljunctionfile
# IJfile
# genoplinkfile
#chromosome   "chr#"
# exonindex
# exonindex_end

library("data.table")
args = commandArgs(TRUE)

exoninforfile = args[1]
alljunctionfile = args[2]
IJfile= args[3]
genoplinkfile= args[4]
chrom=args[5]
exonindex = as.numeric(args[6])
exonindex_end = as.numeric(args[7])
phenofile.IJ = IJfile
phenofile = alljunctionfile
genofile = genoplinkfile
studyname = args[8]
genotypename = args[9]
PSItype = args[10]
print (c(exonindex, exonindex_end))

################################
#  read in the exon information

exon.infor <- read.table(exoninforfile,header=T,sep="\t",stringsAsFactors=F) # Note: notice sep
exonIDs <- exon.infor[,1] 
chrstr <- exon.infor[, 4]
chr<- sub("chr","",chrstr)
SNPstartpos=exon.infor[,6]-200000
SNPendpos=exon.infor[,7]+200000
exonnames <- exonIDs 

# select the right SE (exons in the current chromosome)
selectchr=exon.infor[chrstr==chrom,]
if (exonindex_end > nrow(selectchr)) exonindex_end=nrow(selectchr)       #get the range of exons in the current batch
selected.se=as.vector(as.matrix(selectchr[exonindex:exonindex_end,1]))   #get the ID of selected exons
selected.index=c(1:nrow(exon.infor))[exonIDs %in% selected.se]        #get the row number of selected exons

TMPASSODIR <- paste0("Glimmps_each_exon_cis_", studyname)
library(lme4) 

source(paste(c(strsplit(genotypename,split="/")[[1]][1:11],"scripts/GLiMMPS_functions.R"),collapse="/"))

###############################################################################################################################################
# read in phenotype expression levels in plink format

a <- fread( phenofile ,header=T,sep="\t") # Note: check sep      (total read count)
#a <- fread( "GTEx_76_JC.raw.input.SE.txt.DP.txt" ,header=T,sep="\t") # Note: check sep
print (dim(a))
allreads.data <- data.frame(a)
nsample <- dim(allreads.data)[1] 
nexons <-dim(allreads.data)[2] -2

allreads.matrix <- allreads.data[,seq(1,nexons)+2 ]    #remove the first two columns
IDs.pheno <- as.character(allreads.data[,1])     									# I made change here

b <- fread( phenofile.IJ ,header=T,sep="\t") # Note: check sep
print (dim(b))
#b <- fread( "GTEx_76_JC.raw.input.SE.txt.IC.txt" ,header=T,sep="\t") # Note: check sep
IJ.data <- data.frame(b)
IJ.matrix <- IJ.data[,seq(1,nexons)+2 ]

### create output directory ##
if (!file.exists( TMPASSODIR )) {system(paste ("mkdir",TMPASSODIR)) }

##############
# read in genotype information

#.map file for the current chromosome
#c <-fread(paste("/u/scratch2/p/panyang/sQTL/genotype/",genotypename,"/",chrom,".",genofile, ".map", sep=""), header=F, sep="\t")    # Note: check sep (Yang's code)
c <-fread(paste(genotypename,"/",chrom,".",genofile, ".map", sep=""), header=F, sep="\t")    # Note: check sep (My code)
map <- data.frame(c)
names(map) <- c("Chr","SNPID","cM","Pos")
print(dim(c))

#.raw file for the current chromosome
#d <- fread(paste("/u/scratch2/p/panyang/sQTL/genotype/",genotypename,"/",chrom, ".", genofile, ".map_tpose.raw", sep=""), header=T,sep="\t",  na.strings=c("NA"))     # (Yang's code)
d <- fread(paste(genotypename,"/",chrom, ".", genofile, ".map_tpose.raw", sep=""), header=T,sep="\t",  na.strings=c("NA"))     # (My code)
#d <- fread("chr22.CEU_allsnp_tpose.raw", header=T,sep="\t",na.strings=c("NA")) # Note: sep	
#d <- fread("chr21.Genotype_450_maf0.01.map_tpose.raw", header=T,sep="\t",na.strings=c("NA"))
genoplink <- data.frame(d)
genoplink=genoplink[,-1] 
print(dim(d))

IDs.geno <- names(genoplink)           #  I changed here 

IDs.common <- IDs.geno[is.element(IDs.geno, IDs.pheno)]     #get the samples with genotype information
nsnps <- dim(genoplink)[1] -5                              #  I changed here 
sub.geno <- match(IDs.common , IDs.geno)          #the position of samples in IDs.geno that are also in IDs.common

geno <- as.matrix(genoplink[seq(1, nsnps)+5, sub.geno])    #  I changed here (we basically removed the header and first column from d)

### only take those individuals with genotypes from the count table, and sort the phenotype to be the same order as in the genotype matrix
sub <- match(IDs.common , IDs.pheno)

###############
#  start association for psi~ SNP for every SNP. 


for (i in selected.index){       #for each selected exon
  if (file.exists( paste(TMPASSODIR,"/",exonIDs[i],".asso",sep=""))) next;
  n<- allreads.matrix[sub,i]
  y <- IJ.matrix[sub,i]
  
  snp.filter=map[,1]==chr[i] & map[,4]>SNPstartpos[i] & map[,4]<SNPendpos[i]   #get the SNPs that are on the same chromosome with the exons & within the 200kb range
  if (length(snp.filter[snp.filter])==0){print (paste("nocis", exonIDs[i], sep=" ")); next;}   #if there is no snp passes the filter
  print (paste("start", exonIDs[i], sep=" "))       #if there are snps pass the filter
  #pvals.glmm <- rep(NA,nsnps)
  pvals.lm <- rep(NA,nsnps)
  betas <-  rep(NA,nsnps) 
  for ( gi in  c(1: nsnps)[snp.filter]) {     #for each snp that passes the filter   
    ### data association
    
    #########################
    #missing genotype filter#
    #########################
    #originalSNP=as.vector(geno[gi,])
    #validSNP=originalSNP[which(!is.na(y/n))]        #get the SNP value of samples with corresponding PSI value not equal to NA (if the PSI is NA, then it still doesn't help to have non-missing genotype information)
    ##we want to filter out SNPs with only one non-missing value in any genotype group
    #geno0=length(which(validSNP==0))    #number of samples with genotype 0
    #geno1=length(which(validSNP==1))
    #geno2=length(which(validSNP==2))
    #samplesize=c(geno0,geno1,geno2)
    #samplesize=samplesize[which(samplesize>0)]        
    #if (min(samplesize)==1){         #the sample size can be 0 but if it is not 0, we require it to be greater than 1
    #  next      #if the SNP doesn't meet this requirement, we will skip this SNP
    #}
    
    onedata <- list(n=n,y=y,SNP=as.vector(geno[gi,]) )
    #results.glm <- glm.sQTL ( onedata )
    
    #pvals.glm[gi] <- results.glm$pval
    
    #results.quasi <- glmquasi.sQTL ( onedata )
    #pvals.glmquasi[gi] <- results.quasi$pval
    
    ###############
    # GLiMMPS method 
    #print (onedata)
    #results.glmm <- glmm.sQTL ( onedata )
    #pvals.glmm[gi] <- results.glmm$pval
    #betas[gi] <- results.glmm$betas[2]
    ############
    results.lm <- lm.sQTL ( onedata , psitype=PSItype)
    pvals.lm[gi] <- results.lm$pval
    betas[gi] <- results.lm$betas[2]
    #results.glmmWald <- glmmWald.sQTL ( onedata )
    #pvals.glmmWald[gi] <- results.glmmWald$pval
    
  }
  
  tmpout <- cbind(map[snp.filter,c(1,2,4)], formatC(pvals.lm[snp.filter],format="e",digits=3), round(betas[snp.filter],3))
  
  colnames(tmpout) <- c("Chr","SNPID","Pos","pvals.lm","Beta")
  
  
  
  #########################
  #missing genotype filter#
  #########################
  ##remove rows with NA p value and NA beta (those are SNPs with low sample size and they should be filtered out)
  #missingrow=intersect(which(gsub(" ", "", as.character(tmpout[,"pvals.lm"]), fixed = TRUE)=="NA"),which(is.na(gsub(" ", "", as.character(tmpout[,"Beta"]), fixed = TRUE))))
  #if (length(missingrow)>0){
  #  tmpout=tmpout[-missingrow,]
  #}
  
  
  
  #remove rows with NaN p value and 0 beta (those are sQTL events with same PSI in different genotype groups)
  flatrow=intersect(which(gsub(" ", "", as.character(tmpout[,"pvals.lm"]), fixed = TRUE)=="NaN"),which(gsub(" ", "", as.character(tmpout[,"Beta"]), fixed = TRUE)==0))
  if (length(flatrow)>0){
    tmpout=tmpout[-flatrow,]
  }
  
  write.table(tmpout  , paste(TMPASSODIR,"/",exonIDs[i],".asso",sep=""),row.names=F,col.names=T, quote=F,sep="\t" )
  
}
######################
cat("Finished association!\n")

