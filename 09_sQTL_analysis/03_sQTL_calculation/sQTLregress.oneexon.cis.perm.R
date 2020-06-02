#!/usr/bin/env Rscript

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
# permutation id

library("data.table")
args = commandArgs(TRUE)

exoninforfile = args[1]
alljunctionfile = args[2]
IJfile= args[3]
genoplinkfile= args[4]
chrom=args[5]
exonindex = as.numeric(args[6])
exonindex_end = as.numeric(args[7])
perm_id=as.numeric(args[8])
phenofile.IJ = IJfile
phenofile = alljunctionfile
genofile = genoplinkfile
studyname = args[9]
genotypename = args[10]
PSItype = args[11]
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

# select the right se
selectchr=exon.infor[chrstr==chrom,]
if (exonindex_end > nrow(selectchr)) exonindex_end=nrow(selectchr)
selected.se=as.vector(as.matrix(selectchr[exonindex:exonindex_end,1]))
selected.index=c(1:nrow(exon.infor))[exonIDs %in% selected.se]

permdir=paste("./Glimmps_each_exon_cis_",studyname,"_perm",perm_id, sep="")
TMPASSODIR <-paste("./Glimmps_each_exon_cis_",studyname,"_perm",perm_id,"/", chrom, sep="")

library(lme4) 

source(paste(c(strsplit(genotypename,split="/")[[1]][1:11],"scripts/GLiMMPS_functions.R"),collapse="/"))

###############################################################################################################################################
# read in phenotype expression levels in plink format
allreads.data <- data.frame(fread( phenofile ,header=T,sep="\t"))
nsample <- dim(allreads.data)[1] 
nexons <-dim(allreads.data)[2] -2

allreads.matrix <- allreads.data[,seq(1,nexons)+2 ]
IDs.pheno <- as.character(allreads.data[,1])      									# I made change here

IJ.data <-data.frame(fread( phenofile.IJ ,header=T,sep="\t"))
IJ.matrix <- IJ.data[,seq(1,nexons)+2 ]

### create output directory ##

if (!file.exists( permdir )) {system(paste ("mkdir", permdir)) }
if (!file.exists( TMPASSODIR )) {system(paste ("mkdir",TMPASSODIR)) }

  
# do the association for all SNPs near the target exon

#cmmd <- paste("mkdir ",TMPASSODIR,"/chr",chr,sep="")
#if (!file.exists( paste(TMPASSODIR,"/chr",chr,sep="")) ) {system(cmmd)}

##############
# read in genotype information

#map<-data.frame(fread(paste("/u/scratch2/p/panyang/sQTL/genotype/",genotypename,"/",chrom,".",genofile, ".map", sep=""), header=F, sep="\t"))       #  I changed here (Yang's code)
map<-data.frame(fread(paste(genotypename,"/",chrom,".",genofile, ".map", sep=""), header=F, sep="\t"))       # (My code)
names(map) <- c("Chr","SNPID","cM","Pos")
print(dim(map))

#genoplink <- data.frame(fread(paste("/u/scratch2/p/panyang/sQTL/genotype/",genotypename,"/","perm", perm_id, "/", chrom, ".", genofile,".map_tpose_perm", perm_id, ".raw", sep=""), header=T,sep="\t",  na.strings=c("NA")))     # (Yang's code) 			 
#genoplink <- data.frame(fread(paste(genotypename,"/", chrom, ".", genofile,".map_tpose_perm", perm_id, ".raw", sep=""), header=T,sep="\t",  na.strings=c("NA")))     # (My code, wrong permutation) 			 
genoplink <- data.frame(fread(paste(genotypename,"/",chrom, ".", genofile, ".map_tpose.raw", sep=""), header=T,sep="\t",  na.strings=c("NA")))     # (My code, new permutation method, the permutation happens below)

genoplink=genoplink[,-1] 

print (dim(genoplink))
 
IDs.geno <- names(genoplink)           #  I changed here 
IDs.common <- IDs.geno[is.element(IDs.geno, IDs.pheno)]
nsnps <- dim(genoplink)[1] -5                              #  I changed here 
sub.geno <- match(IDs.common , IDs.geno)
geno <- as.matrix(genoplink[seq(1, nsnps)+5, sub.geno])    #  I changed here 

### only take those individuals with genotypes, and sort the phenotype to be the same order as in the genotype matrix


sub <- match(IDs.common , IDs.pheno)
#IDs.pheno[sub]


###############
#  start association for psi~ SNP for every SNP. 


for (i in selected.index){

	#if (file.exists( paste(TMPASSODIR,"/",exonIDs[i],".asso",sep=""))) next;
	
    n<- allreads.matrix[sub,i]
    y <- IJ.matrix[sub,i]
    
 	snp.filter=map[,1]==chr[i] & map[,4]>SNPstartpos[i] & map[,4]<SNPendpos[i]
 	if (length(snp.filter[snp.filter])==0){print (paste("nocis", exonIDs[i], sep=" ")); next;}
 	print (paste("start", exonIDs[i], sep=" "))
 	#pvals.glmm <- rep(NA,nsnps)
 	pvals.lm <- rep(NA,nsnps)
  betas <-  rep(NA,nsnps) 
  for ( gi in  c(1: nsnps)[snp.filter]) {         
    ### data association
    
    #permutation
    set.seed(perm_id)        #set the seed so that we can reproduce the result
    originalSNP=as.vector(geno[gi,])    #get the original SNP value
    nonNASNP=originalSNP[which(!is.na(originalSNP))]   #for samples with non missing SNP values, get their SNP values
    perm_nonNASNP=sample(nonNASNP)      #permute only the non missing values
    permutedSNP=originalSNP
    permutedSNP[which(!is.na(permutedSNP))]=perm_nonNASNP     #give those permuted non missing values back to the same group of samples (so that missing values will still be missing and non missing values will still be non missing)
    

    #########################
    #missing genotype filter#
    #########################
    #validSNP=permutedSNP[which(!is.na(y/n))]        #get the SNP value of samples with corresponding PSI value not equal to NA (if the PSI is NA, then it still doesn't help to have non-missing genotype information)
    ##we want to filter out SNPs with only one non-missing value in any genotype group
    #geno0=length(which(validSNP==0))    #number of samples with genotype 0
    #geno1=length(which(validSNP==1))
    #geno2=length(which(validSNP==2))
    #samplesize=c(geno0,geno1,geno2)
    #samplesize=samplesize[which(samplesize>0)]        
    #if (min(samplesize)==1){         #the sample size can be 0 but if it is not 0, we require it to be greater than 1
    #  next      #if the SNP doesn't meet this requirement, we will skip this SNP
    #}
    
    
    #onedata <- list(n=n,y=y,SNP=  as.vector(geno[gi,]) )
    onedata <- list(n=n,y=y,SNP=permutedSNP)
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
  #remove rows with NA p value and NA beta (those are SNPs with low sample size and they should be filtered out)
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

