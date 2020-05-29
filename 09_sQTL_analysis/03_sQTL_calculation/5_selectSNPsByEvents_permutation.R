#The purpose of this code is to select (significant) top SNP (with smallest p value. 
#if there are ties, we choose the one that is closest to the splice site)
# for each SE and output sQTL list

args <- commandArgs(TRUE)
FDRcutoff=0.1
numperm=5   #number of permutations

###1. read in exon information###
exoninfo=read.table(args[2],sep="\t",header=T)
rownames(exoninfo)=exoninfo[,"ID"]

###2. for each exon, we find the SNP with the smallest p value###
fin_list= Sys.glob(args[1])   #file name of Glimmpse result for each exon
topsnpinfo=matrix(NA,length(fin_list),8)
for (i in 1:length(fin_list)){        #for each exon
  file=fin_list[i]
  se_id=strsplit(strsplit(file,split="/")[[1]][2],split="\\.")[[1]][1]
  topsnpinfo[i,1]=se_id
  glimmpse=read.table(file,sep="\t",header=TRUE)
  
  #calculate distance between all the SNPs and the target exon
  dis=-(glimmpse[,"Pos"]-((exoninfo[se_id,6]+exoninfo[se_id,7])/2))      #middle of exon body (we choose the reverse because this agrees with our python code, i.e., the distance=middle of exon-SNP position rather than SNP position-middle of exon)
  absdis=abs(dis)
  glimmpse=cbind(glimmpse,absdis,dis)
  #find the SNP with the smallest p value (if ties, find the one with smallest absdis)
  topsnp=as.character(as.matrix(glimmpse[order(glimmpse[,"pvals.lm"],glimmpse[,"absdis"],decreasing=FALSE),][1,]))
  topsnpinfo[i,2:8]=topsnp
}
#output the top SNP information for all the exons
fout=paste('selected.SNPs.glmm.',strsplit(args[1],split="/")[[1]][1],".permutation.txt",sep="")    #top SNPs for each exon
write.table(topsnpinfo,fout,sep="\t",row.names=F,col.names=F,quote=F)

###3. for each exon, read in the permuted p value of the corresponding snp###
permpvalue=matrix(NA,length(fin_list),(2+numperm))
permpvalue[,1]=topsnpinfo[,1]    #exon ID
permpvalue[,2]=topsnpinfo[,3]    #snp ID   
for (i in 1:length(fin_list)){    #for each exon
  file=fin_list[i]
  se_id=strsplit(strsplit(file,split="/")[[1]][2],split="\\.")[[1]][1]
  chr=topsnpinfo[i,2]
  #read in the permuted p value for the current exon-SNP pair
  for (j in 1:numperm){
    newfile=paste(paste(strsplit(file,split="/")[[1]][1],"_perm",j,sep=""),paste("chr",chr,sep=""),strsplit(file,split="/")[[1]][2],sep="/")
    permglimmpse=read.table(newfile,sep="\t",header=TRUE)
    
    permpvalue[i,(j+2)]=min(as.numeric(as.character(permglimmpse[,"pvals.lm"])),na.rm=T)
    #the permuted p value for each exon in each permutation is still the smallest p value across all SNPs (NOT the p value of the exon-snp pair in the observed data)
  }
}

###4. calculate empirical FDR and output significant sQTLs###
observed=sort(as.numeric(topsnpinfo[,5]))       #observed p values in ascending order
permuted=sort(as.numeric(permpvalue[,(2+1:numperm)]))   #permuted p values in ascending order

#calculate FDR
pvaluelist=unique(observed)      #we use observed p values as p value cutoff (permuted p values are not used here)
fdr=rep(NA,length(pvaluelist))                   #fdr corresponding to each p value cutoff
for (i in 1:length(pvaluelist)){
  pcut=pvaluelist[i]                            #for each p value cutoff
  FP=sum(permuted<=pcut)        #for each cutoff, how many false positives (permuted p values smaller than cutoff)
  FPTP=sum(observed<=pcut)      #for each cutoff, how many true positives + false positive (observed p values smaller than cutoff)
  fdr[i]=(FP/numperm)/FPTP      #let's say we have M exons, then total number of observed p values is M and total number of permuted p values is N*M so we need to devide FP by N to normalize the ratio
}
fdr[which(fdr>1)]=1
#fdr greater than 1 means the top permuted p values are actually smaller than the top observed ones (in normal cases, FPTP should be greater than FP but if the null distribution is on the left of alternative distribution, then it basically means everything is false positive)

#get the sQTLs 
sig=which(fdr<=FDRcutoff)
pos=sig[length(sig)]       #this is the position of the last p value that gives a fdr smaller than the cutoff
if (pos>=1){    #if we have at least one sQTL with FDR < cutoff
  pvaluecutoff=pvaluelist[pos]     #get the corresponding p value
  print(pvaluecutoff)
  write.table(pvaluecutoff,"permutation_p_value_cutoff.txt",row.names=F,col.names=F,quote=F)
  sigsQTLinfo=topsnpinfo[which(as.numeric(topsnpinfo[,5])<=pvaluecutoff),]       #get the sQTLs that pass the permutation test
  sigsQTLinfo=sigsQTLinfo[order(as.numeric(sigsQTLinfo[,5]),decreasing=F),]      #order them based on their p values in ascending order
  
  #output information of those sQTLs
  fout_sqtl=paste('selected.sQTL.glmm.',strsplit(args[1],split="/")[[1]][1],".permutation.txt",sep="")      #significant top SNPs (i.e., sQTL) for each exon
  write.table(sigsQTLinfo,fout_sqtl,sep="\t",row.names=F,col.names=F,quote=F)
}else{
  print("no sQTL passes permutation test")       #this line doesn't work for some reason
}

