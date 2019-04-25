args <- commandArgs(TRUE)
br=as.character(args[1])      #brain region

setwd("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/scripts")
brainregion=as.character(as.matrix(read.table("brainregionlist.txt",sep="\t")))
rootinput="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/sQTL_run/logit/JC/"
rootsplicinginput='/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/input_splicing/logit/JC/'
outputpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.2_sQTL_SNP_annotation/result/fraction_of_exon_within_300bp_modified"


splicetypelist=c("SE","A3SS","A5SS")

for (splicetype in splicetypelist){
  ##########################
  #generate p value cutoffs#
  ##########################
  allpvalue=c()
  for (i in 1:length(brainregion)[1]){
    inputpath=paste(rootinput,splicetype,"/",brainregion[i],sep="")
    #read in the top SNP result of all the exons
    setwd(inputpath)
    exonSNP=read.table(paste("selected.SNPs.glmm.Glimmps_each_exon_cis_",brainregion[i],".permutation.txt",sep=""),sep="\t")    #the results based on pvalue and permutation are the same
    
    allpvalue=c(allpvalue,exonSNP[,5])
  }
  numpoint=100
  start=floor(quantile(-log10(allpvalue),c(0.01,0.99))[1])        #this is a more reasonable range since it still gives us enough sQTL events to calculate the proportion at every p value cutoff
  end=ceiling(quantile(-log10(allpvalue),c(0.01,0.99))[2])
  xaxis=seq(from=start,to=end,length.out=numpoint)
  pcut=10^(-xaxis)   #p value cutoffs
  
  #############################################################
  #get the p value and distance info of all exons and all SNPs#
  #############################################################
  checkdis=function(exoninfo,exonID,SNPpos,splicetype){
    strand=as.character(exoninfo[exonID,"strand"])
    if (splicetype=="SE"){
      exonstart=as.numeric(exoninfo[exonID,6])
      exonend=as.numeric(exoninfo[exonID,7])
      if (strand=="+"){
        SS5start=exonend-3
        SS5end=exonend+6
        SS3start=exonstart-20
        SS3end=exonstart+3
      }
      if (strand=="-"){
        SS5start=exonstart-6
        SS5end=exonstart+3
        SS3start=exonend-3
        SS3end=exonend+20
      }
    }
    if (splicetype=="SE"){
        dis=min(abs(SNPpos-SS5start),abs(SNPpos-SS5end),abs(SNPpos-SS3start),abs(SNPpos-SS3end))  
    }
    
    if (splicetype=="A5SS"){
      exonstart=as.numeric(exoninfo[exonID,6])     #we choose the longer exon
      exonend=as.numeric(exoninfo[exonID,7])
      if (strand=="+"){
        ass1=as.numeric(exoninfo[exonID,7])
        ass1_start=ass1-3
        ass1_end=ass1+6
        ass2=as.numeric(exoninfo[exonID,9])
        ass2_start=ass2-3
        ass2_end=ass2+6
      }
      if (strand=="-"){
        ass1=as.numeric(exoninfo[exonID,6])
        ass1_start=ass1-6
        ass1_end=ass1+3
        ass2=as.numeric(exoninfo[exonID,8])
        ass2_start=ass2-6
        ass2_end=ass2+3
      }
    }
    if (splicetype=="A3SS"){
      exonstart=as.numeric(exoninfo[exonID,6])    #we choose the longer exon
      exonend=as.numeric(exoninfo[exonID,7])
      if (strand=="+"){
        ass1=as.numeric(exoninfo[exonID,6])
        ass1_start=ass1-20
        ass1_end=ass1+3
        ass2=as.numeric(exoninfo[exonID,8])
        ass2_start=ass2-20
        ass2_end=ass2+3
      }
      if (strand=="-"){
        ass1=as.numeric(exoninfo[exonID,7])
        ass1_start=ass1-3
        ass1_end=ass1+20
        ass2=as.numeric(exoninfo[exonID,9])
        ass2_start=ass2-3
        ass2_end=ass2+20
      }
    }   
    if (splicetype=="A3SS" || splicetype=="A5SS"){
        dis=min(abs(SNPpos-ass1_start),abs(SNPpos-ass1_end),abs(SNPpos-ass2_start),abs(SNPpos-ass2_end)) 
    }
    return(dis)
  }
  
  
  result2plot=list()
  inputpath=paste(rootinput,splicetype,"/",br,sep="")
  #read in the top SNP result of all the exons
  setwd(inputpath)
  exonSNP=read.table(paste("selected.SNPs.glmm.Glimmps_each_exon_cis_",br,".permutation.txt",sep=""),sep="\t")    #the results based on pvalue and permutation are the same
  
  #read in exon information
  inputsplicing=paste(rootsplicinginput,splicetype,"/exon_info.fromGTF.",splicetype,".txt",sep="")
  exoninfo=read.table(inputsplicing,sep="\t",header=T)
  rownames(exoninfo)=exoninfo[,"ID"]
  
  #for each exon, we get all the SNPs information (position and p value)
  for (e in 1:dim(exonSNP)[1]){
    print(e)
    exon=as.character(exonSNP[e,1])
    result2plot[[exon]]=list()
    exonpath=paste(inputpath,"/","Glimmps_each_exon_cis_",br,sep="")
    setwd(exonpath)
    glimmpse=read.table(paste(exon,"asso",sep="."),sep="\t",header=TRUE)
    #store the p value info of each SNP
    result2plot[[exon]][["pvalue"]]=glimmpse[,"pvals.lm"]
    #for each SNP, calculate distance to the nearest SS
    distance=sapply(as.numeric(glimmpse[,"Pos"]),checkdis,exoninfo=exoninfo,exonID=exon,splicetype)
    result2plot[[exon]][["distance"]]=distance
  }
  
  #############################################
  #calculate proportion at each p value cutoff#
  #############################################
  proportion=rep(NA,numpoint)
  denominator=rep(NA,numpoint)
  numerator=rep(NA,numpoint)
  for (p in 1:length(pcut)){
    #total number of significant exon
    numexon=0  
    #calculate the number of exons with significant SNP within 300bp
    numexonwith300=0
    for (exon in names(result2plot)){       #for each exon
      sigdis=result2plot[[exon]]$distance[which(result2plot[[exon]]$pvalue<=pcut[p])]   #get the distance of all the significant SNPs
      if (length(sigdis)>0){      #if this exon has significant SNPs
        numexon=numexon+1
        if (min(sigdis)<=300){       #if there is at least one significant SNP with distance within 300bp
          numexonwith300=numexonwith300+1        #we consider this exon as an exon with SNP within 300bp
        }
      }
    }
    proportion[p]=numexonwith300/numexon
    denominator[p]=numexon
    numerator[p]=numexonwith300
  }
  proportionmatrix=cbind(cbind(cbind(proportion,xaxis),denominator),numerator)
  colnames(proportionmatrix)=c("proportion","xaxis","denominator","numerator")
  setwd(outputpath)
  write.table(proportionmatrix,paste(splicetype,"_",br,"_proportion.txt",sep=""),sep="\t")
}

