#The purpose of this code is for each sQTL exon, we go through each brain region and find the SNP with the smallest p value
#then we collect that p value, beta value, location of that SNP and ID of that SNP

outputpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/12_region_ubiquitous_vs_specific_splicing_event_heatmap_radar_plot"
splicetype="SE"
PSItype="logit"
type="pvalue"
counttype="JC"
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

findpos=function(exoninfo,shortexonID,snppos,splicetype){         #classify SNPs into five categories
  oneexoninfo=as.character(as.matrix(exoninfo[shortexonID,]))
  strand=oneexoninfo[5]
  exonstart=as.numeric(oneexoninfo[6])
  exonend=as.numeric(oneexoninfo[7])
  SNPpos=as.numeric(snppos)
  
  if (splicetype=="SE"){
    if (strand=="+"){
      #wide splice site
      SS5start=exonend-3
      SS5end=exonend+6
      SS3start=exonstart-20
      SS3end=exonstart+3
      #narrow splice site
      diSS5start=exonend+1
      diSS5end=exonend+2
      diSS3start=exonstart-1
      diSS3end=exonstart
    }
    if (strand=="-"){
      #wide splice site
      SS5start=exonstart-6
      SS5end=exonstart+3
      SS3start=exonend-3
      SS3end=exonend+20
      #narrow splice site
      diSS5start=exonstart-1
      diSS5end=exonstart
      diSS3start=exonend+1
      diSS3end=exonend+2
    }
    newexonstart=exonstart
    newexonend=exonend
    
    dis=min(abs(SNPpos-SS5start),abs(SNPpos-SS5end),abs(SNPpos-SS3start),abs(SNPpos-SS3end)) 
    
    if (SNPpos>=diSS5start && SNPpos<=diSS5end){
      label="dinucleotide"    #5'SS dinucleotide
      return(label)
    }else if (SNPpos>=diSS3start && SNPpos<=diSS3end){
      label="dinucleotide"    #3'SS dinucleotide
      return(label)
    }else if (SNPpos>=SS5start && SNPpos<=SS5end){
      label="SS"      #5'SS
      return(label)
    }else if (SNPpos>=SS3start && SNPpos<=SS3end){
      label="SS"      #3'SS
      return(label)
    }else if (SNPpos>=newexonstart && SNPpos<=newexonend){
      label="exon"
      return(label)
    }else if (dis<=300){
      label="<=300bp"
      return(label)
    }else{
      label=">300bp"
      return(label)
    }
  }
  
  if (splicetype=="A3SS"){
    if (strand=="+"){
      ass1=as.numeric(oneexoninfo[6])
      ass1_start=ass1-20
      ass1_end=ass1+3
      ass2=as.numeric(oneexoninfo[8])
      ass2_start=ass2-20
      ass2_end=ass2+3
      #narrow splice site
      diass1start=ass1-1
      diass1end=ass1
      diass2start=ass2-1
      diass2end=ass2
    }
    if (strand=="-"){
      ass1=as.numeric(oneexoninfo[7])
      ass1_start=ass1-3
      ass1_end=ass1+20
      ass2=as.numeric(oneexoninfo[9])
      ass2_start=ass2-3
      ass2_end=ass2+20
      #narrow splice site
      diass1start=ass1+1
      diass1end=ass1+2
      diass2start=ass2+1
      diass2end=ass2+2
    }
    
    dis=min(abs(SNPpos-ass1_start),abs(SNPpos-ass1_end),abs(SNPpos-ass2_start),abs(SNPpos-ass2_end)) 
    
    if (SNPpos>=diass1start && SNPpos<=diass1end){
      label="dinucleotide"
      return(label)
    }else if (SNPpos>=diass2start && SNPpos<=diass2end){
      label="dinucleotide"
      return(label)
    }else if (SNPpos>=ass1_start && SNPpos<=ass1_end){ 
      label="SS"
      return(label)
    }else if (SNPpos>=ass2_start && SNPpos<=ass2_end){
      label="SS"
      return(label)
    }else if (SNPpos>=exonstart && SNPpos<=exonend){
      label="exon"
      return(label)
    }else if (dis<=300){
      label="<=300bp"
      return(label)
    }else{
      label=">300bp"
      return(label)
    }
  }
  
  if (splicetype=="A5SS"){
    if (strand=="+"){
      ass1=as.numeric(oneexoninfo[7])
      ass1_start=ass1-3
      ass1_end=ass1+6
      ass2=as.numeric(oneexoninfo[9])
      ass2_start=ass2-3
      ass2_end=ass2+6
      #narrow splice site
      diass1start=ass1+1
      diass1end=ass1+2
      diass2start=ass2+1
      diass2end=ass2+2
    }
    if (strand=="-"){
      ass1=as.numeric(oneexoninfo[6])
      ass1_start=ass1-6
      ass1_end=ass1+3
      ass2=as.numeric(oneexoninfo[8])
      ass2_start=ass2-6
      ass2_end=ass2+3
      #narrow splice site
      diass1start=ass1-1
      diass1end=ass1
      diass2start=ass2-1
      diass2end=ass2
    }
    
    dis=min(abs(SNPpos-ass1_start),abs(SNPpos-ass1_end),abs(SNPpos-ass2_start),abs(SNPpos-ass2_end)) 
    
    if (SNPpos>=diass1start && SNPpos<=diass1end){
      label="dinucleotide"
      return(label)
    }else if (SNPpos>=diass2start && SNPpos<=diass2end){
      label="dinucleotide"
      return(label)
    }else if (SNPpos>=ass1_start && SNPpos<=ass1_end){ 
      label="SS"
      return(label)
    }else if (SNPpos>=ass2_start && SNPpos<=ass2_end){
      label="SS"
      return(label)
    }else if (SNPpos>=exonstart && SNPpos<=exonend){
      label="exon"
      return(label)
    }else if (dis<=300){
      label="<=100bp"
      return(label)
    }else{
      label=">300bp"
      return(label)
    }
  }
}

pick_category=function(list_of_category,order_of_importance){      
  #given a list of categories and the order we want to use, this function will return one category from the list based on the order of importance
  reorderedlist=list_of_category[order(match(list_of_category, order_of_importance))]
  return(reorderedlist)
}

importance_order=c("dinucleotide", "SS", "exon", "<=300bp", ">300bp")

#filter out SNPs which are not single base 
SNPfilter=function(x){
  allele0=strsplit(x,split="_")[[1]][3]     #reference allele, the value is 0
  allele1=strsplit(x,split="_")[[1]][4]     #alternative allele, the value is 1
  if (nchar(allele0)==1 && nchar(allele1)==1){
    return("Single")
  }else{
    return("Not_single")
  }
}


getrsID=function(x,rsIDtable){
  shortsnpID=paste(strsplit(x,split="_")[[1]][1],strsplit(x,split="_")[[1]][2],sep="_")
  if (shortsnpID %in% rsIDtable[,1]){
    rsID=as.character(rsIDtable[which(rsIDtable[,1] %in% shortsnpID)[1],2])
  }else{
    rsID=x
  }
  return(rsID)
}


#############################
#Collect all the information#
#############################
#read in the exon information table
exoninfopath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/input_splicing",
                   PSItype,counttype,splicetype,sep="/")
setwd(exoninfopath)
exoninfo=read.table(paste("exon_info.fromGTF.",splicetype,".txt",sep=""),sep="\t",header=T)
rownames(exoninfo)=exoninfo[,"ID"]

#read in the lookup table
setwd(outputpath)
lookup=read.table("Genotype_swap_lookup_table.txt",sep="\t",header=T)

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


#generate matrix for the information we need, i.e., p value, beta, SNP ID, rsID, location of SNP
betamatrix300=betamatrix200k=pvmatrix300=pvmatrix200k=lcmatrix300=lcmatrix200k=topsnpmatrix300=topsnpmatrix200k=rsmatrix300=rsmatrix200k=matrix(NA,length(sQTLexon),length(brainregionlist))
rownames(betamatrix300)=rownames(betamatrix200k)=rownames(pvmatrix300)=rownames(pvmatrix200k)=rownames(lcmatrix300)=rownames(lcmatrix200k)=rownames(topsnpmatrix300)=rownames(topsnpmatrix200k)=rownames(rsmatrix300)=rownames(rsmatrix200k)=sQTLexon
colnames(betamatrix300)=colnames(betamatrix200k)=colnames(pvmatrix300)=colnames(pvmatrix200k)=colnames(lcmatrix300)=colnames(lcmatrix200k)=colnames(topsnpmatrix300)=colnames(topsnpmatrix200k)=colnames(rsmatrix300)=colnames(rsmatrix200k)=brainregionlist

#fill in the five matrix
for (i in 1:length(sQTLexon)){
  exonfullID=sQTLexon[i]
  exonshortID=exonIDconversion[which(exonIDconversion[,"fullID"]==exonfullID),"shortID"]
  print(i)
  
  #go to each brain region
  for (br in 1:length(brainregionlist)){
    currentBR=brainregionlist[br]
    
    #read in the rsID conversion table
    inputrspath=paste(sqtlrootinput,currentBR,sep="/")
    setwd(inputrspath)
    rstable=read.table(paste("selected.sQTL.glmm.Glimmps_each_exon_cis_",currentBR,".pvalue.txt.rsIDmap.txt",sep=""),sep="\t")
    #rownames(rstable)=rstable[,1]
    
    #read in the exon result in the current brain region
    inputpath=paste(sqtlrootinput,currentBR,paste("Glimmps_each_exon_cis_",currentBR,sep=""),sep="/")
    setwd(inputpath)
    exonresult=read.table(paste(exonshortID,".asso",sep=""),sep="\t",header=T)
    
    #remove SNPs that are not single base
    singlebase=sapply(as.matrix(exonresult[,"SNPID"]),SNPfilter)
    exonresult=exonresult[which(singlebase=="Single"),]
    
    #calculate the location of each SNP
    location=sapply(as.matrix(exonresult[,"Pos"]), findpos, exoninfo=exoninfo, shortexonID=exonshortID, splicetype=splicetype)
    exonresult=cbind(exonresult,location)
    
    #######################################
    #for all the SNPs (within 200kb range)#
    #######################################
    #find the row with the smallest p value, largest absolute beta value, and closest distance to the splice site
    exonresult2sort=exonresult
    exonresult2sort[,"Beta"]=abs(exonresult2sort[,"Beta"])
    exonresult2sort[,"location"]=as.character(exonresult2sort[,"location"])
    exonresult2sort[which(exonresult2sort[,"location"]=="dinucleotide"),"location"]=1
    exonresult2sort[which(exonresult2sort[,"location"]=="SS"),"location"]=2
    exonresult2sort[which(exonresult2sort[,"location"]=="exon"),"location"]=3
    exonresult2sort[which(exonresult2sort[,"location"]=="<=300bp"),"location"]=4
    exonresult2sort[which(exonresult2sort[,"location"]==">300bp"),"location"]=5
    #exonresult2sort=exonresult2sort[order(exonresult2sort[,"pvals.lm"],exonresult2sort[,"Beta"],exonresult2sort[,"location"],decreasing=FALSE),]
    exonresult2sort=exonresult2sort[with(exonresult2sort, order(pvals.lm,-Beta,location)), ]    #p value in increasing order, beta in decreasing order (when tied p value), location in increasing order (when tied both p value and beta)
    
    row=which(exonresult[,"SNPID"]==exonresult2sort[1,"SNPID"])
    pvalue=exonresult[row,"pvals.lm"]
    beta=exonresult[row,"Beta"]
    snploci=as.character(exonresult[row,"location"])
    snpID=as.character(exonresult[row,"SNPID"])
    rsID=getrsID(snpID,rstable)
    
    pvmatrix200k[i,br]=pvalue
    lcmatrix200k[i,br]=snploci
    topsnpmatrix200k[i,br]=snpID
    rsmatrix200k[i,br]=rsID
    #check if the SNP genotype is swapped or not
    if (as.character(lookup[snpID,"label"])=="Correct"){
      betamatrix200k[i,br]=beta
    }else{
      betamatrix200k[i,br]=-beta
    }
    
    #row=which(exonresult[,"pvals.lm"]==min(exonresult[,"pvals.lm"]))       #the SNP(s) with the smallest p value
    #if (length(row>1)){            #if we have multiple SNPs with the same p value, we choose the one with the largest effect size
    #  rowbeta=which(abs(exonresult300[row,"Beta"])==max(abs(exonresult300[row,"Beta"])))
    #  if (length(rowbeta)>1){      #if we have multiple SNPs with the same p value and same beta, we choose the one with the 
    #    
    #  }
      
    #  subrow=exonresult[row,][which(abs(exonresult[row,"Beta"])==max(abs(exonresult[row,"Beta"]))),]
    #  pvalue=subrow[,"pvals.lm"]
    #  beta=subrow[,"Beta"]
    #  snploci=as.character(subrow[,"location"])
    #  snpID=as.character(subrow[,"SNPID"])
    #  rsID=getrsID(snpID,rstable)
    #}else{
    #}
    
    #############################
    #for SNPs within 300bp range#
    #############################
    exonresult300=subset(exonresult,exonresult[,"location"]!=">300bp")
    if (dim(exonresult300)[1]>0){           #if we have SNPs within 300bp
      exonresult2sort300=subset(exonresult2sort,exonresult2sort[,"location"]!=5)
      
      row=which(exonresult300[,"SNPID"]==exonresult2sort300[1,"SNPID"])
      pvalue=exonresult300[row,"pvals.lm"]
      beta=exonresult300[row,"Beta"]
      snploci=as.character(exonresult300[row,"location"])
      snpID=as.character(exonresult300[row,"SNPID"])
      rsID=getrsID(snpID,rstable)
      
      pvmatrix300[i,br]=pvalue
      lcmatrix300[i,br]=snploci
      topsnpmatrix300[i,br]=snpID
      rsmatrix300[i,br]=rsID
      #check if the SNP genotype is swapped or not
      if (as.character(lookup[snpID,"label"])=="Correct"){
        betamatrix300[i,br]=beta
      }else{
        betamatrix300[i,br]=-beta
      }
    }
  }
}

outputmatrix300bp=paste(outputpath,"300bp_highp",sep="/")
command=paste("mkdir -p",outputmatrix300bp)
system(command)
setwd(outputmatrix300bp)
write.table(pvmatrix300,"pvmatrix300_highp.txt",sep="\t")
write.table(lcmatrix300,"lcmatrix300_highp.txt",sep="\t")
write.table(topsnpmatrix300,"topsnpmatrix300_highp.txt",sep="\t")
write.table(rsmatrix300,"rsmatrix300_highp.txt",sep="\t")
write.table(betamatrix300,"betamatrix300_highp.txt",sep="\t")

outputmatrix200kb=paste(outputpath,"200kb_highp",sep="/")
command=paste("mkdir -p",outputmatrix200kb)
system(command)
setwd(outputmatrix200kb)
write.table(pvmatrix200k,"pvmatrix200k_highp.txt",sep="\t")
write.table(lcmatrix200k,"lcmatrix200k_highp.txt",sep="\t")
write.table(topsnpmatrix200k,"topsnpmatrix200k_highp.txt",sep="\t")
write.table(rsmatrix200k,"rsmatrix200k_highp.txt",sep="\t")
write.table(betamatrix200k,"betamatrix200k_highp.txt",sep="\t")



#for each sQTL exon, calculate the number of significant brain regions, add that as an additional column
if (type=="pvalue"){
  cutoff=10^-5
}

num.sig.region300=num.sig.region200kb=dis_2_SS300=dis_2_SS200k=rep(NA,length(sQTLexon))
for (i in 1:length(sQTLexon)){
  num.sig.region300[i]=sum(pvmatrix300[i,]<=cutoff)
  num.sig.region200kb[i]=sum(pvmatrix200k[i,]<=cutoff)
  dis_2_SS300[i]=pick_category(lcmatrix300[i,],importance_order)[1]
  dis_2_SS200k[i]=pick_category(lcmatrix200k[i,],importance_order)[1]
}

setwd(outputmatrix300bp)
write.table(cbind(num.sig.region300,dis_2_SS300),"numsigregion.dis2SS.300bp_highp.txt",sep="\t")
setwd(outputmatrix200kb)
write.table(cbind(num.sig.region200kb,dis_2_SS200k),"numsigregion.dis2SS.200kb_highp.txt",sep="\t")






