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

splicetype="SE"
PSItype="logit"
type="pvalue"
counttype="JC"
sqtlrootinput=paste("/path/to/sQTL_run/logit/JC",
                    splicetype,sep="/")
outputpath="/path/to/12_region_ubiquitous_vs_specific_splicing_event_heatmap_radar_plot/Region_ubiquitous_vs_specific_splicing_event_heatmap_radar_plot_step_3_pick_SNP_and_collect_values_for_heatmap_more_stringent_version/single_job_run/result/200kb"
rootinput="/path/to/12_region_ubiquitous_vs_specific_splicing_event_heatmap_radar_plot"

pick_category=function(list_of_category,order_of_importance){      
  #given a list of categories and the order we want to use, this function will return one category from the list based on the order of importance
  reorderedlist=list_of_category[order(match(list_of_category, order_of_importance))]
  return(reorderedlist)
}

importance_order=c("dinucleotide", "SS", "exon", "<=300bp", ">300bp")

collect.p.beta=function(exonshortID,snpID){     #function to get the p value and beta value of a given exon-snp pair in each brain region
  bpe=ppe=matrix(NA,1,length(brainregionlist))
  for (br in 1:length(brainregionlist)){
    currentBR=brainregionlist[br]
    inputpath=paste(sqtlrootinput,currentBR,paste("Glimmps_each_exon_cis_",currentBR,sep=""),sep="/")
    setwd(inputpath)
    exonresult=read.table(paste(exonshortID,".asso",sep=""),sep="\t",header=T)
    rownum=which(exonresult[,"SNPID"]==snpID)
    bpe[1,br]=exonresult[rownum,"Beta"]
    ppe[1,br]=exonresult[rownum,"pvals.lm"]
  }
  return(list(beta=bpe,pval=ppe))
}

#####################################################
#read in the p value/beta/rsID/snpID/location matrix#
#####################################################
setwd(paste(rootinput,"200kb_highp",sep="/"))
pvmatrix200k=read.table("pvmatrix200k_highp.txt",sep="\t",header=T,check.names=F)
lcmatrix200k=read.table("lcmatrix200k_highp.txt",sep="\t",header=T,check.names=F)
topsnpmatrix200k=read.table("topsnpmatrix200k_highp.txt",sep="\t",header=T,check.names=F)

betapickexon=pvaluepickexon=matrix(NA,1,length(brainregionlist))   #p value and beta value of the picked SNP for each exon

###################
#Other preparation#
###################
exoninfopath=paste("/path/to/input_splicing",
                   PSItype,counttype,splicetype,sep="/")
setwd(exoninfopath)
exoninfo=read.table(paste("exon_info.fromGTF.",splicetype,".txt",sep=""),sep="\t",header=T)
rownames(exoninfo)=exoninfo[,"ID"]

#read in the lookup table
setwd(rootinput)
lookup=read.table("Genotype_swap_lookup_table.txt",sep="\t",header=T)

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

#get the permutation p value cutoff for each brain region
cutoff=matrix(NA,length(brainregionlist),1)
rownames(cutoff)=brainregionlist
for (br in 1:length(brainregionlist)){
  currentbr=brainregionlist[br]
  setwd(paste(sqtlrootinput,"/",currentbr,sep=""))
  cutoff[br,1]=as.numeric(as.matrix(read.table("permutation_p_value_cutoff.txt",sep="\t")))
}

#output the result
output.p=paste(outputpath,"pvalue",sep="/")
command=paste("mkdir -p",output.p)
system(command)
output.b=paste(outputpath,"beta",sep="/")
command=paste("mkdir -p",output.b)
system(command)

#####################################################################################################
#for the current exon, pick the significant SNP(s) that are the top one in most of the brain regions#
#####################################################################################################
#exonnum <- as.numeric(Sys.getenv("SGE_TASK_ID"))
for (exonnum in 1:length(sQTLexon)){
  print(exonnum)
  if (sum(pvmatrix200k[exonnum,]<cutoff)>0){      
    sigsnplist=as.character(as.matrix(topsnpmatrix200k[exonnum,which(pvmatrix200k[exonnum,]<cutoff)]))   #get the list of significant SNPs (top SNP doesn't need to be significant)
    
    maxtopnum=max(table(sigsnplist))
    exon2pickstep1=names(which(table(sigsnplist)==maxtopnum))
    exonfullID=rownames(topsnpmatrix200k)[exonnum]
    exonshortID=exonIDconversion[which(exonIDconversion[,"fullID"]==exonfullID),"shortID"]
    
    if (length(exon2pickstep1)==1){              #if we only have one exon after the first round of selection
      #collect the p value and beta value of this exon
      result=collect.p.beta(exonshortID,exon2pickstep1)
      betapickexon=as.matrix(result$beta)
      pvaluepickexon=as.matrix(result$pval)
      rownames(betapickexon)=rownames(pvaluepickexon)=paste(exonfullID,exon2pickstep1,sep="~")
      
      #check if the SNP genotype is swapped or not
      if (as.character(lookup[exon2pickstep1,"label"])=="Correct"){
        betapickexon=betapickexon
      }else{
        betapickexon=-betapickexon
      }
      
      setwd(output.p)
      write.table(pvaluepickexon,paste(exonfullID,"~",exon2pickstep1,"_pvaluepickexon.txt",sep=""),sep="\t",col.names=F)
      
      setwd(output.b)
      write.table(betapickexon,paste(exonfullID,"~",exon2pickstep1,"_betapickexon.txt",sep=""),sep="\t",col.names=F)
    }else{
      ########################################################################################
      #If we have ties, find the significant SNP(s) that is significant in more brain regions#
      ########################################################################################
      signum=rep(NA,length(exon2pickstep1))     #the number of significant regions of each SNP
      for (e in 1:length(exon2pickstep1)){
        signum[e]=sum(collect.p.beta(exonshortID,exon2pickstep1[e])$pval<as.numeric(cutoff))
      }
      maxsignum=max(signum)
      exon2pickstep2=exon2pickstep1[which(signum==maxsignum)]
      
      if (length(exon2pickstep2)==1){                  #if we only have one exon after the second round of selection
        #collect the p value and beta value of this exon
        result=collect.p.beta(exonshortID,exon2pickstep2)
        betapickexon=as.matrix(result$beta)
        pvaluepickexon=as.matrix(result$pval)
        rownames(betapickexon)=rownames(pvaluepickexon)=paste(exonfullID,exon2pickstep2,sep="~")
        
        #check if the SNP genotype is swapped or not
        if (as.character(lookup[exon2pickstep2,"label"])=="Correct"){
          betapickexon=betapickexon
        }else{
          betapickexon=-betapickexon
        }
        
        setwd(output.p)
        write.table(pvaluepickexon,paste(exonfullID,"~",exon2pickstep2,"_pvaluepickexon.txt",sep=""),sep="\t",col.names=F)
        
        setwd(output.b)
        write.table(betapickexon,paste(exonfullID,"~",exon2pickstep2,"_betapickexon.txt",sep=""),sep="\t",col.names=F)
      }else{
        ##########################################################################
        #If we have ties, find the significant SNP(s) that is closest to the exon#
        ##########################################################################
        lclist=rep(NA,length(exon2pickstep2))
        for (e in 1:length(exon2pickstep2)){
          lclist[e]=as.character(as.matrix(lcmatrix200k[exonnum,which(topsnpmatrix200k[exonnum,]==exon2pickstep2[e])][1]))
        }
        closestlc=pick_category(lclist,importance_order)[1]
        exon2pickstep3=exon2pickstep2[which(lclist==closestlc)[1]]
        
        #collect the p value and beta value of this exon
        result=collect.p.beta(exonshortID,exon2pickstep3)
        betapickexon=as.matrix(result$beta)
        pvaluepickexon=as.matrix(result$pval)
        rownames(betapickexon)=rownames(pvaluepickexon)=paste(exonfullID,exon2pickstep3,sep="~")
        
        if (as.character(lookup[exon2pickstep3,"label"])=="Correct"){
          betapickexon=betapickexon
        }else{
          betapickexon=-betapickexon
        }
        
        setwd(output.p)
        write.table(pvaluepickexon,paste(exonfullID,"~",exon2pickstep3,"_pvaluepickexon.txt",sep=""),sep="\t",col.names=F)
        
        setwd(output.b)
        write.table(betapickexon,paste(exonfullID,"~",exon2pickstep3,"_betapickexon.txt",sep=""),sep="\t",col.names=F)
      }
    }
  }
}
  


