splicetypelist=c("SE","A3SS","A5SS")
typelist=c("pvalue","permutation")
rootinput="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/summary/logit/JC"
RBPdb="/u/nobackup/yxing/PROJECT/yidazhan/research/software/deepbind/db/db.tsv"
rootsqtl="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/sQTL_run/logit/JC"
rootsplicinginput='/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/input_splicing/logit/JC/'

#get all sQTL exons
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

#get all RBPs
RBPtable=read.table(RBPdb,sep="\t",header=T)
subRBPtable=subset(RBPtable,RBPtable[,"Species"]=="Homo sapiens")
subRBPtable=subset(subRBPtable,subRBPtable[,"Type"]=="RBP")

#change RBP name (the name of some RBPs in the DeepBind table is different with their name in the gene expression table. They are the same RBP but use different names)
subRBPtable=as.matrix(subRBPtable)
subRBPtable[which(subRBPtable[,"Protein"]=="BRUNOL4"),"Protein"]="CELF4"
subRBPtable[which(subRBPtable[,"Protein"]=="HuR"),"Protein"]="ELAVL1"
subRBPtable[which(subRBPtable[,"Protein"]=="STAR-PAP"),"Protein"]="TUT1"
subRBPtable[which(subRBPtable[,"Protein"]=="BRUNOL5"),"Protein"]="CELF5"
subRBPtable[which(subRBPtable[,"Protein"]=="BRUNOL6"),"Protein"]="CELF6"
subRBPtable[which(subRBPtable[,"Protein"]=="HNRPLL"),"Protein"]="HNRNPLL"

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
    if (SNPpos>=exonstart && SNPpos<=exonend){       #if the SNP is on exon body, the distance is 0
      dis=0
    }else{
    dis=min(abs(SNPpos-SS5start),abs(SNPpos-SS5end),abs(SNPpos-SS3start),abs(SNPpos-SS3end))  
    }
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
    if (SNPpos>=exonstart && SNPpos<=exonend){       #if the SNP is on exon body, the distance is 0
      dis=0
    }else{
    dis=min(abs(SNPpos-ass1_start),abs(SNPpos-ass1_end),abs(SNPpos-ass2_start),abs(SNPpos-ass2_end)) 
    }
  }
  return(dis)
}


for (splicetype in splicetypelist){
  for (type in typelist){
    rootoutput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.3_motif_analysis/individual_exon_binding_peaks_scan/DeepBind/deepbind_input/",splicetype,"/",type,sep="")
    command=paste("mkdir -p",rootoutput)
    system(command)
    
    #get all the sQTL exons
    sQTLexon=c()
    inputpath=paste(rootinput,splicetype,sep="/")
    setwd(inputpath)
    for (i in 1:length(brainregionlist)){
      temp=read.table(paste("all_info_exon_sQTL_GWAS_disease_in_",brainregionlist[i],"_logit_JC_",splicetype,"_",type,".txt",sep=""),sep="\t",header=T)
      sQTLexon=c(sQTLexon,as.character(temp[,"exon_full_ID"]))
    }
    sQTLexon=unique(sQTLexon)
    
    #exon information
    inputsplicing=paste(rootsplicinginput,splicetype,"/exon_info.fromGTF.",splicetype,".txt",sep="")
    exoninfo=read.table(inputsplicing,sep="\t",header=T)
    rownames(exoninfo)=exoninfo[,"ID"]
    
    #generate jobs for each sQTL exon - RBP - significant SNP trio
    joblist=matrix(,nrow=0,ncol=5)
    for (e in 1:length(sQTLexon)){
      shortID=paste("SE",strsplit(sQTLexon[e],split="\\|")[[1]][1],sep="_")
      
      #get all significant SNPs within 300bp
      for (b in 1:length(brainregionlist)){
        br=brainregionlist[b]
        sqtlpath=paste(rootsqtl,splicetype,br,paste("Glimmps_each_exon_cis_",br,sep=""),sep="/")
        setwd(sqtlpath)
        sqtlexon=read.table(paste(shortID,".asso",sep=""),sep="\t",header=T)
        
        #get the p value cutoff
        if (type=="pvalue"){
          cutoff=10^-5
        }
        if (type=="permutation"){
          setwd(paste(rootsqtl,splicetype,br,sep="/"))
          cutoff=as.numeric(as.matrix(read.table("permutation_p_value_cutoff.txt",sep="\t")))
        }
        
        #get significant SNPs
        sigsqtlexon=subset(sqtlexon,sqtlexon[,"pvals.lm"]<=cutoff)
        
        #get SNPs within 300bp
        distance=sapply(as.numeric(sigsqtlexon[,"Pos"]),checkdis,exoninfo=exoninfo,exonID=shortID,splicetype)
        sigsqtlexon=sigsqtlexon[which(distance<=300),]
        
        if (dim(sigsqtlexon)[1]>0){      #if there are significant SNPs within 300bp
          for (s in 1:dim(sigsqtlexon)[1]){         
            snpid=as.character(sigsqtlexon[s,"SNPID"])
            #for each RBP
            for (r in 1:dim(subRBPtable)[1]){
              rbpid=paste(subRBPtable[r,"ID"],subRBPtable[r,"Protein"],sep="_")
              
              #put the current sQTL exon - SNP - RBP trio into the matrix
              newline=c(shortID,snpid,rbpid,br,paste(shortID,snpid,rbpid,sep="_"))
              joblist=rbind(joblist,newline)
            }
          }
        }
      }
    }
    #remove duplications
    colnames(joblist)=c("Exon","SNP","RBP","Brain.region","label")
    rownames(joblist)=joblist[,"label"]
    temp=aggregate(joblist[,"Brain.region"], list(joblist[,"label"]), function(x) toString(paste0(unique(x))))   #merge brain region information for the same exon-SNP-RBP trio)
    uniquejoblist=joblist[temp[,1],]
    uniquejoblist[,"Brain.region"]=temp[,2]
    
    setwd(rootoutput)
    write.table(uniquejoblist,paste(splicetype,"_",type,"_uniquejoblist.txt",sep=""),sep="\t")
    
    print(paste(splicetype,type,dim(uniquejoblist)[1]))
  }
}



#output the exon ID conversion table (short ID to full ID)
for (splicetype in splicetypelist){
  for (type in typelist){
    rootoutput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.3_motif_analysis/individual_exon_binding_peaks_scan/DeepBind/deepbind_input/",splicetype,"/",type,sep="")
    command=paste("mkdir -p",rootoutput)
    system(command)
    
    #get all the sQTL exons
    sQTLexon=c()
    inputpath=paste(rootinput,splicetype,sep="/")
    setwd(inputpath)
    for (i in 1:length(brainregionlist)){
      temp=read.table(paste("all_info_exon_sQTL_GWAS_disease_in_",brainregionlist[i],"_logit_JC_",splicetype,"_",type,".txt",sep=""),sep="\t",header=T)
      sQTLexon=c(sQTLexon,as.character(temp[,"exon_full_ID"]))
    }
    sQTLexon=unique(sQTLexon)

    #get the short ID
    shortIDlist=rep(NA,length(sQTLexon))
    for (e in 1:length(sQTLexon)){
      shortIDlist[e]=paste("SE",strsplit(sQTLexon[e],split="\\|")[[1]][1],sep="_")
    }
    
    setwd(rootoutput)
    write.table(cbind(shortIDlist,sQTLexon),paste(splicetype,"_",type,"_short_long_exon_ID_conversion.txt",sep=""),sep="\t",row.names=F,col.names=F,quote=F)
  }
}


