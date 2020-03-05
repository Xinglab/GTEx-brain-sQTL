#Purpose of this code:
#For all the sQTLs we identified, we focus on the sQTLs with GWAS association
#for the associated GWAS traits, we download the summary statistics of those GWAS studies

######################################################Search and collect GWAS summary statistics##################################################
outputpath="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/0_search_and_download_summary_statistics"
#########################
#read in the sQTL result#
#########################
splicetype="SE"
PSItype="logit"
counttype="JC"
type="pvalue"
rootinput="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/summary"

setwd(paste(rootinput,PSItype,counttype,splicetype,sep="/"))
sqtl=read.table(paste("sQTL_summary_all_brainregion_",PSItype,"_",counttype,"_",splicetype,"_",type,"_sorted.txt",sep=""),sep="\t",header=T,quote="")

################################
#read in the GWAS catalog table#
################################
GWASdb="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project/analysis/1_GLMM/sQTL/GWAS_databases/gwas_catalog_v1.0.1-associations_e89_r2017-07-31.tsv"
GWAS=read.table(GWASdb,sep="\t",header=T,fill=TRUE,check.names=F,quote="",comment.char="")

disease_key_word=toupper(c("Alzheimer","Amyotrophic lateral sclerosis","Parkinson","frontotemporal dementia","Huntington",
                           "epilepsy","autism","schizophrenia","bipolar","depression",
                           "attention deficit hyperactivity disorder","glio","multiple sclerosis","narcolepsy","stroke"))       #change all the names to upper case for comparison

###################################################
#read in the GWAS catalog summary statistics table#
###################################################
#this table contains all the GWAS studies with summary statistics on GWAS catalog
setwd("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project/analysis/1_GLMM/sQTL/GWAS_databases/Summary_statistics")
gwasstudy=read.csv("GWAS_catalog_summary_statistics_study_table.csv",header=T)
gwaslink=read.csv("GWAS_catalog_summary_statistics_FTP_link.csv",header=T)

rownames(gwasstudy)=paste(gsub(" ", "", gwasstudy[,"First.Author"], fixed = TRUE),gwasstudy[,"PubMed.ID"],gwasstudy[,"Study.accession"],sep="_")    #generate unique identifier of each study

####################################################
#1. get sQTLs with GWAS association
sqtl_select=subset(sqtl,!is.na(sqtl[,"disease.ontology"]))

#2. get their GWAS summary statistics information for each sQTL 
GWAS_SS=matrix(NA,dim(sqtl_select)[1],6)
colnames(GWAS_SS)=c("Neurological.asso",         #1: related to neurological disorders, 0: not related
                    "First.Author","PubMed.ID","Study.accession","FTP.link",
                    "downloadable")              #1: with at least one summary statistics to download, 0: no downloadable summary statistics

for (i in 1:dim(sqtl_select)[1]){        #for each sQTL
  #check if the event is related to neurological disorders
  DO=toupper(sqtl_select[i,"disease.ontology"])
  if (grepl(paste(disease_key_word, collapse="|"), DO)){     #if the disease ontology contains any of the 10 disease names
    GWAS_SS[i,"Neurological.asso"]=1
  }else{
    GWAS_SS[i,"Neurological.asso"]=0
  } 
  
  #for each sQTL event, get the GWAS SNP rsID
  rsID=as.character(sqtl_select[i,"GWAS.rsID"])
  
  #get the GWAS study pubmed ID(s) and study accession number(s) of the GWAS snp
  info=GWAS[which(GWAS[,"SNPS"] %in% rsID),]
  pubmedID=""
  studyaccess=""
  firstauthor=""
  ftplink=""
  
  downloadable=0       #we assume the event has no GWAS study with summary statistics
  
  for (j in 1:dim(info)[1]){       #there may be multiple GWAS studies related to the GWAS SNP
    temppubmedID=info[j,"PUBMEDID"]
    tempstudyaccess=as.character(info[j,"STUDY ACCESSION"])
    
    #for some GWAS entries, the first author information in GWAS catalog table is different from the summary statistics table
    if (temppubmedID %in% gwasstudy[,"PubMed.ID"] && tempstudyaccess %in% gwasstudy[,"Study.accession"]){
      #we use the first author informaiton in the summary statistics table (if the study has summary statistics) since it is directly used in the FTP link
      a=which(gwasstudy[,"PubMed.ID"]==temppubmedID)
      b=which(gwasstudy[,"Study.accession"]==tempstudyaccess)
      tempfirstauthor=gsub(" ", "", as.character(gwasstudy[intersect(a,b),"First.Author"]), fixed = TRUE) 
      #tempfirstauthor=gsub(" ", "", as.character(gwasstudy[b,"First.Author"]), fixed = TRUE) 
    }else{
      #if the study has no summary statistics, we have to use the first author information in the GWAS catalog table, which will not affect downstream analysis anyway
      tempfirstauthor=gsub(" ", "", as.character(info[j,"FIRST AUTHOR"]), fixed = TRUE) 
    }
    
    uniqueID=paste(tempfirstauthor,temppubmedID,tempstudyaccess,sep="_")
    
    #for each study, get the GWAS summary statistics information related to this GWAS study
    if (uniqueID %in% rownames(gwasstudy)){       #if the GWAS study has GWAS summary statistics
      #get the FTP link
      search=sapply(as.character(gwaslink[,1]),grepl,pattern=uniqueID)
      tempftplink=as.character(gwaslink[which(search),1])
      
      downloadable=1           #if at least one of the GWAS studies has summary statistics, this event is considered downloadable
    }else{
      tempftplink="no_link"
    }
    
    pubmedID=paste(pubmedID,temppubmedID,sep=",")
    studyaccess=paste(studyaccess,tempstudyaccess,sep=",")
    firstauthor=paste(firstauthor,tempfirstauthor,sep=",")
    ftplink=paste(ftplink,tempftplink,sep=",")
  }

  GWAS_SS[i,"PubMed.ID"]=pubmedID
  GWAS_SS[i,"Study.accession"]=studyaccess
  GWAS_SS[i,"First.Author"]=firstauthor
  GWAS_SS[i,"FTP.link"]=ftplink
  GWAS_SS[i,"downloadable"]=downloadable
}


setwd(outputpath)
result=cbind(sqtl_select,GWAS_SS)
write.table(result,paste("sQTL_GWAS_summary_statistics_",splicetype,"_",PSItype,"_",counttype,"_",type,".txt",sep=""),sep="\t",row.names=F)

#get the subset with downloadable summary statistics
gwasftp=subset(result,result[,"downloadable"]==1)

#get the subset with downloadable summary statistics && neurological disorder related 
gwasftpnd=subset(gwasftp,gwasftp[,"Neurological.asso"]==1)
setwd(outputpath)
write.table(gwasftpnd,paste("sQTL_GWAS_summary_statistics_downloadable_disease_related_",splicetype,"_",PSItype,"_",counttype,"_",type,".txt",sep=""),sep="\t",row.names=F)







