#Purpose of this code:
#1. we check if we have downloaded all the downloadable GWAS summary statistics or not

inputpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/0_search_and_download_summary_statistics"
rootoutputpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/GWAS_summary_statistics"
setwd(inputpath)
GWAS_SS=read.table("sQTL_GWAS_summary_statistics_SE_logit_JC_pvalue.txt",sep="\t",header=T)

#get the subset with downloadable summary statistics
gwasftp=subset(GWAS_SS,GWAS_SS[,"downloadable"]==1)

######################################################check if we have downloaded all the summary statistics available######################################################
allftplink=c()     #all the unique FTP links
alluniqueID=c()    #all the unique study ID
for (i in 1:dim(gwasftp)[1]){
  authorlist=strsplit(as.character(gwasftp[i,"First.Author"]),split=",")[[1]][-1]
  pubmedidlist=strsplit(as.character(gwasftp[i,"PubMed.ID"]),split=",")[[1]][-1]
  studyaccesslist=strsplit(as.character(gwasftp[i,"Study.accession"]),split=",")[[1]][-1]
  ftplinklist=strsplit(as.character(gwasftp[i,"FTP.link"]),split=",")[[1]][-1]
  
  allftplink=c(allftplink,ftplinklist)
  
  for (j in 1:length(ftplinklist)){
    if (ftplinklist[j]=="no_link"){      #for sQTLs with multiple GWAS studies, not all of them have summary statistics to download
      donothing="donothing"
    }else{
      author=authorlist[j]
      pubmedid=pubmedidlist[j]
      studyaccess=studyaccesslist[j]
      ftplink=ftplinklist[j]
      
      label=paste(author,pubmedid,studyaccess,sep="_")
      alluniqueID=c(alluniqueID,label)
    }
  }
}

#total number of GWAS summary statistics to download:
length(unique(alluniqueID))
#36
length(unique(allftplink))
#37      #this one also includes "no_link" as a FTP link so it is 37, not 36







