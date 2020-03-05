#Purpose of this code:
#1. for all the zipped files, we decompress them

inputpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/0_search_and_download_summary_statistics"
rootoutputpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/GWAS_summary_statistics"
setwd(inputpath)
GWAS_SS=read.table("sQTL_GWAS_summary_statistics_SE_logit_JC_pvalue.txt",sep="\t",header=T)

#get the subset with downloadable summary statistics
gwasftp=subset(GWAS_SS,GWAS_SS[,"downloadable"]==1)

alluniqueID=c()    #all the unique study ID
for (i in 1:dim(gwasftp)[1]){
  authorlist=strsplit(as.character(gwasftp[i,"First.Author"]),split=",")[[1]][-1]
  pubmedidlist=strsplit(as.character(gwasftp[i,"PubMed.ID"]),split=",")[[1]][-1]
  studyaccesslist=strsplit(as.character(gwasftp[i,"Study.accession"]),split=",")[[1]][-1]
  ftplinklist=strsplit(as.character(gwasftp[i,"FTP.link"]),split=",")[[1]][-1]
  
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


foldernamelist=unique(alluniqueID)
for (i in 1:length(foldernamelist)){
  path=paste(rootoutputpath,foldernamelist[i],sep="/")
  setwd(path)
  #check if there is any zip file
  ziplist=system("ls *.zip",intern=T)
  if (length(ziplist)==0){        #if there is no .zip file
    donothing="donothing"
  }else{
    for (j in 1:length(ziplist)){
      command=paste("unzip",ziplist[j])
      system(command)
    }
  }
  
  #check if there is any gz file
  gzlist=system("ls *.gz",intern=T)
  if (length(gzlist)==0){        #if there is no .gz file
    donothing="donothing"
  }else{
    for (k in 1:length(gzlist)){
      command=paste("gunzip",gzlist[k])
      system(command)
    }
  }
  
  
  
  #check if there is the harmonised folder
  if (dir.exists("harmonised")){
    subpath=paste(path,"harmonised",sep="/")
    setwd(subpath)
    
    #check if there is any zip file
    ziplist=system("ls *.zip",intern=T)
    if (length(ziplist)==0){        #if there is no .zip file
      donothing="donothing"
    }else{
      for (j in 1:length(ziplist)){
        command=paste("unzip",ziplist[j])
        system(command)
      }
    }
    
    #check if there is any gz file
    gzlist=system("ls *.gz",intern=T)
    if (length(gzlist)==0){        #if there is no .gz file
      donothing="donothing"
    }else{
      for (k in 1:length(gzlist)){
        command=paste("gunzip",gzlist[k])
        system(command)
      }
    }
  }
}




