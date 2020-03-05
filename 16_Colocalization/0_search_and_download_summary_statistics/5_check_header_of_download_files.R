#Purpose of this code:
#for all harmonized summary statistics, we checked the header of them to see:
#(1) are they consistent with each other
#(2) do they have sample size information

inputpath="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/0_search_and_download_summary_statistics"
rootoutputpath="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/GWAS_summary_statistics"
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

  #check if there is the harmonised folder
  if (dir.exists("harmonised")){
    subpath=paste(path,"harmonised",sep="/")
    setwd(subpath)
    
    #get the header information
    filelist=system("ls",intern=T)
    filename=filelist[grep("h.tsv",filelist)]

    command=paste("head -n 1",filename)
    header=system(command,intern=T)
  
    print(foldernamelist[i])
    print(gsub("\t","   ",header))
    print("     ")
    print("     ")
  }
}




