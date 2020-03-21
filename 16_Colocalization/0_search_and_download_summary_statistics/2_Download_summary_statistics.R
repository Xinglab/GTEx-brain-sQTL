#Purpose of this code:
#For all the sQTLs we identified, we focus on the sQTLs with GWAS association
#for the associated GWAS traits, we download the summary statistics of those GWAS studies

######################################################Download GWAS summary statistics##################################################
#############################
#Download summary statistics#
#############################
inputpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/0_search_and_download_summary_statistics"
rootoutputpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/GWAS_summary_statistics"
setwd(inputpath)
GWAS_SS=read.table("sQTL_GWAS_summary_statistics_SE_logit_JC_pvalue.txt",sep="\t",header=T)

#get the subset with downloadable summary statistics
gwasftp=subset(GWAS_SS,GWAS_SS[,"downloadable"]==1)

#download the summary statistics
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
      outputpath=paste(rootoutputpath,label,sep="/")
      #check if the folder is already there
      if (dir.exists(outputpath)){             #if the folder is already there, it means we have already downloaded the summary statistics of this study
        donothing="donothing"
      }else{
        command=paste("mkdir -p",outputpath)
        system(command)
        setwd(outputpath)
        write.table(ftplink,paste(label,"ftplink.txt",sep="_"),sep="",row.names=F,col.names=F,quote=F)
        
        #download the data   
        command=paste("wget -r",ftplink)
        system(command)
        
        #move the data to the output folder
        command=paste("mv",paste(outputpath,strsplit(ftplink,split="ftp://")[[1]][2],"*",sep="/"),outputpath)
        system(command)
        #the reason that we need to do this is because in a wget command, the link will also be interpreted as a path
        #for example, the downloaded file will be stored in outputpath/ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/LockeAE_25673413_GCST002783/
        #rather than just outputpath
        
        #remove the empty folder
        system("rm -r ftp.ebi.ac.uk")
      }
    }
  }
}









