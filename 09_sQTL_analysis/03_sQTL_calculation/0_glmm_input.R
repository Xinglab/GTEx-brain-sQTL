#The purpose of this code is to generate sQTL input for all types of alternative splicing#
splicetype="SE"      #type of alternative splicing, e.g., SE, A3SS, A5SS, MXE, IR
counttype="JC"       #JCEC (junction count + exon body count) or JC (junction count only)
PSItype="logit"      #logit (logit transformed PSI value) or original (original PSI from 0 to 1)

###################
#read in PSI table#
###################
if (PSItype=="original"){   
  inputpath="/input/path/for/original/PSI/value"
  inputname=paste(splicetype,"_",counttype,"_PSI_filter.txt",sep="")
}
if (PSItype=="logit"){     #logit PSI with correction
  inputpath="/input/path/for/logit/transformed/PSI/value"
  inputname=paste(splicetype,"_",counttype,"_logit_corrected_PSI.txt",sep="")
}
exoninfopath="/path/to/rMATS/post/step/output"
exoninfoname=paste("fromGTF.",splicetype,".txt",sep="")
inputbrainpath="./03_Get_sample_annotation/example_output"
inputbrainname="brain_region_table_brain.txt"


setwd(inputpath)
PSI=read.table(inputname,sep="\t",header=T)
PSIdenominator=matrix(1,dim(PSI)[1],dim(PSI)[2])        #the all 1 matrix used as total count
rownames(PSIdenominator)=rownames(PSI)
colnames(PSIdenominator)=colnames(PSI)
PSIdenominator[is.na(PSI)]=0      #change missing value to 0 to mimic the behavior of count data
PSI[is.na(PSI)]=0           #change missing value to 0 to mimic the behavior of count data
norIC=PSI
nortotalRC=PSIdenominator

setwd(exoninfopath)
exoninfo=read.table(exoninfoname,sep="\t",header=T)
nospace_exoninfo=apply(exoninfo,2,function(x)gsub('\\s+', '',x)) #remove all the white spaces in exoninfo
exonname=apply(nospace_exoninfo,1,paste,collapse="|")
exoninfo=exoninfo[which(exonname %in% rownames(PSI)),]


####################
#read in annotation#
####################
#read in brain region information
setwd(inputbrainpath)
oribrainregion=read.table(inputbrainname,sep="\t",header=T)
BRlist=unique(as.character(as.matrix(oribrainregion)))


#################
#generate output#
#################
outputpath="/output/path"
command=paste("mkdir -p ",outputpath,sep="")
system(command)
setwd(outputpath)

#output all SRR IDs into one file
write.table(colnames(PSI),"SRR_ID.txt",sep="\t",row.names=F,col.names=F,quote=F)

#modify exon ID in exoninfo
for (i in 1:dim(PSI)[1]){
  temp=strsplit(rownames(PSI)[i],split="\\|")[[1]]
  if (exoninfo[i,"ID"]==temp[1]){   #make sure that the order of exoninfo and order of PSI is the same
    temp[1]=paste("SE_",temp[1],sep="")
    exoninfo[i,"ID"]=temp[1]
  }else{
    print("order doesn't match")
  }
}
write.table(exoninfo,paste("exon_info.fromGTF.",splicetype,".txt",sep=""),row.names=F,sep="\t",quote=F)

#generate the count table for each brain region
for (BR in BRlist){
  tempBR=gsub("(","",BR,fixed=TRUE)
  tempBR=gsub(")","",tempBR,fixed=TRUE)
  subfolder=gsub(" ", "", tempBR, fixed = TRUE)
  command=paste("mkdir",subfolder)
  system(command)
  
  setwd(paste(outputpath,"/",subfolder,sep=""))
  subnorIC=norIC[,which(as.matrix(oribrainregion) %in% BR)]
  subnortotalRC=nortotalRC[,which(as.matrix(oribrainregion) %in% BR)]
  
  #output SRR ID
  write.table(colnames(subnorIC),"SRR_ID.txt",sep="\t",row.names=F,col.names=F,quote=F)
  
  #1. change the exon name
  rownames(subnorIC)=exoninfo[,"ID"]
  rownames(subnortotalRC)=exoninfo[,"ID"]
  #2. transpose
  subnorIC=t(subnorIC)
  subnortotalRC=t(subnortotalRC)
  #3. add one column in the beginning
  IID=rep(1,dim(subnorIC)[1])
  subnorIC=cbind(IID,subnorIC)
  subnortotalRC=cbind(IID,subnortotalRC)
  #3. remove row name (SRR ID) and put that as a column
  FID=rownames(subnorIC)
  subnorIC=cbind(FID,subnorIC)
  subnortotalRC=cbind(FID,subnortotalRC)
  
  #output
  write.table(subnorIC,"GTEx_brain_IC.txt",sep="\t",quote=F,row.names=F)
  write.table(subnortotalRC,"GTEx_brain_totalRC.txt",sep="\t",quote=F,row.names=F)
  
  setwd(outputpath)  #go back to the root folder
}


