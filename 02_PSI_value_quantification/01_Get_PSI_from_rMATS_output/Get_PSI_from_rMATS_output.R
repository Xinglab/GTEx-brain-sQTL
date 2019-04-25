splicetype="SE"      #type of alternative splicing, e.g., SE, A3SS, A5SS, MXE, IR
counttype="JC"       #JCEC (junction count + exon body count) or JC (junction count only)

##################
#Input parameters#
##################
inputpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/data/brain_rMATS_post"     

inputcount=paste(counttype,".raw.input.",splicetype,".txt",sep="")
inputexon=paste("fromGTF.",splicetype,".txt",sep="")
annotationpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/document/V7_annotation"
annotationname="GTEx_v7_Annotations_SubjectPhenotypesDS.txt"          #phenotype information
IDconversionpath=annotationpath
IDconversionname="gtex_v7_brain.csv"               
rootoutput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/1_get_matrix_from_rMATS/result/",counttype,sep="")

command=paste("mkdir -p ",rootoutput,sep="")    #create output folder if not exist
system(command)

b1path=inputpath
b1name="sample.txt"

######################
#read in count result#
######################
setwd(inputpath)
RCtable=read.table(inputcount,sep="\t",header=T)
samplesize=length(strsplit(as.character(RCtable[1,"IJC_SAMPLE_1"]),split=",")[[1]])

#get exon name as row name
setwd(inputpath)
exoninfo=read.table(inputexon,sep="\t",header=T)
nospace_exoninfo=apply(exoninfo,2,function(x)gsub('\\s+', '',x)) #remove all the white spaces in exoninfo
exonname=apply(nospace_exoninfo,1,paste,collapse="|")

#get SRR ID as column name
setwd(b1path)
raw_SRR=read.table(b1name,sep="\t")   
SRRlist=rep(NA,dim(raw_SRR)[1])
for(i in 1:dim(raw_SRR)[1]){
  SRRlist[i]=strsplit(strsplit(strsplit(as.character(raw_SRR[i,1]),split="/")[[1]][12],split="_")[[1]][1],split=".bam")[[1]][1]
}

#generate tables for output
oriICtable=matrix(NA,dim(RCtable)[1],samplesize)    #original inclusion count
rownames(oriICtable)=exonname
colnames(oriICtable)=SRRlist      #we use SRR ID because it is unique to each sample

oriSCtable=matrix(NA,dim(RCtable)[1],samplesize)      #original skipping count
rownames(oriSCtable)=exonname
colnames(oriSCtable)=SRRlist

oritotalRCtable=matrix(NA,dim(RCtable)[1],samplesize)     #original total count
rownames(oritotalRCtable)=exonname
colnames(oritotalRCtable)=SRRlist

PSItable=matrix(NA,dim(RCtable)[1],samplesize)
rownames(PSItable)=exonname
colnames(PSItable)=SRRlist

#################
#fill the tables#
#################
get_psi<-function(x){  
  #get count
  inclulength=as.numeric(as.character(x[,"IncFormLen"]))
  skilength=as.numeric(as.character(x[,"SkipFormLen"]))
  ratio=inclulength/skilength
  inclucount=as.numeric(strsplit(as.character(x[,"IJC_SAMPLE_1"]),split=",")[[1]])   #original inclusion count
  skicount=as.numeric(strsplit(as.character(x[,"SJC_SAMPLE_1"]),split=",")[[1]])     #original skipping count
  #calculate psi
  PSI=(inclucount/inclulength)/(inclucount/inclulength+skicount/skilength)
  return(list(inclucount,   #original inclusion count
              skicount,     #original skipping count
              inclucount+skicount,     #original total count
              PSI,           #inclusion level
              inclucount,         #normalized inclusion count
              skicount*ratio,     #normalized skipping count (we use skipping count*ratio instead of inclusion count/ratio)
              c(inclucount+skicount*ratio)))   #normalized total count
}

for (i in 1:dim(RCtable)[1]){
  temp=get_psi(RCtable[i,])
  oriICtable[i,]=temp[[1]]
  oriSCtable[i,]=temp[[2]]
  oritotalRCtable[i,]=temp[[3]]
  PSItable[i,]=temp[[4]]
}


###########
#Filtering#
###########
IC=oriICtable
totalRC=oritotalRCtable
PSI=PSItable

#we filter out exons with average PSI<5% or >95% and average total read count <10
avepsi=apply(PSI,1,mean,na.rm=T)
averc=apply(totalRC,1,mean,na.rm=T)    #we use original total read count

rccutoff=10
bool=intersect(intersect(which(avepsi>0.05),which(avepsi<0.95)),which(averc>=rccutoff))
PSI_filter=PSI[bool,]   

#since PCA cannot deal with missing data, we further remove exons with more than 5% missing value and then impute the rest
sumna=function(x){
  temp=sum(is.na(x))
  return(temp)
}
sumNA=apply(PSI_filter,1,sumna)
missing_cutoff=0.05
PSI_filter=PSI_filter[which((sumNA/dim(PSI_filter)[2])<missing_cutoff),]

#we also filter out exons with max(PSI)-min(PSI)<=5%
deltapsi=function(x){
  temp=max(x,na.rm=T)-min(x,na.rm=T)   #if we use this, we will have 12132 exons left
  return(temp)
}
delta_psi=apply(PSI_filter,1,deltapsi)
PSI_filter=PSI_filter[which(delta_psi>0.05),]

#output result
outputpath=rootoutput
setwd(outputpath)

IC_filter=IC[rownames(PSI_filter),]
totalRC_filter=totalRC[rownames(PSI_filter),]
write.table(PSI_filter,paste(splicetype,"_",counttype,"_PSI_filter.txt",sep=""),sep="\t")
write.table(IC_filter,paste(splicetype,"_",counttype,"_IC_filter.txt",sep=""),sep="\t")
write.table(totalRC_filter,paste(splicetype,"_",counttype,"_totalRC_filter.txt",sep=""),sep="\t")




