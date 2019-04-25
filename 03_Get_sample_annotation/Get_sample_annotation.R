outputfolder="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/2_get_phenotype_information/result"

#get SRR ID of all samples
setwd("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/1_get_matrix_from_rMATS/result/JC")
data=read.table("A3SS_JC_IC_filter.txt",sep="\t",header=T)
SRRlist=colnames(data)

#read in annotation information
#for age and gender
annotationpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/GTEx_V7_document/V7_annotation"
annotationname="GTEx_v7_Annotations_SubjectPhenotypesDS.txt"          #phenotype information
setwd(annotationpath)
annotation=read.table(annotationname,sep="\t",header=T)
rownames(annotation)=annotation[,"SUBJID"]

#for brain region
IDconversionpath=annotationpath
IDconversionname="gtex_v7_brain.csv"   
IDconversion=read.csv(IDconversionname,header=T)

#for batch effect
phenotypepath=annotationpath
phenotypename="GTEx_v7_Annotations_SampleAttributesDS.txt"
setwd(phenotypepath)
phenotype=read.table(phenotypename,sep="\t",header=T,quote="",fill=TRUE)
rownames(phenotype)=phenotype[,"SAMPID"]
  
############################
#get annotation information#
############################
samplesize=dim(data)[2]
gendertable=matrix(NA,1,samplesize)
rownames(gendertable)="gender"
colnames(gendertable)=SRRlist

agetable=matrix(NA,1,samplesize)
rownames(agetable)="age"
colnames(agetable)=SRRlist

#ID conversion
sampleID=rep(NA,length(SRRlist))
for (i in 1:length(sampleID)){
  sampleID[i]=as.character(IDconversion[which(IDconversion[,"Run_s"]==SRRlist[i]),"submitted_subject_id_s"])[1]   #there are duplications in the annotation file
}


for (j in 1:length(sampleID)){
  #gender information
  gendertable[1,j]=annotation[sampleID[j],"SEX"]    #we use sampleID here because there is no SRRID information in the annotation file. But this is ok because SRRID and sampleID are in the same order
  #age information
  agetable[1,j]=as.character(annotation[sampleID[j],"AGE"])
}

setwd(outputfolder)
sampleID_SRRID=cbind(SRRlist,sampleID)
rownames(sampleID_SRRID)=SRRlist
colnames(sampleID_SRRID)=c("SRRID","sampleID")

write.table(sampleID_SRRID,"sampleID_SRRID_brain.txt",sep="\t")
write.table(gendertable,"gendertable_brain.txt",sep="\t",quote=F)
write.table(agetable,"agetable_brain.txt",sep="\t",quote=F)

##############################
#get brain region information#
##############################
brain_region=as.matrix(agetable)
rownames(brain_region)="brain_region"

for (i in 1:dim(brain_region)[2]){
  SRRID=colnames(brain_region)[i]
  brain_region[1,i]=as.character(IDconversion[which(IDconversion[,"Run_s"]==SRRID),"body_site_s"])[1]   #there are duplications in the annotation file
}
write.table(brain_region,"brain_region_table_brain.txt",sep="\t")

##############################
#get batch effect information#
##############################
brain_sample_ID=rep(NA,length(SRRlist))
for (i in 1:length(SRRlist)){
  brain_sample_ID[i]=as.character(IDconversion[which(IDconversion[,"Run_s"]==SRRlist[i]),"Sample_Name_s"])[1]   #there are duplications in the annotation file
}

#get batch effect
#SMCENTER: Code for BSS collection site
#SMNABTCH: Nucleic Acid Isolation Batch ID, batch when DNA/RNA was isolated and extracted from a sample
#SMNABTCHT: Type of nucleic acid isolation batch
#SMNABTCHD: Date of nucleic acid isolation batch
#SMGEBTCH: Genotype or Expression Batch ID
#SMGEBTCHD: Date of genotype or expression batch
#SMGEBTCHT: Type of genotype or expression batch
batch_name=c("SMCENTER","SMNABTCH","SMNABTCHT","SMNABTCHD","SMGEBTCH","SMGEBTCHD","SMGEBTCHT")
batch=phenotype[brain_sample_ID,batch_name]
batch=t(batch)
colnames(batch)=SRRlist

write.table(batch,"batch_effect_table_brain.txt",sep="\t")

