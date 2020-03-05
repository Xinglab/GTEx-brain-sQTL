#Purpose of this code:
#for the sQTL exon we want to do colocalization test on, we want to calculate the MAF of all the SNPs tested in the sQTL calculation of this exon
#So this code will submit each sQTL exon-brain region combination as a job, for each job, we calculate the MAF and other information

###############################################################get the exon-brain region combination####################################################################################

job <- as.numeric(Sys.getenv("SGE_TASK_ID"))
print(paste("job ID:",job))

#flashscratch version (when I run the analysis, the nobackup/yxing folder is under migration so this part of the analysis was run under the scratch folder)
inputpath="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/0_search_and_download_summary_statistics"
rootoutputpath="/u/flashscratch/y/yidazhan/GTEx_V7_analysis/1_MAF_calculation/result"
exoninfopath="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/input_splicing/logit/JC/SE"
codefolder="/u/flashscratch/y/yidazhan/GTEx_V7_analysis/1_MAF_calculation"

setwd(inputpath)
GWAS_SS=read.table("sQTL_GWAS_summary_statistics_SE_logit_JC_pvalue.txt",sep="\t",header=T)

#exon full ID - short ID conversion
setwd(exoninfopath)
exonfullinfo=read.table("exon_info.fromGTF.SE.txt",sep="\t",header=T)
rownames(exonfullinfo)=exonfullinfo[,"ID"]

#brain region name conversion
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
formalbrainregionlist=c("Amygdala",
                        "Anterior cingulate cortex BA24",
                        "Caudate basal ganglia",
                        "Cerebellar Hemisphere",
                        "Cerebellum",
                        "Cortex",
                        "Frontal Cortex BA9",
                        "Hippocampus",
                        "Hypothalamus",
                        "Nucleus accumbens basal ganglia",
                        "Putamen basal ganglia",
                        "Spinal cord cervical c-1",
                        "Substantia nigra")


#get the subset with downloadable summary statistics
gwasftp=subset(GWAS_SS,GWAS_SS[,"downloadable"]==1)

#get the list of unique exons
exonlist=as.character(unique(gwasftp[,"exon_short_ID"]))

#get the significant brain regions of each exon (for the same exon, the SNP could be different so the significant brain regions could also be different, here we take the union for each exon)
uniqueallsigregion=rep(NA,length(exonlist))
for (i in 1:length(exonlist)){
  exonshortID=exonlist[i]
  allsigregion=as.character(gwasftp[which(gwasftp[,"exon_short_ID"] %in% exonshortID),"significant.region"])
  if (length(allsigregion)==1){           #if the exon only appears once
    uniqueallsigregion[i]=allsigregion
  }else{                                  #if the exon appears multiple times (either because multiple GWAS studies or different sQTL SNPs)
    uniqueallsigregion[i]=paste(unique(strsplit(paste(allsigregion,collapse=", "),split=", ")[[1]]),collapse=", ")
  }
}


#run the calculation for the current exon-brain region combination
if ((job %% length(brainregionlist))==0){
  i=(job %/% length(brainregionlist))
  br=length(brainregionlist)
}else{
  i=(job %/% length(brainregionlist))+1
  br=job %% length(brainregionlist)
}


exonshortID=exonlist[i]

#get the full exon ID
exonfullID=paste(c(strsplit(exonshortID,split="_")[[1]][2],as.matrix(exonfullinfo[exonshortID,])[-1]),collapse="|")

#generate input folder for that exon
exonpath=paste(rootoutputpath,gsub("\\|",",",exonfullID),sep="/")
command=paste("mkdir -p",exonpath)
system(command)

sig.region=formalbrainregionlist         #we carry out the calculation on all brain regions, not just the significant brain regions
currentBR=brainregionlist[which(formalbrainregionlist %in% sig.region[br])]

exonbrpath=paste(exonpath,currentBR,sep="/")
command=paste("mkdir -p",exonbrpath)
system(command)

###############################################################start the calculation for each combination####################################################################################

brainregion=currentBR
testedexon=exonfullID
#brainregion="Brain-Cerebellum"
#testedexon="189436|ENSG00000186868.15_3|MAPT|chr17|+|44051750|44051837|44049224|44049311|44055740|44055806"

#sampleinputroot="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/input_splicing/logit/JC/SE"
#sQTLruninputroot="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/sQTL_run/logit/JC/SE"
#rootsqtl="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/sQTL_run"
#genoplinkprefix="Genotype_V7_plink_maf0.05"
#rootoutput="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/1_MAF_calculation"
#genotype="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project/data/raw_data/files/genotype/V7_whole_exon_sequencing/GTEx_Analysis_2016-01-15_v7_WGS_652ind_VarID_Lookup_Table.txt"

#flashscratch version
sampleinputroot="/u/project/yxing-read-only/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/input_splicing/logit/JC/SE"
sQTLruninputroot="/u/project/yxing-read-only/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/sQTL_run/logit/JC/SE"
rootsqtl="/u/project/yxing-read-only/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/sQTL_run"
genoplinkprefix="Genotype_V7_plink_maf0.05"
rootoutput="/u/flashscratch/y/yidazhan/GTEx_V7_analysis/1_MAF_calculation"
genotype="/u/project/yxing-read-only/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project/data/raw_data/files/genotype/V7_whole_exon_sequencing/GTEx_Analysis_2016-01-15_v7_WGS_652ind_VarID_Lookup_Table.txt"


PSItype="logit"
counttype="JC"
splicetype="SE"

library("data.table")

##############################################################################
#step 1: get all the SNPs tested for the given exon in the given brain region#
##############################################################################
shortexonID=paste("SE",strsplit(testedexon,split="\\|")[[1]][1],sep="_")
sQTLruninput=paste(sQTLruninputroot,brainregion,paste("Glimmps_each_exon_cis_",brainregion,sep=""),sep="/")
setwd(sQTLruninput)
sQTLresult=read.table(paste(shortexonID,".asso",sep=""),sep="\t",header=T)


################################################################################
#Step 2: calculate MAF and sample size (with genotype information) for each SNP#
################################################################################
#get the sample ID of samples in the current brain region
setwd(paste(strsplit(rootsqtl,split="sQTL_run")[[1]],"/input_splicing/",PSItype,"/",counttype,"/",splicetype,"/",brainregion,sep=""))
IDs.pheno=as.character(as.matrix(read.table("SRR_ID.txt",sep="\t")))

#get the genotype information of all the SNPs in the current brain region
temp=strsplit(testedexon,split="\\|")[[1]]
chrom=temp[4]      #the genotype information is by chromosome
genotypename = paste(rootsqtl,"/",PSItype,"/",counttype,"/",splicetype,"/",brainregion,sep="")
genoplinkfile= genoplinkprefix
genofile = genoplinkfile
coordinate=paste(temp[4],":",temp[6],"-",temp[7],sep="")
genesymbol=temp[3]
#.map file for the current chromosome
c <-fread(paste(genotypename,"/",chrom,".",genofile, ".map", sep=""), header=F, sep="\t")    # Note: check sep (My code)
map <- data.frame(c)
names(map) <- c("Chr","SNPID","cM","Pos")
#.raw file for the current chromosome
d <- fread(paste(genotypename,"/",chrom, ".", genofile, ".map_tpose.raw", sep=""), header=T,sep="\t",  na.strings=c("NA"))     # (My code)
genoplink <- data.frame(d)
genoplink=genoplink[,-1] 
IDs.geno <- names(genoplink)           #sample ID of samples with genotype information in the current brain region
#IDs.common <- IDs.geno[is.element(IDs.geno, IDs.pheno)]     #get the samples with genotype information
IDs.common <- intersect(IDs.geno,IDs.pheno)
nsnps <- dim(genoplink)[1] -5                              #remove the first a few lines
#sub.geno <- match(IDs.common , IDs.geno)          #the position of samples in IDs.geno that are also in IDs.common
#geno <- as.matrix(genoplink[seq(1, nsnps)+5, sub.geno])    #  I changed here (we basically removed the header and first column from d)
geno <- as.matrix(genoplink[seq(1, nsnps)+5, IDs.common])   #this is the genotype table for all SNPs

#calculate MAF and sample size for each SNP
MAF.SS=matrix(NA,dim(sQTLresult)[1],2)
colnames(MAF.SS)=c("MAF","Sample.size")
rowname=rep(NA,dim(MAF.SS)[1])
rsID=rep(NA,dim(MAF.SS)[1])

for (i in 1:dim(sQTLresult)[1]){     #for each SNP
  print(i)
  exongenoname=sQTLresult[i,"SNPID"]
  #get the rsID
  command=paste("grep ","\"",exongenoname,"\""," ",genotype,sep="")
  IDconversion=system(command,intern=T)
  if (strsplit(IDconversion,split="\t")[[1]][6]!="."){       #if we have the rsID based on RS_ID_dbSNP142_GRCh37p13
    rsID[i]=strsplit(IDconversion,split="\t")[[1]][6]
  }else{                                                     #if we don't have, we try RS_ID_dbSNP147_GRCh37p13
    rsID[i]=strsplit(IDconversion,split="\t")[[1]][7]
  }
  
  
  #calculate MAF and sample size
  exongeno=geno[which(map[,"SNPID"]==as.character(exongenoname)),]
  rowname[i]=exongenoname
  ratio=table(exongeno)/sum(table(exongeno))
  
  if ("0" %in% names(ratio)){           #if we have samples with genotype as 0
    ratio0=ratio[which(names(ratio)=="0")]
  }else{
    ratio0=0                            #if we don't have samples with genotype as 0
  }
  
  if ("1" %in% names(ratio)){           #if we have samples with genotype as 1
    ratio1=ratio[which(names(ratio)=="1")]
  }else{
    ratio1=0                            #if we don't have samples with genotype as 1
  }
  
  if ("2" %in% names(ratio)){           #if we have samples with genotype as 2
    ratio2=ratio[which(names(ratio)=="2")]
  }else{
    ratio2=0                            #if we don't have samples with genotype as 2
  }
  
  ratioA=ratio0+ratio1/2
  ratioa=ratio2+ratio1/2
  MAF.SS[i,"MAF"]=min(ratioA,ratioa)
  
  MAF.SS[i,"Sample.size"]=sum(!is.na(exongeno))
}

rownames(MAF.SS)=rowname

output=cbind(MAF.SS,sQTLresult,rsID)
outputpath=paste(rootoutput,"/result/",gsub("\\|",",",testedexon),"/",brainregion,sep="")
#command=paste("mkdir -p",outputpath)
#system(command)

setwd(outputpath)
write.table(output,"Colocalization_input.txt",sep="\t")


print("finished successfully!")











