job <- as.numeric(Sys.getenv("SGE_TASK_ID"))
print(paste("job ID:",job))


gwascatalogfile="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/0_search_and_download_summary_statistics/sQTL_GWAS_summary_statistics_SE_logit_JC_pvalue.txt"
colocresultinputpath="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/3_Colocalization_analysis"     #colocalization analysis result
rootcolocinputpath="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/2_Colocalization_input/result"   #input of colocalization analysis
rootoutputpath="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/4_Colocalization_plot/result"
locuszoom="/u/project/yxing/PROJECT/yidazhan/research/software/locuszoom/bin/locuszoom"

#read in the colocalization result
setwd(colocresultinputpath)
colocresult=read.table("colocalization_result_summary.txt",sep="\t",header=T)

#read in the GWAS catalog information
gwascatalog=read.table(gwascatalogfile,sep="\t",header=T)

getposfromSNPID=function(x){           #return the genomic coordinate from SNP ID
  pos=strsplit(x,split="_")[[1]][2]
  return(pos)
}

extend=100     #extend the plot for 100bp on each side
cutoff=0.75                                                     #PP4 cutoff for significant result


testlabel=colnames(colocresult)[job]                              #there are different ways to do the colocalization analysis, we need to choose the result of one of them
siglist=subset(colocresult,colocresult[,testlabel]>=cutoff)

print(paste("number of plot",dim(siglist)[1],sep=" : "))

if (dim(siglist)[1]>0){
  #for the significant result, collect the input information for plot and make the plot
  for (i in 1:dim(siglist)[1]){
    
    print(paste(paste(testlabel,i,sep=" : ")," out of ",dim(siglist)[1],sep=""))
    
    temp=strsplit(rownames(siglist)[i],split="~")[[1]]
    exonfullID=temp[1]
    gwasstudy=temp[2]
    currentBR=temp[3]
    
    #get the input file for colocalization analysis
    inputpath=paste(rootcolocinputpath,exonfullID,gwasstudy,currentBR,sep="/")
    setwd(inputpath)
    colocinput=read.table("colocalization_input.txt",sep="\t",header=T)
    
    #generate output folder and input file for locuszoom plot
    outputpath=paste(rootoutputpath,testlabel,rownames(siglist)[i],sep="/")
    command=paste("mkdir -p",outputpath)
    system(command)
    
    setwd(outputpath)
    input4plot.gwas=colocinput[,c("sQTL.rsID","GWAS.p_value")]
    input4plot.sqtl=colocinput[,c("sQTL.rsID","sQTL.pvals.lm")]
    colnames(input4plot.gwas)=colnames(input4plot.sqtl)=c("MarkerName","P-value")
    write.table(input4plot.gwas,"input4plot_gwas.txt",sep="\t",row.names=F,quote=F)
    write.table(input4plot.sqtl,"input4plot_sqtl.txt",sep="\t",row.names=F,quote=F)
    
    #collect other information for plot
    #1. range of the plot
    poslist=as.numeric(sapply(as.matrix(colocinput[rownames(input4plot.gwas),"sQTL.SNPID"]),getposfromSNPID))
    #2. chromosome information
    chr=strsplit(exonfullID,split=",")[[1]][4]
    #3. other information (disease relationship, GWAS terms, etc)
    genesymbol=strsplit(exonfullID,split=",")[[1]][3]
    exonID=strsplit(exonfullID,split=",")[[1]][2]
    firstauthor=strsplit(gwasstudy,split="_")[[1]][1]
    pubmedid=strsplit(gwasstudy,split="_")[[1]][2]
    studyaccess=strsplit(gwasstudy,split="_")[[1]][3]
    rownum=intersect(intersect(intersect(intersect(which(gwascatalog[,"gene.symbol"] %in% genesymbol),which(gwascatalog[,"exon.ID"] %in% exonID)),
                                         which(grepl(firstauthor, gwascatalog[,"First.Author"]))),
                               which(grepl(pubmedid, gwascatalog[,"PubMed.ID"]))),
                     which(grepl(studyaccess,gwascatalog[,"Study.accession"])))
    gwasinfo=gwascatalog[rownum,]
    write.table(gwasinfo,"GWAS_information_of_the_current_event.txt",sep="\t")
    
    write.table(siglist[i,],"Colocalization_information_of_the_current_event.txt",sep="\t")
    
    
    #make the plot for GWAS result
    command=paste(locuszoom,"--metal input4plot_gwas.txt",
                  "--pop EUR --build hg19 --source 1000G_March2012",
                  "--chr",strsplit(chr,split="chr")[[1]][2],"--start",range(poslist)[1]-extend,"--end",range(poslist)[2]+extend,
                  "--gwas-cat whole-cat_significant-only",
                  "--plotonly",
                  "--prefix GWAS",
                  paste("signifLine=",-log10(10^-5),sep=""),
                  "signifLineColor=blue")
    system(command)
    
    #make the plot for sQTL result
    command=paste(locuszoom,"--metal input4plot_sqtl.txt",
                  "--pop EUR --build hg19 --source 1000G_March2012",
                  "--chr",strsplit(chr,split="chr")[[1]][2],"--start",range(poslist)[1]-extend,"--end",range(poslist)[2]+extend,
                  "--gwas-cat whole-cat_significant-only",
                  "--plotonly",
                  "--prefix sQTL",
                  paste("signifLine=",-log10(10^-5),sep=""),
                  "signifLineColor=blue")
    system(command)
    
    print(paste(i,"plot successfully!"))
  }
}




##############################
#code to add customized label#
##############################
#testlabel="sQTL.p.MAF_GWAS.original.p.MAF_sQTL.MAF"
#i=11

#label=matrix(,0,3)
#colnames(label)=c("snp",	"string",	"color")
##add sQTL SNP if it is within the range
##sQTLrsID=rownames(input4plot.sqtl)[which(input4plot.sqtl[,"P-value"]==min(input4plot.sqtl[,"P-value"]))]
#sQTLrsID=as.character(gwasinfo[which(gwasinfo[,"sQTL.significant.region"]=="Cortex"),"sQTL.rsID"])      #this sQTL SNP has the same p value with the default picked on (there are a couple of SNPs with the same smallest p value)
#newrow=c(sQTLrsID,"sQTL","blue")
#label=rbind(label,newrow)
#setwd(outputpath)
#write.table(label,"locuszoom_SNP_label.txt",sep="\t",quote=F,row.names=F)

#command=paste(locuszoom,"--metal input4plot_sqtl.txt --denote-markers-file locuszoom_SNP_label.txt",
#              "--pop EUR --build hg19 --source 1000G_March2012",
#              "--chr",strsplit(chr,split="chr")[[1]][2],"--start",range(poslist)[1]-extend,"--end",range(poslist)[2]+extend,
#              "--gwas-cat whole-cat_significant-only",
#              "--plotonly",
#              "--prefix sQTL",
#              paste("signifLine=",-log10(10^-5),sep=""),
#              "signifLineColor=blue")
#system(command)


