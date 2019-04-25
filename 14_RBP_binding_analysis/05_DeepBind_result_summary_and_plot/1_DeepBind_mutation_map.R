job <- as.numeric(Sys.getenv("SGE_TASK_ID"))
args <- commandArgs(TRUE)
splicetype="SE"      #type of alternative splicing, e.g., SE, A3SS, A5SS, MXE, IR
type="pvalue"       #pvalue or permutation
windowsize=20

library(NMF)
library(RColorBrewer)
require(ggplot2)
require(ggseqlogo)
require(cowplot)
#library(reshape)
require(scales)
library(reshape2)
library(gridExtra)
library(grid)
library(ComplexHeatmap)
library(circlize)
library(lemon)

#get all sQTL exons
rootinput="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/summary/logit/JC"
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
sQTLexon=c()
inputpath=paste(rootinput,splicetype,sep="/")
setwd(inputpath)
for (i in 1:length(brainregionlist)){
  temp=read.table(paste("all_info_exon_sQTL_GWAS_disease_in_",brainregionlist[i],"_logit_JC_",splicetype,"_",type,".txt",sep=""),sep="\t",header=T)
  sQTLexon=c(sQTLexon,as.character(temp[,"exon_full_ID"]))
}
sQTLexon=unique(sQTLexon)
shortIDlist=rep(NA,length(sQTLexon))
for (e in 1:length(sQTLexon)){
  shortIDlist[e]=paste("SE",strsplit(sQTLexon[e],split="\\|")[[1]][1],sep="_")
}

#read in the RBP information
RBPdb="/u/nobackup/yxing/PROJECT/yidazhan/research/software/deepbind/db/db.tsv"
RBPtable=read.table(RBPdb,sep="\t",header=T)
subRBPtable=subset(RBPtable,RBPtable[,"Species"]=="Homo sapiens")
subRBPtable=subset(subRBPtable,subRBPtable[,"Type"]=="RBP")

#read in the joblist
rootoutput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.3_motif_analysis/individual_exon_binding_peaks_scan/DeepBind/deepbind_input/",splicetype,"/",type,sep="")
setwd(rootoutput)
uniquejoblist=read.table(paste(splicetype,"_",type,"_uniquejoblist.txt",sep=""),sep="\t")

#get the information for the current job
shortID=as.character(uniquejoblist[job,"Exon"])
br=gsub(" ", "", as.character(uniquejoblist[job,"Brain.region"]), fixed = TRUE)
snpid=as.character(uniquejoblist[job,"SNP"])
rbpid=as.character(uniquejoblist[job,"RBP"])

#generate output folder
output=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.3_motif_analysis/individual_exon_binding_peaks_scan/DeepBind/deepbind_plot/",splicetype,"/",type,sep="")
outputpath=paste(output,paste(strsplit(sQTLexon[which(shortIDlist %in% shortID)],split="\\|")[[1]],collapse=","),paste(snpid,rbpid,sep="~"),sep="/")

#############################
#read in the deepbind result#
#############################
input=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.3_motif_analysis/individual_exon_binding_peaks_scan/DeepBind/deepbind_output/",splicetype,"/",type,sep="")
inputpath=paste(input,paste(strsplit(sQTLexon[which(shortIDlist %in% shortID)],split="\\|")[[1]],collapse=","),paste(snpid,rbpid,sep="~"),sep="/")
setwd(inputpath)

dbscore=read.table("dbscore.txt",header=T,sep="\t",check.names=F)

#get reference and alternative allele
allele0=strsplit(snpid,split="_")[[1]][3]     #reference allele, the value is 0
allele1=strsplit(snpid,split="_")[[1]][4]     #alternative allele, the value is 1
if (allele0==colnames(dbscore)[windowsize+1]){          #the SNP is on the same strand with the exon sequence
  ref=allele0
  alt=allele1
}else{                                        #the SNP is on the opposite strand of the exon sequence
  ref=chartr("ATGC","TACG",allele0)           #reverse complement
  alt=chartr("ATGC","TACG",allele1)
}


######
#plot#
######
if (nchar(ref)==1 && nchar(alt)==1 && dbscore[alt,c(windowsize+1)]!=0){          #we only plot for single nucleotide mutation with binding changing ability    
  
  #1. heatmap: rowname=c("> A","> C","> G","> T"), no colname, no legend (keep at least one legend for the final figure), no clustering dendrogram
  data2plot=as.matrix(dbscore)
  rownames(data2plot)=paste(">",rownames(data2plot))
  
  start=range(as.numeric(as.matrix(data2plot)))[1]
  end=range(as.numeric(as.matrix(data2plot)))[2]
  
  originalsequence=colnames(dbscore)
  colnames(data2plot)=seq(1:length(originalsequence))    #we need to change ATCG to 1-41 since the column names cannot be duplicated
  data2plot=melt(data2plot)
  colnames(data2plot)=c("row.name","col.name","value")
  temp=strsplit(sQTLexon[which(shortIDlist==shortID)],split="\\|")[[1]]
  title=paste(paste(temp[3:7],collapse="_"),as.character(uniquejoblist[job,"label"]),sep="\n")
  
  if (start<0 && end>0){         #if we have both positive and negative value
    pals=rev(brewer.pal(11,"RdBu"))      #11 colors
    #we generate the corresponding values of the 11 colors
    numbin=5
    stepsize_minus=abs(start)/numbin
    stepsize_plus=abs(end)/numbin
    breakslist=c(start,start+stepsize_minus,start+stepsize_minus*2,start+stepsize_minus*3,start+stepsize_minus*4,
                 0,
                 end-stepsize_plus*4,end-stepsize_plus*3,end-stepsize_plus*2,end-stepsize_plus,end)
    vals=rescale(breakslist)      #this is the position of the 11 colors (relative position)
    hm <- ggplot(data = data2plot, aes(x = factor(col.name), y = row.name, fill = value)) + geom_tile() + 
      scale_fill_gradientn(name="", colours=pals, values=vals, breaks = breakslist) 
  }else if (start>=0 && end>=0){      #if all the values are positive
    hm <- ggplot(data = data2plot, aes(x = factor(col.name), y = row.name, fill = value)) + geom_tile() + 
      scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0) 
  }else if (start<=0 && end<=0){      #if all the values are negative
    hm <- ggplot(data = data2plot, aes(x = factor(col.name), y = row.name, fill = value)) + geom_tile() + 
      scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0) 
  }
  sequence=colnames(dbscore)
  sequence[windowsize+1]=tolower(sequence[windowsize+1])       #we label the sQTL SNP in lower case letter
  #The order of edges for plot.margin is unit(c(top, right, bottom, left), units) 
  hm_margins <- c(-1.5,1,-0.5,1)      
  hm <- hm + theme_gray() + 
    coord_equal() + #force the cell to be square
    scale_x_discrete(breaks=as.character(seq(1:length(originalsequence))),     #we change the column name back to ATCG
                     labels=sequence) +            
    theme(#legend.position = "left", legend.direction = "vertical",
          #legend.title = element_text(size = 15), 
          #legend.key.size = unit(1,"cm"),
          #legend.text = element_text(size = 7), 
          #legend.text = element_blank(),
          #panel.grid.major = element_blank(), 
          #panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position="none",
          plot.margin=unit(hm_margins, "cm"),
          #axis.text.x = element_text(hjust = 1),
          text = element_text(size=14, face="bold"), 
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))
  #print(hm)
  #dev.off()
  
  
  #2. sequence logo: column-wise sum (change negative numbers -> positive), the height is propotional to the column-wise sum
  #pdf("sequence logo.pdf")
  sequencelogoheight=apply(abs(dbscore),2,sum)
  custom_mat = matrix(0,dim(dbscore)[1],dim(dbscore)[2])
  rownames(custom_mat)=rownames(dbscore)
  colnames(custom_mat)=colnames(dbscore)
  
  for (i in 1:dim(custom_mat)[2]){          
    custom_mat[colnames(dbscore)[i],i]=sequencelogoheight[i]
  }
  # Generate sequence logo
  plogo_margins <- c(0.5,0.42,-1.3,1.15)        
  plogo=ggseqlogo(custom_mat, method='custom', seq_type='rna',stack_width = 0.88,    #stack_width is the white space between letters, smaller values mean larger space (but this won't change the overall width of the figure)
                  font="helvetica_bold") + ylab('') + ggtitle(title) +
    theme(axis.text.x = element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          plot.margin=unit(plogo_margins, "cm"))
  #print(plogo)
  #dev.off()
  
  #put everything together
  setwd(outputpath)
  pdf(paste(uniquejoblist[job,"label"],"_modified.pdf",sep=""),height=3, width=6)
  grid.arrange(plogo, hm, ncol=1)
  dev.off()
}




