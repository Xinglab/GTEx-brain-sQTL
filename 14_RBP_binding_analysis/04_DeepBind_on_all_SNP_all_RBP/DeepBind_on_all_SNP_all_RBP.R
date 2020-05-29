job <- as.numeric(Sys.getenv("SGE_TASK_ID"))
args <- commandArgs(TRUE)
splicetype=args[1]
type=args[2]

windowsize=20

#get all sQTL exons
rootinput="/path/to/summary/logit/JC"
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
RBPdb="/path/to/deepbind/db/db.tsv"
RBPtable=read.table(RBPdb,sep="\t",header=T)
subRBPtable=subset(RBPtable,RBPtable[,"Species"]=="Homo sapiens")
subRBPtable=subset(subRBPtable,subRBPtable[,"Type"]=="RBP")

#read in the joblist
rootoutput=paste("/path/to/DeepBind/deepbind_input/",splicetype,"/",type,sep="")
setwd(rootoutput)
uniquejoblist=read.table(paste(splicetype,"_",type,"_uniquejoblist.txt",sep=""),sep="\t")

#get the information for the current job
shortID=as.character(uniquejoblist[job,"Exon"])
br=gsub(" ", "", as.character(uniquejoblist[job,"Brain.region"]), fixed = TRUE)
snpid=as.character(uniquejoblist[job,"SNP"])
rbpid=as.character(uniquejoblist[job,"RBP"])

#get genome sequence
functionfile="/path/to/rbp_map_functions.py"
fasta="/path/to/hg19/FASTA/file"
library(reticulate)
source_python(functionfile)
genome=read_genome(fasta)
library(Biostrings)

#generate output folder
output=paste("/output/path/DeepBind/deepbind_output/",splicetype,"/",type,sep="")
path=paste(output,paste(strsplit(sQTLexon[which(shortIDlist %in% shortID)],split="\\|")[[1]],collapse=","),paste(snpid,rbpid,sep="~"),sep="/")
command=paste("mkdir -p",path)
system(command)
setwd(path)

#output the list of significant brain regions
write.table(br,"sig_brain_region.txt",sep="\t",quote=F,row.names=F,col.names=F)

##########################
#run deepbind calculation#
##########################
setwd(path)

#generate .id file for the current RBP
r=which(subRBPtable[,"ID"] %in% strsplit(rbpid,split="_")[[1]][1])
idfile=paste(subRBPtable[r,"ID"],"#",subRBPtable[r,"Protein"], paste("(",subRBPtable[r,"Experiment"],")",sep=""),sep=" ")
write.table(idfile,"RBP.ids",sep="\t",row.names=F,col.names=F,quote=F)

#get the position of the SNP
chr=paste("chr",strsplit(snpid,split="_")[[1]][1],sep="")
pos=as.numeric(strsplit(snpid,split="_")[[1]][2])
strand=strsplit(sQTLexon[which(shortIDlist %in% shortID)],split="\\|")[[1]][5]

#expand from the SNP to two sides for 20 bp (41 bp in total)
region_start=pos-windowsize
region_end=pos+windowsize

#generate a 4x41 matrix for score + score just for wildtype sequence + score just for mutant sequence
dbscore=dbscore_WT=dbscore_WT0=dbscore_WT00=dbscore_MUT=dbscore_MUT0=dbscore_MUT00=countmatrix=matrix(0,nrow=4,ncol=2*windowsize+1)
rownames(dbscore)=rownames(dbscore_WT)=rownames(dbscore_WT0)=rownames(dbscore_WT00)=rownames(dbscore_MUT)=rownames(dbscore_MUT0)=rownames(dbscore_MUT00)=c("A","T","C","G")
colnames(dbscore)=colnames(dbscore_WT)=colnames(dbscore_WT0)=colnames(dbscore_WT00)=colnames(dbscore_MUT)=colnames(dbscore_MUT0)=colnames(dbscore_MUT00)=strsplit(fetch_seq(genome,chr,region_start-1,region_end,strand),split="")[[1]]
n=2*windowsize+1
m=windowsize

#calculate motif score matrix within a given sliding window
windowscore_calculator=function(scorematrix,scorematrix_WT,scorematrix_WT0,scorematrix_WT00,scorematrix_MUT,scorematrix_MUT0,scorematrix_MUT00,wildtypesequence,psscore,windownumber){
  substitute=c("A","T","C","G")
  for (n in 1:nchar(wildtypesequence)){   #for each position in the sequence
    for (t in 1:length(substitute)){      #for each nucleotide that we can substitute into
      #replace the nth nucleotide with the tth nucleotide
      mutseq=as.character(replaceLetterAt(DNAString(wildtypesequence),n,substitute[t]))
      #put the mutant sequence into a .seq file
      mutseqfilename=paste("mutant_seq_",windownumber,"_",n,"_",substitute[t],".seq",sep="")
      write.table(mutseq,mutseqfilename,sep="\t",col.names=F,row.names=F,quote=F)
      #run deepbind to get the score p(s_hat)
      command=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/software/deepbind/deepbind RBP.ids", mutseqfilename)
      ps_hat=as.numeric(system(command,intern=T)[2])
      
      #calculate score (original score)
      scorematrix_WT[t,n]=psscore
      scorematrix_MUT[t,n]=ps_hat
      scorematrix[t,n]=(ps_hat-psscore)*max(0,psscore,ps_hat)
      
      #0: we only keep those positions that contribute to the deepbind score (at least one positive score)
      if (max(0,psscore,ps_hat)>0){       #if we have at least on positive score
        scorematrix_WT0[t,n]=psscore
        scorematrix_MUT0[t,n]=ps_hat
      }else{                              #if we don't have positive scores
        scorematrix_WT0[t,n]=0
        scorematrix_MUT0[t,n]=0
      }
      
      #00: we only keep positive scores (negative means no binding, so we just use 0)
      scorematrix_WT00[t,n]=max(0,psscore)
      scorematrix_MUT00[t,n]=max(0,ps_hat)
    }
  }
  rownames(scorematrix)=rownames(scorematrix_WT)=rownames(scorematrix_WT0)=rownames(scorematrix_WT00)=rownames(scorematrix_MUT)=rownames(scorematrix_MUT0)=rownames(scorematrix_MUT00)=substitute
  colnames(scorematrix)=colnames(scorematrix_WT)=colnames(scorematrix_WT0)=colnames(scorematrix_WT00)=colnames(scorematrix_MUT)=colnames(scorematrix_MUT0)=colnames(scorematrix_MUT00)=strsplit(wildtypesequence,split="")[[1]]
  return(list(SM=scorematrix,SM_WT=scorematrix_WT,SM_WT0=scorematrix_WT0,SM_WT00=scorematrix_WT00,SM_MUT=scorematrix_MUT,SM_MUT0=scorematrix_MUT0,SM_MUT00=scorematrix_MUT00))
}

#sliding window
for (i in 1:(n+m-1)){
  print(i)
  window_start=region_start-windowsize+i
  window_end=window_start+windowsize-1
  #generate a 4x20 matrix
  windowscore=matrix(0,nrow=4,ncol=windowsize)
  windowscore_WT=windowscore_WT0=windowscore_WT00=matrix(0,nrow=4,ncol=windowsize)
  windowscore_MUT=windowscore_MUT0=windowscore_MUT00=matrix(0,nrow=4,ncol=windowsize)
  
  #get the wildtype sequence
  wtseq=fetch_seq(genome,chr,window_start-1,window_end,strand)   #start position is not returned using this function so we need to -1 to include the start position
  #put the sequence into a .seq file
  write.table(wtseq,"wildtype_seq.seq",sep="\t",col.names=F,row.names=F,quote=F)
  
  #run deepbind to get the score p(s)
  command="/u/nobackup/yxing/PROJECT/yidazhan/research/software/deepbind/deepbind RBP.ids wildtype_seq.seq"
  ps=as.numeric(system(command,intern=T)[2])
  #fill the 4x20 matrix
  tempscore=windowscore_calculator(windowscore,windowscore_WT,windowscore_WT0,windowscore_WT00,windowscore_MUT,windowscore_MUT0,windowscore_MUT00,wtseq,ps,i)
  windowscore=tempscore$SM
  windowscore_WT=tempscore$SM_WT
  windowscore_WT0=tempscore$SM_WT0
  windowscore_WT00=tempscore$SM_WT00
  windowscore_MUT=tempscore$SM_MUT
  windowscore_MUT0=tempscore$SM_MUT0
  windowscore_MUT00=tempscore$SM_MUT00
  
  #after filling the whole 4x20 matrix, we add the matrix back to the 4x41 matrix
  if (strand=="+"){
    overlap_column_window=which(seq(window_start,window_end) %in% intersect(seq(region_start,region_end),seq(window_start,window_end)))
    overlap_column_db=which(seq(region_start,region_end) %in% intersect(seq(region_start,region_end),seq(window_start,window_end)))
  }
  if (strand=="-"){
    overlap_column_window=rev(dim(windowscore)[2]-which(seq(window_start,window_end) %in% intersect(seq(region_start,region_end),seq(window_start,window_end)))+1)
    overlap_column_db=rev(dim(dbscore)[2]-which(seq(region_start,region_end) %in% intersect(seq(region_start,region_end),seq(window_start,window_end)))+1)
  }
  dbscore[,overlap_column_db]=dbscore[,overlap_column_db]+windowscore[,overlap_column_window]
  
  dbscore_WT[,overlap_column_db]=dbscore_WT[,overlap_column_db]+windowscore_WT[,overlap_column_window]
  dbscore_WT0[,overlap_column_db]=dbscore_WT0[,overlap_column_db]+windowscore_WT0[,overlap_column_window]
  dbscore_WT00[,overlap_column_db]=dbscore_WT00[,overlap_column_db]+windowscore_WT00[,overlap_column_window]
  
  dbscore_MUT[,overlap_column_db]=dbscore_MUT[,overlap_column_db]+windowscore_MUT[,overlap_column_window]
  dbscore_MUT0[,overlap_column_db]=dbscore_MUT0[,overlap_column_db]+windowscore_MUT0[,overlap_column_window]
  dbscore_MUT00[,overlap_column_db]=dbscore_MUT00[,overlap_column_db]+windowscore_MUT00[,overlap_column_window]
  
  countmatrix[,overlap_column_db]=countmatrix[,overlap_column_db]+1
}

#after filling the 4x41 matrix: calculate the average score: devide all the numbers in the matrix by window size 20
dbscore=dbscore/windowsize
write.table(dbscore,"dbscore.txt",sep="\t")

dbscore_WT=dbscore_WT/windowsize
write.table(dbscore_WT,"dbscore_WT.txt",sep="\t")
dbscore_WT0=dbscore_WT0/windowsize
write.table(dbscore_WT0,"dbscore_WT0.txt",sep="\t")
dbscore_WT00=dbscore_WT00/windowsize
write.table(dbscore_WT00,"dbscore_WT00.txt",sep="\t")

dbscore_MUT=dbscore_MUT/windowsize
write.table(dbscore_MUT,"dbscore_MUT.txt",sep="\t")
dbscore_MUT0=dbscore_MUT0/windowsize
write.table(dbscore_MUT0,"dbscore_MUT0.txt",sep="\t")
dbscore_MUT00=dbscore_MUT00/windowsize
write.table(dbscore_MUT00,"dbscore_MUT00.txt",sep="\t")

#remove all the mutant sequences
command="rm mutant_seq_*.seq"
system(command)

      
