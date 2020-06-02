library("plotrix")
library(ggsci)
library(scales)
library(ggplot2)
theme_ydz <- theme(#plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
  plot.title = element_text(size = rel(1.3), vjust = 2, hjust = 0.5, lineheight = 0.8),
  
  # axis
  axis.title.x = element_text(face="bold", size=16),
  axis.title.y = element_text(face="bold", size=16, angle=90),
  axis.text = element_text(size = rel(1.1)),
  axis.text.x = element_text(hjust = 0.5, vjust = 0, size=14),
  axis.text.y = element_text(vjust = 0.5, hjust = 0, size=14),
  axis.line = element_line(colour = "black"),
  
  # ticks
  axis.ticks = element_line(colour = 'black'),
  
  # legend
  legend.title=element_blank(),
  legend.text = element_text(size = rel(1.1),face="bold"),
  legend.position = "right",
  legend.direction = 'vertical',
  legend.background = element_blank(),
  legend.key = element_rect(fill = NA, colour = NA),
  
  # strip
  strip.text=element_text(size = rel(1.3)),
  
  # panel
  panel.background= element_blank(), 

  aspect.ratio=1,
  complete = T)

####################
#read in annotation#
####################
phenotypepath="./03_Get_sample_annotation/example_output"
setwd(phenotypepath)
oriage=read.table("agetable_brain.txt",sep="\t",header=T)           
origender=read.table("gendertable_brain.txt",sep="\t",header=T)          #current gender is numeric data, not categorical data. It needs to be changed to categorical data for further model fitting
tempgender=origender
origender[1,which(as.numeric(as.matrix(tempgender)[1,])==1)]="male"
origender[1,which(as.numeric(as.matrix(tempgender)[1,])==2)]="female"
oribrainregion=read.table("brain_region_table_brain.txt",sep="\t",header=T)


oribatch=read.table("batch_effect_table_brain.txt",sep="\t",header=T)
rownames(oribatch)=c("collection site",
                  "Nucleic Acid Isolation Batch",
                  "Type of nucleic acid isolation batch",
                  "Date of nucleic acid isolation batch",
                  "Genotype or Expression Batch ID",
                  "Date of genotype or expression batch",
                  "Type of genotype or expression batch")
#SMCENTER: Code for BSS collection site
#SMNABTCH: Nucleic Acid Isolation Batch ID, batch when DNA/RNA was isolated and extracted from a sample
#SMNABTCHT: Type of nucleic acid isolation batch
#SMNABTCHD: Date of nucleic acid isolation batch
#SMGEBTCH: Genotype or Expression Batch ID
#SMGEBTCHD: Date of genotype or expression batch
#SMGEBTCHT: Type of genotype or expression batch

#population information
setwd("/input/path/to/GTEx/phenotype/result")
phenotypename="phs000424.v7.pht002742.v7.p2.c1.GTEx_Subject_Phenotypes.GRU.txt"
phenotype=read.table(phenotypename,sep="\t",skip=10,header=T,fill=TRUE,quote="")
setwd("/input/path/to/GTEx/brain/tissue/sample/annotation")
annotation=read.csv("gtex_v7_brain.csv",header=T)

setwd("/input/to/GTEx/processed/gene/expression/data")
IDconversion=read.table("SRRID_2_GTEXID.txt",sep="\t",header=T)
rownames(IDconversion)=IDconversion[,"sampleID"]

######################
#read in t-SNE result#
######################
rootinput="/input/path/to/gene/expression/tSNE/result"
expressiontypelist=c("allgene","allRBP","allSF")
for (expressiontype in expressiontypelist){
  inputpath=paste(rootinput,expressiontype,sep="/")
  setwd(inputpath)
  sampleplot=read.table(paste(expressiontype,"_tsne_sample_log_expr.txt",sep=""),sep="\t")
  
  
  #get the phenotype information of samples in the expr matrix (not all 1409 sampes are in the expr data)
  age=oriage[,as.character(IDconversion[rownames(sampleplot),"SRRlist"])]
  gender=origender[,as.character(IDconversion[rownames(sampleplot),"SRRlist"])]
  brainregion=oribrainregion[,as.character(IDconversion[rownames(sampleplot),"SRRlist"])]
  batch=oribatch[,as.character(IDconversion[rownames(sampleplot),"SRRlist"])]
  
  #change brain region name into shorter ones
  BRconversion=matrix(NA,length(unique(as.character(as.matrix(brainregion)))),2)
  BRconversion[,1]=unique(as.character(as.matrix(brainregion)))
  BRconversion[,2]=c("Amygdala",
                     "Anterior cingulate cortex",
                     "Caudate",
                     "Cerebellar Hemisphere",
                     "Cerebellum",
                     "Cortex",
                     "Frontal Cortex",
                     "Hippocampus",
                     "Hypothalamus",
                     "Nucleus accumbens",
                     "Putamen",
                     "Spinal cord",
                     "Substantia nigra")                 
  rownames(BRconversion)=BRconversion[,1]
  colnames(BRconversion)=c("original","new")
  short_BR=rep(NA,dim(brainregion)[2])
  for(i in 1:dim(brainregion)[2]){
    short_BR[i]=BRconversion[as.character(brainregion[1,i]),"new"]
  }
  short_BR=t(as.data.frame(short_BR))
  
  
  race=rep(NA,dim(gender)[2])
  for (r in 1:dim(gender)[2]){
    SRRID=colnames(gender)[r]
    sampleID=annotation[which(annotation[,"Run_s"]==SRRID),"submitted_subject_id_s"]
    race[r]=phenotype[which(phenotype[,"SUBJID"]==as.character(sampleID)),"RACE"]
  }
  race[which(race==1)]="Asian"
  race[which(race==2)]="African American"
  race[which(race==3)]="White"
  race[which(race==4)]="American Indian"
  race[which(race==98)]="Not reported"
  race[which(race==99)]="Unknown"
  
  
  ######
  #plot#
  ######
  outputpath=inputpath
  setwd(outputpath)
  
  filename=paste(expressiontype,"_t-SNE_plot.pdf",sep="")
  pdf(filename,width = 8, height = 6)
  #age
  COLOR="age.group"    #this will control the legend name
  temp=as.data.frame(sampleplot)
  age.group=as.factor(as.character(as.matrix(age[1,])))
  p.age=ggplot(temp,aes_string(x=temp$V1, y=temp$V2, color=COLOR)) + 
    geom_point(size=2) + 
    scale_color_npg() + 
    theme_ydz + 
    ggtitle(expressiontype) + 
    labs(x ="\ntSNE-1", y = "tSNE-2\n")
  print(p.age)
  
  #gender
  COLOR="gender.group"    #this will control the legend name
  temp=as.data.frame(sampleplot)
  gender.group=as.factor(as.character(as.matrix(gender[1,])))
  p.gender=ggplot(temp,aes_string(x=temp$V1, y=temp$V2, color=COLOR)) + 
    geom_point(size=2) + 
    scale_color_npg() + 
    theme_ydz + 
    ggtitle(expressiontype) + 
    labs(x ="\ntSNE-1", y = "tSNE-2\n")
  print(p.gender)
  
  #brain region    
  cbPalette <- c("#E53537",     #GTEx Artery color/Amygdala
                 "#EEA760",     #GTEx adipose tissue color/ACC
                 "#C9E5C3",     #GTEx pituitary color/Caudate
                 "#9ACB3C",     #GTEx lung color/CH
                 "#028B45",     #GTEx thyroid color/Cerebellum
                 "#FFD923",     #GTEx nerve color/cortex
                 "#F0EC68",     #GTEx brain color/frontal cortex
                 "#6D67AF",     #GTEx muscle color/hippocampus
                 "#7F388D",     #GTEx Heart color/Hypothalamus
                 "#2999D5",     #GTEx Skin color/Nucleus accumbens
                 "#4EC3C7",     #GTEx breast color/Putamen
                 "#CD8ABC",     #GTEx LCL color/spinal cord
                 "#F287B7")     #GTEx blood color/Substantia nigra
  COLOR="brain.region"    #this will control the legend name
  temp=as.data.frame(sampleplot)
  brain.region=as.factor(as.character(as.matrix(short_BR[1,])))
  p.brain=ggplot(temp,aes_string(x=temp$V1, y=temp$V2, color=COLOR)) + 
    geom_point(size=2) + 
    scale_color_manual(values=cbPalette) + 
    theme_ydz + 
    ggtitle(expressiontype) +
    labs(x ="\ntSNE-1", y = "tSNE-2\n")
  print(p.brain)
  
  #batch effect
  for (b in 1:dim(batch)[1]){
    batch_effect=batch[b,]
    COLOR=gsub(" ", ".", rownames(batch)[b], fixed = TRUE)    #this will control the legend name (this is the name of the batch effect, there cannot be space in the name)
    temp=as.data.frame(sampleplot)
    assign(paste(COLOR,"",sep=""),as.factor(as.character(as.matrix(batch_effect))))
    p=ggplot(temp,aes_string(x=temp$V1, y=temp$V2, color=COLOR)) + 
      geom_point(size=2) + 
      theme_ydz + 
      ggtitle(rownames(batch)[b]) +
      labs(x ="\ntSNE-1", y = "tSNE-2\n") + 
      guides(colour=FALSE)  #remove legend because there are too many levels for some of the batch effect
    print(p)
  }
  
  #race
  COLOR="race.group"    #this will control the legend name
  temp=as.data.frame(sampleplot)
  race.group=as.factor(as.character(as.matrix(race)))
  p.race=ggplot(temp,aes_string(x=temp$V1, y=temp$V2, color=COLOR)) + 
    geom_point(size=2) + 
    scale_color_npg() + 
    theme_ydz + 
    ggtitle(expressiontype) + 
    labs(x ="\ntSNE-1", y = "tSNE-2\n")
  print(p.race)
  
  dev.off()
}

