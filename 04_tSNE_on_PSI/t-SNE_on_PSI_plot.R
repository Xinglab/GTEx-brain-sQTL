####################################
#read in the annotation information#
####################################
label="logit"
phenotypepath="./03_Get_sample_annotation/example_output"
setwd(phenotypepath)
age=read.table("agetable_brain.txt",sep="\t",header=T)           
gender=read.table("gendertable_brain.txt",sep="\t",header=T)          #current gender is numeric data, not categorical data. It needs to be changed to categorical data for further model fitting
tempgender=gender
gender[1,which(as.numeric(as.matrix(tempgender)[1,])==1)]="male"
gender[1,which(as.numeric(as.matrix(tempgender)[1,])==2)]="female"

brainregion=read.table("brain_region_table_brain.txt",sep="\t",header=T)
batch=read.table("batch_effect_table_brain.txt",sep="\t",header=T)
rownames(batch)=c("collection site",
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
phenotypeinput="./04_tSNE_on_PSI/example_input"
setwd(phenotypeinput)
phenotypename="phs000424.v7.pht002742.v7.p2.c1.GTEx_Subject_Phenotypes.GRU.txt"
phenotype=read.table(phenotypename,sep="\t",skip=10,header=T,fill=TRUE,quote="")

annotationinput="./03_Get_sample_annotation/example_input"
setwd(annotationinput)
annotation=read.csv("gtex_v7_brain.csv",header=T)
race=rep(NA,dim(gender)[2])
for (i in 1:dim(gender)[2]){
  SRRID=colnames(gender)[i]
  sampleID=annotation[which(annotation[,"Run_s"]==SRRID),"submitted_subject_id_s"]
  race[i]=phenotype[which(phenotype[,"SUBJID"]==as.character(sampleID)),"RACE"]
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
library("plotrix")
library(ggsci)
library(scales)
library(ggplot2)
theme_ydz <- theme(#plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
                     plot.title = element_text(size = rel(1.3), vjust = 2, hjust = 0.5, lineheight = 0.8),
                     
                     # axis
                     #axis.title = element_text(size = rel(1.2)),
                     #axis.title = element_blank(),
                     axis.title.x = element_text(face="bold", size=16),
                     axis.title.y = element_text(face="bold", size=16, angle=90),
                     axis.text = element_text(size = rel(1.1)),
                     axis.text.x = element_text(hjust = 0.5, vjust = 0, size=14),
                     axis.text.y = element_text(vjust = 0.5, hjust = 0, size=14),
                     axis.line = element_line(colour = "black"),
                     
                     # ticks
                     axis.ticks = element_line(colour = 'black'),
                     
                     # legend
                     legend.position="none",
                     
                     # strip
                     # strip.background=element_rect(fill=NA, color=NA),
                     strip.text=element_text(size = rel(1.3)),
                     
                     # panel
                     #panel.background = element_rect(fill = NA, colour = "black", size = 0.7),
                     panel.background= element_blank(), 
                     #panel.border = element_rect(colour = "grey", fill=NA, size=1),
                     # panel.grid.major = element_line(colour = "grey", size = 0.05),
                     
                     aspect.ratio=1,
                     complete = T)

#change brain region name into shorter ones
BRconversion=matrix(NA,length(unique(as.character(as.matrix(brainregion)))),2)
BRconversion[,1]=unique(as.character(as.matrix(brainregion)))
BRconversion[,2]=c("Cerebellum","Anterior cingulate cortex","Hypothalamus","Frontal Cortex",
                   "Caudate","Amygdala","Spinal cord","Nucleus accumbens","Substantia nigra",
                   "Cortex","Cerebellar Hemisphere","Putamen","Hippocampus")                 
rownames(BRconversion)=BRconversion[,1]
colnames(BRconversion)=c("original","new")
short_BR=rep(NA,dim(brainregion)[2])
for(i in 1:dim(brainregion)[2]){
  short_BR[i]=BRconversion[as.character(brainregion[1,i]),"new"]
}
short_BR=t(as.data.frame(short_BR))

splicelist=c("SE","MXE","A3SS","A5SS","RI")
countlist=c("JC","JCEC")

for (i in 1:length(splicelist)){
  for (j in 1:length(countlist)){
    splicetype=splicelist[i]
    counttype=countlist[j]
    
    ######################
    #read in t-SNE result#
    ######################
    inputpath="./04_tSNE_on_PSI/example_output"
    setwd(inputpath)
    sampleplot=read.table(paste(splicetype,"_",counttype,"_tsne_sample_",label,"_PSI.txt",sep=""),sep="\t",header=T)
    
    ######
    #plot#
    ######
    outputpath="./04_tSNE_on_PSI/example_output"
    setwd(outputpath)
    
    filename=paste(label,"_",splicetype,"_",counttype,"_t-SNE_plot.pdf",sep="")
    pdf(filename,width = 8, height = 6)
    #age
    COLOR="age.group"    #this will control the legend name
    temp=as.data.frame(sampleplot)
    age.group=as.factor(as.character(as.matrix(age[1,])))
    p.age=ggplot(temp,aes_string(x=temp$V1, y=temp$V2, color=COLOR)) + 
      geom_point(size=2) + 
      scale_color_npg() + 
      theme_ydz + 
      ggtitle(paste(splicetype,counttype,"\n")) + 
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
      ggtitle(paste(splicetype,counttype,"\n")) + 
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
      ggtitle(paste(splicetype,counttype,"\n")) +
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
        ggtitle(paste(paste(splicetype,counttype),rownames(batch)[b],sep="\n")) +
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
      ggtitle(paste(splicetype,counttype,"\n")) + 
      labs(x ="\ntSNE-1", y = "tSNE-2\n")
    print(p.race)
    
    dev.off()
    
    
    
    #we put brain region, age and gender information into the same plot
    filename=paste(label,"_",splicetype,"_",counttype,"_t-SNE_plot_BR_AGE_GENDER.pdf",sep="")
    pdf(filename,width = 8, height = 6)
    temp=as.data.frame(sampleplot)
    COLOR="brain.region"    #this will control the legend name
    brain.region=as.factor(as.character(as.matrix(short_BR[1,])))
    p=ggplot(temp,aes_string(x=temp$V1, y=temp$V2, color = COLOR)) + 
      geom_point(size=1.5, aes(shape=age.group), stroke = 1) + 
      scale_color_manual(values=cbPalette) + 
      scale_shape(solid = FALSE) + 
      theme_ydz + 
      ggtitle(paste(splicetype,counttype,"\n")) + 
      labs(x ="\ntSNE-1", y = "tSNE-2\n")
    print(p)
    dev.off()
  }
}



