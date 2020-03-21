#The purpose of this code is to show the same information in a circular bar plot instead of a CDF or violin plot
#This code only has plot section. The data comes from the result of SNP_to_SS_distance_vs_brain_region_specificity - 5 categories.R
#reference: https://www.r-graph-gallery.com/297-circular-barplot-with-groups/
#The code needs to run using module load R/3.4.2

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

splicetypelist=c("SE","A3SS","A5SS")
typelist=c("pvalue","permutation")
PSItype="logit"
counttype="JC"

outputpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.2_sQTL_SNP_annotation/result/SNP_to_SS_distance_vs_brain_region_specificity"
summaryinput="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/summary"
sQTLinput="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/sQTL_run"

library(ggplot2)
#library(reshape)
library(tidyverse)
library(viridis)

pick_category=function(list_of_category,order_of_importance){      
  #given a list of categories and the order we want to use, this function will return one category from the list based on the order of importance
  reorderedlist=list_of_category[order(match(list_of_category, order_of_importance))]
  return(reorderedlist)
}

importance_order=c("dinucleotide", "SS", "exon", "<=300bp", ">300bp")

for (splicetype in splicetypelist){
  for (type in typelist){
    exoninfopath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/input_splicing",
                       PSItype,counttype,splicetype,sep="/")
    setwd(exoninfopath)
    exoninfo=read.table(paste("exon_info.fromGTF.",splicetype,".txt",sep=""),sep="\t",header=T)
    rownames(exoninfo)=exoninfo[,"ID"]
    exonshorterID=rownames(exoninfo)
    
    #1. get all sQTL exons
    sQTLexon=c()
    inputpath=paste(summaryinput,PSItype,counttype,splicetype,sep="/")
    setwd(inputpath)
    for (i in 1:length(brainregionlist)){
      temp=read.table(paste("all_info_exon_sQTL_GWAS_disease_in_",brainregionlist[i],"_",PSItype,"_",counttype,"_",splicetype,"_",type,".txt",sep=""),sep="\t",header=T)
      sQTLexon=c(sQTLexon,as.character(temp[,"exon_full_ID"]))
    }
    sQTLexon=unique(sQTLexon)
    
    #2. get a more stringent list of sQTL events
    cutofftype="FDR10"                                                                 #######change########
    if (cutofftype=="FDR1" || cutofftype=="FDR5"){
      setwd(paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/12_region_ubiquitous_vs_specific_splicing_event_heatmap_radar_plot/Region_ubiquitous_vs_specific_splicing_event_heatmap_radar_plot_step_5_pick_SNP_and_collect_values_for_stacked_bar_plot_different_cutoffs/2_select_top_SNP/result/",cutofftype,sep=""))
    }
    if (cutofftype=="FDR10"){
      setwd("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/12_region_ubiquitous_vs_specific_splicing_event_heatmap_radar_plot/Region_ubiquitous_vs_specific_splicing_event_heatmap_radar_plot_step_3_pick_SNP_and_collect_values_for_heatmap_more_stringent_version/single_job_run/result/200kb")
    }
    pvmatrix200k=read.table("pvmatrix_200kb.txt",sep="\t",header=T)
    newsQTLexon=rep(NA,dim(pvmatrix200k)[1])
    for (i in 1:length(newsQTLexon)){
      newsQTLexon[i]=strsplit(rownames(pvmatrix200k)[i],split="~")[[1]][1]
    }
    
    #3. get other information based on the new list of sQTL events
    setwd("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.2_sQTL_SNP_annotation/result/SNP_to_SS_distance_vs_brain_region_specificity/closest_significant_SNP_categorized_distance/5_categories")
    result=read.table(paste(splicetype,"_",type,"_exonclass.txt",sep=""),sep="\t",header=T,check.names=F)
    result=result[newsQTLexon,]
    
    exonclass=as.matrix(result[,1:13])
    num.sig.region=as.numeric(as.matrix(result[,14]))

    dis_2_SS=rep(NA,length(newsQTLexon))      #the best category for each sQTL exon
    for (e in 1:length(newsQTLexon)){
      dis_2_SS[e]=pick_category(exonclass[e,],importance_order)[1]
    }
    
    #4. make the plot
    distance=c(rep("dinucleotide",13),rep("SS",13), rep("exon",13), rep("<=300bp",13), rep(">300bp",13))
    numbr=rep(c(1:13),5)
    value=numevent=rep(NA,length(numbr))
    for (i in 1:length(distance)){
      category=distance[i]
      number.br=numbr[i]
      #absolute value
      numevent[i]=length(which(num.sig.region[which(dis_2_SS==category)]==number.br))
      #proportion
      value[i]=(length(which(num.sig.region[which(dis_2_SS==category)]==number.br))/length(which(dis_2_SS==category)))*100
    }
    
    #data2plot=cbind(paste(distance,numbr,sep="_"),distance,value)
    data2plot=cbind(numbr,distance,value,numevent)
    df=as.data.frame(data2plot)
    colnames(df)=c("num.sig.region","distance","num.events.prop","num.events")
    df$distance=factor(df$distance, levels = importance_order) 
    df$num.events.prop=as.numeric(as.character(df$num.events.prop))
    df$num.events=as.numeric(as.character(df$num.events))
    
    ###circular bar plot###
    # Set a number of 'empty bar' to add at the end of each group
    empty_bar=3
    to_add = data.frame( matrix(NA, empty_bar*nlevels(df$distance), ncol(df)) )
    colnames(to_add) = colnames(df)
    to_add$distance=rep(levels(df$distance), each=empty_bar)
    df=rbind(df, to_add)
    df=df %>% arrange(distance)
    df$id=seq(1, nrow(df))
    
    # Get the name and the y position of each label
    label_data=df
    number_of_bar=nrow(label_data)
    angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
    label_data$hjust<-ifelse( angle < -90, 1, 0)
    label_data$angle<-ifelse(angle < -90, angle+180, angle)
    
    # prepare a data frame for base lines
    base_data=df %>% 
      group_by(distance) %>% 
      summarize(start=min(id), end=max(id) - empty_bar) %>% 
      rowwise() %>% 
      mutate(title=mean(c(start, end)))
    
    # prepare a data frame for grid (scales)
    grid_data = base_data
    grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
    grid_data$start = grid_data$start - 1
    #grid_data=grid_data[-1,]   #this removes the circular grid line on the proportion axis 
    
    # Make the plot
    setwd(outputpath)
    pdf(paste(splicetype,type,"SNP_to_SS_distance_vs_brain_region_specificity_circular_bar_plot_",cutofftype,".pdf"),height=8,width=8)
    p = ggplot(df, aes(x=as.factor(id), y=num.events.prop, fill=distance)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
      
      geom_bar(aes(x=as.factor(id), y=num.events.prop, fill=distance), stat="identity", alpha=0.5) +
      
      scale_fill_manual(values=c("#FDC086", "#BEAED4", "#0F0C73","#26908E","#626565")) +              #color of each group of bars
      
      # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
      geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
      geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
      geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
      geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
      
      # Add text showing the value of each 100/75/50/25 lines
      annotate("text", x = rep(max(df$id),4), y = c(20, 40, 60, 80), label = c("20%", "40%", "60%", "80%") , color="black", size=3 , angle=0, fontface="bold", hjust=1) +
      
      geom_bar(aes(x=as.factor(id), y=num.events.prop, fill=distance), stat="identity", alpha=0.5) +
      ylim(-100,120) +
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(rep(-1,4), "cm") 
      ) +
      coord_polar() + 
      #geom_text(data=label_data, aes(x=id, y=num.events.prop+10, label=num.sig.region, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=3, angle= label_data$angle, inherit.aes = FALSE ) +     #text label of each bar (num.sig.region on top of each bar)
      geom_text(data=label_data, aes(x=id, y=-15, label=num.sig.region, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=3, angle= label_data$angle, inherit.aes = FALSE ) +     #text label of each bar (num.sig.region below each bar)
      geom_text(data=label_data, aes(x=id, y=num.events.prop+10, label=num.events, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=3, angle= label_data$angle, inherit.aes = FALSE ) +     #text label of each bar (num.events on top of each bar)
      
      # Add base line information
      geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +     #the black curve below each group of bars
      #geom_text(data=base_data, aes(x = title, y = -18, label=distance), colour = rep("black",5), alpha=rep(0.8,5), size=rep(4,5), fontface=rep("bold",5), inherit.aes = FALSE)      #the text label for each group of bars (below each group)
      geom_text(data=base_data, aes(x = title, y = 100, label=distance), colour = rep("black",5), alpha=rep(0.8,5), size=rep(4,5), fontface=rep("bold",5), inherit.aes = FALSE)      #the text label for each group of bars (outside each group)
    
    print(p)
    dev.off()
  }
}








