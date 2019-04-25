#!/bin/sh

###
# tells the cluster to run the job from the current working directory
#$ -cwd

# tells the cluster what info to inform the user through email
#$ -m a

# tells the cluster the user's email
#$ -M yidazhan@mail

# tells the cluster to use 1024MB of memory, and 1 hours of run time
#$ -l h_data=4G,h_rt=3:00:00

#$ -V
 
#$ -e /dev/null           #no log file
#$ -o /dev/null           #no log file   

# tells the cluster to run on 100 machines     
#$ -t 1-15198:1                                 
###

splicetype="A5SS"  ###change###
type="pvalue"    ###change###

# this is most important part, it tells the machine that I want to use R
source /u/local/Modules/default/init/modules.sh
module load R/3.4.0

# run the R code
code_folder="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.3_motif_analysis/individual_exon_binding_peaks_scan/DeepBind"
Rscript $code_folder/1_DeepBind_mutation_map.R $splicetype $type

#SE: 75001-103428:1
#A3SS: 1-16320:1
#A5SS: 1-15198:1
