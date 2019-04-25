#!/bin/sh

###
# tells the cluster to run the job from the current working directory
#$ -cwd

# tells the cluster what info to inform the user through email
#$ -m beas

# tells the cluster the user's email
#$ -M yidazhan@mail

# tells the cluster to use 1024MB of memory, and 1 hours of run time
#$ -l h_data=8G,h_rt=16:00:00

#$ -e ./10_sQTL_result_plot_SE_A3SS_A5SS/
#$ -o ./10_sQTL_result_plot_SE_A3SS_A5SS/
###

code_folder="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/scripts"

# this is most important part, it tells the machine that I want to use R
source /u/local/Modules/default/init/modules.sh
module load R/3.2.3

# run the R code
Rscript $code_folder/2_sQTL_result_all_event_summary.R
