#!/bin/sh

###
# tells the cluster to run the job from the current working directory
#$ -cwd

# tells the cluster what info to inform the user through email
#$ -m beas

# tells the cluster the user's email
#$ -M yidazhan@mail

# tells the cluster to use 1024MB of memory, and 1 hours of run time
#$ -l h_data=12G,h_rt=24:00:00,highp

###

code_folder="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/12_region_ubiquitous_vs_specific_splicing_event_heatmap_radar_plot"

# this is most important part, it tells the machine that I want to use R
source /u/local/Modules/default/init/modules.sh
module load R/3.2.3

# run the R code
Rscript $code_folder/Region_ubiquitous_vs_specific_splicing_event_heatmap_radar_plot_step_2_collect_values_for_heatmap_highp.R
