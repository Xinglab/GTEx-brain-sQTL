#!/bin/sh

###
# tells the cluster to run the job from the current working directory
#$ -cwd

# tells the cluster what info to inform the user through email
#$ -m beas

# tells the cluster the user's email
#$ -M yidazhan@mail

# tells the cluster to use 1024MB of memory, and 1 hours of run time
#$ -l h_data=18G,h_rt=24:00:00

#$ -V

#$ -e ./submit_12_deepbind_region_dependent_sQTL_RBP_plot_part_2_log/
#$ -o ./submit_12_deepbind_region_dependent_sQTL_RBP_plot_part_2_log/

# tells the cluster to run on 100 machines     
#$ -t 1-23:1                                ###change###


# this is most important part, it tells the machine that I want to use R
source /u/local/Modules/default/init/modules.sh
module load R/3.2.3

# run the R code
code_folder="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.3_motif_analysis/individual_exon_binding_peaks_scan/DeepBind"


splicetype="A5SS"     ###change###
type="pvalue"       ###change###
chunksize=100     
Rscript $code_folder/2_DeepBind_region_dependent_sQTL_RBP_plot.R $splicetype $type $chunksize

#SE: 20028 sigpairs   (1-201:1)
#A3SS: 2524 sigpairs  (1-26:1)
#A5SS: 2256 sigpairs  (1-23:1)
