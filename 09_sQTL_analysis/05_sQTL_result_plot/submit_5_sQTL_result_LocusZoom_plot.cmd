#!/bin/sh

###
# tells the cluster to run the job from the current working directory
#$ -cwd

#pass all the environment variables to the script (this is to guarantee that we are using the correct version of python, see below)
#$ -V

# tells the cluster what info to inform the user through email
#$ -m beas

# tells the cluster the user's email
#$ -M yidazhan@mail

# tells the cluster to use 1024MB of memory, and 1 hours of run time
#$ -l h_data=12G,h_rt=48:00:00,highp

#$ -e ./12_sQTL_LocusZoom_plot_log/
#$ -o ./12_sQTL_LocusZoom_plot_log/
###

splicetype="A5SS"        ###change###

counttype="JC"
PSItype="logit"
#PSItype="logit" means the PSI value is after logit transformation + confounder conrrection
#PSItype="original" means the PSI value is before logit transformation and no confounder correction. In this case, we will do logit transformation later in the pipeline
#the idea here is to run the sQTL analysis on data before + after confounder correction so eventually they are both on logit scale


#for logit PSI
splicing_input="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/5_correction_for_technical_confounders/result/"$counttype
#for original PSI
#splicing_input="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/1_get_matrix_from_rMATS/result/"$counttype


summaryinput="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/summary/"$PSItype/$counttype/$splicetype
rootsqtl="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/sQTL_run"
outputpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/plot/"$PSItype/$counttype/$splicetype
IDconvertscript="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/scripts/5_sQTL_result_LocusZoom_plot_rsID_conversion.py"
locuszoom="/u/nobackup/yxing/PROJECT/yidazhan/research/software/locuszoom/bin/locuszoom"

mkdir -p $outputpath

# this is most important part, it tells the machine that I want to use R
source /u/local/Modules/default/init/modules.sh
module load R

# run the R code
code_folder="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/scripts"

type="pvalue"
Rscript $code_folder/5_sQTL_result_LocusZoom_plot.R $splicetype $counttype $PSItype $splicing_input $summaryinput $rootsqtl $outputpath $type $IDconvertscript $locuszoom

type="permutation"
Rscript $code_folder/5_sQTL_result_LocusZoom_plot.R $splicetype $counttype $PSItype $splicing_input $summaryinput $rootsqtl $outputpath $type $IDconvertscript $locuszoom
