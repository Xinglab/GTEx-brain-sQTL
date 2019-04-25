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
#$ -l h_data=8G,h_rt=8:00:00

#$ -e ./11_2_sQTL_sashimi_plot_log/
#$ -o ./11_2_sQTL_sashimi_plot_log/
###


###pay attention: index all the BAM files first before running sashimi plot!!! This only needs to be done once for each BAM file###
###Before running the code, annotate out miniconda in ~/.bashrc to make sure that we only have python2 (rmats2sashimiplot cannot run using python3)


#output the python version to make sure we are using the right version (python2 instead of python3)
which python

splicetype="SE"        ###change###
echo $splicetype

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
totalcountinput="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/1_get_matrix_from_rMATS/result/"$counttype
GWASdbpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project/analysis/1_GLMM/sQTL/GWAS_databases"
GWASdbname="gwas_catalog_v1.0.1-associations_e89_r2017-07-31.tsv"
genoplinkprefix="Genotype_V7_plink_maf0.05"
LDpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project/analysis/1_GLMM/sQTL/LD"


# this is most important part, it tells the machine that I want to use R
source /u/local/Modules/default/init/modules.sh
module load R/3.2.3

# run the R code
code_folder="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/scripts"

type="pvalue"
Rscript $code_folder/4_sQTL_result_sashimi_plot.R $splicetype $counttype $PSItype $splicing_input $summaryinput $rootsqtl $outputpath $type $totalcountinput $GWASdbpath $GWASdbname $genoplinkprefix $LDpath

type="permutation"
Rscript $code_folder/4_sQTL_result_sashimi_plot.R $splicetype $counttype $PSItype $splicing_input $summaryinput $rootsqtl $outputpath $type $totalcountinput $GWASdbpath $GWASdbname $genoplinkprefix $LDpath
