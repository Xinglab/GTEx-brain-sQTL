#!/bin/sh

###
# tells the cluster to run the job from the current working directory
#$ -cwd

# tells the cluster what info to inform the user through email
#$ -m a

# tells the cluster the user's email
#$ -M yidazhan@mail

# tells the cluster to use 1024MB of memory, and 1 hours of run time
#$ -l h_data=4G,h_rt=2:00:00

# tells the cluster to run on 100 machines     
#$ -t 1-819:1                                #44 unique sQTL exon times 13 brain regions times number of gwas studies

#$ -e ./log/
#$ -o ./log/
###


# this is most important part, it tells the machine that I want to use R
source /u/local/Modules/default/init/modules.sh
module load R/3.2.3

# run the R code
code_folder="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/3_Colocalization_analysis"

Rscript $code_folder/3.1_Colocalization_analysis_job_array.R
