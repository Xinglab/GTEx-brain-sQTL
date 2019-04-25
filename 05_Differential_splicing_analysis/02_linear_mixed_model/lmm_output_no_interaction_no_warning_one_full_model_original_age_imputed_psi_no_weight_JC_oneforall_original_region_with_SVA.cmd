#!/bin/sh

###
# tells the cluster to run the job from the current working directory
#$ -cwd

# tells the cluster what info to inform the user through email
#$ -m a

# tells the cluster the user's email
#$ -M yidazhan@mail

# tells the cluster to use 1024MB of memory, and 1 hours of run time
#$ -l h_data=8G,h_rt=23:00:00

#$ -e ./log/
#$ -o ./log/
###

splicetype=$1
code_folder=$2

# this is most important part, it tells the machine that I want to use R
source /u/local/Modules/default/init/modules.sh
module load R

# run the R code
Rscript $code_folder/lmm_output_no_interaction_no_warning_one_full_model_original_age_imputed_psi_no_weight_JC_oneforall_original_region_with_SVA.R $splicetype $code_folder