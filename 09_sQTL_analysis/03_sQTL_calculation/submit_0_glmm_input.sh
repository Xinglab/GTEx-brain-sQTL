#!/bin/sh

###
# tells the cluster to run the job from the current working directory
#$ -cwd

# tells the cluster what info to inform the user through email
#$ -m beas

# tells the cluster the user's email
#$ -M yidazhan@mail

# tells the cluster to use 1024MB of memory, and 1 hours of run time
#$ -l h_data=30G,h_rt=12:00:00

#$ -e ./input_generation_log/
#$ -o ./input_generation_log/
###


#example: 
#qsub -V -N A3SS_original_JC submit_glmm_input.sh A3SS original JC
#qsub -V -N A5SS_original_JC submit_glmm_input.sh A5SS original JC
#qsub -V -N MXE_original_JC submit_glmm_input.sh MXE original JC
#qsub -V -N RI_original_JC submit_glmm_input.sh RI original JC
#qsub -V -N SE_original_JC submit_glmm_input.sh SE original JC

#qsub -V -N A3SS_logit_JC submit_glmm_input.sh A3SS logit JC
#qsub -V -N A5SS_logit_JC submit_glmm_input.sh A5SS logit JC
#qsub -V -N MXE_logit_JC submit_glmm_input.sh MXE logit JC
#qsub -V -N RI_logit_JC submit_glmm_input.sh RI logit JC
#qsub -V -N SE_logit_JC submit_glmm_input.sh SE logit JC

#qsub -V -N A3SS_original_JCEC submit_glmm_input.sh A3SS original JCEC
#qsub -V -N A5SS_original_JCEC submit_glmm_input.sh A5SS original JCEC
#qsub -V -N MXE_original_JCEC submit_glmm_input.sh MXE original JCEC
#qsub -V -N RI_original_JCEC submit_glmm_input.sh RI original JCEC
#qsub -V -N SE_original_JCEC submit_glmm_input.sh SE original JCEC

#qsub -V -N A3SS_logit_JCEC submit_glmm_input.sh A3SS logit JCEC
#qsub -V -N A5SS_logit_JCEC submit_glmm_input.sh A5SS logit JCEC
#qsub -V -N MXE_logit_JCEC submit_glmm_input.sh MXE logit JCEC
#qsub -V -N RI_logit_JCEC submit_glmm_input.sh RI logit JCEC
#qsub -V -N SE_logit_JCEC submit_glmm_input.sh SE logit JCEC


splicetype=$1
PSItype=$2
counttype=$3

# this is most important part, it tells the machine that I want to use R
source /u/local/Modules/default/init/modules.sh
module load R

# run the R code
code_folder="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/scripts"
Rscript $code_folder/0_glmm_input.R $splicetype $PSItype $counttype
