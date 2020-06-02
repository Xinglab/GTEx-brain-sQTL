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

splicetype=$1
PSItype=$2
counttype=$3

# this is most important part, it tells the machine that I want to use R
source /u/local/Modules/default/init/modules.sh
module load R

# run the R code
code_folder="/path/to/0_glmm_input.R"
Rscript $code_folder/0_glmm_input.R $splicetype $PSItype $counttype
