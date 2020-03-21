#!/bin/sh

###
# tells the cluster to run the job from the current working directory
#$ -cwd

# tells the cluster what info to inform the user through email
#$ -m a

# tells the cluster the user's email
#$ -M yidazhan@mail

# tells the cluster to use 1024MB of memory, and 1 hours of run time
#$ -l h_data=8G,h_rt=12:00:00

# tells the cluster to run on 100 machines     
#$ -t 1-819:1                                #44 unique sQTL exon times 13 brain regions times number of gwas studies

#$ -e ./log/
#$ -o ./log/
###


# this is most important part, it tells the machine that I want to use R
source /u/local/Modules/default/init/modules.sh
module load R

# run the R code
code_folder="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/2_Colocalization_input"

Rscript $code_folder/2_Colocalization_input_job_array.R
