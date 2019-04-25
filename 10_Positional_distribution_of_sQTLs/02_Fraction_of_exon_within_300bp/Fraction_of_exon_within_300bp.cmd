#!/bin/sh

###
# tells the cluster to run the job from the current working directory
#$ -cwd

# tells the cluster what info to inform the user through email
#$ -m beas

# tells the cluster the user's email
#$ -M yidazhan@mail

# tells the cluster to use 1024MB of memory, and 1 hours of run time
#$ -l h_data=8G,h_rt=23:00:00

#$ -e ./log/
#$ -o ./log/
###


br=$1
code_folder=$2

echo $br
echo $code_folder

# this is most important part, it tells the machine that I want to use R
source /u/local/Modules/default/init/modules.sh
module load R

# run the R code
Rscript $code_folder/Fraction_of_exon_within_300bp.R $br
