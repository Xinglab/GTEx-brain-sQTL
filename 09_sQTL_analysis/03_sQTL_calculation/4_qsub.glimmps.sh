#!/bin/bash 
#$ -S /bin/bash                   
#$ -R y                                   
#$ -l h_data=14G,h_rt=23:00:00
#$ -V       
#$ -cwd 
#$ -j y       
#$ -m a



export s=`sed -n ${SGE_TASK_ID}p job.list`     #this code just link the job number (1-n:1) with job ID and each row in job.list (the first row in job.list will have job ID=1, e.g., glimmps.o449015.1)
echo $s
#e.g. chr1/chr1.glimmps.0.sh
bash $s         #run the shell script
