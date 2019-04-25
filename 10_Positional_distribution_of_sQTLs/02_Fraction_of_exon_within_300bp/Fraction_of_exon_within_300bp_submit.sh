#!/bin/sh

#$ -cwd 
#$ -m beas
#$ -M yidazhan@mail
#$ -l h_data=2G,h_rt=1:00:00
#$ -e ./log/
#$ -o ./log/

code_folder="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.2_sQTL_SNP_annotation/Figure4A_fraction_of_exon_within_300bp_modified"

cd $code_folder

brainregionlist="brainregionlist.txt"

while read myline || [ -n "$myline" ]                                    
do                                                                       
	jobname=$myline
	#echo $jobname
	/u/systems/UGE8.0.1/bin/lx-amd64/qsub -N $jobname $code_folder/Fraction_of_exon_within_300bp.cmd $myline $code_folder
done < $brainregionlist               #the input file needs to be read like this. By doing this, what happened in the while loop can be used outside the loop                                            
    
    


