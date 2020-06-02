#!/bin/sh

#$ -cwd 
#$ -m beas
#$ -M yidazhan@mail
#$ -l h_data=2G,h_rt=1:00:00
#$ -e ./part1_log/
#$ -o ./part1_log/

code_folder="/path/to/sQTL/scripts"

splicetype="SE"        ###change###
counttype="JC"
PSItype="logit"
#PSItype="logit" means the PSI value is after logit transformation + confounder conrrection
#PSItype="original" means the PSI value is before logit transformation and no confounder correction. In this case, we will do logit transformation later in the pipeline
#the idea here is to run the sQTL analysis on data before + after confounder correction so eventually they are both on logit scale

splicing_input="/path/to/input_splicing/"$PSItype/$counttype/$splicetype
rootoutput="/path/to/sQTL_run/"$PSItype/$counttype/$splicetype 


cd $code_folder

brainregionlist="brainregionlist.txt"

while read myline || [ -n "$myline" ]                                    
do                                                                       
	jobname="p1_"$PSItype"_"$counttype"_"$splicetype"_"$myline
	#echo $jobname
	/u/systems/UGE8.0.1/bin/lx-amd64/qsub -N $jobname $code_folder/pipeline_part1.sh $myline $splicetype $code_folder $splicing_input $rootoutput $counttype $PSItype
done < $brainregionlist               #the input file needs to be read like this. By doing this, what happened in the while loop can be used outside the loop                                            
    
    


