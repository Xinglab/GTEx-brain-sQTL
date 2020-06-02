module load plink
inputpath="/input/path/to/genotype/result/from/1_generate_binary_genotype_file.sh"
outputpath="/output/path"
outputprefix="output_prefix"    #from 1_generate_binary_genotype_file.sh

cd $outputpath
prunedfile=$outputprefix"_pruned"

#on unpruned data:
plink --bfile $inputpath/$outputprefix --pca 635 --out $outputpath/$outputprefix
#on pruned data:
plink --bfile $inputpath/prunedfile --pca 635 --out $outputpath/prunedfile
