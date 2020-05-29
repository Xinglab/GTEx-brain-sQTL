#get binary data from VCF
module load plink

genotypepath="/path/to/GTEx/genotype/file"
genotypename="genotype_file_name.vcf"
outputprefix="output_prefix"
outputpath="/output_path"
plink --vcf $genotypepath/$genotypename --maf 0.05 --make-bed --out $outputpath/$outputprefix

#pruning
cd $outputpath
plink --allow-no-sex --bfile $outputprefix --indep-pairwise 200 100 0.2 --out $outputpath/$outputprefix"_pruning"

#extract data from pruned result
module load plink
cd $outputpath
plink --allow-no-sex --bfile $outputprefix --extract $outputpath/$outputprefix"_pruning"".prune.in" \
--make-bed --out $outputprefix"_pruned"
