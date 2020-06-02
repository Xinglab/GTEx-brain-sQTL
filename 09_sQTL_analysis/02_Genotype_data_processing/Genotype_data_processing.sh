#Purpose of this code:
#generate .ped and .map file from vcd file 

module load vcftools
output_prefix="output_prefix"
vcftools --vcf GTEx_genotype_data.vcf --maf 0.05 --out $output_prefix --plink
