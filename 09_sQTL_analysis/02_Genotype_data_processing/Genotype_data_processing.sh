#Purpose of this code:
#generate .ped and .map file from vcd file 

module load vcftools
cd /u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project/data/raw_data/files/genotype/V7_whole_exon_sequencing
vcftools --vcf GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC_maf_0.01.vcf --maf 0.05 --out Genotype_V7_vcftools_maf0.05 --plink
