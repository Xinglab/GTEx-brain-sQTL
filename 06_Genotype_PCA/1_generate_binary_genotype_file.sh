#reference: http://www.bioinf.wits.ac.za/courses/gwas/cupop.pdf

#get binary data from VCF
module load plink
cd /u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project/data/raw_data/files/genotype/V7_whole_exon_sequencing
outputpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/4.6_data_inspection_PCA_on_genotype/result"
plink --vcf GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC_maf_0.01.vcf --maf 0.05 --make-bed --out $outputpath/Genotype_V7_plink_binary_maf0.05

#pruning
#see https://www.cog-genomics.org/plink/1.9/ld and 
#2.6 Population stratification analysis section in the supplementary of Genetic effects on gene expression across human tissues
module load plink
outputpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/4.6_data_inspection_PCA_on_genotype/result"
cd $outputpath
plink --allow-no-sex --bfile Genotype_V7_plink_binary_maf0.05 --indep-pairwise 200 100 0.2 --out $outputpath/Genotype_V7_plink_binary_pruning_maf0.05

#extract data from pruned result
module load plink
outputpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/4.6_data_inspection_PCA_on_genotype/result"
cd $outputpath
plink --allow-no-sex --bfile Genotype_V7_plink_binary_maf0.05 --extract $outputpath/Genotype_V7_plink_binary_pruning_maf0.05.prune.in \
--make-bed --out Genotype_V7_plink_binary_pruned_maf0.05
