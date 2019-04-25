#tutorial see http://www.bioinf.wits.ac.za/courses/gwas/cupop.pdf + https://www.cog-genomics.org/plink/1.9/strat#pca
#this works. The --pca option is the same as in EIGENSTRAT
#see https://www.biostars.org/p/242370/
#the two evidence is from https://www.cog-genomics.org/plink/1.9/strat#pca and http://cnsgenomics.com/software/gcta/#PCA
module load plink
inputpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/4.6_data_inspection_PCA_on_genotype/result"
outputpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/4.6_data_inspection_PCA_on_genotype/result/PCA"
cd $outputpath
#on unpruned data:
plink --bfile $inputpath/Genotype_V7_plink_binary_maf0.05 --pca 635 --out $outputpath/Genotype_V7_plink_binary_maf0.05
#on pruned data:
plink --bfile $inputpath/Genotype_V7_plink_binary_pruned_maf0.05 --pca 635 --out $outputpath/Genotype_V7_plink_binary_pruned_maf0.05
#we want to output all PCs and since we have 635 donors, we output 635 PCs which are all the PCs