#variable1: full path + file name of plink output
#e.g. /u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project/analysis/1_GLMM/sQTL/sQTL_run/Brain-Cortex_copy/selected_sQTL_LD_GWAS/sig_SNPS_gwas.ld

#The purpose of this code is to annotate all the SNPs that are in high LD with the sQTLs we identified (only SNPs annotated in GWAS catalog AND in high LD with sQTLs are outputted)

#example input:
#sys.argv=['nothing','/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project/analysis/1_GLMM/sQTL/sQTL_run/Brain-Cortex_copy/selected_sQTL_LD_GWAS/sig_SNPS_gwas.ld']

## inputs gwas catalog and plink ld output 
import sys,os,re
from collections import defaultdict

def read_in_plink(plink_fn):     #the purpose of this function is to pair our SNP with SNPs in high LD with our SNP and get the r square of this pair
  plink_in = open(plink_fn,'r')
  firstline = True
  ld={}                      #dictionary with key=SNP pair and key=r square value of this pair
  p = defaultdict(list)      #dictionary with value=SNP in our result and key=SNPs that are in high LD with our SNPs
  for line in plink_in:       
    if firstline:          #skip the first line
      firstline = False
      continue
    fields = line.rstrip().split()
    p[fields[5]].append(fields[2])     #fields[2] is the SNP in our result and fields[5] is the SNPs that are in high LD with our SNPs
    ld[fields[5]+'\t'+fields[2]]=fields[6]
  plink_in.close()
  return p,ld

def main():
  gwas_fn = sys.argv[2]
  plink_fn = sys.argv[1]
  plink_dict,ld = read_in_plink(plink_fn)
  fout=open(plink_fn+'.LD.result.txt','w')
  gwas_dict = {}            #dictionary with key=SNPs that are in high LD with our SNP and value=GWAS catelog information of that SNP
  firstline = True
  gwas_in = open(gwas_fn)
  for line in gwas_in:
    if firstline:        #skip the first line
      firstline = False
      continue
    fields = line.rstrip().split('\t')
    if fields[21] in plink_dict:       #if the SNP in GWAS catalog is one of the SNPs that are in high LD with our SNP
    	if fields[21] in gwas_dict:      #if we have already have this SNP in our result (one SNP can appear in different rows in the GWAS catalog BUT the annotation can be different, so we keep them all)
    		gwas_dict[fields[21]] = [m+';'+n for m,n in zip(gwas_dict[fields[21]],fields)]
    	if fields[21] not in gwas_dict:
    		gwas_dict[fields[21]] = fields   #we add the information of that SNP into our result
  gwas_in.close()
  for snp in gwas_dict:      #for each SNP that is (1) in GWAS catalog and (2) in high LD with our SNP
    for s in plink_dict[snp]:     #get the SNPs in our result that pairs with snp from the previous step
      r_2=''
      if s+'\t'+snp in ld:    #if the pair is in this order
         r_2=ld[s+'\t'+snp]   #get the r square
      else:                   #if the pair is in the reverse order
         r_2=ld[snp+'\t'+s]   #get the r square
      #print s, '\t', snp, '\t',gwas_dict[snp][7], '\t',gwas_dict[snp][34],'\t',gwas_dict[snp][27], '\t',gwas_dict[snp][1],'\t',r_2
      fout.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(s, snp,gwas_dict[snp][7],gwas_dict[snp][34],gwas_dict[snp][27],gwas_dict[snp][1],r_2))

if __name__ == '__main__':
  main()
