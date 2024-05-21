# This file tracks:
#           1) Hand selected SNPs depending on the phenotype
#           2) Top 20 significant SNPs in both GWAS and GIFT for each test




# loop through the .csv from head to tail.
# if a SNP pval is lower than in the current list (which will start with all 1's) then replace it with next lowest SNP (1 in this case)
# repeat until lowest p-values are in the list
# grep out the lines that have those values in them and write them to a SNP table "top_20_SNP_<SAMPLE_NUM>_ID_ID_NUM_HERE>"
# grep out the lines containing hand picked SNP and put those in a SNP table "Focused_SNP_<SAMPLE_NUM>_ID_ID_NUM_HERE>"
# Pass these data onto R script which will plot them? (Need to alter the R script in the physics GWAS GIFT program to do this too....)