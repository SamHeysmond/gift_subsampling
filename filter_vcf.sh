#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=35g
#SBATCH --time=06:00:00
#SBATCH --job-name=filter_vcf
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbysh17@nottingham.ac.uk

##################################################
#REVISED_CODE ####################################
##################################################

##################################################
## FILTER CHECKLIST ##############################
##################################################
# 1 - CHROMOSOMES 1-5 ONLY (NO CHLORO OR MITO)
# 2 - MAF 
# 3 - MAX MISS 
# 4 - MIN Q 
# 5 - BIALLELIC (MAX ALLELES 2)
# 6 - MIN GQ 
# 7 - MIN MEAN DEPTH 
# 8 - MAX MEAN DEPTH
# 9 - MIN DEPTH 
# 10 - MAX DEPTH 
# 11 - SINGLETONS = NONE ALLOWED
# 12 - LD (DISABLED) 
# 13 - IMPUTED MISSING SNPS = TRUE
##################################################
##################################################

##################################################
#### S T A R T  F I L E S ########################

# 1 - raw.vcf.gz (downloaded)
# 2 - beagle.06Aug24.a91.jar (downloaded)

##################################################
##################################################

source ~/.bashrc

cd /gpfs01/home/mbysh17/core_files

##################################################
#### REMOVING RESIDUAL FILES #####################
##################################################
echo "removing residual files . . ."

# rm raw.vcf.gz.tbi
# rm new_chr_only.vcf
rm quality_filtered.vcf
rm new_no_singletons.vcf
rm out.singletons
rm LD_PRUNED.vcf
rm IMPUTED.vcf.gz
rm IMPUTED.vcf
rm IMPUTED.log
rm TAIR10_chr_all.fas.fai
# rm FINAL.vcf 

##################################################
##################################################

##################################################
#### SETTING FILTER VARIABLES ####################
##################################################

MAX_ALLELES=2
MAF=0.01
MISS=0.98
QUAL=30
G_QUAL=25
MIN_DEPTH=15
MAX_DEPTH=30

##################################################
##################################################

echo "Filtering VCF . . . "

##################################################
#### INDEX THE RAW VCF FILE ######################
##################################################

# conda activate bcf_env
# #zip and tabix the file if necessary
# echo "Zipping and indexing files where needed"
# #bgzip 1001genomes_snp_biallelic_only_ACGTN.vcf
# #tabix 1001genomes_snp_biallelic_only_ACGTN.vcf
# # tabix 1001genomes_snp-short-indel_with_tair10_only_ACGTN.vcf.gz
# tabix raw.vcf.gz

##################################################
##################################################

##################################################
#### FILTERING FOR 1 #############################
##################################################

# View only chromosomes 1-5
# echo "selecting chromosomes to view"
# bcftools view -r 1,2,3,4,5 -Oz -o new_chr_only.vcf.gz raw.vcf.gz  
# # doesnt make a gz file so need to make one myself
# echo "indexing file . . ."
# tabix -p vcf new_chr_only.vcf.gz
# conda deactivate 
##################################################
##################################################

##################################################
#### FILTERING FOR 2-10 ##########################
##################################################

# # activate vcftool environment
conda activate vcftools_env

# Filter for depth and quality of SNPs
echo "Filtering with vcftools . . . "
# vcftools --gzvcf new_chr_only.vcf.gz \
# --remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL \
# --max-alleles $MAX_ALLELES --minGQ $G_QUAL \
# --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
# --minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout > quality_filtered.vcf

# min DP and max DP caused a problem so just using min Mean max mean
vcftools --gzvcf new_chr_only.vcf.gz \
--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL \
--max-alleles $MAX_ALLELES --minGQ $G_QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--recode --stdout > quality_filtered.vcf

conda deactivate
##################################################
##################################################

##################################################
### FILTERING FOR 11 #############################
##################################################
conda activate vcftools_env

echo "Finding singletons . . . "
# finds singletons and puts them in a file called (out.singletons)
#vcftools --singletons --vcf depth_quality_no_indels_recoded.recode.vcf
vcftools --singletons --vcf quality_filtered.vcf
# this will generate out.singletons

echo "Removing singletons . . . "
#vcftools --vcf depth_quality_no_indels_recoded.recode.vcf --exclude-positions out.singletons --recode --recode-INFO-all --out no_singletons
#vcftools --vcf depth_quality_no_indels_recoded.recode.vcf --exclude-positions out.singletons --recode --recode-INFO-all --out no_singletons --stdout > no_singletons.vcf
vcftools --vcf quality_filtered.vcf --exclude-positions out.singletons --recode --recode-INFO-all --out new_no_singletons --stdout > new_no_singletons.vcf

# exit vcftools environment
conda deactivate

##################################################
##################################################

##################################################
### FILTERING FOR 12 #############################
##################################################
# enter bcftools environment
# conda activate bcf_env

# #filtering for LD
# echo "Filtering for LD"
# #bcftools +prune -m 0.25 -w 1000 no_singletons.vcf -Ov -o LD_pruned.vcf
# bcftools +prune -m 0.4 -w 10kb new_no_singletons.vcf -Ov -o LD_PRUNED.vcf

# conda deactivate

##################################################
##################################################

##################################################
### FILTERING FOR 13 #############################
##################################################

echo "Imputing SNPs . . . "

conda activate java_env

# java -jar beagle.06Aug24.a91.jar gt=LD_PRUNED.vcf out=IMPUTED
java -Xmx39g -jar beagle.06Aug24.a91.jar gt=new_no_singletons.vcf out=IMPUTED

# outputs IMPUTED.vcf.gz

conda deactivate

# re-headder the vcf to regain lost info from beagle #####
conda activate bcf_env

#unzip since it compresses it
bgzip -d IMPUTED.vcf.gz

# create fasta index file
echo "creating fasta index file"
samtools faidx TAIR10_chr_all.fas

# replace the chr of the chromosome in the .fas.fai file to nothing
sed -i -e 's/Chr//g' TAIR10_chr_all.fas.fai

# reheadder the vcf
echo "Re-headering the VCF since beagle loses information"
bcftools reheader --fai TAIR10_chr_all.fas.fai -o FINAL.vcf IMPUTED.vcf

#this LD pruned and imputed vcf can then be passed into the subsample script to..
#.. pick out the subsample population we want.

##################################################
##################################################

echo "Script finished!"

#end of script
