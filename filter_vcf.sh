#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=16g
#SBATCH --time=24:00:00
#SBATCH --job-name=filter_vcf
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbysh17@nottingham.ac.uk

##################################################
#REVISED_CODE ####################################
##################################################

source ~/.bashrc

cd /gpfs01/home/mbysh17/core_files

echo "Filtering VCF"

conda activate bcf_env

#zip and tabix the file if necessary
echo "Zipping and indexing files where needed"
#bgzip 1001genomes_snp_biallelic_only_ACGTN.vcf
#tabix 1001genomes_snp_biallelic_only_ACGTN.vcf
# tabix 1001genomes_snp-short-indel_with_tair10_only_ACGTN.vcf.gz
tabix raw.vcf.gz

# View only chromosomes 1-5
echo "selecting chromosomes to view"

#bcftools view -r 1,2,3,4,5 1001genomes_snp_biallelic_only_ACGTN.vcf.vcf.gz > chr_only.vcf
#bcftools view -r 1,2,3,4,5 1001genomes_snp-short-indel_with_tair10_only_ACGTN.vcf.gz > new_chr_only.vcf
bcftools view -r 1,2,3,4,5 raw.vcf.gz > new_chr_only.vcf
# doesnt make a gz file

# Filter for allele count of 30 and only for biallelic locations
# AC of 30 liable to change but will need to alter the python GIFT code too
echo "Filtering for allele count and biallelic sites"
#bcftools view --min-ac=30 --max-alleles 2 chr_only.vcf > min_ac_30_biallelic.vcf
bcftools view --min-ac=30 --max-alleles 2 new_chr_only.vcf > new_min_ac_30_biallelic.vcf

# Add in tags for vcftools to use (this may not be needed?)
#echo "Adding in tags"
#bcftools +fill-tags min_ac_30_biallelic.vcf -Ov --output tagged.vcf --t AN,AC

# exit the bcftool environment
conda deactivate

# # activate vcftool environment
conda activate vcftools_env
# # Filter for depth and quality of SNPs
# echo "Filtering for depth and quality"
# #vcftools --vcf tagged.vcf --minDP 3 --minGQ 25 --remove-indels --recode --recode-INFO-all --out depth_quality_no_indels_recoded
# vcftools --vcf new_min_ac_30_biallelic.vcf --minDP 3 --minGQ 25 --remove-indels --recode --recode-INFO-all --out new_depth_quality_no_indels_recoded --stdout > new_depth_quality_no_indels_recoded.vcf

#Filter for missingness in the vcf (5% threshold for Gemma)
# do this AFTER selecting subsample population -> move to python code ?
# max missing of 10 or 9 might be too harsh?. It results in removing any SNPs that have  10 or 9 missing GTs
# added filter for maf of 0.05 but liable to change to 0.1 if too harsh
vcftools --vcf new_depth_quality_no_indels_recoded.vcf --max-missing-count 50 --maf 0.05 --recode --recode-INFO-all --out new_remove_missing --stdout > new_remove_missing.vcf

echo "Finding singletons"
# finds singletons and puts them in a file called (out.singletons)
#vcftools --singletons --vcf depth_quality_no_indels_recoded.recode.vcf
vcftools --singletons --vcf new_remove_missing.vcf
# this will generate out.singletons

echo "Removing singletons"
#vcftools --vcf depth_quality_no_indels_recoded.recode.vcf --exclude-positions out.singletons --recode --recode-INFO-all --out no_singletons
#vcftools --vcf depth_quality_no_indels_recoded.recode.vcf --exclude-positions out.singletons --recode --recode-INFO-all --out no_singletons --stdout > no_singletons.vcf
vcftools --vcf new_remove_missing.vcf --exclude-positions out.singletons --recode --recode-INFO-all --out new_no_singletons --stdout > new_no_singletons.vcf

# exit vcftools environment
conda deactivate

# enter bcftools environment
conda activate bcf_env

#filtering for LD
echo "Filtering for LD"
#bcftools +prune -m 0.25 -w 1000 no_singletons.vcf -Ov -o LD_pruned.vcf
bcftools +prune -m 0.4 -w 10kb new_no_singletons.vcf -Ov -o LD_PRUNED.vcf

conda deactivate

#imput missing SNP data with beagle
echo "Imputing SNPs"

conda activate java_env

java -jar beagle.06Aug24.a91.jar gt=LD_PRUNED.vcf out=IMPUTED

conda deactivate


# re-headder the vcf to regain lost info from beagle
conda activate bcf_env

#unzip since it compresses it
# remove imputed.vcf in case its already there
rm IMPUTED.vcf
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

echo "Script finished!"

#end of script
