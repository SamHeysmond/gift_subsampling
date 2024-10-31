#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=25g
#SBATCH --time=10:00:00
#SBATCH --job-name=filter_vcf
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbysh17@nottingham.ac.uk

##################################################
#### S T A R T  F I L E S ########################
# chromosomes 1-5
# new_chr_only.vcf

##################################################

source ~/.bashrc

cd /gpfs01/home/mbysh17/core_files

mkdir subset_info
##################################################
#### REMOVING RESIDUAL FILES #####################
##################################################
echo "removing residual files"

##################################################
##################################################

##################################################
#### GATHERING VCF STATISTICS ####################
##################################################

# OUT=~/vcftools/cichlid_subset
SUBSET_VCF=~/core_files/new_chr_only.vcf.gz
OUT=~/core_files/subset_info

# compress and set aside vcf to read from
echo "Compressing file"
conda activate bcf_env
bgzip new_chr_only.vcf
bcftools index new_chr_only.vcf.gz
conda deactivate

conda activate vcftools_env

# calc allele freq
echo "Calculating allele freq"
vcftools --gzvcf $SUBSET_VCF --freq2 --out $OUT --max-alleles 2

# calc mean depth per individual
echo "Calculating mean depth per individual"
vcftools --gzvcf $SUBSET_VCF --depth --out $OUT

# calc mean depth per site
echo "Calculating mean depth per site"
vcftools --gzvcf $SUBSET_VCF --site-mean-depth --out $OUT

# calc site quality
echo "Calculating site quality"
vcftools --gzvcf $SUBSET_VCF --site-quality --out $OUT

# calc proportion of missing data per individual
echo "Calculating proportion of missing data per individual"
vcftools --gzvcf $SUBSET_VCF --missing-indv --out $OUT

# calc proportion of missing data per site
echo "Calculating proportion of missing data per site"
vcftools --gzvcf $SUBSET_VCF --missing-site --out $OUT

# calc heterozygosity and inbreeding coefficient per individual
echo "Calculating heterozygosity and inbreeding coefficient per individual"
vcftools --gzvcf $SUBSET_VCF --het --out $OUT


##################################################
##################################################

echo "Script finished!"

#end of script
