#!/bin/bash
#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=20g
#SBATCH --time=04:00:00
#SBATCH --job-name=filter_vcf
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbysh17@nottingham.ac.uk

#change to core directory
cd /gpfs01/home/mbysh17/core_files
# source conda environments
source ~/.bashrc


conda activate bcftools_env

#bcftools +fill-tags 1001genomes_snp_biallelic_only_ACGTN.vcf  -Ov --output output_1.vcf -- -t AN,AC



# filter by allele numer
#bcftools filter -i 'AN < 30' output_1.vcf -Ov -o output_2.vcf

# filter by allele count (maybe try AC<15 if needed?)
#bcftools filter -i 'AC < 30' output_1.vcf -Ov -o output_test.vcf

# works but tiny
#bcftools view -q 0.05:minor 1001genomes_snp_biallelic_only_ACGTN.vcf | bcftools filter -i 'AC>=30 & FMT/DP>=10 & FMT/DP<=50 & F_MISSING<0.01' -Ov -o output_test.vcf

# changed to info/DP -> tiny vcf with nothing in it
#bcftools view -q 0.05:minor 1001genomes_snp_biallelic_only_ACGTN.vcf | bcftools filter -i 'AC>=30 & INFO/DP>=10 & INFO/DP<=50 & F_MISSING<0.01' -Ov -o output_test2.vcf

# removed upper boundary on depth  AND changed back to FMT
bcftools view -q 0.05:minor 1001genomes_snp_biallelic_only_ACGTN.vcf | bcftools filter -i 'AC>=30 & FMT/DP>=10 & F_MISSING<0.01' -Ov -o output_test3.vcf


conda deactivate

echo "bcftools job done"

echo "filter script finished"

# end of filter script