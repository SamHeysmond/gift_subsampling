#!/bin/bash
#SBATCH --partition=hmemq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=250g
#SBATCH --time=10:00:00
#SBATCH --job-name=filter_snp
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbysh17@nottingham.ac.uk

#change to core directory
cd /gpfs01/home/mbysh17/core_files

# source conda environments
source ~/.bashrc

# activate gatk environment
conda activate gatk_env

# create FASTA sequence dictionary file
gatk CreateSequenceDictionary -R TAIR10_chr_all.fas

conda deactivate 

echo "fasta dictionary created"

# activate samtools env
conda activate samtools_env

#create fasta index file
samtools faidx TAIR10_chr_all.fas

conda deactivate

echo "fasta index created"

conda activate bcftools_env

# add in the tags for AN and AC into the vcf -> creating a new vcf
bcftools +fill-tags 1001genomes_snp_biallelic_only_ACGTN.vcf  -Ov --output output_1.vcf -- -t AN,AC

# filter by allele count to match filtering in GIFT python script
# this will give all SNPs that contain less than 30 AC
bcftools filter -i 'AC < 30' output_1.vcf -Ov -o output_2.vcf

conda deactivate

echo "bcftools job done"

conda activate gatk_env

# Convert the VCF to a table with the columns "CHROM" and "POS"
gatk VariantsToTable -V output_2.vcf -F CHROM -F POS -O output_3.table

conda deactivate

echo "done with filtering for SNPs"

# filter some csvs with this data in python (no longer needed)

# conda activate gift_env

# #change to home directory
# cd /gpfs01/home/mbysh17/

# # filter out the SNPs with less than 30 AC from the compiled csv files
# python batch_files/filter_snps.py

# conda deactivate

# echo "filter script finished"

# end of filter script