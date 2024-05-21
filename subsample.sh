#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16g
#SBATCH --time=10:00:00
#SBATCH --job-name=subsample_testing
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbysh17@nottingham.ac.uk

#change to home directory
cd /gpfs01/home/mbysh17
# source conda environments
source ~/.bashrc

#activate environment for subsampling
conda activate subsample_env

#obtain list of all samples in the vcf file
bcftools query -l core_files/1001genomes_snp_biallelic_only_ACGTN.vcf > core_files/all_vcf_samples.txt

#run python script to make subsampled VCF and phenotype file (.csv format)
#selecting 500 samples from the ~1000 samples list and choosing Mo98 phenotype

for ((i=200;i<=1000; i+=200));do

echo "current i val is: ${i}"
python3 batch_files/subsample.py -p core_files/master_list.csv -s core_files/all_vcf_samples.txt -n ${i} -t leaf_ionome_Mo98 -op core_files/subsampled_phenotype_${i}.csv -og core_files/subsampled_genotype_${i}.vcf

done
conda deactivate
# you should now have the files needed to run GWAS and GIFT on a random subsample of the main data

echo "Script finished!"
#end of script
