#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8g
#SBATCH --time=1:00:00
#SBATCH --job-name=batch_template
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbysh17@nottingham.ac.uk

cd /gpfs01/home/mbysh17
cp -r /gpfs01/home/sbzsmb/Students_2024/SNP_Matrix core_files/
cp -r /gpfs01/home/sbzsmb/Students_2024/1001genomes_snp_biallelic_only_ACGTN.vcf core_files/
cp -r /gpfs01/home/sbzsmb/Students_2024/master_list.csv core_files/
cp -r /gpfs01/home/sbzsmb/Students_2024/physics_GWAS_OOP.py core_files/
cp -r /gpfs01/home/sbzsmb/Students_2024/physics_gwas_1.sh core_files/
cp -r /gpfs01/home/sbzsmb/Students_2024/leaf_ionome_Sr88_GWAS.sh core_files/

echo "Script finished!"
#end of script
