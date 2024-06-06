#!/bin/bash
#SBATCH --job-name=physics_GWAS_1
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8g
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs01/home/sbzsmb/OandE/%x.out
#SBATCH --error=/gpfs01/home/sbzsmb/OandE/%x.err

source $HOME/.bash_profile
conda activate ngs_pipe_env_gatk4

date
python3 /gpfs01/home/sbzsmb/Scripts/physics_GWAS_OOP.py -v /gpfs01/home/sbzsmb/Data/1001genomes_snp_biallelic_chrom2_only_ACGTN.vcf -f /gpfs01/home/sbzsmb/Data/master_list.csv -p leaf_ionome_Mo98 -o /gpfs01/home/sbzsmb/Data/Physics_GWAS/Mo_whole_genome_metrics.csv
date
