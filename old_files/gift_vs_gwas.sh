#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g
#SBATCH --time=72:00:00
#SBATCH --job-name=gift_vs_gwas
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbysh17@nottingham.ac.uk

# this script will use the subsampled vcf and phenotype data generated
# to compare gift and gwas by running them both

#change directory to my user home directory
cd /gpfs01/home/mbysh17

#source conda environments
source ~/.bashrc

#start for loop (subsampling at stages of 200 up to 1000 samples at most
for ((i=200;i<=1000;i+=200));do

#activate conda environment for the GIFT method
conda activate gift_env

#run GIFT on the subsampled data
python3 core_files/physics_GWAS_OOP.py -v core_files/subsampled_genotype_${i}.vcf -f core_files/subsampled_phenotype_${i}.csv -p leaf_ionome_Mo98 -o output_files/Mo_whole_genome_metrics_${i}.csv

# exit gift environment
conda deactivate

#activate gwas environment
conda activate gwas_env

#run gwas on subsampled data
pygwas run -t none -a amm -g core_files/SNP_Matrix/full_imputed_SNP_MATRIX -k core_files/SNP_Matrix/full_imputed_SNP_MATRIX/kinship_ibs_binary_mac5.h5py -o output_files/leaf_ionome_Mo98_GWAS_${i}.csv core_files/subsampled_phenotype_${i}.csv

# exit gwas environment
conda deactivate

#exit for loop
done

echo "Script finished!"
#end of script
