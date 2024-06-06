#!/bin/bash 
#SBATCH --job-name=leaf_ionome_Sr88_global 
#SBATCH --partition=shortq 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --mem=2g 
#SBATCH --time=01:30:00 
#SBATCH --output=/gpfs01/home/sbzsmb/OandE/%x.out 
#SBATCH --error=/gpfs01/home/sbzsmb/OandE/%x.err 
#SBATCH --parsable 
echo RUNNING ON `hostname` 

source $HOME/.bash_profile 
conda activate p2_env 


pygwas run -t none -a amm -g /gpfs01/home/sbzsmb/Data/SNP_Matrix/full_imputed_SNP_MATRIX -k /gpfs01/home/sbzsmb/Data/SNP_Matrix/full_imputed_SNP_MATRIX/kinship_ibs_binary_mac5.h5py -o /gpfs01/home/sbzsmb/Running_GWAS/GWAS_Results/Global/leaf_ionome_Sr88_GWAS.csv /gpfs01/home/sbzsmb/Running_GWAS/Phenotypes/leaf_ionome_Sr88.csv