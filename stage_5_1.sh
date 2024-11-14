#!/bin/bash
#SBATCH --partition=hmemq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=51
#SBATCH --mem=650g
#SBATCH --time=10:00:00
#SBATCH --job-name=Stage_5_1
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbysh17@nottingham.ac.uk
#===============================
# This script runs all stage 5 scripts to filter ...
# ... genotypes to only the positions given
# It also creates R scripts to calculating ...
# ... theta paths.
#===============================
# mem down from 800
# cpu down from 51
# hmemq to defq

echo "start of stage 5 script"

# source conda environments
source ~/.bashrc

#change to home directory
cd /gpfs01/home/mbysh17/

######################################################
#### S T A G E : 5

# rm -r  core_files/filtered_genotype_tracker
# rm -r  core_files/theta_paths
rm -r batch_files/stage_5_prerun
rm -r batch_files/stage_5_prorun
rm -r output_files/significant_SNP_data

# mkdir core_files/filtered_genotype_tracker
# mkdir core_files/theta_paths
mkdir batch_files/stage_5_prerun
mkdir batch_files/stage_5_prorun
mkdir output_files/significant_SNP_data

conda activate gift_env

python batch_files/stage_5_filter_significant_at_max_subsample.py

conda deactivate

# end of script