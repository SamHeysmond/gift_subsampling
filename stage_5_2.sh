#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2g
#SBATCH --time=10:00:00
#SBATCH --job-name=Stage_5_2
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbysh17@nottingham.ac.uk
#===============================
# This script runs all stage 4 and 5 scripts 
#===============================
echo "start of stage 5 script"

# source conda environments
source ~/.bashrc

#change to home directory
cd /gpfs01/home/mbysh17/

######################################################
#### S T A G E : 5

rm -r  core_files/theta_paths

mkdir core_files/theta_paths


conda activate gift_env

python batch_files/stage_5_run_and_monitor.py

conda deactivate

# end of script