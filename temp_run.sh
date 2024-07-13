#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4g
#SBATCH --time=72:00:00
#SBATCH --job-name=batch_template
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbysh17@nottingham.ac.uk

#change to core directory
cd /gpfs01/home/mbysh17/core_files
# source conda environments
source ~/.bashrc

#change to home directory
cd /gpfs01/home/mbysh17/

conda activate python3_env

python batch_files/GWAS_run_and_monitor.py

conda deactivate

echo "Script finished!"
#end of script
