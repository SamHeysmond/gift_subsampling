#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10g
#SBATCH --time=1:00:00
#SBATCH --job-name=GWAS_manhattan_plot
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbysh17@nottingham.ac.uk

source ~/.bashrc

conda activate r_env

Rscript make_manhattan_gwas.R

conda deactivate

echo "Script finished!"
#end of script
