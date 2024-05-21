#!/bin/bash
#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=3g
#SBATCH --time=00:10:00
#SBATCH --job-name=clear_files
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbysh17@nottingham.ac.uk

cd /gpfs01/home/mbysh17/

# remove the subsample files which denote the exact samples used per run
rm -rf core_files/subsample_text_files/*

# remove all results from the output folder to save space
# but keep the directory 
rm -rf output_files/*

# clears the slurmoandE file
rm -rf slurmOandE/*

echo "Script finished!"
#end of script
