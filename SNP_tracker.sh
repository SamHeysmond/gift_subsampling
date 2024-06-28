#!/bin/bash
#SBATCH --partition=hmemq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=1350g
#SBATCH --time=7-00:00:00
#SBATCH --job-name=SNP_tracker_sh
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbysh17@nottingham.ac.uk

#change to home directory
cd /gpfs01/home/mbysh17

# source conda environment
source ~/.bashrc

# clear current directory 
rm -rf output_files/R_DATA/

# Set up directory
mkdir output_files/R_DATA/

# enter python 3 environment
conda activate gift_env

# Make the R scripts needed for second stage of analysis
# this script also makes the batch file to run each R script it makes
#python3 batch_files/SNP_tracker.py -d /gpfs01/home/mbysh17/output_files/
python3 batch_files/SNP_tracker.py

# exit conda environment
conda deactivate

echo "end of SNP_tracker.sh"
#end of script