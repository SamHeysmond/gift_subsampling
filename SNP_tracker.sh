#!/bin/bash
#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=90g
#SBATCH --time=01:00:00
#SBATCH --job-name=SNP_tracker_sh
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbysh17@nottingham.ac.uk

# changed over to short q
# changed mem down from 1350GB to 90GB
# changed CPUS down from 32 to 16
# changed time down from 7 days to 1 hour

#change to home directory
cd /gpfs01/home/mbysh17

# source conda environment
source ~/.bashrc

# clear summary data directory if it already exists
rm -rf output_files/R_DATA/

# Set up directory for all summary data to go into 
mkdir output_files/R_DATA/

# enter python 3 environment
conda activate gift_env

# concatonate files needed for second stage of analysis
python3 batch_files/SNP_tracker.py

# exit conda environment
conda deactivate

echo "end of SNP_tracker.sh"
#end of script