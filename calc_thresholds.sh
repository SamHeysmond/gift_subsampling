#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=15g
#SBATCH --time=01:00:00
#SBATCH --job-name=calc_thres
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbysh17@nottingham.ac.uk

#change to core directory
cd /gpfs01/home/mbysh17/
# source conda environments
source ~/.bashrc

# activate environment for python use
conda activate gift_env

# run thresholds calculator script
# python batch_files/calc_thresholds.py
python batch_files/calc_thresholds.py -subsampleFile core_files/subsample_numbers_list.txt


# exit conda environment
conda deactivate 

echo "thres bash script done"
# end of file