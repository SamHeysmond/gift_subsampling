#!/bin/bash
#SBATCH --job-name=run_and_monitor
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1g
#SBATCH --time=6-00:00:00
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbysh17@nottingham.ac.uk

source ~/.bashrc
conda activate python3_env

# To change directory etc, modify the python file run_and_monitor.py
python3 run_and_monitor.py

conda deactivate
#end of script
