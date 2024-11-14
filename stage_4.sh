#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=6g
#SBATCH --time=02:00:00
#SBATCH --job-name=Stage_4
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbysh17@nottingham.ac.uk
#===============================
# This script runs all stage 4 and 5 scripts 
#===============================
echo "start of stage 4 script"

# source conda environments
source ~/.bashrc

#change to home directory
cd /gpfs01/home/mbysh17/

######################################################
#### S T A G E : 4

# clear and set up directories required for this
rm -r batch_files/stage_4_prerun
rm -r batch_files/stage_4_prorun
rm -r batch_files/stage_4_R_scripts
rm -r output_files/STAGE_4_RESULTS

mkdir batch_files/stage_4_prerun
mkdir batch_files/stage_4_prorun
mkdir batch_files/stage_4_R_scripts
mkdir output_files/STAGE_4_RESULTS

conda activate python3_env

# makes the R scripts and batch files that accompany and run them in tandem
python batch_files/stage_4.py

# runs the batch files which run the R scripts
python batch_files/stage_4_run_and_monitor.py

conda deactivate

# end of script