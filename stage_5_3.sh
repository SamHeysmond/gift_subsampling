#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16g
#SBATCH --time=02:00:00
#SBATCH --job-name=Stage_5_3
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbysh17@nottingham.ac.uk
#===============================
# This script runs all stage 5 R scripts to make...
# ... the plots for the paths
#===============================
echo "start of stage 5 script"

# source conda environments
source ~/.bashrc

#change to home directory
cd /gpfs01/home/mbysh17/

######################################################
#### S T A G E : 5

rm -r  core_files/combined_position_data
rm -r  batch_files/stage_5_3_prerun
rm -r  batch_files/stage_5_3_prorun
rm -r  output_files/stage_5_results
rm -r  batch_files/stage_5_3_Rscripts/

mkdir core_files/combined_position_data
mkdir batch_files/stage_5_3_prerun
mkdir batch_files/stage_5_3_prorun
mkdir output_files/stage_5_results
mkdir batch_files/stage_5_3_Rscripts/


conda activate gift_env

python batch_files/stage_5_3.py

python batch_files/stage_5_3_ram.py

conda deactivate

# end of script