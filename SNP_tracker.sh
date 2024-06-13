#!/bin/bash
#SBATCH --partition=hmemq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=1200g
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


# clear current directory output if exists
rm -rf batch_files/parallel_stage2/
rm -rf batch_files/completed_parallel_stage2/
rm -rf output_files/SNP_tracker_R_scripts/
rm -rf output_files/summary_plots/IDEA1/
rm -rf output_files/summary_plots/IDEA2/
rm -rf output_files/summary_plots/IDEA3/
rm -rf output_files/R_DATA/

# Set up directory to put the R scripts into
mkdir batch_files/parallel_stage2/
mkdir batch_files/completed_parallel_stage2/
mkdir output_files/SNP_tracker_R_scripts/
mkdir output_files/summary_plots/
mkdir output_files/summary_plots/IDEA1/
mkdir output_files/summary_plots/IDEA2/
mkdir output_files/summary_plots/IDEA3/
mkdir output_files/R_DATA/

# enter python 3 environment
conda activate gift_env

# Make the R scripts needed for second stage of analysis
# this script also makes the batch file to run each R script it makes
python3 batch_files/SNP_tracker.py -d /gpfs01/home/mbysh17/output_files/

# run the modified "SNP" run and modifier to execute the R scripts in parallel
python3 -u batch_files/Rscript_run_and_monitor.py

# exit conda environment
conda deactivate

echo "end of SNP_tracker.sh"
#end of script