#!/bin/bash
#SBATCH --partition=hmemq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=300g
#SBATCH --time=24:00:00
#SBATCH --job-name=filter_GWAS
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

#clear directory of R parallel and completed
rm -r /gpfs01/home/mbysh17/batch_files/R_parallel/
rm -r /gpfs01/home/mbysh17/batch_files/completed_R_parallel/
rm -r /gpfs01/home/mbysh17/output_files/csv_before_filter/

# make the directories for the subruns
mkdir /gpfs01/home/mbysh17/batch_files/R_parallel/
mkdir /gpfs01/home/mbysh17/batch_files/completed_R_parallel/
mkdir /gpfs01/home/mbysh17/output_files/csv_before_filter/

# move all csvs before filter to a temp folder
cp /gpfs01/home/mbysh17/output_files/*.csv /gpfs01/home/mbysh17/output_files/csv_before_filter/

echo "entering python script"

conda activate modin_env

# run script to fitler data and create R scripts
python batch_files/filter_snps_GWAS.py

echo "filter GWAS script finished, running the R scripts"

# run the script that will run each R script in R_parallel
python batch_files/GWAS_run_and_monitor.py

conda deactivate

# end of filter script