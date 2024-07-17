#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=4g
#SBATCH --time=01:00:00
#SBATCH --job-name=subsample_setup
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbysh17@nottingham.ac.uk

#change to home directory
cd /gpfs01/home/mbysh17

# source conda environments
source ~/.bashrc

#activate environment for subsampling
conda activate subsample_env

#obtain list of all sample IDs in the vcf file
bcftools query -l core_files/1001genomes_snp_biallelic_only_ACGTN.vcf > core_files/all_vcf_samples.txt

#deactivate conda environment
conda deactivate

# make directory for the subsample text files to go in
# but remove the previous one if it exists
rm -rf core_files/subsample_text_files
mkdir core_files/subsample_text_files

# make directory for parallel tasks to go in
# but remove the previous one if it exists
rm -rf batch_files/parallel
mkdir batch_files/parallel

# clear and remake the running/completed folder for jobs
rm -rf batch_files/completed_parallel
mkdir batch_files/completed_parallel

# make job list file in core files
# first remove old JOB_LIST if it exists
rm core_files/JOB_LIST.csv
# write in the header for the csv
echo "JOB_ID,SUBSAMPLE_N,PHENOTYPE" > core_files/JOB_LIST.csv

# remove leftover csv and vcf files in core files if any exist from previous runs
rm core_files/subsampled_phenotype_*.csv
rm core_files/subsampled_genotype_*.vcf

# Make the subsample scripts ========================================
# enter python environment
conda activate python3_env

# run the python script which will make all of the scripts needed
python3 batch_files/subsample_script_setup.py

echo "Subsample_setup.sh finished!"

#end of script
