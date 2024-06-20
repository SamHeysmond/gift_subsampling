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

#obtain list of all samples in the vcf file
bcftools query -l core_files/1001genomes_snp_biallelic_only_ACGTN.vcf > core_files/all_vcf_samples.txt

#deactivate conda environment
conda deactivate

# maybe add a way to make output_files/ directory and remove old one?

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
# first remove old JOB_LIST
# but remove the previous one if it exists
rm core_files/JOB_LIST.csv
echo "JOB_ID,SUBSAMPLE_N,PHENOTYPE" > core_files/JOB_LIST.csv

# remove leftover csv and vcf files in core files
rm core_files/subsampled_phenotype_*.csv
rm core_files/subsampled_genotype_*.vcf

# Make the subsample scripts ========================================
# enter python environment
conda activate python3_env

# run the python script which will make all of the scripts needed
python3 batch_files/subsample_script_setup.py

echo "Subsample_setup.sh finished!"

#end of script
