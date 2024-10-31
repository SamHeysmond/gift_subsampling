#!/bin/bash
#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3g
#SBATCH --time=00:10:00
#SBATCH --job-name=subsample_setup
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbysh17@nottingham.ac.uk

##################################################
#REVISED_CODE ####################################
##################################################

#change to home directory
cd /gpfs01/home/mbysh17

# source conda environments
source ~/.bashrc

#activate environment for subsampling
conda activate bcf_env

#obtain list of all sample IDs in the vcf file
bcftools query -l core_files/FINAL.vcf > core_files/all_vcf_samples.txt

#deactivate conda environment
conda deactivate

# make directory for the subsample text files to go in
# but remove the previous one if it exists
rm -rf core_files/subsample_text_files
mkdir core_files/subsample_text_files

# make directory for the subsample data (genotype and phenotype)
# but remove the previous one if it exists
rm -rf core_files/subsampled_data
mkdir core_files/subsampled_data

# clear and remake the genotype tracker folder
# this stores the genotype of all SNPs for each sample in a given GIFT run
rm -rf core_files/genotype_tracker
mkdir core_files/genotype_tracker

# make directory for parallel tasks to go in
# but remove the previous one if it exists
rm -rf batch_files/parallel
mkdir batch_files/parallel

# clear and remake the running/completed folder for jobs
rm -rf batch_files/processing_parallel
mkdir batch_files/processing_parallel

rm -rf batch_files/completed_parallel
mkdir batch_files/completed_parallel

# make job list file in core files
# first remove old JOB_LIST if it exists
rm core_files/JOB_LIST.csv
# write in the header for the csv
echo "JOB_ID,SUBSAMPLE_N,PHENOTYPE" > core_files/JOB_LIST.csv

# Make the subsample scripts ========================================

# enter python environment
conda activate python3_env

# run the python script which will make all of the individual scripts needed
python3 batch_files/subsample_script_setup.py -subsampleFile core_files/subsample_numbers_list.txt

echo "Subsample_setup.sh finished!"
#end of script
