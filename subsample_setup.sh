#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4g
#SBATCH --time=00:10:00
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

# make directory for the subsample text files to go in
mkdir core_files/subsample_text_files/

# make directory for parallel tasks to go in
mkdir batch_files/parallel

#multiply the subsample sh script 100 times for 100 scripts to be made
# CHANGED THIS TO 10 WHILE TESTING TO SPEED THINGS UP!
for ((n=1;n<=10;n++));do

echo "current n val is: ${n}"

# copy the template for subsample and run into the parallel folder
cp batch_files/subsample_and_run.sh batch_files/parallel/subsample_and_run_${n}.sh

#exit for loop
done

# make job list file in core files
echo "JOB_ID,SUBSAMPLE_N,PHENOTYPE" > core_files/JOB_LIST.csv

echo "Subsample_setup.sh finished!"
#end of script
