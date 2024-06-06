#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8g
#SBATCH --time=08:00:00
#SBATCH --job-name=subrun
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbysh17@nottingham.ac.uk

#change to home directory
cd /gpfs01/home/mbysh17
# source conda environments
source ~/.bashrc

# in case any environment is active
conda deactivate

#activate environment for subsampling
conda activate subsample_env

# run python script to make subsampled VCF and phenotype file (.csv format)
# selecting 500 samples from the ~1000 samples list and choosing Mo98 phenotype

#set i value to number of samples desired
i=200

phenotype=leaf_ionome_Mo98

echo "current i val is: ${i}"

#append the JOB ID and subsample number (i) to a file
# !!! this file will need REMOVING at end of pipeline (maybe after the make_manhattan plots R script) !!!!!!
echo "${SLURM_JOB_ID},${i},${phenotype}" >> core_files/JOB_LIST.csv 

# run the subsample script to get the appropriate number of subsamples
python3 batch_files/subsample.py -p core_files/master_list.csv -s core_files/all_vcf_samples.txt -n ${i} -t ${phenotype} -op core_files/subsampled_phenotype_${i}_${SLURM_JOB_ID}.csv -og core_files/subsampled_genotype_${i}_${SLURM_JOB_ID}.vcf -ri ${SLURM_JOB_ID}

#deactivate conda environment
conda deactivate

# you should now have the files needed to run GWAS and GIFT on a random subsample of the main data
echo "==================================="
echo "Subsampling finished for job: ${SLURM_JOB_ID}, moving on to gift_vs_gwas stage."
echo "==================================="

# the following code will use the subsampled vcf and phenotype data 
# to compare gift and gwas by running them both

#activate conda environment for the GIFT method
conda activate gift_env

#run GIFT on the subsampled data
python3 core_files/physics_GWAS_OOP.py -v core_files/subsampled_genotype_${i}_${SLURM_JOB_ID}.vcf -f core_files/subsampled_phenotype_${i}_${SLURM_JOB_ID}.csv -p ${phenotype} -o output_files/${phenotype}_whole_genome_metrics_${i}_${SLURM_JOB_ID}.csv -id ${SLURM_JOB_ID}

# exit gift environment
conda deactivate

echo "GIFT for: ${SLURM_JOB_ID} finished"

# remove the subsampled vcf file to save on storage (TEMP PAUSED)
# rm core_files/subsampled_genotype_${i}_${SLURM_JOB_ID}.vcf

#activate gwas environment
conda activate gwas_env

#run gwas on subsampled data
pygwas run -t none -a amm -g core_files/SNP_Matrix/full_imputed_SNP_MATRIX -k core_files/SNP_Matrix/full_imputed_SNP_MATRIX/kinship_ibs_binary_mac5.h5py -o output_files/${phenotype}_GWAS_${i}_${SLURM_JOB_ID}.csv core_files/subsampled_phenotype_${i}_${SLURM_JOB_ID}.csv

# exit gwas environment
conda deactivate

# remove the subsampled phenotype file to save on storage (TEMP PAUSED)
# rm core_files/subsampled_phenotype_${i}_${SLURM_JOB_ID}.csv

echo "GWAS for: ${SLURM_JOB_ID} finished"
#end of script
