#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=10g
#SBATCH --time=24:00:00
#SBATCH --job-name=GWAS_protocol
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbysh17@nottingham.ac.uk

##################################################
#REVISED_CODE ####################################
##################################################

source ~/.bashrc

cd /gpfs01/home/mbysh17/ProjectDS_Revised

conda activate bcf_env

echo "Compressing vcf"

bgzip no_singletons.recode.vcf

echo "Indexing VCF"

tabix no_singletons.recode.vcf.gz

echo "Running GWAS pipeline"

conda deactivate

#conda activate ??

# optional covariate file for third input
#bash run_gwas_gemma.sh phenotype.tsv vcf_file.gz > log.txt


echo "Script finished!"

#end of script
