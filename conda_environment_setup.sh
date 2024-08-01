#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20g
#SBATCH --time=24:00:00
#SBATCH --job-name=env_setup
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbysh17@nottingham.ac.uk

# This script sets up all necessary conda environments for the analysis

# standard cd command across all batch scripts
cd /gpfs01/home/mbysh17

# source conda from user home
source ~/.bashrc

# set up channels
conda config --add channels r
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels fvogt257

# Environment for running gift (keep)
conda create -n gift_env -y
conda activate gift_env
conda install python=3.10 -y
conda install r-base -y
conda install bioconductor-biocinstaller -y # (check)
conda install r-tidyverse -y
conda install r-htmlwidgets -y # (check)
conda install r-manhattanly -y # (check)
conda install r-ggrepel -y
conda install r-ggplot2 -y
conda install r-dplyr -y
conda install r-ggstatsplot -y 
pip3 install pandas
conda deactivate

# modin install tricky (due to dependencies)- try one of the below in order
# try the ray install one first but other lines may be necessary if ray fails
pip install "modin[ray]"
# pip install modin --pre
# pip install "modin[all]" 
# pip install -U modin
conda deactivate

# Python 3 environment 
conda create -n python3_env -y
conda activate python3_env
conda install python=3.10 -y
pip3 install pandas
conda deactivate

# Gwas environment for running GWAS software
conda create -n gwas_env -y
conda activate gwas_env
conda install python=2.7 -y
pip install PyGWAS
conda deactivate

# R environment
conda create -n r_env -y
conda activate r_env
conda install r-base -y
conda install bioconductor-biocinstaller -y
conda install r-qqman
conda install r-tidyverse -y
conda install r-dplyr -y
conda install r-ggrepel -y
conda deactivate

# Gatk Environment
conda create -n gatk_env -y
conda activate gatk_env
conda install -c bioconda gatk4 -y
conda install python -y
conda deactivate

# for post-project testing vcf prune environment
conda create -n bcf_env bcftools tabix -y
conda deactivate

# environment for subsampling
conda create -n subsample_env -y
conda activate subsample_env
conda install bcftools -y
conda install python -y
pip3 install pandas
conda deactivate

# samtools environment
conda create -n samtools_env 
conda activate samtools_env
conda install samtools -y
conda deactivate

# Bcftools environment
conda create -n bcftools_env
conda activate bcftools_env
conda install bcftools -y
conda deactivate

# Bedtools environment
conda create -n bedtools_env -y
conda activate bedtools_env
conda install bedtools -y
conda deactivate

echo "Script finished!"
#end of script
