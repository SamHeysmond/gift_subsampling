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

# standard cd command across all batch scripts
cd /gpfs01/home/mbysh17

# source conda from user home
source ~/.bashrc

# # set up channels
# conda config --add channels r
# conda config --add channels bioconda
# conda config --add channels conda-forge
# conda config --add channels fvogt257

# # #conda create -n gift_env -y
# conda activate gift_env
# # conda install python=3.10 -y
# # conda install r-base -y
# # conda install bioconductor-biocinstaller -y
# # conda install r-tidyverse -y
# # conda install r-htmlwidgets -y
# # conda install r-manhattanly -y
# # conda install r-ggrepel -y
# conda install r-ggplot2 -y
# conda install r-dplyr -y
# conda install r-ggstatsplot -y # need to cite this (type into R to find the citation again)
# # pip3 install pandas
# conda deactivate

# conda create -n vcf2gwas_env -y
# conda activate vcf2gwas_env
# conda install vcf2gwas -y
# conda deactivate


# # modin install tricky- try one of the below
# # try the ray install one first but other lines may be necessary
# pip install "modin[ray]"
# # pip install modin --pre
# # pip install "modin[all]" 
# # pip install -U modin
# # pip3 install "modin[all]"

# conda deactivate

# conda create -n python3_env -y
# conda activate python3_env
# conda install python=3.10 -y
# pip3 install pandas
# conda deactivate

# conda create -n gwas_env -y
# conda activate gwas_env
# conda install python=2.7 -y
# pip install PyGWAS
# conda deactivate

# conda create -n r_env -y
# conda activate r_env
# conda install r-base -y
# conda install bioconductor-biocinstaller -y
# conda install r-qqman
# conda install r-tidyverse -y
# conda install r-dplyr -y
# conda install r-ggrepel -y
# conda deactivate

# conda create -n gwas_env2 -y
# conda activate gwas_env2
# conda install python=2.7
# pip install PyGWAS
# conda deactivate

# conda create -n gatk_env -y
# conda activate gatk_env
# conda install -c bioconda gatk4 -y
# conda install python -y
# conda deactivate

# conda create -n bcf_env bcftools tabix -y
# conda deactivate

# conda create -n subsample_env -y
# conda activate subsample_env
# conda install bcftools -y
# conda install python -y
# pip3 install pandas
# conda deactivate


# create samtools environment
# conda create -n samtools_env 
# conda activate samtools_env
# conda install samtools -y
# conda deactivate

# create bcftools environment
# conda create -n bcftools_env
# conda activate bcftools_env
# conda install bcftools -y
# conda deactivate

# create bedtools environment
conda create -n bedtools_env -y
conda activate bedtools_env
conda install bedtools -y
# conda deactivate


echo "Script finished!"
#end of script
