#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8g
#SBATCH --time=01:00:00
#SBATCH --job-name=env_setup
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbysh17@nottingham.ac.uk
cd /gpfs01/home/mbysh17

source ~/.bashrc
#conda config --add channels r
#conda config --add channels bioconda
#conda config --add channels conda-forge

#conda create -n gift_env -y
conda activate gift_env
#conda install python=3.12
#conda install r-base -y
#conda install bioconductor-biocinstaller -y
#conda install r-tidyverse -y
#conda install r-htmlwidgets -y
#conda install r-manhattanly -y
#conda install r-ggrepel -y
conda install -c conda-forge modin-all
#pip3 install pandas
#pip3 install "modin[all]"
conda deactivate

#conda create -n python3_env -y
#conda activate python3_env
#conda install python=3.12 -y
#pip3 install pandas
#conda deactivate

#conda create -n gwas_env -y
#conda activate gwas_env
#conda install python=2.7 -y
#pip install PyGWAS
#conda deactivate

#conda create -n r_env -y
#conda activate r_env
#conda install r-base -y
#conda install bioconductor-biocinstaller -y
#conda install r-qqman
#conda install r-tidyverse -y
#conda install r-dplyr -y
#conda install r-ggrepel -y
#conda deactivate

#conda create -n gwas_env2 -y
#conda activate gwas_env2
#conda install python=2.7
#pip install PyGWAS
#echo "testing command for pygwas"
#pygwas -h
#echo "end of testing command for pygwas"
#conda deactivate

#conda create -n gatk_env -y
#conda activate gatk_env
#conda install -c bioconda gatk4 -y
#conda install python -y
#conda deactivate

#conda create -n bcf_env bcftools tabix -y
#conda deactivate

#conda create -n subsample_env -y
#conda activate subsample_env
#conda install bcftools -y
#conda install python -y
#pip3 install pandas
#conda deactivate

echo "Script finished!"
#end of script
