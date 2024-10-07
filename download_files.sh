#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16g
#SBATCH --time=24:00:00
#SBATCH --job-name=download_files
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbysh17@nottingham.ac.uk

##################################################
#REVISED_CODE ####################################
##################################################


source ~/.bashrc
cd /gpfs01/home/mbysh17/core_files/

conda activate axel_env

# echo "Downloading vcf"
# #wget https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_with_tair10_only_ACGTN.vcf.gz 
# #axel -n 16 https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_with_tair10_only_ACGTN.vcf.gz --output=1001genomes_snp-short-indel_with_tair10_only_ACGTN.vcf.gz
# axel -n 16 https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz --output=raw.vcf.gz

# # use option -q for quiet if it takes up too much space in error file


# echo "downloading beagle"
# axel -n 16 https://faculty.washington.edu/browning/beagle/beagle.06Aug24.a91.jar --output=beagle.06Aug24.a91.jar

echo "downloading fasta reference file"
# axel -n 16 https://www.arabidopsis.org/download/file?path=Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz --output=TAIR10_chr_all.fas.gz

axel -n 16 https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz --output=TAIR10_chr_all.fas.gz


echo "decompressing fasta reference file"
bgzip -d TAIR10_chr_all.fas.gz


echo "Downloading phenotype file"
# download phenotype file (need to do this live)
axel -n 16 https://ffionexplorer.nottingham.ac.uk/ionmap/session/73c35af63087870deee190afb8c8b089/download/downloadData1?w= --output=selected_Ion_data.csv


axel -n 16 https://github.com/genetics-statistics/GEMMA/archive/refs/heads/master.zip

unzip gwas_gemma-master.zip

conda deactivate
echo "Script finished!"
#end of script
