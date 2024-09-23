#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16g
#SBATCH --time=24:00:00
#SBATCH --job-name=download_vcf
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out
#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err
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

echo "downloading beagle"
axel -n 16 https://faculty.washington.edu/browning/beagle/beagle.06Aug24.a91.jar --output=beagle.06Aug24.a91.jar

echo "downloading fasta reference file"
axel -n 16 https://www.arabidopsis.org/download/file?path=Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz --output=TAIR10_chr_all.fas.gz

echo "decompressing fasta reference file"
bgzip -d TAIR10_chr_all.fas.gz


#echo "Downloading phenotype file"
# download phenotype file
#axel -n 16 https://ffionexplorer.nottingham.ac.uk/ionmap/session/d57e687b98c81cbc0d54670a13741db1/download/downloadData1?w=

conda deactivate
echo "Script finished!"
#end of script
