# SCRIPT INFO
# This script will make the dozens of subruns for each...
# PHENOTYPE, SUBSAMPLE NUMBER, and COPY of that subsample num (100 copies each)
# so for 2 phenotypes that's -> 2x 5 x 100 files = 1000 batch files...whew


# open the list of phenotypes to analyse
# contains format (no header)
# leaf_ionome_Mo98
# leaf_ionome_Na23
# .... 
phenotypes_input=open("core_files/phenotypes_list.txt","r")

# alter this to change the number of subsamples per run
# keep in mind the if statements below will also need modifying too if so
subsample_list=[200,400,600,800,1000]

#loop through each phenotype with the following settings
for phenotype in phenotypes_input:
    phenotype = phenotype.replace('\n','')
    print ("Current phenotype is: ", phenotype)
    # loop for each different amount of samples
    for subsample_num in subsample_list:

        # make 100 copies of each file 
        # (alter this to change the number of tests to be done)
        for copynum in range(1,100):
            bash_script_output=open("batch_files/parallel/subrun_"+
                                    str(phenotype)+
                                    "_"+
                                    str(subsample_num)+
                                    "_"+
                                    str(copynum)+
                                    ".sh"
                                    ,"w"
                                    )
            # SLURM variables
            bash_script_output.write(f'#!/bin/bash\n')
            bash_script_output.write(f'#SBATCH --partition=defq\n')
            bash_script_output.write(f'#SBATCH --nodes=1\n')
            bash_script_output.write(f'#SBATCH --ntasks=1\n')
            bash_script_output.write(f'#SBATCH --cpus-per-task=4\n')
            # give appropriate memory and time according to the job's subsample number... 
            # ...(bigger sample -> more resources)
            if subsample_num==int(200):
                bash_script_output.write(f'#SBATCH --mem=6g\n')
                bash_script_output.write(f'#SBATCH --time=06:00:00\n')
            elif subsample_num==int(400):
                bash_script_output.write(f'#SBATCH --mem=8g\n')
                bash_script_output.write(f'#SBATCH --time=12:00:00\n')
            elif subsample_num==int(600):
                bash_script_output.write(f'#SBATCH --mem=8g\n')
                bash_script_output.write(f'#SBATCH --time=18:00:00\n')
            elif subsample_num==int(800):
                bash_script_output.write(f'#SBATCH --mem=8g\n')
                bash_script_output.write(f'#SBATCH --time=24:00:00\n')
            elif subsample_num==int(1000):
                bash_script_output.write(f'#SBATCH --mem=10g\n')
                bash_script_output.write(f'#SBATCH --time=30:00:00\n')
            bash_script_output.write(f'#SBATCH --job-name=subrun\n')
            bash_script_output.write(f'#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out\n')
            bash_script_output.write(f'#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err\n')
            bash_script_output.write(f'#SBATCH --mail-type=ALL\n')
            bash_script_output.write(f'#SBATCH --mail-user=mbysh17@nottingham.ac.uk\n')
            bash_script_output.write(f'#===============================\n')
            bash_script_output.write(f'#change to home directory\n')
            bash_script_output.write(f'cd /gpfs01/home/mbysh17\n')
            bash_script_output.write(f'# source conda environments\n')
            bash_script_output.write(f'source ~/.bashrc\n')
            bash_script_output.write(f'# in case any environment is active, deactivate conda env\n')
            bash_script_output.write(f'conda deactivate\n')
            bash_script_output.write(f'#activate environment for subsampling\n')
            bash_script_output.write(f'conda activate subsample_env\n')

            # variables
            bash_script_output.write(f'#set variables for this script\n')
            bash_script_output.write(f'i={subsample_num}\n')
            bash_script_output.write(f'phenotype={phenotype}\n')

            #job ID
            bash_script_output.write(f'#append the JOB ID and subsample number (i) to a file\n')
            bash_script_output.write('echo "${SLURM_JOB_ID},${i},${phenotype}" >> core_files/JOB_LIST.csv \n')
            bash_script_output.write(f'# run the subsample script to get the appropriate number of subsamples\n')
            bash_script_output.write('python3 batch_files/subsample.py -p core_files/master_list.csv -s core_files/all_vcf_samples.txt -n ${i} -t ${phenotype} -op core_files/subsampled_phenotype_${i}_${SLURM_JOB_ID}.csv -og core_files/subsampled_genotype_${i}_${SLURM_JOB_ID}.vcf -ri ${SLURM_JOB_ID}\n')
            bash_script_output.write(f'#deactivate conda environment\n')
            bash_script_output.write(f'conda deactivate\n')
            bash_script_output.write(f'# you should now have the files needed to run GWAS and GIFT on a random subsample of the main data\n')
            bash_script_output.write('echo "Subsampling finished for job: ${SLURM_JOB_ID}, moving on to gift_vs_gwas stage."\n')
            bash_script_output.write(f'# the following code will use the subsampled vcf and phenotype data \n')
            bash_script_output.write(f'# to compare gift and gwas by running them both\n')
            bash_script_output.write(f'#activate conda environment for the GIFT method\n')
            bash_script_output.write(f'conda activate gift_env\n')
            bash_script_output.write(f'#run GIFT on the subsampled data\n')
            bash_script_output.write('python3 -W ignore core_files/physics_GWAS_OOP.py -v core_files/subsampled_genotype_${i}_${SLURM_JOB_ID}.vcf -f core_files/subsampled_phenotype_${i}_${SLURM_JOB_ID}.csv -p ${phenotype} -o output_files/${phenotype}_whole_genome_metrics_${i}_${SLURM_JOB_ID}.csv -id ${SLURM_JOB_ID}\n')
            bash_script_output.write(f'# exit gift environment\n')
            bash_script_output.write(f'conda deactivate\n')
            bash_script_output.write('echo "GIFT for: ${SLURM_JOB_ID} finished"\n')
            bash_script_output.write(f'# remove the subsampled vcf file to save on storage\n')
            bash_script_output.write('rm core_files/subsampled_genotype_${i}_${SLURM_JOB_ID}.vcf\n')
            bash_script_output.write(f'#activate gwas environment\n')
            bash_script_output.write(f'conda activate gwas_env\n')
            bash_script_output.write(f'#run gwas on subsampled data\n')
            bash_script_output.write('pygwas run -t none -a amm -g core_files/SNP_Matrix/full_imputed_SNP_MATRIX -k core_files/SNP_Matrix/full_imputed_SNP_MATRIX/kinship_ibs_binary_mac5.h5py -o output_files/${phenotype}_GWAS_${i}_${SLURM_JOB_ID}.csv core_files/subsampled_phenotype_${i}_${SLURM_JOB_ID}.csv\n')
            bash_script_output.write(f'# exit gwas environment\n')
            bash_script_output.write(f'conda deactivate\n')
            bash_script_output.write(f'# remove the subsampled phenotype file to save on storage \n')
            bash_script_output.write('rm core_files/subsampled_phenotype_${i}_${SLURM_JOB_ID}.csv\n')
            bash_script_output.write('echo "GWAS for: ${SLURM_JOB_ID} finished"\n')

            # calling for the python script to first MAKE the R script for GWAS manhattan plots
            # to make manhattan plots from the GWAS data
            bash_script_output.write(f'conda activate python3_env\n')
            bash_script_output.write('python3 batch_files/make_r_scripts.py -id ${SLURM_JOB_ID} -i ${i} -p ${phenotype} -o output_files/\n')
            bash_script_output.write(f'conda deactivate\n')
            bash_script_output.write(f'conda activate r_env\n')
            bash_script_output.write('Rscript output_files/${SLURM_JOB_ID}_${i}_${phenotype}.R\n')
            bash_script_output.write(f'conda deactivate\n')
            bash_script_output.write(f'#end of script')
            bash_script_output.close()
print("====================\n")
print("subsample_script_setup.py has finished!")
print("====================\n")
# end of script