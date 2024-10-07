# SCRIPT INFO
# This script will make the dozens of subruns for each...
# METHOD, PHENOTYPE, SUBSAMPLE NUMBER, and COPY of that subsample num (100 copies each)
# so for 2 phenotypes that's -> 2x 2x 5 x 100 files = 2000 batch files...whew

# open the list of phenotypes to analyse
# contains the following format (no header):
# leaf_ionome_Mo98
#Mo98
# leaf_ionome_Na23
#Na23
# .... 
phenotypes_input=open("core_files/phenotypes_list.txt","r")

# alter this to change the number of subsamples per run
# keep in mind the if statements below will also need modifying too if so
subsample_list=[200,400,600,800,999]

#loop through each phenotype with the following settings
for phenotype in phenotypes_input:
    phenotype = phenotype.replace('\n','')
    switch = True
    # loop for each different amount of samples
    for subsample_num in subsample_list:

        # make 100 copies of each file 
        # (alter this to change the number of tests to be done)
        for copynum in range(1,101): #increased to 101 to get 100 copies
            if subsample_num==int(999) and switch==False:
                pass
            else:
                # puts all the runs in a folder called "batch_files/parallel/...."
                output_f=open("batch_files/parallel/subrun_"+
                                        str(phenotype)+
                                        "_"+
                                        str(subsample_num)+
                                        "_"+
                                        str(copynum)+
                                        ".sh"
                                        ,"w"
                                        )
                # SLURM variables
                output_f.write(f'#!/bin/bash\n')
                output_f.write(f'#SBATCH --partition=shortq\n')
                output_f.write(f'#SBATCH --nodes=1\n')
                output_f.write(f'#SBATCH --ntasks=1\n')
                output_f.write(f'#SBATCH --cpus-per-task=2\n')
                # give appropriate memory and time according to the job's subsample number... 
                # ...(bigger sample -> more resources)
                if subsample_num==int(200):
                    output_f.write(f'#SBATCH --mem=4g\n')
                    output_f.write(f'#SBATCH --time=01:00:00\n')
                elif subsample_num==int(400):
                    output_f.write(f'#SBATCH --mem=5g\n')
                    output_f.write(f'#SBATCH --time=01:00:00\n')
                elif subsample_num==int(600):
                    output_f.write(f'#SBATCH --mem=6g\n')
                    output_f.write(f'#SBATCH --time=01:00:00\n')
                elif subsample_num==int(800):
                    output_f.write(f'#SBATCH --mem=6g\n')
                    output_f.write(f'#SBATCH --time=01:00:00\n')
                elif subsample_num==int(999):
                    # only need to do 1 script for 999 since all the same samples (max sample num)
                    switch = False
                    output_f.write(f'#SBATCH --mem=6g\n')
                    output_f.write(f'#SBATCH --time=02:00:00\n')
                    
                output_f.write(f'#SBATCH --job-name=subrun\n')
                output_f.write(f'#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out\n')
                output_f.write(f'#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err\n')
                output_f.write(f'#SBATCH --mail-type=ALL\n')
                output_f.write(f'#SBATCH --mail-user=mbysh17@nottingham.ac.uk\n')
                output_f.write(f'#===============================\n')
                output_f.write(f'#change to home directory\n')
                output_f.write(f'cd /gpfs01/home/mbysh17\n')
                output_f.write(f'# source conda environments\n')
                output_f.write(f'source ~/.bashrc\n')
                output_f.write(f'\n')
                output_f.write(f'# in case any environment is active, deactivate conda env\n')
                output_f.write(f'conda deactivate\n')
                output_f.write(f'#activate environment for subsampling\n')
                output_f.write(f'conda activate subsample_env\n')
                output_f.write(f'\n')

                # variables for each subrun
                output_f.write(f'#set variables for this script\n')
                output_f.write(f'i={subsample_num}\n')
                output_f.write(f'phenotype={phenotype}\n') 
                output_f.write(f'\n')
                # should contain "processing_parallel" in the file path
                output_f.write("original_location=$(scontrol show job $SLURM_JOB_ID | awk -F= '/Command=/{print $2}')\n")
                output_f.write(f'\n')

                #job ID updating to the csv
                output_f.write(f'#append the JOB ID and subsample number (i) to a file\n')
                output_f.write('echo "${SLURM_JOB_ID},${i},${phenotype}" >> core_files/JOB_LIST.csv \n')
                output_f.write(f'\n')

                # main code to run the subsampling
                output_f.write(f'# run the subsample script to get the appropriate number of subsamples\n')
                output_f.write('python3 batch_files/subsample.py -p core_files/leaf_phenotype.csv -s core_files/all_vcf_samples.txt -n ${i} -t ${phenotype} -op core_files/subsampled_data/subsampled_phenotype_${i}_${SLURM_JOB_ID}.csv -og core_files/subsampled_data/subsampled_genotype_${i}_${SLURM_JOB_ID}.vcf -ri ${SLURM_JOB_ID}\n')
                output_f.write(f'\n')

                output_f.write(f'#deactivate conda environment\n')
                output_f.write(f'conda deactivate\n')
                output_f.write(f'\n')
                # you should now have the files needed to run GWAS and GIFT on a random subsample population

                #running gift and gwas for a given subsample population
                output_f.write('echo "Subsampling finished for job: ${SLURM_JOB_ID}, moving on to gift_vs_gwas stage." \n')
                output_f.write(f'\n')
                # the following code will use the subsampled vcf and phenotype data to compare gift and gwas by running them both

                #################################################
                #### PERFORM GIFT#######

                #activate conda environment for the GIFT method
                output_f.write(f'conda activate gift_env\n')
                output_f.write(f'\n')
                output_f.write(f'echo "Running GIFT_CORE.py" \n')
                output_f.write(f'\n')

                # Use ignore command to prevent the warning from stacking up in the slurm output folder.
                output_f.write('python3 -W ignore batch_files/GIFT_CORE.py -v core_files/subsampled_data/subsampled_genotype_${i}_${SLURM_JOB_ID}.vcf -f core_files/subsampled_data/subsampled_phenotype_${i}_${SLURM_JOB_ID}.csv -p ${phenotype} -o output_files/${phenotype}_whole_genome_metrics_${i}_${SLURM_JOB_ID}.csv -id ${SLURM_JOB_ID} -s ${i}\n')
                output_f.write(f'\n')

                # exit gift environment for python script
                output_f.write(f'conda deactivate\n')

                # enter the R environment to continue the GIFT calculations
                output_f.write(f'\n')
                output_f.write(f'conda activate r_env\n')
                output_f.write(f'\n')
                output_f.write(f'echo "Calculating PSNP8 with gift_testing_giota.R" \n')
                output_f.write(f'\n')
                output_f.write('Rscript batch_files/gift_testing_giota.R core_files/subsampled_data/subsampled_phenotype_${i}_${SLURM_JOB_ID}.csv core_files/genotype_tracker/${SLURM_JOB_ID}_genotypes.csv output_files/leaf_ionome_${phenotype}_whole_genome_metrics_${i}_${SLURM_JOB_ID}.csv\n')
                output_f.write(f'\n')


                output_f.write('echo "GIFT for: ${SLURM_JOB_ID} finished" \n')
                output_f.write(f'conda deactivate\n')
                output_f.write(f'\n')

                ########################################################
                ######################################################################

                # remove the subsampled vcf file to save on storage\n')
                #bash_script_output.write('rm core_files/subsampled_data/subsampled_genotype_${i}_${SLURM_JOB_ID}.vcf\n')
                #bash_script_output.write(f'\n')

                #################################################
                #### PERFORM GWAS #######
                # obtain the list of phenotype and accession IDs separately for the GWAS pipeline
                output_f.write(f'\n')
                output_f.write(f'echo "Preparing files for GWAS" \n')
                output_f.write(f'\n')
                output_f.write(f"cd core_files/subsampled_data/\n")
                output_f.write("tail -n+2 subsampled_phenotype_${i}_${SLURM_JOB_ID}.csv | awk -F ',' '{print $1}' > ${SLURM_JOB_ID}_accession_list.txt\n")
                output_f.write("tail -n+2 subsampled_phenotype_${i}_${SLURM_JOB_ID}.csv | awk -F ',' '{print $2}' > ${SLURM_JOB_ID}_phenotype_list.tsv\n")
                output_f.write(f"cd ~\n")
                output_f.write(f'\n')

                    # activate GWAS environment
                output_f.write(f"conda activate gwas_pipeline \n")
                output_f.write(f'\n')
                    # works but misses some SNPs i think due to relationship matrix? thinks they are pop_strat related?
                #output_f.write("bash gwas_gemma-master/run_gwas_gemma.sh core_files/subsampled_data/${SLURM_JOB_ID}_phenotype_list.tsv core_files/subsampled_data/subsampled_genotype_${i}_${SLURM_JOB_ID}.vcf\n")

                    # Generate .ped and .map files (for relationship matrix if needed)
                output_f.write('vcftools --vcf core_files/subsampled_data/subsampled_genotype_${i}_${SLURM_JOB_ID}.vcf --plink --out core_files/subsampled_data/subsampled_genotype_${i}_${SLURM_JOB_ID}\n')
                output_f.write(f'\n')

                    # generate bed,bim and bam files for linear model
                output_f.write('plink --file core_files/subsampled_data/subsampled_genotype_${i}_${SLURM_JOB_ID} --make-bed --out core_files/subsampled_data/subsampled_genotype_${i}_${SLURM_JOB_ID}\n')
                output_f.write(f'\n')

                    # modify the fam file to remove the -9 and add in the phenotype 
                        # maybe find a way to do this without the weird modification?
                output_f.write("cut -d' ' -f1,2,3,4,5 core_files/subsampled_data/subsampled_genotype_${i}_${SLURM_JOB_ID}.fam > core_files/subsampled_data/subsampled_genotype_${i}_${SLURM_JOB_ID}_modified.fam\n")
                output_f.write(f'\n')
                output_f.write("paste -d ' ' core_files/subsampled_data/subsampled_genotype_${i}_${SLURM_JOB_ID}_modified.fam core_files/subsampled_data/${SLURM_JOB_ID}_phenotype_list.tsv > core_files/subsampled_data/subsampled_genotype_${i}_${SLURM_JOB_ID}.fam\n")
                output_f.write(f'\n')

                    # perform the gemma association test with linear model only
                output_f.write('gemma -bfile core_files/subsampled_data/subsampled_genotype_${i}_${SLURM_JOB_ID} -lm 2 -o ${SLURM_JOB_ID}\n')
                output_f.write(f'\n')
                output_f.write(f"conda deactivate \n")
                output_f.write(f'\n')
                output_f.write(f'\n')

                # convert the output of the GWAS into a format usable for my code
                    # activate python 3 environment
                output_f.write(f'\n')
                output_f.write(f'echo "Converting GWAS format" \n')
                output_f.write(f'\n')
                output_f.write(f'conda activate python3_env\n')
                output_f.write(f'\n')
                    # convert the GWAS data to correct format
                #output_f.write('python3 batch_files/convert_GWAS_format.py -i output/${SLURM_JOB_ID}_phenotype_list.assoc.clean.txt -o output_files/leaf_ionome_${phenotype}_GWAS_${i}_${SLURM_JOB_ID}.csv\n')
                output_f.write('python3 batch_files/convert_GWAS_format.py -i output/${SLURM_JOB_ID}.assoc.txt -o output_files/leaf_ionome_${phenotype}_GWAS_${i}_${SLURM_JOB_ID}.csv\n')
                output_f.write(f'\n')
                output_f.write(f'conda deactivate\n')
                output_f.write(f'\n')
                output_f.write(f'\n')
                

                #//////////////////////////////////////////////////////////////////////////////////
                # DEPRICATED CODE /////////////////////////////////////////
                #activate gwas environment
                #bash_script_output.write(f'conda activate gwas_env\n')

                #run gwas on subsampled data
                #bash_script_output.write('pygwas run -t none -a amm -g core_files/SNP_Matrix/full_imputed_SNP_MATRIX -k core_files/SNP_Matrix/full_imputed_SNP_MATRIX/kinship_ibs_binary_mac5.h5py -o output_files/${phenotype}_GWAS_${i}_${SLURM_JOB_ID}.csv core_files/subsampled_phenotype_${i}_${SLURM_JOB_ID}.csv\n')
                
                # exit gwas environment
                #bash_script_output.write(f'conda deactivate\n')
                #//////////////////////////////////////////////////////////////////////////////////
                #//////////////////////////////////////////////////////////////////////////////////

                # # remove the subsampled phenotype file to save on storage
                # bash_script_output.write('rm core_files/subsampled_phenotype_${i}_${SLURM_JOB_ID}.csv\n')
                # bash_script_output.write('echo "GWAS for: ${SLURM_JOB_ID} finished" \n')
                ########################################################
                ######################################################################

                ######## TEMP DISABLED WHILE WORKING ON GIFT FIRST ####################
                ####################################################################
                # # calling for the python script to make the R script for GWAS and GIFT manhattan plots
                output_f.write(f'conda activate python3_env\n')
                output_f.write('python3 batch_files/make_r_scripts.py -m GWAS -id ${SLURM_JOB_ID} -i ${i} -p ${phenotype} -o output_files/\n')
                output_f.write('python3 batch_files/make_r_scripts.py -m GIFT -id ${SLURM_JOB_ID} -i ${i} -p ${phenotype} -o output_files/\n')
                
                # # # updated to NOT run the R script to save space (image generation prevented)
                # # bash_script_output.write(f'conda deactivate\n')
                # # bash_script_output.write(f'conda activate r_env\n')
                # # #bash_script_output.write('Rscript output_files/${SLURM_JOB_ID}_${i}_${phenotype}.R\n')

                output_f.write(f'conda deactivate\n')

                ######################################################################
                #### REMOVE NON-RESULT FILES TO SAVE ON STORAGE #######

                output_f.write('rm output/${SLURM_JOB_ID}.assoc.txt\n')
                output_f.write('rm output/${SLURM_JOB_ID}.log.txt\n')
                output_f.write('\n')
                output_f.write('rm core_files/subsample_text_files/subsamples_${i}_${SLURM_JOB_ID}.txt\n')
                output_f.write('rm core_files/genotype_tracker/${SLURM_JOB_ID}_genotypes.csv\n')
                output_f.write('\n')
                    # instead of just removing all, i could compress the vcf and phenotype files
                output_f.write('rm core_files/subsampled_data/subsampled_genotype_${i}_${SLURM_JOB_ID}.*\n')
                output_f.write('rm core_files/subsampled_data/${SLURM_JOB_ID}_accession_list.txt\n')
                output_f.write('rm core_files/subsampled_data/${SLURM_JOB_ID}_phenotype_list.tsv\n')
                output_f.write('rm core_files/subsampled_data/subsampled_genotype_${i}_${SLURM_JOB_ID}_modified.fam\n')
                output_f.write('\n')
                output_f.write('\n')

                #######################################################
                ######################################################################
                    # original location is set to processing_parallel
                output_f.write('mv ${original_location} batch_files/completed_parallel/\n')
                output_f.write('\n')
                output_f.write('echo "Moved file to completed folder" \n')
                output_f.write('\n')
                output_f.write('\n')

                output_f.write(f'#end of script')
                output_f.write(f'\n')
                output_f.close()

print("====================\n")
print("subsample_script_setup.py has finished!")
print("====================\n")
# end of script