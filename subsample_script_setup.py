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

# packages
import argparse
import math
import pandas
import os 

def fetch_subsample_numbers_list(subsample_numbers_list_file):

    subsample_num_list=[]

    temp_file=open(f"{subsample_numbers_list_file}","r")

    for line in temp_file:

        subsample_number=line.replace('\n','')

        subsample_num_list.append(subsample_number)
    
    print("Subsample number list fetched: ",flush=True)
    print(subsample_num_list,flush=True)

    return subsample_num_list

parser=argparse.ArgumentParser(description="W.I.P")

parser.add_argument('-subsampleFile', 
                    type=str, 
                    metavar='location of subsample file', 
                    required=True, 
                    help='location of the file containing list of all sample numbers'
                    )

# stores input data and parses them
args= parser.parse_args() 

phenotypes_input=open("core_files/phenotypes_list.txt","r")

# set the subsample lists to what is in the file
subsample_list=fetch_subsample_numbers_list(args.subsampleFile)

cpu_num = 4
cores=cpu_num-1


# STEPS FOR MATRYOSHKA CODE
# set the data sets for each subsample number by 
    # 1) subsampling from each one and saving to files like subsample.py would do
        # see output 1,2 and 3 in subsample.py
    # 2) then set the output of these to same format as subsample .py but using the ID name for the run (gives unique run IDs
        # e.g. using a mv command but with the variable name to the SLURM job ID in named format
    # 3) ensure copynum range is 1 to 2 so that theres only 1 copy at each subsample level
        # to get repeats i will need to run the WHOLE set of scripts multiple times (e.g. 10 times over)

#1)
#set input file
sample_list_from_vcf = pandas.read_csv('core_files/all_vcf_samples.txt' , sep=',', usecols=[0], names=['Accession_ID'], header=None)

print("sample_list_from_vcf : ",flush=True)
print(sample_list_from_vcf,flush=True)

# open phenotype data file (set to Mo98 for this)
sample_list_from_phenotype_file = pandas.read_csv("core_files/leaf_phenotype.csv", sep=',', usecols=['Accession_ID',"Mo98"])
print("sample_list_from_phenotype_file : ",flush=True)
print(sample_list_from_phenotype_file,flush=True)

# merge data where Accession_ID is in both the vcf and phenotype files
consensus_ID_dataframe = (sample_list_from_vcf.reset_index(drop=True)[["Accession_ID"]].merge(sample_list_from_phenotype_file.reset_index(drop=True), on=['Accession_ID'], how='inner',left_index=False,right_index=False))

print("consensus_ID_dataframe_with_phenotypes : ",flush=True)
print(consensus_ID_dataframe,flush=True)

print(f"Length of consensus dataframe (should be 1000 or 1003 ish): {len(consensus_ID_dataframe)}", flush=True)

# now subsample for MATRYOSHKA for each number in subsample num list
for subsample_number in reversed(subsample_list):

    # if its the max number it'll take from consensus, not from itself
    if subsample_number==max(subsample_list):
        # take from consensus
        subsampled_dataframe = consensus_ID_dataframe.sample(n=int(subsample_number))

    else:
        # take from previous (subsampling from itself)
        subsampled_dataframe=subsampled_dataframe.sample(n=int(subsample_number))

    #### 
    
    print("Subsampled_dataframe before reset index: ",flush=True)
    print(subsampled_dataframe,flush=True)

    subsampled_dataframe.reset_index(drop=True, inplace=True)
    print("Subsampled_dataframe after reset index: ",flush=True)
    print(subsampled_dataframe,flush=True)

    # open file for writing all the IDs that were subsampled for this particular run
    subsampled_IDs_path = (f'core_files/subsample_text_files/subsamples_{subsample_number}_1234.txt')

    print("subsampled_dataframe data types before sorting: ",flush=True)
    print(subsampled_dataframe.dtypes,flush=True)

    # convert phenotype to numeric and sort in ascending order
    subsampled_dataframe[[f'Mo98']] = subsampled_dataframe[[f'Mo98']].apply(pandas.to_numeric)
    subsampled_dataframe = subsampled_dataframe.sort_values(by='Mo98',ascending=True)

    print("subsampled_dataframe after sorting: ",flush=True)
    print(subsampled_dataframe,flush=True)

    print("subsampled_dataframe data types after sorting: ",flush=True)
    print(subsampled_dataframe.dtypes,flush=True)

    # OUTPUT 1 
        # ordered subsample list of just the IDs 
    subsampled_dataframe.to_csv(subsampled_IDs_path,columns=['Accession_ID'],index=False, header=False)
        # NEEDS MV COMMAND BELOW

    # OUTPUT 2
        #create subsampled vcf
        # set the ID to 1234-> this will be changed in the single run scripts that are temporarily active
    os.system(f'bcftools view --samples-file core_files/subsample_text_files/subsamples_{subsample_number}_1234.txt core_files/FINAL.vcf > core_files/subsampled_data/subsampled_genotype_{subsample_number}_1234.vcf')
        # NEEDS MV COMMAND BELOW

    # OUTPUT 3
        # subsample of the accessions and phenotype data together
    subsampled_dataframe.to_csv(f"core_files/subsampled_data/subsampled_phenotype_{subsample_number}_1234.csv",header=True,index=False)
        # NEEDS MV COMMAND BELOW


#loop through each phenotype with the following settings
for phenotype in phenotypes_input:
    phenotype = phenotype.replace('\n','')
    switch = True

    # loop for each different amount of samples
    for subsample_num in subsample_list:

        # give appropriate memory and time according to the job's subsample number... 
        # ...(bigger sample -> more resources) adjust scaling when necessary
        # memory_limit=math.ceil(2*(int(subsample_num)/100) + (20))
        # future update subsample number list into dictionary with memory as key/value alongside subsample
            # or as data frame

        # 200 -> 22 // 600 -> 30 // 1000 -> 38 (GB)
        # allow 90% of the total memory limit (in case of spills etc?)

        hours_given=math.ceil(int(subsample_num)/100)
        # 200 -> 4 hours
        # make 100 copies of each file 
        # for copynum in range(1,101):
        for copynum in range(1,2): # temp changed
            if subsample_num=='999' and switch==False:
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
                output_f.write(f'#SBATCH --partition=hmemq\n')
                output_f.write(f'#SBATCH --nodes=1\n')
                output_f.write(f'#SBATCH --ntasks=1\n')
                output_f.write(f'#SBATCH --cpus-per-task={cpu_num}\n')
                
                if subsample_num=='999':
                    # only need to do 1 script for 999 since all the same samples (max sample num)
                    switch = False
                    output_f.write(f'#SBATCH --mem=40g\n')
                    memory_limit=40
                if subsample_num=='800':
                    output_f.write(f'#SBATCH --mem=40g\n')
                    memory_limit=40
                if subsample_num=='600':
                    output_f.write(f'#SBATCH --mem=32g\n')
                    memory_limit=32
                if subsample_num=='400':
                    output_f.write(f'#SBATCH --mem=28g\n')
                    memory_limit=28
                if subsample_num=='200':
                    output_f.write(f'#SBATCH --mem=24g\n')
                    memory_limit=24
                if subsample_num=='100':
                    output_f.write(f'#SBATCH --mem=20g\n')
                    memory_limit=20
                if subsample_num=='50':
                    output_f.write(f'#SBATCH --mem=16g\n')
                    memory_limit=16

                core_mem_limit = math.ceil(0.9*memory_limit)
                    
                output_f.write(f'#SBATCH --time={hours_given}:00:00\n')             
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
                output_f.write(f'core_mem_limit={core_mem_limit}\n') 
                output_f.write(f'\n')
                # should contain "processing_parallel" in the file path
                output_f.write("original_location=$(scontrol show job $SLURM_JOB_ID | awk -F= '/Command=/{print $2}')\n")
                output_f.write(f'\n')

                # SPACING
                output_f.write(f'\n')
                output_f.write(f'\n')
                # SPACING

                # temp for MATRYOSHKA
                output_f.write(f'echo "Moving files for MATRYOSHKA"\n')
                # MV output 1
                output_f.write('mv core_files/subsample_text_files/subsamples_'+str(subsample_num)+'_1234.txt core_files/subsample_text_files/subsamples_'+str(subsample_num)+'_${SLURM_JOB_ID}.txt\n')
                output_f.write(f'\n')

                # MV output 2
                output_f.write('mv core_files/subsampled_data/subsampled_genotype_'+str(subsample_num)+'_1234.vcf core_files/subsampled_data/subsampled_genotype_'+str(subsample_num)+'_${SLURM_JOB_ID}.vcf\n')
                output_f.write(f'\n')

                # MV output 3
                output_f.write('mv core_files/subsampled_data/subsampled_phenotype_'+str(subsample_num)+'_1234.csv core_files/subsampled_data/subsampled_phenotype_'+str(subsample_num)+'_${SLURM_JOB_ID}.csv\n')
                output_f.write(f'\n')

                # SPACING
                output_f.write(f'\n')
                output_f.write(f'\n')
                # SPACING

                #job ID updating to the csv
                output_f.write(f'#append the JOB ID and subsample number (i) to a file\n')
                output_f.write('echo "${SLURM_JOB_ID},${i},${phenotype}" >> core_files/JOB_LIST.csv \n')
                output_f.write(f'echo "Subsample number: {subsample_num}"  \n')
                output_f.write(f'echo "memory assigned: {memory_limit}"  \n')
                output_f.write(f'echo "GIFT core mem lim: {core_mem_limit}"  \n')
                output_f.write(f'\n')

                # main code to run the subsampling (TEMP DISABLED FOR MATRYOSHKA)
                output_f.write(f'# run the subsample script to get the appropriate number of subsamples\n')
                # output_f.write('python3 batch_files/subsample.py -p core_files/leaf_phenotype.csv \\')
                # output_f.write(f'-s core_files/all_vcf_samples.txt \\')
                # output_f.write(f'-n $i \\')
                # output_f.write(f'-t $phenotype \\')
                # output_f.write('-op core_files/subsampled_data/subsampled_phenotype_${i}_${SLURM_JOB_ID}.csv \\')
                # output_f.write('-og core_files/subsampled_data/subsampled_genotype_${i}_${SLURM_JOB_ID}.vcf \\')
                # output_f.write(f'-ri $SLURM_JOB_ID ')
                # output_f.write(f'\n')

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
                # output_f.write('python3 -W ignore batch_files/GIFT_CORE.py -v core_files/subsampled_data/subsampled_genotype_${i}_${SLURM_JOB_ID}.vcf -f core_files/subsampled_data/subsampled_phenotype_${i}_${SLURM_JOB_ID}.csv -p ${phenotype} -o output_files/${phenotype}_whole_genome_metrics_${i}_${SLURM_JOB_ID}.csv -id ${SLURM_JOB_ID} -s ${i} -m ${core_mem_limit}\n')
                
                output_f.write(f'python3 -W ignore batch_files/GIFT_CORE.py -v core_files/subsampled_data/subsampled_genotype_{subsample_num}_$SLURM_JOB_ID.vcf \\')
                output_f.write(f'\n')
                output_f.write(f'-f core_files/subsampled_data/subsampled_phenotype_{subsample_num}_$SLURM_JOB_ID.csv \\')
                output_f.write(f'\n')
                output_f.write(f'-p {phenotype} \\')
                output_f.write(f'\n')
                # depricated
                output_f.write(f'-o output_files/{phenotype}_whole_genome_metrics_{subsample_num}_$SLURM_JOB_ID.csv \\')
                output_f.write(f'\n')
                output_f.write(f'-id $SLURM_JOB_ID \\')
                output_f.write(f'\n')
                # aRGS.S depricated
                output_f.write(f'-s {subsample_num} \\')
                output_f.write(f'\n')
                output_f.write(f'-m {core_mem_limit} ')
                output_f.write(f'\n')
                output_f.write(f'\n')

                # exit gift environment for python script
                output_f.write(f'conda deactivate\n')

                # enter the R environment to continue the GIFT calculations
                output_f.write(f'\n')
                output_f.write(f'conda activate r_env\n')
                output_f.write(f'\n')
                output_f.write(f'echo "Calculating PSNP8 with gift_testing_giota.R" \n')
                output_f.write(f'\n')

                # output_f.write(f'Rscript batch_files/gift_testing_giota.R core_files/subsampled_data/subsampled_phenotype_{subsample_num}_$SLURM_JOB_ID.csv core_files/genotype_tracker/$SLURM_JOB_ID_genotypes.csv output_files/leaf_ionome_{phenotype}_whole_genome_metrics_{subsample_num}_$SLURM_JOB_ID.csv {cores} \n')
                    # works
                #output_f.write(f'Rscript batch_files/gift_testing_giota.R core_files/subsampled_data/subsampled_phenotype_{subsample_num}_$SLURM_JOB_ID.csv core_files/genotype_tracker/genotypes_$SLURM_JOB_ID.csv output_files/leaf_ionome_{phenotype}_whole_genome_metrics_{subsample_num}_$SLURM_JOB_ID.csv \n')
                output_f.write(f'Rscript batch_files/gift_testing_giota.R core_files/subsampled_data/subsampled_phenotype_{subsample_num}_$SLURM_JOB_ID.csv \\')
                output_f.write(f'core_files/genotype_tracker/genotypes_$SLURM_JOB_ID.csv \\')
                output_f.write(f'output_files/leaf_ionome_{phenotype}_whole_genome_metrics_{subsample_num}_$SLURM_JOB_ID.csv ')
                output_f.write(f'\n')
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
                #output_f.write('gemma -bfile core_files/subsampled_data/subsampled_genotype_${i}_${SLURM_JOB_ID} -lm 2 -o ${SLURM_JOB_ID}\n')
                # maf is 0.05 since 0.01 of 1000 is 10 and 0.05 of 200 is 10
                #output_f.write('gemma -bfile core_files/subsampled_data/subsampled_genotype_${i}_${SLURM_JOB_ID} -lm 2 -miss 0.1 -maf 0.05 -o ${SLURM_JOB_ID}\n')
                    # -miss 1 or 0 not sure yet -> aim is to allow ALL through analysis since i filtered earlier
                output_f.write('gemma -bfile core_files/subsampled_data/subsampled_genotype_${i}_${SLURM_JOB_ID} -lm 2 -miss 1 -maf 0 -hwe 0 -r2 1 -o ${SLURM_JOB_ID}\n')


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
                
                # # remove the subsampled phenotype file to save on storage
                # bash_script_output.write('rm core_files/subsampled_phenotype_${i}_${SLURM_JOB_ID}.csv\n')
                # bash_script_output.write('echo "GWAS for: ${SLURM_JOB_ID} finished" \n')
                ########################################################
                ######################################################################

                ######## RAW RUN MANHATTAN PLOTS ####################
                ####################################################################
                # # calling for the python script to make the R script for GWAS and GIFT manhattan plots
                output_f.write(f'conda activate python3_env\n')
                output_f.write('python3 batch_files/make_r_scripts.py -m GWAS -id ${SLURM_JOB_ID} -i ${i} -p ${phenotype} -o output_files/\n')
                output_f.write('python3 batch_files/make_r_scripts.py -m GIFT -id ${SLURM_JOB_ID} -i ${i} -p ${phenotype} -o output_files/\n')
                output_f.write(f'conda deactivate\n')
                # # # updated to run R scripts
                output_f.write(f'\n')
                output_f.write('echo "Running R scripts" \n')
                output_f.write(f'conda activate r_env\n')
                # temp enambled
                output_f.write('Rscript output_files/GWAS_${phenotype}_${i}_${SLURM_JOB_ID}.R\n')
                output_f.write('Rscript output_files/GIFT_${phenotype}_${i}_${SLURM_JOB_ID}.R\n')
                output_f.write(f'conda deactivate\n')
                output_f.write(f'\n')
                ######################################################################
                #### REMOVE NON-RESULT FILES TO SAVE ON STORAGE #######
                    ### PAUSED FOR TESTING
                output_f.write('\n')
                output_f.write('echo "Removing used files" \n')
                output_f.write('\n')
                output_f.write('rm output/${SLURM_JOB_ID}.assoc.txt\n')
                output_f.write('rm output/${SLURM_JOB_ID}.log.txt\n')
                output_f.write('\n')
                     # keep for now (so comment out to disable the removal)
                # output_f.write('rm core_files/subsample_text_files/subsamples_${i}_${SLURM_JOB_ID}.txt\n')
                # output_f.write('rm core_files/genotype_tracker/genotypes_${SLURM_JOB_ID}.csv\n')

                ### WRITE IN A REMOVE FOR subsampled_data/subsampled_phenotpe_{i}_{SLURM_JOB_ID}.csv

                output_f.write('\n')
                    # instead of just removing all, i could compress the vcf and phenotype files
                output_f.write('rm core_files/subsampled_data/subsampled_genotype_${i}_${SLURM_JOB_ID}.*\n')

                    # keep for now (so comment out to disable the removal)
                # output_f.write('rm core_files/subsampled_data/${SLURM_JOB_ID}_accession_list.txt\n')
                # output_f.write('rm core_files/subsampled_data/${SLURM_JOB_ID}_phenotype_list.tsv\n')
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

# end of script