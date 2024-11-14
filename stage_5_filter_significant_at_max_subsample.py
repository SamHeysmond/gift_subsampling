
# packages
# import modin.pandas as pandas
# import ray
import pandas
import math
import os

# testing
from multiprocessing import Pool

# functions


# 4 CPUS and 30 GB for now OUT OF (32 CPUs and 70GB )
# os.environ["MODIN_CPUS"] = "4"

# memory_for_ray= 30000000000
# ray.init(_plasma_directory="/tmp", object_store_memory=memory_for_ray,num_cpus=4) # setting to disable out of core in Ray and obj store mem increase
		

# change these functions to import from a main py file later
def fetch_subsample_numbers_list(subsample_numbers_list_file):

    subsample_num_list=[]

    temp_file=open(f"{subsample_numbers_list_file}","r")

    for line in temp_file:

        subsample_number=line.replace('\n','')

        subsample_num_list.append(subsample_number)
    
    print("Subsample number list fetched: ",flush=True)
    print(subsample_num_list,flush=True)

    return subsample_num_list

def fetch_phenotype_list(phenotype_list_file):
    phenotypes_list=[]
    temp_file=open(f"{phenotype_list_file}","r")

    for line in temp_file:

        subsample_number=line.replace('\n','')

        phenotypes_list.append(subsample_number)
    
    print("Phenotype list fetched: ",flush=True)
    print(phenotypes_list,flush=True)

    return phenotypes_list

def log_function(value):
    return -(math.log(value,10))

# def f(x):
    #return x * x
def gather_function(job_ID):
    #does this for job_ID in list_of_job_ids_to_gather:

        #temp disable
    temp_df = pandas.read_csv(f'{PATH_TO_MAIN}core_files/genotype_tracker/genotypes_{job_ID}.csv')
    #print("temp df read: ",flush=True)
    # print(temp_df,flush=True)

    # filter based on the position data gathered
        # not sure if i usedd this or not
    #temp_df_filtered = temp_df.filter(items=position_columns_filter)
    # works but gets slow above 400
        # temp_df_filtered = temp_df.loc[:, position_columns_filter]
    # testing drop method instead (WORKS- temp disable)
    temp_df_filtered=temp_df.drop(columns=[col for col in temp_df.columns if col not in position_columns_filter])

    # print("temp df FILTERED: ",flush=True)

    # save the filtered data to a folder
    # temp disable
    temp_df_filtered.to_csv(f"{PATH_TO_MAIN}core_files/filtered_genotype_tracker/filtered_genotypes_{job_ID}_{subsample_number}.csv",header=True,index=False)

    ### BATCH FILE 1 CREATION
    # saved to batch_files/stage_5_prerun
    batch1_out = open(f"{PATH_TO_MAIN}batch_files/stage_5_prerun/batch1_{job_ID}.sh","w")

    # necessary start to the file
    batch1_out.write(f'#!/bin/bash\n')
    batch1_out.write(f'#SBATCH --partition=defq\n')
    batch1_out.write(f'#SBATCH --nodes=1\n')
    batch1_out.write(f'#SBATCH --ntasks=1\n')
    batch1_out.write(f'#SBATCH --cpus-per-task=4\n')
    batch1_out.write(f'#SBATCH --mem=8g\n')
    batch1_out.write(f'#SBATCH --time=01:00:00\n')
    batch1_out.write(f'#SBATCH --job-name=S5_batch1\n')
    batch1_out.write(f'#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out\n')
    batch1_out.write(f'#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err\n')
    batch1_out.write(f'#SBATCH --mail-type=ALL\n')
    batch1_out.write(f'#SBATCH --mail-user=mbysh17@nottingham.ac.uk\n')
    batch1_out.write(f'#===============================\n')
    batch1_out.write(f'echo "start OF IDEA 5 subrun"\n')
    batch1_out.write(f'#change to home directory\n')
    batch1_out.write(f'cd /gpfs01/home/mbysh17\n')
    batch1_out.write(f'# source conda environments\n')
    batch1_out.write(f'source ~/.bashrc\n')
    batch1_out.write(f'conda deactivate\n')
    batch1_out.write(f'conda activate r_env\n') # in gift_env but may change to R env if needed
    # run the R script
    # batch1_out.write(f'Rscript batch_files/gift_testing_giota.R core_files/subsampled_data/subsampled_phenotype_{subsample_num}_$SLURM_JOB_ID.csv \\')
    # batch1_out.write(f'core_files/genotype_tracker/genotypes_$SLURM_JOB_ID.csv \\')
    # batch1_out.write(f'output_files/leaf_ionome_{phenotype}_whole_genome_metrics_{subsample_num}_$SLURM_JOB_ID.csv ')

    batch1_out.write(f'Rscript batch_files/get_paths.R core_files/subsampled_data/subsampled_phenotype_{subsample_number}_{job_ID}.csv \\')
    batch1_out.write(f'core_files/filtered_genotype_tracker/filtered_genotypes_{job_ID}_{subsample_number}.csv \\')
    batch1_out.write(f'core_files/theta_paths/delta_theta_paths_{subsample_number}_{job_ID}.csv ')

    batch1_out.write(f'conda deactivate\n')
    batch1_out.write(f'echo "END OF STAGE 5 subrun batch 1"\n')
    batch1_out.write(f'# end of file\n')
    batch1_out.close()


# set constant for file path
PATH_TO_MAIN = "/gpfs01/home/mbysh17/"

# fetch data
subsample_numbers_list = fetch_subsample_numbers_list(PATH_TO_MAIN+"core_files/subsample_numbers_list.txt")
phenotype_list=fetch_phenotype_list(PATH_TO_MAIN+"core_files/phenotypes_list.txt")

# get last in the list (maximum value)
    # may use "max" function instead in future since the list may be ordered wrong
max_subsample_num = max(subsample_numbers_list)
print(f"Max subsample number = {max_subsample_num}",flush=True)

# obtain threshold data
threshold_data =pandas.read_csv(f'{PATH_TO_MAIN}output_files/R_DATA/THRESHOLDS.csv')

print(f"Threshold data: ",flush=True)
print(threshold_data,flush=True)

# must repeat whole process for each phenotype
for phenotype in phenotype_list:

    print(f"phenotype =  {phenotype}",flush=True)

    # fetch csv for the _ALL file for current phenotype 
    # {phenotype}_{GIFT/GWAS}_{subsample_num}_ALL.csv
    phenotype_data_at_max_subsample= pandas.read_csv(f'{PATH_TO_MAIN}output_files/R_DATA/{phenotype}_GIFT_{max_subsample_num}_ALL.csv')

    # subset the threshold data to get the threshold at max subsamples
    subsetted_threshold_data=threshold_data[(threshold_data['PHENOTYPE']==phenotype) & 
                        (threshold_data['SUBSAMPLE_NUM']==int(max_subsample_num)) & 
                        (threshold_data['PVAL_TYPE']=="AVERAGE_PSNP8") &
                        (threshold_data['THRESHOLD_TYPE']=="NNNH")
                        ]

    print("Subsetted threshold data",flush=True)
    print(subsetted_threshold_data,flush=True)

    # obtain threshold value specifically
    subsetted_threshold_value = subsetted_threshold_data.iloc[0]["THRESHOLD_VALUE"]

    # apply log scale to the data
    phenotype_data_at_max_subsample["AVERAGE_PSNP8"]= phenotype_data_at_max_subsample["AVERAGE_PSNP8"].apply(log_function)

    # filter to keep only results that lie above the threshold
    phenotype_data_at_max_subsample=phenotype_data_at_max_subsample.loc[phenotype_data_at_max_subsample['AVERAGE_PSNP8']>=subsetted_threshold_value]

    # take the top 20 of these
    # cropped_pheno_data = phenotype_data_at_max_subsample.head(20).copy()
    # print("Cropped pheno data:", flush=True)
    # print(cropped_pheno_data,flush=True)

    #instead take top 20 ordered by value
        # order by PVALUE (highest first since they have been log transformed already)
    phenotype_data_at_max_subsample.sort_values(by="AVERAGE_PSNP8",ascending=False,inplace=True)
    cropped_pheno_data = phenotype_data_at_max_subsample.head(20).copy()
    print("Cropped pheno data:", flush=True)
    print(cropped_pheno_data,flush=True)

    # concatonate chromosome and position columns
    cropped_pheno_data.loc[:, 'POSITION_DATA'] = cropped_pheno_data['CHR'].astype(str) + ":" + cropped_pheno_data["POS"].astype(str) 


    # get the positions data from this dataframe so that the position data follows format CHROMOSOME:BASEPAIRS e.g. 1:3030 
        # save to a list for filtering the temp_df soon

    position_columns_filter=[]
    position_columns_filter=list(cropped_pheno_data['POSITION_DATA'])
    # print("position_columns_filter ",flush=True)
    # print(position_columns_filter,flush=True)

    # output the position data to csv
    cropped_pheno_data.to_csv(f'{PATH_TO_MAIN}core_files/{phenotype}_position_data_keep.csv',header=True,index=False)

    # for each of the subsample levels- gather then filter all the genotype data and put the filtered data into new folder
    for subsample_number in subsample_numbers_list:

        # get and process list of all significant SNPs for current subsample level if there are any
        all_snps_at_subsample_level= pandas.read_csv(f'{PATH_TO_MAIN}output_files/R_DATA/{phenotype}_GIFT_{subsample_number}_ALL.csv')

        # subset the threshold data to get the threshold at max subsamples
        subsetted_threshold_data=threshold_data[(threshold_data['PHENOTYPE']==phenotype) & 
                        (threshold_data['SUBSAMPLE_NUM']==int(subsample_number)) & 
                        (threshold_data['PVAL_TYPE']=="AVERAGE_PSNP8") &
                        (threshold_data['THRESHOLD_TYPE']=="NNNH")
                        ]

        # obtain threshold value specifically
        subsetted_threshold_value = subsetted_threshold_data.iloc[0]["THRESHOLD_VALUE"]

        # apply log scale to the data
        all_snps_at_subsample_level["AVERAGE_PSNP8"]= all_snps_at_subsample_level["AVERAGE_PSNP8"].apply(log_function)

        # filter to keep only results that lie above the threshold
        significant_snps=all_snps_at_subsample_level.loc[all_snps_at_subsample_level['AVERAGE_PSNP8']>=subsetted_threshold_value].copy()
        
        if len(significant_snps)==0:
            # input default values into row 0
            significant_snps.loc[0]=[0,0,subsample_number,0,0,0,0,0]
            # add the new column but give it value of 0
            significant_snps['POSITION_DATA'] = 0

        else:
            #format a new column for position info comparison
            significant_snps.loc[:, 'POSITION_DATA'] = significant_snps['CHR'].astype(str) + "_" + cropped_pheno_data["POS"].astype(str) 

        # save to csv for stage5_3.py to read later
        significant_snps.to_csv(f"{PATH_TO_MAIN}output_files/significant_SNP_data/GIFT_SIGNIFICANT_{subsample_number}.csv",header=True,index=False)

        # temp disable
        # if int(subsample_number)!=400 and int(subsample_number)!=600 and int(subsample_number)!=999:
        #     # skip for these two subsample numbers
        #     print(f"Skipped level: {subsample_number}",flush=True)
        #     continue

        # get JOB list .csv which contains the ID and subsample num with phenotype for each run
        job_list_df = pandas.read_csv(f'{PATH_TO_MAIN}core_files/JOB_LIST.csv',header=0)

        print("Job list df",flush=True)
        print(job_list_df)

        # filter it for the specific subsample number and phenotype in question
        filtered_job_list_df = job_list_df[(job_list_df['SUBSAMPLE_N'].astype(int)==int(subsample_number)) & 
                                            (job_list_df['PHENOTYPE'].astype(str)==str(phenotype))].copy()
        
        print("filtered_job_list_df",flush=True)
        print(filtered_job_list_df,flush=True)

        list_of_job_ids_to_gather = list(filtered_job_list_df["JOB_ID"])
        print("list of job ids...: ",flush=True)
        print(list_of_job_ids_to_gather,flush=True)

        del job_list_df
        del filtered_job_list_df
        # gather list of files that follow the naming protocol 
            # each ID within SUBSAMPLE_LEVEL = 400 etc
        if int(subsample_number)<=200:
            print("Entering <=200 section",flush=True)
            with Pool(50) as pool: # changed from 50
                pool.map(gather_function, list_of_job_ids_to_gather)

        elif int(subsample_number)==400:
            print("Entering 400 section",flush=True)
            with Pool(50) as pool:
                pool.map(gather_function, list_of_job_ids_to_gather)

        elif int(subsample_number)==600:
            print("Entering 600 section",flush=True)
            with Pool(25) as pool: # reduced to 20 from 25
                pool.map(gather_function, list_of_job_ids_to_gather)

        elif int(subsample_number)==int(max_subsample_num):
            print("Entering and max section",flush=True)

            import modin.pandas as pandas
            import ray
            # use modin version and set up ray instance etc
            # (25 CPUS and 105 GB) for now OUT OF (51 CPUs and 500GB )
            os.environ["MODIN_CPUS"] = "25" #reduced to 20 from 25

            memory_for_ray= 150000000000
            ray.init(_plasma_directory="/tmp", object_store_memory=memory_for_ray,num_cpus=25) # setting to disable out of core in Ray and obj store mem increase

            for job_id in list_of_job_ids_to_gather:
                gather_function(job_id)

        print(f"Pool stuff done for subsample {subsample_number}",flush=True)

        

print("stage_5_filter_significant_at_max_subsample.py finished",flush=True)

# end of script