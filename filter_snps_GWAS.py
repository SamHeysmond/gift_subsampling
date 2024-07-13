import argparse, os, math
import concurrent.futures
import pandas

# so i can see all the columns when testing with print
pandas.set_option('display.max_columns',None)
pandas.options.display.max_columns=None


# function to process csv file
def process_csv(csv_f):
    print("Iteration started",flush=True)
    print(f"reading file: {PATH_TO_MAIN}output_files/{csv_f}",flush=True)
    # print(f"phenotype: {phenotype}",flush=True)
    # print(f"subsample_number: {subsample_number}",flush=True)

    base_df = pandas.read_csv(f"{PATH_TO_MAIN}output_files/{csv_f}")
    
    # outer will bring in everything- if it overlaps itll be "both", otherwise itll be "left only" or "right only"
    # we want to keep the "left only" data
    print("Merge command...",flush=True)
    merged_df=pandas.merge(base_df,removal_df,how="outer",on=['chromosomes','positions'],indicator=True)

    print(merged_df.head(),flush=True)

    # Delete rows where data is in both 
    print("Drop command 1...",flush=True)
    merged_df = merged_df.drop(merged_df[merged_df['_merge'] == 'both'].index)

    # Delete rows where data is right only
    print("Drop command 2...",flush=True)
    merged_df = merged_df.drop(merged_df[merged_df['_merge'] == 'right_only'].index)

    print(merged_df.head(),flush=True)

    # remove merge column; no longer needed
    final_df=merged_df.drop(columns=['_merge'])
    final_df=final_df.dropna()
    print(final_df.head(),flush=True)

    # save to FILTERED file
    print("saving to csv...",flush=True)
    final_df.to_csv(f"{PATH_TO_MAIN}output_files/{csv_f}",header=True,index=False)

    print("Iteration finished",flush=True)


# placeholder till parseargs will work
# will implement an argument that inputs home user directory automatically
# -> or lets user decide
PATH_TO_MAIN = "/gpfs01/home/mbysh17/"

print("inside python script",flush=True)

# df containing snp locations to be removed
removal_df = pandas.read_csv(f"{PATH_TO_MAIN}core_files/output_3.table",sep='\t')

# rename columns to fit the GWAS output csv header
removal_df.rename(columns={"CHROM":'chromosomes',"POS":"positions"},inplace=True)

print("head of removal df:", removal_df.head(),flush=True)

csv_files=[]

for file in os.listdir(PATH_TO_MAIN+"output_files/csv_before_filter"):
    if file.endswith(".csv") and file.__contains__("T20")==False:
        #print(file,": ++++++++++++ ADDED ++++++++++++ !")
        csv_files.append(file)
    else:
        #print(file,": ////////// SKIPPED //////// !!")
        pass


print("////////////////////////////////////////////////",flush=True)
print("Some of csv files list",csv_files[0:5],flush=True)
print("////////////////////////////////////////////////",flush=True)


GWAS_files=[]

# split csv files into categories (PHENO, METHOD, SUBNUM)
# PROGRESS METER
print("////////////////////////////////////////////////",flush=True)
print("Fetching GWAS files...",flush=True)
print("////////////////////////////////////////////////",flush=True)
for csv_file in csv_files:
    csv_file_name = csv_file.split("_")
    #print(f"Current csv file to be split is... {csv_file}")
    #print(f"position 6 of its name is... {csv_file[6]}")

    # Mo98 GWAS lists
    if csv_file_name[2]=="Mo98" and csv_file_name[3]=='GWAS' :
        GWAS_files.append(csv_file)

    # ---
    # Na23 GWAS lists
    elif csv_file_name[2]=="Na23" and csv_file_name[3]=='GWAS':
        GWAS_files.append(csv_file) 

# PROGRESS METER
print("////////////////////////////////////////////////",flush=True)
print("Finished fetching GWAS files!",flush=True)
print("////////////////////////////////////////////////",flush=True)

print("////////////////////////////////////////////////",flush=True)
print("Some of GWAS files list",GWAS_files[0:5],flush=True)
print("////////////////////////////////////////////////",flush=True)



# filtering the GWAS csv files 
# note: this can be  turned into multiprocessing OR multithreading if you'd like

# Na23_GWAS_200_ALL.csv
# only need to do this for GWAS as already done in the GIFT code


# can change to ThreadPoolExecutor instead if uts I/O bound
# use processpoolexecutor if its CPU bound
with concurrent.futures.ThreadPoolExecutor() as executor:

    # runs the function with all items in the list as inputs
    # do so in parallel
    executor.map(process_csv, GWAS_files)

R_files=[]
# fetch all the R scripts
for file in os.listdir(PATH_TO_MAIN+"output_files"):
    if file.endswith(".R"):
        file_split=file.split("_")
        if file_split[3]=="whole":
            pass
        else:
            #print(file,": ++++++++++++ ADDED ++++++++++++ !",flush=True)
            R_files.append(file)
    else:
        #print(file,": ////////// SKIPPED //////// !!")
        pass

for R_file in R_files:
    # location where the batch scripts will be written to

    R_file_name = R_file.replace(".R","")
    R_batch=open(PATH_TO_MAIN+"batch_files/R_parallel/"+str(R_files.index(R_file))+"_"+str(R_file_name)+".sh","w")
    # necessary start to the file
    R_batch.write(f'#!/bin/bash\n')
    R_batch.write(f'#SBATCH --partition=defq\n')
    R_batch.write(f'#SBATCH --nodes=1\n')
    R_batch.write(f'#SBATCH --ntasks=1\n')
    R_batch.write(f'#SBATCH --cpus-per-task=3\n')
    R_batch.write(f'#SBATCH --mem=7g\n')
    R_batch.write(f'#SBATCH --time=1:00:00\n')
    R_batch.write(f'#SBATCH --job-name=R_subrun_filtered_rerun\n')
    R_batch.write(f'#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out\n')
    R_batch.write(f'#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err\n')
    R_batch.write(f'#SBATCH --mail-type=ALL\n')
    R_batch.write(f'#SBATCH --mail-user=mbysh17@nottingham.ac.uk\n')
    R_batch.write(f'#===============================\n')
    R_batch.write(f'echo "start of GWAS filtered rerun script"\n')
    R_batch.write(f'#change to home directory\n')
    R_batch.write(f'cd /gpfs01/home/mbysh17\n')
    R_batch.write(f'# source conda environments\n')
    R_batch.write(f'source ~/.bashrc\n')
    R_batch.write(f'conda deactivate\n')
    R_batch.write(f'conda activate r_env\n')
    R_batch.write(f'# R SCRIPT FOR GWAS FILTERED RERUN\n')
    R_batch.write(f'Rscript {PATH_TO_MAIN}output_files/{R_file_name}.R\n')
    R_batch.write(f'conda deactivate\n')
    R_batch.write(f'echo "End of GWAS filtered rerun script"\n')
    R_batch.write(f'conda deactivate\n')
    R_batch.write(f'# END OF FILE\n')
    R_batch.close()
    # end of function

print("Python filtering done",flush=True)