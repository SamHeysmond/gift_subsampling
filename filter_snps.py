import modin.pandas as pandas
import ray


ray.init(_plasma_directory="/tmp", object_store_memory=100000000000) # setting to disable out of core in Ray and obj store mem increase
# obj store memory currently at 100GB


# so i can see all the columns when testing with print
pandas.set_option('display.max_columns',None)
pandas.options.display.max_columns=None


# placeholder till parseargs will work
# will implement an argument that inputs home user directory automatically
# -> or lets user decide
PATH_TO_MAIN = "/gpfs01/home/mbysh17/"



subsample_num_list=[200,400,600,800,1000] # can later update this to read from earlier scripts or something

phenotype_list=["Mo98","Na23"] # can later update this to read from the phenotype text file


method_type_list=["GIFT","GWAS"]

# df containing snp locations to be removed
removal_df = pandas.read_csv(f"{PATH_TO_MAIN}core_files/output_3.table",sep='\t')
removal_df.rename(columns={"CHROM":'CHR'},inplace=True)

# filtering the average csv files 
# note: this can be  turned into multiprocessing OR multithreading if you'd like

# Na23_GWAS_200_ALL.csv
# only need to do this for GWAS as already done in the GIFT code

for phenotype in phenotype_list:

    for subsample_number in subsample_num_list:

        # Example variable inputs are as follows: 
        # example 1 Mo98,200,AVERAGE_P
        # example 2 Mo98,200,AVERAGE_PSNP4
        # ...
        # example x Mo98,400,AVERAGE_P
        # ....
        print("Iteration started",flush=True)
        print(f"phenotype: {phenotype}",flush=True)
        print(f"subsample_number: {subsample_number}",flush=True)

        base_df = pandas.read_csv(f"{PATH_TO_MAIN}output_files/R_DATA/{phenotype}_GWAS_{subsample_number}_ALL.csv")
        
        # outer will bring in everything- if it overlaps itll be "both", otherwise itll be "left only" or "right only"
        # we want to keep the "left only" data
        merged_df=pandas.merge(base_df,removal_df,how="outer",on=['CHR','POS'],indicator=True)

        print(merged_df.head,flush=True)

        # Delete rows where data is in both 
        merged_df = merged_df.drop(merged_df[merged_df['_merge'] == 'both'].index)

        # Delete rows where data is right only
        merged_df = merged_df.drop(merged_df[merged_df['_merge'] == 'right_only'].index)

        print(merged_df.head,flush=True)

        # remove merge column; no longer needed
        final_df=merged_df.drop(columns=['_merge'])
        final_df=final_df.dropna()
        print(final_df.head,flush=True)

        # save to FILTERED file
        final_df.to_csv(f"{PATH_TO_MAIN}output_files/R_DATA_FILTERED/{phenotype}_GWAS_{subsample_number}_ALL.csv",header=True,index=False)

        print("Iteration finished",flush=True)
        

print("Python filtering done",flush=True)