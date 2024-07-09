# import modin.pandas as pandas
# import ray

import pandas
from math import log10

# ray.init(_plasma_directory="/tmp", object_store_memory=10000000000) # setting to disable out of core in Ray and obj store mem increase
# obj store memory currently at 10GB


# so i can see all the columns when testing with print
pandas.set_option('display.max_columns',None)
pandas.options.display.max_columns=None

# IDEA 3 (updated R)
def calc_thresholds(phenotype,subsample_number,pval_type,threshold_df):
    print("Entered FUNCTION: calc thresholds",flush=True)
    
    threshold_df= pandas.read_csv(f'{PATH_TO_MAIN}output_files/R_DATA/THRESHOLDS.csv')

    # Read the appropriate CSV file
    if pval_type=="AVERAGE_ABS_THETA": 
        pass
        ### No need to calc threshold with theta (no current way to do this)

    else:
          # for GWAS data....
        if pval_type=="AVERAGE_P":

            #csv_df = pandas.read_csv(f"{PATH_TO_MAIN}output_files/R_DATA_FILTERED/{phenotype}_GWAS_{subsample_number}_ALL.csv")

            # CHANGED back to R_DATA since thats where the filtered GWAS data goes now
            csv_df = pandas.read_csv(f"{PATH_TO_MAIN}output_files/R_DATA/{phenotype}_GWAS_{subsample_number}_ALL.csv")


        # for GIFT data (PSNP4 and 5 specifically)
        else:

            # fetch the data of the csv for the current phenotype and current method (GIFT)
            csv_df = pandas.read_csv(f"{PATH_TO_MAIN}output_files/R_DATA/{phenotype}_GIFT_{subsample_number}_ALL.csv")
            # calculate threshold if it isnt abs theta (since abs theta has no threshold yet)

        print(f"Calculating for pval type: {pval_type}",flush=True)

        #Calculate the BHY threshold
        m = len(csv_df)
        m=float(m)

        #R_out.write(f'gwasResults <- csv_data[order(csv_data${pval_type}),]\n')
        s=1.0

        pvals_list = csv_df[f'{pval_type}']

        # sort in order of ascension 
        pvals_list=sorted(pvals_list)

        del csv_df

        print("Pval list looks like: ", pvals_list[0:5],flush=True)

        for i,p in enumerate(pvals_list):
            i+=1

            if i>2: # different to GIFT python code...

                s=s+1/(i-1)

            # 0.05 = the FHD threshold
            thes_pval = ((i + 1.0) / m) * 0.05 / s
            if p> thes_pval:
                break

        bhy_thres = thes_pval
        transformed_bhy_thres = (-log10(thes_pval))

        bonferroni_thres = -log10(0.05/int(subsample_number))

        # export these values to the dataframe
        # first make new row for the data
        new_row_BF = pandas.Series({'PHENOTYPE':phenotype,
                    'SUBSAMPLE_NUM':subsample_number,
                    'PVAL_TYPE':pval_type,
                    'THRESHOLD_TYPE':'BF',
                    'THRESHOLD_VALUE':bonferroni_thres
                    })
        

        
        new_row_BHY= pandas.Series({'PHENOTYPE':phenotype,
                    'SUBSAMPLE_NUM':subsample_number,
                    'PVAL_TYPE':pval_type,
                    'THRESHOLD_TYPE':'BHY',
                    'THRESHOLD_VALUE':transformed_bhy_thres
                    })
        
        
        threshold_df=pandas.concat([threshold_df,new_row_BF.to_frame().T],ignore_index=True)
        threshold_df=pandas.concat([threshold_df,new_row_BHY.to_frame().T],ignore_index=True)
        threshold_df.to_csv(f'{PATH_TO_MAIN}output_files/R_DATA/THRESHOLDS.csv',header=True,index=False)
    # end of function




#packages
PATH_TO_MAIN = "/gpfs01/home/mbysh17/"

subsample_num_list=[200,400,600,800,1000] # can later update this to read from earlier scripts or something

phenotype_list=["Mo98","Na23"] # can later update this to read from the phenotype text file

pvals=["AVERAGE_P","AVERAGE_PSNP4","AVERAGE_PSNP5","AVERAGE_ABS_THETA"] # these should always be the same 4 types

# make dataframe that will house all the threshold data
threshold_df = pandas.DataFrame(columns=['PHENOTYPE', 'SUBSAMPLE_NUM', 'PVAL_TYPE','THRESHOLD_TYPE','THRESHOLD_VALUE'])
threshold_df.to_csv(f'{PATH_TO_MAIN}output_files/R_DATA/THRESHOLDS.csv',header=True,index=False)

# R script creation and running
for phenotype in phenotype_list:

    for subsample_number in subsample_num_list:

        for pval_type in pvals:

            # Example variable inputs are as follows: 
            # example 1 Mo98,200,AVERAGE_P
            # example 2 Mo98,200,AVERAGE_PSNP4
            # ...
            # example x Mo98,400,AVERAGE_P
            # ....
            calc_thresholds(phenotype,subsample_number,pval_type,threshold_df)

            

# export the threshold table to R_DATA



print("Threshold Calcs finished",flush=True)