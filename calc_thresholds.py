# import packages
import pandas
from math import log10

# so i can see all the columns when testing with print
pandas.set_option('display.max_columns',None)
pandas.options.display.max_columns=None

# set constant for file path
PATH_TO_MAIN = "/gpfs01/home/mbysh17/"

# Function to calculate threshold data
def calc_thresholds(phenotype,subsample_number,pval_type,threshold_df):

    # Read the appropriate CSV file
    threshold_df= pandas.read_csv(f'{PATH_TO_MAIN}output_files/R_DATA/THRESHOLDS.csv')
   
    if pval_type=="AVERAGE_ABS_THETA": 
        pass
        ### No need to calc threshold with theta (no current way to do this)
    else:
          # for GWAS data....
        if pval_type=="AVERAGE_P":

            # read the compiled GWAS data file from R_DATA folder
            csv_df = pandas.read_csv(f"{PATH_TO_MAIN}output_files/R_DATA/{phenotype}_GWAS_{subsample_number}_ALL.csv")


        # for GIFT data (PSNP4 and 5 specifically)
        else:

            # fetch the compiled data of the csv for the current phenotype and current method (GIFT)
            csv_df = pandas.read_csv(f"{PATH_TO_MAIN}output_files/R_DATA/{phenotype}_GIFT_{subsample_number}_ALL.csv")

        ##########################################
        #### Calculate the BHY threshold
        
        m = len(csv_df)
        m=float(m)
        s=1.0
        pvals_list = csv_df[f'{pval_type}']

        # sort in order of ascension 
        pvals_list=sorted(pvals_list)

        # save on memory by clearing the csv dataframe 
        del csv_df

        for i,p in enumerate(pvals_list):

            if i>2: # different to GIFT python code...

                s=s+1/(i-1)

            # 0.05 = the FHD threshold
            thes_pval = ((i + 1.0) / m) * 0.05 / s

            if p> thes_pval:

                break

        bhy_thres = thes_pval

        transformed_bhy_thres = (-log10(thes_pval))

        ############################################
        #### calculate the bonferroni threshold

        bonferroni_thres = -log10(0.05/int(subsample_number))

        # export the above values to the dataframe

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
        
        # concatonate the threshold value information to the dataframe and export as CSV file
        threshold_df=pandas.concat([threshold_df,new_row_BF.to_frame().T],ignore_index=True)
        threshold_df=pandas.concat([threshold_df,new_row_BHY.to_frame().T],ignore_index=True)
        threshold_df.to_csv(f'{PATH_TO_MAIN}output_files/R_DATA/THRESHOLDS.csv',header=True,index=False)

    # end of function


# Function to calculate threshold data (updated for BY threshold)
def calc_thresholds(phenotype,subsample_number,pval_type,threshold_df):

    # Read the appropriate CSV file
    threshold_df= pandas.read_csv(f'{PATH_TO_MAIN}output_files/R_DATA/THRESHOLDS.csv')
   
    if pval_type=="AVERAGE_ABS_THETA": 
        pass
        ### No need to calc threshold with theta (no current way to do this)
    else:
          # for GWAS data....
        if pval_type=="AVERAGE_P":

            # read the compiled GWAS data file from R_DATA folder
            csv_df = pandas.read_csv(f"{PATH_TO_MAIN}output_files/R_DATA/{phenotype}_GWAS_{subsample_number}_ALL.csv")


        # for GIFT data (PSNP4 and 5 specifically)
        else:

            # fetch the compiled data of the csv for the current phenotype and current method (GIFT)
            csv_df = pandas.read_csv(f"{PATH_TO_MAIN}output_files/R_DATA/{phenotype}_GIFT_{subsample_number}_ALL.csv")

        ##########################################
        #### Calculate the BHY threshold
        # obtain pvalues 
        pvals = csv_df[f'{pval_type}']

        # save on memory by clearing the csv dataframe 
        del csv_df

        # set up the constants
        alpha=0.05 # this is the fdr_threshold we want to control
        k = float(len(pvals))
        BY_threshold = None
        BY_threshold_rank= None

        # cumulative sum of K needed for BY procedure
        k_sum=0
        for value in range(1,int(k)+1):
            k_sum += 1/value

        # calculate different version of alpha for BY
        alpha_prime = (alpha/k_sum)
        #print(f"alpha: {alpha} /// k: {k} /// k_sum: {k_sum} /// alpha_prime: {alpha_prime}")

        # sort pvalues into ascending order
        s_pvals = sorted(pvals)
        # print("Sorted pvals are as follows")
        # print(s_pvals[0:5])

        for i, p in enumerate(s_pvals):

            i+=1 # must start at "rank 1"

            # check my values
            #print(f"Current i: {i} /// current p : {p} ///")

            # calculate adjusted value
            #BH calculation
            # adjusted_p = alpha*(i/k) 

            # BY calculation
            adjusted_p = alpha_prime*(i/k)

            # check my values
            #print(f"adjusted_p: {adjusted_p}")

            if p<adjusted_p:

                #print(f"{p} < {adjusted_p}")

                BY_threshold = adjusted_p
                BY_threshold_rank = i

                pass
            else:
                #print(f"{p} >{adjusted_p}: LOOP BROKEN!")
                break
        if BY_threshold != None:
            transformed_bhy_thres = (-log10(BY_threshold))
        else:
            transformed_bhy_thres = 0

        ############################################
        #### calculate the bonferroni threshold

        bonferroni_thres = -log10(0.05/int(subsample_number))

        # export the above values to the dataframe

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
        
        # concatonate the threshold value information to the dataframe and export as CSV file
        threshold_df=pandas.concat([threshold_df,new_row_BF.to_frame().T],ignore_index=True)
        threshold_df=pandas.concat([threshold_df,new_row_BHY.to_frame().T],ignore_index=True)
        threshold_df.to_csv(f'{PATH_TO_MAIN}output_files/R_DATA/THRESHOLDS.csv',header=True,index=False)

    # end of function





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

# end of script