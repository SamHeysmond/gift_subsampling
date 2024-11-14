# import packages
import pandas
from math import log10
import numpy
import statsmodels.stats.multitest
import argparse

# so i can see all the columns when testing with print
pandas.set_option('display.max_columns',None)
pandas.options.display.max_columns=None

# set constant for file path
PATH_TO_MAIN = "/gpfs01/home/mbysh17/"

# TEMP EDIT FOR SINGLE RUNS
apply_correction="NO"

# new function
def calc_thresholds_main(phenotype,subsample_number,pval_type,threshold_df,apply_correction):
    # Read the appropriate CSV file
    threshold_df= pandas.read_csv(f'{PATH_TO_MAIN}output_files/R_DATA/THRESHOLDS.csv')

        # for GWAS data....
    if pval_type=="AVERAGE_P" or pval_type=="pvals":

        # read the compiled GWAS data file from R_DATA folder
        # csv_df = pandas.read_csv(f"{PATH_TO_MAIN}output_files/R_DATA/{phenotype}_GWAS_{subsample_number}_ALL.csv")
        csv_df = pandas.read_csv(f"{PATH_TO_MAIN}output_files/R_DATA/{phenotype}_GWAS_{subsample_number}_UNCORRECTED.csv")

    # for GIFT data (PSNP4 and 5 specifically)
    elif pval_type=="AVERAGE_PSNP8" or pval_type=="PSNP8":

        # fetch the compiled data of the csv for the current phenotype and current method (GIFT)
        # csv_df = pandas.read_csv(f"{PATH_TO_MAIN}output_files/R_DATA/{phenotype}_GIFT_{subsample_number}_ALL.csv")
        csv_df = pandas.read_csv(f"{PATH_TO_MAIN}output_files/R_DATA/{phenotype}_GIFT_{subsample_number}_UNCORRECTED.csv")

    ##########################################
    #### 
    # obtain pvalues 
    pvals = list(csv_df[f'{pval_type}'])

    # sort pvalues into ascending order (smallest to biggest)
    s_pvals = sorted(pvals)

        # GWAS steps
    if pval_type == "AVERAGE_P" or pval_type=="pvals":

        # obtain 99% and 95% NULL HYPOTHESIS 
        NNBF= numpy.percentile(s_pvals,1)
        NFBF= numpy.percentile(s_pvals,5)

        # LOG TRANSFORM THEM
        NNBF=-log10(NNBF)
        NFBF=-log10(NFBF)

        if apply_correction=="YES":

            print("Applying correction...",flush=True)
            input_list_pvals=list(csv_df[f'{pval_type}'])

            # defaults to BH. BY is too strong here.
            rejected,corrected_pvals= statsmodels.stats.multitest.fdrcorrection(input_list_pvals,
                                                                                        alpha=0.05)

            # save corrected pvalues to the dataframe
            csv_df[f'{pval_type}'] = corrected_pvals

            # save the dataframe as a corrected version of itself
            csv_df.to_csv(f"{PATH_TO_MAIN}output_files/R_DATA/{phenotype}_GWAS_{subsample_number}_ALL.csv",header=True,index=False)
        else:
            print("NO correction...",flush=True)
        new_row_NN = pandas.Series({'PHENOTYPE':phenotype,
                    'SUBSAMPLE_NUM':subsample_number,
                    'PVAL_TYPE':pval_type,
                    'THRESHOLD_TYPE':'NNBF',
                    'THRESHOLD_VALUE':NNBF
                    })
        
        new_row_NF= pandas.Series({'PHENOTYPE':phenotype,
                    'SUBSAMPLE_NUM':subsample_number,
                    'PVAL_TYPE':pval_type,
                    'THRESHOLD_TYPE':'NFBF',
                    'THRESHOLD_VALUE':NFBF
                    })

        # GIFT steps
    elif pval_type=="AVERAGE_PSNP8" or pval_type=="PSNP8":

        # obtain 99% and 95% quantile values (before correction)
        # NNNH= numpy.percentile(s_pvals,99)
        # NFNH= numpy.percentile(s_pvals,95)
        # to do this, must obtain the SMALLEST 1% and 5% values (these will be biggest after transformed)
        NNNH= numpy.percentile(s_pvals,1)
        NFNH= numpy.percentile(s_pvals,5)

        print("Null hypotheses before log tranformation",flush=True)
        print(f"NF: {NFNH} ",flush=True)
        print(f"NN: {NNNH} ",flush=True)

        NNNH=-log10(NNNH)
        NFNH=-log10(NFNH)
        
        print("Null hypotheses after log tranformation",flush=True)
        print(f"NF: {NFNH} ",flush=True)
        print(f"NN: {NNNH} ",flush=True)

        if apply_correction=="YES":
            print("Applying correction...",flush=True)
            # correct all the pvalues using BH or BY (im using BY this time)                                                                      
            input_list_pvals=list(csv_df[f'{pval_type}'])

            rejected,corrected_pvals= statsmodels.stats.multitest.fdrcorrection(input_list_pvals,
                                                                                alpha=0.05)
 
            csv_df[f'{pval_type}'] = corrected_pvals

            # return corrected dataframe (temp disabled)
            csv_df.to_csv(f"{PATH_TO_MAIN}output_files/R_DATA/{phenotype}_GIFT_{subsample_number}_ALL.csv",header=True,index=False)
        else:
            print("NO correction...",flush=True)

        # export threshold values to the dataframe
            # first make new row for the data
        new_row_NN = pandas.Series({'PHENOTYPE':phenotype,
                    'SUBSAMPLE_NUM':subsample_number,
                    'PVAL_TYPE':pval_type,
                    'THRESHOLD_TYPE':'NNNH',
                    'THRESHOLD_VALUE':NNNH
                    })
        
        new_row_NF= pandas.Series({'PHENOTYPE':phenotype,
                    'SUBSAMPLE_NUM':subsample_number,
                    'PVAL_TYPE':pval_type,
                    'THRESHOLD_TYPE':'NFNH',
                    'THRESHOLD_VALUE':NFNH
                    })

    
    # concatonate the threshold value information to the dataframe and export as CSV file
        # threshold_df=pandas.concat([threshold_df,new_row_BF.to_frame().T],ignore_index=True)
        # threshold_df=pandas.concat([threshold_df,new_row_BHY.to_frame().T],ignore_index=True)
    threshold_df=pandas.concat([threshold_df,new_row_NN.to_frame().T],ignore_index=True)
    threshold_df=pandas.concat([threshold_df,new_row_NF.to_frame().T],ignore_index=True)
    threshold_df.to_csv(f'{PATH_TO_MAIN}output_files/R_DATA/THRESHOLDS.csv',header=True,index=False)

    # end of function

def fetch_phenotype_list(phenotype_list_file):
    phenotypes_list=[]
    temp_file=open(f"{phenotype_list_file}","r")

    for line in temp_file:

        subsample_number=line.replace('\n','')

        phenotypes_list.append(subsample_number)
    
    print("Phenotype list fetched: ",flush=True)
    print(phenotypes_list,flush=True)

    return phenotypes_list


# subsample_num_list=[200,400,600,800,1000] # can later update this to read from earlier scripts or something
    # new
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


# subsample_num_list=[200,400,600,800,999]
subsample_num_list=fetch_subsample_numbers_list(args.subsampleFile)

# phenotype_list=["Mo98","Na23"] # can later update this to read from the phenotype text file
phenotype_list=fetch_phenotype_list(PATH_TO_MAIN+"core_files/phenotypes_list.txt")

# pvals=["AVERAGE_P","AVERAGE_PSNP4","AVERAGE_PSNP5","AVERAGE_ABS_THETA"] # these should always be the same 4 types
    # new
# temp swapped
pvals=["AVERAGE_P","AVERAGE_PSNP8"] 
# pvals=["pvals","PSNP8"] 

# make dataframe that will house all the threshold data
threshold_df = pandas.DataFrame(columns=['PHENOTYPE', 'SUBSAMPLE_NUM', 'PVAL_TYPE','THRESHOLD_TYPE','THRESHOLD_VALUE'])
threshold_df.to_csv(f'{PATH_TO_MAIN}output_files/R_DATA/THRESHOLDS.csv',header=True,index=False)

# R script creation and running
for phenotype in phenotype_list:

    for subsample_number in subsample_num_list:

        for pval_type in pvals:
            # Example variable inputs are as follows: 
            # example 1 Mo98,200,AVERAGE_P
            # example 2 Mo98,200,AVERAGE_PSNP8
            # calc_thresholds(phenotype,subsample_number,pval_type,threshold_df)
            calc_thresholds_main(phenotype,subsample_number,pval_type,threshold_df,apply_correction)

# export the threshold table to R_DATA

print("Thresholds calculted and data corrected. script finished",flush=True)

## OLD CODE TEMP STORE
        # new_row_BF = pandas.Series({'PHENOTYPE':phenotype,
        #             'SUBSAMPLE_NUM':subsample_number,
        #             'PVAL_TYPE':pval_type,
        #             'THRESHOLD_TYPE':'BF',
        #             'THRESHOLD_VALUE':bonferroni_thres
        #             })
        
        # new_row_BHY= pandas.Series({'PHENOTYPE':phenotype,
        #             'SUBSAMPLE_NUM':subsample_number,
        #             'PVAL_TYPE':pval_type,
        #             'THRESHOLD_TYPE':'BHY',
        #             'THRESHOLD_VALUE':transformed_bhy_thres
        #             })

# end of script