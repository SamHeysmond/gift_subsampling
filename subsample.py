# Required files:
# 1) input vcf of all genotypes (used in bcftools command towards the end)
# 2) input phenotype file (.csv format)
# Output files:
# 1) A subsampled vcf e.g. core_files/subsampled_500.vcf 
# 2) A sumsampled phenotype file (.csv format) e.g. core_files/subsampled_phenotype_500.csv
# 3) Subsampled list of subsample ids (.txt format) e.g. core_files/subsamples/subsamples_500.txt'

#packages
import argparse
import random
import os
import pandas

parser=argparse.ArgumentParser(description="subsamples a given number of individuals from master_list.csv")

parser.add_argument('-p', 
                    type=str, 
                    metavar='master list of all phenotype data', 
                    required=True, 
                    help='Input for the master list of phenotype data'
                    )
parser.add_argument('-s', 
                    type=str, 
                    metavar='list of all samples from the vcf in .txt format', 
                    required=True, 
                    help='list of all samples from the vcf in .txt format'
                    )
parser.add_argument('-n', 
                    type=str, 
                    metavar='subsample_number', 
                    required=True, 
                    help='How many samples you want to subsample from the master_list.csv'
                    )
parser.add_argument('-t', 
                    type=str, 
                    metavar='trait', 
                    required=True, 
                    help='Which trait do we focus on from the csv file'
                    )
parser.add_argument('-op', 
                    type=str, 
                    metavar='output phenotype', 
                    required=True, 
                    help='Output csv of subsampled individuals phenotype (only 1 phenotype). Format: Accession_ID,Na23/Mo98'
                    )
parser.add_argument('-og', 
                    type=str, 
                    metavar='output genotype', 
                    required=True, 
                    help='Output of vcf only containing the subsampled IDs'
                    )
parser.add_argument('-ri', 
                    type=str, 
                    metavar='run ID', 
                    required=True, 
                    help='For batch useage. Inputs the ID to make the file unique'
                    )

# stores input data and parses them
args= parser.parse_args() 

#set input file
sample_list_from_vcf = pandas.read_csv(args.s, sep=',', usecols=[0], names=['Accession_ID'], header=None)

print("sample_list_from_vcf : ",flush=True)
print(sample_list_from_vcf,flush=True)

# open phenotype data file
input_file_phenotype_data = open(args.p,'r')
sample_list_from_phenotype_file = pandas.read_csv(args.p, sep=',', usecols=['Accession_ID',args.t])
print("sample_list_from_phenotype_file : ",flush=True)
print(sample_list_from_phenotype_file,flush=True)

# merge data where Accession_ID is in both the vcf and phenotype files
consensus_ID_dataframe = (sample_list_from_vcf.reset_index(drop=True)[["Accession_ID"]].merge(sample_list_from_phenotype_file.reset_index(drop=True), on=['Accession_ID'], how='inner',left_index=False,right_index=False))

### TESTING
    # outer will bring in everything- if it overlaps itll be "both", otherwise itll be "left only" or "right only"
    # we want to keep the "left only" data
merged_df=pandas.merge(sample_list_from_phenotype_file,sample_list_from_vcf,how="outer",on=['Accession_ID'],indicator=True)

print("- - - - - - - - - - - - - - - - - - ",flush= True)
print("Merged dataframe initially:",flush= True)
print(merged_df,flush= True)
print("- - - - - - - - - - - - - - - - - - ",flush= True)

    # Delete rows where data is in both 
merged_df = merged_df.drop(merged_df[merged_df['_merge'] == 'left_only'].index)

print("- - - - - - - - - - - - - - - - - - ",flush= True)
print("Merged dataframe after removing left only",flush= True)
print(merged_df,flush= True)
print("- - - - - - - - - - - - - - - - - - ",flush= True)

    # Delete rows where data is right only
merged_df = merged_df.drop(merged_df[merged_df['_merge'] == 'right_only'].index)

print("- - - - - - - - - - - - - - - - - - ",flush= True)
print("Merged dataframe after removing right only",flush= True)
print(merged_df,flush= True)
print("- - - - - - - - - - - - - - - - - - ",flush= True)

    # remove merge column; no longer needed
final_df=merged_df.drop(columns=['_merge'])
final_df=final_df.dropna()

print("- - - - - - - - - - - - - - - - - - ",flush= True)
print("Merged dataframe after removing merge column",flush= True)
print(merged_df,flush= True)
print("- - - - - - - - - - - - - - - - - - ",flush= True)

####



print("consensus_ID_dataframe_with_phenotypes : ",flush=True)
print(consensus_ID_dataframe,flush=True)

print(f"Length of consensus dataframe (should be 1000 or 1003 ish): {len(consensus_ID_dataframe)}", flush=True)

# now randomly select a list of ID's for the subsample number e.g. 200 or 400 etc
# subsampled_dataframe = consensus_dataframe.sample(n=int(args.n))
# consensus_dataframe = consensus_dataframe.reset_index()  # make sure indexes pair with number of rows
subsampled_dataframe = consensus_ID_dataframe.sample(n=int(args.n))

print("Subsampled_dataframe before reset index: ",flush=True)
print(subsampled_dataframe,flush=True)

subsampled_dataframe.reset_index(drop=True, inplace=True)
print("Subsampled_dataframe after reset index: ",flush=True)
print(subsampled_dataframe,flush=True)

# open file for writing all the IDs that were subsampled for this particular run
subsampled_IDs_path = 'core_files/subsample_text_files/subsamples_'+str(args.n)+'_'+str(args.ri)+'.txt'

print("subsampled_dataframe data types before sorting: ",flush=True)
print(subsampled_dataframe.dtypes,flush=True)

# subsampled_dataframe = subsampled_dataframe.drop(columns='index')
subsampled_dataframe[[f'{args.t}']] = subsampled_dataframe[[f'{args.t}']].apply(pandas.to_numeric)
subsampled_dataframe = subsampled_dataframe.sort_values(by=args.t,ascending=True)

print("subsampled_dataframe after sorting: ",flush=True)
print(subsampled_dataframe,flush=True)

print("subsampled_dataframe data types after sorting: ",flush=True)
print(subsampled_dataframe.dtypes,flush=True)


# OUTPUT 1 
    # ordered subsample list of just the IDs 
subsampled_dataframe.to_csv(subsampled_IDs_path,columns=['Accession_ID'],index=False, header=False)

# OUTPUT 2
    # subsampled vcf
#now use bcftools to subsample from the vcf to make a subsampled VCF of only what is in the picked bin
os.system('bcftools view --samples-file core_files/subsample_text_files/subsamples_'+str(args.n)+'_'+str(args.ri)+'.txt core_files/FINAL.vcf > '+args.og)
    # output subsampled vcf (output:1) 

# OUTPUT 3
    # subsample of the accessions and phenotype data together
subsampled_dataframe.to_csv(args.op,header=True,index=False)

# testing print
print("End of subsample script reached",flush=True)
# end of file