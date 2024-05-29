## Script description
# This script will subsample
# Required files:
# 1) input vcf of all genotypes
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
                    help='Output csv of subsampled individuals phenotype (only 1 phenotype)'
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

#create empty samples list for storing the VCF sample IDs
samples_list_vcf =[]

# Main steps:
# 1) read all sample ID's from samples.txt file (sourced from main vcf) into a list

#set input file
input_file_samples_vcf=open(args.s,'r')
#read each line and append to the list
for line in input_file_samples_vcf:
    #line = read_next_line(input_file_samples_vcf)
    cleanline = line.replace('\n','')
    samples_list_vcf.append(str(cleanline))
#print("Samples list vcf: ", samples_list_vcf)

# open phenotype data file
input_file_phenotype_data = open(args.p,'r')

row=0
for line in input_file_phenotype_data:
    line = line.replace('\n','')
    line=line.split(',')
    if row == 0:
        # find index of the header for the trait/phenotype we are looking into
        phenotype_index = line.index(args.t)
        #create dataframe with the correct two headers (ID,phenotype/trait data)
        consensus_dataframe = pandas.DataFrame(columns=[line[0], line[phenotype_index]])

    else: #otherwise look to add data to the dataframe
        # if it has NA in the specific phenotype column
        if "NA" in str(line[phenotype_index]).upper():
            print("NA detected, ommitting this line")

        else: #add new row to the consensus dataframe
            new_row=pandas.Series({"1001_Genomes_ID":line[0],str(args.t):line[phenotype_index]})
            consensus_dataframe=pandas.concat([consensus_dataframe, new_row.to_frame().T], ignore_index=True)

    row +=1

input_file_phenotype_data.close()

print("Consensus_dataframe BEFORE checking for consensus")
print(consensus_dataframe)

# 3) make a consensus list of all samples that exist in both the vcf list and the excel list
# using the excel list as a basis (in the form of a pandas dataframe)
consensus_dataframe = consensus_dataframe.reset_index()  # make sure indexes pair with number of rows

for index, row in consensus_dataframe.iterrows():

    # check if each of the rows of the consensus dataframe are in the list of IDs from the vcf
    # if not -> remove the row from the dataframe 
    if str(row["1001_Genomes_ID"]) in samples_list_vcf:
        # do nothing since this ID is in both VCF and phenotype data files so can stay!
        pass 

    else: #remove the row
        print("Sample from phenotype file not found in vcf, removing ID : ", row["1001_Genomes_ID"])
        consensus_dataframe = consensus_dataframe.drop([index])

print("Consensus_dataframe AFTER checking for consensus")
print(consensus_dataframe)

# now randomly select a list of ID's for the subsample number e.g. 200 or 400 etc
subsampled_dataframe = consensus_dataframe.sample(n=int(args.n))
consensus_dataframe = consensus_dataframe.reset_index()  # make sure indexes pair with number of rows
print("Subsampled dataframe")
print(subsampled_dataframe)

subsampled_IDs=open('core_files/subsample_text_files/subsamples_'+str(args.n)+'_'+str(args.ri)+'.txt','w')

# open the output phenotype file
output_file_subsampled_phenotype = open(args.op,"w")

#consensus_dataframe = consensus_dataframe.reset_index()  # make sure indexes pair with number of rows
for index, row in subsampled_dataframe.iterrows():

    #write in each ID from the subsampled dataframe
    subsampled_IDs.write(str(row["1001_Genomes_ID"])+'\n')

    #write in each ID and phenotype information from the subsampled dataframe to a subsampled phenotype file
    output_file_subsampled_phenotype.write(str(row["1001_Genomes_ID"])+','+str(row[str(args.t)])+'\n') 

# close files vvv

# Output subsample ID text list (output:2) finished being made!
subsampled_IDs.close()

# Output phenotype csv (output:3) finished being made!
output_file_subsampled_phenotype.close()

#now use bcftools to subsample from the vcf to make a subsampled VCF of only what is in the picked bin
os.system('bcftools view --samples-file core_files/subsample_text_files/subsamples_'+str(args.n)+'_'+str(args.ri)+'.txt core_files/1001genomes_snp_biallelic_only_ACGTN.vcf > '+args.og)
# output subsampled vcf (output:1) finished being made! 

# testing print
print("End of subsample script reached")
# end of file