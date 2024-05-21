## Script description
# This script will subsample
# Required files:
# 1) input vcf of all genotypes
# 2) input phenotype file (.csv format)
# Output files:
# 1) A subsampled vcf e.g. core_files/subsampled_500.vcf 
# 2) A sumsampled phenotype file (.csv format) e.g. core_files/subsampled_phenotype_500.csv
# 3) Subsampled list of subsample ids (.txt format) e.g. core_files/subsamples_500.txt'

#packages
import argparse
import random
import os

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

#create empty sampels list
samples_list_vcf =[]
samples_list_pheno =[]
consensus_sample_list=[]
picked_bin=[] 

# Main steps:
# 1) read all sample ID's from samples.txt file (sourced from main vcf) into a list

#set input file
input_file_samples_vcf=open(args.s,'r')
#read each line and append to the list
for line in input_file_samples_vcf:
    #line = read_next_line(input_file_samples_vcf)
    cleanline = line.replace('\n','')
    samples_list_vcf.append(cleanline)
print("Samples list vcf: ", samples_list_vcf)

# 2) read all sample ID's from the phenotype spreadsheet to make a list (REPEAT CODE?)
input_file_phenotype_data = open(args.p,'r')
for line in input_file_phenotype_data:
    #line = read_next_line(input_file_phenotype_data)
    cleanline=line.split(',')
    #only add the first column of data since this is the sample ID
    samples_list_pheno.append(cleanline[0])
# remove first item as this is just the title
# print("Samples list pheno: ", samples_list_pheno)
samples_list_pheno.pop(0)
# print("Samples list pheno without first item: ", samples_list_pheno)
input_file_samples_vcf.close()
input_file_phenotype_data.close()

# 3) make a consensus list of all samples that exist in both the vcf list and the excel list
# using the vcf list as a basis
for sample1 in samples_list_vcf:
    if sample1 in samples_list_pheno:
        consensus_sample_list.append(sample1)
    #else:
        #print("Sample from vcf not found in phenotype file : ", sample1)


# based on the number of subsampling needed, try to make a subsampled list of that number
# from the consensus list of samples that exist in both files
# goes from 0 to num of samples-1 e.g. 0 to 499 (for n = 500)

for n in range(int(args.n)):
    #print("Attempt (n):",n)
    new_number_picked = False

    while new_number_picked == False:
        #start off assuming it is true that the next ID has not yet been picked
        new_number_picked = True

        #pick random number from 0 to length of consensus samples list
        random_id_index = random.randint(0,len(consensus_sample_list)-1)
        '''
        #check to see if this sample has been picked previously
        if consensus_sample_list[random_id_index] in picked_bin:
            #switches to false if it finds any instance of this ID already picked
            new_number_picked = False
            consensus_sample_list.pop(random_id_index)
        '''
        # if the value contains an "NA", remove it and dont pick it
        if "NA" in str(consensus_sample_list[random_id_index]).upper():
            consensus_sample_list.pop(random_id_index)

            # set the flag to false so we loop again
            new_number_picked = False

    #add the sample to the picked bin
    picked_bin.append(consensus_sample_list[random_id_index])   

    #remove the value you just picked from the consensus sample list so you cant pick it again
    consensus_sample_list.pop(random_id_index)

#finish by printing set of IDs picked 
#print("Picked bin:", picked_bin)
subsampled_IDs=open('core_files/subsample_text_files/subsamples_'+str(args.n)+'_'+str(args.ri)+'.txt','w')
#write the subsampled data into a new text file
for ID in picked_bin:
    subsampled_IDs.write(str(ID)+'\n')
subsampled_IDs.close()

#now use bcftools to subsample from the vcf to make a subsampled VCF of only what is in the picked bin
os.system('bcftools view --samples-file core_files/subsample_text_files/subsamples_'+str(args.n)+'_'+str(args.ri)+'.txt core_files/1001genomes_snp_biallelic_only_ACGTN.vcf > '+args.og)

#Sift through the spreadsheet file to only include the phenotype and samples from the picked_bin
input_file_phenotype_data = open(args.p,'r')
output_file_subsampled_phenotype = open(args.op,"w")
for index, line in enumerate(input_file_phenotype_data):
    cleanline=line.split(',')
    #for the header
    if index == 0:
        #print(cleanline)
        #find index of the phenotype header 
        Phenotype_index = cleanline.index(args.t)
        #writes only the correct headers to the new csv file
        output_file_subsampled_phenotype.write(str(cleanline[0])+','+str(cleanline[Phenotype_index])+'\n') 
    else:
        #check if the ID matches to one of the picked IDs (only writing it in if so)
        if cleanline[0] in picked_bin:
            output_file_subsampled_phenotype.write(str(cleanline[0])+','+str(cleanline[Phenotype_index])+'\n')

# closing phenotype data inputs and outputs
input_file_phenotype_data.close()
output_file_subsampled_phenotype.close()

# testing print
print("End of subsample script reached")
# end of file