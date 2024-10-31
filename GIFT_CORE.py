# Packages

import argparse, random, os
from math import sqrt, erfc, log10
import numpy as np

import modin.pandas as pandas
import ray



def ordered_list(phenotypes_file, phenotype):
	string_based_list=[]
      
	#Import phenotype file as an array and sort by the size of the second column ()
	pheno_array=pandas.read_csv(phenotypes_file)
	pheno_array=pheno_array.sort_values(by=[phenotype]) # Note: sorted smallest first to largest last
      
	# Make a list of the sample names, ordered by value
	pheno_order=pheno_array.iloc[:,0].tolist()
      
	#? Need to filter out identical values here so I don't have to do it manually
	print("pheno_order:",flush=True)
	print(pheno_order[0:5],flush=True)
      
	# convert the integer values into a string
	for integer_value in pheno_order:
		string_based_list.append(str(integer_value))
	return string_based_list #P# Jon 2.1
	

def get_position(chromosome,position):
	position_data=str(chromosome)+":"+str(position)
	return position_data

# Arguments for the script
parser = argparse.ArgumentParser(description="Sums all DPs (depths) in a vcf.")
parser.add_argument('-v', type=str, metavar='input_vcf', required=True, help='The input vcf file.')
parser.add_argument('-f', type=str, metavar='phenotypes_file', required=True, help='The input phenotype file. It is a .csv file with a one line headder and two columns: individual and phenotype. Phenotype is a numeric value and they must be in order of size.')
parser.add_argument('-p', type=str, metavar='phenotype', required=True, help='The phenotype - exactly as written in the phenotype file headder.')
parser.add_argument('-o', type=str, metavar='output_file', required=True, help='The output file.')
# ### start of sam edit (4) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# add in an ID tag to allow for proper file naming of outputs and tracking of run IDs
parser.add_argument('-id', type=str, metavar='run_ID', required=True, help='Run ID of the batch file.')
parser.add_argument('-s', type=str, metavar='subsample_num', required=True, help='Subsample number used in this test for threshold calculation.')
parser.add_argument('-m', type=str, metavar='memory', required=True, help='Memory (in GB to limit modin (pandas) with.)')
# ### end of sam edit (4) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
args = parser.parse_args()

# ray.init(_plasma_directory="/tmp", object_store_memory=700000000000) # setting to disable out of core in Ray and obj store mem increase
	# obj store memory currently at args.m defined
# memory_for_ray = int(args.m)*int(10000000000)
# Example:
	# 30000000000 = 30GB
multiplier= 10**9
memory_for_ray = int(args.m)*int(multiplier)
print("Memory for ray (BYTES) is... ",flush=True)
print(memory_for_ray,flush=True)

print("Memory for ray (GB) is... ",flush=True)
print(int(memory_for_ray)/multiplier,flush=True)

os.environ["MODIN_CPUS"] = "3"

ray.init(_plasma_directory="/tmp", object_store_memory=memory_for_ray,num_cpus=3) # setting to disable out of core in Ray and obj store mem increase
		
#test print
print("genotype_tracker_df initialising: ",  flush=True)

num_header = 0
with open(args.v) as file:
    for line in file.readlines():
        if line.startswith("##"):
            num_header += 1
        else:
            break

# read in the vcf file
vcf_df=pandas.read_csv(args.v,sep='\t',dtype='string',skiprows=num_header)
vcf_df=vcf_df.rename({'#CHROM':'CHROM'},axis=1)
vcf_df.drop(['ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'],axis=1,inplace=True)

# vcf_df.insert(loc=2, column='new col', )
vcf_df.insert(0,"position_data",value=['NA' for i in range(vcf_df.shape[0])])
vcf_df['position_data']=vcf_df.apply(lambda x: get_position(x.CHROM,x.POS),axis =1)

vcf_df.drop(['CHROM','POS',],axis=1,inplace=True)

	# dont use regex=true as it fails to work
vcf_df.replace({"0|0":"-1","1|1":"1"}, inplace=True)

# print("Vcf dataframe: ",flush=True)
# print(vcf_df.iloc[:,0:5],flush=True)

print("====================================",flush=True)
print("====================================",flush=True)

print("Transposing dataframe",flush=True)
genotype_df=vcf_df.transpose()



print("====================================",flush=True)
print("====================================",flush=True)
# save some memory (i hope)
print("Deleting old vcf dataframe",flush=True)
del vcf_df

print("====================================",flush=True)
print("====================================",flush=True)

# print("transposed df: ",flush=True)
# print(genotype_df.iloc[:,0:5],flush=True)

print("====================================",flush=True)
print("====================================",flush=True)


# print("Current columns are: ",flush=True)
col_list=genotype_df.columns.tolist()
# print(col_list[0:5],flush=True)

print("====================================",flush=True)
print("====================================",flush=True)

print("redoing column headers ",flush=True)
#set first row (positions) as header
genotype_df.columns = genotype_df.iloc[0]
genotype_df = genotype_df[1:]


# print(genotype_df.iloc[:,0:5],flush=True)
print("====================================",flush=True)
print("====================================",flush=True)

# print("new columns are: ",flush=True)
col_list=genotype_df.columns.tolist()
# print(col_list[0:5],flush=True)

# print("====================================",flush=True)
# print("====================================",flush=True)

# 	# this seemingly does nothing
# print("renaming first column ",flush=True)
# try:
# 	genotype_df.rename(columns={"position_data":"ACCESSION_ID"},inplace=True)
# except:
#       print("Could not change column name ",flush=True)
# print(genotype_df.iloc[:,0:5],flush=True)

# print("====================================",flush=True)
# print("====================================",flush=True)

ordered_phenotype_index=ordered_list(args.f, args.p)


print("Fetching ordered index to implement",flush=True)
print("====================================",flush=True)

# print(ordered_phenotype_index[0:5],flush=True)

print(f"Ordered index length: {len(ordered_phenotype_index)}",flush=True)
print("====================================",flush=True)

# print("resetting index ",flush=True)
# genotype_df = genotype_df.reset_index()  # make sure indexes pair with number of rows
# print(genotype_df.iloc[:,0:5],flush=True)

# reorder rows
print("reindexing the file",flush=True)
genotype_df=genotype_df.reindex(ordered_phenotype_index)

# print(genotype_df,flush=True)

print("====================================",flush=True)
print("====================================",flush=True)

print("Exporting the genotype dataframe",flush=True)
print("Path is: /gpfs01/home/mbysh17/core_files/genotype_tracker/"+args.id+"_genotypes.csv",flush=True)
genotype_df.to_csv("/gpfs01/home/mbysh17/core_files/genotype_tracker/genotypes_"+args.id+".csv",header=True,index=False)

# apply genotype function to all of the subsample columns e.g. 0|0 becomes -1 and 1|1 becomes 1

# combine CHROM and POS into 1 column for position data

# Then do either method A or B

# METHOD A
	# transpose the dataframe so columbs become rows and rows become columns
    # re-order the rows (accessions) so that they match the order of phenotype (with smallest at top)
	# output the df

# METHOD B
    # append the posidion data column as a list of headers for the output df
	# Following the order of phenotypes-take the column and add it as a row to an output dataframe (genotype_tracker_df)
		# do this for each accession

# output the genotype tracker to a csv file
# genotype_tracker_df.to_csv("core_files/genotype_tracker/"+args.id+"_genotypes.csv",header=True,index=False)

########################################################
# end of file