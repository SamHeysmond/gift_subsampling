# SCRIPT INFO 
# This script will write the custom R scripts used...
# ... to generate manhattan plots based on GWAS data...
# for each JOB_ID run!


import argparse

parser=argparse.ArgumentParser(description="makes R scripts for a given run")
parser.add_argument('-id', 
                    type=str, 
                    metavar='master list of all phenotype data', 
                    required=True, 
                    help='Input for the master list of phenotype data'
                    )
parser.add_argument('-i', 
                    type=str, 
                    metavar='input subsample num', 
                    required=True, 
                    help='tracks current subsample num e.g. 200 or 400'
                    )
parser.add_argument('-p', 
                    type=str, 
                    metavar='input phenotype', 
                    required=True, 
                    help='tracks the current phenotype e.g. leaf_ionome_Mo98'
                    )
parser.add_argument('-o', 
                    type=str, 
                    metavar='output directory for script', 
                    required=True, 
                    help='where the specific R script from this code will be sent to'
                    )

# stores input data and parses them
args= parser.parse_args() 

Rscript_output=open(str(args.o)+str(args.id)+"_"+str(args.i)+"_"+str(args.p)+".R",)

Rscript_output.write(f'\n')

Rscript_output.close()