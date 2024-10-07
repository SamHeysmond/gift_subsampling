#packages
import argparse
import pandas

parser=argparse.ArgumentParser(description="")

parser.add_argument('-i', 
                    type=str, 
                    metavar='GWAS association results', 
                    required=True, 
                    help='give an input from the GWAS pipeline output (ends in .assoc.clean.txt)'
                    )
parser.add_argument('-o', 
                    type=str, 
                    metavar='', 
                    required=True, 
                    help='CSV output for future code usage'
                    )


# stores input data and parses them
args= parser.parse_args() 


#set input file to a dataframe (it should be able to read the tsv disguised as a txt)
# use with the pipeline
# GWAS_data=pandas.read_csv(args.i, sep='\t', usecols=['CHR','BP','P'])

# use from own code
GWAS_data=pandas.read_csv(args.i, sep='\t', usecols=['chr','ps','p_lrt'])

print("Dataframe before reformat",flush=True)
GWAS_data.head()

# rename the columns to fit the format I need (Saves having to do it in later python script)
GWAS_data.rename(columns={
                "chr":"CHROM",
                "ps":"POS",
                "p_lrt":"pvals"
                },inplace=True)

# use this for now until i reverse the need for it later
# GWAS_data.rename(columns={
#                 "CHR":"chromosomes",
#                 "BP":"positions",
#                 "P":"pvals"
#                 },inplace=True)

print("Dataframe after reformat",flush=True)
GWAS_data.head()

#output format
#CHROM, POS, PVAL
GWAS_data.to_csv(args.o,sep=',', index=False)

# testing print
print("GWAS data successfully reformatted to csv",flush=True)
# end of file