# This file tracks:
#           1) Hand selected SNPs depending on the phenotype
#           2) Top 20 significant SNPs in both GWAS and GIFT for each test




# loop through the .csv from head to tail.

# if a SNP pval is lower than in the current list (which will start with all 1's) then replace it with next lowest SNP (1 in this case)

# repeat until lowest p-values are in the list

# grep out the lines that have those values in them and write them to a SNP table "top_20_SNP_<SAMPLE_NUM>_ID_ID_NUM_HERE>"

# grep out the lines containing hand picked SNP and put those in a SNP table "Focused_SNP_<SAMPLE_NUM>_ID_ID_NUM_HERE>"

# Pass these data onto R script which will plot them? (Need to alter the R script in the physics GWAS GIFT program to do this too....)

# ===================

# looping through all .csv made for given ID and given method (gift/gwas)

#packages
import argparse

parser=argparse.ArgumentParser(description="subsamples a given number of individuals from master_list.csv")

parser.add_argument('-jl', 
                    type=str, 
                    metavar='master list of all phenotype data', 
                    required=True, 
                    help='Input for the master list of phenotype data'
                    )
parser.add_argument('-d', # this is not the output where this program WITES to, its what it TAKES from
                    type=str, 
                    metavar='directory for output folder', 
                    required=True, 
                    help='Lets the program read the csv files made from GWAS and GIFT in output'
                    )
parser.add_argument('-o', 
                    type=str, 
                    metavar='output directory of result from SNP_tracker.py', 
                    required=True, 
                    help='The directory you want the csv and figures to be put in'
                    )

# stores input data and parses them
args= parser.parse_args() 

# opening necessary files
input_jobs_list=open(args.jl,'r')

output_SNP_tracker_csv=open(args.o+"SNP_tracker_tally.csv","w")

# for graphs (WIP)
#output_SNP_tracker_Rscript=.....

# set up list of phenotypes sorted

stored_phenotypes=[] #will contain things like leaf_ionome_Mo98....leaf_ionome_Rad50 etc etc
cumulative_run_number=0 # out of 100 per run (or n_runs)
current_phenotype="" #string of the current phenotype being tracked e.g. leaf_ionome_Mo98
cumulative_GWAS_significance=0 # out of 100 (or N) GWAS parallel runs how many show a SNP as significant
cumulative_GIFT_significance=0 # out of 100 (or N) GIFT parallel runs how many show a SNP as significant

N_GWAS_tests=0 #may combine these since they SHOULD be the same ALWAYS
N_GIFT_tests=0

## REMINDER OF JOB LIST FORMAT
# JOB_ID,SUBSAMPLE_N,PHENOTYPE
# 521856,200,leaf_ionome_Mo98
# 521857,200,leaf_ionome_Mo98
# …… ^^ the above lines show two RUNS, each with 200 SAMPLES from the vcf
# 521900,400,leaf_ionome_Mo98
#  …..
# 60000,200,leaf_ionome_Rad50

for line in input_jobs_list:
    # split on the comma to make it into a list
    clean_line=line.split(",")

    # read the ID of the job
    current_job_id=clean_line[0]

    # check if the PHENOTYPE is different from whats in the list (should be yes if just started too)
    if clean_line[2] not in stored_phenotypes:

        # if we arent at the first sample(/run?)  THIS FEELS WRONG -NEEDS CHANGING
        if cumulative_run_number!=0:

            #write in the current information to the output file since we're 
            # moving onto another phenotype now
            output_SNP_tracker_csv.write(str(current_phenotype)+
                                        str(current_SNP_chromosome)+
                                        str(current_SNP_position)+
                                        str(cumulative_run_number)+
                                        str(cumulative_GWAS_significance)+
                                        str(N_GWAS_tests)+
                                        str(cumulative_GIFT_significance)+
                                        str(N_GIFT_tests)
                                        )
        else: 
            # add the new phenotype to the list
            stored_phenotypes.append(str(clean_line[2]))


    # find the T20 (absolute theta) GIFT snps for that ID 
    # REMINDER OF FORMAT 
    # CHROM,POS,PVAL
    # 3,123413,0.00001
    # 4,1233141,0.0003
    # 3,5435133,0.0312
    input_T20_absolute_theta=open(args.d+current_job_id+"_T20_absolute_theta.csv","r")

    # loop through T20 absolute theta SNPs for 1/100 Job IDs (for 200 samples) -> it would be 1/run_num
    for abs_theta_line in input_T20_pSNP4:
        clean_abs_theta_line=abs_theta_line.split(",")

        # store the chromosome and position 
        current_SNP_chromosome=clean_abs_theta_line[0]
        current_SNP_position=clean_abs_theta_line[1]
        

    # find the T20 (pSNP4) GIFT snps for that ID 
    input_T20_pSNP4=open(args.d+current_job_id+"_T20_pSNP4.csv","r")

    # find the T20 (absolute theta) GIFT snps for that ID 
    input_T20_pSNP5=open(args.d+current_job_id+"_T20_pSNP5.csv","r")

    # find the HAND PICKED GIFT snps for that ID (PENDING UPDATE)
    # input_T20_HAND_PICKED=open(args.d+current_job_id+"_T20_HAND_PICKED.csv","r")

    # find the T20 GWAS snps for that ID

    # find the HAND PICKED GWAS snps for that ID (PENDING UPDATE)

    # run the numbers

    # write to the csv

# output csv format example

# Phenotype,CHR,POS,Subsample_N,N_Times_Significant_GWAS,N_GWAS_tests,N_Times_Significant_GIFT,N_GIFT_tests
# leaf_ionome_Mo98,2,123,200,85,100,99,100


# output figure ideas
# idea 1)
# Graph: bar chart
# X axis = SNP position (GWAS/GIFT share same position for each position with two adjacent bars)
# Y axis = % of times deemed significant

# idea 2) 
# graph: box and whisker plot
# x Axis = GWAS_average_T20(or custom SNP)_detection_rate VS GIFTaverage_T20(or custom SNP)_detection_rate 
# Y Axis = detection rate (%)
# Purpose: allows to view spread of detection rate of the T20 (or custom) SNPs to see if GIFT detects better

