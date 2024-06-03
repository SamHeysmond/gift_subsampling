# File version: 03_06_2024

# This file:
#   1) Tracks Hand selected SNPs depending on the phenotype  (X)
#   2) Tracks Top 20 significant SNPs in both GWAS and GIFT for each test (X)
#   3) Stores all SNPs for GWAS and GIFT, averaging them and making a csv for R to use (yes)




# loop through the .csv from head to tail.

# if a SNP pval is lower than in the current list (which will start with all 1's) then replace it with next lowest SNP (1 in this case)

# repeat until lowest p-values are in the list

# grep out the lines that have those values in them and write them to a SNP table "top_20_SNP_<SAMPLE_NUM>_ID_ID_NUM_HERE>"

# grep out the lines containing hand picked SNP and put those in a SNP table "Focused_SNP_<SAMPLE_NUM>_ID_ID_NUM_HERE>"

# Pass these data onto R script which will plot them? (Need to alter the R script in the physics GWAS GIFT program to do this too....)

# ===================

# looping through all .csv made for given ID and given method (gift/gwas)

#packages
import argparse, pandas, os, math

parser=argparse.ArgumentParser(description="subsamples a given number of individuals from master_list.csv")

parser.add_argument('-jl', 
                    type=str, 
                    metavar='master list of all phenotype data', 
                    required=True, 
                    help='Input for the master list of phenotype data'
                    )
parser.add_argument('-d', # this is NOT the output where this program WITES to, its what it TAKES from
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

#output_SNP_tracker_csv=open(args.o+"SNP_tracker_tally.csv","w")


# Phenotype,CHR,POS,Subsample_N,N_Times_Significant_GWAS,N_GWAS_tests,N_Sig_GIFT_Absolute_theta,N_Sig_GIFT_pSNP4,N_Sig_GIFT_pSNP5,N_GIFT_tests
# leaf_ionome_Mo98,2,123,200,85,100,99,79,80,100
# might remove the N_SIGS_GIFT_ABSOLUTE THETA UNTIL I GET A THRESHOLD TO USE!
SNP_dataframe = pandas.DataFrame(columns=['PHENOTYPE',
                                          'CHR',
                                          'POS',
                                          'SUBSAMPLE_NUM',
                                          'N_SIGS_GWAS',
                                          'N_GWAS_TESTS',
                                          'N_SIGS_GIFT_ABSOLUTE_THETA',
                                          'N_SIGS_GIFT_PSNP4',
                                          'N_SIGS_GIFT_PSNP5',
                                          'N_GIFT_TESTS'])

# this is across all subsample numbers
General_T20_SNPs_dataframe = pandas.DataFrame(columns=['PHENOTYPE',
                                          'CHR',
                                          'POS',
                                          'N_SIGS_GWAS',
                                          'N_GWAS_TESTS',
                                          'N_SIGS_GIFT_ABSOLUTE_THETA',
                                          'N_SIGS_GIFT_PSNP4',
                                          'N_SIGS_GIFT_PSNP5',
                                          'N_GIFT_TESTS'])

# variables for positive and negative control
positive_control_upper_bound = 0
positive_control_chromosome = 0
negative_control_upper_bound = 0
negative_control_chromosome = 0


# for graphs (WIP)
#output_SNP_tracker_Rscript=.....

# set up list of phenotypes sorted
stored_phenotypes=[] #will contain things like leaf_ionome_Mo98....leaf_ionome_Na23 etc etc
cumulative_run_number=0 # out of 100 per run (or n_runs)
current_phenotype="" #string of the current phenotype being tracked e.g. leaf_ionome_Mo98
cumulative_GWAS_significance=0 # out of 100 (or N) GWAS parallel runs how many show a SNP as significant
cumulative_GIFT_significance=0 # out of 100 (or N) GIFT parallel runs how many show a SNP as significant
N_GWAS_tests=0 #may combine these since they SHOULD be the same ALWAYS
N_GIFT_tests=0

# set up positive and negative control dataframes 
# each SNP is uniquely identified by CHR, POS and SUBSAMPLE NUMBER
positive_control_Mo98_df = pandas.DataFrame(columns=[
                                                    'CHR', 
                                                    'POS',
                                                    'SUBSAMPLE_NUM',
                                                    'N_SIGS_GWAS',
                                                    'N_GWAS_TESTS',
                                                    'N_SIGS_GIFT_ABSOLUTE_THETA',
                                                    'N_SIGS_GIFT_PSNP4',
                                                    'N_SIGS_GIFT_PSNP5',
                                                    'N_GIFT_TESTS'
                                                    ])

# dataframe for ALL SNPs for IDEA 3 (FOR Mo98)
Mo98_ALL_SNPS_GIFT_df = pandas.DataFrame(columns=[
                                                    'CHR', 
                                                    'POS',
                                                    'SUBSAMPLE_NUM',  # for each subsamp number 200-1000
                                                    'TOTAL_PSNP4',
                                                    'TOTAL_PSNP5',
                                                    'TOTAL_ABS_THETA',
                                                    'TIMES_APPEARED', # COULD be 90 or 99 or 100 or 10 who knows
                                                    'TOTAL_GIFT'   #SHOULD be 100
                                                    ])

Mo98_ALL_SNPS_GWAS_df = pandas.DataFrame(columns=[
                                                    'CHR', 
                                                    'POS',
                                                    'SUBSAMPLE_NUM',  # for each subsamp number 200-1000
                                                    'TOTAL_P',
                                                    'TIMES_APPEARED', # COULD be 90 or 99 or 100 or 10 who knows
                                                    'TOTAL_GWAS'   #SHOULD be 100
                                                    ])

Na23_ALL_SNPS_GIFT_df = pandas.DataFrame(columns=[
                                                    'CHR', 
                                                    'POS',
                                                    'SUBSAMPLE_NUM',  # for each subsamp number 200-1000
                                                    'TOTAL_PSNP4',
                                                    'TOTAL_PSNP5',
                                                    'TOTAL_ABS_THETA',
                                                    'TIMES_APPEARED', # COULD be 90 or 99 or 100 or 10 who knows
                                                    'TOTAL_GIFT'   #SHOULD be 100
                                                    ])

Na23_ALL_SNPS_GWAS_df = pandas.DataFrame(columns=[
                                                    'CHR', 
                                                    'POS',
                                                    'SUBSAMPLE_NUM',  # for each subsamp number 200-1000
                                                    'TOTAL_P',
                                                    'TIMES_APPEARED', # COULD be 90 or 99 or 100 or 10 who knows
                                                    'TOTAL_GWAS'   #SHOULD be 100
                                                    ])



# further dataframes can then look at these SNPs and compare their n_sig values between GIFT and GWAS at DIFFERENT subsample nums

## REMINDER OF JOB LIST FORMAT
# JOB_ID,SUBSAMPLE_N,PHENOTYPE
# 521856,200,leaf_ionome_Mo98
# 521857,200,leaf_ionome_Mo98
# …… ^^ the above lines show two RUNS, each with 200 SAMPLES from the vcf
# 521900,400,leaf_ionome_Mo98
#  …..
# 60000,200,leaf_ionome_Rad50

# EACH JOB ID IS CONNECTED TO THE FOLLOWING THINGS
# GWAS T20 SNPS
# GIFT ABSOLUTE THETA T20 SNPS
# GIFT PSNP4 T20 SNPS
# GIFT PSNP5 T20 SNPS
# HANDPICKED SNPS? (WIP?)
# looping through jobs list file


# GWAS csv name = output_files/${phenotype}_GWAS_${i}_${SLURM_JOB_ID}.csv
# GIFT csv name = output_files/${phenotype}_whole_genome_metrics_${i}_${SLURM_JOB_ID}.csv
# phenotype name examples:
# leaf_ionome_Mo98
# leaf_ionome_Na23

# Reminder of CSV format (GIFT) NAME leaf_ionome_Mo98_whole_genome_metrics_600_732692.csv
# CHROM,POS,largest_theta,smallest_theta,absolute_theta,theta_range,largest_relative_theta,smallest_relative_theta,absolute_relative_theta,range_relative_theta,min_p,mean_p,log_mean_p,bigest_theta_p,pSNP4,pSNP5
# 1,73,6.285714285714285,-4.224489795918359,6.285714285714285,10.510204081632644,0.3845193508114856,-0.25842696629213435,0.3845193508114856,0.6429463171036199,4.185363300872138e-05,nan,nan,0.00015134587316541535,0.15259768149369662,0.6580333017260325

# from the above i want: CHROM,POS,absolute_theta,pSNP4,pSNP5

# reminder of csv format (GWAS) NAME leaf_ionome_Mo98_GWAS_600_732692.csv
# chromosomes,positions,pvals,mafs,macs,genotype_var_perc
# 1,55,0.621946624343,0.0016694490818,1,0.000407516956329

def GATHER_ALL_GIFT_SNPS(current_phenotype_df, ALL_GIFT_SNPS_df,Total_GIFT,current_subsample_num):

    # increment the GIFT total counter
    Total_GIFT+=1

    # loop through all rows of the current csv file being analysed
    for current_index, current_row in current_phenotype_df.iterrows():

        #set the assumption that the current SNP in the current csv file is NOT already in the total SNP dataframe
        in_all_SNPS=False

        # loop through the total SNP dataframe
        for ALL_SNPS_index, ALL_SNPS_row in ALL_GIFT_SNPS_df.iterrows():
                
            # check if the current SNP for the current subsample matches anything in the ALL_SNPS dataframe for the following parameters.
            if int(current_row['POS']) == int(ALL_SNPS_row['POS']) and int(current_row['CHR']) == int(ALL_SNPS_row['CHR']) and int(current_subsample_num) == int(ALL_SNPS_row['SUBSAMPLE_NUM']):
                
                # set flag to true
                in_all_SNPS=True

                # update the row that its on (no need to update CHR, POS or SUBSAMPLE_NUM -> just the rest)
                ALL_GIFT_SNPS_df.at[ALL_SNPS_index,'TOTAL_PSNP4']= float(ALL_SNPS_row['TOTAL_PSNP4']) + float(current_row['pSNP4'])
                ALL_GIFT_SNPS_df.at[ALL_SNPS_index,'TOTAL_PSNP5']= float(ALL_SNPS_row['TOTAL_PSNP5']) + float(current_row['pSNP5'])
                ALL_GIFT_SNPS_df.at[ALL_SNPS_index,'TOTAL_ABS_THETA']= float(ALL_SNPS_row['TOTAL_ABS_THETA']) + float(current_row['absolute_theta'])
                ALL_GIFT_SNPS_df.at[ALL_SNPS_index,'TIMES_APPEARED']= int(ALL_SNPS_row['TIMES_APPEARED']) + 1
                ALL_GIFT_SNPS_df.at[ALL_SNPS_index,'TOTAL_GIFT']= int(ALL_SNPS_row['TOTAL_GIFT']) + 1

        # if it doesnt find this SNP in the main dataframe, then add this SNP to the main dataframe for all SNPs
        if in_all_SNPS==False:
            
                # add new data to dataframe
                # add in a line for the current line of data in the calculation (NOT SURE IF NEED TO CONVERT TO INTEGER OR NOT HERE)
            new_row=pandas.Series({
                                    "CHR":current_row['CHR'],
                                    "POS":current_row['POS'],
                                    "SUBSAMPLE_NUM":int(current_subsample_num),
                                    "TOTAL_PSNP4":current_row['pSNP4'],
                                    "TOTAL_PSNP5":current_row['pSNP5'],
                                    "TOTAL_ABS_THETA":current_row['absolute_theta'],
                                    "TIMES_APPEARED":"1",
                                    "TOTAL_GIFT":Total_GIFT
                                    })
            
            # join this new row (for this SNP) to the main dataframe
            ALL_GIFT_SNPS_df =pandas.concat([ALL_GIFT_SNPS_df, new_row.to_frame().T], ignore_index=True)

            # re-sort the df so CHR and POS go in ascending order -> it looks nicer that way.
            ALL_GIFT_SNPS_df.sort_values(by=["CHR","POS"], axis=0, ascending=[True,True],inplace=True, na_position='first')

    print("GIFT GATHER function done")
    return ALL_GIFT_SNPS_df, Total_GIFT

def GATHER_ALL_GWAS_SNPS(current_phenotype_df,ALL_GWAS_SNPS_df,Total_GWAS,current_subsample_num):
    
    # increment the GWAS total counter
    Total_GWAS+=1

    # loop through all rows of the current csv file being analysed
    for current_index, current_row in current_phenotype_df.iterrows():

        #set the assumption that the current SNP in the current csv file is NOT already in the total SNP dataframe
        in_all_SNPS=False

        # loop through total SNP dataframe
        for ALL_SNPS_index, ALL_SNPS_row in ALL_GWAS_SNPS_df.iterrows():
                
            #check if the current SNP for the current subsample matches anything in the ALL_SNPS dataframe
            if int(current_row['positions']) == int(ALL_SNPS_row['POS']) and int(current_row['chromosomes']) == int(ALL_SNPS_row['CHR']) and int(current_subsample_num) == int(ALL_SNPS_row['SUBSAMPLE_NUM']):
                
                # set flag to true
                in_all_SNPS=True

                # update the row that its on (no need to update CHR, POS or SUBSAMPLE_NUM -> just the rest)
                ALL_GWAS_SNPS_df.at[ALL_SNPS_index,'TOTAL_P']= float(ALL_SNPS_row['TOTAL_P']) + float(current_row['pvals'])
                ALL_GWAS_SNPS_df.at[ALL_SNPS_index,'TIMES_APPEARED']= int(ALL_SNPS_row['TIMES_APPEARED']) + 1
                ALL_GWAS_SNPS_df.at[ALL_SNPS_index,'TOTAL_GWAS']= int(ALL_SNPS_row['TOTAL_GWAS']) + 1

        # if it doesnt find this SNP in the main dataframe, then add it to the main dataframe for all SNPs
        if in_all_SNPS==False:
                # add new data to dataframe
                # add in a line for the current line of data in the calculation (# NOT SURE IF NEED TO CONVERT TO INTEGER OR NOT HERE)
            new_row=pandas.Series({
                                    "CHR":current_row['chromosomes'],
                                    "POS":current_row['positions'],
                                    "SUBSAMPLE_NUM":int(current_subsample_num),
                                    "TOTAL_P":current_row['pvals'],
                                    "TIMES_APPEARED":"1", # not sure if string or int here
                                    "TOTAL_GWAS":Total_GWAS
                                    })
            
            # add the new row (current SNP) to the main SNP dataframe
            ALL_GWAS_SNPS_df =pandas.concat([ALL_GWAS_SNPS_df, new_row.to_frame().T], ignore_index=True)

            # re-sort the df so CHR and POS go in ascending order -> it looks nicer that way.
            ALL_GWAS_SNPS_df.sort_values(by=["CHR","POS"], axis=0, ascending=[True,True],inplace=True, na_position='first')
    print("GWAS GATHER function finished")
    return ALL_GWAS_SNPS_df, Total_GWAS

def CALCULATE_AVERAGE_SNPS_GIFT(ALL_GIFT_SNPS_df):
    ALL_GIFT_SNPS_df["AVERAGE_PSNP4"] = ALL_GIFT_SNPS_df["TOTAL_PSNP4"] / ALL_GIFT_SNPS_df["TIMES_APPEARED"]
    ALL_GIFT_SNPS_df["AVERAGE_PSNP5"] = ALL_GIFT_SNPS_df["TOTAL_PSNP5"] / ALL_GIFT_SNPS_df["TIMES_APPEARED"]
    ALL_GIFT_SNPS_df["AVERAGE_ABS_THETA"] =ALL_GIFT_SNPS_df["TOTAL_ABS_THETA"] / ALL_GIFT_SNPS_df["TIMES_APPEARED"]

    return ALL_GIFT_SNPS_df

def CALCULATE_AVERAGE_SNPS_GWAS(ALL_GWAS_SNPS_df):
    ALL_GWAS_SNPS_df["AVERAGE_P"] = ALL_GWAS_SNPS_df["TOTAL_P"] / ALL_GWAS_SNPS_df["TIMES_APPEARED"]

    return ALL_GWAS_SNPS_df

def MAKE_R_SCRIPT(phenotype,subsample_number,pval_type):
    # make R script for each P value type (pSNP4, pSNP5, abs theta) -> check phys gwas for things like abs theta incase different
    R_out=open("output_files/"+str(phenotype)+"_"+str(subsample_number)+"_"+str(pval_type)+"_AVG_MANHATTAN.R","w")
    R_out.write(f'#R script for making manhattan plots with ggplot\n')
    R_out.write(f'library("tidyverse")\n')
    R_out.write(f'library("ggrepel")\n')   
    R_out.write(f'\n')

    # for GWAS data
    if pval_type=="TOTAL_P":
        R_out.write(f'GWAS_ALL_SNPS_DATA<-read.csv("output_files/{phenotype}_ALL_SNPS_GWAS.csv",header=true)\n')  
        R_out.write(f'GWAS_SUBSAMPLE_{subsample_number}_SNPS_DATA <- subset(GWAS_ALL_SNPS_DATA,SUBSAMPLE_NUM=={subsample_number})\n') 
    
    # for GIFT data
    else:
        R_out.write(f'GIFT_ALL_SNPS_DATA<-read.csv("output_files/{phenotype}_ALL_SNPS_GWAS.csv",header=true)\n')  
        R_out.write(f'{pval_type}_SUBSAMPLE_{subsample_number}_SNPS_DATA <- subset(GWAS_ALL_SNPS_DATA,SUBSAMPLE_NUM=={subsample_number})\n') 
    
    # cumulative calculations
    R_out.write(f'don <- {pval_type}_SUBSAMPLE_{subsample_number}_SNPS_DATA %>%\n')
    R_out.write(f'     # Compute CHR size\n')
    R_out.write(f'     group_by(CHR) %>% \n')
    R_out.write(f'     summarise(chr_len=max(POS)) %>%\n')
    R_out.write(f'     # Calculate cumulative position of each CHR\n')
    R_out.write(f'     mutate(tot=cumsum(chr_len)-chr_len) %>%\n')
    R_out.write(f'     select(-chr_len) %>%\n')
    R_out.write(f'     # Add this info to the initial dataset\n')
    R_out.write(f'     left_join({pval_type}_SUBSAMPLE_{subsample_number}_SNPS_DATA, ., by=c("CHR"="CHR")) %>%\n')
    R_out.write(f'     # Add a cumulative position of each SNP\n')
    R_out.write(f'     arrange(CHR, POS) %>%\n')
    R_out.write(f'     mutate( BPcum=POS+tot) \n')
    R_out.write(f'axisdf = don %>%\n')
    R_out.write(f'     group_by(CHR) %>%\n')
    R_out.write(f'     summarize(center=( max(BPcum) + min(BPcum) ) / 2 )\n')

    # set y limit for the graph (not sure if this changes anything here though- just copied from other script)
    if pval_type=="TOTAL_ABS_THETA": #test
        R_out.write(f'ylim <- abs(floor(min({pval_type}_SUBSAMPLE_{subsample_number}_SNPS_DATA${pval_type}))) +1\n')
    else:
        R_out.write(f'ylim <- abs(floor(log10(min({pval_type}_SUBSAMPLE_{subsample_number}_SNPS_DATA${pval_type})))) +1\n')

    R_out.write(f'#open png\n')
    R_out.write(f'png("output_files/{phenotype}_{subsample_number}_{pval_type}"_AVG_MANHATTAN.png", bg = "white", width = 9.75, height = 3.25, units = "in", res = 1200, pointsize = 4)\n')

    # make the plot
    if pval_type=="TOTAL_ABS_THETA":
        R_out.write(f'ggplot(don, aes(x=BPcum, y=({pval_type}), color=as_factor(CHR))) +\n')
    else:
        R_out.write(f'ggplot(don, aes(x=BPcum, y=-log10({pval_type}), color=as_factor(CHR))) +\n')
    R_out.write(f'     # Show all points\n')
    R_out.write(f'     geom_point(alpha=0.5) +\n')
    R_out.write(f'     # custom X axis:\n')
    R_out.write(f'     scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +\n')
    # TEMP PAUSE FOR TESTING
    #Rscript_output.write(f'     scale_y_continuous(expand = c(0, 0) ) + # remove space between plot area and x axis\n')
    # label axis
    R_out.write(f'     # my axis labels\n')
    if pval_type=="TOTAL_ABS_THETA":
        R_out.write(f'     labs(y= "({pval_type})", x = "chromosome position")+\n')
    else:
        R_out.write(f'     labs(y= "-log10({pval_type})", x = "chromosome position")+\n')
    # add a theme
    R_out.write(f'     # Custom the theme:\n')
    R_out.write(f'     theme_minimal() +\n')
    R_out.write(f'     guides(colour="none")\n')
    R_out.write(f'     theme(\n')
    R_out.write(f'       panel.border = element_blank(),\n')
    R_out.write(f'       panel.grid.major.x = element_blank(),\n')
    R_out.write(f'       panel.grid.minor.x = element_blank(),\n')
    R_out.write(f'       axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)\n')
    R_out.write(f'      )\n')
    R_out.write(f'dev.off()\n')
    R_out.write(f'#END OF SCRIPT\n')
    R_out.close()

    os.system("output_files/"+str(phenotype)+"_"+str(subsample_number)+"_"+str(pval_type)+"_AVG_MANHATTAN.R")
    return

Total_GIFT_Mo98 = 0
Total_GWAS_Mo98 = 0

Total_GIFT_Na23 = 0
Total_GWAS_Na23 = 0
csv_files = os.listdir(args.d)

for csv_file in csv_files:

    # split on the underscore to make lists
    csv_file=csv_file.split("_")

    #check for the specific phenotype based on the naming convention
    # AND split code into GIFT and GWAS specific code

    if csv_files[2] == "Mo98" and csv_files[3]=='whole':# GIFT code vvvvv
        
        #Total_GIFT+=1 (in function instead)
        Current_Mo98_dataframe=pandas.read_csv(csv_file)

        Mo98_ALL_SNPS_GIFT_df,Total_GIFT_Mo98=GATHER_ALL_GIFT_SNPS(Current_Mo98_dataframe,Mo98_ALL_SNPS_GIFT_df,Total_GIFT_Mo98,int(csv_files[6]))

        ''' THRESHOLD CODE WIP
        # calculate thresholds (1) and (2) for each p value (pSNP4 and pSNP5 (and abs theta?))
        # threshold (1.1) pSNP4 BYH threshold
        m = len(Current_Mo98_dataframe)
        Current_Mo98_dataframe.sort_values(by=["pSNP4"], axis=0, ascending=True,inplace=True, na_position='first')
        Current_Mo98_dataframe.reset_index()  # make sure indexes pair with number of rows
        s = 1
        i = 0
        for current_index, current_row in Current_Mo98_dataframe.iterrows():
            p = float(current_row["pSNP4"])
            i += 1
            if i>1:
                    s = s + 1/(i-1)
            thres_pval = ((i+1)/m) * 0.05 / s
            if (p>thres_pval):
                break
        BHY_pSNP4_thres = -(math.log(thres_pval, 10))
        
        # threshold (1.2) pSNP4 BT threshold
        BT_pSNP4 = 0.05 / (len(Current_Mo98_dataframe) * int(csv_files[6])) # multiply by subsample num which is in file name
        BF_THRES_pSNP4 = -(math.log(BT_pSNP4,10))

         # threshold (2.1) pSNP5 BYH threshold
        m = len(Current_Mo98_dataframe)
        Current_Mo98_dataframe.sort_values(by=["pSNP5"], axis=0, ascending=True,inplace=True, na_position='first')
        Current_Mo98_dataframe.reset_index()  # make sure indexes pair with number of rows
        s = 1
        i = 0
        for current_index, current_row in Current_Mo98_dataframe.iterrows():
            p = float(current_row["pSNP5"])
            i += 1
            if i>1:
                    s = s + 1/(i-1)
            thres_pval = ((i+1)/m) * 0.05 / s
            if (p>thres_pval):
                break
        BHY_pSNP5_thres = -(math.log(thres_pval, 10))
        
        # threshold (2.2) pSNP5 BT threshold
        BT_pSNP5 = 0.05 / (len(Current_Mo98_dataframe) * int(csv_files[6])) # multiply by subsample num which is in file name
        BF_THRES_pSNP5 = -(math.log(BT_pSNP5,10))

         # threshold (3.1) absolute theta threshold?
         # awaiting confirmation on how to threshold this value
         # place holder bonferroni correction
        BT_abstheta = 0.05 / (len(Current_Mo98_dataframe) * int(csv_files[6])) # multiply by subsample num which is in file name
        BF_THRES_abstheta = BT_abstheta

        ########
        ## EDIT - Will need to calculate distance from peak to the mot1 gene and add this to the boundary!!
        #######################
        ## code for that probably goes here

        ###################
        # now cycle through the csv and find positive control SNPs in the MOT1 gene
        # MOT1 gene location boundaries 
        positive_control_chromosome = 2
        positive_control_lower_bound_position = 10933005
        positive_control_upper_bound_position = 10934604

        for current_index, current_row in Current_Mo98_dataframe.iterrows():
             THRESHOLD CODE WIP
            if positive_control_lower_bound_position<=int(current_row["POS"])<=positive_control_upper_bound_position:
                # then check if its significant or not for each type of pval
                if float(current_row["pSNP4"])>= BHY_pSNP4_thres:
                    #print("This SNP is significant for pSNP4 under BHY")
                    # check if this SNP with its CHR POS AND SUBSAMPLE_NUM is already in the dataframe
                    already_there = False
                    for inner_index, inner_row in positive_control_Mo98_df.iterrows():
                        if int(current_row['POS']) == int(inner_row['POS']) and int(current_row['CHR']) == int(inner_row['CHR']) and int(csv_files[6]) == int(inner_row['SUBSAMPLE_NUM']):
                                already_there = True
                if float(current_row["pSNP4"])>= BF_THRES_pSNP4:
                    #print("This SNP is significant for pSNP4 under BF")
                    pass
        '''
            
    elif csv_files[2]=="Mo98" and csv_files[3]=='GWAS':# GWAS code vvv      
        #Total_GWAS+=1
        Current_Mo98_dataframe=pandas.read_csv(csv_file)
        Mo98_ALL_SNPS_GWAS_df,Total_GWAS_Mo98 = GATHER_ALL_GWAS_SNPS(Current_Mo98_dataframe,Mo98_ALL_SNPS_GWAS_df,Total_GWAS_Mo98, int(csv_files[6]))
        
        ''' THRESHOLD CODE WIP
        # calculate thresholds (1) and (2) for GWAS p values

        # threshold (1.1) pSNP4 BYH threshold
        m = len(Current_Mo98_dataframe)
        Current_Mo98_dataframe.sort_values(by=["pvals"], axis=0, ascending=True,inplace=True, na_position='first')
        Current_Mo98_dataframe.reset_index()  # make sure indexes pair with number of rows
        s = 1
        i = 0
        for current_index, current_row in Current_Mo98_dataframe.iterrows():
            p = float(current_row["pval"])
            i += 1
            if i>1:
                    s = s + 1/(i-1)
            thres_pval = ((i+1)/m) * 0.05 / s
            if (p>thres_pval):
                break
        BHY_p_thres = -(math.log(thres_pval, 10))
        
        # threshold (1.2) pSNP4 BT threshold
        BT_p = 0.05 / (len(Current_Mo98_dataframe) * int(csv_files[6])) # multiply by subsample num which is in file name
        BF_THRES_p = -(math.log(BT_p,10))

        ########
        ## EDIT - Will need to calculate distance from peak to the mot1 gene and add this to the boundary!!
        #######################
        ## code for that probably goes here

        ###################
        # now cycle through the csv and find positive control SNPs in the MOT1 gene
        # MOT1 gene location boundaries 

        THRESHOLD CODE WIP
        positive_control_chromosome = 2
        positive_control_lower_bound_position = 10933005
        positive_control_upper_bound_position = 10934604

        for current_index, current_row in Current_Mo98_dataframe.iterrows():

            if positive_control_lower_bound_position<=int(current_row["POS"])<=positive_control_upper_bound_position:
                # then check if its significant or not for each type of pval
                if float(current_row["pSNP4"])>= BHY_pSNP4_thres:
                    print("This SNP is significant for pSNP4 under BHY")
                    # check if this SNP with its CHR POS AND SUBSAMPLE_NUM is already in the dataframe
                    already_there = False
                    for inner_index, inner_row in positive_control_Mo98_df.iterrows():
                        if int(current_row['POS']) == int(inner_row['POS']) and int(current_row['CHR']) == int(inner_row['CHR']) and int(csv_files[6]) == int(inner_row['SUBSAMPLE_NUM']):
                                already_there = True
                if float(current_row["pSNP4"])>= BF_THRES_pSNP4:
                    print("This SNP is significant for pSNP4 under BF")

        '''

    elif csv_files[2] == "Na23" and csv_files[3]=='whole':
        Current_Na23_dataframe=pandas.read_csv(csv_file) 
        Na23_ALL_SNPS_GIFT_df,Total_GIFT_Na23=GATHER_ALL_GIFT_SNPS(Current_Na23_dataframe,Na23_ALL_SNPS_GIFT_df,Total_GIFT_Na23,int(csv_files[6]))

    elif csv_files[2] == "Na23" and csv_files[3]=='GWAS':
        #Total_GWAS+=1
        Current_Na23_dataframe=pandas.read_csv(csv_file)
        Na23_ALL_SNPS_GWAS_df,Total_GWAS_Na23 = GATHER_ALL_GWAS_SNPS(Current_Na23_dataframe,Na23_ALL_SNPS_GWAS_df,Total_GWAS_Na23, int(csv_files[6]))

    else:
        print("failure to detect phenotype and GWAS/GIFT method")

# Now go through each of the total SNP dataframes and calculate the average P values for each SNP
Mo98_ALL_SNPS_GIFT_df = CALCULATE_AVERAGE_SNPS_GIFT(Mo98_ALL_SNPS_GIFT_df)
Mo98_ALL_SNPS_GWAS_df = CALCULATE_AVERAGE_SNPS_GWAS(Mo98_ALL_SNPS_GWAS_df)

Na23_ALL_SNPS_GIFT_df = CALCULATE_AVERAGE_SNPS_GIFT(Na23_ALL_SNPS_GIFT_df)
Na23_ALL_SNPS_GWAS_df = CALCULATE_AVERAGE_SNPS_GWAS(Na23_ALL_SNPS_GWAS_df)

# write results to csv files for an R script to run them
Mo98_ALL_SNPS_GIFT_df.to_csv("output_files/Mo98_ALL_SNPS_GIFT.csv", index=False)
Mo98_ALL_SNPS_GWAS_df.to_csv("output_files/Mo98_ALL_SNPS_GWAS.csv", index=False)

Na23_ALL_SNPS_GIFT_df.to_csv("output_files/Na23_ALL_SNPS_GIFT.csv", index=False)
Na23_ALL_SNPS_GWAS_df.to_csv("output_files/Na23_ALL_SNPS_GWAS.csv", index=False)

subsample_num_list=[200,400,600,800,1000]
phenotype_list=["Mo98","Na23"]
pvals=["TOTAL_P","TOTAL_PSNP4","TOTAL_PSNP5","TOTAL_ABS_THETA"]
# R script creation and running
for phenotype in phenotype_list:

    pass

    for subsample_number in subsample_num_list:

        for pval_type in pvals:

            MAKE_R_SCRIPT(phenotype,subsample_number,pval_type)
          
            
            
            pass












'''
for line in input_jobs_list:
    # split on the comma to make it into a list
    clean_line=line.split(",")

    # read the ID of the job
    current_job_id=clean_line[0]

    # check if the PHENOTYPE is different from whats in the list (should be yes if just started too)
    if clean_line[2] not in stored_phenotypes:
        # current phenotype not seen yet so either new phenotype or begining
        # now simply write in the GWAS T20 SNP list first
        input_GWAS_T20=open(args.d+clean_line[2]+"_GWAS_T20_SNPS_"+clean_line[1]+"_"+clean_line[0]+".csv","r")
        for GWAS_T20_SNP in input_GWAS_T20:


            new_row=pandas.Series({"CHROM":CHROM,"POS":POS,"PVAL":absolute_theta})
            SNP_dataframe =pandas.concat([SNP_dataframe, new_row.to_frame().T], ignore_index=True)

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

# Phenotype,CHR,POS,Subsample_N,N_Times_Significant_GWAS,N_GWAS_tests,N_Sig_GIFT_Absolute_theta,N_Sig_GIFT_pSNP4,N_Sig_GIFT_pSNP5,N_GIFT_tests
# leaf_ionome_Mo98,2,123,200,85,100,99,79,80,100


# the following graphs are for:
# Graph: bar chart
# OF ONE PHENOTYPE
# OF ONE SAMPLE NUMBER
# EXAMPLE Pheno = leaf_ionome_Mo96, Subsample_N = 200
# Can also include hand picked Snps

# output figure ideas
# idea 1)

# X axis = SNP position (GWAS/GIFT share same chr + position with two adjacent bars)
# Y axis = % of times deemed significant

# idea 2) 
# graph: box and whisker plot
# x Axis = GWAS_average_T20(or custom SNP)_detection_rate VS GIFTaverage_T20(or custom SNP)_detection_rate 
# Y Axis = detection rate (%)
# Purpose: allows to view spread of detection rate of the T20 (or custom) SNPs to see if GIFT detects better

'''
