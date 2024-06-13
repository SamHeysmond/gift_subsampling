# TESTING OF PARALLELISATION

#packages
import pandas, os, math
import time
import concurrent.futures

# start time
t1=time.perf_counter()

# so i can see all the columns when testing with print
pandas.set_option('display.max_columns',None)
pandas.options.display.max_columns=None

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

# set up positive and negative control dataframes (might combine later for efficiency and space)
# each SNP is uniquely identified by CHR, POS and SUBSAMPLE NUMBER
positive_control_Mo98_df = pandas.DataFrame(columns=[
                                                    'CHR', 
                                                    'POS',
                                                    'PVAL_TYPE',
                                                    'SUBSAMPLE_NUM',
                                                    'VALUE'
                                                    ])

negative_control_Mo98_df = pandas.DataFrame(columns=[
                                                    'CHR', 
                                                    'POS',
                                                    'PVAL_TYPE',
                                                    'SUBSAMPLE_NUM',
                                                    'VALUE'
                                                    ])

positive_control_Na23_df = pandas.DataFrame(columns=[
                                                    'CHR', 
                                                    'POS',
                                                    'PVAL_TYPE',
                                                    'SUBSAMPLE_NUM',
                                                    'VALUE'
                                                    ])

negative_control_Na23_df = pandas.DataFrame(columns=[
                                                    'CHR', 
                                                    'POS',
                                                    'PVAL_TYPE',
                                                    'SUBSAMPLE_NUM',
                                                    'VALUE'
                                                    ])

#IDEA 3 dataframe for ALL SNPs for (for Mo98 GIFT)
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

#IDEA 3 dataframe for ALL SNPs for (for Mo98 GWAS)
Mo98_ALL_SNPS_GWAS_df = pandas.DataFrame(columns=[
                                                    'CHR', 
                                                    'POS',
                                                    'SUBSAMPLE_NUM',  # for each subsamp number 200-1000
                                                    'TOTAL_P',
                                                    'TIMES_APPEARED', # COULD be 90 or 99 or 100 or 10 who knows
                                                    'TOTAL_GWAS'   #SHOULD be 100
                                                    ])

#IDEA 3 dataframe for ALL SNPs for (for Na23 GIFT)
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

#IDEA 3 dataframe for ALL SNPs for (for Na23 GWAS)
Na23_ALL_SNPS_GWAS_df = pandas.DataFrame(columns=[
                                                    'CHR', 
                                                    'POS',
                                                    'SUBSAMPLE_NUM',  # for each subsamp number 200-1000
                                                    'TOTAL_P',
                                                    'TIMES_APPEARED', # COULD be 90 or 99 or 100 or 10 who knows
                                                    'TOTAL_GWAS'   #SHOULD be 100
                                                    ])

# mark for removal?
T20_SNPS_GIFT_AND_GWAS = pandas.DataFrame(columns=[
                                                    'CHR', 
                                                    'POS',
                                                    'SUBSAMPLE_NUM',  # for each subsamp number 200-1000
                                                    'TOTAL_PSNP4',
                                                    'TOTAL_PSNP5',
                                                    'TOTAL_ABS_THETA',
                                                    'TIMES_APPEARED', # COULD be 90 or 99 or 100 or 10 who knows
                                                    'TOTAL_GIFT'   #SHOULD be 100
                                                    ])

# Reminder of CSV format (GIFT) NAME leaf_ionome_Mo98_whole_genome_metrics_600_732692.csv
# CHROM,POS,largest_theta,smallest_theta,absolute_theta,theta_range,largest_relative_theta,smallest_relative_theta,absolute_relative_theta,range_relative_theta,min_p,mean_p,log_mean_p,bigest_theta_p,pSNP4,pSNP5
# 1,73,6.285714285714285,-4.224489795918359,6.285714285714285,10.510204081632644,0.3845193508114856,-0.25842696629213435,0.3845193508114856,0.6429463171036199,4.185363300872138e-05,nan,nan,0.00015134587316541535,0.15259768149369662,0.6580333017260325

# from the above i want: CHROM,POS,absolute_theta,pSNP4,pSNP5

# reminder of csv format (GWAS) NAME leaf_ionome_Mo98_GWAS_600_732692.csv
# chromosomes,positions,pvals,mafs,macs,genotype_var_perc
# 1,55,0.621946624343,0.0016694490818,1,0.000407516956329

# IDEA 3.1
# combined in the sense that it stacks/combines snps from same location and works for GWAS and GIFT
# TRACE 2
def IDEA_3_GATHER_ALL_SNPS_COMBINED(GWAS_or_GIFT,current_phenotype_df, GWAS_OR_GIFT_ALL_SNPS_df,Total_GIFT_or_GWAS,current_subsample_num): 
    print("Entered FUNCTION: IDEA_3_GATHER_ALL_SNPS_COMBINED",flush=True)

    if GWAS_or_GIFT == "GWAS":
        CHR = "chromosomes"
        POS = "positions"
        # increment the GWAS total counter
        Total_GIFT_or_GWAS+=1

        #############################
        #### NEW CODE 2
        ##########################

        # convert format of current csv to new format 
        temp_dataframe = current_phenotype_df.copy()

        #   changing the names of current columns
        temp_dataframe.rename(columns={CHR:'CHR',POS:'POS','pvals':'TOTAL_P'}, inplace=True)

        # keep specific columns needed
        temp_dataframe=temp_dataframe[['CHR','POS','TOTAL_P']]

        # adding in two new columns at specific index values
        temp_dataframe["TIMES_APPEARED"] = 1
        temp_dataframe["TOTAL_GWAS"] = Total_GIFT_or_GWAS
        temp_dataframe.insert(2,"SUBSAMPLE_NUM",current_subsample_num)
        
  
        # concatonate to the cumulative dataframe
        #print("Temp GWAS dataframe looks like: ",flush=True)
        #print(temp_dataframe.head(),flush=True)

        #GWAS_OR_GIFT_ALL_SNPS_df =pandas.concat([GWAS_OR_GIFT_ALL_SNPS_df,temp_dataframe], ignore_index=True)
        # testing this
        GWAS_OR_GIFT_ALL_SNPS_df=pandas.concat([GWAS_OR_GIFT_ALL_SNPS_df, temp_dataframe],ignore_index=True)

        # group by and sort the data
        GWAS_OR_GIFT_ALL_SNPS_df = GWAS_OR_GIFT_ALL_SNPS_df.groupby(['CHR','POS','SUBSAMPLE_NUM'],as_index=False).agg({'TOTAL_P':'sum','TIMES_APPEARED':'sum','TOTAL_GWAS':'max'})
        #############################
        ##########################

    elif GWAS_or_GIFT == "GIFT":
        CHR="CHROM"
        POS ="POS"

        # increment the GIFT total counter
        Total_GIFT_or_GWAS+=1

        #############################
        #### NEW CODE 2
        ##########################
        # convert format of current csv to new format 
        temp_dataframe = current_phenotype_df.copy()
        #   changing the names of current columns
        temp_dataframe.rename(columns={CHR:'CHR',POS:'POS','pSNP4':'TOTAL_PSNP4','pSNP5':'TOTAL_PSNP5','absolute_theta':'TOTAL_ABS_THETA'}, inplace=True)
        # keep specific columns needed
        temp_dataframe=temp_dataframe[['CHR','POS','TOTAL_PSNP4','TOTAL_PSNP5','TOTAL_ABS_THETA']]
        
        # adding in two new columns at specific index values
        temp_dataframe.insert(2,"SUBSAMPLE_NUM",current_subsample_num)
        temp_dataframe["TIMES_APPEARED"] = 1
        temp_dataframe["TOTAL_GIFT"] = Total_GIFT_or_GWAS
        # concatonate to the cumulative dataframe
        
        #GWAS_OR_GIFT_ALL_SNPS_df =pandas.concat([temp_dataframe,GWAS_OR_GIFT_ALL_SNPS_df], ignore_index=True)
        # testing the below code line
        GWAS_OR_GIFT_ALL_SNPS_df=pandas.concat([GWAS_OR_GIFT_ALL_SNPS_df,temp_dataframe],ignore_index=True)
        
        # group by and sort the data
        GWAS_OR_GIFT_ALL_SNPS_df = GWAS_OR_GIFT_ALL_SNPS_df.groupby(['CHR','POS','SUBSAMPLE_NUM'],as_index=False).agg({'TOTAL_PSNP4':'sum','TOTAL_PSNP5':'sum','TOTAL_ABS_THETA':'sum','TIMES_APPEARED':'sum','TOTAL_GIFT':'max'})
        #############################
        ##########################

    print("IDEA 3 GIFT GATHER function done",flush=True)

    return GWAS_OR_GIFT_ALL_SNPS_df, Total_GIFT_or_GWAS
    
    # end of function code

# IDEA 3
def IDEA_3_CALCULATE_AVERAGE_SNPS_GIFT(ALL_GIFT_SNPS_df):
    print("Entered FUNCTION: IDEA_3_CALCULATE_AVERAGE_SNPS_GIFT",flush=True)

    ALL_GIFT_SNPS_df["AVERAGE_PSNP4"] = ALL_GIFT_SNPS_df["TOTAL_PSNP4"] / ALL_GIFT_SNPS_df["TIMES_APPEARED"]
    ALL_GIFT_SNPS_df["AVERAGE_PSNP5"] = ALL_GIFT_SNPS_df["TOTAL_PSNP5"] / ALL_GIFT_SNPS_df["TIMES_APPEARED"]
    ALL_GIFT_SNPS_df["AVERAGE_ABS_THETA"] =ALL_GIFT_SNPS_df["TOTAL_ABS_THETA"] / ALL_GIFT_SNPS_df["TIMES_APPEARED"]

    return ALL_GIFT_SNPS_df

# IDEA 3
def IDEA_3_CALCULATE_AVERAGE_SNPS_GWAS(ALL_GWAS_SNPS_df):
    print("Entered FUNCTION: IDEA_3_CALCULATE_AVERAGE_SNPS_GWAS",flush=True)
    ALL_GWAS_SNPS_df["AVERAGE_P"] = ALL_GWAS_SNPS_df["TOTAL_P"] / ALL_GWAS_SNPS_df["TIMES_APPEARED"]

    return ALL_GWAS_SNPS_df

# IDEA 3
def IDEA_3_R_AND_BATCH(phenotype,subsample_number,pval_type):
    print("Entered FUNCTION: IDEA_3_R_AND_BATCH",flush=True)
    #########################################
    ### MAKE R SCRIPT ##################
    #############################

    # make R script for each P value type (pSNP4, pSNP5, abs theta, GWAS_P)
    R_out=open("/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Github/gift_subsampling/test_env/output_files/SNP_tracker_R_scripts/"+str(phenotype)+"_"+str(subsample_number)+"_"+str(pval_type)+"_AVG_MANHATTAN.R","w")
    R_out.write(f'#R script for making manhattan plots with ggplot\n')
    R_out.write(f'library("tidyverse")\n')
    #R_out.write(f'library("ggplot2")\n')
    #R_out.write(f'library("dplyr")\n')
    # R_out.write(f'library("ggrepel")\n')   # Dont need ggrepel for now
    R_out.write(f'print("Start of IDEA 3 R")\n')
    R_out.write(f'\n')

    # for GWAS data....
    if pval_type=="TOTAL_P":

        # fetch the data of the csv for the current phenotype and current method (GWAS)
        R_out.write(f'GWAS_ALL_SNPS_DATA<-read.csv("output_files/R_DATA/{phenotype}_ALL_SNPS_GWAS.csv",header=TRUE)\n')  
        
        # subsample this dataset to the current subsample number (e.g. 400)
        R_out.write(f'{pval_type}_SUBSAMPLE_{subsample_number}_SNPS_DATA <- subset(GWAS_ALL_SNPS_DATA,SUBSAMPLE_NUM=={subsample_number})\n') 
    
    # for GIFT data
    else:

        # fetch the data of the csv for the current phenotype and current method (GIFT)
        R_out.write(f'GIFT_ALL_SNPS_DATA<-read.csv("output_files/R_DATA/{phenotype}_ALL_SNPS_GIFT.csv",header=TRUE)\n')  

        # subsample this dataset to the current subsample number (e.g. 400)
        # dont worry about pval type separation ,that bit comes soon.
        R_out.write(f'{pval_type}_SUBSAMPLE_{subsample_number}_SNPS_DATA <- subset(GIFT_ALL_SNPS_DATA,SUBSAMPLE_NUM=={subsample_number})\n') 
    
    # cumulative calculations
    R_out.write(f'mydata <- {pval_type}_SUBSAMPLE_{subsample_number}_SNPS_DATA %>%\n')
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
    R_out.write(f'axisdf = mydata %>%\n')
    R_out.write(f'     group_by(CHR) %>%\n')
    R_out.write(f'     summarize(center=( max(BPcum) + min(BPcum) ) / 2 )\n')

    # set y limit for the graph (not sure if this changes anything here though)
    # y limit depends on the pval type in question
    if pval_type=="TOTAL_ABS_THETA": # testing without log10
        R_out.write(f'ylim <- abs(floor(min({pval_type}_SUBSAMPLE_{subsample_number}_SNPS_DATA${pval_type}))) +1\n')
    else:
        R_out.write(f'ylim <- abs(floor(log10(min({pval_type}_SUBSAMPLE_{subsample_number}_SNPS_DATA${pval_type})))) +1\n')

    R_out.write(f'#open png\n')
    # Send output to IDEA 3 summary plot folder
    R_out.write(f'png("output_files/summary_plots/IDEA3/{phenotype}_{subsample_number}_{pval_type}_AVG_MANHATTAN.png", bg = "white", width = 9.75, height = 3.25, units = "in", res = 1200, pointsize = 4)\n')

    # make the plot
    # this is where the plot depends on the pval type. 
    # for absolute theta, there is no -log10 since it isnt quite a "pvalue"
    if pval_type=="TOTAL_ABS_THETA":
        R_out.write(f'ggplot(mydata, aes(x=BPcum, y=({pval_type}), color=as_factor(CHR))) +\n')
    else:
        R_out.write(f'ggplot(mydata, aes(x=BPcum, y=(-log10({pval_type})), color = as_factor(CHR))) +\n')
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
    R_out.write(f'print("End of IDEA 3 R script")\n')
    R_out.write(f'#END OF SCRIPT\n')
    R_out.close()

    ###################################
    ### MAKE BATCH SCRIPT ######
    #####################

    #  make a batch file for the script and put it into batch_files/parallel_stage2/

    # location where the batch scripts will be written to
    R_batch=open("/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Github/gift_subsampling/test_env/batch_files/parallel_stage2/"+str(phenotype)+"_"+str(subsample_number)+"_"+str(pval_type)+"_AVG_MANHATTAN.sh","w")
    
    # necessary start to the file
    R_batch.write(f'#!/bin/bash\n')
    R_batch.write(f'#SBATCH --partition=defq\n')
    R_batch.write(f'#SBATCH --nodes=1\n')
    R_batch.write(f'#SBATCH --ntasks=1\n')
    R_batch.write(f'#SBATCH --cpus-per-task=3\n')
    R_batch.write(f'#SBATCH --mem=8g\n')
    R_batch.write(f'#SBATCH --time=1:00:00\n')
    R_batch.write(f'#SBATCH --job-name=R_subrun\n')
    R_batch.write(f'#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out\n')
    R_batch.write(f'#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err\n')
    R_batch.write(f'#SBATCH --mail-type=ALL\n')
    R_batch.write(f'#SBATCH --mail-user=mbysh17@nottingham.ac.uk\n')
    R_batch.write(f'#===============================\n')
    R_batch.write(f'echo "start of IDEA 3 batch script"\n')
    R_batch.write(f'#change to home directory\n')
    R_batch.write(f'cd /gpfs01/home/mbysh17\n')
    R_batch.write(f'# source conda environments\n')
    R_batch.write(f'source ~/.bashrc\n')
    R_batch.write(f'conda deactivate\n')
    R_batch.write(f'conda activate gift_env\n')
    R_batch.write(f'# R SCRIPT FOR (IDEA 3) AVG MANHATTAN PLOT\n')
    R_batch.write(f'Rscript output_files/SNP_tracker_R_scripts/{phenotype}_{subsample_number}_{pval_type}_AVG_MANHATTAN.R\n')
    R_batch.write(f'conda deactivate\n')
    R_batch.write(f'echo "End of IDEA 3 batch script"\n')
    R_batch.write(f'conda deactivate\n')
    R_batch.write(f'# END OF FILE\n')
    R_batch.close()
    # end of function

# IDEA 1.1
def GET_T20_LOCATIONS_AT_1000(all_snps_dataframe,pval_type): #parameters of: phenotype, subsample number, and pval type (implied method)
    print("Entered FUNCTION: IDEA_1_GET_T20_LOCATIONS_AT_1000",flush=True)
    # need to pass in...
    # (for each phenotype) GWAS ALL SNPS dataframe  (x2)
    # (for each phenotype) GIFT ALL SNPS dataframe  (x2)
    
    # 1) sort GWAS all SNPs by 1000 subsample and the PVAL (so biggest is at the top)
    #Mo98_ALL_SNPS_GWAS_df.sort_values(by=["SUBSAMPLE_NUM","AVERAGE_P"], axis=0, ascending=[False,True],inplace=True, na_position='first')
    all_snps_dataframe.sort_values(by=["SUBSAMPLE_NUM",pval_type], axis=0, ascending=[False,True],inplace=True, na_position='first')

    # 2) take top 20 of these averaged GWAS values
    # example format:
    # CHR    POS    SUBSAMPLE_NUM    TOTAL_P    TIMES_APPEARED    TOTAL_GWAS    AVERAGE_P
    # 1      342    1000             7          98                100           0.05
    # 4      787    1000             12         99                100           0.90
    current_pval_T20_df=all_snps_dataframe.head(20)

    # keep ONLY the location information
    current_pval_T20_df=current_pval_T20_df[['CHR','POS']]

    # return the dataframe of the top 20 snps at 1000 subsamples
    return current_pval_T20_df

# IDEA 1.3.1
def IDEA_1_ACCUMULATE_T20_SNP_DATA(
        current_dataframe,
        GWAS_P_locations_dataframe,
        PSNP4_locations_dataframe,
        PSNP5_locations_dataframe,
        ABS_THETA_locations_dataframe,
        cumulative_t20_dataframe,
        GWAS_or_GIFT,
        subsample_level
        ):
    print("Entered FUNCTION: IDEA_1_ACCUMULATE_T20_SNP_DATA",flush=True)
    # if its a GWAS csv then we're just looking at GWAS_P values to add in

    if GWAS_or_GIFT=="GWAS":
        CHR = "chromosomes"
        POS = "positions"
        current_dataframe.rename(columns={CHR:'CHR',POS:'POS'}, inplace=True)

    elif GWAS_or_GIFT == "GIFT":
        CHR="CHROM"
        POS ="POS"
        current_dataframe.rename(columns={CHR:'CHR',POS:'POS'}, inplace=True)

    ##############
    ###### NEW CODE
    ##############
    if GWAS_or_GIFT=="GWAS":
        # take the CHR and POS from the locations dataframe as a "df_key"
        df_key=GWAS_P_locations_dataframe.loc[:,["CHR","POS"]]

        #then use inner join to produce records ONLY where chr and pos match whats in the dataframe
        df_out_1 = current_dataframe.merge(
            df_key,
            how='inner',
            left_on=['CHR','POS'],
            right_on=['CHR','POS']
        )

        # keep these columns
        df_out_2= df_out_1.loc[:,["CHR","POS","pvals"]]

        # add in subsample num and pval type column data
        df_out_2["SUBSAMPLE_NUM"] = subsample_level
        df_out_2.insert(2,"PVAL_TYPE","AVERAGE_P")

        # concatonate it with the main t20 dataframe

        #cumulative_t20_dataframe =pandas.concat([df_out_2,cumulative_t20_dataframe], ignore_index=True)
        # testing below code line
        cumulative_t20_dataframe=pandas.concat([cumulative_t20_dataframe,df_out_2],ignore_index=True)

    elif GWAS_or_GIFT=="GIFT":
        ####
        ## PSNP4
        ####
        # take the CHR and POS from the locations dataframe as a "df_key"
        df_key=PSNP4_locations_dataframe.loc[:,["CHR","POS"]]

        #then use inner join to produce records ONLY where chr and pos match whats in the dataframe
        df_out_1 = current_dataframe.merge(
            df_key,
            how='inner',
            left_on=['CHR','POS'],
            right_on=['CHR','POS']
        )

        # keep these columns
        df_out_2= df_out_1.loc[:,["CHR","POS","pSNP4"]]

        #process the table to concatonate it to the main t20 cumulative dataframe
        df_out_2.rename(columns={'pSNP4':"PSNP4"})

        # melt down the dataframe for pos and neg
        df_out_2=pandas.melt(df_out_2, id_vars=['CHR','POS'],var_name="PVAL_TYPE",value_name="VALUE")

        # add in the subsample level
        df_out_2["SUBSAMPLE_NUM"] = subsample_level

        # CONCATONATE
        #cumulative_t20_dataframe =pandas.concat([cumulative_t20_dataframe,df_out_2], ignore_index=True)
        cumulative_t20_dataframe=pandas.concat([cumulative_t20_dataframe,df_out_2],ignore_index=True)

        # repeat for other two types of p val but concatonate to df_out_2 instead of equals to it
        ####
        ## PSNP5
        ####
        # take the CHR and POS from the locations dataframe as a "df_key"
        df_key=PSNP5_locations_dataframe.loc[:,["CHR","POS"]]

        #then use inner join to produce records ONLY where chr and pos match whats in the dataframe
        df_out_1 = current_dataframe.merge(
            df_key,
            how='inner',
            left_on=['CHR','POS'],
            right_on=['CHR','POS']
        )
   
        # keep these columns
        df_out_2= df_out_1.loc[:,["CHR","POS","pSNP5"]]

        #process the table to concatonate it to the main t20 cumulative dataframe
        df_out_2.rename(columns={'pSNP5':"PSNP5"})

        # melt down the dataframe for pos and neg
        df_out_2=pandas.melt(df_out_2, id_vars=['CHR','POS'],var_name="PVAL_TYPE",value_name="VALUE")

        # add in the subsample level
        df_out_2["SUBSAMPLE_NUM"] = subsample_level

        # CONCATONATE
        #cumulative_t20_dataframe =pandas.concat([cumulative_t20_dataframe,df_out_2], ignore_index=True)
        cumulative_t20_dataframe=pandas.concat([cumulative_t20_dataframe,df_out_2],ignore_index=True)
        ####
        ## absolute_theta
        ####

           # take the CHR and POS from the locations dataframe as a "df_key"
        df_key=ABS_THETA_locations_dataframe.loc[:,["CHR","POS"]]

        #then use inner join to produce records ONLY where chr and pos match whats in the dataframe
        df_out_1 = current_dataframe.merge(
            df_key,
            how='inner',
            left_on=['CHR','POS'],
            right_on=['CHR','POS']
        )
   
        # keep these columns
        df_out_2= df_out_1.loc[:,["CHR","POS","absolute_theta"]]

        #process the table to concatonate it to the main t20 cumulative dataframe
        df_out_2.rename(columns={'absolute_theta':"ABS_THETA"})

        # melt down the dataframe for pos and neg
        df_out_2=pandas.melt(df_out_2, id_vars=['CHR','POS'],var_name="PVAL_TYPE",value_name="VALUE")

        # add in the subsample level
        df_out_2["SUBSAMPLE_NUM"] = subsample_level

        # testing to see what the df_out_2 looks like with absolute theta in it
        #print("Temp t20 gwas dataframe looks like: ",flush=True)
        #print(df_out_2.head(),flush=True)

        # CONCATONATE
        #cumulative_t20_dataframe =pandas.concat([cumulative_t20_dataframe,df_out_2], ignore_index=True)
        cumulative_t20_dataframe=pandas.concat([cumulative_t20_dataframe,df_out_2],ignore_index=True)
    ####
    ####
    ####
    return cumulative_t20_dataframe

# IDEA 1.4.1
def IDEA_1_MAKE_R_SCRIPT(
        phenotype,
        cumulative_t20_dataframe,
        ):
    print("Entered FUNCTION: IDEA_1_MAKE_R_SCRIPT",flush=True)
    cumulative_t20_dataframe.to_csv("/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Github/gift_subsampling/test_env/output_files/R_DATA/"+str(phenotype)+"_T20_TRACKED.csv", index = False)

    #####################################
    ####### MAKING THE R SCRIPT
    ######################################

    # write the R script to pair with the SNP and its data
    CURRENT_SNP_R_SCRIPT=open("/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Github/gift_subsampling/test_env/output_files/SNP_tracker_R_scripts/"+str(phenotype)+"_T20_TRACKED.R","w")
    CURRENT_SNP_R_SCRIPT.write(f'#R script for making box plots with ggplot\n')
    #CURRENT_SNP_R_SCRIPT.write(f'library("ggplot2")\n')
    CURRENT_SNP_R_SCRIPT.write(f'library("tidyverse")\n')
    CURRENT_SNP_R_SCRIPT.write(f'print("start of IDEA1 R script")\n')
    CURRENT_SNP_R_SCRIPT.write(f'T20_TRACKED_DATA<- read.csv("output_files/R_DATA/{phenotype}_T20_TRACKED.csv", header= TRUE, sep=",")\n')
    CURRENT_SNP_R_SCRIPT.write(f'\n')
    CURRENT_SNP_R_SCRIPT.write(f'#open png\n')
    CURRENT_SNP_R_SCRIPT.write(f'png("output_files/summary_plots/IDEA1/{phenotype}_T20_TRACKED.png", bg = "white", width = 9.75, height = 3.25, units = "in", res = 1200, pointsize = 4)\n')
    CURRENT_SNP_R_SCRIPT.write(f'ggplot(T20_TRACKED_DATA, aes(x=PVAL_TYPE, y=VALUE, fill=SUBSAMPLE_NUM)) +\n')
    CURRENT_SNP_R_SCRIPT.write(f'   geom_boxplot()\n')
    CURRENT_SNP_R_SCRIPT.write(f'# no facet_wrap for this code\n')
    CURRENT_SNP_R_SCRIPT.write(f'\n')
    CURRENT_SNP_R_SCRIPT.write(f'dev.off()\n')
    CURRENT_SNP_R_SCRIPT.write(f'print("End of IDEA1 R script")\n')
    CURRENT_SNP_R_SCRIPT.write(f'# END OF R SCRIPT')
    CURRENT_SNP_R_SCRIPT.close()

# IDEA 1.5.1
def IDEA_1_MAKE_BASH_SCRIPT(
        phenotype,
        ):
    print("Entered FUNCTION: IDEA_1_MAKE_BASH_SCRIPT",flush=True)
    #####################################
    ####### MAKING THE BATCH SCRIPT
    ######################################

    # location where the batch scripts will be written to
    CURRENT_SNP_BATCH=open("/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Github/gift_subsampling/test_env/batch_files/parallel_stage2/"+str(phenotype)+"_T20_TRACKED.sh","w")
    
    # necessary start to the file
    CURRENT_SNP_BATCH.write(f'#!/bin/bash\n')
    CURRENT_SNP_BATCH.write(f'#SBATCH --partition=defq\n')
    CURRENT_SNP_BATCH.write(f'#SBATCH --nodes=1\n')
    CURRENT_SNP_BATCH.write(f'#SBATCH --ntasks=1\n')
    CURRENT_SNP_BATCH.write(f'#SBATCH --cpus-per-task=3\n')
    CURRENT_SNP_BATCH.write(f'#SBATCH --mem=8g\n')
    CURRENT_SNP_BATCH.write(f'#SBATCH --time=1:00:00\n')
    CURRENT_SNP_BATCH.write(f'#SBATCH --job-name=R_subrun\n')
    CURRENT_SNP_BATCH.write(f'#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out\n')
    CURRENT_SNP_BATCH.write(f'#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err\n')
    CURRENT_SNP_BATCH.write(f'#SBATCH --mail-type=ALL\n')
    CURRENT_SNP_BATCH.write(f'#SBATCH --mail-user=mbysh17@nottingham.ac.uk\n')
    CURRENT_SNP_BATCH.write(f'#===============================\n')
    CURRENT_SNP_BATCH.write(f'echo "start OF IDEA 1 batch script" \n')
    CURRENT_SNP_BATCH.write(f'#change to home directory\n')
    CURRENT_SNP_BATCH.write(f'cd /gpfs01/home/mbysh17\n')
    CURRENT_SNP_BATCH.write(f'# source conda environments\n')
    CURRENT_SNP_BATCH.write(f'source ~/.bashrc\n')
    CURRENT_SNP_BATCH.write(f'conda deactivate\n')
    CURRENT_SNP_BATCH.write(f'conda activate gift_env\n')
    CURRENT_SNP_BATCH.write(f'# R SCRIPT FOR (IDEA 1) BOXPLOT\n')
    CURRENT_SNP_BATCH.write(f'Rscript output_files/SNP_tracker_R_scripts/{phenotype}_T20_TRACKED.R"\n')
    CURRENT_SNP_BATCH.write(f'conda deactivate\n')
    CURRENT_SNP_BATCH.write(f'echo "END OF IDEA 1 batch script" \n')
    CURRENT_SNP_BATCH.write(f'# end of script')
    CURRENT_SNP_BATCH.close()

# IDEA 2.1.1
def IDEA_2_CONTROL_CHECK(current_phenotype_df,
            positive_control_df,
            positive_control_chromosome,
            positive_control_LB,
            positive_control_UB,
            negative_control_df,
            negative_control_chromosome,
            negative_control_LB,
            negative_control_UB,
            GWAS_or_GIFT,
            current_subsample_num
            ):
    print("Entered FUNCTION: IDEA_2",flush=True)
    #print("State of GWAS or GIFT: ",GWAS_or_GIFT,flush=True)
    # list used for splitting 1 line from a GIFT csv into separate lines for each pval type (in the list)
    list_of_pval_types = ["pSNP4","pSNP5","absolute_theta"]

    if GWAS_or_GIFT == "GWAS":
        CHR = "chromosomes"
        POS = "positions"
    elif GWAS_or_GIFT == "GIFT": 
        CHR="CHROM"
        POS="POS"

    ####
    ### New code
    ####

    # check if a given row in the dataframe fits the current guidelines
    temp_positive_control_df=current_phenotype_df.copy()
    temp_negative_control_df=current_phenotype_df.copy()


    temp_positive_control_df = temp_positive_control_df[(temp_positive_control_df[CHR]==positive_control_chromosome) & (temp_positive_control_df[POS]>= positive_control_LB) & (temp_positive_control_df[POS]<= positive_control_UB)]
    temp_negative_control_df = temp_negative_control_df[(temp_negative_control_df[CHR]==negative_control_chromosome) & (temp_negative_control_df[POS]>= negative_control_LB) & (temp_negative_control_df[POS]<= negative_control_UB)]

    # melt the data to change columns into separate rows depending on GIFT or GWAS method
    # only melt for GIFT, not needed for GWAS
    if GWAS_or_GIFT == "GWAS":
            
        #   changing the names of current columns
        temp_positive_control_df.rename(columns={CHR:'CHR',POS:'POS','pvals':'VALUE'}, inplace=True)
        temp_negative_control_df.rename(columns={CHR:'CHR',POS:'POS','pvals':'VALUE'}, inplace=True)

        # keep specific columns needed
        temp_positive_control_df=temp_positive_control_df[['CHR','POS','VALUE']]
        temp_negative_control_df=temp_negative_control_df[['CHR','POS','VALUE']]
        
        #add column for pval type and subsample number
        temp_positive_control_df["SUBSAMPLE_NUM"] = current_subsample_num             
        temp_positive_control_df.insert(2,"PVAL_TYPE","GWAS_P")

        temp_negative_control_df["SUBSAMPLE_NUM"] = current_subsample_num 
        temp_negative_control_df.insert(2,"PVAL_TYPE","GWAS_P")


    elif GWAS_or_GIFT == "GIFT": 

        #change column names
        temp_positive_control_df.rename(columns={CHR:'CHR',POS:'POS','absolute_theta':'ABS_THETA','pSNP4':'PSNP4','pSNP5':'PSNP5'}, inplace=True)
        temp_negative_control_df.rename(columns={CHR:'CHR',POS:'POS','absolute_theta':'ABS_THETA','pSNP4':'PSNP4','pSNP5':'PSNP5'}, inplace=True)

        # keep specific columns needed
        temp_positive_control_df=temp_positive_control_df[['CHR','POS','ABS_THETA','PSNP4','PSNP5']]
        temp_negative_control_df=temp_negative_control_df[['CHR','POS','ABS_THETA','PSNP4','PSNP5']]

        # melt down the dataframe for pos and neg
        temp_negative_control_df=pandas.melt(temp_negative_control_df, id_vars=['CHR','POS'],var_name="PVAL_TYPE",value_name="VALUE")
        temp_positive_control_df=pandas.melt(temp_positive_control_df, id_vars=['CHR','POS'],var_name="PVAL_TYPE",value_name="VALUE")

        # add in the necessary columns for subsample number  
        temp_positive_control_df["SUBSAMPLE_NUM"] = current_subsample_num
        temp_negative_control_df["SUBSAMPLE_NUM"] = current_subsample_num

    #concatonate either or both dataframes if they arent empty
    if len(temp_positive_control_df)>0:

        #positive_control_df =pandas.concat([positive_control_df,temp_positive_control_df], ignore_index=True)
        positive_control_df=pandas.concat([positive_control_df,temp_positive_control_df],ignore_index=True)

    if len(temp_negative_control_df)>0:

        #negative_control_df =pandas.concat([negative_control_df,temp_negative_control_df], ignore_index=True)
        negative_control_df=pandas.concat([negative_control_df,temp_negative_control_df],ignore_index=True)
    
    print("IDEA_2 function done",flush=True)
    
    ####
    ####
    ####
    # END OF FUNCTION
    # return positive and negative control dataframes
    return positive_control_df,negative_control_df

# IDEA 2.3.1
def IDEA_2_MAKE_R_AND_BASH_SCRIPT(
                        phenotype,
                        control_dataframe,
                        positive_or_negative
                        ):
    print("Entered FUNCTION: IDEA_2_MAKE_R_AND_BASH_SCRIPT",flush=True)
    control_dataframe.to_csv("/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Github/gift_subsampling/test_env/output_files/R_DATA/"+str(phenotype)+"_"+str(positive_or_negative)+"_control.csv", index = False)

    #####################################
    ####### MAKING THE R SCRIPT
    ######################################

    # write the R script to pair with the SNP and its data
    CURRENT_SNP_R_SCRIPT=open("/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Github/gift_subsampling/test_env/output_files/SNP_tracker_R_scripts/"+str(phenotype)+"_"+str(positive_or_negative)+"_control.R","w")
    CURRENT_SNP_R_SCRIPT.write(f'#R script for making box plots with ggplot\n')
    CURRENT_SNP_R_SCRIPT.write(f'library("tidyverse")\n')
    #CURRENT_SNP_R_SCRIPT.write(f'library("ggplot2")\n')
    CURRENT_SNP_R_SCRIPT.write(f'print("start of IDEA 2 R script")\n')
    CURRENT_SNP_R_SCRIPT.write(f'{positive_or_negative}_control_data<- read.csv("/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Github/gift_subsampling/test_env/output_files/R_DATA/{phenotype}_{positive_or_negative}_control.csv", header= TRUE, sep=",")\n')
    CURRENT_SNP_R_SCRIPT.write(f'\n')
    CURRENT_SNP_R_SCRIPT.write(f'#open png\n')
    CURRENT_SNP_R_SCRIPT.write(f'png("/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Github/gift_subsampling/test_env/output_files/summary_plots/IDEA2/{phenotype}_{positive_or_negative}_control.png", bg = "white", width = 9.75, height = 3.25, units = "in", res = 1200, pointsize = 4)\n')
    CURRENT_SNP_R_SCRIPT.write(f'ggplot({positive_or_negative}_control_data, aes(x=PVAL_TYPE, y=VALUE, fill=SUBSAMPLE_NUM)) +\n')
    CURRENT_SNP_R_SCRIPT.write(f'   geom_boxplot()\n')
    CURRENT_SNP_R_SCRIPT.write(f'# no facet_wrap for this code\n')
    CURRENT_SNP_R_SCRIPT.write(f'\n')
    CURRENT_SNP_R_SCRIPT.write(f'dev.off()\n')
    CURRENT_SNP_R_SCRIPT.write(f'print("End of IDEA 2 R script")\n')
    CURRENT_SNP_R_SCRIPT.write(f'# END OF R SCRIPT')
    CURRENT_SNP_R_SCRIPT.close()

    #####################################
    ####### MAKING THE BATCH SCRIPT
    ######################################

    # location where the batch scripts will be written to
    CURRENT_SNP_BATCH=open("/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Github/gift_subsampling/test_env/batch_files/parallel_stage2/"+str(phenotype)+"_"+str(positive_or_negative)+"_control.sh","w")
    
    # necessary start to the file
    CURRENT_SNP_BATCH.write(f'#!/bin/bash\n')
    CURRENT_SNP_BATCH.write(f'#SBATCH --partition=defq\n')
    CURRENT_SNP_BATCH.write(f'#SBATCH --nodes=1\n')
    CURRENT_SNP_BATCH.write(f'#SBATCH --ntasks=1\n')
    CURRENT_SNP_BATCH.write(f'#SBATCH --cpus-per-task=3\n')
    CURRENT_SNP_BATCH.write(f'#SBATCH --mem=8g\n')
    CURRENT_SNP_BATCH.write(f'#SBATCH --time=1:00:00\n')
    CURRENT_SNP_BATCH.write(f'#SBATCH --job-name=R_subrun\n')
    CURRENT_SNP_BATCH.write(f'#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out\n')
    CURRENT_SNP_BATCH.write(f'#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err\n')
    CURRENT_SNP_BATCH.write(f'#SBATCH --mail-type=ALL\n')
    CURRENT_SNP_BATCH.write(f'#SBATCH --mail-user=mbysh17@nottingham.ac.uk\n')
    CURRENT_SNP_BATCH.write(f'#===============================\n')
    CURRENT_SNP_BATCH.write(f'echo "start OF IDEA 2 btach script"\n')
    CURRENT_SNP_BATCH.write(f'#change to home directory\n')
    CURRENT_SNP_BATCH.write(f'cd /mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Github/gift_subsampling/test_env\n')
    CURRENT_SNP_BATCH.write(f'# source conda environments\n')
    CURRENT_SNP_BATCH.write(f'source ~/.bashrc\n')
    CURRENT_SNP_BATCH.write(f'conda deactivate\n')
    CURRENT_SNP_BATCH.write(f'conda activate gift_env\n')
    CURRENT_SNP_BATCH.write(f'# R SCRIPT FOR (IDEA 2) BOXPLOT\n')
    CURRENT_SNP_BATCH.write(f'Rscript /mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Github/gift_subsampling/test_env/output_files/SNP_tracker_R_scripts/{phenotype}_{positive_or_negative}_control.R\n')
    CURRENT_SNP_BATCH.write(f'conda deactivate\n')
    CURRENT_SNP_BATCH.write(f'echo "END OF IDEA 2 btach script"\n')
    CURRENT_SNP_BATCH.write(f'# end of file\n')
    CURRENT_SNP_BATCH.close()

    # END OF FUNCTION

Total_GIFT_Mo98 = 0
Total_GWAS_Mo98 = 0

Total_GIFT_Na23 = 0
Total_GWAS_Na23 = 0

#csv_files = os.listdir(args.d)
current_wd = os.getcwd()

print("Current wd is: ", current_wd,flush=True)

csv_files=[]

print("csv file list BEFORE sort: ======================================= \n",flush=True)
print(csv_files,flush=True)
print("\n ======================================= \n",flush=True)

for file in os.listdir("/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Github/gift_subsampling/test_env/output_files/"):
    if file.endswith(".csv") and file.__contains__("T20")==False:
        print(file,": ++++++++++++ ADDED ++++++++++++ !")
        csv_files.append(file)
    else:
        print(file,": ////////// SKIPPED //////// !!")


#csv_files = [f for f in os.listdir(str(args.d)) if os.path.isdir(f)]
#csv_files=os.listdir("/gpfs01/home/mbysh17/output_files")


# might need to turn this into a function soon (similar for loop used twice...)
# loop through each csv file in the given directory (which contain GWAS or GIFT data)

csv_file_index =0
csv_file_index_max=int(len(csv_files)-1)

for csv_file in csv_files: 

    csv_file_path=("/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Github/gift_subsampling/test_env/output_files/"+str(csv_file))

    print("////////////////////////////////////////////////",flush=True)
    completion_percentage=float((csv_file_index/csv_file_index_max)*100)
    print("1ST LOOP CSV FILE COMPLETION %: ", completion_percentage,flush=True)
    print(csv_file_path,flush=True)
    print("////////////////////////////////////////////////",flush=True)
    csv_file_index+=1
    # split on the underscore to make lists
    csv_file=csv_file.split("_")

    #check for the specific phenotype based on the naming convention
    # AND split code into GIFT and GWAS specific code
    
    # if statement 1.1
    if csv_file[2] == "Mo98" and csv_file[3]=='whole':
        print("Entered if statement 1.1",flush=True)

        #Total_GIFT+=1 (in function instead)
        Current_Mo98_dataframe=pandas.read_csv(csv_file_path)
    
        GWAS_or_GIFT = "GIFT"
        ####################################################
        # IDEA 3.1 #####################################
        
        Mo98_ALL_SNPS_GIFT_df,Total_GIFT_Mo98=IDEA_3_GATHER_ALL_SNPS_COMBINED(
                                                                GWAS_or_GIFT,
                                                                Current_Mo98_dataframe,
                                                                Mo98_ALL_SNPS_GIFT_df,
                                                                Total_GIFT_Mo98,
                                                                int(csv_file[6])
                                                                )

        ####################################################
        # IDEA 2.1 #####################################

        # MOT1 gene location boundaries 
        positive_control_chromosome = 2
        positive_control_LB = 10933005
        positive_control_UB = 10934604

        # rad50 gene location boundaries
        negative_control_chromosome = 2
        negative_control_LB =13600431
        negative_control_UB =13609104

        positive_control_Mo98_df, negative_control_Mo98_df=IDEA_2_CONTROL_CHECK(Current_Mo98_dataframe,
                                                                positive_control_Mo98_df,
                                                                positive_control_chromosome,
                                                                positive_control_LB,
                                                                positive_control_UB,
                                                                negative_control_Mo98_df,
                                                                negative_control_chromosome,
                                                                negative_control_LB,
                                                                negative_control_UB,
                                                                GWAS_or_GIFT,
                                                                int(csv_file[6]) #SUBSAMPLE NUMBER FOR GIFT FILE
                                                                )


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
        BT_pSNP4 = 0.05 / (len(Current_Mo98_dataframe) * int(csv_file[6])) # multiply by subsample num which is in file name
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
        BT_pSNP5 = 0.05 / (len(Current_Mo98_dataframe) * int(csv_file[6])) # multiply by subsample num which is in file name
        BF_THRES_pSNP5 = -(math.log(BT_pSNP5,10))

         # threshold (3.1) absolute theta threshold?
         # awaiting confirmation on how to threshold this value
         # place holder bonferroni correction
        BT_abstheta = 0.05 / (len(Current_Mo98_dataframe) * int(csv_file[6])) # multiply by subsample num which is in file name
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
            if positive_control_lower_bound_position<=int(current_row[POS])<=positive_control_upper_bound_position:
                # then check if its significant or not for each type of pval
                if float(current_row["pSNP4"])>= BHY_pSNP4_thres:
                    #print("This SNP is significant for pSNP4 under BHY")
                    # check if this SNP with its CHR POS AND SUBSAMPLE_NUM is already in the dataframe
                    already_there = False
                    for inner_index, inner_row in positive_control_Mo98_df.iterrows():
                        if int(current_row[POS]) == int(inner_row['POS']) and int(current_row[CHR]) == int(inner_row['CHR']) and int(csv_file[6]) == int(inner_row['SUBSAMPLE_NUM']):
                                already_there = True
                if float(current_row["pSNP4"])>= BF_THRES_pSNP4:
                    #print("This SNP is significant for pSNP4 under BF")
                    pass
        '''

    # if statement 1.2     
    elif csv_file[2]=="Mo98" and csv_file[3]=='GWAS' and csv_file[4]!="T20":# GWAS code vvv      
        
        print("Entered if statement 1.2",flush=True)
        #Total_GWAS+=1
        Current_Mo98_dataframe=pandas.read_csv(csv_file_path)

        GWAS_or_GIFT="GWAS"

        Mo98_ALL_SNPS_GWAS_df,Total_GWAS_Mo98 = IDEA_3_GATHER_ALL_SNPS_COMBINED(GWAS_or_GIFT,Current_Mo98_dataframe,Mo98_ALL_SNPS_GWAS_df,Total_GWAS_Mo98, int(csv_file[4]))

        # MOT1 gene location boundaries 
        positive_control_chromosome = 2
        positive_control_LB = 10933005
        positive_control_UB = 10934604

        # rad50 gene location boundaries
        negative_control_chromosome = 2
        negative_control_LB =13600431
        negative_control_UB =13609104

        positive_control_Mo98_df, negative_control_Mo98_df= IDEA_2_CONTROL_CHECK(Current_Mo98_dataframe,
                                                                    positive_control_Mo98_df,
                                                                    positive_control_chromosome,
                                                                    positive_control_LB,
                                                                    positive_control_UB,
                                                                    negative_control_Mo98_df,
                                                                    negative_control_chromosome,
                                                                    negative_control_LB,
                                                                    negative_control_UB,
                                                                    GWAS_or_GIFT,
                                                                    int(csv_file[4]) #SUBSAMPLE NUMBER FOR GIFT FILE
                                                                    )
        
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
        BT_p = 0.05 / (len(Current_Mo98_dataframe) * int(csv_file[6])) # multiply by subsample num which is in file name
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

            if positive_control_lower_bound_position<=int(current_row[POS])<=positive_control_upper_bound_position:
                # then check if its significant or not for each type of pval
                if float(current_row["pSNP4"])>= BHY_pSNP4_thres:
                    print("This SNP is significant for pSNP4 under BHY")
                    # check if this SNP with its CHR POS AND SUBSAMPLE_NUM is already in the dataframe
                    already_there = False
                    for inner_index, inner_row in positive_control_Mo98_df.iterrows():
                        if int(current_row[POS]) == int(inner_row['POS']) and int(current_row[CHR]) == int(inner_row['CHR']) and int(csv_file[6]) == int(inner_row['SUBSAMPLE_NUM']):
                                already_there = True
                if float(current_row["pSNP4"])>= BF_THRES_pSNP4:
                    print("This SNP is significant for pSNP4 under BF")

        '''

    # if statement 1.3
    elif csv_file[2] == "Na23" and csv_file[3]=='whole':

        print("Entered if statement 1.3",flush=True)
        Current_Na23_dataframe=pandas.read_csv(csv_file_path) 

        GWAS_or_GIFT="GIFT"

        # TRACE 1
        Na23_ALL_SNPS_GIFT_df,Total_GIFT_Na23=IDEA_3_GATHER_ALL_SNPS_COMBINED(GWAS_or_GIFT,
                                                                              Current_Na23_dataframe,
                                                                              Na23_ALL_SNPS_GIFT_df,
                                                                              Total_GIFT_Na23,
                                                                              int(csv_file[6]))

        #chr4:6,391,854-6,395,922
        
        # HKT1 gene location boundaries 
        positive_control_chromosome = 4
        positive_control_LB = 6391854
        positive_control_UB = 6395922

        # rad50 gene location boundaries
        negative_control_chromosome = 2
        negative_control_LB =13600431
        negative_control_UB =13609104

    
        positive_control_Na23_df, negative_control_Na23_df=IDEA_2_CONTROL_CHECK(Current_Na23_dataframe,
                                                                    positive_control_Na23_df,
                                                                    positive_control_chromosome,
                                                                    positive_control_LB,
                                                                    positive_control_UB,
                                                                    negative_control_Na23_df,
                                                                    negative_control_chromosome,
                                                                    negative_control_LB,
                                                                    negative_control_UB,
                                                                    GWAS_or_GIFT,
                                                                    int(csv_file[6]) 
                                                                    )

    # if statement 1.4
    elif csv_file[2] == "Na23" and csv_file[3]=='GWAS' and csv_file[4]!="T20": #
        print("Entered if statement 1.4",flush=True)
        
        #Total_GWAS+=1
        Current_Na23_dataframe=pandas.read_csv(csv_file_path)

        GWAS_or_GIFT="GWAS"

        Na23_ALL_SNPS_GWAS_df,Total_GWAS_Na23 = IDEA_3_GATHER_ALL_SNPS_COMBINED(GWAS_or_GIFT,
                                                                                Current_Na23_dataframe,
                                                                                Na23_ALL_SNPS_GWAS_df,
                                                                                Total_GWAS_Na23, 
                                                                                int(csv_file[4]))

        # HKT1 gene location boundaries 
        positive_control_chromosome = 4
        positive_control_LB = 6391854
        positive_control_UB = 6395922

        #MUT1 homologue (MLH1) gene location boundaries
        # chr4:5,816,942-5,821,066
        negative_control_chromosome = 4
        negative_control_LB =5816942
        negative_control_UB =5821066
        
        positive_control_Na23_df, negative_control_Na23_df= IDEA_2_CONTROL_CHECK(Current_Na23_dataframe,
                                                                    positive_control_Na23_df,
                                                                    positive_control_chromosome,
                                                                    positive_control_LB,
                                                                    positive_control_UB,
                                                                    negative_control_Na23_df,
                                                                    negative_control_chromosome,
                                                                    negative_control_LB,
                                                                    negative_control_UB,
                                                                    GWAS_or_GIFT,
                                                                    int(csv_file[4]) #SUBSAMPLE NUMBER FOR GIFT FILE
                                                                    )



# end time
t2 = time.perf_counter()

print(f"FINISHED FIRST LOOP IN {t2-t1} seconds!")


####################################################
# IDEA 2.2 #####################################
############################################
IDEA_2_MAKE_R_AND_BASH_SCRIPT("Mo98", positive_control_Mo98_df,"positive")
IDEA_2_MAKE_R_AND_BASH_SCRIPT("Mo98", negative_control_Mo98_df,"negative")
IDEA_2_MAKE_R_AND_BASH_SCRIPT("Na23", positive_control_Na23_df,"positive")
IDEA_2_MAKE_R_AND_BASH_SCRIPT("Na23", negative_control_Na23_df,"negative")

print("IDEA 2 FINISHED",flush=True)
############################################
# IDEA 2.2 #####################################
####################################################


####################################################
# IDEA 3.2 #####################################
############################################
# Now go through each of the total SNP dataframes and calculate the average P values for each SNP
Mo98_ALL_SNPS_GIFT_df = IDEA_3_CALCULATE_AVERAGE_SNPS_GIFT(Mo98_ALL_SNPS_GIFT_df)
Mo98_ALL_SNPS_GWAS_df = IDEA_3_CALCULATE_AVERAGE_SNPS_GWAS(Mo98_ALL_SNPS_GWAS_df)

Na23_ALL_SNPS_GIFT_df = IDEA_3_CALCULATE_AVERAGE_SNPS_GIFT(Na23_ALL_SNPS_GIFT_df)
Na23_ALL_SNPS_GWAS_df = IDEA_3_CALCULATE_AVERAGE_SNPS_GWAS(Na23_ALL_SNPS_GWAS_df)

############################################
# IDEA 3.2 #####################################
####################################################

####################################################
# IDEA 3.3 #####################################
############################################

# testing
print("Idea 3 dataframes look like this: \n",Mo98_ALL_SNPS_GIFT_df,"\n",Mo98_ALL_SNPS_GWAS_df,"\n", Na23_ALL_SNPS_GIFT_df,"\n",Na23_ALL_SNPS_GWAS_df,flush=True)

# write results to csv files for an R script to run them
Mo98_ALL_SNPS_GIFT_df.to_csv("/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Github/gift_subsampling/test_env/output_files/R_DATA/Mo98_ALL_SNPS_GIFT.csv", index=False)
Mo98_ALL_SNPS_GWAS_df.to_csv("/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Github/gift_subsampling/test_env/output_files/R_DATA/Mo98_ALL_SNPS_GWAS.csv", index=False)

Na23_ALL_SNPS_GIFT_df.to_csv("/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Github/gift_subsampling/test_env/output_files/R_DATA/Na23_ALL_SNPS_GIFT.csv", index=False)
Na23_ALL_SNPS_GWAS_df.to_csv("/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Github/gift_subsampling/test_env/output_files/R_DATA/Na23_ALL_SNPS_GWAS.csv", index=False)

############################################
# IDEA 3.3 #####################################
####################################################

####################################################
# IDEA 3.4 #####################################
############################################
subsample_num_list=[200,400,600,800,1000] # can later update this to read from earlier scripts or something

phenotype_list=["Mo98","Na23"] # can later update this to read from the phenotype text file

pvals=["TOTAL_P","TOTAL_PSNP4","TOTAL_PSNP5","TOTAL_ABS_THETA"] # these should always be the same 4 types
# R script creation and running
for phenotype in phenotype_list:

    for subsample_number in subsample_num_list:

        for pval_type in pvals:

            # Example variable inputs are as follows: 
            # example 1 Mo98,200,TOTAL_P
            # example 2 Mo98,200,TOTAL_PSNP4
            # ...
            # example x Mo98,400,TOTAL_P
            # ....
            IDEA_3_R_AND_BATCH(phenotype,subsample_number,pval_type)

print("IDEA 3 FINISHED",flush=True)
############################################
# IDEA 3.4 #####################################
####################################################


####################################################
# IDEA 1.1 ######################################
average_pvals_list=["AVERAGE_P","AVERAGE_PSNP4","AVERAGE_PSNP5","AVERAGE_ABS_THETA"]
for phenotype in phenotype_list:

    for pval_type in average_pvals_list:

        if phenotype=="Mo98":

            if pval_type == "AVERAGE_P":

                # get T20 SNPs locations for GWAS_P for 1000 subsample
                Mo98_GWAS_P_T20_SNP_locations_df = GET_T20_LOCATIONS_AT_1000(Mo98_ALL_SNPS_GWAS_df,pval_type) 

            elif pval_type == "AVERAGE_PSNP4":

                Mo98_PSNP4_T20_SNP_locations_df = GET_T20_LOCATIONS_AT_1000(Mo98_ALL_SNPS_GIFT_df,pval_type) 

            elif pval_type == "AVERAGE_PSNP5":

                Mo98_PSNP5_T20_SNP_locations_df = GET_T20_LOCATIONS_AT_1000(Mo98_ALL_SNPS_GIFT_df,pval_type) 

            elif pval_type == "AVERAGE_ABS_THETA":

                Mo98_ABS_THETA_T20_SNP_locations_df = GET_T20_LOCATIONS_AT_1000(Mo98_ALL_SNPS_GIFT_df,pval_type) 

        elif phenotype == "Na23":
            
            if pval_type == "AVERAGE_P":

                # get T20 SNPs locations for GWAS_P for 1000 subsample
                Na23_GWAS_P_T20_SNP_locations_df = GET_T20_LOCATIONS_AT_1000(Na23_ALL_SNPS_GWAS_df,pval_type) 

            elif pval_type == "AVERAGE_PSNP4":

                Na23_PSNP4_T20_SNP_locations_df = GET_T20_LOCATIONS_AT_1000(Na23_ALL_SNPS_GIFT_df,pval_type) 

            elif pval_type == "AVERAGE_PSNP5":

                Na23_PSNP5_T20_SNP_locations_df = GET_T20_LOCATIONS_AT_1000(Na23_ALL_SNPS_GIFT_df,pval_type) 

            elif pval_type == "AVERAGE_ABS_THETA":
                
                Na23_ABS_THETA_T20_SNP_locations_df = GET_T20_LOCATIONS_AT_1000(Na23_ALL_SNPS_GIFT_df,pval_type) 

# might need to turn csv reading into a function...

####################################################
# IDEA 1.2 ######################################

# set up cumulative T20 dataframe for each phenotype

Mo98_cumulative_t20_dataframe = pandas.DataFrame(columns=[
                                            'CHR',
                                            'POS',
                                            'PVAL_TYPE',
                                            'SUBSAMPLE_NUM',
                                            'VALUE',
                                            ])

Na23_cumulative_t20_dataframe = pandas.DataFrame(columns=[
                                            'CHR',
                                            'POS',
                                            'PVAL_TYPE',
                                            'SUBSAMPLE_NUM',
                                            'VALUE',
                                            ])

####################################################
# IDEA 1.3.X ######################################
csv_file_index=0

for csv_file in csv_files: 

    csv_file_path=("/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Github/gift_subsampling/test_env/output_files/"+str(csv_file))

    # split on the underscore to make lists
    csv_file=csv_file.split("_")

    print("////////////////////////////////////////////////",flush=True)
    completion_percentage=float((csv_file_index/csv_file_index_max)*100)
    print("2ND LOOP CSV FILE COMPLETION %: ", completion_percentage,flush=True)
    print(csv_file_path,flush=True)
    print("////////////////////////////////////////////////",flush=True)
    csv_file_index+=1

    #check for the specific phenotype based on the naming convention
    # AND split code into GIFT and GWAS specific code
    

    if csv_file[2] == "Mo98" and csv_file[3]=='whole':# GIFT code vvvvv
        
        Current_Mo98_dataframe=pandas.read_csv(csv_file_path)
       
        # set method and level variables
        GWAS_or_GIFT = "GIFT"
        subsample_level = int(csv_file[6]) 

        # testing
        '''
        this_csv="_".join(csv_file)
        print("Currently reading file: ",this_csv,flush=True)
        print(Current_Mo98_dataframe.head(),flush=True)
        '''
        # update the cumulative t20 dataframe
        Mo98_cumulative_t20_dataframe = IDEA_1_ACCUMULATE_T20_SNP_DATA(
            Current_Mo98_dataframe,
            Mo98_GWAS_P_T20_SNP_locations_df,
            Mo98_PSNP4_T20_SNP_locations_df, # passed in but not needed
            Mo98_PSNP5_T20_SNP_locations_df, # passed in but not needed
            Mo98_ABS_THETA_T20_SNP_locations_df, # passed in but not needed 
            Mo98_cumulative_t20_dataframe,
            GWAS_or_GIFT,
            subsample_level
            )

    elif csv_file[2]=="Mo98" and csv_file[3]=='GWAS' and csv_file[4]!="T20":# GWAS code vvv      
        
        Current_Mo98_dataframe=pandas.read_csv(csv_file_path)
        
        # set method and level variables
        GWAS_or_GIFT = "GWAS"
        subsample_level = int(csv_file[4])

        # testing
        '''
        this_csv="_".join(csv_file)
        print("Currently reading file: ",this_csv,flush=True)
        print(Current_Mo98_dataframe.head(),flush=True)
        '''

        # update the cumulative t20 dataframe
        Mo98_cumulative_t20_dataframe = IDEA_1_ACCUMULATE_T20_SNP_DATA(
            Current_Mo98_dataframe,
            Mo98_GWAS_P_T20_SNP_locations_df,
            Mo98_PSNP4_T20_SNP_locations_df, # passed in but not needed
            Mo98_PSNP5_T20_SNP_locations_df, # passed in but not needed
            Mo98_ABS_THETA_T20_SNP_locations_df, # passed in but not needed 
            Mo98_cumulative_t20_dataframe,
            GWAS_or_GIFT,
            subsample_level
            )

    elif csv_file[2] == "Na23" and csv_file[3]=='whole':

        Current_Na23_dataframe=pandas.read_csv(csv_file_path) 

        # set method and level variables
        GWAS_or_GIFT = "GIFT"
        subsample_level = int(csv_file[6])

        # testing
        '''
        this_csv="_".join(csv_file)
        print("Currently reading file: ",this_csv,flush=True)
        print(Current_Na23_dataframe.head(),flush=True)
        '''
        # update the cumulative t20 dataframe
        Na23_cumulative_t20_dataframe = IDEA_1_ACCUMULATE_T20_SNP_DATA(
            Current_Na23_dataframe,
            Na23_GWAS_P_T20_SNP_locations_df,
            Na23_PSNP4_T20_SNP_locations_df, # passed in but not needed
            Na23_PSNP5_T20_SNP_locations_df, # passed in but not needed
            Na23_ABS_THETA_T20_SNP_locations_df, # passed in but not needed 
            Na23_cumulative_t20_dataframe,
            GWAS_or_GIFT,
            subsample_level
            )

    elif csv_file[2] == "Na23" and csv_file[3]=='GWAS' and csv_file[4]!="T20": 

        Current_Na23_dataframe=pandas.read_csv(csv_file_path) 

        # set method and level variables
        GWAS_or_GIFT = "GWAS"
        subsample_level = int(csv_file[4])

        # testing
        '''
        this_csv="_".join(csv_file)
        print("Currently reading file: ",this_csv,flush=True)
        print(Current_Na23_dataframe.head(),flush=True)
        '''

        # update the cumulative t20 dataframe
        Na23_cumulative_t20_dataframe = IDEA_1_ACCUMULATE_T20_SNP_DATA(
            Current_Na23_dataframe,
            Na23_GWAS_P_T20_SNP_locations_df,
            Na23_PSNP4_T20_SNP_locations_df, # passed in but not needed
            Na23_PSNP5_T20_SNP_locations_df, # passed in but not needed
            Na23_ABS_THETA_T20_SNP_locations_df, # passed in but not needed 
            Na23_cumulative_t20_dataframe,
            GWAS_or_GIFT,
            subsample_level
            )

####################################################
# IDEA 1.4 ######################################
IDEA_1_MAKE_R_SCRIPT("Mo98",Mo98_cumulative_t20_dataframe)
IDEA_1_MAKE_R_SCRIPT("Na23",Na23_cumulative_t20_dataframe)

####################################################
# IDEA 1.5 ######################################
IDEA_1_MAKE_BASH_SCRIPT("Mo98")
IDEA_1_MAKE_BASH_SCRIPT("Na23")

print("END OF IDEA 1",flush=True)

print("END OF SNP_TRACKER.PY",flush=True)
# end of file