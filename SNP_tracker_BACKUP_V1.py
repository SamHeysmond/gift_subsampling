# File version: 03_06_2024

# This file:
#   1) Tracks all values of top 20 significant PVAL types for each pval (PSNP4, GWAS_P etc...) 
#      across multiple sample levels to put into a box plot
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
import argparse, os, math
import concurrent.futures
import modin.pandas as pandas
import ray

ray.init(_plasma_directory="/tmp", object_store_memory=600000000000) # setting to disable out of core in Ray and obj store mem increase
# obj store memory currently at 600GB


# so i can see all the columns when testing with print
pandas.set_option('display.max_columns',None)
pandas.options.display.max_columns=None


# placeholder till parseargs will work
# will implement an argument that inputs home user directory automatically
# -> or lets user decide
PATH_TO_MAIN = "/gpfs01/home/mbysh17/"


# Reminder of CSV format (GIFT) NAME leaf_ionome_Mo98_whole_genome_metrics_600_732692.csv
# CHROM,POS,largest_theta,smallest_theta,absolute_theta,theta_range,largest_relative_theta,smallest_relative_theta,absolute_relative_theta,range_relative_theta,min_p,mean_p,log_mean_p,bigest_theta_p,pSNP4,pSNP5
# 1,73,6.285714285714285,-4.224489795918359,6.285714285714285,10.510204081632644,0.3845193508114856,-0.25842696629213435,0.3845193508114856,0.6429463171036199,4.185363300872138e-05,nan,nan,0.00015134587316541535,0.15259768149369662,0.6580333017260325

# reminder of csv format (GWAS) NAME leaf_ionome_Mo98_GWAS_600_732692.csv
# chromosomes,positions,pvals,mafs,macs,genotype_var_perc
# 1,55,0.621946624343,0.0016694490818,1,0.000407516956329

# IDEA 3.1
# combined in the sense that it stacks/combines snps from same location and works for GWAS and GIFT
# TRACE 2

# IDEA 3 NEW (CHECKED!)
def IDEA_3_process_all_snps_file(df_to_process,GIFT_or_GWAS,TOTAL_GIFT_OR_GWAS,subsample_num):
    print("////////////////////////////////////////////////",flush=True)
    print("Entered IDEA3 Function",flush=True)
    print("////////////////////////////////////////////////",flush=True)

    if GIFT_or_GWAS == "GWAS":
        CHR="chromosomes"
        POS="positions"

        df_to_process.rename(columns={CHR:'CHR',POS:'POS','pvals':'TOTAL_P'},inplace=True)

        df_to_process["TIMES_APPEARED"] = 1
        df_to_process["TOTAL_GWAS"]=TOTAL_GIFT_OR_GWAS # should be 100 at EACH subsample level e.g. at 800 sample -> 100 total GWAS
        df_to_process.insert(2,"SUBSAMPLE_NUM",subsample_num)

        df_to_process = df_to_process.groupby(['CHR','POS','SUBSAMPLE_NUM'],as_index=False).agg({'TOTAL_P':'sum','TIMES_APPEARED':'sum','TOTAL_GWAS':'max'})
        
        df_to_process["AVERAGE_P"] = df_to_process["TOTAL_P"] / df_to_process["TIMES_APPEARED"]
    
    elif GIFT_or_GWAS=="GIFT":

        CHR ="CHROM"
        POS="POS"

        df_to_process.rename(columns={CHR:'CHR',POS:'POS','pSNP4':'TOTAL_PSNP4','pSNP5':'TOTAL_PSNP5','absolute_theta':'TOTAL_ABS_THETA'}, inplace=True)
        
        # adding in two new columns at specific index values
        df_to_process["TIMES_APPEARED"] = 1
        df_to_process["TOTAL_GIFT"] = TOTAL_GIFT_OR_GWAS
        df_to_process.insert(2,"SUBSAMPLE_NUM",subsample_num)
        
        # group by and sort the data (THIS TAKES A WHILE AND LOTS OF MEMORY)
        df_to_process = df_to_process.groupby(['CHR','POS','SUBSAMPLE_NUM'],as_index=False).agg({'TOTAL_PSNP4':'sum','TOTAL_PSNP5':'sum','TOTAL_ABS_THETA':'sum','TIMES_APPEARED':'sum','TOTAL_GIFT':'max'})
            
        df_to_process["AVERAGE_PSNP4"] = df_to_process["TOTAL_PSNP4"] / df_to_process["TIMES_APPEARED"]
        df_to_process["AVERAGE_PSNP5"] = df_to_process["TOTAL_PSNP5"] / df_to_process["TIMES_APPEARED"]
        df_to_process["AVERAGE_ABS_THETA"] =df_to_process["TOTAL_ABS_THETA"] / df_to_process["TIMES_APPEARED"]

    return df_to_process

# IDEA 3 (CHECKED! X2)
def IDEA_3_R_AND_BATCH(phenotype,subsample_number,pval_type):
    print("Entered FUNCTION: IDEA_3_R_AND_BATCH",flush=True)
    #########################################
    ### MAKE R SCRIPT ##################
    #############################

    # make R script for each P value type (pSNP4, pSNP5, abs theta, GWAS_P)
    R_out=open(PATH_TO_MAIN+"output_files/SNP_tracker_R_scripts/"+str(phenotype)+"_"+str(subsample_number)+"_"+str(pval_type)+"_MANHATTAN.R","w")
    R_out.write(f'#R script for making manhattan plots with ggplot\n')
    R_out.write(f'library("tidyverse")\n')
    #R_out.write(f'library("ggplot2")\n')
    #R_out.write(f'library("dplyr")\n')
    # R_out.write(f'library("ggrepel")\n')   # Dont need ggrepel for now
    R_out.write(f'print("Start of IDEA 3 R")\n')
    R_out.write(f'\n')

    # for GWAS data....
    if pval_type=="AVERAGE_P":

        # fetch the data of the csv for the current phenotype and current method (GWAS)
        R_out.write(f'{pval_type}_SUBSAMPLE_{subsample_number}_SNPS_DATA<-read.csv("{PATH_TO_MAIN}output_files/R_DATA/{phenotype}_GWAS_{subsample_number}_ALL.csv",header=TRUE)\n')  
        
    # for GIFT data
    else:

        # fetch the data of the csv for the current phenotype and current method (GIFT)
        R_out.write(f'{pval_type}_SUBSAMPLE_{subsample_number}_SNPS_DATA<-read.csv("{PATH_TO_MAIN}output_files/R_DATA/{phenotype}_GIFT_{subsample_number}_ALL.csv",header=TRUE)\n')  

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
    if pval_type=="AVG_ABS_THETA": # testing without log10
        R_out.write(f'ylim <- abs(floor(min(mydata))) +1\n')
    else:
        R_out.write(f'ylim <- abs(floor(log10(min(mydata)))) +1\n')

    R_out.write(f'#open png\n')
    # Send output to IDEA 3 summary plot folder
    R_out.write(f'png("{PATH_TO_MAIN}output_files/summary_plots/IDEA3/{phenotype}_{subsample_number}_{pval_type}_MANHATTAN.png", bg = "white", width = 9.75, height = 3.25, units = "in", res = 1200, pointsize = 4)\n')

    # make the plot
    # this is where the plot depends on the pval type. 
    # for absolute theta, there is no -log10 since it isnt quite a "pvalue"
    if pval_type=="AVG_ABS_THETA":
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
    R_batch=open(PATH_TO_MAIN+"batch_files/parallel_stage2/"+str(phenotype)+"_"+str(subsample_number)+"_"+str(pval_type)+"_MANHATTAN.sh","w")
    
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
    R_batch.write(f'Rscript {PATH_TO_MAIN}output_files/SNP_tracker_R_scripts/{phenotype}_{subsample_number}_{pval_type}_MANHATTAN.R\n')
    R_batch.write(f'conda deactivate\n')
    R_batch.write(f'echo "End of IDEA 3 batch script"\n')
    R_batch.write(f'conda deactivate\n')
    R_batch.write(f'# END OF FILE\n')
    R_batch.close()
    # end of function

# IDEA 1.1 (CHECKED! x2)
def GET_T20_LOCATIONS_AT_1000(dataframe_name,pval_type,phenotype): #parameters of: phenotype, subsample number, and pval type (implied method)
    print("Entered FUNCTION: IDEA_1_GET_T20_LOCATIONS_AT_1000",flush=True)

    # read from csv
    temp_all_snps_dataframe = pandas.read_csv(PATH_TO_MAIN+"output_files/R_DATA/"+dataframe_name)
    # 1) sort GWAS all SNPs by 1000 subsample and the PVAL (so biggest is at the top)
    #Mo98_ALL_SNPS_GWAS_df.sort_values(by=["SUBSAMPLE_NUM","AVERAGE_P"], axis=0, ascending=[False,True],inplace=True, na_position='first')
    if pval_type=="AVERAGE_ABS_THETA": # sort to get highest at top
         temp_all_snps_dataframe.sort_values(by=pval_type, axis=0, ascending=False,inplace=True, na_position='first')
    
    # otherwise sort to lowest (for pval)
    else:
        temp_all_snps_dataframe.sort_values(by=pval_type, axis=0, ascending=True,inplace=True, na_position='first')

    # 2) take top 20 of these averaged GWAS values
    # example format:
    # CHR    POS    SUBSAMPLE_NUM    TOTAL_P    TIMES_APPEARED    TOTAL_GWAS    AVERAGE_P
    # 1      342    1000             7          98                100           0.05
    # 4      787    1000             12         99                100           0.90
    current_pval_T20_df=temp_all_snps_dataframe.head(20)

    # delete the temp variable
    del temp_all_snps_dataframe

    # keep ONLY the location information in the t20 variable
    current_pval_T20_df=current_pval_T20_df[['CHR','POS']]

    # Write dataframe to csv
    current_pval_T20_df_path = PATH_TO_MAIN+"output_files/R_DATA/"+phenotype+"_"+pval_type+"_T20_LOCATIONS.csv"
    current_pval_T20_df.to_csv(current_pval_T20_df_path,header=True,index=False)

    # delete the variable to save memory
    del current_pval_T20_df

    # return the name of the dataframe (as a path)
    return current_pval_T20_df_path

# IDEA 1.3.1 (NEEDS FIX)
def IDEA_1_ACCUMULATE_T20_SNP_DATA(current_dataframe_main,
                                    GWAS_P_locations_dataframe_path,
                                    PSNP4_locations_dataframe_path,
                                    PSNP5_locations_dataframe_path,
                                    ABS_THETA_locations_dataframe_path,
                                    cumulative_t20_dataframe_path,
                                    GWAS_or_GIFT,
                                    subsample_level
                                    ):
    print("Entered FUNCTION: IDEA_1_ACCUMULATE_T20_SNP_DATA",flush=True)

    # read the current cumulative t20 dataframe
    try:
        cumulative_t20_dataframe=pandas.read_csv(cumulative_t20_dataframe_path)
    except:
        print("Current t20 dataframe could not be read (not yet made or error occured). This message should only appear twice")
        cumulative_t20_dataframe=pandas.DataFrame(columns=[
                                                    'CHR', 
                                                    'POS',
                                                    'PVAL_TYPE',  # for each subsamp number 200-1000
                                                    'SUBSAMPLE_NUM',
                                                    'VALUE' #there was an extra comma here- oops
                                                    ])

   

    if GWAS_or_GIFT=="GWAS":
        # read in current dataframe (but only certain columns)
        current_dataframe=current_dataframe_main[['chromosomes','positions','pvals']].copy()
        # read in the location dataframe
        locations_dataframe=pandas.read_csv(GWAS_P_locations_dataframe_path)

        current_dataframe.rename(columns={'chromosomes':'CHR','positions':'POS','pvals':'VALUE'}, inplace=True)

        # Should have columns
        #   CHR, POS, VALUE
        print("GWAS current_dataframe with column name change",flush=True)
        print(current_dataframe.head(),flush=True)

        df_out =(locations_dataframe.reset_index(drop=True)[["CHR", "POS"]].merge(current_dataframe.reset_index(drop=True), on=["CHR", "POS"], how="inner",left_index=False, right_index=False,))

        # insert the columns for current PVAL type and subsample number
        #   insert the pval type e.g. PSNP4 (all the way down this dataframe)
        df_out.insert(2,"PVAL_TYPE","GWAS_P")

        #   insert the subsample level (all the way down this dataframe)
        df_out.insert(3,"SUBSAMPLE_NUM",subsample_level)

        # Should have columns
        #   CHR, POS, PVAL_TYPE, SUBSAMPLE_NUM, VALUE
        print("GWAS df_out head",flush=True)
        print(df_out.head(),flush=True)

        # concatonate it with the main t20 dataframe
        cumulative_t20_dataframe=pandas.concat([cumulative_t20_dataframe,df_out])

        print("GWAS cumulative t20 dataframe head:",flush=True)
        print(cumulative_t20_dataframe.head(),flush=True)

    elif GWAS_or_GIFT=="GIFT":
        
        GIFT_locations_dataframes=[PSNP4_locations_dataframe_path,PSNP5_locations_dataframe_path,ABS_THETA_locations_dataframe_path]
        GIFT_column_to_change=['pSNP4','pSNP5','absolute_theta']
        GIFT_column_to_keep=['PSNP4','PSNP5',"ABS_THETA"]
        my_index = 0

        for my_index in range(0,3): # 0 , 1 , 2  STOP

            # take a copy of the curent dataframe instead (only keeping the correct pval type each iteration)
            current_dataframe=current_dataframe_main[['CHROM','POS',GIFT_column_to_change[my_index]]].copy()
            
            # fetch appropriate locations dataframe
            locations_dataframe = pandas.read_csv(GIFT_locations_dataframes[my_index])

            # rename chromosome column in the current dataframe (copied from main)
            current_dataframe.rename(columns={'CHROM':'CHR',GIFT_column_to_change[my_index]:"VALUE"}, inplace=True)

            # MERGE where the main dataframe contains locations of the locations dataframe (based on CHR and POS)
            df_out =(locations_dataframe.reset_index(drop=True)[["CHR", "POS"]].merge(current_dataframe.reset_index(drop=True), on=["CHR", "POS"], how="inner",left_index=False, right_index=False))

            print(f"Header of GIFT merged for {GIFT_column_to_keep[my_index]}",flush=True)
            # Should have following columns: (example: pSNP4)
            # CHR, POS, VALUE
            # 1,   24,   0.00213
            # ..,   .. ,  ......
            print(df_out.head(),flush=True)

            # insert the columns for current PVAL type and subsample number
            #   insert the pval type e.g. PSNP4 (all the way down this dataframe)
            df_out.insert(2,"PVAL_TYPE",GIFT_column_to_keep[my_index])

            #   insert the subsample level (all the way down this dataframe)
            df_out.insert(3,"SUBSAMPLE_NUM",subsample_level)

            print(f"GIFT merged with PVAL_TYPE and SUBSAMPLE_NUM columns",flush=True)
            # Should have following columns: (example: pSNP4)
            # CHR, POS, PVAL_TYPE, SUBSAMPLE_NUM,VALUE
            # 1,   24,   PSNP4,      200     ,   0.00213
            # ..,   .. ,  ......
            print(df_out.head(),flush=True)

            # concat the current pval data to main dataframe e.g. PSNP4 stuff
            cumulative_t20_dataframe=pandas.concat([cumulative_t20_dataframe,df_out])

            # testing 
            print("Cumulative dataframe head",flush=True)
            print(cumulative_t20_dataframe.head(),flush=True)
    
            # clear variables to save on memory
            del df_out

    #write the cumulative_t20 dataframe to csv
    cumulative_t20_dataframe.to_csv(cumulative_t20_dataframe_path,header=True,index=False)

    # save on memory by deleting the dataframe variable
    del cumulative_t20_dataframe

# IDEA 1.4.1 (CHECKED!)
def IDEA_1_MAKE_R_SCRIPT(
        phenotype,
        cumulative_t20_dataframe_name,
        ):
    print("Entered FUNCTION: IDEA_1_MAKE_R_SCRIPT",flush=True)

    #####################################
    ####### MAKING THE R SCRIPT
    ######################################

    cumulative_t20_dataframe_name = cumulative_t20_dataframe_name.replace(".csv","")

    # write the R script to pair with the SNP and its data
    CURRENT_SNP_R_SCRIPT=open(PATH_TO_MAIN+"output_files/SNP_tracker_R_scripts/"+str(cumulative_t20_dataframe_name)+".R","w")
    CURRENT_SNP_R_SCRIPT.write(f'#R script for making box plots with ggplot\n')
    #CURRENT_SNP_R_SCRIPT.write(f'library("ggplot2")\n')
    CURRENT_SNP_R_SCRIPT.write(f'library("tidyverse")\n')
    CURRENT_SNP_R_SCRIPT.write(f'print("start of IDEA1 R script")\n')
    CURRENT_SNP_R_SCRIPT.write(f'T20_TRACKED_DATA<- read.csv("{PATH_TO_MAIN}output_files/R_DATA/{cumulative_t20_dataframe_name}.csv", header= TRUE, sep=",")\n')
    CURRENT_SNP_R_SCRIPT.write(f'# Convert the subsample number to a FACTOR variable\n')
    CURRENT_SNP_R_SCRIPT.write(f'T20_TRACKED_DATA$SUBSAMPLE_NUM<-factor(T20_TRACKED_DATA$SUBSAMPLE_NUM)\n')
    CURRENT_SNP_R_SCRIPT.write(f'\n')
    CURRENT_SNP_R_SCRIPT.write(f'#open png\n')
    CURRENT_SNP_R_SCRIPT.write(f'png("{PATH_TO_MAIN}output_files/summary_plots/IDEA1/{cumulative_t20_dataframe_name}.png", bg = "white", width = 9.75, height = 3.25, units = "in", res = 1200, pointsize = 4)\n')
    CURRENT_SNP_R_SCRIPT.write(f'ggplot(T20_TRACKED_DATA, aes(x=PVAL_TYPE, y=VALUE, fill=SUBSAMPLE_NUM)) +\n')
    CURRENT_SNP_R_SCRIPT.write(f'   geom_boxplot() +\n')
    CURRENT_SNP_R_SCRIPT.write(f'   facet_wrap(~PVAL_TYPE, scale="free")\n')
    CURRENT_SNP_R_SCRIPT.write(f'\n')
    CURRENT_SNP_R_SCRIPT.write(f'dev.off()\n')
    CURRENT_SNP_R_SCRIPT.write(f'print("End of IDEA1 R script")\n')
    CURRENT_SNP_R_SCRIPT.write(f'# END OF R SCRIPT')
    CURRENT_SNP_R_SCRIPT.close()

# IDEA 1.5.1 (CHECKED!) fixed rscript run error
def IDEA_1_MAKE_BASH_SCRIPT(
        phenotype,
        cumulative_t20_dataframe_name
        ):
    print("Entered FUNCTION: IDEA_1_MAKE_BASH_SCRIPT",flush=True)
    #####################################
    ####### MAKING THE BATCH SCRIPT
    ######################################
    cumulative_t20_dataframe_name = cumulative_t20_dataframe_name.replace(".csv","")
    # location where the batch scripts will be written to
    CURRENT_SNP_BATCH=open(PATH_TO_MAIN+"batch_files/parallel_stage2/"+str(cumulative_t20_dataframe_name)+".sh","w")
    
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
    CURRENT_SNP_BATCH.write(f'Rscript {PATH_TO_MAIN}output_files/SNP_tracker_R_scripts/{cumulative_t20_dataframe_name}.R\n')
    CURRENT_SNP_BATCH.write(f'conda deactivate\n')
    CURRENT_SNP_BATCH.write(f'echo "END OF IDEA 1 batch script" \n')
    CURRENT_SNP_BATCH.write(f'# end of script')
    CURRENT_SNP_BATCH.close()

# IDEA 2.1.1 (CHECKED!)
def IDEA_2_CONTROL_CHECK(current_df,
            positive_control_df_name,
            negative_control_df_name,
            this_phenotype,
            GIFT_or_GWAS,
            subsample_number
            ):
    print("Entered FUNCTION: IDEA_2",flush=True)
    #print("State of GWAS or GIFT: ",GWAS_or_GIFT,flush=True)
    # list used for splitting 1 line from a GIFT csv into separate lines for each pval type (in the list)
    if this_phenotype=="Mo98":
        #MO98 data
        #MOT1 gene location boundaries 
        positive_control_chromosome = 2
        positive_control_LB = 10933005
        positive_control_UB = 10934604

        # rad50 gene location boundaries
        negative_control_chromosome = 2
        negative_control_LB =13600431
        negative_control_UB =13609104

    elif this_phenotype=="Na23":
        # Na23 stuff 
        # HKT1 gene location boundaries 
        positive_control_chromosome = 4
        positive_control_LB = 6391854
        positive_control_UB = 6395922

        # rad50 gene location boundaries
        negative_control_chromosome = 2
        negative_control_LB =13600431
        negative_control_UB =13609104

    if GIFT_or_GWAS == "GWAS":
        CHR = "chromosomes"
        POS = "positions"

    elif GIFT_or_GWAS == "GIFT": 
        CHR="CHROM"
        POS="POS"

    # copy current df into temp pos and neg dfs 
    temp_positive_control_df=current_df.copy()
    temp_negative_control_df=current_df.copy()

    # select rows that fit the positive or negative control regions
    temp_positive_control_df = temp_positive_control_df[(temp_positive_control_df[CHR]==positive_control_chromosome) & (temp_positive_control_df[POS]>= positive_control_LB) & (temp_positive_control_df[POS]<= positive_control_UB)]
    temp_negative_control_df = temp_negative_control_df[(temp_negative_control_df[CHR]==negative_control_chromosome) & (temp_negative_control_df[POS]>= negative_control_LB) & (temp_negative_control_df[POS]<= negative_control_UB)]

    if GIFT_or_GWAS == "GWAS":
            
        # changing the names of current columns
        temp_positive_control_df.rename(columns={CHR:'CHR',POS:'POS','pvals':'VALUE'}, inplace=True)
        temp_negative_control_df.rename(columns={CHR:'CHR',POS:'POS','pvals':'VALUE'}, inplace=True)
        
        # COLS : CHR, POS, VALUE

        #add column for pval type and subsample number
        temp_positive_control_df.insert(2,"PVAL_TYPE","GWAS_P")
        temp_negative_control_df.insert(2,"PVAL_TYPE","GWAS_P")
        # COLS : CHR, POS, PVAL_TYPE, VALUE
        
        temp_positive_control_df.insert(3,"SUBSAMPLE_NUM",subsample_number) 
        temp_negative_control_df.insert(3,"SUBSAMPLE_NUM",subsample_number) 
        # COLS : CHR, POS, PVAL_TYPE, SUBSAMPLE_NUM,VALUE

        # END COLS : CHR, POS, PVAL_TYPE, SUBSAMPLE_NUM, VALUE
    elif GIFT_or_GWAS == "GIFT": 

        #change column names
        temp_positive_control_df.rename(columns={CHR:'CHR',POS:'POS','absolute_theta':'ABS_THETA','pSNP4':'PSNP4','pSNP5':'PSNP5'}, inplace=True)
        temp_negative_control_df.rename(columns={CHR:'CHR',POS:'POS','absolute_theta':'ABS_THETA','pSNP4':'PSNP4','pSNP5':'PSNP5'}, inplace=True)

        # melt down the dataframe for pos and neg
        # should combine the column titles ABS_THETA, PSNP4 etc... into a column with their respective values under new column "VALUE"
        temp_negative_control_df=pandas.melt(temp_negative_control_df, id_vars=['CHR','POS'],var_name="PVAL_TYPE",value_name="VALUE")
        temp_positive_control_df=pandas.melt(temp_positive_control_df, id_vars=['CHR','POS'],var_name="PVAL_TYPE",value_name="VALUE")

        # COLS : CHR, POS, PVAL_TYPE, VALUE

        # add in the necessary columns for subsample number  
        temp_positive_control_df.insert(3,"SUBSAMPLE_NUM",subsample_number)
        temp_negative_control_df.insert(3,"SUBSAMPLE_NUM",subsample_number)

        # COLS : CHR, POS, PVAL_TYPE, SUBSAMPLE_NUM, VALUE

    #concatonate either or both dataframes if they arent empty
    if len(temp_positive_control_df)>0:
        
        #positive_control_df_name =  e.g. Mo98_GWAS_positive_control.csv
        try:
            positive_control_df=pandas.read_csv(PATH_TO_MAIN+"output_files/R_DATA/"+positive_control_df_name)

            print("Successfully read the negative control df",flush=True)

            positive_control_df=pandas.concat([positive_control_df,temp_positive_control_df],ignore_index=True)
            
        except:
            print("Tried to read positive control dataframe but it didnt exist...",flush=True)
            print("Setting current temp positive control to dataframe...",flush=True)
            positive_control_df=temp_positive_control_df
               

        positive_control_df.to_csv(PATH_TO_MAIN+"output_files/R_DATA/"+positive_control_df_name,header=True,index=False)

        del positive_control_df

    if len(temp_negative_control_df)>0:

        try:
            negative_control_df=pandas.read_csv(PATH_TO_MAIN+"output_files/R_DATA/"+negative_control_df_name)
            negative_control_df=pandas.concat([negative_control_df,temp_negative_control_df],ignore_index=True)
            negative_control_df=temp_negative_control_df
            print("Successfully read the negative control df",flush=True)

        except:
            print("Tried to read negative control dataframe but it didnt exist...",flush=True)
            print("Setting current temp positive control to dataframe...",flush=True)
            negative_control_df=temp_negative_control_df

        negative_control_df.to_csv(PATH_TO_MAIN+"output_files/R_DATA/"+negative_control_df_name,header=True,index=False)

        del negative_control_df
    
    print("IDEA_2 function done",flush=True)
    
    # END OF FUNCTION

# IDEA 2.3.1 (CHECKED!x2)
def IDEA_2_MAKE_R_AND_BASH_SCRIPT(
                        phenotype,
                        control_dataframe_name,
                        positive_or_negative
                        ):
    print("Entered FUNCTION: IDEA_2_MAKE_R_AND_BASH_SCRIPT",flush=True)

    #####################################
    ####### MAKING THE R SCRIPT
    ######################################
    control_dataframe_name=control_dataframe_name.replace(".csv","")

    # write the R script to pair with the SNP and its data
    CURRENT_SNP_R_SCRIPT=open(PATH_TO_MAIN+"output_files/SNP_tracker_R_scripts/"+control_dataframe_name+".R","w")
    CURRENT_SNP_R_SCRIPT.write(f'#R script for making box plots with ggplot\n')
    CURRENT_SNP_R_SCRIPT.write(f'library("tidyverse")\n')
    #CURRENT_SNP_R_SCRIPT.write(f'library("ggplot2")\n')
    CURRENT_SNP_R_SCRIPT.write(f'print("start of IDEA 2 R script")\n')
    CURRENT_SNP_R_SCRIPT.write(f'{control_dataframe_name}_data<- read.csv("{PATH_TO_MAIN}output_files/R_DATA/{control_dataframe_name}.csv", header= TRUE, sep=",")\n')
    CURRENT_SNP_R_SCRIPT.write(f'\n')
    CURRENT_SNP_R_SCRIPT.write(f'# Convert the subsample number to a FACTOR variable\n')
    CURRENT_SNP_R_SCRIPT.write(f'{control_dataframe_name}_data$SUBSAMPLE_NUM<-factor({control_dataframe_name}_data$SUBSAMPLE_NUM)\n')
    CURRENT_SNP_R_SCRIPT.write(f'\n')
    CURRENT_SNP_R_SCRIPT.write(f'#open png\n')
    CURRENT_SNP_R_SCRIPT.write(f'png("{PATH_TO_MAIN}output_files/summary_plots/IDEA2/{control_dataframe_name}.png", bg = "white", width = 9.75, height = 3.25, units = "in", res = 1200, pointsize = 4)\n')
    CURRENT_SNP_R_SCRIPT.write(f'ggplot({control_dataframe_name}_data, aes(x=PVAL_TYPE, y=VALUE, fill=SUBSAMPLE_NUM)) +\n')
    CURRENT_SNP_R_SCRIPT.write(f'   geom_boxplot() +\n')
    CURRENT_SNP_R_SCRIPT.write(f'   facet_wrap(~PVAL_TYPE, scale="free")\n')
    CURRENT_SNP_R_SCRIPT.write(f'dev.off()\n')
    CURRENT_SNP_R_SCRIPT.write(f'print("End of IDEA 2 R script")\n')
    CURRENT_SNP_R_SCRIPT.write(f'# END OF R SCRIPT')
    CURRENT_SNP_R_SCRIPT.close()

    #####################################
    ####### MAKING THE BATCH SCRIPT
    ######################################

    # location where the batch scripts will be written to
    CURRENT_SNP_BATCH=open(PATH_TO_MAIN+"batch_files/parallel_stage2/"+control_dataframe_name+".sh","w")
    
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
    CURRENT_SNP_BATCH.write(f'echo "start OF IDEA 2 batch script"\n')
    CURRENT_SNP_BATCH.write(f'#change to home directory\n')
    CURRENT_SNP_BATCH.write(f'cd /gpfs01/home/mbysh17\n')
    CURRENT_SNP_BATCH.write(f'# source conda environments\n')
    CURRENT_SNP_BATCH.write(f'source ~/.bashrc\n')
    CURRENT_SNP_BATCH.write(f'conda deactivate\n')
    CURRENT_SNP_BATCH.write(f'conda activate gift_env\n') # in gift_env but may change to R env if needed
    CURRENT_SNP_BATCH.write(f'# R SCRIPT FOR (IDEA 2) BOXPLOT\n')
    CURRENT_SNP_BATCH.write(f'Rscript {PATH_TO_MAIN}output_files/SNP_tracker_R_scripts/{control_dataframe_name}.R\n')
    CURRENT_SNP_BATCH.write(f'conda deactivate\n')
    CURRENT_SNP_BATCH.write(f'echo "END OF IDEA 2 batch script"\n')
    CURRENT_SNP_BATCH.write(f'# end of file\n')
    CURRENT_SNP_BATCH.close()

    # END OF FUNCTION



# initialise variables

Total_GIFT_Mo98 = 0
Total_GWAS_Mo98 = 0
Total_GIFT_Na23 = 0
Total_GWAS_Na23 = 0

csv_files=[]

#print("csv file list BEFORE sort: ======================================= \n",flush=True)
#print(csv_files,flush=True)
#print("\n ======================================= \n",flush=True)

for file in os.listdir(PATH_TO_MAIN+"output_files"):
    if file.endswith(".csv") and file.__contains__("T20")==False:
        #print(file,": ++++++++++++ ADDED ++++++++++++ !")
        csv_files.append(file)
    else:
        #print(file,": ////////// SKIPPED //////// !!")
        pass

# might need to turn for loops into a function soon (similar for loop used twice...)
#   could then MAYBE multithread the loops? depends...


csv_file_index =0
csv_file_index_max=int(len(csv_files)-1)

# PHENOTYPE_MEHTOD_SUBNUM_files=[]
Mo98_GIFT_200_files =[]
Mo98_GIFT_400_files =[]
Mo98_GIFT_600_files =[]
Mo98_GIFT_800_files =[]
Mo98_GIFT_1000_files =[]

Mo98_GWAS_200_files =[]
Mo98_GWAS_400_files =[]
Mo98_GWAS_600_files =[]
Mo98_GWAS_800_files =[]
Mo98_GWAS_1000_files =[]

Na23_GIFT_200_files =[]
Na23_GIFT_400_files =[]
Na23_GIFT_600_files =[]
Na23_GIFT_800_files =[]
Na23_GIFT_1000_files =[]

Na23_GWAS_200_files =[]
Na23_GWAS_400_files =[]
Na23_GWAS_600_files =[]
Na23_GWAS_800_files =[]
Na23_GWAS_1000_files =[]

# split csv files into categories (PHENO, METHOD, SUBNUM)
# PROGRESS METER
print("////////////////////////////////////////////////",flush=True)
print("Splitting the file names into groups...",flush=True)
print("////////////////////////////////////////////////",flush=True)
for csv_file in csv_files:
    csv_file_name = csv_file.split("_")
    #print(f"Current csv file to be split is... {csv_file}")
    #print(f"position 6 of its name is... {csv_file[6]}")
    # Mo98 GWAS lists
    if csv_file_name[2] == "Mo98" and csv_file_name[3]=='whole': 
        if csv_file_name[6]=="200":
            Mo98_GIFT_200_files.append(csv_file)
            #print(f"Appended {csv_file} to the Mo98_GIFT_200 files list!")

        elif csv_file_name[6]=="400":
            Mo98_GIFT_400_files.append(csv_file)

        elif csv_file_name[6]=="600":
            Mo98_GIFT_600_files.append(csv_file)

        elif csv_file_name[6]=="800":
            Mo98_GIFT_800_files.append(csv_file)

        elif csv_file_name[6]=="1000":
            Mo98_GIFT_1000_files.append(csv_file)

        # REMOVED THE !=T20 as there shouldnt be any T20 files in here anymore
    # Mo98 GWAS lists
    elif csv_file_name[2]=="Mo98" and csv_file_name[3]=='GWAS' :
        if csv_file_name[4]=="200":# GWAS code vvv  
            Mo98_GWAS_200_files.append(csv_file)

        elif csv_file_name[4]=="400":# GWAS code vvv  
            Mo98_GWAS_400_files.append(csv_file)  

        elif csv_file_name[4]=="600":# GWAS code vvv  
            Mo98_GWAS_600_files.append(csv_file)

        elif csv_file_name[4]=="800":# GWAS code vvv  
            Mo98_GWAS_800_files.append(csv_file) 

        elif csv_file_name[4]=="1000":# GWAS code vvv  
            Mo98_GWAS_1000_files.append(csv_file)

    # ---
    # Na23 GIFT lists
    elif csv_file_name[2] == "Na23" and csv_file_name[3]=='whole': 
        if csv_file_name[6]=="200":
            Na23_GIFT_200_files.append(csv_file)

        elif csv_file_name[6]=="400":
            Na23_GIFT_400_files.append(csv_file)

        elif csv_file_name[6]=="600":
            Na23_GIFT_600_files.append(csv_file)

        elif csv_file_name[6]=="800":
            Na23_GIFT_800_files.append(csv_file)

        elif csv_file_name[6]=="1000":
            Na23_GIFT_1000_files.append(csv_file)

    # Na23 GWAS lists
    elif csv_file_name[2]=="Na23" and csv_file_name[3]=='GWAS':
        if csv_file_name[4]=="200":# GWAS code vvv  
            Na23_GWAS_200_files.append(csv_file)

        elif csv_file_name[4]=="400":# GWAS code vvv  
            Na23_GWAS_400_files.append(csv_file) 

        elif csv_file_name[4]=="600":# GWAS code vvv  
            Na23_GWAS_600_files.append(csv_file)  

        elif csv_file_name[4]=="800":# GWAS code vvv  
            Na23_GWAS_800_files.append(csv_file) 

        elif csv_file_name[4]=="1000":# GWAS code vvv  
            Na23_GWAS_1000_files.append(csv_file)  

# PROGRESS METER
print("////////////////////////////////////////////////",flush=True)
print("Finished splitting into groups!",flush=True)
print("////////////////////////////////////////////////",flush=True)

# PROGRESS METER
print("////////////////////////////////////////////////",flush=True)
print("TEST Na23 GWAS 200!"+Na23_GWAS_200_files[1],flush=True)
print("////////////////////////////////////////////////",flush=True)

# PROGRESS METER
print("////////////////////////////////////////////////",flush=True)
print("TEST Na23 GWAS 200!"+Na23_GWAS_400_files[1],flush=True)
print("////////////////////////////////////////////////",flush=True)

# PROGRESS METER
print("////////////////////////////////////////////////",flush=True)
print("TEST Na23 GWAS 200!"+Na23_GWAS_600_files[1],flush=True)
print("////////////////////////////////////////////////",flush=True)

# PROGRESS METER
print("////////////////////////////////////////////////",flush=True)
print("TEST Na23 GWAS 200!"+Na23_GWAS_800_files[1],flush=True)
print("////////////////////////////////////////////////",flush=True)

# PROGRESS METER
print("////////////////////////////////////////////////",flush=True)
print("TEST Na23 GWAS 200!"+Na23_GWAS_1000_files[1],flush=True)
print("////////////////////////////////////////////////",flush=True)


# concatonate each into one big file -> x20 csv files for 2 phenotypes (10 per phenotype)

######################################################
####################################
### IDEA 3
# store each of the lists of the files into one big list to loop through
list_of_list_of_files=[Mo98_GIFT_200_files,
                       Mo98_GIFT_400_files,
                       Mo98_GIFT_600_files,
                       Mo98_GIFT_800_files,
                       Mo98_GIFT_1000_files,
                       Mo98_GWAS_200_files,
                       Mo98_GWAS_400_files,
                       Mo98_GWAS_600_files,
                       Mo98_GWAS_800_files,
                       Mo98_GWAS_1000_files,
                       Na23_GIFT_200_files,
                       Na23_GIFT_400_files,
                       Na23_GIFT_600_files,
                       Na23_GIFT_800_files,
                       Na23_GIFT_1000_files,
                       Na23_GWAS_200_files,
                       Na23_GWAS_400_files,
                       Na23_GWAS_600_files,
                       Na23_GWAS_800_files,
                       Na23_GWAS_1000_files]

### #######################
### LOOP 1 (IDEA 1 AND 2)
###
# Reminder of CSV format (GIFT) NAME leaf_ionome_Mo98_whole_genome_metrics_600_732692.csv
# CHROM,POS,largest_theta,smallest_theta,absolute_theta,theta_range,largest_relative_theta,smallest_relative_theta,absolute_relative_theta,range_relative_theta,min_p,mean_p,log_mean_p,bigest_theta_p,pSNP4,pSNP5
# 1,73,6.285714285714285,-4.224489795918359,6.285714285714285,10.510204081632644,0.3845193508114856,-0.25842696629213435,0.3845193508114856,0.6429463171036199,4.185363300872138e-05,nan,nan,0.00015134587316541535,0.15259768149369662,0.6580333017260325

# reminder of csv format (GWAS) NAME leaf_ionome_Mo98_GWAS_600_732692.csv
# chromosomes,positions,pvals,mafs,macs,genotype_var_perc
# 1,55,0.621946624343,0.0016694490818,1,0.000407516956329
for list_of_files in list_of_list_of_files:

    # read first item to determine what it is
    first_item_name=list_of_files[1].split("_")

    if first_item_name[3]=="whole":
        GIFT_or_GWAS="GIFT"
        subsample_number = first_item_name[6]

    elif first_item_name[3]=="GWAS":
        GIFT_or_GWAS="GWAS"
        subsample_number = first_item_name[4]
   
   #first_item_name[2] == phenotype

    # PROGRESS METER
    print("////////////////////////////////////////////////",flush=True)
    print(f"Concatonating {first_item_name[2]}_{GIFT_or_GWAS}_{subsample_number}",flush=True)
    print("////////////////////////////////////////////////",flush=True)

    if GIFT_or_GWAS=="GIFT":

        # concatonate all the files in the given list
        this_df=pandas.concat([pandas.read_csv(PATH_TO_MAIN+"output_files/"+csv_file,sep=",",na_filter=False,usecols=["CHROM","POS","absolute_theta","pSNP4","pSNP5"]) for csv_file in list_of_files],ignore_index=True)

        # once concatonated, check for pos + neg control values in this function
        IDEA_2_CONTROL_CHECK(this_df,
                             str(first_item_name[2])+"_positive_control.csv", # positive control df
                             str(first_item_name[2])+"_negative_control.csv", # negative control df
                             first_item_name[2], # phenotype
                             GIFT_or_GWAS, # GIFT or gwas
                             subsample_number # subsample num
                             )

        # then process the dataframe to combine locations and calculate averages
        this_df=IDEA_3_process_all_snps_file(this_df, GIFT_or_GWAS,len(list_of_files),subsample_number)

    elif GIFT_or_GWAS=="GWAS":

        this_df=pandas.concat([pandas.read_csv(PATH_TO_MAIN+"output_files/"+csv_file,sep=",",na_filter=False,usecols=["chromosomes","positions","pvals"]) for csv_file in list_of_files],ignore_index=True)

        # IDEA 2 function
        IDEA_2_CONTROL_CHECK(this_df,
                             str(first_item_name[2])+"_positive_control.csv", # positive control df
                             str(first_item_name[2])+"_negative_control.csv", # negative control df
                             first_item_name[2], # phenotype
                             GIFT_or_GWAS, # GIFT or gwas
                             subsample_number # subsample num
                             )

        this_df=IDEA_3_process_all_snps_file(this_df, GIFT_or_GWAS,len(list_of_files),subsample_number)

    # save to csv
    this_df.to_csv(PATH_TO_MAIN+"output_files/R_DATA/"+first_item_name[2]+"_"+GIFT_or_GWAS+"_"+subsample_number+"_ALL.csv",header=True,index=False)

    # delete df to save mem
    del this_df

##########
#########################
#########################################

####################################################
# IDEA 2.2 #####################################
############################################
IDEA_2_MAKE_R_AND_BASH_SCRIPT("Mo98", "Mo98_positive_control.csv","positive")
IDEA_2_MAKE_R_AND_BASH_SCRIPT("Mo98", "Mo98_negative_control.csv","negative")
IDEA_2_MAKE_R_AND_BASH_SCRIPT("Na23", "Na23_positive_control.csv","positive")
IDEA_2_MAKE_R_AND_BASH_SCRIPT("Na23", "Na23_negative_control.csv","negative")

print("IDEA 2 FINISHED",flush=True)
############################################
# IDEA 2.2 #####################################
####################################################


####################################################
# IDEA 3.4 #####################################
############################################
subsample_num_list=[200,400,600,800,1000] # can later update this to read from earlier scripts or something

phenotype_list=["Mo98","Na23"] # can later update this to read from the phenotype text file

pvals=["AVERAGE_P","AVERAGE_PSNP4","AVERAGE_PSNP5","AVERAGE_ABS_THETA"] # these should always be the same 4 types
# R script creation and running
for phenotype in phenotype_list:

    for subsample_number in subsample_num_list:

        for pval_type in pvals:

            # Example variable inputs are as follows: 
            # example 1 Mo98,200,AVERAGE_P
            # example 2 Mo98,200,AVERAGE_PSNP4
            # ...
            # example x Mo98,400,AVERAGE_P
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

                # get T20 SNPs locations for GWAS_P for 1000 subsample -> RETURNS THE NAME
                Mo98_GWAS_P_T20_SNP_locations_path = GET_T20_LOCATIONS_AT_1000("Mo98_GWAS_1000_ALL.csv",pval_type,phenotype) 

            elif pval_type == "AVERAGE_PSNP4":

                Mo98_PSNP4_T20_SNP_locations_path = GET_T20_LOCATIONS_AT_1000("Mo98_GIFT_1000_ALL.csv",pval_type,phenotype) 

            elif pval_type == "AVERAGE_PSNP5":

                Mo98_PSNP5_T20_SNP_locations_path = GET_T20_LOCATIONS_AT_1000("Mo98_GIFT_1000_ALL.csv",pval_type,phenotype) 

            elif pval_type == "AVERAGE_ABS_THETA":

                Mo98_ABS_THETA_T20_SNP_locations_path = GET_T20_LOCATIONS_AT_1000("Mo98_GIFT_1000_ALL.csv",pval_type,phenotype) 

        elif phenotype == "Na23":
            
            if pval_type == "AVERAGE_P":

                # get T20 SNPs locations for GWAS_P for 1000 subsample
                Na23_GWAS_P_T20_SNP_locations_path = GET_T20_LOCATIONS_AT_1000("Na23_GWAS_1000_ALL.csv",pval_type,phenotype) 

            elif pval_type == "AVERAGE_PSNP4":

                Na23_PSNP4_T20_SNP_locations_path = GET_T20_LOCATIONS_AT_1000("Na23_GIFT_1000_ALL.csv",pval_type,phenotype) 

            elif pval_type == "AVERAGE_PSNP5":

                Na23_PSNP5_T20_SNP_locations_path = GET_T20_LOCATIONS_AT_1000("Na23_GIFT_1000_ALL.csv",pval_type,phenotype) 

            elif pval_type == "AVERAGE_ABS_THETA":
                
                Na23_ABS_THETA_T20_SNP_locations_path = GET_T20_LOCATIONS_AT_1000("Na23_GIFT_1000_ALL.csv",pval_type,phenotype) 

####################################################
# IDEA 1.2 ######################################

# set up cumulative T20 dataframe paths for each phenotype

Mo98_cumulative_t20_dataframe_path = PATH_TO_MAIN+"output_files/R_DATA/Mo98_cumulative_t20_dataframe.csv"

Na23_cumulative_t20_dataframe_path = PATH_TO_MAIN+"output_files/R_DATA/Na23_cumulative_t20_dataframe.csv"


####################################################
# IDEA 1.3.X ######################################
csv_file_index=0


#1) Concat all the files (same way as before) but WITHOUT the process_files funciton 
#           (as dont want averages -> want raw values)

for list_of_files in list_of_list_of_files:
    # read first item to determine what it is
    first_item_name=list_of_files[1].split("_")
    if first_item_name[3]=="whole":
        GIFT_or_GWAS="GIFT"
        subsample_number = first_item_name[6]

    elif first_item_name[3]=="GWAS":
        GIFT_or_GWAS="GWAS"
        subsample_number = first_item_name[4]
   
    this_phenotype=first_item_name[2]

    # PROGRESS METER
    print("////////////////////////////////////////////////",flush=True)
    print(f"Concatonating {first_item_name[2]}_{GIFT_or_GWAS}_{subsample_number}",flush=True)
    print("////////////////////////////////////////////////",flush=True)

    # concat files based on gift or gwas
    if GIFT_or_GWAS=="GIFT":

        this_df=pandas.concat([pandas.read_csv(PATH_TO_MAIN+"output_files/"+csv_file,sep=",",na_filter=False,usecols=["CHROM","POS","absolute_theta","pSNP4","pSNP5"]) for csv_file in list_of_files],ignore_index=True)

    elif GIFT_or_GWAS=="GWAS":

        this_df=pandas.concat([pandas.read_csv(PATH_TO_MAIN+"output_files/"+csv_file,sep=",",na_filter=False,usecols=["chromosomes","positions","pvals"]) for csv_file in list_of_files],ignore_index=True)

    # Testing
    print("////////////////////////////////////////////////",flush=True)
    print(f"Head of current df",flush=True)
    print(this_df.head(),flush=True)
    print("////////////////////////////////////////////////",flush=True)



    # accumulate snps based on phenotype and gift or gwas
    if this_phenotype=="Mo98":
        IDEA_1_ACCUMULATE_T20_SNP_DATA(this_df,
                                        Mo98_GWAS_P_T20_SNP_locations_path,
                                        Mo98_PSNP4_T20_SNP_locations_path,
                                        Mo98_PSNP5_T20_SNP_locations_path,
                                        Mo98_ABS_THETA_T20_SNP_locations_path,
                                        Mo98_cumulative_t20_dataframe_path, 
                                        GIFT_or_GWAS,
                                        subsample_number
                                        )
            
    elif this_phenotype=="Na23":
        IDEA_1_ACCUMULATE_T20_SNP_DATA(this_df,
                                        Na23_GWAS_P_T20_SNP_locations_path,
                                        Na23_PSNP4_T20_SNP_locations_path,
                                        Na23_PSNP5_T20_SNP_locations_path,
                                        Na23_ABS_THETA_T20_SNP_locations_path,
                                        Na23_cumulative_t20_dataframe_path, 
                                        GIFT_or_GWAS,
                                        subsample_number
                                    )

    # save to csv (disabled for now- dont think i need to store it)
    #this_df.to_csv(PATH_TO_MAIN+"output_files/R_DATA/"+first_item_name[2]+"_"+GIFT_or_GWAS+"_"+subsample_number+"_IDEA_1.csv",header=True,index=False)

    # delete df to save mem
    del this_df

####################################################
# IDEA 1.4 ######################################
IDEA_1_MAKE_R_SCRIPT("Mo98","Mo98_cumulative_t20_dataframe.csv")
IDEA_1_MAKE_R_SCRIPT("Na23","Na23_cumulative_t20_dataframe.csv")

####################################################
# IDEA 1.5 ######################################
IDEA_1_MAKE_BASH_SCRIPT("Mo98","Mo98_cumulative_t20_dataframe.csv")
IDEA_1_MAKE_BASH_SCRIPT("Na23","Na23_cumulative_t20_dataframe.csv")

print("END OF IDEA 1",flush=True)

print("END OF SNP_TRACKER.PY",flush=True)
# end of file