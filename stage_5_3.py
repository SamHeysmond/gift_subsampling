# packages

import pandas
import os

# change these functions to import from a main py file later
def fetch_subsample_numbers_list(subsample_numbers_list_file):

    subsample_num_list=[]

    temp_file=open(f"{subsample_numbers_list_file}","r")

    for line in temp_file:

        subsample_number=line.replace('\n','')

        subsample_num_list.append(subsample_number)
    
    print("Subsample number list fetched: ",flush=True)
    print(subsample_num_list,flush=True)

    return subsample_num_list

def fetch_phenotype_list(phenotype_list_file):
    phenotypes_list=[]
    temp_file=open(f"{phenotype_list_file}","r")

    for line in temp_file:

        subsample_number=line.replace('\n','')

        phenotypes_list.append(subsample_number)
    
    print("Phenotype list fetched: ",flush=True)
    print(phenotypes_list,flush=True)

    return phenotypes_list

def write_R_script_and_shell(subsample_number,position_info):
    # saved to batch_files/stage_5_prerun
    R_out = open(f"{PATH_TO_MAIN}batch_files/stage_5_3_Rscripts/{subsample_number}_{position_info}.R","w")
    R_out.write(f'# start of file \n')
    R_out.write(f'library(ggplot2)\n')
    R_out.write(f'library(tidyverse)\n')
    R_out.write(f' \n')
    R_out.write(f' \n')
    R_out.write(f'ThPaths_Data<- read.csv("{PATH_TO_MAIN}/core_files/combined_position_data/{subsample_num}_{position_info}.csv")\n')
    
    # R_out.write(f'ThPaths_Data$AverageValues<-rowMeans(ThPaths_Data)\n')
    # R_out.write(f'ThPaths_Data$AverageValues<-ThPaths_Data$AverageValues/{subsample_number}\n')
    # R_out.write(f'ThPaths_Data$Row_SD <- apply(ThPaths_Data[,-ncol(ThPaths_Data)], 1, sd)\n')
    # R_out.write(f'ThPaths_Data$Row_SD <- (ThPaths_Data$Row_SD)/{subsample_number}\n')

    # R_out.write(f'ThPaths_Data$Index <- seq_len(nrow(ThPaths_Data)) \n')
    # R_out.write(f'ThPaths_Data$Index <- ThPaths_Data$Index/{subsample_number} \n')

    # R_out.write(f'ThPaths_Data$Index <- ThPaths_Data$Index/{subsample_number}\n')

    # R_out.write(f'ThPaths_Data_long <- ThPaths_Data %>%\n')
    # R_out.write(f'  pivot_longer(cols = -Index,names_to = "Line", values_to = "Value") \n')
    # R_out.write(f' \n')

    # R_out.write(f'png("{PATH_TO_MAIN}output_files/stage_5_results/{subsample_number}_{position_info}.png") \n')
    # R_out.write(f'ggplot(ThPaths_Data_long, aes(x=Index/{subsample_number}, y=Value/{subsample_number}, color=Line)) +\n')
    # R_out.write(f'  geom_line(alpha=0.5) +\n')
    # R_out.write(f'  geom_point() +\n')
    # R_out.write(f'  theme(legend.position = "none") \n')

    #### SECOND PLOT
    # R_out.write(f'ggplot(ThPaths_Data, aes(x=Index/{subsample_number}, y=AverageValues/{subsample_number}, color=Line)) +\n')
    # R_out.write(f'  geom_line(alpha=0.5) +\n')
    # R_out.write(f'  geom_point() +\n')
    # R_out.write(f'  theme(legend.position = "none") \n')

    ### THIRD PLOT
    # R_out.write(f'ggplot(ThPaths_Data, aes(x=Index, y=AverageValues)) +\n')
    # R_out.write(f'  geom_line() +\n')
    # R_out.write(f'  geom_point() +\n')
    # R_out.write(f'  geom_errorbar(aes(ymin=AverageValues-Row_SD, ymax = AverageValues+Row_SD, color="red"))\n')
    # R_out.write(f'   \n')

    ### FOURTH PLOT
    R_out.write(f'ThPaths_Data$Index <- seq_len(nrow(ThPaths_Data)) \n')
    R_out.write(f'ThPaths_Data$Index <- ThPaths_Data$Index/{subsample_number} \n')
    
    R_out.write(f'ThPaths_Data_long <- ThPaths_Data %>%\n')
    R_out.write(f'  pivot_longer(cols = -Index, names_to = "Line", values_to = "Value") \n')
    R_out.write(f' \n')
    if int(subsample_num)==999:
        R_out.write(f'ThPaths_Data <- ThPaths_Data %>% \n')
        R_out.write(f'    pivot_longer(cols = -c(Index),names_to = "Line", values_to = "Value") \n')
    elif int(subsample_num)!=999 and Matryoshka==True:
        R_out.write(f'ThPaths_Data <- ThPaths_Data %>% \n')
        R_out.write(f'    pivot_longer(cols = -c(Index),names_to = "Line", values_to = "Value") \n')
    elif int(subsample_num)!=999 and Matryoshka==False:
        R_out.write(f'ThPaths_Data$AverageValues<-apply(ThPaths_Data[, !(names(ThPaths_Data) %in% c("AverageValues","Index"))], 1, mean) \n')
        R_out.write(f'ThPaths_Data$Row_SD <- apply(ThPaths_Data[, !(names(ThPaths_Data) %in% c("AverageValues","Index"))], 1, sd)\n')
        R_out.write(f'ThPaths_Data <- ThPaths_Data %>% \n')
        R_out.write(f'    pivot_longer(cols = -c(Index,AverageValues, Row_SD),names_to = "Line", values_to = "Value") \n')


    R_out.write(f'png("{PATH_TO_MAIN}output_files/stage_5_results/{subsample_number}_{position_info}.png",res=75, height=900, width = 900) \n')

    if int(subsample_num)==999:
        R_out.write(f'ggplot(ThPaths_Data, aes(x=Index/{subsample_number} , y=Value/{subsample_number} , color="black")) +\n')
        R_out.write(f'  geom_line(color="black")+\n')
        R_out.write(f'   geom_point(color="black")+\n')
        R_out.write(f'   theme(axis.text=element_text(size=30),\n')
        R_out.write(f'      axis.title=element_text(size=30,face="bold"),\n')
        if position_info in list_of_currently_significant_snps:
            R_out.write(f'      plot.background = element_rect(fill = alpha("green",0.5)))\n')
        else:
            R_out.write(f'      plot.background = element_rect(fill = alpha("red",0.5)))\n')
        R_out.write(f'   \n')
        R_out.write(f'   \n')
    
    elif int(subsample_num)!=999 and Matryoshka==True:
        R_out.write(f'ggplot(ThPaths_Data, aes(x=Index/{subsample_number} , y=Value/{subsample_number} , color="black")) +\n')
        R_out.write(f'  geom_line(color="black")+\n')
        R_out.write(f'   geom_point(color="black")+\n')
        R_out.write(f'   theme(axis.text=element_text(size=30),\n')
        R_out.write(f'      axis.title=element_text(size=30,face="bold"),\n')
        if position_info in list_of_currently_significant_snps:
            R_out.write(f'      plot.background = element_rect(fill = alpha("green",0.5)))\n')
        else:
            R_out.write(f'      plot.background = element_rect(fill = alpha("red",0.5)))\n')
        R_out.write(f'   \n')
        R_out.write(f'   \n')
        
    elif int(subsample_num)!=999 and Matryoshka==False:

        R_out.write(f'ggplot(ThPaths_Data_long, aes(x=Index/{subsample_number}, y=Value/{subsample_number}, color=Line)) +\n')
        R_out.write(f'  geom_line() +\n')
        R_out.write(f'  geom_point() +\n')
        R_out.write(f'  theme(legend.position = "none") +\n')
        R_out.write(f'  geom_errorbar(data=ThPaths_Data,aes(x=Index/{subsample_number}, ymin=(AverageValues-Row_SD)/{subsample_number}, ymax = (AverageValues+Row_SD)/{subsample_number}), color="black")+\n')
        R_out.write(f'  geom_line(data=ThPaths_Data, aes(x=Index/{subsample_number},y=AverageValues/{subsample_number}),color="red") +\n')
        if int(subsample_num)>=400:
            R_out.write(f'  geom_point(data=ThPaths_Data,aes(x=Index/{subsample_number},y=AverageValues/{subsample_number}),color="red") +\n')
        else:
            R_out.write(f'  geom_point(data=ThPaths_Data,aes(x=Index/{subsample_number},y=AverageValues/{subsample_number}),color="black") +\n')
        R_out.write(f'  theme(axis.text=element_text(size=30), \n')
        R_out.write(f'      axis.title=element_text(size=30,face="bold"),\n')
        
        
        significant_snp_data=pandas.read_csv(f"{PATH_TO_MAIN}output_files/significant_SNP_data/GIFT_SIGNIFICANT_{subsample_number}.csv")

        # obtain list of significant SNPs at this subsample level
        list_of_currently_significant_snps=list(significant_snp_data['POSITION_DATA'])

        # check if position info is in the currently significant SNPs for GIFT 
        print(f"Current position info to check for: {position_info}",flush=True)
        if position_info in list_of_currently_significant_snps:
            R_out.write(f'      plot.background = element_rect(fill = alpha("green",0.5)))\n')
        else:
            R_out.write(f'      plot.background = element_rect(fill = alpha("red",0.5)))\n')
            
        R_out.write(f'   \n')
        R_out.write(f'   \n')
    # R_out.write(f'ggplot(ThPaths_Data, aes(x=Index, y=AverageValues)) +  \n')
    # R_out.write(f'  geom_line() +\n')
    # R_out.write(f'  geom_point() +\n')
    # R_out.write(f'  geom_errorbar(aes(ymin=AverageValues-Row_SD, ymax = AverageValues+Row_SD, color="red"))\n')
    # R_out.write(f'   \n')
    R_out.write(f'dev.off() \n')
    R_out.write(f' \n')
    R_out.write(f'# end of script\n')
    R_out.close()

    shell_out=open(f"{PATH_TO_MAIN}batch_files/stage_5_3_prerun/{subsample_number}_{position_info}.sh","w")
    shell_out.write(f'#!/bin/bash\n')
    shell_out.write(f'#SBATCH --partition=defq\n')
    shell_out.write(f'#SBATCH --nodes=1\n')
    shell_out.write(f'#SBATCH --ntasks=1\n')
    shell_out.write(f'#SBATCH --cpus-per-task=2\n')
    shell_out.write(f'#SBATCH --mem=6g\n')
    shell_out.write(f'#SBATCH --time=01:00:00\n')
    shell_out.write(f'#SBATCH --job-name=S5_3_batch\n')
    shell_out.write(f'#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out\n')
    shell_out.write(f'#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err\n')
    shell_out.write(f'#SBATCH --mail-type=ALL\n')
    shell_out.write(f'#SBATCH --mail-user=mbysh17@nottingham.ac.uk\n')
    shell_out.write(f'#===============================\n')
    shell_out.write(f'echo "start OF IDEA 5_3 subrun"\n')
    shell_out.write(f'#change to home directory\n')
    shell_out.write(f'cd /gpfs01/home/mbysh17\n')
    shell_out.write(f'# source conda environments\n')
    shell_out.write(f'source ~/.bashrc\n')
    shell_out.write(f'conda deactivate\n')
    shell_out.write(f'conda activate r_env\n') # in gift_env but may change to R env if needed
    shell_out.write(f'Rscript {PATH_TO_MAIN}batch_files/stage_5_3_Rscripts/{subsample_number}_{position_info}.R\n')
    shell_out.write(f' \n')
    shell_out.write(f' \n')
    shell_out.write(f'# end of file \n')
    shell_out.close()

# set constant for file path
PATH_TO_MAIN = "/gpfs01/home/mbysh17/"

# temp for matryoshka
Matryoshka=True

subsample_numbers_list = fetch_subsample_numbers_list(PATH_TO_MAIN+"core_files/subsample_numbers_list.txt")
phenotype_list=fetch_phenotype_list(PATH_TO_MAIN+"core_files/phenotypes_list.txt")

for phenotype in phenotype_list:
    # Mo98

    # get positions to keep for this pehnotype
    positions_df=pandas.read_csv(f'{PATH_TO_MAIN}core_files/{phenotype}_position_data_keep.csv')
    positions_list=[]
    positions_list=list(positions_df['POSITION_DATA'])
    print("positions list: ",flush=True)
    print(positions_list,flush=True)

    for subsample_num in subsample_numbers_list:
        # 50, 100, 200, 400, 600, 999

        current_list_of_files = []
        # gather all names of files that exist for this subsample number (should be 100 each except max subsample)
        for file in os.listdir(PATH_TO_MAIN+"core_files/theta_paths"):

            if file.endswith(".csv") and file.__contains__(f"_{subsample_num}_")==True:
                # adds all the csv files needed
                current_list_of_files.append(file)
                # name format: delta_theta_paths_400_1620462.csv
            else:
                pass
        
        print("current list of files: ",flush=True)
        print(*current_list_of_files,sep='\n')
        
        for position in positions_list:
            # make empty dataframe for the current position
            current_position_df=pandas.DataFrame()

            # open each file one by one and make new column using the ID as a header and ..
                # naming the dataframe by the position as CHROM_POS.csv format
            for current_dataset in current_list_of_files:
                
                temp_df=pandas.read_csv(f'{PATH_TO_MAIN}core_files/theta_paths/{current_dataset}')

                # print("Head of temp_df",flush=True)
                # print(temp_df.head(10))

                dataset_name = current_dataset.split("_")
                current_ID = dataset_name[4].replace(".csv","") 
                # print(f"current ID: {current_ID}",flush=True)


                # make new column with ID as the header (future update remove the .wp or something...)
                current_position_df[f'{current_ID}'] = temp_df[f'{position}.Wp{current_ID}']

                del temp_df
                        
            # write the position dataframe to file
            reformatted_position=position.replace(":","_")
            current_position_df.to_csv(f'{PATH_TO_MAIN}core_files/combined_position_data/{subsample_num}_{reformatted_position}.csv',header=True,index=False)
        

            write_R_script_and_shell(subsample_num,reformatted_position)
            
            del current_position_df

print("End of stage_5_3.py",flush=True)

# end of script