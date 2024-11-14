# to do notes
    # make some global variables for repeated variables inside functions?

#packages
import argparse

PATH_TO_MAIN = "/gpfs01/home/mbysh17/"


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

# IDEA 3 RSCRIPT UPDATED (THRESHOLDS - testing)
def IDEA_3_R_AND_BATCH(phenotype,subsample_number,pval_type):
    print("Entered FUNCTION: IDEA_3_R_AND_BATCH",flush=True)
    #########################################
    ### MAKE R SCRIPT ##################
    #############################

    # make R script for each P value type (pSNP4, pSNP5, abs theta, GWAS_P)
    R_out=open(PATH_TO_MAIN+"output_files/SNP_tracker_R_scripts/"+str(phenotype)+"_"+str(subsample_number)+"_"+str(pval_type)+"_MANHATTAN.R","w")
    R_out.write(f'#R script for making average manhattan plots with ggplot\n')
    R_out.write(f'library("tidyverse")\n')
    # R_out.write(f'library("ggrepel")\n')   # Dont need ggrepel for now
    R_out.write(f'print("Start of IDEA 3 R")\n')
    R_out.write(f'\n')

    #########################################
    ### need to repeat the below steps to do a before (uncorrected) and after (_all.csv) so that i can compare the graphs and use the _all for further analysis only ##################
    #############################

    # for GWAS data....
    if pval_type=="AVERAGE_P":

        # fetch the data of the csv for the current phenotype and current method (GWAS)
        R_out.write(f'csv_data<-read.csv("{PATH_TO_MAIN}output_files/R_DATA/{phenotype}_GWAS_{subsample_number}_ALL.csv",header=TRUE)\n')

    # for GIFT data
    else:

        # fetch the data of the csv for the current phenotype and current method (GIFT)
        R_out.write(f'csv_data<-read.csv("{PATH_TO_MAIN}output_files/R_DATA/{phenotype}_GIFT_{subsample_number}_ALL.csv",header=TRUE)\n')  

    # cumulative calculations
    R_out.write(f'mydata <- csv_data %>%\n')
    R_out.write(f'     # Compute CHR size\n')
    R_out.write(f'     group_by(CHR) %>% \n')
    R_out.write(f'     summarise(chr_len=max(POS)) %>%\n')
    R_out.write(f'     # Calculate cumulative position of each CHR\n')
    R_out.write(f'     mutate(tot=cumsum(chr_len)-chr_len) %>%\n')
    R_out.write(f'     select(-chr_len) %>%\n')
    R_out.write(f'     # Add this info to the initial dataset\n')
    R_out.write(f'     left_join(csv_data, ., by=c("CHR"="CHR")) %>%\n')
    R_out.write(f'     # Add a cumulative position of each SNP\n')
    R_out.write(f'     arrange(CHR, POS) %>%\n')
    R_out.write(f'     mutate( BPcum=POS+tot) \n')
    R_out.write(f'axisdf = mydata %>%\n')
    R_out.write(f'     group_by(CHR) %>%\n')
    R_out.write(f'     summarize(center=( max(BPcum) + min(BPcum) ) / 2 )\n')

    # set y limit for the graph (not sure if this changes anything here though)
    # y limit depends on the pval type in question
    R_out.write(f'ylim <- abs(floor(log10(min(mydata)))) +1\n')

    # import and calculate with threshold data
    R_out.write(f'threshold_data_csv<-read.csv("{PATH_TO_MAIN}output_files/R_DATA/THRESHOLDS.csv", header= TRUE, sep=",")\n')
    R_out.write(f'# subset for the current phenotype and pval type\n')

    if pval_type=="AVERAGE_P":
            # GWAS THRESHOLDS IMPORTED
        R_out.write(f'NN_THRESHOLD_DATA<-subset(threshold_data_csv,PHENOTYPE=="{phenotype}" & PVAL_TYPE=="{pval_type}" & SUBSAMPLE_NUM=={subsample_number} & THRESHOLD_TYPE=="NNBF")\n')
        R_out.write(f'NF_THRESHOLD_DATA<-subset(threshold_data_csv,PHENOTYPE=="{phenotype}" & PVAL_TYPE=="{pval_type}" & SUBSAMPLE_NUM=={subsample_number} & THRESHOLD_TYPE=="NFBF")\n')

    else:
            # GIFT THRESHOLDS IMPORTED
        R_out.write(f'NN_THRESHOLD_DATA<-subset(threshold_data_csv,PHENOTYPE=="{phenotype}" & PVAL_TYPE=="{pval_type}" & SUBSAMPLE_NUM=={subsample_number} & THRESHOLD_TYPE=="NNNH")\n')
        R_out.write(f'NF_THRESHOLD_DATA<-subset(threshold_data_csv,PHENOTYPE=="{phenotype}" & PVAL_TYPE=="{pval_type}" & SUBSAMPLE_NUM=={subsample_number} & THRESHOLD_TYPE=="NFNH")\n')

    R_out.write(f'#subset for each type of threshold\n')
    R_out.write(f'NN_THRESHOLD<-NN_THRESHOLD_DATA$THRESHOLD_VALUE\n')
    R_out.write(f'NF_THRESHOLD<-NF_THRESHOLD_DATA$THRESHOLD_VALUE\n')
    R_out.write(f'print("NN_THRESHOLD value")\n')
    R_out.write(f'print(NN_THRESHOLD)\n')
    R_out.write(f'print("NF_THRESHOLD value")\n')
    R_out.write(f'print(NF_THRESHOLD)\n')
    R_out.write(f'\n')


    # R_out.write('if (BY_THRESHOLD == 0){\n')
    # R_out.write(f'BY_THRESHOLD=max(-log10(csv_data${pval_type}))\n')
    # R_out.write('}\n')

        ## look into how many points are above the threshold
    R_out.write(f'#calculate amount of points above the 95th threshold as a percentage of all points\n')
    R_out.write(f'percent_sig<-length(which(-log10(csv_data${pval_type})>NF_THRESHOLD))/length(csv_data${pval_type})\n')
    R_out.write(f'percent_sig<-percent_sig*100\n')
        # new (changed rounding to 4 sf or dp? not sure which)
    R_out.write(f'percent_sig<-round(percent_sig,4)\n')

    R_out.write(f'#open png\n')
    # Send output to IDEA 3 summary plot folder
    R_out.write(f'png("{PATH_TO_MAIN}output_files/summary_plots/IDEA3/{phenotype}_{subsample_number}_{pval_type}_MANHATTAN.png", bg = "white", width = 7.5, height = 2.5, units = "in", res = 1200, pointsize = 3)\n')

    # make the plot
        # this is where the plot depends on the pval type. 
    R_out.write(f'ggplot(mydata, aes(x=BPcum, y=(-log10({pval_type})), color = as_factor(CHR))) +\n')
    R_out.write(f'     # Show all points\n')
    R_out.write(f'     geom_point(alpha=0.5) +\n')
    R_out.write(f'     # custom X axis:\n')
    R_out.write(f'     scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +\n')
    # TEMP PAUSE FOR TESTING
    #Rscript_output.write(f'     scale_y_continuous(expand = c(0, 0) ) + # remove space between plot area and x axis\n')
    # label axis
    R_out.write(f'     # my axis labels\n')

    R_out.write(f'     labs(y= "-log10({pval_type})", x = "chromosome position")+\n')

    #only add in thresholds if the pval type isnt abs theta

    R_out.write(f'     # threshold line (95th)\n')
    R_out.write(f'     geom_hline(yintercept=NF_THRESHOLD, linetype="dashed", color = "red")+\n')
    R_out.write(f'     #threshold line (99th)\n')
    R_out.write(f'     geom_hline(yintercept=NN_THRESHOLD, linetype="dashed", color = "blue")+\n')
        # find way to automate the x position
    R_out.write(f'     annotate("text",x=20000000,y=NF_THRESHOLD-1,size=2.5,label=paste0("% above 95th%: ",percent_sig))+\n')


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
    R_batch.write(f'#SBATCH --job-name=R_subrun_IDEA3\n')
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

# IDEA 1.4.1 (MODIFIED) (thresholds - testing)
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
    CURRENT_SNP_R_SCRIPT.write(f'library("ggplot2")\n')
    CURRENT_SNP_R_SCRIPT.write(f'library("dplyr")\n')
    CURRENT_SNP_R_SCRIPT.write(f'library("ggstatsplot")\n')
    #CURRENT_SNP_R_SCRIPT.write(f'library("tidyverse")\n')
    CURRENT_SNP_R_SCRIPT.write(f'print("start of IDEA1 R script")\n')
    CURRENT_SNP_R_SCRIPT.write(f'T20_TRACKED_DATA<- read.csv("{PATH_TO_MAIN}output_files/R_DATA/{cumulative_t20_dataframe_name}.csv", header= TRUE, sep=",")\n')
    CURRENT_SNP_R_SCRIPT.write(f'\n')
    CURRENT_SNP_R_SCRIPT.write(f'PVAL_TYPES_LIST = c("GWAS_P","PSNP8")\n')
    CURRENT_SNP_R_SCRIPT.write(f'SUBSAMPLE_NUM_LIST = c({my_subsample_number_list})\n')
    CURRENT_SNP_R_SCRIPT.write(f'\n')
    CURRENT_SNP_R_SCRIPT.write(f'my_cols=c("CHR","POS","PVAL_TYPE","SUBSAMPLE_NUM","VALUE")\n')
    CURRENT_SNP_R_SCRIPT.write(f'# Make empty dataframe to append to later\n')
    CURRENT_SNP_R_SCRIPT.write(f'final_df = data.frame(matrix(nrow=0,ncol=length(my_cols)))\n')
    CURRENT_SNP_R_SCRIPT.write(f'colnames(final_df)=my_cols\n')
    CURRENT_SNP_R_SCRIPT.write(f'\n')
    CURRENT_SNP_R_SCRIPT.write('for (pval_type in PVAL_TYPES_LIST){\n')
    CURRENT_SNP_R_SCRIPT.write(f'   print(pval_type)\n')
    CURRENT_SNP_R_SCRIPT.write('    for (subsample_num in SUBSAMPLE_NUM_LIST){\n')
    CURRENT_SNP_R_SCRIPT.write(f'   new_df<-subset(T20_TRACKED_DATA,PVAL_TYPE==pval_type & SUBSAMPLE_NUM==subsample_num)\n')
    # only aggregate if theres something in the new df
    CURRENT_SNP_R_SCRIPT.write('    if(nrow(new_df)>0){\n')
    CURRENT_SNP_R_SCRIPT.write(f'   #average the values out (gets 20 averages)\n')
    CURRENT_SNP_R_SCRIPT.write(f'   new_df=aggregate(VALUE~CHR+POS,new_df,mean)\n')
    CURRENT_SNP_R_SCRIPT.write(f'   # add back the lost columns\n')
    CURRENT_SNP_R_SCRIPT.write(f'   new_df$PVAL_TYPE<-pval_type\n')
    CURRENT_SNP_R_SCRIPT.write(f'   new_df$SUBSAMPLE_NUM<-subsample_num\n')
   

    # CURRENT_SNP_R_SCRIPT.write('   if (pval_type=="ABS_THETA"){\n')
    # CURRENT_SNP_R_SCRIPT.write(f'       z_scores<-(new_df$VALUE - mean(new_df$VALUE)) / sd(new_df$VALUE)\n')
    # CURRENT_SNP_R_SCRIPT.write(f'       outliers<-abs(z_scores)>3\n')
    # CURRENT_SNP_R_SCRIPT.write(f'       filtered_new_df<-new_df[!outliers,]\n')
    # CURRENT_SNP_R_SCRIPT.write('   }else{\n')
    # CURRENT_SNP_R_SCRIPT.write(f'       z_scores<-(new_df$VALUE - mean(new_df$VALUE))/ sd(new_df$VALUE)\n')
    # CURRENT_SNP_R_SCRIPT.write(f'       outliers<-abs(z_scores)>3\n')
    # CURRENT_SNP_R_SCRIPT.write(f'       filtered_new_df<-new_df[!outliers,]\n')
        # new since avoiding outlier filtering

    CURRENT_SNP_R_SCRIPT.write(f'       filtered_new_df<-new_df\n')
    CURRENT_SNP_R_SCRIPT.write(f'       # apply -log10 function\n')
    CURRENT_SNP_R_SCRIPT.write(f'       filtered_new_df<-filtered_new_df%>%\n')
    CURRENT_SNP_R_SCRIPT.write(f'           mutate(VALUE=-log10(VALUE))\n')
    CURRENT_SNP_R_SCRIPT.write(f'       # change name to -log10(pvaltype)\n')
    CURRENT_SNP_R_SCRIPT.write(f'       filtered_new_df<-filtered_new_df%>%\n')
    CURRENT_SNP_R_SCRIPT.write(f'           mutate(PVAL_TYPE=paste("-log10(",PVAL_TYPE,")"))\n')
    CURRENT_SNP_R_SCRIPT.write(f'   final_df<-rbind(final_df,filtered_new_df)\n')
    CURRENT_SNP_R_SCRIPT.write('        }\n')
    CURRENT_SNP_R_SCRIPT.write('   }\n')
    CURRENT_SNP_R_SCRIPT.write('}\n')
    CURRENT_SNP_R_SCRIPT.write(f'\n')
    CURRENT_SNP_R_SCRIPT.write(f'# Convert the subsample number to a FACTOR variable\n')
    CURRENT_SNP_R_SCRIPT.write(f'final_df$SUBSAMPLE_NUM<-factor(final_df$SUBSAMPLE_NUM)\n')
    CURRENT_SNP_R_SCRIPT.write(f'\n')
    CURRENT_SNP_R_SCRIPT.write(f'threshold_data_csv<-read.csv("{PATH_TO_MAIN}output_files/R_DATA/THRESHOLDS.csv", header= TRUE, sep=",")\n')

    # new graph plot with NO facet mode and define axis
    # my_pval_list=["GWAS_P","PSNP4","PSNP5","ABS_THETA"]
    # av_list=["AVERAGE_P","AVERAGE_PSNP4","AVERAGE_PSNP5"] # "AVERAGE_ABS_THETA" Not needed here
        # new
    my_pval_list=["GWAS_P","PSNP8",]
    av_list=["AVERAGE_P","AVERAGE_PSNP8"] 

    pval_index = 0
    # specific y lims for pvals
    if phenotype=="Mo98":
        # have this change depending on a multiple of the highest pvalue
        ylim=80

    elif phenotype=="Na23":
        # have this change depending on a multiple of the highest pvalue
        ylim=75

    for item in my_pval_list:


        # import and calculate with threshold data
        CURRENT_SNP_R_SCRIPT.write(f'# subset for the current phenotype and pval type\n')
        CURRENT_SNP_R_SCRIPT.write(f'threshold_data<-subset(threshold_data_csv,PHENOTYPE=="{phenotype}" & PVAL_TYPE=="{av_list[pval_index]}")\n')
       
        CURRENT_SNP_R_SCRIPT.write(f'#subset for each type of threshold\n')
        # CURRENT_SNP_R_SCRIPT.write(f'#Bonferroni (BF)\n')
        # CURRENT_SNP_R_SCRIPT.write(f'BF_threshold_data<-subset(threshold_data,THRESHOLD_TYPE=="BF" )\n')
        # CURRENT_SNP_R_SCRIPT.write(f'BF_threshold_data$SUBSAMPLE_NUM<-factor(BF_threshold_data$SUBSAMPLE_NUM)\n')
        # CURRENT_SNP_R_SCRIPT.write(f'\n')
        # CURRENT_SNP_R_SCRIPT.write(f'#(BHY)\n')
        # CURRENT_SNP_R_SCRIPT.write(f'BHY_threshold_data<-subset(threshold_data,THRESHOLD_TYPE=="BHY" )\n')
        # CURRENT_SNP_R_SCRIPT.write(f'BHY_threshold_data$SUBSAMPLE_NUM<-factor(BHY_threshold_data$SUBSAMPLE_NUM)\n')
        # CURRENT_SNP_R_SCRIPT.write(f'\n')
            # im going to use...
            # GWAS
                # 95% AND 99% BF 
            # GIFT
                # 95% and 99% NH (null hypothesis)

# replace all instances of "BF_threshold_data" for (NF)95%_threshold_data (NF_threshold_data)
# replace all instances of "BHY_threshold_data" for (NN)95% _threshold_data (NN_threshold_data)
        if item=="GWAS_P":
            CURRENT_SNP_R_SCRIPT.write(f'#95th percentile data (bonferroni)\n')
            CURRENT_SNP_R_SCRIPT.write(f'NF_threshold_data<-subset(threshold_data,THRESHOLD_TYPE=="NFBF" )\n')
            CURRENT_SNP_R_SCRIPT.write(f'NF_threshold_data$SUBSAMPLE_NUM<-factor(NF_threshold_data$SUBSAMPLE_NUM)\n')
            CURRENT_SNP_R_SCRIPT.write(f'\n')
            CURRENT_SNP_R_SCRIPT.write(f'#99th percentile data (bonferroni)\n')
            CURRENT_SNP_R_SCRIPT.write(f'NN_threshold_data<-subset(threshold_data,THRESHOLD_TYPE=="NNBF" )\n')
            CURRENT_SNP_R_SCRIPT.write(f'NN_threshold_data$SUBSAMPLE_NUM<-factor(NN_threshold_data$SUBSAMPLE_NUM)\n')
            CURRENT_SNP_R_SCRIPT.write(f'\n')

        elif item =="PSNP8":
            CURRENT_SNP_R_SCRIPT.write(f'#95th percentile data (Null hypothesis data)\n')
            CURRENT_SNP_R_SCRIPT.write(f'NF_threshold_data<-subset(threshold_data,THRESHOLD_TYPE=="NFNH" )\n')
            CURRENT_SNP_R_SCRIPT.write(f'NF_threshold_data$SUBSAMPLE_NUM<-factor(NF_threshold_data$SUBSAMPLE_NUM)\n')
            CURRENT_SNP_R_SCRIPT.write(f'\n')
            CURRENT_SNP_R_SCRIPT.write(f'#99th percentile data (Null hypothesis data)\n')
            CURRENT_SNP_R_SCRIPT.write(f'NN_threshold_data<-subset(threshold_data,THRESHOLD_TYPE=="NNNH" )\n')
            CURRENT_SNP_R_SCRIPT.write(f'NN_threshold_data$SUBSAMPLE_NUM<-factor(NN_threshold_data$SUBSAMPLE_NUM)\n')
            CURRENT_SNP_R_SCRIPT.write(f'\n')



        CURRENT_SNP_R_SCRIPT.write(f'# dataframe for 95% threshold info\n')
        CURRENT_SNP_R_SCRIPT.write(f'threshold_lines1 <- data.frame(\n')
        CURRENT_SNP_R_SCRIPT.write(f'   SUBSAMPLE_NUM = NF_threshold_data$SUBSAMPLE_NUM,\n')
        CURRENT_SNP_R_SCRIPT.write(f'   x = as.numeric(NF_threshold_data$SUBSAMPLE_NUM) - 0.25,\n')
        CURRENT_SNP_R_SCRIPT.write(f'   xend = as.numeric(NF_threshold_data$SUBSAMPLE_NUM) + 0.25,\n')
        CURRENT_SNP_R_SCRIPT.write(f'   y = NF_threshold_data$THRESHOLD_VALUE,\n')
        CURRENT_SNP_R_SCRIPT.write(f'   yend = NF_threshold_data$THRESHOLD_VALUE,\n')
        CURRENT_SNP_R_SCRIPT.write(f'   threshold_type = "95% Threshold")\n')
        CURRENT_SNP_R_SCRIPT.write(f'\n')
        CURRENT_SNP_R_SCRIPT.write(f'# dataframe for the 99% threshold\n')
        CURRENT_SNP_R_SCRIPT.write(f'threshold_lines2 <- data.frame(\n')
        CURRENT_SNP_R_SCRIPT.write(f'   SUBSAMPLE_NUM = NN_threshold_data$SUBSAMPLE_NUM,\n')
        CURRENT_SNP_R_SCRIPT.write(f'   x = as.numeric(NN_threshold_data$SUBSAMPLE_NUM) - 0.25,\n')
        CURRENT_SNP_R_SCRIPT.write(f'   xend = as.numeric(NN_threshold_data$SUBSAMPLE_NUM) + 0.25,\n')
        CURRENT_SNP_R_SCRIPT.write(f'   y = NN_threshold_data$THRESHOLD_VALUE,\n')
        CURRENT_SNP_R_SCRIPT.write(f'   yend = NN_threshold_data$THRESHOLD_VALUE,\n')
        CURRENT_SNP_R_SCRIPT.write(f'   threshold_type = "99% Threshold")\n')
        CURRENT_SNP_R_SCRIPT.write(f'\n')




        CURRENT_SNP_R_SCRIPT.write(f'# combine the threshold dataframes\n')
        CURRENT_SNP_R_SCRIPT.write(f'threshold_lines <- rbind(threshold_lines1, threshold_lines2)\n')
        CURRENT_SNP_R_SCRIPT.write(f'\n')

        # main data graph
        CURRENT_SNP_R_SCRIPT.write(f'#open png\n')
        # CURRENT_SNP_R_SCRIPT.write(f'png("{PATH_TO_MAIN}output_files/summary_plots/IDEA1/{cumulative_t20_dataframe_name}_{item}.png", bg = "white", width = 5.75, height = 8.25, units = "in", res = 1200, pointsize = 4)\n')
            # new
        CURRENT_SNP_R_SCRIPT.write(f'png("{PATH_TO_MAIN}output_files/summary_plots/IDEA1/{cumulative_t20_dataframe_name}_{item}.png", bg = "white", width = 5, height = 7, units = "in", res = 1200, pointsize = 3)\n')
        CURRENT_SNP_R_SCRIPT.write(f'this_graph_data<-subset(final_df,PVAL_TYPE=="-log10( {item} )")\n')
        CURRENT_SNP_R_SCRIPT.write(f'myplot<-ggplot(this_graph_data, aes(x=SUBSAMPLE_NUM, y=VALUE, fill=SUBSAMPLE_NUM))+\n')
        CURRENT_SNP_R_SCRIPT.write(f'   guides(fill = guide_legend(override.aes = list(size = 7)))+\n')
        CURRENT_SNP_R_SCRIPT.write(f'   geom_boxplot(outlier.shape = NA) +\n')
        CURRENT_SNP_R_SCRIPT.write(f'   geom_point(size=0.5,color="black",alpha = 1,position = position_jitterdodge(jitter.width = 0.1),show.legend=FALSE)+\n')
        CURRENT_SNP_R_SCRIPT.write(f'   scale_y_continuous(limits=c(0,{ylim}))+\n')
        CURRENT_SNP_R_SCRIPT.write(f'   geom_segment(data = threshold_lines,\n')
        CURRENT_SNP_R_SCRIPT.write(f'       aes(x = x, xend = xend, y = y, yend = yend, color = threshold_type), \n')
        CURRENT_SNP_R_SCRIPT.write(f'       linetype = "dashed", show.legend = TRUE) +\n')
        CURRENT_SNP_R_SCRIPT.write(f'   scale_color_manual(name = "Thresholds",\n')
        CURRENT_SNP_R_SCRIPT.write(f'       values = c("95% Threshold" = "red", "99% Threshold" = "blue")) +\n')
        CURRENT_SNP_R_SCRIPT.write(f'   scale_fill_manual(name = "Subsample numbers",\n')
        CURRENT_SNP_R_SCRIPT.write(f'       values = c("200" = "#e04e43", "400" = "#e0892b", "600" = "#d1d145", "800" = "#5fb547", "999" = "#51dee0")) +\n')
        CURRENT_SNP_R_SCRIPT.write(f'   guides(fill = guide_legend(order = 1), color = guide_legend(order = 2)) +\n')
        CURRENT_SNP_R_SCRIPT.write(f'   theme(legend.position = "right",\n')
        CURRENT_SNP_R_SCRIPT.write(f'       legend.text = element_text(size = 13),\n')
        CURRENT_SNP_R_SCRIPT.write(f'       legend.title = element_text(size = 13),\n')
        CURRENT_SNP_R_SCRIPT.write(f'       legend.key.size = unit(0.8,"cm"))\n')
        CURRENT_SNP_R_SCRIPT.write(f'print(myplot)\n')
        CURRENT_SNP_R_SCRIPT.write(f'\n')
        CURRENT_SNP_R_SCRIPT.write(f'dev.off()\n')
        CURRENT_SNP_R_SCRIPT.write(f'\n')

        # stats test (KW)
            # find way to make text bigger for this one (so more readable)
        CURRENT_SNP_R_SCRIPT.write(f'png("{PATH_TO_MAIN}output_files/summary_plots/IDEA1/{cumulative_t20_dataframe_name}_{item}_KW_TEST.png", bg = "white", width = 6.25, height = 8.25, units = "in", res = 1200, pointsize = 4)\n')
        CURRENT_SNP_R_SCRIPT.write(f'ggbetweenstats(\n')
        CURRENT_SNP_R_SCRIPT.write(f'   data=this_graph_data,\n')
        CURRENT_SNP_R_SCRIPT.write(f'   x = SUBSAMPLE_NUM,\n')
        CURRENT_SNP_R_SCRIPT.write(f'   y=VALUE,\n')
        CURRENT_SNP_R_SCRIPT.write(f'   type = "nonparametric",\n')
        CURRENT_SNP_R_SCRIPT.write(f'   p.adjust.method = "bonferroni",\n')
        CURRENT_SNP_R_SCRIPT.write(f'   pairwise.display = "all")\n')
        CURRENT_SNP_R_SCRIPT.write(f'dev.off()\n')
        CURRENT_SNP_R_SCRIPT.write(f'\n')

        # increment the pval_index 
        pval_index+=1


    CURRENT_SNP_R_SCRIPT.write(f'print("End of IDEA1 R script")\n')
    CURRENT_SNP_R_SCRIPT.write(f'# END OF R SCRIPT')
    CURRENT_SNP_R_SCRIPT.close()

# IDEA 1.5.1 (MODIFIED)
def IDEA_1_MAKE_BASH_SCRIPT(cumulative_t20_dataframe_name):
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
    CURRENT_SNP_BATCH.write(f'#SBATCH --job-name=R_subrun_IDEA1\n')
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

# IDEA 2.3.1 (MODIFIED)
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
    # #install.packages("ggstatsplot") through conda first
    #CURRENT_SNP_R_SCRIPT.write(f'library("tidyverse")\n')
    CURRENT_SNP_R_SCRIPT.write(f'library("ggplot2")\n')
    CURRENT_SNP_R_SCRIPT.write(f'library("dplyr")\n')
    CURRENT_SNP_R_SCRIPT.write(f'library("ggstatsplot")\n')
    CURRENT_SNP_R_SCRIPT.write(f'print("start of IDEA 2 R script")\n')
    CURRENT_SNP_R_SCRIPT.write(f'csv_data<- read.csv("{PATH_TO_MAIN}output_files/R_DATA/{control_dataframe_name}.csv", header= TRUE, sep=",")\n')
    CURRENT_SNP_R_SCRIPT.write(f'\n')
    # CURRENT_SNP_R_SCRIPT.write(f'PVAL_TYPES_LIST = c("GWAS_P","PSNP4","PSNP5","ABS_THETA")\n')
    # CURRENT_SNP_R_SCRIPT.write(f'SUBSAMPLE_NUM_LIST = c(200,400,600,800,1000)\n')
        # new
    CURRENT_SNP_R_SCRIPT.write(f'PVAL_TYPES_LIST = c("GWAS_P","PSNP8")\n')
    CURRENT_SNP_R_SCRIPT.write(f'SUBSAMPLE_NUM_LIST = c(200,400,600,800,999)\n')
    CURRENT_SNP_R_SCRIPT.write(f'\n')
    CURRENT_SNP_R_SCRIPT.write(f'my_cols=c("CHR","POS","PVAL_TYPE","SUBSAMPLE_NUM","VALUE")\n')
    CURRENT_SNP_R_SCRIPT.write(f'final_df = data.frame(matrix(nrow=0,ncol=length(my_cols)))\n')
    CURRENT_SNP_R_SCRIPT.write(f'colnames(final_df)=my_cols\n')
    CURRENT_SNP_R_SCRIPT.write(f'\n')
    # CURRENT_SNP_R_SCRIPT.write('for (pval_type in PVAL_TYPES_LIST){\n')
    # CURRENT_SNP_R_SCRIPT.write('    for (subsample_num in SUBSAMPLE_NUM_LIST){\n')
    # CURRENT_SNP_R_SCRIPT.write(f'           new_df<-subset(csv_data,PVAL_TYPE==pval_type & SUBSAMPLE_NUM==subsample_num)\n')
    # CURRENT_SNP_R_SCRIPT.write('           if(nrow(new_df)>0){\n')
    # CURRENT_SNP_R_SCRIPT.write(f'           #average the values out if at the same loci\n')
    # CURRENT_SNP_R_SCRIPT.write(f'           new_df=aggregate(VALUE~CHR+POS,new_df,mean)\n')
    # CURRENT_SNP_R_SCRIPT.write(f'           # add back the lost columns\n')
    # CURRENT_SNP_R_SCRIPT.write(f'           new_df$PVAL_TYPE<-pval_type\n')
    # CURRENT_SNP_R_SCRIPT.write(f'           new_df$SUBSAMPLE_NUM<-subsample_num\n')
    # CURRENT_SNP_R_SCRIPT.write('           }\n')
    # CURRENT_SNP_R_SCRIPT.write(f'\n')
    # CURRENT_SNP_R_SCRIPT.write('            if (pval_type=="ABS_THETA"){\n')
    # CURRENT_SNP_R_SCRIPT.write(f'               z_scores<-(new_df$VALUE - mean(new_df$VALUE)) / sd(new_df$VALUE)\n')
    # CURRENT_SNP_R_SCRIPT.write(f'               outliers<-abs(z_scores)>3\n')
    # CURRENT_SNP_R_SCRIPT.write(f'               filtered_new_df<-new_df[!outliers,]\n')
    # CURRENT_SNP_R_SCRIPT.write('           }else{\n')
    # CURRENT_SNP_R_SCRIPT.write(f'               z_scores<-(new_df$VALUE - mean(new_df$VALUE))/ sd(new_df$VALUE)\n')
    #     # maybe change how outliers work -> might not remove them anymore in future? or at least make supplementary file with them in 
    # CURRENT_SNP_R_SCRIPT.write(f'               outliers<-abs(z_scores)>3\n')
    # CURRENT_SNP_R_SCRIPT.write(f'               filtered_new_df<-new_df[!outliers,]\n')
    # CURRENT_SNP_R_SCRIPT.write(f'               # apply -log10 function\n')
    # CURRENT_SNP_R_SCRIPT.write(f'               filtered_new_df<-filtered_new_df%>%\n')
    # CURRENT_SNP_R_SCRIPT.write(f'                   mutate(VALUE=-log10(VALUE))\n')
    # CURRENT_SNP_R_SCRIPT.write(f'                # change name to -log10(pvaltype)\n')
    # CURRENT_SNP_R_SCRIPT.write(f'               filtered_new_df<-filtered_new_df%>%\n')
    # CURRENT_SNP_R_SCRIPT.write(f'                   mutate(PVAL_TYPE=paste("-log10(",PVAL_TYPE,")"))\n')
    # CURRENT_SNP_R_SCRIPT.write('               }\n')
    # CURRENT_SNP_R_SCRIPT.write(f'               final_df<-rbind(final_df,filtered_new_df)\n')
    # CURRENT_SNP_R_SCRIPT.write('       }\n')
    # CURRENT_SNP_R_SCRIPT.write('   }\n')
    # CURRENT_SNP_R_SCRIPT.write(f'\n')

        # new (removed absolute theta check)
    CURRENT_SNP_R_SCRIPT.write('for (pval_type in PVAL_TYPES_LIST){\n')
    CURRENT_SNP_R_SCRIPT.write('    for (subsample_num in SUBSAMPLE_NUM_LIST){\n')
    CURRENT_SNP_R_SCRIPT.write(f'           new_df<-subset(csv_data,PVAL_TYPE==pval_type & SUBSAMPLE_NUM==subsample_num)\n')
    CURRENT_SNP_R_SCRIPT.write('           if(nrow(new_df)>0){\n')
    CURRENT_SNP_R_SCRIPT.write(f'           #average the values out if at the same loci\n')
    CURRENT_SNP_R_SCRIPT.write(f'           new_df=aggregate(VALUE~CHR+POS,new_df,mean)\n')
    CURRENT_SNP_R_SCRIPT.write(f'           # add back the lost columns\n')
    CURRENT_SNP_R_SCRIPT.write(f'           new_df$PVAL_TYPE<-pval_type\n')
    CURRENT_SNP_R_SCRIPT.write(f'           new_df$SUBSAMPLE_NUM<-subsample_num\n')
    # CURRENT_SNP_R_SCRIPT.write(f'           z_scores<-(new_df$VALUE - mean(new_df$VALUE))/ sd(new_df$VALUE)\n')
    #     # maybe change how outliers work -> might not remove them anymore in future? or at least make supplementary file with them in 
    # CURRENT_SNP_R_SCRIPT.write(f'           outliers<-abs(z_scores)>3\n')
    # CURRENT_SNP_R_SCRIPT.write(f'           filtered_new_df<-new_df[!outliers,]\n')
    CURRENT_SNP_R_SCRIPT.write(f'           filtered_new_df<-new_df\n')
    CURRENT_SNP_R_SCRIPT.write(f'           # apply -log10 function\n')
    CURRENT_SNP_R_SCRIPT.write(f'           filtered_new_df<-filtered_new_df%>%\n')
    CURRENT_SNP_R_SCRIPT.write(f'               mutate(VALUE=-log10(VALUE))\n')
    CURRENT_SNP_R_SCRIPT.write(f'            # change name to -log10(pvaltype)\n')
    CURRENT_SNP_R_SCRIPT.write(f'           filtered_new_df<-filtered_new_df%>%\n')
    CURRENT_SNP_R_SCRIPT.write(f'               mutate(PVAL_TYPE=paste("-log10(",PVAL_TYPE,")"))\n')
    CURRENT_SNP_R_SCRIPT.write(f'           final_df<-rbind(final_df,filtered_new_df)\n')
    CURRENT_SNP_R_SCRIPT.write('            }\n')
    CURRENT_SNP_R_SCRIPT.write('       }\n')
    CURRENT_SNP_R_SCRIPT.write('   }\n')
    CURRENT_SNP_R_SCRIPT.write(f'\n')
    CURRENT_SNP_R_SCRIPT.write(f'# Convert the subsample number to a FACTOR variable\n')
    CURRENT_SNP_R_SCRIPT.write(f'final_df$SUBSAMPLE_NUM<-factor(final_df$SUBSAMPLE_NUM)\n')
    CURRENT_SNP_R_SCRIPT.write(f'\n')

        # # new graph plot with NO facet mode and define axis
        # my_pval_list=["GWAS_P","PSNP4","PSNP5","ABS_THETA"]

        # # list of titles of the averages (needed for analysis of threshold)
        # av_list=["AVERAGE_P","AVERAGE_PSNP4","AVERAGE_PSNP5"] # "AVERAGE_ABS_THETA" Not needed here

    # new
        # new graph plot with NO facet mode and define axis
    my_pval_list=["GWAS_P","PSNP8"]

    # list of titles of the averages (needed for analysis of threshold)
    av_list=["AVERAGE_P","AVERAGE_PSNP8"]

    # writes in the threshold data once to save space
    CURRENT_SNP_R_SCRIPT.write(f'threshold_data_csv<-read.csv("{PATH_TO_MAIN}output_files/R_DATA/THRESHOLDS.csv", header= TRUE, sep=",")\n')

    #specific y limit to make the data spread out as much as possible
        # change these to be automated when possible
    if phenotype=="Mo98":

        if positive_or_negative=="positive":
            ylim=75

        elif positive_or_negative=="negative":
            ylim=25
            
    elif phenotype=="Na23":

        if positive_or_negative=="positive":
            ylim=50

        elif positive_or_negative=="negative":
            ylim=20
    pval_index = 0

    for pval_type in my_pval_list:

        # no thresholds implemented for ABS_THETA, so just draws in the boxplot as needed
        # if pval_type == "ABS_THETA":
        #     # main data graph
        #     CURRENT_SNP_R_SCRIPT.write(f'#open png\n')
        #     CURRENT_SNP_R_SCRIPT.write(f'png("{PATH_TO_MAIN}output_files/summary_plots/IDEA2/{control_dataframe_name}_{pval_type}.png", bg = "white", width = 5.75, height = 8.25, units = "in", res = 1200, pointsize = 4)\n')
        #     CURRENT_SNP_R_SCRIPT.write(f'this_graph_data<-subset(final_df,PVAL_TYPE=="{pval_type}")\n')
        #     CURRENT_SNP_R_SCRIPT.write(f'myplot<-ggplot(this_graph_data, aes(x=SUBSAMPLE_NUM, y=VALUE, fill=SUBSAMPLE_NUM))+ \n')
        #     CURRENT_SNP_R_SCRIPT.write(f'   guides(fill = guide_legend(override.aes = list(size = 7)))+\n')
        #     CURRENT_SNP_R_SCRIPT.write(f'   geom_boxplot(outlier.shape = NA) +\n')
        #     CURRENT_SNP_R_SCRIPT.write(f'   geom_point(size=0.5,color="black",alpha = 1,position = position_jitterdodge(jitter.width = 0.1),show.legend=FALSE)+\n')
        #     CURRENT_SNP_R_SCRIPT.write(f'   scale_y_continuous(limits=c(0,125))+\n')
        #     CURRENT_SNP_R_SCRIPT.write(f'   scale_fill_manual(name = "Subsample numbers",\n')
        #     CURRENT_SNP_R_SCRIPT.write(f'       values = c("200" = "#e04e43", "400" = "#e0892b", "600" = "#d1d145", "800" = "#5fb547", "1000" = "#51dee0"))+\n')
        #     CURRENT_SNP_R_SCRIPT.write(f'    guides(fill = guide_legend(order = 1))+\n')
        #     CURRENT_SNP_R_SCRIPT.write(f'   theme(legend.position = "right",\n')
        #     CURRENT_SNP_R_SCRIPT.write(f'       legend.text = element_text(size = 13),\n')
        #     CURRENT_SNP_R_SCRIPT.write(f'       legend.title = element_text(size = 13),\n')
        #     CURRENT_SNP_R_SCRIPT.write(f'       legend.key.size = unit(0.8,"cm"))\n')
        #     CURRENT_SNP_R_SCRIPT.write(f'print(myplot)\n')
        #     CURRENT_SNP_R_SCRIPT.write(f'dev.off()\n')
        #     CURRENT_SNP_R_SCRIPT.write(f'\n')
        #     # KW stats test
        #     CURRENT_SNP_R_SCRIPT.write(f'png("{PATH_TO_MAIN}output_files/summary_plots/IDEA2/{control_dataframe_name}_{pval_type}_KW_TEST.png", bg = "white", width = 6.25, height = 8.25, units = "in", res = 1200, pointsize = 4)\n')
        #     CURRENT_SNP_R_SCRIPT.write(f'ggbetweenstats(\n')
        #     CURRENT_SNP_R_SCRIPT.write(f'   data=this_graph_data,\n')
        #     CURRENT_SNP_R_SCRIPT.write(f'   x = SUBSAMPLE_NUM,\n')
        #     CURRENT_SNP_R_SCRIPT.write(f'   y=VALUE,\n')
        #     CURRENT_SNP_R_SCRIPT.write(f'   type = "nonparametric",\n')
        #     CURRENT_SNP_R_SCRIPT.write(f'   p.adjust.method = "bonferroni",\n')
        #     CURRENT_SNP_R_SCRIPT.write(f'   pairwise.display = "all")\n')
        #     CURRENT_SNP_R_SCRIPT.write(f'dev.off()\n')

        # else:
            # end of if-else (DEPRECATED)

        # import threshold data
        CURRENT_SNP_R_SCRIPT.write(f'# subset for the current phenotype and pval type\n')
        CURRENT_SNP_R_SCRIPT.write(f'threshold_data<-subset(threshold_data_csv,PHENOTYPE=="{phenotype}" & PVAL_TYPE=="{av_list[pval_index]}")\n')

        # CURRENT_SNP_R_SCRIPT.write(f'#subset for each type of threshold\n')
        # CURRENT_SNP_R_SCRIPT.write(f'#Bonferroni (BF)\n')
        # CURRENT_SNP_R_SCRIPT.write(f'BF_threshold_data<-subset(threshold_data,THRESHOLD_TYPE=="BF" )\n')
        # CURRENT_SNP_R_SCRIPT.write(f'BF_threshold_data$SUBSAMPLE_NUM<-factor(BF_threshold_data$SUBSAMPLE_NUM)\n')
        # CURRENT_SNP_R_SCRIPT.write(f'\n')
        # CURRENT_SNP_R_SCRIPT.write(f'#(BHY)\n')
        # CURRENT_SNP_R_SCRIPT.write(f'BHY_threshold_data<-subset(threshold_data,THRESHOLD_TYPE=="BHY" )\n')
        # CURRENT_SNP_R_SCRIPT.write(f'BHY_threshold_data$SUBSAMPLE_NUM<-factor(BHY_threshold_data$SUBSAMPLE_NUM)\n')
        # CURRENT_SNP_R_SCRIPT.write(f'\n')
        # CURRENT_SNP_R_SCRIPT.write(f'# dataframe for BF threshold info\n')
        # CURRENT_SNP_R_SCRIPT.write(f'threshold_lines1 <- data.frame(\n')
        # CURRENT_SNP_R_SCRIPT.write(f'   SUBSAMPLE_NUM = BF_threshold_data$SUBSAMPLE_NUM,\n')
        # CURRENT_SNP_R_SCRIPT.write(f'   x = as.numeric(BF_threshold_data$SUBSAMPLE_NUM) - 0.25,\n')
        # CURRENT_SNP_R_SCRIPT.write(f'   xend = as.numeric(BF_threshold_data$SUBSAMPLE_NUM) + 0.25,\n')
        # CURRENT_SNP_R_SCRIPT.write(f'   y = BF_threshold_data$THRESHOLD_VALUE,\n')
        # CURRENT_SNP_R_SCRIPT.write(f'   yend = BF_threshold_data$THRESHOLD_VALUE,\n')
        # CURRENT_SNP_R_SCRIPT.write(f'   threshold_type = "BF Threshold")\n')
        # CURRENT_SNP_R_SCRIPT.write(f'\n')
        # CURRENT_SNP_R_SCRIPT.write(f'# dataframe for the BHY threshold\n')
        # CURRENT_SNP_R_SCRIPT.write(f'threshold_lines2 <- data.frame(\n')
        # CURRENT_SNP_R_SCRIPT.write(f'   SUBSAMPLE_NUM = BHY_threshold_data$SUBSAMPLE_NUM,\n')
        # CURRENT_SNP_R_SCRIPT.write(f'   x = as.numeric(BHY_threshold_data$SUBSAMPLE_NUM) - 0.25,\n')
        # CURRENT_SNP_R_SCRIPT.write(f'   xend = as.numeric(BHY_threshold_data$SUBSAMPLE_NUM) + 0.25,\n')
        # CURRENT_SNP_R_SCRIPT.write(f'   y = BHY_threshold_data$THRESHOLD_VALUE,\n')
        # CURRENT_SNP_R_SCRIPT.write(f'   yend = BHY_threshold_data$THRESHOLD_VALUE,\n')
        # CURRENT_SNP_R_SCRIPT.write(f'   threshold_type = "BHY Threshold")\n')
        # CURRENT_SNP_R_SCRIPT.write(f'\n')

        #
        # replace all instances of "BF_threshold_data" for (NF)95%_threshold_data (NF_threshold_data)
        # replace all instances of "BHY_threshold_data" for (NN)95% _threshold_data (NN_threshold_data)
        if pval_type == "GWAS_P":
            CURRENT_SNP_R_SCRIPT.write(f'#95% Bonferroni \n')
            CURRENT_SNP_R_SCRIPT.write(f'NF_threshold_data<-subset(threshold_data,THRESHOLD_TYPE=="NFBF" )\n')
            CURRENT_SNP_R_SCRIPT.write(f'NF_threshold_data$SUBSAMPLE_NUM<-factor(NF_threshold_data$SUBSAMPLE_NUM)\n')
            CURRENT_SNP_R_SCRIPT.write(f'\n')
            CURRENT_SNP_R_SCRIPT.write(f'# 99% bonferroni \n')
            CURRENT_SNP_R_SCRIPT.write(f'NN_threshold_data<-subset(threshold_data,THRESHOLD_TYPE=="NNBF" )\n')
            CURRENT_SNP_R_SCRIPT.write(f'NN_threshold_data$SUBSAMPLE_NUM<-factor(NN_threshold_data$SUBSAMPLE_NUM)\n')
            CURRENT_SNP_R_SCRIPT.write(f'\n')

        elif pval_type == "PSNP8":
            CURRENT_SNP_R_SCRIPT.write(f'#95% null hypothesis\n')
            CURRENT_SNP_R_SCRIPT.write(f'NF_threshold_data<-subset(threshold_data,THRESHOLD_TYPE=="NFNH" )\n')
            CURRENT_SNP_R_SCRIPT.write(f'NF_threshold_data$SUBSAMPLE_NUM<-factor(NF_threshold_data$SUBSAMPLE_NUM)\n')
            CURRENT_SNP_R_SCRIPT.write(f'\n')
            CURRENT_SNP_R_SCRIPT.write(f'#99% null hypothesis\n')
            CURRENT_SNP_R_SCRIPT.write(f'NN_threshold_data<-subset(threshold_data,THRESHOLD_TYPE=="NNNH" )\n')
            CURRENT_SNP_R_SCRIPT.write(f'NN_threshold_data$SUBSAMPLE_NUM<-factor(NN_threshold_data$SUBSAMPLE_NUM)\n')
            CURRENT_SNP_R_SCRIPT.write(f'\n')


        CURRENT_SNP_R_SCRIPT.write(f'# dataframe for NF threshold info\n')
        CURRENT_SNP_R_SCRIPT.write(f'threshold_lines1 <- data.frame(\n')
        CURRENT_SNP_R_SCRIPT.write(f'   SUBSAMPLE_NUM = NF_threshold_data$SUBSAMPLE_NUM,\n')
        CURRENT_SNP_R_SCRIPT.write(f'   x = as.numeric(NF_threshold_data$SUBSAMPLE_NUM) - 0.25,\n')
        CURRENT_SNP_R_SCRIPT.write(f'   xend = as.numeric(NF_threshold_data$SUBSAMPLE_NUM) + 0.25,\n')
        CURRENT_SNP_R_SCRIPT.write(f'   y = NF_threshold_data$THRESHOLD_VALUE,\n')
        CURRENT_SNP_R_SCRIPT.write(f'   yend = NF_threshold_data$THRESHOLD_VALUE,\n')
        CURRENT_SNP_R_SCRIPT.write(f'   threshold_type = "95% Threshold")\n')
        CURRENT_SNP_R_SCRIPT.write(f'\n')
        CURRENT_SNP_R_SCRIPT.write(f'# dataframe for the BHY threshold\n')
        CURRENT_SNP_R_SCRIPT.write(f'threshold_lines2 <- data.frame(\n')
        CURRENT_SNP_R_SCRIPT.write(f'   SUBSAMPLE_NUM = NN_threshold_data$SUBSAMPLE_NUM,\n')
        CURRENT_SNP_R_SCRIPT.write(f'   x = as.numeric(NN_threshold_data$SUBSAMPLE_NUM) - 0.25,\n')
        CURRENT_SNP_R_SCRIPT.write(f'   xend = as.numeric(NN_threshold_data$SUBSAMPLE_NUM) + 0.25,\n')
        CURRENT_SNP_R_SCRIPT.write(f'   y = NN_threshold_data$THRESHOLD_VALUE,\n')
        CURRENT_SNP_R_SCRIPT.write(f'   yend = NN_threshold_data$THRESHOLD_VALUE,\n')
        CURRENT_SNP_R_SCRIPT.write(f'   threshold_type = "99% Threshold")\n')
        CURRENT_SNP_R_SCRIPT.write(f'\n')





        CURRENT_SNP_R_SCRIPT.write(f'# combine the threshold dataframes\n')
        CURRENT_SNP_R_SCRIPT.write(f'threshold_lines <- rbind(threshold_lines1, threshold_lines2)\n')
        CURRENT_SNP_R_SCRIPT.write(f'\n')
        # Main data plot
        CURRENT_SNP_R_SCRIPT.write(f'#open png\n')
        #CURRENT_SNP_R_SCRIPT.write(f'png("{PATH_TO_MAIN}output_files/summary_plots/IDEA2/{control_dataframe_name}_{pval_type}.png", bg = "white", width = 5.75, height = 8.25, units = "in", res = 1200, pointsize = 4)\n')
            # new
        CURRENT_SNP_R_SCRIPT.write(f'png("{PATH_TO_MAIN}output_files/summary_plots/IDEA2/{control_dataframe_name}_{pval_type}.png", bg = "white", width = 5, height = 7, units = "in", res = 1200, pointsize = 3)\n')
        CURRENT_SNP_R_SCRIPT.write(f'this_graph_data<-subset(final_df,PVAL_TYPE=="-log10( {pval_type} )")\n')
        CURRENT_SNP_R_SCRIPT.write(f'\n')
        CURRENT_SNP_R_SCRIPT.write(f'myplot<-ggplot(this_graph_data, aes(x=SUBSAMPLE_NUM, y=VALUE, fill=SUBSAMPLE_NUM))+ \n')
        CURRENT_SNP_R_SCRIPT.write(f'   guides(fill = guide_legend(override.aes = list(size = 7)))+\n')
        CURRENT_SNP_R_SCRIPT.write(f'   geom_boxplot(outlier.shape = NA) +\n')
        CURRENT_SNP_R_SCRIPT.write(f'   geom_point(size=0.5,color="black",alpha = 1,position = position_jitterdodge(jitter.width = 0.1),show.legend=FALSE)+\n')
        CURRENT_SNP_R_SCRIPT.write(f'   scale_y_continuous(limits=c(0,{ylim}))+\n')
        CURRENT_SNP_R_SCRIPT.write(f'   geom_segment(data = threshold_lines,\n')
        CURRENT_SNP_R_SCRIPT.write(f'       aes(x = x, xend = xend, y = y, yend = yend, color = threshold_type), \n')
        CURRENT_SNP_R_SCRIPT.write(f'       linetype = "dashed", show.legend = TRUE) +\n')
        CURRENT_SNP_R_SCRIPT.write(f'   scale_color_manual(name = "Thresholds",\n')
        CURRENT_SNP_R_SCRIPT.write(f'       values = c("95% Threshold" = "red", "99% Threshold" = "blue")) +\n')
        CURRENT_SNP_R_SCRIPT.write(f'   scale_fill_manual(name = "Subsample numbers",\n')
        CURRENT_SNP_R_SCRIPT.write(f'       values = c("200" = "#e04e43", "400" = "#e0892b", "600" = "#d1d145", "800" = "#5fb547", "999" = "#51dee0")) +\n')
        CURRENT_SNP_R_SCRIPT.write(f'   guides(fill = guide_legend(order = 1), color = guide_legend(order = 2)) +\n')
        CURRENT_SNP_R_SCRIPT.write(f'   theme(legend.position = "right",\n')
        CURRENT_SNP_R_SCRIPT.write(f'       legend.text = element_text(size = 13),\n')
        CURRENT_SNP_R_SCRIPT.write(f'       legend.title = element_text(size = 13),\n')
        CURRENT_SNP_R_SCRIPT.write(f'       legend.key.size = unit(0.8,"cm"))\n')
        CURRENT_SNP_R_SCRIPT.write(f'print(myplot)\n')
        CURRENT_SNP_R_SCRIPT.write(f'dev.off()\n')
        # Stats test plot (Kruskal wallace)
        CURRENT_SNP_R_SCRIPT.write(f'\n')
        CURRENT_SNP_R_SCRIPT.write(f'png("{PATH_TO_MAIN}output_files/summary_plots/IDEA2/{control_dataframe_name}_{pval_type}_KW_TEST.png", bg = "white", width = 6.25, height = 8.25, units = "in", res = 1200, pointsize = 4)\n')
        CURRENT_SNP_R_SCRIPT.write(f'ggbetweenstats(\n')
        CURRENT_SNP_R_SCRIPT.write(f'   data=this_graph_data,\n')
        CURRENT_SNP_R_SCRIPT.write(f'   x = SUBSAMPLE_NUM,\n')
        CURRENT_SNP_R_SCRIPT.write(f'   y=VALUE,\n')
        CURRENT_SNP_R_SCRIPT.write(f'   type = "nonparametric",\n')
        CURRENT_SNP_R_SCRIPT.write(f'   p.adjust.method = "bonferroni",\n')
        CURRENT_SNP_R_SCRIPT.write(f'   pairwise.display = "all")\n')
        CURRENT_SNP_R_SCRIPT.write(f'dev.off()\n')

        # increment the pval_index
        pval_index+=1


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
    CURRENT_SNP_BATCH.write(f'#SBATCH --job-name=R_subrun_IDEA2\n')
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


parser=argparse.ArgumentParser(description="W.I.P")

parser.add_argument('-subsampleFile', 
                    type=str, 
                    metavar='location of subsample file', 
                    required=True, 
                    help='location of the file containing list of all sample numbers'
                    )

# stores input data and parses them
args= parser.parse_args() 

subsample_num_list=fetch_subsample_numbers_list(args.subsampleFile)
my_subsample_number_list=""
for item in subsample_num_list:
    my_subsample_number_list += str(item)
    my_subsample_number_list += ","

my_subsample_number_list = my_subsample_number_list[:-1]

############################################
IDEA_2_MAKE_R_AND_BASH_SCRIPT("Mo98", "Mo98_positive_control.csv","positive")
IDEA_2_MAKE_R_AND_BASH_SCRIPT("Mo98", "Mo98_negative_control.csv","negative")
IDEA_2_MAKE_R_AND_BASH_SCRIPT("Na23", "Na23_positive_control.csv","positive")
IDEA_2_MAKE_R_AND_BASH_SCRIPT("Na23", "Na23_negative_control.csv","negative")

print("IDEA 2 FINISHED",flush=True)

#subsample_num_list=[200,400,600,800,1000] # can later update this to read from earlier scripts or something
    # new
# subsample_num_list=[200,400,600,800,999] # can later update this to read from earlier scripts or something

# subsample_num_list=[50,100,400,600,999] # can later update this to read from earlier scripts or something

# phenotype_list=["Mo98","Na23"] # can later update this to read from the phenotype text file

phenotype_list=fetch_phenotype_list(PATH_TO_MAIN+"core_files/phenotypes_list.txt")

# pvals=["AVERAGE_P","AVERAGE_PSNP4","AVERAGE_PSNP5","AVERAGE_ABS_THETA"] # these should always be the same 4 types
    # new
pvals=["AVERAGE_P","AVERAGE_PSNP8"] # these should always be the same 4 types

# R script creation and running
for phenotype in phenotype_list:

    for subsample_number in subsample_num_list:

        for pval_type in pvals:

            # Example variable inputs are as follows: 
            # example 1 Mo98,200,AVERAGE_P
            # example 2 Mo98,200,AVERAGE_PSNP8
            # ...
            # example x Mo98,400,AVERAGE_P
            # ....
            IDEA_3_R_AND_BATCH(phenotype,subsample_number,pval_type)

print("IDEA 3 FINISHED",flush=True)

# set up cumulative T20 dataframe paths for each phenotype
Mo98_cumulative_t20_dataframe_path = PATH_TO_MAIN+"output_files/R_DATA/Mo98_cumulative_t20_dataframe.csv"
Na23_cumulative_t20_dataframe_path = PATH_TO_MAIN+"output_files/R_DATA/Na23_cumulative_t20_dataframe.csv"


# make the R scripts for IDEA 1
IDEA_1_MAKE_R_SCRIPT("Mo98","Mo98_cumulative_t20_dataframe.csv")
IDEA_1_MAKE_R_SCRIPT("Na23","Na23_cumulative_t20_dataframe.csv")

# make the bash scripts for IDEA 1
# IDEA_1_MAKE_BASH_SCRIPT("Mo98","Mo98_cumulative_t20_dataframe.csv")
# IDEA_1_MAKE_BASH_SCRIPT("Na23","Na23_cumulative_t20_dataframe.csv")
    # new without phenotype
IDEA_1_MAKE_BASH_SCRIPT("Mo98_cumulative_t20_dataframe.csv")
IDEA_1_MAKE_BASH_SCRIPT("Na23_cumulative_t20_dataframe.csv")

print("END OF IDEA 1",flush=True)

print("END OF SNP_TRACKER_R_AND_BASH_MAKER",flush=True)
# end of file