### SCRIPT INFO
# Makes the R and bash scripts to get
    # 1) boxplot of standard deviations across all subsample levels
    # 2) individual manhattan plots of std dev results for each sub level

## packages
import pandas

### Functions

def fetch_phenotype_list(phenotype_list_file):
    phenotypes_list=[]
    temp_file=open(f"{phenotype_list_file}","r")

    for line in temp_file:

        subsample_number=line.replace('\n','')

        phenotypes_list.append(subsample_number)
    
    print("Phenotype list fetched: ",flush=True)
    print(phenotypes_list,flush=True)

    return phenotypes_list


def fetch_subsample_numbers_list(subsample_numbers_list_file):

    subsample_num_list=[]

    temp_file=open(f"{subsample_numbers_list_file}","r")

    for line in temp_file:

        subsample_number=line.replace('\n','')

        subsample_num_list.append(subsample_number)
    
    print("Subsample number list fetched: ",flush=True)
    print(subsample_num_list,flush=True)

    return subsample_num_list

# constants
# set constant for file path
PATH_TO_MAIN = "/gpfs01/home/mbysh17/"

# fetch phenotype and subsample number lists
phenotype_list=fetch_phenotype_list(PATH_TO_MAIN+"core_files/phenotypes_list.txt")
subsample_numbers_list=fetch_subsample_numbers_list(PATH_TO_MAIN+"core_files/subsample_numbers_list.txt")
method_list=["GIFT","GWAS"]

# file locations
R_script_out_directory=(PATH_TO_MAIN+"batch_files/stage_4_R_scripts/")
Shell_out_directory=(PATH_TO_MAIN+"batch_files/stage_4_prerun/")
Stdv_results_directory=(PATH_TO_MAIN+"output_files/R_DATA/std_dev_results.csv")
Plot_output_directory=(PATH_TO_MAIN+"output_files/STAGE_4_RESULTS/")

# plot the standard deviations summary as a box plot for all subsamples
R_out=open(str(R_script_out_directory)+"boxplot_stdv_plots.R","w")
R_out.write(f' \n')
R_out.write(f'library(tidyverse) \n')
R_out.write(f'std_dev_results<-read.csv("{Stdv_results_directory}", header= TRUE, sep=",") \n')
R_out.write(f'std_dev_results$PVAL_TYPE <- as.factor(std_dev_results$PVAL_TYPE) \n')
R_out.write(f'std_dev_results$SUBSAMPLE_NUM <- as.factor(std_dev_results$SUBSAMPLE_NUM) \n')
R_out.write(f' \n')
R_out.write(f'png("{Plot_output_directory}boxplot.png", unit="px", width=1080, height=720,res=150) \n')
#R_out.write(f'png("{Plot_output_directory}{method_type}_{subsample_number}_MANHATTAN.png", bg = "white", width = 7.5, height = 2.5, units = "in", res = 1200, pointsize = 3) \n')

R_out.write(f'ggplot(std_dev_results, aes(x=SUBSAMPLE_NUM, y=pvalue_stdv, fill=PVAL_TYPE)) + \n')
R_out.write(f' geom_boxplot()\n')
R_out.write(f'dev.off() \n')
R_out.write(f' \n')
R_out.close()

# plot the standard deviations as a manhattan plot for each subsample
    # update later for multiple phenotypes
phenotype="Mo98"

for subsample_number in subsample_numbers_list:
    R_out=open(str(R_script_out_directory)+str(subsample_number)+"_stdv_plots.R","w")
    R_out.write(f'# STAGE 4 R SCRIPT \n')
    R_out.write(f'library(tidyverse) \n')
    R_out.write(f' \n')

        # may not need the quotes here: keep an eye out for this
    R_out.write(f'std_dev_results<-read.csv("{Stdv_results_directory}", header= TRUE, sep=",") \n')
    R_out.write(f' \n')
    R_out.write(f'std_dev_results$PVAL_TYPE <- as.factor(std_dev_results$PVAL_TYPE) \n')
    R_out.write(f'std_dev_results$SUBSAMPLE_NUM <- as.factor(std_dev_results$SUBSAMPLE_NUM) \n')
    R_out.write(f' \n')
    R_out.write(f'subsample_dataset_GIFT<-subset(std_dev_results,SUBSAMPLE_NUM=={subsample_number} & PVAL_TYPE=="PSNP8")\n')
    R_out.write(f'subsample_dataset_GWAS<-subset(std_dev_results,SUBSAMPLE_NUM=={subsample_number} & PVAL_TYPE=="GWAS_P")\n')
    R_out.write(f' \n')



    for method_type in method_list:
            
        # obtain threshold data
        threshold_data =pandas.read_csv(f'{PATH_TO_MAIN}output_files/R_DATA/THRESHOLDS.csv')

        if method_type=="GWAS":      
            average_column_name="AVERAGE_P"

            # subset the threshold data to get the threshold at current subsamples:
            subsetted_threshold_data=threshold_data[(threshold_data['PHENOTYPE']==str(phenotype)) & 
                            (threshold_data['SUBSAMPLE_NUM']==int(subsample_number)) & 
                            (threshold_data['PVAL_TYPE']=="AVERAGE_P") &
                            (threshold_data['THRESHOLD_TYPE']=="NNBF")
                            ]
            
            # obtain threshold value specifically
            subsetted_threshold_value = subsetted_threshold_data.iloc[0]["THRESHOLD_VALUE"]

        elif method_type == "GIFT":
            average_column_name="AVERAGE_PSNP8"

            # subset the threshold data to get the threshold at current subsamples:
            subsetted_threshold_data=threshold_data[(threshold_data['PHENOTYPE']==str(phenotype)) & 
                            (threshold_data['SUBSAMPLE_NUM']==int(subsample_number)) & 
                            (threshold_data['PVAL_TYPE']=="AVERAGE_PSNP8") &
                            (threshold_data['THRESHOLD_TYPE']=="NNNH")
                            ]
            
            # obtain threshold value specifically
            subsetted_threshold_value = subsetted_threshold_data.iloc[0]["THRESHOLD_VALUE"]

        ### HIGHLIGHT OF SIGNIFICANT SNPS (remove annotate)
        # get data to annotate with
            # read file of SNPs for the given subsample number
            # {phenotype}_{GIFT/GWAS}_{subsample_num}_ALL.csv
            # update later to change based on phenotype
        R_out.write(f'ALL_SNPS<-read.csv("output_files/R_DATA/{phenotype}_{method_type}_{subsample_number}_ALL.csv",header=TRUE)\n')

        # log transform the appropriate column of data 
        R_out.write(f'ALL_SNPS${average_column_name} <- (-log10(ALL_SNPS${average_column_name}))\n')

        # filter to keep only those above the threshold
        R_out.write(f'SIGNIFICANT_SNPS <- ALL_SNPS %>% \n') 
        R_out.write(f'  filter({average_column_name}>={subsetted_threshold_value}) \n') 

        R_out.write(f' \n')
        R_out.write(f'my_data <- subsample_dataset_{method_type} %>% \n')
        R_out.write(f'  group_by(CHROM) %>%  \n')
        R_out.write(f'  summarise(chr_len=max(POS)) %>% \n')
        R_out.write(f'  mutate(tot=cumsum(chr_len)-chr_len) %>% \n')
        R_out.write(f'  select(-chr_len) %>% \n')
        R_out.write(f'  left_join(subsample_dataset_{method_type}, ., by=c("CHROM"="CHROM")) %>% \n')
        R_out.write(f'  arrange(CHROM,POS) %>% \n')
        R_out.write(f'  mutate( BPcum=POS+tot) \n')
        R_out.write(f' \n')
        R_out.write(f'axisdf = my_data %>% \n')
        R_out.write(f'  group_by(CHROM) %>% \n')
        R_out.write(f'  summarize(center=( max(BPcum) + min(BPcum) ) / 2 ) \n')
        R_out.write(f' \n')
        
        # make identifier columns for each dataframe 
        R_out.write(f'my_data<- my_data %>% \n')
        R_out.write(f'  mutate(identifier=paste(CHROM,POS,sep="_")) \n')
        R_out.write(f' \n')
        R_out.write(f'SIGNIFICANT_SNPS<- SIGNIFICANT_SNPS %>% \n')
        R_out.write(f'  mutate(identifier=paste(CHR,POS,sep="_")) \n')
        R_out.write(f' \n')

        # cross reference significant snps to the _ALL results to determine which need to be highlighted with black dots
        R_out.write(f'my_data<- my_data%>% \n')
        R_out.write(f'  mutate(is_highlight= identifier %in% SIGNIFICANT_SNPS$identifier) \n')

        R_out.write(f'png("{Plot_output_directory}{method_type}_{subsample_number}_MANHATTAN.png", bg = "white", width = 7.5, height = 2.5, units = "in", res = 1200, pointsize = 3) \n')
        R_out.write(f' \n')  
        R_out.write(f'ggplot(my_data, aes(x=BPcum, y=pvalue_stdv, color=as_factor(CHROM)))+ \n')
        R_out.write(f'  geom_point(alpha=0.5)+ \n')
        R_out.write(f'  scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center)+ \n')
        R_out.write(f'  labs(x = "position", color = "Chromosomes") + \n')
        # add highlighting overlay if there are any to highlight
        R_out.write(f'  geom_point(data = subset(my_data, is_highlight==TRUE), color="black", alpha=0.5) +\n')
        
        if subsample_number == "600":
            R_out.write(f'  ylim(0,7) \n')
        elif subsample_number == "400":
            R_out.write(f'  ylim(0,6) \n')
        else:         
            R_out.write(f'  ylim(0,5) \n')
        R_out.write(f'dev.off() \n')
        R_out.write(f' \n')

    R_out.write(f' \n')
    R_out.write(f' \n')
    R_out.write(f'# end of R script \n')
    R_out.close()

    # Create shell script to run the manhattan plots of standard deviation

    Shell_out=open(str(Shell_out_directory)+str(subsample_number)+"_stdv_plots.sh","w")
    # necessary start to the file
    Shell_out.write(f'#!/bin/bash\n')
    Shell_out.write(f'#SBATCH --partition=defq\n')
    Shell_out.write(f'#SBATCH --nodes=1\n')
    Shell_out.write(f'#SBATCH --ntasks=1\n')
    Shell_out.write(f'#SBATCH --cpus-per-task=3\n')
    Shell_out.write(f'#SBATCH --mem=8g\n')
    Shell_out.write(f'#SBATCH --time=1:00:00\n')
    Shell_out.write(f'#SBATCH --job-name=S4_SUB_MANHATTAN_{subsample_number}\n')
    Shell_out.write(f'#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out\n')
    Shell_out.write(f'#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err\n')
    Shell_out.write(f'#SBATCH --mail-type=ALL\n')
    Shell_out.write(f'#SBATCH --mail-user=mbysh17@nottingham.ac.uk\n')
    Shell_out.write(f'#===============================\n')
    Shell_out.write(f'echo "start SUBRUN_S4_{subsample_number}"\n')
    Shell_out.write(f'#change to home directory\n')
    Shell_out.write(f'cd /gpfs01/home/mbysh17\n')
    Shell_out.write(f'# source conda environments\n')
    Shell_out.write(f'source ~/.bashrc\n')
    Shell_out.write(f'conda deactivate\n')
    Shell_out.write(f'conda activate r_env\n')
    Shell_out.write(f'# Run the right R script \n')
    # str(R_script_out_directory)+str(subsample_number)+"_stdv_plots.R
    Shell_out.write(f'Rscript {R_script_out_directory}{subsample_number}_stdv_plots.R\n')
    Shell_out.write(f'conda deactivate\n')
    Shell_out.write(f'echo "END OF STAGE 4 SUBRUN"\n')
    Shell_out.write(f'# end of file\n')
    Shell_out.close()

# batch file to get box plots
    # (str(R_script_out_directory)+"boxplot_stdv_plots.R","w")
Shell_out=open(str(Shell_out_directory)+"boxplot_stdv_plots.sh","w")
# necessary start to the file
Shell_out.write(f'#!/bin/bash\n')
Shell_out.write(f'#SBATCH --partition=defq\n')
Shell_out.write(f'#SBATCH --nodes=1\n')
Shell_out.write(f'#SBATCH --ntasks=1\n')
Shell_out.write(f'#SBATCH --cpus-per-task=3\n')
Shell_out.write(f'#SBATCH --mem=8g\n')
Shell_out.write(f'#SBATCH --time=1:00:00\n')
Shell_out.write(f'#SBATCH --job-name=S4_SUB_boxplot\n')
Shell_out.write(f'#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out\n')
Shell_out.write(f'#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err\n')
Shell_out.write(f'#SBATCH --mail-type=ALL\n')
Shell_out.write(f'#SBATCH --mail-user=mbysh17@nottingham.ac.uk\n')
Shell_out.write(f'#===============================\n')
Shell_out.write(f'echo "start SUBRUN_S4_{subsample_number}"\n')
Shell_out.write(f'#change to home directory\n')
Shell_out.write(f'cd /gpfs01/home/mbysh17\n')
Shell_out.write(f'# source conda environments\n')
Shell_out.write(f'source ~/.bashrc\n')
Shell_out.write(f'conda deactivate\n')
Shell_out.write(f'conda activate r_env\n')
Shell_out.write(f'# Run the right R script \n')
# str(R_script_out_directory)+str(subsample_number)+"_stdv_plots.R
Shell_out.write(f'Rscript {R_script_out_directory}boxplot_stdv_plots.R\n')
Shell_out.write(f'conda deactivate\n')
Shell_out.write(f'echo "END OF STAGE 4 BOXPLOT"\n')
Shell_out.write(f'# end of file\n')
Shell_out.close()


print("stage_4.py script finished!",flush=True)

# end of script
