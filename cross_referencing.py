import pandas
import os
import math

PATH_TO_MAIN = "/gpfs01/home/mbysh17/"

# add all excel files to lists
# csv_files_Arabidopsis_Thaliana=[
#     "Mo98_GIFT_1000_ALL.csv",
#     "Mo98_GIFT_200_ALL.csv",
#     "Na23_GIFT_1000_ALL.csv",
#     "Na23_GIFT_200_ALL.csv",
#     "Mo98_GWAS_1000_ALL.csv",
#     "Mo98_GWAS_200_ALL.csv",
#     "Na23_GWAS_1000_ALL.csv",
#     "Na23_GWAS_200_ALL.csv"
#     ]


csv_files_Arabidopsis_Thaliana=[]

for file in os.listdir(PATH_TO_MAIN+"output_files/R_DATA/"):
    if file.endswith(".csv") and file.__contains__("ALL")==True:
        csv_files_Arabidopsis_Thaliana.append(file)
        print(file,": ++++++++++++ ADDED ++++++++++++ !", flush = True)
    else:
        pass


#Begin writing the batch script that will do the cross referencing

 # location where the batch scripts will be written to
my_batch=open(PATH_TO_MAIN+"batch_files/cross_reference_script.sh","w")

# necessary start to the file
my_batch.write(f'#!/bin/bash\n')
my_batch.write(f'#SBATCH --partition=defq\n')
my_batch.write(f'#SBATCH --nodes=1\n')
my_batch.write(f'#SBATCH --ntasks=1\n')
my_batch.write(f'#SBATCH --cpus-per-task=4\n')
my_batch.write(f'#SBATCH --mem=10g\n')
my_batch.write(f'#SBATCH --time=24:00:00\n')
my_batch.write(f'#SBATCH --job-name=cross_reference_script\n')
my_batch.write(f'#SBATCH --output=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.out\n')
my_batch.write(f'#SBATCH --error=/gpfs01/home/mbysh17/slurmOandE/slurm-%x-%j.err\n')
my_batch.write(f'#SBATCH --mail-type=ALL\n')
my_batch.write(f'#SBATCH --mail-user=mbysh17@nottingham.ac.uk\n')
my_batch.write(f'#===============================\n')
my_batch.write(f'echo "start cross_reference_script" \n')
my_batch.write(f'#change to home directory\n')
my_batch.write(f'cd /gpfs01/home/mbysh17\n')
my_batch.write(f'# source conda environments\n')
my_batch.write(f'source ~/.bashrc\n')
my_batch.write(f'echo "START OF crossreferencing script" \n')
my_batch.write(f'# activate env\n')
my_batch.write(f'conda activate bedtools_env\n')
# read in the threshold file
threshold_data =pandas.read_csv(f'{PATH_TO_MAIN}output_files/R_DATA/THRESHOLDS.csv')

for csv_file in csv_files_Arabidopsis_Thaliana:

    print("Start of for loop:",flush=True)

    csv_file_name = csv_file.split("_")
    if csv_file_name[1]=="GIFT":
        current_method = "GIFT"
    elif csv_file_name[1]=="GWAS":
        current_method = "GWAS"

    current_phenotype=csv_file_name[0]
    print(f"current phenotype: {current_phenotype}")

    current_subsample_num=int(csv_file_name[2])

    print(f"current subsample_num: {current_subsample_num}")

    plant_data=pandas.read_csv(f'{PATH_TO_MAIN}output_files/R_DATA/{csv_file}')

    # filter the data based on the right threshold (subsample level, phenotype, pval type)
        # i will use 1000, 200 Mo98, Na23, AVERAGE_P, PSNP5 for now

    def log_function(value):
        return -(math.log(value,10))

    if current_method=="GWAS":

        print("Threshold data before filter",flush=True)
        print(threshold_data,flush=True)

        BHY_thresh_data=threshold_data[(threshold_data['PHENOTYPE']==current_phenotype) & 
                                  (threshold_data['SUBSAMPLE_NUM']==current_subsample_num) & 
                                  (threshold_data['PVAL_TYPE']=="AVERAGE_P") &
                                  (threshold_data['THRESHOLD_TYPE']=="BHY")
                                  ]

        print("Threshold data after filter",flush=True)
        print(BHY_thresh_data,flush=True)

        BHY_thresh_value = BHY_thresh_data.iloc[0]["THRESHOLD_VALUE"]

        print("Threshold value)",flush=True)
        print(BHY_thresh_value,flush=True)
         
        print("=====================================",flush=True)
        print("Head of GWAS databaste (top 10) before log10 ",flush=True)
        print(plant_data.head(10))


        plant_data["AVERAGE_P"] = plant_data["AVERAGE_P"].apply(log_function)


        print("=====================================",flush=True)
        print("Head of GWAS databaste (top 10) after log10 ",flush=True)
        print(plant_data.head(10))

        plant_data=plant_data.loc[plant_data['AVERAGE_P']>=BHY_thresh_value]

        print("=====================================",flush=True)
        print("Head of GWAS databaste (top 10) after threshold filtering",flush=True)
        print(plant_data.head(10))

        # sort in descending order (so biggest are first after log transformation)
        plant_data=plant_data.sort_values(by="AVERAGE_P",ascending=False)

        print("=====================================",flush=True)
        print("Head of GIFT databaste (top 10) after sorting",flush=True)
        print(plant_data.head(10))

    elif current_method == "GIFT":

        print("Threshold data before filter",flush=True)
        print(threshold_data,flush=True)

        BHY_thresh_data=threshold_data[(threshold_data['PHENOTYPE']==current_phenotype) & 
                                (threshold_data['SUBSAMPLE_NUM']==current_subsample_num) & 
                                (threshold_data['PVAL_TYPE']=="AVERAGE_PSNP5") &
                                (threshold_data['THRESHOLD_TYPE']=="BHY")
                                ]
        print("Threshold data after filter",flush=True)
        print(BHY_thresh_data,flush=True)

        BHY_thresh_value = BHY_thresh_data.iloc[0]["THRESHOLD_VALUE"]

        print("Threshold value)",flush=True)
        print(BHY_thresh_value,flush=True)

        print("=====================================",flush=True)
        print("Head of GIFT databaste (top 10) before log10 ",flush=True)
        print(plant_data.head(10))

        plant_data["AVERAGE_PSNP5"] = plant_data["AVERAGE_PSNP5"].apply(log_function)

        print("=====================================",flush=True)
        print("Head of GIFT databaste (top 10) after log10 ",flush=True)
        print(plant_data.head(10))

        plant_data=plant_data.loc[plant_data['AVERAGE_PSNP5']>=BHY_thresh_value]

        print("=====================================",flush=True)
        print("Head of GIFT databaste (top 10) after threshold filtering",flush=True)
        print(plant_data.head(10))

        # sort in descending order (so biggest are first after log transformation)
        plant_data=plant_data.sort_values(by="AVERAGE_PSNP5",ascending=False)

        print("=====================================",flush=True)
        print("Head of GIFT databaste (top 10) after sorting",flush=True)
        print(plant_data.head(10))

   
    # rename headers if necessary
    #plant_data.rename(columns={""},inplace=True)

    # cut down the data after filtering
    plant_data=plant_data[['CHR','POS']]

    # duplicate the position so we haev a "start" and "end" position for overlap searching
    plant_data["END"] = plant_data["POS"]

    # RENAME columns to fit bed format
    plant_data.rename(columns={'CHR':'chromosome','POS':'START'},inplace=True)


    #add "Chr" to the first column

    plant_data['chromosome']="Chr"+plant_data['chromosome'].astype(str)

    # account for linkage disequlibrium by expanding a region around each SNP
        # NEED TO CONFIRM LD

    # plant_data['START']=plant_data['START'] - 100000

    # plant_data['END']=plant_data['END'] + 100000

    # write to bed file
    bedfile_name=PATH_TO_MAIN+str("output_files/GENES_DATA/")+csv_file.replace(".csv",".bed")
    intersect_result_file_name=PATH_TO_MAIN+str("output_files/GENES_DATA/Intersect_results_")+csv_file.replace(".csv",".txt")
    final_intersect_result_file_name=PATH_TO_MAIN+str("output_files/GENES_DATA/FINAL_Intersect_results_")+csv_file.replace(".csv",".txt")
    plant_data.to_csv(bedfile_name,index=False,header=False,sep='\t')
    plant_data_gff3_location=str(PATH_TO_MAIN+"core_files/TAIR10_GFF3_genes.gff")


    # careful with -b since this file is the one loaded into memory (so if mem overloaded swap the files around)
    my_batch.write(f'bedtools intersect -a {bedfile_name} -b {plant_data_gff3_location} -wb | grep "gene" -w | sort --unique 1> {intersect_result_file_name}\n')
    
    # for only genes that have been named
    #my_batch.write("awk -F\'\"\' \'{if ($5==\"; gene_name \") print $6}\' "+intersect_result_file_name+" | sort -u > gene_names_only"+intersect_result_file_name+"\n")
    
    # backslash="\\"
    # # for genes with gene ids (including those unnamed genes)
    # my_batch.write("awk \'{if ($6==\"gene\") print $13}\' "+intersect_result_file_name+" | sort -u | tr -d "+backslash+"\"\; > gene_IDs_only"+intersect_result_file_name+"\n")

     # for genes with gene ids (including those unnamed genes)
    #my_batch.write("awk '{if ($3==\"gene\") print $9}' "+intersect_result_file_name+" | awk -F \'[=;]\' \'{print $2}\' | sort -u > "+intersect_result_file_name+"\n")
    #my_batch.write("awk '{if ($3==\"gene\") print $9}' "+intersect_result_file_name+" | sort -u > "+intersect_result_file_name+"\n")
    
    # works byut needs filtering
    #my_batch.write("awk '{print $12}' "+intersect_result_file_name+" | sort -u > "+intersect_result_file_name+"\n")

    #my_batch.write("awk -F '[=;]' '{print $2}' "+intersect_result_file_name+"  > "+final_intersect_result_file_name+"\n")

    my_batch.write("awk -F '[=;]' '{print $2}' "+intersect_result_file_name+" | sort --unique 1>  "+final_intersect_result_file_name+"\n")

    my_batch.write(f'\n')

    print("End of for loop:",flush=True)

my_batch.write(f'conda deactivate \n')
my_batch.write(f'echo "END OF cross_reference_script" \n')
my_batch.write(f'# end of script')
my_batch.close()

print("End of cross_referencing python script",flush=True)
# End of script

# this works
#bedtools intersect -a Sheep_MapResults_rlnaB12.bed -b Ovis_aries.Oar_v3.1.112.chr.gff3 -wb | grep "ID=gene" | sort --unique 1> latest_test.txt

# latest command
# bedtools intersect -a Equine_MapResults_rHeight.bed -b Equus_caballus.EquCab3.0.112.gtf -wb | grep "gene" -w | sort --unique 1> latest_test.txt

# THIS WORKS to get gene names only
# awk -F'"' '{if ($5=="; gene_name ") print $6}' latest_test.txt | sort -u > gene_names_only.txt

# THIS WORKS to get gene IDs only
# awk '{if ($6=="gene") print $13}' intersect_resultsMapResults_Height.txt | sort -u | tr -d \"\;


# refine
# cat latest_test.txt | grep "gene_name .*[A-Za-z0-9]." -o > only_genes.txt

## get ONLY the top X

## transform to bed format ()


## write to bed file (must be tab delimited and contain CHR, START, END  columns)

# # pip install python-calamine
# horse_residual_height_data=pandas.read_excel("test_env/MapResults_rHeight.xls", index_col=None, engine="calamine")
# print(horse_residual_height_data.head())

# horse_residual_height_data['mlog10_pSNP8'] = horse_residual_height_data['mlog10_pSNP8'].fillna(9999)

# horse_residual_height_data = horse_residual_height_data.sort_values(by=['mlog10_pSNP8'],ascending=True)
# print("Head")
# print(horse_residual_height_data.head(100))

# pos_top=50
# pos_bot=70
# print(f'Viewing positions {pos_top} - {pos_bot}')
# print(horse_residual_height_data.iloc[pos_top:pos_bot,])

# print("Tail")
# horse_residual_height_data = horse_residual_height_data.sort_values(by=['mlog10_pSNP8'],ascending=False)
# print(horse_residual_height_data.head(100))

