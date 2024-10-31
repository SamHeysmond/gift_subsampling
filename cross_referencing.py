import pandas
import os
import math

PATH_TO_MAIN = "/gpfs01/home/mbysh17/"

csv_files_Arabidopsis_Thaliana=[]

def log_function(value):
    return -(math.log(value,10))

def add_to_end(data):
    return (data + 1000)

def take_from_start(data):
    return (data-1000)

# add all the compiled files to a list (ending in _ALL.csv)
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
my_batch.write(f'#SBATCH --partition=shortq\n')
my_batch.write(f'#SBATCH --nodes=1\n')
my_batch.write(f'#SBATCH --ntasks=1\n')
my_batch.write(f'#SBATCH --cpus-per-task=2\n')
my_batch.write(f'#SBATCH --mem=8g\n')
my_batch.write(f'#SBATCH --time=02:00:00\n')
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



    if current_method=="GWAS":

        subsetted_threshold_data=threshold_data[(threshold_data['PHENOTYPE']==current_phenotype) & 
                                  (threshold_data['SUBSAMPLE_NUM']==current_subsample_num) & 
                                  (threshold_data['PVAL_TYPE']=="AVERAGE_P") &
                                  (threshold_data['THRESHOLD_TYPE']=="NNBF")
                                  ]


        Threshold_Value = subsetted_threshold_data.iloc[0]["THRESHOLD_VALUE"]

        # if float(Threshold_Value)==0.0:
        #     print("Erroneous threshold encountered, setting to 9999",flush=True)
        #     Threshold_Value = 9999

        
        print("Threshold value)",flush=True)
        print(Threshold_Value,flush=True)
         
        plant_data["AVERAGE_P"] = plant_data["AVERAGE_P"].apply(log_function)

        plant_data=plant_data.loc[plant_data['AVERAGE_P']>=Threshold_Value]

        print("GWAS data before sort",flush=True)
        print(plant_data.head(10))

        # sort in descending order (so biggest are first after log transformation)
        plant_data=plant_data.sort_values(by="AVERAGE_P",ascending=False)
        print("GWAS data after sort",flush=True)
        print(plant_data.head(10))

        # now take the top 200 SNPs (temp disabled)
        # cropped_plant_data=plant_data.head(200)

        # print("(Head of cropped plant data)",flush=True)
        # print(cropped_plant_data.head(10))

    elif current_method == "GIFT":

        print("Threshold data before filter",flush=True)
        print(threshold_data,flush=True)

        subsetted_threshold_data=threshold_data[(threshold_data['PHENOTYPE']==current_phenotype) & 
                                (threshold_data['SUBSAMPLE_NUM']==current_subsample_num) & 
                                (threshold_data['PVAL_TYPE']=="AVERAGE_PSNP8") &
                                (threshold_data['THRESHOLD_TYPE']=="NNNH")
                                ]
        print("Threshold data after filter",flush=True)
        print(subsetted_threshold_data,flush=True)

        Threshold_Value = subsetted_threshold_data.iloc[0]["THRESHOLD_VALUE"]

        # if float(Threshold_Value)==0.0:
        #     print("Erroneous threshold encountered, setting to 9999",flush=True)
        #     Threshold_Value = 9999

        print("Threshold value)",flush=True)
        print(Threshold_Value,flush=True)


        plant_data["AVERAGE_PSNP8"] = plant_data["AVERAGE_PSNP8"].apply(log_function)


        plant_data=plant_data.loc[plant_data['AVERAGE_PSNP8']>=Threshold_Value]
        
        print("PSNP8 data before sort",flush=True)
        print(plant_data.head(10))

        # sort in descending order (so biggest are first after log transformation)
        plant_data=plant_data.sort_values(by="AVERAGE_PSNP8",ascending=False)
        print("PSNP8 data after sort",flush=True)
        print(plant_data.head(10))

        # now take the top 200 SNPs (temp disabled)
        # cropped_plant_data=plant_data.head(200)
   
        # print("(Head of cropped plant data)",flush=True)
        # print(plant_data.head(10))

    # rename headers if necessary
    #plant_data.rename(columns={""},inplace=True)

    # cut down the data after filtering
    plant_data=plant_data[['CHR','POS']]

    # duplicate the position so we haev a "start" and "end" position for overlap searching
        # maybe could add 1000 base pairs either side here?
        # if a number was less than 0 though, would need to be set to 0 e.g. -545 -> 0
    plant_data["END"] = plant_data["POS"]



    print("(Head of plant data before +/-1000 changes)",flush=True)
    print(plant_data.head(10))

    # trying to add 1000 to start and end here

    # RENAME columns to fit bed format
    plant_data.rename(columns={'CHR':'chromosome','POS':'START'},inplace=True)

    plant_data["END"]=plant_data["END"].apply(add_to_end)
    plant_data["START"]=plant_data["START"].apply(take_from_start)

    print("(Head of plant data after +/-1000 changes)",flush=True)
    print(plant_data.head(10))
    
    #add "Chr" to the first column
    plant_data['chromosome']="Chr"+plant_data['chromosome'].astype(str)

    # write to bed file
    bedfile_name=PATH_TO_MAIN+str("output_files/GENES_DATA/")+csv_file.replace(".csv",".bed")
    intersect_result_file_name=PATH_TO_MAIN+str("output_files/GENES_DATA/Intersect_results_")+csv_file.replace(".csv",".txt")
    final_intersect_result_file_name=PATH_TO_MAIN+str("output_files/GENES_DATA/FINAL_Intersect_results_")+csv_file.replace(".csv",".txt")
    plant_data.to_csv(bedfile_name,index=False,header=False,sep='\t')
    plant_data_gff3_location=str(PATH_TO_MAIN+"core_files/TAIR10_GFF3_genes.gff")


    # careful with -b since this file is the one loaded into memory (so if mem overloaded swap the files around)
    my_batch.write(f'bedtools intersect -a {bedfile_name} -b {plant_data_gff3_location} -wb | grep "gene" -w | sort --unique 1> {intersect_result_file_name}\n')
    
    my_batch.write("awk -F '[=;]' '{print $2}' "+intersect_result_file_name+" | sort --unique 1>  "+final_intersect_result_file_name+"\n")
    
    #my_batch.write('echo "EndOfFile" >>' + final_intersect_result_file_name+' \n')


    my_batch.write(f'\n')

    print("End of for loop:",flush=True)

my_batch.write(f'conda deactivate \n')
my_batch.write(f'echo "END OF cross_reference_script" \n')

my_batch.write(f'echo "greenlight" > core_files/light.txt \n')
my_batch.write(f'# end of script')
my_batch.close()

print("End of cross_referencing python script",flush=True)

# write red light to file
f=open(f"{PATH_TO_MAIN}/core_files/light.txt","w")
f.write("redlight")
f.close()

# End of script