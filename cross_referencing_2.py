import pandas ,os 


PATH_TO_MAIN = "/gpfs01/home/mbysh17/"

GIFT_gene_results_files=[]
GWAS_gene_results_files=[]


Mo98_gene_tracker_df = pandas.DataFrame(columns=['DATA_SOURCE'])
Na23_gene_tracker_df = pandas.DataFrame(columns=['DATA_SOURCE'])


for file in os.listdir(PATH_TO_MAIN+"output_files/GENES_DATA"):
    if file.endswith(".txt") and file.__contains__("FINAL") and file.__contains__("GIFT")==True:
        GIFT_gene_results_files.append(file)
        print(file,": ++++++++++++ ADDED ++++++++++++ !", flush = True)
    else:
        pass



for file in os.listdir(PATH_TO_MAIN+"output_files/GENES_DATA"):
    if file.endswith(".txt") and file.__contains__("FINAL") and file.__contains__("GWAS")==True:
        GWAS_gene_results_files.append(file)
        print(file,": ++++++++++++ ADDED ++++++++++++ !", flush = True)
    else:
        pass


Mo98_genes_data = pandas.read_csv(f'{PATH_TO_MAIN}core_files/Mo_genes_data.csv')
Mo98_gene_list=Mo98_genes_data["Gene_ID"].tolist()

Na23_genes_data =  pandas.read_csv(f'{PATH_TO_MAIN}core_files/Na_genes_data.csv')
Na23_gene_list=Na23_genes_data["Gene_ID"].tolist()


print(Mo98_genes_data.head())

for gene_ID in Mo98_gene_list:

    # adds each Gene_ID of interest as a new column
    Mo98_gene_tracker_df[f'{gene_ID}'] = ''


for gene_ID in Na23_gene_list:

    # adds each Gene_ID of interest as a new column
    Na23_gene_tracker_df[f'{gene_ID}'] = ''



print(Mo98_gene_tracker_df)

current_Mo98_index =0  
current_Na23_index =0


GIFT_gene_results_files = sorted(GIFT_gene_results_files, reverse=True, key= lambda ele:int(ele.split("_")[5]))

GWAS_gene_results_files= sorted(GWAS_gene_results_files, reverse=True, key= lambda ele:int(ele.split("_")[5]))


for gene_list in GIFT_gene_results_files:

    print(f"Gene list: {gene_list}",flush=True)

    gene_list_name = gene_list.split("_")
    
    print(f"Gene list name: {gene_list_name}",flush=True)

    current_phenotype = gene_list_name[3]
    current_method = gene_list_name[4]
    current_subsample_num = gene_list_name[5]

    gene_list_name_modified = current_phenotype+"_"+current_method+"_"+current_subsample_num

    current_gene_list = pandas.DataFrame(columns=["Genes"])
    current_gene_list["Genes"] = pandas.read_csv(f'{PATH_TO_MAIN}output_files/GENES_DATA/{gene_list}', header=None)

    current_gene_list=current_gene_list["Genes"].tolist()    

    print(current_gene_list[0:4])

    if current_phenotype=="Mo98":
        # first append empty row
        Mo98_gene_tracker_df.loc[len(Mo98_gene_tracker_df)] = pandas.Series(dtype='string')

        # add in name of source to source column
        Mo98_gene_tracker_df.loc[current_Mo98_index,'DATA_SOURCE']=str(gene_list_name_modified)


        for gene_ID in Mo98_gene_list:

            if gene_ID in current_gene_list:

                #print(f"Found gene: {gene_ID} in file: {gene_list}!")

                # add in info if the gene ID was found in the file or not
                Mo98_gene_tracker_df.loc[current_Mo98_index,gene_ID]="YES"
            else:
                Mo98_gene_tracker_df.loc[current_Mo98_index,gene_ID]="NO"
        
        current_Mo98_index+=1

        print("Head of Mo98 dataframe tracker")

        print(Mo98_gene_tracker_df,flush=True)


    elif current_phenotype == "Na23":
        # first append empty row
        Na23_gene_tracker_df.loc[len(Na23_gene_tracker_df)] = pandas.Series(dtype='string')



        # add in name of source to source column
        Na23_gene_tracker_df.loc[current_Na23_index,'DATA_SOURCE']=str(gene_list_name_modified)


        for gene_ID in Na23_gene_list:

            if gene_ID in current_gene_list:

                #print(f"Found gene: {gene_ID} in file: {gene_list}!")

                # add in info if the gene ID was found in the file or not
                Na23_gene_tracker_df.loc[current_Na23_index,gene_ID]="YES"
            else:
                Na23_gene_tracker_df.loc[current_Na23_index,gene_ID]="NO"

        current_Na23_index+=1


for gene_list in GWAS_gene_results_files:

    print(f"Gene list: {gene_list}",flush=True)

    gene_list_name = gene_list.split("_")
    
    print(f"Gene list name: {gene_list_name}",flush=True)

    current_phenotype = gene_list_name[3]
    current_method = gene_list_name[4]
    current_subsample_num = gene_list_name[5]

    gene_list_name_modified = current_phenotype+"_"+current_method+"_"+current_subsample_num

    current_gene_list = pandas.DataFrame(columns=["Genes"])
    current_gene_list["Genes"] = pandas.read_csv(f'{PATH_TO_MAIN}output_files/GENES_DATA/{gene_list}', header=None)

    current_gene_list=current_gene_list["Genes"].tolist()    

    print(current_gene_list[0:4])

    if current_phenotype=="Mo98":
        # first append empty row
        Mo98_gene_tracker_df.loc[len(Mo98_gene_tracker_df)] = pandas.Series(dtype='string')

        # add in name of source to source column
        Mo98_gene_tracker_df.loc[current_Mo98_index,'DATA_SOURCE']=str(gene_list_name_modified)


        for gene_ID in Mo98_gene_list:

            if gene_ID in current_gene_list:

                #print(f"Found gene: {gene_ID} in file: {gene_list}!")

                # add in info if the gene ID was found in the file or not
                Mo98_gene_tracker_df.loc[current_Mo98_index,gene_ID]="YES"
            else:
                Mo98_gene_tracker_df.loc[current_Mo98_index,gene_ID]="NO"
        
        current_Mo98_index+=1

        print("Head of Mo98 dataframe tracker")

        print(Mo98_gene_tracker_df,flush=True)


    elif current_phenotype == "Na23":
        # first append empty row
        Na23_gene_tracker_df.loc[len(Na23_gene_tracker_df)] = pandas.Series(dtype='string')



        # add in name of source to source column
        Na23_gene_tracker_df.loc[current_Na23_index,'DATA_SOURCE']=str(gene_list_name_modified)


        for gene_ID in Na23_gene_list:

            if gene_ID in current_gene_list:

                #print(f"Found gene: {gene_ID} in file: {gene_list}!")

                # add in info if the gene ID was found in the file or not
                Na23_gene_tracker_df.loc[current_Na23_index,gene_ID]="YES"
            else:
                Na23_gene_tracker_df.loc[current_Na23_index,gene_ID]="NO"

        current_Na23_index+=1



# add gene names to the headers
for row in range(0,len(Mo98_genes_data)):

    current_ID = str(Mo98_genes_data["Gene_ID"][row])
    new_header = str(current_ID) + "(" + str(Mo98_genes_data["Gene_Name"][row]) + ")"

    Mo98_gene_tracker_df = Mo98_gene_tracker_df.rename(columns = {current_ID:new_header})

for row in range(0,len(Na23_gene_list)):

    current_ID = str(Na23_genes_data["Gene_ID"][row])
    new_header = str(current_ID) + "(" + str(Na23_genes_data["Gene_Name"][row]) + ")"
    
    Na23_gene_tracker_df = Na23_gene_tracker_df.rename(columns = {current_ID:new_header})

# save both tracker dataframes to csv in the GENES_DATA file

Mo98_gene_tracker_df.to_csv(f"{PATH_TO_MAIN}output_files/GENES_DATA/Mo98_Gene_Tracker.csv",header=True,index=False)
Na23_gene_tracker_df.to_csv(f"{PATH_TO_MAIN}output_files/GENES_DATA/Na23_Gene_Tracker.csv",header=True,index=False)

print("End of cross_referencing_2 python script",flush=True)

