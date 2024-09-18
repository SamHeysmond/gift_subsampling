#import modin.pandas as pandas
#import ray
import pandas

# set ray memory to 20GB
#ray.init(_plasma_directory="/tmp", object_store_memory=20000000000)

pandas.set_option('display.max_columns',None)
pandas.options.display.max_columns=None

PATH_TO_MAIN = "/gpfs01/home/mbysh17/"

phenotype_df = pandas.read_csv(PATH_TO_MAIN+"core_files/selected_Ion_data.csv",index_col=[0])

print("Start",flush=True)
print(phenotype_df,flush=True)


# Filter out unwanted columns
drop_cols=[1,2,3,4,5,7]

phenotype_df.drop(phenotype_df.columns[drop_cols],axis=1,inplace=True)

print("After drop cols",flush=True)
print(phenotype_df,flush=True)

#Filter out seed ionome data
leaf_df = phenotype_df[phenotype_df['material'].str.contains("leaf")]
print("After seed filter",flush=True)
print(leaf_df,flush=True)

leaf_df.to_csv(PATH_TO_MAIN+"core_files/leaf_phenotype.csv",header=True,index=False)

#Filter out for specific metals
Mo98_df = leaf_df[['Accession_ID','Mo98']]
Na23_df = leaf_df[['Accession_ID','Na23']]

print("After metal filter",flush=True)
print(Mo98_df,flush=True)
print(Na23_df,flush=True)

Mo98_df.to_csv(PATH_TO_MAIN+"core_files/Mo98_phenotype.csv",header=True,index=False)
Na23_df.to_csv(PATH_TO_MAIN+"core_files/Na23_phenotype.csv",header=True,index=False)

#end of script
