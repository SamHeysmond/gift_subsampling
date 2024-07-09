# GIFT vs GWAS subsampling and testing
This README.md file gives a summary of what each of the files on the github contains. Any references for local files are also detailed in their description. 
For referencing and information on external software, websites etc, please see "External_References_And_Info.md"
# Files and what they contain
> ## Workflow_Timeline.md
+ Introduction
   + Briefly introduces problem at hand 
   + Aims of analysis
   + Expected outcomes (for each stage)
+ Data
    + Input and output data for each stage of analysis
    + Any details about this data that are needed to be known e.g. if input or output needs naming format etc
+ Code
  + The steps where code was used is denoted in this file along with any additional files or requirements.
+ References
    + To data that is needed (references of software and papers etc are in Referenecs_and_info.md)
---
# SCRIPTS : PREP

> ## 1 conda_environment_setup.sh
+ Installation and setup of all environments used
    + Dependendcies
    + Version numbers
+ Author: Sam Heysmond

# SCRIPTS : STAGE 1

> ## 2 subsample_setup.sh
+ Sets up folders and files needed for analysis
+ Removes previous files in those locations (refreshing the data)
  + (This can be disabled by commenting out the rm commands and echo command for JOB_LIST.csv)
+ Runs: subsample_script_setup.py (SCRIPT 3)
+ STAGE: 1
+ Author: Sam Heysmond

> ## 3 subsample_script_setup.py
+ Creates all the batch files which will subsample the vcf data and perform GIFT and GWAS
  + STAGE: 1
  + Creates: SCRIPT 6
+ Author: Sam Heysmond

> ## 4 subsample.py
+ The script used to subsample the vcf and phenotype file ID's to match the subsample number and be present in both genotype and phenotype files.
+ Author: Sam Heysmond

> ## 4 run_and_monitor.sh
+ Runs: SCRIPT 5 (run_and_monitor.py)
  + this is required to run to manage the batch scripts for the entirety of the subsampling and testing of GIFT and GWAS. So, give this file max run time.
+ Author: Sam Heysmond

> ## 5 run_and_monitor.py
+ Manages all of the batch scripts in "batch_files/parallel", running as many in parallel as possible. When the batch file has been put in the que, it moves it to a different folder "batch_files/completed_parallel"
+ Runs: SCRIPT 6 (subrun_{phenotype}_{subsampleNum}\_{copynum}.sh)
+ Author: Sam Heysmond

> ## 6 subrun_{phenotype}_{subsampleNum}\_{copynum}.sh
+ Self-made script (made from SCRIPT 3)
+ Each of these scripts runs stage 1 of analysis which is the running of GWAS and GIFT on subsampled data.
+ Currently doesnt run SCRIPT 10 ({slurm_job_id}_{i}_{phenotype}.R) because GWAS filtering needs to happen first
+ Runs:  SCRIPT 4 (SUBSAMPLE.PY), script 7 (PHYSICS_GWAS_OOP.py), SCRIPT 8 (pygwas) ,SCRIPT 9 (make_r_scripts)
+ Author: Sam Heysmond

> ## 7 physics_GWAS_OOP.py
+ The core script for GIFT analysis, this script will run a GIFT analysis on one set of input data per run.
+ This script has an inbuilt R script maker for the results so it also creates manhattan plot graphics.
+ The code has been updated to calculate and plot top 20 SNP data for each type of SNP E.G. PSNP4
+ Authors: Assistant Prof S Bray, Prof J Wattis, Prof C Rauch
+ Edited by: Sam Heysmond

> ## 8 pygwas
+ The core script for GWAS analysis, this script will run a GWAS analysis on one set of input data per run.
+ Reference: PENDING

> ## 9 make_r_scripts.py
+ Creates: SCRIPT 10 ({slurm_job_id}_{i}_{phenotype}.R)
+ Author: Sam Heysmond

> ## 10 {slurm_job_id}_{i}_{phenotype}.R
+ The individual R scripts that make manhattan plots for the GWAS data (since GIFT code includes the creation of an R script).
+ Calculates and plots the top 20 SNPs for the current dataset
+ ((THIS SCRIPT WILL NEED RUNNING LATER AFTER FILTERING GWAS DATA))
+ Author: Sam Heysmond


##### FILTER STUFF GOES HERE
> ## 13
+ The
+ Author: 

# SCRIPTS : STAGE 2

#### will need to re-run all the below data

> ## 11 SNP_tracker.sh
+ The main shell script which runs stage 2 of analysis which is gathering data on the SNPs and creating summary plots to compare both GIFT and GWAS results.
+ Runs: SCRIPT 12 (SNP_tracker.py)
+ Author: Sam Heysmond
  
> ## 12 SNP_tracker.py
+ The python script which concatonates the GWAS and GIFT data and processes it to satisfied summary plot/IDEAS 1-3 as described in the workflow timeline.
+ Author: Sam Heysmond

> ## 13 SNP_tracker_R_and_BASH_maker.sh
+ Runs: SCRIPT 14 (SNP_tracker_R_and_BASH_maker.py), SCRIPT 17 (Rscript_run_and_monitor.py)
+ Author: 

> ## 14 SNP_tracker_R_and_BASH_maker.py
+ Creates:  SCRIPT 15({phenotype}_cumulative_t20_dataframe.sh) SCRIPT 16 ({phenotype}_cumulative_t20_dataframe.R), SCRIPT 17 ({phenotype}\_{positive/negative}\_control.sh), SCRIPT 18 ({phenotype}_{positive/negative}_control.R), SCRIPT 19 ({phenotype}\_{subsample_number}_{pval_type}_MANHATTAN.sh),SCRIPT 20 ({phenotype}\_{subsample_number}_{pval_type}_MANHATTAN.R)
+ Author: 

> ## 15 {phenotype}_cumulative_t20_dataframe.sh
+ Runs: SCRIPT 16 ({phenotype}_cumulative_t20_dataframe.R)
+ Author: 

> ## 16 {phenotype}_cumulative_t20_dataframe.R
+ Top 20 SNPs at 1000 subsample for each SNP type (in both GWAS and GIFT) tracked through 1000 to 200 subsample level.
+ Author: 
  
> ## 17 {phenotype}_{positive/negative}_control.sh
+ Runs: SCRIPT 18 ({phenotype}_{positive/negative}_control.R) 
+ Author: 

> ## 18 {phenotype}_{positive/negative}_control.R
+ R script for IDEA 2

> ## 19 {phenotype}_{subsample_number}_{pval_type}_MANHATTAN.sh
+ Runs: SCRIPT 20 ({phenotype}_{subsample_number}_{pval_type}_MANHATTAN.R)
+ Author: 
  
> ## 20 {phenotype}_{subsample_number}_{pval_type}_MANHATTAN.R
+ The
+ Author: 

> ## 13
+ The

> ## 21 Rscript_run_and_monitor.py
+ Runs: SCRIPT 15 ({phenotype}_cumulative_t20_dataframe.sh), SCRIPT 17 ({phenotype}\_{positive/negative}\_control.sh) , SCRIPT 19 ({phenotype}_{subsample_number}_{pval_type}_MANHATTAN.sh)
+ Author: 


+ Author: 
  
> ## 13
+ The
+ Author: 
  
> ## 13
+ The
+ Author: 

> ## 13
+ The
+ Author: 


# SCRIPTS : STAGE 3
  




> ## 13
+ The
+ Author: 
  
> ## 13
+ The
+ Author: 

> ## 13
+ The
+ Author: 
  
> ## 13
+ The
+ Author: 
  
> ## 13
+ The
+ Author: 

> ## 13
+ The
+ Author: 

# DATA FILES (INPUT)
> ## 1 SNP_Matrix
+ (Get info on this)
+ Used as input for the GWAS software "pygwas"
+ Reference: PENDING

> ## 2 1001genomes_snp_biallelic_only_ACGTN.vcf
+ Genotype information from the 1001 genomes project on biallelic SNPs data from Arabidopsis Thaliana
+ Reference: PENDING

> ## 3 master_list.csv
+ Leaf and seed ionome phenotype data from 1001 Genomes database on Arabidopsis Thaliana
+ Reference: PENDING

> ## 4 phenotypes_list.txt
+ List of which phenotype columns to analyse. The scripts will use this to know what phenotypes to look into.
+ Reference: N/A

# DATA FILES (OUTPUT)
> ## 1 JOB_LIST.csv
+ Output list of all subrun jobs for stage 1 of analysis. This is partly designed for error checking to ensure all jobs were properly completed, but in future versions of scripts it may be used to help with analysis and segregation of data.
+ Reference: N/A

> ## 2 all_vcf_samples.txt
+ List of all IDs of samples from the vcf folder. Some scripts may use this data to check that the vcf and master list (phenotypes) contain the same sample IDs before analysis can begin.
Reference: N/A

> ## 3 subsamples_{subsample_num}_{JOB_ID}.txt
+ List of IDs of all the subsamples for a particular job run. This is useful for error checking but is also used in the program to compare between the VCF and phenotype data for the sample IDs that have been selected.
---
> ## References_and_info.md
  + References
      + Scientific papers describing the tool/algorithm
      + Github links to tools
      + Scripts and software that were used
---