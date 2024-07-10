# GIFT vs GWAS subsampling and testing
This README.md file gives a summary of what each of the files on the github contains. Any references for local files are also detailed in their description. 
For referencing and information on external software, websites etc, please see "External_References_And_Info.md"
# Files and what they contain / do
> ## Workflow_Timeline.md
+ Introduction
   + Briefly introduces problem at hand 
   + Aims of analysis
   + Expected outcomes (for each stage)
+ Data info
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

> ## 2 filtering_vcf.sh 
+ Creates a table "core_files/output_3.table" which contains all SNPs that need to be taken out of GWAS data to match the data of GIFT
+ Author: Sam Heysmond

# SCRIPTS : STAGE 1

> ## 3 subsample_setup.sh
+ Sets up folders and files needed for analysis
+ Removes previous files in those locations (refreshing the data)
  + (This can be disabled by commenting out the rm commands and echo command for JOB_LIST.csv)
+ Runs: subsample_script_setup.py (SCRIPT 4)
+ Author: Sam Heysmond

> ## 4 subsample_script_setup.py
+ Creates all the batch files which will subsample the vcf data and perform GIFT and GWAS
  + Creates: SCRIPT 8 (subrun_{phenotype}_{subsampleNum}\_{copynum}.sh)
+ Author: Sam Heysmond

> ## 5 subsample.py
+ The script used to subsample the vcf and phenotype file ID's to match the subsample number and be present in both genotype and phenotype files.
+ Author: Sam Heysmond

> ## 6 run_and_monitor.sh
+ Runs: SCRIPT 7 (run_and_monitor.py)
  + this is required to run to manage the batch scripts for the entirety of the subsampling and testing of GIFT and GWAS. So, give this file max run time.
+ Author: Sam Heysmond

> ## 7 run_and_monitor.py
+ Manages all of the batch scripts in "batch_files/parallel", running as many in parallel as possible. When the batch file has been put in the que, it moves it to a different folder "batch_files/completed_parallel"
+ Runs: SCRIPT 8 (subrun_{phenotype}_{subsampleNum}\_{copynum}.sh)
+ Author: Sam Heysmond

> ## 8 subrun_{phenotype}_{subsampleNum}\_{copynum}.sh
+ Self-made script (made from SCRIPT 4)
+ Each of these scripts runs stage 1 of analysis which is the running of GWAS and GIFT on subsampled data.
+ Currently doesnt run SCRIPT 12 ({slurm_job_id}_{subsample_num}_{phenotype}.R) because GWAS filtering needs to happen first
+ Runs:  SCRIPT 5 (subsample.py), SCRIPT 9 (PHYSICS_GWAS_OOP.py), SCRIPT 10 (pygwas), SCRIPT 11 (make_r_scripts)
+ Author: Sam Heysmond

> ## 9 physics_GWAS_OOP.py
+ The core script for GIFT analysis, this script will run a GIFT analysis on one set of input data per run.
+ This script has an inbuilt R script maker for the results so it also creates manhattan plot graphics.
+ The code has been updated to calculate and plot top 20 SNP data for each type of SNP E.G. PSNP4
+ Figure: 1A-Z 
+ Authors: Assistant Prof S Bray, Prof J Wattis, Prof C Rauch
+ Edited by: Sam Heysmond

> ## 10 pygwas
+ The core script for GWAS analysis, this script will run a GWAS analysis on one set of input data per run.
+ Reference: PENDING

> ## 11 make_r_scripts.py
+ Creates: SCRIPT 12 ({slurm_job_id}_{subsample_num}_{phenotype}.R)
+ Author: Sam Heysmond

> ## 12 {slurm_job_id}_{subsample_num}_{phenotype}.R
+ The individual R scripts that make manhattan plots for the GWAS data (since GIFT code includes the creation of an R script).
+ Calculates and plots the top 20 SNPs for the current dataset
+ This script gets run later after GWAS is filtered by SCRIPT 16 ({job_ID}_{subsample_number}_{phenotype}.sh)  
+ Figure: 1A-Z 
+ Author: Sam Heysmond

> ## 13 filtering_GWAS.sh
+ Copies all csv files in "output_files/" into a folder "output_files/csv_before_filter"
+ Runs: SCRIPT 14 (filter_snps_GWAS.py), SCRIPT 15 (GWAS_run_and_monitor.py)
+ Author: 

> ## 14 filter_snps_GWAS.py
+ Reads all csv files in "output_files/csv_before_filter" 
+ Reads core_files/output_3.table
+ Filters only the GWAS csv files from this list (uses multithreading to process multiple csvs at once)
  + Outputs the files back into "output_files/" overwriting the original csv (GWAS ONLY)
+ Makes SCRIPT X ({job_ID}_{subsample_number}_{phenotype}.sh) a bash script (which will run the R script for GWAS)
  + Puts it into R_parallel
+ Author: Sam Heysmond

> ## 15 GWAS_run_and_monitor.py
+ Runs all scripts inside R_parallel 
  + Puts them into completed_R_parallel when done
+ Author: Sam Heysmond

> ## 16 {job_ID}_{subsample_number}_{phenotype}.sh
+ The script to run R scripts for GWAS generated in stage 1
+ Figure: 1A-Z
+ Author: Sam Heysmond

# SCRIPTS : STAGE 2

> ## 17 SNP_tracker.sh
+ The main shell script which runs stage 2 of analysis which is gathering data on the SNPs and creating summary plots to compare both GIFT and GWAS results.
+ Runs: SCRIPT 18 (SNP_tracker.py)
+ Author: Sam Heysmond
  
> ## 18 SNP_tracker.py
+ The python script which concatonates the GWAS and GIFT data and processes it to satisfied summary plot/IDEAS 1-3 as described in the workflow timeline.
+ Outputs all the data to output_files/R_DATA
+ Author: Sam Heysmond

> ## 19 calc_thresholds.sh
+ Runs: SCRIPT 20 (calc_thresholds.py)
+ Author: 

> ## 20 calc_thresholds.py
+ looks inside of R_DATA (For gift) and R_DATA_FILTERED (for GWAS)
  + CHANGE THIS AFTER RE-RUNNING THE SCRIPTS SINCE IT WILL ALL BE IN R_DATA
+ Calculates thresholds for each subsample level, pval type and phenotype
+ Makes the file R_DATA/THRESHOLDS.csv
+ Author: Sam Heysmond

> ## 21 SNP_tracker_R_and_BASH_maker.sh
+ Runs: SCRIPT 22 (SNP_tracker_R_and_BASH_maker.py), SCRIPT 29 (Rscript_run_and_monitor.py)
+ Author: Sam Heysmond

> ## 22 SNP_tracker_R_and_BASH_maker.py
+ Creates:  SCRIPT 23({phenotype}_cumulative_t20_dataframe.sh) SCRIPT 24 ({phenotype}_cumulative_t20_dataframe.R), SCRIPT 25 ({phenotype}\_{positive/negative}\_control.sh), SCRIPT 26 ({phenotype}_{positive/negative}_control.R), SCRIPT 27 ({phenotype}\_{subsample_number}_{pval_type}_MANHATTAN.sh),SCRIPT 28 ({phenotype}\_{subsample_number}_{pval_type}_MANHATTAN.R)
+ WARNING: current version calculates threshold within the code for IDEA3, but IDEA1,2 use the file THRESHOLDS.csv
  + I will try to update this in future to make use of the THRESHOLDS.csv file to speed up runtime.
+ Author: Sam Heysmond

> ## 23 {phenotype}_cumulative_t20_dataframe.sh
+ Runs: SCRIPT 24 ({phenotype}_cumulative_t20_dataframe.R)
+ Author: 

> ## 24 {phenotype}_cumulative_t20_dataframe.R
+ Top 20 SNPs at 1000 subsample for each SNP type (in both GWAS and GIFT) tracked through 1000 to 200 subsample level.
+ Figure: 2A-Z , 3A-Z 
+ Author: Sam Heysmond
  
> ## 25 {phenotype}_{positive/negative}_control.sh
+ Runs: SCRIPT 26 ({phenotype}_{positive/negative}_control.R) 
+ Author: Sam Heysmond

> ## 26 {phenotype}_{positive/negative}_control.R
+ R script which tracks SNPs that fall into positive or negative control regions across phenotypes, subsample numbers and pval types
+ Figure: 2A-Z , 3A-Z 
+ Author: Sam Heysmond

> ## 27 {phenotype}_{subsample_number}_{pval_type}_MANHATTAN.sh
+ Runs: SCRIPT 28 ({phenotype}_{subsample_number}_{pval_type}_MANHATTAN.R)
+ Author: Sam Heysmond
  
> ## 28 {phenotype}_{subsample_number}_{pval_type}_MANHATTAN.R
+ The R script that makes average manhattan plots for each subsample level, pval type and phenotype.
+ Figure: 2A-Z 
+ Author: Sam Heysmond

> ## 29 Rscript_run_and_monitor.py
+ Runs: SCRIPT 15 ({phenotype}_cumulative_t20_dataframe.sh), SCRIPT 17 ({phenotype}\_{positive/negative}\_control.sh) , SCRIPT 19 ({phenotype}_{subsample_number}_{pval_type}_MANHATTAN.sh)
+ Author: 

# SCRIPTS : STAGE 3

> ## 30 Stage_3.sh
+ Runs: SCRIPT 31 (Zoom_In.py) , SCRIPT 31 (cross_referencing.py)
+ Author: Sam Heysmond
  
> ## 31 Zoom_In.py
+ Figure: 4A-Z 
+ Author: Sam Heysmond

> ## 31 cross_referencing.py
+ Figure: 5 / Table 1 
+ Author: Sam Heysmond


# PENDING POSITION

> ## X filter_snps.py (dont need?)
+ Looks into the R_DATA file and filters out the SNPs in GWAS that GIFT already removed
+ Outputs the filtered GWAS data into R_DATA_FILTERED
+ THIS NEEDS CHANGING TO LOOK AT ALL FILES IN "output_files" FOLDER INSTEAD! OR replace for other script combo(filtering_GWAS.sh/filter_snps_GWAS.py) <- ive replaced with these files
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