# GIFT vs GWAS subsampling and testing
This README.md file gives a summary of what each of the files on the github contains. Any references for local files are also detailed in their description. 
For referencing and information on external software, websites etc, please see [Other_references.md](Other_references.md)
---
## Key:
> {variable_name} = A varible that is assigned many values e.g. {subsample_num} may be 200,400 or 1000

> Figure contribution: 1A-Z = The script helps create the specified figure in final writeup (A-Z is placeholder as subfigure labels can vary within this range)
---
# Files and a breif description of what they contain / do
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
    + To data that is needed 
      + references of tools are in Referenecs_and_info.md
      + references of scripts/files are in README.md
---
> ## Other_references.md
  + References
      + Scientific papers describing the tool/algorithm/website (where possible)
        + links to the website if it is a site
      + Github links to tools
---
# SCRIPTS : PREP STAGE

> ## 0 get_files.sh
+ Fetching local files necessary for the analysis
+ Author: Sam Heysmond

> ## 1 conda_environment_setup.sh
+ Installation and setup of all environments used
    + Dependendcies
    + Version numbers
+ Author: Sam Heysmond

> ## 2 filtering_vcf.sh
+ Creating a table of all SNPs below AC of 30 from the main vcf
+ Author: Sam Heysmond

---
# SCRIPTS : STAGE 1

> ## 3 subsample_setup.sh
+ Sets up folders and files needed for subsample analysis
  + But first removes previous files in those locations to ensure no data remains to start with
    + (This can be disabled by commenting out the rm commands and echo command for JOB_LIST.csv)
+ Runs: subsample_script_setup.py (SCRIPT 4)
+ Author: Sam Heysmond

> ## 4 subsample_script_setup.py
+ Creates all the batch files which will subsample the vcf data and perform GIFT and GWAS analysis
  + Creates: SCRIPT 8 (subrun_{phenotype}_{subsampleNum}\_{copynum}.sh)
+ Author: Sam Heysmond

> ## 5 subsample.py
+ The script used to subsample the vcf and phenotype file ID's which match in both genotype and phenotype files (for a given subsample number e.g. 400).
+ Author: Sam Heysmond

> ## 6 run_and_monitor.sh
+ Runs: SCRIPT 7 (run_and_monitor.py)
  + this is required to run to manage the batch scripts for the entirety of the subsampling and testing of GIFT and GWAS; give this file max run time the server allows (in this case it was 7 days).
+ Author: Sam Heysmond

> ## 7 run_and_monitor.py
+ Manages all of the batch scripts in "batch_files/parallel", running as many in parallel as possible. When the batch file has been put in the que, it moves it to a different folder "batch_files/completed_parallel"
+ Runs: SCRIPT 8 (subrun_{phenotype}_{subsampleNum}\_{copynum}.sh)
+ Author: Sam Heysmond

> ## 8 subrun_{phenotype}_{subsampleNum}\_{copynum}.sh
+ Each of these scripts runs stage 1 of analysis which is the running of GWAS and GIFT on subsampled data.
+ Currently doesnt run SCRIPT 12 ({slurm_job_id}_{subsample_num}_{phenotype}.R) to generate plots because GWAS filtering needs to happen first.
+ Created by: SCRIPT 4
+ Runs:  SCRIPT 5 (subsample.py), SCRIPT 9 (PHYSICS_GWAS_OOP.py), SCRIPT 10 (pygwas), SCRIPT 11 (make_r_scripts)
+ Author: Sam Heysmond

> ## 9 physics_GWAS_OOP.py
+ The core script for GIFT analysis, this script will run a GIFT analysis on one set of input data per run.
+ This script has an inbuilt R script maker for the results so it also creates manhattan plot graphics.
+ The code has been updated to calculate and plot top 20 SNP data for each type of SNP E.G. PSNP4
+ Figure contribution: 1A-Z 
+ Authors: Assistant Prof S Bray, Prof J Wattis, Prof C Rauch
+ Edited by: Sam Heysmond

> ## 10 pygwas
+ The core script for GWAS analysis, this script will run a GWAS analysis on one set of input data per run.
+ Author: Ümit Seren
+ Github: [PyGWAS github](https://github.com/timeu/PyGWAS/tree/master)

> ## 11 make_r_scripts.py
+ Creates: SCRIPT 12 ({slurm_job_id}_{subsample_num}_{phenotype}.R)
+ Author: Sam Heysmond

> ## 12 {slurm_job_id}_{subsample_num}_{phenotype}.R
+ The individual R scripts that make manhattan plots for the GWAS data (since GIFT code includes the creation of an R script).
+ Calculates and plots the top 20 SNPs for the current dataset (GIFT code also has this embedded within)
+ These scripts are run after GWAS is filtered by SCRIPT 16 ({job_ID}_{subsample_number}_{phenotype}.sh)  
+ Figure contribution: 1A-Z 
+ Author: Sam Heysmond

> ## 13 filtering_GWAS.sh
+ Copies all csv files in "output_files/" into a folder "output_files/csv_before_filter"
+ Runs: SCRIPT 14 (filter_snps_GWAS.py), SCRIPT 15 (GWAS_run_and_monitor.py)
+ Author: Sam Heysmond

> ## 14 filter_snps_GWAS.py
+ Reads all csv files in "output_files/csv_before_filter" 
+ Reads core_files/output_3.table
+ Filters only the GWAS csv files from this list (uses multithreading to process multiple csvs at once)
  + Outputs the files back into "output_files/" overwriting the original csv (GWAS ONLY)
+ Makes SCRIPT 16 ({job_ID}_{subsample_number}_{phenotype}.sh) a bash script (which will run the R script for GWAS manhattan plots)
  + Puts it into R_parallel ready for running by SCRIPT 16
+ Author: Sam Heysmond

> ## 15 GWAS_run_and_monitor.py
+ Runs all scripts inside R_parallel 
  + Puts them into completed_R_parallel when done
+ Author: Assistant Prof S Bray
+ Edited by: Sam Heysmond

> ## 16 {job_ID}_{subsample_number}_{phenotype}.sh
+ The script to run R scripts for GWAS generated in stage 1
+ Figure contribution: 1A-Z
+ Author: Sam Heysmond

# SCRIPTS : STAGE 2

> ## 17 SNP_tracker.sh
+ The main shell script for stage 2 of analysis which is gathering data on the SNPs and creating summary plots to compare both GIFT and GWAS results.
+ Runs: SCRIPT 18 (SNP_tracker.py)
+ Author: Sam Heysmond
  
> ## 18 SNP_tracker.py
+ The python script which concatonates the GWAS and GIFT data, and processes it to satisfy summary plot/IDEAS 1-3 as described in the [workflow timeline](Workflow_Timeline_GIFT_vs_GWAS.md).
+ Outputs all the data to output_files/R_DATA
+ Author: Sam Heysmond

> ## 19 calc_thresholds.sh
+ Runs: SCRIPT 20 (calc_thresholds.py)
+ Author: Sam Heysmond

> ## 20 calc_thresholds.py
+ Looks inside of R_DATA, gathering the contatenated files with averages of all GIFT and GWAS runs.
+ Calculates thresholds for each subsample level, pval type and phenotype
+ Makes the file R_DATA/THRESHOLDS.csv
+ Author: Sam Heysmond

> ## 21 SNP_tracker_R_and_BASH_maker.sh
+ Runs: SCRIPT 22 (SNP_tracker_R_and_BASH_maker.py), SCRIPT 29 (Rscript_run_and_monitor.py)
+ Author: Sam Heysmond

> ## 22 SNP_tracker_R_and_BASH_maker.py
+ Creates:  SCRIPT 23({phenotype}\_cumulative\_t20\_dataframe.sh) SCRIPT 24 ({phenotype}\_cumulative\_t20\_dataframe.R), SCRIPT 25 ({phenotype}\_{positive/negative}\_control.sh), SCRIPT 26 ({phenotype}_{positive/negative}_control.R), SCRIPT 27 ({phenotype}\_{subsample_number}\_{pval_type}\_MANHATTAN.sh), SCRIPT 28 ({phenotype}\_{subsample_number}\_{pval_type}_MANHATTAN.R)
+ N.B: current version calculates threshold within the code for IDEA3, but IDEA1,2 use the file THRESHOLDS.csv
  + I will try to update this in future to make use of the THRESHOLDS.csv file to speed up runtime and stop re-calculation of them in other scripts.
+ Author: Sam Heysmond

> ## 23 {phenotype}_cumulative_t20_dataframe.sh
+ Runs: SCRIPT 24 ({phenotype}_cumulative_t20_dataframe.R)
+ Author: Sam Heysmond

> ## 24 {phenotype}_cumulative_t20_dataframe.R
+ Top 20 SNPs at 1000 subsample for each SNP type (in both GWAS and GIFT) tracked through 1000 to 200 subsample level.
+ Performs KW test on differences of the boxplots
+ Figure contribution: 2A-Z , 3A-Z 
+ Author: Sam Heysmond
  
> ## 25 {phenotype}_{positive/negative}_control.sh
+ Runs: SCRIPT 26 ({phenotype}_{positive/negative}_control.R) 
+ Author: Sam Heysmond

> ## 26 {phenotype}_{positive/negative}_control.R
+ R script which tracks SNPs that fall into positive or negative control regions (genes) across phenotypes, subsample numbers and pval types
+ Also runs KW test on the boxplots generated
+ Figure contribution: 2A-Z , 3A-Z 
+ Author: Sam Heysmond

> ## 27 {phenotype}_{subsample_number}_{pval_type}_MANHATTAN.sh
+ Runs: SCRIPT 28 ({phenotype}\_{subsample_number}_{pval_type}_MANHATTAN.R)
+ Author: Sam Heysmond
  
> ## 28 {phenotype}_{subsample_number}_{pval_type}_MANHATTAN.R
+ The R script that makes average manhattan plots for each subsample level, pval type and phenotype.
+ Figure contribution: 2A-Z 
+ Author: Sam Heysmond

> ## 29 Rscript_run_and_monitor.py
+ Runs: SCRIPT 15 ({phenotype}\_cumulative_t20_dataframe.sh), SCRIPT 17 ({phenotype}\_{positive/negative}\_control.sh) , SCRIPT 19 ({phenotype}\_{subsample_number}_{pval_type}_MANHATTAN.sh)
+ Author: Assistant Prof S Bray
+ Edited by: Sam Heysmond

# SCRIPTS : STAGE 3

> ## 30 Stage_3.sh
+ Runs: SCRIPT 31 (Zoom_In.py) , SCRIPT 31 (cross_referencing.py)
+ Author: Sam Heysmond
  
> ## 31 Zoom_In.py
+ Creates R scripts that make a Zoomed in manhattan plot around the peaks for Na and Mo phenotypes
+ Figure contribution: 4A-Z 
+ Author: Sam Heysmond

> ## 32 cross_referencing.py
+ Cross references significant SNPs with Arabidopsis thaliana genes by making a shell script to do so.
+ Figure contribution: Table 1 
+ Author: Sam Heysmond

> ## 33 cross_referencing_2.py
+ Cross references genes that were picked up by the previous script (SCRIPT 32) against dataframes of genes of interest to track when and where they were detected as significant.
+ Figure contribution: Table 1 
+ Author: Sam Heysmond
  
---
# DATA FILES (INPUT)
### These are files that we start with at the beginning of the analysis
> ## 1 SNP_Matrix
+ Used as input for the GWAS software "pygwas"
+ Reference: [1001 Genomes Project database](https://1001genomes.org/data/GMI-MPI/releases/v3.1/)

> ## 2 1001genomes_snp_biallelic_only_ACGTN.vcf
+ Genotype information from the 1001 genomes project on biallelic SNPs data from Arabidopsis Thaliana
  + This file was modified by Assistant Prof S Bray and created from "1001genomes_snp-short-indel_only_ACGTN.vcf" in the database.
+ Reference: [1001 Genomes Project database](https://1001genomes.org/data/GMI-MPI/releases/v3.1/)
+ Edited by: Assistant prof S Bray

> ## 3 master_list.csv
+ Leaf and seed ionome phenotype data from 1001 Genomes database on Arabidopsis Thaliana
+ Reference: [1001 Genomes Project database](https://1001genomes.org/data/GMI-MPI/releases/v3.1/)

> ## 4 phenotypes_list.txt
+ List of which phenotype columns to analyse. The scripts will use this to know what phenotypes to look into.
+ Reference: N/A

> ## 5 TAIR10_GFF_genes.gff
+ GFF3 file filtered to only contain gene information.
+ Reference: Assistant prof S Bray

> ## 6 Mo_genes_data.csv
+ List of gene ID and names of interest for genes related to Mo phenotype
+ Reference: N/A

> ## 7 Na_genes_data.csv
+ List of gene ID and names of interest for genes related to Na phenotype
+ Reference: N/A

---
# DATA FILES (OUTPUT)
### These are data files that are output during the analysis. These files may also be input for other scripts but are classed initially as output as we did not begin with these files.
> ## 1 JOB_LIST.csv
+ Output list of all subrun jobs for stage 1 of analysis. This is partly designed for error checking to ensure all jobs were properly completed, but in future versions of scripts it may be used to help with analysis and segregation of data.
+ Reference: N/A

> ## 2 all_vcf_samples.txt
+ List of all IDs of samples from the vcf folder. Some scripts may use this data to check that the vcf and master list (phenotypes) contain the same sample IDs before analysis can begin.
Reference: N/A

> ## 3 subsamples_{subsample_num}_{JOB_ID}.txt
+ List of IDs of all the subsamples for a particular job run. This is useful for error checking but is also used in the program to compare between the VCF and phenotype data for the sample IDs that have been selected.
  
> ## 4 output_1.vcf
+ Modified version of the 1001 genomes VCF that includes the tags AC and AN
+ Reference: N/A

> ## 5 output_2.vcf
+ Filtered version of output_1.vcf which only includes the sites that have less than the filtered AC (in this case AC<30)
+ Reference: N/A

> ## 6 output_3.table
+ List of all the SNP positions that had an AC<30 in the VCF (which is all the snps in output_2.vcf)
+ This list is used to filter out SNPs from GWAS output files.
+ Reference: N/A

> ## 7 {phenotype}_GWAS_MANHATTAN_{subsample_number}_{ID}.png
+ GWAS manhattan plot for all subsample numbers and both phenotypes
+ Figure: 1A-Z
+ Reference: N/A

> ## 8 {phenotype}_whole_genome_metrics_{subsample_num}_{ID}_{pval_type}.png
+ GIFT manhattan plot for all subsample numbers and both phenotypes
+ Figure: 1A-Z
+ Reference: N/A

> ## 9 {phenotype}_GWAS_MANHATTAN_{subsample_number}_{ID}.csv
+ GWAS manhattan plot data
+ Reference: N/A

> ## 10 {phenotype}_whole_genome_metrics_{subsample_num}_{ID}_{pval_type}.csv
+ GIFT manhattan plot data
+ Reference: N/A

> ## 11 {phenotype}_GWAS_T20_SNPS_{subsample_num}_{ID}.csv
+ The position and value of the T20 snps in GWAS
+ Reference: N/A

> ## 12 {ID}_T20_{pval_type}.csv
+ The position and value of the T20 snps for a specific pval type and job ID in GIFT
+ Reference: N/A

> ## 13 {phenotype}\_{GIFT/GWAS}_{subsample_num}_ALL.csv
+ Association results from the GIFT/GWAS output which has been averaged across all tests for the subsample number in the name.
+ Reference: N/A

> ## 14 {phenotype}_{positive/negative}_control.csv
+ All SNPs for all phenotypes, subsample numbers and pval types that lie within the positive and negative control regions
+ Reference: N/A

> ## 15 {phenotype}_cumulative_t20_dataframe.csv
+ The top 20 SNPs from each pval type at subsample 1000 tracked through to subsample 200
+ Reference: N/A

> ## 16 {phenotype}\_AVERAGE_{pval_type}_T20_LOCATIONS.csv
+ chromosome and position only of the top 20 snps
+ Reference: N/A

> ## 17 THRESHOLDS.csv
+ Bonferroni (BF) and Benjamini-Hochberg-Yekutieli thresholds for each pval type, phenotype and subsample level
  + Based on the average data for each subsample level
+ Reference: N/A

> ## 18 {phenotype}\_cumulative_t20_dataframe_{pval_type}.png
+ Box plot of the top 20 SNPs for specified phenotype and pval type across all subsample levels
+ Figure: 2A-Z
+ Reference: N/A

> ## 19 {phenotype}_cumulative_t20_dataframe_{pval_type}_KW_TEST.png
+ Kruskal wallace test between the boxes on the boxplot for OUTPUT 18
+ Figure: 3A-Z
+ Reference: N/A

> ## 20 {phenotype}\_{positive/negative}\_control_{pval_type}.png
+ Box plots of all the values of SNPs for a specified pval type and phenotype across all subsample levels
+ Figure: 2A-Z
+ Reference: N/A

> ## 21 {phenotype}\_{positive/negative}\_control_{pval_type}_KW_TEST.png
+ Kuskal Wallace test on OUTPUT 20 boxplots
+ Figure: 3A-Z
+ Reference: N/A

> ## 22 {phenotype}\_{subsample_num}\_AVERAGE_{pval_type}_MANHATTAN.png
+ Manhattan plots based on the averages data for each phenotype, subsample number and pval type combination
+ Figure: 2A-Z
+ Reference: N/A

> ## 23 {phenotype}\_{GIFT/GWAS}\_AVERAGE_{pval_type}_ZOOM.png
+ A zoomed in manhattan plot around the peak and positive control genes
+ Figure: 4A-Z
+ Reference: N/A

> ## 24 {phenotype}\_{GIFT/GWAS}_{subsample_num}_ALL.bed
+ The chromosome and position of all SNPs for a given combination of phenotype, method and subsample number that are above the threshold
  + For GIFT representation PSNP5 was used
  + For GWAS representation P was used 
+ Reference: N/A

> ## 25 Intersect_results{phenotype}_{GIFT/GWAS}_{subsample_num}_ALL.txt
+ Information from the gff3 file where OUTPUT 24 crosses over with data in the GFF file
+ Reference: N/A

> ## 26 FINAL_Intersect_results_{phenotype}_{GIFT/GWAS}_{subsample_num}_ALL.txt
+ List of gene IDs filtered out from OUTPUT 25
+ Reference: N/A

> ## 27 {phenotype}_Gene_Tracker.csv
+ Detection summary of genes of interest for each phenotype
+ Reference: N/A

> ## 28 slurm_error.err
+ All the errors output by the Ada server
+ Reference: N/A


> ## 29 slurm_out.out
+ All the output given by the Ada server
  + This is usually things we intend to be output such as "print('Script finished',flush=True)" in python scripts.
+ Reference: N/A

---
