# gift subsampling and testing
This README.md file gives a summary of what each of the files on the github contains 
For referencing information, please see "References_and_info.md"
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
# SCRIPTS

> ## 1 conda_environment_setup.sh
+ Installation and setup of all environments used
    + Dependendcies
    + Version numbers


> ## 2 subsample_setup.sh
+ Sets up folders and files needed for analysis
+ Removes previous files in those locations (refreshing the data)
  + (This can be disabled by commenting out the rm commands and echo command for JOB_LIST.csv)
+ Runs: subsample_script_setup.py


> ## 3 subsample_script_setup.py
+ Creates all the batch files which will subsample the vcf data and perform GIFT and GWAS

> ## 4 run_and_monitor.sh
+ Runs: run_and_monitor.py
  + this is required to run to manage the batch scripts for the entirety of the subsampling and testing of GIFT and GWAS. So, give this file max run time.

> ## 5 run_and_monitor.py
+ Manages all of the batch scripts in "batch_files/parallel", running as many in parallel as possible. When the batch file has been put in the que, it moves it to a different folder "batch_files/completed_parallel"
+ Runs: subrun_{phenotype}_{subsampleNum}_{copynum}.sh

> ## 6 subrun_{phenotype}_{subsampleNum}_{copynum}.sh
+ Self-made script
+ I

> ## 7 physics_GWAS_OOP.py
+ The core script for GIFT analysis, this script will run a GIFT analysis on one set of input data per run.
+ This script has an inbuilt R script maker for the results so it also creates manhattan plot graphics
+ The code has been updated to calculate and plot top 20 SNP data for each type of SNP E.G. PSNP4

> ## 8 pygwas
+ The core script for GWAS analysis, this script will run a GWAS analysis on one set of input data per run.

> ## 9 make_r_scripts.py
+ Creates manhattan plots from the GWAS data
+ Calculates and plots the top 20 SNPs for the current dataset
+ These R scripts are run during the subrun code (SCRIPT 6)

> ## {slurm_job_id}_{i}_{phenotype}.R
+ I

> ## SNP_tracker.sh
+ I
  
> ## SNP_tracker.py
+ I



# DATA FILES (INPUT)


# DATA FILES (OUTPUT)



> ## Windows_Method.sh 
+ Installation steps
    + Dependendcies
    + Version numbers
+ Main code for step 2 in Workflow_Timeline.md
---
> ## References_and_info.md
  + References
      + Scientific papers describing the tool/algorithm
      + Github links to tools
      + Scripts and software that were used
---
> ## color_h.py
+ For Authorship and information see "References_and_info.md"
+ Code used to create a hyrophobicity scale on the aligned proteins

> ## pymol_code_Yasmin.txt
+ Code for creating figure images in pymol
+ Instructions for obtaining protein charge information from APBS website
+ For Authorship and information see "References_and_info.md"
---