# Introduction
GWAS has been the go-to method for analysing SNP data for years but its underlying methods may lead to lack of biological sensitivity due to relying on averages to calculate significance. I have been tasked with comparing GIFT, a new method for analysing SNP data, to GWAS. In this analysis I aim to test how well GIFT performs compared to GWAS through downsampling of data, specifically to see how well GIFT and GWAS hold up when sample sizes are low. I also aim to investigate some of the positive detected peaks to validate the detection of GIFT. In summary, I begin with some VCF (genotype) and csv (masterlist) data on Arabidopsis Thaliana and hope to end up with some figures that will demonstrate the effectiveness of GIFT and GWAS at lower sample sizes.

---
# Timeline
The timeline is split into 3 main stages:
### 1) Subsampling the data and running GWAS vs GIFT (on Arabidopsis Thaliana data)
### 2) Concatonating the results from stage 1 and comparing GIFT vs GWAS on the basis of 3 main ideas (see stage 2 for more on these). This is done on  Arabidopsis Thaliana data.
### 3) Investigating significant results from GIFT and GWAS and performing GO analysis

### Indexing will be denoted as such: Stage.Step.StepSection e.g. 1.3.2
---
># 1.1 Set up workspace
##  1.1.1 Instructions
+ We begin by setting up our workplace on the HPC. First we organise the file structure (see File_structure.md) and then downloading and installing conda (following the instructions on the official website for Linux install https://conda.io/projects/conda/en/latest/user-guide/install/linux.html) so we can set up conda environments for the analysis. See "External_References_And_Info.md" for conda website referencing.

## 1.1.2 - Files (input)
+ Miniconda3-latest-Linux-x86_64.sh (downloaded from conda website)

## 1.1.3 - software versions needed
+ conda(miniconda) (latest)
+ conda_environment_setup.sh (script)

## 1.1.4 Expected outcome
+ A file structure and conda environments ready to begin analysis.

---

> #  1.2 Create scripts for stage 1
##  1.2.1 -Instructions
+ Enter the phenotype columns (from master_list.csv) that you wish to analyse into the file phenotypes_list.txt. Then run the script subsample_setup.sh using the command sbatch subsample_setup.sh from the batch_files folder.

## 1.2.2 - Files (input)
+ file 1) phenotypes_list.txt
+ file 2) 1001genomes_snp_biallelic_only_ACGTN.vcf
+ file 3) master_list.csv

## 1.2.3 - Files (output)
+ file 1) JOB_LIST.csv
+ file 2) subrun_{phenotype}_{subsampleNum}\_{copynum}.sh e.g. subrun_Mo98_400_3.sh

## 1.2.4 - software versions needed
+ subsample_script_setup.py (script)
+ subsample_env (conda environment)
+ python3_env (conda environment)

## 1.2.5 Expected outcome
+ A folder at ~/batch_files/parallel/ which contains all the generated scripts needed for stage 1 of analysis.

---
> #  1.3 Execute the scripts for stage 1 
##  1.2.1 -Instructions
+ Run the script run_and_monitor.sh using the command "sbatch run_and_monitor.sh" from the folder "batch_files/".

## 1.2.2 - Files (input)
+ file 1) subrun_{phenotype}_{subsampleNum}\_{copynum}.sh e.g. subrun_Mo98_400_3.sh. 
  + In our case we have 1000 individual subrun files (500 for Mo98 and 500 for Na23)
  + These are input and run automatically by the run_and_monitor.py script.

## 1.2.3 - Files (output)
+ file 1) JOB_LIST.csv update per subrun, giving information on the JOB_ID, the subsample number and phenotype.
+ file 2) subsamples_{subsample_num}_{JOB_ID}.txt
+ file 3) 4 Manhattan plots for each subrun (1 for GWAS and 3 for GIFT) e.g. leaf_ionome_Mo98_GWAS_MANHATTAN_1000_812222.png
+ file 4) 4 R scripts per subrun that make each of the manhattan plots 
+ file 5) 4 csv files per subrun that contain data on the top 20 most significant SNPs
+ file 6) 2 csv files per subrun that contain all the SNP data analysed (1 for GIFT, 1 for GWAS)

## 1.2.4 - software versions needed
+ run_and_monitor.py 
  + python3_env (conda environment)
+ subrun_{phenotype}_{subsampleNum}\_{copynum}.sh e.g. subrun_Mo98_400_3.sh (background script)
  + subsample_env (conda environment)
  + gift_env (conda environment)
  + gwas_env (conda environment)
  + python3_env (conda environment)
  + r_env (conda environment)

## 1.2.5 Expected outcome
+ For each phenotype, expect at least 3 days of runtime. In our case we used Leaf ionome of Mo98 and Na23 and this took around 6-7 days to process on the HPC. You should obtain all the SNP analysis results as denoted in step 1.2.3 and this should conclude stage 1 of analysis.
---

> #  2.1 Generate summary plots for GIFT vs GWAS 
##  1.2.1 -Instructions
+ Run the SNP_tracker.sh script from batch files using "sbatch SNP_tracker.sh". This script will concatonate the information gained during stage 1 and generate individual batch and R scripts for the script "Rscript_run_and_monitor.py" to manage and execute.

## 1.2.2 - Files (input)
+ directory 1) "output_files". All files within this are passed in for the script to handle
  + This is currently set to a specific directory in the code but in later version will be adaptable in automatically detecting the file structure.

## 1.2.3 - Files (output)
+ The total number of files is denoted below; we looked into two phenotypes so half these files are for one phenotype each
+ file 1) IDEA 1 R script: 2x
+ file 1) IDEA 1 batch file: 2x
+ file 1) IDEA 1 figures: 2x
+ file 1) IDEA 2 R script: 4x
+ file 1) IDEA 2 batch file: 4x
+ file 1) IDEA 2 figures: 4x
+ file 1) IDEA 3 R script: 40x
+ file 1) IDEA 3 batch file: 40x
+ file 1) IDEA 3 figures: 40x

## 1.2.4 - software versions needed
+ SNP_tracker.py (concatting data and generating subrun scripts)
  + gift_env
+ Rscript_run_and_monitor.py (managing all the subrun scripts)
  + gift_env
+ IDEA 1 script e.g. Mo98_cumulative_t20_dataframe.sh (background subrun script)
  + gift_env
+ IDEA 2 script e.g. Na23_negative_control.sh (background subrun script)
  + gift_env
+ IDEA 3 script e.g. Mo98_800_AVERAGE_P_MANHATTAN.sh (background subrun script)
  + gift_env
+ 

## 1.2.5 Expected outcome
+ The figures generated (in 1.2.3) will match to the ideas as denoted in 1.2.3. Idea1: Track the top 20 most significant SNPs from GWAS and GIFT across the FULL dataset, using a box plot for each downsample. Idea2: Track snps within the boundaries of positive and negative control regions per phenotype. Idea3: Create an average manhattan plot for each down sampling level for each SNP type e.g. Abs theta, PSNP4 etc. These figures should be placed in the correct folder "summary_plots" and relative IDEA folder e.g. Idea 1 figures get placed in "IDEA1" folder as described in "file_structure.md".
---























