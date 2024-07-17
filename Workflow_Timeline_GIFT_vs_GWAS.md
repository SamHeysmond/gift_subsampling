# Introduction
GWAS has been the go-to method for analysing SNP data for years but its underlying methods may lead to lack of biological sensitivity due to relying on averages to calculate significance. I have been tasked with comparing GIFT, a new method for analysing SNP data, to GWAS. In this analysis I aim to test how well GIFT performs compared to GWAS through downsampling of data, specifically to see how well GIFT and GWAS hold up when sample sizes are low. I also aim to investigate some of the positive detected peaks to validate the detection of GIFT. In summary, I begin with some VCF (genotype) and csv (masterlist) data on Arabidopsis Thaliana and hope to end up with some figures that will demonstrate the effectiveness of GIFT and GWAS at lower sample sizes.

---
# Timeline
The timeline is split into 3 main stages:
### 1) Subsampling the data and running GWAS vs GIFT (on Arabidopsis Thaliana data)
### 2) Concatonating the results from stage 1 and comparing GIFT vs GWAS on the basis of 3 main ideas (see stage 2 for more on these). 
### 3) Investigating significant results from GIFT and GWAS by zooming into the peaks and performing GO analysis to see which genes are picked up by each method.

### Indexing will be denoted as such: Stage.Step.Section e.g. 1.3.2
---
># 0.1.0 Set up workspace
##  0.1.1 Instructions
+ We begin by setting up our workplace on the HPC. First we organise the file structure see [file_structure.md](file_structure.md) and then downloading and installing conda (following the instructions on the official website for Linux install https://conda.io/projects/conda/en/latest/user-guide/install/linux.html) so we can set up conda environments for the analysis. See [Other_references.md](Other_references.md) for conda website referencing.

## 0.1.2 - Files (input)
+ Miniconda3-latest-Linux-x86_64.sh (downloaded from conda website)
+ SNP_Matrix core_files/
+ 1001genomes_snp_biallelic_only_ACGTN.vcf
+ master_list.csv
+ physics_GWAS_OOP.py
+ physics_gwas_1.sh
+ leaf_ionome_Sr88_GWAS.sh

## 0.1.3 - Files (output)
+ output_1.vcf (temp file of main vcf with AN and AC labels)
+ output_2.vcf (temp file of a VCF with all SNPs that have AC<30>)
+ output_3.table (SNP location data from output_2.vcf)

## 0.1.4 - scripts needed
+ get_files.sh 
+ conda_environment_setup.sh 
+ filtering_vcf.sh

## 0.1.5 Expected outcome
+ We obtained a file structure containing the core files needed and conda environments ready to begin analysis. We also get a file "output_3.table" which has all the locations of SNPs with AC<30 . The latter two input files were discarded once they were read as they were examples of how to run some of the software and so were not needed after understanding them.

---

> #  1.1.0 Create scripts for stage 1
##  1.1.1 -Instructions
+ We entered the phenotyped columns (from master_list.csv) that we wished to analyse into the file phenotypes_list.txt and placed it into the core_files directory. Then ran the script subsample_setup.sh using the command sbatch subsample_setup.sh from the batch_files folder.

## 1.1.2 - Files (input)
+ file 1) phenotypes_list.txt
+ file 2) 1001genomes_snp_biallelic_only_ACGTN.vcf
+ file 3) master_list.csv

## 1.1.3 - Files (output)
+ file 1) JOB_LIST.csv
+ file 2) subrun_{phenotype}_{subsampleNum}\_{copynum}.sh e.g. subrun_Mo98_400_3.sh

## 1.1.4 - software versions and scripts
+ subsample_script_setup.py (script)

## 1.1.5 Expected outcome
+ A folder at ~/batch_files/parallel/ which contains all the generated scripts needed for stage 1 of analysis ready to be run. A JOB_LIST csv file prepared to store jobs when they get run shortly after this step. 

---
> #  1.3.0 Execute the scripts for stage 1 
##  1.3.1 -Instructions
+ We ran the script run_and_monitor.sh using the command "sbatch run_and_monitor.sh" from the folder "batch_files/".

## 1.3.2 - Files (input)
+ Folder 1)batch_files/parallel/ containing all the subrun scripts BEFORE being run e.g. subrun_{phenotype}_{subsampleNum}\_{copynum}.sh (subrun_Mo98_400_3.sh.)
  + In our case we have 2000 individual subrun files (1000 for Mo98 and 1000 for Na23)
  + These are input and run automatically by the run_and_monitor.py script.

## 1.3.3 - Files (output)
+ Folder 1) batch_files/completed_parallel/ containing all the subrun scripts AFTER being run e.g. subrun_{phenotype}_{subsampleNum}\_{copynum}.sh (subrun_Mo98_400_3.sh.)
+ file 1) core_files/all_vcf_samples.txt
+ file 2) core_files/subsample_test_files/subsamples_{subsample_num}_{ID}
+ file 2) core_files/subsampled_phenotype_${subsample_num}_${SLURM_JOB_ID}.csv
+ file 3) core_files/subsampled_genotype_${subsample_num}_${SLURM_JOB_ID}.vcf
+ file 4) JOB_LIST.csv updated per subrun, giving information on the JOB_ID, the subsample number and phenotype.
+ file 5) subsamples_{subsample_num}_{JOB_ID}.txt
+ file 6) 2 R scripts per subrun that will make each of the manhattan plots for GIFT and GWAS
+ file 7) 4 csv files per subrun that contain data on the top 20 most significant SNPs
+ file 8) 2 csv files per subrun that contain all the SNP data analysed (1 for GIFT, 1 for GWAS)
+ file 9) Manhattan plots for GWAS: named as such ({phenotype}_GIFT_MANHATTAN_{subsample_num}_{ID}.png)

## 1.3.4 - scripts needed
+ run_and_monitor.sh
+ run_and_monitor.py 
+ subrun_{phenotype}_{subsampleNum}\_{copynum}.sh e.g. subrun_Mo98_400_3.sh (background script)
+ subsample.py
+ physics_GWAS_OOP.py
+ pygwas
+ make_r_scripts.py


## 1.3.5 Expected outcome
+ A list of all sample IDs present inside the vcf file (required for the subsample_script_setup.py).
+ Temporary files of all subsampled IDs for a particular JOB_ID for a subrun.
+ Temporary files of subsampled phenotype (in csv format) and genotype (in vcf format) used for individual subruns.
+ CSV information for GIFT and GWAS of the results of each subrun along with some information on the TOP20 SNPs for GIFT 
  + GWAS has its top 20 SNPs calculated later in R scripts when they are run.
+ Manhattan plots for GIFT only 
  + GWAS manhattan plots are created later in R scripts after filtering
+ For each phenotype, expect at least 3 days of runtime on the current version. In our case we used Leaf ionome of Mo98 and Na23 and this took around 6-7 days to process on the HPC. 

---

> # 1.4.0 Post GWAS filtering
## -Instructions
+ Because the GWAS input file was a SNP matrix and not a VCF we still needed to filter it for SNPs that fall under AC of 30 on all the csvs that were output from pygwas. To do this we ran the below .sh scripts to leave us with filtered csvs of our GWAS data.

## 1.4.1 Files (input)
+ Folder 1) output_files/ All of the csv files from pygwas that exist in this folder
+ File 1) core_files/output_3.table


## 1.4.2 Files (output)
+ Folder 1) batch_files/R_parallel/ all R scripts for GWAS befor being run
+ Folder 2) batch_files/completed_R_parallel/ all R scripts for GWAS after being run
+ File 1) output_files/ All filtered GWAS csv files
+ File 2) GWAS manhattan plots {phenotype}_GWAS_MANHATTAN_{subsample_num}_{ID}.png

## 1.4.3 scripts needed
+ filtering_GWAS.sh
+ filter_snps_GWAS.py
+ filtering_GWAS.sh
+ GWAS_run_and_monitor.py
+ {job_ID}_{subsample_number}_{phenotype}.sh (background script)

## 1.4.4 Expected outcome
+ We expected all the GWAS csv files to be filtered and replaced in the folder "output_files" for further analysis. And we obtained all images of our GWAS manhattan plots. This concluded the main steps for STAGE 1.
---

> #  2.1.0 Concatonate GWAS and GIFT data
##  2.1.1 -Instructions
+ We ran the SNP_tracker.sh script from batch files using "sbatch SNP_tracker.sh".

## 2.1.2 - Files (input)
+ Folder 1) "output_files/". All GWAS and GIFT csv files are passed in 
  + This is currently set to a specific directory in the code but in later version will be adaptable in automatically detecting the file structure.

## 2.1.3 - Files (output)
+ IDEA 1: cumulative t20 dataframes across all subsamples e.g. Mo98_cumulative_t20_dataframe.csv
+ IDEA 2: positive and negative control data for each phenotype e.g. Mo98_positive_control.csv
+ IDEA 3: Combined information for all runs of a phenotype, subsample and method combo e.g. Mo98_GWAS_200_ALL.csv

## 2.1.4 - scripts needed
+ SNP_tracker.sh
+ SNP_tracker.py (concatting data)


## 2.1.5 Expected outcome
+ Combines information to fit the format of three ideas
  +  IDEA1: Track the top 20 most significant SNPs from GWAS and GIFT across the FULL dataset, for each downsample. 
  +  IDEA2: Track snps within the boundaries of positive and negative control regions per phenotype. 
  +  IDEA3: Combine data from all runs for a given phenotype, subsampel and method combo
  
---

> # 2.2.0 Calculate threshold values for figures
## -Instructions
+ We used a script to look into the data created previously and calculate two thresholds (bonferroni and Benjamini-Hochberg-Yekutieli)

## 2.2.1 Files (input)
+ Folder 1) R_DATA/ all the averages files for each combo of subrun format: {phenotype}\_{GIFT/GWAS}_{subsample_num}_ALL.csv

## 2.2.2 Files (output)
+ File 1) THRESHOLDS.csv

## 2.2.3 scripts needed
+ calc_thresholds.sh
+ calc_thresholds.py

## 2.2.4 Expected outcome
+ We obtained a file containing all threshold information for each average of the run combinations, for instance for each phenotype, subsample number and pvalue type there exists a BF threshold anda BHY threshold value.

> # 2.3.0 Create figures for STAGE 2
## -Instructions
+ Next we needed to make all the R scritps for the concatonated data and run each of them to get our figures. To do this we ran "sbatch SNP_tracker_R_and_BASH_maker.sh" from batch_files directory.

## 2.3.1 Files (input)
+ R_DATA/ All files within here that were made in the last two steps of STAGE2
  + Each of the scripts (numbred 23-28) will call for their own required files for the IDEAS they contribute to

## 2.3.2 Files (output)
+ Folder 1) output_files/summary_plots/IDEA1/ All IDEA 1 figures
+ Folder 2) output_files/summary_plots/IDEA2/ All IDEA 2 figures
+ Folder 3) output_files/summary_plots/IDEA3/ All IDEA 3 figures

## 2.3.3 scripts needed
+ SNP_tracker_R_and_BASH_maker.sh
+ SNP_tracker_R_and_BASH_maker.py
+ Rscript_run_and_monitor.py
+ background scripts
  + SCRIPTS 23-28 which are made by SNP_tracker_R_and_BASH_maker.py and run by Rscript_run_and_monitor.py
  
## 2.3.4 Expected outcome
+ We created three folders with figures for each of the IDEAS aforementioned, ending our STAGE 2 of analysis.

---

> # 3.1.0 Investigating the peaks at 1000 subsamples and crossreferencing SNPs with genes

## 3.1.1 Instructions
+ For our final stage of analysis we began by running "sbatch Stage_3.sh" from the batch_files folder to run the scripts in the below scripts section.

## 3.1.2 Files (input)
+ TAIR10_GFF_genes.gff
+ Data of average values from 1000 subsample levels obtained previously in STAGE 2, in the folder R_DATA/
  + Format: {phenotype}_{method_type}_1000_ALL.csv
+ THRESHOLDS.csv
  + Previously calculated thresholds file

## 3.1.3 Files (output)
+ Folder 1) stage_3_scripts/ folder with the R scripts to make the ZOOM plots
  + format : {phenotype}\_{GIFT/GWAS}\_AVERAGE_{pval_type}_ZOOM.R
+ Folder 2) summary_plots/stage_3/ folder fille with ZOOM plots
  + format : {phenotype}_{GIFT/GWAS}_AVERAGE_{pval_type}_ZOOM.png
+ Folder 3) GO_DATA
  + Bed files for top 1000 significant SNPs 
    + format : {phenotype}_{GIFT/GWAS}_{subsample_num}_ALL.bed
  + Intersect results between the GFF and bed files
    + format : Intersect_results{phenotype}_{GIFT/GWAS}_{subsample_num}_ALL.txt
  + FINAL intersect results that have been filtered to only give unique gene IDs
    + format : FINAL_Intersect_results_{phenotype}_{GIFT/GWAS}_{subsample_num}_ALL.txt

## 3.1.4 scripts needed
+ Stage_3.sh
+ Zoom_In.py
+ cross_referencing.py

## 3.1.5  Expected outcome
+ We gathered zoomed in manhattan plots containing information on the peak and key positive control genes in the area and the SNP positions around the peak relative to thresholds. We also obtained lists of gene IDs that the top 1000 SNPs (from select combinations of subsamples, methods and phentoypes) had crossed over with.

> # 3.2.0 GO analysis
## 3.2.1 -Instructions
+ To investigate which functions the significant genes were involved in we finally performed a GO analysis using DAVID.

## 3.2.2 Files (input)
+ File 1) Intersect results from the previous step in the folder GO_DATA/
  + format : FINAL_Intersect_results_{phenotype}_{GIFT/GWAS}_{subsample_num}_ALL.txt

## 3.2.3 Files (output)
+ N/A we wrote select outputs from the website manually into a table

## 3.2.4 websites needed
+ DAVID
+ PANTHERS

## 3.2.5 Expected outcome
+ We found descriptions of biological processes, cellular components and molecular function for genes that were significantly related to one another in the GO analysis.
























