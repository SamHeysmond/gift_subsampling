file structure for the analysis:
---
```
User Home
|----batch_files
|   |----get_files.sh
|   |----conda_environment_setup.sh
|   |----make_r_scripts.py
|   |----Rscript_run_and_monitor.py
|   |----run_and_monitor.py
|   |----run_and_monitor.sh
|   |----SNP_tracker.py
|   |----SNP_tracker.sh
|   |----subsample.py
|   |----subsample_script_setup.py
|   |----subsample_setup.sh
|   |----filtering_vcf.sh
|   |----GWAS_run_and_monitor.py
|   |----filter_snps_GWAS.py
|   |----filtering_GWAS.sh
|   |----calc_thresholds.sh
|   |----calt_thresholds.py
|   |----SNP_tracker_R_and_BASH_maker.py
|   |----SNP_tracker_R_and_BASH_maker.sh
|   |----Zoom_In.py
|   |----Stage_3.sh
|   |----cross_referencing.py
|   |----cross_reference_script.sh
|   |----parallel/
|       |----Stage_1_subrun_scripts e.g. subrun_leaf_ionome_Na23_1000_100.sh (before complete)
|   |----completed_parallel
|       |---- Stage_1_subrun_scripts e.g. subrun_leaf_ionome_Na23_1000_100.sh (after complete)
|   |----parallel_stage2/
|       |----IDEA 1 Script (before run) e.g. Mo98_cumulative_t20_dataframe.sh
|       |----IDEA_2_Script (before run) e.g. Na23_negative_control.sh
|       |----IDEA_3_Script (before run) e.g. Mo98_800_AVERAGE_P_MANHATTAN.sh
|   |----completed_parallel_stage2/
|       |----IDEA 1 Script (after run) e.g. Mo98_cumulative_t20_dataframe.sh
|       |----IDEA 2 Script (after run) e.g. Na23_negative_control.sh
|       |----IDEA 3 Script (after run) e.g. Mo98_800_AVERAGE_P_MANHATTAN.sh
|   |----R_parallel/
|       |----{number}_{id}_{subsample_number}_{phenotype}.R (GWAS R script)
|   |----completed_R_parallel/
|       |----{number}_{id}_{subsample_number}_{phenotype}.R (GWAS R script)
|
|
|----core_files/
|   |----SNP_Matrix (folder)
|   |----1001genomes_snp_biallelic_only_ACGTN.vcf
|   |----master_list.csv
|   |----phenotypes_list.txt
|   |----physics_GWAS_OOP.py
|   |----all_vcf_samples.txt
|   |----JOB_LIST.csv
|   |----TAIR10_GFF3_genes.gff
|   |----output_1.vcf
|   |----output_2.vcf
|   |----output_3.table
|   |----Mo_genes_data.csv
|   |----Na_genes_data.csv
|   |----subsample_text_files/
|       |----subsamples_600_812188.txt (example)
|
|
|----output_files/
|   |----{phenotype}_{GIFT/GWAS}_MANHATTAN_{subsample_num}_{ID}.png 
|   |----{ID}_{subsample_num}_{phenotype}.R
|   |----leaf_ionome_Mo98_whole_genome_metrics_1000_812222.csv (csv example: GIFT data)
|   |----leaf_ionome_Mo98_GWAS_1000_812222.csv (csv example: GWAS data)
|   |----{phenotype}_{GIFT/GWAS}_T20_SNPS_{subsample_num}_{ID}.csv (T20 GWAS csv data)
|   |----R_DATA/
|       |----{phenotype}_{GIFT/GWAS}_{subsample_num}_ALL.csv
|       |----{phenotype}_{positive/negative}_control.csv
|       |----{phenotype}_cumulative_t20_dataframe.csv
|       |----{phenotype}_AVERAGE_{pval_type}_T20_LOCATIONS.csv
|       |----THRESHOLDS.csv
|
|   |----SNP_tracker_R_scripts/
|       |----{phenotype}_{subsample_num}_AVERAGE_{pval_type}_MANHATTAN.R
|       |----{phenotype}_{positive/negative}_control.R
|       |----{phenotype}_cumulative_t20_dataframe.R
|
|   |----summary_plots/
|       |----IDEA1
|           |----{phenotype}_cumulative_t20_dataframe_{pval_type}.png
|           |----{phenotype}_cumulative_t20_dataframe_{pval_type}_KW_TEST.png
|
|       |----IDEA2
|           |----{phenotype}_{positive/negative}_control_{pval_type}.png
|           |----{phenotype}_{positive/negative}_control_{pval_type}_KW_TEST.png
|        
|       |----IDEA3
|           |----{phenotype}_{subsample_num}_AVERAGE_{pval_type}_MANHATTAN.png
|
|       |----stage_3
|           |----{phenotype}_{GIFT/GWAS}_AVERAGE_{pval_type}_ZOOM.png
|
|   |----stage_3_scripts/
|       |----{phenotype}_{GIFT/GWAS}_AVERAGE_{pval_type}_ZOOM.R
|
|   |----GENES_DATA/
|       |----{phenotype}_{GIFT/GWAS}_{subsample_num}_ALL.bed
|       |----Intersect_results{phenotype}_{GIFT/GWAS}_{subsample_num}_ALL.txt
|       |----FINAL_Intersect_results_{phenotype}_{GIFT/GWAS}_{subsample_num}_ALL.txt
|       |----Na23_Gene_Tracker.csv 
|       |----Mo98_Gene_Tracker.csv
|
|
|----slurm0andE
|   |---- slurm_error.err
|   |---- slurm_out.out
```


