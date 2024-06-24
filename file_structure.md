file structure for the analysis:
---
```
User Home
|----batch_files
|   |----conda_environment_setup.sh
|   |----make_r_scripts.py
|   |----Rscript_run_and_monitor.py
|   |----run_and_monitor.py
|   |----run_and_monitor.sh
|   |----SNP_tracker.py
|   |----SNP_tracker.sh
|   |----subsample_script_setup.py
|   |----subsample_setup.sh
|   |----parallel/
|       |----Stage_1_subrun_scripts e.g. subrun_leaf_ionome_Na23_1000_100.sh (before complete)
|   |----completed_parallel
|       |---- Stage_1_subrun_scripts e.g. subrun_leaf_ionome_Na23_1000_100.sh (after complete)
|   |----parallel_stage2/
|       |----IDEA_1_Script_before_run e.g. Mo98_cumulative_t20_dataframe.sh
|       |----IDEA_2_Script_before_run e.g. Na23_negative_control.sh
|       |----IDEA_3_Script_before_run e.g. Mo98_800_AVERAGE_P_MANHATTAN.sh
|   |----completed_parallel_stage2/
|       |----IDEA_1_Script_after_run e.g. Mo98_cumulative_t20_dataframe.sh
|       |----IDEA_2_Script_after_run e.g. Na23_negative_control.sh
|       |----IDEA_3_Script_after_run e.g. Mo98_800_AVERAGE_P_MANHATTAN.sh
|
|----core_files/
|   |----1001genomes_snp_biallelic_only_ACGTN.vcf
|   |----master_list.csv
|   |----phenotypes_list.txt
|   |----physics_GWAS_OOP.py
|   |----all_vcf_samples.txt
|   |----JOB_LIST.csv
|   |----subsample_text_files/
|       |----subsamples_600_812188.txt (example)
|
|----output_files/
|   |----leaf_ionome_Mo98_GWAS_MANHATTAN_1000_812222.png (manhattan plot example)
|   |----812222_1000_leaf_ionome_Mo98.R (R script example)
|   |----leaf_ionome_Mo98_whole_genome_metrics_1000_812222.csv (csv example: GIFT data)
|   |----leaf_ionome_Mo98_GWAS_1000_812222.csv (csv example: GWAS data)
|   |----leaf_ionome_Mo98_GWAS_T20_SNPS_1000_812222.csv(T20 GWAS csv data)
|   |----R_DATA/
|   |----SNP_tracker_R_scripts/
|   |----summary_plots/
|
|+ file 3) 4 Manhattan plots for each subrun (1 for GWAS and 3 for GIFT) e.g. leaf_ionome_Mo98_GWAS_MANHATTAN_1000_812222.png
+ file 4) 4 R scripts per subrun that make each of the manhattan plots 
+ file 5) 4 csv files per subrun that contain data on the top 20 most significant SNPs
+ file 6) 2 csv files per subrun that contain all the SNP data analysed (1 for GIFT, 1 for GWAS)
|
|
|
|----slurm0andE
|   |---- slurm_error.err
|   |---- slurm_out.out
|
|
|
```


==

packages/button
├── lib
│   ├── button.d.ts
│   ├── button.js
│   ├── button.js.map
│   ├── button.stories.d.ts
│   ├── button.stories.js
│   ├── button.stories.js.map
│   ├── index.d.ts
│   ├── index.js
│   └── index.js.map
├── package.json
├── src
│   ├── button.stories.tsx
│   ├── button.tsx
│   └── index.ts
└── tsconfig.json