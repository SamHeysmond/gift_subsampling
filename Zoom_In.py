
# directory to my home dir.
PATH_TO_MAIN = "/gpfs01/home/mbysh17/"

# loop for 
    # 1000 SUBSAMPLE
    # all pval type
    # both phenotypes (Mo98 and Na23)

# set up variables to loop through (the metrics and phenotypes)
pval_type_list=['AVERAGE_P','AVERAGE_PSNP4','AVERAGE_PSNP5','AVERAGE_ABS_THETA']
phenotype_list=["Mo98","Na23"]

# loop through all phenotypes and pval metrics
for phenotype in phenotype_list:

    for pval_type in pval_type_list:

        if pval_type=="AVERAGE_P":
            method_type="GWAS"

        else:
            method_type="GIFT"

        R_out=open(f"{PATH_TO_MAIN}output_files/stage_3_scripts/{phenotype}_{method_type}_{pval_type}_Zoom.R","w")
        R_out.write(f'#R script for zoomed in manhattan plot\n')
        R_out.write(f'library("ggplot2")\n')
        R_out.write(f'print("Start of {method_type} Zoom.R")\n')
        R_out.write(f'\n')
        if method_type=="GIFT":
            R_out.write(f'gwasResults<-read.csv("{PATH_TO_MAIN}output_files/R_DATA/{phenotype}_{method_type}_1000_ALL.csv",header=TRUE)\n')
        elif method_type=="GWAS":
            R_out.write(f'gwasResults<-read.csv("{PATH_TO_MAIN}output_files/R_DATA/{phenotype}_{method_type}_1000_ALL.csv",header=TRUE)\n')
        R_out.write(f'\n')

        # import threshold information from threshold file
        R_out.write(f'threshold_data_csv<-read.csv("{PATH_TO_MAIN}output_files/R_DATA/THRESHOLDS.csv",header=TRUE,sep=",")\n')
        R_out.write(f'threshold_data<-subset(threshold_data_csv,PHENOTYPE==\"{phenotype}\" & PVAL_TYPE==\"{pval_type}\" & SUBSAMPLE_NUM==\"1000\")\n')
        R_out.write(f'\n')
        R_out.write(f'BF_threshold_data<-subset(threshold_data,THRESHOLD_TYPE=="BF" )\n')
        R_out.write(f'BHY_threshold_data<-subset(threshold_data,THRESHOLD_TYPE=="BHY" )\n')
        R_out.write(f'\n')
        R_out.write(f'#filter results to the range desired\n')

        # Positive control gene boundary for Mo98 phenotype
        if phenotype=="Mo98":
            chr_of_interest = 2
            gene_start=10933005
            gene_end = 10934604
            pos_control_gene_name="MOT1"
            reverse_gene = True
            
        # Positive control gene boundary for Na23 phenotype
        elif phenotype=="Na23":
            chr_of_interest = 4
            gene_start=6391854
            gene_end = 6395922
            pos_control_gene_name="HKT1"
            reverse_gene = False
        
        # import boundaries for the gene of interest that we set above
        R_out.write(f'chr_of_interest <-{chr_of_interest}\n')
        R_out.write(f'gene_start = {gene_start}\n')
        R_out.write(f'gene_end = {gene_end}\n')
        R_out.write(f'\n')
        R_out.write(f'bp_start = gene_start - (10*(gene_end-gene_start))\n')
        R_out.write(f'bp_end = gene_end + (10*(gene_end-gene_start))\n')
        R_out.write(f'\n')
        R_out.write(f'threshold_label_xval = (0.25*(bp_end-bp_start)) + bp_start\n')
        R_out.write(f'\n')
        R_out.write(f'df_filtered<-gwasResults[gwasResults$CHR == chr_of_interest & gwasResults$POS >= bp_start & gwasResults$POS <= bp_end, ]\n')
        R_out.write(f'\n')

        if pval_type =="AVERAGE_ABS_THETA":
            # get maximum of raw abs theta value for the peak 
            R_out.write(f'max_snp_val<-gwasResults[which.max(gwasResults${pval_type}),]\n')  

        else:
            # get minimum of raw pval,psnp4,psnp5 value for the peak 
            R_out.write(f'max_snp_val<-gwasResults[which.min(gwasResults${pval_type}),]\n')

        R_out.write(f'\n')
        R_out.write(f'png("{PATH_TO_MAIN}output_files/summary_plots/stage_3/{phenotype}_{method_type}_{pval_type}_ZOOM.png", bg = "white", width = 9.75, height = 3.25, units = "in", res = 1200, pointsize = 4)\n')
        R_out.write(f'\n')

        if pval_type =="AVERAGE_ABS_THETA":
            R_out.write(f'ggplot(df_filtered, aes(x = POS, y = ({pval_type}))) +\n')  

        else:
            R_out.write(f'ggplot(df_filtered, aes(x = POS, y = -log10({pval_type}))) +\n')

        R_out.write(f'  geom_point() +\n')
        R_out.write(f'  theme_minimal() +\n')
        R_out.write(f'  labs(\n')
        R_out.write(f'          title = paste("Zoom in of chromosome:", chr_of_interest, "from", bp_start, "to", bp_end),\n')

        if pval_type =="AVERAGE_ABS_THETA":
            R_out.write(f'          x = "BP", y = "AVERAGE_ABS_THETA")+\n')

            # check direction of gene for direction of arrow
            if reverse_gene == True:
                R_out.write(f'    geom_segment(x = {gene_end}, y = max(df_filtered${pval_type})+5,\n')
                R_out.write(f'    xend ={gene_start} , yend = max(df_filtered${pval_type})+5,\n')

            else:
                R_out.write(f'    geom_segment(x = {gene_start}, y = max(df_filtered${pval_type})+5,\n')
                R_out.write(f'    xend = {gene_end}, yend = max(df_filtered${pval_type})+5,\n')

        else:
            R_out.write(f'  x = "BP", y = "-log10({pval_type})")+\n')

            if reverse_gene == True:
                R_out.write(f'    geom_segment(x = {gene_end}, y = (max(-log10(df_filtered${pval_type}))+5),\n')
                R_out.write(f'      xend = {gene_start}, yend = max(-log10(df_filtered${pval_type}))+5,\n')

            else:
                R_out.write(f'    geom_segment(x = {gene_start}, y = (max(-log10(df_filtered${pval_type}))+5),\n')
                R_out.write(f'      xend = {gene_end}, yend = max(-log10(df_filtered${pval_type}))+5,\n')

        R_out.write(f'      lineend = "round", linejoin = "round",linewidth = 1, arrow = arrow(length = unit(0.3, "inches")),\n')
        R_out.write(f'      colour = "#EC7014")+\n')
        R_out.write(f'\n')
        R_out.write(f'annotate("text",x = (gene_start+gene_end)/2,\n')

        # dont log10 the abs theta values
        if pval_type =="AVERAGE_ABS_THETA":
            R_out.write(f'y = max(df_filtered${pval_type}) + 5,\n')

        else:
            R_out.write(f'y = max(-log10(df_filtered${pval_type})) + 5,\n')
            
        R_out.write(f'label = paste("{pos_control_gene_name}"),color = "blue",hjust = 0.5)+\n')
        R_out.write(f'\n')
        R_out.write(f'annotate("rect", xmin = gene_start, xmax = gene_end, \n')
        R_out.write(f'\n')

        if pval_type =="AVERAGE_ABS_THETA":
            R_out.write(f'ymin = 0, ymax = max(df_filtered${pval_type})+5,\n')

        else:
             R_out.write(f'ymin = 0, ymax = max(-log10(df_filtered${pval_type}))+5,\n')
       
        R_out.write(f'alpha = .2)+\n')
        R_out.write(f'\n')
        R_out.write(f'  geom_vline(aes(xintercept=max_snp_val$POS), color ="red")+\n')
        R_out.write(f'\n')
        R_out.write(f'geom_text(aes(x=max_snp_val$POS -500, label="PEAK", \n')

        if pval_type =="AVERAGE_ABS_THETA":
            R_out.write(f'  y=max_snp_val${pval_type} -5),colour="red", angle=90)\n')
        else:
            R_out.write(f'  y=-log10(max_snp_val${pval_type})-5),colour="red", angle=90)+\n') 
            R_out.write(f' geom_hline(aes(yintercept=BF_threshold_data$THRESHOLD_VALUE, linetype="BF"),\n')
            R_out.write(f'       col = "red")+\n')
            R_out.write(f'\n')
            R_out.write(f' geom_hline(aes(yintercept=BHY_threshold_data$THRESHOLD_VALUE, linetype="BHY"),\n')
            R_out.write(f'      col = "blue")+\n')
            R_out.write(f'\n')
            R_out.write(f'scale_linetype_manual(name = "Thresholds", values = c(2, 2), \n')
            R_out.write(f'      guide = guide_legend(override.aes = list(color = c("red", "blue"))))\n')
        R_out.write(f'\n')
        R_out.write(f'dev.off()\n')
        R_out.close()

# end of script