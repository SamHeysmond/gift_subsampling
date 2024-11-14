# SCRIPT INFO 
# This script will write the custom R scripts for manhattan plots..
#   .. based on GWAS data, for each JOB_ID run!

import argparse

parser=argparse.ArgumentParser(description="makes R scripts for a given run")

# create the arguments for inputs for this script
parser.add_argument('-m', 
                    type=str, 
                    metavar='method used', 
                    required=True, 
                    help='GIFT or GWAS put after the flag'
                    )
parser.add_argument('-id', 
                    type=str, 
                    metavar='SLURM_JOB_ID', 
                    required=True, 
                    help='Input for the SLURM job ID'
                    )
parser.add_argument('-i', 
                    type=str, 
                    metavar='input subsample num', 
                    required=True, 
                    help='tracks current subsample num e.g. 200 or 400'
                    )
parser.add_argument('-p', 
                    type=str, 
                    metavar='input phenotype', 
                    required=True, 
                    help='tracks the current phenotype e.g. Mo98 or Na23'
                    )
parser.add_argument('-o', 
                    type=str, 
                    metavar='output directory for script', 
                    required=True, 
                    help='where the specific R script from this code will be sent to'
                    )


# stores input data and parses them
args= parser.parse_args() 

if args.m == "GIFT":
    print("GIFT method selected, creating R script",flush=True)

    # this is GIFT but its the file format currently in use (will change later)
    method_variable = "whole_genome_metrics"
    p_value_name="PSNP8"

elif args.m == "GWAS":
    print("GWAS method selected, creating R script",flush=True)
    method_variable = "GWAS"
    p_value_name="pvals"

# make the R file which goes into the given directory (args.o)
Rscript_output=open(str(args.o)+str(args.m)+"_"+str(args.p)+"_"+str(args.i)+"_"+str(args.id)+".R","w")

Rscript_output.write(f'#R script for making manhattan plots with ggplot\n')
Rscript_output.write(f'# info: This script will loop for every job ID to create manhattan plots ..\n')
Rscript_output.write(f'# ... from GWAS or GIFT data\n')
Rscript_output.write(f'# library installs\n')
Rscript_output.write(f'library("tidyverse")\n')
Rscript_output.write(f'library("ggrepel")\n')

# Calculate top 20 SNPs for the given run of GIFT or GWAS
Rscript_output.write(f'Results<-read.csv("output_files/leaf_ionome_{args.p}_{method_variable}_{args.i}_{args.id}.csv",header=TRUE)\n')
Rscript_output.write(f'# get top 20 SNPs from the data\n')
Rscript_output.write(f'T20_SNPS <- Results %>%\n')
Rscript_output.write(f'         select(CHROM,POS,{p_value_name}) %>%\n')
Rscript_output.write(f'         #sort by pvalues lowest first\n')
Rscript_output.write(f'         arrange({p_value_name}) %>%\n')
Rscript_output.write(f'         # take the top 20 values\n')
Rscript_output.write(f'         slice_head(n=20)\n')
Rscript_output.write(f'# write the top 20 snps for this ID to a csv file\n')
Rscript_output.write(f'write.csv(T20_SNPS,"output_files/T20_SNPS_{args.m}_{args.p}_{args.i}_{args.id}.csv",row.names = FALSE)\n')

# threshold calculation for pvals
# change this to instead correct all pvalues by BH or BY procedure in future update
Rscript_output.write(f'# Calculate the Benjamini yekutieli threshold\n')
Rscript_output.write(f'm <- nrow(Results)\n')
Rscript_output.write(f'Results <- Results[order(Results${p_value_name}),]\n')

# Rscript_output.write(f's <- 1.0\n')
# Rscript_output.write(f'i <- 0\n')
# Rscript_output.write('for (p in Results$'+str(p_value_name)+') {\n')
# Rscript_output.write(f'  p\n')
# Rscript_output.write(f'  i <- i+1\n')
# Rscript_output.write('  if (i > 1) {\n')
# Rscript_output.write(f'    s <- s + 1.0/(i-1)\n')
# Rscript_output.write('  }\n')
# Rscript_output.write(f'  thes_pval <- ((i + 1.0) / m) * 0.05 / s\n')
# Rscript_output.write('  if (p > thes_pval) {break\n')
# Rscript_output.write('  }\n')
# Rscript_output.write('}\n')
# Rscript_output.write(f'thes_pval_original <- thes_pval\n')
# Rscript_output.write(f'by_thres <- -log10(thes_pval)\n')
# Rscript_output.write(f'# calculate bonferroni_threshold\n')

# # Should the *1135 change depending on subsample number? if so -> e.g. 200 samples would be *200
# Rscript_output.write(f'bt <- 0.05 / (nrow(Results)*{args.i})\n')
# Rscript_output.write(f'bf_thres <- -log10(bt)\n')

if args.m == "GIFT":
    # 99TH and 95th percentile
    Rscript_output.write(f'NNNH<-quantile(Results$PSNP8,probs=0.01)\n')
    Rscript_output.write(f'NNNH<-unname(NNNH)\n')
    Rscript_output.write(f'NNNH<- -log10(NNNH)\n')
    Rscript_output.write(f'\n')
    Rscript_output.write(f'NFNH<-quantile(Results$PSNP8,probs=0.05)\n')
    Rscript_output.write(f'NFNH<-unname(NFNH)\n')
    Rscript_output.write(f'NFNH<- -log10(NFNH)\n')

# obtain 99% and 95% bonferroni

elif args.m == "GWAS":
    # 99th and 95th bonferoni
    # test apply fix here with Results$pvals
        # not giving anything ...
    # Rscript_output.write(f'NNBF<- -log10(0.01/(nrow(Results$pvals)))\n')
    # Rscript_output.write(f'\n')
    # Rscript_output.write(f'NFBF<- -log10(0.05/(nrow(Results$pvals)))\n')
    # Rscript_output.write(f'\n')
    Rscript_output.write(f'NNNH<-quantile(Results$pvals,probs=0.01)\n')
    Rscript_output.write(f'NNNH<-unname(NNNH)\n')
    Rscript_output.write(f'NNNH<- -log10(NNNH)\n')
    Rscript_output.write(f'\n')
    Rscript_output.write(f'NFNH<-quantile(Results$pvals,probs=0.05)\n')
    Rscript_output.write(f'NFNH<-unname(NFNH)\n')
    Rscript_output.write(f'NFNH<- -log10(NFNH)\n')

# annotation prep stuff
Rscript_output.write(f'#make new columns where the default highlight and annotation is no\n')
Rscript_output.write(f'Results<-Results %>%\n')
Rscript_output.write(f'         mutate(is_highlight = "no")\n')
Rscript_output.write(f'Results<-Results %>%\n')
Rscript_output.write(f'         mutate(is_annotate = "no") \n')
Rscript_output.write(f'# if column CHR and POS match up to a row in T20SNPSthat has same CHR and POS...\n')
Rscript_output.write(f'# ... then mutate the gwas result to flag for highlight and annotation.\n')
Rscript_output.write('for (T20_index in 1:nrow(T20_SNPS)) {\n')
Rscript_output.write('     for (gwas_index in 1:nrow(Results)) {\n')
Rscript_output.write(f'         if (Results$CHROM[gwas_index]== T20_SNPS$CHROM[T20_index] & \n')
Rscript_output.write('          Results$POS[gwas_index]== T20_SNPS$POS[T20_index]){\n')
Rscript_output.write(f'         Results$is_highlight[gwas_index]="yes"\n')
Rscript_output.write(f'         Results$is_annotate[gwas_index]="yes"\n')
Rscript_output.write('          } \n')
Rscript_output.write('      }\n')
Rscript_output.write('}\n')

# cumulative calculations
Rscript_output.write(f'my_data <- Results %>%\n')
Rscript_output.write(f'     # Compute chromosome size\n')
Rscript_output.write(f'     group_by(CHROM) %>% \n')
Rscript_output.write(f'     summarise(chr_len=max(POS)) %>%\n')
Rscript_output.write(f'     # Calculate cumulative position of each chromosome\n')
Rscript_output.write(f'     mutate(tot=cumsum(chr_len)-chr_len) %>%\n')
Rscript_output.write(f'     select(-chr_len) %>%\n')
Rscript_output.write(f'     # Add this info to the initial dataset\n')
Rscript_output.write(f'     left_join(Results, ., by=c("CHROM"="CHROM")) %>%\n')
Rscript_output.write(f'     # Add a cumulative position of each SNP\n')
Rscript_output.write(f'     arrange(CHROM,POS) %>%\n')
Rscript_output.write(f'     mutate( BPcum=POS+tot) \n')
Rscript_output.write(f'axisdf = my_data %>%\n')
Rscript_output.write(f'     group_by(CHROM) %>%\n')
Rscript_output.write(f'     summarize(center=( max(BPcum) + min(BPcum) ) / 2 )\n')
# set y limit for the graph (not sure if this changes anything here though- just copied from other script)
Rscript_output.write(f'ylim <- abs(floor(log10(min(Results${p_value_name})))) +5\n')

Rscript_output.write(f'#open png\n')
# Rscript_output.write(f'png("output_files/{args.p}_GWAS_MANHATTAN_{args.i}_{args.id}.png", bg = "white", width = 9.75, height = 3.25, units = "in", res = 1200, pointsize = 4)\n')
Rscript_output.write(f'png("output_files/{args.m}_{args.p}_{args.i}_{args.id}_MANHATTAN.png", bg = "white", width = 7.5, height = 2.5, units = "in", res = 1200, pointsize = 3)\n')
# make the plot
Rscript_output.write(f'ggplot(my_data, aes(x=BPcum, y=-log10({p_value_name}), color=as_factor(CHROM))) +\n')
Rscript_output.write(f'     # Show all points\n')
Rscript_output.write(f'     geom_point(alpha=0.5) +\n')
# commented out
#Rscript_output.write(f'     geom_point( aes(color=as.factor(CHROM)), alpha=0.5, size=1.3) +\n')

#Rscript_output.write(f'     scale_color_manual(values = rep(c("cyan", "blue"), 22 )) +\n')

# write in the threshold lines
#if args.m == "GIFT":
Rscript_output.write(f'     geom_hline(yintercept = NNNH, color = "red", linetype = "dashed") +\n')
Rscript_output.write(f'     geom_hline(yintercept = NFNH, color = "blue", linetype = "dashed") +\n')
# else:
#     Rscript_output.write(f'     geom_hline(yintercept = NNBF, color = "red", linetype = "dashed") +\n')
#     Rscript_output.write(f'     geom_hline(yintercept = NFBF, color = "blue", linetype = "dashed") +\n')

# Rscript_output.write(f'     geom_hline(yintercept = bf_thres, color = "red", linetype = "dashed") +\n')
# Rscript_output.write(f'     geom_hline(yintercept = by_thres, color = "blue", linetype = "dashed") +\n')

Rscript_output.write(f'     # custom X axis:\n')
Rscript_output.write(f'     scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center ) +\n')

# TEMP PAUSE FOR TESTING
#Rscript_output.write(f'     scale_y_continuous(expand = c(0, 0) ) + # remove space between plot area and x axis\n')

# Highlighting and labelling based on earlier calculations
Rscript_output.write(f'     # add highlighting and labels\n')
Rscript_output.write(f'     geom_point(data = subset(my_data, is_highlight=="yes"), color="orange", size=1.3) +\n')
Rscript_output.write(f'     geom_label_repel(data=subset(my_data, is_annotate=="yes"),\n')
Rscript_output.write(f'         xlim = c(-Inf,Inf), # added -> make room for labels\n')
# paused y lim determiner here for now
#Rscript_output.write(f'         ylim=c(-Inf,Inf), # added -> make room for labels\n')
Rscript_output.write(f'         min.segment.length = 0, # added -> draw lines\n')
Rscript_output.write(f'         max.overlaps = Inf, #added to see what happens\n')
Rscript_output.write(f'         force = 5,\n')
Rscript_output.write(f'         force_pull=0.01,\n')
Rscript_output.write(f'         max.time=2,\n')
Rscript_output.write(f'         aes(label=POS),\n')
#Rscript_output.write(f'         size=1.5) +\n')
Rscript_output.write(f'         size=2.3) +\n')

# label axis
Rscript_output.write(f'     # my axis labels\n')
Rscript_output.write(f'     labs(y= "-log10({p_value_name})", x = "chromosome position")+\n')

# add a theme (may remove in updates)
Rscript_output.write(f'     # Custom the theme:\n')
Rscript_output.write(f'     theme_minimal() +\n')
Rscript_output.write(f'     guides(colour="none")\n')
Rscript_output.write(f'     theme(\n')
Rscript_output.write(f'       panel.border = element_blank(),\n')
Rscript_output.write(f'       panel.grid.major.x = element_blank(),\n')
Rscript_output.write(f'       panel.grid.minor.x = element_blank(),\n')
Rscript_output.write(f'       axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)\n')
Rscript_output.write(f'      )\n')
Rscript_output.write(f'dev.off()\n')
Rscript_output.write(f'#END OF SCRIPT\n')
Rscript_output.close()
print("=======================================================",flush=True)
print("R script for raw data manhattan plot has been created",flush=True)
print("=======================================================", flush=True)
# END OF SCRIPT