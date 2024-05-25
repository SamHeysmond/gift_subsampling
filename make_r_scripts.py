# SCRIPT INFO 
# This script will write the custom R scripts used...
# ... to generate manhattan plots based on GWAS data...
# for each JOB_ID run!


import argparse

parser=argparse.ArgumentParser(description="makes R scripts for a given run")
parser.add_argument('-id', 
                    type=str, 
                    metavar='master list of all phenotype data', 
                    required=True, 
                    help='Input for the master list of phenotype data'
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
                    help='tracks the current phenotype e.g. leaf_ionome_Mo98'
                    )
parser.add_argument('-o', 
                    type=str, 
                    metavar='output directory for script', 
                    required=True, 
                    help='where the specific R script from this code will be sent to'
                    )

# stores input data and parses them
args= parser.parse_args() 

Rscript_output=open(str(args.o)+str(args.id)+"_"+str(args.i)+"_"+str(args.p)+".R",)

Rscript_output.write(f'#R script for making manhattan plots with ggplot\n')
Rscript_output.write(f'# info: This script will loop for every job ID to create manhattan plots\n')
Rscript_output.write(f'# from GWAS data produced via pygwas (other script)\n')
Rscript_output.write(f'# library installs\n')
Rscript_output.write(f'library("tidyverse")\n')
# mark following line for removal in future update
Rscript_output.write(f'#JOB_LIST<- read.csv("/gpfs01/home/mbysh17/core_files/JOB_LIST.csv",header=TRUE)\n')
Rscript_output.write('gwasResults<-read.csv("output_files/{args.p}_GWAS_{args.i}_{args.id}.csv",header=TRUE)\n')
Rscript_output.write(f'# get top 20 SNPs from the GWAS data\n')
Rscript_output.write(f'T20_SNPS <- gwasResults %>%\n')
Rscript_output.write(f'         select(chromosomes,positions,pvals) %>%\n')
Rscript_output.write(f'         #sort by pvalues lowest first\n')
Rscript_output.write(f'         arrange(pvals) %>%\n')
Rscript_output.write(f'         # take the top 20 values\n')
Rscript_output.write(f'         slice_head(n=20)\n')
Rscript_output.write(f'# write the top 20 snps for this ID to a csv file\n')
Rscript_output.write('write.csv("output_files/{args.p}_GWAS_T20_SNPS_{args.i}_{args.id}.csv",row.names = FALSE)\n')
Rscript_output.write(f'#make new columns where the default highlight and annotation is no\n')
Rscript_output.write(f'gwasResults<-gwasResults %>%\n')
Rscript_output.write(f'         mutate(is_highlight = "no")\n')
Rscript_output.write(f'gwasResults<-gwasResults %>%\n')
Rscript_output.write(f'         mutate(is_annotate = "no") \n')
Rscript_output.write(f'# if column CHR and POS match up to a row in T20SNPSthat has same CHR and POS...\n')
Rscript_output.write(f'# ... then mutate the gwas result to flag for highlight and annotation.\n')
Rscript_output.write('for (T20_index in 1:nrow(T20_SNPS)) {\n')
Rscript_output.write('     for (gwas_index in 1:nrow(gwasResults)) {\n')
Rscript_output.write(f'         if (gwasResults$chromosomes[gwas_index]== T20_SNPS$chromosomes[T20_index] & \n')
Rscript_output.write('          gwasResults$positions[gwas_index]== T20_SNPS$positions[T20_index]){\n')
Rscript_output.write(f'         gwasResults$is_highlight[gwas_index]="yes"\n')
Rscript_output.write(f'         gwasResults$is_annotate[gwas_index]="yes"\n')
Rscript_output.write('          } \n')
Rscript_output.write('      }\n')
Rscript_output.write('}\n')
Rscript_output.write(f'don <- gwasResults %>%\n')
Rscript_output.write(f'     # Compute chromosome size\n')
Rscript_output.write(f'     group_by(chromosomes) %>% \n')
Rscript_output.write(f'     summarise(chr_len=max(positions)) %>%\n')
Rscript_output.write(f'     # Calculate cumulative position of each chromosome\n')
Rscript_output.write(f'     mutate(tot=cumsum(chr_len)-chr_len) %>%\n')
Rscript_output.write(f'     select(-chr_len) %>%\n')
Rscript_output.write(f'     # Add this info to the initial dataset\n')
Rscript_output.write(f'     left_join(gwasResults, ., by=c("chromosomes"="chromosomes")) %>%\n')
Rscript_output.write(f'     # Add a cumulative position of each SNP\n')
Rscript_output.write(f'     arrange(chromosomes, positions) %>%\n')
Rscript_output.write(f'     mutate( BPcum=positions+tot) \n')
Rscript_output.write(f'axisdf = don %>%\n')
Rscript_output.write(f'     group_by(chromosomes) %>%\n')
Rscript_output.write(f'     summarize(center=( max(BPcum) + min(BPcum) ) / 2 )\n')
Rscript_output.write(f'#open png\n')
Rscript_output.write('png("output_files/{args.p}_GWAS_MANHATTAN_{args.i}_{args.id}.png", bg = "white", width = 9.75, height = 4, units = "in", res = 1200, pointsize = 4)\n')
Rscript_output.write(f'ggplot(don, aes(x=BPcum, y=-log10(pvals))) +\n')
Rscript_output.write(f'     # Show all points\n')
Rscript_output.write(f'     geom_point( aes(color=as.factor(chromosomes)), alpha=0.8, size=1.3) +\n')
Rscript_output.write(f'     scale_color_manual(values = rep(c("cyan", "blue"), 22 )) +\n')
Rscript_output.write(f'     # custom X axis:\n')
Rscript_output.write(f'     scale_x_continuous( label = axisdf$chromosomes, breaks= axisdf$center ) +\n')
Rscript_output.write(f'     scale_y_continuous(expand = c(0, 0) ) + # remove space between plot area and x axis\n')
Rscript_output.write(f'     # add highlighting and labels\n')
Rscript_output.write(f'     geom_point(data = subset(don, is_highlight=="yes"), color="orange", size=1.3) +\n')
Rscript_output.write(f'     geom_label_repel(data=subset(don, is_annotate=="yes"),\n')
Rscript_output.write(f'         xlim = c(-Inf,Inf), # added -> make room for labels\n')
Rscript_output.write(f'         ylim=c(-Inf,Inf), # added -> make room for labels\n')
Rscript_output.write(f'         min.segment.length = 0, # added -> draw lines\n')
Rscript_output.write(f'         max.overlaps = Inf, #added to see what happens\n')
Rscript_output.write(f'         force = 5,\n')
Rscript_output.write(f'         force_pull=0.01,\n')
Rscript_output.write(f'         max.time=2,\n')
Rscript_output.write(f'         aes(label=positions),\n')
Rscript_output.write(f'         size=1.5) +\n')
Rscript_output.write(f'     # my axis labels\n')
Rscript_output.write(f'     labs(y= "-log10(pvalue)", x = "chromosome position")+\n')
Rscript_output.write(f'     # Custom the theme:\n')
Rscript_output.write(f'     theme_bw() +\n')
Rscript_output.write(f'     theme(\n')
Rscript_output.write(f'         legend.position="none",\n')
Rscript_output.write(f'         panel.border = element_blank(),\n')
Rscript_output.write(f'         panel.grid.major.x = element_blank(),\n')
Rscript_output.write(f'         panel.grid.minor.x = element_blank()\n')
Rscript_output.write(f'         )\n')
Rscript_output.write(f'dev.off()\n')
Rscript_output.write(f'#END OF SCRIPT\n')
Rscript_output.close()
print("============================================")
print("R script for manhattan plot has been created")
print("============================================")
# END OF SCRIPT