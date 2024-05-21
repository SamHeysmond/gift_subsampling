#R script for making manhattan plots with ggplot

# info: This script will loop for every job ID to create manhattan plots
# from GWAS data produced via pygwas (other script)

# library installs
library("tidyverse")
library("ggrepel")



# reads in the list of all job IDs with their respective subsample number
# format (ID, Subsample_N, phenotype) (this file gets written to by subsample_and_run.sh)
JOB_LIST<- read.csv('/gpfs01/home/mbysh17/core_files/JOB_LIST.csv')

# store each column as a list
jobs<-c(JOB_LIST$JOB_ID)
sample_n<-c(JOB_LIST$SAMPLE_N)
phenotype<c(JOB_LIST$PHENOTYPE)

#for loop start
for (i in jobs){
  
#load the data from input of GWAS result
filename <- paste0('/gpfs01/home/mbysh17/output_files/',phenotype[i],'_GWAS_',sample_n[i],'_',jobs[i],'.csv')
gwasResults<-read.csv(filename)

# local testing code commented out
#gwasResults<-read.csv('leaf_ionome_Mo98_GWAS_200_503066.csv')

# get top 20 SNPs from the GWAS data
T20_SNPS <- gwasResults %>%
              select(positions,pvals) %>%
              #sory by pvalues lowest first
              arrange(pvals) %>%
              select(positions)%>%
              slice_head(n=20)
T20_SNPS<- as.list(T20_SNPS)
T20_SNPS

don <- gwasResults %>%
  # Compute chromosome size
  group_by(chromosomes) %>% 
  summarise(chr_len=max(positions)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwasResults, ., by=c("chromosomes"="chromosomes")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chromosomes, positions) %>%
  mutate( BPcum=positions+tot) %>%

  # Add highlight and annotation information
  mutate( is_highlight=ifelse(positions %in% T20_SNPS, "yes", "no")) %>%
  mutate( is_annotate=ifelse(-log10(pvals)>4, "yes", "no")) 

  axisdf = don %>%
  group_by(chromosomes) %>%
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

png("MY_TEST_PLOT.png", bg = "white", width = 9.75, height = 3.25, units = "in", res = 1200, pointsize = 4)
  
ggplot(don, aes(x=BPcum, y=-log10(pvals))) +

    # Show all points
    geom_point( aes(color=as.factor(chromosomes)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("cyan", "blue"), 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$chromosomes, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    
    # Add highlighted points
    geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
    
    # Add label using ggrepel to avoid overlapping
    geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=positions), size=2) +
  
    # my axis labels
    labs(y= "-log10(pvalue)", x = "chromosome position")+
  
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
dev.off()


# for loop end!
}

# OLD CODE BELOW (MAY REMOVE SOON) THAT USED QQman/plot thing
# read the GWAS result file
#data = read.table(filename,sep=',',header=TRUE)
#data

# read the data for T20 snp5 to store information on top 20 SNPs as calculated in gift python script
# T20_SNP5_LOCATION<-paste0('output_files/',i,'_T20_pSNP5.csv')
# T20_SNP5_TABLE<- read.csv=(T20_SNP5_LOCATION)
#SNPs of interest
# POI=c(T20_SNP5_TABLE$POS)

# subset the data and take out top 20 most significant SNPs based on p value
# !!! DO THIS LATER? !!!!

#open png to write to by first making the name
# output_name <- paste0 ('/gpfs01/home/mbysh17/output_files/manhattan_',phenotype[i],'_GWAS_',sample_n[i],'_',jobs[i],'.png')
# png(output_name)

#creating custom SNP column to substitute in for now
# length_ <- dim(data)[1]
# length_
# data$SNP <- paste('SNP',1:length_)

#make the manhattan plot

# manhattan(data, chr="chromosomes",bp="positions",p="pvals", snp="SNP",logp=TRUE,col=c("#1772E0","#17D3E0"))

#close pdf file
# dev.off()



#end of R script
