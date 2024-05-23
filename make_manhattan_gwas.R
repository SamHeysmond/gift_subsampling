#R script for making manhattan plots with ggplot

# info: This script will loop for every job ID to create manhattan plots
# from GWAS data produced via pygwas (other script)

# library installs
library("tidyverse")
library("ggrepel")

# reads in the list of all job IDs with their respective subsample number
# format (ID, Subsample_N, phenotype) (this file gets written to by subsample_and_run.sh)
JOB_LIST<- read.csv('/gpfs01/home/mbysh17/core_files/JOB_LIST.csv',header=TRUE)

# store each column as a list
#jobs<-c(JOB_LIST$JOB_ID)
#sample_n<-c(JOB_LIST$SAMPLE_N)
#phenotype<c(JOB_LIST$PHENOTYPE)

#for loop start
for (i in 1:nrow(JOB_LIST)){
    
  #load the data from input of GWAS result
  
  filename <- paste0('/gpfs01/home/mbysh17/output_files/',JOB_LIST$PHENOTYPE[i],'_GWAS_',JOB_LIST$SUBSAMPLE_N[i],'_',JOB_LIST$JOB_ID[i],'.csv')
  gwasResults<-read.csv(filename,header=TRUE)
  
  # get top 20 SNPs from the GWAS data
  T20_SNPS <- gwasResults %>%
                select(chromosomes,positions,pvals) %>%
                #sort by pvalues lowest first
                arrange(pvals) %>%
                # take the top 20 values
                slice_head(n=20)
  
  T20_SNPS_FILENAME <- paste0('/gpfs01/home/mbysh17/output_files/',JOB_LIST$PHENOTYPE[i],'_GWAS_T20_SNPS_',JOB_LIST$SUBSAMPLE_N[i],'_',JOB_LIST$JOB_ID[i],'.csv')
  # write the top 20 snps for this ID to a csv file
  write.csv(T20_SNPS,T20_SNPS_FILENAME,row.names = FALSE)
  
  #make new columns where the default highlight and annotation is no
  gwasResults<-gwasResults %>%
    mutate(is_highlight = "no") 
  
  gwasResults<-gwasResults %>%
    mutate(is_annotate = "no") 
  
# if column CHR and POS match up to a row in T20SNPSthat has same CHR and POS...
# ... then mutate the gwas result to flag for highlight and annotation.
  for (T20_index in 1:nrow(T20_SNPS)) {
    for (gwas_index in 1:nrow(gwasResults)) {
      if (gwasResults$chromosomes[gwas_index]== T20_SNPS$chromosomes[T20_index] & 
          gwasResults$positions[gwas_index]== T20_SNPS$positions[T20_index]){
          gwasResults$is_highlight[gwas_index]="yes"
          gwasResults$is_annotate[gwas_index]="yes"
      } 
    }
  }
  
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
    mutate( BPcum=positions+tot) 
  
  axisdf = don %>%
  group_by(chromosomes) %>%
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  pngname <- paste0('/gpfs01/home/mbysh17/output_files/',JOB_LIST$PHENOTYPE[i],'_GWAS_MANHATTAN_',JOB_LIST$SUBSAMPLE_N[i],'_',JOB_LIST$JOB_ID[i],'.png')
  
  png(pngname, bg = "white", width = 9.75, height = 3.25, units = "in", res = 1200, pointsize = 4)
  
  ggplot(don, aes(x=BPcum, y=-log10(pvals))) +

    # Show all points
    geom_point( aes(color=as.factor(chromosomes)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("cyan", "blue"), 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$chromosomes, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) + # remove space between plot area and x axis
    
    # add highlighting and labels
    geom_point(data = subset(don, is_highlight=="yes"), color="orange", size=2) +
    geom_label_repel(data=subset(don, is_annotate=="yes"),
                     xlim = c(-Inf,Inf), # added -> make room for labels
                     ylim=c(-Inf,Inf), # added -> make room for labels
                     min.segment.length = 0, # added -> draw lines
                     max.overlaps = Inf, #added to see what happens
                     aes(label=positions), 
                     size=2) +
    
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
