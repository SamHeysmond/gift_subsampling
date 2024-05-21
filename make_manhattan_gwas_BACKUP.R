#R script for making manhattan plots
library(qqman)

subsample_num<- c(600)

for (i in subsample_num){
filename <- paste0('/gpfs01/home/mbysh17/output_files/leaf_ionome_Mo98_GWAS_',i,'.csv')
#load the data
data = read.table(filename,sep=',',header=TRUE)
data

#open png to write to
output_name <- paste0 ('/gpfs01/home/mbysh17/output_files/manhattan_plots_',i,'.png')
png(output_name)

#creating custom SNP column to substitute in for now
length_ <- dim(data)[1]
length_
data$SNP <- paste('SNP',1:length_)

#make the manhattan plot
manhattan(data, chr="chromosomes",bp="positions",p="pvals", snp="SNP",logp=TRUE,col=c("#1772E0","#17D3E0"))

#close pdf file
dev.off()
}

#end of R script
