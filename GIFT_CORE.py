#! /usr/bin/python3
# Testing command:
# cd ~
# python3 /Users/sian_bray/Dropbox/Bray/Research/Physics_GWAS/Scripts/physics_GWAS_OOP_V2.py -v /Users/sian_bray/Dropbox/Bray/Research/Physics_GWAS/Data/MOT1_Test/Mot1_biallelic_10kb_AF_filtered.vcf -f /Users/sian_bray/Dropbox/Salt/2_Data/Phenotype_Data/master_list.csv -p leaf_ionome_Mo98 -o ~/physics_gwas_test.csv

# Real commands:
# cd /Users/sian_bray/Dropbox/Bray/Research/Physics_GWAS/Data/Testing
# python3 /Users/sian_bray/Dropbox/Bray/Research/Physics_GWAS/Scripts/physics_GWAS_OOP_V2.py -v /Users/sian_bray/Dropbox/Bray/Research/Physics_GWAS/Data/VCFs/1001genomes_snp-biallelic_only_ACGTN.vcf -f /Users/sian_bray/Dropbox/Salt/2_Data/Phenotype_Data/master_list.csv -p leaf_ionome_Mo98 -o /Users/sian_bray/Dropbox/Bray/Research/Physics_GWAS/Data/Testing/Mo_whole_genome_metrics.csv
# python3 /Users/sian_bray/Dropbox/Bray/Research/Physics_GWAS/Scripts/physics_GWAS_OOP_V2.py -v /Users/sian_bray/Dropbox/Bray/Research/Physics_GWAS/Data/VCFs/1001genomes_snp-biallelic_only_ACGTN.vcf -f /Users/sian_bray/Dropbox/Salt/2_Data/Phenotype_Data/master_list.csv -p leaf_ionome_Na23 -o /Users/sian_bray/Dropbox/Bray/Research/Physics_GWAS/Data/Testing/Na_whole_genome_metrics.csv

# python3 /Users/sian_bray/Dropbox/Bray/Research/Physics_GWAS/Scripts/physics_GWAS_OOP_V2.py -v /Users/sian_bray/Dropbox/Bray/Research/Physics_GWAS/Data/VCFs/1001genomes_snp-biallelic_only_ACGTN.vcf -f /Users/sian_bray/Dropbox/Salt/2_Data/Phenotype_Data/master_list.csv -p leaf_ionome_Mo98 -o /Users/sian_bray/Dropbox/Bray/Research/Physics_GWAS/Data/Testing/Mo_p_value.csv
# python3 /Users/sian_bray/Dropbox/Bray/Research/Physics_GWAS/Scripts/physics_GWAS_OOP_V2.py -v /Users/sian_bray/Dropbox/Bray/Research/Physics_GWAS/Data/VCFs/1001genomes_snp-biallelic_only_ACGTN.vcf -f /Users/sian_bray/Dropbox/Salt/2_Data/Phenotype_Data/master_list.csv -p leaf_ionome_Na23 -o /Users/sian_bray/Dropbox/Bray/Research/Physics_GWAS/Data/Testing/Na_p_value.csv

# ## = comments for Jon
# ### = responses from Jon
# #? = to do
#P# = publication reference

import argparse, random, pandas, os
from math import sqrt, erfc, log10
import numpy as np

# Creates a list of individual names ordered by phenotype value.
# Input is a csv file with two columns; individual and phenotype value.
# A headder is expected.
# Output is a list of individuals.
def ordered_list(phenotypes_file, phenotype):

	# Import phenotype file as an array and sort by the size of the second column ()
	pheno_array=pandas.read_csv(phenotypes_file)
	pheno_array=pheno_array.sort_values(by=[phenotype]) # Note: sorted smallest first to largest last
	# Make a list of the sample names, ordered by value
	pheno_order=pheno_array.iloc[:,0].tolist()
	#? Need to filter out identical values here so I don't have to do it manually
	print("pheno_order:",flush=True)
	print(pheno_order,flush=True)
	return pheno_order #P# Jon 2.1

# Takes a vcf file and generates a list of the headder items for later use.
def read_vcf_head(vcf_file):

	file=open(vcf_file, 'r')
	for line in file:
		if '#CHROM' in line:
			newline = line.replace('\n','')
			vcf_order=newline.split('\t')
			#vcf_order=line.split('\t')
			return vcf_order
	file.close()

# Calculate multiple test corrections - from pygwas mtcorr.py

# Implements Benjamini-Hochberg FDR threshold (1995)
# Input is a list of p-values and a threshold
def get_bh_thres(pvals, fdr_thres=0.05):

    m = len(pvals)
    s_pvals = sorted(pvals) 
    for i, p in enumerate(s_pvals):
        thes_pval = ((i + 1.0) / float(m)) * fdr_thres
        if p > thes_pval:
            break
        
    return {'thes_pval':thes_pval, 'thres_i':i}
        
# Implements the Benjamini-Hochberg-Yekutieli procedure (2001).
# Assumes arbitrary dependence between variables.
# Input is a list of p-values and a threshold
def get_bhy_thres(pvals, fdr_thres=0.05):

    m = len(pvals)
    m_float = float(m)
    s_pvals = sorted(pvals)
    s = 1.0
    for i, p in enumerate(s_pvals):

		# is it >2 because used enumerate here? or is it diff to R script for other reason?
        if i > 2:
            s = s + 1.0/(i-1)
        thes_pval = ((i + 1.0) / m_float) * fdr_thres / s
        if p > thes_pval:
            break
    return {'thes_pval':thes_pval, 'thres_i':i}


# Create a class 'field'
# This class takes each VCF line and processes it (turns it into an ordered list of genotype values)
# You can then perform different methods on each line
class field:
	def __init__(self, line, pheno_order, vcf_order):
		self.line=line
		self.pheno_order=pheno_order
		self.vcf_order=vcf_order

		# process the line
		line=line.replace('\n', '')
		line=line.split('\t')
		genotype_values=[] # Empty list for the -1, 0 and 1 values

		# Save and split the format field.
		line_format=line[8]
		line_format=line_format.split(':')

### SAM EDIT ############################
##### Store the chromosome and position together as a string variable "chr:pos" 
# this variable will be the column header for the genotype tracker
		current_chromosome = str(line[0])
		current_position = str(line[1])
		# e.g. 2:142 (chromosome 2: position 142)
		position_variable = current_chromosome+":"+current_position
####################################################################################

		# For each individual in the phenotype file, locate the genotype and give it a value between -1 and +1 (biallelic diploid data only).
		# In this case the individuals are those left behand after subsampling
		# this entire for loop is called on each position in the vcf (i think)
		# pheno_order is a list of sample IDs ordered smallest to largest in phenotype
		for individual in pheno_order:

			individual=str(individual)

			try: # Get the location index of the ID of the individual in the vcf. e.g. index 24
				# starts with the smallest individual for that phenotype since going through pheno_order
				index=vcf_order.index(individual)
			except ValueError: # if the individual is not in the vcf
				print(f"individual {individual} was NOT in the VCF ",flush=True)
				continue

			current_column=line[index] # Get the column that contains genotypes for this individual.
			current_column=current_column.split(':') # Split into fields.
			current_genotypes=current_column[line_format.index('GT')] # Find the genotypes for this individual.

			#  This changes unphased symbol to phased symbol
			current_genotypes=current_genotypes.replace('/','|') # Make sure you split phased and/or unphased data.
			
			current_genotypes=current_genotypes.split('|')
			current_value=[]
			total=0

			# Make a list of all the genotype values for that individual

			# example looking at column 10020 on row (line) (most confident):
			# 1|1 ->  sum =2  total = 2 -> current_value =

			for things in current_genotypes:
				if things == "0" or things == "1":
					current_value.append(int(things))
					total+=1
				else:
					print(f"Things was missing in genotypes {current_genotypes}",flush=True)

			if current_value != []: # If there is any genotype data (sum of an empty list is 0)
				current_value=sum(current_value)
				current_value=current_value/total
				# current_value=round(current_value * 2) / 2 # round to the nearest 0.5 in case of pop-level genotyping or polyploids

				# does it still round somehow?

				# case of dominant and recessive alleles example:
				# AA = +1, Aa = 0, aa = -1
					# should only ever be +1 or -1 when dealing with biallelic sites?
				# not sure on how it lands at either 0, 0.5 or 1
					# does it only look at one column at a time? i.e. only ever looks at 0|0 , 0|1 or 1|1
				if current_value == 0:
					genotype_values.append(-1)
				if current_value == 0.5:
					genotype_values.append(0)
				if current_value == 1:
					genotype_values.append(1)
			else:
				print(f"Value was MISSING for individual: {individual}",flush=True)

				# TEMP FOR TESTING
				print(f"Temp value given:0 {individual}",flush=True)
				genotype_values.append(0)
				
### SAM EDIT ############################
##### Once done adding genotypes for a given SNP for each individual, output this to a table

		# test print
		#print(f"Genotype values: {genotype_values}",flush=True)
		#print(f"Genotype values length: {len(genotype_values)}",flush=True)

		#append the genotype values as a new column with the position in the vcf as the header (instead of snp3 etc)
		genotype_tracker_df[f'{position_variable}'] = genotype_values
		
		#test print
		#print("genotype_tracker_df updated: ",  flush=True)
		#print(genotype_tracker_df, flush=True)

####################################################################################
		self.ordered_states=genotype_values #P# Jon 2.2
		
		# obtain length of all the values (the -1's, 0's and 1's)
		N = len(genotype_values)

		# obtain counts for each of the types of values
		N_plus = genotype_values.count(1)
		N_minus = genotype_values.count(-1)
		N_zero = genotype_values.count(0)

		self.N=N
		self.N_plus=N_plus
		self.N_minus=N_minus
		self.N_zero=N_zero

		self.line=line
		#? can maybe remove everything but ordered states when finished for efficency

	# Filter here (doesnt this count for 30 not 15 since its 15 min on - and on + states?)
	# changed this to be more lenient
	def sense_check(self, min_SNPs=15): 
		ordered_states = self.ordered_states
		if ordered_states: # Check that there is something in ordered_states

			# counts number of +1 states
			count_plus=ordered_states.count(1)

			# counts number of -1 states
			count_minus=ordered_states.count(-1)
			if count_plus >= min_SNPs and count_minus >= min_SNPs:
				return True
		return False # If you didn't manage to return True

	# Calculate the cumulative sum of the ordered genotype states
	def cum_sum(self): 

		return list(np.cumsum(self.ordered_states))

	# Calculate the straight path i.e. a straigh line from the start to the end of the cumulative sum
	def straight_path(self): 
		ordered_states = self.ordered_states
		final_dest=sum(ordered_states)
		per_step=final_dest/len(ordered_states)
		straight_line=[]
		for step, num in enumerate(ordered_states, 1):

			#samq makes a straight line by going from point A to point B in equal steps?
			straight_line.append(per_step*step)
			
		return straight_line #P# Jon 2.7

	# Calculate a random path
	def rand_sum(self): 

		genotype_values=self.ordered_states
		random.shuffle(genotype_values)
		return np.cumsum(genotype_values)

	# Calculate theta for each position in the ordered states #? Could do to make this re-usable for values that are not the ordered state i.e. a given random path
	# SAM EDIT: Changed from false to True!
	def calc_theta(self, J_2_8=True): 

		# Create the stright line (i.e. from 0 to to sum of the values):
		# straight line is the (total sum / steps) * the step
		#final_dest=sum(self.ordered_states)
		#try:
		#	per_step=final_dest/len(self.ordered_states)
		#except ZeroDivisionError:
		#	per_step=0
		#cumulative_sum=list(np.cumsum(self.ordered_states))
		#straight_line=[]
		#for step, num in enumerate(cumulative_sum, 1):
		#	straight_line.append(per_step*step)

		# Creates the straight_line which is the same as the straight path (see straight_path() function above)
		straight_line=self.straight_path()
		cumulative_sum=self.cum_sum()

		N=self.N
		N_plus=self.N_plus
		N_minus=self.N_minus
		N_zero=self.N_zero

		#P# Jon 2.8
		if J_2_8 == True:
			theta_plus = []
			theta_minus = []
			W_plus_j = 0 # number of plus states at j
			W_minus_j = 0 # number of minuss states at j
			for j, state in enumerate(cumulative_sum):
				# j + 1 to make 1-based
				if state == 1: 
					W_plus_j += 1
				if state == -1:
					W_minus_j += 1
				theta_plus.append(W_plus_j - (((j+1)*N_plus)/N))
				theta_minus.append(W_minus_j - (((j+1)*N_minus)/N))

#################################################
#### SAM EDIT FOR UPDATE TO PSNP8 ###############
				#theta 0 (random case) (same as for thetaJ) 
				#second half of 253 - 254
				theta_0 = (((j+1)*N_plus)/N) - (((j+1)*N_minus)/N)
				# theta_zero can be calculated from these two, see Jon's paper

#################################################
#### SAM EDIT FOR UPDATE TO PSNP8 ###############

			# for Theta J
			# w+ - w-
			theta_j = theta_plus - theta_minus

			# for delta theta
				# difference between Theta J and Theta 0 (random case)
			delta_theta = theta_j - theta_0
			
#################################################
#################################################
		# Calculate theta
		theta=[]
		for count, val in enumerate(cumulative_sum):
			# append the value of the difference from the straight line for each count value i.e. each position on the line of the cumulative sum
			theta.append(val-straight_line[count])

		# Calculate relative theta

		# largest and smallest raw theta 
		# samq where the line differes from the straight line the most and the least?)
		largest_theta=max(theta)
		smallest_theta=min(theta)

		# Calculate relative theta
		# newtheta_plus =oldtheta_plus times N / ( Nplus times (Nzero + Nminus ) )
		# should give a newtheta between -1 and +1
		largest_relative_theta = largest_theta * N / (N_plus * (N_zero + N_minus))
		smallest_relative_theta = smallest_theta * N / (N_plus * (N_zero + N_minus))

		# generate largest, smallest, range, relative and absolute theta
		if abs(largest_theta) > abs(smallest_theta):
			absolute_theta=abs(largest_theta)
			absolute_relative_theta=abs(largest_relative_theta)
		else:
			absolute_theta=abs(smallest_theta)
			absolute_relative_theta=abs(smallest_relative_theta)

		theta_range=largest_theta-smallest_theta
		range_relative_theta=largest_relative_theta-smallest_relative_theta

		#head: largest_theta,smallest_theta,absolute_theta,theta_range,largest_relative_theta,smallest_relative_theta,absolute_relative_theta,range_relative_theta
		return largest_theta,smallest_theta,absolute_theta,theta_range,largest_relative_theta,smallest_relative_theta,absolute_relative_theta,range_relative_theta

	def calc_pvals(self): ## Jon - this is where I calculate the P-values

		# N = Total states
		# theta_j = k - j * w_0_plus
		# j = first j elements in list
		# N_plus = number of plus states
		# N_minus = number of minus states
		# k = number of plus states in the first j elements
		# w_0_plus = frequency of plus states
		# j * w_0_plus = expected + states at j

		ordered_states = self.ordered_states

		# Calculate number of states: N, N_plus and N_minus
		N, N_plus, N_minus, N_zero = self.N, self.N_plus, self.N_minus, self.N_zero ## Jon - I use the function on line 111 to calculate all the N's

		# Create empty tuple of P values
		p_vals = []
		bigest_theta = None
		bigest_theta_j = None

		# Calculate straight line and cumulative sum
		cumsum = self.cum_sum() ## Jon - I use the function on line 128 to calculate the cumulative sum
		straight_line = self.straight_path() ## Jon - I use the function on line 132 to calculate the straight path

		# To keep a theta list
		theta = []
		# For the first j elements in the ordered list
		# samq ordered states that occur when the phenotype is ordered from smallest to biggest?
		for j, state in enumerate(ordered_states, 1):

			# Calculate theta at j, j - 1 to make it zero-based
			theta_j = cumsum[j-1] - straight_line[j-1]
			theta.append(theta_j)

			if bigest_theta == None:
				bigest_theta = theta_j
				bigest_theta_j = int(j-1)

			if theta_j > bigest_theta: ### Jon's response:  (lines 244-246 - we can ignore - my ideas/calculations regarding the biggest value of theta have been superseded by the pSNPmin, pSNP2/pSNP3/pSNP4, / pSNP5 )
				bigest_theta = theta_j
				bigest_theta_j = int(j-1)

			# Calculate p at j
			# But not for the last line due to theta_j and N-j both being zero!
			if theta_j == N:
				p_vals.append(1)
			else:
				p = erfc((N*sqrt(N)*theta_j)/(sqrt(2*j*(N-j)*N_plus*N_minus))) ## Jon - I calculate the raw p value here
				p_vals.append(p)

		# temporary nan or 0 as last value fix:
		#p_vals[-1]=1 ## Jon - the last value almost always becomes 0 or 'nan' so I am just removing it for now ### Jon's response: trying to evaluate the last item in the list will give 'erfc(0/0)', so setting this to 1 is the correct outcome (or just ignore it, and do later averaging only over elements  1,2,3, .... j .... N-1)

		# Calculate values to return:
		min_p = min(p_vals) # smallest p ## Jon - taking the lowest overall p value
		mean_p = sum(p_vals) / len(p_vals) # average p ## Jon - taking the average p value

		# -log10 average ## Jon - taking the -log10 average p value
		log_p_vals=[]
		for no in p_vals:
			try:
				log_p_vals.append(-log10(no))
			except ValueError: ## Jon - this was an earlier escape for any 'nan' values
				log_p_vals.append(0)
		log_mean_p = sum(log_p_vals) / len(log_p_vals)

		# Greatest theta p_value
		bigest_theta_p = p_vals[bigest_theta_j] ## Jon - p value for the largest theta, bigest_theta_j is the index position of the biggest theta value
		# sigma for pSNP4, j - 1 to make it zero-based ## Jon - below are the calculations fos pSNP4 and pSNP5, x += x means x = x+x, so e.g. 5 += 5 is 10
		sigma_4 = 0
		for j, state in enumerate(ordered_states, 1):
			if j < N:
				sigma_4 += abs(theta[j-1]) / sqrt(j*(N-j))

		# pSNP4
		Z = sigma_4 * sqrt(N)/sqrt(2*N_plus*N_minus)
		pSNP4 = erfc(Z)

		# sigma for pSNP5, j - 1 to make it zero-based
		sigma_5 = 0
		for j, state in enumerate(ordered_states, 1):
			if j < N:
				sigma_5 += (theta[j-1])/ sqrt(j*(N-j)) 
		
		# pSNP5
		Z = sigma_5 * sqrt(N)/sqrt(2*N_plus*N_minus)
		pSNP5 = erfc(abs(Z))

		# return the smallest P value
		return min_p, mean_p, log_mean_p, bigest_theta_p, pSNP4, pSNP5


# Make R plots, requires the R libraries; tidyverse, htmlwidgets and manhattanly
def R_plots(output_file, metric='largest_relative_theta', p_value=True): #?# write if p_value == True to -log10 and bf/bhy correct
	R_out=open(f'{output_file[:-4]}_{metric}.R', 'w')
	png_out=f'{output_file[:-4]}_{metric}.png'
	html_out=f'{output_file[:-4]}_{metric}.html'
	# write the R Script
	R_out.write(f'library("tidyverse")\n')
	R_out.write(f'library("htmlwidgets")\n')
	R_out.write(f'library("ggrepel")\n')
	R_out.write(f'#library("manhattanly")\n')

	# does this need header and separator flags?
	R_out.write(f'GWAS_result1 <- read.csv("{output_file}")\n')

	# Issues with NA encountered- TEMP PAUSE THIS COMMAND
	#R_out.write(f'GWAS_result1[GWAS_result1==Inf] <- NA\n')
	#R_out.write(f'GWAS_result1[GWAS_result1==-Inf] <- NA\n')
	#R_out.write(f'GWAS_result<-GWAS_result1[complete.cases(GWAS_result1),]\n')
	R_out.write(f'GWAS_result <-GWAS_result1\n')

	# 'absolute_theta'
	# ### start of sam edit (1) <<<<<<<<<<<<<<<<<<<<< # moved this stuff under p_value ==True only

	#check if metric is absolute theta or a p value (psnp4 or 5)
	if metric == 'absolute_theta':
		# set p value flag to FALSE
		p_value = False

		# reads the T20 snp file based on which metric is being graphed
		snps_R_variable=str("T20_Absolute_theta_SNPs")
		R_out.write(f'T20_Absolute_theta_SNPs <- read.csv("output_files/'+args.id+'_T20_absolute_theta.csv", header= TRUE, sep=",")\n')
	
	elif metric =="pSNP4":
		snps_R_variable=str("T20_pSNP4_SNPs")
		R_out.write(f'T20_pSNP4_SNPs <- read.csv("output_files/'+args.id+'_T20_pSNP4.csv", header= TRUE, sep=",")\n')
		p_value = True

	elif metric =="pSNP5":
		snps_R_variable=str("T20_pSNP5_SNPs")
		R_out.write(f'T20_pSNP5_SNPs <- read.csv("output_files/'+args.id+'_T20_pSNP5.csv", header= TRUE, sep=",")\n')
		p_value = True
	### ADD IN CODE TO READ A CSV FOR CUSTOM SNPS OF INTEREST MUCH LIKE THE ABOVE BUT SEPARATE!!
		#code here for any SNPs in particular (but will require a file to be made first for it to read from)
	# ### end of sam edit (1) >>>>>>>>>>>>>>>>>>>>>>>>>>>

	if p_value == True:
		
		R_out.write(f'# Calculate the BHY threshold\n')
		R_out.write(f'm <- nrow(GWAS_result)\n')
		R_out.write(f'GWAS_result <- GWAS_result[order(GWAS_result${metric}),]\n')
		R_out.write(f's <- 1.0\n')
		R_out.write(f'i <- 0\n')
		R_out.write('for (p in GWAS_result$'+metric+') {\n')
		R_out.write(f'  p\n')
		R_out.write(f'  i <- i+1\n')
		R_out.write('  if (i > 1) {\n')
		R_out.write(f'    s <- s + 1.0/(i-1)\n')
		R_out.write('  }\n')
		R_out.write(f'  thes_pval <- ((i + 1.0) / m) * 0.05 / s\n')
		R_out.write('  if (p > thes_pval) {break\n')
		R_out.write('  }\n')
		R_out.write('}\n')
		R_out.write(f'thes_pval_original <- thes_pval\n')
		# ### start of sam edit (1.2) <<<<<<<<<<<<<<<<<<<<< 
		R_out.write(f'bhy_thres <- -log10(thes_pval)\n')
		R_out.write(f'# calculate bonferroni_threshold\n')

		# Should the *1135 change depending on subsample number? if so -> e.g. 200 samples would be *200
		R_out.write(f'bt <- 0.05 / (nrow(GWAS_result)*{args.s}) # times max number of tests per p-value (subsample number)\n')
		R_out.write(f'bf_thres <- -log10(bt)\n')
		# ### END of sam edit (1.2) <<<<<<<<<<<<<<<<<<<<< 
	R_out.write(f'data_cum <- GWAS_result %>% \n')
	R_out.write(f'  group_by(CHROM) %>% \n')
	R_out.write(f'  summarise(max_bp = max(POS)) %>% \n')
	R_out.write(f'  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% \n')
	R_out.write(f'  select(CHROM, bp_add)\n')
	R_out.write(f'GWAS_result <- GWAS_result %>% \n')
	R_out.write(f'  inner_join(data_cum, by = "CHROM") %>% \n')
	R_out.write(f'  mutate(bp_cum = POS + bp_add)\n')
	# ### start of sam edit (2) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	# first add column for annotation and highlighting, setting each to no (default)
	R_out.write(f'GWAS_result<- GWAS_result %>% \n')
	R_out.write(f'  mutate(is_highlight = "no")\n')
	R_out.write(f'GWAS_result<- GWAS_result %>% \n')
	R_out.write(f'  mutate(is_annotate = "no")\n')
	# if column CHR and POS match up to a row in POI that has same CHR and POS...
	# then mutate the gwas result to highlight and annotate the data
	R_out.write('for (T20_index in 1:nrow('+snps_R_variable+')){ \n') 
	R_out.write('	for(gwas_index in 1:nrow(GWAS_result)) { \n')
	R_out.write('		if (GWAS_result$CHROM[gwas_index]=='+snps_R_variable+'$CHROM[T20_index] & \n')
	R_out.write('			GWAS_result$POS[gwas_index]== '+snps_R_variable+'$POS[T20_index]){ \n')
	R_out.write('			GWAS_result$is_highlight[gwas_index]="yes" \n')
	R_out.write('			GWAS_result$is_annotate[gwas_index]="yes" \n')
	R_out.write('		} \n')
	R_out.write('	} \n')
	R_out.write('} \n')
	# ### end of sam edit (2) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	R_out.write(f'axis_set <- GWAS_result %>% \n')
	R_out.write(f'  group_by(CHROM) %>% \n')
	R_out.write(f'  summarize(center = mean(bp_cum))\n')
	R_out.write(f'ylim <- abs(floor(log10(min(GWAS_result${metric})))) +1\n')
	R_out.write(f'png("{png_out}", bg = "white", width = 9.75, height = 3.25, units = "in", res = 1200, pointsize = 4)\n')
	
	if p_value == True: 
		# ### START of sam edit (2.2) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
		R_out.write(f'manhplot <- ggplot(GWAS_result, aes(x = bp_cum, y = (-log10({metric})), # size = 1, \n')
		# ### ENDof sam edit (2.2) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	if p_value == False:
		R_out.write(f'manhplot <- ggplot(GWAS_result, aes(x = bp_cum, y = ({metric}), # size = 1, \n')

	R_out.write(f'                                  color = as_factor(CHROM))) +\n')
	R_out.write(f'  geom_point(alpha = 0.5) +\n')
	
	if p_value == True:
		R_out.write(f'geom_hline(yintercept = bf_thres, color = "red", linetype = "dashed") +\n')
		R_out.write(f'geom_hline(yintercept = bhy_thres, color = "blue", linetype = "dashed") +\n')
	
	R_out.write(f'  scale_x_continuous(label = axis_set$CHROM, breaks = axis_set$center) +\n')
	# ### start of sam edit (3) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	#highlight in orange
	R_out.write(f'  geom_point(data = subset(GWAS_result, is_highlight=="yes"), color="orange", size=2) +\n')
	# add a label using ggrepel
	R_out.write(f'  geom_label_repel(data=subset(GWAS_result, is_annotate=="yes"),\n')
	R_out.write(f'					xlim=c(-Inf,Inf), \n')
	R_out.write(f'					ylim=c(-Inf,Inf), \n')
	R_out.write(f'					min.segment.length=0, \n')
	R_out.write(f'					max.overlaps = Inf, \n')
	R_out.write(f' 					aes(label=POS), \n')
	R_out.write(f'					size=2) + \n')
	# alt idea using geom label but havent tried it yet
	#R_out.write(f'  geom_label(data = subset(GWAS_result, is_highlight=="yes")) +\n')
	# ### end of sam edit (3)>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	R_out.write(f'  labs(x = NULL, \n')
	R_out.write(f'       y = "{metric}") + \n')
	R_out.write(f'  theme_minimal() +\n')
	R_out.write(f'  guides(colour="none")\n')
	R_out.write(f'  theme(\n')
	R_out.write(f'    panel.border = element_blank(),\n')
	R_out.write(f'    panel.grid.major.x = element_blank(),\n')
	R_out.write(f'    panel.grid.minor.x = element_blank(),\n')
	R_out.write(f'    axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)\n')
	R_out.write(f'  )\n')
	R_out.write(f'print(manhplot)\n')
	R_out.write(f'dev.off()\n')
	R_out.write(f'#attach(GWAS_result)\n')
	R_out.write(f'#GWAS_result <- GWAS_result[order(-{metric}),]\n')
	R_out.write(f'#detach(GWAS_result)\n')
	R_out.write(f'#top_perc <- (nrow(GWAS_result)/100) * 0.01\n')
	R_out.write(f'#top_perc <-round(top_perc, digits = 0)\n')
	R_out.write(f'#GWAS_result <- head(GWAS_result,top_perc)\n')
	R_out.write(f'#names(GWAS_result)[names(GWAS_result) == "CHROM"] <- "CHR"\n')
	R_out.write(f'#GWAS_result$P <- (0.0000001/GWAS_result${metric})\n')
	R_out.write(f'#names(GWAS_result)[names(GWAS_result) == "POS"] <- "BP"\n')
	R_out.write(f'#html_file <- manhattanly(GWAS_result, annotation1 = "CHR", annotation2 = "BP")\n')
	R_out.write(f'#saveWidget(html_file, "{html_out}", selfcontained = T)\n')
	# Close the files!
	R_out.close()
	# Run the R code!
	to_run=f'Rscript {output_file[:-4]}_{metric}.R'
	print(to_run)
	os.system(to_run)

# Run the code!
if __name__ == '__main__':

	# File input
	parser = argparse.ArgumentParser(description="Sums all DPs (depths) in a vcf.")
	parser.add_argument('-v', type=str, metavar='input_vcf', required=True, help='The input vcf file.')
	parser.add_argument('-f', type=str, metavar='phenotypes_file', required=True, help='The input phenotype file. It is a .csv file with a one line headder and two columns: individual and phenotype. Phenotype is a numeric value and they must be in order of size.')
	parser.add_argument('-p', type=str, metavar='phenotype', required=True, help='The phenotype - exactly as written in the phenotype file headder.')
	parser.add_argument('-o', type=str, metavar='output_file', required=True, help='The output file.')
	# ### start of sam edit (4) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	# add in an ID tag to allow for proper file naming of outputs and tracking of run IDs
	parser.add_argument('-id', type=str, metavar='run_ID', required=True, help='Run ID of the batch file.')
	parser.add_argument('-s', type=str, metavar='subsample_num', required=True, help='Subsample number used in this test for threshold calculation.')
	# ### end of sam edit (4) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	args = parser.parse_args()

	# Generate ordered phenotypes and headder lists

	#list of the sample names, ordered by phenotype smallest to largest
	ordered_pheno=ordered_list(args.f, args.p)


### SAM EDIT ############################
##### create the dataframe with the order of accessions in orderof phenotype
	genotype_tracker_df = pandas.DataFrame(columns=['Accession_ID'])
	# sets initial column equal to the accession IDs of individuals

	#genotype_tracker_df.iloc[:,'Accession_ID'] = ordered_pheno
	genotype_tracker_df['Accession_ID'] = ordered_pheno

	#test print
	print("genotype_tracker_df initialising: ",  flush=True)
	print(genotype_tracker_df, flush=True)
########################################################

	#header of VCF including things like POS, ID and each sample e.g. 10013
	headder=read_vcf_head(args.v)
	print(f"header order: {headder}",flush=True)
	
	# open vcf file and output file
	vcf=open(args.v, 'r')

	# tempararily disabled
	# output_file=open(args.o, 'w+')
	# # Write the headder
	# #output_file.write('CHROM,POS,largest_theta,smallest_theta,absolute_theta,theta_range,largest_relative_theta,smallest_relative_theta,absolute_relative_theta,range_relative_theta\n')
	# #output_file.write('CHROM,POS,min_p,mean_p,log_mean_p,bigest_theta_p,pSNP4,pSNP5\n')
	# output_file.write('CHROM,POS,largest_theta,smallest_theta,absolute_theta,theta_range,largest_relative_theta,smallest_relative_theta,absolute_relative_theta,range_relative_theta,min_p,mean_p,log_mean_p,bigest_theta_p,pSNP4,pSNP5\n')
 
# ### start of sam edit (5) <<<<<<<<<<<<<<<<<<<< (disabled)
	# add in tracker csv file (Format: CHROM, POS, PVAL)
	# these will store and track the top 20 SNPs for each output e.g. SNP4 ...

	# T20_absolute_theta_output=open("output_files/"+args.id+"_T20_absolute_theta.csv","w")
	# # write in header CHROM,POS,PVAL
	# T20_absolute_theta_output.write("CHROM,POS,PVAL"+"\n")
	# T20_absolute_theta_output.close()

	# T20_pSNP4_output=open("output_files/"+args.id+"_T20_pSNP4.csv","w")
	# #write in header CHROM,POS,PVAL
	# T20_pSNP4_output.write("CHROM,POS,PVAL"+"\n")
	# T20_pSNP4_output.close()

	# T20_pSNP5_output=open("output_files/"+args.id+"_T20_pSNP5.csv","w")
	# #write in header CHROM,POS,PVAL
	# T20_pSNP5_output.write("CHROM,POS,PVAL"+"\n")
	# T20_pSNP5_output.close()
# ### end of sam edit (5) <<<<<<<<<<<<<<<<<<<<

	# calculate what I want
	for line in vcf:
		if '#' not in line:
			#line is line of VCF -> takes all information from all samples in the vcf (itll only consider samples in ordered phenotypes of course)
			# ordered_pheno is  orederd phenotypes and headder lists from the input phenotype and input phenotype file.
			# header is the header of the vcf
			x=field(line, ordered_pheno, headder) 
			if int(x.line[1]) % 10000 == 0:
				print(x.line[1])  # remove this print?

			# ensure that theres at least 15 + and 15 - states
			# disabled this if statement for now using 1==2
			if x.sense_check() and 1==2:
				largest_theta,smallest_theta,absolute_theta,theta_range,largest_relative_theta,smallest_relative_theta,absolute_relative_theta,range_relative_theta = x.calc_theta()
				#output_file.write(f'{CHROM},{POS},{largest_theta},{smallest_theta},{absolute_theta},{theta_range},{largest_relative_theta},{smallest_relative_theta},{absolute_relative_theta},{range_relative_theta}\n')
				CHROM,POS = x.line[0], x.line[1]
				min_p, mean_p, log_mean_p, bigest_theta_p, pSNP4, pSNP5 = x.calc_pvals()
				#output_file.write(f'{CHROM},{POS},{min_p},{mean_p},{log_mean_p},{bigest_theta_p},{pSNP4},{pSNP5}\n')
				output_file.write(f'{CHROM},{POS},{largest_theta},{smallest_theta},{absolute_theta},{theta_range},{largest_relative_theta},{smallest_relative_theta},{absolute_relative_theta},{range_relative_theta},{min_p},{mean_p},{log_mean_p},{bigest_theta_p},{pSNP4},{pSNP5}\n')
				
# ### start of sam edit (6) <<<<<<<<<<<<<<<<<<<<
				# open up the csv to import into dataframes
				dataFrame_absolute_theta = pandas.read_csv("output_files/"+args.id+"_T20_absolute_theta.csv")
				dataFrame_pSNP4 = pandas.read_csv("output_files/"+args.id+"_T20_pSNP4.csv")
				dataFrame_pSNP5 = pandas.read_csv("output_files/"+args.id+"_T20_pSNP5.csv")

				# check the T20 files to see how many lines they have
				n_abs_theta_lines= len(dataFrame_absolute_theta)
				n_pSNP4_lines = len(dataFrame_pSNP4)
				n_pSNP5_lines = len(dataFrame_pSNP5)
				# old code below
				#with open(r"output_files/"+args.id+"_T20_pSNP4.csv", 'r') as fp:
					#n_abs_theta_lines = len(fp.readlines())

				# ===============================================================
				# ABSOLUTE_THETA SECTION
				# =============================================================
				# if theres room for more SNPs (i.e. 20 limit not reached) 
				# then add the current absolute theta value
				# <20 since pandas doesnt count the header?
				if n_abs_theta_lines<20:
					# add in a line for the current line of data in the calculation
					new_row=pandas.Series({"CHROM":CHROM,"POS":POS,"PVAL":absolute_theta})
					dataFrame_absolute_theta =pandas.concat([dataFrame_absolute_theta, new_row.to_frame().T], ignore_index=True)

					# re-sort the csv so lowest pval is at the top and biggest pval is last row
					dataFrame_absolute_theta.sort_values(by=["PVAL"], axis=0, ascending=False,inplace=True, na_position='first')

				# if no room in list -> check if abs theta is bigger than current lowest abs theta in t20 list
				# by checking if its larger than last pval in file (should be the highest)
				else:

					# source highest pval from dataframe (should be last index and pval col)
					lowest_T20_abs_theta=dataFrame_absolute_theta.iloc[-1]["PVAL"]
					if float(lowest_T20_abs_theta)<float(absolute_theta):
						# replace smallest value in list (last item in the list)
						# first remove the smallest abs theta (last item in the dictionary)
						dataFrame_absolute_theta=dataFrame_absolute_theta.drop(dataFrame_absolute_theta.index[-1])
						
						# add in a line for the current line of data in the calculation
						new_row=pandas.Series({"CHROM":CHROM,"POS":POS,"PVAL":absolute_theta})
						dataFrame_absolute_theta =pandas.concat([dataFrame_absolute_theta, new_row.to_frame().T], ignore_index=True)

						# re-sort the csv so highest abs theta is at the top and lowest is at the bottom
						dataFrame_absolute_theta.sort_values(by=["PVAL"], axis=0, ascending=False,inplace=True, na_position='first')
				# ===============================================================
				# PSNP4 SECTION
				# =============================================================
				# if theres room for more SNPs (i.e. 20 limit not reached) 
				# then add the current pSNP4
				# <20 since pandas doesnt count the header?
				if n_pSNP4_lines<20:
					# add in a line for the current line of data in the calculation
					new_row=pandas.Series({"CHROM":CHROM,"POS":POS,"PVAL":pSNP4})

					dataFrame_pSNP4 =pandas.concat([dataFrame_pSNP4, new_row.to_frame().T], ignore_index=True)

					# re-sort the csv so lowest pval is at the top and biggest pval is last row
					dataFrame_pSNP4.sort_values(by=["PVAL"], axis=0, ascending=True,inplace=True, na_position='first')

				# if no room in list -> check if its smaller than the biggest value in list
				# by checking if its larger than last pval in file (should be the highest)
				else:

					# source highest pval from dataframe (should be last index and pval col)
					highest_pval=dataFrame_pSNP4.iloc[-1]["PVAL"]
					if float(highest_pval)>float(pSNP4):
						# replace biggest value in list (last item in the list)
						# first remove the biggest item (last item in the dictionary)
						dataFrame_pSNP4=dataFrame_pSNP4.drop(dataFrame_pSNP4.index[-1])
						
						# add in a line for the current line of data in the calculation
						new_row=pandas.Series({"CHROM":CHROM,"POS":POS,"PVAL":pSNP4})
						dataFrame_pSNP4 =pandas.concat([dataFrame_pSNP4, new_row.to_frame().T], ignore_index=True)

						# re-sort the csv so lowest pval is at the top and biggest pval is last row
						dataFrame_pSNP4.sort_values(by=["PVAL"], axis=0, ascending=True,inplace=True, na_position='first')
				# ===============================================================
				# PSNP5 SECTION
				# =============================================================

				# if theres room for more SNPs (i.e. 20 limit not reached) 
				# then add the current pSNP5
				# <20 since pandas doesnt count the header?
				if n_pSNP5_lines<20:
					# add in a line for the current line of data in the calculation
					new_row=pandas.Series({"CHROM":CHROM,"POS":POS,"PVAL":pSNP5})

					dataFrame_pSNP5 =pandas.concat([dataFrame_pSNP5, new_row.to_frame().T], ignore_index=True)

					# re-sort the csv so lowest pval is at the top and biggest pval is last row
					dataFrame_pSNP5.sort_values(by=["PVAL"], axis=0, ascending=True,inplace=True, na_position='first')

				# if no room in list -> check if its smaller than the biggest value in list
				# by checking if its larger than last pval in file (should be the highest)
				else:

					# source highest pval from dataframe (should be last index and pval col)
					highest_pval=dataFrame_pSNP5.iloc[-1]["PVAL"]
					if float(highest_pval)>float(pSNP5):
						# replace biggest value in list (last item in the list)
						# first remove the biggest item (last item in the dictionary)
						dataFrame_pSNP5=dataFrame_pSNP5.drop(dataFrame_pSNP5.index[-1])
						
						# add in a line for the current line of data in the calculation
						new_row=pandas.Series({"CHROM":CHROM,"POS":POS,"PVAL":pSNP5})
						dataFrame_pSNP5 =pandas.concat([dataFrame_pSNP5, new_row.to_frame().T], ignore_index=True)

						# re-sort the csv so lowest pval is at the top and biggest pval is last row
						dataFrame_pSNP5.sort_values(by=["PVAL"], axis=0, ascending=True,inplace=True, na_position='first')

				# write all dataframes back to csv 
				dataFrame_absolute_theta.to_csv("output_files/"+args.id+"_T20_absolute_theta.csv", index=False)
				dataFrame_pSNP4.to_csv("output_files/"+args.id+"_T20_pSNP4.csv", index=False)
				dataFrame_pSNP5.to_csv("output_files/"+args.id+"_T20_pSNP5.csv", index=False)

				# COMMENTS FOR AN IDEA OF CUSTOM SNPS OF INTEREST (NOT T_20)
					# Check for SNPs of interest (by CHROM and POS since rsID isnt tracked here or in GWAS script)
					# open file with positions of interest (POIs): POI_file=open("core_files/POI_SNPs.txt", 'r')
					# loop through it and store all 20 positions as a list: for line in ....
					# check if position is in POI_list : if so, write the full output line to a file : "POI_SNPs.csv"
				# ^^^^ this is not yet implemented, only T_20 tracking is implemented as of now

# ### end of sam edit (6) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> (disabled)
	#output_file.close()

### SAM EDIT ############################
	# output the genotype tracker to a csv file
	genotype_tracker_df.to_csv("core_files/genotype_tracker/"+args.id+"_genotypes.csv",header=True,index=False)

########################################################


	# Plot stuff
	#headder items to plot: largest_theta,smallest_theta,absolute_theta,theta_range,largest_relative_theta,smallest_relative_theta,absolute_relative_theta,range_relative_theta
	#R_plots(args.o, metric='largest_theta')
	#R_plots(args.o, metric='smallest_theta')
	
	#temp disabled
	#R_plots(args.o, metric='absolute_theta')
	
	#R_plots(args.o, metric='theta_range')
	#R_plots(args.o, metric='max_theta_plus')
	#R_plots(args.o, metric='max_theta_minus')
	#R_plots(args.o, metric='largest_relative_theta')
	#R_plots(args.o, metric='smallest_relative_theta')
	#R_plots(args.o, metric='absolute_relative_theta')
	#R_plots(args.o, metric='range_relative_theta')
	#R_plots(args.o, metric='min_p')
	#R_plots(args.o, metric='mean_p')
	#R_plots(args.o, metric='log_mean_p')
	#R_plots(args.o, metric='bigest_theta_p')

	#temp disabled
	# R_plots(args.o, metric='pSNP4')
	# R_plots(args.o, metric='pSNP5')
# end of file