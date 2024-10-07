# modified to accept arguments
args<-commandArgs(trailingOnly=TRUE)

print(args)

# read the data
# Pheno <- read.csv("Pheno.csv")
# Geno <- read.csv("Geno.csv")

# read the data (modified)
Pheno <- read.csv(args[1]) # phenotype file 
# e.g.core_files/subsampled_data/subsampled_phenotype_200_1366612.csv

Geno <- read.csv(args[2]) # genotype file 

#remove accession ID column
Geno <-Geno[-1]

# print(Geno)

# e.g. core_files/genotype_tracker/1366612_genotypes.csv
# Function to calculate the number of each microstate for all SNPs
# Input variables:
#   Geno: genotype dataset (row: individuals, col: SNPs)
# Output variables:
#   Nmpz: dataframe of number of microstates for each SNP

Function_Nmpz <- function(Geno){
  
  Nm <- data.frame(Nm = apply(Geno, 2, function(x){sum(x==-1, na.rm = TRUE)}))
  Np <- data.frame(Np = apply(Geno, 2, function(x){sum(x==1, na.rm = TRUE)})) 
  Nz <- data.frame(Nz = apply(Geno, 2, function(x){sum(x==0, na.rm = TRUE)}))
  Nmpz <- cbind(Nm,Np,Nz)
  
  return(Nmpz)
}

# Call Function_Nmpz
Nmpz <- Function_Nmpz(Geno)
# Function to order the individuals according to their phenotype and create the sorted datasets
# Input variables:
#   Geno: geno dataset (row: individuals, col: SNPs)
#   Pheno: pheno dataset (row: individuals, col: id, phenotype)
# Output variables:
#   List of sorted geno (row: individuals, col: SNPs) and pheno datasets (row: individuals, col: phenotype)

Fun_sortData <- function(Geno, Pheno){
  GenoSorted <- data.frame(matrix(NA, ncol = ncol(Geno), nrow = nrow(Pheno)))
  
  pheno.ordered <- sort(Pheno[,2], decreasing = FALSE, index.return = TRUE)
  PhenoSorted <- data.frame(IID = (Pheno[,1])[pheno.ordered$ix], trait = pheno.ordered$x)
  
  GenoSorted <- Geno[pheno.ordered$ix,]
  
  return(list(Geno = GenoSorted, Pheno = PhenoSorted))
}

# Call Fun_sortData
sortData <- Fun_sortData(Geno, Pheno)
# Function to compute the cumulative distributions of the ordered and random configurations
# Input variables:
#   Geno: sorted microstates dataset (row: individuals, col: SNPs)
#   Nmpz: dataframe of number of microstates for each SNP (rows)
# Output variables:
#   List of Wp, Wm, Wz, W0p, W0m, W0z: Cumulative distributions of +1,-1 respectively for ordered and random configurations

Fun_CumDistr <- function(Geno, Nmpz){
  library(tidyr) # tibble, pivot_wider
  library(dplyr)
  library(zoo) # na.locf
  
  CumDistrW <- as.data.frame(apply(Geno, 2, function(x){
    # Convert vector of microstates to data frame
    df <- tibble(position_j = seq_along(x), value = x)
    uniq.val <- as.character(unlist(unique(df[!is.na(df[,2]),2])))
    
    # Create cumulative sum for each unique value
    result <- df %>%
      group_by(value) %>%
      mutate(cumulative_sum = cumsum(value == value)) %>%
      ungroup()
    
    # Pivot the result to have one column per unique value
    result <- pivot_wider(result, names_from = value, values_from = cumulative_sum, values_fill = 0)
    
    # Replace zeros with NA values from the index onward
    result <- apply(result[,uniq.val], 2, function(x){
      first_non_zero_index <- which(x != 0)[1]
      x[first_non_zero_index:length(x)][x[first_non_zero_index:length(x)] == 0] <- NA
      return(x)
    })
    
    # Remove rows with the missing genotype
    indToRemove <- which(is.na(x))
    if (length(indToRemove)!=0){result <- result[-indToRemove,]}
    
    # Use na.locf() to fill NA values with the last non-NA value
    result <- as.data.frame(na.locf(result))
    Wm <- unlist(result[,which(uniq.val=="-1")])
    Wp <- unlist(result[,which(uniq.val=="1")])
    Wz <- unlist(result[,which(uniq.val=="0")])
    
    # If have less than three states, replace the null vectors with zeros
    if (is.null(Wm)) {Wm <- rep(0, times = nrow(df))}
    if (is.null(Wp)) {Wp <- rep(0, times = nrow(df))}
    if (is.null(Wz)) {Wz <- rep(0, times = nrow(df))}
    
    # Make the vectors having the same length by filling last positions with NA 
    if (length(Wm)!=nrow(df)) {Wm <- c(Wm, rep(NA, nrow(df)-length(Wm)))}
    if (length(Wp)!=nrow(df)) {Wp <- c(Wp, rep(NA, nrow(df)-length(Wp)))}
    if (length(Wz)!=nrow(df)) {Wz <- c(Wz, rep(NA, nrow(df)-length(Wz)))}
    
    return(list(Wp = Wp, Wm = Wm, Wz = Wz))
  }))
  
  Wp <- CumDistrW[,grep("Wp$", names(CumDistrW), value = TRUE)]
  Wm <- CumDistrW[,grep("Wm$", names(CumDistrW), value = TRUE)]
  Wz <- CumDistrW[,grep("Wz$", names(CumDistrW), value = TRUE)]
  
  CumDistrW0 <- apply(Nmpz, 1, function(Nmpz){
    Nmpz <- unlist(Nmpz)
    N <- sum(Nmpz)
    W0p <- c(1:N) * Nmpz[2] / N
    W0m <- c(1:N) * Nmpz[1] / N
    W0z <- c(1:N) * Nmpz[3] / N
    return(data.frame(W0p = W0p, W0m = W0m, W0z = W0z))
  })
  
  W0p <- as.data.frame(do.call(qpcR:::cbind.na, lapply(CumDistrW0, `[[`, "W0p")))
  W0m <- as.data.frame(do.call(qpcR:::cbind.na, lapply(CumDistrW0, `[[`, "W0m")))
  W0z <- as.data.frame(do.call(qpcR:::cbind.na, lapply(CumDistrW0, `[[`, "W0z")))
  
  return(list(Wp = Wp, Wm = Wm, Wz = Wz, W0p = W0p, W0m = W0m, W0z = W0z))
}


# Call Fun_CumDistr
CumDistr <- Fun_CumDistr(sortData$Geno,Nmpz)
# Function to calculate Theta-paths
# Input variables:
#   CumDistr: List consisting of the Cumulative Distributions Wp, Wm, Wz, W0p, W0m, W0z
# Output variables:
#   List of Thp = Wp - W0p,
#           Thm = Wm - W0m,
#           Thj = Wp - Wm (the phenotype-responding genetic path),
#           Th0j = W0p - W0m (the default genetic path),
#           DThj = Thj - Th0j

Fun_ThPaths <- function(CumDistr){
  
  N <- nrow(CumDistr$Wp)
  num_SNPs <- ncol(CumDistr$Wp)
  
  Thp <- data.frame(matrix(NA, nrow = N, ncol = num_SNPs))
  Thm <- data.frame(matrix(NA, nrow = N, ncol = num_SNPs))
  Thj <- data.frame(matrix(NA, nrow = N, ncol = num_SNPs))
  Th0j <- data.frame(matrix(NA, nrow = N, ncol = num_SNPs))
  DThj <- data.frame(matrix(NA, nrow = N, ncol = num_SNPs))
  
  Thp <- CumDistr$Wp - CumDistr$W0p
  Thm <- CumDistr$Wm -CumDistr$W0m
  Thj <- CumDistr$Wp - CumDistr$Wm
  Th0j <- CumDistr$W0p - CumDistr$W0m  
  DThj <- Thj - Th0j
  
  return(list(Thp = Thp, Thm = Thm, Thj = Thj, Th0j = Th0j, DThj = DThj))
}


# Call Fun_CumDistr
ThPaths <- Fun_ThPaths(CumDistr)
# Function to calculate the p-value
# Input variables:
#   Geno: sorted geno dataset (row: individuals, col: SNPs)
#   Nmpz: dataframe of number of microstates for each SNP (col1 = N_-, col2 = N_+, col3 = N_0)
#   ThPaths: List of Theta-paths (output of Fun_ThPaths)
# Output variables:
#   pSNP8: vector consists of pSNP8 values

Fun_pSNP8 <- function(Geno, Nmpz, ThPaths){
  pSNP8 <- c()
  
  num_SNPs <- ncol(Geno)
  DThj <- as.data.frame(ThPaths$DThj)
  
  phi <- apply(DThj, 2, function(x){max(x, na.rm = TRUE)})
  phi_tilde <- apply(DThj, 2, function(x){min(x, na.rm = TRUE)})
  
  for (i in 1:num_SNPs) {
    
    no.NA <- !is.na(Geno[,i])
    N <- sum(no.NA)
    
    if (sum(Nmpz[i,]>=3)==3) { #check if we have three states - compute the pSNP8 for three-states
      e <- exp( -(8*N*(phi[i]-phi_tilde[i])^2) / (N*(N-Nmpz[i,3])-(Nmpz[i,2]-Nmpz[i,1])^2) )
      
      pSNP8[i] <- ( (1/2)-(1/pi)*atan(sqrt(-phi[i]*phi_tilde[i])/((phi[i]-phi_tilde[i])*sqrt(2))) ) *e
      
    } #end if condition for three states
    if (sum(Nmpz[i,]>=3)==2) {  # length(uniq.states)==2
      uniq.states <- data.frame(table(Geno[no.NA,i]))
      ind_geno <- which(uniq.states$Freq>=3)
      if (length(intersect(c(-1,0),uniq.states[ind_geno,1]))==2){ #check if we have -1 and 0 states to use Thjm and Thjz
        e <- exp( -(2*N*(phi[i]-phi_tilde[i])^2) / (uniq.states$Freq[ind_geno[1]]*uniq.states$Freq[ind_geno[2]]) )
        
        pSNP8[i] <- ( (1/2)-(1/pi)*atan(sqrt(-2*phi[i]*phi_tilde[i])/(2*(phi[i]-phi_tilde[i]))) ) *e
      } #end if condition for two states: -1, 0
      
      if (length(intersect(c(+1,0),uniq.states[ind_geno,1]))==2){ #check if we have +1 and 0 states to use Thjp and Thjz
        e <- exp( -(2*N*(phi[i]-phi_tilde[i])^2) / (uniq.states$Freq[ind_geno[1]]*uniq.states$Freq[ind_geno[2]]) )
        
        pSNP8[i] <- ( (1/2)-(1/pi)*atan(sqrt(-2*phi[i]*phi_tilde[i])/(2*(phi[i]-phi_tilde[i]))) ) *e
      } #end if condition for two states: +1, 0
      
      else if (length(intersect(c(-1,+1),uniq.states[ind_geno,1]))==2){ #check if we have -1 and +1 states to use Thjm and Thjp
        e <- exp( -(N*(phi[i]-phi_tilde[i])^2) / (2*uniq.states$Freq[ind_geno[1]]*uniq.states$Freq[ind_geno[2]]) )
        
        pSNP8[i] <- ( (1/2)-(1/pi)*atan(sqrt(-2*phi[i]*phi_tilde[i])/(2*(phi[i]-phi_tilde[i]))) ) *e
      } #end if condition for two states: -1, +1
    } #end if condition for two states
    
    if (sum(Nmpz[i,]>=3)<2) { #check if we have only one state
      pSNP8[i] <- 1 # set this to one so after taking the logarithm we have zero (non-significant)
    }
    
  }
  
  return(pSNP8)
  
}

# Call Fun_SNP8
pSNP8 <- Fun_pSNP8(sortData$Geno, Nmpz, ThPaths)
# print the resutls
# -log10(pSNP8)
##   [1] 1.4345565 2.3454476 1.5697092 2.5239315 3.1935865 1.6160982 2.5742162
##   [8] 3.1256929 2.9821764 1.9326580 1.6793189 3.0262339 2.4061109 1.6486148
##  [15] 1.5497902 3.0382329 1.7020914 3.2935544 1.3761550 1.3354568 2.4547571
##  [22] 1.1228967 2.2266280 2.7440259 1.2148204 1.9734982 1.9475189 2.1434411
##  [29] 1.3919811 2.1936408 2.5317805 1.3188666 1.7082395 1.9450845 2.1775905
##  [36] 1.7300716 1.2846930 1.0549732 2.8971978 1.8222817 1.9476803 1.9127795
##  [43] 1.6736846 2.7095063 2.0195736 0.8211802 1.3038139 1.7330951 2.2028351
##  [50] 2.1987472 1.2562649 2.2974478 1.6003969 2.4449933 2.2070828 1.3465240
##  [57] 1.1036016 1.8353216 1.3929625 2.5440730 2.0228890 0.7589827 1.1150197
##  [64] 1.8376851 1.9618578 1.0709342 1.7790705 1.6563086 1.2526108 1.7698584
##  [71] 2.8232969 1.4887981 2.0604609 1.5271504 2.0585893 2.2422488 3.7996934
##  [78] 1.9326931 3.3241667 1.9402297 1.8718083 1.6985661 2.6009021 1.5239820
##  [85] 2.1025233 1.6055342 2.2033212 2.2920782 1.7581886 2.9413408 1.9618478
##  [92] 2.2291531 2.0062640 1.6339028 1.1921031 1.4439298 2.2549471 1.9962322
##  [99] 1.5363508 1.9746386

positions_list <-colnames(Geno)

positions_list_without_X<- gsub("X","",positions_list)

split_positions<-strsplit(positions_list_without_X, "\\.")

chromosomes<-sapply(split_positions,function(x) as.numeric(x[1]))

SNP_positions<-sapply(split_positions,function(x) as.numeric(x[2]))

# result_df = data.frame(CHROM=c(chromosomes),
#                        POS=c(SNP_positions),
#                        transformed_PSNP8=c(-log10(pSNP8))
#                        )

result_df = data.frame(CHROM=c(chromosomes),
                       POS=c(SNP_positions),
                       PSNP8=c(pSNP8))

write.csv(result_df,file=args[3],row.names=FALSE,quote=FALSE)