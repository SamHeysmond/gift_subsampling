## Packages
library(data.table)
library(parallel)
library(foreach)
library(doParallel)

### INPUT FROM WIP CODE####

args<-commandArgs(trailingOnly=TRUE)

print(args)

print("Reading phenotype data . . . ")
Pheno <- fread(args[1]) # phenotype file 
# e.g.core_files/subsampled_data/subsampled_phenotype_200_1494720.csv

print("Reading genotype data . . .")
Geno = fread(args[2]) # genotype file 
# e.g. 1494720_genotypes.csv
#############

# read the data
# Pheno <- fread("Pheno.csv")
# Geno <- fread("Geno.csv")

# Geno<-Geno[,1:50]

# Function to calculate the number of each microstate for all SNPs
# Input variables:
#   Geno: genotype dataset (row: individuals, col: SNPs)
# Output variables:
#   Nmpz: dataframe of number of microstates for each SNP

Function_Nmpz <- function(Geno){
  
  Nm <- apply(Geno, 2, function(x){sum(x==-1, na.rm = TRUE)})
  Np <- apply(Geno, 2, function(x){sum(x==1, na.rm = TRUE)})
  Nz <- apply(Geno, 2, function(x){sum(x==0, na.rm = TRUE)})
  Nmpz <- cbind(Nm,Np,Nz)
  
  return(Nmpz)
}

# Call Function_Nmpz
print("Calling function Nmpz . . .")
Nmpz <- Function_Nmpz(Geno)

# Function to compute the cumulative distributions of the ordered and random configurations
# Input variables:
#   Geno: sorted microstates dataset (row: individuals, col: SNPs)
#   Nmpz: dataframe of number of microstates for each SNP (rows)
# Output variables:
#   List of Wp, Wm, Wz, W0p, W0m, W0z: Cumulative distributions of +1,-1 respectively for ordered and random configurations

# try defining globally
CalcCumDistrW <-function(x){
  list(Wp = cumsum(x == 1),
        Wz = cumsum(x == 0),
        Wm = cumsum(x == -1))
}


Fun_CumDistr <- function(Geno, Nmpz){
  library(tidyr) # tibble, pivot_wider
  library(dplyr)
  library(zoo) # na.locf

# Geno, 2 means margin of 2 which means for each column i.e. each position(SNP)
  # CalcCumDistrW <-function(x){
  #   list(Wp = cumsum(x == 1),
  #        Wz = cumsum(x == 0),
  #        Wm = cumsum(x == -1))
  # }
  
  # Use lapply to apply the function to all columns and flatten the result in one step
  print("Calling function CumDistr W function . . .")
    # this works on small data set but slow on bit data set
  #CumDistrW <- Geno[, unlist(lapply(.SD, CalcCumDistrW), recursive = FALSE)]
    # multi core version

    # this works (sorta? but do this without mcl apply)
#  CumDistrW <- Geno[, unlist(mclapply(.SD, CalcCumDistrW,mc.cores=6), recursive = FALSE)]

  options(mc.cores = 4)
  cl <- makeCluster(4)
  clusterEvalQ(cl, Sys.info())
  # Check what each worker is doing and log outputs
  clusterEvalQ(cl, {
  log_file <- file("worker_log.txt", "a")
  writeLines(paste("Starting job on worker", Sys.getpid()), log_file)
   writeLines(paste("sys info", Sys.getpid()), log_file)
  close(log_file)
  })
  clusterExport(cl,"CalcCumDistrW")
    # give Geno as a list due to data overhead issue causing the workers to idle
  CumDistrW <- Geno[, unlist(parLapply(cl, as.list(Geno) , function(x){
    list(Wp = cumsum(x == 1),
        Wz = cumsum(x == 0),
        Wm = cumsum(x == -1))
  }), recursive = FALSE)]

  # Rename columns
  print("Renaming columns . . .")
  stopCluster(cl)
  new_colnames <- unlist(lapply(names(Geno), function(col) paste0(col, ".", c("Wp", "Wz", "Wm"))))
  setnames(CumDistrW, new_colnames)
  
  # Subset the data table by selecting columns whose names end with the desired suffix
  Wp <- CumDistrW[, grep(paste0("Wp", "$"), colnames(CumDistrW), value = TRUE), with = FALSE]
  Wm <- CumDistrW[, grep(paste0("Wm", "$"), colnames(CumDistrW), value = TRUE), with = FALSE]
  Wz <- CumDistrW[, grep(paste0("Wz", "$"), colnames(CumDistrW), value = TRUE), with = FALSE]

  print("Calling function CumDistr W0 function . . .")
  CumDistrW0 <- apply(Nmpz, 1, function(Nmpz){
    Nmpz <- unlist(Nmpz)
      # total number of plus minus and zero states (should be accession number)
    N <- sum(Nmpz)
    
    W0p <- c(1:N) * Nmpz[2] / N
    W0m <- c(1:N) * Nmpz[1] / N
    W0z <- c(1:N) * Nmpz[3] / N
    # save(W0p,W0m,W0z, file="temp_store6.RData")
    return(data.frame(W0p = W0p, W0m = W0m, W0z = W0z))
  })
  
  
  # W0p <- as.data.frame(do.call(qpcR:::cbind.na, lapply(CumDistrW0, `[[`, "W0p")))
  # W0m <- as.data.frame(do.call(qpcR:::cbind.na, lapply(CumDistrW0, `[[`, "W0m")))
  # W0z <- as.data.frame(do.call(qpcR:::cbind.na, lapply(CumDistrW0, `[[`, "W0z")))
  
  W0p <- as.data.table(do.call(qpcR:::cbind.na, lapply(CumDistrW0, `[[`, "W0p")))
  W0m <- as.data.table(do.call(qpcR:::cbind.na, lapply(CumDistrW0, `[[`, "W0m")))
  W0z <- as.data.table(do.call(qpcR:::cbind.na, lapply(CumDistrW0, `[[`, "W0z")))
  
  return(list(Wp = Wp, Wm = Wm, Wz = Wz, W0p = W0p, W0m = W0m, W0z = W0z))
  # i instead want to return a list of data tables for each of the accessions (margin 1)
}

# Call Fun_CumDistr
# CumDistr <- Fun_CumDistr(sortData$Geno,Nmpz)
print("Calling function CumDistr . . .")
CumDistr <- Fun_CumDistr(Geno,Nmpz)
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
print("Calling function ThPaths . . .")
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
    
    # no.NA <- !is.na(Geno[,i])
    # N <- sum(no.NA)
    N <- nrow(Geno)
    if (sum(Nmpz[i,]>=3)==3) { #check if we have three states - compute the pSNP8 for three-states
      e <- exp( -(8*N*(phi[i]-phi_tilde[i])^2) / (N*(N-Nmpz[i,3])-(Nmpz[i,2]-Nmpz[i,1])^2) )
      
      pSNP8[i] <- ( (1/2)-(1/pi)*atan(sqrt(-phi[i]*phi_tilde[i])/((phi[i]-phi_tilde[i])*sqrt(2))) ) *e
      
    } #end if condition for three states
    if (sum(Nmpz[i,]>=3)==2) {  # length(uniq.states)==2
      # states_table<-Geno[, .N, by=eval(noquote(c(colnames(Geno)[i])))]
      # states_table<- data.table(Var1=names(table(Geno$noquote(colnames(Geno)[i]))),
      #                           Freq= as.integer(table(Geno$noquote(colnames(Geno)[i]))))
      # 
      uniq.states <- data.frame(table(Geno[,..i]))
      
      # uniq.states <- data.frame(table(Geno[,i]))
      

      
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
# pSNP8 <- Fun_pSNP8(sortData$Geno, Nmpz, ThPaths)
print("Calling function pSNP8 . . .")
pSNP8 <- Fun_pSNP8(Geno, Nmpz, ThPaths)
# print the resutls
# -log10(pSNP8)
##   [1] 1.4345565 2.3454476 1.5697092 2.5239315 3.1935865 1.6160982 2.5742162
##   [8] 3.1256929 2.9821764 1.9326580 1.6793189 3.0262339 2.4061109 1.6486148

print("Tidying up . . .")

positions_list <-colnames(Geno)

positions_list_without_X<- gsub("X","",positions_list)

split_positions<-strsplit(positions_list_without_X, ":")

chromosomes<-sapply(split_positions,function(x) as.numeric(x[1]))

SNP_positions<-sapply(split_positions,function(x) as.numeric(x[2]))

# print("positions_list")
# print(positions_list)

# print("positons_list_without_X")
# print(positions_list_without_X)
# print("split_positions")
# print(split_positions)

# print("chromosomes")
# print(chromosomes)

# print("SNP_positions")
# print(SNP_positions)

# result_df = data.frame(CHROM=c(chromosomes),
#                        POS=c(SNP_positions),
#                        transformed_PSNP8=c(-log10(pSNP8))
#                        )

result_df = data.frame(CHROM=c(chromosomes),
                       POS=c(SNP_positions),
                       PSNP8=c(pSNP8))

print("Exporting data . . .")
write.csv(result_df,file=args[3],row.names=FALSE,quote=FALSE)

                  