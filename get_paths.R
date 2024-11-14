

library(data.table)
library(parallel)
library(foreach)
library(doParallel)

# read the data

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
        Wm = cumsum(x == -1))
      # Wz = cumsum(x == 0),
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
  # Parameters for retry
  max_attempts <- 5  # Maximum number of retry attempts
  retry_delay <- 2   # Delay in seconds between attempts
  attempt <- 1       # Start with the first attempt
  success <- FALSE   # Flag to indicate successful cluster creation

  # Retry loop
  while (attempt <= max_attempts && !success) {
    try({
    options(mc.cores = 3)
    cl <- makeCluster(3)
    clusterExport(cl,"CalcCumDistrW")
    success<-TRUE
    message("Cluster created successfully on attempt ", attempt)
     }, silent = TRUE)

    if (!success) {
    message("Attempt ", attempt, " failed. Retrying in ", retry_delay, " seconds...")
    Sys.sleep(retry_delay)  # Wait before retrying
    attempt <- attempt + 1  # Increment the attempt count
    }
  }
  
  # Check if the cluster was successfully created
  if (!success) {
    stop("Failed to create cluster after ", max_attempts, " attempts.")
  }
    # give Geno as a list due to data overhead issue causing the workers to idle
  CumDistrW <- Geno[, unlist(parLapply(cl, as.list(Geno) , function(x){
    list(Wp = cumsum(x == 1),
        Wm = cumsum(x == -1))
      # Wz = cumsum(x == 0),
  }), recursive = FALSE)]

  # Rename columns
  print("Renaming columns . . .")
  stopCluster(cl)
  # new_colnames <- unlist(lapply(names(Geno), function(col) paste0(col, ".", c("Wp", "Wz", "Wm"))))
  new_colnames <- unlist(lapply(names(Geno), function(col) paste0(col, ".", c("Wp","Wm"))))
  setnames(CumDistrW, new_colnames)
  
  # Subset the data table by selecting columns whose names end with the desired suffix
  Wp <- CumDistrW[, grep(paste0("Wp", "$"), colnames(CumDistrW), value = TRUE), with = FALSE]
  Wm <- CumDistrW[, grep(paste0("Wm", "$"), colnames(CumDistrW), value = TRUE), with = FALSE]

    # may not need this
  # Wz <- CumDistrW[, grep(paste0("Wz", "$"), colnames(CumDistrW), value = TRUE), with = FALSE]

  print("Calling function CumDistr W0 function . . .")
  CumDistrW0 <- apply(Nmpz, 1, function(Nmpz){
    Nmpz <- unlist(Nmpz)
      # total number of plus minus and zero states (should be accession number)
    N <- sum(Nmpz)
    
    W0p <- c(1:N) * Nmpz[2] / N
    W0m <- c(1:N) * Nmpz[1] / N

    # might not need this
    # W0z <- c(1:N) * Nmpz[3] / N
    # save(W0p,W0m,W0z, file="temp_store6.RData")
    #return(data.frame(W0p = W0p, W0m = W0m, W0z = W0z))
      # works
    # return(data.table(W0p = W0p, W0m = W0m, W0z = W0z))
    return(data.table(W0p = W0p, W0m = W0m))
  })
  
  
  # W0p <- as.data.frame(do.call(qpcR:::cbind.na, lapply(CumDistrW0, `[[`, "W0p")))
  # W0m <- as.data.frame(do.call(qpcR:::cbind.na, lapply(CumDistrW0, `[[`, "W0m")))
  # W0z <- as.data.frame(do.call(qpcR:::cbind.na, lapply(CumDistrW0, `[[`, "W0z")))
  
  W0p <- as.data.table(do.call(qpcR:::cbind.na, lapply(CumDistrW0, `[[`, "W0p")))
  W0m <- as.data.table(do.call(qpcR:::cbind.na, lapply(CumDistrW0, `[[`, "W0m")))
  # W0z <- as.data.table(do.call(qpcR:::cbind.na, lapply(CumDistrW0, `[[`, "W0z")))
  
  #return(list(Wp = Wp, Wm = Wm, Wz = Wz, W0p = W0p, W0m = W0m, W0z = W0z))
  return(list(Wp = Wp, Wm = Wm, W0p = W0p, W0m = W0m))

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
  
  # works but slightly more memory?
  # Thp <- data.table(matrix(NA, nrow = N, ncol = num_SNPs))
  # Thm <- data.table(matrix(NA, nrow = N, ncol = num_SNPs))
  Thj <- data.table(matrix(NA, nrow = N, ncol = num_SNPs))
  Th0j <- data.table(matrix(NA, nrow = N, ncol = num_SNPs))
  DThj <- data.table(matrix(NA, nrow = N, ncol = num_SNPs))
  
    #not used
  # Thp <- CumDistr$Wp - CumDistr$W0p
  #   # not used
  # Thm <- CumDistr$Wm -CumDistr$W0m
    # used
  Thj <- CumDistr$Wp - CumDistr$Wm
    # Used
  Th0j <- CumDistr$W0p - CumDistr$W0m  
    # used
  DThj <- Thj - Th0j
  
  # return(list(Thp = Thp, Thm = Thm, Thj = Thj, Th0j = Th0j, DThj = DThj))
  return(list(Thj = Thj, Th0j = Th0j, DThj = DThj))
}


# Call Fun_CumDistr
print("Calling function ThPaths . . .")
ThPaths <- Fun_ThPaths(CumDistr)

#write.csv(result_df,file=args[3],row.names=FALSE,quote=FALSE)
  # works

delta_df <-ThPaths$DThj 
split_string<-strsplit(args[3], "_")[[1]]
id_number<- split_string[7]
cleaned_id_number <- gsub("\\.csv$", "", id_number)
colnames(delta_df) <- paste0(colnames(delta_df), cleaned_id_number)
# write.csv(delta_df,file="Delta_Theta_Paths.csv",row.names=FALSE,quote=FALSE)
write.csv(delta_df,file=args[3],row.names=FALSE,quote=FALSE)

# end of script