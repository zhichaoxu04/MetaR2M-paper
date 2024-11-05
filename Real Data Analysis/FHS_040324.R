# .libPaths("/home/zxu7/R/ubuntu/4.2.0")
# ----- Load the packages
library(GMMAT)
library(SIS)
library(tidyverse)
library(HDMT)
library(tidyr)
cat("Load Packages Completed", "\n")

# ------ Load the data from local or HPC
# load("S:/Rotation/PW/META/RDA/FourStudies_011223.RData")
load("/rsrch3/scratch/biostatistics/zxu7/Rotation/PW/META/RDA/FourStudies_011223.RData")  # Change Rdata
load("/rsrch3/scratch/biostatistics/zxu7/Rotation/PW/META/RDA/FHS_Gene_TranscriptionID_annotation.RData")  # Change Rdata
source("/rsrch3/scratch/biostatistics/zxu7/Rotation/PW/META/RDA/040324/RDA_Infer_101023.R") # Change Function 
print("Load Data Completed")
outputDir <- "/rsrch3/scratch/biostatistics/zxu7/Rotation/PW/META/RDA/040324/Result/"

# input
args <- commandArgs(trailingOnly=TRUE)
aID1 <- as.numeric(args[1])   # DATASET
aID2 <- as.numeric(args[2])   # Random Seed

# aID1 <- 1
# aID2 <- 2023



seed <- aID2
iter.max <- 3
OutExpDat <- tidyr::expand_grid(data.frame(outcome = c("SBP_adj", "HDL"),
                                           exposure = c("Age", "Sex")),
                                method = c("iSIS"),
                                FDR = c(TRUE,FALSE),
                                DataNum = c(1:4))
outcome <- OutExpDat$outcome[aID1]
exposure <- OutExpDat$exposure[aID1]
FDR <- OutExpDat$FDR[aID1]
method <- OutExpDat$method[aID1]
DataNum <- OutExpDat$DataNum[aID1]

if(FDR == TRUE){
  FDR_out <- "FDR"
}else if (FDR == FALSE){
  FDR_out <- "NOFDR"
}

cat("Random seed= ", seed, "\n")
cat("Iter.max = ", iter.max, "\n")
cat("outcome is ", outcome, "\n")
cat("exposure is ", exposure, "\n")

if(DataNum == 1){
  data <- Microarray_OFF %>% 
    select(shareid, Age, Sex, Currsmk, Alcohol2, BMI, cohort, # X and Cov
           all_of(outcome), 
           starts_with("ID_")) %>%  # Med
    filter(complete.cases(.)) 
  
} else if(DataNum == 2){
  data <- Microarray_GEN %>% 
    select(shareid, Age, Sex, Currsmk, Alcohol2, BMI, cohort, # X and Cov
           all_of(outcome), 
           starts_with("ID_")) %>%  # Med
    filter(complete.cases(.)) 
  
} else if(DataNum == 3){
  data <- RNAseq_OFF %>% 
    select(shareid, ENSG00000223972.5:last_col()) %>% 
    select(where(~sum(.) != 0)) %>% 
    left_join(RNAseq_OFF %>% select(shareid, Age, Sex, Currsmk, Alcohol2, BMI, cohort, # X and Cov
                                     all_of(outcome)), by = "shareid") %>% 
    relocate(shareid, Age, Sex, Currsmk, Alcohol2, BMI, cohort, # X and Cov
             all_of(outcome)) %>% 
    filter(complete.cases(.)) 
} else if(DataNum == 4){
  data <- RNAseq_GEN %>% 
    select(shareid, ENSG00000223972.5:last_col()) %>% 
    select(where(~sum(.) != 0)) %>% 
    left_join(RNAseq_GEN %>% select(shareid, Age, Sex, Currsmk, Alcohol2, BMI, cohort, # X and Cov
                                    all_of(outcome)), by = "shareid") %>% 
    relocate(shareid, Age, Sex, Currsmk, Alcohol2, BMI, cohort, # X and Cov
             all_of(outcome)) %>% 
    filter(complete.cases(.)) 
  
} 

dataset <- c("Microarray_OFF", "Microarray_GEN", "RNAseq_OFF", "RNAseq_GEN")
cat("Dataset is ", dataset[DataNum], "\n")

# ----------------------------------- 
#         Analysis
# -----------------------------------


cat("N Row of data = ", dim(data)[1], "\n") 
cat("N Col of data = ", dim(data)[2], "\n") 
cat("Last colnames of dataset:", tail(colnames(data)), "\n")
cat("First 9 colnames of dataset:", head(colnames(data), 9), "\n")

# ----- 

n <- nrow(data)   # Sample Size
cat("Sample Size = ", n, "\n") 

set.seed(seed)
idx1 <- sample(1:n, ceiling(n/2), replace = FALSE)

col1 <- data %>% 
  slice(idx1) %>% 
  select(-cohort) %>%  
  select(where(~sum(.) != 0)) %>%
  colnames()
col2 <- data %>% 
  slice(-idx1) %>% 
  select(-cohort) %>%  
  select(where(~sum(.) != 0)) %>%
  colnames()
colFinal <- col1[col1 %in% col2]

data <- data %>% select(all_of(colFinal))

X <- as.numeric(data %>% select(all_of(exposure)) %>% pull())   # X is the exposure

if(DataNum %in% c(1,2)){
  M <- as.matrix(data %>% select(starts_with("ID_")))
}else if (DataNum %in% c(3,4)){
  M <- as.matrix(data %>% select(starts_with("E")))
}


cat("First 3 meditors name is:", colnames(M)[1:3], "\n")

d <- ncol(M)   
cat("number of M is:", d, "\n")


Covar <- data %>% select(Age:BMI, -all_of(exposure))  # Other 3 variables
cat("Covariates names are:", colnames(Covar), "\n")

Y <- as.numeric(data %>% select(all_of(outcome)) %>% pull()) # Y is the outcome

result <- R2(Y = Y, M = M, Covar = Covar, X = X, d = d, n = n, 
             iter.max = iter.max, nsis = NULL, first.half = FALSE, seed = seed,
             FDR = FDR, FDRCutoff = 0.2, method=method, idx1 = idx1)
result[["Geneset1"]] <- colnames(M)[result$select1]
result[["Geneset2"]] <- colnames(M)[result$select2]

result$output

setwd(outputDir)
save("result", file = paste0(outputDir, "FHS_", outcome, "_", exposure, "_", FDR_out, "_", dataset[DataNum], "_", method, ".RData"))


