# ----- Load the packages
library(GMMAT)
library(SIS)
library(tidyverse)
library(HDMT)
library(tidyr)
cat("Load Packages Completed", "\n")

outputDir <- "/rsrch3/scratch/biostatistics/zxu7/Rotation/PW/META/RDA/072724/Result/"

# input
args <- commandArgs(trailingOnly=TRUE)
aID1 <- as.numeric(args[1])   
aID2 <- as.numeric(args[2])   
aID3 <- as.numeric(args[3])
aID4 <- as.numeric(args[4])
aID5 <- as.numeric(args[5])

# ------ Load the data from local or HPC

if(aID2 == 1){
  outcome <- "SBP_adj"
  exposure <- "Age"
}else(aID2 == 2){
  outcome <- "HDL"
  exposure <- "Sex"
}

cat("Outcome is: ", outcome, "\n Exposure is: ", exposure, "\n")

if(aID1 == 1){
  Out1 <- "FHS"
  # load("/Users/xu/Library/CloudStorage/OneDrive-InsideMDAnderson/MDACC/PW/Meta_analysis/RDA/rData/FHS_FourStudies_072824.rData")
  load("/rsrch3/scratch/biostatistics/zxu7/Rotation/PW/META/RDA/FHS_FourStudies_072824.rData")
  
  if(aID3 == 1){
    DatName <- "Offsrping-Micro"
    data <- Offspring_Micro %>% 
      dplyr::select(shareid, Age, Sex, Currsmk, Alcohol2, BMI, Study, PC1:PC10,
                    all_of(outcome), 
                    all_of(exposure), 
                    starts_with("ID_")) %>%  # Med
      filter(complete.cases(.)) 
    
  } else if(aID3 == 2){
    DatName <- "GEN3-Micro"
    data <- GEN3_Micro %>% 
      dplyr::select(shareid, Age, Sex, Currsmk, Alcohol2, BMI, Study, PC1:PC10,
                    all_of(outcome), 
                    all_of(exposure), 
                    starts_with("ID_")) %>%  # Med
      filter(complete.cases(.)) 
    
  } else if(aID3 == 3){
    DatName <- "Offsrping-RNAseq"
    data <- Offspring_RNAseq %>% 
      dplyr::select(shareid, Age, Sex, Currsmk, Alcohol2, BMI, Study, PC1:PC10,
                    all_of(outcome), 
                    all_of(exposure), 
                    ENSG00000223972.5:last_col()) %>%  # Med
      filter(complete.cases(.)) 
    
  } else if(aID3 == 4){
    DatName <- "GEN3-RNAseq"
    data <- GEN3_RNAseq %>% 
      dplyr::select(shareid, Age, Sex, Currsmk, Alcohol2, BMI, Study, PC1:PC10,
                    all_of(outcome), 
                    all_of(exposure), 
                    ENSG00000223972.5:last_col()) %>%  # Med
      filter(complete.cases(.)) 
    
  } 
}else if(aID1 == 2){
  Out1 <- "MESA"
  load("/rsrch3/scratch/biostatistics/zxu7/Rotation/PW/META/RDA/MESA_FourStudies_072824.rData")
  
  if(aID3 == 1){
    DatName <- "MESA1"
    data <- MESA_1 %>% 
      dplyr::rename(shareid = sidno, 
                    Age = age1c, Sex = gender1, Currsmk = cursmk1, Alcohol2 = curalc1, BMI = bmi1c, 
                    HDL = hdl1, SBP_adj = sbp1c_adj) %>% 
      dplyr::select(shareid, Age, Sex, Currsmk, Alcohol2, BMI, Study, PC1:PC10,
                    all_of(outcome), 
                    all_of(exposure), 
                    starts_with("E")) %>%  # Med
      filter(complete.cases(.)) 
    
  } else if(aID3 == 2){
    DatName <- "MESA2"
    data <- MESA_2 %>% 
      dplyr::rename(shareid = sidno, 
                    Age = age1c, Sex = gender1, Currsmk = cursmk1, Alcohol2 = curalc1, BMI = bmi1c, 
                    HDL = hdl1, SBP_adj = sbp1c_adj) %>% 
      dplyr::select(shareid, Age, Sex, Currsmk, Alcohol2, BMI, Study, PC1:PC10,
                    all_of(outcome), 
                    all_of(exposure), 
                    starts_with("E")) %>%  # Med
      filter(complete.cases(.)) 
    
  } else if(aID3 == 3){
    DatName <- "MESA3"
    data <- MESA_3 %>% 
      dplyr::rename(shareid = sidno, 
                    Age = age1c, Sex = gender1, Currsmk = cursmk1, Alcohol2 = curalc1, BMI = bmi1c, 
                    HDL = hdl1, SBP_adj = sbp1c_adj) %>% 
      dplyr::select(shareid, Age, Sex, Currsmk, Alcohol2, BMI, Study, PC1:PC10,
                    all_of(outcome), 
                    all_of(exposure), 
                    starts_with("E")) %>%  # Med
      filter(complete.cases(.)) 
    
  } else if(aID3 == 4){
    DatName <- "MESA4"
    data <- MESA_4 %>% 
      dplyr::rename(shareid = sidno, 
                    Age = age1c, Sex = gender1, Currsmk = cursmk1, Alcohol2 = curalc1, BMI = bmi1c, 
                    HDL = hdl1, SBP_adj = sbp1c_adj) %>% 
      dplyr::select(shareid, Age, Sex, Currsmk, Alcohol2, BMI, Study, PC1:PC10,
                    all_of(outcome), 
                    all_of(exposure), 
                    starts_with("E")) %>%  # Med
      filter(complete.cases(.)) 
    
  } else if(aID3 == 5){
    DatName <- "MESA-ALL"
    data <- MESA_all %>% 
      dplyr::rename(shareid = sidno, 
                    Age = age1c, Sex = gender1, Currsmk = cursmk1, Alcohol2 = curalc1, BMI = bmi1c, 
                    HDL = hdl1, SBP_adj = sbp1c_adj) %>% 
      dplyr::select(shareid, Age, Sex, Currsmk, Alcohol2, BMI, Study, PC1:PC10,
                    all_of(outcome), 
                    all_of(exposure), 
                    starts_with("E")) %>%  # Med
      filter(complete.cases(.)) 
    
  } 
}else if(aID1 == 3){
  Out1 <- "MEGA"
  load("/rsrch3/scratch/biostatistics/zxu7/Rotation/PW/META/RDA/MEGA_Dat_072824.rData")
  DatName <- "MEGA"
  data <- MEGA_dat_PCs %>% 
    dplyr::select(shareid, Age, Sex, Currsmk, Alcohol2, BMI, Study, PC1:PC10,
                  all_of(outcome), 
                  all_of(exposure), 
                  TTLL10:last_col()) %>%  # Med
    filter(complete.cases(.)) 
}

cat("Load data completed \n")

source("/rsrch3/scratch/biostatistics/zxu7/Rotation/PW/META/RDA/072724/RDA_Infer_101023.R") 
cat("Load function completed \n")


cat("Dataset is ", DatName, "\n")
cat("# of coloums data = ", dim(data)[2], "\n") 


# ----- 
n <- nrow(data)   # Sample Size
cat("Sample Size = ", n, "\n") 

if(aID4 == 1){
  Out4 <- "Point"
  data <- data
}else if(aID4 == 2){
  set.seed(aID5)
  Out4 <- "Bootstrap"
  data <- data %>% 
    slice_sample(n = n, replace = TRUE) 
}

set.seed(aID5)
idx1 <- sample(1:n, ceiling(n/2), replace = FALSE)

col1 <- data %>% 
  slice(idx1) %>% 
  select(-Study) %>%  
  select(where(~sum(.) != 0)) %>%
  colnames()
col2 <- data %>% 
  slice(-idx1) %>% 
  select(-Study) %>%  
  select(where(~sum(.) != 0)) %>%
  colnames()
colFinal <- col1[col1 %in% col2]

data <- data %>% select(all_of(colFinal))

X <- as.numeric(data %>% select(all_of(exposure)) %>% pull()) 

if(DatName %in% c("Offsrping-Micro", "GEN3-Micro")){
  M <- as.matrix(data %>% dplyr::select(starts_with("ID_")))
}else if (DatName %in% c("MEGA")){
  M <- as.matrix(data %>% dplyr::select(TTLL10:last_col()))
}else{
  M <- as.matrix(data %>% dplyr::select(starts_with("E")))
}

d <- ncol(M)   
cat("# of M is:", d, "\n")

Covar <- data %>% select(Age:PC10, -all_of(exposure))  # Other 3 variables
cat("Covariates names are:", colnames(Covar), "\n")

Y <- as.numeric(data %>% select(all_of(outcome)) %>% pull()) # Y is the outcome

result <- R2(Y = Y, M = M, Covar = Covar, X = X, d = d, n = n, 
             iter.max = 3, nsis = NULL, first.half = FALSE, seed = 2024,
             FDR = FALSE, FDRCutoff = 0.2, method="iSIS", idx1 = idx1)
result[["Geneset1"]] <- colnames(M)[result$select1]
result[["Geneset2"]] <- colnames(M)[result$select2]

result$output

setwd(outputDir)
save("result", file = paste0(outputDir, Out1,"_",outcome, "_", exposure, "_", DatName, "_", Out4, "_", aID5, ".RData"))

cat("Job Done! \n")
























