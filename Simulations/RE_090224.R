suppressPackageStartupMessages({
  library(tidyverse)
  library(SIS)
  library(mvmeta)
  library(foreach)
  library(doSNOW)
  library(parallel)
  library(iterators)
  library(snow)
  library(ggplot2)
  library(metafor)
  library(MASS)
  library(meta)
  })

# source("S:/Rotation/PW/META/Sim/031824/Functions_121422.R")
source("/rsrch3/scratch/biostatistics/zxu7/Rotation/PW/META/Sim/090224/Functions_090224.R")

# input
args <- commandArgs(trailingOnly=TRUE)
aID1 <- as.numeric(args[1])  # Num of cores to run
aID2 <- as.numeric(args[2])  # Num of Iterations
aID3 <- as.numeric(args[3])  # Num of Sample Size
aID4 <- as.numeric(args[4])  # Proportion & Combination of mediators of Split Sample size for meta
aID5 <- as.numeric(args[5])  # Random seed for alpha and beta
aID6 <- as.numeric(args[6])  # Correlation Structure
aID7 <- as.numeric(args[7])  # RE or FE

# aID1 <- 1
# aID2 <- 2
# aID3 <- 1500
# aID4 <- 1
# aID5 <- 1
# aID6 <- 3
# aID7 <- 1

aID7 <- 2
if(aID7 == 1){
  FixedEffect <- TRUE
  FixedEffectOut <- "FE"
}else if (aID7 == 2){
  FixedEffect <- FALSE
  FixedEffectOut <- "RE"
}
cat("Random Effect or Fixed Effect: ", FixedEffectOut, "\n")

Pooled <- FALSE # Whether perform pooled data analysis
cat("Whether to run pooled analysis: ", Pooled, "\n")

# SampleSize must be 3*P
# splitProp <- list(rep(1/5, 5), rep(1/8, 8), rep(1/10, 10), rep(1/16, 16), rep(1/20, 20))

splitProp <- list(rep(1/5, 5), rep(1/8, 8), rep(1/10, 10), rep(1/16, 16) , rep(1/20, 20))
Meds <- list(c(150, 150, 150, 1050))

Combin <- expand.grid(splitProp, Meds) %>% arrange(Var1)
cat("This job is running for: ", as.character(Combin[aID4, ]), "\n")
cat("Total Sample Size is: ", aID3, "\n")

TrueM <- unlist(Combin$Var2[aID4])[1]
p1 <- unlist(Combin$Var2[aID4])[2] 
p2 <- unlist(Combin$Var2[aID4])[3]
p3 <- unlist(Combin$Var2[aID4])[4]
p <- TrueM + p1 + p2 + p3

CorM <- aID6
cat("Correlation is: ", CorM, "\n")

if(CorM == 3){
  corMat <- diag(p)
  corMatDat <- matrix(rnorm((TrueM+p1)^2, 0, 0.1), nrow = (TrueM+p1), ncol = (TrueM+p1))
  corMatDat[upper.tri(corMatDat)] <- t(corMatDat)[upper.tri(corMatDat)]
  diag(corMatDat) <- 1
  corMatDat <- Matrix::nearPD(corMatDat, corr = TRUE, base.matrix = TRUE, maxit = 300)$mat
  corMat[1:(TrueM+p1), 1:(TrueM+p1)] <- corMatDat
  diag(corMat) <- 1
}else{
  corMat <- NA
}


# Start parallel computing
cl <- makeCluster(20, type = "SOCK", outfile = "clusterout.txt")
doSNOW::registerDoSNOW(cl)

TrueR2Sim <- 1000
TrueR2List <- foreach(i=1:TrueR2Sim, .combine='c', .inorder = F) %dopar% {
  GetTrueR2(i, res.sd = 2, dist = "chi")
}

TrueR2 <- mean(as.numeric(TrueR2List))
TrueR2SD <- sd(as.numeric(TrueR2List))

cat("True R2 based on", TrueR2Sim, "reps for this setting is: ", TrueR2, "\n")
cat("SD of True R2 based on", TrueR2Sim, "reps for this setting is: ", TrueR2SD, "\n")

# Stop parallel computing
stopCluster(cl)
  


cat("Start Simulations", "\n")

# Start parallel computing
cl <- makeCluster(aID1, type = "SOCK", outfile = "clusterout.txt")
doSNOW::registerDoSNOW(cl)

result <- foreach(rep = 1:aID2, .combine = "rbind", .inorder = F, .errorhandling = "pass") %dopar% {
  
    startTime2 <- Sys.time()
    # splitSize = aID3*rep(1, length(unlist(Combin$Var1[aID4])))
    splitSize = aID3*unlist(Combin$Var1[aID4])
    
    Study_Q = length(unlist(Combin$Var1[aID4]))
    
    nStudies <- length(splitSize)
    resultDF_S <- data.frame(Study = seq_len(nStudies), SampleSize = splitSize, EstR2_S_tem = NA, EstSD_S_tem = NA, TrueR2 = NA)
    
    for(study_it in 1:nStudies){
      
      ABDat <- genAB(seed=study_it*1000+rep, 
                     TrueM = TrueM, 
                     p1 = p1, 
                     p2 = p2, 
                     p3 = p3)
      alpha <- ABDat$alpha
      beta <- ABDat$beta
      
      set.seed(rep)
      dat <- genDat(N=splitSize[study_it], alpha=alpha, beta=beta, r = 3, res.sd = 2, dist = "chi", CorM = CorM, 
                    TrueM = TrueM, 
                    p1 = p1, 
                    p2 = p2, 
                    p3 = p3, 
                    corMat = corMat, TrueR2only=FALSE)
      # X[[study_it]] <- dat$x
      # Y[[study_it]] <- dat$Y
      # M[[study_it]] <- dat$M
      resultDF_S$TrueR2[study_it] <- dat$TrueR2
      
      res_S <- infer_r2(dat$x, dat$M, dat$Y, iter.max=3, penalty="MCP", tune="bic")
      resultDF_S$EstR2_S_tem[study_it] <- res_S$r2_est
      resultDF_S$EstSD_S_tem[study_it] <- res_S$std_asym
      
    }
    
    # Define a function to perform meta-analysis with various methods and calculate metrics
    
    methods <- c("FE", "DL", "HE","HSk", "SJ", "ML", "REML", "EB", "PM", "PMM")
    variance <- resultDF_S$EstSD_S_tem^2
    true_effect <- TrueR2
    
    # Perform analysis
    results <- perform_analysis(resultDF_S$EstR2_S_tem, variance, true_effect, methods)
    
    # Concatenate all results into a single vector with unique names
    all_results_vector <- c()
    for (method in names(results)) {
      # Prefix each result component with the method name for clarity
      method_results <- setNames(results[[method]], paste(method, names(results[[method]]), sep = "_"))
      all_results_vector <- c(all_results_vector, method_results)
    }
    
    TrueR2_S <- paste(round(resultDF_S$TrueR2, 2), collapse = "_")
    
    # Pooled Effects
    if(Pooled == TRUE){
      startTime1 <- Sys.time()
      X_O <- as.numeric(do.call("c", X)) 
      Y_O <- as.numeric(do.call("c", Y)) 
      M_O <- as.matrix(do.call("rbind", M))
      rm(X,Y,M)
      
      res_O <- infer_r2(X_O, M_O, Y_O,  iter.max=3, penalty="MCP", tune="bic")
      EstR2_O <- res_O$r2_est
      AsymSE_O <- res_O$std_asym
      Bias_O <- EstR2_O - TrueR2
      CP_O <- as.numeric(ifelse(abs(Bias_O) <= qnorm(0.975) * AsymSE_O, 1, 0))
      endTime1 <- Sys.time()
      TimeUsed1 <- difftime(endTime1, startTime1, units = "secs")
      
    }else{
      startTime1 <- Sys.time()
      EstR2_O <- -999
      AsymSE_O <- -999
      Bias_O <- -999
      CP_O <- -999
      endTime1 <- Sys.time()
      TimeUsed1 <- difftime(endTime1, startTime1, units = "secs")
      
    }

  
  return(c(TrueR2=TrueR2, Study_Q=Study_Q, EstR2_O=EstR2_O, AsymSE_O=AsymSE_O, Bias_O=Bias_O, CP_O=CP_O,
           TrueR2_S=TrueR2_S, CorM=CorM, REFE = FixedEffectOut,
           all_results_vector,
           SampleSize=aID3, Seed=aID5, 
           SplitProp=as.character(paste(round(unlist(Combin$Var1[aID4]),3), collapse = "_")), 
           Meds=as.character(paste(round(unlist(Combin$Var2[aID4]),3), collapse = "_"))))

}

# Stop parallel computing
stopCluster(cl)

result <- as.data.frame(result)
outputDir <- "/rsrch3/scratch/biostatistics/zxu7/Rotation/PW/META/Sim/090224/Result/Skew/" ### Modify
outRoot <- "PopMean_"
setwd(outputDir)
# save(result, file = paste0(outRoot, aID1, "_", aID2, "_", aID3, "_", aID4, "_", aID5, "_", aID6, "_", FixedEffectOut, ".rData"))
write.table(result, paste0(outRoot, aID1, "_", aID2, "_", aID3, "_", aID4, "_", aID5, "_", aID6, "_", FixedEffectOut, ".txt"), append=F, quote=F, row.names=F, col.names=T)




