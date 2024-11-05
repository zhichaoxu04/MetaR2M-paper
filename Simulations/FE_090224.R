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
source("/rsrch3/scratch/biostatistics/zxu7/Rotation/PW/META/Sim/090224/Functions_090324.R")
cat("Source Code finished \n")

# input
args <- commandArgs(trailingOnly=TRUE)
aID1 <- as.numeric(args[1])  # Num of cores to run
aID2 <- as.numeric(args[2])  # Num of Iterations
aID3 <- as.numeric(args[3])  # Num of Sample Size
aID4 <- as.numeric(args[4])  # Scenario
aID5 <- as.numeric(args[5])  # Random seed for alpha and beta
aID6 <- as.numeric(args[6])  # Correlation Structure

# aID1 <- 1
# aID2 <- 2
# aID3 <- 300
# aID4 <- 1
# aID5 <- 1
# aID6 <- 3

aID7 <- 1
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

splitProp <- list(c(1), 
                  c(1/2, 1/2), c(1/3, 2/3), 
                  rep(1/3, 3), c(1/4, 1/4, 1/2), 
                  rep(1/4, 4), 
                  rep(1/5, 5))
# Meds <- list(c(150,0,0,1350), c(150,0,150,1200), c(150,150,0,1200), c(150,150,150,1050), c(15, 0, 0, 1485), c(5, 0, 0, 1495))
Meds <- list(c(0, 50, 50, 1400))

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

ABDat <- genAB(seed=aID5, 
               TrueM = TrueM, 
               p1 = p1, 
               p2 = p2, 
               p3 = p3)
alpha <- ABDat$alpha
beta <- ABDat$beta

if(TrueM == 0){
  TrueR2 <- 0
}else{
  TrueR2 <- GetTrueR2(aID5, res.sd = 1, dist = "normal")
}

cat("Start Simulations", "\n")

# Start parallel computing
cl <- makeCluster(aID1, type = "SOCK", outfile = "clusterout.txt")
doSNOW::registerDoSNOW(cl)

result <- foreach(rep = 1:aID2, .combine = "rbind", .inorder = F, .errorhandling = "pass") %dopar% {
  set.seed(rep)
  dat <- genDat(N=aID3, alpha=alpha, beta=beta, r = 3, res.sd = 1, dist = "normal", CorM = CorM, 
                TrueM = TrueM, 
                p1 = p1, 
                p2 = p2, 
                p3 = p3, 
                corMat = corMat, TrueR2only=FALSE)
  X <- dat$x
  Y <- dat$Y
  M <- dat$M

  splitSize = aID3*unlist(Combin$Var1[aID4])
  Study_Q = length(unlist(Combin$Var1[aID4]))
  nStudies <- length(splitSize)
  resultDF_S <- data.frame(Study = seq_len(nStudies), SampleSize = splitSize, 
                           EstR2_S_tem = NA, EstSD_S_tem = NA, TrueR2 = TrueR2)
  
  for(study_it in 1:nStudies){
    group <- rep(c(1:nStudies), times = splitSize)
    X_tem <- dat$x[group == study_it]
    Y_tem <- dat$Y[group == study_it]
    M_tem <- dat$M[group == study_it, ]
    
    res_S <- infer_r2(X_tem, M_tem, Y_tem, iter.max=3, penalty="MCP", tune="bic")
    resultDF_S$EstR2_S_tem[study_it] <- res_S$r2_est
    resultDF_S$EstSD_S_tem[study_it] <- res_S$std_asym
    
  }
  
  analysis <- metafor::rma(yi=resultDF_S$EstR2_S_tem, vi=resultDF_S$EstSD_S_tem^2, method = "FE")
  coef <- as.numeric(analysis$beta)
  se <- as.numeric(analysis$se)
  bias <- coef - TrueR2
  cp <- as.numeric(abs(bias) <= qnorm(0.975) * se)
  
  TrueR2_S <- paste(round(resultDF_S$TrueR2, 2), collapse = "_")
  
  
  return(c(TrueR2=TrueR2, Study_Q=Study_Q, 
           TrueR2_S=TrueR2_S, CorM=CorM, REFE = FixedEffectOut,
           FE_Estimate=coef, FE_SE=se, FE_Bias=bias, FE_CP=cp,
           SampleSize=aID3, Seed=aID5, 
           SplitProp=as.character(paste(round(unlist(Combin$Var1[aID4]),3), collapse = "_")), 
           Meds=as.character(paste(round(unlist(Combin$Var2[aID4]),3), collapse = "_"))))
  
}

# Stop parallel computing
stopCluster(cl)

result <- as.data.frame(result)
result

outputDir <- "/rsrch3/scratch/biostatistics/zxu7/Rotation/PW/META/Sim/090224/Result/NoM/" ### Modify
outRoot <- "PopMean_"
setwd(outputDir)
# save(result, file = paste0(outRoot, aID1, "_", aID2, "_", aID3, "_", aID4, "_", aID5, "_", aID6, "_", FixedEffectOut, ".rData"))
write.table(result, paste0(outRoot, aID1, "_", aID2, "_", aID3, "_", aID4, "_", aID5, "_", aID6, "_", FixedEffectOut, ".txt"), append=F, quote=F, row.names=F, col.names=T)


