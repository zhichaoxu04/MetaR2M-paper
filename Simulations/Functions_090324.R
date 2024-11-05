library(tidyverse)
library(SIS)
library(mvmeta)
library(meta)
library(aod) # For wald test

genAB <- function(seed, TrueM = 50, p1 = 150, p2 = 100, p3 = 300){
  set.seed(seed) 
  
  alpha <- c(rnorm(TrueM, mean = 0, sd = 1.5), 
             rep(0, p1),
             rnorm(p2, mean = 0, sd = 1.5), 
             rep(0, p3))
  
  beta <- c(rnorm(TrueM, mean = 0, sd = 1.5), 
            rnorm(p1, mean = 0, sd = 1.5), 
            rep(0, p2 + p3))
  return(list(alpha=alpha, beta=beta))
}

CalTrueR2 <- function(alpha, beta, r=3, res.sd=1){
  vary <- (r + sum(alpha * beta))^2 + sum(beta^2) + res.sd^2
  Rsq.YM.true <- as.vector(((r + sum(alpha * beta))^2 - r^2/(1 + sum(alpha[beta!=0]^2)) + sum(beta[alpha!=0]^2))/vary)
  Rsq.YX.true <- as.vector((r + sum(alpha * beta))^2/vary)
  Rsq.YMX.true <- as.vector(((r + sum(alpha * beta))^2 + sum(beta[alpha!=0]^2))/vary)
  TrueR2 <- Rsq.YM.true + Rsq.YX.true - Rsq.YMX.true
  SOS.true <- TrueR2/Rsq.YX.true
  return(TrueR2=as.numeric(TrueR2))
}


genDat <- function(N, alpha, beta, r = 3, res.sd = 1, CorM = c(1,2,3), TrueM, p1 ,p2 ,p3, corMat, dist = "normal", TrueR2only=TRUE){
  
  if(!length(alpha) == length(beta)){
    stop("The length of Alpha and Beta is different!")
  }
  
  p <- TrueM+p1+p2+p3
  ck <- 0
  
  N_large <- 100000
  x_large <- rnorm(N_large) # Generate the same x
  M_large <- matrix(0, nrow = N_large, ncol = p)
  
  if (CorM == 1){
    for (i in 1:p){M_large[, i] <- alpha[i] * x_large + rnorm(N_large, 0, 1)}
    
  }else if(CorM == 2){
    for (i in 1:500){M_large[, i] <- alpha[i] * x_large + rnorm(N_large, 0, 1)}
    for (i in 501:1000){M_large[, i] <- 2*M_large[, i-500] + rnorm(N_large, 0, 1)}
    for (i in 1001:p){M_large[, i] <- -M_large[, 1] + rnorm(N_large, 0, 1)}
    
  }else if(CorM == 3){
    Residual_large <- MASS::mvrnorm(N_large, mu = rep(0, p), Sigma = corMat)
    for (i in 1:p){M_large[, i] <- alpha[i] * x_large + Residual_large[, i]}
    
  }
  
  if(dist == "normal"){
    Y_large <- r * x_large + as.vector(beta %*% t(M_large)) + rnorm(N_large, 0, res.sd)
  }else if(dist == "chi"){
    Y_large <- r * x_large + as.vector(beta %*% t(M_large)) + as.numeric(scale(rchisq(N_large, df = res.sd)))
  }
  
  
  olsYXM.large <- lm(Y_large ~ cbind(x_large, M_large[, 1:(TrueM+p1)]))
  VYXM_large <- var(olsYXM.large$residuals)
  olsYM.large <- lm(Y_large ~ M_large[, 1:(TrueM+p1)])
  VYM_large <- var(olsYM.large$residuals)
  olsYX.large <- lm(Y_large ~ x_large)
  VYX_large <- var(olsYX.large$residuals)
  VY_large <- var(Y_large)
  
  R45.true <- 1 - (VYX_large + VYM_large - VYXM_large)/VY_large
  Rsq.YM.true <- 1 - VYM_large/VY_large
  Rsq.YX.true <- 1 - VYX_large/VY_large
  Rsq.YMX.true <- 1 - VYXM_large/VY_large
  SOS.true <- R45.true/Rsq.YX.true
  
  if(TrueR2only == TRUE){
    x=NA
    M=NA
    Y=NA
  }else{
    x <- rnorm(N, 0, 1) # Generate the same x
    M <- matrix(0, nrow = N, ncol = p)
    
    if (CorM == 1){
      for (i in 1:p){M[, i] <- alpha[i] * x + rnorm(N, 0, 1)}
      
    }else if(CorM == 2){
      for (i in 1:500){M[, i] <- alpha[i] * x + rnorm(N, 0, 1)}
      for (i in 501:1000){M[, i] <- 2*M[, i-500] + rnorm(N, 0, 1)}
      for (i in 1001:p){M[, i] <- -M[, 1] + rnorm(N, 0, 1)}
      
    }else if(CorM == 3){
      Residual <- MASS::mvrnorm(N, mu = rep(0, p), Sigma = corMat)
      for (i in 1:p){M[, i] <- alpha[i] * x + Residual[, i]}
      
    }
    
    if(dist == "normal"){
      Y <- r * x + as.vector(beta %*% t(M)) + rnorm(N, 0, res.sd)
    }else if(dist == "chi"){
      Y <- r * x + as.vector(beta %*% t(M)) + as.numeric(scale(rchisq(N, df = res.sd)))
    }
    
    
    M <- scale(M, center = T, scale = T) # Standardized M
  }
  
  
  
  return(list(alpha=alpha, beta=beta, x=as.numeric(x), M=as.matrix(M), Y=as.numeric(Y), TrueR2=as.numeric(R45.true)))
  
}

GetTrueR2 <- function(seed, res.sd = 1, dist = "normal"){
  ABDat <- genAB(seed=seed, 
                 TrueM = TrueM, 
                 p1 = p1, 
                 p2 = p2, 
                 p3 = p3)
  alpha <- ABDat$alpha
  beta <- ABDat$beta
  
  dat <- genDat(N=1, alpha=alpha, beta=beta, r = 3, res.sd = res.sd, CorM = CorM, dist=dist,
                TrueM = TrueM, 
                p1 = p1, 
                p2 = p2, 
                p3 = p3, 
                corMat = corMat, TrueR2only=TRUE)
  return(dat$TrueR2)
}


select_variable <- function(W, Y, iter.max = 3, penalty = c("SCAD", "MCP", "lasso"), tune = c("bic", "ebic", "aic", "cv")) {
  m <- tryCatch(SIS::SIS(x = W, y = Y, 
                         family = "gaussian", tune = tune, seed = 1, 
                         penalty = penalty, 
                         iter.max = iter.max,
                         nsis = NULL),  
                error = function(c) NA)
    select_vars <- union(m$ix, ncol(W))

  if(!is.na(m)[1]){
    return(list(select_vars = as.numeric(select_vars), penalty=penalty,tune=tune, iter.max=iter.max))
  }else{
    return(list(select_vars = NA, penalty=penalty,tune=tune, iter.max=iter.max))
  }
  
}

crossfit <- function(X, M, Y, idx, iter.max=3, penalty="MCP", tune="bic") {
  d <- ncol(M)
  W <- cbind(M, X)
  
  select_result <- select_variable(W = W[idx, ], Y = Y[idx], iter.max=iter.max, penalty=penalty, tune=tune)
  select_vars <- select_result$select_vars

  if (!is.na(select_vars[1])){
    select_vars <- union(select_result$select_vars, (d+1))
    
    m_yx <- lm(Y[-idx] ~ X[-idx])
    err_yx <- m_yx$residuals
    
    m_yw <- lm(Y[-idx] ~ W[-idx, select_vars])
    err_yw <- m_yw$residuals

    m_yz <- lm(Y[-idx] ~ W[-idx, select_vars[select_vars != (d + 1)]])
    err_yz <- m_yz$residuals
    
    err_y <- Y[-idx] - mean(Y[-idx])
    
    v_yx <- mean(err_yx^2)
    v_yw <- mean(err_yw^2)
    v_yz <- mean(err_yz^2)
    v_y <- mean(err_y^2)
    
    err <- cbind(err_yx^2,err_yz^2,err_yw^2,err_y^2)
    A <- cov(err)
    
  } else {
    m_yx <- lm(Y[-idx] ~ X[-idx])
    err_yx <- m_yx$residuals
    
    m_yw <- lm(Y[-idx] ~ X[-idx]) # No M
    err_yw <- m_yw$residuals
    
    m_yz <- lm(Y[-idx] ~ 1) # W is X when no M selected
    err_yz <- m_yz$residuals
    
    err_y <- Y[-idx] - mean(Y[-idx])
    
    v_yx <- mean(err_yx^2)
    v_yw <- mean(err_yw^2)
    v_yz <- mean(err_yz^2)
    v_y <- mean(err_y^2)
    
    err <- cbind(err_yx^2,err_yz^2,err_yw^2,err_y^2)
    A <- cov(err)
  }
  
  return(list(v_yw = v_yw, v_yz = v_yz, v_yx = v_yx, v_y = v_y, A = A))
}

infer_r2 <- function(X, M, Y, iter.max=3, penalty="MCP", tune="bic") {
  n <- nrow(M)
  
  # indices of observations in a subsample.
  idx <- sample(1:n, ceiling(n / 2), replace = FALSE)
  
  m1 <- crossfit(X, M, Y, idx, iter.max, penalty, tune)
  m2 <- crossfit(X, M, Y, -idx, iter.max, penalty, tune)
  
  A <- 0.5 * (m1$A + m2$A)
  v_yw <- 0.5 * (m1$v_yw + m2$v_yw)
  v_yz <- 0.5 * (m1$v_yz + m2$v_yz)
  v_yx <- 0.5 * (m1$v_yx + m2$v_yx)
  v_y <- var(Y)
  
  a <- c(-1/v_y, -1/v_y, 1/v_y, (v_yx+v_yz-v_yw)/v_y^2)
  v <- t(a) %*% A %*% a
    
  r2_est <- 1.0 - (v_yx + v_yz - v_yw) / v_y
  
  SOS <- r2_est / (1-v_yx/v_y)
  
  list(r2_est = r2_est, v_y = v_y, 
       v = as.numeric(v), 
       std_asym = as.numeric(sqrt(v) / sqrt(n)), 
       SOS=SOS)
}



Sim <- function(N=1000, splitSize=c(400, 200, 100, 300), seed=14){
  
  ABDat <- genAB(seed=8, TrueM = 50, p1 = 100, p2 = 100, p3 = 300)
  alpha <- ABDat$alpha
  beta <- ABDat$beta
  
  set.seed(seed)
  dat <- genDat(N=N, alpha=alpha, beta=beta, r = 3, res.sd = 1)
  X <- dat$x
  Y <- dat$Y
  M <- dat$M
  TrueR2 <- dat$TrueR2
  
  res_O <- infer_r2(X, M, Y,  iter.max=3, penalty="MCP", tune="bic")
  EstR2_O <- res_O$r2_est
  EstSD_O <- res_O$std_asym
  Bias_O <- EstR2_O - TrueR2
  CP_O <- as.numeric(ifelse(abs(Bias_O) <= qnorm(0.975) * EstSD_O, 1, 0))
  
  if(sum(splitSize) != N){
    stop("The sum of the split size does NOT match the total sample size")
  }else{
    nStudies <- length(splitSize)
    resultDF_S <- data.frame(Study = seq_len(nStudies), SampleSize = splitSize, EstR2_S_tem = NA, EstSD_S_tem = NA)
    for(study_it in 1:nStudies){
      group <- rep(c(1:nStudies), times = splitSize)
      X_tem <- dat$x[group == study_it]
      Y_tem <- dat$Y[group == study_it]
      M_tem <- dat$M[group == study_it, ]
      
      res_S <- infer_r2(X_tem, M_tem, Y_tem, iter.max=3, penalty="MCP", tune="bic")
      resultDF_S$EstR2_S_tem[study_it] <- res_S$r2_est
      resultDF_S$EstSD_S_tem[study_it] <- res_S$std_asym
      
    }
    
    maFixed <- mvmeta(resultDF_S$EstR2_S_tem ~ 1, S=resultDF_S$EstSD_S_tem^2, method="fixed")
    EstR2_S <- as.numeric(coef(maFixed))
    Bias_S <- as.numeric(coef(maFixed)) - TrueR2
    EstSD_S <- as.numeric(sqrt(maFixed$vcov))
    CP_S <- as.numeric(ifelse(abs(Bias_S) <= qnorm(0.975) * EstSD_S, 1, 0))
  }
  
  return(c(TrueR2, EstR2_O, EstSD_O, Bias_O, CP_O, EstR2_S, EstSD_S, Bias_S, CP_S))
  
}

perform_analysis <- function(effect_size, variance, true_effect, methods) {
  results <- list()
  
  for (method in methods) {
    if (method %in% c("FE", "DL", "HE","HSk", "SJ", "ML", "REML", "EB", "PM", "GENQ", "PMM", "GENQM")) {
      analysis <- metafor::rma(effect_size, variance, method = method)
      coef <- as.numeric(analysis$beta)
      se <- as.numeric(analysis$se)
    } else {
      next
    }
    
    bias <- coef - true_effect
    cp <- as.numeric(abs(bias) <= qnorm(0.975) * se)
    
    # Store results as a named vector
    results_vector <- c(Estimate = coef, SE = se, Bias = bias, CP = cp)
    names(results_vector) <- c("Estimate", "SE", "Bias", "CP")
    
    # Add the named vector to the results list
    results[[method]] <- results_vector
  }
  
  return(results)
}




