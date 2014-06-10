
#################################
### This file contains all    ###
### functions for simulations ###
#################################

library(compiler)
enableJIT(3)

#this function is used for fast at risk set calculation
cumsumRev <- function(a) rev(cumsum(rev(a)))

#this function checks if the working directory is ivsurv and changes it to simulations
setwd2sim <- function() {
  wd <- getwd();if (substr(wd, nchar(wd) - 5,nchar(wd)) == "ivsurv"){setwd(paste(getwd(), "/simulations", sep = ""))};getwd()
}

loadPackages <- function() {
  packages <- c("survival", "ggplot2", "reshape2", "foreach", 
                "Matrix", "SimCorMultRes", "grid", "BB")
  print("Loading...")
  for (i in packages) {library(i, character.only=T); print(i)}
  packages
}

runParallel <- function(ncores = NULL){
  switch(Sys.info()[["sysname"]], 
         "Windows" = {
           library(doSNOW)
           num.cores <- ifelse(is.null(ncores), getDoParWorkers(), ncores)
           cl <- makeCluster(max(1, num.cores - 1)); registerDoSNOW(cl)
         },
         "Linux" = {
           library(doMC); registerDoMC()
           if (length(grep("bigmem", Sys.info()[["nodename"]])) > 0) {
             if(is.null(ncores)) {
               options(cores = 15)
             } else {options(cores = ncores)}
           } else { options(cores = max(getDoParWorkers() - 1, 1))}
           cl <- NULL
         })
  cl
}


norta <- function(N = N, cor.matrix = cor.matrix, conf.corr.X = 0, instrument.strength = 0.8, beta0 = 0, beta1 = 1) {
  # generate correlated random variables with the IV assumptions using the NORTA method
  require(Matrix)
  require(SimCorMultRes)
  if (!is.numeric(N) | N < 1) 
    stop("'N' must be greater than or equal to one")
  N <- as.integer(N)
  ans <- rsmvnorm(R = N, cor.matrix = cor.matrix)
  probs <- pnorm(ans)
  U <- ans[,1] #confounder
  Z <- rnorm(N)
  X1 <- instrument.strength * Z + sqrt(1 - instrument.strength^2) * rnorm(N)
  X <- conf.corr.X * U + sqrt(1 - conf.corr.X^2) * X1
  ans <- cbind(Z, X, qexp(probs[,2], rate = exp(1 * (beta0 + beta1 * X))), U)
  colnames(ans) <- c("Z", "X", "Y", "U")
  data.frame(ans)
}

genIVData <- function(N = N, Z2XCoef, U2XCoef, U2YCoef, beta0 = 0, beta1 = 1,
                      survival.distribution = c("exponential", "normal"),
                      confounding.function = c("linear", "exponential", "square"),
                      break2sls = FALSE, break.method = c("collider", "error"), error.amount = 0.01) {
  if (!is.numeric(N) | N < 1) 
    stop("'N' must be greater than or equal to one")
  surv.dist <- match.arg(survival.distribution)
  confounding.function <- match.arg(confounding.function)
  break.method <- match.arg(break.method)
  
  
  
  N <- as.integer(N)
  
  
  if (break2sls) {
    if (break.method == "collider") {
      uy.confounder <- rnorm(N)
      U <- rnorm(N, sd = 1) + U2YCoef * uy.confounder
      err <- rnorm(N) + U2YCoef * uy.confounder
    } else if (break.method == "error") {
      U <- rnorm(N, sd = 1)
      err <- rnorm(N)
    }
  } else {
    U <- rnorm(N, sd = 1)
    err <- rnorm(N)
  }
  
  
  #U <- U / sd(U)
  Z <- rnorm(N, sd = 1)
  #Z <- Z / sd(Z)
  if (break2sls & break.method == "error") {
    X <- switch(confounding.function,
                linear = (Z2XCoef * Z) + (U2XCoef * U) + error.amount * err + rnorm(N, sd = 1),
                exponential = (Z2XCoef * Z) + (U2XCoef * exp(U)) + error.amount * err + rnorm(N, sd = 1),
                square = (Z2XCoef * Z) + (U2XCoef * U^2) + error.amount * err + rnorm(N, sd = 1))
  } else {
    X <- switch(confounding.function,
                linear = (Z2XCoef * Z) + (U2XCoef * U) + rnorm(N, sd = 1),
                exponential = (Z2XCoef * Z) + (U2XCoef * exp(U)) + rnorm(N, sd = 1),
                square = (Z2XCoef * Z) + (U2XCoef * U^2) + rnorm(N, sd = 1))
  }
  #X <- X / sd(X)
  
  if (!break2sls) {
    Y <- switch(surv.dist,
                exponential = rexp(N, rate = exp(-(beta0 + beta1 * X + U2YCoef * U))),
                normal = exp(beta0 + beta1 * X + U2YCoef * U + err))
  } else {
    if (break.method == "error") {
      Y <- exp(beta0 + beta1 * X + U2YCoef * U + err)
    } else {
      Y <- exp(beta0 + beta1 * X + err)
    }
  }

  ans <- cbind(Z, X, Y, U)
  colnames(ans) <- c("Z", "X", "Y", "U")
  data.frame(ans)
}

genIVDataOld <- function(N = N, Z2XCoef, U2XCoef, U2YCoef, beta0 = 0, beta1 = 1,
                      survival.distribution = c("exponential", "normal")) {
  if (!is.numeric(N) | N < 1) 
    stop("'N' must be greater than or equal to one")
  surv.dist <- match.arg(survival.distribution)
  N <- as.integer(N)
  U <- rnorm(N, sd = 1)
  #U <- U / sd(U)
  Z <- rnorm(N, sd = 1)
  #Z <- Z / sd(Z)
  X <- (Z2XCoef * Z) + (U2XCoef * U) + rnorm(N, sd = 1)
  #X <- X / sd(X)
  Y <- switch(surv.dist,
              exponential = rexp(N, rate = exp(-(beta0 + beta1 * X + U2YCoef * U))),
              normal = exp(beta0 + beta1 * X + U2YCoef * U + rnorm(N)))
  
  ans <- cbind(Z, X, Y, U)
  colnames(ans) <- c("Z", "X", "Y", "U")
  data.frame(ans)
}


genMultivarIVData <- function(N = N, Z2XCoef, U2XCoef, U2YCoef, beta, num.confounded,
                              survival.distribution = c("exponential", "normal"),
                              intercept = F, break2sls = FALSE) {
  if (!is.numeric(N) | N < 1) 
    stop("'N' must be greater than or equal to one")
  surv.dist <- match.arg(survival.distribution)
  num.vars <- length(beta)
  N <- as.integer(N)
  U <- matrix(rnorm(N * num.confounded), ncol = num.confounded)
  Z <- matrix(rnorm(N * num.confounded), ncol = num.confounded)

  if (intercept) {
    no.int.beta <- beta[-1]
    Z.append <- U.append <- array(0, dim = c(N,(num.vars - 1)))
    Z.append[, 1:num.confounded] <- Z
    U.append[, 1:num.confounded] <- U
    rand.mat <- matrix(rnorm(N * (num.vars - 1)), ncol = (num.vars - 1))
    X <- Z.append * Z2XCoef + U.append * U2XCoef + rand.mat
    #X2Y <- t(t(X) * no.int.beta)
    beta <- no.int.beta
    colnames(X) <- paste("X", 1:length(no.int.beta), sep = "")
  } else {
    Z.append <- U.append <- array(0, dim = c(N, num.vars))
    Z.append[, 1:num.confounded] <- Z
    U.append[, 1:num.confounded] <- U
    rand.mat <- matrix(rnorm(N * num.vars), ncol = num.vars)
    X <- Z.append * Z2XCoef + U.append * U2XCoef + rand.mat
    #X2Y <- t(t(X) * beta)
    colnames(X) <- paste("X", 1:length(beta), sep = "")
  }

  Y <- switch(surv.dist,
              exponential = rexp(N, rate = exp(-(X %*% beta + rowSums(U)))),
              normal = exp(X %*% beta + rowSums(U) + rnorm(N, sd = 0.2)))
  
  ret <- list(Z = Z, X = X, Y = Y, U = U)
  ret
}

roughResidSDEstAFTIV <- function(dat) {
  # returns rough estimate of sd of residuals
  # for smoothed rank estimator for AFT-IV
  lm.fit <- lm(X ~ Z, data = dat)
  dat$Xhat <- lm.fit$fitted.values
  lm.fit.2 <- lm(t ~ Xhat, data = dat)
  sd(resid(lm.fit.2))
}

roughResidSDEstAFT <- function(dat) {
  # returns rough estimate of sd of residuals
  # for smoothed rank estimator for AFT
  if (all(dat$t >= 0)) {dat$t <- log(dat$t)}
  lm.fit <- lm(t ~ X, data = dat)
  sd(resid(lm.fit))
}


simIVSurvivalData <- function(sample.size, conf.corr.X = 0.0, conf.corr.Y, instrument.strength, 
                              lambda, beta0, beta1, verbose = F, norta = F,
                              survival.distribution = c("exponential", "normal"),
                              confounding.function = c("linear", "exponential", "square"),
                              break2sls = FALSE, break.method = c("collider", "error"),
                              error.amount = 0.01) {
  #conf.corr.X == confounder correlation with X
  #conf.corr.Y == confounder correlation with Y
  
  surv.dist <- match.arg(survival.distribution)
  confounding.function <- match.arg(confounding.function)
  break.method <- match.arg(break.method)
  
  if (norta){
    #correlation matrix for U, Y
    sigma <- cbind(c(1,conf.corr.Y), c(conf.corr.Y, 1))
    #generate IV random variables
    vars <- norta(sample.size, sigma, conf.corr.X = conf.corr.X, instrument.strength = instrument.strength, beta0, beta1)
  } else {
    #generate IV random variables
    vars <- genIVData(sample.size, Z2XCoef = instrument.strength, U2XCoef = conf.corr.X, U2YCoef = conf.corr.Y, 
                      beta0, beta1, confounding.function = confounding.function,
                      survival.distribution = surv.dist, break2sls = break2sls, break.method = break.method,
                      error.amount = error.amount)
  }
  
  #if specified, print sample correlation between Z, X, U, and Y
  if (verbose == T) {print (cor(vars))}
  
  X <- vars$X #covariate
  Z <- vars$Z #instrument
  #failure variable
  Fail.time <- vars$Y
  #generate independent censoring data
  Cen.time <- rexp(sample.size, rate = lambda)
  #failure indicator
  delta <- 1 * (Cen.time >= Fail.time)
  # you can also use X<- pmin(Fail.time, Cen.time)
  t <- apply(data.frame(Fail.time, Cen.time), 1, min)
  #return variables in data.frame
  data.simu <- data.frame(t, delta, X, Z)
  #store correlations as attribute
  attr(data.simu, "cor") <- cor(vars)
  data.simu
}

simIVMultivarSurvivalData <- function(sample.size, conf.corr.X = 0.0, conf.corr.Y = 0, instrument.strength, 
                                      lambda, beta, survival.distribution = c("exponential", "normal"),
                                      num.confounded, intercept = F, break2sls = FALSE) {
  #conf.corr.X == confounder correlation with X
  #conf.corr.Y == confounder correlation with Y
  
  surv.dist <- match.arg(survival.distribution)
  
  #generate IV random variables
  vars <- genMultivarIVData(N = sample.size, Z2XCoef = instrument.strength, U2XCoef = conf.corr.X, 
                            U2YCoef = conf.corr.Y, beta, num.confounded = num.confounded,
                            survival.distribution = surv.dist, intercept = intercept)
  
  X <- vars$X #covariate
  Z <- vars$Z #instrument
  # failure variable
  Fail.time <- vars$Y
  # generate independent censoring data
  Cen.time <- rexp(sample.size, rate = lambda)
  # failure indicator
  delta <- 1 * (Cen.time >= Fail.time)
  # you can also use X<- pmin(Fail.time, Cen.time)
  t <- apply(data.frame(Fail.time, Cen.time), 1, min)
  log.t <- log(t)
  if (intercept) {
    log.t <- log.t + beta[1]
    t <- exp(log.t)
  }
  # return variables in data.frame
  survival <- data.frame(t, delta, log.t)
  
  # return a list of the survival aspects of 
  # the data and the covariate matrix
  ret <- list(survival = survival, X = X, Z = Z)
  dat <- data.frame(Z, X[,1:ncol(Z)], t, vars$U)
  names(dat) <- c( "Z", "X","Y", "U")
  attr(ret, "cor") <- cor(dat)
  class(ret) <- "survival.data"
  ret
}


SimIVDataCompareEstimators <- function(type, n.sims, sample.size, conf.corr.X = 0.0, conf.corr.Y = 0.0, instrument.strength,
                                       lambda, beta0, beta1, seed = NULL, norta=F, 
                                       survival.distribution = c("exponential", "normal"), break2sls = FALSE,
                                       break.method = c("collider", "error"), error.amount = 0.01){
  # This function simulates ('n.sims'-times) survival data with a confounding variable U and an instrument Z
  # and estimates beta using the regular AFT estimating equation and also using the IV estimating equation
  # proposed by Professor Yu. It stores the results in vectors and returns a list containing these vectors.
  
  #set seed if supplied by user
  if (!is.null(seed)) {set.seed(seed)}
  
  #check to make sure user specified allowed estimating equations
  types <- c("AFT", "AFT-IV", "AFT-2SLS", "AFT-IPCW", "AFT-2SLS-xhat")
  funcs <- c("vEvalAFTScore", "vEvalAFTivScore", "vEvalAFT2SLSScore", "vEvalAFTivIPCWScore", "vEvalAFT2SLSxhatScore")
  for (i in length(type)) {if (!is.element(type[i], types)) {stop("'type' must only contain 'AFT', 'AFT-IV',' AFT-2SLS' or 'AFT-IPCW'")}}
  
  beta.store <- array(0, dim = c(n.sims, length(type)))
  l <- 0
  pct.censored <- NULL
  num.errors <- tot.cors <- 0
  pb <- txtProgressBar(min = 0, max = n.sims, style = 3)
  while (l < n.sims){
    possibleError <- tryCatch({
      #####  GENERATE DATA  #########
      Data.simu <- simIVSurvivalData(sample.size, conf.corr.X = conf.corr.X, conf.corr.Y = conf.corr.Y, 
                                     instrument.strength = instrument.strength, lambda = lambda, 
                                     beta0 = beta0, beta1=beta1, norta = F,
                                     survival.distribution = survival.distribution, 
                                     break2sls = break2sls, break.method = break.method,
                                     error.amount = error.amount)

      beta.tmp <- array(0, dim = length(type))
      
      #solve each estimating equation for beta1
      for (e in 1:length(type)){
        #return correct estimating equation function
        est.eqn <- match.fun(funcs[[match(type[e], types)]])
        
        #solve for beta using bisection method
        beta.tmp[e] <- uniroot(est.eqn, interval = c(beta1 - 3, beta1 + 3), tol = 0.001, "data.simu" = Data.simu)$root
        
      }
      
    }, error=function(e) e)
    if(inherits(possibleError, "error")){
      num.errors <- num.errors + 1
      next
    } else {
      l <- l + 1
      #store results on successful attempt
      beta.store[l, ] <- beta.tmp
      setTxtProgressBar(pb, l)
      pct.censored[l] <- mean(1 - Data.simu$delta)
      #add correlation matrix so we can take mean
      tot.cors <- tot.cors + attr(Data.simu, "cor")
    }
  }
  avg.cors <- tot.cors / n.sims
  close(pb)
  print (paste("Number of errors:", num.errors))
  res <- list()
  for (i in 1:ncol(beta.store)){res[[i]] <- beta.store[, i]}
  names(res) <- type
  attr(res, "truth") <- beta1
  attr(res, "pct.censored") <- mean(pct.censored)
  attr(res, "avg.cor") <- avg.cors
  class(res) <- "AFTsim"
  res
}



SimIVDataCompareEstimatorsMultivar <- function(type, n.sims, sample.size, conf.corr.X = 0.0, conf.corr.Y = 0.0, instrument.strength,
                                               lambda, beta, seed = NULL, norta=F, intercept = F,
                                               survival.distribution = c("exponential", "normal"), break2sls = FALSE){
  # This function simulates ('n.sims'-times) survival data with a confounding variable U and an instrument Z
  # and estimates beta using the regular AFT estimating equation and also using the IV estimating equation
  # proposed by Professor Yu. It stores the results in vectors and returns a list containing these vectors.
  
  #set seed if supplied by user
  if (!is.null(seed)) {set.seed(seed)}
  
  #check to make sure user specified allowed estimating equations
  types <- c("AFT", "AFT-IV", "AFT-2SLS", "AFT-IPCW")
  funcs <- c("AFTScorePre", "AFTivScorePre", "AFT2SLSScorePre", "AFTivIPCWScorePre")
  funcs.sm <- c("AFTScoreSmoothPre", "AFTivScoreSmoothPre", "AFT2SLSScoreSmoothPre")
  for (i in length(type)) {if (!is.element(type[i], types)) {stop("'type' must only contain 'AFT', 'AFT-IV',' AFT-2SLS' or 'AFT-IPCW'")}}
  
  beta.store <- array(0, dim = c(n.sims, length(type)))
  l <- 0
  pct.censored <- NULL
  num.errors <- tot.cors <- 0
  pb <- txtProgressBar(min = 0, max = n.sims, style = 3)
  while (l < n.sims){
    possibleError <- tryCatch({
      #####  GENERATE DATA  #########
      Data.simu <- simIVMultivarSurvivalData(sample.size, conf.corr.X = conf.corr.X, conf.corr.Y = conf.corr.Y, 
                                             instrument.strength = instrument.strength, lambda = lambda, 
                                             beta = beta, intercept = intercept,
                                             survival.distribution = survival.distribution,
                                             num.confounded = 1, break2sls = break2sls)
      if (intercept) {
        Data.simu$X <- cbind(rep(1, nrow(Data.simu$X)), Data.simu$X)
        colnames(Data.simu$X)[1] <- "Intercept"
      }
      beta.tmp <- array(0, dim = length(type))
      
      #solve each estimating equation for beta1
      for (e in 1:length(type)){
        #return correct estimating equation function
        est.eqn <- match.fun(funcs[[match(type[e], types)]])
        est.eqn.sm <- match.fun(funcs.sm[[match(type[e], types)]])
        attr(est.eqn, "name") <- funcs[[match(type[e], types)]]
        attr(est.eqn.sm, "name") <- funcs.sm[[match(type[e], types)]]
        dfsane.tol <- 0.000001
        if (type[e] == "AFT-IPCW") {dfsane.tol <- 8}
        ssf <- 1e10
        
        ct <- 0
        if (e == 1) {
          init.par <- NULL
        } else {
          init.par <- est$par
          init.par <- NULL
        }

        #solve for beta using deriv-free spectral method
        est <- repFitAFT(tol = dfsane.tol, data = Data.simu, est.eqn = est.eqn, est.eqn.sm = est.eqn.sm,
                         instrument.names = "prop_endo", confounded.x.names = paste("X", 1:1, sep=""), 
                         fit.method = "dfsane", init.par = init.par, final.fit = F,
                         method = c(2), control = list(tol = dfsane.tol, maxit = 1000, trace = F, M = c(500)), quiet = T)
        
        conf.x.idx <- match(paste("X", 1:1, sep=""), names(est$par))
        beta.tmp[e] <- est$par[conf.x.idx]
        
      }
      
    }, error=function(e) e)
    if(inherits(possibleError, "error")){
      num.errors <- num.errors + 1
      print (possibleError)
      next
    } else {
      l <- l + 1
      #store results on successful attempt
      beta.store[l, ] <- beta.tmp
      setTxtProgressBar(pb, l)
      pct.censored[l] <- mean(1 - Data.simu$survival$delta)
      #add correlation matrix so we can take mean
      tot.cors <- tot.cors + attr(Data.simu, "cor")
    }
  }
  avg.cors <- tot.cors / n.sims
  close(pb)
  print (paste("Number of errors:", num.errors))
  res <- list()
  for (i in 1:ncol(beta.store)){res[[i]] <- beta.store[, i]}
  names(res) <- type
  attr(res, "truth") <- beta[1]
  attr(res, "pct.censored") <- mean(pct.censored)
  attr(res, "avg.cor") <- avg.cors
  class(res) <- "AFTsim"
  res
}


summary.AFTsim <- function(res) {
  # summary function for the output of the 
  # 'SimIVDataCompareEstimators' function
  stopifnot(class(res) == "AFTsim")
  results <- data.frame(array(0, dim = c(length(res), 10)))
  n.data <- length(res[[1]])
  results[,1] <- names(res)
  colnames(results) <- c("Estimator", "Mean", 
                         "LCI", "UCI", "sd", "Q2.5", "Q5", "Med", "Q95", "Q97.5")
  for (i in 1:length(res)){
    results[i, 2] <- mean(res[[i]])
    results[i, 5] <- sd <- sd(res[[i]])
    results[i, 3:4] <- c(mean(res[[i]]) - 1.96 * (sd / sqrt(n.data)), 
                         mean(res[[i]]) + 1.96 * (sd / sqrt(n.data)))
    results[i, 6:10] <- quantile(res[[i]], probs = c(0.025, 0.05, 0.5, 0.95, 0.975))
  }
  results[,2:10] <- round(results[,2:10], 4)
  print(paste("Average Pct Censored:", 100 * attr(res, "pct.censored")))
  attr(results, "cor.U.Y") <- attr(res, "avg.cor")[3,4]
  attr(results, "cor.U.X") <- attr(res, "avg.cor")[2,4]
  attr(results, "cor.Z.X") <- attr(res, "avg.cor")[1,2]
  print(sprintf("Cor(U,Y): %g  ||  Cor(U,X): %g  ||  Cor(Z,X): %g", 
          round(attr(res, "avg.cor")[3,4], 4), round(attr(res, "avg.cor")[2,4], 4), 
                round(attr(res, "avg.cor")[1,2], 4)))
  results
}


createGrid <- function(U2X.range, U2Y.range, Z2X.range, n.data, sample.size, lambda.range){
  stopifnot(is.numeric(U2X.range) & is.numeric(U2Y.range) & is.numeric(Z2X.range) 
            & is.numeric(lambda.range))
  if (!all(n.data %% 1 == 0) | !all(sample.size %% 1 == 0)) {stop("n.data and sample.size must be integers")}
  grid <- expand.grid(U2X.range, U2Y.range, Z2X.range, n.data, sample.size, lambda.range)
  colnames(grid) <- c("Conf.Corr.X", "Conf.Corr.Y", "Instrument.Strength", "n.data", "sample.size", "lambda")
  attr(grid, "grid") <- "sim.grid"
  grid
}

simulateGrid <- function(est.eqns, grid, beta, seed = NULL, 
                         survival.distribution = c("exponential", "normal"),
                         confounding.function = c("linear", "exponential", "square"),
                         break2sls = FALSE, break.method = c("collider", "error"),
                         error.amount = 0.01) {
  if (is.null(attr(grid, "grid"))) {stop("Grid must be created with createGrid function")}
  
  results <- array(0, dim = c(length(est.eqns), nrow(grid), (ncol(grid) + 12)))
  dimnames(results)[[1]] <- est.eqns
  dimnames(results)[[3]] <- c(colnames(grid), "Cor.U.Y", "Cor.U.X", "Cor.Z.X", "Mean", 
                              "LCI", "UCI", "sd", "Q2.5", "Q5", "Med", "Q95", "Q97.5")
  
  #simulate situation when there IS a confounder present for
  #varying levels of correlation between U and X and U and Y
  if (length(beta) == 1) {
    raw.results <- foreach(i = 1:nrow(grid), .packages = packages) %dopar% {
      print(sprintf("Simulation %g / %g", i, nrow(grid)))
      res <- SimIVDataCompareEstimators(type = est.eqns, n.sims = grid$n.data[i], 
                                        sample.size = grid$sample.size[i], 
                                        conf.corr.X = grid$Conf.Corr.X[i], 
                                        conf.corr.Y = grid$Conf.Corr.Y[i], 
                                        instrument.strength = grid$Instrument.Strength[i], 
                                        lambda = grid$lambda[i], beta0 = 0, beta1 = beta, seed = seed,
                                        survival.distribution = survival.distribution, break2sls = break2sls,
                                        break.method = break.method, error.amount = error.amount)
      for (j in 1:ncol(grid)) {attr(res, colnames(grid)[j]) <- grid[i, j]}
      res
    }
  } else {
    raw.results <- foreach(i = 1:nrow(grid), .packages = packages) %dopar% {
      print(sprintf("Simulation %g / %g", i, nrow(grid)))
      res <- SimIVDataCompareEstimatorsMultivar(type = est.eqns, n.sims = grid$n.data[i], 
                                                sample.size = grid$sample.size[i], 
                                                conf.corr.X = grid$Conf.Corr.X[i], 
                                                conf.corr.Y = grid$Conf.Corr.Y[i], 
                                                instrument.strength = grid$Instrument.Strength[i], 
                                                lambda = grid$lambda[i], beta = beta, seed = seed,
                                                survival.distribution = survival.distribution)
      for (j in 1:ncol(grid)) {attr(res, colnames(grid)[j]) <- grid[i, j]}
      res
    }
  }
  sum.results <- lapply(raw.results, function(res) {
    cors <- c(attr(res, "avg.cor")[3,4], attr(res, "avg.cor")[2,4], attr(res, "avg.cor")[1,2])
    cors <- matrix(rep(cors,length(est.eqns)), ncol = length(cors), byrow = T)
    cbind(cors, as.matrix(summary(res)[,2:10]))
  })
  for (i in 1:length(sum.results)) {results[, i, (ncol(grid) + 1):(ncol(grid) + 12)] <- sum.results[[i]]}
  for (i in 1:dim(results)[1]){results[i, , (1:ncol(grid))] <- as.matrix(grid)}
  return.results <- list(raw = raw.results, summary.array = results)
  class(return.results) <- "simulation.grid"
  return.results
}




fitAFT <- function(data, est.eqn = NULL, instrument.names, confounded.x.names, 
                   init.par = NULL, init.method = c("lm", "bisection"),
                   fit.method = c("dfsane", "multiStart", "nleqslv", "sane"), ...) {
  
  require(BB); require(nleqslv)
  fit.method <- match.arg(fit.method)
  if (fit.method == "dfsane" & is.null(est.eqn)) {
    stop("Please supply est.eqn if using fit.method: df.sane")
  }
  init.method <- match.arg(init.method)
  #instr.loc <- match(instrument.names, colnames(data$X))
  conf.x.loc <- match(confounded.x.names, colnames(data$X))
  ZXmat <- data$X

  if (attr(est.eqn, "name") == "AFT2SLSScorePre" | attr(est.eqn, "name") == "AFT2SLSScoreSmoothPre"
      | attr(est.eqn, "name") == "AFT2SLSScoreAllPre" | attr(est.eqn, "name") == "AFT2SLSScoreSmoothAllPre") {
    #if 2SLS is used, replace Z with Xhat 
    Z <- as.matrix(data$Z)
    if (attr(est.eqn, "name") == "AFT2SLSScoreAllPre" | attr(est.eqn, "name") == "AFT2SLSScoreSmoothAllPre") {
      for (v in 1:ncol(data$X)) {
        data$X[,v] <- lm(data$X[,v] ~ Z)$fitted.values
      }
    } else {
      for (v in 1:length(conf.x.loc)) {
        dat.lm <- data.frame(response = data$X[,conf.x.loc[v]], predictor = Z)
        Xhat <- lm(response ~ predictor, data = dat.lm)$fitted.values
        data$X[,conf.x.loc[v]] <- Xhat
      }
    }
    if (attr(est.eqn, "name") == "AFT2SLSScoreSmoothPre" | attr(est.eqn, "name") == "AFT2SLSScoreSmoothAllPre") {
      est.eqn <- AFTScoreSmoothPre
      attr(est.eqn, "name") <- "AFTScoreSmoothPre"
    } else {
      est.eqn <- AFTScorePre
      attr(est.eqn, "name") <- "AFTScorePre"
    }

  } else {
    ZXmat[,conf.x.loc] <- data$Z
  }

  
  if (is.null(init.par)) {
    #Xhat <- lm(data$X[,conf.x.loc] ~ data$Z)$fitted.values
    XXhat <- as.matrix(data$X)
    #XXhat[,conf.x.loc] <- Xhat
    init.par <- lm(data$survival$log.t ~ XXhat-1)$coefficients
    init.par[which(is.na(init.par))] <- 0
    names(init.par) <- colnames(data$X)
    #num.vars <- length(init.par)
    #names(init.par)[conf.x.loc] <- instrument.names
    if (init.method == "bisection") {
      interval <- list()
      for (i in 1:num.vars) {
        interval[[i]] <- c(init.par[i] - 1, init.par[i] + 1)
      }
      init.par <- coordinateBisection(est.eqn, interval = interval, num.vars = num.vars, 
                                      max.iter = 500, bisection.tol = 0.25,
                                      survival=data$survival, X = as.matrix(data$X), ZXmat = as.matrix(ZXmat))
      print (init.par)
      print (sum(est.eqn(init.par, survival=data$survival, X = as.matrix(data$X), ZXmat = as.matrix(ZXmat)) ^ 2))
    }
    names(init.par) <- colnames(data$X)
  }
  

  
  if (fit.method == "dfsane") {
    #Derivative-Free Spectral Approach for solving nonlinear systems of equations
    #from CRAN package 'BB'
    if (attr(est.eqn, "name") == "AFTivIPCWScorePre") {
      #generate function G_c() for ICPW 
      GC <- genKMCensoringFunc(data$survival)
      df <- BBsolve(par = init.par, fn = est.eqn, ...,
                    survival=data$survival, X = as.matrix(data$X), ZXmat = as.matrix(ZXmat), GC = GC)
      fval <- est.eqn(beta = df$par, survival=data$survival, X = as.matrix(data$X), ZXmat = as.matrix(ZXmat), GC = GC)
    } else {
      df <- BBsolve(par = init.par, fn = est.eqn, ...,
                    survival=data$survival, X = as.matrix(data$X), ZXmat = as.matrix(ZXmat))
      fval <- est.eqn(beta = df$par, survival=data$survival, X = as.matrix(data$X), ZXmat = as.matrix(ZXmat))
    }
    ret <- list(par = df$par, fval = fval, iter = df$iter)
  } else if (fit.method == "nleqslv" | fit.method == "sane") {
    if (!is.null(est.eqn)) {
      if (attr(est.eqn, "name") != "AFTScoreSmoothPre" & attr(est.eqn, "name") != "AFTivScoreSmoothPre") {
        
        warning(paste(paste("Arg: est.eqn =", attr(est.eqn, "name")), 
                                            "not be used. Used AFTivScoreSmoothPre instead"))
        est.eqn <- AFTivScoreSmoothPre
      }
    } else {
      warning("est.eqn not supplied. Used AFTivScoreSmoothPre")
      est.eqn <- AFTivScoreSmoothPre
    }
    if (fit.method == "nleqslv") {
      df <- nleqslv(x = init.par, fn = est.eqn, ..., 
                    survival=data$survival, X = as.matrix(data$X), ZXmat = as.matrix(ZXmat), tau = 0.01)
      ret <- list(par = df$x, fval = df$fvec, iter = df$iter)
      #print(df)
    } else if (fit.method == "sane") {
      df <- sane(par = init.par, fn = est.eqn, ...,
                 survival=data$survival, X = as.matrix(data$X), ZXmat = as.matrix(ZXmat), tau = 0.01)
      fval <- est.eqn(beta = df$par, survival=data$survival, X = as.matrix(data$X), ZXmat = as.matrix(ZXmat), tau = 0.01)
      ret <- list(par = df$par, fval = fval, iter = df$iter)
    }
  } else if (fit.method == "multiStart") {
    #Derivative-Free Spectral Approach for solving nonlinear systems of equations
    #from CRAN package 'BB'
    n.starts <- 10
    p <- length(init.par)
    p0 <- matrix(rnorm(n.starts * p), n.starts, p) 
    p0 <- rbind(p0, init.par)
    if (attr(est.eqn, "name") == "AFTivIPCWScorePre") {
      #generate function G_c() for ICPW 
      GC <- genKMCensoringFunc(data$survival)
      df <- multiStart(par = p0, fn = est.eqn, ..., action = "solve",
                       survival=data$survival, X = as.matrix(data$X), ZXmat = as.matrix(ZXmat), GC = GC)
      fval <- est.eqn(beta = df$par, survival=data$survival, X = as.matrix(data$X), ZXmat = as.matrix(ZXmat), GC = GC)
    } else {
      df <- multiStart(par = p0, fn = est.eqn, ..., action = "solve",
                       survival=data$survival, X = as.matrix(data$X), ZXmat = as.matrix(ZXmat))
      #take only the best ones and fit again
      good.init.idx <- which(df$fvalue < 500)
      p0 <- df$par[good.init.idx,] + matrix(rnorm(length(good.init.idx) * p, sd = 5e-6), length(good.init.idx), p) 
      
      df <- multiStart(par = p0, fn = est.eqn, ..., action = "solve",
                       survival=data$survival, X = as.matrix(data$X), ZXmat = as.matrix(ZXmat))
      
      good.init.idx <- which(df$fvalue < 100)
      df$par <- df$par[good.init.idx,]
      df$fvalue <- df$fvalue[good.init.idx]
      colnames(p0) <- colnames(data$X)
      fval <- est.eqn(beta = df$par, survival=data$survival, X = as.matrix(data$X), ZXmat = as.matrix(ZXmat))
    }
    ret <- list(par = df$par, fval = fval, iter = df$iter)
  }
  ret$sum.sq.fval <- sum(ret$fval^2)
  if (init.method == "multiStart") {ret <- df}
  ret$call <- match.call()
  class(ret) <- "aft.fit"
  ret
}

repFitAFT <- function(tol = 5, maxit = 25, data, est.eqn = NULL, est.eqn.sm = NULL, instrument.names, confounded.x.names,
                      init.par = NULL, init.method = c("lm", "bisection"), final.fit = T,
                      fit.method = c("dfsane", "multiStart", "nleqslv", "sane"), ...) {
  # repeatedly calls fitAFT until convergence. decreasing noise 
  # is added to coefficients at each call to encourage solution
  ct <- 0
  ssf <- best.ssf <- 1e10
  while (ssf > tol & ct <= maxit) {
    ct <- ct + 1
    old.ssf <- ssf
    
    #solve for beta using deriv-free spectral method
    est <- fitAFT(data = data, est.eqn = est.eqn, 
                  instrument.names = instrument.names, confounded.x.names = confounded.x.names, 
                  fit.method = fit.method, init.par = init.par, ...)
    ssf <- est$sum.sq.fval
    sd <- 5e-6 * min(sqrt(best.ssf), 5e2)
    if (ssf < best.ssf) {
      best.est <- est
      best.ssf <- ssf
      init.par <- best.est$par + rnorm(length(est$par), sd = sd)
    } else {
      init.par <- best.est$par + rnorm(length(est$par), sd = sd)
    }
    print (sprintf("Current ssf: %g   Best ssf: %g, sd: %g", ssf, best.ssf, sd))
  }
  if (final.fit) {
    best.est <- fitAFT(data = data, est.eqn = est.eqn.sm, 
                       instrument.names = instrument.names, confounded.x.names = confounded.x.names, 
                       fit.method = "nleqslv", init.par = init.par, global = "dbldog", method = "Broyden", 
                       control = list(ftol=1e-3, btol=1e-4, trace=1))
  }
  best.est$call <- match.call()
  best.est
}


bootstrapCI <- function(aft.fit, data, n.bootstraps = 999, percent, packages = NULL, func.names) {
  require(foreach)
  stopifnot(class(data) == "survival.data")
  stopifnot(class(aft.fit) == "aft.fit")
  
  n.bootstraps <- as.integer(n.bootstraps)
  call <- aft.fit$call
  #initialize with coefficient estimates from initial
  #fit to speed up fitting process
  
  newmaxit <- 32
  call[[match("maxit", names(call))]] <- as.name("newmaxit")
  coefficients <- aft.fit$par
  if (!is.null(call[[match("init.par", names(call))]])) {
    call[[match("init.par", names(call))]] <- as.name("coefficients")
  } else {
    call$init.par <- as.name("coefficients")
  }
  boot.estimates <- foreach(i = 1:n.bootstraps, .packages = packages, .combine = cbind, 
                            .export = c(func.names, "newmaxit", "coefficients")) %dopar% {
    print(sprintf("Bootstrap %g / %g", i, n.bootstraps))
    set.seed(i+123)
    #resample data
    samp.idx <- sample(1:nrow(data$X), nrow(data$X), replace = T)
    data.samp <- data
    data.samp$X <- data.samp$X[samp.idx,]
    data.samp$Z <- data.samp$Z[samp.idx]
    data.samp$survival <- data.samp$survival[samp.idx,]
    call[[match("data", names(call))]] <- as.name("data.samp")
    #fit coefficients with resampled data
    gc()
    boot.par <- eval(call)$par
    boot.par
  }
  
  ret.mat <- data.frame(array(0, dim = c(length(coefficients), 5)))
  names(ret.mat) <- c("Name", "Beta", "ci.lower", "ci.upper", "se")
  
  for (i in 1:length(coefficients)) {
    boot.est <- boot.estimates[i,]
    alpha.o.2 <- (1 - percent) / 2
    basic.quants <- quantile(boot.est, probs = c(1 - alpha.o.2, alpha.o.2))
    basic.CIs <- 2 * coefficients[i] - basic.quants
    se.hat <- sd(boot.est)
    ret.mat[i,2:ncol(ret.mat)] <- c(coefficients[i], basic.CIs, se.hat)
  }
  ret.mat$Name <- names(coefficients)
  ret = list(results = ret.mat, boot.est = boot.estimates)
  class(ret) <- "bootstrap.results"
  ret
}

print.bootstrap.results <- function(bsr, true.beta = NULL) {
  res <- bsr$results
  if (!is.null(true.beta)) {
    res$True.Beta <- true.beta
    res <- res[,c(1,ncol(res),2:(ncol(res) - 1))]
  }
  res[,2:ncol(res)] <- round(res[,2:ncol(res)], 4)
  print(res)
}

