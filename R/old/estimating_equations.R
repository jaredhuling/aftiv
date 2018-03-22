
#################################
### This file contains all    ###
### functions for est eqns    ###
#################################





svBootstrapOld <- function(beta, esteqn, B, nobs, GC = NULL, VarFunc = NULL, ...) 
{
  p  <- length(beta)
  Mb <- matrix(rexp(nobs * B, rate = 1), ncol = B)
  
  is.ipcw <- attr(esteqn, "name") == "evalAFTivIPCWScorePrec" | 
    attr(esteqn, "name") == "AFTivIPCWScorePre"
  
  if (is.ipcw) 
  {
    #time.t <- list(...)[["data.simu"]]$t
    #Un1 <- t(apply( Mb, 2, function(x) esteqn(beta = beta, multiplier.wts = x, GC = GC, ...) ))
    
    
    #print("ipcw var1")
    #print(var(t(Un1) ))
    
    Un1 <- array(NA, dim = c(B, p) )
    if (attr(esteqn, "name") == "evalAFTivIPCWScorePrec") 
    {
      data <- list(...)[["data.simu"]]
      
      for (i in 1:B) 
      {
        samp.idx <- sample.int(nobs, nobs, replace = TRUE)
        GC.boot  <- genKMCensoringFunc(data[samp.idx,])
        Un1[i,]  <- esteqn(beta      = beta, 
                           GC        = GC.boot, 
                           data.simu = data[samp.idx,])
      }
    } else 
    {
      list.dat   <- list(...)
      survival   <- list.dat[["survival"]]
      X          <- list.dat[["X"]]
      ZXmat      <- list.dat[["ZXmat"]]
      conf.x.loc <- list.dat[["conf.x.loc"]]
      
      for (i in 1:B) 
      {
        samp.idx <- sample.int(nobs, nobs, replace = TRUE)
        GC.boot  <- genKMCensoringFunc(survival[samp.idx,])
        Un1[i,]  <- esteqn(beta       = beta, 
                           GC         = GC.boot, 
                           survival   = survival[samp.idx,],
                           X          = X[samp.idx,], 
                           ZXmat      = ZXmat[samp.idx,], 
                           conf.x.loc = conf.x.loc)
      }
    }
    #varG <- VarFunc(time.t)
    #varG <- sum(varG[!is.nan(varG)]^2) * nobs * nobs
  } else 
  {
    Un1 <- t(apply( Mb, 2, function(x) esteqn(beta = beta, multiplier.wts = x, ...) ))
    if (p == 1) 
    {
      Un1 <- t(Un1)
    }
  }
  
  
  if (is.ipcw) 
  {
    V <- var(Un1)
    #VG <- varG * diag(ncol(V))
    #V <- V + VG
  } else 
  {
    V <- var(Un1)
  }
  print(attr(esteqn, "name"))
  #print(V)
  #print(solve(V))
  ##############Zb <- mvrnorm(n = B, rep(0, p), Sigma = solve(V))
  
  if (p == 1) 
  {
    rt <- as.numeric(sqrt(V))
    if (is.ipcw) 
    {
      cat("V:", V, "\n")
      print(rt)
      print(1 / rt)
      Zb <- (1/rt) * (matrix(4 * rbinom(p * B, 1, 0.5) - 2, ncol = p))
      #Zb <- (1/rt) * (matrix(runif(p * B, min=-2, max=2), ncol = p))
      #Zb <- (1/rt) * (matrix(rt(p * B, df = 2), ncol = p))
      Zb <- (2/(rt * (log(nobs)^(0.525)))  ) * (matrix(rnorm(p * B), ncol = p))
    } else 
    {
      cat("V:", V, "\n")
      Zb <- (1/rt) * (matrix(rnorm(p * B), ncol = p))
    }
    #Zb <- rt * (matrix(2 * rnorm(p * B) - 1, ncol = p))
  } else 
  {
    eeg <- eigen(V)
    rt  <- solve(eeg$vectors %*% diag(sqrt(eeg$values) ) %*% solve(eeg$vectors))
    if (is.ipcw) 
    {
      #Zb <- 2 * t(rt %*% t(2 * matrix(rbinom(p * B, 1, 0.5), ncol = p) - 1))/(log(nobs)^(0.525))
      Zb <- t(rt %*% t(3 * matrix(rbinom(p * B, 1, 0.5), ncol = p) - 1.5))/(log(log(nobs)) )
      #Zb <- 2 * t(rt %*% t(2 * matrix(rbinom(p * B, 1, 0.5), ncol = p) - 1))
    } else 
    {
      Zb <- t(rt %*% t(2 * matrix(rbinom(p * B, 1, 0.5), ncol = p) - 1))
    }
  }
  
  if (is.ipcw) 
  {
    Un2 <- t(apply( Zb, 1, function(x) esteqn(beta = beta + x / sqrt(nobs),  GC = GC, ...) ))
    
  } else 
  {
    Un2 <- t(apply( Zb, 1, function(x) esteqn(beta = beta + x / sqrt(nobs), ...) ))
  }
  
  if (p == 1) 
  {
    Un2 <- t(Un2)
  }
  
  sigma.hat <- solve(var(Un2))
  if (is.ipcw) 
  {
    sigma.hat <- sigma.hat
    se.hat    <- sqrt(diag(sigma.hat)) / sqrt(nobs)
  } else 
  {
    se.hat    <- sqrt(diag(sigma.hat)) / sqrt(nobs)
  }
  cat("The est: ", beta, "\n")
  print("The variance: ")
  print(se.hat)
  
  list(se.hat = se.hat, var = sigma.hat, V = V)  
}






evalAFTScore <- function(data.simu, beta, multiplier.wts = NULL)
{ 
  #the score function for the AFT model
  
  #transform to log-time
  data.simu$t <- log(data.simu$t)
  
  #store the T_i - bX_i term (error)
  data.simu <- data.frame(data.simu, err = data.simu$t - beta * data.simu$X)
  
  #sort according to error size ####observed failure time 
  #data.simu <- data.simu[order(data.simu$X),] 
  ord.idx <- order(data.simu$err)
  data.simu <- data.simu[ord.idx,] 
  
  #the denominator of the at-risk comparison term  
  data.simu <- data.frame(data.simu, at.risk.terms = nrow(data.simu):1)
  
  if (!is.null(multiplier.wts)) {
    multiplier.wts <- multiplier.wts[ord.idx]
    #the numerator of the at-risk comparison term  
    data.simu <- data.frame(data.simu, at.risk.X.terms = cumsumRev(data.simu$X * multiplier.wts))
    
    #return the score   
    return(sum(data.simu$delta * (data.simu$at.risk.terms * data.simu$X * multiplier.wts - data.simu$at.risk.X.terms)) / sqrt(nrow(data.simu)))
  } else {
    data.simu <- data.frame(data.simu, at.risk.X.terms = cumsumRev(data.simu$X))
    
    #return the score   
    return(sum(data.simu$delta * (data.simu$at.risk.terms * data.simu$X - data.simu$at.risk.X.terms)) / sqrt(nrow(data.simu)))
  }
  
}
vEvalAFTScore <- Vectorize(evalAFTScore, vectorize.args = "beta")
attr(vEvalAFTScore, "name") <- "evalAFTScore"

evalAFTivScore <- function(data.simu, beta, multiplier.wts = NULL)
{ 
  #the score function for the AFT model with instrumental variables
  
  #transform to log-time
  data.simu$t <- log(data.simu$t)
  
  #store the T_i - bX_i term (error)
  data.simu$err = data.simu$t - beta * data.simu$X
  
  #sort according to error size ####observed failure time 
  #data.simu <- data.simu[order(data.simu$X),]   
  ord.idx <- order(data.simu$err)
  data.simu <- data.simu[ord.idx,] 
  
  #the denominator of the at-risk comparison term  
  data.simu$at.risk.terms <- nrow(data.simu):1
  
  if (!is.null(multiplier.wts)) {
    multiplier.wts <- multiplier.wts[ord.idx]
    #the numerator of the at-risk comparison term  
    data.simu$at.risk.Z.terms <- cumsumRev(data.simu$Z * multiplier.wts)
    
    #return the score   
    return(sum(data.simu$delta * (data.simu$at.risk.terms * data.simu$Z * multiplier.wts - data.simu$at.risk.Z.terms)) / sqrt(nrow(data.simu)))
  } else {
    #the numerator of the at-risk comparison term  
    data.simu$at.risk.Z.terms <- cumsumRev(data.simu$Z)
    
    #return the score   
    return(sum(data.simu$delta * (data.simu$at.risk.terms * data.simu$Z - data.simu$at.risk.Z.terms)) / sqrt(nrow(data.simu)))
  }
  
}
vEvalAFTivScore <- Vectorize(evalAFTivScore, vectorize.args = "beta")
attr(vEvalAFTivScore, "name") <- "evalAFTivScore"

evalAFTivIPCWScore <- function(data.simu, beta, multiplier.wts = NULL)
{ 
  #the score function for the AFT model with instrumental variables and Inverse Probability Censoring Weighting
  
  #generate function G_c() for ICPW
  GC <- genKMCensoringFunc(data.simu)
  
  data.simu$t.original <- data.simu$t
  
  #transform to log-time
  data.simu$t <- log(data.simu$t)
  
  #store the T_i - bX_i term (error)
  data.simu$err <- data.simu$t - beta * data.simu$X
  
  #store beta * X
  data.simu$bX <- beta * data.simu$X
  
  #creates G_c(T_i) term in the IPCW estimating equation
  data.simu$GCT <- GC(exp(data.simu$t))
  
  #sort according to error size ####observed failure time 
  #data.simu <- data.simu[order(data.simu$X),]   
  ord.idx <- order(data.simu$err)
  data.simu <- data.simu[ord.idx,] 
  
  #create indicator to set to zero terms where GCT == 0 and 
  #set to 1 so no dividing by zero occurs
  zero.indicator <- 1 * (data.simu$GCT != 0)
  data.simu$GCT[which(data.simu$GCT == 0)] <- 1
  
  if (!is.null(multiplier.wts)) {
    multiplier.wts <- multiplier.wts[ord.idx]
  }
  
  #the denominator of the at-risk comparison term  
  data.simu <- genIPCWNumDenom(data.simu, GC, multiplier.wts)
  
  #system.time(data.simu1 <- genIPCWNumDenom2(data.simu, GC))
  #system.time(data.simu2 <- genIPCWNumDenom(data.simu, GC))
  #system.time(data.simu2 <- genIPCWNumDenomSlower(data.simu, GC))
  #all(data.simu1$ICPW.at.risk.terms == data.simu2$ICPW.at.risk.terms)
  
  if (is.null(multiplier.wts)) {
    #return the score   
    return(sum(zero.indicator * (data.simu$delta / data.simu$GCT) * 
          (data.simu$IPCW.at.risk.terms * data.simu$Z - data.simu$IPCW.at.risk.Z.terms)) / (sqrt(nrow(data.simu))))
  } else {
    #return the score   
    return(sum(zero.indicator * (data.simu$delta / data.simu$GCT) * 
          (data.simu$IPCW.at.risk.terms * data.simu$Z * multiplier.wts - data.simu$IPCW.at.risk.Z.terms)) / (sqrt(nrow(data.simu))))
  }
}
vEvalAFTivIPCWScore <- Vectorize(evalAFTivIPCWScore, vectorize.args = "beta")
attr(vEvalAFTivIPCWScore, "name") <- "evalAFTivIPCWScore"


evalAFTivIPCWScoreSubset <- function(data.simu, beta)
{ 
  #the score function for the AFT model with instrumental variables and Inverse Probability Censoring Weighting
  
  data.simu$t.original <- data.simu$t
  
  #transform to log-time
  data.simu$t <- log(data.simu$t)
  
  #store beta * X
  data.simu$bX <- beta * data.simu$X
  
  bxquant <- quantile(data.simu$bX, probs = 0.95)
  
  #data.simu <- data.simu[which(data.simu$bX < bxquant),]
  
  
  #store the T_i - bX_i term (error)
  data.simu$err <- data.simu$t - beta * data.simu$X
  
  errquant <- quantile(data.simu$err, probs = 0.99)
  
  #data.simu <- data.simu[which(data.simu$err < errquant),]

  #generate function G_c() for ICPW
  GC <- genKMCensoringFunc(data.simu)
  #GC <- dexp
  
  #creates G_c(T_i) term in the IPCW estimating equation
  data.simu$GCT <- GC(data.simu$t.original)
  
  #sort according to error size ####observed failure time 
  #data.simu <- data.simu[order(data.simu$X),]   
  data.simu <- data.simu[order(data.simu$err),] 
  
  #create indicator to set to zero terms where GCT == 0 and 
  #set to 1 so no dividing by zero occurs
  zero.indicator <- 1 * (data.simu$GCT != 0)
  data.simu$GCT[which(data.simu$GCT == 0)] <- 1
  
  #the denominator of the at-risk comparison term  
  data.simu <- genIPCWNumDenom(data.simu, GC)
  #system.time(data.simu1 <- genIPCWNumDenom2(data.simu, GC))
  #system.time(data.simu2 <- genIPCWNumDenom(data.simu, GC))
  #system.time(data.simu2 <- genIPCWNumDenomSlower(data.simu, GC))
  #all(data.simu1$ICPW.at.risk.terms == data.simu2$ICPW.at.risk.terms)
  
  data.simu$IPCW.at.risk.terms <- ifelse(data.simu$IPCW.at.risk.terms == 0, 1, data.simu$IPCW.at.risk.terms)
  
  #return the score   
  sum(zero.indicator * (data.simu$delta / data.simu$GCT) * 
        (data.simu$Z - data.simu$IPCW.at.risk.Z.terms / data.simu$IPCW.at.risk.terms)) / sqrt(nrow(data.simu))
}
vEvalAFTivIPCWScoreSubset <- Vectorize(evalAFTivIPCWScoreSubset, vectorize.args = "beta")


genIPCWNumDenom <- cmpfun(function(dat, GC.func, multiplier.wts = NULL){
  #dat is a data.frame
  #GC.func is a function
  num <- denom <- array(0, dim = c(nrow(dat),1))
  bX <- dat$bX
  if (is.null(multiplier.wts)) {
    Z <- dat$Z
  } else {
    #bX <- dat$bX * multiplier.wts
    Z <- dat$Z * multiplier.wts
  }
  
  err <- dat$err
  len <- nrow(dat)
  for (i in 1:len){
    err.i <- err[i]
    #num.tmp <- denom.tmp <- 0
    #for (j in i:nrow(dat)){
    #  ipcw <- GC.func(dat$bX[j] + err.i)
    #  num.tmp <- num.tmp + dat$Z[j] / ipcw
    #  denom.tmp <- denom.tmp + 1 / ipcw
    #}
    ind.zero <- F
    ipcw <- GC.func(exp(bX[i:len] + err.i))
    if (all(ipcw == 0)){
      ipcw <- rep(0.0001, length(ipcw))
      ind.zero <- T
    }
    #ind.to.zero <- 1 * (ipcw != 0) 
    ind.to.zero <- 1
    ipcw[which(ipcw == 0)] <- min(ipcw[which(ipcw != 0)]) / 5
    if (!ind.zero){
      num[i] <- sum((Z[i:len] / ipcw) * ind.to.zero)
      denom[i] <- sum((1 / ipcw) * ind.to.zero)
    } else {num[i] <- denom[i] <- 1e5}
  }
  dat$IPCW.at.risk.Z.terms <- num
  dat$IPCW.at.risk.terms <- denom
  dat
})

genIPCWNumDenomSubset <- cmpfun(function(dat, GC.func){
  #dat is a data.frame
  #GC.func is a function
  num <- denom <- array(0, dim = c(nrow(dat),1))
  bX <- dat$bX
  Z <- dat$Z
  err <- dat$err
  len <- nrow(dat)
  for (i in 1:len){
    err.i <- err[i]
    #num.tmp <- denom.tmp <- 0
    #for (j in i:nrow(dat)){
    #  ipcw <- GC.func(dat$bX[j] + err.i)
    #  num.tmp <- num.tmp + dat$Z[j] / ipcw
    #  denom.tmp <- denom.tmp + 1 / ipcw
    #}
    ind.zero <- F
    func.eval <- exp(bX[i:len] + err.i)
    #ipcw.quant <- quantile(func.eval, probs = 0.95)
    ipcw <- GC.func(func.eval)
    if (all(ipcw == 0)){
      ipcw <- rep(0.01, length(ipcw))
      ind.zero <- T
    }
    ipcw[which(ipcw == 0)] <- min(ipcw[which(ipcw != 0)]) / 2
    ipcw.quant <- quantile(ipcw, probs = 0.15)
    #ipcw <- ifelse(ipcw > 0.2, ipcw, 0.2)
    #ipcw.to.zero <- 1 * (ipcw <= ipcw.quant)
    #ipcw.to.zero <- 1 * (ipcw != 0)
    if (!ind.zero){
      num[i] <- sum((Z[i:len] / ipcw))
      denom[i] <- sum((1 / ipcw))
    } else {num[i] <- denom[i] <- 0}
  }
  dat$IPCW.at.risk.Z.terms <- num
  dat$IPCW.at.risk.terms <- denom
  dat
})



genIPCWNumDenom2 <- cmpfun(function(dat, GC.func){
  #dat is a data.frame
  #GC.func is a function
  num <- denom <- array(0, dim = c(nrow(dat),1))
  unlist(lapply(1:nrow(dat), function(i) {
    err.i <- dat$err[i]
    ind.zero <- F
    ipcw <- GC.func(exp(dat$bX[i:nrow(dat)] + err.i))
    if (all(ipcw == 0)){
      ipcw <- rep(0.01, length(ipcw))
      ind.zero <- T
    }
    ipcw[which(ipcw == 0)] <- min(ipcw[which(ipcw != 0)]) / 2
    if (!ind.zero){
      num[i] <- sum(dat$Z[i:nrow(dat)] / ipcw)
      denom[i] <- sum(1 / ipcw)
    } else {num[i] <- denom[i] <- 0}
  }))
  dat$IPCW.at.risk.Z.terms <- num
  dat$IPCW.at.risk.terms <- denom
  dat
})




genIPCWNumDenomSlower <- cmpfun(function(dat, GC.func){
  #dat is a data.frame
  #GC.func is a function
  res <- unlist(lapply(dat$err, function(x) {
    idx <- which(dat$err >= x)
    ipcw <- GC.func(exp(dat$bX[idx] + x))
    if (all(ipcw == 0)){ipcw <- rep(0.001, length(ipcw))}
    ipcw[which(ipcw == 0)] <- min(ipcw[which(ipcw != 0)]) / 2
    c(sum(dat$Z[idx] / ipcw), sum(1 / ipcw))
  }))
  idx <- seq(1,length(res) - 1,by=2)
  dat$IPCW.at.risk.Z.terms <- res[idx]
  dat$ICPW.at.risk.terms <- res[idx + 1]
  dat
})

# genKMCensoringFunc <- function(data, cox = FALSE, X = NULL){
#   #returns the G_c function for inverse probability censored weighting
#   #by obtaining the Kaplan-Meier estimate for censoring and returning
#   #a step function of the estimate
#   require(survival)
#   data$delta.G <- 1 - data$delta
#   if (!cox)
#   {
#     if (is.null(data$t.original)) {
#       if (is.null(data$t)) {
#         cens.fit <- survfit(Surv(exp(log.t), delta.G) ~ 1, data = data)
#         cens.dist <- data.frame(c(0, summary(cens.fit)$time), c(1, summary(cens.fit)$surv))
#       } else {
#         cens.fit <- survfit(Surv(t, delta.G) ~ 1, data = data)
#         cens.dist <- data.frame(c(min(data$t), summary(cens.fit)$time), c(1, summary(cens.fit)$surv))
#       }
#     } else {
#       cens.fit <- survfit(Surv(t.original, delta.G) ~ 1, data = data)
#       cens.dist <- data.frame(c(min(data$t.original), summary(cens.fit)$time), c(1, summary(cens.fit)$surv))
#     }
#     
#     tmpFunc <- stepfun(cens.dist[,1], c(1,cens.dist[,2]), f=0)
#     attr(tmpFunc, "cox") <- FALSE
#   } else 
#   {
#     if(is.null(X)) stop("Must provide covariates X to be modeled if cox = TRUE")
#     
#     if (is.null(data$t.original)) {
#       if (is.null(data$t)) 
#       {
#         cph <- coxph(Surv(exp(data$log.t), data$delta.G) ~ X)
#         bh <- basehaz(cph)
#         baseline.survival <- exp(-bh$hazard)
#         cens.dist.cph <- data.frame(c(0, bh$time), 
#                                     c(1, baseline.survival))
#       } else 
#       {
#         cph <- coxph(Surv(data$t, data$delta.G) ~ X)
#         bh <- basehaz(cph)
#         baseline.survival <- exp(-bh$hazard)
#         cens.dist.cph <- data.frame(c(min(data$t), bh$time), 
#                                     c(1, baseline.survival))
#       }
#     } else {
#       cph <- coxph(Surv(data$t.original, data$delta.G) ~ X)
#       bh <- basehaz(cph)
#       baseline.survival <- exp(-bh$hazard)
#       cens.dist.cph <- data.frame(c(min(data$t.original), bh$time), 
#                                   c(1, baseline.survival))
#     }
#     
#     
#     # S_0(t) = exp(-Lambda(t))
#     # S(t | X) = S_0(t) ^ exp(x'beta)
#     
#     tmpFunc.cph <- stepfun(cens.dist.cph[,1], c(1,cens.dist.cph[,2]), f = 0)
#     
#     coef.cph <- coef(cph)
#     names(coef.cph) <- NULL
#     
#     ## return function to compute estimated not-censoring
#     ## probability conditional on x (covariates)
#     tmpFunc <- function(t, x)
#     {
#       if (is.matrix(x))
#       {
#         return( tmpFunc.cph(t) ^ exp(drop(x %*% coef.cph)) )
#       } else 
#       {
#         return( tmpFunc.cph(t) ^ exp(   sum(x * coef.cph)) )
#       }
#     }
#     attr(tmpFunc, "cox") <- TRUE
#   }
# 
#   tmpFunc
# }


AFTivIPCWScore <- function(beta, data.simu)
{ 
  if (class(data.simu) != "survival.data") {stop("Need to use data generated by simIVMultivarSurvivalData")}
  
  #the score function for the AFT model with instrumental variables and Inverse Probability Censoring Weighting
  survival <- data.simu$survival
  X <- data.simu$X
  Z <- data.simu$Z
  n <- nrow(X)
  nvars <- ncol(X)
  
  #generate function G_c() for ICPW
  GC <- genKMCensoringFunc(survival)
  
  #transform to log-time
  survival$t <- log(survival$t)
  
  
  #creates G_c(T_i) term in the IPCW estimating equation
  survival$GCT <- GC(exp(survival$t))
  
  #store the T_i - bX_i term (error)
  survival$err = survival$t - X %*% beta
  
  survival$bX <- X %*% beta
  
  #sort according to error size ####observed failure time 
  #data.simu <- data.simu[order(data.simu$X),]  
  order.idx <- order(survival$err)
  survival <- survival[order.idx,] 
  X <- as.matrix(X[order.idx,])
  Z <- as.matrix(Z[order.idx,]) 
  
  
  #create indicator to set to zero terms where GCT == 0 and 
  #set to 1 so no dividing by zero occurs
  zero.indicator <- 1 * (survival$GCT != 0)
  survival$GCT[which(survival$GCT == 0)] <- 1
  
  #first col is as.risk.terms, remaining are at.risk.Z.terms
  at.risk.mat <- genIPCWNumDenomMultivar(survival, Z, GC)
  
  #generate empty vector to return eventually
  ret.vec <- numeric(nvars)
  for (i in 1:nvars) {
    #return the score   
    ret.vec[i] <- sum(zero.indicator * (survival$delta / survival$GCT) * (at.risk.mat[,1] * Z[,i] - at.risk.mat[,i + 1])) / sqrt(n)
  }
  
  ret.vec

}

AFTivIPCWScorePreKM <- function(beta, data.simu, GC)
{ 
  if (class(data.simu) != "survival.data") {stop("Need to use data generated by simIVMultivarSurvivalData")}
  
  #the score function for the AFT model with instrumental variables and Inverse Probability Censoring Weighting
  survival <- data.simu$survival
  X <- data.simu$X
  Z <- data.simu$Z
  n <- nrow(X)
  nvars <- ncol(X)
  
  #transform to log-time
  tt <- log(survival$t)
  
  #creates G_c(T_i) term in the IPCW estimating equation
  GCT <- GC( survival$t )
  
  bX <- X %*% beta
  
  #store the T_i - bX_i term (error)
  err = survival$t - bX
  
  
  
  #sort according to error size ####observed failure time 
  #data.simu <- data.simu[order(data.simu$X),]  
  order.idx <- order(err)
  survival <- survival[order.idx,]
  tt <- tt[order.idx]
  bX <- bX[order.idx]
  err <- err[order.idx]
  X <- as.matrix(X[order.idx,])
  Z <- as.matrix(Z[order.idx,]) 
  
  
  #create indicator to set to zero terms where GCT == 0 and 
  #set to 1 so no dividing by zero occurs
  zero.indicator <- 1 * (GCT != 0)
  GCT[which(GCT == 0)] <- 1
  
  #first col is as.risk.terms, remaining are at.risk.Z.terms
  at.risk.mat <- genIPCWNumDenomMultivar2(bX, Z, err, GC)
  
  #generate empty vector to return eventually
  ret.vec <- numeric(nvars)
  for (i in 1:nvars) {
    #return the score   
    ret.vec[i] <- sum(zero.indicator * (survival$delta / GCT) * (at.risk.mat[,1] * Z[,i] - at.risk.mat[,i + 1])) / sqrt(n)
  }
  
  ret.vec
  
}




absAFTivScorePre <- function(beta, ...){sum(AFTivScorePre(beta, ...)^2)}



AFTScoreSmoothPreVec <- function(beta, survival, X, ZXmat, tau = 1e-3)
{ 
  #faster version for smaller datasets. for large datasets
  #this takes up too much memory
  #the score function for the AFT model
  
  n <- nrow(X)
  nvars <- ncol(X)
  
  #transform to log-time
  if (is.null(survival$log.t)) { 
    log.t <- log(survival$t)
  } else (log.t <- survival$log.t)
  
  delta <- survival$delta
  
  #store the T_i - bX_i term (error)
  err = as.vector(log.t - X %*% beta)
  
  #the denominator of the at-risk comparison term  
  out.sig <- sigmoid(outer(err, err, "-"), tau)
  at.risk.terms <- colSums(out.sig)
  
  # #generate empty vector to return eventually
  #   ret.vec <- numeric(nvars)
  #   for (i in 1:nvars) {
  #     #the numerator of the at-risk comparison term 
  #     at.risk.X.terms <- unlist(lapply(err, function(x) sum(X[,i] * sigmoid((err - x), tau))))
  #     
  #     #return the score   
  #     ret.vec[i] <- sum(delta * (at.risk.terms * X[,i] - at.risk.X.terms)) / sqrt(n)
  #   }
  ret.vec <- apply(X, 2, function(xi) {
    #at.risk.X.terms <- rowSums(sweep(out.sig, MARGIN = 2, xi, '*'))
    at.risk.X.terms <- colSums(out.sig * xi)
    sum(delta * (at.risk.terms * xi - at.risk.X.terms)) / sqrt(n)
  })
  
  ret.vec
}



absAFTivScoreSmoothPre <- function(...){sum(AFTivScoreSmoothPre(...)^2)}

evalAFT2SLSScore <- function(data.simu, beta, multiplier.wts = NULL)
{ 
  #the score function for the AFT model with instrumental variables
  
  #transform to log-time
  data.simu$t <- log(data.simu$t)
  
  #predict X with Z (1st stage of 2SLS process) to get Xhat
  data.simu$X.hat <- lm(X ~ Z, data = data.simu)$fitted.values
  
  #store the T_i - bX_i term (error)
  data.simu$err = data.simu$t - beta * data.simu$X.hat
  
  #sort according to error size ####observed failure time 
  #data.simu <- data.simu[order(data.simu$X),]   
  ord.idx <- order(data.simu$err)
  data.simu <- data.simu[ord.idx,] 
  
  #the denominator of the at-risk comparison term  
  data.simu$at.risk.terms <- nrow(data.simu):1
  if (is.null(multiplier.wts)) {
    
    #the numerator of the at-risk comparison term  
    data.simu$t.risk.Z.terms <- cumsumRev(data.simu$X.hat)
    
    #return the score   
    return(sum(data.simu$delta * (data.simu$at.risk.terms * data.simu$X.hat - data.simu$t.risk.Z.terms)) / sqrt(nrow(data.simu)))
    
  } else {
    multiplier.wts <- multiplier.wts[ord.idx]
    
    #the numerator of the at-risk comparison term  
    data.simu$t.risk.Z.terms <- cumsumRev(data.simu$X.hat * multiplier.wts)
    
    #return the score   
    return(sum(data.simu$delta * (data.simu$at.risk.terms * data.simu$X.hat * multiplier.wts - data.simu$t.risk.Z.terms)) / sqrt(nrow(data.simu)))
    
  }
  
}
vEvalAFT2SLSScore <- Vectorize(evalAFT2SLSScore, vectorize.args = "beta")
attr(vEvalAFT2SLSScore, "name") <- "evalAFT2SLSScore"

evalAFTivScoreSmooth <- function(data.simu, beta)
{ 
  #the score function for the AFT model with instrumental variables
  
  #transform to log-time
  data.simu$t <- log(data.simu$t)
  
  #store the T_i - bX_i term (error)
  data.simu$err = data.simu$t - beta * data.simu$X
  
  #sort according to error size ####observed failure time 
  #data.simu <- data.simu[order(data.simu$X),]   
  data.simu <- data.simu[order(data.simu$err),] 
  
  #get a rough estimate of sd for smoothed method
  stdev <- roughResidSDEstAFTIV(data.simu)
  
  #the denominator of the at-risk comparison term  
  data.simu$at.risk.terms <- unlist(lapply(data.simu[,"err"], function(x) 
    sum(pnorm(nrow(data.simu)^(0.26) * (data.simu$err - x) / stdev))))
  
  #the numerator of the at-risk comparison term  
  data.simu$at.risk.X.terms <- unlist(lapply(data.simu[,"err"], function(x) 
    sum(data.simu$Z * pnorm(nrow(data.simu)^(0.26) * (data.simu$err - x) / stdev))))
  
  #return the score   
  sum(data.simu$delta * (data.simu$at.risk.terms * data.simu$Z - data.simu$at.risk.Z.terms)) / sqrt(nrow(data.simu))
}
vEvalAFTivScoreSmooth <- Vectorize(evalAFTivScoreSmooth, vectorize.args = "beta")

evalAFTScoreSmooth <- function(beta, data.simu, stdev)
{ 
  #the score function for the AFT model
  
  n <- nrow(data.simu)
  
  #transform to log-time
  data.simu$t <- log(data.simu$t)
  
  #store the T_i - bX_i term (error)
  data.simu <- data.frame(data.simu, err = data.simu$t - beta * data.simu$X)
  
  #sort according to error size ####observed failure time 
  #data.simu <- data.simu[order(data.simu$X),]   
  data.simu <- data.simu[order(data.simu$err),] 
  
  #the denominator of the at-risk comparison term  
  data.simu$at.risk.terms <- unlist(lapply(data.simu[,"err"], function(x) 
    sum(pnorm(n^(0.26) * (data.simu$err - x) / stdev))))
  
  #the numerator of the at-risk comparison term  
  data.simu$at.risk.X.terms <- unlist(lapply(data.simu[,"err"], function(x) 
    sum(data.simu$X * pnorm(n^(0.26) * (data.simu$err - x) / stdev))))
  
  #return the score   
  sum(data.simu$delta * (data.simu$at.risk.terms * data.simu$X - data.simu$at.risk.X.terms)) / sqrt(n)
}
vEvalAFTScoreSmooth <- Vectorize(evalAFTScoreSmooth, vectorize.args = "beta")

AFTScoreSmooth <- function(beta, data.simu, stdev)
{ 
  if (class(data.simu) != "survival.data") {stop("Need to use data generated by simIVMultivarSurvivalData")}
  
  #the score function for the AFT model
  
  survival <- data.simu$survival
  X <- data.simu$X
  n <- nrow(X)
  nvars <- ncol(X)
  
  #transform to log-time
  survival$t <- log(survival$t)
  
  #store the T_i - bX_i term (error)
  survival$err = survival$t - X %*% beta
  
  #sort according to error size ####observed failure time 
  #data.simu <- data.simu[order(data.simu$X),]  
  order.idx <- order(survival$err)
  survival <- survival[order.idx,] 
  X <- X[order.idx,] 
  
  #the denominator of the at-risk comparison term  
  survival$at.risk.terms <- unlist(lapply(survival$err, function(x) 
    sum(pnorm(n^(0.26) * (survival$err - x) / stdev))))
  
  #generate empty vector to return eventually
  ret.vec <- numeric(nvars)
  for (i in 1:nvars) {
    #the numerator of the at-risk comparison term 
    at.risk.X.terms <- unlist(lapply(survival$err, function(x) 
        sum(X[,i] * pnorm(n^(0.26) * (survival$err - x) / stdev))))
    
    #return the score   
    ret.vec[i] <- sum(survival$delta * (survival$at.risk.terms * X[,i] - at.risk.X.terms)) / sqrt(n)
  }
  
  ret.vec
}



jacobianAFTSmooth <- function(beta, data.simu, stdev)
{ 
  if (class(data.simu) != "survival.data") {stop("Need to use data generated by simIVMultivarSurvivalData")}
  
  #the score function for the AFT model
  
  survival <- data.simu$survival
  X <- data.simu$X
  n <- nrow(X)
  nvars <- ncol(X)
  
  #transform to log-time
  survival$t <- log(survival$t)
  
  #store the T_i - bX_i term (error)
  survival$err = survival$t - X %*% beta
  
  #sort according to error size ####observed failure time 
  #data.simu <- data.simu[order(data.simu$X),]  
  order.idx <- order(survival$err)
  survival <- survival[order.idx,] 
  X <- as.matrix(X[order.idx,]) 
  
  
  X_err <- cbind(X, survival$err)
  #generate empty vector to return eventually
  jacobian <- as.matrix(array(0, dim = c(nvars, nvars)))
  for (l in 1:nvars) {
    for (m in 1:nvars) {
      #the numerator of the at-risk comparison term 
      x.sq.terms <- apply(X_err[,c(l, m, nvars + 1)], 1, function(x) 
        sum((X[,l] - x[1]) * (X[,m]  - x[2]) * dnorm(n^(0.26) * (x[3] - survival$err) / stdev)))
      
      #return the J(l,m) element   
      jacobian[l, m] <- ((n^(0.26)) / stdev) * sum((survival$delta * (x.sq.terms))) / sqrt(n)
    }
  }
  
  jacobian
}

evalGradAFT <- function(beta, data.simu, stdev)
{ 
  #the score function for the AFT model
  
  #transform to log-time
  data.simu$t <- log(data.simu$t)
  
  #store the T_i - bX_i term (error)
  data.simu <- data.frame(data.simu, err = data.simu$t - beta * data.simu$X)
  
  #sort according to error size ####observed failure time 
  #data.simu <- data.simu[order(data.simu$X),]   
  data.simu <- data.simu[order(data.simu$err),] 
  
  #the numerator of the at-risk comparison term  
  data.simu$x.sq.terms <- apply(data.simu[,c("X","err")], 1, function(x) 
    sum((data.simu$X - x[1]) * (data.simu$X - x[1]) * dnorm(nrow(data.simu)^(0.26) * (x[2] - data.simu$err) / stdev)))

  #return the score   
  as.matrix(((nrow(data.simu)^(0.26)) / stdev) * sum((data.simu$delta * (data.simu$x.sq.terms))) / sqrt(nrow(data.simu)))
}
vEvalGradAFT <- Vectorize(evalGradAFT, vectorize.args = "beta")


AFTScore <- function(beta, data.simu)
{ 
  if (class(data.simu) != "survival.data") {stop("Need to use data generated by simIVMultivarSurvivalData")}
  
  #the score function for the AFT model
  
  survival <- data.simu$survival
  X <- data.simu$X
  n <- nrow(X)
  nvars <- ncol(X)
  
  #transform to log-time
  survival$t <- log(survival$t)
  
  #store the T_i - bX_i term (error)
  survival$err = survival$t - X %*% beta
  
  #sort according to error size ####observed failure time 
  #data.simu <- data.simu[order(data.simu$X),]  
  order.idx <- order(survival$err)
  survival <- survival[order.idx,] 
  X <- as.matrix(X[order.idx,])
  
  #the denominator of the at-risk comparison term  
  at.risk.terms <- n:1
  
  ret.vec <- numeric(nvars)
  for (i in 1:nvars) {
    #the numerator of the at-risk comparison term  
    at.risk.X.terms <- cumsumRev(X[,i])
    
    #return the score   
    ret.vec[i] <- sum(survival$delta * (at.risk.terms * X[,i] - at.risk.X.terms)) / sqrt(n)
  }  
  ret.vec
}

absAFTScore <- function(beta, data.simu)
{ 
  if (class(data.simu) != "survival.data") {stop("Need to use data generated by simIVMultivarSurvivalData")}
  
  #the score function for the AFT model
  
  survival <- data.simu$survival
  X <- data.simu$X
  n <- nrow(X)
  nvars <- ncol(X)
  
  #transform to log-time
  survival$t <- log(survival$t)
  
  #store the T_i - bX_i term (error)
  survival$err = survival$t - X %*% beta
  
  #sort according to error size ####observed failure time 
  #data.simu <- data.simu[order(data.simu$X),]  
  order.idx <- order(survival$err)
  survival <- survival[order.idx,] 
  X <- as.matrix(X[order.idx,])
  
  #the denominator of the at-risk comparison term  
  at.risk.terms <- n:1
  
  ret.vec <- numeric(nvars)
  for (i in 1:nvars) {
    #the numerator of the at-risk comparison term  
    at.risk.X.terms <- cumsumRev(X[,i])
    
    #return the score   
    ret.vec[i] <- sum(survival$delta * (at.risk.terms * X[,i] - at.risk.X.terms)) / sqrt(n)
  }  
  sum(abs(ret.vec))
}


AFTivScore <- function(beta, data.simu)
{ 
  if (class(data.simu) != "survival.data") {stop("Need to use data generated by simIVMultivarSurvivalData")}
  
  #the score function for the AFT model
  
  survival <- data.simu$survival
  X <- data.simu$X
  Z <- data.simu$Z
  n <- nrow(X)
  nvars <- ncol(X)
  
  #transform to log-time
  survival$t <- log(survival$t)
  
  #store the T_i - bX_i term (error)
  survival$err = survival$t - X %*% beta
  
  #sort according to error size ####observed failure time 
  #data.simu <- data.simu[order(data.simu$X),]  
  order.idx <- order(survival$err)
  survival <- survival[order.idx,] 
  X <- as.matrix(X[order.idx,])
  Z <- as.matrix(Z[order.idx,]) 
  
  #the denominator of the at-risk comparison term  
  at.risk.terms <- n:1
  
  ret.vec <- numeric(nvars)
  for (i in 1:nvars) {
    #the numerator of the at-risk comparison term  
    at.risk.Z.terms <- cumsumRev(Z[,i])
    
    #return the score   
    ret.vec[i] <- sum(survival$delta * (at.risk.terms * Z[,i] - at.risk.Z.terms)) / sqrt(n)
  }  
  ret.vec
}
vAFTivScore <- Vectorize(AFTivScore, vectorize.args = "beta")



absAFTivScore <- function(beta, data.simu)
{ 
  if (class(data.simu) != "survival.data") {stop("Need to use data generated by simIVMultivarSurvivalData")}
  
  #the score function for the AFT model
  
  survival <- data.simu$survival
  X <- data.simu$X
  Z <- data.simu$Z
  n <- nrow(X)
  nvars <- ncol(X)
  
  #transform to log-time
  survival$t <- log(survival$t)
  
  #store the T_i - bX_i term (error)
  survival$err = survival$t - X %*% beta
  
  #sort according to error size ####observed failure time 
  #data.simu <- data.simu[order(data.simu$X),]  
  order.idx <- order(survival$err)
  survival <- survival[order.idx,] 
  X <- as.matrix(X[order.idx,])
  Z <- as.matrix(Z[order.idx,]) 
  
  #the denominator of the at-risk comparison term  
  at.risk.terms <- n:1
  
  ret.vec <- numeric(nvars)
  for (i in 1:nvars) {
    #the numerator of the at-risk comparison term  
    at.risk.Z.terms <- cumsumRev(Z[,i])
    
    #return the score   
    ret.vec[i] <- sum(survival$delta * (at.risk.terms * Z[,i] - at.risk.Z.terms)) / sqrt(n)
  }  
  sum(abs(ret.vec))
}


AFTivScoreSmoothNorm <- function(beta, data.simu, stdev)
{ 
  if (class(data.simu) != "survival.data") {stop("Need to use data generated by simIVMultivarSurvivalData")}
  
  #the score function for the AFT model
  
  survival <- data.simu$survival
  X <- data.simu$X
  Z <- data.simu$Z
  n <- nrow(X)
  nvars <- ncol(X)
  
  #transform to log-time
  survival$t <- log(survival$t)
  
  #store the T_i - bX_i term (error)
  survival$err = survival$t - X %*% beta
  
  #sort according to error size ####observed failure time 
  #data.simu <- data.simu[order(data.simu$X),]  
  order.idx <- order(survival$err)
  survival <- survival[order.idx,] 
  X <- as.matrix(X[order.idx,])
  Z <- as.matrix(Z[order.idx,]) 
  
  #the denominator of the at-risk comparison term  
  survival$at.risk.terms <- unlist(lapply(survival$err, function(x) 
    sum(pnorm(n^(0.26) * (survival$err - x) / stdev))))
  
  #generate empty vector to return eventually
  ret.vec <- numeric(nvars)
  for (i in 1:nvars) {
    #the numerator of the at-risk comparison term 
    at.risk.X.terms <- unlist(lapply(survival$err, function(x) 
      sum(Z[,i] * pnorm(n^(0.26) * (survival$err - x) / stdev))))
    
    #return the score   
    ret.vec[i] <- sum(survival$delta * (survival$at.risk.terms * Z[,i] - at.risk.X.terms)) / sqrt(n)
  }
  
  ret.vec
}



AFTivScoreSmooth <- function(beta, data.simu)
{ 
  if (class(data.simu) != "survival.data") {stop("Need to use data generated by simIVMultivarSurvivalData")}
  
  #the score function for the AFT model
  
  survival <- data.simu$survival
  X <- as.matrix(data.simu$X)
  Z <- as.matrix(data.simu$Z)
  n <- nrow(X)
  nvars <- ncol(X)
  
  #transform to log-time
  survival$t <- log(survival$t)
  
  #store the T_i - bX_i term (error)
  survival$err = survival$t - X %*% beta
  
  #sort according to error size ####observed failure time 
  #data.simu <- data.simu[order(data.simu$X),]  
  order.idx <- order(survival$err)
  survival <- survival[order.idx,] 
  X <- as.matrix(X[order.idx,])
  Z <- as.matrix(Z[order.idx,]) 
  
  #the denominator of the at-risk comparison term  
  survival$at.risk.terms <- unlist(lapply(survival$err, function(x) 
    sum(sigmoid((survival$err - x), 0.5) )))
  
  #generate empty vector to return eventually
  ret.vec <- numeric(nvars)
  for (i in 1:nvars) {
    #the numerator of the at-risk comparison term 
    at.risk.X.terms <- unlist(lapply(survival$err, function(x) 
      sum(Z[,i] * sigmoid((survival$err - x), 0.5))))
    
    #return the score   
    ret.vec[i] <- sum(survival$delta * (survival$at.risk.terms * Z[,i] - at.risk.X.terms)) / sqrt(n)
  }
  
  ret.vec
}


jacobianAFTivSmooth <- function(beta, data.simu, stdev)
{ 
  if (class(data.simu) != "survival.data") {stop("Need to use data generated by simIVMultivarSurvivalData")}
  
  #the score function for the AFT model
  
  survival <- data.simu$survival
  X <- data.simu$X
  Z <- data.simu$Z
  n <- nrow(X)
  nvars <- ncol(X)
  
  #transform to log-time
  survival$t <- log(survival$t)
  
  #store the T_i - bX_i term (error)
  survival$err = survival$t - X %*% beta
  
  #sort according to error size ####observed failure time 
  #data.simu <- data.simu[order(data.simu$X),]  
  order.idx <- order(survival$err)
  survival <- survival[order.idx,] 
  Z <- as.matrix(Z[order.idx,])
  X <- as.matrix(X[order.idx,]) 
  
  
  X_err <- cbind(X, Z, survival$err)
  #generate empty vector to return eventually
  jacobian <- as.matrix(array(0, dim = c(nvars, nvars)))
  for (l in 1:nvars) {
    for (m in 1:nvars) {
      #the numerator of the at-risk comparison term 
      x.sq.terms <- apply(X_err[,c(l, m, nvars + l, nvars + m, 2 * nvars + 1)], 1, function(x) 
        sum((Z[,l] - x[3]) * (X[,m]  - x[2]) * dnorm(n^(0.26) * (x[5] - survival$err) / stdev)))
      
      #return the J(l,m) element   
      jacobian[l, m] <- ((n^(0.26)) / stdev) * sum((survival$delta * (x.sq.terms))) / sqrt(n)
    }
  }
  
  jacobian
}







jacobianAFTivIPCWScoreSmooth <- function(beta, data.simu, stdev)
{ 
  if (class(data.simu) != "survival.data") {stop("Need to use data generated by simIVMultivarSurvivalData")}
  
  #the score function for the AFT model
  
  survival <- data.simu$survival
  X <- data.simu$X
  Z <- data.simu$Z
  n <- nrow(X)
  nvars <- ncol(X)
  
  #generate function G_c() for ICPW
  GC <- genKMCensoringFunc(survival)
  
  #transform to log-time
  survival$t <- log(survival$t)
  
  
  #creates G_c(T_i) term in the IPCW estimating equation
  survival$GCT <- GC(exp(survival$t))
  
  #store the T_i - bX_i term (error)
  survival$err = survival$t - X %*% beta
  
  survival$bX <- X %*% beta
  
  #sort according to error size ####observed failure time 
  #data.simu <- data.simu[order(data.simu$X),]  
  order.idx <- order(survival$err)
  survival <- survival[order.idx,] 
  X <- as.matrix(X[order.idx,])
  Z <- as.matrix(Z[order.idx,]) 
  
  
  #create indicator to set to zero terms where GCT == 0 and 
  #set to 1 so no dividing by zero occurs
  zero.indicator <- 1 * (survival$GCT != 0)
  survival$GCT[which(survival$GCT == 0)] <- 1
  
  
  
  X_err <- cbind(X, Z, survival$err, survival$bX)
  #generate empty vector to return eventually
  jacobian <- as.matrix(array(0, dim = c(nvars, nvars)))
  for (l in 1:nvars) {
    for (m in 1:nvars) {
      #the numerator of the at-risk comparison term 
      x.sq.terms <- apply(X_err[,c(l, m, nvars + l, nvars + m, 2 * nvars + 1, 2 * nvars + 2)], 1, function(x) {
        denom <- GC(exp(x[6] + survival$err))
        ind.2.drop <- 1 * (denom != 0)
        denom[which(denom == 0)] <- 1 #to prevent dividing by 0
        div.denom <- ind.2.drop * 1 / denom
        sum((Z[,l] - x[3]) * (X[,m]  - x[2]) * dnorm(n^(0.26) * (x[5] - survival$err) / stdev) * div.denom)
        })
      
      #return the J(l,m) element   
      jacobian[l, m] <- ((n^(0.26)) / stdev) * sum(((survival$delta / survival$GCT) * (x.sq.terms))) / sqrt(n)
    }
  }
  
  jacobian
  
}













  ###############################################
  ###      Regular Estimating Equations       ###
  ###############################################



# AFTScorePre <- function(beta, survival, X, ZXmat)
# { 
#   #the score function for the AFT model
#   
#   n <- nrow(X)
#   nvars <- ncol(X)
#   
#   #transform to log-time
#   if (is.null(survival$log.t)) {survival$log.t <- log(survival$t)}
#   
#   #store the T_i - bX_i term (error)
#   survival$err = survival$log.t - X %*% beta
#   
#   #sort according to error size ####observed failure time 
#   #data.simu <- data.simu[order(data.simu$X),]  
#   order.idx <- order(survival$err)
#   survival <- survival[order.idx,] 
#   X <- as.matrix(X[order.idx,]) 
#   
#   #the denominator of the at-risk comparison term  
#   at.risk.terms <- n:1 / n
#   
#   ret.vec <- numeric(nvars)
#   for (i in 1:nvars) {
#     #the numerator of the at-risk comparison term  
#     at.risk.Z.terms <- cumsumRev(X[,i]) / n
# 
#     #return the score   
#     ret.vec[i] <- sum(survival$delta * (at.risk.terms * X[,i] - at.risk.Z.terms)) / sqrt(n)
#   }  
#   ret.vec
# }
# 
# attr(AFTScorePre, "name") <- "AFTScorePre"
# 
# 
# 
# AFTivScorePre <- function(beta, survival, X, ZXmat)
# { 
#   #the score function for the AFT model
#   
#   n <- nrow(X)
#   nvars <- ncol(X)
#   
#   #transform to log-time
#   if (is.null(survival$log.t)) {survival$log.t <- log(survival$t)}
#   
#   #store the T_i - bX_i term (error)
#   survival$err = survival$log.t - X %*% beta
#   
#   #sort according to error size ####observed failure time 
#   #data.simu <- data.simu[order(data.simu$X),]  
#   order.idx <- order(survival$err)
#   survival <- survival[order.idx,] 
#   X <- as.matrix(X[order.idx,])
#   ZXmat <- as.matrix(ZXmat[order.idx,]) 
#   
#   #the denominator of the at-risk comparison term  
#   at.risk.terms <- n:1 / n
#   
#   ret.vec <- numeric(nvars)
#   for (i in 1:nvars) {
#     #the numerator of the at-risk comparison term  
#     at.risk.Z.terms <- cumsumRev(ZXmat[,i]) / n
#     
#     #return the score   
#     ret.vec[i] <- sum(survival$delta * (at.risk.terms * ZXmat[,i] - at.risk.Z.terms)) / sqrt(n)
#   }  
#   ret.vec
# }
# attr(AFTivScorePre, "name") <- "AFTivScorePre"
# 
# AFTivIPCWScorePre <- function(beta, survival, X, ZXmat, GC)
# { 
#   
#   #the score function for the AFT model with instrumental variables and Inverse Probability Censoring Weighting
#   n <- nrow(X)
#   nvars <- ncol(X)
#   
#   #transform to log-time
#   if (is.null(survival$log.t)) {survival$log.t <- log(survival$t)}
#   
#   #creates G_c(T_i) term in the IPCW estimating equation
#   survival$GCT <- GC(survival$t)
#   
#   #store the T_i - bX_i term (error)
#   survival$err = survival$log.t - X %*% beta
#   
#   survival$bX <- X %*% beta
#   
#   #sort according to error size ####observed failure time 
#   #data.simu <- data.simu[order(data.simu$X),]  
#   order.idx <- order(survival$err)
#   survival <- survival[order.idx,] 
#   X <- as.matrix(X[order.idx,])
#   #Z <- as.matrix(Z[order.idx,])
#   ZXmat <- as.matrix(ZXmat[order.idx,])
#   
#   #create indicator to set to zero terms where GCT == 0 and 
#   #set to 1 so no dividing by zero occurs
#   zero.indicator <- 1 * (survival$GCT != 0)
#   survival$GCT[which(survival$GCT == 0)] <- 1
#   
#   #first col is as.risk.terms, remaining are at.risk.Z.terms
#   at.risk.mat <- genIPCWNumDenomMultivar(survival, ZXmat, GC) / n
#   
#   #generate empty vector to return eventually
#   ret.vec <- numeric(nvars)
#   for (i in 1:nvars) {
#     #return the score   
#     ret.vec[i] <- sum(zero.indicator * (survival$delta / survival$GCT) * (at.risk.mat[,1] * ZXmat[,i] - at.risk.mat[,i + 1])) / sqrt(n)
#   }
#   
#   ret.vec
#   
# }
# 
# attr(AFTivIPCWScorePre, "name") <- "AFTivIPCWScorePre"
# 
# 
# 
# genIPCWNumDenomMultivar2 <- cmpfun(function(bX, Z, err, GC.func){
#   #dat is a data.frame
#   #GC.func is a function
#   num.vars <- ncol(Z)
#   Z <- as.matrix(Z)
#   n.obs <- nrow(bX)
#   
#   at.risk.list <- lapply(1:n.obs, function(i) {
#     err.i <- err[i]
#     ind.zero <- F
#     ipcw <- GC.func(exp(bX[i:n.obs] + err.i))
#     
#     if (all(ipcw == 0)){
#       len.i <- length(ipcw)
#       print(len.i)
#       ipcw <- rep(0.01, len.i)
#       ind.zero <- T
#     }
#     #ipcw[which(ipcw == 0)] <- min(ipcw[which(ipcw != 0)]) / 2
#     ret.vec <- array(0, dim=num.vars+1)
#     if (!ind.zero){
#       ret.vec[1] <- sum(1 / ipcw)
#       ret.vec[2:(num.vars+1)] <- if(i != n.obs ){ 
#         colSums(Z[i:n.obs,] / ipcw)
#       } else {
#         Z[n.obs,] / ipcw
#       }
#         
#       return (ret.vec)
#     } else {
#       return (ret.vec)
#     }
#   })
#   do.call(rbind, at.risk.list)
# })
# 
# AFTivIPCWScorePre2 <- function(beta, survival, X, ZXmat, GC)
# { 
#   
#   #the score function for the AFT model with instrumental variables and Inverse Probability Censoring Weighting
#   n <- nrow(X)
#   nvars <- ncol(X)
#   
#   #transform to log-time
#   if (is.null(survival$log.t)) {
#     log.t <- log(survival$t)
#   } else {log.t <- survival$log.t}
#   
#   t <- survival$t
#   
#   #creates G_c(T_i) term in the IPCW estimating equation
#   GCT <- GC(t)
#   
#   bX <- X %*% beta
#   
#   #store the T_i - bX_i term (error)
#   err = log.t - bX
#   
#   delta <- survival$delta
#   
#   
#   
#   #sort according to error size ####observed failure time 
#   #data.simu <- data.simu[order(data.simu$X),]  
#   order.idx <- order(err)
#   #survival <- survival[order.idx,] 
#   GCT <- GCT[order.idx]
#   err <- err[order.idx]
#   delta <- delta[order.idx]
#   X <- as.matrix(X[order.idx,])
#   bX <- as.matrix(bX[order.idx,])
#   ZXmat <- as.matrix(ZXmat[order.idx,])
#   
#   
#   
#   #create indicator to set to zero terms where GCT == 0 and 
#   #set to 1 so no dividing by zero occurs
#   zero.indicator <- 1 * (GCT != 0)
#   GCT[which(GCT == 0)] <- 1
#   
#   #first col is as.risk.terms, remaining are at.risk.Z.terms
#   at.risk.mat <- genIPCWNumDenomMultivar2(bX = bX, err = err, Z = ZXmat, GC)
# 
#   ret.vec <- colSums(zero.indicator * (delta / GCT) * (at.risk.mat[,1] * ZXmat - at.risk.mat[,2:(nvars + 1)])) / sqrt(n)
#   
# }
# 
#  ###############################################
#  ###     Smoothed Estimating Equations       ###
#  ###############################################
# 
# AFTivScoreSmoothPre <- function(beta, survival, X, ZXmat, tau = 1e-3)
# { 
#   
#   #the score function for the AFT model
#   
#   n <- nrow(X)
#   nvars <- ncol(X)
#   
#   #transform to log-time
#   if (is.null(survival$log.t)) {
#     log.t <- log(survival$t)
#   } else (log.t <- survival$log.t)
#   
#   delta <- survival$delta
#   
#   if (n < 10000) {
#     #store the T_i - bX_i term (error)
#     err = as.vector(log.t - X %*% beta)
#     
#     #the denominator of the at-risk comparison term  
#     out.sig <- sigmoid(outer(err, err, "-"), tau)
#     at.risk.terms <- colSums(out.sig)
#   } else {
#     #store the T_i - bX_i term (error)
#     err = log.t - X %*% beta
#     
#     #the denominator of the at-risk comparison term  
#     at.risk.terms <- unlist(lapply(err, function(x) sum(sigmoid((err - x), tau) )))
#   }
#   #generate empty vector to return eventually
#   #   ret.vec <- numeric(nvars)
#   #   for (i in 1:nvars) {
#   #     #the numerator of the at-risk comparison term 
#   #     at.risk.X.terms <- unlist(lapply(survival$err, function(x) 
#   #       sum(ZXmat[,i] * sigmoid((survival$err - x), tau))))
#   #     
#   #     #return the score   
#   #     ret.vec[i] <- sum(survival$delta * (survival$at.risk.terms * ZXmat[,i] - at.risk.X.terms)) / sqrt(n)
#   #   }
#   
#   
#   if (n < 10000) {
#     ret.vec <- apply(ZXmat, 2, function(ZX) {
#       at.risk.X.terms <- colSums(out.sig * ZX)
#       sum(delta * (at.risk.terms * ZX - at.risk.X.terms)) / sqrt(n)
#     })
#   } else {
#     ret.vec <- apply(ZXmat, 2, function(ZX) {
#       at.risk.X.terms <- unlist(lapply(err, function(x) sum(ZX * sigmoid((err - x), tau))))
#       sum(delta * (at.risk.terms * ZX - at.risk.X.terms)) / sqrt(n)
#     })
#   }
#   
#   ret.vec
# }
# 
# attr(AFTivScoreSmoothPre, "name") <- "AFTivScoreSmoothPre"
# 
# genIPCWNumDenomMultivar <- cmpfun(function(dat, Z, GC.func){
#   #dat is a data.frame
#   #GC.func is a function
#   num.vars <- ncol(Z)
#   num <- denom <- array(0, dim = c(nrow(dat),1))
#   at.risk.list <- lapply(1:nrow(dat), function(i) {
#     err.i <- dat$err[i]
#     ind.zero <- F
#     ipcw <- GC.func(exp(dat$bX[i:nrow(dat)] + err.i))
#     if (all(ipcw == 0)){
#       ipcw <- rep(0.01, length(ipcw))
#       ind.zero <- T
#     }
#     ipcw[which(ipcw == 0)] <- min(ipcw[which(ipcw != 0)]) / 2
#     ret.vec <- array(0, dim=num.vars+1)
#     if (!ind.zero){
#       ret.vec[1] <- sum(1 / ipcw)
#       for (j in 1:num.vars) {
#         ret.vec[j+1] <- sum(Z[i:nrow(dat),j] / ipcw)
#       }
#       return (ret.vec)
#     } else {
#       return (ret.vec)
#     }
#   })
#   do.call(rbind, at.risk.list)
# })
# 
# AFTScoreSmoothPre <- function(beta, survival, X, ZXmat, tau = 1e-3)
# { 
#   
#   #the score function for the AFT model
#   
#   n <- nrow(X)
#   nvars <- ncol(X)
#   
#   #transform to log-time
#   if (is.null(survival$log.t)) { 
#     log.t <- log(survival$t)
#   } else (log.t <- survival$log.t)
#   
#   delta <- survival$delta
#   
#   if (n < 10000) {
#     #store the T_i - bX_i term (error)
#     err = as.vector(log.t - X %*% beta)
#     
#     #the denominator of the at-risk comparison term  
#     out.sig <- sigmoid(outer(err, err, "-"), tau)
#     at.risk.terms <- colSums(out.sig) / n
#   } else {
#     #store the T_i - bX_i term (error)
#     err = log.t - X %*% beta
#     
#     #the denominator of the at-risk comparison term  
#     at.risk.terms <- unlist(lapply(err, function(x) sum(sigmoid((err - x), tau) ))) / n
#   }
#   # #generate empty vector to return eventually
#   #   ret.vec <- numeric(nvars)
#   #   for (i in 1:nvars) {
#   #     #the numerator of the at-risk comparison term 
#   #     at.risk.X.terms <- unlist(lapply(err, function(x) sum(X[,i] * sigmoid((err - x), tau))))
#   #     
#   #     #return the score   
#   #     ret.vec[i] <- sum(delta * (at.risk.terms * X[,i] - at.risk.X.terms)) / sqrt(n)
#   #   }
#   if (n < 10000) {
#     ret.vec <- apply(X, 2, function(xi) {
#       at.risk.X.terms <- colSums(out.sig * xi) / n
#       sum(delta * (at.risk.terms * xi - at.risk.X.terms)) / sqrt(n)
#     })
#   } else {
#     ret.vec <- apply(X, 2, function(xi) {
#       at.risk.X.terms <- unlist(lapply(err, function(x) sum(xi * sigmoid((err - x), tau)))) / n
#       sum(delta * (at.risk.terms * xi - at.risk.X.terms)) / sqrt(n)
#     })
#   }
#   
#   ret.vec
# }
# 
# attr(AFTScoreSmoothPre, "name") <- "AFTScoreSmoothPre"


