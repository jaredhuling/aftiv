

lsBootstrap <- function(beta, esteqn, B, nobs, GC = NULL, VarFunc = NULL, n.riskFunc = NULL, ...) 
{
  p <- length(beta)

  
  is.ipcw <- attr(esteqn, "name") == "evalAFTivIPCWScorePrec" | 
    attr(esteqn, "name") == "AFTivIPCWScorePre"
  
  Mb <- matrix(rexp(nobs * B, rate = 1), ncol = B)
  if (is.ipcw) 
  {
    #Un2 <- t(apply( Mb, 2, function(x) esteqn(beta = beta, multiplier.wts = x, GC = GC, ...) ))
    
    #varG <- VarFunc(time.t)
    #n.risk <- n.riskFunc(time.t)
    #idx.nnan <- !is.nan(varG)
    #varG <- (sum(delta[idx.nnan]*varG[idx.nnan]^2) ) * nobs #* sum(delta[idx.nnan]==1) #* nobs #sum(delta[idx.nnan]==1) #) #* sum(delta[idx.nnan]==1) # #
    ##* sqrt(nobs) # * nobs # * nobs
    ##print(varG)
    ##varG <- 0
    #print("varg1")
    #print(varG)
    
    #biases.sum <- rep(NA, B)
    #for (i in 1:B) {
    #  samp.idx <- sample.int(length(time.t), length(time.t), replace = TRUE)
    #  sfb <- survfit(Surv(time.t[samp.idx], delta[samp.idx]) ~ 1)
    #  serv <- summary(sfb)$surv 
    #  
    #  dfsurv = data.frame(c(0, summary(sfb)$time), c(serv[1], serv)); 
    #  tmpFuncSurv <- stepfun(dfsurv[,1], c(dfsurv[1,2],dfsurv[,2]), f=0);
    #  biases <- tmpFuncSurv(time.t[samp.idx][which(delta[samp.idx]==1)]) #- 1+pexp(obs.time[which(status==1)], 1)
    #  #biases <- tmpFuncSurv(obs.time)
    #  #biases[which(is.nan(biases))] <- cur.vars[which(is.nan(biases))-5]
    #  biases.sum[i] <- sum(biases)
    #  #rmeans[i] <- summary(sfb, rmean=TRUE)$table[["*rmean"]]
    #}
    
    #varG <- var(biases.sum)
    #print("varg2")
    #print(varG)
    
    Un2 <- array(NA, dim = c(B, p) )
    if (attr(esteqn, "name") == "evalAFTivIPCWScorePrec") 
    {
      data <- list(...)[["data.simu"]]
      
      for (i in 1:B) 
      {
        samp.idx <- sample.int(nobs, nobs, replace = TRUE)
        GC.boot  <- genKMCensoringFunc(data[samp.idx,])
        Un2[i,]  <- esteqn(beta = beta, GC = GC.boot, data.simu = data[samp.idx,])
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
        Un2[i,]  <- esteqn(beta       = beta, 
                           GC         = GC.boot, 
                           survival   = survival[samp.idx,],
                           X          = X[samp.idx,], 
                           ZXmat      = ZXmat[samp.idx,], 
                           conf.x.loc = conf.x.loc,
                           multiplier.wts = Mb[,i])
      }
    }
    
  } else 
  {
    Un2 <- t(apply( Mb, 2, function(x) esteqn(beta = beta, multiplier.wts = x, ...) ))
    if (p == 1) 
    {
      Un2 <- t(Un2)
    }
  }
  V <- var(Un2)
 
  if (p == 1)
  {
    rt <- as.numeric(sqrt(V))
    if (is.ipcw) 
    {
      cat("V:", V, "\n")
      Zb <- (1/rt) * (matrix(rnorm(p * B), ncol = p))
    } else 
    {
      Zb <- (1/rt) * (matrix(rnorm(p * B), ncol = p))
    }
    #Zb <- rt * (matrix(2 * rnorm(p * B) - 1, ncol = p))
  } else 
  {
    if (is.ipcw) 
    {
      # eeg <- eigen(V)
      # rt  <- solve(eeg$vectors %*% diag(sqrt(eeg$values) ) %*% solve(eeg$vectors))
      # Zb  <- t(rt %*% t(matrix(rnorm(p * B), ncol = p)))
      Zb <- matrix(rnorm(p * B), ncol = p)
    } else  
    {
      Zb <- matrix(rnorm(p * B), ncol = p)
    }
  }
  

  #Zb <- matrix(rnorm(p * B, sd = nobs^0.35), ncol = p)
  
  
  
  
  if (is.ipcw) 
  {
    #time.t <- list(...)[["data.simu"]]
    #delta  <- data$delta
    #time.t <- data$t
    Un1    <- t(apply( Zb, 1, function(x) esteqn(beta = beta + x / sqrt(nobs), GC = GC, ...) ))
  } else 
  {
    Un1    <- t(apply( Zb, 1, function(x) esteqn(beta = beta + x / sqrt(nobs), ...) ))
  }
  
  if (p == 1) 
  {
    Un1 <- t(Un1)
  }
  AA     <- lm.fit(y = Un1, x = Zb)$coefficients
  AA.inv <- solve(AA)
  
  
  variance <- AA.inv %*% V %*% t(AA.inv)
    
  
  se.hat <- sqrt(diag(variance)) / sqrt(nobs)
  
  #print(attr(esteqn, "name"))
  #cat("The est: ", beta, "\n")
  #print("The se: ")
  #print(se.hat)
  
  list(se.hat = se.hat, 
       var    = variance, 
       A      = AA, 
       V      = V)
}

svBootstrap <- function(beta, esteqn, B, nobs, 
                        GC = NULL, VarFunc = NULL, 
                        var.method = c("binomial", "normal"), ...) 
{
  p  <- length(beta)
  Mb <- matrix(rexp(nobs * B, rate = 1), ncol = B)
  var.method <- match.arg(var.method)
  
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
                           conf.x.loc = conf.x.loc,
                           multiplier.wts = Mb[,i])
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

    Zb <- (1/rt) * (matrix(rnorm(p * B), ncol = p))
    #Zb <- rt * (matrix(2 * rnorm(p * B) - 1, ncol = p))
  } else 
  {
    if (var.method == "binomial")
    {
      eeg <- eigen(V)
      rt  <- solve(eeg$vectors %*% diag(sqrt(eeg$values) ) %*% solve(eeg$vectors))
      if (is.ipcw) 
      {
        #Zb <- t(rt %*% t(3 * matrix(rbinom(p * B, 1, 0.5), ncol = p) - 1.5)) / ( log(log(nobs)) )
        Zb <- t(rt %*% t(2 * matrix(rbinom(p * B, 1, 0.5), ncol = p) - 1))
      } else 
      {
        Zb <- t(rt %*% t(2 * matrix(rbinom(p * B, 1, 0.5), ncol = p) - 1))
      }
    } else
    {
      Zb <- mvrnorm(B, rep(0, p), Sigma = solve(V))
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
  
  list(se.hat = se.hat, var = sigma.hat, V = V)  
}

bootstrapVarUnivar <- function(est, 
                               est.eqn,
                               data, 
                               method     = c("ls", "sv", "full.bootstrap"), 
                               B          = 1000L, 
                               GC         = NULL, 
                               VarFunc    = NULL, 
                               n.riskFunc = NULL) 
{
  # takes a fitted aft.fit object and computes estimated variance using
  # a bootstrap approach. supports 2 fast multiplier boostrap techniques
  # in addition to traditional boostrap estimate
  
  method <- match.arg(method)
  
  #stopifnot(class(data) == "survival.data")  
  
  #conf.x.loc <- match(fitted.obj$confounded.x.names, colnames(data$X))
  #ZXmat <- data$X
  
  if (is.null(nrow(data$X))) 
  {
    nobs <- length(data$X)
  } else 
  {
    nobs <- nrow(data$X)
  }
    
  
  if (method == "ls") 
  {
    if (attr(est.eqn, "name") == "evalAFTivIPCWScorePrec") 
    {
      bs <- lsBootstrap(beta = est, 
                        esteqn = est.eqn, 
                        B = B, 
                        nobs = nobs, 
                        GC = GC, 
                        VarFunc = VarFunc, 
                        n.riskFunc = n.riskFunc, 
                        data.simu = data)
    } else {
      bs <- lsBootstrap(beta = est, 
                        esteqn = est.eqn, 
                        B = B, 
                        nobs = nobs, 
                        GC = NULL, 
                        VarFunc = NULL, 
                        n.riskFunc = NULL, 
                        data.simu = data)
    }
    
  } else if (method == "sv") 
  {
    if (attr(est.eqn, "name") == "evalAFTivIPCWScorePrec") 
    {
      bs <- svBootstrap(beta = est, 
                        esteqn = est.eqn, 
                        B = B, 
                        nobs = nobs, 
                        data.simu = data, 
                        GC = GC, 
                        VarFunc = VarFunc)
    } else 
    {
      bs <- svBootstrap(beta = est, 
                        esteqn = est.eqn, 
                        B = B, 
                        nobs = nobs, 
                        data.simu = data)
    }
    
  } else {
    stop("method not supported yet")
  }
  
  bs
}

bootstrapVar <- function(fitted.obj, data, 
                         method = c("ls", "sv", "full.bootstrap"), 
                         B = 1000L) 
{
  # takes a fitted aft.fit object and computes estimated variance using
  # a bootstrap approach. supports 2 fast multiplier boostrap techniques
  # in addition to traditional boostrap estimate
  
  method <- match.arg(method)
  stopifnot(class(data) == "survival.data")
  stopifnot(class(fitted.obj) == "aft.fit")
  
  
  
  conf.x.loc <- match(fitted.obj$confounded.x.names, colnames(data$X))
  ZXmat <- data$X
  
  if (fitted.obj$final.fit) 
  {
    est.eqn <- fitted.obj$est.eqn.sm
  } else 
  {
    est.eqn <- fitted.obj$est.eqn
  }

  if (attr(est.eqn, "name") == "AFT2SLSScorePre" | attr(est.eqn, "name") == "AFT2SLSScoreSmoothPre"
      | attr(est.eqn, "name") == "AFT2SLSScoreAllPre" | attr(est.eqn, "name") == "AFT2SLSScoreSmoothAllPre") 
  {
    #if 2SLS is used, replace Z with Xhat 
    Z <- as.matrix(data$Z)
    if (attr(est.eqn, "name") == "AFT2SLSScoreAllPre" | attr(est.eqn, "name") == "AFT2SLSScoreSmoothAllPre") 
    {
      for (v in 1:ncol(data$X)) 
      {
        data$X[,v] <- lm(data$X[,v] ~ Z)$fitted.values
      }
    } else 
    {
      for (v in 1:length(conf.x.loc)) 
      {
        dat.lm <- data.frame(response = data$X[,conf.x.loc[v]], predictor = Z)
        Xhat   <- lm(response ~ predictor, data = dat.lm)$fitted.values
        data$X[,conf.x.loc[v]] <- Xhat
      }
    }
    if (attr(est.eqn, "name") == "AFT2SLSScoreSmoothPre" | attr(est.eqn, "name") == "AFT2SLSScoreSmoothAllPre") 
    {
      est.eqn               <- AFTScoreSmoothPre
      attr(est.eqn, "name") <- "AFTScoreSmoothPre"
    } else 
    {
      est.eqn               <- AFTScorePre
      attr(est.eqn, "name") <- "AFTScorePre"
    }
    

    
  } else 
  {
    ZXmat[,conf.x.loc] <- data$Z
  }
  
  if (method == "ls") 
  {
    if (attr(est.eqn, "name") == "AFTivIPCWScorePre") 
    {
      #generate function G_c() for ICPW 
      GC <- genKMCensoringFunc(data$survival)
      bs <- lsBootstrap(beta = fitted.obj$par, 
                        esteqn = est.eqn, 
                        B = B, 
                        nobs = fitted.obj$nobs, 
                        GC = GC, 
                        survival = data$survival,
                        X = data$X, 
                        ZXmat = ZXmat, 
                        conf.x.loc = conf.x.loc)
    } else 
    {
      bs <- lsBootstrap(beta = fitted.obj$par, 
                        esteqn = est.eqn, 
                        B = B, 
                        nobs = fitted.obj$nobs, 
                        survival = data$survival,
                        X = data$X, 
                        ZXmat = ZXmat)
    }

  } else if (method == "sv") 
  {
    if (attr(est.eqn, "name") == "AFTivIPCWScorePre") 
    {
      #generate function G_c() for ICPW 
      GC <- genKMCensoringFunc(data$survival)
      bs <- svBootstrap(beta = fitted.obj$par, 
                        esteqn = est.eqn, 
                        B = B, 
                        nobs = fitted.obj$nobs, 
                        GC = GC, 
                        survival = data$survival,
                        X = data$X, 
                        ZXmat = ZXmat, 
                        conf.x.loc = conf.x.loc)
    } else 
    {
      bs <- svBootstrap(beta = fitted.obj$par, 
                        esteqn = est.eqn, 
                        B = B, 
                        nobs = fitted.obj$nobs, 
                        survival = data$survival,
                        X = data$X, 
                        ZXmat = ZXmat)
    }
  } else {
    stop("method not supported yet")
  }
  
  bs
}



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



bootstrapVarEstEqn <- function(est, data, 
                               est.eqn, 
                               method = c("ls", "sv", "full.bootstrap"), 
                               B = 1000L) 
{
  # takes a fitted aft.fit object and computes estimated variance using
  # a bootstrap approach. supports 2 fast multiplier boostrap techniques
  # in addition to traditional boostrap estimate
  
  method <- match.arg(method)
  stopifnot(class(data) == "survival.data")
  stopifnot(class(fitted.obj) == "aft.fit")
  
  
  
  conf.x.loc <- match(fitted.obj$confounded.x.names, colnames(data$X))
  ZXmat      <- data$X

  
  if (attr(est.eqn, "name") == "AFT2SLSScorePre" | attr(est.eqn, "name") == "AFT2SLSScoreSmoothPre"
      | attr(est.eqn, "name") == "AFT2SLSScoreAllPre" | attr(est.eqn, "name") == "AFT2SLSScoreSmoothAllPre") 
  {
    #if 2SLS is used, replace Z with Xhat 
    Z <- as.matrix(data$Z)
    if (attr(est.eqn, "name") == "AFT2SLSScoreAllPre" | attr(est.eqn, "name") == "AFT2SLSScoreSmoothAllPre") 
    {
      for (v in 1:ncol(data$X)) 
      {
        data$X[,v] <- lm(data$X[,v] ~ Z)$fitted.values
      }
    } else 
    {
      for (v in 1:length(conf.x.loc)) 
      {
        dat.lm <- data.frame(response = data$X[,conf.x.loc[v]], predictor = Z)
        Xhat   <- lm(response ~ predictor, data = dat.lm)$fitted.values
        data$X[,conf.x.loc[v]] <- Xhat
      }
    }
    if (attr(est.eqn, "name") == "AFT2SLSScoreSmoothPre" | attr(est.eqn, "name") == "AFT2SLSScoreSmoothAllPre") 
    {
      est.eqn <- AFTScoreSmoothPre
      attr(est.eqn, "name") <- "AFTScoreSmoothPre"
    } else 
    {
      est.eqn <- AFTScorePre
      attr(est.eqn, "name") <- "AFTScorePre"
    }
    
    
    
  } else 
  {
    ZXmat[,conf.x.loc] <- data$Z
  }
  
  if (method == "ls") 
  {
    if (attr(est.eqn, "name") == "AFTivIPCWScorePre") 
    {
      #generate function G_c() for ICPW 
      GC <- genKMCensoringFunc(data$survival)
      bs <- lsBootstrap(beta = fitted.obj$par, 
                        esteqn = est.eqn, 
                        B = B,
                        nobs = fitted.obj$nobs, 
                        GC = GC, 
                        survival = data$survival,
                        X = data$X, 
                        ZXmat = ZXmat)
    } else 
    {
      bs <- lsBootstrap(beta = fitted.obj$par, 
                        esteqn = est.eqn, 
                        B = B, 
                        nobs = fitted.obj$nobs, 
                        survival = data$survival,
                        X = data$X, 
                        ZXmat = ZXmat)
    }
    
  } else if (method == "sv") 
  {
    if (attr(est.eqn, "name") == "AFTivIPCWScorePre") 
    {
      #generate function G_c() for ICPW 
      GC <- genKMCensoringFunc(data$survival)
      bs <- svBootstrap(beta = fitted.obj$par, 
                        esteqn = est.eqn, 
                        B = B, 
                        nobs = fitted.obj$nobs, 
                        GC = GC, 
                        survival = data$survival,
                        X = data$X, 
                        ZXmat = ZXmat)
    } else 
    {
      bs <- svBootstrap(beta = fitted.obj$par, 
                        esteqn = est.eqn, 
                        B = B, 
                        nobs = fitted.obj$nobs, 
                        survival = data$survival,
                        X = data$X, 
                        ZXmat = ZXmat)
    }
  } else 
  {
    stop("method not supported yet")
  }
  
  bs
}




bootstrapCI <- function(aft.fit, data, n.bootstraps = 999, percent, packages = NULL, func.names) 
{
  # brute force bootstrap approach
  stopifnot(class(data) == "survival.data")
  stopifnot(class(aft.fit) == "aft.fit")
  
  n.bootstraps <- as.integer(n.bootstraps)
  call         <- aft.fit$call
  #initialize with coefficient estimates from initial
  #fit to speed up fitting process
  
  newmaxit <- 32
  call[[match("maxit", names(call))]] <- as.name("newmaxit")
  coefficients <- aft.fit$par
  if (!is.null(call[[match("init.par", names(call))]])) 
  {
    call[[match("init.par", names(call))]] <- as.name("coefficients")
  } else 
  {
    call$init.par <- as.name("coefficients")
  }
  boot.estimates <- foreach(i = 1:n.bootstraps, .packages = packages, .combine = cbind, 
                            .export = c(func.names, "newmaxit", "coefficients")) %dopar% {
                              print(sprintf("Bootstrap %g / %g", i, n.bootstraps))
                              set.seed(i + 123)
                              #resample data
                              samp.idx           <- sample(1:nrow(data$X), nrow(data$X), replace = T)
                              data.samp          <- data
                              data.samp$X        <- data.samp$X[samp.idx,]
                              data.samp$Z        <- data.samp$Z[samp.idx]
                              data.samp$survival <- data.samp$survival[samp.idx,]
                              call[[match("data", names(call))]] <- as.name("data.samp")
                              #fit coefficients with resampled data
                              eval(call)$par
                            }
  
  ret.mat        <- data.frame(array(0, dim = c(length(coefficients), 5)))
  names(ret.mat) <- c("Name", "Beta", "ci.lower", "ci.upper", "se")
  
  for (i in 1:length(coefficients)) 
  {
    boot.est     <- boot.estimates[i,]
    alpha.o.2    <- (1 - percent) / 2
    basic.quants <- quantile(boot.est, probs = c(1 - alpha.o.2, alpha.o.2))
    basic.CIs    <- 2 * coefficients[i] - basic.quants
    se.hat       <- sd(boot.est)
    ret.mat[i,2:ncol(ret.mat)] <- c(coefficients[i], basic.CIs, se.hat)
  }
  ret.mat$Name <- names(coefficients)
  ret          <- list(results = ret.mat, boot.est = boot.estimates)
  class(ret)   <- "bootstrap.results"
  ret
}

print.bootstrap.results <- function(bsr, true.beta = NULL) 
{
  res <- bsr$results
  if (!is.null(true.beta)) 
  {
    res$True.Beta <- true.beta
    res <- res[,c(1,ncol(res),2:(ncol(res) - 1))]
  }
  res[,2:ncol(res)] <- round(res[,2:ncol(res)], 4)
  print(res)
}


