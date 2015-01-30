

lsBootstrap <- function(beta, esteqn, B, nobs, ...) {
  p <- length(beta)
  Zb <- matrix(rnorm(p * B), ncol = p)
  Un1 <- t(apply( Zb, 1, function(x) esteqn(beta = beta + x / sqrt(nobs), ...) ))
  AA <- lm.fit(y = Un1, x = Zb)$coefficients

  Mb <- matrix(rexp(nobs * B, rate = 1), ncol = B)
  
  Un2 <- t(apply( Mb, 2, function(x) esteqn(beta = beta, multiplier.wts = x, ...) ))
  
  V <- var(Un2)
  
  AA.inv <- solve(AA)
  
  variance <- AA.inv %*% V %*% t(AA.inv)
  
  se.hat <- sqrt(diag(variance)) / sqrt(nobs)
  
  list(se.hat = se.hat, var = variance, A = AA, V = V)
}

svBootstrap <- function(beta, esteqn, B, nobs, ...) {
  p <- length(beta)
  
  Mb <- matrix(rexp(nobs * B, rate = 1), ncol = B)
  
  Un1 <- t(apply( Mb, 2, function(x) esteqn(beta = beta, multiplier.wts = x, ...) ))
  
  V <- var(Un1)
  
  Zb <- mvrnorm(n = B, rep(0, p), Sigma = solve(V))
  
  Un2 <- t(apply( Zb, 1, function(x) esteqn(beta = beta + x / sqrt(nobs), ...) ))
  
  sigma.hat <- solve(var(Un2))
  
  se.hat <- sqrt(diag(sigma.hat)) / sqrt(nobs)
  
  list(se.hat = se.hat, var = sigma.hat, V = V)  
}


bootstrapVar <- function(fitted.obj, data, 
                         method = c("ls", "sv", "full.bootstrap"), 
                         B = 1000L) {
  # takes a fitted aft.fit object and computes estimated variance using
  # a bootstrap approach. supports 2 fast multiplier boostrap techniques
  # in addition to traditional boostrap estimate
  
  method <- match.arg(method)
  stopifnot(class(data) == "survival.data")
  stopifnot(class(fitted.obj) == "aft.fit")
  
  conf.x.loc <- match(fitted.obj$confounded.x.names, colnames(data$X))
  ZXmat <- data$X
  
  if (fitted.obj$final.fit) {
    est.eqn <- fitted.obj$est.eqn.sm
  } else {
    est.eqn <- fitted.obj$est.eqn
  }
  
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
    
    if (attr(est.eqn, "name") == "AFTivIPCWScorePre") {
      #generate function G_c() for ICPW 
      GC <- genKMCensoringFunc(data$survival)
      
    }
    
  } else {
    ZXmat[,conf.x.loc] <- data$Z
  }
  
  if (method == "ls") {
    bs <- lsBootstrap(beta = fitted.obj$par, esteqn = est.eqn, 
                      B = B, nobs = fitted.obj$nobs)
  } else if (method == "sv") {
    bs <- svBootstrap(beta = fitted.obj$par, esteqn = est.eqn, 
                      B = B, nobs = fitted.obj$nobs)
  } else {
    stop("method not supported yet")
  }
  
  bs
}
