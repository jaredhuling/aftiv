

lsBootstrap <- function(beta, esteqn, B, nobs, ...) {
  p <- length(beta)
  Zb <- matrix(rnorm(p * B), ncol = p)
  Un1 <- t(apply( Zb, 1, function(x) esteqn(beta = beta + x / sqrt(nobs), ...) ))
  AA <- lm.fit(y = Un1, x = Zb)$coefficients

  Mb <- matrix(rexp(nobs * B, rate = 1), ncol = B)
  
  Un2 <- t(apply( Mb, 2, function(x) esteqn(beta = beta, multiplier.wts = x, ...) ))
  
  V <- var(Un2)
  
  AA.inv <- solve(AA)
  
  list(var = AA.inv %*% V %*% t(AA.inv), A = AA, V = V)
}

svBootstrap <- function(beta, esteqn, B, nobs, ...) {
  p <- length(beta)
  
  Mb <- matrix(rexp(nobs * B, rate = 1), ncol = B)
  
  Un1 <- t(apply( Mb, 2, function(x) esteqn(beta = beta, multiplier.wts = x, ...) ))
  
  V <- var(Un1)
  
  Zb <- mvrnorm(n = B, rep(0, p), Sigma = solve(V))
  
  Un2 <- t(apply( Zb, 1, function(x) esteqn(beta = beta + x / sqrt(nobs), ...) ))
  
  sigma.hat <- var(Un2)
  
  list(var = solve(sigma.hat), V = V)  
}


bootstrapVar <- function(fitted.obj, data, est.eqn) {
  
  stopifnot(class(data) == "survival.data")
  stopifnot(class(aft.fit) == "aft.fit")
  
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
    
    if (attr(est.eqn, "name") == "AFTivIPCWScorePre") {
      #generate function G_c() for ICPW 
      GC <- genKMCensoringFunc(data$survival)
      
    }
    
  } else {
    ZXmat[,conf.x.loc] <- data$Z
  }
}
