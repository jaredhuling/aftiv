

AFTScorePre.cpp <- function(beta, survival, X, ZXmat) 
{
  
  if (is.null(survival$log.t)) 
  {
    log.t <- log(survival$t)
  } else 
  {
    log.t <- survival$log.t
  }
  
  err       <- log.t - X %*% beta
  order.idx <- order(err)
  X         <- as.matrix(X[order.idx,]) 
  delta     <- survival$delta
  delta     <- delta[order.idx]
  
  .Call("AFTScorePre", XX = X, delta_vec = delta, PACKAGE = "aftiv")
}

AFTivScorePre.cpp <- function(beta, survival, X, ZXmat) 
{
  
  if (is.null(survival$log.t)) 
  {
    log.t <- log(survival$t)
  } else 
  {
    log.t <- survival$log.t
  }
  
  err       <- log.t - X %*% beta
  order.idx <- order(err)
  ZXmat     <- as.matrix(ZXmat[order.idx,]) 
  delta     <- survival$delta
  delta     <- delta[order.idx]
  
  .Call("AFTScorePre", XX = ZXmat, delta_vec = delta, PACKAGE = "aftiv")
}

