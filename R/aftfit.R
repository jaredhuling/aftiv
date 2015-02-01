

aftfit <- function(formula, 
                   data, 
                   instrument, 
                   confounded.x.names, 
                   method = c("AFT", "AFT-IV", "AFT-2SLS", "AFT-IPCW"),
                   smoothed = FALSE,
                   weights = NULL, 
                   bootstrap = FALSE,
                   boot.method = c("ls", "sv", "full.bootstrap"),
                   B = 1000L,
                   na.action, 
                   init = NULL, 
                   return.data=FALSE,
                   tol = 1e-5,
                   maxit = 100,
                   verbose = 0,
                   BB.control = NULL,
                   ...) {
  
  method <- match.arg(method, several.ok = TRUE)
  boot.method <- match.arg(boot.method)
  Call <- match.call()
  
  # create a call to model.frame() that contains the formula (required)
  #  and any other of the relevant optional arguments
  # then evaluate it in the proper frame
  indx <- match(c("formula", "data", "weights", "subset", "na.action"),
                names(Call), nomatch=0) 
  if (indx[1] ==0) stop("A formula argument is required")
  temp <- Call[c(1,indx)]  # only keep the arguments we wanted
  temp[[1]] <- as.name('model.frame')  # change the function called
  
  temp$formula <- if(missing(data)) terms(formula)
                  else terms(formula, data=data)
  
  mf <- eval(temp, parent.frame())
  if (nrow(mf) ==0) stop("No (non-missing) observations")
  Terms <- terms(mf)
  
  Y <- model.extract(mf, "response")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
  type <- attr(Y, "type")
  if (type!='right' && type!='counting')
    stop(paste("AFT model doesn't support \"", type,
               "\" survival data", sep=''))
  weights <- model.weights(mf)
  nobs <- nrow(Y)   #remember this before any time transforms
  
  types <- c("AFT", "AFT-IV", "AFT-2SLS", "AFT-IPCW")
  funcs <- c("AFTScorePre", "AFTivScorePre", "AFT2SLSScorePre", "AFTivIPCWScorePre")
  funcs.sm <- c("AFTScoreSmoothPre", "AFTivScoreSmoothPre", "AFT2SLSScoreSmoothPre")
  
  contrast.arg <- NULL  #due to shared code with model.matrix.coxph
  
  X <- model.matrix(Terms, mf, contrasts=contrast.arg)
  #print(head(X))
  # drop the intercept after the fact, and also drop strata if necessary
  
  Xatt <- attributes(X) 
  adrop = 0
  xdrop <- Xatt$assign %in% adrop 
  #xdrop <- Xatt$assign %in% adrop  #columns to drop (always the intercept)
  X <- X[, !xdrop, drop=FALSE]
  attr(X, "assign") <- Xatt$assign[!xdrop]
  if (any(adrop>0)) attr(X, "contrasts") <- Xatt$contrasts[-adrop]
  else attr(X, "contrasts") <- Xatt$contrasts
  attr(X, "contrasts") <- Xatt$contrasts
  offset <- model.offset(mf)
  if (is.null(offset) | all(offset==0)) offset <- rep(0., nrow(mf))
  
  assign <- attrassign(X, Terms)
  contr.save <- attr(X, "contrasts")
  nvars <- ncol(X)
  if (missing(init)) init <- rep(0, nvars)
  
  surv.dat <- vector(mode = "list", length = 3)
  names(surv.dat) <- c("survival", "X", "Z")
  surv.dat[[1]] <- data.frame(Y[,1], Y[,2])

  colnames(surv.dat[[1]]) <- c("log.t", "delta")
  surv.dat[[2]] <- X
  surv.dat[[3]] <- instrument
  attr(surv.dat, "class") <- "survival.data"
  
  beta <- se.hat <- array(0, dim = c(length(type), nvars))
  fit.objects <- bootstrap.objects <- vector(length = length(type), mode = "list")
  
  quiet <- TRUE
  trace <- FALSE
  if (verbose) {
    trace <- TRUE
    if (verbose <= 1) {
      quiet <- FALSE
    }
  }
  
  if (is.null(BB.control)) {
    BB.ns <- TRUE
  } else {
    BB.ns <- FALSE
  }

  #solve each estimating equation for beta
  for (e in 1:length(method)){
    #return correct estimating equation function
    est.eqn <- match.fun(funcs[[match(method[e], types)]])
    est.eqn.sm <- match.fun(funcs.sm[[match(method[e], types)]])
    attr(est.eqn, "name") <- funcs[[match(method[e], types)]]
    attr(est.eqn.sm, "name") <- funcs.sm[[match(method[e], types)]]
    
    if (BB.ns) {
      dfsane.tol <- 0.000001
      if (method[e] == "AFT-IPCW") {dfsane.tol <- 8}
      BB.control <- list(tol = dfsane.tol, maxit = 1000, 
                         trace = trace, M = c(500))
    }
    ssf <- 1e10
    
    ct <- 0
    if (e == 1) {
      init.par <- init
    } else {
      init.par <- est$par
      init.par <- NULL
    }
    
    if (verbose) {
      cat("Fitting model", e)
    }
    
    #solve for beta using deriv-free spectral method
    est <- repFitAFT(tol = dfsane.tol, data = surv.dat, est.eqn = est.eqn, est.eqn.sm = est.eqn.sm,
                     instrument.names = NULL, confounded.x.names = confounded.x.names, 
                     fit.method = "dfsane", init.par = init.par, final.fit = smoothed,
                     method = c(2), control = BB.control, quiet = quiet)
    
    if (bootstrap) {
      # bootstrap estimate
      bootstrap.objects[[e]] <- bootstrapVar(aft.repfit, dat.sim, B = B, method = boot.method)
      se.hat[e, ] <- bootstrap.objects[[e]]$se.hat
    }
    
    beta[e, ] <- est$par
    fit.objects[[e]] <- est
    
    
  }
  list(beta = beta, se = se.hat, fit.objects = fit.objects, 
       bootstrap.objects = bootstrap.objects,
       B = B, Call = Call)
}