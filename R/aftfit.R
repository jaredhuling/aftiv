

aftfit <- function(formula, 
                   data, 
                   instrument, 
                   confounded.x.names, 
                   method              = c("AFT", "AFT-IV", "AFT-2SLS", "AFT-IPCW"),
                   smoothed            = FALSE,
                   weights             = NULL, 
                   bootstrap           = FALSE,
                   boot.method         = c("ls", "sv", "full.bootstrap"),
                   B                   = 1000L,
                   dependent.censoring = FALSE,
                   na.action, 
                   init                = NULL, 
                   return.data         = FALSE,
                   tol                 = 1e-5,
                   maxit               = 100L,
                   verbose             = 0,
                   BB.control          = NULL,
                   ...) {
  
  method      <- match.arg(method, several.ok = TRUE)
  boot.method <- match.arg(boot.method)
  Call        <- match.call()

  # create a call to model.frame() that contains the formula (required)
  #  and any other of the relevant optional arguments
  # then evaluate it in the proper frame
  indx <- match(c("formula", "data", "weights", "subset", "na.action"),
                names(Call), nomatch=0) 
  if (indx[1] == 0) stop("A formula argument is required")
  temp      <- Call[c(1,indx)]  # only keep the arguments we wanted
  temp[[1]] <- as.name('model.frame')  # change the function called
  
  temp$formula <- if(missing(data)) terms(formula)
                  else terms(formula, data=data)
  
  mf    <- eval(temp, parent.frame())
  if (nrow(mf) == 0) stop("No (non-missing) observations")
  Terms <- terms(mf)
  
  Y <- model.extract(mf, "response")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
  type <- attr(Y, "type")
  if (type!='right' && type!='counting')
    stop(paste("AFT model doesn't support \"", type,
               "\" survival data", sep=''))
  weights <- model.weights(mf)
  nobs    <- nrow(Y)   #remember this before any time transforms
  
  types    <- c("AFT", "AFT-IV", "AFT-2SLS", "AFT-IPCW")
  funcs    <- c("AFTScorePre", "AFTivScorePre", "AFT2SLSScorePre", "AFTivIPCWScorePre")
  funcs.sm <- c("AFTScoreSmoothPre", "AFTivScoreSmoothPre", "AFT2SLSScoreSmoothPre", "AFTivIPCWScoreSmooth")
  
  contrast.arg <- NULL  #due to shared code with model.matrix.coxph
  
  X <- model.matrix(Terms, mf, contrasts=contrast.arg)
  #print(head(X))
  # drop the intercept after the fact, and also drop strata if necessary
  
  Xatt  <- attributes(X) 
  adrop <- 0
  xdrop <- Xatt$assign %in% adrop 
  #xdrop <- Xatt$assign %in% adrop  #columns to drop (always the intercept)
  X <- X[, !xdrop, drop=FALSE]
  attr(X, "assign") <- Xatt$assign[!xdrop]
  if (any(adrop>0)) attr(X, "contrasts") <- Xatt$contrasts[-adrop]
  else attr(X, "contrasts") <- Xatt$contrasts
  attr(X, "contrasts") <- Xatt$contrasts
  offset <- model.offset(mf)
  if (is.null(offset) | all(offset==0)) offset <- rep(0., nrow(mf))
  
  assign     <- attrassign(X, Terms)
  contr.save <- attr(X, "contrasts")
  nvars      <- ncol(X)
  if (missing(init)) init <- rep(0, nvars)
  
  surv.dat        <- vector(mode = "list", length = 3)
  names(surv.dat) <- c("survival", "X", "Z")
  surv.dat[[1]]   <- data.frame(Y[,1], Y[,2])
  
  if (is.null(dim(instrument))) 
  {
    instrument <- matrix(instrument, ncol = 1)
  }
  
  stopifnot(nrow(X) == nrow(instrument))

  colnames(surv.dat[[1]]) <- c("log.t", "delta")
  surv.dat[[2]]           <- X
  surv.dat[[3]]           <- instrument
  attr(surv.dat, "class") <- "survival.data"
  
  beta               <- se.hat                   <- array(0, dim = c(length(method), nvars))
  rownames(beta)     <- rownames(se.hat)         <- method
  colnames(beta)     <- colnames(se.hat)         <- colnames(X)
  fit.objects        <- bootstrap.objects        <- vector(length = length(method), mode = "list")
  names(fit.objects) <- names(bootstrap.objects) <- method
  
  quiet <- TRUE
  trace <- FALSE
  if (verbose > 1) 
  {
    trace <- TRUE
    quiet <- FALSE
  }
  
  if (is.null(BB.control)) 
  {
    BB.ns <- TRUE
  } else 
  {
    BB.ns <- FALSE
  }

  # solve each estimating equation for beta
  for (e in 1:length(method))
  {
    # return correct estimating equation function
    est.eqn                  <- match.fun(funcs[[match(method[e], types)]])
    attr(est.eqn, "name")    <- funcs[[match(method[e], types)]]
    
    # and return its smooth counterpart
    est.eqn.sm               <- match.fun(funcs.sm[[match(method[e], types)]])
    attr(est.eqn.sm, "name") <- funcs.sm[[match(method[e], types)]]

    if (BB.ns) 
    {
      dfsane.tol <- 0.000001
      if (method[e] == "AFT-IPCW") 
      {
        dfsane.tol <- 8
      }
      BB.control <- list(tol   = dfsane.tol, 
                         maxit = 1000, 
                         trace = trace, 
                         M     = c(500))
    }
    ssf <- 1e10
    
    ct <- 0
    if (e == 1) 
    {
      init.par <- init
    } else 
    {
      init.par <- est$par
      #init.par <- NULL
    }
    
    if (verbose) 
    {
      cat("Fitting model", e, "\n")
    }
    
    # solve for beta using deriv-free spectral method
    est <- repFitAFT(tol                 = dfsane.tol, 
                     data                = surv.dat, 
                     est.eqn             = est.eqn, 
                     est.eqn.sm          = est.eqn.sm,
                     instrument.names    = NULL, 
                     confounded.x.names  = confounded.x.names, 
                     maxit               = maxit, 
                     dependent.censoring = dependent.censoring,
                     fit.method          = "dfsane", 
                     init.par            = init.par, 
                     final.fit           = smoothed,
                     method              = c(2), 
                     control             = BB.control, 
                     quiet               = quiet)
    
    if (bootstrap) 
    {
      if (verbose) 
      {
        cat("Bootstrapping model", e, "\n")
      }
      
      # bootstrap estimate
      bootstrap.objects[[e]] <- bootstrapVar(est, surv.dat, B = B, method = boot.method)
      se.hat[e, ]            <- bootstrap.objects[[e]]$se.hat
    }
    
    beta[e, ]        <- est$par
    fit.objects[[e]] <- est
    
    
  }
  ret <- list(beta              = beta, 
              se                = se.hat, 
              fit.objects       = fit.objects, 
              bootstrap.objects = bootstrap.objects,
              B                 = B, 
              Call              = Call)
  class(ret) <- "aftfits"
  ret
}


repFitAFT <- function(tol                 = 5, 
                      maxit               = 25, 
                      data, 
                      est.eqn             = NULL, 
                      est.eqn.sm          = NULL, 
                      instrument.names, 
                      confounded.x.names,
                      init.par            = NULL, 
                      init.method         = c("lm", "bisection"), 
                      final.fit           = TRUE,
                      dependent.censoring = FALSE,
                      fit.method          = c("dfsane", "multiStart", "nleqslv", "sane"), 
                      ...) 
{
  # repeatedly calls fitAFT until convergence. decreasing noise 
  # is added to coefficients at each call to encourage solution
  ct  <- 0
  ssf <- best.ssf <- 1e10
  while (ssf > tol & ct <= maxit) 
  {
    ct <- ct + 1
    old.ssf <- ssf
    
    #solve for beta using deriv-free spectral method
    est <- fitAFT(data                = data, 
                  est.eqn             = est.eqn, 
                  instrument.names    = instrument.names, 
                  confounded.x.names  = confounded.x.names, 
                  dependent.censoring = dependent.censoring,
                  fit.method          = fit.method, 
                  init.par            = init.par, 
                  ...)
    ssf <- est$sum.sq.fval
    sd  <- 5e-3 * min(sqrt(best.ssf), 50)
    if (ssf < best.ssf) 
    {
      best.est <- est
      best.ssf <- ssf
      init.par <- best.est$par + rnorm(length(est$par), sd = sd)
    } else 
    {
      init.par <- best.est$par + rnorm(length(est$par), sd = sd)
    }
    print (sprintf("Current ssf: %g   Best ssf: %g, sd: %g", ssf, best.ssf, sd))
  }
  if (final.fit) 
  {
    best.est <- fitAFT(data                = data,
                       est.eqn             = est.eqn.sm, 
                       instrument.names    = instrument.names, 
                       confounded.x.names  = confounded.x.names, 
                       fit.method          = "nleqslv", 
                       init.par            = init.par, 
                       global              = "dbldog", 
                       method              = "Broyden", 
                       dependent.censoring = dependent.censoring,
                       control             = list(ftol  = 1e-3, 
                                                  btol  = 1e-4, 
                                                  trace = 1))
  }
  best.est$est.eqn        <- est.eqn
  best.est$est.eqn.smooth <- est.eqn.sm
  best.est$final.fit      <- final.fit
  best.est$call           <- match.call()
  best.est
}



fitAFT <- function(data, 
                   est.eqn             = NULL, 
                   instrument.names, 
                   confounded.x.names, 
                   init.par            = NULL, 
                   init.method         = c("lm", "bisection"),
                   fit.method          = c("dfsane", "multiStart", "nleqslv", "sane"), 
                   dependent.censoring = FALSE, 
                   ...) 
{
  fit.method <- match.arg(fit.method)
  if (fit.method == "dfsane" & is.null(est.eqn)) 
  {
    stop("Please supply est.eqn if using fit.method: df.sane")
  }
  GC          <- NULL
  init.method <- match.arg(init.method)
  #instr.loc  <- match(instrument.names, colnames(data$X))
  conf.x.loc  <- match(confounded.x.names, colnames(data$X))
  ZXmat       <- data$X
  
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
  
  
  if (is.null(init.par)) 
  {
    #Xhat <- lm(data$X[,conf.x.loc] ~ data$Z)$fitted.values
    XXhat <- as.matrix(data$X)
    #XXhat[,conf.x.loc] <- Xhat
    init.par <- lm(data$survival$log.t ~ XXhat-1)$coefficients
    init.par[which(is.na(init.par))] <- 0
    names(init.par) <- colnames(data$X)
    #num.vars <- length(init.par)
    #names(init.par)[conf.x.loc] <- instrument.names
    if (init.method == "bisection") 
    {
      interval <- list()
      for (i in 1:num.vars) 
      {
        interval[[i]] <- c(init.par[i] - 1, init.par[i] + 1)
      }
      init.par <- coordinateBisection(est.eqn, 
                                      interval = interval, 
                                      num.vars = num.vars, 
                                      max.iter = 500, 
                                      bisection.tol = 0.25,
                                      survival = data$survival, 
                                      X = as.matrix(data$X), 
                                      ZXmat = as.matrix(ZXmat))
      print (init.par)
      print (sum(est.eqn(init.par, survival=data$survival, X = as.matrix(data$X), ZXmat = as.matrix(ZXmat)) ^ 2))
    }
    names(init.par) <- colnames(data$X)
  }
  
  
  
  if (fit.method == "dfsane") 
  {
    #Derivative-Free Spectral Approach for solving nonlinear systems of equations
    #from CRAN package 'BB'
    if (attr(est.eqn, "name") == "AFTivIPCWScorePre") 
    {
      #generate function G_c() for ICPW 
      if (dependent.censoring)
      {
        GC <- genKMCensoringFunc(data$survival, cox = TRUE, X = as.matrix(cbind(data$X, data$Z)))
      } else 
      {
        GC <- genKMCensoringFunc(data$survival)
      }
      df   <- BBsolve(par        = init.par, 
                      fn         = est.eqn, 
                      ...,
                      survival   = data$survival, 
                      X          = as.matrix(data$X), 
                      ZXmat      = as.matrix(ZXmat), 
                      Z          = data$Z, 
                      GC         = GC, 
                      conf.x.loc = conf.x.loc)
      
      fval <- est.eqn(beta       = df$par, 
                      survival   = data$survival, 
                      X          = as.matrix(data$X), 
                      ZXmat      = as.matrix(ZXmat), 
                      Z          = data$Z, 
                      GC         = GC,
                      conf.x.loc = conf.x.loc)
    } else 
    {
      df   <- BBsolve(par      = init.par, 
                      fn       = est.eqn, 
                      ...,
                      survival = data$survival, 
                      X        = as.matrix(data$X), 
                      ZXmat    = as.matrix(ZXmat))
      
      fval <- est.eqn(beta     = df$par, 
                      survival = data$survival, 
                      X        = as.matrix(data$X), 
                      ZXmat    = as.matrix(ZXmat))
    }
    ret <- list(par                = df$par, 
                fval               = fval, 
                iter               = df$iter,
                nobs               = nrow(data$X), 
                nvars              = length(df$par),
                GC                 = GC, 
                instrument.names   = instrument.names, 
                confounded.x.names = confounded.x.names)
  } else if (fit.method == "nleqslv" | fit.method == "sane") 
  {
    if (!is.null(est.eqn)) 
    {
      if (attr(est.eqn, "name") != "AFTScoreSmoothPre" & attr(est.eqn, "name") != "AFTivScoreSmoothPre") 
      {
        
        warning(paste(paste("Arg: est.eqn =", attr(est.eqn, "name")), 
                      "not be used. Used AFTivScoreSmoothPre instead"))
        est.eqn <- AFTivScoreSmoothPre
      }
    } else 
    {
      warning("est.eqn not supplied. Used AFTivScoreSmoothPre")
      est.eqn <- AFTivScoreSmoothPre
    }
    if (fit.method == "nleqslv") 
    {
      df  <- nleqslv(x        = init.par, 
                     fn       = est.eqn, 
                     ..., 
                     survival = data$survival, 
                     X        = as.matrix(data$X), 
                     ZXmat    = as.matrix(ZXmat), 
                     tau      = 0.01)
      
      ret <- list(par                = df$x, 
                  fval               = df$fvec, 
                  iter               = df$iter,
                  nobs               = nrow(data$X), 
                  nvars              = length(df$par),
                  GC                 = GC, 
                  instrument.names   = instrument.names, 
                  confounded.x.names = confounded.x.names)
      #print(df)
    } else if (fit.method == "sane") 
    {
      df   <- sane(par      = init.par, 
                   fn       = est.eqn, 
                   ...,
                   survival = data$survival, 
                   X        = as.matrix(data$X), 
                   ZXmat    = as.matrix(ZXmat), 
                   tau      = 0.01)
      
      fval <- est.eqn(beta     = df$par, 
                      survival = data$survival, 
                      X        = as.matrix(data$X), 
                      ZXmat    = as.matrix(ZXmat),
                      tau      = 0.01)
      
      ret  <- list(par                = df$par, 
                   fval               = fval, 
                   iter               = df$iter,
                   nobs               = nrow(data$X), 
                   nvars              = length(df$par),
                   GC                 = GC, 
                   instrument.names   = instrument.names, 
                   confounded.x.names = confounded.x.names)
    }
  } else if (fit.method == "multiStart") 
  {
    #Derivative-Free Spectral Approach for solving nonlinear systems of equations
    #from CRAN package 'BB'
    n.starts <- 10
    p        <- length(init.par)
    p0       <- matrix(rnorm(n.starts * p), n.starts, p) 
    p0       <- rbind(p0, init.par)
    if (attr(est.eqn, "name") == "AFTivIPCWScorePre") 
    {
      #generate function G_c() for ICPW 
      if (dependent.censoring)
      {
        GC <- genKMCensoringFunc(data$survival, cox = TRUE, X = as.matrix(cbind(data$X, data$Z)))
      } else 
      {
        GC <- genKMCensoringFunc(data$survival)
      }
      df   <- multiStart(par      = p0, 
                         fn       = est.eqn, 
                         ..., 
                         action   = "solve",
                         survival = data$survival, 
                         X        = as.matrix(data$X), 
                         ZXmat    = as.matrix(ZXmat), 
                         Z        = as.matrix(data$Z), 
                         GC       = GC)
      
      fval <- est.eqn(beta     = df$par, 
                      survival = data$survival, 
                      X        = as.matrix(data$X), 
                      ZXmat    = as.matrix(ZXmat), 
                      GC       = GC)
    } else 
    {
      df <- multiStart(par      = p0,  
                       fn       = est.eqn, 
                       ..., 
                       action   = "solve",
                       survival = data$survival, 
                       X        = as.matrix(data$X), 
                       ZXmat    = as.matrix(ZXmat))
      #take only the best ones and fit again
      good.init.idx <- which(df$fvalue < 500)
      p0 <- df$par[good.init.idx,] + matrix(rnorm(length(good.init.idx) * p, sd = 5e-6), length(good.init.idx), p) 
      
      df <- multiStart(par      = p0, 
                       fn       = est.eqn, 
                       ..., 
                       action   = "solve",
                       survival = data$survival, 
                       X        = as.matrix(data$X), 
                       ZXmat    = as.matrix(ZXmat))
      
      good.init.idx <- which(df$fvalue < 100)
      df$par        <- df$par[good.init.idx,]
      df$fvalue     <- df$fvalue[good.init.idx]
      colnames(p0)  <- colnames(data$X)
      fval <- est.eqn(beta     = df$par, 
                      survival = data$survival, 
                      X        = as.matrix(data$X), 
                      ZXmat    = as.matrix(ZXmat))
    }
    ret <- list(par = df$par, fval = fval, iter = df$iter, 
                nobs = nrow(data$X), nvars = length(df$par),
                GC = GC, 
                instrument.names = instrument.names, 
                confounded.x.names = confounded.x.names)
  }
  ret$sum.sq.fval <- sum(ret$fval^2)
  if (init.method == "multiStart") 
  {
    ret <- df
  }
  ret$call   <- match.call()
  class(ret) <- "aft.fit"
  ret
}







bootstrap.model <- function(model, 
                            formula, 
                            data, 
                            instrument, 
                            confounded.x.names, 
                            method      = c("AFT", "AFT-IV", "AFT-2SLS", "AFT-IPCW"),
                            smoothed    = FALSE,
                            weights     = NULL, 
                            bootstrap   = FALSE,
                            boot.method = c("ls", "sv", "full.bootstrap"),
                            B           = 100L,
                            na.action, 
                            init        = NULL, 
                            return.data = FALSE,
                            tol         = 1e-5,
                            maxit       = 100,
                            verbose     = 0,
                            BB.control  = NULL,
                            ...)
{
  
  method      <- match.arg(method, several.ok = TRUE)
  boot.method <- match.arg(boot.method)
  Call        <- match.call()
  
  # create a call to model.frame() that contains the formula (required)
  #  and any other of the relevant optional arguments
  # then evaluate it in the proper frame
  indx <- match(c("formula", "data", "weights", "subset", "na.action"),
                names(Call), nomatch = 0) 
  if (indx[1] == 0) stop("A formula argument is required")
  temp      <- Call[c(1,indx)]  # only keep the arguments we wanted
  temp[[1]] <- as.name('model.frame')  # change the function called
  
  temp$formula <- if(missing(data)) terms(formula)
  else terms(formula, data=data)
  
  mf <- eval(temp, parent.frame())
  if (nrow(mf) == 0) stop("No (non-missing) observations")
  Terms <- terms(mf)
  
  Y <- model.extract(mf, "response")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
  
  type <- attr(Y, "type")
  
  if (type!='right' && type!='counting')
  {
    stop(paste("AFT model doesn't support \"", type,
               "\" survival data", sep=''))
  }
  weights  <- model.weights(mf)
  nobs     <- nrow(Y)   #remember this before any time transforms
  
  types    <- c("AFT", "AFT-IV", "AFT-2SLS", "AFT-IPCW")
  funcs    <- c("AFTScorePre", "AFTivScorePre", "AFT2SLSScorePre", "AFTivIPCWScorePre")
  funcs.sm <- c("AFTScoreSmoothPre", "AFTivScoreSmoothPre", "AFT2SLSScoreSmoothPre", "AFTivIPCWScoreSmooth")
  
  contrast.arg <- NULL  #due to shared code with model.matrix.coxph
  
  X <- model.matrix(Terms, mf, contrasts=contrast.arg)
  #print(head(X))
  # drop the intercept after the fact, and also drop strata if necessary
  
  Xatt  <- attributes(X) 
  adrop <- 0
  xdrop <- Xatt$assign %in% adrop 
  #xdrop <- Xatt$assign %in% adrop  #columns to drop (always the intercept)
  X <- X[, !xdrop, drop=FALSE]
  attr(X, "assign") <- Xatt$assign[!xdrop]
  if (any(adrop>0)) attr(X, "contrasts") <- Xatt$contrasts[-adrop]
  else attr(X, "contrasts") <- Xatt$contrasts
  attr(X, "contrasts") <- Xatt$contrasts
  offset <- model.offset(mf)
  if (is.null(offset) | all(offset==0)) offset <- rep(0., nrow(mf))
  
  assign     <- attrassign(X, Terms)
  contr.save <- attr(X, "contrasts")
  nvars      <- ncol(X)
  if (missing(init)) init <- rep(0, nvars)
  
  surv.dat        <- vector(mode = "list", length = 3)
  names(surv.dat) <- c("survival", "X", "Z")
  surv.dat[[1]]   <- data.frame(Y[,1], Y[,2])
  
  if (is.null(dim(instrument))) 
  {
    instrument <- matrix(instrument, ncol = 1)
  }
  
  stopifnot(nrow(X) == nrow(instrument))
  
  colnames(surv.dat[[1]]) <- c("log.t", "delta")
  surv.dat[[2]]           <- X
  surv.dat[[3]]           <- instrument
  attr(surv.dat, "class") <- "survival.data"
  
  beta               <- se.hat                   <- array(0, dim = c(length(method), nvars))
  rownames(beta)     <- rownames(se.hat)         <- method
  colnames(beta)     <- colnames(se.hat)         <- colnames(X)
  fit.objects        <- bootstrap.objects        <- vector(length = length(method), mode = "list")
  names(fit.objects) <- names(bootstrap.objects) <- method
  
  quiet <- TRUE
  trace <- FALSE
  
  if (verbose > 1) 
  {
    trace <- TRUE
    quiet <- FALSE
  }
  
  if (is.null(BB.control)) 
  {
    BB.ns <- TRUE
  } else 
  {
    BB.ns <- FALSE
  }
  
  #solve each estimating equation for beta
  for (e in 1:length(method))
  {
    # return correct estimating equation function
    est.eqn                  <- match.fun(funcs[[match(method[e], types)]])
    attr(est.eqn, "name")    <- funcs[[match(method[e], types)]]
    
    # and return the smooth counterpart
    est.eqn.sm               <- match.fun(funcs.sm[[match(method[e], types)]])
    attr(est.eqn.sm, "name") <- funcs.sm[[match(method[e], types)]]
    
    
    
    ct <- 0
    
    if (bootstrap) 
    {
      if (verbose) 
      {
        cat("Bootstrapping model", e, "\n")
      }
      # bootstrap estimate
      bootstrap.objects[[e]] <- bootstrapVar(model$fit.objects[[e]], surv.dat, B = B, method = boot.method)
      se.hat[e, ]            <- bootstrap.objects[[e]]$se.hat
    }
    
    #beta[e, ] <- est$par
    #fit.objects[[e]] <- est
    
    
  }
  ret <- model
  ret$bootstrap.objects <- bootstrap.objects
  ret$se <- se.hat
  class(ret) <- "aftiv.boot.fit"
  ret
}
