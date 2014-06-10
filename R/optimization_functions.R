library(compiler)
enableJIT(3)


stableNewtonRaphson <- function(f.sm, f.ind = NULL, grad.func = NULL, num.vars, 
                                init.method = c("lm", "bisection"),
                                interval = NULL, bisection.tol = 0.25, bisection.steps = 500, 
                                stop.tol = 1e-05, eps = 1e-06, ...) {
  init.method <- match.arg(init.method)
  if (init.method == "bisection"){
    if (is.null(interval)){
      if (is.null(data.simu$survival$log.t)) {
        data.simu$survival$log.t <- log(data.simu$survival$t)
      }
      Xhat <- lm(data.simu$X ~ data.simu$Z)$fitted.values
      est <- lm(data.simu$survival$log.t ~ Xhat-1)$coefficients
      est <- lm(data.simu$survival$log.t ~ data.simu$X-1)$coefficients
      interval <- list()
      for (i in 1:length(est)) {
        interval[[i]] <- c(est[i] - 3, est[i] + 3)
      }
    }
    if (!is.null(f.ind)) {
      #get stable initial estimate for NR with modified coordinate bisection method
      est <- coordinateBisection(f.ind, interval = interval, num.vars = num.vars, bisection.tol = bisection.tol,
                                 max.iter = bisection.steps, data.simu = dat)

    } else {
      #get stable initial estimate for NR with modified coordinate bisection method
      est <- coordinateBisection(f.sm, interval = interval, num.vars = num.vars, bisection.tol = bisection.tol,
                                 max.iter = bisection.steps, data.simu = dat, stdev = stdev)
    }
  } else if(init.method == "lm"){
    if (is.null(data.simu$survival$log.t)) {
      data.simu$survival$log.t <- log(data.simu$survival$t)
    }
    Xhat <- lm(data.simu$X ~ data.simu$Z)$fitted.values
    est <- lm(data.simu$survival$log.t ~ Xhat-1)$coefficients
  }
  
  grad.est <- gradientDescent(est, f = f.sm, grad.func  = grad.func, 
                              data.simu = dat, stdev = stdev, nmax=2)
  #grad.est$final
  #grad.est <- optim(est, fn=absAFTivScore, method = "Nelder-Mead", data.simu=dat, control=list(maxit=500))
  #initialize NR with coordinate bisection estimates
  nr.est <- dampedNewtonRaphsonRoot(grad.est$final, f = f.sm, grad.func = grad.func, stop.tol = 1e-05, eps = 1e-06, ...)
  #nr.est
  list(bisection.est = est, nstep = nr.est$nstep, final.est = nr.est$final, func.val = nr.est$funcval)
}

stableBroyden <- function(f.sm, f.ind = NULL, grad.func = NULL, num.vars, 
                          init.method = c("lm", "bisection", "par"), init.par = NULL,
                          interval = NULL, bisection.tol = 0.25, bisection.steps = 500, 
                          stop.tol = 1e-05, eps = 1e-06, ...) {
  init.method <- match.arg(init.method)
  if (init.method == "bisection"){
    if (is.null(interval)){
      if (is.null(data.simu$survival$log.t)) {
        data.simu$survival$log.t <- log(data.simu$survival$t)
      }
      Xhat <- lm(data.simu$X ~ data.simu$Z)$fitted.values
      est <- lm(data.simu$survival$log.t ~ Xhat-1)$coefficients
      est <- lm(data.simu$survival$log.t ~ data.simu$X-1)$coefficients
      interval <- list()
      for (i in 1:length(est)) {
        interval[[i]] <- c(est[i] - 3, est[i] + 3)
      }
    }
    if (!is.null(f.ind)) {
      #get stable initial estimate for NR with modified coordinate bisection method
      est <- coordinateBisection(f.ind, interval = interval, num.vars = num.vars, bisection.tol = bisection.tol,
                                 max.iter = bisection.steps, data.simu = dat)
      
    } else {
      #get stable initial estimate for NR with modified coordinate bisection method
      est <- coordinateBisection(f.sm, interval = interval, num.vars = num.vars, bisection.tol = bisection.tol,
                                 max.iter = bisection.steps, data.simu = dat, stdev = stdev)
    }
  } else if(init.method == "lm") {
    if (is.null(data.simu$survival$log.t)) {
      data.simu$survival$log.t <- log(data.simu$survival$t)
    }
    Xhat <- lm(data.simu$X ~ data.simu$Z)$fitted.values
    est <- lm(data.simu$survival$log.t ~ Xhat-1)$coefficients
  } else if(init.method == "par") {
    stopifnot(!is.null(init.par))
    est <- init.par
  }
  

  #initialize NR with lm or bisection estimates
  broyden.est <- nleqslv(est, fn = f.sm, method="Broyden", global="gline", ...)
  #nr.est
  list(init.est = est, nstep = broyden.est$iter, final.est = broyden.est$x, func.val = broyden.est$fvec)
}

newtonRaphsonRoot <- function(initial.par, f, nmax = 25, stop.tol = 1e-05, eps = 1e-06, grad.func = NULL, ...) {
  
  if(is.null(grad.func)) {
    grad.func = function(x) gradient(f, x=true.beta, ...)
  }
  ctr = 0
  new.par = initial.par
  old.par = initial.par - 1
  while(ctr < nmax & sqrt(sum((new.par - old.par)^2)) > stop.tol) {
    old.par = new.par
    new.par = old.par - solve(grad.func(old.par, ...)) %*% f(old.par, ...)
    ctr = ctr + 1
  }
  list(nstep = ctr, initial = initial.par, final = new.par,
       funcval = f(new.par, ...))
}

dampedNewtonRaphsonRoot <- function(initial.par, f, nmax = 25, stop.tol = 1e-05, eps = 1e-06, grad.func = NULL, ...) {
  
  if(is.null(grad.func)) {
    grad.func = function(x) gradient(f, x=true.beta, ...)
  }
  ctr = 0
  new.par = initial.par
  old.par = initial.par - 1
  while(ctr < nmax & sqrt(sum((new.par - old.par)^2)) > stop.tol) {
    old.par = new.par
    delta = solve(grad.func(old.par, ...)) %*% f(old.par, ...)
    new.par = old.par - delta
    k <- 0
    while(sum(f(new.par, ...)^2) > sum(f(old.par, ...)^2)) {
      k <- k + 1
      old.par = new.par
      delta = delta / 2
      new.par = old.par - delta
      print (new.par)
      if (k > 15) {
        new.par <- new.par + rnorm(length(new.par), sd = 0.001)
      }
    }
    ctr = ctr + 1
  }
  list(nstep = ctr, initial = initial.par, final = new.par,
       funcval = f(new.par, ...))
}



gradientDescent <- function(initial.par, f, nmax = 25, stop.tol = 1e-05,
                            eps = 1e-06, grad.func = NULL, ...) {
  # function to implement gradient-descent. A line search 
  # is performed to select step size to prevent over-shooting
  
  #if(is.null(grad.func))
  #  grad.func = function(x) Gradmat(x, infcn, eps)
  steps = NULL
  possible.steps <- c(1e-7,1e-6,1e-5)
  last.step <- length(possible.steps)
  new.par = initial.par
  old.par = new.par - 1
  ctr = 0
  while(ctr < nmax & sqrt(sum((new.par - old.par)^2)) > stop.tol) {
    ctr = ctr + 1
    old.par = new.par
    J <- grad.func(old.par, ...)
    f.k <- f(old.par, ...)
    h.k <- -t(J) %*% f.k
    f.x.h <- array(0, dim = length(possible.steps))
    #find step size that minimizes f.x.h
    for (i in 1:length(possible.steps)) {
      eval <- old.par + possible.steps[i] * h.k
      f.x.h.tmp <- f(eval, ...)
      f.x.h[i] <- t(f.x.h.tmp) %*% f.x.h.tmp
    }
    last.step <- which.min(f.x.h)
    cur.step <- possible.steps[last.step]
    steps[ctr] <- cur.step
    ### update estimate
    new.par = old.par + cur.step * h.k
  }
  list(nstep = ctr, initial = initial.par, final = new.par,
       funcval = f(new.par, ...), steps = steps)
}

coordinateBisection <- function(f, interval, num.vars, max.iter=25, bisection.tol = 0.25, ...) {
  if (length(interval) != num.vars) {
    intervals <- rep(list(interval), num.vars)
  } else {intervals <- interval}
  for (i in 1:max.iter) {
    if (i == 1) {
      old.est <- unlist(lapply(intervals, function(x) sum(x) / 2))
    } else {old.est <- cur.est}
    for (j in 1:num.vars) {
      mid.pt <- sum(intervals[[j]])/2
      func.lower <- func.mid.pt <- unlist(lapply(intervals, function(x) x[1]))
      func.higher <- unlist(lapply(intervals, function(x) x[2]))
      func.mid.pt[j] <- mid.pt
      f.vals <- array(0, dim = 3)
      f.vals[1] <- f(func.lower, ...)[j]
      f.vals[2] <- f(func.mid.pt, ...)[j]
      f.vals[3] <- f(func.higher, ...)[j]
      pts <- c(intervals[[j]][1], func.mid.pt[j], intervals[[j]][2])
      if (length(unique(sign(f.vals))) == 1) {
        f.vals.sm <- f.vals[-which.min(abs(f.vals))]
        pts.sm <- pts[-which.min(abs(f.vals))]
        val.1 <- pts[which.min(abs(f.vals))]
        val.2 <- pts.sm[which.min(abs(f.vals.sm))]
        closest.func.val <- f.vals[which.min(abs(f.vals))]
        
        mid.func.val <- f.vals.sm[which.min(abs(f.vals.sm))]
        new.val <- 0.85 * (val.1 - closest.func.val * ((val.1 - val.2) / (closest.func.val - mid.func.val)))
        intervals[[j]][1] <- new.val
        intervals[[j]][2] <- val.2
      } else {
        if (sign(f.vals[2]) == sign(f.vals[1])) {
          intervals[[j]][1] <- mid.pt
        } else {
          intervals[[j]][2] <- mid.pt
        }
      }
    }
    cur.est <- unlist(lapply(intervals, function(x) sum(x) / 2))
    print (intervals)
    if (all(abs(cur.est - old.est) < bisection.tol)) {break}
  }
  estimates <- unlist(lapply(intervals, function(x) sum(x) / 2))
  estimates
}


coordinateSecant <- function(f, interval, num.vars, max.iter=25, tol = 0.25, ...) {
  beta.i <- list()
  beta.i[[1]] <- rep(interval[2], num.vars)
  beta.i[[2]] <- rep(interval[1], num.vars)
  for (i in 1:max.iter) {
    vals.1 <- beta.i[[i+1]]
    vals.2 <- beta.i[[i]]
    
    new.val <- numeric(num.vars)
    for (j in 1:num.vars) {

      val.1 <- vals.1[j]
      val.2 <- vals.2[j]
      f.vals.1 <- f(vals.1, ...) #n - 1
      f.vals.2 <- f(vals.2, ...) #n - 2
      closer.func.val <- f.vals.1[j]
      farther.func.val <- f.vals.2[j]
      
      new.val[j] <- val.1 - closer.func.val * ((val.1 - val.2) / (closer.func.val - farther.func.val))
      new.val[j] <- 0.75 * new.val[j]
      vals.2[j] <- vals.1[j] 
      vals.1[j] <- new.val[j]
    }
    beta.i[[i+2]] <- new.val
    print (new.val)
    if (all(abs(new.val - beta.i[[i+1]]) < tol)) {break}
  }
  estimates <- beta.i[[i+2]]
  estimates
}
