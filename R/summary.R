summary.aftfits <- function(object,  conf.int = 0.95, 
                            models = 1:length(object$fit.objects), scale = 1, 
                            return.summary = FALSE, ...) {
  aft <- object
  beta <- aft$beta
  if (is.null(aft$beta)) {   # Null model
    return(object)  #The summary method is the same as print in this case
  }
  nabeta <- !(is.na(beta))          #non-missing coefs
  beta2 <- beta[nabeta]
  if(is.null(beta) | is.null(aft$se))
    stop("Input is not valid")
  se <- aft$se
  n.fits <- length(models)
  
  rval.list <- vector(mode = "list", length = n.fits)
  names(rval.list) <- 1:n.fits
  
  k <- 0
  for (i in models) {
    k <- k + 1
    rval <- list(call=aft$Call, #fail=aft$fail, na.action=aft$na.action,
                 n=aft$fit.objects[[i]]$nobs)
    #if (!is.null(aft$nevent)) rval$nevent <- aft$nevent
    
    
    tmp <- cbind(beta[i,], se[i,], beta[i,]/se[i,],
                 1 - pchisq((beta[i,] / se[i,])^2, 1))
    dimnames(tmp) <- list(names(beta[i,]), c("coef", "se(coef)", "z", "Pr(>|z|)"))
    
  
    rval$beta <- tmp
    
    if (conf.int) {
      z <- qnorm((1 + conf.int)/2, 0, 1)
      beta[i,] <- beta[i,] * scale
      se[i,] <- se[i,] * scale
      tmp <- cbind(beta[i,], 
                   se[i,], 
                   beta[i,] - z * se[i,],
                   beta[i,] + z * se[i,])
      dimnames(tmp) <- list(names(beta[i,]), c("coef",
                                               "se(coef)",
                                               paste("lower .", round(100 * conf.int, 2), sep = ""),
                                               paste("upper .", round(100 * conf.int, 2), sep = "")))
      rval$conf.int <- tmp
    }
    if (is.R()) class(rval)    <-"summary.aftfit"
    else        oldClass(rval) <- "summary.aftfit"
    cat("*********************", "\n\n")
    cat("Method:", names(aft$fit.objects)[i], "\n\n")
    print(rval)
    names(rval.list)[k] <- names(aft$fit.objects)[i]
    rval.list[[k]] <- rval
  }
  
    df <- length(beta2)
  #logtest <- -2 * (aft$loglik[1] - aft$loglik[2])
  #rval$logtest <- c(test=logtest,
  #                  df=df,
  #                  pvalue=1 - pchisq(logtest, df))
  #rval$sctest <- c(test=aft$score,
  #                 df=df,
  #                 pvalue=1 - pchisq(aft$score, df))
  #rval$rsq<-c(rsq=1-exp(-logtest/aft$n),
  #            maxrsq=1-exp(2*aft$loglik[1]/aft$n))
  #rval$waldtest<-c(test=as.vector(round(aft$wald.test, 2)),
  #                 df=df,
  #                 pvalue=1 - pchisq(as.vector(aft$wald.test), df))
  #if (!is.null(aft$rscore))
  #  rval$robscore<-c(test=aft$rscore,
  #                   df=df,
  #                   pvalue=1 - pchisq(aft$rscore, df))
  #rval$used.robust<-!is.null(aft$naive.var)
  
  #if (!is.null(aft$concordance)) {
  #  if (is.matrix(aft$concordance)) temp <- colSums(aft$concordance)
  #  else temp <- aft$concordance
  #  rval$concordance <- c("concordance"= (temp[1] + temp[3]/2)/
  #                          sum(temp[1:3]), "se"= temp[5]/(2*sum(temp[1:3])))
  #}
  
  
  if (return.summary) 
  {
    return(rval.list)
  }
}


print.summary.aftfit <-
  function(x, digits = max(getOption('digits')-3, 3),  
           signif.stars = getOption("show.signif.stars"), ...) {
    if (!is.null(x$call)) {
      cat("Call:\n")
      dput(x$call)
      cat("\n")
    }
    #if (!is.null(x$fail)) {
    #  cat(" Coxreg failed.", x$fail, "\n")
    #  return()
    #}
    #savedig <- options(digits = digits)
    #on.exit(options(savedig))
    
    #omit <- x$na.action
    cat("  n=", x$n)
    #if (!is.null(x$nevent)) cat(", number of events=", x$nevent, "\n")
    #else cat("\n")
    #if (length(omit))
    #  cat("   (", naprint(omit), ")\n", sep="")
    
    #if (nrow(x$coef)==0) {   # Null model
    #  cat ("   Null model\n")
    #  return()
    #}
    
    
    if(!is.null(x$coefficients)) {
      cat("\n")
      if (is.R()) printCoefmat(x$coefficients, digits=digits,
                               signif.stars=signif.stars, ...)
      else prmatrix(x$coefficients)
    }
    if(!is.null(x$conf.int)) {
      cat("\n")
      #prmatrix(cbind(x$conf.int, x$beta[,3:4]))
      printCoefmat(cbind(x$conf.int, x$beta[,3:4]), digits = digits, signif.stars = signif.stars,
                   na.print = "NA", ...)
    } else {
      cat("\n")
      #prmatrix(x$beta)
      printCoefmat(x$beta, digits = digits, signif.stars = signif.stars,
                   na.print = "NA", ...)
    }
    cat("\n")
    
    if (!is.null(x$concordance)) {
      cat("Concordance=", format(round(x$concordance[1],3)),
          " (se =", format(round(x$concordance[2], 3)),")\n")
    }
    #cat("Rsquare=", format(round(x$rsq["rsq"],3)),
    #    "  (max possible=", format(round(x$rsq["maxrsq"],3)),
    #    ")\n" )
    
    #cat("Likelihood ratio test= ", format(round(x$logtest["test"], 2)), "  on ",
    #    x$logtest["df"], " df,", "   p=", format(x$logtest["pvalue"]),
    #    "\n", sep = "")
    #cat("Wald test            = ", format(round(x$waldtest["test"], 2)), "  on ",
    #    x$waldtest["df"], " df,", "   p=", format(x$waldtest["pvalue"]),
    #    "\n", sep = "")
    #cat("Score (logrank) test = ", format(round(x$sctest["test"], 2)), "  on ",
    #    x$sctest["df"]," df,", "   p=", format(x$sctest["pvalue"]), sep ="")
    #if (is.null(x$robscore))
    #  cat("\n\n")
    #else cat(",   Robust = ", format(round(x$robscore["test"], 2)), 
    #         "  p=", format(x$robscore["pvalue"]), "\n\n", sep="")   
    
    #if (x$used.robust)
    #  cat("  (Note: the likelihood ratio and score tests",
    #      "assume independence of\n     observations within a cluster,",
    #      "the Wald and robust score tests do not).\n")
    invisible()
  }

