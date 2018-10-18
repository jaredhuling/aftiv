

#this function is used for fast at risk set calculation
cumsumRev <- function(a) rev(cumsum(rev(a)))

sigmoid <- function(x, tau) 
{
  1 / (1 + exp(-x / tau))
}

#this function checks if the working directory is ivsurv and changes it to simulations
setwd2sim <- function() 
{
  wd <- getwd();if (substr(wd, nchar(wd) - 5,nchar(wd)) == "ivsurv"){setwd(paste(getwd(), "/simulations", sep = ""))};getwd()
}

loadPackages <- function(load = TRUE) 
{
  packages <- c("survival", "ggplot2", "reshape2", "foreach", 
                "Matrix", "SimCorMultRes", "grid", "BB")
  if (load)
  {
    print("Loading...")
    for (i in packages) {library(i, character.only = TRUE); print(i)}
  }
  packages
}

#' @export
runParallel <- function(ncores = NULL)
{
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
