
<!-- README.md is generated from README.Rmd. Please edit that file -->

# aftiv

## Installation

You can install the development version of aftiv from
[github](https://github.com/jaredhuling/aftiv) with:

``` r
devtools::install_github("jaredhuling/aftiv")
```

## Example

This is a basic example which shows you how to fit a semiparametric AFT
model with instrumental variable estimation:

``` r
library(aftiv)
```

``` r

## simulate data
set.seed(1)
true.beta <- c(1,-0.5,-0.5,0)
dat <- simIVMultivarSurvivalData(500,1,1,-1,1,true.beta,num.confounded = 1,
                                 cens.distribution = "lognormal",
                                 confounding.function = "exp")

## delta is event indicator, log.t is log of the observed time
## X are the covariates, the first of which is the exposure of interest, the 
## rest are covariates to adjust for
df <- data.frame(dat$survival[c("delta", "log.t")], dat$X)

## Z is the instrument, related to the first variable in X
Z <- dat$Z

system.time(aftf <- aftfit(Surv(log.t, delta) ~ ., data = df, 
                           instrument = Z, 
                           confounded.x.names = "X1", # name of the exposure of interest
                           method = c("AFT",       # naive, unadjusted (biased) estimator
                                      "AFT-2SLS",  # 2-stage approach that relies on IV model
                                      "AFT-IV",    # incorrect approach
                                      "AFT-IPCW"), # proposed approach of Huling, et al
                           boot.method = "ls", 
                           B = 200L, ## number of bootstrap iterations
                           bootstrap = TRUE)) ## use bootstrap for Conf Intervals
#> [1] "Current ssf: 3.4443e-08   Best ssf: 3.4443e-08, sd: 0.1"
#> [1] "Current ssf: 6.62061e-08   Best ssf: 6.62061e-08, sd: 0.1"
#> [1] "Current ssf: 3.04507e-08   Best ssf: 3.04507e-08, sd: 0.1"
#> [1] "Current ssf: 1.93511   Best ssf: 1.93511, sd: 0.1"
#>    user  system elapsed 
#>  18.587   1.265  20.198
```

Investigate results:

``` r
summary(aftf)
#> ********************* 
#> 
#> Method: AFT 
#> 
#> Call:
#> aftfit(formula = Surv(log.t, delta) ~ ., data = df, instrument = Z, 
#>     confounded.x.names = "X1", method = c("AFT", "AFT-2SLS", 
#>         "AFT-IV", "AFT-IPCW"), bootstrap = TRUE, boot.method = "ls", 
#>     B = 200L)
#> 
#>   n= 500
#>         coef  se(coef) lower .95 upper .95      z Pr(>|z|)    
#> X1  1.119493  0.053212  1.015200  1.223787 21.038  < 2e-16 ***
#> X2 -0.374554  0.098035 -0.566698 -0.182409 -3.821 0.000133 ***
#> X3 -0.336222  0.104813 -0.541652 -0.130792 -3.208 0.001337 ** 
#> X4 -0.001018  0.106226 -0.209217  0.207181 -0.010 0.992355    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> ********************* 
#> 
#> Method: AFT-2SLS 
#> 
#> Call:
#> aftfit(formula = Surv(log.t, delta) ~ ., data = df, instrument = Z, 
#>     confounded.x.names = "X1", method = c("AFT", "AFT-2SLS", 
#>         "AFT-IV", "AFT-IPCW"), bootstrap = TRUE, boot.method = "ls", 
#>     B = 200L)
#> 
#>   n= 500
#>         coef  se(coef) lower .95 upper .95      z Pr(>|z|)    
#> X1  1.037481  0.161432  0.721080  1.353882  6.427  1.3e-10 ***
#> X2 -0.501683  0.195723 -0.885292 -0.118074 -2.563   0.0104 *  
#> X3 -0.351659  0.182110 -0.708589  0.005271 -1.931   0.0535 .  
#> X4 -0.003302  0.186932 -0.369681  0.363078 -0.018   0.9859    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> ********************* 
#> 
#> Method: AFT-IV 
#> 
#> Call:
#> aftfit(formula = Surv(log.t, delta) ~ ., data = df, instrument = Z, 
#>     confounded.x.names = "X1", method = c("AFT", "AFT-2SLS", 
#>         "AFT-IV", "AFT-IPCW"), bootstrap = TRUE, boot.method = "ls", 
#>     B = 200L)
#> 
#>   n= 500
#>        coef se(coef) lower .95 upper .95      z Pr(>|z|)    
#> X1  0.91527  0.04839   0.82042   1.01012 18.913  < 2e-16 ***
#> X2 -0.35576  0.09620  -0.54430  -0.16722 -3.698 0.000217 ***
#> X3 -0.30964  0.10819  -0.52168  -0.09760 -2.862 0.004208 ** 
#> X4 -0.01589  0.10406  -0.21985   0.18806 -0.153 0.878597    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> ********************* 
#> 
#> Method: AFT-IPCW 
#> 
#> Call:
#> aftfit(formula = Surv(log.t, delta) ~ ., data = df, instrument = Z, 
#>     confounded.x.names = "X1", method = c("AFT", "AFT-2SLS", 
#>         "AFT-IV", "AFT-IPCW"), bootstrap = TRUE, boot.method = "ls", 
#>     B = 200L)
#> 
#>   n= 500
#>        coef se(coef) lower .95 upper .95      z Pr(>|z|)    
#> X1  1.02465  0.07786   0.87204   1.17726 13.160   <2e-16 ***
#> X2 -0.10566  0.12148  -0.34376   0.13243 -0.870   0.3844    
#> X3 -0.01731  0.13382  -0.27959   0.24497 -0.129   0.8971    
#> X4  0.28338  0.12154   0.04516   0.52159  2.331   0.0197 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
