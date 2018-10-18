
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
true.beta <- c(1,-0.5,-0.5,0.25,0,0,-0.75,0.75)
dat <- simIVMultivarSurvivalData(500,1,1,-1,1,true.beta,num.confounded = 1,
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
#> [1] "Current ssf: 1.4995e-06   Best ssf: 1.4995e-06, sd: 0.1"
#> [1] "Current ssf: 6.9141e-07   Best ssf: 6.9141e-07, sd: 2.44908e-06"
#> [1] "Current ssf: 2.56856e-07   Best ssf: 2.56856e-07, sd: 0.1"
#> [1] "Current ssf: 9.77557e-07   Best ssf: 9.77557e-07, sd: 0.1"
#> [1] "Current ssf: 0.108804   Best ssf: 0.108804, sd: 0.1"
#>    user  system elapsed 
#>  20.323   1.679  22.075
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
#>        coef se(coef) lower .95 upper .95      z Pr(>|z|)    
#> X1  1.09747  0.04677   1.00580   1.18914 23.466  < 2e-16 ***
#> X2 -0.46356  0.09811  -0.65584  -0.27127 -4.725 2.30e-06 ***
#> X3 -0.35479  0.10184  -0.55440  -0.15518 -3.484 0.000495 ***
#> X4  0.05797  0.12096  -0.17911   0.29504  0.479 0.631779    
#> X5 -0.02851  0.08339  -0.19196   0.13494 -0.342 0.732440    
#> X6  0.18180  0.11646  -0.04645   0.41005  1.561 0.118506    
#> X7 -0.72951  0.10480  -0.93492  -0.52410 -6.961 3.38e-12 ***
#> X8  0.69238  0.14131   0.41542   0.96933  4.900 9.59e-07 ***
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
#>        coef se(coef) lower .95 upper .95      z Pr(>|z|)    
#> X1  1.13557  0.16621   0.80981   1.46134  6.832 8.36e-12 ***
#> X2 -0.61031  0.17777  -0.95874  -0.26188 -3.433 0.000597 ***
#> X3 -0.32688  0.20328  -0.72531   0.07154 -1.608 0.107832    
#> X4  0.13173  0.18370  -0.22832   0.49177  0.717 0.473334    
#> X5 -0.03477  0.17683  -0.38135   0.31181 -0.197 0.844106    
#> X6  0.28027  0.24526  -0.20043   0.76097  1.143 0.253142    
#> X7 -0.65419  0.19761  -1.04150  -0.26688 -3.310 0.000931 ***
#> X8  0.81627  0.24276   0.34047   1.29206  3.362 0.000772 ***
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
#> X1  1.01349  0.05289   0.90984   1.11715 19.164  < 2e-16 ***
#> X2 -0.46464  0.09680  -0.65437  -0.27492 -4.800 1.59e-06 ***
#> X3 -0.33461  0.12385  -0.57735  -0.09187 -2.702   0.0069 ** 
#> X4  0.05092  0.11734  -0.17906   0.28090  0.434   0.6643    
#> X5 -0.03376  0.08862  -0.20745   0.13993 -0.381   0.7032    
#> X6  0.18014  0.12616  -0.06713   0.42740  1.428   0.1533    
#> X7 -0.68517  0.11376  -0.90813  -0.46221 -6.023 1.71e-09 ***
#> X8  0.69186  0.13769   0.42199   0.96172  5.025 5.04e-07 ***
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
#> X1  0.95980  0.11552   0.73338   1.18622  8.308  < 2e-16 ***
#> X2 -0.47636  0.06357  -0.60095  -0.35177 -7.494 6.69e-14 ***
#> X3 -0.30442  0.09209  -0.48490  -0.12393 -3.306 0.000947 ***
#> X4  0.14009  0.11519  -0.08568   0.36586  1.216 0.223927    
#> X5  0.08107  0.07087  -0.05783   0.21998  1.144 0.252633    
#> X6  0.04462  0.08193  -0.11597   0.20520  0.545 0.586051    
#> X7 -0.55813  0.10900  -0.77177  -0.34449 -5.120 3.05e-07 ***
#> X8  0.64701  0.11152   0.42843   0.86559  5.802 6.57e-09 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
