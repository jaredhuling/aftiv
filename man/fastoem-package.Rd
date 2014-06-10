\name{fastoem-package}
\alias{fastoem-package}
\alias{fastoem}
\docType{package}
\title{
Orthogonalizing EM algorithm for penalized least squares problems.
}
\description{
oemfit fits least squares problems using Orthogonalizing EM algorithm. oemfit supports penalized least squares problems, including LASSO, ELASTIC-NET, etc.
}
\details{
\tabular{ll}{
Package: \tab fastoem\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-12-26\cr
License: \tab  GPL (>= 2)\cr
Depends: \tab Rcpp, RcppEigen, Matrix, svd, igraph\cr
}
oemfitFast is the main OEM function. It takes x and y as arguments for data and produces the regularization path over a grid of values for lambda. It accepts x as either a 'matrix' or 'sparseMatrix'
}
\author{
Jared Huling

Maintainer: Jared Huling <jaredhuling@gmail.com>

}
\references{
Xiong, S., Dai, B. and Qian, P. Z. G. (2011), \emph{Orthogonalizing Penalized Regression}, http://www.
stat.wisc.edu/~zhiguang/OEM.pdf

}

\keyword{ package }

\examples{
n.obs <- 1e5
n.vars <- 150

true.beta <- rnorm(n.vars)

x <- matrix(rnorm(n.obs * n.vars), n.obs, n.vars)
y <- rnorm(n.obs, sd = 3) + x \%*\% true.beta

fit <- oemfitFast(x = x, y = y, penalty = "lasso")
}
