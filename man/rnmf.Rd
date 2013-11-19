\docType{package}
\name{rnmf}
\alias{rnmf}
\alias{rnmf-package}
\title{rnmf}
\usage{
  rnmf(A, k = 5, alpha = 0, beta = 0, maxit = 50,
    tol = 0.01, trim = FALSE, ini.H = NULL, ini.W = NULL,
    ini.zeta = NULL, my.seed = NULL, my.mode = "cell",
    quiet = FALSE, nreg = 1, p1 = NA, n1 = NA)
}
\arguments{
  \item{A}{The input non-negative matrix to be factorized.
  \code{A}}

  \item{k}{The number of columns of W \code{k}}

  \item{alpha}{Non-negative number. Controls ||W||_2^2
  \code{alpha}}

  \item{beta}{Non-negative number. Controls sparsity of H
  \code{beta}}

  \item{maxit}{Positive integer. Maximum number of
  iterations.\code{maxit}}

  \item{tol}{Tolerance for convergence \code{tol}}

  \item{trim}{Logical. Trim or not?\code{trim}}

  \item{ini.H}{Initial H \code{ini.H}}

  \item{ini.W}{Initial W \code{ini.W}}

  \item{ini.zeta}{Initial zeta \code{ini.zeta}}

  \item{my.seed}{If missing initial W, use this seed to
  generate it \code{my.seed}}

  \item{my.mode}{Which variation to use? one of the
  following strings 'cell','col', 'smooth', 'row',
  'rowsmooth' \code{my.mode}}

  \item{quiet}{Display a short report on screen after the
  function is run? \code{quiet}}

  \item{nreg}{number of runs in "cell" variation.
  \code{nreg}}

  \item{p1}{Don't ask questions. \code{p1}}

  \item{n1}{Don't ask questions. \code{n1}}
}
\value{
  output W and H.
}
\description{
  rnmf

  Performs robust non-negative matrix factorizaion on a
  non-negative matrix A. A~=WH, where A is a p by n matrix,
  W is a p by k matrix and H is a k by n matrix. (This
  parts to be finished)
}
\examples{
# R code here showing how your function works (To be finished)
}
\keyword{keywords}

