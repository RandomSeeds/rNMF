\docType{package}
\name{rnmf}
\alias{rnmf}
\alias{rnmf-package}
\title{rnmf}
\usage{
  rnmf(A, k = 5, alpha = 0, beta = 0, maxit = 50,
    tol = 0.01, trim = FALSE, ini.H = NULL, ini.W = NULL,
    ini.zeta = NULL, my.seed = NULL, variation = "cell",
    quiet = FALSE, nreg = 1, p1 = NA, n1 = NA)
}
\arguments{
  \item{A}{A non-negative numerical matrix. The input
  non-negative matrix to be decomposed into W \%*\% H.}

  \item{k}{A non-negative integer. The number of columns of
  W.}

  \item{alpha}{Non-negative number. Controls ||W||_2^2.}

  \item{beta}{Non-negative number. Controls sparsity of H.}

  \item{maxit}{Positive integer. Maximum number of
  iterations.}

  \item{tol}{A number between 0 and 1. Tolerance for
  convergence. We suggest a number smaller than 0.01.}

  \item{trim}{A number between 0 and 1. Percentage of
  information to be trimmed. If trim = 0, the procedure
  will be .}

  \item{ini.H}{A numeric matrix. Initial right matrix H in
  A ~= W \%*\% H. Must have the same number of columns as
  A, and k rows.}

  \item{ini.W}{A numeric matrix. Initial left matrix W in A
  ~= W \%*\% H. Must have the same number of rows as A, and
  k columns.}

  \item{ini.zeta}{A logical matrix. Must have the same
  dimension as A. A "FALSE" in ini.zeta corresponds to an
  entry marked as to be trimmed. The number of "FALSE" in
  ini.zeta must be <= round(trim * p * n).}

  \item{my.seed}{Numeric. The seed for initialization of W
  or H. Ignored if ini.H or ini.W are not NULL.}

  \item{variation}{Character string. Which variation to
  use? The options are: 'cell','col', 'smooth', 'row',
  'rowsmooth'}

  \item{quiet}{Logical. Display a short report on screen
  after the function is run?}

  \item{nreg}{number of runs in "cell" variation. Not
  tested. DON'T change its default value. (and don't ask
  questions).}

  \item{p1}{Don't ask questions.}

  \item{n1}{Don't ask questions.}
}
\value{
  A list of the following items: "W": Left matrix of the
  decomposition. "H": Right matrix of the decomposition.
  "fit": W \%*\% H. "trimmed1": A list of locations of
  trimmed cells after stage one in each iteration. (to be
  updated) "trimmed2": A list of locations of trimmed cells
  after stage two in each iteration. (to be updated)
  "niter": Number of iterations performed.
}
\description{
  rnmf

  Performs robust penalized non-negative matrix
  factorizaion on a non-negative matrix A to obtain W and
  H, such that A ~= W %*% H, where A is a p by n matrix, W
  is a p by k matrix and H is a k by n matrix. Outliers in
  A are detected and trimmed (by cells, columns or rows,
  see details). Here we assume each row of a represent a
  feature/variable, and each column of A is an
  observation/sample. The objective function is ||A -
  W\%*\%H||_{trimmed} + alpha * ||W||_2^2 + beta *
  sum|H[:,j]|^2. The minimization process is a generalized
  alternating least square method.
}
\examples{
## Load a clean single simulated tumor image.
data("Tumor")
## Add 5\% corruptions.
Tumor.corrupted = Tumor
Tumor.corrupted[sample(1:4900, round(0.05 * 4900), replace = FALSE)] = 1
## Do rnmf with different settings
res.rnmf1 = rnmf(Tumor.corrupted, trim = FALSE)
res.rnmf2 = rnmf(Tumor.corrupted, tol = 0.001, trim = 0.06)
res.rnmf3 = rnmf(Tumor.corrupted, k = 10, beta = 0.1, tol = 0.001, trim = 0.06, my.seed = 123, variation = "smooth")
par(mfrow = c(2,2))
image(Tumor.corrupted, main = "Corrupted")
image(res.rnmf1$fit, main = "rnmf (no trimming) fit")
image(res.rnmf2$fit, main = "rnmf (smooth) fit 2")
image(res.rnmf3$fit, main = "rnmf (smooth) fit 3")
}
\keyword{Non-negative}
\keyword{Robustness,}
\keyword{Sparsity,}
\keyword{Trimming}
\keyword{factorization,}
\keyword{matrix}

