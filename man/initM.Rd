\name{initM}
\alias{initM}
\title{Randomly initialize a matrix.}
\usage{
  initM(large, nrow, ncol, small = 0, my.seed = NULL)
}
\arguments{
  \item{large}{The upper limit of the uniform
  distribution.}

  \item{nrow}{The number of rows of the matrix to be
  generated.}

  \item{ncol}{The number of columns of the matrix to be
  generated.}

  \item{small}{The lower limit of the uniform
  distribution.}

  \item{my.seed}{The seed to be used for the initialization
  of the matrix.}
}
\value{
  A nrow by ncol matrix.
}
\description{
  Randomly initialize a matrix with uniform distribution.
}
\examples{
initM(1,10,10)
}

