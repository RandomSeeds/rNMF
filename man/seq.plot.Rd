\name{seq.plot}
\alias{seq.plot}
\title{Generates files containing sequential plotting of results of rNMF.}
\usage{
  seq.plot(res, width = 1280, height = 720,
    mylayout = c(2, 5), force = FALSE)
}
\arguments{
  \item{res}{A list, result from the 'rnmf()' function.}

  \item{width}{Width of the page(s).}

  \item{height}{Height of the page(s).}

  \item{mylayout}{Layouts of the graphs on pages. For
  example mylayout = c(4,5) would put 20 graphs on each
  page, in 4 by 5 matrices}

  \item{force}{logic. Force to produce all pages.}
}
\value{
  Nothing. However files will be created.
}
\description{
  The files contain images of rNMF reconstructions with
  subsets of basis vectors in W (and corresponding vectors
  in H). Let k be the number of columns of W. 'seq.plot'
  will generates 2k graphs in one or more png files. The
  first k graphs are W[,u] %*% H[u,], u = 1,..,k. The next
  k graphs are W[,(u1,u2)] %*% H[(u1,u2),], (u1,u2) =
  (1,2),...,(k-1,k), (k,1). The basis vectors in W are
  ordered by the range of values in them.
}
\examples{
## Load a clean single simulated tumor image.
data("Tumor")
## Add 5\% corruptions.
Tumor.corrupted = Tumor
set.seed(1)
Tumor.corrupted[sample(1:4900, round(0.05 * 4900), replace = FALSE)] = 1
## Do rnmf with different settings
res = rnmf(A = Tumor.corrupted, trim = 0.05, variation = "cell", my.seed = 1)
seq.plot(res, mylayout = c(2,5))
}

