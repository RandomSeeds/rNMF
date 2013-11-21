\name{smooth2}
\alias{smooth2}
\title{Smooth a cell with the average of nearby non-outliers.}
\usage{
  smooth2(i, image.f, p, n)
}
\arguments{
  \item{i}{The row number.}

  \item{image.f}{a data.frame.}

  \item{p}{The number of rows of the image matrix.}

  \item{n}{The number of columns of the image matrix.}
}
\value{
  The value that is filled into the missing cell.
}
\description{
  "smooth2" is used in "smoothing2". It takes input i, the
  index of the row in image.f that needs to be smoothed. It
  returns the mean of the closest non-outlier pixels.
  Searching begins with a disk with radius one. The radius
  increases if the search fails. image.f is a data frame
  that contains x coordinates, y coordinates, pixel values,
  logic flags "outs". NOTE: The function needs to be
  improved by changing the search region from a disk to a
  square.
}
\examples{
## Load a clean single simulated tumor image.
}

