## Setup
## install.packages(c("devtools", "roxygen2", "knitr"))

## install.packages("roxygen2_4.0.0.tgz", repos = NULL)
## Load the libraries
library("devtools")
library("roxygen2")
library("knitr")
## Create the package directory
## create("rnmf")

document("rnmf")
tools::showNonASCII(readLines("./rnmf/man/rnmf.Rd")) 
## Build the package (a tar ball)
system("R CMD build rnmf")
system("R CMD check rnmf")

remove.packages("rnmf")
install.packages("rnmf_0.5.tar.gz", repos = NULL, type = "source")
library(rnmf)
?rnmf
example(rnmf)

