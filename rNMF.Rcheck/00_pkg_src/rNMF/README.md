rNMF
====

This is a development repo for the robust nonnegative matrix factorization (rNMF) R-package. rNMF decomposes a high dimensional nonnegative data set with potential corruptions to a product of two low rank matrices with a separate outlier set. 
Package vignette: http://cran.r-project.org/web/packages/rNMF/vignettes/rNMF.html

### Installation
To install and load the published package from CRAN, run the following in R:

```
install.packages("rNMF")
library(rNMF)
```

To install the develop version, run the following in R:

```
install.packages("devtools")
library(devtools)
install_github("RandomSeeds/rNMF")
library(rNMF)
```

### Use the package
The vignette of the package is at

http://cran.r-project.org/web/packages/rNMF/vignettes/rNMF.html

To see the help of the main function, run:

```
?rnmf
```
To see an example, do:

```
example(rnmf)
```
