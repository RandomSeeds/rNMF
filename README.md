rNMF
====

This is a development repo for the robust nonnegative matrix factorization(rNMF) R-package.

Package Vignette: http://hal.case.edu/~yifan/rNMF.html

To install and load the package, do the following.

In R console, install dependent packages:

```
install.packages(c("nnls", "knitr"))
```
Download "rnmf_0.5.0.tar.gz". In R console, run the following line to install the package locally:

```
install.packages("/Path to the package/rnmf_0.5.0.tar.gz", repos = NULL, type = "source")
```
RESTART R. 

Load the package:

```
library("rnmf")
```

To see the help of the main function, run

```
?rnmf
```
To see an example, do

```
example(rnmf)
```
