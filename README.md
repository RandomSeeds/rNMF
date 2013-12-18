rnmf
====

robust NMF (under construction)

To install, and load the package, do the following.

In R console, install dependent packages (the dependent packages cannot be automatically installed because “rnmf” is local):

```
install.packages(c("nnls", "plyr"))
```
In R console, run the following:

```
install.packages("/Path to the package/rnmf_0.1.0.1tar.gz", repos = NULL, type = "source")
```
RESTART R. On my machine, the help files including examples do not show correctly before R is restarted.

Run

```
library("rnmf")
```
To see the help of the main function, do

```
?rnmf
```
To see an example, do

```
example(rnmf)
```
