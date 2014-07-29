# rnmf

## robust NMF

### For users
To install, and load the package, do the following.

Download 'rnmf_0.5tar.gz' (that's all you need) to your local machine.

In R console, run:

```
install.packages("/Path to the package/rnmf_0.5tar.gz", repos = NULL, type = "source")
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

### For developers
  * Main functions are in '/R'. Data is in '/data'. 
  * 'DESCRIPTION' and '/R/rnmf-package.r' contain descriptions for the package.
  * '/R/rnmf.R' and '/R/seq_plot.R' are the two main functions of the package. 
  * Use 'build.R' to build the package.