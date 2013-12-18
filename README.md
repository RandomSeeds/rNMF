rnmf
====

robust NMF (under construction)

To install, and load the package, do the following.


0. In R console, install dependent packages (the dependent packages cannot be automatically installed because “rnmf” is local):

> install.packages(c("nnls", "plyr"))

1. In R console, run
> install.packages("/Path to the package/rnmf_0.1.0.tar.gz", repos = NULL, type = "source")
change "Path to the package" to the path to the “rnmf_0.1.0.tar.gz” file.

1.5. RESTART R. On my machine, the help files including examples do not show correctly before R is restarted.

2. Run
> library("rnmf")

3. To see the help of the main function, do
> ?rnmf

4. To see an example, do
> example(rnmf)
