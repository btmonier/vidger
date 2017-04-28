# ggviseq
Make rapid visualizations of RNAseq data in R

Installation
------------

The easiest way to obtain this package is to install `devtools` and pull the package contents directly from GitHub.

``` r
# Development version from GitHub
# install.packages('devtools')

devtools::install_github('btmonier/ggviseq')
```

New Functions
-------------
* `vsBoxplot()`
* `vsDEGmat()`
* `vsFourWay()`
* `vsMAPlot()`
* `vsScatterMatrix()`
* `vsScatterPlot()`
* `vsVolcano()`
* `vsVolcanoMatrix()`

Loading test data
-----------------
I have loaded three objects from three seperate analysis types:

* `df.cuff` A `cuffdiff` output file.
* `df.deseq` A `DESeq2` object class.
* `df.edger` An `edgeR` object class.

Getting help
------------
For additional information on these functions, please see the given documentation in the `ggvseq` package by adding the `?` help operator before any of the given functions in this package or by using the `help()` function. 
